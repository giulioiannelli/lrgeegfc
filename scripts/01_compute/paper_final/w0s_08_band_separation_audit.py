#!/usr/bin/env python3
"""Is low_gamma's near-clearing a genuine near-miss, or an artefact of the gate?

Part A's grid puts low_gamma one place behind beta on the trace axis and left
that unexamined. Sibling lane W0-A independently reports that adding low_gamma
to the required-null set shrinks the band-selectivity window to 0.51 octaves,
and specifically that low_gamma *leaks* at the low end of the plateau
(f = 0.07, 0.10) and is null only for f >= 0.14. Two lanes therefore point at
the same weakness in the project's second result, and nobody owns it. This
script owns it.

The call was pre-registered in scope report section 12.5, before any of these
numbers were computed:

* **genuine near-miss** -- the bands are simply not separable at n = 10 -- if
  low_gamma passes calibration, its axis cluster survives leave-one-patient-out,
  its margin is not carried by a minority of patients, AND the *paired*
  per-patient difference ``beta - low_gamma`` fails to clear;
* **gate artefact** if it fails calibration, or its clearing collapses under
  LOO, or the paired difference clears.

Two tests do the work, and neither has been run on these bands before.

1. **Paired band difference.** Reading two marginal q-values side by side is
   not a test of whether the bands differ. The difference is tested directly:
   per patient, per scale, ``margin_beta - margin_low_gamma``, gated with the
   locked sign-flip cluster test in both directions. A cohort that cannot
   separate the bands and a cohort in which they genuinely differ produce the
   same pair of marginal q-values, and only the paired test tells them apart.
2. **Knob-resolved decomposition.** The headline is a median over the plateau
   f in {0.07, 0.10, 0.14, 0.20}. If low_gamma's cluster is carried by the low-f
   end and gone at high f, the integrated number is averaging a knob artefact
   rather than four equally good readings -- and that would reproduce W0-A's
   finding from a different statistic on the same substrate.

Every number is recomputed from the cells on disk. Nothing is read back from a
prior CSV or a prior report.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.metrics.cohort_gate import (
    DESCRIPTIVE_ALPHA,
    axis_cluster_gate,
    axis_reversal_gate,
    calibrate_from_surrogates,
    effective_tests,
)
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from w0s_00_artifacts import load_cells                       # noqa: E402

OUT = Path(os.environ.get(
    "W0S_OUT_ANALYSIS", ROOT / "data" / "paper_final" / "lane_s_scale"))
READOUT = "T"
N_PERM = int(os.environ.get("W0S_NPERM", "10000"))
N_CAL = 200
PAIR = ("beta", "low_gamma")


#: A cell is "calibrated" when its measured false-positive rate at nominal 0.05
#: does not exceed this -- the same tolerance ``gate_grid`` applies.
FPR_TOL = 0.10


def _cluster(M, seed=7):
    return axis_cluster_gate(M, n_perm=N_PERM, rng=np.random.default_rng(seed))


def _calibration_over_axis(S, n_draws=N_CAL, seed=5):
    """Held-out-realization FPR at every scale of the axis, summarised.

    :func:`calibrate_from_surrogates` calibrates ONE cell, i.e. one scale. A
    band's axis is 16 cells, and a band that is well calibrated on average but
    badly calibrated where its cluster sits is not calibrated for the claim
    being made. The worst scale is therefore reported alongside the median.
    """
    fpr = np.array([calibrate_from_surrogates(
        S[:, :, j], n_draws=n_draws,
        rng=np.random.default_rng(seed + j))["fpr_0.05"]
        for j in range(S.shape[2])], float)
    ok = np.isfinite(fpr)
    return dict(fpr_05_med=float(np.nanmedian(fpr)),
                fpr_05_max=float(np.nanmax(fpr)) if ok.any() else np.nan,
                n_scales_calibrated=int(np.sum(fpr[ok] <= FPR_TOL)),
                n_scales_checked=int(ok.sum()),
                calibrated=bool(ok.any() and np.nanmax(fpr) <= FPR_TOL))


def main() -> None:
    cs = load_cells()
    OUT.mkdir(parents=True, exist_ok=True)
    s = np.asarray(cs.s, float)
    print(f"[w0s-bandsep] readout {READOUT} | {len(cs.bands)} bands x "
          f"{cs.n_scales} scales | R={cs.R} | fracs={list(cs.fracs)}", flush=True)

    # -- 1. knob-integrated per-band picture, recomputed ------------------- #
    rows = []
    M = {}
    for b in cs.bands:
        M[b], labs = cs.margins(b, READOUT)
        ac = _cluster(M[b])
        O, S, _ = cs.cell(b, READOUT)
        cal = _calibration_over_axis(S)
        # how concentrated is the cohort margin: patients positive at the
        # cluster's own scales, and the single largest patient's share of the
        # cohort mean. A margin carried by a minority is not a cohort effect.
        sel = (slice(None), slice(ac["cluster"][0], ac["cluster"][1] + 1)
               if ac["cluster"] is not None else slice(None))
        Mc = M[b][sel]
        per_pat = np.nanmean(Mc, axis=1)
        pos = int((per_pat > 0).sum())
        top_share = (float(np.nanmax(per_pat) / np.nansum(per_pat))
                     if np.nansum(per_pat) > 0 else np.nan)
        loo_p = [_cluster(np.delete(M[b], k, axis=0), seed=100 + k)["p"]
                 for k in range(M[b].shape[0])]
        rows.append(dict(
            band=b, n_pat=M[b].shape[0],
            margin_med=float(np.nanmedian(M[b])),
            margin_mean=float(np.nanmean(M[b])),
            cluster_p=ac["p"], cluster_mass=ac["mass"],
            cluster_lo=(None if ac["cluster"] is None else float(s[ac["cluster"][0]])),
            cluster_hi=(None if ac["cluster"] is None else float(s[ac["cluster"][1]])),
            cluster_width_pts=(0 if ac["cluster"] is None
                               else ac["cluster"][1] - ac["cluster"][0] + 1),
            **cal,
            n_pat_positive_in_cluster=pos,
            top_patient_share=top_share,
            loo_p_max=float(np.max(loo_p)), loo_p_min=float(np.min(loo_p)),
            loo_worst_patient=labs[int(np.argmax(loo_p))],
            loo_n_flips=int(sum(p >= DESCRIPTIVE_ALPHA for p in loo_p)),
            n_eff_pr=effective_tests(M[b])["n_eff_pr"]))
        print(f"  {b:11s} cluster p={ac['p']:.4f} margin={rows[-1]['margin_med']:+.4f} "
              f"loo_p_max={rows[-1]['loo_p_max']:.4f} "
              f"calibrated={rows[-1]['calibrated']}", flush=True)
    band = pd.DataFrame(rows)
    band["cluster_q"] = np.asarray(bh_fdr(band.cluster_p.to_numpy()), float)
    band.to_csv(OUT / "band_separation_bands.csv", index=False)

    # -- 2. paired band difference ---------------------------------------- #
    a, bnd = PAIR
    D = M[a] - M[bnd]
    rv = axis_reversal_gate(D, s, n_perm=N_PERM, rng=np.random.default_rng(21))
    per_scale = []
    for j in range(cs.n_scales):
        d = D[:, j]
        d = d[np.isfinite(d) & (d != 0)]
        pg = float(wilcoxon(d, alternative="greater")[1]) if d.size >= 5 else np.nan
        pl = float(wilcoxon(d, alternative="less")[1]) if d.size >= 5 else np.nan
        per_scale.append(dict(s=float(s[j]), d_med=float(np.nanmedian(D[:, j])),
                              n_pos=int((D[:, j] > 0).sum()),
                              p_beta_gt=pg, p_beta_lt=pl))
    ps = pd.DataFrame(per_scale)
    ps["q_beta_gt"] = np.asarray(bh_fdr(ps.p_beta_gt.to_numpy()), float)
    ps["q_beta_lt"] = np.asarray(bh_fdr(ps.p_beta_lt.to_numpy()), float)
    ps.insert(0, "pair", f"{a}_minus_{bnd}")
    ps.to_csv(OUT / "band_separation_paired_per_scale.csv", index=False)
    print(f"  paired {a}-{bnd}: cluster p(beta>low_g)={rv['p_pos']:.4f} "
          f"p(low_g>beta)={rv['p_neg']:.4f} | per-scale min q "
          f"{np.nanmin(ps.q_beta_gt):.3f}/{np.nanmin(ps.q_beta_lt):.3f}", flush=True)

    # -- 3. knob-resolved decomposition ----------------------------------- #
    krows = []
    for fi, f in enumerate(cs.fracs):
        for b in cs.bands:
            O, S, _ = cs.cell_at(b, READOUT, fi)
            # patient_margin, applied scale by scale: O is (K, nS), S is
            # (K, R, nS), and the locked margin is obs minus the patient's own
            # surrogate median at that same scale.
            Mf = O - np.nanmedian(S, axis=1)
            ac = _cluster(Mf, seed=200 + fi)
            krows.append(dict(frac=float(f), band=b, cluster_p=ac["p"],
                              cluster_mass=ac["mass"],
                              margin_med=float(np.nanmedian(Mf)),
                              cluster_width_pts=(0 if ac["cluster"] is None else
                                                 ac["cluster"][1] - ac["cluster"][0] + 1)))
    kn = pd.DataFrame(krows)
    kn["cluster_q"] = np.nan
    for f in kn.frac.unique():
        sel = kn.frac == f
        kn.loc[sel, "cluster_q"] = np.asarray(
            bh_fdr(kn.loc[sel, "cluster_p"].to_numpy()), float)
    kn.to_csv(OUT / "band_separation_by_fraction.csv", index=False)
    piv = kn.pivot(index="band", columns="frac", values="cluster_q").round(4)
    print("  cluster q by plateau fraction (BH within fraction, 6 bands):",
          flush=True)
    print(piv.to_string(), flush=True)

    # -- 3b. do the two bands separate ALONG THE KNOB? --------------------- #
    # Exploratory, and labelled as such: it composes the two pre-registered
    # tests rather than adding a third object. The knob decomposition can show
    # two bands responding to the sparsification fraction in opposite
    # directions while the knob-integrated paired test sees nothing, because
    # the median over fractions is exactly what averages that away. One test,
    # not a search: each patient's margin is regressed on log f (averaged over
    # scales), and the two bands' slopes are compared pairwise.
    def _knob_slope(bd):
        sl = np.full(len(cs.have(bd)), np.nan)
        prof = np.stack([  # (nF, K) mean margin over the scale axis
            np.nanmean(cs.cell_at(bd, READOUT, fi)[0]
                       - np.nanmedian(cs.cell_at(bd, READOUT, fi)[1], axis=1),
                       axis=1)
            for fi in range(len(cs.fracs))])
        lf = np.log(np.asarray(cs.fracs, float))
        for k in range(prof.shape[1]):
            y = prof[:, k]
            if np.isfinite(y).all():
                sl[k] = np.polyfit(lf, y, 1)[0]
        return sl, prof

    sl_a, prof_a = _knob_slope(a)
    sl_b, prof_b = _knob_slope(bnd)
    dks = sl_a - sl_b
    dks = dks[np.isfinite(dks)]
    p_knob = (float(wilcoxon(dks, alternative="greater")[1])
              if dks.size >= 5 else np.nan)
    pd.DataFrame(dict(
        patient=cs.have(a), knob_slope_beta=sl_a, knob_slope_low_gamma=sl_b,
        diff=sl_a - sl_b)).to_csv(OUT / "band_separation_knob_slopes.csv",
                                  index=False)
    kf = []
    for fi, f in enumerate(cs.fracs):
        Da = (cs.cell_at(a, READOUT, fi)[0]
              - np.nanmedian(cs.cell_at(a, READOUT, fi)[1], axis=1))
        Db = (cs.cell_at(bnd, READOUT, fi)[0]
              - np.nanmedian(cs.cell_at(bnd, READOUT, fi)[1], axis=1))
        acf = _cluster(Da - Db, seed=300 + fi)
        kf.append(dict(frac=float(f), p_beta_gt_lowg=acf["p"],
                       mass=acf["mass"],
                       d_margin_med=float(np.nanmedian(Da - Db))))
    kfd = pd.DataFrame(kf)
    kfd["q"] = np.asarray(bh_fdr(kfd.p_beta_gt_lowg.to_numpy()), float)
    kfd.to_csv(OUT / "band_separation_paired_by_fraction.csv", index=False)
    print(f"  knob slope of margin (per patient, vs log f): "
          f"beta {np.nanmean(sl_a):+.4f} low_gamma {np.nanmean(sl_b):+.4f} "
          f"| paired p(beta rises more) = {p_knob:.4f}", flush=True)
    print("  paired beta-low_gamma cluster gate, per fraction:", flush=True)
    print(kfd.round(4).to_string(index=False), flush=True)

    # -- 4. the pre-registered call --------------------------------------- #
    lg = band[band.band == bnd].iloc[0]
    be = band[band.band == a].iloc[0]
    paired_clears = bool(min(rv["p_pos"], rv["p_neg"]) < DESCRIPTIVE_ALPHA)
    loo_collapses = bool(lg.loo_p_max >= DESCRIPTIVE_ALPHA)
    uncalibrated = not bool(lg.calibrated)
    minority = bool(lg.n_pat_positive_in_cluster < 6)
    verdict = ("gate_artefact" if (paired_clears or loo_collapses or uncalibrated)
               else "genuine_near_miss")
    lo_f = kn[(kn.band == bnd) & (kn.frac <= 0.10)].cluster_p.to_numpy()
    hi_f = kn[(kn.band == bnd) & (kn.frac >= 0.14)].cluster_p.to_numpy()
    call = dict(
        pair=f"{a}_vs_{bnd}", verdict=verdict,
        low_gamma_cluster_p=float(lg.cluster_p), low_gamma_cluster_q=float(lg.cluster_q),
        beta_cluster_p=float(be.cluster_p), beta_cluster_q=float(be.cluster_q),
        low_gamma_margin_med=float(lg.margin_med), beta_margin_med=float(be.margin_med),
        paired_p_beta_gt=float(rv["p_pos"]), paired_p_lowg_gt=float(rv["p_neg"]),
        paired_clears=paired_clears,
        low_gamma_calibrated=bool(lg.calibrated),
        low_gamma_fpr_05_max=float(lg.fpr_05_max),
        low_gamma_fpr_05_med=float(lg.fpr_05_med),
        beta_calibrated=bool(be.calibrated), beta_fpr_05_max=float(be.fpr_05_max),
        low_gamma_loo_p_max=float(lg.loo_p_max), loo_collapses=loo_collapses,
        low_gamma_n_pat_positive=int(lg.n_pat_positive_in_cluster),
        minority_carried=minority,
        low_gamma_p_lowf_max=float(np.max(lo_f)) if lo_f.size else np.nan,
        low_gamma_p_highf_min=float(np.min(hi_f)) if hi_f.size else np.nan,
        knob_dependent=bool(lo_f.size and hi_f.size
                            and np.max(lo_f) < DESCRIPTIVE_ALPHA
                            <= np.min(hi_f)),
        knob_slope_beta=float(np.nanmean(sl_a)),
        knob_slope_low_gamma=float(np.nanmean(sl_b)),
        knob_slope_paired_p=p_knob,
        paired_by_fraction_min_q=float(np.nanmin(kfd["q"])))
    pd.DataFrame([call]).to_csv(OUT / "band_separation_verdict.csv", index=False)
    for k, v in call.items():
        print(f"    {k:26s} {v}", flush=True)
    print(f"[w0s-bandsep] wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
