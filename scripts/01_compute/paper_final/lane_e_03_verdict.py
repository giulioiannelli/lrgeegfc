#!/usr/bin/env python3
r"""Lane E / E1+E2b+E4 -- gates, the role-swap separability test, and the scale axis.

Consumes the knob-integrated grid of ``lane_e_01_encinf_knob_grid.py`` and emits
every number the lane report quotes. Nothing here loads FC or draws a surrogate;
the substrate, the null and the functionals are all upstream.

WHAT IS COMPUTED

1. **Knob integration.** Per (patient, band, scale, functional) the margin is
   formed at each plateau fraction against that fraction's own surrogate median,
   then integrated as the MEDIAN over the four fractions, with the across-
   fraction spread reported. A single-fraction number is never reported alone.
   Because one shuffled dense matrix feeds all four fractions, the fraction
   median of realization ``r`` is itself a well-defined knob-integrated
   surrogate realization, so the per-scale gate is run on those.

2. **The primary gate**: one sign-flip cluster-mass test over the whole 16-point
   scale axis per (band, functional) -- the axis is worth ~1.5 independent tests,
   so per-scale correction over-charges by roughly tenfold. The declared primary
   family is the SIX BANDS WITHIN ONE FUNCTIONAL; the 24-cell family (4
   functionals x 6 bands) is the conservative secondary. Both are reported.

3. **The role-swap test (C2, separability).** ``task_learn`` and ``task_test``
   exchanged on the observed graphs. The construction, the split-half baseline
   arms, the durations and the graphs are identical; only the semantic
   assignment moves. Under the null that the two task blocks are exchangeable,
   ``d = T(real) - T(swapped)`` is symmetric about zero per patient, so the
   cohort test is an exact sign-flip test needing no surrogate at all. Two
   algebraic identities make this readable and are verified numerically here:
   ``T_infspec(swap) = -T_infspec`` exactly, and ``T_learn(swap) = T_test``
   exactly, so the swap test on ``T_learn`` IS the paired contrast between
   encoding persistence and whole-task persistence.

4. **Scale.** ``n_eff`` of the scale axis per functional, the per-patient
   margin-profile centroid in ``log s`` paired within patient, and ``N_eff(s)``
   / ``m(s)`` so no scale is described in words without a number.

Outputs: data/paper_final/lane_e_encinf/verdict/*.csv
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PATIENTS_4PHASE
from lrg_eegfc.utils.fc.heat_multiscale import scale_resolution
from lrg_eegfc.utils.metrics.cohort_gate import (
    axis_cluster_gate,
    cohort_margin_gate,
    effective_tests,
)
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr, boot_ci_mean
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.substrate import CANONICAL, canonical_eig, canonical_scale_grid

ROOT = setup_script_env()
BASE = Path(os.environ.get("LANE_E_OUT",
                           ROOT / "data" / "paper_final" / "lane_e_encinf"))
GRID, OUT = BASE / "grid", BASE / "verdict"
FUNCS = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")
N_PERM = int(os.environ.get("LANE_E_NPERM", "10000"))
SEED = 20260831


# --------------------------------------------------------------------------- #
def load_grid(pats, bands):
    """Per band: per-fraction and knob-integrated observed / swapped / margin arrays.

    ``margin_frac`` is ``(K, nF, nS, nK)``; ``margin`` is its median over
    fractions, ``spread`` its across-fraction range. Because one shuffled dense
    matrix feeds all fractions, ``surr_knob`` -- the fraction median of
    realization ``r`` -- is itself a knob-integrated surrogate realization.
    """
    s, fracs, store = None, None, {}
    for band in bands:
        O, W, MF, OK, SK, keep = [], [], [], [], [], []
        for p in pats:
            f = GRID / "cells" / f"{p}__{band}.npz"
            if not f.exists():
                continue
            d = np.load(f)
            s, fracs = d["s"], d["fracs"]
            obs, swap, surr = d["obs"], d["swap"], d["surr"]   # (nF,nS,nK)/(nF,R,nS,nK)
            O.append(obs)
            W.append(swap)
            MF.append(obs - np.nanmedian(surr, axis=1))        # (nF,nS,nK)
            OK.append(np.nanmedian(obs, axis=0))
            SK.append(np.nanmedian(surr, axis=0))              # (R,nS,nK)
            keep.append(p)
        MF = np.array(MF)
        store[band] = dict(patients=keep, obs=np.array(O), swap=np.array(W),
                           margin_frac=MF, obs_knob=np.array(OK),
                           surr_knob=np.array(SK),
                           margin=np.nanmedian(MF, axis=1) if MF.size else MF,
                           spread=(np.nanmax(MF, axis=1) - np.nanmin(MF, axis=1))
                           if MF.size else MF)
    return s, fracs, store


def _cluster(M, rng, sign=1.0):
    r = axis_cluster_gate(sign * np.asarray(M, float), n_perm=N_PERM, rng=rng)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", default="")
    a = ap.parse_args()
    bands = a.bands.split(",") if a.bands else list(BRAIN_BANDS_NAMES)
    pats = list(PATIENTS_4PHASE)
    OUT.mkdir(parents=True, exist_ok=True)
    s, fr, store = load_grid(pats, bands)
    if s is None:
        raise SystemExit(f"no cells under {GRID/'cells'}")
    nS = s.size

    # ---------------- 0. algebraic identities of the swap, verified ---------- #
    idn = []
    for band in bands:
        st = store[band]
        if st["obs"].size == 0:
            continue
        o, w = st["obs"], st["swap"]
        idn.append(dict(
            band=band,
            max_abs_dev_infspec_antisym=float(np.nanmax(np.abs(w[..., 2] + o[..., 2]))),
            max_abs_dev_learnswap_eq_test=float(np.nanmax(np.abs(w[..., 1] - o[..., 0]))),
            max_abs_dev_testswap_eq_learn=float(np.nanmax(np.abs(w[..., 0] - o[..., 1]))),
        ))
    pd.DataFrame(idn).to_csv(OUT / "swap_identities.csv", index=False)

    # ---------------- 1. primary gate: per-band axis cluster ---------------- #
    rows = []
    rng = np.random.default_rng(SEED)
    for band in bands:
        st = store[band]
        if not st["patients"]:
            continue
        for k, fn in enumerate(FUNCS):
            M = st["margin"][:, :, k]
            r = _cluster(M, rng)
            ne = effective_tests(M)
            mean, lo, hi = boot_ci_mean(np.nanmean(M, axis=1), B=10000,
                                        rng=np.random.default_rng(SEED + k))
            rows.append(dict(
                band=band, func=fn, n_pat=M.shape[0], p_cluster=r["p"],
                mass=r["mass"], cluster_lo=(r["cluster"] or (np.nan, np.nan))[0],
                cluster_hi=(r["cluster"] or (np.nan, np.nan))[1],
                s_lo=s[r["cluster"][0]] if r["cluster"] else np.nan,
                s_hi=s[r["cluster"][1]] if r["cluster"] else np.nan,
                margin_med=float(np.nanmedian(M)),
                margin_mean=mean, ci_lo=lo, ci_hi=hi,
                knob_spread_med=float(np.nanmedian(st["spread"][:, :, k])),
                n_eff_pr=ne.get("n_eff_pr"), n_eff_cn=ne.get("n_eff_cn"),
                mean_r=ne.get("mean_offdiag"),
            ))
    cl = pd.DataFrame(rows)
    for fn in FUNCS:                      # PRIMARY family: 6 bands within functional
        m = cl.func == fn
        cl.loc[m, "q_within_func"] = bh_fdr(cl.loc[m, "p_cluster"].values)
    cl["q_all24"] = bh_fdr(cl["p_cluster"].values)          # conservative secondary
    cl.to_csv(OUT / "cluster_gate.csv", index=False)

    # ---------------- 2. LOO of the cluster verdict ------------------------- #
    loo = []
    for band in bands:
        st = store[band]
        if not st["patients"]:
            continue
        for k, fn in enumerate(FUNCS):
            M = st["margin"][:, :, k]
            for i, p in enumerate(st["patients"]):
                keep = np.ones(M.shape[0], bool); keep[i] = False
                r = _cluster(M[keep], np.random.default_rng(SEED + 7 * i))
                loo.append(dict(band=band, func=fn, dropped=p, p_cluster=r["p"]))
    pd.DataFrame(loo).to_csv(OUT / "cluster_gate_loo.csv", index=False)

    # ---------------- 3. per-scale grid (secondary, reported in full) ------- #
    grid = []
    for band in bands:
        st = store[band]
        if not st["patients"]:
            continue
        for k, fn in enumerate(FUNCS):
            for j in range(nS):
                g = cohort_margin_gate(st["obs_knob"][:, j, k],
                                       st["surr_knob"][:, :, j, k],
                                       labels=st["patients"], full=True,
                                       expect_n=len(pats),
                                       rng=np.random.default_rng(SEED + j))
                grid.append(dict(band=band, func=fn, s=float(s[j]),
                                 **{kk: v for kk, v in g.items() if kk != "loo_p"}))
    gd = pd.DataFrame(grid)
    for fn in FUNCS:
        m = gd.func == fn
        gd.loc[m, "q_grid_within_func"] = bh_fdr(gd.loc[m, "p"].values)
    gd.to_csv(OUT / "per_scale_grid.csv", index=False)

    # ---------------- 4. per-fraction gate (knob robustness) ---------------- #
    frac_rows = []
    for band in bands:
        st = store[band]
        if not st["patients"]:
            continue
        for fi, f in enumerate(fr):
            Mf = st["margin_frac"][:, fi]                         # (K,nS,nK)
            for k, fn in enumerate(FUNCS):
                r = _cluster(Mf[:, :, k], np.random.default_rng(SEED + fi))
                frac_rows.append(dict(band=band, func=fn, frac=float(f),
                                      p_cluster=r["p"],
                                      margin_med=float(np.nanmedian(Mf[:, :, k]))))
    pd.DataFrame(frac_rows).to_csv(OUT / "per_fraction_gate.csv", index=False)

    # ---------------- 5. the role-swap separability test -------------------- #
    swp = []
    for band in bands:
        st = store[band]
        if not st["patients"]:
            continue
        d_all = np.nanmedian(st["obs"] - st["swap"], axis=1)      # knob median (K,nS,nK)
        for k, fn in enumerate(FUNCS):
            d = d_all[:, :, k]
            rp = _cluster(d, np.random.default_rng(SEED + 11), +1.0)
            rn = _cluster(d, np.random.default_rng(SEED + 11), -1.0)
            # per-scale exact sign-flip is equivalent to the signed-rank on d
            ps = [float(wilcoxon(d[:, j], alternative="greater").pvalue)
                  if np.any(d[:, j] != 0) else np.nan for j in range(nS)]
            swp.append(dict(band=band, func=fn, n_pat=d.shape[0],
                            p_cluster_pos=rp["p"], p_cluster_neg=rn["p"],
                            d_med=float(np.nanmedian(d)),
                            n_scales_p05_pos=int(np.nansum(np.array(ps) < 0.05)),
                            best_p_scale=float(np.nanmin(ps)) if np.isfinite(ps).any() else np.nan))
    sw = pd.DataFrame(swp)
    for fn in FUNCS:
        m = sw.func == fn
        sw.loc[m, "q_pos_within_func"] = bh_fdr(sw.loc[m, "p_cluster_pos"].values)
    sw.to_csv(OUT / "swap_test.csv", index=False)

    # ---------------- 5b. duration-ratio regression on the swap statistic --- #
    # task_test is LONGER than task_learn in all 10 patients (ratio 0.40-0.74),
    # so the swap null's exchangeability assumption is violated in a cohort-
    # consistent way. Recordings are never truncated in this project; the
    # sanctioned control is a duration-ratio regression, which asks whether the
    # swap statistic would still be non-zero at a ratio of 1.
    dur_p = ROOT / "data" / "paper_final" / "w0b_nulls" / "durations.json"
    if dur_p.exists():
        durs = json.loads(dur_p.read_text())
        drows = []
        for band in bands:
            st = store[band]
            if not st["patients"]:
                continue
            x = np.array([np.log10(durs[p]["task_learn"] / durs[p]["task_test"])
                          for p in st["patients"]])
            d_all = np.nanmedian(st["obs"] - st["swap"], axis=1)
            for k, fn in enumerate(FUNCS):
                y = np.nanmean(d_all[:, :, k], axis=1)     # scale-averaged, descriptive
                ok = np.isfinite(x) & np.isfinite(y)
                if ok.sum() < 5:
                    continue
                xx, yy = x[ok], y[ok]
                sl, ic = np.polyfit(xx, yy, 1)
                rr = np.corrcoef(xx, yy)[0, 1]
                brng = np.random.default_rng(SEED)
                bi = np.array([np.polyfit(xx[i], yy[i], 1)[1]
                               for i in (brng.integers(0, xx.size, (2000, xx.size)))])
                drows.append(dict(band=band, func=fn, n=int(ok.sum()), slope=float(sl),
                                  intercept_at_ratio1=float(ic),
                                  intercept_ci_lo=float(np.percentile(bi, 2.5)),
                                  intercept_ci_hi=float(np.percentile(bi, 97.5)),
                                  r_dur_vs_d=float(rr),
                                  mean_d=float(np.nanmean(y[ok]))))
        pd.DataFrame(drows).to_csv(OUT / "swap_duration_regression.csv", index=False)

    # ---------------- 6. scale position: centroids, paired ------------------ #
    cen = []
    logs = np.log10(s)
    for band in bands:
        st = store[band]
        if not st["patients"]:
            continue
        C = np.full((len(st["patients"]), len(FUNCS)), np.nan)
        for i in range(len(st["patients"])):
            for k in range(len(FUNCS)):
                w = np.clip(st["margin"][i, :, k], 0, None)
                if np.nansum(w) > 0:
                    C[i, k] = float(np.nansum(w * logs) / np.nansum(w))
        for ka, kb in ((1, 3), (1, 2), (0, 3)):
            ok = np.isfinite(C[:, ka]) & np.isfinite(C[:, kb])
            if ok.sum() < 5:
                continue
            dd = C[ok, kb] - C[ok, ka]
            mean, lo, hi = boot_ci_mean(dd, B=10000, rng=np.random.default_rng(SEED))
            cen.append(dict(band=band, func_a=FUNCS[ka], func_b=FUNCS[kb], n=int(ok.sum()),
                            centroid_s_a=float(10 ** np.nanmedian(C[ok, ka])),
                            centroid_s_b=float(10 ** np.nanmedian(C[ok, kb])),
                            mean_shift_log10=mean, ci_lo=lo, ci_hi=hi,
                            shift_factor=float(10 ** mean),
                            p_two=float(wilcoxon(dd).pvalue) if np.any(dd != 0) else np.nan))
    ce = pd.DataFrame(cen)
    if not ce.empty:
        ce["q_two"] = bh_fdr(ce["p_two"].values)
    ce.to_csv(OUT / "scale_centroids.csv", index=False)

    # ---------------- 6b. per-patient profiles (fluctuations are signal) ----- #
    pp = []
    for band in bands:
        st = store[band]
        for i, p in enumerate(st["patients"]):
            for k, fn in enumerate(FUNCS):
                for j in range(nS):
                    pp.append(dict(band=band, patient=p, func=fn, s=float(s[j]),
                                   margin=float(st["margin"][i, j, k]),
                                   obs=float(st["obs_knob"][i, j, k]),
                                   swap=float(np.nanmedian(st["swap"][i, :, j, k])),
                                   knob_spread=float(st["spread"][i, j, k])))
    pd.DataFrame(pp).to_csv(OUT / "per_patient_profiles.csv", index=False)

    # ---------------- 7. scale units ---------------------------------------- #
    un = []
    for band in bands:
        for p in store[band]["patients"]:
            try:
                ev, V = canonical_eig(p, "A", band)
            except Exception:                                      # noqa: BLE001
                continue
            r = scale_resolution(ev, V, s)
            for j in range(nS):
                un.append(dict(band=band, patient=p, s=float(s[j]),
                               n_eff=float(r["n_eff"][j]), m_comm=float(r["m_comm"][j]),
                               N=int(r["N"])))
    pd.DataFrame(un).to_csv(OUT / "scale_units.csv", index=False)

    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="lane_e_verdict", substrate=CANONICAL.label(),
        fracs=[float(x) for x in fr], n_perm=N_PERM, bands=bands, cohort=pats,
        primary_family="6 bands within one functional (axis-cluster gate)",
        secondary_family="24 cells (4 functionals x 6 bands)",
        s_grid=[float(x) for x in s], seed=SEED), indent=2))

    # ---------------- console summary --------------------------------------- #
    pd.set_option("display.width", 200)
    print("\n=== swap algebraic identities (must be ~0) ===")
    print(pd.DataFrame(idn).to_string(index=False))
    print("\n=== PRIMARY gate: axis-cluster, q within functional ===")
    print(cl[["band", "func", "p_cluster", "q_within_func", "q_all24",
              "margin_med", "knob_spread_med", "s_lo", "s_hi", "n_eff_pr"]]
          .to_string(index=False))
    print("\n=== role-swap separability test ===")
    print(sw.to_string(index=False))
    print("\n=== scale centroids ===")
    print(ce.to_string(index=False) if not ce.empty else "(none)")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
