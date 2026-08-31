#!/usr/bin/env python3
"""S0 criteria, S3 fine grid, S4 verdict -- does the trace vary across scales?

Everything here runs through the locked gate
(:mod:`lrg_eegfc.utils.metrics.cohort_gate`) unmodified: margins never raw
values, per-scale never best-scale, both multiplicity families side by side,
leave-one-patient-out of the verdict, and held-out-realization calibration with
uncalibrated cells withheld rather than caveated.

Stages
------
1. **Gate** every readout on the (band x scale) grid, whole-grid BH within a
   readout, plus the per-band sign-flip axis-cluster test (BH over 6 bands) --
   the family the locked contract uses, since a 32-point sweep is not 32 tests.
2. **Cross-scale independence**, observed and null-referenced. ``effective_tests``
   on the per-patient margin matrix gives the observed ``n_eff``; the same
   statistic on held-out surrogate realizations of the *same readout* gives the
   noise floor. A readout that decorrelates its scales only as much as its own
   noise does has bought nothing, which is the trap the pre-registration exists
   to catch.
3. **Calibration** per cell, `refuse_uncalibrated=True`.
4. **The three pre-registered criteria** (scope report §5.2), evaluated verbatim
   and reported for every construction attempted, including the failures.
5. **S3**: the scale axis re-read on its two halves -- the near-local regime
   ``s < 3`` (18 of 32 grid points; the locked contract grid starts at s = 1 and
   never samples below it) against the spatially near-global regime ``s >= 3``.

Nothing is tested against zero anywhere.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.metrics.cohort_gate import (
    axis_cluster_gate,
    calibrate_from_surrogates,
    effective_tests,
    gate_grid,
)
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from w0s_00_artifacts import load_cells                      # noqa: E402

OUT = Path(os.environ.get(
    "W0S_OUT_ANALYSIS", ROOT / "data" / "paper_final" / "lane_s_scale"))

#: The readouts carried through the full contract (bootstrap CI, LOO,
#: calibration). The stratum-resolved surface is gated too, but exploratorily.
HEADLINE = ("T", "Ccon_diag", "Tloc_diag", "Qcon_diag", "Qloc_diag", "Thei")
N_PERM = 10_000
N_CAL_DRAWS = 100
N_NEFF_DRAWS = 100
EXPECT_N = 10
FINE_MAX = 3.0


def _shape_stats(M: np.ndarray, s: np.ndarray) -> dict:
    """Two direct readings of "does the effect vary along the scale axis?".

    ``slope_abs``
        cohort mean of each patient's own ``|Spearman(margin, log s)|`` -- how
        monotonically scale-dependent an individual profile is. Defined for
        every patient, unlike a coefficient of variation, which needs a positive
        mean and silently drops patients.
    ``shape_agreement``
        mean off-diagonal correlation between patients' **z-scored** profiles --
        do patients agree on *where* along the axis the effect is larger? This
        is what separates the two readings a flat cohort curve admits: a genuine
        within-patient invariance, and patients each having a differently-shaped
        profile that averages flat. A Friedman cannot tell those apart, which is
        why it is not used here.

    Both are referenced to held-out surrogate draws by the caller. Neither is
    ever compared to zero.
    """
    ok = np.isfinite(M).all(axis=0)
    X, ls = M[:, ok], np.log(np.asarray(s)[ok])
    out = dict(slope_abs=np.nan, slope_signed=np.nan, shape_agreement=np.nan,
               n_finite=int(ok.sum()))
    if X.shape[1] < 4 or X.shape[0] < 3:
        return out
    rho = np.array([spearmanr(X[k], ls)[0] for k in range(X.shape[0])])
    rho = rho[np.isfinite(rho)]
    if rho.size:
        out["slope_abs"] = float(np.mean(np.abs(rho)))
        out["slope_signed"] = float(np.mean(rho))
    sd = X.std(axis=1, ddof=1)
    if np.all(sd > 0):
        Z = (X - X.mean(axis=1, keepdims=True)) / sd[:, None]
        C = np.corrcoef(Z)
        iu = np.triu_indices(C.shape[0], 1)
        out["shape_agreement"] = float(np.nanmean(C[iu]))
    return out


def _null_draws(cs, band: str, readout: str, s: np.ndarray, n_draws: int) -> dict:
    """``n_eff_pr``, ``slope_abs`` and ``shape_agreement`` on held-out surrogate margins.

    Under the null the observation is exchangeable with its own surrogates, so
    these are draws from each statistic's null distribution *at this readout's
    own noise level*. That last clause is the point: a noisier statistic
    decorrelates its scales and inflates ``n_eff`` for free, so only this
    comparison can tell extra information from extra noise.
    """
    _, S, _ = cs.cell(band, readout)
    n = min(n_draws, S.shape[1])
    out = {k: np.full(n, np.nan) for k in ("n_eff_pr", "slope_abs",
                                           "shape_agreement")}
    for i in range(n):
        Mn = cs.heldout_margins(band, readout, i)
        out["n_eff_pr"][i] = effective_tests(Mn)["n_eff_pr"]
        sh = _shape_stats(Mn, s)
        out["slope_abs"][i] = sh["slope_abs"]
        out["shape_agreement"][i] = sh["shape_agreement"]
    return out


def _upper_p(obs: float, null: np.ndarray) -> float:
    n = null[np.isfinite(null)]
    if not n.size or not np.isfinite(obs):
        return np.nan
    return float((1 + int((n >= obs).sum())) / (n.size + 1))


def main() -> None:
    cs = load_cells()
    OUT.mkdir(parents=True, exist_ok=True)
    s = cs.s
    fine = s < FINE_MAX
    print(f"[w0s-verdict] {len(cs.readouts)} readouts x {len(cs.bands)} bands "
          f"x {cs.n_scales} scales | fine sub-grid s<{FINE_MAX}: "
          f"{int(fine.sum())} points", flush=True)

    # ------------------------------------------------------------------ #
    # 1. per-cell gate, all readouts
    # ------------------------------------------------------------------ #
    grids = []
    for m in cs.readouts:
        cells, labs0 = [], None
        for b in cs.bands:
            O, S, labs = cs.cell(b, m)
            labs0 = labs
            for j in range(cs.n_scales):
                cells.append(({"readout": m, "band": b, "s": float(s[j])},
                              O[:, j], S[:, :, j]))
        full = m in HEADLINE
        df = gate_grid(cells, labels=labs0, full=full, n_boot=2000,
                       expect_n=EXPECT_N,
                       calibrate=full, calibration_draws=N_CAL_DRAWS,
                       refuse_uncalibrated=True)
        grids.append(df)
        print(f"  gated {m:12s} cleared(q<.05) {int((df.q < .05).sum()):4d}"
              f"/{len(df)}", flush=True)
    gate = pd.concat(grids, ignore_index=True)
    gate.to_csv(OUT / "gate_grid.csv", index=False)

    # ------------------------------------------------------------------ #
    # 2. axis-cluster gate + cross-scale independence (observed and null)
    # ------------------------------------------------------------------ #
    rows = []
    for m in cs.readouts:
        for b in cs.bands:
            M, _ = cs.margins(b, m)
            e = effective_tests(M)
            ac = axis_cluster_gate(M, n_perm=N_PERM,
                                   rng=np.random.default_rng(7))
            ac_f = axis_cluster_gate(M[:, fine], n_perm=N_PERM,
                                     rng=np.random.default_rng(7))
            ac_c = axis_cluster_gate(M[:, ~fine], n_perm=N_PERM,
                                     rng=np.random.default_rng(7))
            nd = _null_draws(cs, b, m, s, N_NEFF_DRAWS)
            nn = nd["n_eff_pr"][np.isfinite(nd["n_eff_pr"])]
            sh = _shape_stats(M, s)
            e_f = effective_tests(M[:, fine])
            rows.append(dict(
                readout=m, band=b,
                n_eff_pr=e["n_eff_pr"], n_eff_cn=e["n_eff_cn"],
                mean_offdiag=e["mean_offdiag"], n_finite=e["n_finite"],
                n_eff_pr_fine=e_f["n_eff_pr"], mean_offdiag_fine=e_f["mean_offdiag"],
                n_eff_null_med=float(np.median(nn)) if nn.size else np.nan,
                n_eff_null_p95=float(np.percentile(nn, 95)) if nn.size else np.nan,
                p_neff_above_null=_upper_p(e["n_eff_pr"], nd["n_eff_pr"]),
                slope_abs=sh["slope_abs"], slope_signed=sh["slope_signed"],
                slope_abs_null_med=float(np.nanmedian(nd["slope_abs"])),
                p_slope_above_null=_upper_p(sh["slope_abs"], nd["slope_abs"]),
                shape_agreement=sh["shape_agreement"],
                shape_agreement_null_med=float(
                    np.nanmedian(nd["shape_agreement"])),
                p_shape_above_null=_upper_p(sh["shape_agreement"],
                                            nd["shape_agreement"]),
                cluster_p=ac["p"], cluster_mass=ac["mass"],
                cluster_lo=(None if ac["cluster"] is None
                            else float(s[ac["cluster"][0]])),
                cluster_hi=(None if ac["cluster"] is None
                            else float(s[ac["cluster"][1]])),
                cluster_p_fine=ac_f["p"], cluster_p_coarse=ac_c["p"],
                margin_med=float(np.nanmedian(M)),
                knob_iqr_med=float(np.nanmedian(cs.knob_spread(b, m))),
            ))
        print(f"  axis {m:12s} done", flush=True)
    ax = pd.DataFrame(rows)
    for col, out in (("cluster_p", "cluster_q"), ("cluster_p_fine", "cluster_q_fine"),
                     ("cluster_p_coarse", "cluster_q_coarse")):
        ax[out] = np.nan
        for m in cs.readouts:
            sel = ax.readout == m
            ok = sel & ax[col].notna()
            if ok.any():
                ax.loc[ok, out] = bh_fdr(ax.loc[ok, col].to_numpy())
    ax.to_csv(OUT / "axis_and_neff.csv", index=False)

    # ------------------------------------------------------------------ #
    # 3. the three pre-registered criteria (scope report §5.2)
    # ------------------------------------------------------------------ #
    inc = ax[ax.readout == "T"].set_index("band")
    crit = []
    for m in cs.readouts:
        sub = ax[ax.readout == m].set_index("band")
        # (a) same phenomenon
        rho_beta = np.nan
        for b in ("beta",):
            Mi, _ = cs.margins(b, "T")
            Mm, _ = cs.margins(b, m)
            ai = np.nanmean(Mi, axis=1)
            am = np.nanmean(Mm, axis=1)
            ok = np.isfinite(ai) & np.isfinite(am)
            if ok.sum() >= 5:
                rho_beta = float(spearmanr(ai[ok], am[ok])[0])
        beta_clears = bool(sub.loc["beta", "cluster_q"] < 0.05)
        theta_null = bool(not (sub.loc["theta", "cluster_q"] < 0.05))
        crit_a = bool(np.isfinite(rho_beta) and rho_beta >= 0.5
                      and beta_clears and theta_null)
        # (b) more independent information, above its own noise
        ratio = (sub["n_eff_pr"] / inc["n_eff_pr"]).to_numpy()
        n_abs = int(np.nansum(ratio >= 2.0))
        n_null = int(np.nansum(sub["p_neff_above_null"].to_numpy() < 0.05))
        crit_b = bool(n_abs >= 3 and n_null >= 3)
        # (c) calibrated
        g = gate[gate.readout == m]
        if "fpr_05" in g and g["fpr_05"].notna().any():
            frac_cal = float((g["fpr_05"] <= 0.10).mean())
        else:
            frac_cal = np.nan
        crit_c = bool(np.isfinite(frac_cal) and frac_cal >= 0.95)
        crit.append(dict(
            readout=m, rho_aggregate_vs_incumbent_beta=rho_beta,
            beta_clears=beta_clears, theta_null=theta_null, criterion_a=crit_a,
            n_bands_neff_2x=n_abs, n_bands_neff_above_null=n_null,
            n_eff_pr_median=float(np.nanmedian(sub["n_eff_pr"])),
            n_eff_null_med_median=float(np.nanmedian(sub["n_eff_null_med"])),
            criterion_b=crit_b, frac_cells_calibrated=frac_cal,
            criterion_c=crit_c,
            passes_all=bool(crit_a and crit_b and crit_c)))
    dc = pd.DataFrame(crit)
    dc.to_csv(OUT / "criteria.csv", index=False)

    # ------------------------------------------------------------------ #
    # 4. S3 -- the scale axis on its two halves, incumbent statistic
    # ------------------------------------------------------------------ #
    s3 = ax[ax.readout == "T"][
        ["band", "n_eff_pr", "mean_offdiag", "n_eff_pr_fine",
         "mean_offdiag_fine", "cluster_p", "cluster_q", "cluster_p_fine",
         "cluster_q_fine", "cluster_p_coarse", "cluster_q_coarse",
         "cluster_lo", "cluster_hi"]].copy()
    s3.to_csv(OUT / "s3_fine_vs_coarse.csv", index=False)

    # per-scale incumbent profile, the shape of the curve is itself a result
    prof = gate[gate.readout == "T"][
        ["band", "s", "margin_med", "p", "q", "n_pat", "frac_pos",
         "rank_biserial", "ci_lo", "ci_hi", "loo_p_max", "fpr_05"]]
    prof.to_csv(OUT / "incumbent_scale_profile.csv", index=False)

    # ------------------------------------------------------------------ #
    print("\n-- criteria (scope report §5.2), every construction attempted --")
    print(dc[["readout", "rho_aggregate_vs_incumbent_beta", "criterion_a",
              "n_eff_pr_median", "n_eff_null_med_median",
              "n_bands_neff_2x", "n_bands_neff_above_null", "criterion_b",
              "frac_cells_calibrated", "criterion_c", "passes_all"]]
          .to_string(index=False))
    print("\n-- S3: the axis on its two halves (incumbent) --")
    print(s3.to_string(index=False))
    print("\n-- does the effect vary along the axis? (null-referenced) --")
    print(ax[ax.readout.isin(HEADLINE)][
        ["readout", "band", "slope_abs", "slope_abs_null_med",
         "p_slope_above_null", "shape_agreement",
         "shape_agreement_null_med", "p_shape_above_null"]]
        .to_string(index=False))
    print(f"\n[w0s-verdict] -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
