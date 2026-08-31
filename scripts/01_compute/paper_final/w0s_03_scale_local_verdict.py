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


def _neff_null(cs, band: str, readout: str, n_draws: int) -> np.ndarray:
    """``n_eff_pr`` of held-out surrogate margin matrices -- the readout's noise floor."""
    _, S, _ = cs.cell(band, readout)
    R = S.shape[1]
    out = np.full(min(n_draws, R), np.nan)
    for i in range(out.size):
        out[i] = effective_tests(cs.heldout_margins(band, readout, i))["n_eff_pr"]
    return out


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
            nn = _neff_null(cs, b, m, N_NEFF_DRAWS)
            nn = nn[np.isfinite(nn)]
            p_neff = (float((1 + int((nn >= e["n_eff_pr"]).sum())) / (nn.size + 1))
                      if nn.size and np.isfinite(e["n_eff_pr"]) else np.nan)
            e_f = effective_tests(M[:, fine])
            rows.append(dict(
                readout=m, band=b,
                n_eff_pr=e["n_eff_pr"], n_eff_cn=e["n_eff_cn"],
                mean_offdiag=e["mean_offdiag"], n_finite=e["n_finite"],
                n_eff_pr_fine=e_f["n_eff_pr"], mean_offdiag_fine=e_f["mean_offdiag"],
                n_eff_null_med=float(np.median(nn)) if nn.size else np.nan,
                n_eff_null_p95=float(np.percentile(nn, 95)) if nn.size else np.nan,
                p_neff_above_null=p_neff,
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
                rho_beta = float(spearmanr(ai[ok], am[ok]).statistic)
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
    print(f"\n[w0s-verdict] -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
