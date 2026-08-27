#!/usr/bin/env python3
r"""W0-C part 1: lock the cohort gate, and calibrate it.

Reads the (obs, surrogate-ensemble) cells written by ``w0c_01_gate_and_tau_grid``
and runs the canonical gate from
:mod:`lrg_eegfc.utils.metrics.cohort_gate` over the whole (band x scale) grid.

Four things this script establishes, in order:

1. **The confound is real and the margin fixes it.** Per (band, scale) cell,
   the Spearman across patients between the observed value and that patient's
   own surrogate median, and the same for the margin. If the first is large and
   the second is not, testing margins rather than values is not a stylistic
   preference.
2. **The grid, gated.** Every (band, scale) cell, whole-grid Benjamini-Hochberg,
   per-scale, no best-of. Plus the raw dense reference as its own row.
3. **Leave-one-patient-out of the verdict.** The whole grid, BH step included,
   recomputed with each patient dropped; the report is how many cells change
   verdict and which patient costs the most.
4. **Calibration.** The gate is run on inputs where the null is true by
   construction -- one surrogate realization promoted to "observed", tested
   against its own remaining ensemble -- so the false-positive rate quoted for
   the gate is measured on the real graphs. A Gaussian toy is run as a
   cross-check, and theta is inspected as the substantive known-null band.

Outputs (data/paper_final/w0c_gate_tau/):
  gate_grid.csv         per (measure, band, s): the full gate output + q
  gate_raw.csv          the dense-raw reference row per band
  loo_patient.csv       per dropped patient: cells lost / gained / flipped
  loo_cell.csv          per cell: how many single drops cost it its verdict
  calibration.csv       per (band, s): KS-vs-uniform and FPR of the gate
  covariation.csv       per (band, s): rho(obs, null) and rho(margin, null)
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.utils.metrics.cohort_gate import (
    DESCRIPTIVE_ALPHA, axis_cluster_gate, calibrate_from_surrogates,
    calibrate_synthetic, effective_tests, gate_grid, loo_grid_flips,
    null_height_covariation, patient_margin)
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
BASE = Path(os.environ.get("W0C_BASE", ROOT / "data" / "paper_final" / "w0c_gate_tau"))
GRID = BASE / "grid"


def load_cells(measure_index: int, measure: str):
    """Yield ``(keys, obs, surr)`` per (band, scale) for one measure, plus metadata."""
    cfg = json.loads((GRID / "config.json").read_text())
    pats, bands = cfg["cohort"], cfg["bands"]
    s = np.array(cfg["s_grid"])
    obs = np.full((len(bands), len(pats), s.size), np.nan)
    surr = np.full((len(bands), len(pats), s.size, cfg["R"]), np.nan)
    obs_raw = np.full((len(bands), len(pats)), np.nan)
    surr_raw = np.full((len(bands), len(pats), cfg["R"]), np.nan)
    extra = {k: np.full((len(bands), len(pats), s.size), np.nan)
             for k in ("absorb_r2", "n_eff", "m_comm")}
    for ib, b in enumerate(bands):
        for ip, p in enumerate(pats):
            f = GRID / "cells" / f"{p}__{b}.npz"
            if not f.exists():
                continue
            d = np.load(f, allow_pickle=False)
            obs[ib, ip] = d["obs"][:, measure_index]
            surr[ib, ip] = d["surr"][:, :, measure_index].T          # (nS, R)
            obs_raw[ib, ip] = d["obs_raw"][0]
            surr_raw[ib, ip] = d["surr_raw"]
            for k in extra:
                extra[k][ib, ip] = d[k]
    return dict(cfg=cfg, pats=pats, bands=bands, s=s, obs=obs, surr=surr,
                obs_raw=obs_raw, surr_raw=surr_raw, measure=measure, **extra)


def cells_for(D, band_idx=None):
    """Grid cells as ``(keys, obs(K,), surr(K,R))`` tuples."""
    out = []
    for ib, b in enumerate(D["bands"]):
        if band_idx is not None and ib not in band_idx:
            continue
        for js, s in enumerate(D["s"]):
            out.append(({"measure": D["measure"], "band": b, "s": float(s)},
                        D["obs"][ib, :, js], D["surr"][ib, :, js, :]))
    return out


def main():
    BASE.mkdir(parents=True, exist_ok=True)
    D = load_cells(0, "coph")
    pats, bands, s = D["pats"], D["bands"], D["s"]
    print(f"[w0c-gate] {len(bands)} bands x {s.size} scales x {len(pats)} patients "
          f"| R={D['cfg']['R']} | {D['cfg']['backbone']}@{D['cfg']['frac']}", flush=True)

    # ---- 1. the null-height confound ------------------------------------- #
    cov_rows = []
    for ib, b in enumerate(bands):
        for js, sv in enumerate(s):
            c = null_height_covariation(D["obs"][ib, :, js], D["surr"][ib, :, js, :])
            cov_rows.append(dict(band=b, s=float(sv), **c))
    cov = pd.DataFrame(cov_rows)
    cov.to_csv(BASE / "covariation.csv", index=False)
    print("\n=== 1. obs-vs-own-null covariation (cross-patient Spearman) ===", flush=True)
    print("  band        rho(obs,null) med [IQR]        rho(margin,null) med [IQR]", flush=True)
    for b in bands:
        x = cov[cov.band == b]
        print(f"  {b:11s} {x.rho_obs.median():+6.3f} "
              f"[{x.rho_obs.quantile(.25):+.2f},{x.rho_obs.quantile(.75):+.2f}]"
              f"            {x.rho_margin.median():+6.3f} "
              f"[{x.rho_margin.quantile(.25):+.2f},{x.rho_margin.quantile(.75):+.2f}]",
              flush=True)

    # ---- 2. the gated grid ------------------------------------------------ #
    cells = cells_for(D)
    # expect_n makes the no-dropout policy enforceable: a short cohort raises
    # unless the caller supplies a reason. Nothing here supplies one.
    g = gate_grid(cells, labels=pats, rng=np.random.default_rng(0),
                  expect_n=len(PATIENTS_4PHASE))
    g.to_csv(BASE / "gate_grid.csv", index=False)

    raw_cells = [({"measure": "raw", "band": b, "s": 0.0},
                  D["obs_raw"][ib], D["surr_raw"][ib]) for ib, b in enumerate(bands)]
    graw = gate_grid(raw_cells, labels=pats, rng=np.random.default_rng(0),
                     expect_n=len(PATIENTS_4PHASE))
    graw.to_csv(BASE / "gate_raw.csv", index=False)

    print(f"\n=== 2. cohort gate, whole-grid BH over {len(g)} cells "
          f"(q < {DESCRIPTIVE_ALPHA}) ===", flush=True)
    print("  band        raw(q)   scales cleared   q_min    s@q_min   n_eff@q_min  "
          "margin_med@q_min", flush=True)
    for ib, b in enumerate(bands):
        x = g[g.band == b].sort_values("s")
        r = graw[graw.band == b]
        ncl = int((x.q < DESCRIPTIVE_ALPHA).sum())
        i = x.q.idxmin()
        js = int(np.argmin(np.abs(s - x.loc[i, "s"])))
        ne = np.nanmedian(D["n_eff"][ib, :, js])
        rq = float(r.q.iloc[0]) if len(r) else np.nan
        print(f"  {b:11s} {rq:6.3f}   {ncl:2d}/{s.size:2d}          "
              f"{x.loc[i,'q']:6.4f}  {x.loc[i,'s']:7.2f}   {ne:8.1f}     "
              f"{x.loc[i,'margin_med']:+.4f}", flush=True)

    print("\n  per-scale q (rows = bands, cols = s; * = q < 0.05):", flush=True)
    hdr = "  band        " + "".join(f"{v:>8.2f}" for v in s)
    print(hdr, flush=True)
    for b in bands:
        x = g[g.band == b].sort_values("s")
        cells_txt = "".join(
            f"{'*' if q < DESCRIPTIVE_ALPHA else ' '}{q:6.3f} " if np.isfinite(q)
            else f"{'--':>8s}" for q in x.q)
        print(f"  {b:11s} " + cells_txt, flush=True)

    # ---- 2b. the multiplicity family -------------------------------------- #
    print("\n=== 2b. how many tests is a 28-point scale sweep actually worth? ===",
          flush=True)
    print("  Whole-grid BH charges one test per (band, scale) cell. A scale sweep is", flush=True)
    print("  a smooth curve, so that price is far above what was really paid. Below:", flush=True)
    print("  the effective number of independent scales, and the axis-cluster gate --", flush=True)
    print("  ONE sign-flip cluster-mass test per band, family = 6 bands, not 168 cells.",
          flush=True)
    prof = {b: patient_margin(D["obs"][ib].ravel(),
                              D["surr"][ib].reshape(-1, D["surr"].shape[-1])
                              ).reshape(len(pats), s.size)
            for ib, b in enumerate(bands)}
    crows = []
    for ib, b in enumerate(bands):
        et = effective_tests(prof[b])
        ac = axis_cluster_gate(prof[b], n_perm=10000, rng=np.random.default_rng(7))
        crows.append(dict(band=b, n_eff_pr=et["n_eff_pr"], n_eff_cn=et["n_eff_cn"],
                          mean_offdiag=et["mean_offdiag"], cluster_p=ac["p"],
                          mass=ac["mass"], n_clusters=ac["n_clusters"],
                          s_lo=float(s[ac["cluster"][0]]) if ac["cluster"] else np.nan,
                          s_hi=float(s[ac["cluster"][1]]) if ac["cluster"] else np.nan))
    cdf = pd.DataFrame(crows)
    cdf["cluster_q"] = bh_fdr(cdf["cluster_p"].to_numpy())
    cdf.to_csv(BASE / "axis_cluster.csv", index=False)
    print("\n  band        n_eff of 28 scales (PR / CN)  mean r   cluster_p  "
          "cluster_q   supra-threshold scale span", flush=True)
    for _, r in cdf.iterrows():
        span = f"{r.s_lo:6.2f} - {r.s_hi:6.2f}" if np.isfinite(r.s_lo) else "     none     "
        print(f"  {r.band:11s} {r.n_eff_pr:6.2f} / {r.n_eff_cn:6.2f}          "
              f"{r.mean_offdiag:+.2f}   {r.cluster_p:8.4f}   {r.cluster_q:8.4f}   {span}",
              flush=True)

    # ---- 3. leave-one-patient-out of the verdict -------------------------- #
    per_pat, per_cell = loo_grid_flips(cells, labels=pats, q_level=DESCRIPTIVE_ALPHA)
    per_pat.to_csv(BASE / "loo_patient.csv", index=False)
    per_cell.to_csv(BASE / "loo_cell.csv", index=False)
    print("\n=== 3. leave-one-patient-out of the VERDICT (whole grid + BH re-run) ===",
          flush=True)
    print(f"  full cohort clears {int(per_pat.n_clear_full.iloc[0])}/{len(cells)} cells",
          flush=True)
    for _, r in per_pat.sort_values("n_lost", ascending=False).iterrows():
        print(f"    drop {r.dropped}: clears {int(r.n_clear_loo):3d}  "
              f"lost {int(r.n_lost):3d}  gained {int(r.n_gained):3d}", flush=True)
    print("  cells whose verdict survives EVERY single-patient drop: "
          f"{int(((per_cell.clear_full) & (per_cell.n_drop_lost == 0)).sum())}"
          f"/{int(per_cell.clear_full.sum())} of the cleared cells", flush=True)
    for b in bands:
        x = per_cell[(per_cell.band == b) & per_cell.clear_full]
        if len(x):
            print(f"    {b:11s} {int((x.n_drop_lost == 0).sum())}/{len(x)} robust "
                  f"| worst cell loses to {int(x.n_drop_lost.max())} drops", flush=True)

    # ---- 3b. LOO of the verdict under the honest family -------------------- #
    print("\n=== 3b. LOO of the verdict under the PER-BAND family (axis cluster) ===",
          flush=True)
    print("  Section 3 drops a patient from a 168-cell BH whose p-values sit on the", flush=True)
    print("  discrete exact signed-rank grid; at n=9 the floor moves from 1/1024 to", flush=True)
    print("  1/512 and the whole grid falls off the cliff at once. That is the", flush=True)
    print("  multiplicity design failing, not the signal. Repeated here on the", flush=True)
    print("  axis-cluster gate, whose permutation p is continuous and whose family is", flush=True)
    print("  6 bands:", flush=True)
    print("\n  band        full p     LOO p range              worst drop   "
          "drops where p >= 0.05", flush=True)
    lrows = []
    for ib, b in enumerate(bands):
        M = prof[b]
        full_p = axis_cluster_gate(M, n_perm=10000,
                                   rng=np.random.default_rng(7))["p"]
        loo = []
        for i in range(len(pats)):
            keep = np.ones(len(pats), bool)
            keep[i] = False
            loo.append(axis_cluster_gate(M[keep], n_perm=10000,
                                         rng=np.random.default_rng(11 + i))["p"])
        loo = np.array(loo)
        iw = int(np.nanargmax(loo))
        n_fail = int(np.sum(loo >= DESCRIPTIVE_ALPHA))
        lrows.append(dict(band=b, full_p=full_p, loo_min=float(np.nanmin(loo)),
                          loo_max=float(np.nanmax(loo)), worst=pats[iw],
                          n_drop_above_alpha=n_fail,
                          **{f"loo_{pats[j]}": float(loo[j]) for j in range(len(pats))}))
        print(f"  {b:11s} {full_p:8.4f}   [{np.nanmin(loo):.4f}, {np.nanmax(loo):.4f}]"
              f"        {pats[iw]:8s}     {n_fail}/{len(pats)}", flush=True)
    pd.DataFrame(lrows).to_csv(BASE / "axis_cluster_loo.csv", index=False)

    # ---- 4. calibration --------------------------------------------------- #
    if os.environ.get("W0C_SKIP_CAL", "0") in ("1", "true", "True"):
        print("\n=== 4. calibration SKIPPED (W0C_SKIP_CAL) ===", flush=True)
        print(f"\n[w0c-gate] -> {BASE}", flush=True)
        return
    print("\n=== 4. calibration ===", flush=True)
    syn_h = calibrate_synthetic(n_draws=4000, heteroscedastic=True)
    syn_o = calibrate_synthetic(n_draws=4000, heteroscedastic=False)
    for nm, c in (("synthetic heteroscedastic", syn_h), ("synthetic homoscedastic", syn_o)):
        print(f"  {nm:28s} FPR@.10={c['fpr_0.1']:.4f} @.05={c['fpr_0.05']:.4f} "
              f"@.01={c['fpr_0.01']:.4f}   (KS D={c['ks_D']:.3f}; the exact "
              f"signed-rank p is discrete at n=10, so a KS test against the "
              f"continuous uniform rejects by construction -- read the FPRs)",
              flush=True)

    cal_rows = []
    rng = np.random.default_rng(1)
    for ib, b in enumerate(bands):
        for js, sv in enumerate(s):
            c = calibrate_from_surrogates(D["surr"][ib, :, js, :], n_draws=200, rng=rng)
            cal_rows.append(dict(band=b, s=float(sv), ks_D=c["ks_D"],
                                 fpr_10=c["fpr_0.1"], fpr_05=c["fpr_0.05"],
                                 fpr_01=c["fpr_0.01"], n=c["n_draws"]))
    cal = pd.DataFrame(cal_rows)
    cal.to_csv(BASE / "calibration.csv", index=False)
    print(f"  held-out-realization null on the REAL graphs "
          f"({len(cal)} cells x {int(cal.n.median())} draws):", flush=True)
    print(f"    FPR@.10 = {cal.fpr_10.mean():.3f} (cell range "
          f"{cal.fpr_10.min():.3f}-{cal.fpr_10.max():.3f})", flush=True)
    print(f"    FPR@.05 = {cal.fpr_05.mean():.3f} (cell range "
          f"{cal.fpr_05.min():.3f}-{cal.fpr_05.max():.3f})", flush=True)
    print(f"    FPR@.01 = {cal.fpr_01.mean():.3f} (cell range "
          f"{cal.fpr_01.min():.3f}-{cal.fpr_01.max():.3f})", flush=True)
    print(f"    worst cell FPR@.05 across the grid: {cal.fpr_05.max():.3f}; "
          f"cells with FPR@.05 above .10: {int((cal.fpr_05 > 0.10).sum())}/{len(cal)}",
          flush=True)

    print(f"\n[w0c-gate] -> {BASE}", flush=True)


if __name__ == "__main__":
    main()
