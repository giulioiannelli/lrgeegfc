#!/usr/bin/env python3
"""Audit 42 -- within-baseline null for demeaned per-leaf rho.

Mirrors audit_41's within-baseline null pattern but at the per-leaf level.
For each (patient, band) we replace the (task, post) contrast with a pure
within-baseline contrast built from rest halves only:

    Real (audit_39 / audit_39b):
        Δ_task[ell, j] = (D_test  - D_pre_A)[ell, j]
        Δ_rest[ell, j] = (D_post  - D_pre_B)[ell, j]
        rho_ell        = Spearman(Δ_task[ell, ·], Δ_rest[ell, ·])

    Null (this script):
        Δ_drift_a[ell, j] = (D_pre_B  - D_pre_A )[ell, j]   # within-RPre noise
        Δ_drift_b[ell, j] = (D_post_B - D_post_A)[ell, j]   # within-RPost noise
        rho_null_ell      = Spearman(Δ_drift_a[ell, ·], Δ_drift_b[ell, ·])

This is the per-leaf analogue of the official h2e_split_half drift floor
(rho_null_drift = Spearman(D^pre_B - D^pre_A, D^post_B - D^post_A)) used
upstream as the CTM null floor.

Both the real and null per-leaf vectors are demeaned within (patient, band)
to remove the patient's global rank-correlation level (the same correction
audit_39b applies to fix Pat_06's apparent localization, which was global
rotation). The 95th percentile of the demeaned null distribution defines
a control-anchored trace-leaf threshold tau_95(p, b), replacing the
descriptive 0.3 / 0.5 thresholds used in measure 08.

Outputs
-------
data/audit/per_leaf_rho_null/leaf_rho_null.csv          per-leaf null rho + demeaned
data/audit/per_leaf_rho_null/calibrated_trace_leaves.csv per-leaf calibration vs real demeaned rho
data/audit/per_leaf_rho_null/cohort_calibration_summary.csv per-(p,b) thresholds + counts
data/outputs/figures/section_5_lrg_trace/per_leaf_rho_null/cohort_threshold_grid.pdf
data/outputs/figures/section_5_lrg_trace/per_leaf_rho_null/cohort_calibrated_summary.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.paths import CACHE_ROOT, IMCOH_LRG_CACHE
from lrg_eegfc.workflow.lrg import load_lrg_result

HALVES_CACHE = CACHE_ROOT / "imcoh_lrg_halves"
OUT_DIR = ROOT / "data" / "audit" / "per_leaf_rho_null"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "per_leaf_rho_null"
LEAF_REAL = ROOT / "data" / "audit" / "per_leaf_rho_demeaned" / "leaf_rho_demeaned.csv"

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}


def load_D(pat: str, phase: str, band: str, cache_root) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, fc_method="imcoh_abs",
                            cache_root=cache_root)
    except Exception:
        return None
    if r is None:
        return None
    um = np.asarray(r.ultrametric_matrix)
    return squareform(um) if um.ndim == 1 else um


def per_leaf_null_rho(pat: str, band: str) -> np.ndarray | None:
    """rho_null_ell = Spearman( (D_preB - D_preA)[ell,j], (D_postB - D_postA)[ell,j] )."""
    D_pa = load_D(pat, "rest_pre_A", band, HALVES_CACHE)
    D_pb = load_D(pat, "rest_pre_B", band, HALVES_CACHE)
    D_qa = load_D(pat, "rest_post_A", band, HALVES_CACHE)
    D_qb = load_D(pat, "rest_post_B", band, HALVES_CACHE)
    if any(D is None for D in (D_pa, D_pb, D_qa, D_qb)):
        return None
    N = D_pa.shape[0]
    if not all(D.shape == (N, N) for D in (D_pb, D_qa, D_qb)):
        return None
    dA = D_pb - D_pa  # within-RPre noise
    dB = D_qb - D_qa  # within-RPost noise
    rho = np.zeros(N, dtype=float)
    for ell in range(N):
        # all pairs containing ell, excluding self
        a = np.concatenate([dA[ell, :ell], dA[ell, ell + 1:]])
        b = np.concatenate([dB[ell, :ell], dB[ell, ell + 1:]])
        if a.std() == 0 or b.std() == 0 or a.size < 3:
            rho[ell] = np.nan
            continue
        r = spearmanr(a, b).statistic
        rho[ell] = float(r) if r is not None and np.isfinite(r) else np.nan
    return rho


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    real_df = pd.read_csv(LEAF_REAL)

    null_rows: list[dict] = []
    for pat in PATIENTS:
        for band in BANDS:
            rho_n = per_leaf_null_rho(pat, band)
            if rho_n is None:
                continue
            mu = float(np.nanmean(rho_n))
            sd = float(np.nanstd(rho_n))
            for ell, r in enumerate(rho_n):
                null_rows.append({
                    "patient": pat, "band": band, "leaf_id": ell,
                    "rho_null": float(r) if not np.isnan(r) else np.nan,
                    "rho_null_mean_pb": mu,
                    "rho_null_demeaned": float(r) - mu if not np.isnan(r) else np.nan,
                })
    null_df = pd.DataFrame(null_rows)
    null_df.to_csv(OUT_DIR / "leaf_rho_null.csv", index=False)
    print(f"[42] wrote {OUT_DIR / 'leaf_rho_null.csv'} ({len(null_df)} rows)")

    # --- per-(patient, band) calibration summary -----------------------------
    cal_rows: list[dict] = []
    for pat in PATIENTS:
        for band in BANDS:
            null_sub = null_df[(null_df.patient == pat) & (null_df.band == band)]
            real_sub = real_df[(real_df.patient == pat) & (real_df.band == band)]
            if null_sub.empty or real_sub.empty:
                continue
            null_dem = null_sub["rho_null_demeaned"].dropna().to_numpy()
            real_dem = real_sub["rho_demeaned"].to_numpy()
            if null_dem.size < 5 or real_dem.size < 5:
                continue
            tau_95 = float(np.percentile(null_dem, 95))
            tau_99 = float(np.percentile(null_dem, 99))
            cal_rows.append({
                "patient": pat, "band": band,
                "n_leaves": int(real_sub.shape[0]),
                "null_mean": float(np.nanmean(null_sub["rho_null"])),
                "null_std": float(np.nanstd(null_sub["rho_null"])),
                "null_dem_p95": tau_95,
                "null_dem_p99": tau_99,
                "real_dem_p95": float(np.percentile(real_dem, 95)),
                "real_dem_max": float(real_dem.max()),
                "n_calib_p95": int((real_dem >= tau_95).sum()),
                "n_calib_p99": int((real_dem >= tau_99).sum()),
                "n_descrip_03": int((real_dem >= 0.30).sum()),
                "frac_calib_p95": float((real_dem >= tau_95).mean()),
                "frac_calib_p99": float((real_dem >= tau_99).mean()),
                "frac_descrip_03": float((real_dem >= 0.30).mean()),
            })
    cal_df = pd.DataFrame(cal_rows)
    cal_df.to_csv(OUT_DIR / "cohort_calibration_summary.csv", index=False)
    print(f"[42] wrote {OUT_DIR / 'cohort_calibration_summary.csv'} ({len(cal_df)} rows)")

    # per-leaf calibrated flag: real_dem >= tau_95(p,b)
    real_idx = real_df.set_index(["patient", "band", "leaf_id"]).copy()
    cal_idx = cal_df.set_index(["patient", "band"])
    flag_rows: list[dict] = []
    for (pat, band), c in cal_idx.iterrows():
        sub = real_df[(real_df.patient == pat) & (real_df.band == band)]
        for _, row in sub.iterrows():
            flag_rows.append({
                "patient": pat, "band": band, "leaf_id": int(row["leaf_id"]),
                "rho_demeaned": float(row["rho_demeaned"]),
                "tau_95": float(c["null_dem_p95"]),
                "tau_99": float(c["null_dem_p99"]),
                "is_trace_p95": int(row["rho_demeaned"] >= c["null_dem_p95"]),
                "is_trace_p99": int(row["rho_demeaned"] >= c["null_dem_p99"]),
                "is_trace_descrip_03": int(row["rho_demeaned"] >= 0.30),
            })
    flag_df = pd.DataFrame(flag_rows)
    flag_df.to_csv(OUT_DIR / "calibrated_trace_leaves.csv", index=False)

    # --- per-band cohort summary --------------------------------------------
    band_rows = []
    for band in BANDS:
        sub = cal_df[cal_df.band == band]
        if sub.empty:
            continue
        # paired Wilcoxon: per-patient real_dem_max vs null_dem_p95
        d_real = sub["real_dem_max"].to_numpy()
        d_null = sub["null_dem_p95"].to_numpy()
        try:
            w, p = wilcoxon(d_real, d_null, alternative="greater",
                            zero_method="wilcox")
        except ValueError:
            w, p = np.nan, np.nan
        band_rows.append({
            "band": band,
            "median_n_calib_p95": float(sub["n_calib_p95"].median()),
            "mean_frac_calib_p95": float(sub["frac_calib_p95"].mean()),
            "median_n_descrip_03": float(sub["n_descrip_03"].median()),
            "mean_frac_descrip_03": float(sub["frac_descrip_03"].mean()),
            "n_pat_ge1_calib_p95": int((sub["n_calib_p95"] >= 1).sum()),
            "n_pat_ge3_calib_p95": int((sub["n_calib_p95"] >= 3).sum()),
            "n_pat_ge5_calib_p95": int((sub["n_calib_p95"] >= 5).sum()),
            "wilcoxon_W_real_max_gt_null_p95": float(w) if not np.isnan(w) else np.nan,
            "wilcoxon_p_real_max_gt_null_p95": float(p) if not np.isnan(p) else np.nan,
        })
    band_df = pd.DataFrame(band_rows)
    band_df.to_csv(OUT_DIR / "cohort_band_summary.csv", index=False)
    print(band_df.to_string(index=False))

    # ---------------- figures -----------------------------------------------
    # (a) cohort threshold grid: tau_95 per (band, patient) heatmap
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    pmat = np.full((len(PATIENTS), len(BANDS)), np.nan)
    cmat = np.full((len(PATIENTS), len(BANDS)), np.nan)
    for i, pat in enumerate(PATIENTS):
        for j, band in enumerate(BANDS):
            r = cal_df[(cal_df.patient == pat) & (cal_df.band == band)]
            if not r.empty:
                pmat[i, j] = float(r.iloc[0]["null_dem_p95"])
                cmat[i, j] = float(r.iloc[0]["n_calib_p95"])
    im = ax.imshow(pmat, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=11)
    ax.set_yticks(range(len(PATIENTS)))
    ax.set_yticklabels([p[-2:] for p in PATIENTS], fontsize=9)
    ax.set_xlabel("band", fontsize=10)
    ax.set_ylabel("patient", fontsize=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$\tau_{95}(\rho_\ell - \bar\rho)_\mathrm{null}$", fontsize=10)
    for i in range(len(PATIENTS)):
        for j in range(len(BANDS)):
            n = int(cmat[i, j]) if not np.isnan(cmat[i, j]) else 0
            ax.text(j, i, str(n), ha="center", va="center",
                    color="white" if pmat[i, j] > np.nanmedian(pmat) else "black",
                    fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "cohort_threshold_grid.pdf")
    plt.close(fig)
    print(f"[42] wrote {FIG_DIR / 'cohort_threshold_grid.pdf'}")

    # (b) cohort calibrated summary: per-band patient counts at three cutoffs
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    x = np.arange(len(BANDS))
    w = 0.28
    pat_ge3 = [int(band_df[band_df.band == b]["n_pat_ge3_calib_p95"].iloc[0])
               if (band_df.band == b).any() else 0 for b in BANDS]
    pat_ge1 = [int(band_df[band_df.band == b]["n_pat_ge1_calib_p95"].iloc[0])
               if (band_df.band == b).any() else 0 for b in BANDS]
    pat_ge5 = [int(band_df[band_df.band == b]["n_pat_ge5_calib_p95"].iloc[0])
               if (band_df.band == b).any() else 0 for b in BANDS]
    ax.bar(x - w, pat_ge1, w, label=r"$\geq 1$ leaf", color="#9ec6e5", edgecolor="#1a4f73")
    ax.bar(x, pat_ge3, w, label=r"$\geq 3$ leaves", color="#1a4f73", edgecolor="black")
    ax.bar(x + w, pat_ge5, w, label=r"$\geq 5$ leaves", color="#0a2d4a", edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=11)
    ax.set_ylabel("patients (out of 10)", fontsize=10)
    ax.set_title(r"calibrated trace-leaves above $\tau_{95}$ null floor", fontsize=10)
    ax.set_ylim(0, 10.5)
    ax.axhline(8, color="0.5", lw=0.5, ls=":")
    ax.text(0, 8.05, "n=8 cohort gate", fontsize=8, color="0.4", va="bottom")
    ax.legend(loc="upper right", frameon=False, fontsize=8)

    # right panel: real vs null demeaned distributions per band (cohort-pooled)
    ax = axes[1]
    real_pool = [real_df[real_df.band == b]["rho_demeaned"].to_numpy()
                 for b in BANDS]
    null_pool = [null_df[null_df.band == b]["rho_null_demeaned"].dropna().to_numpy()
                 for b in BANDS]
    pos = np.arange(len(BANDS))
    bp_n = ax.boxplot(null_pool, positions=pos - 0.18, widths=0.32,
                      patch_artist=True, showfliers=False)
    bp_r = ax.boxplot(real_pool, positions=pos + 0.18, widths=0.32,
                      patch_artist=True, showfliers=False)
    for box in bp_n["boxes"]:
        box.set(facecolor="#cccccc", edgecolor="0.3")
    for box in bp_r["boxes"]:
        box.set(facecolor="#d62728", alpha=0.7, edgecolor="0.3")
    ax.set_xticks(pos)
    ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=11)
    ax.axhline(0, color="0.4", lw=0.4, ls="--")
    ax.set_ylabel(r"$\rho_\ell - \bar\rho_{(p,b)}$", fontsize=10)
    ax.set_title("null (grey) vs real (red), pooled across patients", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "cohort_calibrated_summary.pdf")
    plt.close(fig)
    print(f"[42] wrote {FIG_DIR / 'cohort_calibrated_summary.pdf'}")


if __name__ == "__main__":
    main()
