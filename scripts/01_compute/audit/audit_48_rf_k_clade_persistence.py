#!/usr/bin/env python3
"""Audit 48 — RF(k) clade-persistence triangle (Section 5.4 multiscale companion).

For each (patient, band) cell at tau = 1 / lambda_max, cut the dendrogram at
every meaningful k and count the fraction of clades from one phase that
survive into the next under Jaccard >= theta. The triangle scalar at k:

    T_RF(k; p, b) = P(taskt, rsPost; k) - P(rsPre, taskt; k)

with sign convention T_RF > 0 = trace direction (more clades persist from
task into rest_post than from rest_pre into task at this cut level). Cohort
scalar is the mean across k in [2, N(p) - 5]. Within-baseline null replaces
taskt with a second rest_pre half:

    T_RF_null(k; p, b) = P(rsPre_B, rsPost_A; k) - P(rsPre_A, rsPre_B; k)

Two thresholds reported: theta in {0.70, 0.85}.

Outputs
-------
data/reports/section_5_lrg_trace/14_rf_clade_persistence/
  tables/
    Td_per_patient_per_band.csv  per (patient, band, threshold) row
    per_k_long.csv               per (patient, band, threshold, k) row
    cohort_summary.csv           per (band, threshold) row
    pat03_dropout.csv            Pat_03 dropout under absolute + null gates
  figures/
    rf_k_cohort_heatmap.pdf      cohort n_trace(k)/10 heatmap

Scope: .agents/guides/task-persistence-investigation/2026-05-06_rf-k-multiscale-measure.md
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr, rank_biserial
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrgsglib.plotlib.colorbars import imshow_colorbar_caxdivider

PATIENTS = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
BANDS = list(BRAIN_BANDS_NAMES)
THETAS = (0.70, 0.85)
HALVES = CACHE_ROOT / "imcoh_lrg_halves"

OUT_TABLES = ROOT / "data" / "reports" / "section_5_lrg_trace" / "14_rf_clade_persistence" / "tables"
OUT_FIGS = ROOT / "data" / "reports" / "section_5_lrg_trace" / "14_rf_clade_persistence" / "figures"


# ---------------------------------------------------------------------------
# Clade persistence via contingency table (vectorised)
# ---------------------------------------------------------------------------

def best_jaccards(labels_A: np.ndarray, labels_B: np.ndarray) -> np.ndarray:
    """For each cluster in labels_A, return its best Jaccard match in labels_B.

    Vectorised: contingency table M[i, j] = |C_A^i intersect C_B^j|, sizes via
    row/column sums, J[i, j] = M / (|C_A^i| + |C_B^j| - M); take the row max.
    """
    uniq_a, idx_a = np.unique(labels_A, return_inverse=True)
    uniq_b, idx_b = np.unique(labels_B, return_inverse=True)
    nA = len(uniq_a)
    nB = len(uniq_b)
    M = np.zeros((nA, nB), dtype=np.int64)
    np.add.at(M, (idx_a, idx_b), 1)
    sizes_a = M.sum(axis=1)
    sizes_b = M.sum(axis=0)
    union = sizes_a[:, None] + sizes_b[None, :] - M
    with np.errstate(divide="ignore", invalid="ignore"):
        J = np.where(union > 0, M / union, 0.0)
    return J.max(axis=1)


def t_rf_grid(Z_a: np.ndarray, Z_mid: np.ndarray, Z_b: np.ndarray,
              k_grid: np.ndarray, thetas: tuple[float, ...]) -> dict[float, np.ndarray]:
    """T_RF(k) = P(Z_mid, Z_b; k) - P(Z_a, Z_mid; k) for every theta in thetas."""
    out = {theta: np.empty(len(k_grid), dtype=float) for theta in thetas}
    for i, k in enumerate(k_grid):
        L_a = fcluster(Z_a, t=int(k), criterion="maxclust")
        L_m = fcluster(Z_mid, t=int(k), criterion="maxclust")
        L_b = fcluster(Z_b, t=int(k), criterion="maxclust")
        bj_am = best_jaccards(L_a, L_m)
        bj_mb = best_jaccards(L_m, L_b)
        for theta in thetas:
            P_am = float(np.sum(bj_am >= theta) / len(bj_am)) if len(bj_am) else float("nan")
            P_mb = float(np.sum(bj_mb >= theta) / len(bj_mb)) if len(bj_mb) else float("nan")
            out[theta][i] = P_mb - P_am
    return out


def k_range_for(n_leaves: int) -> np.ndarray:
    """Meaningful cut levels: k = 2..(n_leaves - 5)."""
    return np.arange(2, max(2, n_leaves - 5) + 1, dtype=int)


# ---------------------------------------------------------------------------
# LRG load
# ---------------------------------------------------------------------------

def load_phase(pat: str, band: str, phase: str, half: bool = False):
    if half:
        return load_lrg_result(pat, phase, band, "imcoh_abs", cache_root=HALVES)
    return load_lrg_result(pat, phase, band, "imcoh_abs")


# ---------------------------------------------------------------------------
# Wilcoxon helpers (one-sided "greater" because trace = positive here)
# ---------------------------------------------------------------------------

def wilcoxon_one_sided_greater(x: np.ndarray, y: np.ndarray | None = None) -> float:
    try:
        if y is None:
            return float(wilcoxon(x, alternative="greater", zero_method="wilcox").pvalue)
        return float(wilcoxon(x, y, alternative="greater", zero_method="wilcox").pvalue)
    except ValueError:
        return float("nan")


def rb_paired_greater(x: np.ndarray, y: np.ndarray | None = None) -> float:
    try:
        diff = x if y is None else (np.asarray(x) - np.asarray(y))
        diff = diff[diff != 0]
        if len(diff) == 0:
            return float("nan")
        return float(rank_biserial(diff))
    except Exception:
        return float("nan")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    OUT_TABLES.mkdir(parents=True, exist_ok=True)
    OUT_FIGS.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict] = []
    long_rows: list[dict] = []

    for pat in PATIENTS:
        print(f"=== {pat} ===")
        for band in BANDS:
            res_pre = load_phase(pat, band, "rest_pre")
            res_tt = load_phase(pat, band, "task_test")
            res_post = load_phase(pat, band, "rest_post")
            if any(r is None for r in (res_pre, res_tt, res_post)):
                print(f"  {band}: real LRG cache miss, skip")
                continue
            Z_pre, Z_tt, Z_post = (
                res_pre.linkage_matrix, res_tt.linkage_matrix, res_post.linkage_matrix,
            )
            n_leaves = int(Z_pre.shape[0]) + 1
            k_grid = k_range_for(n_leaves)

            res_pa = load_phase(pat, band, "rest_pre_A", half=True)
            res_pb = load_phase(pat, band, "rest_pre_B", half=True)
            res_post_a = load_phase(pat, band, "rest_post_A", half=True)
            has_null = all(r is not None for r in (res_pa, res_pb, res_post_a))

            t_real = t_rf_grid(Z_pre, Z_tt, Z_post, k_grid, THETAS)
            if has_null:
                Z_pa = res_pa.linkage_matrix
                Z_pb = res_pb.linkage_matrix
                Z_pst_a = res_post_a.linkage_matrix
                k_grid_null = k_range_for(int(Z_pa.shape[0]) + 1)
                t_null = t_rf_grid(Z_pa, Z_pb, Z_pst_a, k_grid_null, THETAS)
            else:
                t_null = None

            for theta in THETAS:
                t_real_arr = t_real[theta]
                t_real_cohort = float(np.nanmean(t_real_arr))
                if has_null:
                    t_null_cohort = float(np.nanmean(t_null[theta]))
                else:
                    t_null_cohort = float("nan")
                summary_rows.append({
                    "patient": pat,
                    "band": band,
                    "threshold": theta,
                    "n_leaves": n_leaves,
                    "k_min": int(k_grid[0]),
                    "k_max": int(k_grid[-1]),
                    "T_RF_cohort": t_real_cohort,
                    "T_RF_null_cohort": t_null_cohort,
                    "T_RF_per_k_json": json.dumps([round(float(v), 6) for v in t_real_arr]),
                })
                for k_i, k in enumerate(k_grid):
                    long_rows.append({
                        "patient": pat,
                        "band": band,
                        "threshold": theta,
                        "k": int(k),
                        "T_RF": float(t_real_arr[k_i]),
                    })
            print(f"  {band}: n_leaves={n_leaves} k_max={k_grid[-1]} null={'y' if has_null else 'n'}")

    summary_df = pd.DataFrame(summary_rows)
    long_df = pd.DataFrame(long_rows)
    summary_df.to_csv(OUT_TABLES / "Td_per_patient_per_band.csv", index=False)
    long_df.to_csv(OUT_TABLES / "per_k_long.csv", index=False)

    # ---- cohort summary ----
    coh_rows = []
    for band in BANDS:
        for theta in THETAS:
            sub = summary_df[(summary_df["band"] == band) & (summary_df["threshold"] == theta)]
            real = sub["T_RF_cohort"].values.astype(float)
            null = sub["T_RF_null_cohort"].values.astype(float)
            n = len(sub)
            n_trace = int(np.sum(real > 0))
            p_abs = wilcoxon_one_sided_greater(real)
            rb_abs = rb_paired_greater(real)
            valid_null = ~np.isnan(null)
            if n > 0 and valid_null.sum() == n:
                p_null = wilcoxon_one_sided_greater(real, null)
                rb_null = rb_paired_greater(real, null)
                n_above_null = int(np.sum(real > null))
            else:
                p_null, rb_null, n_above_null = float("nan"), float("nan"), -1
            sub_drop = sub[sub["patient"] != "Pat_03"]
            real_drop = sub_drop["T_RF_cohort"].values.astype(float)
            p_drop = wilcoxon_one_sided_greater(real_drop)
            coh_rows.append({
                "band": band,
                "threshold": theta,
                "n_pat": n,
                "n_trace": n_trace,
                "T_RF_median": float(np.median(real)) if n else float("nan"),
                "T_RF_null_median": float(np.nanmedian(null)) if n else float("nan"),
                "wilcoxon_p_abs": p_abs,
                "rb_abs": rb_abs,
                "n_above_null": n_above_null,
                "wilcoxon_p_null": p_null,
                "rb_null": rb_null,
                "wilcoxon_p_drop": p_drop,
            })
    coh = pd.DataFrame(coh_rows)
    # Within-probe BH-FDR at m=6 (6 bands), separately per (theta, layer)
    coh["q_abs_within_m6"] = np.nan
    coh["q_null_within_m6"] = np.nan
    for theta in THETAS:
        mask = coh["threshold"] == theta
        for col, qcol in [("wilcoxon_p_abs", "q_abs_within_m6"),
                          ("wilcoxon_p_null", "q_null_within_m6")]:
            p = coh.loc[mask, col].values.astype(float)
            valid = ~np.isnan(p)
            if valid.sum() == 0:
                continue
            q_full = np.full(p.shape, np.nan)
            q_full[valid] = bh_fdr(p[valid])
            coh.loc[mask, qcol] = q_full
    coh.to_csv(OUT_TABLES / "cohort_summary.csv", index=False)

    # ---- Pat_03 dropout long-form ----
    drop_rows = []
    for band in BANDS:
        for theta in THETAS:
            sub_full = summary_df[(summary_df["band"] == band) & (summary_df["threshold"] == theta)]
            sub_drop = sub_full[sub_full["patient"] != "Pat_03"]
            real_full = sub_full["T_RF_cohort"].values.astype(float)
            real_drop = sub_drop["T_RF_cohort"].values.astype(float)
            null_full = sub_full["T_RF_null_cohort"].values.astype(float)
            null_drop = sub_drop["T_RF_null_cohort"].values.astype(float)
            null_full_ok = not np.any(np.isnan(null_full))
            null_drop_ok = not np.any(np.isnan(null_drop))
            drop_rows.append({
                "band": band,
                "threshold": theta,
                "n_full": len(sub_full),
                "n_drop": len(sub_drop),
                "n_trace_full": int(np.sum(real_full > 0)),
                "n_trace_drop": int(np.sum(real_drop > 0)),
                "p_abs_full": wilcoxon_one_sided_greater(real_full),
                "p_abs_drop": wilcoxon_one_sided_greater(real_drop),
                "p_null_full": wilcoxon_one_sided_greater(real_full, null_full)
                                if null_full_ok else float("nan"),
                "p_null_drop": wilcoxon_one_sided_greater(real_drop, null_drop)
                                if null_drop_ok else float("nan"),
            })
    pd.DataFrame(drop_rows).to_csv(OUT_TABLES / "pat03_dropout.csv", index=False)

    # ---- cohort heatmap n_trace(k)/10 over (band, k), one panel per theta ----
    if not summary_df.empty:
        min_kmax = int(summary_df["k_max"].min())
        common_k = np.arange(2, min_kmax + 1)
        fig, axes = plt.subplots(len(THETAS), 1,
                                 figsize=(8.4, 3.6 * len(THETAS)),
                                 sharex=True)
        if len(THETAS) == 1:
            axes = [axes]
        for ax, theta in zip(axes, THETAS):
            H = np.full((len(BANDS), len(common_k)), np.nan)
            for bi, band in enumerate(BANDS):
                sub = long_df[(long_df["band"] == band) & (long_df["threshold"] == theta)]
                for ki, k in enumerate(common_k):
                    vals = sub[sub["k"] == k]["T_RF"].values.astype(float)
                    if len(vals) > 0:
                        H[bi, ki] = float(np.sum(vals > 0)) / len(vals)
            im = ax.imshow(
                H, aspect="auto", cmap="RdBu_r", vmin=0, vmax=1, origin="lower",
                extent=[common_k[0] - 0.5, common_k[-1] + 0.5,
                        -0.5, len(BANDS) - 0.5],
                interpolation="nearest",
            )
            ax.set_yticks(np.arange(len(BANDS)))
            ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS])
            ax.set_title(rf"$\theta = {theta:.2f}$  —  cohort $n_{{trace}}(k)/10$")
            ax.set_xlabel("$k$ (cut level, # clusters)")
            imshow_colorbar_caxdivider(im, ax, size="2.5%", pad=0.04)
        fig.tight_layout()
        fig.savefig(OUT_FIGS / "rf_k_cohort_heatmap.pdf")
        plt.close(fig)

    # ---- console summary ----
    print("\n=== cohort summary ===")
    print(coh.to_string(index=False))
    print(f"\n[audit_48] tables -> {OUT_TABLES}")
    print(f"          figure  -> {OUT_FIGS / 'rf_k_cohort_heatmap.pdf'}")


if __name__ == "__main__":
    main()
