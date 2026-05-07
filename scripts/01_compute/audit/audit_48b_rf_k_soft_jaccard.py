#!/usr/bin/env python3
"""Audit 48b — Soft (graded) RF(k): mean best Jaccard, no threshold gate.

For each (patient, band) at tau = 1 / lambda_max, cut the dendrogram at every
meaningful k. For each clade in C^phi_A(k), take its best Jaccard match in
C^phi_B(k); report the MEAN over clades (continuous, in [0, 1]):

    P_soft(phi_A, phi_B; k) = (1 / |C^phi_A(k)|) *
                              sum_{C in C^phi_A(k)} max_{C' in C^phi_B(k)} J(C, C')

Triangle scalar at k:

    T_RFsoft(k; p, b) = P_soft(taskt, rsPost; k) - P_soft(rsPre, taskt; k)

Sign convention: positive = trace direction.

Cohort scalar = mean over k in [2, N - 5]. This addresses RF(theta)'s
threshold-rejection of near-same clades: a clade with best-Jaccard 0.65
contributes 0.65 to soft RF instead of 0 (under theta=0.70).

Outputs (under data/reports/section_5_lrg_trace/14_rf_clade_persistence/):
  tables/Td_per_patient_per_band_soft.csv
  tables/per_k_long_soft.csv
  tables/cohort_summary_soft.csv
  tables/pat03_dropout_soft.csv
  figures/rf_k_soft_cohort_heatmap.pdf

Scope: .agents/guides/task-persistence-investigation/2026-05-06_rf-k-multiscale-measure.md
       (graded variant added 2026-05-07 in response to RF threshold mismatch
       at beta vs KC lambda=0).
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
HALVES = CACHE_ROOT / "imcoh_lrg_halves"

OUT_TABLES = ROOT / "data" / "reports" / "section_5_lrg_trace" / "14_rf_clade_persistence" / "tables"
OUT_FIGS = ROOT / "data" / "reports" / "section_5_lrg_trace" / "14_rf_clade_persistence" / "figures"


# ---------------------------------------------------------------------------
# Soft persistence via vectorised contingency table
# ---------------------------------------------------------------------------

def best_jaccards(labels_A: np.ndarray, labels_B: np.ndarray) -> np.ndarray:
    """For each cluster in labels_A, return best Jaccard match in labels_B."""
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


def t_rfsoft_grid(Z_a: np.ndarray, Z_mid: np.ndarray, Z_b: np.ndarray,
                  k_grid: np.ndarray) -> np.ndarray:
    """T_RFsoft(k) = mean(best_J(Z_mid->Z_b)) - mean(best_J(Z_a->Z_mid))."""
    out = np.empty(len(k_grid), dtype=float)
    for i, k in enumerate(k_grid):
        L_a = fcluster(Z_a, t=int(k), criterion="maxclust")
        L_m = fcluster(Z_mid, t=int(k), criterion="maxclust")
        L_b = fcluster(Z_b, t=int(k), criterion="maxclust")
        bj_am = best_jaccards(L_a, L_m)
        bj_mb = best_jaccards(L_m, L_b)
        out[i] = float(bj_mb.mean()) - float(bj_am.mean())
    return out


def k_range_for(n_leaves: int) -> np.ndarray:
    return np.arange(2, max(2, n_leaves - 5) + 1, dtype=int)


def load_phase(pat: str, band: str, phase: str, half: bool = False):
    if half:
        return load_lrg_result(pat, phase, band, "imcoh_abs", cache_root=HALVES)
    return load_lrg_result(pat, phase, band, "imcoh_abs")


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
            t_real = t_rfsoft_grid(Z_pre, Z_tt, Z_post, k_grid)
            t_real_cohort = float(np.nanmean(t_real))

            res_pa = load_phase(pat, band, "rest_pre_A", half=True)
            res_pb = load_phase(pat, band, "rest_pre_B", half=True)
            res_post_a = load_phase(pat, band, "rest_post_A", half=True)
            has_null = all(r is not None for r in (res_pa, res_pb, res_post_a))
            if has_null:
                Z_pa = res_pa.linkage_matrix
                Z_pb = res_pb.linkage_matrix
                Z_pst_a = res_post_a.linkage_matrix
                k_grid_null = k_range_for(int(Z_pa.shape[0]) + 1)
                t_null = t_rfsoft_grid(Z_pa, Z_pb, Z_pst_a, k_grid_null)
                t_null_cohort = float(np.nanmean(t_null))
            else:
                t_null_cohort = float("nan")

            summary_rows.append({
                "patient": pat,
                "band": band,
                "n_leaves": n_leaves,
                "k_min": int(k_grid[0]),
                "k_max": int(k_grid[-1]),
                "T_RFsoft_cohort": t_real_cohort,
                "T_RFsoft_null_cohort": t_null_cohort,
                "T_RFsoft_per_k_json": json.dumps([round(float(v), 6) for v in t_real]),
            })
            for k_i, k in enumerate(k_grid):
                long_rows.append({
                    "patient": pat,
                    "band": band,
                    "k": int(k),
                    "T_RFsoft": float(t_real[k_i]),
                })
            print(f"  {band}: n_leaves={n_leaves} T={t_real_cohort:+.4f} null={t_null_cohort:+.4f}")

    summary_df = pd.DataFrame(summary_rows)
    long_df = pd.DataFrame(long_rows)
    summary_df.to_csv(OUT_TABLES / "Td_per_patient_per_band_soft.csv", index=False)
    long_df.to_csv(OUT_TABLES / "per_k_long_soft.csv", index=False)

    # ---- cohort summary ----
    coh_rows = []
    for band in BANDS:
        sub = summary_df[summary_df["band"] == band]
        real = sub["T_RFsoft_cohort"].values.astype(float)
        null = sub["T_RFsoft_null_cohort"].values.astype(float)
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
        real_drop = sub_drop["T_RFsoft_cohort"].values.astype(float)
        null_drop = sub_drop["T_RFsoft_null_cohort"].values.astype(float)
        p_drop = wilcoxon_one_sided_greater(real_drop)
        p_drop_null = wilcoxon_one_sided_greater(real_drop, null_drop) \
            if not np.any(np.isnan(null_drop)) else float("nan")
        coh_rows.append({
            "band": band,
            "n_pat": n,
            "n_trace": n_trace,
            "T_RFsoft_median": float(np.median(real)) if n else float("nan"),
            "T_RFsoft_null_median": float(np.nanmedian(null)) if n else float("nan"),
            "wilcoxon_p_abs": p_abs,
            "rb_abs": rb_abs,
            "n_above_null": n_above_null,
            "wilcoxon_p_null": p_null,
            "rb_null": rb_null,
            "wilcoxon_p_drop": p_drop,
            "wilcoxon_p_drop_null": p_drop_null,
        })
    coh = pd.DataFrame(coh_rows)
    p_abs_arr = coh["wilcoxon_p_abs"].values.astype(float)
    p_null_arr = coh["wilcoxon_p_null"].values.astype(float)
    coh["q_abs_within_m6"] = bh_fdr(p_abs_arr) if not np.any(np.isnan(p_abs_arr)) else np.full_like(p_abs_arr, np.nan)
    valid_pn = ~np.isnan(p_null_arr)
    q_null_full = np.full(len(coh), np.nan)
    if valid_pn.sum() > 0:
        q_null_full[valid_pn] = bh_fdr(p_null_arr[valid_pn])
    coh["q_null_within_m6"] = q_null_full
    coh.to_csv(OUT_TABLES / "cohort_summary_soft.csv", index=False)

    # ---- Pat_03 dropout long-form ----
    pd.DataFrame([{
        "band": r["band"],
        "n_full": r["n_pat"],
        "p_abs_full": r["wilcoxon_p_abs"],
        "p_abs_drop": r["wilcoxon_p_drop"],
        "p_null_full": r["wilcoxon_p_null"],
        "p_null_drop": r["wilcoxon_p_drop_null"],
    } for r in coh_rows]).to_csv(OUT_TABLES / "pat03_dropout_soft.csv", index=False)

    # ---- cohort heatmap n_trace(k)/10 ----
    if not summary_df.empty:
        min_kmax = int(summary_df["k_max"].min())
        common_k = np.arange(2, min_kmax + 1)
        H = np.full((len(BANDS), len(common_k)), np.nan)
        for bi, band in enumerate(BANDS):
            sub = long_df[long_df["band"] == band]
            for ki, k in enumerate(common_k):
                vals = sub[sub["k"] == k]["T_RFsoft"].values.astype(float)
                if len(vals) > 0:
                    H[bi, ki] = float(np.sum(vals > 0)) / len(vals)
        fig, ax = plt.subplots(1, 1, figsize=(8.4, 3.8))
        im = ax.imshow(
            H, aspect="auto", cmap="RdBu_r", vmin=0, vmax=1, origin="lower",
            extent=[common_k[0] - 0.5, common_k[-1] + 0.5, -0.5, len(BANDS) - 0.5],
            interpolation="nearest",
        )
        ax.set_yticks(np.arange(len(BANDS)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS])
        ax.set_title(r"soft RF$(k)$  —  cohort $n_{trace}(k)/10$  (continuous mean best-Jaccard)")
        ax.set_xlabel("$k$ (cut level, # clusters)")
        imshow_colorbar_caxdivider(im, ax, size="2.5%", pad=0.04)
        fig.tight_layout()
        fig.savefig(OUT_FIGS / "rf_k_soft_cohort_heatmap.pdf")
        plt.close(fig)

    print("\n=== cohort summary (soft RF) ===")
    print(coh.to_string(index=False))
    print(f"\n[audit_48b] tables -> {OUT_TABLES}")
    print(f"           figure  -> {OUT_FIGS / 'rf_k_soft_cohort_heatmap.pdf'}")


if __name__ == "__main__":
    main()
