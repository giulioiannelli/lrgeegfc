#!/usr/bin/env python3
"""Systematic exploration of reorganization metrics across patients and bands.

For each metric, computes a Reorganization Index (RI) per (patient, band):
  RI = D_cross / D_within
where:
  D_within = mean(D(rest_pre,rest_post), D(task_learn,task_test))   [same-condition pairs]
  D_cross  = mean of 4 cross-condition distances             [rest vs task]

RI > 1  →  task phases differ from rest phases  →  reorganization
RI ~ 1  →  no systematic difference              →  persistence

Then checks cross-patient consistency: for each band, do ALL patients agree
on whether RI > 1 or RI ~ 1?

Output: summary table + heatmap figure.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from lrg_eegfc.config import BRAIN_BANDS_NAMES, PHASE_LABELS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrgsglib.utils.basic.linalg import (
    tree_baker_gamma,
    tree_cophenetic_correlation,
    tree_fowlkes_mallows_index,
    tree_robinson_foulds_distance,
    ultrametric_rank_correlation,
    ultrametric_scaled_distance,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
from lrg_eegfc.config.const import PATIENTS_4PHASE
PATIENTS = PATIENTS_4PHASE
BANDS = BRAIN_BANDS_NAMES
PHASES = list(PHASE_LABELS)
from lrg_eegfc.config.paths import LRG_CACHE

# Phase grouping for reorganization index
REST_PHASES = ["rest_pre", "rest_post"]
TASK_PHASES = ["task_learn", "task_test"]
WITHIN_PAIRS = [("rest_pre", "rest_post"), ("task_learn", "task_test")]
CROSS_PAIRS = [
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("rest_post", "task_learn"), ("rest_post", "task_test"),
]

OUTPUT_DIR = FIGURES_ROOT / "metric_exploration"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ===================================================================
# Metric functions: each takes (Z1, Z2, U1_sq, U2_sq, n_nodes) → float
# Returns a DISTANCE (higher = more different)
# ===================================================================

def metric_multiscale_ari(Z1, Z2, U1, U2, n):
    """Integrated (1 - ARI) across normalized threshold."""
    thresholds = np.linspace(0.01, 1.0, 100)
    max_h1, max_h2 = Z1[:, 2].max(), Z2[:, 2].max()
    ari_vals = []
    for t in thresholds:
        l1 = fcluster(Z1, t=t * max_h1, criterion="distance")
        l2 = fcluster(Z2, t=t * max_h2, criterion="distance")
        ari_vals.append(adjusted_rand_score(l1, l2))
    ari_vals = np.array(ari_vals)
    return float(np.trapz(1 - ari_vals, thresholds) / (thresholds[-1] - thresholds[0]))


def metric_multiscale_nmi(Z1, Z2, U1, U2, n):
    """Integrated (1 - NMI) across normalized threshold."""
    thresholds = np.linspace(0.01, 1.0, 100)
    max_h1, max_h2 = Z1[:, 2].max(), Z2[:, 2].max()
    nmi_vals = []
    for t in thresholds:
        l1 = fcluster(Z1, t=t * max_h1, criterion="distance")
        l2 = fcluster(Z2, t=t * max_h2, criterion="distance")
        nmi_vals.append(normalized_mutual_info_score(l1, l2))
    nmi_vals = np.array(nmi_vals)
    return float(np.trapz(1 - nmi_vals, thresholds) / (thresholds[-1] - thresholds[0]))


def metric_multiscale_fm(Z1, Z2, U1, U2, n):
    """Integrated (1 - FM) across k (number of clusters)."""
    from sklearn.metrics.cluster import contingency_matrix
    k_values = np.unique(np.geomspace(2, max(2, n // 2), 60).astype(int))
    fm_vals = []
    for k in k_values:
        l1 = fcluster(Z1, k, criterion="maxclust")
        l2 = fcluster(Z2, k, criterion="maxclust")
        C = contingency_matrix(l1, l2)
        tp = (C * (C - 1)).sum() / 2
        sr = (C.sum(axis=1) * (C.sum(axis=1) - 1)).sum() / 2
        sc = (C.sum(axis=0) * (C.sum(axis=0) - 1)).sum() / 2
        fm = float(tp / np.sqrt(sr * sc)) if sr > 0 and sc > 0 else 0.0
        fm_vals.append(fm)
    fm_vals = np.array(fm_vals)
    log_k = np.log2(k_values)
    return float(np.trapz(1 - fm_vals, log_k) / (log_k[-1] - log_k[0]))


def metric_robinson_foulds(Z1, Z2, U1, U2, n):
    return tree_robinson_foulds_distance(Z1, Z2, normalized=True)


def metric_cophenetic_dist(Z1, Z2, U1, U2, n):
    """1 - cophenetic correlation (convert similarity → distance)."""
    return 1.0 - tree_cophenetic_correlation(Z1, Z2)


def metric_baker_gamma_dist(Z1, Z2, U1, U2, n):
    """1 - Baker's Gamma (convert similarity → distance)."""
    return 1.0 - tree_baker_gamma(Z1, Z2)


def metric_ultrametric_scaled(Z1, Z2, U1, U2, n):
    """Ultrametric scaled distance (log, normalized). Range [0,1]."""
    return ultrametric_scaled_distance(U1, U2, scale="log", normalize=True)


def metric_ultrametric_rank_dist(Z1, Z2, U1, U2, n):
    """1 - Spearman rank correlation of ultrametric entries."""
    corr = ultrametric_rank_correlation(U1, U2, method="spearman")
    return 1.0 - (corr if corr is not None else 0.0)


def metric_log_cophenetic_dist(Z1, Z2, U1, U2, n):
    """1 - Pearson correlation of log-transformed cophenetic distances."""
    c1 = cophenet(Z1)
    c2 = cophenet(Z2)
    log_c1 = np.log(np.clip(c1, 1e-12, None))
    log_c2 = np.log(np.clip(c2, 1e-12, None))
    corr, _ = pearsonr(log_c1, log_c2)
    return 1.0 - corr


def metric_log_cophenetic_spearman_dist(Z1, Z2, U1, U2, n):
    """1 - Spearman correlation of log cophenetic distances."""
    c1 = cophenet(Z1)
    c2 = cophenet(Z2)
    log_c1 = np.log(np.clip(c1, 1e-12, None))
    log_c2 = np.log(np.clip(c2, 1e-12, None))
    corr, _ = spearmanr(log_c1, log_c2)
    return 1.0 - corr


def metric_multiscale_ari_weighted(Z1, Z2, U1, U2, n):
    """Multiscale ARI but weighted by information content.

    Weight each scale by k*(N-k)/N^2 where k = number of clusters,
    emphasizing intermediate partitions (not trivial 1-cluster or N-clusters).
    """
    thresholds = np.linspace(0.01, 1.0, 100)
    max_h1, max_h2 = Z1[:, 2].max(), Z2[:, 2].max()
    weighted_vals = []
    weights = []
    for t in thresholds:
        l1 = fcluster(Z1, t=t * max_h1, criterion="distance")
        l2 = fcluster(Z2, t=t * max_h2, criterion="distance")
        k1 = len(np.unique(l1))
        k2 = len(np.unique(l2))
        k_avg = (k1 + k2) / 2
        # Weight: peaks at k=N/2, zero at k=1 and k=N
        w = k_avg * (n - k_avg) / (n * n)
        ari = adjusted_rand_score(l1, l2)
        weighted_vals.append((1 - ari) * w)
        weights.append(w)
    total_w = np.trapz(weights, thresholds)
    if total_w == 0:
        return 0.0
    return float(np.trapz(weighted_vals, thresholds) / total_w)


def metric_variation_of_info(Z1, Z2, U1, U2, n):
    """Integrated Variation of Information across normalized threshold."""
    thresholds = np.linspace(0.01, 1.0, 100)
    max_h1, max_h2 = Z1[:, 2].max(), Z2[:, 2].max()
    vi_vals = []
    for t in thresholds:
        l1 = fcluster(Z1, t=t * max_h1, criterion="distance")
        l2 = fcluster(Z2, t=t * max_h2, criterion="distance")
        # VI = H(l1|l2) + H(l2|l1) = H(l1) + H(l2) - 2*I(l1,l2)
        # Normalize by log(n)
        from sklearn.metrics import mutual_info_score
        h1 = _entropy(l1)
        h2 = _entropy(l2)
        mi = mutual_info_score(l1, l2)
        vi = (h1 + h2 - 2 * mi) / max(np.log(n), 1e-12)
        vi_vals.append(vi)
    vi_vals = np.array(vi_vals)
    return float(np.trapz(vi_vals, thresholds) / (thresholds[-1] - thresholds[0]))


def _entropy(labels):
    """Shannon entropy of a label vector."""
    _, counts = np.unique(labels, return_counts=True)
    probs = counts / counts.sum()
    return -np.sum(probs * np.log(probs + 1e-12))


# ===================================================================
# All metrics registry
# ===================================================================
METRICS = {
    "MS-ARI": metric_multiscale_ari,
    "MS-NMI": metric_multiscale_nmi,
    "MS-FM": metric_multiscale_fm,
    "MS-ARI-W": metric_multiscale_ari_weighted,
    "MS-VI": metric_variation_of_info,
    "RF": metric_robinson_foulds,
    "1-Coph": metric_cophenetic_dist,
    "1-Baker": metric_baker_gamma_dist,
    "U-Scaled": metric_ultrametric_scaled,
    "1-U-Rank": metric_ultrametric_rank_dist,
    "1-LogCoph": metric_log_cophenetic_dist,
    "1-LogCophS": metric_log_cophenetic_spearman_dist,
}


# ===================================================================
# Analysis
# ===================================================================
def load_all_lrg():
    """Load all LRG results, return dict keyed by (patient, phase, band)."""
    data = {}
    for pat in PATIENTS:
        for band in BANDS:
            for phase in PHASES:
                lrg = load_lrg_result(pat, phase, band, "msc", cache_root=LRG_CACHE)
                if lrg is not None:
                    # Store square ultrametric matrix
                    u = lrg.ultrametric_matrix
                    u_sq = squareform(u) if u.ndim == 1 else u
                    data[(pat, phase, band)] = {
                        "Z": lrg.linkage_matrix,
                        "U": u_sq,
                        "n": lrg.n_nodes,
                    }
    return data


def compute_pairwise_distances(data, metric_fn, patient, band):
    """Compute 4x4 distance matrix for one patient/band using metric_fn."""
    n_phases = len(PHASES)
    D = np.zeros((n_phases, n_phases))
    for i, pi in enumerate(PHASES):
        for j, pj in enumerate(PHASES):
            if i >= j:
                continue
            key_i = (patient, pi, band)
            key_j = (patient, pj, band)
            if key_i not in data or key_j not in data:
                D[i, j] = D[j, i] = np.nan
                continue
            di, dj = data[key_i], data[key_j]
            val = metric_fn(di["Z"], dj["Z"], di["U"], dj["U"], di["n"])
            D[i, j] = D[j, i] = val
    return D


def reorganization_index(D):
    """Compute RI from a 4x4 phase distance matrix.

    RI = D_cross / D_within where:
      within = mean(D[rest_pre,rest_post], D[task_learn,task_test])
      cross  = mean of 4 rest-vs-task distances
    """
    phase_idx = {p: i for i, p in enumerate(PHASES)}

    d_within = np.mean([
        D[phase_idx[a], phase_idx[b]] for a, b in WITHIN_PAIRS
    ])
    d_cross = np.mean([
        D[phase_idx[a], phase_idx[b]] for a, b in CROSS_PAIRS
    ])

    if d_within == 0 or np.isnan(d_within):
        return np.nan
    return d_cross / d_within


def main():
    print("Loading all LRG results...")
    data = load_all_lrg()
    print(f"  Loaded {len(data)} entries")

    # Compute RI for all (metric, patient, band)
    results = []
    for metric_name, metric_fn in METRICS.items():
        print(f"\nComputing {metric_name}...")
        for pat in PATIENTS:
            for band in BANDS:
                D = compute_pairwise_distances(data, metric_fn, pat, band)
                ri = reorganization_index(D)

                # Also store individual distances for inspection
                phase_idx = {p: i for i, p in enumerate(PHASES)}
                d_rest = D[phase_idx["rest_pre"], phase_idx["rest_post"]]
                d_task = D[phase_idx["task_learn"], phase_idx["task_test"]]
                d_cross_vals = [D[phase_idx[a], phase_idx[b]] for a, b in CROSS_PAIRS]
                d_cross = np.mean(d_cross_vals)

                results.append({
                    "metric": metric_name,
                    "patient": pat,
                    "band": band,
                    "RI": ri,
                    "D_within": (d_rest + d_task) / 2,
                    "D_cross": d_cross,
                    "D_rest": d_rest,
                    "D_task": d_task,
                })
                print(f"  {pat}/{band}: RI={ri:.3f}  "
                      f"(rest={d_rest:.4f} task={d_task:.4f} cross={d_cross:.4f})")

    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_DIR / "metric_exploration_results.csv", index=False)

    # ===================================================================
    # Analysis: Which metric gives most consistent cross-patient pattern?
    # ===================================================================
    print("\n" + "=" * 80)
    print("CROSS-PATIENT CONSISTENCY ANALYSIS")
    print("=" * 80)

    # For each metric, compute per-band statistics across patients
    consistency_scores = []
    for metric_name in METRICS:
        mdf = df[df["metric"] == metric_name]
        print(f"\n--- {metric_name} ---")
        print(f"{'Band':<12} {'Mean RI':>8} {'Std RI':>8} {'Min RI':>8} {'Max RI':>8} "
              f"{'All>1':>6} {'All<1':>6} {'CV':>8}")

        band_stats = []
        for band in BANDS:
            bdf = mdf[mdf["band"] == band]
            ris = bdf["RI"].values
            mean_ri = np.nanmean(ris)
            std_ri = np.nanstd(ris)
            cv = std_ri / abs(mean_ri) if abs(mean_ri) > 1e-6 else np.inf
            all_above = np.all(ris > 1.0)
            all_below = np.all(ris < 1.0)
            print(f"{band:<12} {mean_ri:>8.3f} {std_ri:>8.3f} {np.min(ris):>8.3f} "
                  f"{np.max(ris):>8.3f} {str(all_above):>6} {str(all_below):>6} {cv:>8.3f}")
            band_stats.append({
                "band": band,
                "mean_ri": mean_ri,
                "std_ri": std_ri,
                "cv": cv,
                "all_agree_reorg": all_above,
                "all_agree_persist": all_below,
                "n_agree": sum(ris > 1) if mean_ri > 1 else sum(ris < 1),
            })

        bsdf = pd.DataFrame(band_stats)
        # Consistency score: how many bands have unanimous agreement?
        n_unanimous = bsdf["all_agree_reorg"].sum() + bsdf["all_agree_persist"].sum()
        # Also: mean CV across bands (lower = more consistent)
        mean_cv = bsdf["cv"].mean()
        # Combined: number of unanimous bands, tiebreak by low CV
        consistency_scores.append({
            "metric": metric_name,
            "n_unanimous_bands": n_unanimous,
            "mean_cv": mean_cv,
            "n_reorg_bands": bsdf["all_agree_reorg"].sum(),
            "n_persist_bands": bsdf["all_agree_persist"].sum(),
        })

    csdf = pd.DataFrame(consistency_scores).sort_values(
        ["n_unanimous_bands", "mean_cv"], ascending=[False, True]
    )

    print("\n" + "=" * 80)
    print("METRIC RANKING (by cross-patient consistency)")
    print("=" * 80)
    print(f"{'Metric':<14} {'Unanimous':>10} {'Reorg':>7} {'Persist':>8} {'Mean CV':>8}")
    print("-" * 50)
    for _, row in csdf.iterrows():
        print(f"{row['metric']:<14} {row['n_unanimous_bands']:>10.0f} "
              f"{row['n_reorg_bands']:>7.0f} {row['n_persist_bands']:>8.0f} "
              f"{row['mean_cv']:>8.3f}")

    # ===================================================================
    # Figure: RI heatmaps for top metrics
    # ===================================================================
    top_metrics = csdf["metric"].values[:6]  # Top 6

    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    axes = axes.ravel()

    for idx, metric_name in enumerate(top_metrics):
        ax = axes[idx]
        mdf = df[df["metric"] == metric_name]

        # Build patient x band RI matrix
        ri_matrix = np.zeros((len(PATIENTS), len(BANDS)))
        for i, pat in enumerate(PATIENTS):
            for j, band in enumerate(BANDS):
                row = mdf[(mdf["patient"] == pat) & (mdf["band"] == band)]
                ri_matrix[i, j] = row["RI"].values[0] if len(row) > 0 else np.nan

        im = ax.imshow(ri_matrix, cmap="RdBu_r", vmin=0.5, vmax=1.5, aspect="auto")
        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS],
                           rotation=45, ha="right", fontsize=10)
        ax.set_yticks(range(len(PATIENTS)))
        ax.set_yticklabels(PATIENTS, fontsize=10)
        ax.set_title(metric_name, fontsize=13, fontweight="bold")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        # Annotate
        for i in range(len(PATIENTS)):
            for j in range(len(BANDS)):
                val = ri_matrix[i, j]
                color = "white" if abs(val - 1.0) > 0.25 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=8, fontweight="bold", color=color)

    fig.suptitle(
        "Reorganization Index (RI) by Patient and Band\n"
        "RI > 1: reorganization (red)  |  RI < 1: persistence (blue)  |  RI = 1: no pattern",
        fontsize=14, fontweight="bold", y=0.98,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(OUTPUT_DIR / "metric_consistency_heatmap.pdf", dpi=150, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "metric_consistency_heatmap.png", dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved heatmap: {OUTPUT_DIR / 'metric_consistency_heatmap.png'}")

    # ===================================================================
    # Figure 2: Per-band RI with patient dots for best metric
    # ===================================================================
    best_metric = csdf.iloc[0]["metric"]
    mdf = df[df["metric"] == best_metric]

    fig, ax = plt.subplots(figsize=(12, 6))
    x_positions = np.arange(len(BANDS))
    colors_pat = plt.colormaps["Set2"](np.linspace(0, 1, len(PATIENTS)))

    for i, pat in enumerate(PATIENTS):
        ris = [mdf[(mdf["patient"] == pat) & (mdf["band"] == b)]["RI"].values[0]
               for b in BANDS]
        ax.scatter(x_positions + (i - 1.5) * 0.08, ris, s=80, color=colors_pat[i],
                   label=pat, zorder=5, edgecolors="black", linewidths=0.5)

    # Band means
    for j, band in enumerate(BANDS):
        vals = mdf[mdf["band"] == band]["RI"].values
        ax.plot([j - 0.2, j + 0.2], [np.mean(vals)] * 2, color="black",
                linewidth=3, zorder=6)

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1.5, alpha=0.7)
    ax.set_xticks(x_positions)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS], fontsize=12)
    ax.set_ylabel("Reorganization Index", fontsize=12)
    ax.set_title(f"Best metric: {best_metric}\n"
                 f"RI = D(rest↔task) / D(within-condition)", fontsize=14, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)

    fig.savefig(OUTPUT_DIR / "best_metric_ri_by_band.pdf", dpi=150, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "best_metric_ri_by_band.png", dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUTPUT_DIR / 'best_metric_ri_by_band.png'}")


if __name__ == "__main__":
    main()
