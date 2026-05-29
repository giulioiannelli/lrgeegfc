#!/usr/bin/env python3
"""Sweep of partition-based dendrogram comparison metrics.

For two dendrograms Z1 and Z2 (same patient, same band, different phases):
1. Cut each at k clusters using fcluster(Z, k, criterion='maxclust')
2. Compare the resulting partitions using ARI, NMI, and VI
3. Try k = 2, 3, 4, 5, 6, 8, 10, 15, 20
4. Also compute aggregates: mean, max, weighted (1/k)

Hypotheses tested:
  H1: sim(task_learn, task_test) highest among all pairs
  H2: sim(task_test, rest_post) > sim(rest_pre, rest_post)
  H3: within-type sim > cross-type sim
  H4: Reorganization strength decreases with frequency
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
)

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import LRG_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
from lrg_eegfc.config.const import PATIENTS_4PHASE
PATIENTS_4PH = PATIENTS_4PHASE
PATIENTS_ALL = list(PATIENT_PHASES.keys())
BANDS = BRAIN_BANDS_NAMES

K_VALUES = [2, 3, 4, 5, 6, 8, 10, 15, 20]

REST = {"rest_pre", "rest_post"}
TASK = {"task_learn", "task_test"}

# Phase pair labels for 4-phase patients (6 pairs)
PAIR_LABELS_4PH = [
    ("rest_pre", "task_learn"),
    ("rest_pre", "task_test"),
    ("rest_pre", "rest_post"),
    ("task_learn", "task_test"),
    ("task_learn", "rest_post"),
    ("task_test", "rest_post"),
]

PAIR_SHORT = {
    ("rest_pre", "task_learn"): "Pre-TL",
    ("rest_pre", "task_test"): "Pre-TT",
    ("rest_pre", "rest_post"): "Pre-Post",
    ("task_learn", "task_test"): "TL-TT",
    ("task_learn", "rest_post"): "TL-Post",
    ("task_test", "rest_post"): "TT-Post",
}

OUT_DIR = FIGURES_ROOT / "metric_exploration" / "partition_multiscale"


def classify_pair(p1: str, p2: str) -> str:
    s = {p1, p2}
    if s <= REST or s <= TASK:
        return "within"
    return "cross"


# ---------------------------------------------------------------------------
# Variation of Information
# ---------------------------------------------------------------------------
def variation_of_information(labels1: np.ndarray, labels2: np.ndarray) -> float:
    """Compute Variation of Information between two clusterings.

    VI(C1, C2) = H(C1) + H(C2) - 2*I(C1, C2)
    Lower = more similar. Range [0, log(n)].
    """
    n = len(labels1)
    if n == 0:
        return 0.0

    # Build contingency
    clusters1 = np.unique(labels1)
    clusters2 = np.unique(labels2)

    # Entropy of each partition
    def entropy(labels):
        _, counts = np.unique(labels, return_counts=True)
        p = counts / n
        return -np.sum(p * np.log(p + 1e-30))

    h1 = entropy(labels1)
    h2 = entropy(labels2)

    # Mutual information
    mi = 0.0
    for c1 in clusters1:
        mask1 = labels1 == c1
        n1 = mask1.sum()
        for c2 in clusters2:
            mask2 = labels2 == c2
            nij = (mask1 & mask2).sum()
            if nij > 0:
                mi += (nij / n) * np.log((nij * n) / (n1 * mask2.sum()) + 1e-30)

    vi = h1 + h2 - 2 * mi
    return max(vi, 0.0)


# ---------------------------------------------------------------------------
# Natural k from dendrogram
# ---------------------------------------------------------------------------
def find_natural_k(Z: np.ndarray) -> int:
    """Find natural number of clusters from largest gap in merge heights."""
    heights = Z[:, 2]
    if len(heights) < 2:
        return 2
    gaps = np.diff(heights)
    natural_k = len(heights) - np.argmax(gaps)
    return max(natural_k, 2)


# ---------------------------------------------------------------------------
# Load all LRG data
# ---------------------------------------------------------------------------
def load_all() -> dict:
    """Load linkage matrices for all (patient, phase, band)."""
    data = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for ph in phases:
                lrg = load_lrg_result(pat, ph, band, "msc", cache_root=LRG_CACHE)
                if lrg is None:
                    print(f"  WARNING: Missing {pat} {ph} {band}")
                    continue
                data[(pat, ph, band)] = {
                    "Z": lrg.linkage_matrix,
                    "n_nodes": lrg.n_nodes,
                }
    return data


# ---------------------------------------------------------------------------
# Compute partition metrics for a pair
# ---------------------------------------------------------------------------
def compute_pair_metrics(
    Z1: np.ndarray, Z2: np.ndarray, n1: int, n2: int
) -> dict:
    """Compute ARI, NMI, VI at each k, plus aggregates."""
    results = {}

    # Check node count mismatch
    if n1 != n2:
        # Use minimum - both dendrograms have indices 0..n-1
        # We restrict to the first min(n1,n2) nodes
        n = min(n1, n2)
        print(f"    WARNING: n_nodes mismatch ({n1} vs {n2}), using first {n} nodes")
    else:
        n = n1

    # Per-k metrics
    ari_vals = {}
    nmi_vals = {}
    vi_vals = {}

    for k in K_VALUES:
        if k >= n:
            continue  # Skip if k >= number of nodes

        labels1 = fcluster(Z1, k, criterion="maxclust")
        labels2 = fcluster(Z2, k, criterion="maxclust")

        # If node counts differ, restrict to shared nodes
        if n1 != n2:
            labels1 = labels1[:n]
            labels2 = labels2[:n]

        ari = adjusted_rand_score(labels1, labels2)
        nmi = normalized_mutual_info_score(labels1, labels2)
        vi = variation_of_information(labels1, labels2)

        ari_vals[k] = ari
        nmi_vals[k] = nmi
        vi_vals[k] = vi

        results[f"ARI_k{k}"] = ari
        results[f"NMI_k{k}"] = nmi
        results[f"VI_k{k}"] = vi

    # Aggregates
    if ari_vals:
        ks = sorted(ari_vals.keys())
        ari_arr = np.array([ari_vals[k] for k in ks])
        nmi_arr = np.array([nmi_vals[k] for k in ks])
        vi_arr = np.array([vi_vals[k] for k in ks])
        weights = np.array([1.0 / k for k in ks])
        weights /= weights.sum()

        results["mean_ARI"] = np.mean(ari_arr)
        results["max_ARI"] = np.max(ari_arr)
        results["weighted_ARI"] = np.average(ari_arr, weights=weights)

        results["mean_NMI"] = np.mean(nmi_arr)
        results["max_NMI"] = np.max(nmi_arr)
        results["weighted_NMI"] = np.average(nmi_arr, weights=weights)

        results["mean_VI"] = np.mean(vi_arr)
        results["min_VI"] = np.min(vi_arr)
        results["weighted_VI"] = np.average(vi_arr, weights=weights)

    # Natural k
    nat_k1 = find_natural_k(Z1)
    nat_k2 = find_natural_k(Z2)
    results["natural_k_1"] = nat_k1
    results["natural_k_2"] = nat_k2

    # ARI/NMI at the average natural k (rounded)
    avg_nat_k = int(round((nat_k1 + nat_k2) / 2))
    avg_nat_k = max(2, min(avg_nat_k, n - 1))

    labels1_nat = fcluster(Z1, avg_nat_k, criterion="maxclust")
    labels2_nat = fcluster(Z2, avg_nat_k, criterion="maxclust")
    if n1 != n2:
        labels1_nat = labels1_nat[:n]
        labels2_nat = labels2_nat[:n]

    results["ARI_natural_k"] = adjusted_rand_score(labels1_nat, labels2_nat)
    results["NMI_natural_k"] = normalized_mutual_info_score(labels1_nat, labels2_nat)
    results["VI_natural_k"] = variation_of_information(labels1_nat, labels2_nat)
    results["avg_natural_k"] = avg_nat_k

    return results


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------
def compute_all_pairs(data: dict) -> pd.DataFrame:
    """Compute partition metrics for all phase pairs."""
    rows = []
    for pat, phases in PATIENT_PHASES.items():
        avail = [ph for ph in phases if any((pat, ph, b) in data for b in BANDS)]
        for band in BANDS:
            avail_band = [ph for ph in phases if (pat, ph, band) in data]
            for p1, p2 in combinations(avail_band, 2):
                d1 = data[(pat, p1, band)]
                d2 = data[(pat, p2, band)]
                metrics = compute_pair_metrics(
                    d1["Z"], d2["Z"], d1["n_nodes"], d2["n_nodes"]
                )
                row = {
                    "patient": pat,
                    "band": band,
                    "phase1": p1,
                    "phase2": p2,
                    "pair": f"{p1}-{p2}",
                    "pair_type": classify_pair(p1, p2),
                }
                row.update(metrics)
                rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Hypothesis testing
# ---------------------------------------------------------------------------
def test_hypotheses(df: pd.DataFrame) -> dict:
    """Test H1-H4 for each metric variant.

    Returns dict: metric_name -> {band -> {hypothesis -> {result, details}}}
    """
    # Identify all metric columns (ARI_k*, NMI_k*, VI_k*, aggregates, natural)
    metric_cols = [c for c in df.columns if c.startswith(("ARI_", "NMI_", "VI_",
                                                          "mean_", "max_", "min_",
                                                          "weighted_"))]
    # Remove non-metric columns
    metric_cols = [c for c in metric_cols if c not in
                   ("natural_k_1", "natural_k_2", "avg_natural_k")]

    results = {}

    for metric in metric_cols:
        # Determine if higher = more similar (ARI, NMI) or lower = more similar (VI)
        is_similarity = not metric.startswith(("VI_", "min_VI", "mean_VI", "weighted_VI"))

        band_results = {}
        for band in BANDS:
            band_df = df[(df["band"] == band)].copy()
            # Only use 4-phase patients for hypothesis testing
            band_4ph = band_df[band_df["patient"].isin(PATIENTS_4PH)]

            hyp = {}

            # --- H1: sim(task_learn, task_test) highest ---
            h1_pass_count = 0
            h1_total = 0
            for pat in PATIENTS_4PH:
                pat_df = band_4ph[band_4ph["patient"] == pat]
                if pat_df.empty:
                    continue
                tl_tt = pat_df[
                    (pat_df["phase1"] == "task_learn") & (pat_df["phase2"] == "task_test")
                ][metric].values
                if len(tl_tt) == 0:
                    continue
                tl_tt_val = tl_tt[0]
                other_vals = pat_df[
                    ~((pat_df["phase1"] == "task_learn") & (pat_df["phase2"] == "task_test"))
                ][metric].values

                if len(other_vals) == 0:
                    continue
                h1_total += 1
                if is_similarity:
                    if tl_tt_val >= np.max(other_vals):
                        h1_pass_count += 1
                else:
                    if tl_tt_val <= np.min(other_vals):
                        h1_pass_count += 1

            hyp["H1_TL_TT_highest"] = {
                "pass": h1_pass_count,
                "total": h1_total,
                "unanimous": h1_pass_count == h1_total and h1_total > 0,
                "relaxed": h1_pass_count >= h1_total - 1 and h1_total > 0,
            }

            # --- H2: sim(task_test, rest_post) > sim(rest_pre, rest_post) ---
            h2_pass_count = 0
            h2_total = 0
            for pat in PATIENTS_4PH:
                pat_df = band_4ph[band_4ph["patient"] == pat]
                tt_post = pat_df[
                    (pat_df["phase1"] == "task_test") & (pat_df["phase2"] == "rest_post")
                ][metric].values
                pre_post = pat_df[
                    (pat_df["phase1"] == "rest_pre") & (pat_df["phase2"] == "rest_post")
                ][metric].values
                if len(tt_post) == 0 or len(pre_post) == 0:
                    continue
                h2_total += 1
                if is_similarity:
                    if tt_post[0] > pre_post[0]:
                        h2_pass_count += 1
                else:
                    if tt_post[0] < pre_post[0]:
                        h2_pass_count += 1

            hyp["H2_TT_Post_gt_Pre_Post"] = {
                "pass": h2_pass_count,
                "total": h2_total,
                "unanimous": h2_pass_count == h2_total and h2_total > 0,
                "relaxed": h2_pass_count >= h2_total - 1 and h2_total > 0,
            }

            # --- H3: within-type sim > cross-type sim ---
            h3_pass_count = 0
            h3_total = 0
            for pat in PATIENTS_4PH:
                pat_df = band_4ph[band_4ph["patient"] == pat]
                within = pat_df[pat_df["pair_type"] == "within"][metric].values
                cross = pat_df[pat_df["pair_type"] == "cross"][metric].values
                if len(within) == 0 or len(cross) == 0:
                    continue
                h3_total += 1
                if is_similarity:
                    if np.mean(within) > np.mean(cross):
                        h3_pass_count += 1
                else:
                    if np.mean(within) < np.mean(cross):
                        h3_pass_count += 1

            hyp["H3_within_gt_cross"] = {
                "pass": h3_pass_count,
                "total": h3_total,
                "unanimous": h3_pass_count == h3_total and h3_total > 0,
                "relaxed": h3_pass_count >= h3_total - 1 and h3_total > 0,
            }

            band_results[band] = hyp

        # --- H4: Reorganization strength decreases with frequency ---
        # Compute reorg strength per band (mean within - mean cross) averaged across patients
        reorg_by_band = []
        for band in BANDS:
            band_4ph = df[(df["band"] == band) & (df["patient"].isin(PATIENTS_4PH))]
            within_vals = band_4ph[band_4ph["pair_type"] == "within"][metric].values
            cross_vals = band_4ph[band_4ph["pair_type"] == "cross"][metric].values
            if len(within_vals) > 0 and len(cross_vals) > 0:
                if is_similarity:
                    gap = np.mean(within_vals) - np.mean(cross_vals)
                else:
                    gap = np.mean(cross_vals) - np.mean(within_vals)
                reorg_by_band.append(gap)
            else:
                reorg_by_band.append(np.nan)

        # Check if monotonically decreasing (or mostly decreasing)
        valid_gaps = [(i, g) for i, g in enumerate(reorg_by_band) if not np.isnan(g)]
        if len(valid_gaps) >= 4:
            gaps_arr = [g for _, g in valid_gaps]
            # Count decreasing steps
            n_decreasing = sum(
                gaps_arr[i] > gaps_arr[i + 1] for i in range(len(gaps_arr) - 1)
            )
            n_steps = len(gaps_arr) - 1
            h4_monotone = n_decreasing == n_steps
            h4_relaxed = n_decreasing >= n_steps - 1
        else:
            h4_monotone = False
            h4_relaxed = False

        # Add H4 to results (band-independent)
        for band in BANDS:
            if band in band_results:
                band_results[band]["H4_freq_decrease"] = {
                    "monotone": h4_monotone,
                    "relaxed": h4_relaxed,
                    "reorg_by_band": reorg_by_band,
                }

        results[metric] = band_results

    return results


# ---------------------------------------------------------------------------
# Z-score normalization within (patient, band)
# ---------------------------------------------------------------------------
def zscore_normalize(df: pd.DataFrame, metric_col: str) -> pd.DataFrame:
    """Add z-scored column for a metric within each (patient, band)."""
    col_z = f"{metric_col}_z"
    df[col_z] = np.nan
    for (pat, band), grp in df.groupby(["patient", "band"]):
        vals = grp[metric_col].values
        if len(vals) > 1 and np.std(vals) > 1e-10:
            z = (vals - np.mean(vals)) / np.std(vals)
            df.loc[grp.index, col_z] = z
        else:
            df.loc[grp.index, col_z] = 0.0
    return df


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_summary_heatmap(hyp_results: dict, out_path: Path):
    """Heatmap: metric x band, showing hypothesis pass counts."""
    metrics = sorted(hyp_results.keys())
    hypotheses = ["H1_TL_TT_highest", "H2_TT_Post_gt_Pre_Post", "H3_within_gt_cross"]

    fig, axes = plt.subplots(1, 3, figsize=(20, max(8, len(metrics) * 0.35)))

    for ax_idx, hyp_name in enumerate(hypotheses):
        ax = axes[ax_idx]
        matrix = np.zeros((len(metrics), len(BANDS)))

        for i, metric in enumerate(metrics):
            for j, band in enumerate(BANDS):
                if band in hyp_results[metric]:
                    h = hyp_results[metric][band].get(hyp_name, {})
                    total = h.get("total", 0)
                    passed = h.get("pass", 0)
                    if total > 0:
                        matrix[i, j] = passed / total
                    else:
                        matrix[i, j] = np.nan

        im = ax.imshow(matrix, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS], fontsize=8)
        ax.set_yticks(range(len(metrics)))
        ax.set_yticklabels(metrics, fontsize=7)
        ax.set_title(hyp_name.replace("_", " "), fontsize=10, fontweight="bold")

        # Annotate cells
        for i in range(len(metrics)):
            for j in range(len(BANDS)):
                val = matrix[i, j]
                if not np.isnan(val):
                    color = "white" if val < 0.3 or val > 0.7 else "black"
                    ax.text(j, i, f"{val:.0%}", ha="center", va="center",
                            fontsize=6, color=color)

    fig.colorbar(im, ax=axes, shrink=0.6, label="Fraction of patients passing")
    fig.suptitle("Partition Metric Hypothesis Testing\n(fraction of 4-phase patients passing)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


def plot_ari_profiles(df: pd.DataFrame, best_metric: str, out_path: Path):
    """Per-pair profiles across bands for the best metric (like logCosine profiles)."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for idx, pat in enumerate(PATIENTS_4PH):
        ax = axes[idx]
        pat_df = df[df["patient"] == pat]

        for p1, p2 in PAIR_LABELS_4PH:
            pair_df = pat_df[
                (pat_df["phase1"] == p1) & (pat_df["phase2"] == p2)
            ]
            if pair_df.empty:
                continue
            vals = []
            for band in BANDS:
                band_row = pair_df[pair_df["band"] == band]
                if not band_row.empty and best_metric in band_row.columns:
                    vals.append(band_row[best_metric].values[0])
                else:
                    vals.append(np.nan)

            pair_type = classify_pair(p1, p2)
            ls = "-" if pair_type == "within" else "--"
            lw = 2.0 if pair_type == "within" else 1.5
            label = PAIR_SHORT.get((p1, p2), f"{p1}-{p2}")
            ax.plot(range(len(BANDS)), vals, marker="o", ls=ls, lw=lw, label=label)

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS], fontsize=9)
        ax.set_title(pat, fontsize=11, fontweight="bold")
        ax.set_ylabel(best_metric, fontsize=9)
        ax.legend(fontsize=7, ncol=2)
        ax.grid(True, alpha=0.3)

    fig.suptitle(f"Partition Metric Profiles: {best_metric}", fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


def plot_scale_comparison(df: pd.DataFrame, out_path: Path):
    """How ARI changes with k for each phase pair (averaged across 4-phase patients)."""
    ari_cols = [f"ARI_k{k}" for k in K_VALUES if f"ARI_k{k}" in df.columns]
    k_avail = [k for k in K_VALUES if f"ARI_k{k}" in df.columns]

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    for idx, band in enumerate(BANDS):
        ax = axes[idx // 3, idx % 3]
        band_df = df[(df["band"] == band) & (df["patient"].isin(PATIENTS_4PH))]

        for p1, p2 in PAIR_LABELS_4PH:
            pair_data = band_df[
                (band_df["phase1"] == p1) & (band_df["phase2"] == p2)
            ]
            if pair_data.empty:
                continue

            # Average across patients
            means = []
            sems = []
            for k in k_avail:
                col = f"ARI_k{k}"
                vals = pair_data[col].dropna().values
                means.append(np.mean(vals) if len(vals) > 0 else np.nan)
                sems.append(np.std(vals) / np.sqrt(len(vals)) if len(vals) > 1 else 0)

            pair_type = classify_pair(p1, p2)
            ls = "-" if pair_type == "within" else "--"
            lw = 2.0 if pair_type == "within" else 1.5
            label = PAIR_SHORT.get((p1, p2), f"{p1}-{p2}")

            ax.errorbar(k_avail, means, yerr=sems, marker="o", ls=ls, lw=lw,
                        label=label, capsize=3, markersize=4)

        ax.set_xlabel("k (number of clusters)", fontsize=9)
        ax.set_ylabel("ARI", fontsize=9)
        ax.set_title(f"{BRAIN_BAND_TEX_DICT.get(band, band)}", fontsize=11, fontweight="bold")
        ax.legend(fontsize=6, ncol=2)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color="gray", lw=0.5, ls=":")
        ax.set_xscale("log")
        ax.set_xticks(k_avail)
        ax.set_xticklabels([str(k) for k in k_avail], fontsize=7)

    fig.suptitle("ARI vs number of clusters k\n(mean +/- SEM across 4-phase patients)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


def plot_natural_k_analysis(df: pd.DataFrame, data: dict, out_path: Path):
    """Analyze natural k distribution across (patient, band, phase)."""
    # Compute natural_k for every (patient, phase, band)
    nat_k_rows = []
    for (pat, ph, band), d in data.items():
        nk = find_natural_k(d["Z"])
        nat_k_rows.append({
            "patient": pat,
            "phase": ph,
            "band": band,
            "natural_k": nk,
            "n_nodes": d["n_nodes"],
        })
    nat_k_df = pd.DataFrame(nat_k_rows)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel 1: Distribution of natural_k
    ax = axes[0, 0]
    all_nk = nat_k_df["natural_k"].values
    ax.hist(all_nk, bins=range(2, max(all_nk) + 2), edgecolor="black", alpha=0.7)
    ax.set_xlabel("Natural k", fontsize=10)
    ax.set_ylabel("Count", fontsize=10)
    ax.set_title("Distribution of natural k across all conditions", fontsize=10, fontweight="bold")
    ax.axvline(np.median(all_nk), color="red", ls="--", lw=1.5, label=f"Median={np.median(all_nk):.0f}")
    ax.axvline(np.mean(all_nk), color="blue", ls="--", lw=1.5, label=f"Mean={np.mean(all_nk):.1f}")
    ax.legend(fontsize=8)

    # Panel 2: Natural k by band
    ax = axes[0, 1]
    band_data = []
    for band in BANDS:
        vals = nat_k_df[nat_k_df["band"] == band]["natural_k"].values
        band_data.append(vals)
    bp = ax.boxplot(band_data, labels=[BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS],
                    patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("lightblue")
    ax.set_ylabel("Natural k", fontsize=10)
    ax.set_title("Natural k by frequency band", fontsize=10, fontweight="bold")

    # Panel 3: Natural k by phase
    ax = axes[1, 0]
    phase_order = ["rest_pre", "task_learn", "task_test", "rest_post"]
    phase_data = []
    for ph in phase_order:
        vals = nat_k_df[nat_k_df["phase"] == ph]["natural_k"].values
        phase_data.append(vals)
    bp = ax.boxplot(phase_data, labels=phase_order, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("lightyellow")
    ax.set_ylabel("Natural k", fontsize=10)
    ax.set_title("Natural k by phase", fontsize=10, fontweight="bold")

    # Panel 4: Natural k by patient
    ax = axes[1, 1]
    pat_data = []
    for pat in PATIENTS_ALL:
        vals = nat_k_df[nat_k_df["patient"] == pat]["natural_k"].values
        pat_data.append(vals)
    bp = ax.boxplot(pat_data, labels=PATIENTS_ALL, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("lightgreen")
    ax.set_ylabel("Natural k", fontsize=10)
    ax.set_title("Natural k by patient", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", rotation=45)

    fig.suptitle("Natural k Analysis (largest gap in merge heights)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")

    return nat_k_df


# ---------------------------------------------------------------------------
# Consistency scoring (for identifying best metric)
# ---------------------------------------------------------------------------
def score_metrics_consistency(df: pd.DataFrame) -> pd.DataFrame:
    """Score each metric by sign consistency of the within-cross gap.

    For each (metric, patient, band):
      gap = mean(within_pairs) - mean(cross_pairs)
    Check if gap sign is consistent across patients.
    """
    metric_cols = [c for c in df.columns if c.startswith(("ARI_", "NMI_", "VI_",
                                                          "mean_", "max_", "min_",
                                                          "weighted_"))]
    metric_cols = [c for c in metric_cols if c not in
                   ("natural_k_1", "natural_k_2", "avg_natural_k")]

    score_rows = []
    for metric in metric_cols:
        is_similarity = not metric.startswith(("VI_", "min_VI", "mean_VI", "weighted_VI"))

        band_unanimous = 0
        band_relaxed = 0
        total_bands = 0

        for band in BANDS:
            pat_gaps = []
            for pat in PATIENTS_4PH:
                pat_band = df[(df["patient"] == pat) & (df["band"] == band)]
                within = pat_band[pat_band["pair_type"] == "within"][metric].dropna().values
                cross = pat_band[pat_band["pair_type"] == "cross"][metric].dropna().values
                if len(within) > 0 and len(cross) > 0:
                    if is_similarity:
                        gap = np.mean(within) - np.mean(cross)
                    else:
                        gap = np.mean(cross) - np.mean(within)  # Higher VI = less similar
                    pat_gaps.append(gap)

            if len(pat_gaps) >= 3:
                total_bands += 1
                n_pos = sum(g > 0 for g in pat_gaps)
                n_neg = sum(g < 0 for g in pat_gaps)
                if n_pos == len(pat_gaps) or n_neg == len(pat_gaps):
                    band_unanimous += 1
                if n_pos >= len(pat_gaps) - 1 or n_neg >= len(pat_gaps) - 1:
                    band_relaxed += 1

        score_rows.append({
            "metric": metric,
            "unanimous_bands": band_unanimous,
            "relaxed_bands": band_relaxed,
            "total_bands": total_bands,
            "is_similarity": is_similarity,
        })

    score_df = pd.DataFrame(score_rows)
    score_df = score_df.sort_values(
        ["unanimous_bands", "relaxed_bands"], ascending=[False, False]
    ).reset_index(drop=True)
    return score_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("PARTITION-BASED DENDROGRAM COMPARISON SWEEP")
    print("=" * 70)

    print("\n[1/6] Loading LRG data...")
    data = load_all()
    print(f"  Loaded {len(data)} entries")

    print("\n[2/6] Computing partition metrics for all pairs...")
    df = compute_all_pairs(data)
    print(f"  Computed {len(df)} pair comparisons")

    # Save raw results
    csv_path = OUT_DIR / "results.csv"
    df.to_csv(csv_path, index=False)
    print(f"  Saved raw results to {csv_path}")

    print("\n[3/6] Testing hypotheses...")
    hyp_results = test_hypotheses(df)

    # Print hypothesis summary
    print("\n  --- Hypothesis Summary (per-metric, across bands) ---")
    metric_hyp_scores = []
    for metric, band_hyps in hyp_results.items():
        h1_u = sum(1 for b in BANDS if b in band_hyps
                   and band_hyps[b].get("H1_TL_TT_highest", {}).get("unanimous", False))
        h1_r = sum(1 for b in BANDS if b in band_hyps
                   and band_hyps[b].get("H1_TL_TT_highest", {}).get("relaxed", False))
        h2_u = sum(1 for b in BANDS if b in band_hyps
                   and band_hyps[b].get("H2_TT_Post_gt_Pre_Post", {}).get("unanimous", False))
        h2_r = sum(1 for b in BANDS if b in band_hyps
                   and band_hyps[b].get("H2_TT_Post_gt_Pre_Post", {}).get("relaxed", False))
        h3_u = sum(1 for b in BANDS if b in band_hyps
                   and band_hyps[b].get("H3_within_gt_cross", {}).get("unanimous", False))
        h3_r = sum(1 for b in BANDS if b in band_hyps
                   and band_hyps[b].get("H3_within_gt_cross", {}).get("relaxed", False))
        total_u = h1_u + h2_u + h3_u
        total_r = h1_r + h2_r + h3_r

        metric_hyp_scores.append({
            "metric": metric,
            "H1_unan": h1_u, "H1_relax": h1_r,
            "H2_unan": h2_u, "H2_relax": h2_r,
            "H3_unan": h3_u, "H3_relax": h3_r,
            "total_unan": total_u, "total_relax": total_r,
        })

    hyp_df = pd.DataFrame(metric_hyp_scores).sort_values(
        ["total_unan", "total_relax", "H3_unan"], ascending=[False, False, False]
    ).reset_index(drop=True)
    print(hyp_df.to_string(index=False))

    print("\n[4/6] Scoring metrics by within-cross gap consistency...")
    consistency_df = score_metrics_consistency(df)
    print(consistency_df.head(20).to_string(index=False))

    # Identify best metric
    best_metric = consistency_df.iloc[0]["metric"]
    print(f"\n  BEST metric by gap consistency: {best_metric}")

    # Also consider hypothesis-based best
    if not hyp_df.empty:
        hyp_best = hyp_df.iloc[0]["metric"]
        print(f"  BEST metric by hypothesis score: {hyp_best}")

    print("\n[5/6] Generating figures...")

    # Figure 1: Summary heatmap
    plot_summary_heatmap(hyp_results, OUT_DIR / "summary_heatmap.pdf")

    # Figure 2: ARI profiles for best metric
    plot_ari_profiles(df, best_metric, OUT_DIR / "ari_profiles.pdf")

    # Also plot for hyp_best if different
    if not hyp_df.empty and hyp_best != best_metric:
        plot_ari_profiles(df, hyp_best, OUT_DIR / f"ari_profiles_{hyp_best}.pdf")

    # Figure 3: Scale comparison (ARI vs k)
    plot_scale_comparison(df, OUT_DIR / "scale_comparison.pdf")

    # Figure 4: Natural k analysis
    nat_k_df = plot_natural_k_analysis(df, data, OUT_DIR / "natural_k_analysis.pdf")

    # Save natural_k data
    nat_k_df.to_csv(OUT_DIR / "natural_k.csv", index=False)

    # Save hypothesis scores
    hyp_df.to_csv(OUT_DIR / "hypothesis_scores.csv", index=False)

    # Save consistency scores
    consistency_df.to_csv(OUT_DIR / "consistency_scores.csv", index=False)

    print("\n[6/6] Detailed results for best metrics...")

    # Print per-pair, per-band values for the best metric
    print(f"\n  --- {best_metric}: Per-pair mean across 4-phase patients ---")
    print(f"  {'Pair':<12}", end="")
    for band in BANDS:
        print(f"  {BRAIN_BAND_TEX_DICT.get(band, band):>10}", end="")
    print()

    for p1, p2 in PAIR_LABELS_4PH:
        label = PAIR_SHORT.get((p1, p2), f"{p1}-{p2}")
        print(f"  {label:<12}", end="")
        for band in BANDS:
            pair_band = df[
                (df["patient"].isin(PATIENTS_4PH))
                & (df["band"] == band)
                & (df["phase1"] == p1)
                & (df["phase2"] == p2)
            ]
            if not pair_band.empty and best_metric in pair_band.columns:
                val = pair_band[best_metric].mean()
                print(f"  {val:>10.4f}", end="")
            else:
                print(f"  {'N/A':>10}", end="")
        print()

    # Print within-cross gap for the best metric
    print(f"\n  --- {best_metric}: Within-cross gap per (patient, band) ---")
    is_sim = not best_metric.startswith(("VI_", "min_VI", "mean_VI", "weighted_VI"))
    print(f"  {'Band':<12}", end="")
    for pat in PATIENTS_4PH:
        print(f"  {pat:>10}", end="")
    print(f"  {'Mean':>10}  {'Sign':>6}")

    for band in BANDS:
        print(f"  {band:<12}", end="")
        gaps = []
        for pat in PATIENTS_4PH:
            pb = df[(df["patient"] == pat) & (df["band"] == band)]
            w = pb[pb["pair_type"] == "within"][best_metric].dropna().values
            c = pb[pb["pair_type"] == "cross"][best_metric].dropna().values
            if len(w) > 0 and len(c) > 0:
                gap = np.mean(w) - np.mean(c) if is_sim else np.mean(c) - np.mean(w)
                gaps.append(gap)
                print(f"  {gap:>+10.4f}", end="")
            else:
                print(f"  {'N/A':>10}", end="")
        if gaps:
            mg = np.mean(gaps)
            sign = "ALL+" if all(g > 0 for g in gaps) else (
                "ALL-" if all(g < 0 for g in gaps) else "MIXED")
            print(f"  {mg:>+10.4f}  {sign:>6}")
        else:
            print()

    # Additional: Print z-scored version
    print(f"\n  --- Z-scored {best_metric}: Mean per pair across patients ---")
    df = zscore_normalize(df, best_metric)
    zcol = f"{best_metric}_z"
    print(f"  {'Pair':<12}  {'Type':<6}", end="")
    for band in BANDS:
        print(f"  {band:>10}", end="")
    print()

    for p1, p2 in PAIR_LABELS_4PH:
        label = PAIR_SHORT.get((p1, p2), f"{p1}-{p2}")
        ptype = classify_pair(p1, p2)[:3]
        print(f"  {label:<12}  {ptype:<6}", end="")
        for band in BANDS:
            pair_band = df[
                (df["patient"].isin(PATIENTS_4PH))
                & (df["band"] == band)
                & (df["phase1"] == p1)
                & (df["phase2"] == p2)
            ]
            if not pair_band.empty and zcol in pair_band.columns:
                val = pair_band[zcol].mean()
                print(f"  {val:>+10.3f}", end="")
            else:
                print(f"  {'N/A':>10}", end="")
        print()

    print("\n" + "=" * 70)
    print("DONE. All outputs in:", OUT_DIR)
    print("=" * 70)


if __name__ == "__main__":
    main()
