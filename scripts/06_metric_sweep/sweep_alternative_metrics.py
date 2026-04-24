#!/usr/bin/env python3
"""Sweep of alternative dendrogram comparison metrics for H2 (task trace).

H2: Does rest_post resemble task_test more than rest_pre?
  - H2a: dist(Pre,Post) - dist(TT,Post) > 0  (distance metrics)
          sim(TT,Post) - sim(Pre,Post) > 0     (similarity metrics)
  - H2b: dist(Pre,TT) - dist(TT,Post) > 0     (distance metrics)
          sim(TT,Post) - sim(Pre,TT) > 0       (similarity metrics)

For each metric and band, count unanimity across 4 patients.

Metrics implemented:
  1. Cophenetic Pearson correlation
  2. Cophenetic Spearman correlation
  3. Baker's Gamma
  4. Normalized L1 on ultrametric
  5. Normalized L2 on ultrametric
  6. Top-k merge agreement (k=5, 10, 20)
  7. Weighted partition ARI (coarse-weighted and fine-weighted)
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, cophenet
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import adjusted_rand_score

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import LRG_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PATIENTS = PATIENTS_4PHASE
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
BANDS = BRAIN_BANDS_NAMES
OUT_DIR = FIGURES_ROOT / "metric_exploration" / "alternative_metrics"

# Phase pairs for H2 contrasts
PAIR_LABELS = [
    ("rest_pre", "task_learn"),
    ("rest_pre", "task_test"),
    ("rest_pre", "rest_post"),
    ("task_learn", "task_test"),
    ("task_learn", "rest_post"),
    ("task_test", "rest_post"),
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_all() -> dict:
    """Load linkage + ultrametric for all (patient, phase, band)."""
    data = {}
    for pat in PATIENTS:
        for band in BANDS:
            for ph in PHASES:
                lrg = load_lrg_result(pat, ph, band, "msc", cache_root=LRG_CACHE)
                if lrg is None:
                    print(f"  WARNING: Missing {pat} {ph} {band}")
                    continue
                Z = lrg.linkage_matrix
                U = lrg.ultrametric_matrix  # condensed 1D
                c_condensed = cophenet(Z)   # cophenetic distances, condensed 1D
                data[(pat, ph, band)] = {
                    "Z": Z,
                    "U": U,
                    "c_cond": c_condensed,
                    "n_nodes": lrg.n_nodes,
                }
    return data


# ---------------------------------------------------------------------------
# Metric 1: Cophenetic Pearson correlation (similarity, higher = more similar)
# ---------------------------------------------------------------------------
def cophenetic_pearson(d1, d2):
    """Pearson correlation between condensed cophenetic vectors."""
    r, _ = pearsonr(d1["c_cond"], d2["c_cond"])
    return r


# ---------------------------------------------------------------------------
# Metric 2: Cophenetic Spearman correlation (similarity)
# ---------------------------------------------------------------------------
def cophenetic_spearman(d1, d2):
    """Spearman rank correlation of cophenetic vectors."""
    r, _ = spearmanr(d1["c_cond"], d2["c_cond"])
    return r


# ---------------------------------------------------------------------------
# Metric 3: Baker's Gamma (similarity)
# ---------------------------------------------------------------------------
def bakers_gamma(d1, d2):
    """Baker's Gamma: Spearman rank correlation of cophenetic distances.

    For each pair (i,j), both dendrograms assign a cophenetic height.
    Baker's Gamma is the Spearman correlation of these heights.
    Equivalent to cophenetic Spearman when using ultrametric vectors
    (which already are the cophenetic heights from the linkage).
    We use the ultrametric vectors directly for a cleaner implementation.
    """
    r, _ = spearmanr(d1["U"], d2["U"])
    return r


# ---------------------------------------------------------------------------
# Metric 4: Normalized L1 on ultrametric (distance, lower = more similar)
# ---------------------------------------------------------------------------
def norm_l1_ultrametric(d1, d2):
    """L1 distance between min-max normalized ultrametric vectors. Range [0,1]."""
    u1 = d1["U"].copy()
    u2 = d2["U"].copy()
    # Min-max normalize to [0,1]
    r1 = u1.max() - u1.min()
    r2 = u2.max() - u2.min()
    if r1 > 0:
        u1 = (u1 - u1.min()) / r1
    if r2 > 0:
        u2 = (u2 - u2.min()) / r2
    return np.mean(np.abs(u1 - u2))


# ---------------------------------------------------------------------------
# Metric 5: Normalized L2 on ultrametric (distance, lower = more similar)
# ---------------------------------------------------------------------------
def norm_l2_ultrametric(d1, d2):
    """L2 distance between min-max normalized ultrametric vectors."""
    u1 = d1["U"].copy()
    u2 = d2["U"].copy()
    r1 = u1.max() - u1.min()
    r2 = u2.max() - u2.min()
    if r1 > 0:
        u1 = (u1 - u1.min()) / r1
    if r2 > 0:
        u2 = (u2 - u2.min()) / r2
    return np.sqrt(np.mean((u1 - u2) ** 2))


# ---------------------------------------------------------------------------
# Metric 6: Top-k merge agreement (similarity, higher = more similar)
# ---------------------------------------------------------------------------
def _top_k_merge_agreement(Z1, Z2, k):
    """Fraction of top-k merges that agree between two dendrograms.

    A "merge" is defined by the pair of cluster indices being merged.
    We look at the last k merges (the coarsest structural decisions).
    Agreement means both dendrograms merge the same set of leaf nodes
    at those top levels.
    """
    n1 = Z1.shape[0] + 1
    n2 = Z2.shape[0] + 1

    def get_top_k_merge_sets(Z, n, k_val):
        """For the top k_val merges, get the set of leaves in each merged pair."""
        # Build membership for each node
        membership = {i: frozenset([i]) for i in range(n)}
        for idx, row in enumerate(Z):
            left, right = int(row[0]), int(row[1])
            new_id = n + idx
            membership[new_id] = membership[left] | membership[right]

        # Top k merges = last k rows
        merge_pairs = set()
        for row in Z[-k_val:]:
            left, right = int(row[0]), int(row[1])
            pair = frozenset([membership[left], membership[right]])
            merge_pairs.add(pair)
        return merge_pairs

    actual_k = min(k, Z1.shape[0], Z2.shape[0])
    if actual_k == 0:
        return 0.0

    merges1 = get_top_k_merge_sets(Z1, n1, actual_k)
    merges2 = get_top_k_merge_sets(Z2, n2, actual_k)

    overlap = len(merges1 & merges2)
    return overlap / actual_k


def top_k5_merge(d1, d2):
    return _top_k_merge_agreement(d1["Z"], d2["Z"], 5)


def top_k10_merge(d1, d2):
    return _top_k_merge_agreement(d1["Z"], d2["Z"], 10)


def top_k20_merge(d1, d2):
    return _top_k_merge_agreement(d1["Z"], d2["Z"], 20)


# ---------------------------------------------------------------------------
# Metric 7: Weighted partition ARI (similarity, higher = more similar)
# ---------------------------------------------------------------------------
def _variation_of_information(labels1, labels2):
    """Compute VI between two clusterings."""
    n = len(labels1)
    if n == 0:
        return 0.0

    def ent(labels):
        _, counts = np.unique(labels, return_counts=True)
        p = counts / n
        return -np.sum(p * np.log(p + 1e-30))

    h1, h2 = ent(labels1), ent(labels2)
    clusters1 = np.unique(labels1)
    clusters2 = np.unique(labels2)
    mi = 0.0
    for c1 in clusters1:
        m1 = labels1 == c1
        n1 = m1.sum()
        for c2 in clusters2:
            m2 = labels2 == c2
            nij = (m1 & m2).sum()
            if nij > 0:
                mi += (nij / n) * np.log((nij * n) / (n1 * m2.sum()) + 1e-30)
    return max(h1 + h2 - 2 * mi, 0.0)


K_VALUES = [2, 3, 4, 5, 6, 8, 10, 15, 20]


def weighted_ari_coarse(d1, d2):
    """Weighted ARI emphasizing coarse structure (weight = 1/k)."""
    Z1, Z2 = d1["Z"], d2["Z"]
    n = min(d1["n_nodes"], d2["n_nodes"])
    weights, ari_vals = [], []
    for k in K_VALUES:
        if k >= n:
            continue
        l1 = fcluster(Z1, k, criterion="maxclust")[:n]
        l2 = fcluster(Z2, k, criterion="maxclust")[:n]
        ari_vals.append(adjusted_rand_score(l1, l2))
        weights.append(1.0 / k)
    if not ari_vals:
        return None
    weights = np.array(weights)
    weights /= weights.sum()
    return np.average(ari_vals, weights=weights)


def weighted_ari_fine(d1, d2):
    """Weighted ARI emphasizing fine structure (weight = k)."""
    Z1, Z2 = d1["Z"], d2["Z"]
    n = min(d1["n_nodes"], d2["n_nodes"])
    weights, ari_vals = [], []
    for k in K_VALUES:
        if k >= n:
            continue
        l1 = fcluster(Z1, k, criterion="maxclust")[:n]
        l2 = fcluster(Z2, k, criterion="maxclust")[:n]
        ari_vals.append(adjusted_rand_score(l1, l2))
        weights.append(float(k))
    if not ari_vals:
        return None
    weights = np.array(weights)
    weights /= weights.sum()
    return np.average(ari_vals, weights=weights)


def weighted_vi_coarse(d1, d2):
    """Weighted VI emphasizing coarse structure (weight = 1/k). Distance metric."""
    Z1, Z2 = d1["Z"], d2["Z"]
    n = min(d1["n_nodes"], d2["n_nodes"])
    weights, vi_vals = [], []
    for k in K_VALUES:
        if k >= n:
            continue
        l1 = fcluster(Z1, k, criterion="maxclust")[:n]
        l2 = fcluster(Z2, k, criterion="maxclust")[:n]
        vi_vals.append(_variation_of_information(l1, l2))
        weights.append(1.0 / k)
    if not vi_vals:
        return None
    weights = np.array(weights)
    weights /= weights.sum()
    return np.average(vi_vals, weights=weights)


def weighted_vi_fine(d1, d2):
    """Weighted VI emphasizing fine structure (weight = k). Distance metric."""
    Z1, Z2 = d1["Z"], d2["Z"]
    n = min(d1["n_nodes"], d2["n_nodes"])
    weights, vi_vals = [], []
    for k in K_VALUES:
        if k >= n:
            continue
        l1 = fcluster(Z1, k, criterion="maxclust")[:n]
        l2 = fcluster(Z2, k, criterion="maxclust")[:n]
        vi_vals.append(_variation_of_information(l1, l2))
        weights.append(float(k))
    if not vi_vals:
        return None
    weights = np.array(weights)
    weights /= weights.sum()
    return np.average(vi_vals, weights=weights)


def mean_ari(d1, d2):
    """Unweighted mean ARI across k values. (similarity)."""
    Z1, Z2 = d1["Z"], d2["Z"]
    n = min(d1["n_nodes"], d2["n_nodes"])
    ari_vals = []
    for k in K_VALUES:
        if k >= n:
            continue
        l1 = fcluster(Z1, k, criterion="maxclust")[:n]
        l2 = fcluster(Z2, k, criterion="maxclust")[:n]
        ari_vals.append(adjusted_rand_score(l1, l2))
    return np.mean(ari_vals) if ari_vals else None


def mean_vi(d1, d2):
    """Unweighted mean VI across k values. (distance)."""
    Z1, Z2 = d1["Z"], d2["Z"]
    n = min(d1["n_nodes"], d2["n_nodes"])
    vi_vals = []
    for k in K_VALUES:
        if k >= n:
            continue
        l1 = fcluster(Z1, k, criterion="maxclust")[:n]
        l2 = fcluster(Z2, k, criterion="maxclust")[:n]
        vi_vals.append(_variation_of_information(l1, l2))
    return np.mean(vi_vals) if vi_vals else None


# ---------------------------------------------------------------------------
# Metric registry: name -> (function, is_similarity)
#   is_similarity=True:  higher = more similar
#   is_similarity=False: lower  = more similar (distance)
# ---------------------------------------------------------------------------
METRICS = {
    "CophPearson":       (cophenetic_pearson,   True),
    "CophSpearman":      (cophenetic_spearman,  True),
    "BakersGamma":       (bakers_gamma,         True),
    "NormL1_Ult":        (norm_l1_ultrametric,  False),
    "NormL2_Ult":        (norm_l2_ultrametric,  False),
    "TopK5_Merge":       (top_k5_merge,         True),
    "TopK10_Merge":      (top_k10_merge,        True),
    "TopK20_Merge":      (top_k20_merge,        True),
    "WeightARI_Coarse":  (weighted_ari_coarse,  True),
    "WeightARI_Fine":    (weighted_ari_fine,     True),
    "WeightVI_Coarse":   (weighted_vi_coarse,   False),
    "WeightVI_Fine":     (weighted_vi_fine,      False),
    "MeanARI":           (mean_ari,             True),
    "MeanVI":            (mean_vi,              False),
}


# ---------------------------------------------------------------------------
# Compute pairwise metric values
# ---------------------------------------------------------------------------
def compute_pairwise(data: dict) -> pd.DataFrame:
    """Compute all metrics for all (patient, band, phase-pair) combos."""
    rows = []
    for pat in PATIENTS:
        for band in BANDS:
            for p1, p2 in PAIR_LABELS:
                key1 = (pat, p1, band)
                key2 = (pat, p2, band)
                if key1 not in data or key2 not in data:
                    continue
                d1 = data[key1]
                d2 = data[key2]
                row = {
                    "patient": pat,
                    "band": band,
                    "phase1": p1,
                    "phase2": p2,
                    "pair": f"{p1}-{p2}",
                }
                for name, (fn, _) in METRICS.items():
                    try:
                        val = fn(d1, d2)
                        row[name] = val if val is not None else np.nan
                    except Exception as e:
                        print(f"    WARN: {name} failed for {pat} {band} {p1}-{p2}: {e}")
                        row[name] = np.nan
                rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# H2 contrast computation
# ---------------------------------------------------------------------------
def compute_h2_contrasts(df: pd.DataFrame) -> pd.DataFrame:
    """For each metric, patient, band compute H2a and H2b contrasts.

    H2a (task persists into rest_post):
      similarity: sim(TT,Post) - sim(Pre,Post)   positive = task persists
      distance:   dist(Pre,Post) - dist(TT,Post)  positive = task persists

    H2b (rest_post closer to task_test than rest_pre is):
      similarity: sim(TT,Post) - sim(Pre,TT)     positive = task persists
      distance:   dist(Pre,TT) - dist(TT,Post)    positive = task persists
    """
    rows = []
    for pat in PATIENTS:
        for band in BANDS:
            pat_band = df[(df["patient"] == pat) & (df["band"] == band)]
            # Get required pairs
            pre_post = pat_band[(pat_band["phase1"] == "rest_pre") &
                                (pat_band["phase2"] == "rest_post")]
            tt_post = pat_band[(pat_band["phase1"] == "task_test") &
                               (pat_band["phase2"] == "rest_post")]
            pre_tt = pat_band[(pat_band["phase1"] == "rest_pre") &
                              (pat_band["phase2"] == "task_test")]

            if pre_post.empty or tt_post.empty or pre_tt.empty:
                continue

            for name, (_, is_sim) in METRICS.items():
                v_pre_post = pre_post[name].values[0]
                v_tt_post = tt_post[name].values[0]
                v_pre_tt = pre_tt[name].values[0]

                if np.isnan(v_pre_post) or np.isnan(v_tt_post) or np.isnan(v_pre_tt):
                    continue

                if is_sim:
                    h2a = v_tt_post - v_pre_post   # positive = task persists
                    h2b = v_tt_post - v_pre_tt      # positive = task persists
                else:
                    h2a = v_pre_post - v_tt_post   # positive = task persists
                    h2b = v_pre_tt - v_tt_post      # positive = task persists

                rows.append({
                    "metric": name,
                    "patient": pat,
                    "band": band,
                    "H2a": h2a,
                    "H2b": h2b,
                    "is_similarity": is_sim,
                    # Raw values for diagnostics
                    "v_pre_post": v_pre_post,
                    "v_tt_post": v_tt_post,
                    "v_pre_tt": v_pre_tt,
                })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Unanimity scoring
# ---------------------------------------------------------------------------
def score_unanimity(h2_df: pd.DataFrame) -> pd.DataFrame:
    """For each (metric, band), count how many patients have positive H2a/H2b.

    Unanimity = all 4 patients positive.
    """
    rows = []
    for name in METRICS:
        for band in BANDS:
            mb = h2_df[(h2_df["metric"] == name) & (h2_df["band"] == band)]
            if mb.empty:
                continue

            n_pats = len(mb)
            h2a_pos = (mb["H2a"] > 0).sum()
            h2b_pos = (mb["H2b"] > 0).sum()
            h2a_unan = h2a_pos == n_pats
            h2b_unan = h2b_pos == n_pats
            h2a_mean = mb["H2a"].mean()
            h2b_mean = mb["H2b"].mean()

            rows.append({
                "metric": name,
                "band": band,
                "n_patients": n_pats,
                "H2a_positive": h2a_pos,
                "H2b_positive": h2b_pos,
                "H2a_unanimous": h2a_unan,
                "H2b_unanimous": h2b_unan,
                "H2a_mean": h2a_mean,
                "H2b_mean": h2b_mean,
            })
    return pd.DataFrame(rows)


def summarize_metrics(unanimity_df: pd.DataFrame) -> pd.DataFrame:
    """Summarize: for each metric, how many bands have H2a/H2b unanimity."""
    rows = []
    for name in METRICS:
        mdf = unanimity_df[unanimity_df["metric"] == name]
        h2a_unan_bands = mdf["H2a_unanimous"].sum()
        h2b_unan_bands = mdf["H2b_unanimous"].sum()
        total_bands = len(mdf)
        # Which bands are unanimous?
        h2a_bands = list(mdf[mdf["H2a_unanimous"]]["band"])
        h2b_bands = list(mdf[mdf["H2b_unanimous"]]["band"])

        rows.append({
            "metric": name,
            "H2a_unan_bands": int(h2a_unan_bands),
            "H2b_unan_bands": int(h2b_unan_bands),
            "H2_combined": int(h2a_unan_bands + h2b_unan_bands),
            "total_bands": total_bands,
            "H2a_which": ", ".join(h2a_bands) if h2a_bands else "-",
            "H2b_which": ", ".join(h2b_bands) if h2b_bands else "-",
        })
    summary = pd.DataFrame(rows).sort_values(
        ["H2_combined", "H2a_unan_bands", "H2b_unan_bands"],
        ascending=[False, False, False]
    ).reset_index(drop=True)
    return summary


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_summary_heatmap(unanimity_df: pd.DataFrame, summary_df: pd.DataFrame,
                         out_dir: Path):
    """Heatmap: rows = metrics (sorted by H2 score), columns = bands.
    Color = H2 unanimity (fraction of patients positive).
    Two panels: H2a and H2b.
    """
    # Sort metrics by combined H2 score
    metric_order = summary_df["metric"].tolist()
    band_order = BANDS

    fig, axes = plt.subplots(1, 2, figsize=(14, max(5, len(metric_order) * 0.4)))

    for ax_idx, (hyp, label) in enumerate([
        ("H2a", "H2a: sim(TT,Post) vs sim(Pre,Post)"),
        ("H2b", "H2b: sim(TT,Post) vs sim(Pre,TT)"),
    ]):
        ax = axes[ax_idx]
        pos_col = f"{hyp}_positive"
        n_col = "n_patients"

        matrix = np.full((len(metric_order), len(band_order)), np.nan)
        for i, metric in enumerate(metric_order):
            for j, band in enumerate(band_order):
                row = unanimity_df[
                    (unanimity_df["metric"] == metric) &
                    (unanimity_df["band"] == band)
                ]
                if not row.empty:
                    frac = row[pos_col].values[0] / row[n_col].values[0]
                    matrix[i, j] = frac

        im = ax.imshow(matrix, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)

        ax.set_xticks(range(len(band_order)))
        ax.set_xticklabels(
            [BRAIN_BAND_TEX_DICT.get(b, b) for b in band_order],
            fontsize=10,
        )
        ax.set_yticks(range(len(metric_order)))
        ax.set_yticklabels(metric_order, fontsize=9)
        ax.set_title(label, fontsize=11, fontweight="bold")

        # Annotate
        for i in range(len(metric_order)):
            for j in range(len(band_order)):
                val = matrix[i, j]
                if np.isfinite(val):
                    n_pos = int(round(val * 4))
                    text = f"{n_pos}/4"
                    color = "white" if val < 0.3 or val > 0.8 else "black"
                    fontweight = "bold" if val == 1.0 else "normal"
                    ax.text(j, i, text, ha="center", va="center",
                            fontsize=8, color=color, fontweight=fontweight)

    fig.colorbar(im, ax=axes, shrink=0.6, label="Fraction patients with positive H2")
    fig.suptitle(
        "H2 Task Trace Unanimity: Alternative Dendrogram Metrics\n"
        "(positive = rest_post resembles task_test more than baseline)",
        fontsize=13, fontweight="bold",
    )
    plt.tight_layout()
    path = out_dir / "h2_unanimity_heatmap.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


def plot_h2_contrast_profiles(h2_df: pd.DataFrame, out_dir: Path):
    """Per-metric, per-patient H2a contrast across bands (line profiles)."""
    metric_names = list(METRICS.keys())
    n_metrics = len(metric_names)
    n_cols = 4
    n_rows = (n_metrics + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.2 * n_rows),
                             squeeze=False)

    for idx, name in enumerate(metric_names):
        ax = axes[idx // n_cols, idx % n_cols]
        mdf = h2_df[h2_df["metric"] == name]

        for pat in PATIENTS:
            pat_data = mdf[mdf["patient"] == pat]
            if pat_data.empty:
                continue
            vals = []
            for band in BANDS:
                row = pat_data[pat_data["band"] == band]
                vals.append(row["H2a"].values[0] if not row.empty else np.nan)
            ax.plot(range(len(BANDS)), vals, marker="o", markersize=4,
                    label=pat, linewidth=1.5)

        ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels(
            [BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS],
            fontsize=8,
        )
        ax.set_title(name, fontsize=10, fontweight="bold")
        ax.set_ylabel("H2a contrast", fontsize=8)
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(fontsize=7, ncol=2)

    # Remove empty axes
    for idx in range(n_metrics, n_rows * n_cols):
        axes[idx // n_cols, idx % n_cols].set_visible(False)

    fig.suptitle(
        "H2a Contrast per Patient and Band\n"
        "(positive = rest_post closer to task_test than to rest_pre)",
        fontsize=13, fontweight="bold",
    )
    plt.tight_layout()
    path = out_dir / "h2a_contrast_profiles.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


def plot_h2b_contrast_profiles(h2_df: pd.DataFrame, out_dir: Path):
    """Per-metric, per-patient H2b contrast across bands (line profiles)."""
    metric_names = list(METRICS.keys())
    n_metrics = len(metric_names)
    n_cols = 4
    n_rows = (n_metrics + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.2 * n_rows),
                             squeeze=False)

    for idx, name in enumerate(metric_names):
        ax = axes[idx // n_cols, idx % n_cols]
        mdf = h2_df[h2_df["metric"] == name]

        for pat in PATIENTS:
            pat_data = mdf[mdf["patient"] == pat]
            if pat_data.empty:
                continue
            vals = []
            for band in BANDS:
                row = pat_data[pat_data["band"] == band]
                vals.append(row["H2b"].values[0] if not row.empty else np.nan)
            ax.plot(range(len(BANDS)), vals, marker="o", markersize=4,
                    label=pat, linewidth=1.5)

        ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels(
            [BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS],
            fontsize=8,
        )
        ax.set_title(name, fontsize=10, fontweight="bold")
        ax.set_ylabel("H2b contrast", fontsize=8)
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(fontsize=7, ncol=2)

    for idx in range(n_metrics, n_rows * n_cols):
        axes[idx // n_cols, idx % n_cols].set_visible(False)

    fig.suptitle(
        "H2b Contrast per Patient and Band\n"
        "(positive = rest_post closer to task_test than rest_pre is)",
        fontsize=13, fontweight="bold",
    )
    plt.tight_layout()
    path = out_dir / "h2b_contrast_profiles.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


def plot_per_patient_detail(h2_df: pd.DataFrame, unanimity_df: pd.DataFrame,
                            out_dir: Path):
    """Detailed heatmap: for each metric, show per-patient H2a sign for every band."""
    metric_order = list(METRICS.keys())

    fig, ax = plt.subplots(figsize=(12, max(5, len(metric_order) * 0.45)))

    # Build matrix: rows = metrics, columns = (band, patient) pairs
    col_labels = []
    for band in BANDS:
        for pat in PATIENTS:
            col_labels.append(f"{BRAIN_BAND_TEX_DICT.get(band, band)}\n{pat[-2:]}")

    matrix = np.full((len(metric_order), len(col_labels)), np.nan)
    for i, metric in enumerate(metric_order):
        for j_band, band in enumerate(BANDS):
            for j_pat, pat in enumerate(PATIENTS):
                col_idx = j_band * len(PATIENTS) + j_pat
                row = h2_df[
                    (h2_df["metric"] == metric) &
                    (h2_df["band"] == band) &
                    (h2_df["patient"] == pat)
                ]
                if not row.empty:
                    h2a = row["H2a"].values[0]
                    matrix[i, col_idx] = 1.0 if h2a > 0 else 0.0

    im = ax.imshow(matrix, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=7, rotation=0)
    ax.set_yticks(range(len(metric_order)))
    ax.set_yticklabels(metric_order, fontsize=9)

    # Annotate
    for i in range(len(metric_order)):
        for j in range(len(col_labels)):
            val = matrix[i, j]
            if np.isfinite(val):
                text = "+" if val == 1.0 else "-"
                color = "white" if val < 0.5 else "black"
                ax.text(j, i, text, ha="center", va="center",
                        fontsize=7, color=color, fontweight="bold")

    # Draw vertical separators between bands
    for k in range(1, len(BANDS)):
        ax.axvline(k * len(PATIENTS) - 0.5, color="black", lw=1.0)

    ax.set_title(
        "H2a Sign per Patient: + = task trace, - = no trace\n"
        "(green = positive, red = negative)",
        fontsize=11, fontweight="bold",
    )
    plt.tight_layout()
    path = out_dir / "h2a_per_patient_detail.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


def plot_raw_values_diagnostic(df: pd.DataFrame, out_dir: Path):
    """Show raw metric values for the 3 key pairs (Pre-Post, TT-Post, Pre-TT),
    averaged across patients, to diagnose sensitivity."""
    key_pairs = [
        ("rest_pre", "rest_post", "Pre-Post"),
        ("task_test", "rest_post", "TT-Post"),
        ("rest_pre", "task_test", "Pre-TT"),
    ]
    metric_names = list(METRICS.keys())
    n_metrics = len(metric_names)
    n_cols = 4
    n_rows = (n_metrics + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5 * n_cols, 3.2 * n_rows),
                             squeeze=False)

    for idx, name in enumerate(metric_names):
        ax = axes[idx // n_cols, idx % n_cols]

        for p1, p2, label in key_pairs:
            pair_data = df[(df["phase1"] == p1) & (df["phase2"] == p2)]
            means = []
            sems = []
            for band in BANDS:
                vals = pair_data[pair_data["band"] == band][name].dropna().values
                means.append(np.mean(vals) if len(vals) > 0 else np.nan)
                sems.append(np.std(vals) / np.sqrt(len(vals)) if len(vals) > 1 else 0)

            ax.errorbar(range(len(BANDS)), means, yerr=sems,
                        marker="o", markersize=4, label=label,
                        capsize=3, linewidth=1.5)

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels(
            [BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS], fontsize=8
        )
        ax.set_title(name, fontsize=10, fontweight="bold")
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(fontsize=7)

    for idx in range(n_metrics, n_rows * n_cols):
        axes[idx // n_cols, idx % n_cols].set_visible(False)

    fig.suptitle(
        "Raw Metric Values: Key Phase Pairs (mean +/- SEM across patients)",
        fontsize=13, fontweight="bold",
    )
    plt.tight_layout()
    path = out_dir / "raw_values_diagnostic.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("ALTERNATIVE DENDROGRAM METRICS - H2 (Task Trace) SWEEP")
    print("=" * 70)

    # Load data
    print("\n[1/5] Loading LRG data...")
    data = load_all()
    print(f"  Loaded {len(data)} entries ({len(PATIENTS)} patients x "
          f"{len(PHASES)} phases x {len(BANDS)} bands)")

    # Compute pairwise metrics
    print("\n[2/5] Computing all metrics for all phase pairs...")
    df = compute_pairwise(data)
    print(f"  Computed {len(df)} pair comparisons across {len(METRICS)} metrics")

    # Save raw results
    csv_path = OUT_DIR / "pairwise_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"  Saved to {csv_path}")

    # Compute H2 contrasts
    print("\n[3/5] Computing H2 contrasts...")
    h2_df = compute_h2_contrasts(df)
    h2_csv_path = OUT_DIR / "h2_contrasts.csv"
    h2_df.to_csv(h2_csv_path, index=False)
    print(f"  Saved to {h2_csv_path}")

    # Score unanimity
    print("\n[4/5] Scoring unanimity...")
    unanimity_df = score_unanimity(h2_df)
    unanimity_csv = OUT_DIR / "h2_unanimity.csv"
    unanimity_df.to_csv(unanimity_csv, index=False)

    summary_df = summarize_metrics(unanimity_df)
    summary_csv = OUT_DIR / "h2_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"  Saved unanimity to {unanimity_csv}")
    print(f"  Saved summary to {summary_csv}")

    # Print summary table
    print("\n" + "=" * 90)
    print("H2 UNANIMITY SUMMARY (sorted by combined H2a+H2b unanimous bands)")
    print("=" * 90)
    print(f"{'Rank':<5} {'Metric':<20} {'H2a':>4} {'H2b':>4} {'Comb':>5}  "
          f"{'H2a bands':<30} {'H2b bands'}")
    print("-" * 95)
    for i, row in summary_df.iterrows():
        print(f"{i+1:<5} {row['metric']:<20} {row['H2a_unan_bands']:>4}/6 "
              f"{row['H2b_unan_bands']:>4}/6 {row['H2_combined']:>5}/12  "
              f"{row['H2a_which']:<30} {row['H2b_which']}")

    # Detailed per-patient, per-band results for top 3
    print("\n" + "=" * 90)
    print("TOP METRICS - Per-patient H2a detail")
    print("=" * 90)
    for _, srow in summary_df.head(5).iterrows():
        name = srow["metric"]
        print(f"\n--- {name} (H2a: {srow['H2a_unan_bands']}/6, "
              f"H2b: {srow['H2b_unan_bands']}/6) ---")
        print(f"  {'Band':<12}", end="")
        for pat in PATIENTS:
            print(f"  {pat:>10}", end="")
        print(f"  {'Unan':>6}")

        for band in BANDS:
            print(f"  {band:<12}", end="")
            vals = []
            for pat in PATIENTS:
                row = h2_df[
                    (h2_df["metric"] == name) &
                    (h2_df["band"] == band) &
                    (h2_df["patient"] == pat)
                ]
                if not row.empty:
                    v = row["H2a"].values[0]
                    vals.append(v)
                    print(f"  {v:>+10.5f}", end="")
                else:
                    print(f"  {'N/A':>10}", end="")
            unan = all(v > 0 for v in vals) if vals else False
            sign = "YES" if unan else "no"
            print(f"  {sign:>6}")

    # Plots
    print("\n[5/5] Generating figures...")
    plot_summary_heatmap(unanimity_df, summary_df, OUT_DIR)
    plot_h2_contrast_profiles(h2_df, OUT_DIR)
    plot_h2b_contrast_profiles(h2_df, OUT_DIR)
    plot_per_patient_detail(h2_df, unanimity_df, OUT_DIR)
    plot_raw_values_diagnostic(df, OUT_DIR)

    # Final summary
    print("\n" + "=" * 90)
    print("FINDINGS")
    print("=" * 90)

    # Best metric for H2a
    best_h2a = summary_df.sort_values("H2a_unan_bands", ascending=False).iloc[0]
    print(f"\n  Best metric for H2a: {best_h2a['metric']} "
          f"({best_h2a['H2a_unan_bands']}/6 bands unanimous)")
    print(f"    Unanimous bands: {best_h2a['H2a_which']}")

    # Best metric for H2b
    best_h2b = summary_df.sort_values("H2b_unan_bands", ascending=False).iloc[0]
    print(f"\n  Best metric for H2b: {best_h2b['metric']} "
          f"({best_h2b['H2b_unan_bands']}/6 bands unanimous)")
    print(f"    Unanimous bands: {best_h2b['H2b_which']}")

    # Best combined
    best_comb = summary_df.iloc[0]
    print(f"\n  Best combined: {best_comb['metric']} "
          f"(H2a: {best_comb['H2a_unan_bands']}/6, H2b: {best_comb['H2b_unan_bands']}/6)")

    # Check which bands are consistently found across metrics
    print("\n  Band consistency across all metrics:")
    for band in BANDS:
        band_unan = unanimity_df[unanimity_df["band"] == band]
        h2a_count = band_unan["H2a_unanimous"].sum()
        h2b_count = band_unan["H2b_unanimous"].sum()
        print(f"    {band:<12}: H2a unanimous in {h2a_count:>2}/{len(METRICS)} metrics, "
              f"H2b in {h2b_count:>2}/{len(METRICS)} metrics")

    print(f"\n  All outputs in: {OUT_DIR}")
    print("=" * 90)


if __name__ == "__main__":
    main()
