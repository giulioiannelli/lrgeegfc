#!/usr/bin/env python3
"""Massive sweep of TREE-STRUCTURAL similarity metrics.

Goal: find a metric where the within-vs-cross gap has CONSISTENT SIGN
across ALL patients for each band.

All metrics compare dendrogram structure (cophenetic distances, merge order,
tree topology) — NOT flat partitions like ARI.

For each metric, computes:
  similarity(tree_a, tree_b) → float in some range

Then for each (patient, band):
  gap = mean_similarity(within_pairs) - mean_similarity(cross_pairs)
  gap > 0 → within-condition trees more similar → REORGANIZATION

The winning metric is the one where gap sign is unanimous across patients
for the most bands.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, spearmanr, kendalltau, wasserstein_distance

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import LRG_CACHE
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrgsglib.utils.basic.linalg import (
    tree_baker_gamma,
    tree_cophenetic_correlation,
    tree_robinson_foulds_distance,
    ultrametric_rank_correlation,
    ultrametric_scaled_distance,
)

# ---------------------------------------------------------------------------
PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
PATIENTS = list(PATIENT_PHASES.keys())
BANDS = BRAIN_BANDS_NAMES

REST = {"rest_pre", "rest_post"}
TASK = {"task_learn", "task_test"}


def classify_pair(p1, p2):
    s = {p1, p2}
    if s <= REST or s <= TASK:
        return "within"
    return "cross"


# ===================================================================
# Pre-compute all needed data
# ===================================================================
def load_all():
    """Load linkage + cophenetic for all (patient, phase, band)."""
    data = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for ph in phases:
                lrg = load_lrg_result(pat, ph, band, "msc", cache_root=LRG_CACHE)
                if lrg is None:
                    continue
                Z = lrg.linkage_matrix
                c_condensed = cophenet(Z)
                c_sq = squareform(c_condensed)
                U = lrg.ultrametric_matrix
                U_sq = squareform(U) if U.ndim == 1 else U
                data[(pat, ph, band)] = {
                    "Z": Z,
                    "c_cond": c_condensed,
                    "c_sq": c_sq,
                    "U_sq": U_sq,
                    "n": lrg.n_nodes,
                }
    return data


# ===================================================================
# Metric definitions — all return SIMILARITY (higher = more similar)
# ===================================================================

def m_cophenetic_pearson(d1, d2):
    """Pearson correlation of cophenetic distances."""
    return pearsonr(d1["c_cond"], d2["c_cond"])[0]


def m_cophenetic_spearman(d1, d2):
    """Spearman correlation of cophenetic distances."""
    return spearmanr(d1["c_cond"], d2["c_cond"])[0]


def m_cophenetic_kendall(d1, d2):
    """Kendall tau of cophenetic distances."""
    # Subsample for speed (Kendall is O(n^2))
    n = len(d1["c_cond"])
    if n > 5000:
        idx = np.random.RandomState(42).choice(n, 5000, replace=False)
        return kendalltau(d1["c_cond"][idx], d2["c_cond"][idx])[0]
    return kendalltau(d1["c_cond"], d2["c_cond"])[0]


def m_log_cophenetic_pearson(d1, d2):
    """Pearson correlation of LOG cophenetic distances."""
    lc1 = np.log(np.clip(d1["c_cond"], 1e-12, None))
    lc2 = np.log(np.clip(d2["c_cond"], 1e-12, None))
    return pearsonr(lc1, lc2)[0]


def m_log_cophenetic_spearman(d1, d2):
    """Spearman correlation of LOG cophenetic distances."""
    lc1 = np.log(np.clip(d1["c_cond"], 1e-12, None))
    lc2 = np.log(np.clip(d2["c_cond"], 1e-12, None))
    return spearmanr(lc1, lc2)[0]


def m_baker_gamma(d1, d2):
    """Baker's Gamma coefficient."""
    return tree_baker_gamma(d1["Z"], d2["Z"])


def m_robinson_foulds_sim(d1, d2):
    """1 - Robinson-Foulds distance (similarity)."""
    return 1.0 - tree_robinson_foulds_distance(d1["Z"], d2["Z"], normalized=True)


def m_ultrametric_rank_spearman(d1, d2):
    """Spearman rank correlation of ultrametric matrices."""
    return ultrametric_rank_correlation(d1["U_sq"], d2["U_sq"], method="spearman")


def m_ultrametric_scaled_sim(d1, d2):
    """1 - ultrametric scaled distance (log, normalized)."""
    return 1.0 - ultrametric_scaled_distance(d1["U_sq"], d2["U_sq"],
                                              scale="log", normalize=True)


def m_cophenetic_profile_corr(d1, d2):
    """Correlation of node cophenetic profiles.

    For each node i, compute mean cophenetic distance to all others.
    Then correlate these N-dimensional vectors.
    """
    p1 = d1["c_sq"].mean(axis=1)
    p2 = d2["c_sq"].mean(axis=1)
    return pearsonr(p1, p2)[0]


def m_log_cophenetic_profile_corr(d1, d2):
    """Correlation of node LOG cophenetic profiles."""
    lc1 = np.log(np.clip(d1["c_sq"], 1e-12, None))
    lc2 = np.log(np.clip(d2["c_sq"], 1e-12, None))
    p1 = lc1.mean(axis=1)
    p2 = lc2.mean(axis=1)
    return pearsonr(p1, p2)[0]


def m_merge_order_spearman(d1, d2):
    """Spearman correlation of merge orders.

    Rank pairs by cophenetic distance (merge order). Compare rankings.
    Same as cophenetic Spearman but conceptually clearer.
    """
    return spearmanr(d1["c_cond"], d2["c_cond"])[0]


def m_log_ratio_consistency(d1, d2):
    """Consistency of log-ratio of cophenetic distances.

    If trees are similar, log(c1/c2) should have low variance.
    Return: 1 / (1 + std(log(c1/c2))).
    """
    c1 = np.clip(d1["c_cond"], 1e-12, None)
    c2 = np.clip(d2["c_cond"], 1e-12, None)
    log_ratio = np.log(c1 / c2)
    return 1.0 / (1.0 + np.std(log_ratio))


def m_wasserstein_sim(d1, d2):
    """1 / (1 + Wasserstein distance of log cophenetic distributions)."""
    lc1 = np.log(np.clip(d1["c_cond"], 1e-12, None))
    lc2 = np.log(np.clip(d2["c_cond"], 1e-12, None))
    w = wasserstein_distance(lc1, lc2)
    return 1.0 / (1.0 + w)


def m_quantile_sim(d1, d2):
    """1 - RMSE of quantile profiles (in log space).

    Compare distributions at quantiles [0.05, 0.1, ..., 0.95].
    """
    qs = np.linspace(0.05, 0.95, 19)
    lc1 = np.log(np.clip(d1["c_cond"], 1e-12, None))
    lc2 = np.log(np.clip(d2["c_cond"], 1e-12, None))
    q1 = np.quantile(lc1, qs)
    q2 = np.quantile(lc2, qs)
    # Normalize by range
    r = max(np.ptp(q1), np.ptp(q2), 1e-12)
    rmse = np.sqrt(np.mean((q1 - q2) ** 2)) / r
    return 1.0 - min(rmse, 1.0)


def m_top_k_overlap(d1, d2, k_frac=0.1):
    """Jaccard overlap of top-k strongest connections.

    Rank pairs by cophenetic distance (lowest = closest = strongest).
    Check overlap of bottom k pairs.
    """
    n = len(d1["c_cond"])
    k = max(int(n * k_frac), 10)
    idx1 = np.argsort(d1["c_cond"])[:k]
    idx2 = np.argsort(d2["c_cond"])[:k]
    overlap = len(set(idx1) & set(idx2))
    return overlap / k


def m_hierarchical_embedding_sim(d1, d2):
    """Cosine similarity of hierarchical embeddings.

    Use sorted cophenetic distances as embedding vectors.
    """
    s1 = np.sort(d1["c_cond"])
    s2 = np.sort(d2["c_cond"])
    dot = np.dot(s1, s2)
    norm = np.linalg.norm(s1) * np.linalg.norm(s2)
    if norm < 1e-12:
        return 0.0
    return dot / norm


def m_log_hierarchical_embedding_sim(d1, d2):
    """Cosine similarity of sorted LOG cophenetic distances."""
    s1 = np.sort(np.log(np.clip(d1["c_cond"], 1e-12, None)))
    s2 = np.sort(np.log(np.clip(d2["c_cond"], 1e-12, None)))
    dot = np.dot(s1, s2)
    norm = np.linalg.norm(s1) * np.linalg.norm(s2)
    if norm < 1e-12:
        return 0.0
    return dot / norm


# ===================================================================
# Registry
# ===================================================================
METRICS = {
    "CophPearson": m_cophenetic_pearson,
    "CophSpearman": m_cophenetic_spearman,
    "CophKendall": m_cophenetic_kendall,
    "LogCophPearson": m_log_cophenetic_pearson,
    "LogCophSpearman": m_log_cophenetic_spearman,
    "BakerGamma": m_baker_gamma,
    "1-RF": m_robinson_foulds_sim,
    "UltRankSpearman": m_ultrametric_rank_spearman,
    "1-UltScaled": m_ultrametric_scaled_sim,
    "CophProfile": m_cophenetic_profile_corr,
    "LogCophProfile": m_log_cophenetic_profile_corr,
    "LogRatioConsist": m_log_ratio_consistency,
    "WassersteinSim": m_wasserstein_sim,
    "QuantileSim": m_quantile_sim,
    "TopKOverlap": m_top_k_overlap,
    "HierEmbedCos": m_hierarchical_embedding_sim,
    "LogHierEmbedCos": m_log_hierarchical_embedding_sim,
}


# ===================================================================
# Main computation
# ===================================================================
def compute_gaps(data):
    """For each (metric, patient, band), compute within-cross gap.

    Returns dict: (metric_name, patient, band) → gap (float)
    Also returns raw similarities.
    """
    gaps = {}
    sims = {}

    for metric_name, metric_fn in METRICS.items():
        print(f"\n  Computing {metric_name}...")
        for pat in PATIENTS:
            phases = PATIENT_PHASES[pat]
            for band in BANDS:
                available = [ph for ph in phases if (pat, ph, band) in data]
                pairs = list(combinations(available, 2))

                within_sims = []
                cross_sims = []
                for p1, p2 in pairs:
                    d1 = data[(pat, p1, band)]
                    d2 = data[(pat, p2, band)]
                    try:
                        s = metric_fn(d1, d2)
                        if s is None or np.isnan(s):
                            continue
                    except Exception:
                        continue

                    cat = classify_pair(p1, p2)
                    sims[(metric_name, pat, band, (p1, p2))] = s
                    if cat == "within":
                        within_sims.append(s)
                    else:
                        cross_sims.append(s)

                if within_sims and cross_sims:
                    gap = np.mean(within_sims) - np.mean(cross_sims)
                    gaps[(metric_name, pat, band)] = gap

    return gaps, sims


def score_metrics(gaps):
    """Score each metric by cross-patient consistency.

    For each (metric, band), check if ALL patients agree on gap sign.
    """
    results = []

    for metric_name in METRICS:
        band_scores = []
        for band in BANDS:
            pat_gaps = []
            for pat in PATIENTS:
                key = (metric_name, pat, band)
                if key in gaps:
                    pat_gaps.append(gaps[key])

            if len(pat_gaps) < 3:
                continue

            all_pos = all(g > 0 for g in pat_gaps)
            all_neg = all(g < 0 for g in pat_gaps)
            unanimous = all_pos or all_neg
            mean_gap = np.mean(pat_gaps)
            std_gap = np.std(pat_gaps)
            min_gap = min(pat_gaps)
            max_gap = max(pat_gaps)
            cv = std_gap / abs(mean_gap) if abs(mean_gap) > 1e-6 else 999

            band_scores.append({
                "band": band,
                "unanimous": unanimous,
                "direction": "reorg" if all_pos else ("persist" if all_neg else "mixed"),
                "mean_gap": mean_gap,
                "std_gap": std_gap,
                "cv": cv,
                "min_gap": min_gap,
                "max_gap": max_gap,
                "n_pos": sum(g > 0 for g in pat_gaps),
                "n_neg": sum(g < 0 for g in pat_gaps),
                "n_patients": len(pat_gaps),
            })

        n_unanimous = sum(b["unanimous"] for b in band_scores)
        mean_cv = np.mean([b["cv"] for b in band_scores]) if band_scores else 999
        # Bonus: count bands with strong unanimous gap
        n_strong = sum(b["unanimous"] and abs(b["mean_gap"]) > 0.02 for b in band_scores)

        results.append({
            "metric": metric_name,
            "n_unanimous": n_unanimous,
            "n_strong_unanimous": n_strong,
            "mean_cv": mean_cv,
            "band_scores": band_scores,
        })

    # Sort: most unanimous bands, then strongest, then lowest CV
    results.sort(key=lambda r: (-r["n_unanimous"], -r["n_strong_unanimous"], r["mean_cv"]))
    return results


def main():
    print("Loading data...")
    data = load_all()
    print(f"  Loaded {len(data)} entries")

    print("\nComputing all metrics...")
    gaps, sims = compute_gaps(data)

    print("\n\nScoring metrics by cross-patient consistency...")
    results = score_metrics(gaps)

    # Print ranking
    print("\n" + "=" * 90)
    print("METRIC RANKING — by cross-patient consistency")
    print("(unanimous = ALL patients agree on gap sign for that band)")
    print("=" * 90)
    print(f"{'Rank':<5} {'Metric':<20} {'Unanimous':>10} {'Strong':>7} {'MeanCV':>8}")
    print("-" * 55)
    for i, r in enumerate(results):
        print(f"{i+1:<5} {r['metric']:<20} {r['n_unanimous']:>10}/6 "
              f"{r['n_strong_unanimous']:>7}/6 {r['mean_cv']:>8.2f}")

    # Detail for top 5
    print("\n\n" + "=" * 90)
    print("TOP METRICS — detailed band breakdown")
    print("=" * 90)
    for r in results[:5]:
        print(f"\n--- {r['metric']} ({r['n_unanimous']}/6 unanimous) ---")
        print(f"  {'Band':<12} {'Dir':>7} {'Mean':>8} {'Std':>8} {'Min':>8} "
              f"{'Max':>8} {'Pos':>4}/{' Neg':>4} {'Unani':>6}")
        for bs in r["band_scores"]:
            print(f"  {bs['band']:<12} {bs['direction']:>7} {bs['mean_gap']:>+8.4f} "
                  f"{bs['std_gap']:>8.4f} {bs['min_gap']:>+8.4f} {bs['max_gap']:>+8.4f} "
                  f"{bs['n_pos']:>4}/{bs['n_neg']:>4} {'YES' if bs['unanimous'] else 'no':>6}")

    # Per-patient detail for the winner
    if results:
        winner = results[0]["metric"]
        print(f"\n\n{'='*90}")
        print(f"WINNER: {winner} — per-patient gaps")
        print(f"{'='*90}")
        print(f"  {'Band':<12}", end="")
        for pat in PATIENTS:
            print(f"  {pat:>8}", end="")
        print()
        for band in BANDS:
            print(f"  {band:<12}", end="")
            for pat in PATIENTS:
                key = (winner, pat, band)
                if key in gaps:
                    g = gaps[key]
                    print(f"  {g:>+8.4f}", end="")
                else:
                    print(f"  {'N/A':>8}", end="")
            print()


if __name__ == "__main__":
    main()
