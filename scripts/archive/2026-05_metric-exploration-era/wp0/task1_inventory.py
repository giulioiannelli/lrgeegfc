#!/usr/bin/env python3
"""WP0 Task 1 — Generate metric inventory document.

Produces:
  data/wp0_metric_exploration/task1_metric_inventory.md

Run: python scripts/wp0/task1_inventory.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from _common import OUT_ROOT, _METRIC_META  # noqa: E402

OUT = OUT_ROOT / "task1_metric_inventory.md"


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# WP0 Task 1 — Metric Inventory",
        "",
        "Complete catalogue of all comparison metrics used in this work package.",
        "",
        "## Overview",
        "",
        "| # | Name | Category | Distance? | Range | Reference value |",
        "|---|------|----------|-----------|-------|-----------------|",
    ]

    categories = {
        "Ultrametric value-based": [
            "matrix_distance", "scaled_distance", "rank_distance",
            "quantile_rmse", "permutation_robust",
        ],
        "Tree-structural": [
            "tree_robinson_foulds", "tree_cophenetic_corr",
            "tree_baker_gamma", "tree_fowlkes_mallows",
        ],
        "Partition-based": ["ari", "cluster_swap"],
        "Raw FC baselines": [
            "fc_frobenius", "fc_scaled_frobenius",
            "fc_rank_distance", "fc_mean_abs_diff",
        ],
    }

    idx = 1
    for cat, names in categories.items():
        for name in names:
            label, is_dist, rng, ref = _METRIC_META[name]
            dist_str = "Yes" if is_dist else "No (similarity)"
            lines.append(f"| {idx} | `{name}` | {cat} | {dist_str} | {rng} | {ref} |")
            idx += 1

    lines += [
        "",
        f"**Total: {idx - 1} metrics** (11 LRG-based + 4 raw FC baselines)",
        "",
        "---",
        "",
        "## Detailed Descriptions",
        "",
    ]

    # Detailed per-category sections
    details = {
        "matrix_distance": {
            "what": "Element-wise Euclidean distance between upper triangles of two ultrametric (cophenetic) distance matrices.",
            "input": "Two N×N ultrametric matrices U1, U2 (from LRG dendrograms at fixed τ = τ' = 1/λ_max).",
            "formula": "d = ||vec_upper(U1) - vec_upper(U2)||_2",
            "sensitivity": "Dominated by large-scale merges (high ultrametric distances). Unbounded — depends on matrix scale.",
            "issues": "Not normalized; values not directly comparable across patients with different N or different eigenvalue spectra.",
        },
        "scaled_distance": {
            "what": "Euclidean distance after log-transforming and normalizing both ultrametric vectors to unit norm.",
            "input": "Two N×N ultrametric matrices U1, U2.",
            "formula": "v_i = log10(vec_upper(U_i)) / ||log10(vec_upper(U_i))||_2; d = ||v1 - v2||_2",
            "sensitivity": "Scale-invariant. Sensitive to relative ordering of merge distances across all levels.",
            "issues": "Well-behaved. Bounded [0, 1] when normalized. Good candidate for cross-patient comparison.",
        },
        "rank_distance": {
            "what": "One minus Spearman rank correlation between upper-triangle vectors of two ultrametric matrices.",
            "input": "Two N×N ultrametric matrices U1, U2.",
            "formula": "d = 1 - ρ_Spearman(vec_upper(U1), vec_upper(U2))",
            "sensitivity": "Purely rank-based — ignores magnitude, sensitive only to ordering changes in pairwise distances.",
            "issues": "Bounded [0, 2]. Very robust to scale differences. May miss magnitude-only changes.",
        },
        "quantile_rmse": {
            "what": "RMSE between quantile profiles (in log space) of two ultrametric vectors.",
            "input": "Two N×N ultrametric matrices U1, U2. Default quantiles: (0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0).",
            "formula": "d = sqrt(mean((q_log(U1) - q_log(U2))^2)) over 7 quantile points",
            "sensitivity": "Distributional — captures shape differences (tail behavior, central tendency).",
            "issues": "Unbounded. Compressed representation — only 7 numbers summarize the full distribution.",
        },
        "permutation_robust": {
            "what": "Ultrametric matrix distance after reordering leaves to a canonical order (controls for labeling permutations).",
            "input": "Two linkage matrices Z1, Z2.",
            "formula": "Reorder to canonical leaf order, then Euclidean distance of ultrametric vectors.",
            "sensitivity": "Controls for the arbitrary leaf ordering in dendrograms. Same as matrix_distance modulo permutation.",
            "issues": "Unbounded. Computationally heavier than matrix_distance.",
        },
        "tree_robinson_foulds": {
            "what": "Fraction of bipartitions (internal edges) present in one tree but not the other.",
            "input": "Two linkage matrices Z1, Z2.",
            "formula": "RF = |Σ(T1) △ Σ(T2)| / (2(N-3)); normalized to [0,1]",
            "sensitivity": "Purely topological — ignores merge heights, only compares which splits exist.",
            "issues": "Very coarse for small trees. Can saturate (go to 1.0) even for modestly different trees.",
        },
        "tree_cophenetic_corr": {
            "what": "Pearson correlation between cophenetic distance vectors of two trees.",
            "input": "Two linkage matrices Z1, Z2.",
            "formula": "r = Pearson(coph(Z1), coph(Z2))",
            "sensitivity": "Linear relationship between merge distances — captures both topology and relative heights.",
            "issues": "Can be high even when trees differ in absolute scale. Similarity metric (1 = identical).",
        },
        "tree_baker_gamma": {
            "what": "Rank concordance (Goodman-Kruskal gamma) of pairwise distances in two trees.",
            "input": "Two linkage matrices Z1, Z2.",
            "formula": "γ = (C - D) / (C + D) where C = concordant, D = discordant pairs",
            "sensitivity": "Rank-based version of cophenetic correlation. More robust to outlier merge heights.",
            "issues": "Similarity metric (1 = identical). Can be slow for large N (O(N^4) pairs).",
        },
        "tree_fowlkes_mallows": {
            "what": "Geometric mean of precision and recall when comparing flat clusterings from two trees.",
            "input": "Two linkage matrices Z1, Z2, cut at k clusters (default: auto).",
            "formula": "FM = sqrt(precision × recall) where precision/recall measure shared node pairs",
            "sensitivity": "Depends on cut level k. Captures cluster-level agreement, not full tree structure.",
            "issues": "Similarity metric (1 = identical). Sensitive to choice of k. [0, 1].",
        },
        "ari": {
            "what": "Adjusted Rand Index between flat clusterings at the Ψ-optimal dendrogram cut.",
            "input": "Two label vectors from fcluster(Z, optimal_threshold, criterion='distance').",
            "formula": "ARI = (RI - E[RI]) / (max(RI) - E[RI])",
            "sensitivity": "Measures partition agreement at a single (Ψ-optimal) hierarchical level.",
            "issues": "Corrected for chance (0 = random, 1 = identical). Can be negative. Depends on optimal_threshold stability.",
        },
        "cluster_swap": {
            "what": "Fraction of nodes that remain in the same community after optimal label alignment.",
            "input": "Two label vectors (same as ARI), aligned via Hungarian algorithm.",
            "formula": "swap = mean(aligned_labels_a == labels_b)",
            "sensitivity": "Intuitive: fraction of nodes that 'stay put'. Not corrected for chance.",
            "issues": "Similarity metric. [0, 1]. More interpretable than ARI but less statistically rigorous.",
        },
        "fc_frobenius": {
            "what": "Frobenius (L2) distance between upper triangles of raw FC adjacency matrices.",
            "input": "Two N×N MSC matrices W1, W2 (no LRG processing).",
            "formula": "d = ||vec_upper(W1) - vec_upper(W2)||_2",
            "sensitivity": "Dominated by large-magnitude entries. Captures overall FC magnitude changes.",
            "issues": "Unbounded — depends on N and FC magnitude. Not comparable across patients without normalization.",
        },
        "fc_scaled_frobenius": {
            "what": "Frobenius distance normalized by total magnitude of both matrices.",
            "input": "Two N×N MSC matrices W1, W2.",
            "formula": "d = ||v1 - v2||_2 / (||v1||_2 + ||v2||_2)",
            "sensitivity": "Scale-invariant version of Frobenius. Focuses on relative pattern changes.",
            "issues": "Bounded [0, 1]. Good candidate for cross-patient comparison.",
        },
        "fc_rank_distance": {
            "what": "One minus Spearman rank correlation between upper triangles of raw FC matrices.",
            "input": "Two N×N MSC matrices W1, W2.",
            "formula": "d = 1 - ρ_Spearman(vec_upper(W1), vec_upper(W2))",
            "sensitivity": "Rank-based — detects reordering of edge strengths regardless of magnitude.",
            "issues": "Bounded [0, 2]. Comparable across patients (rank-invariant to N and scale).",
        },
        "fc_mean_abs_diff": {
            "what": "Mean absolute difference between upper-triangle entries of raw FC matrices.",
            "input": "Two N×N MSC matrices W1, W2.",
            "formula": "d = mean(|W1_ij - W2_ij|) over upper triangle",
            "sensitivity": "Average per-edge change. Less dominated by outlier edges than Frobenius.",
            "issues": "Unbounded — depends on FC magnitude scale. Not directly comparable across patients.",
        },
    }

    for cat, names in categories.items():
        lines.append(f"### {cat}")
        lines.append("")
        for name in names:
            d = details[name]
            label = _METRIC_META[name][0]
            lines.append(f"#### `{name}` — {label}")
            lines.append("")
            lines.append(f"- **What it compares:** {d['what']}")
            lines.append(f"- **Input:** {d['input']}")
            lines.append(f"- **Formula:** `{d['formula']}`")
            lines.append(f"- **Sensitivity:** {d['sensitivity']}")
            lines.append(f"- **Known issues:** {d['issues']}")
            lines.append("")

    lines += [
        "---",
        "",
        "## Source Locations",
        "",
        "| Category | Source file |",
        "|----------|-----------|",
        "| Ultrametric value-based (5) | `lrgsglib/src/lrgsglib/utils/basic/linalg.py` |",
        "| Tree-structural (4) | `lrgsglib/src/lrgsglib/utils/basic/linalg.py` |",
        "| Partition-based (2) | `src/lrg_eegfc/utils/metrics/reorganization.py` |",
        "| Raw FC baselines (4) | `scripts/wp0/_common.py` (new for WP0) |",
        "| Metric registry | `src/lrg_eegfc/utils/metrics/reorganization.py:build_metric_specs()` |",
        "",
        "## Resolution Axis Mapping",
        "",
        "| Metric | τ-independent? | h-multiscale? | Notes |",
        "|--------|---------------|--------------|-------|",
        "| matrix_distance through permutation_robust | Fixed τ = τ' | No (uses full ultrametric) | Compare D(τ') matrices |",
        "| tree_* metrics | Fixed τ = τ' | No (uses full tree) | Compare dendrogram topology |",
        "| ari, cluster_swap | Fixed τ = τ' | YES — can vary cut level h | Compare partitions at h |",
        "| fc_* metrics | τ-independent | No | Baseline — no LRG involved |",
        "",
    ]

    OUT.write_text("\n".join(lines))
    print(f"Task 1 complete: {OUT}")


if __name__ == "__main__":
    main()
