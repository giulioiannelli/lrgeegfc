#!/usr/bin/env python3
"""Diagnostic figures for 3 outlier cases in mean_VI phase comparison.

Outliers investigated:
  1. Pat_02 theta  — TL<->TT mean_VI=0.877 (worst pair, expected best)
  2. Pat_08 delta  — Pre<->Post mean_VI=0.144 (most similar), TL<->TT=0.527
  3. Pat_03 high_gamma — Pre<->Post mean_VI=0.964 (extreme rest reorganisation)

For each outlier case the script produces a multi-panel PDF:
  Row A: side-by-side dendrograms for all 4 phases (coloured at k=3 and k=6)
  Row B: VI(k) curves for key phase pairs (TL<->TT, Pre<->Post)
  Row C: node-switching count at each k for the same pairs

Output directory:
  data/figures/metric_exploration/outlier_diagnostics/
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import mutual_info_score

# ── Project imports ────────────────────────────────────────────────────
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.paths import LRG_CACHE, FIGURES_ROOT

# ── Constants ──────────────────────────────────────────────────────────
CACHE_ROOT = LRG_CACHE
OUTDIR = FIGURES_ROOT / "metric_exploration" / "outlier_diagnostics"
OUTDIR.mkdir(parents=True, exist_ok=True)

PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHASE_SHORT = {
    "rest_pre": "Pre", "task_learn": "TL", "task_test": "TT", "rest_post": "Post"
}

OUTLIERS = [
    {"patient": "Pat_02", "band": "theta",
     "title": "Pat_02 theta -- TL<->TT mean_VI=0.877 (expected similar)"},
    {"patient": "Pat_08", "band": "delta",
     "title": "Pat_08 delta -- Pre<->Post=0.144 (most similar), TL<->TT=0.527"},
    {"patient": "Pat_03", "band": "high_gamma",
     "title": "Pat_03 high_gamma -- Pre<->Post mean_VI=0.964 (extreme rest reorg)"},
]

K_RANGE = np.arange(2, 21)

# Key pairs to highlight in VI(k) and switching plots
KEY_PAIRS = [
    ("task_learn", "task_test"),
    ("rest_pre", "rest_post"),
]
# Also include all 6 pairs for full context
ALL_PAIRS = list(combinations(PHASES, 2))

PAIR_COLORS = {
    ("task_learn", "task_test"): "#2196F3",   # blue  - within task
    ("rest_pre", "rest_post"):       "#4CAF50",   # green - within rest
    ("rest_pre", "task_learn"):    "#9E9E9E",   # grey
    ("rest_pre", "task_test"):     "#BDBDBD",   # light grey
    ("task_learn", "rest_post"):   "#FF9800",   # orange
    ("task_test", "rest_post"):    "#E91E63",   # pink
}
PAIR_LABELS = {
    ("task_learn", "task_test"): "TL<->TT",
    ("rest_pre", "rest_post"):       "Pre<->Post",
    ("rest_pre", "task_learn"):    "Pre<->TL",
    ("rest_pre", "task_test"):     "Pre<->TT",
    ("task_learn", "rest_post"):   "TL<->Post",
    ("task_test", "rest_post"):    "TT<->Post",
}

# Dendrogram colour palette (for cluster colouring)
CLUSTER_COLORS = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#42d4f4", "#f032e6", "#bfef45", "#fabebe", "#469990",
    "#e6beff", "#9A6324", "#ffe119", "#800000", "#aaffc3",
    "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#000000",
]


# ── Helpers ────────────────────────────────────────────────────────────
def variation_of_information(labels1: np.ndarray, labels2: np.ndarray) -> float:
    """Variation of Information between two label vectors."""
    n = len(labels1)
    if n == 0:
        return 0.0
    h1 = -sum(
        (np.sum(labels1 == c) / n) * np.log(np.sum(labels1 == c) / n)
        for c in np.unique(labels1)
    )
    h2 = -sum(
        (np.sum(labels2 == c) / n) * np.log(np.sum(labels2 == c) / n)
        for c in np.unique(labels2)
    )
    mi = mutual_info_score(labels1, labels2)
    return h1 + h2 - 2 * mi


def get_partitions(linkage_matrix: np.ndarray, k_range: np.ndarray) -> dict:
    """Return {k: labels} for each k in k_range using fcluster."""
    return {k: fcluster(linkage_matrix, t=k, criterion="maxclust") for k in k_range}


def count_node_switches(labels1: np.ndarray, labels2: np.ndarray) -> int:
    """Count how many nodes changed cluster between two aligned label vectors.

    Since cluster IDs are arbitrary, we use a greedy matching to align them
    and then count mismatches.
    """
    from scipy.optimize import linear_sum_assignment

    # Build cost matrix (negative overlap) for Hungarian algorithm
    u1 = np.unique(labels1)
    u2 = np.unique(labels2)
    cost = np.zeros((len(u1), len(u2)))
    for i, c1 in enumerate(u1):
        for j, c2 in enumerate(u2):
            cost[i, j] = -np.sum((labels1 == c1) & (labels2 == c2))
    row_ind, col_ind = linear_sum_assignment(cost)
    # Build mapping
    mapping = {u2[col_ind[i]]: u1[row_ind[i]] for i in range(len(row_ind))}
    # Remap labels2
    labels2_mapped = np.array([mapping.get(l, -1) for l in labels2])
    return int(np.sum(labels1 != labels2_mapped))


def color_dendrogram_by_partition(ax, Z, labels_k, title, k_val, n_nodes):
    """Plot a dendrogram on ax, colouring leaves by their cluster in labels_k."""
    # Compute dendrogram without plotting to get leaf order
    R = dendrogram(Z, no_plot=True, count_sort="descending")
    leaf_order = R["leaves"]

    # Map each cluster to a colour
    unique_clusters = np.unique(labels_k)
    cluster_to_color = {}
    for i, c in enumerate(sorted(unique_clusters)):
        cluster_to_color[c] = CLUSTER_COLORS[i % len(CLUSTER_COLORS)]

    # Build leaf-colour mapping (by original index)
    leaf_color_map = {}
    for idx in range(n_nodes):
        leaf_color_map[idx] = cluster_to_color[labels_k[idx]]

    # Create link-colour function: colour a link by its descendant leaves
    def link_color_func(node_id):
        """Return colour if all descendant leaves belong to same cluster."""
        # Get all leaves under this node
        leaves = _get_leaves(Z, node_id, n_nodes)
        clusters_here = set(labels_k[l] for l in leaves)
        if len(clusters_here) == 1:
            return cluster_to_color[clusters_here.pop()]
        return "#888888"  # mixed: grey

    R2 = dendrogram(
        Z, ax=ax, count_sort="descending",
        link_color_func=link_color_func,
        no_labels=True,
    )

    # Colour the x-tick positions by leaf cluster
    leaf_order2 = R2["leaves"]
    xlbls = ax.get_xticks()

    ax.set_title(title, fontsize=9, fontweight="bold")
    ax.set_ylabel("Ultrametric distance", fontsize=7)
    ax.tick_params(axis="both", labelsize=6)

    # Add coloured bars at bottom to show cluster membership
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    bar_height = (ylim[1] - ylim[0]) * 0.03
    for i, leaf_idx in enumerate(leaf_order2):
        x_pos = 5 + i * 10  # default dendrogram spacing
        color = cluster_to_color[labels_k[leaf_idx]]
        ax.bar(x_pos, bar_height, bottom=ylim[0], width=9,
               color=color, edgecolor="none", alpha=0.7)


def _get_leaves(Z, node_id, n_leaves):
    """Recursively get all leaf indices under a node in the linkage tree."""
    if node_id < n_leaves:
        return [node_id]
    row = int(node_id - n_leaves)
    left = int(Z[row, 0])
    right = int(Z[row, 1])
    return _get_leaves(Z, left, n_leaves) + _get_leaves(Z, right, n_leaves)


# ── Main ───────────────────────────────────────────────────────────────
def make_outlier_figure(case: dict) -> Path:
    """Create one multi-panel diagnostic figure for a single outlier case."""
    patient = case["patient"]
    band = case["band"]
    title = case["title"]

    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")

    # ---- Load LRG results for all 4 phases ----
    results = {}
    partitions = {}
    for phase in PHASES:
        r = load_lrg_result(patient, phase, band, fc_method="msc", cache_root=CACHE_ROOT)
        if r is None:
            print(f"  WARNING: No cache for {patient} {phase} {band}")
            return None
        results[phase] = r
        partitions[phase] = get_partitions(r.linkage_matrix, K_RANGE)
        print(f"  Loaded {phase}: n_nodes={r.n_nodes}, "
              f"optimal_threshold={r.optimal_threshold:.4f}")

    # Check all phases have same number of nodes
    n_nodes_list = [results[p].n_nodes for p in PHASES]
    print(f"  Node counts: {dict(zip(PHASES, n_nodes_list))}")
    # Use minimum for safe indexing
    min_nodes = min(n_nodes_list)

    # ---- Compute VI at each k for all pairs ----
    vi_curves = {}
    switch_curves = {}
    for p1, p2 in ALL_PAIRS:
        vi_vals = []
        switch_vals = []
        for k in K_RANGE:
            labels1 = partitions[p1][k][:min_nodes]
            labels2 = partitions[p2][k][:min_nodes]
            vi_vals.append(variation_of_information(labels1, labels2))
            switch_vals.append(count_node_switches(labels1, labels2))
        vi_curves[(p1, p2)] = np.array(vi_vals)
        switch_curves[(p1, p2)] = np.array(switch_vals)

    # Print summary
    for pair in ALL_PAIRS:
        mean_vi = np.mean(vi_curves[pair])
        label = PAIR_LABELS[pair]
        print(f"  {label}: mean_VI={mean_vi:.3f}")

    # ---- Create figure ----
    # Layout: 4 rows
    #   Row 0: dendrograms (k=3 colouring), 4 columns
    #   Row 1: dendrograms (k=6 colouring), 4 columns
    #   Row 2: VI(k) curves (left) + switch count (right)
    #   Row 3: VI at each k as heatmap for all 6 pairs
    fig = plt.figure(figsize=(20, 22))
    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.98)

    outer_gs = gridspec.GridSpec(4, 1, figure=fig, height_ratios=[1, 1, 1, 0.7],
                                  hspace=0.35, top=0.95, bottom=0.04)

    # ---- Row 0: Dendrograms coloured at k=3 ----
    gs_dend3 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer_gs[0],
                                                 wspace=0.25)
    for i, phase in enumerate(PHASES):
        ax = fig.add_subplot(gs_dend3[i])
        r = results[phase]
        labels_k3 = partitions[phase][3]
        color_dendrogram_by_partition(
            ax, r.linkage_matrix, labels_k3,
            f"{PHASE_SHORT[phase]} (k=3)", k_val=3, n_nodes=r.n_nodes
        )
        if i == 0:
            ax.set_ylabel("Ultrametric dist\n(k=3 colouring)", fontsize=8)

    # ---- Row 1: Dendrograms coloured at k=6 ----
    gs_dend6 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer_gs[1],
                                                 wspace=0.25)
    for i, phase in enumerate(PHASES):
        ax = fig.add_subplot(gs_dend6[i])
        r = results[phase]
        labels_k6 = partitions[phase][6]
        color_dendrogram_by_partition(
            ax, r.linkage_matrix, labels_k6,
            f"{PHASE_SHORT[phase]} (k=6)", k_val=6, n_nodes=r.n_nodes
        )
        if i == 0:
            ax.set_ylabel("Ultrametric dist\n(k=6 colouring)", fontsize=8)

    # ---- Row 2: VI(k) curves + switch counts ----
    gs_curves = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer_gs[2],
                                                  wspace=0.30)

    # --- Left: VI(k) ---
    ax_vi = fig.add_subplot(gs_curves[0])
    for pair in ALL_PAIRS:
        is_key = pair in KEY_PAIRS
        lw = 2.5 if is_key else 1.0
        alpha = 1.0 if is_key else 0.4
        zorder = 10 if is_key else 1
        ax_vi.plot(K_RANGE, vi_curves[pair],
                   color=PAIR_COLORS[pair],
                   linewidth=lw, alpha=alpha, zorder=zorder,
                   label=PAIR_LABELS[pair],
                   marker="o" if is_key else None,
                   markersize=4)

    # Add mean_VI as horizontal lines for key pairs
    for pair in KEY_PAIRS:
        mean_val = np.mean(vi_curves[pair])
        ax_vi.axhline(mean_val, color=PAIR_COLORS[pair], linestyle="--",
                      alpha=0.5, linewidth=1)
        ax_vi.text(K_RANGE[-1] + 0.3, mean_val,
                   f"mean={mean_val:.3f}", color=PAIR_COLORS[pair],
                   fontsize=7, va="center")

    ax_vi.set_xlabel("Number of clusters (k)", fontsize=10)
    ax_vi.set_ylabel("Variation of Information", fontsize=10)
    ax_vi.set_title("VI at each scale k", fontsize=11, fontweight="bold")
    ax_vi.legend(fontsize=7, loc="upper left", ncol=2)
    ax_vi.set_xlim(K_RANGE[0] - 0.5, K_RANGE[-1] + 2)
    ax_vi.grid(True, alpha=0.3)

    # --- Right: Node switches ---
    ax_sw = fig.add_subplot(gs_curves[1])
    for pair in ALL_PAIRS:
        is_key = pair in KEY_PAIRS
        lw = 2.5 if is_key else 1.0
        alpha = 1.0 if is_key else 0.4
        zorder = 10 if is_key else 1
        pct = switch_curves[pair] / min_nodes * 100
        ax_sw.plot(K_RANGE, pct,
                   color=PAIR_COLORS[pair],
                   linewidth=lw, alpha=alpha, zorder=zorder,
                   label=PAIR_LABELS[pair],
                   marker="o" if is_key else None,
                   markersize=4)

    ax_sw.set_xlabel("Number of clusters (k)", fontsize=10)
    ax_sw.set_ylabel("Nodes switching cluster (%)", fontsize=10)
    ax_sw.set_title("Partition instability at each scale", fontsize=11, fontweight="bold")
    ax_sw.legend(fontsize=7, loc="upper left", ncol=2)
    ax_sw.set_xlim(K_RANGE[0] - 0.5, K_RANGE[-1] + 0.5)
    ax_sw.grid(True, alpha=0.3)

    # ---- Row 3: VI heatmap (pairs x k) ----
    gs_hm = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer_gs[3])
    ax_hm = fig.add_subplot(gs_hm[0])

    vi_matrix = np.zeros((len(ALL_PAIRS), len(K_RANGE)))
    for i, pair in enumerate(ALL_PAIRS):
        vi_matrix[i, :] = vi_curves[pair]

    im = ax_hm.imshow(vi_matrix, aspect="auto", cmap="YlOrRd",
                       interpolation="nearest")
    ax_hm.set_xticks(range(len(K_RANGE)))
    ax_hm.set_xticklabels([str(k) for k in K_RANGE], fontsize=8)
    ax_hm.set_yticks(range(len(ALL_PAIRS)))
    ax_hm.set_yticklabels([PAIR_LABELS[p] for p in ALL_PAIRS], fontsize=9)
    ax_hm.set_xlabel("Number of clusters (k)", fontsize=10)
    ax_hm.set_title("VI heatmap: all pairs x scales", fontsize=11, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax_hm, shrink=0.8, label="VI")

    # Annotate cells with values
    for i in range(vi_matrix.shape[0]):
        for j in range(vi_matrix.shape[1]):
            val = vi_matrix[i, j]
            text_color = "white" if val > vi_matrix.max() * 0.6 else "black"
            ax_hm.text(j, i, f"{val:.2f}", ha="center", va="center",
                      fontsize=5.5, color=text_color)

    # ---- Save ----
    out_path = OUTDIR / f"outlier_{patient}_{band}.pdf"
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"\n  Saved: {out_path}")
    return out_path


# ── Run all 3 outlier cases ──────────────────────────────────────────
if __name__ == "__main__":
    for case in OUTLIERS:
        make_outlier_figure(case)
    print("\nDone. All figures saved to:", OUTDIR)
