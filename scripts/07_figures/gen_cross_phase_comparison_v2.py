#!/usr/bin/env python3
"""Generate cross-phase structural comparison figures (v2: FM + RF + Cophenetic).

Produces a 3-row figure per (patient, band):
  Row 1: MSC networks with LRG community partition (4 phases)
  Row 2: Dendrograms with colored branches (log scale)
  Row 3: Multiscale Fowlkes-Mallows | Robinson-Foulds heatmap | Cophenetic corr heatmap

Patients excluded: Pat_06, Pat_07 (missing phases).

Output: data/figures/cross_phase_comparison_v2/{Patient}/fig_cross_phase_{band}.pdf
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.cluster.hierarchy import cophenet, dendrogram, fcluster, to_tree
from scipy.optimize import linear_sum_assignment
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# Project imports
# ---------------------------------------------------------------------------
from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config import BRAIN_BANDS_NAMES, PHASE_LABELS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import LRG_CACHE, MSC_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrgsglib.utils.basic.linalg import (
    tree_cophenetic_correlation,
    tree_fowlkes_mallows_index,
    tree_robinson_foulds_distance,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_08"]
BANDS = BRAIN_BANDS_NAMES
PHASES = list(PHASE_LABELS)

OUTPUT_ROOT = FIGURES_ROOT / "cross_phase_comparison_v2"

N_COMMUNITIES = 8
LAYOUT_K = 0.1
LAYOUT_SEED = 42
EDGE_POWER = 2.0
MAX_WIDTH = 4.0

# 12 distinct colors (no yellow) for cluster matching
CLUSTER_COLORS_HEX = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628",
    "#f781bf", "#17becf", "#8dd3c7", "#bebada", "#fb8072", "#80b1d3",
]


# ===================================================================
# Cluster matching (unchanged from v1)
# ===================================================================
def _get_cluster_assignments(linkage_matrix, n_nodes, n_clusters):
    return fcluster(linkage_matrix, n_clusters, criterion="maxclust")


def _match_clusters_hungarian(labels_ref, labels_target):
    clusters_ref = np.unique(labels_ref)
    clusters_target = np.unique(labels_target)
    n_ref, n_target = len(clusters_ref), len(clusters_target)

    jaccard = np.zeros((n_ref, n_target))
    for i, ca in enumerate(clusters_ref):
        set_a = set(np.where(labels_ref == ca)[0])
        for j, cb in enumerate(clusters_target):
            set_b = set(np.where(labels_target == cb)[0])
            inter = len(set_a & set_b)
            union = len(set_a | set_b)
            jaccard[i, j] = inter / union if union > 0 else 0

    cost = 1 - jaccard
    if n_ref != n_target:
        max_n = max(n_ref, n_target)
        padded = np.ones((max_n, max_n))
        padded[:n_ref, :n_target] = cost
        row_ind, col_ind = linear_sum_assignment(padded)
    else:
        row_ind, col_ind = linear_sum_assignment(cost)

    mapping = {}
    for ti, ri in zip(col_ind[:n_target], row_ind[:n_target]):
        if ti < n_target and ri < n_ref:
            mapping[clusters_target[ti]] = (clusters_ref[ri], jaccard[ri, ti])
        elif ti < n_target:
            mapping[clusters_target[ti]] = (None, 0.0)
    return mapping


def compute_matched_cluster_colors(lrg_results, phases, n_communities):
    ref_phase = phases[0]
    cluster_assignments = {}
    for phase in phases:
        lrg = lrg_results[phase]
        cluster_assignments[phase] = _get_cluster_assignments(
            lrg.linkage_matrix, lrg.n_nodes, n_communities
        )

    ref_labels = cluster_assignments[ref_phase]
    ref_sizes = sorted(
        [(c, np.sum(ref_labels == c)) for c in np.unique(ref_labels)],
        key=lambda x: -x[1],
    )
    ref_c2color = {c: i for i, (c, _) in enumerate(ref_sizes)}

    phase_node_colors = {}
    phase_cluster_color_maps = {}

    for phase in phases:
        labels = cluster_assignments[phase]
        phase_clusters = np.unique(labels)
        if phase == ref_phase:
            c2color = ref_c2color.copy()
        else:
            mapping = _match_clusters_hungarian(ref_labels, labels)
            c2color = {}
            used = set()
            pairs = sorted(
                [(tc, rc, j) for tc, (rc, j) in mapping.items()
                 if rc is not None and j > 0.1],
                key=lambda x: -x[2],
            )
            for tc, rc, _ in pairs:
                idx = ref_c2color[rc]
                if idx not in used:
                    c2color[tc] = idx
                    used.add(idx)
            avail = sorted(set(range(len(CLUSTER_COLORS_HEX))) - used)
            for tc in phase_clusters:
                if tc not in c2color:
                    c2color[tc] = avail.pop(0) if avail else len(used) % len(CLUSTER_COLORS_HEX)
                    used.add(c2color[tc])

        cmap = {cl: CLUSTER_COLORS_HEX[idx] for cl, idx in c2color.items()}
        phase_cluster_color_maps[phase] = cmap
        phase_node_colors[phase] = [
            plt.matplotlib.colors.to_rgba(cmap[labels[i]]) for i in range(len(labels))
        ]
    return phase_node_colors, cluster_assignments, phase_cluster_color_maps


def _get_link_color_func(linkage_matrix, cluster_labels, color_map, n_nodes, threshold):
    tree_root, _ = to_tree(linkage_matrix, rd=True)
    merge_heights = linkage_matrix[:, 2]
    node_clusters: dict[int, set] = {}

    def _collect(node):
        if node.id in node_clusters:
            return node_clusters[node.id]
        if node.is_leaf():
            s = {cluster_labels[node.id]}
        else:
            s = _collect(node.get_left()) | _collect(node.get_right())
        node_clusters[node.id] = s
        return s

    _collect(tree_root)

    def func(node_id):
        if node_id < n_nodes:
            return color_map.get(cluster_labels[node_id], "gray")
        link_idx = node_id - n_nodes
        if merge_heights[link_idx] > threshold:
            return "black"
        clusters = node_clusters.get(node_id, set())
        if len(clusters) == 1:
            return color_map.get(list(clusters)[0], "gray")
        return "black"

    return func


# ===================================================================
# New metrics
# ===================================================================
def compute_multiscale_fm(linkage1, linkage2, n_nodes, k_values=None):
    """Compute Fowlkes-Mallows index at multiple numbers of clusters k."""
    if k_values is None:
        k_values = np.unique(np.geomspace(2, max(2, n_nodes // 2), 60).astype(int))
    fm_vals = []
    for k in k_values:
        l1 = fcluster(linkage1, k, criterion="maxclust")
        l2 = fcluster(linkage2, k, criterion="maxclust")
        # FM index: sqrt(precision * recall) on pair-counting
        n = len(l1)
        tp = fp = fn = 0
        # Vectorized pair counting via contingency
        from sklearn.metrics.cluster import contingency_matrix
        C = contingency_matrix(l1, l2)
        tp = (C * (C - 1)).sum() / 2
        sum_rows = C.sum(axis=1)
        sum_cols = C.sum(axis=0)
        pairs_in_1 = (sum_rows * (sum_rows - 1)).sum() / 2
        pairs_in_2 = (sum_cols * (sum_cols - 1)).sum() / 2
        if pairs_in_1 == 0 or pairs_in_2 == 0:
            fm_vals.append(0.0)
        else:
            fm_vals.append(float(tp / np.sqrt(pairs_in_1 * pairs_in_2)))
    return k_values, np.array(fm_vals)


def compute_fm_integrated_distance(linkage1, linkage2, n_nodes, k_values=None):
    """Integrated FM distance: integral of (1 - FM(k)) over log(k)."""
    k_vals, fm = compute_multiscale_fm(linkage1, linkage2, n_nodes, k_values)
    log_k = np.log2(k_vals)
    d = np.trapz(1 - fm, log_k) / (log_k[-1] - log_k[0])
    return d


# ===================================================================
# Figure
# ===================================================================
def _plot_triangular_heatmap(ax, matrix, labels, title, cmap, vmin, vmax,
                             fmt=".2f", similarity=False):
    """Plot lower-triangular heatmap with annotations."""
    n = len(labels)
    # Mask upper triangle + diagonal
    mask = np.triu(np.ones((n, n), dtype=bool))
    display = np.where(mask, np.nan, matrix)

    im = ax.imshow(display, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.06, shrink=0.85)

    # Annotate cells
    for i in range(n):
        for j in range(n):
            if i > j:
                val = matrix[i, j]
                color = "white" if (val - vmin) / (vmax - vmin + 1e-12) > 0.6 else "black"
                ax.text(j, i, f"{val:{fmt}}", ha="center", va="center",
                        fontsize=11, fontweight="bold", color=color)

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=10)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")


def create_full_comparison_figure(
    lrg_results, msc_matrices, phases, band, patient,
    fm_curves, fm_distances, rf_matrix, coph_matrix,
):
    """3-row figure: networks, dendrograms, FM curves + RF heatmap + Cophenetic heatmap."""
    fig = plt.figure(figsize=(26, 20))

    # Use nested gridspecs for flexible row-3 layout
    gs_top = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[1, 1, 1.3],
                               hspace=0.32, top=0.94, bottom=0.04,
                               left=0.03, right=0.97)

    gs_row0 = gs_top[0].subgridspec(1, 4, wspace=0.15)
    gs_row1 = gs_top[1].subgridspec(1, 4, wspace=0.25)
    gs_row2 = gs_top[2].subgridspec(1, 3, wspace=0.35, width_ratios=[1.4, 1, 1])

    n_phases = len(phases)

    # Shared layout from backbone of first phase
    G_first = nx.from_numpy_array(msc_matrices[phases[0]])
    w_first = np.array([G_first[u][v]["weight"] for u, v in G_first.edges()])
    thr_layout = np.percentile(w_first, 75)
    G_bb = nx.Graph()
    G_bb.add_nodes_from(G_first.nodes())
    for u, v in G_first.edges():
        if G_first[u][v]["weight"] >= thr_layout:
            G_bb.add_edge(u, v, weight=G_first[u][v]["weight"])
    shared_pos = nx.spring_layout(G_bb, seed=LAYOUT_SEED, k=LAYOUT_K, iterations=100)

    # Matched cluster colors
    phase_nc, cluster_asgn, phase_ccm = compute_matched_cluster_colors(
        lrg_results, phases, N_COMMUNITIES
    )
    cmap_edges = plt.colormaps["viridis"]

    # ---- Row 1: networks ----
    for i, phase in enumerate(phases):
        lrg = lrg_results[phase]
        msc = msc_matrices[phase]
        nc = phase_nc[phase]
        labels = cluster_asgn[phase]
        n_cl = len(np.unique(labels))

        ax = fig.add_subplot(gs_row0[0, i])
        G = nx.from_numpy_array(msc)
        edges = list(G.edges(data=True))
        weights = np.array([e[2]["weight"] for e in edges])
        scaled = np.power(weights, EDGE_POWER)
        widths = MAX_WIDTH * scaled
        alphas = np.power(weights, EDGE_POWER)
        ecols = [(*cmap_edges(w)[:3], a) for w, a in zip(weights, alphas)]

        nx.draw_networkx_nodes(
            G, shared_pos, ax=ax, node_size=60, node_color=nc,
            alpha=0.9, edgecolors="white", linewidths=0.3,
        )
        nx.draw_networkx_edges(G, shared_pos, ax=ax, width=widths, edge_color=ecols)
        ax.set_title(f"{phase}\n({n_cl} comm.)", fontsize=13, fontweight="bold")
        ax.axis("off")

    # ---- Row 2: dendrograms ----
    for i, phase in enumerate(phases):
        lrg = lrg_results[phase]
        n_nodes = lrg.n_nodes
        labels = cluster_asgn[phase]
        ccm = phase_ccm[phase]

        merge_h = lrg.linkage_matrix[:, 2]
        cut_idx = n_nodes - N_COMMUNITIES
        if 0 < cut_idx < len(merge_h):
            threshold = (merge_h[cut_idx - 1] + merge_h[cut_idx]) / 2
        elif cut_idx == 0:
            threshold = merge_h[0] / 2
        else:
            threshold = lrg.optimal_threshold

        ax = fig.add_subplot(gs_row1[0, i])
        lcf = _get_link_color_func(lrg.linkage_matrix, labels, ccm, n_nodes, threshold)
        dendrogram(
            lrg.linkage_matrix, ax=ax, orientation="top",
            no_labels=True, link_color_func=lcf,
        )
        ax.axhline(threshold, color="blue", linestyle="--", lw=2)
        heights = lrg.linkage_matrix[:, 2]
        ax.set_yscale("log")
        ax.set_ylim(heights[0] * 0.8, heights[-1] * 1.01)
        ax.set_ylabel(r"$\mathcal{D}$", fontsize=9)
        ax.grid(alpha=0.2)

    # ---- Row 3, panel 1: Multiscale Fowlkes-Mallows ----
    ax_fm = fig.add_subplot(gs_row2[0, 0])
    colors = plt.colormaps["tab10"](np.linspace(0, 1, len(fm_curves)))
    for idx, (pair_key, (k_vals, fm_vals)) in enumerate(fm_curves.items()):
        d = fm_distances[pair_key]
        ax_fm.plot(
            k_vals, fm_vals, color=colors[idx], linewidth=2,
            label=f"{pair_key} (D={d:.3f})",
        )
    ax_fm.set_xlabel("Number of clusters (k)", fontsize=11)
    ax_fm.set_ylabel("Fowlkes-Mallows Index", fontsize=11)
    ax_fm.set_title("Multiscale Fowlkes-Mallows", fontsize=13, fontweight="bold")
    ax_fm.set_xscale("log", base=2)
    ax_fm.legend(fontsize=7, loc="upper right")
    ax_fm.grid(alpha=0.3)
    ax_fm.set_ylim(-0.05, 1.05)

    # ---- Row 3, panel 2: Robinson-Foulds heatmap ----
    ax_rf = fig.add_subplot(gs_row2[0, 1])
    _plot_triangular_heatmap(
        ax_rf, rf_matrix, phases,
        title="Robinson-Foulds Distance",
        cmap="YlOrRd", vmin=0, vmax=1,
    )

    # ---- Row 3, panel 3: Cophenetic correlation heatmap ----
    ax_coph = fig.add_subplot(gs_row2[0, 2])
    _plot_triangular_heatmap(
        ax_coph, coph_matrix, phases,
        title="Cophenetic Correlation",
        cmap="RdYlGn", vmin=0, vmax=1, similarity=True,
    )

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    fig.suptitle(
        f"{patient} — {band_tex} band | Cross-Phase Structural Reorganization",
        fontsize=18, fontweight="bold", y=0.98,
    )
    return fig


# ===================================================================
# Driver
# ===================================================================
def generate_for_patient_band(patient, band, verbose=True):
    if verbose:
        print(f"\n{'='*60}")
        print(f"  {patient} / {band}")
        print(f"{'='*60}")

    # Load LRG + MSC
    lrg_results = {}
    msc_matrices = {}
    for phase in PHASES:
        lrg = load_lrg_result(patient, phase, band, "msc", cache_root=LRG_CACHE)
        msc = load_msc_matrix(patient, phase, band, cache_root=MSC_CACHE, sparsify="none")
        if lrg is None or msc is None:
            print(f"  SKIP: missing data for {phase}")
            return None
        lrg_results[phase] = lrg
        msc_matrices[phase] = msc
        if verbose:
            print(f"  Loaded {phase}: {lrg.n_nodes} nodes")

    n_nodes = lrg_results[PHASES[0]].n_nodes
    n_phases = len(PHASES)

    # Shared k values for FM
    k_values = np.unique(np.geomspace(2, max(2, n_nodes // 2), 60).astype(int))

    # Compute all metrics
    fm_curves = {}
    fm_distances = {}
    rf_matrix = np.zeros((n_phases, n_phases))
    coph_matrix = np.eye(n_phases)  # diagonal = 1

    for i, pi in enumerate(PHASES):
        Zi = lrg_results[pi].linkage_matrix
        for j, pj in enumerate(PHASES):
            if i >= j:
                continue
            Zj = lrg_results[pj].linkage_matrix
            key = f"{pi} vs {pj}"

            # Fowlkes-Mallows curve
            k_vals, fm_vals = compute_multiscale_fm(Zi, Zj, n_nodes, k_values)
            fm_curves[key] = (k_vals, fm_vals)
            fm_distances[key] = compute_fm_integrated_distance(Zi, Zj, n_nodes, k_values)

            # Robinson-Foulds
            rf = tree_robinson_foulds_distance(Zi, Zj, normalized=True)
            rf_matrix[i, j] = rf
            rf_matrix[j, i] = rf

            # Cophenetic correlation
            cc = tree_cophenetic_correlation(Zi, Zj)
            coph_matrix[i, j] = cc
            coph_matrix[j, i] = cc

            if verbose:
                print(f"  {key}: FM_D={fm_distances[key]:.4f}  RF={rf:.4f}  Coph={cc:.4f}")

    # Create figure
    fig = create_full_comparison_figure(
        lrg_results, msc_matrices, PHASES, band, patient,
        fm_curves, fm_distances, rf_matrix, coph_matrix,
    )

    out_dir = OUTPUT_ROOT / patient
    out_dir.mkdir(parents=True, exist_ok=True)

    pdf_path = out_dir / f"fig_cross_phase_{band}.pdf"
    png_path = out_dir / f"fig_cross_phase_{band}.png"
    fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
    fig.savefig(png_path, dpi=100, bbox_inches="tight")
    plt.close(fig)

    if verbose:
        print(f"  Saved: {pdf_path}")
    return pdf_path


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patients", nargs="*", default=None)
    parser.add_argument("--bands", nargs="*", default=None)
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args()

    patients = args.patients or PATIENTS
    bands = args.bands or BANDS
    verbose = not args.quiet

    print("Generating cross-phase comparison figures (v2: FM + RF + Cophenetic)")
    print(f"  Patients: {patients}")
    print(f"  Bands:    {bands}")
    print(f"  Output:   {OUTPUT_ROOT}")

    generated = []
    for patient in patients:
        for band in bands:
            path = generate_for_patient_band(patient, band, verbose=verbose)
            if path is not None:
                generated.append(path)

    print(f"\n{'='*60}")
    print(f"Done. Generated {len(generated)} figures.")
    for p in generated:
        print(f"  {p}")


if __name__ == "__main__":
    main()
