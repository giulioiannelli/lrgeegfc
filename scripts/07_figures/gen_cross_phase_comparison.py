#!/usr/bin/env python3
"""Generate cross-phase structural comparison figures for all patients and bands.

Produces a 3-row figure per (patient, band):
  Row 1: MSC networks with LRG community partition (4 phases)
  Row 2: Dendrograms with colored branches (log scale)
  Row 3: Multiscale ARI curves + MDS phase similarity embedding

Patients excluded: Pat_06, Pat_07 (missing phases).

Output: data/figures/cross_phase_comparison/{Patient}/fig_cross_phase_{band}.pdf
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, to_tree
from scipy.optimize import linear_sum_assignment
from sklearn.manifold import MDS
from sklearn.metrics import adjusted_rand_score

# ---------------------------------------------------------------------------
# Project imports
# ---------------------------------------------------------------------------
from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config import BRAIN_BANDS_NAMES, PHASE_LABELS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import LRG_CACHE, MSC_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.msc import load_msc_matrix

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_08"]
BANDS = BRAIN_BANDS_NAMES  # delta, theta, alpha, beta, low_gamma, high_gamma
PHASES = list(PHASE_LABELS)

OUTPUT_ROOT = FIGURES_ROOT / "cross_phase_comparison"

N_THRESHOLDS = 100
N_COMMUNITIES = 8
LAYOUT_K = 0.1
LAYOUT_SEED = 42
EDGE_POWER = 2.0
MAX_WIDTH = 4.0

# 12 distinct colors (no yellow) for cluster matching
CLUSTER_COLORS_HEX = [
    "#e41a1c",  # red
    "#377eb8",  # blue
    "#4daf4a",  # green
    "#984ea3",  # purple
    "#ff7f00",  # orange
    "#a65628",  # brown
    "#f781bf",  # pink
    "#17becf",  # cyan
    "#8dd3c7",  # teal
    "#bebada",  # lavender
    "#fb8072",  # salmon
    "#80b1d3",  # light blue
]


# ===================================================================
# Cluster matching (Hungarian on Jaccard similarity)
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
    """Return matched node colors, cluster assignments, and color maps."""
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
                [
                    (tc, rc, j)
                    for tc, (rc, j) in mapping.items()
                    if rc is not None and j > 0.1
                ],
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
    tree_root, tree_nodes = to_tree(linkage_matrix, rd=True)
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
# Multiscale ARI
# ===================================================================
def compute_multiscale_ari(linkage1, linkage2, n_nodes, n_thresholds=100):
    max_h1 = linkage1[:, 2].max()
    max_h2 = linkage2[:, 2].max()
    thresholds = np.linspace(0.01, 1.0, n_thresholds)
    ari_vals = []
    for t in thresholds:
        l1 = fcluster(linkage1, t=t * max_h1, criterion="distance")
        l2 = fcluster(linkage2, t=t * max_h2, criterion="distance")
        ari_vals.append(adjusted_rand_score(l1, l2))
    return thresholds, np.array(ari_vals)


def compute_structural_distance(linkage1, linkage2, n_nodes, n_thresholds=100):
    thresholds, ari = compute_multiscale_ari(linkage1, linkage2, n_nodes, n_thresholds)
    d = np.trapz(1 - ari, thresholds) / (thresholds[-1] - thresholds[0])
    return d


# ===================================================================
# Main figure function
# ===================================================================
def create_full_comparison_figure(
    lrg_results,
    msc_matrices,
    phases,
    band,
    patient,
    ari_curves,
    distances,
):
    """3-row figure: networks, dendrograms, ARI + MDS."""
    fig = plt.figure(figsize=(24, 18))
    gs = gridspec.GridSpec(
        3, 4, figure=fig, height_ratios=[1, 1, 1.2], hspace=0.35, wspace=0.25
    )
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

        ax = fig.add_subplot(gs[0, i])
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

        ax = fig.add_subplot(gs[1, i])
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

    # ---- Row 3 left: ARI curves ----
    ax_ari = fig.add_subplot(gs[2, 0:2])
    colors = plt.colormaps["tab10"](np.linspace(0, 1, len(ari_curves)))
    for idx, (pair_key, (thresholds, ari_values)) in enumerate(ari_curves.items()):
        ax_ari.plot(
            thresholds, ari_values, color=colors[idx], linewidth=2,
            label=f"{pair_key} (D={distances[pair_key]:.3f})",
        )
    ax_ari.set_xlabel("Normalized Threshold (h)", fontsize=11)
    ax_ari.set_ylabel("Adjusted Rand Index", fontsize=11)
    ax_ari.set_title("Multiscale Partition Similarity", fontsize=13, fontweight="bold")
    ax_ari.legend(fontsize=8, loc="lower right")
    ax_ari.grid(alpha=0.3)
    ax_ari.set_xlim(0, 1)
    ax_ari.set_ylim(-0.1, 1.05)

    # ---- Row 3 right: MDS embedding ----
    ax_mds = fig.add_subplot(gs[2, 2:4])
    dist_matrix = np.zeros((n_phases, n_phases))
    for i, pi in enumerate(phases):
        for j, pj in enumerate(phases):
            if i < j:
                dist_matrix[i, j] = distances[f"{pi} vs {pj}"]
                dist_matrix[j, i] = dist_matrix[i, j]

    mds = MDS(
        n_components=2, dissimilarity="precomputed",
        random_state=42, normalized_stress="auto",
    )
    coords = mds.fit_transform(dist_matrix)
    cmap_sim = plt.colormaps["RdYlGn"]
    phase_cols = plt.colormaps["Set2"](np.linspace(0, 1, n_phases))

    for i in range(n_phases):
        for j in range(i + 1, n_phases):
            sim = 1 - dist_matrix[i, j]
            ax_mds.plot(
                [coords[i, 0], coords[j, 0]], [coords[i, 1], coords[j, 1]],
                color=cmap_sim(sim), linewidth=max(1.0, sim * 10), zorder=1,
            )
    for i, (phase, c) in enumerate(zip(phases, coords)):
        ax_mds.scatter(c[0], c[1], s=600, c=[phase_cols[i]],
                       edgecolors="black", linewidths=2, zorder=5)
        ax_mds.annotate(phase, c, fontsize=11, fontweight="bold",
                        ha="center", va="center", zorder=10)

    ax_mds.set_title("Phase Structural Similarity\n(MDS from distances)",
                     fontsize=13, fontweight="bold")
    ax_mds.axis("equal")
    ax_mds.axis("off")
    sm = plt.cm.ScalarMappable(cmap=cmap_sim, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax_mds, orientation="vertical", fraction=0.046, pad=0.04)
    cbar.set_label("Similarity (1-D)", fontsize=9)

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    fig.suptitle(
        f"{patient} — {band_tex} band | Cross-Phase Structural Reorganization",
        fontsize=18, fontweight="bold", y=0.99,
    )
    return fig


# ===================================================================
# Driver
# ===================================================================
def generate_for_patient_band(patient, band, verbose=True):
    """Load data, compute ARI, produce figure, save."""
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

    # Compute multiscale ARI for all phase pairs
    n_nodes = lrg_results[PHASES[0]].n_nodes
    ari_curves = {}
    distances = {}
    for pi, pj in combinations(PHASES, 2):
        thresholds, ari_values = compute_multiscale_ari(
            lrg_results[pi].linkage_matrix,
            lrg_results[pj].linkage_matrix,
            n_nodes, N_THRESHOLDS,
        )
        d = compute_structural_distance(
            lrg_results[pi].linkage_matrix,
            lrg_results[pj].linkage_matrix,
            n_nodes, N_THRESHOLDS,
        )
        key = f"{pi} vs {pj}"
        ari_curves[key] = (thresholds, ari_values)
        distances[key] = d
        if verbose:
            print(f"  {key}: D = {d:.4f}")

    # Create figure
    fig = create_full_comparison_figure(
        lrg_results, msc_matrices, PHASES, band, patient,
        ari_curves, distances,
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
    parser.add_argument(
        "--patients", nargs="*", default=None,
        help="Patient IDs (default: all valid)",
    )
    parser.add_argument(
        "--bands", nargs="*", default=None,
        help="Frequency bands (default: all 6)",
    )
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args()

    patients = args.patients or PATIENTS
    bands = args.bands or BANDS
    verbose = not args.quiet

    print(f"Generating cross-phase comparison figures")
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
