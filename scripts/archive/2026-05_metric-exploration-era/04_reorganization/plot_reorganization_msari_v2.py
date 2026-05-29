#!/usr/bin/env python3
"""Improved MS-ARI reorganization analysis — all patients, all scales.

Changes from v1:
  - ARI clipped to [0,1] before integrating → D ∈ [0,1] guaranteed
  - Default integral from h=0.5 (skip noisy small-scale partitions)
  - Three scale stripes: micro [0, 0.33], meso [0.33, 0.66], macro [0.66, 1.0]
  - Show pairwise distances per (patient, band), not a single collapsed RI
  - Key pairs: rest_pre↔rest_post, task_learn↔task_test, task_test↔rest_post
  - Diagnostic panel for Pat_08 theta (outlier investigation)
  - All 6 patients (Pat_06: rest only, Pat_07: 3 phases)
  - PDF only, no PNG

Output: data/figures/metric_exploration/
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, to_tree
from scipy.optimize import linear_sum_assignment
from sklearn.manifold import MDS
from sklearn.metrics import adjusted_rand_score

from lrg_eegfc.config import BRAIN_BANDS_NAMES, PHASE_LABELS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.msc import load_msc_matrix

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
PATIENTS = list(PATIENT_PHASES.keys())
BANDS = BRAIN_BANDS_NAMES
BAND_TEX = [BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS]

USE_CREMA = "--crema" in sys.argv

from lrg_eegfc.config.paths import LRG_CACHE as _LRG_CACHE, LRG_CREMA_CACHE, MSC_CACHE
LRG_CACHE = LRG_CREMA_CACHE if USE_CREMA else _LRG_CACHE
OUTPUT_DIR = FIGURES_ROOT / "crema_metric_exploration" if USE_CREMA else FIGURES_ROOT / "metric_exploration"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

_MSC_LOAD_KW = (
    dict(sparsify="ecm_adaptive", nperseg=4096, ecm_alpha_min=0.01,
         ecm_alpha_max=0.50, ecm_n_ensemble=100, ecm_weight_scale=1000)
    if USE_CREMA else dict(sparsify="none")
)

# Key phase pairs for analysis
KEY_PAIRS = [
    ("rest_pre", "rest_post"),       # rest stability
    ("task_learn", "task_test"),  # task stability
    ("task_test", "rest_post"),     # post-task reorganization
]
KEY_PAIR_LABELS = {
    ("rest_pre", "rest_post"): "rest_pre ↔ rest_post",
    ("task_learn", "task_test"): "task_learn ↔ task_test",
    ("task_test", "rest_post"): "task_test ↔ rest_post",
}

# Scale stripes
SCALE_STRIPES = {
    "micro": (0.0, 0.33),
    "meso": (0.33, 0.66),
    "macro": (0.66, 1.0),
}

# Diagnostic figure constants
N_COMMUNITIES = 8
LAYOUT_K = 0.1
LAYOUT_SEED = 42
EDGE_POWER = 2.0
MAX_WIDTH = 4.0
CLUSTER_COLORS_HEX = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
    "#ff7f00", "#a65628", "#f781bf", "#17becf",
    "#8dd3c7", "#bebada", "#fb8072", "#80b1d3",
]

PAT_COLORS = {
    "Pat_02": "#66c2a5", "Pat_03": "#fc8d62", "Pat_05": "#8da0cb",
    "Pat_06": "#a6d854", "Pat_07": "#ffd92f", "Pat_08": "#e78ac3",
}
PAT_MARKERS = {
    "Pat_02": "o", "Pat_03": "s", "Pat_05": "D",
    "Pat_06": "^", "Pat_07": "v", "Pat_08": "P",
}

N_THRESHOLDS = 200  # Resolution of ARI curve


# ===================================================================
# Core: Multiscale ARI with clipping and flexible bounds
# ===================================================================
def compute_ari_curve(Z1, Z2, n_thresholds=N_THRESHOLDS):
    """Compute ARI(h) for h in [0.01, 1.0].

    Returns (thresholds, ari_values) where ari_values are clipped to [0, 1].
    """
    max_h1 = Z1[:, 2].max()
    max_h2 = Z2[:, 2].max()
    thresholds = np.linspace(0.01, 1.0, n_thresholds)
    ari_vals = np.empty(n_thresholds)
    for i, t in enumerate(thresholds):
        l1 = fcluster(Z1, t=t * max_h1, criterion="distance")
        l2 = fcluster(Z2, t=t * max_h2, criterion="distance")
        ari_vals[i] = adjusted_rand_score(l1, l2)
    return thresholds, ari_vals


def msari_distance(thresholds, ari_vals, h_lo=0.5, h_hi=1.0):
    """Compute MS-ARI distance over [h_lo, h_hi] with ARI clipped to [0,1].

    D = (1/(h_hi - h_lo)) · ∫_{h_lo}^{h_hi} (1 - max(0, ARI(h))) dh

    Returns D ∈ [0, 1].
    """
    mask = (thresholds >= h_lo) & (thresholds <= h_hi)
    t = thresholds[mask]
    a = ari_vals[mask]
    if len(t) < 2:
        return np.nan
    a_clipped = np.clip(a, 0.0, 1.0)
    integrand = 1.0 - a_clipped
    return float(np.trapz(integrand, t) / (t[-1] - t[0]))


# ===================================================================
# Data loading
# ===================================================================
def load_all_lrg():
    """Load all LRG results keyed by (patient, phase, band)."""
    data = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for phase in phases:
                lrg = load_lrg_result(pat, phase, band, "msc", cache_root=LRG_CACHE)
                if lrg is not None:
                    data[(pat, phase, band)] = lrg
    return data


def compute_all_ari_curves(lrg_data):
    """Pre-compute ARI curves for all patient/band/pair combinations.

    Returns dict keyed by (patient, band, (phase_a, phase_b)) → (thresholds, ari_vals).
    """
    curves = {}
    for pat, phases in PATIENT_PHASES.items():
        pairs = list(combinations(phases, 2))
        for band in BANDS:
            for p1, p2 in pairs:
                k1 = (pat, p1, band)
                k2 = (pat, p2, band)
                if k1 in lrg_data and k2 in lrg_data:
                    t, a = compute_ari_curve(
                        lrg_data[k1].linkage_matrix,
                        lrg_data[k2].linkage_matrix,
                    )
                    curves[(pat, band, (p1, p2))] = (t, a)
                    print(f"  {pat}/{band}/{p1}↔{p2}: done")
    return curves


def compute_distances_from_curves(curves, h_lo=0.5, h_hi=1.0):
    """Compute MS-ARI distances from pre-computed ARI curves."""
    dists = {}
    for key, (t, a) in curves.items():
        dists[key] = msari_distance(t, a, h_lo, h_hi)
    return dists


# ===================================================================
# Cluster matching (for diagnostic panel)
# ===================================================================
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
    max_n = max(n_ref, n_target)
    padded = np.ones((max_n, max_n))
    padded[:n_ref, :n_target] = cost
    row_ind, col_ind = linear_sum_assignment(padded)
    mapping = {}
    for ti, ri in zip(col_ind[:n_target], row_ind[:n_target]):
        if ti < n_target and ri < n_ref:
            mapping[clusters_target[ti]] = (clusters_ref[ri], jaccard[ri, ti])
        elif ti < n_target:
            mapping[clusters_target[ti]] = (None, 0.0)
    return mapping


def compute_matched_cluster_colors(lrg_results, phases):
    ref_phase = phases[0]
    cluster_assignments = {}
    for phase in phases:
        lrg = lrg_results[phase]
        cluster_assignments[phase] = fcluster(
            lrg.linkage_matrix, N_COMMUNITIES, criterion="maxclust"
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
    node_clusters = {}

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
# Figure 1: Diagnostic panel — Pat_08 theta
# ===================================================================
def plot_diagnostic_panel(lrg_data, curves, patient="Pat_08", band="theta"):
    """3-row figure: networks + dendrograms + ARI curves."""
    phases = PATIENT_PHASES[patient]
    n_phases = len(phases)

    # Load MSC matrices
    msc_matrices = {}
    for phase in phases:
        msc = load_msc_matrix(patient, phase, band, cache_root=MSC_CACHE, **_MSC_LOAD_KW)
        msc_matrices[phase] = msc

    lrg_results = {ph: lrg_data[(patient, ph, band)] for ph in phases}

    fig = plt.figure(figsize=(6 * n_phases, 18))
    gs = gridspec.GridSpec(3, n_phases, figure=fig, height_ratios=[1, 1, 1.4],
                           hspace=0.30, wspace=0.25)

    # Shared layout
    G_first = nx.from_numpy_array(msc_matrices[phases[0]])
    w_first = np.array([G_first[u][v]["weight"] for u, v in G_first.edges()])
    thr_layout = np.percentile(w_first, 75)
    G_bb = nx.Graph()
    G_bb.add_nodes_from(G_first.nodes())
    for u, v in G_first.edges():
        if G_first[u][v]["weight"] >= thr_layout:
            G_bb.add_edge(u, v, weight=G_first[u][v]["weight"])
    shared_pos = nx.spring_layout(G_bb, seed=LAYOUT_SEED, k=LAYOUT_K, iterations=100)

    # Cluster matching
    phase_nc, cluster_asgn, phase_ccm = compute_matched_cluster_colors(lrg_results, phases)
    cmap_edges = plt.colormaps["viridis"]

    # Row 1: Networks
    for i, phase in enumerate(phases):
        lrg = lrg_results[phase]
        msc = msc_matrices[phase]
        nc = phase_nc[phase]
        n_cl = len(np.unique(cluster_asgn[phase]))
        ax = fig.add_subplot(gs[0, i])
        G = nx.from_numpy_array(msc)
        edges = list(G.edges(data=True))
        weights = np.array([e[2]["weight"] for e in edges])
        scaled = np.power(weights, EDGE_POWER)
        widths = MAX_WIDTH * scaled
        alphas = np.power(weights, EDGE_POWER)
        ecols = [(*cmap_edges(w)[:3], a) for w, a in zip(weights, alphas)]
        nx.draw_networkx_nodes(G, shared_pos, ax=ax, node_size=60, node_color=nc,
                               alpha=0.9, edgecolors="white", linewidths=0.3)
        nx.draw_networkx_edges(G, shared_pos, ax=ax, width=widths, edge_color=ecols)
        ax.set_title(f"{phase}\n({n_cl} comm.)", fontsize=13, fontweight="bold")
        ax.axis("off")

    # Row 2: Dendrograms
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
        dendrogram(lrg.linkage_matrix, ax=ax, orientation="top",
                    no_labels=True, link_color_func=lcf)
        ax.axhline(threshold, color="blue", linestyle="--", lw=2)
        ax.set_yscale("log")
        ax.set_ylim(merge_h[0] * 0.8, merge_h[-1] * 1.01)
        ax.set_ylabel(r"$\mathcal{D}$", fontsize=9)
        ax.grid(alpha=0.2)

    # Row 3: ARI curves (spanning all columns)
    ax_ari = fig.add_subplot(gs[2, :])
    pairs = list(combinations(phases, 2))
    colors = plt.colormaps["tab10"](np.linspace(0, 1, len(pairs)))

    for idx, (p1, p2) in enumerate(pairs):
        key = (patient, band, (p1, p2))
        if key not in curves:
            continue
        t, a = curves[key]
        d_full = msari_distance(t, a, 0.01, 1.0)
        d_h50 = msari_distance(t, a, 0.5, 1.0)
        ax_ari.plot(t, a, color=colors[idx], linewidth=2,
                    label=f"{p1}↔{p2}  D[0,1]={d_full:.3f}  D[.5,1]={d_h50:.3f}")

    # Shade scale stripes
    stripe_colors = {"micro": "#ff0000", "meso": "#00cc00", "macro": "#0000ff"}
    for name, (lo, hi) in SCALE_STRIPES.items():
        ax_ari.axvspan(lo, hi, alpha=0.08, color=stripe_colors[name],
                       label=f"{name} [{lo:.2f}, {hi:.2f}]")

    # Mark h=0.5 cutoff
    ax_ari.axvline(0.5, color="red", linestyle=":", lw=2, alpha=0.7, label="h=0.5 cutoff")
    ax_ari.axhline(0, color="gray", linestyle="-", lw=0.5, alpha=0.5)

    ax_ari.set_xlabel("Normalized threshold (h)", fontsize=12)
    ax_ari.set_ylabel("Adjusted Rand Index", fontsize=12)
    ax_ari.set_title("Multiscale ARI curves — all phase pairs", fontsize=13, fontweight="bold")
    ax_ari.legend(fontsize=8, loc="lower right", ncol=2)
    ax_ari.grid(alpha=0.3)
    ax_ari.set_xlim(0, 1.02)
    ax_ari.set_ylim(-0.15, 1.05)

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    fig.suptitle(
        f"{patient} — {band_tex} | Diagnostic Panel",
        fontsize=16, fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    path = OUTPUT_DIR / f"diagnostic_{patient}_{band}.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


# ===================================================================
# Figure 2: Pairwise distances (h >= 0.5) — 3 key pairs
# ===================================================================
def plot_pairwise_distances(dists_h50):
    """3 panels, one per key pair. Each shows patient × band distances."""
    fig, axes = plt.subplots(1, 3, figsize=(24, 7))
    x = np.arange(len(BANDS))

    for ax_idx, pair in enumerate(KEY_PAIRS):
        ax = axes[ax_idx]
        pair_label = KEY_PAIR_LABELS[pair]

        # Collect patients that have this pair
        offsets_list = []
        for pat in PATIENTS:
            phases = PATIENT_PHASES[pat]
            if pair[0] in phases and pair[1] in phases:
                offsets_list.append(pat)

        n_pats = len(offsets_list)
        if n_pats == 0:
            ax.set_title(f"{pair_label}\n(no data)")
            continue

        offsets = np.linspace(-0.25, 0.25, n_pats)

        for i, pat in enumerate(offsets_list):
            vals = []
            for band in BANDS:
                key = (pat, band, pair)
                vals.append(dists_h50.get(key, np.nan))
            ax.bar(x + offsets[i], vals, width=0.5 / n_pats,
                   color=PAT_COLORS[pat], label=pat,
                   edgecolor="black", linewidth=0.3, alpha=0.85)

        # Band means
        for j, band in enumerate(BANDS):
            band_vals = [dists_h50.get((pat, band, pair), np.nan)
                         for pat in offsets_list]
            band_vals = [v for v in band_vals if not np.isnan(v)]
            if band_vals:
                m = np.mean(band_vals)
                ax.plot([j - 0.3, j + 0.3], [m, m], color="black", lw=2.5, zorder=6)

        ax.set_xticks(x)
        ax.set_xticklabels(BAND_TEX, fontsize=11)
        ax.set_ylabel("MS-ARI distance (h ≥ 0.5)", fontsize=11)
        ax.set_title(pair_label, fontsize=13, fontweight="bold")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(alpha=0.2, axis="y")
        ax.set_ylim(0, 1.05)

    fig.suptitle(
        "Pairwise MS-ARI Distances — key phase pairs (h ≥ 0.5, ARI clipped to [0,1])",
        fontsize=15, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    path = OUTPUT_DIR / "pairwise_distances_h50.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


# ===================================================================
# Figure 3: Scale stripes comparison
# ===================================================================
def plot_scale_stripes(curves):
    """3 rows (stripes) × 6 cols (bands), each cell is a patient × distance bar."""
    stripe_names = list(SCALE_STRIPES.keys())
    n_stripes = len(stripe_names)
    n_bands = len(BANDS)

    fig, axes = plt.subplots(n_stripes, n_bands, figsize=(4 * n_bands, 4 * n_stripes),
                             sharey=True)

    # Only use patients with ≥3 phases for the key pairs
    for si, stripe_name in enumerate(stripe_names):
        h_lo, h_hi = SCALE_STRIPES[stripe_name]

        for bi, band in enumerate(BANDS):
            ax = axes[si, bi]

            # For each key pair, show grouped bars
            bar_width = 0.25
            pair_offsets = [-bar_width, 0, bar_width]
            pair_colors = ["#2166ac", "#4393c3", "#b2182b"]

            for pi, pair in enumerate(KEY_PAIRS):
                pats_with_pair = [p for p in PATIENTS
                                  if pair[0] in PATIENT_PHASES[p]
                                  and pair[1] in PATIENT_PHASES[p]]
                vals = []
                pat_labels = []
                for pat in pats_with_pair:
                    key = (pat, band, pair)
                    if key in curves:
                        t, a = curves[key]
                        d = msari_distance(t, a, h_lo, h_hi)
                        vals.append(d)
                        pat_labels.append(pat)

                if vals:
                    x_pos = np.arange(len(vals))
                    ax.bar(x_pos + pair_offsets[pi], vals, bar_width,
                           color=pair_colors[pi], alpha=0.8,
                           label=KEY_PAIR_LABELS[pair] if bi == 0 else None)

            if si == 0:
                ax.set_title(BAND_TEX[bi], fontsize=12, fontweight="bold")
            if bi == 0:
                ax.set_ylabel(f"{stripe_name}\n[{h_lo:.2f}, {h_hi:.2f}]",
                              fontsize=11, fontweight="bold")
            ax.set_ylim(0, 1.05)
            ax.grid(alpha=0.2, axis="y")

            # X-tick labels
            pats_4phase = [p for p in PATIENTS
                           if "task_test" in PATIENT_PHASES[p]]
            ax.set_xticks(range(len(pats_4phase)))
            ax.set_xticklabels([p.replace("Pat_0", "P") for p in pats_4phase],
                               fontsize=8, rotation=45)

    # Legend from first row
    axes[0, 0].legend(fontsize=7, loc="upper left")

    fig.suptitle(
        "MS-ARI Distance by Scale Stripe — micro / meso / macro\n"
        "(ARI clipped to [0,1], each stripe integrated separately)",
        fontsize=15, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    path = OUTPUT_DIR / "scale_stripes_comparison.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


# ===================================================================
# Figure 4: Per-patient heatmaps (phase × phase, h >= 0.5)
# ===================================================================
def plot_per_patient_heatmaps(dists_h50):
    for pat in PATIENTS:
        phases = PATIENT_PHASES[pat]
        n = len(phases)

        fig, axes_grid = plt.subplots(2, 3, figsize=(20, 12))
        axes_flat = axes_grid.ravel()

        for j, band in enumerate(BANDS):
            ax = axes_flat[j]
            matrix = np.full((n, n), np.nan)
            for i1, p1 in enumerate(phases):
                for i2, p2 in enumerate(phases):
                    if i1 < i2:
                        pair = (p1, p2)
                        val = dists_h50.get((pat, band, pair), np.nan)
                        matrix[i1, i2] = val
                        matrix[i2, i1] = val

            mask = np.eye(n, dtype=bool)
            display = np.where(mask, np.nan, matrix)
            im = ax.imshow(display, cmap="YlOrRd", vmin=0, vmax=1, aspect="equal")
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.06, shrink=0.85)

            for i1 in range(n):
                for i2 in range(n):
                    if i1 != i2:
                        val = matrix[i1, i2]
                        if not np.isnan(val):
                            color = "white" if val > 0.55 else "black"
                            ax.text(i2, i1, f"{val:.2f}", ha="center", va="center",
                                    fontsize=11, fontweight="bold", color=color)

            ax.set_xticks(range(n))
            ax.set_yticks(range(n))
            ax.set_xticklabels(phases, rotation=45, ha="right", fontsize=9)
            ax.set_yticklabels(phases, fontsize=9)
            ax.set_title(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=13, fontweight="bold")

        fig.suptitle(
            f"{pat} — MS-ARI phase distances (h ≥ 0.5, clipped)\n"
            f"({n} phases: {', '.join(phases)})",
            fontsize=14, fontweight="bold", y=1.01,
        )
        plt.tight_layout()
        path = OUTPUT_DIR / f"heatmaps_msari_v2_{pat}.pdf"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {path}")


# ===================================================================
# Results table
# ===================================================================
def build_results_table(curves):
    """Build comprehensive CSV with distances at all scale ranges."""
    rows = []
    for pat in PATIENTS:
        phases = PATIENT_PHASES[pat]
        pairs = list(combinations(phases, 2))
        for band in BANDS:
            for pair in pairs:
                key = (pat, band, pair)
                if key not in curves:
                    continue
                t, a = curves[key]
                row = {
                    "patient": pat,
                    "band": band,
                    "phase_a": pair[0],
                    "phase_b": pair[1],
                    "D_full": msari_distance(t, a, 0.01, 1.0),
                    "D_h50": msari_distance(t, a, 0.5, 1.0),
                }
                for name, (lo, hi) in SCALE_STRIPES.items():
                    row[f"D_{name}"] = msari_distance(t, a, lo, hi)
                rows.append(row)
    return pd.DataFrame(rows)


# ===================================================================
# Main
# ===================================================================
def main():
    print("Loading all LRG results...")
    lrg_data = load_all_lrg()
    print(f"  Loaded {len(lrg_data)} entries")

    print("\nComputing ARI curves (this takes a moment)...")
    curves = compute_all_ari_curves(lrg_data)
    print(f"  Computed {len(curves)} ARI curves")

    print("\nComputing distances (h >= 0.5)...")
    dists_h50 = compute_distances_from_curves(curves, h_lo=0.5, h_hi=1.0)

    # Results CSV
    print("\nBuilding results table...")
    df = build_results_table(curves)
    csv_path = OUTPUT_DIR / "msari_v2_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")

    # Figures
    print("\n--- Diagnostic panel: Pat_08 theta ---")
    plot_diagnostic_panel(lrg_data, curves, "Pat_08", "theta")

    print("\n--- Pairwise distances (h >= 0.5) ---")
    plot_pairwise_distances(dists_h50)

    print("\n--- Scale stripes comparison ---")
    plot_scale_stripes(curves)

    print("\n--- Per-patient heatmaps ---")
    plot_per_patient_heatmaps(dists_h50)

    # Print summary
    print("\n" + "=" * 80)
    print("KEY PAIRWISE DISTANCES (h >= 0.5) — SUMMARY")
    print("=" * 80)
    for pair in KEY_PAIRS:
        label = KEY_PAIR_LABELS[pair]
        print(f"\n  {label}:")
        pats_with = [p for p in PATIENTS
                     if pair[0] in PATIENT_PHASES[p] and pair[1] in PATIENT_PHASES[p]]
        print(f"  {'Band':<12}", end="")
        for pat in pats_with:
            print(f"  {pat:>8}", end="")
        print(f"  {'Mean':>8}  {'Std':>8}")
        print(f"  {'-'*12}", end="")
        for _ in pats_with:
            print(f"  {'--------':>8}", end="")
        print(f"  {'--------':>8}  {'--------':>8}")
        for band in BANDS:
            vals = []
            print(f"  {band:<12}", end="")
            for pat in pats_with:
                d = dists_h50.get((pat, band, pair), np.nan)
                vals.append(d)
                print(f"  {d:>8.3f}", end="")
            clean = [v for v in vals if not np.isnan(v)]
            m = np.mean(clean) if clean else np.nan
            s = np.std(clean) if clean else np.nan
            print(f"  {m:>8.3f}  {s:>8.3f}")

    # Scale stripe summary
    print("\n" + "=" * 80)
    print("SCALE STRIPE COMPARISON — rest_pre↔rest_post (h-range means across patients)")
    print("=" * 80)
    pats_all = [p for p in PATIENTS if "rest_post" in PATIENT_PHASES[p]]
    pair = ("rest_pre", "rest_post")
    print(f"  {'Band':<12}", end="")
    for name in SCALE_STRIPES:
        print(f"  {name:>10}", end="")
    print()
    for band in BANDS:
        print(f"  {band:<12}", end="")
        for name, (lo, hi) in SCALE_STRIPES.items():
            vals = []
            for pat in pats_all:
                key = (pat, band, pair)
                if key in curves:
                    t, a = curves[key]
                    vals.append(msari_distance(t, a, lo, hi))
            m = np.mean(vals) if vals else np.nan
            print(f"  {m:>10.3f}", end="")
        print()


if __name__ == "__main__":
    main()
