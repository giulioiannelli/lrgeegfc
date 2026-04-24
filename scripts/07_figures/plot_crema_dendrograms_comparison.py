#!/usr/bin/env python3
"""Side-by-side dendrogram comparison across phases for CReMa-validated LRG.

For each (patient, band): one wide figure with phases side by side.
Each dendrogram has colored leaf labels, cluster-colored links, log-scale,
and PSI-optimal cut line.

Pass --crema (default) or --dense to choose LRG cache.

Output: data/figures/crema_dendrogram_comparison/{Pat_XX}/
"""
from __future__ import annotations

import sys
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, leaves_list
from scipy.optimize import linear_sum_assignment
from scipy.signal import find_peaks

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import LRG_CACHE as _LRG_CACHE, LRG_CREMA_CACHE, SEEG_DATAPATH, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.visuals.lrg import compute_partition_stability_index, _load_channel_labels

# ---------------------------------------------------------------------------
USE_CREMA = "--dense" not in sys.argv

LRG_CACHE = LRG_CREMA_CACHE if USE_CREMA else _LRG_CACHE
TAG = "crema" if USE_CREMA else "dense"
OUTPUT_ROOT = FIGURES_ROOT / f"{TAG}_dendrogram_comparison"
DATASET_ROOT = SEEG_DATAPATH

PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
BANDS = BRAIN_BANDS_NAMES

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
_DISTINCT_COLORS = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#d63b00", "#f032e6", "#469990", "#9A6324", "#800000",
    "#808000", "#000075", "#e03080", "#1a7820", "#7030c0",
    "#b05010", "#206060", "#c04080", "#305090", "#704020",
]


def get_distinct_colors(n):
    if n <= len(_DISTINCT_COLORS):
        return [matplotlib.colors.to_rgba(c) for c in _DISTINCT_COLORS[:n]]
    colors = []
    for i in range(n):
        hue = (i * 0.618033988749895) % 1.0
        colors.append(matplotlib.colors.hsv_to_rgb([hue, 0.75, 0.85]))
    return colors


def get_subtree_cluster(Z, clusters, n_nodes):
    node_cluster = {}
    for i in range(n_nodes):
        node_cluster[i] = clusters[i]
    for i in range(len(Z)):
        left, right = int(Z[i, 0]), int(Z[i, 1])
        cl = node_cluster.get(left, -1)
        cr = node_cluster.get(right, -1)
        node_cluster[n_nodes + i] = cl if (cl == cr and cl != -1) else -1
    return node_cluster


JACCARD_THRESHOLD = 0.1   # min overlap to reuse a color


def _match_clusters_hungarian(labels_ref, labels_target):
    """Hungarian matching on Jaccard overlap between cluster assignments."""
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
    cost = 1.0 - jaccard
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


def compute_consistent_color_maps(lrg_results, phases, n_cut):
    """Compute cluster_cmap dicts for each phase with consistent colors.

    Uses the first phase as reference.  Subsequent phases are Hungarian-matched.
    """
    palette = get_distinct_colors(max(n_cut, 20))

    # Cluster assignments per phase
    assignments = {}
    for phase in phases:
        lrg = lrg_results[phase]
        N = lrg.n_nodes
        nc = max(2, min(n_cut, N - 1))
        assignments[phase] = fcluster(lrg.linkage_matrix, nc, criterion="maxclust")

    # Reference phase: assign colors by descending cluster size
    ref_phase = phases[0]
    ref_labels = assignments[ref_phase]
    ref_clusters = np.unique(ref_labels)
    ref_sizes = sorted(
        [(c, np.sum(ref_labels == c)) for c in ref_clusters],
        key=lambda x: -x[1],
    )
    ref_c2idx = {c: i for i, (c, _) in enumerate(ref_sizes)}

    phase_cmaps = {}
    for phase in phases:
        labels = assignments[phase]
        phase_clusters = np.unique(labels)
        if phase == ref_phase:
            c2idx = ref_c2idx.copy()
        else:
            mapping = _match_clusters_hungarian(ref_labels, labels)
            c2idx = {}
            used = set()
            # greedily assign by decreasing Jaccard
            pairs = sorted(
                [(tc, rc, j) for tc, (rc, j) in mapping.items()
                 if rc is not None and j > JACCARD_THRESHOLD],
                key=lambda x: -x[2],
            )
            for tc, rc, _ in pairs:
                idx = ref_c2idx[rc]
                if idx not in used:
                    c2idx[tc] = idx
                    used.add(idx)
            # fresh colors for unmatched clusters
            avail = sorted(set(range(len(palette))) - used)
            for tc in phase_clusters:
                if tc not in c2idx:
                    c2idx[tc] = avail.pop(0) if avail else len(used) % len(palette)
                    used.add(c2idx[tc])

        phase_cmaps[phase] = {c: palette[idx] for c, idx in c2idx.items()}

    return phase_cmaps


MAX_PANELS = 4
MIN_N_SPACING = 4
LOG_DROP_MAX = 0.5


def find_prominent_cuts(Z):
    """Find up to MAX_PANELS prominent PSI peaks, sorted by n ascending."""
    psi_vals, n_comms = compute_partition_stability_index(Z)
    if len(psi_vals) < 3:
        return [2]

    peaks_idx, _ = find_peaks(psi_vals, prominence=0.001)
    global_max_idx = np.argmax(psi_vals)
    candidates = np.unique(np.concatenate([[global_max_idx], peaks_idx]))
    candidates = candidates[np.argsort(psi_vals[candidates])[::-1]]

    if len(candidates) == 0:
        return [int(n_comms[global_max_idx])]

    # Select prominent peaks: spaced in n, within log-drop threshold
    prominent = [candidates[0]]
    kept_n = {int(n_comms[candidates[0]])}
    for ci in candidates[1:]:
        if len(prominent) >= MAX_PANELS:
            break
        n_val = int(n_comms[ci])
        if any(abs(n_val - kn) < MIN_N_SPACING for kn in kept_n):
            continue
        log_drop = np.log10(psi_vals[prominent[-1]]) - np.log10(max(psi_vals[ci], 1e-12))
        if log_drop <= LOG_DROP_MAX:
            prominent.append(ci)
            kept_n.add(n_val)

    result = sorted(set(int(n_comms[i]) for i in prominent))
    return result


# ---------------------------------------------------------------------------
# Draw one dendrogram panel
# ---------------------------------------------------------------------------
def draw_dendro_panel(ax, Z, N, n_cut, labels_list, show_ylabel=False,
                      cluster_cmap=None):
    """Draw a single dendrogram with cluster coloring at n_cut."""
    merge_heights = Z[:, 2]
    n_cut = max(2, min(n_cut, N - 1))

    # Cut height
    cut_idx = N - n_cut
    if 0 < cut_idx < len(merge_heights):
        cut_height = (merge_heights[cut_idx - 1] + merge_heights[cut_idx]) / 2
    elif cut_idx == 0:
        cut_height = merge_heights[0] * 0.5
    else:
        cut_height = merge_heights[0] * 0.5

    # Cluster assignments
    clusters = fcluster(Z, n_cut, criterion="maxclust")
    leaf_order = leaves_list(Z)

    # Use provided color map or fall back to left-to-right appearance order
    if cluster_cmap is None:
        palette = get_distinct_colors(n_cut)
        first_seen = {}
        for leaf in leaf_order:
            c = clusters[leaf]
            if c not in first_seen:
                first_seen[c] = len(first_seen)
        cluster_lr = sorted(first_seen.keys(), key=lambda c: first_seen[c])
        cluster_cmap = {c: palette[i] for i, c in enumerate(cluster_lr)}

    # Link color function
    node_cluster = get_subtree_cluster(Z, clusters, N)

    def _lcf(node_id, _nc=node_cluster, _cm=cluster_cmap):
        c = _nc.get(node_id, -1)
        return matplotlib.colors.to_hex(_cm[c]) if c != -1 else "0.5"

    # Draw
    dendrogram(Z, ax=ax, no_labels=True, link_color_func=_lcf)

    # Cut line
    ax.axhline(cut_height, color="0.3", ls="--", lw=1.5)

    # Log scale
    tmin = merge_heights[merge_heights > 0].min() * 0.5
    tmax = merge_heights.max() * 2.0
    ax.set_yscale("log")
    ax.set_ylim(tmin, tmax)

    # Staggered leaf labels
    for i, leaf_idx in enumerate(leaf_order):
        x = 5 + i * 10
        y_frac = -0.02 if i % 2 == 0 else -0.07
        label = labels_list[leaf_idx] if leaf_idx < len(labels_list) else f"Ch{leaf_idx}"
        color = cluster_cmap[clusters[leaf_idx]]
        ax.text(x, y_frac, label, transform=ax.get_xaxis_transform(),
                rotation=90, ha="center", va="top", fontsize=7,
                fontweight="bold", color=color, clip_on=False)

    ax.set_xticks([])
    ax.grid(axis="y", alpha=0.2)
    if show_ylabel:
        ax.set_ylabel(r"$\mathcal{D}\,/\,\mathcal{D}_{\max}$", fontsize=13)


# ---------------------------------------------------------------------------
# Plot figures for a (patient, band): one figure per cut level
# ---------------------------------------------------------------------------
def plot_phase_dendrograms(patient, band, phases, output_dir):
    """Side-by-side dendrograms for each phase at multiple PSI cut levels."""

    # Load all LRG results
    lrg_results = {}
    for phase in phases:
        lrg = load_lrg_result(patient, phase, band, "msc", cache_root=LRG_CACHE)
        if lrg is not None:
            lrg_results[phase] = lrg

    if not lrg_results:
        return []

    available_phases = [p for p in phases if p in lrg_results]
    n_phases = len(available_phases)

    # Load channel labels
    n_nodes = list(lrg_results.values())[0].n_nodes
    labels_list = _load_channel_labels(patient, DATASET_ROOT)
    if not labels_list or len(labels_list) < n_nodes:
        labels_list = [f"Ch{i}" for i in range(n_nodes)]

    # Collect prominent PSI cuts from ALL phases (union), deduplicate
    all_cuts = set()
    for phase in available_phases:
        cuts = find_prominent_cuts(lrg_results[phase].linkage_matrix)
        all_cuts.update(cuts)
    all_cuts = sorted(all_cuts)

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []

    for n_cut in all_cuts:
        # Compute consistent color maps across phases via Hungarian matching
        phase_cmaps = compute_consistent_color_maps(
            lrg_results, available_phases, n_cut
        )

        fig_w = max(n_nodes * 0.14, 14) * n_phases / 2
        fig_h = 8
        fig, axes = plt.subplots(1, n_phases, figsize=(fig_w, fig_h))
        if n_phases == 1:
            axes = [axes]

        for col, phase in enumerate(available_phases):
            ax = axes[col]
            lrg = lrg_results[phase]
            draw_dendro_panel(ax, lrg.linkage_matrix, lrg.n_nodes, n_cut,
                              labels_list, show_ylabel=(col == 0),
                              cluster_cmap=phase_cmaps[phase])
            ax.set_title(f"{phase}", fontsize=13, fontweight="bold")

        fig.suptitle(
            f"{patient} — {band_tex} — n={n_cut} communities ({TAG.upper()} LRG)",
            fontsize=15, fontweight="bold", y=1.02,
        )
        plt.tight_layout()
        out_path = output_dir / f"{band}_n{n_cut}_phase_dendrograms.pdf"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        paths.append(out_path)

    return paths


# ---------------------------------------------------------------------------
def main():
    print(f"Mode: {TAG.upper()} | LRG cache: {LRG_CACHE}")
    total = 0

    for patient, phases in PATIENT_PHASES.items():
        out_dir = OUTPUT_ROOT / patient
        print(f"\n{patient} ({len(phases)} phases)")
        for band in BANDS:
            paths = plot_phase_dendrograms(patient, band, phases, out_dir)
            if paths:
                cuts = [p.stem.split("_n")[1].split("_")[0] for p in paths]
                print(f"  {band}: n={','.join(cuts)}  ({len(paths)} figs)")
                total += len(paths)
            else:
                print(f"  {band}: SKIP (no LRG data)")

    print(f"\nDone: {total} figures in {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
