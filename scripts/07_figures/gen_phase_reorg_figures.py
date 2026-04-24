#!/usr/bin/env python3
"""Section 5.3 — Phase Reorganization Structure (Pat_02).

Figures:
  1. fig_reorg_matrices_Pat_02.pdf        — 3x6 grid of 4x4 reorganization matrices
  2. fig_ultrametric_sidebyside_Pat_02.pdf — raw ultrametric matrices (theta vs high-gamma)
  3. fig_dendrograms_phase_comparison_Pat_02.pdf — phase-resolved dendrograms (3 bands)
  4. fig_node_swaps_Pat_02_theta_*.pdf     — spatial node-swap maps (theta band)

Run: python scripts/gen_phase_reorg_figures.py
"""
from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec

# ── Global font scaling (affects ALL figures) ─────────────────────────
plt.rcParams.update({
    "font.size":        14,
    "axes.titlesize":   15,
    "axes.labelsize":   14,
    "xtick.labelsize":  12,
    "ytick.labelsize":  12,
    "legend.fontsize":  12,
})
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, dendrogram, leaves_list
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import squareform

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PHASE_LABELS,
)
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.visuals.lrg import compute_partition_stability_index, _find_psi_optimal_partition

# ── Config ────────────────────────────────────────────────────────────
PATIENT = "Pat_02"
from lrg_eegfc.config.paths import LRG_CACHE, FIGURES_ROOT, METRIC_CONCORDANCE_CACHE, SEEG_DATAPATH
FC_METHOD = "msc"

OUTPUT_DIR = FIGURES_ROOT / "phase_reorganization" / PATIENT
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DATA_DIR = METRIC_CONCORDANCE_CACHE

PHASES = list(PHASE_LABELS)
BANDS = list(BRAIN_BANDS_NAMES)
BAND_TEX = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]
PHASE_SHORT = {"rest_pre": "rest_pre", "task_learn": "tLearn", "task_test": "tTest", "rest_post": "rest_post"}

SELECTED_METRICS = ["scaled_distance", "tree_robinson_foulds", "tree_fowlkes_mallows"]
METRIC_LABELS = {
    "scaled_distance": "Scaled distance",
    "tree_robinson_foulds": "Robinson\u2013Foulds",
    "tree_fowlkes_mallows": "Fowlkes\u2013Mallows",
}

# Bands for dendrograms (3 bands: high-reorg, mid, low-reorg)
DENDRO_BANDS = ["theta", "beta", "high_gamma"]
DENDRO_BAND_LABELS = [BRAIN_BAND_TEX_DICT[b] + " band" for b in DENDRO_BANDS]

# Community colors
COMM_COLORS = plt.cm.tab10.colors


# ── Helpers ───────────────────────────────────────────────────────────
def align_labels_hungarian(ref_labels, target_labels):
    """Align target labels to reference via maximum overlap (Hungarian)."""
    ref = np.asarray(ref_labels)
    tgt = np.asarray(target_labels)
    u_ref, u_tgt = np.unique(ref), np.unique(tgt)
    dim = max(len(u_ref), len(u_tgt))
    cost = np.zeros((dim, dim))
    for i, r in enumerate(u_ref):
        for j, t in enumerate(u_tgt):
            cost[i, j] = -np.sum((ref == r) & (tgt == t))
    ri, ci = linear_sum_assignment(cost)
    mapping = {}
    for r, c in zip(ri, ci):
        if r < len(u_ref) and c < len(u_tgt):
            mapping[u_tgt[c]] = u_ref[r]
    used = set(mapping.values())
    nxt = max(max(used) + 1, max(u_ref) + 1) if used else 1
    for t in u_tgt:
        if t not in mapping:
            mapping[t] = nxt
            nxt += 1
    return np.array([mapping[l] for l in tgt])


def load_all_lrg():
    """Load all LRG results for PATIENT."""
    results = {}
    for phase in PHASES:
        for band in BANDS:
            res = load_lrg_result(PATIENT, phase, band, FC_METHOD, cache_root=LRG_CACHE)
            if res is not None:
                results[(phase, band)] = res
    return results


def get_optimal_partition(linkage_matrix):
    """Get optimal n* from PSI."""
    psi, n_comm = compute_partition_stability_index(linkage_matrix)
    if len(psi) == 0:
        return 2
    return _find_psi_optimal_partition(psi, n_comm)


def build_node_to_comm(Z, labels):
    """Build node->community mapping for the full dendrogram tree."""
    n_nodes = len(labels)
    node_to_comm = {}
    for i in range(n_nodes):
        node_to_comm[i] = int(labels[i])
    for i in range(len(Z)):
        left_id = int(Z[i, 0])
        right_id = int(Z[i, 1])
        lc = node_to_comm.get(left_id, -1)
        rc = node_to_comm.get(right_id, -1)
        node_to_comm[n_nodes + i] = lc if lc == rc else -1
    return node_to_comm


def make_link_color_func(node_to_comm):
    """Create link color function from node-to-community mapping."""
    def fn(k):
        c = node_to_comm.get(k, -1)
        if c >= 0:
            return mcolors.to_hex(COMM_COLORS[c % len(COMM_COLORS)])
        return "#aaaaaa"
    return fn


# ── Load data ─────────────────────────────────────────────────────────
print(f"Loading LRG results for {PATIENT}...")
lrg = load_all_lrg()
print(f"  {len(lrg)} results loaded")

df_metrics = pd.read_csv(DATA_DIR / f"metric_concordance_{PATIENT}.csv")


# ======================================================================
# Figure 1: 3x6 reorganization matrices
# ======================================================================
print("\nFigure 1: Reorganization matrices...")

fig = plt.figure(figsize=(15, 7.5))
# Manual gridspec: 3 rows, 6 heatmap cols + 1 colorbar col per row
# Use nested gridspec for clean colorbar placement
gs_outer = GridSpec(3, 1, figure=fig, hspace=0.35)

for row_idx, metric in enumerate(SELECTED_METRICS):
    sub = df_metrics[df_metrics["metric"] == metric]
    vmin_row = sub["value"].min()
    vmax_row = sub["value"].max()

    gs_row = gs_outer[row_idx].subgridspec(1, 7, width_ratios=[1, 1, 1, 1, 1, 1, 0.08], wspace=0.15)

    last_im = None
    for col_idx, band in enumerate(BANDS):
        ax = fig.add_subplot(gs_row[0, col_idx])

        mat = np.zeros((4, 4))
        for _, r in sub[sub["band"] == band].iterrows():
            i = PHASES.index(r["phase_A"])
            j = PHASES.index(r["phase_B"])
            mat[i, j] = r["value"]
            mat[j, i] = r["value"]

        im = ax.imshow(mat, cmap="YlOrRd", vmin=vmin_row, vmax=vmax_row,
                        aspect="equal", interpolation="nearest")
        last_im = im

        for i in range(4):
            for j in range(4):
                val = mat[i, j]
                norm_val = (val - vmin_row) / max(vmax_row - vmin_row, 1e-10)
                color = "white" if norm_val > 0.6 else "black"
                fmt = f"{val:.3f}" if val < 1 else f"{val:.1f}"
                ax.text(j, i, fmt, ha="center", va="center",
                       fontsize=10, color=color)

        if row_idx == 2:
            ax.set_xticks(range(4))
            ax.set_xticklabels([PHASE_SHORT[p] for p in PHASES],
                               rotation=35, ha="right")
        else:
            ax.set_xticks([])
        if col_idx == 0:
            ax.set_yticks(range(4))
            ax.set_yticklabels([PHASE_SHORT[p] for p in PHASES])
            ax.set_ylabel(METRIC_LABELS[metric], fontweight="bold")
        else:
            ax.set_yticks([])

        if row_idx == 0:
            ax.set_title(BAND_TEX[col_idx])

    # Colorbar in last column
    cax = fig.add_subplot(gs_row[0, 6])
    cbar = fig.colorbar(last_im, cax=cax)

fig1_path = OUTPUT_DIR / f"fig_reorg_matrices_{PATIENT}.pdf"
fig.savefig(fig1_path, bbox_inches="tight", dpi=300)
plt.close(fig)
print(f"  Saved: {fig1_path}")


# ======================================================================
# Figure 2: Ultrametric matrices side-by-side (theta vs high-gamma)
# ======================================================================
print("\nFigure 2: Ultrametric matrices side-by-side...")

res_ref = lrg[("rest_pre", "theta")]
leaf_order = leaves_list(res_ref.linkage_matrix)

contrast_bands = ["theta", "high_gamma"]
band_labels_2 = [BRAIN_BAND_TEX_DICT[b] + " band" for b in contrast_bands]

fig = plt.figure(figsize=(13, 6.5))
gs = GridSpec(2, 5, figure=fig, width_ratios=[1, 1, 1, 1, 0.05], wspace=0.08, hspace=0.15)

for row_idx, band in enumerate(contrast_bands):
    matrices = []
    for phase in PHASES:
        res = lrg.get((phase, band))
        if res is None:
            matrices.append(None)
            continue
        U = squareform(res.ultrametric_matrix)
        U_ordered = U[np.ix_(leaf_order, leaf_order)]
        matrices.append(U_ordered)

    all_vals = np.concatenate([m[m > 0].ravel() for m in matrices if m is not None])
    vmin_log = np.log10(all_vals.min())
    vmax_log = np.log10(all_vals.max())

    last_im = None
    for col_idx, (phase, mat) in enumerate(zip(PHASES, matrices)):
        ax = fig.add_subplot(gs[row_idx, col_idx])
        if mat is None:
            ax.set_visible(False)
            continue

        mat_log = np.log10(np.where(mat > 0, mat, np.nan))
        im = ax.imshow(mat_log, cmap="magma", vmin=vmin_log, vmax=vmax_log,
                        aspect="equal", interpolation="nearest")
        last_im = im
        ax.set_xticks([])
        ax.set_yticks([])

        if row_idx == 0:
            ax.set_title(phase)
        if col_idx == 0:
            ax.set_ylabel(band_labels_2[row_idx], fontweight="bold")

    # Colorbar
    cax = fig.add_subplot(gs[row_idx, 4])
    cbar = fig.colorbar(last_im, cax=cax)
    cbar.set_label(r"$\log_{10}\, D$")

fig2_path = OUTPUT_DIR / f"fig_ultrametric_sidebyside_{PATIENT}.pdf"
fig.savefig(fig2_path, bbox_inches="tight", dpi=300)
plt.close(fig)
print(f"  Saved: {fig2_path}")


# ======================================================================
# Figure 3: Phase-resolved dendrograms (3 bands, vertical)
# ======================================================================
print("\nFigure 3: Phase-resolved dendrograms...")

fig, axes = plt.subplots(len(DENDRO_BANDS), 4, figsize=(16, 4 * len(DENDRO_BANDS)))

for row_idx, band in enumerate(DENDRO_BANDS):
    partitions = {}
    n_stars = {}
    for phase in PHASES:
        res = lrg.get((phase, band))
        if res is None:
            continue
        n_star = get_optimal_partition(res.linkage_matrix)
        labels = fcluster(res.linkage_matrix, n_star, criterion="maxclust")
        partitions[phase] = labels
        n_stars[phase] = n_star

    if "rest_pre" in partitions:
        ref = partitions["rest_pre"]
        for phase in PHASES[1:]:
            if phase in partitions:
                partitions[phase] = align_labels_hungarian(ref, partitions[phase])

    for col_idx, phase in enumerate(PHASES):
        ax = axes[row_idx, col_idx]
        res = lrg.get((phase, band))
        if res is None or phase not in partitions:
            ax.set_visible(False)
            continue

        labels = partitions[phase]
        n_star = n_stars[phase]
        Z = res.linkage_matrix
        n_nodes = len(labels)
        merge_heights = Z[:, 2]

        # tmin / tmax from merge heights
        tmin = merge_heights[0] * 0.8
        tmax = merge_heights[-1] * 1.05

        # Cut height for n_star clusters (midpoint between adjacent merges)
        cut_idx = n_nodes - n_star
        if 0 < cut_idx < len(merge_heights):
            cut_height = (merge_heights[cut_idx - 1] + merge_heights[cut_idx]) / 2
        else:
            cut_height = merge_heights[0] * 0.5

        # Community-colored links
        node_to_comm = build_node_to_comm(Z, labels)
        link_color_fn = make_link_color_func(node_to_comm)

        # Vertical dendrogram (orientation="right"), truncated
        dn = dendrogram(
            Z, ax=ax, orientation="right",
            truncate_mode="lastp", p=min(40, n_nodes - 1),
            link_color_func=link_color_fn,
            above_threshold_color="#aaaaaa",
            no_labels=True,
            color_threshold=0,
        )

        # Log-scale x-axis with tmin/tmax limits
        ax.set_xscale("log")
        ax.set_xlim(tmin, tmax)

        # Cut height line (vertical for horizontal dendrograms)
        ax.axvline(cut_height, color="k", linestyle="--", linewidth=1.0, alpha=0.7)

        if row_idx == 0:
            ax.set_title(f"{phase}\n$n^*={n_star}$")
        else:
            ax.set_title(f"$n^*={n_star}$")
        if col_idx == 0:
            ax.set_ylabel(DENDRO_BAND_LABELS[row_idx], fontweight="bold")
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

fig.tight_layout()
fig3_path = OUTPUT_DIR / f"fig_dendrograms_phase_comparison_{PATIENT}.pdf"
fig.savefig(fig3_path, bbox_inches="tight", dpi=300)
plt.close(fig)
print(f"  Saved: {fig3_path}")


# ======================================================================
# Figure 4: Node swaps in theta band (brain space)
# ======================================================================
print("\nFigure 4: Node swaps (theta band)...")

try:
    from nilearn import plotting as ni_plot
    from lrg_eegfc.visuals.spatial import load_spatial_metadata, prepare_spatial_coordinates

    DATASET_ROOT = SEEG_DATAPATH
    metadata = load_spatial_metadata(PATIENT, DATASET_ROOT)
    coords_mni = prepare_spatial_coordinates(metadata, scale="mm", to_mni=True)
    N = coords_mni.shape[0]

    # Get theta partitions — use consensus n* (median across phases)
    # so Jaccard measures spatial rearrangement, not granularity change
    theta_nstars_raw = {}
    for phase in PHASES:
        res = lrg.get((phase, "theta"))
        if res is None:
            continue
        theta_nstars_raw[phase] = get_optimal_partition(res.linkage_matrix)

    n_star_consensus = int(np.median(list(theta_nstars_raw.values())))
    print(f"  Per-phase n*: {theta_nstars_raw}")
    print(f"  Consensus n* (median): {n_star_consensus}")

    theta_partitions = {}
    for phase in PHASES:
        res = lrg.get((phase, "theta"))
        if res is None:
            continue
        labels = fcluster(res.linkage_matrix, n_star_consensus, criterion="maxclust")
        theta_partitions[phase] = labels

    if "rest_pre" in theta_partitions:
        ref = theta_partitions["rest_pre"]
        for phase in PHASES[1:]:
            if phase in theta_partitions:
                theta_partitions[phase] = align_labels_hungarian(ref, theta_partitions[phase])

    n_lrg = len(theta_partitions.get("rest_pre", []))
    if n_lrg < N:
        print(f"  Note: LRG giant component {n_lrg} nodes, full {N}")
        coords_plot = coords_mni[:n_lrg]
    else:
        coords_plot = coords_mni

    transitions = [
        ("rest_pre", "task_learn", "rest_pre -> tLearn"),
        ("task_learn", "task_test", "tLearn -> tTest"),
        ("task_test", "rest_post", "tTest -> rest_post"),
    ]

    def node_jaccard_reorg(lab_a, lab_b):
        """Per-node reorganization score via Jaccard distance.

        For each node i, compute J = 1 - |C_A(i) ∩ C_B(i)| / |C_A(i) ∪ C_B(i)|
        where C_A(i) is the set of nodes sharing node i's community in partition A.

        Returns array of shape (n_nodes,) in [0, 1]:
          0 = node's community is perfectly preserved
          1 = node's community is completely reorganized
        """
        n = len(lab_a)
        scores = np.zeros(n)
        for i in range(n):
            ca = set(np.where(lab_a == lab_a[i])[0])
            cb = set(np.where(lab_b == lab_b[i])[0])
            scores[i] = 1.0 - len(ca & cb) / len(ca | cb)
        return scores

    for p1, p2, label in transitions:
        if p1 not in theta_partitions or p2 not in theta_partitions:
            continue

        lab1 = theta_partitions[p1]
        lab2 = theta_partitions[p2]
        reorg = node_jaccard_reorg(lab1, lab2)
        n_swap = int(np.sum(lab1 != lab2))
        frac = n_swap / len(lab1) * 100
        mean_j = np.mean(reorg)

        # Node size proportional to reorganization (bigger = more reorganized)
        node_sizes = 20 + 80 * reorg

        display = ni_plot.plot_markers(
            node_values=reorg,
            node_coords=coords_plot,
            node_cmap="YlOrRd",
            node_vmin=0.0,
            node_vmax=1.0,
            node_size=node_sizes,
            display_mode="lzry",
            title=(f"{label}  |  {n_swap}/{len(lab1)} swapped ({frac:.0f}%)"
                   f"  |  mean Jaccard = {mean_j:.2f}  (n*={n_star_consensus})"
                   f"\nColor = per-node Jaccard distance (0 = stable, 1 = fully reorganized)"),
            colorbar=True,
        )

        out = OUTPUT_DIR / f"fig_node_swaps_{PATIENT}_theta_{p1}_{p2}.pdf"
        display.savefig(str(out), dpi=200)
        display.close()
        print(f"  Saved: {out.name}  (mean Jaccard={mean_j:.3f})")

    plt.close("all")

except ImportError:
    print("  SKIP: nilearn not available")
except Exception as e:
    print(f"  WARN: Figure 4 failed: {e}")
    import traceback
    traceback.print_exc()


# ======================================================================
# Summary CSV
# ======================================================================
print("\nSaving summary CSV...")

rows = []
for metric in SELECTED_METRICS:
    sub = df_metrics[df_metrics["metric"] == metric]
    for band in BANDS:
        band_sub = sub[sub["band"] == band].copy()
        band_sub = band_sub.sort_values("value")
        for rank, (_, r) in enumerate(band_sub.iterrows(), 1):
            rows.append({
                "band": band,
                "metric": metric,
                "phase_A": r["phase_A"],
                "phase_B": r["phase_B"],
                "value": r["value"],
                "rank_within_band": rank,
            })

csv_path = OUTPUT_DIR / f"phase_reorg_summary_{PATIENT}.csv"
pd.DataFrame(rows).to_csv(csv_path, index=False)
print(f"  Saved: {csv_path}")

print(f"\nDone! All figures in {OUTPUT_DIR}")
