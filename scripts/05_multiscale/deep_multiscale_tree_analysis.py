#!/usr/bin/env python3
"""Deep multiscale structural analysis of LRG dendrogram trees.

Analyses WHERE in the tree (fine / meso / coarse scale) dendrograms differ
between recording phases, and how structural divergence relates to
Variation of Information (VI) at each scale.

Outputs
-------
data/figures/multiscale_investigation/tree_structure_deep/
    fig1_scale_band_cophenetic_corr.png
    fig2_cumulative_merge_profiles.png
    fig3_joint_vi_cophenetic.png
    fig4_max_divergence_heatmap.png
    fig5_cumulative_merge_divergence.png
    fig6_delta_patient_gallery.png
    fig7_dissociation_scatter.png
    fig8_triple_overlay.png
    scale_band_correlations.csv
    cumulative_merge_divergence.csv
    vi_profiles.csv
    FINDINGS.md
"""

from __future__ import annotations

import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

# Suppress warnings for clean output
warnings.filterwarnings("ignore")

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D

from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, sem

# --- Project imports ---
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT

# ===========================================================================
# Configuration
# ===========================================================================
PATIENTS = PATIENTS_4PHASE
PHASES = list(PHASE_LABELS)  # rest_pre, task_learn, task_test, rest_post
BANDS = list(BRAIN_BANDS.keys())  # delta ... high_gamma
FC_METHOD = "msc"
N_SCALE_BANDS = 5
K_VALUES = [2, 3, 4, 5, 6, 8, 10, 12, 15, 20]
N_MERGE_GRID = 50
OUTDIR = FIGURES_ROOT / "multiscale_investigation" / "tree_structure_deep"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Phase-pair classification
WITHIN_CONDITION = [("rest_pre", "rest_post"), ("task_learn", "task_test")]
CROSS_CONDITION = [
    ("rest_pre", "task_learn"),
    ("rest_pre", "task_test"),
    ("rest_post", "task_learn"),
    ("rest_post", "task_test"),
]
# Consecutive transitions (for Fig 5)
TRANSITIONS = [
    ("rest_pre", "task_learn"),
    ("task_learn", "task_test"),
    ("task_test", "rest_post"),
]

# ===========================================================================
# Utility functions
# ===========================================================================

from lrg_eegfc.utils.metrics import compute_vi  # canonical implementation


def normalized_vi(labels1, labels2):
    """VI normalized by joint entropy (0 = identical, 1 = maximally different)."""
    n = len(labels1)
    if n == 0:
        return 0.0
    vi = compute_vi(labels1, labels2)
    # Joint entropy
    joint = {}
    for l1, l2 in zip(labels1, labels2):
        k = (l1, l2)
        joint[k] = joint.get(k, 0) + 1
    h_joint = 0.0
    for cnt in joint.values():
        p = cnt / n
        if p > 0:
            h_joint -= p * np.log(p)
    if h_joint == 0:
        return 0.0
    return vi / h_joint


def get_cophenetic_full(linkage_matrix, n_nodes):
    """Get full cophenetic (ultrametric) distance matrix from linkage."""
    # The cophenetic distance between two observations is the height at which
    # they first merge in the dendrogram
    coph = np.zeros((n_nodes, n_nodes))
    # Build from linkage: each row in Z = [idx1, idx2, distance, count]
    # Use scipy-compatible approach
    from scipy.cluster.hierarchy import cophenet as _coph_fn
    # cophenet(Z) returns (correlation_with_Y, condensed_cophenetic_distances)
    # when called without Y, it returns just the condensed cophenetic distances
    coph_condensed = _coph_fn(linkage_matrix)
    # Actually: cophenet(Z) with no Y returns a tuple (c, d) where c=correlation
    # Let me check the return
    if isinstance(coph_condensed, tuple):
        coph_condensed = coph_condensed[1]  # second element is the distances
    return squareform(coph_condensed)


def scale_band_cophenetic_correlation(Z1, Z2, n1, n2, n_bands=N_SCALE_BANDS):
    """Compare cophenetic distances in scale bands.

    Returns
    -------
    band_corrs : array of shape (n_bands,)
        Spearman correlation in each scale band
    band_labels : list of str
        Labels for each band
    """
    # For trees with different numbers of nodes, we need to handle this.
    # When trees have different n_nodes, their cophenetic matrices have
    # different sizes. We use the ultrametric (condensed) distances directly.
    # BUT: we can only compare if n_nodes are the same (same giant component).
    # If different, return NaN.
    if n1 != n2:
        return np.full(n_bands, np.nan), [f"Band {i+1}" for i in range(n_bands)]

    n = n1
    # Get cophenetic distances in condensed form
    from scipy.cluster.hierarchy import cophenet as _coph_fn
    c1 = _coph_fn(Z1)
    c2 = _coph_fn(Z2)
    if isinstance(c1, tuple):
        c1 = c1[1]
    if isinstance(c2, tuple):
        c2 = c2[1]
    # c1, c2 are now condensed distance vectors
    c1 = np.asarray(c1, dtype=float)
    c2 = np.asarray(c2, dtype=float)

    if len(c1) != len(c2):
        return np.full(n_bands, np.nan), [f"Band {i+1}" for i in range(n_bands)]

    # Pool distances to define scale bands
    all_dists = np.concatenate([c1, c2])
    quantiles = np.linspace(0, 100, n_bands + 1)
    boundaries = np.percentile(all_dists, quantiles)

    band_corrs = np.full(n_bands, np.nan)
    band_labels = []
    for i in range(n_bands):
        lo, hi = boundaries[i], boundaries[i + 1]
        # Average distance of each pair across the two trees
        avg_dist = (c1 + c2) / 2.0
        if i < n_bands - 1:
            mask = (avg_dist >= lo) & (avg_dist < hi)
        else:
            mask = (avg_dist >= lo) & (avg_dist <= hi)

        band_labels.append(f"Q{int(quantiles[i])}-{int(quantiles[i+1])}")

        if np.sum(mask) < 3:
            band_corrs[i] = np.nan
            continue

        rho, _ = spearmanr(c1[mask], c2[mask])
        band_corrs[i] = rho

    return band_corrs, band_labels


def cumulative_merge_profile(Z, n_grid=N_MERGE_GRID):
    """Compute cumulative fraction of merges below each height.

    Returns
    -------
    heights : array of shape (n_grid,)
    cum_frac : array of shape (n_grid,)
    """
    merge_heights = Z[:, 2]
    max_h = merge_heights.max()
    min_h = merge_heights.min()
    if max_h == min_h:
        return np.array([max_h]), np.array([1.0])
    heights = np.linspace(min_h, max_h, n_grid)
    n_merges = len(merge_heights)
    cum_frac = np.array([np.sum(merge_heights <= h) / n_merges for h in heights])
    return heights, cum_frac


def cumulative_merge_divergence(Z1, Z2, n_grid=N_MERGE_GRID):
    """L1 distance between cumulative merge profiles on a common height grid.

    Also returns the height of maximum divergence.
    """
    h1 = Z1[:, 2]
    h2 = Z2[:, 2]
    min_h = min(h1.min(), h2.min())
    max_h = max(h1.max(), h2.max())
    if max_h == min_h:
        return 0.0, min_h, np.array([min_h]), np.zeros(1), np.zeros(1)
    heights = np.linspace(min_h, max_h, n_grid)
    n1 = len(h1)
    n2 = len(h2)
    cf1 = np.array([np.sum(h1 <= h) / n1 for h in heights])
    cf2 = np.array([np.sum(h2 <= h) / n2 for h in heights])
    diff = np.abs(cf1 - cf2)
    l1 = np.trapz(diff, heights) / (max_h - min_h)  # normalize by height range
    max_div_idx = np.argmax(diff)
    max_div_height = heights[max_div_idx]
    return l1, max_div_height, heights, cf1, cf2


def vi_profile(Z1, Z2, n1, n2, k_values=K_VALUES):
    """Compute VI at multiple scales k (number of clusters).

    Returns array of VI values, one per k.
    """
    if n1 != n2:
        return np.full(len(k_values), np.nan)

    vis = []
    for k in k_values:
        if k > n1:
            vis.append(np.nan)
            continue
        labels1 = fcluster(Z1, k, criterion="maxclust")
        labels2 = fcluster(Z2, k, criterion="maxclust")
        vis.append(compute_vi(labels1, labels2))
    return np.array(vis)


def nvi_profile(Z1, Z2, n1, n2, k_values=K_VALUES):
    """Compute normalized VI at multiple scales k."""
    if n1 != n2:
        return np.full(len(k_values), np.nan)

    nvis = []
    for k in k_values:
        if k > n1:
            nvis.append(np.nan)
            continue
        labels1 = fcluster(Z1, k, criterion="maxclust")
        labels2 = fcluster(Z2, k, criterion="maxclust")
        nvis.append(normalized_vi(labels1, labels2))
    return np.array(nvis)


def sackin_index(Z, n_nodes):
    """Compute the Sackin imbalance index of a dendrogram.

    Sackin = sum over all leaves of (depth of that leaf).
    Normalized by n * log2(n).
    """
    # Build tree depths from linkage
    depths = np.zeros(2 * n_nodes - 1)
    for i, row in enumerate(Z):
        left, right = int(row[0]), int(row[1])
        parent = n_nodes + i
        depths[parent] = 0  # root has depth 0 logically; we track from merges
    # Actually compute leaf depths by traversing from root
    n_internal = len(Z)
    children = {}
    for i, row in enumerate(Z):
        parent_id = n_nodes + i
        children[parent_id] = (int(row[0]), int(row[1]))
    # BFS from root
    root = n_nodes + n_internal - 1
    depth_map = {root: 0}
    queue = [root]
    while queue:
        node = queue.pop(0)
        if node in children:
            for child in children[node]:
                depth_map[child] = depth_map[node] + 1
                queue.append(child)
    # Sum leaf depths
    leaf_depths = [depth_map.get(i, 0) for i in range(n_nodes)]
    sackin = sum(leaf_depths)
    # Normalize
    norm = n_nodes * np.log2(n_nodes) if n_nodes > 1 else 1.0
    return sackin / norm


def classify_pair(p1, p2):
    """Classify a phase pair as 'within' or 'cross' condition."""
    pair = tuple(sorted([p1, p2]))
    for a, b in WITHIN_CONDITION:
        if pair == tuple(sorted([a, b])):
            return "within"
    return "cross"


# ===========================================================================
# Data loading
# ===========================================================================
print("Loading all LRG results...")
results = {}  # results[patient][band][phase] = LRGResult
n_loaded = 0
n_failed = 0
for patient in PATIENTS:
    results[patient] = {}
    for band in BANDS:
        results[patient][band] = {}
        for phase in PHASES:
            r = load_lrg_result(patient, phase, band, FC_METHOD)
            if r is not None:
                results[patient][band][phase] = r
                n_loaded += 1
            else:
                print(f"  WARNING: Missing {patient}/{band}/{phase}")
                n_failed += 1
print(f"Loaded {n_loaded} results, {n_failed} missing.")

# ===========================================================================
# ANALYSIS 1: Scale-band cophenetic correlations
# ===========================================================================
print("\n=== Analysis 1: Scale-band cophenetic correlations ===")
scale_band_rows = []
for patient in PATIENTS:
    for band in BANDS:
        for p1, p2 in combinations(PHASES, 2):
            r1 = results[patient][band].get(p1)
            r2 = results[patient][band].get(p2)
            if r1 is None or r2 is None:
                continue
            band_corrs, band_labels = scale_band_cophenetic_correlation(
                r1.linkage_matrix, r2.linkage_matrix, r1.n_nodes, r2.n_nodes
            )
            pair_type = classify_pair(p1, p2)
            for i, (bc, bl) in enumerate(zip(band_corrs, band_labels)):
                scale_band_rows.append({
                    "patient": patient,
                    "freq_band": band,
                    "phase1": p1,
                    "phase2": p2,
                    "pair_type": pair_type,
                    "scale_band_idx": i,
                    "scale_band_label": bl,
                    "spearman_corr": bc,
                })

df_scale = pd.DataFrame(scale_band_rows)
df_scale.to_csv(OUTDIR / "scale_band_correlations.csv", index=False)
print(f"  Saved {len(df_scale)} rows to scale_band_correlations.csv")

# ===========================================================================
# ANALYSIS 2: Cumulative merge profiles
# ===========================================================================
print("\n=== Analysis 2: Cumulative merge profiles ===")
merge_div_rows = []
# Store cumulative profiles for Fig 2 (Pat_02, alpha)
cum_profiles_example = {}

for patient in PATIENTS:
    for band in BANDS:
        # Store individual profiles for this patient/band
        for phase in PHASES:
            r = results[patient][band].get(phase)
            if r is None:
                continue
            if patient == "Pat_02" and band == "alpha":
                h, cf = cumulative_merge_profile(r.linkage_matrix)
                cum_profiles_example[phase] = (h, cf)

        for p1, p2 in combinations(PHASES, 2):
            r1 = results[patient][band].get(p1)
            r2 = results[patient][band].get(p2)
            if r1 is None or r2 is None:
                continue
            l1_dist, max_div_h, _, _, _ = cumulative_merge_divergence(
                r1.linkage_matrix, r2.linkage_matrix
            )
            pair_type = classify_pair(p1, p2)
            merge_div_rows.append({
                "patient": patient,
                "freq_band": band,
                "phase1": p1,
                "phase2": p2,
                "pair_type": pair_type,
                "l1_divergence": l1_dist,
                "max_divergence_height": max_div_h,
            })

df_merge = pd.DataFrame(merge_div_rows)
df_merge.to_csv(OUTDIR / "cumulative_merge_divergence.csv", index=False)
print(f"  Saved {len(df_merge)} rows to cumulative_merge_divergence.csv")

# ===========================================================================
# ANALYSIS 3: VI profiles
# ===========================================================================
print("\n=== Analysis 3: VI profiles ===")
vi_rows = []
for patient in PATIENTS:
    for band in BANDS:
        for p1, p2 in combinations(PHASES, 2):
            r1 = results[patient][band].get(p1)
            r2 = results[patient][band].get(p2)
            if r1 is None or r2 is None:
                continue
            vis = vi_profile(r1.linkage_matrix, r2.linkage_matrix, r1.n_nodes, r2.n_nodes)
            nvis = nvi_profile(r1.linkage_matrix, r2.linkage_matrix, r1.n_nodes, r2.n_nodes)
            pair_type = classify_pair(p1, p2)
            for ki, k in enumerate(K_VALUES):
                vi_rows.append({
                    "patient": patient,
                    "freq_band": band,
                    "phase1": p1,
                    "phase2": p2,
                    "pair_type": pair_type,
                    "k": k,
                    "vi": vis[ki],
                    "nvi": nvis[ki],
                })

df_vi = pd.DataFrame(vi_rows)
df_vi.to_csv(OUTDIR / "vi_profiles.csv", index=False)
print(f"  Saved {len(df_vi)} rows to vi_profiles.csv")

# ===========================================================================
# ANALYSIS 4: Sackin index per tree
# ===========================================================================
print("\n=== Analysis 4: Sackin imbalance index ===")
sackin_rows = []
for patient in PATIENTS:
    for band in BANDS:
        for phase in PHASES:
            r = results[patient][band].get(phase)
            if r is None:
                continue
            s = sackin_index(r.linkage_matrix, r.n_nodes)
            sackin_rows.append({
                "patient": patient,
                "freq_band": band,
                "phase": phase,
                "sackin": s,
                "n_nodes": r.n_nodes,
            })
df_sackin = pd.DataFrame(sackin_rows)
df_sackin.to_csv(OUTDIR / "sackin_index.csv", index=False)
print(f"  Saved {len(df_sackin)} rows to sackin_index.csv")


# ===========================================================================
# FIGURE 1: Scale-band cophenetic correlation profiles
# ===========================================================================
print("\n=== Figure 1: Scale-band cophenetic correlation profiles ===")
fig1, axes1 = plt.subplots(2, 3, figsize=(16, 10), sharey=True)
axes1 = axes1.flatten()

for bi, band in enumerate(BANDS):
    ax = axes1[bi]
    df_b = df_scale[df_scale["freq_band"] == band]
    if df_b.empty:
        ax.set_title(BRAIN_BAND_TEX_DICT[band])
        continue

    # Within-condition
    df_within = df_b[df_b["pair_type"] == "within"]
    if not df_within.empty:
        means_w = df_within.groupby("scale_band_idx")["spearman_corr"].mean()
        sems_w = df_within.groupby("scale_band_idx")["spearman_corr"].apply(
            lambda x: sem(x.dropna()) if len(x.dropna()) > 1 else 0
        )
        ax.errorbar(
            means_w.index, means_w.values, yerr=sems_w.values,
            color="steelblue", marker="o", linewidth=2, capsize=3,
            label="Within-condition", zorder=3
        )

    # Cross-condition
    df_cross = df_b[df_b["pair_type"] == "cross"]
    if not df_cross.empty:
        means_c = df_cross.groupby("scale_band_idx")["spearman_corr"].mean()
        sems_c = df_cross.groupby("scale_band_idx")["spearman_corr"].apply(
            lambda x: sem(x.dropna()) if len(x.dropna()) > 1 else 0
        )
        ax.errorbar(
            means_c.index, means_c.values, yerr=sems_c.values,
            color="firebrick", marker="s", linewidth=2, capsize=3,
            label="Cross-condition", zorder=3
        )

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14)
    ax.set_xlabel("Scale band (fine -> coarse)", fontsize=10)
    ax.set_xticks(range(N_SCALE_BANDS))
    # Use the labels from the first available row
    first_labels = df_b.drop_duplicates("scale_band_idx").sort_values("scale_band_idx")["scale_band_label"].tolist()
    if len(first_labels) == N_SCALE_BANDS:
        ax.set_xticklabels(first_labels, fontsize=8)
    ax.axhline(0, color="gray", linestyle="--", alpha=0.5)
    ax.set_ylim(-0.3, 1.05)
    if bi == 0:
        ax.set_ylabel("Spearman correlation", fontsize=12)
    if bi == 0:
        ax.legend(fontsize=9, loc="lower left")

fig1.suptitle(
    "Scale-band cophenetic correlation: WHERE in the tree do phases differ?",
    fontsize=14, fontweight="bold", y=0.98
)
fig1.tight_layout(rect=[0, 0, 1, 0.95])
fig1.savefig(OUTDIR / "fig1_scale_band_cophenetic_corr.png", dpi=200, bbox_inches="tight")
plt.close(fig1)
print("  Saved fig1_scale_band_cophenetic_corr.png")


# ===========================================================================
# FIGURE 2: Cumulative merge profiles (Pat_02, alpha)
# ===========================================================================
print("\n=== Figure 2: Cumulative merge profiles ===")
fig2, ax2 = plt.subplots(figsize=(8, 6))
phase_colors = {
    "rest_pre": "#2166ac", "task_learn": "#d6604d",
    "task_test": "#f4a582", "rest_post": "#4393c3"
}
phase_names = {
    "rest_pre": "Resting Pre", "task_learn": "Task Learn",
    "task_test": "Task Test", "rest_post": "Resting Post"
}
for phase in PHASES:
    if phase in cum_profiles_example:
        h, cf = cum_profiles_example[phase]
        ax2.plot(h, cf, color=phase_colors[phase], linewidth=2.5,
                 label=phase_names[phase])

ax2.set_xlabel("Merge height (ultrametric distance)", fontsize=12)
ax2.set_ylabel("Cumulative fraction of merges", fontsize=12)
ax2.set_title("Cumulative merge profiles -- Pat_02, alpha band", fontsize=13, fontweight="bold")
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)
fig2.tight_layout()
fig2.savefig(OUTDIR / "fig2_cumulative_merge_profiles.png", dpi=200, bbox_inches="tight")
plt.close(fig2)
print("  Saved fig2_cumulative_merge_profiles.png")


# ===========================================================================
# FIGURE 3: Joint VI + cophenetic by scale
# ===========================================================================
print("\n=== Figure 3: Joint VI + cophenetic ===")
fig3, axes3 = plt.subplots(2, 3, figsize=(18, 11))
axes3 = axes3.flatten()

# Map k values to approximate scale band indices for overlay
# k=2 is coarsest (scale band 4), k=20 is finest (scale band 0)
# We approximate: reverse the k ordering to align with scale band

for bi, band in enumerate(BANDS):
    ax = axes3[bi]
    ax2_right = ax.twinx()

    df_b_vi = df_vi[df_vi["freq_band"] == band]
    df_b_sc = df_scale[df_scale["freq_band"] == band]

    if df_b_vi.empty or df_b_sc.empty:
        ax.set_title(BRAIN_BAND_TEX_DICT[band])
        continue

    # VI profile averaged across all patients and phase pairs
    vi_mean = df_b_vi.groupby("k")["vi"].mean()
    vi_sem_vals = df_b_vi.groupby("k")["vi"].apply(
        lambda x: sem(x.dropna()) if len(x.dropna()) > 1 else 0
    )

    # Scale-band corr averaged across all patients and phase pairs
    sc_mean = df_b_sc.groupby("scale_band_idx")["spearman_corr"].mean()
    sc_sem_vals = df_b_sc.groupby("scale_band_idx")["spearman_corr"].apply(
        lambda x: sem(x.dropna()) if len(x.dropna()) > 1 else 0
    )

    # Plot VI (left y-axis)
    l1 = ax.errorbar(
        vi_mean.index, vi_mean.values, yerr=vi_sem_vals.values,
        color="darkviolet", marker="o", linewidth=2, capsize=3,
        label="VI(k)"
    )
    ax.set_xlabel("k (number of clusters)", fontsize=10)
    ax.set_ylabel("VI (nats)", fontsize=10, color="darkviolet")
    ax.tick_params(axis="y", labelcolor="darkviolet")

    # Plot cophenetic correlation (right y-axis)
    # Map scale bands: 0=fine ... 4=coarse
    # X-axis for scale bands uses the same axis but secondary
    l2 = ax2_right.errorbar(
        sc_mean.index, sc_mean.values, yerr=sc_sem_vals.values,
        color="darkorange", marker="s", linewidth=2, capsize=3,
        linestyle="--", label="Cophenetic corr"
    )
    ax2_right.set_ylabel("Spearman corr", fontsize=10, color="darkorange")
    ax2_right.tick_params(axis="y", labelcolor="darkorange")

    # Title and legend
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14)

    # Combined legend
    lines = [l1, l2]
    labels_legend = ["VI(k) - partition change", "Cophenetic corr - structural similarity"]
    ax.legend(lines, labels_legend, fontsize=8, loc="upper left")

fig3.suptitle(
    "Joint analysis: partition divergence (VI) vs structural similarity (cophenetic)\n"
    "Left: VI grows with k; Right: cophenetic corr per scale band",
    fontsize=13, fontweight="bold", y=1.00
)
fig3.tight_layout(rect=[0, 0, 1, 0.95])
fig3.savefig(OUTDIR / "fig3_joint_vi_cophenetic.png", dpi=200, bbox_inches="tight")
plt.close(fig3)
print("  Saved fig3_joint_vi_cophenetic.png")


# ===========================================================================
# FIGURE 4: Scale of maximum divergence heatmaps
# ===========================================================================
print("\n=== Figure 4: Maximum divergence heatmaps ===")
fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(16, 6))

# 4a: Scale band of LOWEST cophenetic correlation (most structural divergence)
max_div_struct = np.full((len(PATIENTS), len(BANDS)), np.nan)
for pi, patient in enumerate(PATIENTS):
    for bi, band in enumerate(BANDS):
        df_pb = df_scale[(df_scale["patient"] == patient) & (df_scale["freq_band"] == band)]
        if df_pb.empty:
            continue
        avg_by_band = df_pb.groupby("scale_band_idx")["spearman_corr"].mean()
        if not avg_by_band.empty:
            max_div_struct[pi, bi] = avg_by_band.idxmin()

im4a = ax4a.imshow(max_div_struct, cmap="viridis", aspect="auto",
                    vmin=0, vmax=N_SCALE_BANDS - 1)
ax4a.set_xticks(range(len(BANDS)))
ax4a.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=11)
ax4a.set_yticks(range(len(PATIENTS)))
ax4a.set_yticklabels(PATIENTS, fontsize=11)
ax4a.set_title("Scale band of lowest cophenetic corr\n(0=fine, 4=coarse)", fontsize=12, fontweight="bold")
plt.colorbar(im4a, ax=ax4a, label="Scale band index")
# Add text annotations
for pi in range(len(PATIENTS)):
    for bi in range(len(BANDS)):
        val = max_div_struct[pi, bi]
        if not np.isnan(val):
            ax4a.text(bi, pi, f"{int(val)}", ha="center", va="center",
                     fontsize=10, fontweight="bold",
                     color="white" if val > 2 else "black")

# 4b: k value of HIGHEST VI (scale of max partition divergence)
max_div_vi = np.full((len(PATIENTS), len(BANDS)), np.nan)
for pi, patient in enumerate(PATIENTS):
    for bi, band in enumerate(BANDS):
        df_pb = df_vi[(df_vi["patient"] == patient) & (df_vi["freq_band"] == band)]
        if df_pb.empty:
            continue
        avg_by_k = df_pb.groupby("k")["vi"].mean()
        if not avg_by_k.empty:
            max_div_vi[pi, bi] = avg_by_k.idxmax()

im4b = ax4b.imshow(max_div_vi, cmap="plasma", aspect="auto")
ax4b.set_xticks(range(len(BANDS)))
ax4b.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=11)
ax4b.set_yticks(range(len(PATIENTS)))
ax4b.set_yticklabels(PATIENTS, fontsize=11)
ax4b.set_title("k of highest VI\n(larger k = finer scale)", fontsize=12, fontweight="bold")
plt.colorbar(im4b, ax=ax4b, label="k (number of clusters)")
for pi in range(len(PATIENTS)):
    for bi in range(len(BANDS)):
        val = max_div_vi[pi, bi]
        if not np.isnan(val):
            ax4b.text(bi, pi, f"{int(val)}", ha="center", va="center",
                     fontsize=10, fontweight="bold",
                     color="white" if val > 10 else "black")

fig4.suptitle(
    "WHERE is maximum divergence? Structural (left) vs partition (right)",
    fontsize=14, fontweight="bold", y=1.02
)
fig4.tight_layout()
fig4.savefig(OUTDIR / "fig4_max_divergence_heatmap.png", dpi=200, bbox_inches="tight")
plt.close(fig4)
print("  Saved fig4_max_divergence_heatmap.png")


# ===========================================================================
# FIGURE 5: Cumulative merge divergence per transition
# ===========================================================================
print("\n=== Figure 5: Cumulative merge divergence per transition ===")
fig5, axes5 = plt.subplots(1, 2, figsize=(14, 6))

# 5a: L1 divergence per transition, by band
transition_labels = [f"{p1} -> {p2}" for p1, p2 in TRANSITIONS]
# Also add rest_pre -> rest_post for reference
ALL_TRANSITIONS = TRANSITIONS + [("rest_pre", "rest_post")]
all_trans_labels = transition_labels + ["rest_pre -> rest_post"]

bar_data = {}
for band in BANDS:
    bar_data[band] = []
    for p1, p2 in ALL_TRANSITIONS:
        vals = []
        for patient in PATIENTS:
            row = df_merge[(df_merge["patient"] == patient) &
                          (df_merge["freq_band"] == band) &
                          (df_merge["phase1"] == p1) &
                          (df_merge["phase2"] == p2)]
            if row.empty:
                # Try reverse order
                row = df_merge[(df_merge["patient"] == patient) &
                              (df_merge["freq_band"] == band) &
                              (df_merge["phase1"] == p2) &
                              (df_merge["phase2"] == p1)]
            if not row.empty:
                vals.append(row.iloc[0]["l1_divergence"])
        bar_data[band].append(np.nanmean(vals) if vals else np.nan)

ax5a = axes5[0]
x = np.arange(len(ALL_TRANSITIONS))
width = 0.12
band_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
for bi, band in enumerate(BANDS):
    ax5a.bar(x + bi * width - 0.3, bar_data[band], width,
             label=BRAIN_BAND_TEX_DICT[band], color=band_colors[bi], alpha=0.85)
ax5a.set_xticks(x)
ax5a.set_xticklabels(all_trans_labels, fontsize=9, rotation=20, ha="right")
ax5a.set_ylabel("Normalized L1 divergence", fontsize=11)
ax5a.set_title("Merge profile divergence per transition", fontsize=12, fontweight="bold")
ax5a.legend(fontsize=8, ncol=3, loc="upper left")
ax5a.grid(axis="y", alpha=0.3)

# 5b: Heatmap of L1 divergence, bands x all phase pairs
all_pairs = list(combinations(PHASES, 2))
pair_labels = [f"{p1[:4]}-{p2[:4]}" for p1, p2 in all_pairs]
heat_data = np.full((len(BANDS), len(all_pairs)), np.nan)
for bi, band in enumerate(BANDS):
    for pi, (p1, p2) in enumerate(all_pairs):
        vals = []
        for patient in PATIENTS:
            row = df_merge[(df_merge["patient"] == patient) &
                          (df_merge["freq_band"] == band) &
                          (((df_merge["phase1"] == p1) & (df_merge["phase2"] == p2)) |
                           ((df_merge["phase1"] == p2) & (df_merge["phase2"] == p1)))]
            if not row.empty:
                vals.append(row.iloc[0]["l1_divergence"])
        heat_data[bi, pi] = np.nanmean(vals) if vals else np.nan

ax5b = axes5[1]
im5 = ax5b.imshow(heat_data, cmap="YlOrRd", aspect="auto")
ax5b.set_xticks(range(len(all_pairs)))
ax5b.set_xticklabels(pair_labels, fontsize=9, rotation=45, ha="right")
ax5b.set_yticks(range(len(BANDS)))
ax5b.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=11)
ax5b.set_title("Merge divergence: bands x phase pairs", fontsize=12, fontweight="bold")
plt.colorbar(im5, ax=ax5b, label="L1 divergence")
for bi in range(len(BANDS)):
    for pi in range(len(all_pairs)):
        val = heat_data[bi, pi]
        if not np.isnan(val):
            ax5b.text(pi, bi, f"{val:.2f}", ha="center", va="center", fontsize=8,
                     color="white" if val > np.nanmax(heat_data) * 0.6 else "black")

fig5.suptitle(
    "Which transitions cause the biggest structural reshuffling?",
    fontsize=14, fontweight="bold", y=1.02
)
fig5.tight_layout()
fig5.savefig(OUTDIR / "fig5_cumulative_merge_divergence.png", dpi=200, bbox_inches="tight")
plt.close(fig5)
print("  Saved fig5_cumulative_merge_divergence.png")


# ===========================================================================
# FIGURE 6: Per-patient gallery for delta band
# ===========================================================================
print("\n=== Figure 6: Delta band patient gallery ===")
fig6, axes6 = plt.subplots(len(PATIENTS), 1, figsize=(10, 4 * len(PATIENTS)), sharex=False)
if len(PATIENTS) == 1:
    axes6 = [axes6]

for pi, patient in enumerate(PATIENTS):
    ax = axes6[pi]
    for phase in PHASES:
        r = results[patient]["delta"].get(phase)
        if r is None:
            continue
        h, cf = cumulative_merge_profile(r.linkage_matrix, n_grid=80)
        ax.plot(h, cf, color=phase_colors[phase], linewidth=2,
                label=phase_names[phase])

    ax.set_ylabel("Cumulative fraction", fontsize=10)
    ax.set_title(f"{patient} - delta band", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, alpha=0.3)

    # Annotate with divergence metrics
    div_texts = []
    for p1, p2 in [("rest_pre", "rest_post"), ("rest_pre", "task_learn")]:
        r1 = results[patient]["delta"].get(p1)
        r2 = results[patient]["delta"].get(p2)
        if r1 is not None and r2 is not None:
            l1d, _, _, _, _ = cumulative_merge_divergence(r1.linkage_matrix, r2.linkage_matrix)
            div_texts.append(f"L1({p1[:4]}-{p2[:4]})={l1d:.3f}")
    if div_texts:
        ax.text(0.02, 0.95, "\n".join(div_texts), transform=ax.transAxes,
                fontsize=8, verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

axes6[-1].set_xlabel("Merge height", fontsize=11)
fig6.suptitle(
    "Delta band: cumulative merge profiles per patient\n(where do the merges happen?)",
    fontsize=14, fontweight="bold", y=1.01
)
fig6.tight_layout()
fig6.savefig(OUTDIR / "fig6_delta_patient_gallery.png", dpi=200, bbox_inches="tight")
plt.close(fig6)
print("  Saved fig6_delta_patient_gallery.png")


# ===========================================================================
# FIGURE 7: Dissociation scatter
# ===========================================================================
print("\n=== Figure 7: Dissociation scatter ===")
fig7, axes7 = plt.subplots(2, 3, figsize=(16, 10))
axes7 = axes7.flatten()

for bi, band in enumerate(BANDS):
    ax = axes7[bi]
    x_vals, y_vals, pt_labels = [], [], []

    for patient in PATIENTS:
        # Scale of max structural divergence (lowest cophenetic corr)
        df_pb_sc = df_scale[(df_scale["patient"] == patient) & (df_scale["freq_band"] == band)]
        if df_pb_sc.empty:
            continue
        avg_sc = df_pb_sc.groupby("scale_band_idx")["spearman_corr"].mean()
        if avg_sc.empty:
            continue
        scale_max_struct_div = avg_sc.idxmin()

        # Scale of max VI (which k has highest VI)
        df_pb_vi = df_vi[(df_vi["patient"] == patient) & (df_vi["freq_band"] == band)]
        if df_pb_vi.empty:
            continue
        avg_vi = df_pb_vi.groupby("k")["vi"].mean()
        if avg_vi.empty:
            continue
        k_max_vi = avg_vi.idxmax()

        # Map k to approximate scale band index:
        # k=2 -> band 4 (coarse), k=20 -> band 0 (fine)
        # Linear interpolation
        k_to_band = N_SCALE_BANDS - 1 - (k_max_vi - min(K_VALUES)) / (max(K_VALUES) - min(K_VALUES)) * (N_SCALE_BANDS - 1)

        x_vals.append(scale_max_struct_div)
        y_vals.append(k_to_band)
        pt_labels.append(patient)

    ax.scatter(x_vals, y_vals, s=100, c=band_colors[bi], edgecolors="black", zorder=3)
    for x, y, lab in zip(x_vals, y_vals, pt_labels):
        ax.annotate(lab, (x, y), fontsize=7, ha="left", va="bottom",
                   xytext=(3, 3), textcoords="offset points")

    # Diagonal reference
    ax.plot([0, N_SCALE_BANDS - 1], [0, N_SCALE_BANDS - 1], "k--", alpha=0.3)
    ax.set_xlim(-0.5, N_SCALE_BANDS - 0.5)
    ax.set_ylim(-0.5, N_SCALE_BANDS - 0.5)
    ax.set_xlabel("Scale of max structural div\n(cophenetic, 0=fine, 4=coarse)", fontsize=9)
    ax.set_ylabel("Scale of max partition div\n(VI, mapped to scale bands)", fontsize=9)
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

fig7.suptitle(
    "Do structural and partition changes happen at the same scale?\n"
    "Diagonal = aligned; off-diagonal = dissociation",
    fontsize=14, fontweight="bold", y=1.02
)
fig7.tight_layout()
fig7.savefig(OUTDIR / "fig7_dissociation_scatter.png", dpi=200, bbox_inches="tight")
plt.close(fig7)
print("  Saved fig7_dissociation_scatter.png")


# ===========================================================================
# FIGURE 8: Triple overlay (VI + 1-cophenetic_corr + Sackin change)
# ===========================================================================
print("\n=== Figure 8: Triple overlay ===")
fig8, axes8 = plt.subplots(2, 3, figsize=(18, 11))
axes8 = axes8.flatten()

for bi, band in enumerate(BANDS):
    ax = axes8[bi]

    # --- VI profile averaged across patients and phase pairs ---
    df_b_vi = df_vi[df_vi["freq_band"] == band]
    if df_b_vi.empty:
        ax.set_title(BRAIN_BAND_TEX_DICT[band])
        continue
    vi_mean = df_b_vi.groupby("k")["vi"].mean()
    vi_sem_v = df_b_vi.groupby("k")["vi"].apply(
        lambda x: sem(x.dropna()) if len(x.dropna()) > 1 else 0
    )

    ax.errorbar(vi_mean.index, vi_mean.values, yerr=vi_sem_v.values,
                color="darkviolet", marker="o", linewidth=2, capsize=3,
                label="VI(k)")
    ax.set_xlabel("k (clusters) / Scale band index", fontsize=10)
    ax.set_ylabel("VI (nats)", fontsize=10, color="darkviolet")
    ax.tick_params(axis="y", labelcolor="darkviolet")

    # --- 1 - cophenetic correlation ---
    ax_r = ax.twinx()
    df_b_sc = df_scale[df_scale["freq_band"] == band]
    if not df_b_sc.empty:
        sc_mean = df_b_sc.groupby("scale_band_idx")["spearman_corr"].mean()
        sc_sem_v = df_b_sc.groupby("scale_band_idx")["spearman_corr"].apply(
            lambda x: sem(x.dropna()) if len(x.dropna()) > 1 else 0
        )
        # 1 - corr = structural divergence
        struct_div_mean = 1.0 - sc_mean
        ax_r.errorbar(sc_mean.index, struct_div_mean.values, yerr=sc_sem_v.values,
                      color="darkorange", marker="s", linewidth=2, capsize=3,
                      linestyle="--", label="1 - coph. corr")
    ax_r.set_ylabel("1 - cophenetic corr", fontsize=10, color="darkorange")
    ax_r.tick_params(axis="y", labelcolor="darkorange")

    # --- Sackin change (as |delta Sackin| averaged across phase pairs) ---
    sackin_changes = []
    for patient in PATIENTS:
        sackin_vals = {}
        for phase in PHASES:
            srow = df_sackin[(df_sackin["patient"] == patient) &
                            (df_sackin["freq_band"] == band) &
                            (df_sackin["phase"] == phase)]
            if not srow.empty:
                sackin_vals[phase] = srow.iloc[0]["sackin"]
        # Compute average absolute Sackin change across phase pairs
        for p1, p2 in combinations(PHASES, 2):
            if p1 in sackin_vals and p2 in sackin_vals:
                sackin_changes.append(abs(sackin_vals[p1] - sackin_vals[p2]))

    if sackin_changes:
        mean_sackin_change = np.mean(sackin_changes)
        # Add as horizontal annotation band
        ax.axhspan(0, 0, alpha=0)  # placeholder
        ax.text(0.98, 0.02, f"|dSackin| = {mean_sackin_change:.3f}",
                transform=ax.transAxes, fontsize=9, ha="right", va="bottom",
                bbox=dict(facecolor="lightgreen", alpha=0.7, boxstyle="round,pad=0.3"))

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14)

    # Combined legend
    lines_a, labels_a = ax.get_legend_handles_labels()
    lines_b, labels_b = ax_r.get_legend_handles_labels()
    ax.legend(lines_a + lines_b, labels_a + labels_b, fontsize=8, loc="upper left")

fig8.suptitle(
    "Triple overlay: partition divergence (VI), structural divergence (1-coph corr), balance change (Sackin)\n"
    "Where do they align and where do they diverge?",
    fontsize=13, fontweight="bold", y=1.01
)
fig8.tight_layout()
fig8.savefig(OUTDIR / "fig8_triple_overlay.png", dpi=200, bbox_inches="tight")
plt.close(fig8)
print("  Saved fig8_triple_overlay.png")


# ===========================================================================
# FINDINGS
# ===========================================================================
print("\n=== Generating FINDINGS.md ===")
findings = []
findings.append("# Deep Multiscale Structural Analysis -- Findings\n")
findings.append(f"Patients: {', '.join(PATIENTS)}")
findings.append(f"Bands: {', '.join(BANDS)}")
findings.append(f"FC method: {FC_METHOD}")
findings.append(f"Scale bands: {N_SCALE_BANDS}")
findings.append(f"k values: {K_VALUES}\n")

# Finding 1: Where is cophenetic correlation lowest?
findings.append("## 1. Scale of maximum structural divergence (lowest cophenetic correlation)\n")
for band in BANDS:
    df_b = df_scale[df_scale["freq_band"] == band]
    if df_b.empty:
        continue
    avg = df_b.groupby("scale_band_idx")["spearman_corr"].mean()
    min_idx = avg.idxmin()
    min_val = avg.min()
    max_idx = avg.idxmax()
    max_val = avg.max()
    findings.append(f"- **{band}**: Lowest corr at scale band {min_idx} (corr={min_val:.3f}), "
                   f"highest at band {max_idx} (corr={max_val:.3f})")

# Finding 2: Within vs cross condition
findings.append("\n## 2. Within-condition vs cross-condition structural similarity\n")
for band in BANDS:
    df_b = df_scale[df_scale["freq_band"] == band]
    if df_b.empty:
        continue
    within_mean = df_b[df_b["pair_type"] == "within"]["spearman_corr"].mean()
    cross_mean = df_b[df_b["pair_type"] == "cross"]["spearman_corr"].mean()
    diff = within_mean - cross_mean
    findings.append(f"- **{band}**: Within={within_mean:.3f}, Cross={cross_mean:.3f}, "
                   f"diff={diff:+.3f}")

# Finding 3: VI peak scale
findings.append("\n## 3. Scale of maximum partition divergence (highest VI)\n")
for band in BANDS:
    df_b = df_vi[df_vi["freq_band"] == band]
    if df_b.empty:
        continue
    avg = df_b.groupby("k")["vi"].mean()
    max_k = avg.idxmax()
    max_vi = avg.max()
    findings.append(f"- **{band}**: Max VI at k={max_k} (VI={max_vi:.3f})")

# Finding 4: Merge profile divergence
findings.append("\n## 4. Merge profile divergence (L1 distance)\n")
for band in BANDS:
    df_b = df_merge[df_merge["freq_band"] == band]
    if df_b.empty:
        continue
    mean_l1 = df_b["l1_divergence"].mean()
    max_l1_row = df_b.loc[df_b["l1_divergence"].idxmax()]
    findings.append(f"- **{band}**: Mean L1={mean_l1:.4f}, "
                   f"max L1={max_l1_row['l1_divergence']:.4f} "
                   f"({max_l1_row['patient']} {max_l1_row['phase1']}-{max_l1_row['phase2']})")

# Finding 5: Sackin index
findings.append("\n## 5. Dendrogram balance (Sackin index)\n")
for band in BANDS:
    df_b = df_sackin[df_sackin["freq_band"] == band]
    if df_b.empty:
        continue
    # Phase-wise means
    phase_means = df_b.groupby("phase")["sackin"].mean()
    findings.append(f"- **{band}**: " + ", ".join(
        f"{p}={v:.3f}" for p, v in phase_means.items()))

# Finding 6: Dissociation summary
findings.append("\n## 6. Structural vs partition divergence dissociation\n")
for band in BANDS:
    # Structural
    df_b_sc = df_scale[df_scale["freq_band"] == band]
    df_b_vi = df_vi[df_vi["freq_band"] == band]
    if df_b_sc.empty or df_b_vi.empty:
        continue
    sc_avg = df_b_sc.groupby("scale_band_idx")["spearman_corr"].mean()
    struct_scale = sc_avg.idxmin()  # 0=fine, 4=coarse
    vi_avg = df_b_vi.groupby("k")["vi"].mean()
    vi_k_max = vi_avg.idxmax()
    # Map k to scale: k=2 -> coarse, k=20 -> fine
    findings.append(
        f"- **{band}**: Max structural div at scale band {struct_scale} "
        f"(0=fine,4=coarse); max VI at k={vi_k_max} "
        f"({'fine' if vi_k_max > 10 else 'meso' if vi_k_max > 5 else 'coarse'} scale)"
    )

# Finding 7: Which transitions cause biggest changes
findings.append("\n## 7. Biggest transitions (merge profile divergence)\n")
for p1, p2 in ALL_TRANSITIONS:
    df_t = df_merge[((df_merge["phase1"] == p1) & (df_merge["phase2"] == p2)) |
                    ((df_merge["phase1"] == p2) & (df_merge["phase2"] == p1))]
    mean_l1 = df_t["l1_divergence"].mean()
    findings.append(f"- {p1} -> {p2}: Mean L1 = {mean_l1:.4f}")

# Write findings
findings_text = "\n".join(findings) + "\n"
with open(OUTDIR / "FINDINGS.md", "w") as f:
    f.write(findings_text)
print("  Saved FINDINGS.md")

print("\n" + "=" * 60)
print("ANALYSIS COMPLETE")
print(f"All outputs in: {OUTDIR}")
print("=" * 60)
