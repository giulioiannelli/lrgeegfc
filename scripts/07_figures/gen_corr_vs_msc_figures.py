#!/usr/bin/env python3
"""Generate figures for Section 2.5: Correlation vs MSC comparison.

Produces 5 PDF figures comparing Pearson correlation and MSC-based FC.

  Fig 1: fig_weight_distributions_corr_vs_msc.pdf  (2x3, overlaid histograms)
  Fig 2: fig_scatter_corr_vs_msc.pdf                (2x3, hexbin scatter)
  Fig 3: fig_adjacency_side_by_side.pdf              (2x6, heatmaps)
  Fig 4: fig_thresholding_sensitivity.pdf            (2x2, metrics vs threshold)
  Fig 5: fig_node_strength_distributions.pdf         (2x3, KDE of node strengths)

Run: python scripts/gen_corr_vs_msc_figures.py
"""
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.stats import gaussian_kde, pearsonr, spearmanr
import networkx as nx

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT,
)
from lrg_eegfc.workflow.corr import compute_corr_matrix
from lrg_eegfc.workflow.msc import load_msc_matrix

# ── Config ────────────────────────────────────────────────────────────
PATIENT = "Pat_02"
PHASE = "rest_pre"
from lrg_eegfc.config.paths import MSC_CACHE, FIGURES_ROOT
NPERSEG = 4096
OUTPUT_DIR = FIGURES_ROOT / "corr_vs_msc" / PATIENT
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

BAND_COLORS = {
    "delta": "#4C72B0", "theta": "#55A868", "alpha": "#C44E52",
    "beta": "#8172B3", "low_gamma": "#CCB974", "high_gamma": "#64B5CD",
}
n_bands = len(BRAIN_BANDS_NAMES)


def upper_tri(mat):
    """Extract upper-triangle values (k=1, excludes diagonal)."""
    return mat[np.triu_indices_from(mat, k=1)]


# ======================================================================
# STEP 0: Compute / load all matrices
# ======================================================================
print("=" * 60)
print("Computing correlation matrices (filter_type='none')...")
print("=" * 60)
corr_raw = {}
for band in BRAIN_BANDS_NAMES:
    result = compute_corr_matrix(
        PATIENT, PHASE, band,
        filter_type="none", zero_diagonal=True,
        filter_order=4, verbose=True,
    )
    corr_raw[band] = result.adjacency_matrix
    print(f"  {band}: shape={result.adjacency_matrix.shape}, mean={result.mean_corr:.4f}")

print("\nComputing correlation matrices (filter_type='abs')...")
corr_abs = {}
for band in BRAIN_BANDS_NAMES:
    result = compute_corr_matrix(
        PATIENT, PHASE, band,
        filter_type="abs", zero_diagonal=True,
        filter_order=4, verbose=True,
    )
    corr_abs[band] = result.adjacency_matrix
    print(f"  {band}: shape={result.adjacency_matrix.shape}, mean={result.mean_corr:.4f}")

print("\nLoading MSC matrices...")
msc = {}
for band in BRAIN_BANDS_NAMES:
    mat = load_msc_matrix(
        PATIENT, PHASE, band,
        cache_root=MSC_CACHE, sparsify="none", n_surrogates=0, nperseg=NPERSEG,
    )
    assert mat is not None, f"MSC cache missing for {band}"
    msc[band] = mat
    print(f"  {band}: shape={mat.shape}, mean={upper_tri(mat).mean():.4f}")

N = msc["beta"].shape[0]
print(f"\nAll matrices loaded. N={N} channels, {N*(N-1)//2} edges.\n")


# ======================================================================
# FIG 1: Weight distributions — correlation vs MSC per band
# ======================================================================
print("Fig 1: weight distributions (corr vs MSC)...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)

bins_corr = np.linspace(-1, 1, 60)
bins_msc = np.linspace(0, 1, 40)

for k, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes.ravel()[k]
    c_vals = upper_tri(corr_raw[band])
    m_vals = upper_tri(msc[band])

    ax.hist(c_vals, bins=bins_corr, density=True, alpha=0.5, color="steelblue",
            edgecolor="none")
    ax.hist(m_vals, bins=bins_msc, density=True, alpha=0.5, color="coral",
            edgecolor="none")

    ax.axvline(0, color="black", ls=":", lw=0.8, alpha=0.5)

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=16, fontweight="bold")
    if k >= 3:
        ax.set_xlabel("Edge weight", fontsize=13)
    if k % 3 == 0:
        ax.set_ylabel("Density", fontsize=13)
    ax.grid(alpha=0.3)

# Shared legend at top — push up to avoid overlapping band titles
legend_handles = [
    Line2D([0], [0], color="steelblue", lw=8, alpha=0.5, label="Correlation"),
    Line2D([0], [0], color="coral", lw=8, alpha=0.5, label="MSC"),
]
fig.legend(handles=legend_handles, loc="upper center", ncol=2, fontsize=14,
           frameon=True, fancybox=True, bbox_to_anchor=(0.5, 1.06))

out = OUTPUT_DIR / "fig_weight_distributions_corr_vs_msc.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# FIG 2: Scatter — correlation vs MSC per band
# ======================================================================
print("Fig 2: scatter corr vs MSC...")
from matplotlib.colors import LogNorm

fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)

hb_last = None
for k, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes.ravel()[k]
    c_vals = upper_tri(corr_raw[band])
    m_vals = upper_tri(msc[band])

    hb = ax.hexbin(c_vals, m_vals, gridsize=80, cmap="inferno", mincnt=1,
                   norm=LogNorm(), rasterized=True)
    hb_last = hb

    ax.axhline(0, color="grey", ls=":", lw=0.8)
    ax.axvline(0, color="grey", ls=":", lw=0.8)
    ax.plot([0, 1], [0, 1], "w--", lw=1, alpha=0.6)

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=18, fontweight="bold")
    x_min = min(c_vals.min(), -0.05)
    ax.set_xlim(x_min, 1)
    ax.set_ylim(0, 1)
    if k >= 3:
        ax.set_xlabel("Correlation", fontsize=16)
    if k % 3 == 0:
        ax.set_ylabel("MSC", fontsize=16)
    ax.tick_params(labelsize=12)

# Shared colorbar
cbar = fig.colorbar(hb_last, ax=axes, shrink=0.6, pad=0.02, aspect=30)
cbar.set_label("Edge count (log scale)", fontsize=14)
cbar.ax.tick_params(labelsize=12)

out = OUTPUT_DIR / "fig_scatter_corr_vs_msc.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# FIG 3: Side-by-side adjacency matrices
# ======================================================================
print("Fig 3: adjacency side-by-side...")
fig, axes = plt.subplots(2, n_bands, figsize=(3.2 * n_bands, 3.0 * 2),
                         constrained_layout=True)

for j, band in enumerate(BRAIN_BANDS_NAMES):
    # Top row: raw correlation
    ax_c = axes[0, j]
    im_c = ax_c.imshow(corr_raw[band], cmap="RdBu_r", vmin=-1, vmax=1,
                       aspect="equal", interpolation="none")
    ax_c.set_xticks([]); ax_c.set_yticks([])
    ax_c.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=20, fontweight="bold")
    if j == 0:
        ax_c.set_ylabel("Correlation", fontsize=20, fontweight="bold")

    # Bottom row: MSC
    ax_m = axes[1, j]
    im_m = ax_m.imshow(msc[band], cmap="magma", vmin=0, vmax=0.6,
                       aspect="equal", interpolation="none")
    ax_m.set_xticks([]); ax_m.set_yticks([])
    if j == 0:
        ax_m.set_ylabel("MSC", fontsize=20, fontweight="bold")

cbar_c = fig.colorbar(im_c, ax=axes[0, :], shrink=0.95, pad=0.015, aspect=25)
cbar_c.ax.tick_params(labelsize=16)
cbar_m = fig.colorbar(im_m, ax=axes[1, :], shrink=0.95, pad=0.015, aspect=25)
cbar_m.ax.tick_params(labelsize=16)

out = OUTPUT_DIR / "fig_adjacency_side_by_side.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# FIG 4: Thresholding sensitivity
# ======================================================================
print("Fig 4: thresholding sensitivity...")

thresholds = np.linspace(0, 0.8, 50)


def compute_threshold_metrics(mat, thresholds):
    """Compute graph metrics at each threshold level."""
    n = mat.shape[0]
    densities, n_comps, largest_fracs, clusterings = [], [], [], []
    for theta in thresholds:
        binary = (mat > theta).astype(float)
        np.fill_diagonal(binary, 0)
        triu = upper_tri(binary)
        densities.append(triu.mean())

        G = nx.from_numpy_array(binary)
        components = list(nx.connected_components(G))
        n_comps.append(len(components))
        largest_fracs.append(max(len(c) for c in components) / n)
        clusterings.append(nx.average_clustering(G))

    return {
        "density": np.array(densities),
        "n_components": np.array(n_comps),
        "largest_component": np.array(largest_fracs),
        "clustering": np.array(clusterings),
    }


# Precompute for all bands
metrics_msc, metrics_corr = {}, {}
for band in BRAIN_BANDS_NAMES:
    print(f"  Computing threshold metrics for {band}...")
    metrics_msc[band] = compute_threshold_metrics(msc[band], thresholds)
    metrics_corr[band] = compute_threshold_metrics(corr_abs[band], thresholds)

metric_panels = [
    ("density", "Edge density"),
    ("n_components", "Connected components"),
    ("largest_component", "Largest component (frac.)"),
    ("clustering", "Clustering coefficient"),
]

fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

for m_idx, (metric, ylabel) in enumerate(metric_panels):
    ax = axes.ravel()[m_idx]
    for band in BRAIN_BANDS_NAMES:
        color = BAND_COLORS[band]
        ax.plot(thresholds, metrics_msc[band][metric],
                color=color, ls="-", lw=1.5)
        ax.plot(thresholds, metrics_corr[band][metric],
                color=color, ls="--", lw=1.5, alpha=0.7)

    ax.set_xlabel(r"Threshold $\theta$", fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(ylabel, fontsize=14, fontweight="bold")
    ax.grid(alpha=0.3)

# Shared legend outside
style_handles = [
    Line2D([0], [0], color="0.4", ls="-", lw=2, label="MSC"),
    Line2D([0], [0], color="0.4", ls="--", lw=2, label="|Corr|"),
]
band_handles = [
    Line2D([0], [0], color=BAND_COLORS[b], ls="-", lw=3,
           label=BRAIN_BAND_TEX_DICT[b])
    for b in BRAIN_BANDS_NAMES
]
fig.legend(handles=style_handles + band_handles, loc="upper center",
           ncol=8, fontsize=12, frameon=True, fancybox=True,
           bbox_to_anchor=(0.5, 1.05))

out = OUTPUT_DIR / "fig_thresholding_sensitivity.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# FIG 5: Node strength distributions
# ======================================================================
print("Fig 5: node strength distributions...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)

for k, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes.ravel()[k]

    # Compute node strengths
    cr = corr_raw[band].copy(); np.fill_diagonal(cr, 0)
    ca = corr_abs[band].copy(); np.fill_diagonal(ca, 0)
    mm = msc[band].copy(); np.fill_diagonal(mm, 0)

    s_raw = cr.sum(axis=1)
    s_abs = ca.sum(axis=1)
    s_msc = mm.sum(axis=1)

    for strengths, label, color in [
        (s_raw, "Corr (raw)", "steelblue"),
        (s_abs, "|Corr|", "darkorange"),
        (s_msc, "MSC", "coral"),
    ]:
        x_grid = np.linspace(strengths.min() - 1, strengths.max() + 1, 200)
        kde = gaussian_kde(strengths, bw_method=0.3)
        ax.plot(x_grid, kde(x_grid), color=color, lw=2)
        ax.fill_between(x_grid, kde(x_grid), alpha=0.1, color=color)

    ax.axvline(0, color="black", ls=":", lw=0.8, alpha=0.5)
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=16, fontweight="bold")
    if k >= 3:
        ax.set_xlabel("Node strength $s_i$", fontsize=13)
    if k % 3 == 0:
        ax.set_ylabel("Density", fontsize=13)
    ax.grid(alpha=0.3)

# Shared legend at top
legend_handles = [
    Line2D([0], [0], color="steelblue", lw=3, label="Corr (raw)"),
    Line2D([0], [0], color="darkorange", lw=3, label="|Corr|"),
    Line2D([0], [0], color="coral", lw=3, label="MSC"),
]
fig.legend(handles=legend_handles, loc="upper center", ncol=3, fontsize=14,
           frameon=True, fancybox=True, bbox_to_anchor=(0.5, 1.06))

out = OUTPUT_DIR / "fig_node_strength_distributions.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


print(f"\nDone! All figures saved to {OUTPUT_DIR}")
