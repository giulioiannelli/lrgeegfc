#!/usr/bin/env python3
"""Generate FC section figures that only need cached MSC matrices (fast).

Produces:
  Section 2.3 — FC Estimation:
    fig_adjacency_heatmaps_{patient}.pdf      (per patient, phases x bands)
    fig_adjacency_heatmaps_all_patients_rsPre.pdf
    fig_surrogate_comparison_{patient}.pdf

  Section 2.4 — Network Characterization:
    fig_weight_distributions_{patient}.pdf     (grid: 4 phases x 6 bands, log-binned + KDE)
    fig_weight_distributions_overlay_{patient}.pdf  (overlay: 6 bands, phases superposed)
    fig_strength_distribution_{patient}.pdf    (boxplots per band per phase)
    fig_network_metrics_summary_{patient}.pdf  (heatmaps: 8 metrics x 6 bands x 4 phases)
    fig_network_graph_{patient}_{phase}.pdf    (spatial PCA layout, 1x6 bands)
    fig_network_graph_spring_{patient}_{phase}.pdf  (spring layout, 1x6 bands)
    fig_network_metrics_cross_patient_{phase}.pdf   (heatmap: patients x bands, 8 metrics)
    fig_metrics_cross_patient_{phase}.pdf      (line plots: patients x bands, 4 key metrics)
    fig_eigenvalue_spectrum_{patient}.pdf       (Laplacian eigenvalues per band)

Run: python scripts/gen_fc_figures_fast.py
"""
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.collections import LineCollection
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA

from lrg_eegfc.config.const import (
    BRAIN_BANDS, BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PHASE_LABELS,
)
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.visuals.spatial import load_spatial_metadata, prepare_spatial_coordinates

# ── Config ────────────────────────────────────────────────────────────
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
PHASE = "rest_pre"
from lrg_eegfc.config.paths import MSC_CACHE, FC_FIG_CACHE, SEEG_DATAPATH, FIGURES_ROOT

COH_CACHE = FC_FIG_CACHE
DATASET_ROOT = SEEG_DATAPATH
NPERSEG = 4096

CMAP_MATRIX = "magma"
BAND_COLORS = {
    "delta": "#4C72B0", "theta": "#55A868", "alpha": "#C44E52",
    "beta": "#8172B3", "low_gamma": "#CCB974", "high_gamma": "#64B5CD",
}
PHASE_COLORS = {
    "rest_pre": "#1f77b4", "task_learn": "#ff7f0e",
    "task_test": "#2ca02c", "rest_post": "#d62728",
}
PHASE_MARKERS = {"rest_pre": "o", "task_learn": "s", "task_test": "D", "rest_post": "^"}


# ── Helpers ───────────────────────────────────────────────────────────
def load_msc_data_for_patient(patient, phases):
    """Load MSC matrices from cache, fallback to coherence tensor."""
    from lrg_eegfc.utils.fc.msc.msc import band_average_msc
    msc_data = {}
    for band in BRAIN_BANDS_NAMES:
        msc_data[band] = {}
        for phase in phases:
            mat = load_msc_matrix(
                patient, phase, band,
                cache_root=MSC_CACHE, sparsify="none", n_surrogates=0, nperseg=NPERSEG,
            )
            if mat is not None:
                msc_data[band][phase] = mat
    # Fallback: if nothing loaded, try coherence tensor
    loaded = sum(1 for b in msc_data for p in msc_data[b])
    if loaded == 0:
        for phase in phases:
            coh_path = COH_CACHE / f"{patient}_{phase}_coh_nperseg-{NPERSEG}.npz"
            if coh_path.exists():
                data = np.load(coh_path)
                W_bands = band_average_msc(data["Coh"], data["freqs"], BRAIN_BANDS)
                for band in BRAIN_BANDS_NAMES:
                    msc_data[band][phase] = W_bands[band]
    return msc_data


def exp_scale(t, vmin, vmax):
    """Exponential interpolation: t in [0,1] -> [vmin, vmax]."""
    return vmin * np.exp(t * np.log(vmax / vmin))


def compute_network_metrics(mat):
    """Compute network metrics from an adjacency matrix."""
    mat = mat.copy()
    np.fill_diagonal(mat, 0)
    triu_idx = np.triu_indices_from(mat, k=1)
    triu_vals = mat[triu_idx]

    strengths = mat.sum(axis=1)
    G = nx.from_numpy_array(mat)
    cc = nx.average_clustering(G, weight="weight")

    si = strengths[triu_idx[0]]
    sj = strengths[triu_idx[1]]
    strength_assort = np.corrcoef(si, sj)[0, 1]

    sorted_w = np.sort(triu_vals)[::-1]
    n_top5 = max(1, int(0.05 * len(sorted_w)))
    weight_conc = sorted_w[:n_top5].sum() / (sorted_w.sum() + 1e-12)

    return {
        "mean_weight": triu_vals.mean(),
        "median_weight": np.median(triu_vals),
        "density_gt005": np.mean(triu_vals > 0.05),
        "density_gt01": np.mean(triu_vals > 0.1),
        "mean_strength": strengths.mean(),
        "std_strength": strengths.std(),
        "clustering_coeff": cc,
        "strength_assort": strength_assort,
        "weight_concentration_5pct": weight_conc,
    }


# ======================================================================
# SEC 2.3, FIG: Adjacency heatmaps (bands x phases) — one per patient
# ======================================================================
n_bands = len(BRAIN_BANDS_NAMES)
n_phases = len(PHASE_LABELS)

for patient in PATIENTS:
    OUTPUT_DIR = FIGURES_ROOT / "report_fc_section" / patient
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Adjacency heatmaps ({patient})...")
    msc_data = load_msc_data_for_patient(patient, PHASE_LABELS)

    fig, axes = plt.subplots(
        n_phases, n_bands, figsize=(3.2 * n_bands, 3.0 * n_phases),
        constrained_layout=True,
    )

    for j, band in enumerate(BRAIN_BANDS_NAMES):
        for i, phase in enumerate(PHASE_LABELS):
            ax = axes[i, j]
            mat = msc_data.get(band, {}).get(phase)
            if mat is None:
                ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                        transform=ax.transAxes, fontsize=14, color="0.5")
                ax.set_xticks([]); ax.set_yticks([])
                continue
            im = ax.imshow(mat, cmap=CMAP_MATRIX, vmin=0, vmax=0.6,
                           aspect="equal", interpolation="none")
            ax.set_xticks([]); ax.set_yticks([])
            if j == 0:
                ax.set_ylabel(phase, fontsize=20, fontweight="bold")
            if i == 0:
                ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=22, fontweight="bold")

    cbar = fig.colorbar(im, ax=axes, shrink=0.95, pad=0.015, aspect=40)
    cbar.ax.tick_params(labelsize=16)

    out = OUTPUT_DIR / f"fig_adjacency_heatmaps_{patient}.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


# ── All-patients rest_pre comparison ─────────────────────────────────────
print("Adjacency heatmaps (all patients, rest_pre)...")
fig, axes = plt.subplots(
    len(PATIENTS), n_bands, figsize=(3.2 * n_bands, 3.0 * len(PATIENTS)),
    constrained_layout=True,
)
for i, patient in enumerate(PATIENTS):
    msc_data = load_msc_data_for_patient(patient, [PHASE])
    for j, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[i, j]
        mat = msc_data.get(band, {}).get(PHASE)
        if mat is None:
            ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                    transform=ax.transAxes, fontsize=14, color="0.5")
            ax.set_xticks([]); ax.set_yticks([])
            continue
        im = ax.imshow(mat, cmap=CMAP_MATRIX, vmin=0, vmax=0.6,
                       aspect="equal", interpolation="none")
        ax.set_xticks([]); ax.set_yticks([])
        if j == 0:
            ax.set_ylabel(patient, fontsize=20, fontweight="bold")
        if i == 0:
            ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=22, fontweight="bold")

cbar = fig.colorbar(im, ax=axes, shrink=0.95, pad=0.015, aspect=40)
cbar.ax.tick_params(labelsize=16)

out = FIGURES_ROOT / "report_fc_section" / "fig_adjacency_heatmaps_all_patients_rsPre.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ── Remaining figures use Pat_02 ──────────────────────────────────────
PATIENT = "Pat_02"
OUTPUT_DIR = FIGURES_ROOT / "report_fc_section" / PATIENT
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
msc_data = load_msc_data_for_patient(PATIENT, PHASE_LABELS)


# ======================================================================
# SEC 2.4, FIG 1: Weight distributions (log-binned + KDE)
# ======================================================================
# --- Grid version: 4 phases x 6 bands ---
print("Weight distributions (grid)...")
bins_log = np.logspace(-3, 0, 30)

fig, axes = plt.subplots(n_phases, n_bands, figsize=(3.5 * n_bands, 3.0 * n_phases),
                         constrained_layout=True)
for i, phase in enumerate(PHASE_LABELS):
    for j, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[i, j]
        mat = msc_data.get(band, {}).get(phase)
        if mat is None:
            ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                    transform=ax.transAxes, fontsize=11, color="0.5")
            ax.set_xscale("log"); ax.set_yscale("log")
            continue
        triu = mat[np.triu_indices_from(mat, k=1)]
        triu_pos = triu[triu > 0]

        # Histogram
        ax.hist(triu_pos, bins=bins_log, color="0.7", edgecolor="white", lw=0.3)

        # KDE in log space — scale to match histogram counts
        if len(triu_pos) > 10:
            log_vals = np.log10(triu_pos)
            kde = gaussian_kde(log_vals, bw_method=0.15)
            x_log = np.linspace(-3, 0, 200)
            x_lin = 10**x_log
            pdf_log = kde(x_log)
            # pdf_log integrates to 1 over log10(x).
            # Histogram bins are uniform in log10 with width delta_log.
            # Expected count per bin = N * delta_log * pdf_log(bin_center).
            delta_log = 3.0 / (len(bins_log) - 1)
            kde_counts = len(triu_pos) * delta_log * pdf_log
            ax.plot(x_lin, kde_counts, color="k", lw=1.5)

        # Median line
        ax.axvline(np.median(triu_pos), color="red", ls="--", lw=1.0)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-3, 1)
        if i == 0:
            ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14, fontweight="bold")
        if j == 0:
            ax.set_ylabel(phase, fontsize=14, fontweight="bold")
        if i < n_phases - 1:
            ax.set_xticklabels([])

out = OUTPUT_DIR / f"fig_weight_distributions_{PATIENT}.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")

# --- Overlay version: 6 bands, phases superposed ---
print("Weight distributions (overlay)...")
fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
axes_flat = axes.ravel()

for k, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes_flat[k]
    for phase in PHASE_LABELS:
        mat = msc_data.get(band, {}).get(phase)
        if mat is None:
            continue
        triu = mat[np.triu_indices_from(mat, k=1)]
        triu_pos = triu[triu > 0]
        if len(triu_pos) < 10:
            continue

        # Faint step histogram
        ax.hist(triu_pos, bins=bins_log, color=PHASE_COLORS[phase],
                alpha=0.15, histtype="stepfilled")

        # KDE — same log-space scaling as grid version
        log_vals = np.log10(triu_pos)
        kde = gaussian_kde(log_vals, bw_method=0.15)
        x_log = np.linspace(-3, 0, 200)
        x_lin = 10**x_log
        pdf_log = kde(x_log)
        delta_log = 3.0 / (len(bins_log) - 1)
        kde_counts = len(triu_pos) * delta_log * pdf_log
        ax.plot(x_lin, kde_counts, color=PHASE_COLORS[phase],
                lw=2.0, label=phase)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-3, 1)
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14, fontweight="bold")
    if k == 0:
        ax.legend(fontsize=9)
    if k >= 3:
        ax.set_xlabel("$A_{ij}$", fontsize=11)
    if k % 3 == 0:
        ax.set_ylabel("Count", fontsize=11)

out = OUTPUT_DIR / f"fig_weight_distributions_overlay_{PATIENT}.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# SEC 2.4, FIG 2: Node strength distribution (boxplots)
# ======================================================================
print("Strength distribution...")
fig, axes = plt.subplots(1, n_bands, figsize=(4 * n_bands, 5), constrained_layout=True)

for k, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes[k]
    box_data = []
    box_labels = []
    for phase in PHASE_LABELS:
        mat = msc_data.get(band, {}).get(phase)
        if mat is None:
            continue
        m = mat.copy()
        np.fill_diagonal(m, 0)
        box_data.append(m.sum(axis=1))
        box_labels.append(phase)

    if box_data:
        bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True, showfliers=False,
                        widths=0.6, medianprops=dict(color="black", lw=1.5))
        for patch, phase in zip(bp["boxes"], box_labels):
            patch.set_facecolor(PHASE_COLORS[phase])
            patch.set_alpha(0.7)

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14, fontweight="bold")
    ax.tick_params(axis="x", rotation=45)
    if k == 0:
        ax.set_ylabel(r"Node strength $s_i$", fontsize=12)
    ax.grid(axis="y", alpha=0.3)

out = OUTPUT_DIR / f"fig_strength_distribution_{PATIENT}.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# SEC 2.4, FIG 3: Network metrics summary (heatmaps)
# ======================================================================
print("Network metrics summary...")
METRIC_NAMES = [
    "mean_weight", "median_weight", "density_gt005", "density_gt01",
    "mean_strength", "std_strength", "clustering_coeff", "strength_assort",
]
METRIC_LABELS = [
    "Mean weight", "Median weight", "Density (>0.05)", "Density (>0.1)",
    "Mean strength", "Std strength", "Clustering coeff.", "Strength assort.",
]

# Compute metrics for all band x phase
records = []
for band in BRAIN_BANDS_NAMES:
    for phase in PHASE_LABELS:
        mat = msc_data.get(band, {}).get(phase)
        if mat is None:
            continue
        metrics = compute_network_metrics(mat)
        metrics["band"] = band
        metrics["phase"] = phase
        records.append(metrics)

df_metrics = pd.DataFrame(records)

fig, axes = plt.subplots(2, 4, figsize=(20, 10), constrained_layout=True)
axes_flat = axes.ravel()

for m_idx, (metric, label) in enumerate(zip(METRIC_NAMES, METRIC_LABELS)):
    ax = axes_flat[m_idx]
    grid = np.full((n_phases, n_bands), np.nan)
    for i, phase in enumerate(PHASE_LABELS):
        for j, band in enumerate(BRAIN_BANDS_NAMES):
            row = df_metrics[(df_metrics["band"] == band) & (df_metrics["phase"] == phase)]
            if len(row) > 0:
                grid[i, j] = row[metric].values[0]

    ax.imshow(grid, cmap="viridis", aspect="auto")
    ax.set_box_aspect(1)

    # Annotate cells
    for i in range(n_phases):
        for j in range(n_bands):
            if not np.isnan(grid[i, j]):
                val = grid[i, j]
                fmt = f"{val:.3f}" if abs(val) < 10 else f"{val:.1f}"
                ax.text(j, i, fmt, ha="center", va="center", fontsize=8,
                        color="white" if val < (np.nanmax(grid) + np.nanmin(grid)) / 2 else "black")

    ax.set_xticks(range(n_bands))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES], fontsize=9)
    ax.set_yticks(range(n_phases))
    ax.set_yticklabels(PHASE_LABELS, fontsize=9, rotation=90, va="center")
    ax.set_title(label, fontsize=11, fontweight="bold")

out = OUTPUT_DIR / f"fig_network_metrics_summary_{PATIENT}.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# SEC 2.4, FIG 4: Network graph — spatial layout (PCA projection)
# ======================================================================
print("Network graph (spatial)...")
metadata = load_spatial_metadata(PATIENT, DATASET_ROOT)
coords_3d = metadata[["x", "y", "z"]].values.astype(float)
valid = ~np.isnan(coords_3d).any(axis=1)

matrices = {}
for band in BRAIN_BANDS_NAMES:
    matrices[band] = load_msc_matrix(PATIENT, PHASE, band, cache_root=MSC_CACHE,
                                      sparsify="none", n_surrogates=0, nperseg=NPERSEG)

N = matrices["beta"].shape[0]
pca = PCA(n_components=2)
coords_2d = pca.fit_transform(coords_3d[valid])
pos_pca = np.full((N, 2), np.nan)
pos_pca[valid] = coords_2d

fig, axes = plt.subplots(1, 6, figsize=(30, 5))

for k, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes[k]
    mat = matrices[band].copy()
    np.fill_diagonal(mat, 0)

    strengths = mat.sum(axis=1)
    s_norm = (strengths - strengths.min()) / (strengths.max() - strengths.min() + 1e-12)

    triu_r, triu_c = np.triu_indices(N, k=1)
    weights = mat[triu_r, triu_c]
    w_min_val, w_max_val = weights.min(), weights.max()
    w_norm = (weights - w_min_val) / (w_max_val - w_min_val + 1e-12)

    widths = exp_scale(w_norm, 0.15, 12.0)
    alphas = exp_scale(w_norm, 0.08, 0.8)

    sort_idx = np.argsort(widths)
    segments, edge_colors, edge_widths = [], [], []
    for idx in sort_idx:
        i, j = triu_r[idx], triu_c[idx]
        if np.isnan(pos_pca[i]).any() or np.isnan(pos_pca[j]).any():
            continue
        segments.append([pos_pca[i], pos_pca[j]])
        edge_colors.append((0.3, 0.3, 0.3, alphas[idx]))
        edge_widths.append(widths[idx])

    lc = LineCollection(segments, colors=edge_colors, linewidths=edge_widths,
                        zorder=1, rasterized=True)
    ax.add_collection(lc)

    node_sizes = 20 + 180 * s_norm
    scatter = ax.scatter(pos_pca[valid, 0], pos_pca[valid, 1],
                         c=strengths[valid], cmap="viridis",
                         s=node_sizes[valid], edgecolors="white", linewidths=0.5,
                         zorder=5, vmin=strengths.min(), vmax=strengths.max())

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=16, fontweight="bold")
    ax.axis("off")
    ax.set_aspect("equal")

    if k == 5:
        cbar = fig.colorbar(scatter, ax=ax, shrink=0.7, pad=0.02)
        cbar.set_label("Node strength", fontsize=11)
        cbar.ax.tick_params(labelsize=9)

fig.subplots_adjust(wspace=0.05)
out = OUTPUT_DIR / f"fig_network_graph_{PATIENT}_{PHASE}.pdf"
fig.savefig(out, bbox_inches="tight", dpi=200)
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# SEC 2.4, FIG 4b: Network graph — spring layout
# ======================================================================
print("Network graph (spring)...")
K_SPRING = {
    "delta": 5, "theta": 6, "alpha": 8,
    "beta": 18, "low_gamma": 25, "high_gamma": 30,
}

fig, axes = plt.subplots(1, 6, figsize=(30, 5))

for k_idx, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes[k_idx]
    mat = matrices[band].copy()
    np.fill_diagonal(mat, 0)

    strengths = mat.sum(axis=1)
    s_norm = (strengths - strengths.min()) / (strengths.max() - strengths.min() + 1e-12)

    G = nx.from_numpy_array(mat)
    k_val = K_SPRING[band] / np.sqrt(N)
    pos = nx.spring_layout(G, k=k_val, iterations=500, seed=42, weight="weight")
    pos_arr = np.array([pos[i] for i in range(N)])

    triu_r, triu_c = np.triu_indices(N, k=1)
    weights = mat[triu_r, triu_c]

    # Rank-based normalization: uniform quantile mapping
    ranks = np.argsort(np.argsort(weights)).astype(float)
    t = ranks / (len(ranks) - 1)

    widths = exp_scale(t, 0.5, 6.0)
    grey_vals = 0.82 - 0.65 * t  # light grey -> near-black

    sort_idx = np.argsort(t)
    segments, edge_colors, edge_widths = [], [], []
    for idx in sort_idx:
        i, j = triu_r[idx], triu_c[idx]
        g = grey_vals[idx]
        segments.append([pos_arr[i], pos_arr[j]])
        edge_colors.append((g, g, g, 1.0))
        edge_widths.append(widths[idx])

    lc = LineCollection(segments, colors=edge_colors, linewidths=edge_widths,
                        zorder=1, rasterized=True)
    ax.add_collection(lc)

    node_sizes = 20 + 180 * s_norm
    scatter = ax.scatter(pos_arr[:, 0], pos_arr[:, 1],
                         c=strengths, cmap="viridis",
                         s=node_sizes, edgecolors="white", linewidths=0.5,
                         zorder=5, vmin=strengths.min(), vmax=strengths.max())

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=16, fontweight="bold")
    ax.axis("off")

    if k_idx == 5:
        cbar = fig.colorbar(scatter, ax=ax, shrink=0.7, pad=0.02)
        cbar.set_label("Node strength", fontsize=11)
        cbar.ax.tick_params(labelsize=9)

    xmin, xmax = pos_arr[:, 0].min(), pos_arr[:, 0].max()
    ymin, ymax = pos_arr[:, 1].min(), pos_arr[:, 1].max()
    margin = 0.08 * max(xmax - xmin, ymax - ymin)
    ax.set_xlim(xmin - margin, xmax + margin)
    ax.set_ylim(ymin - margin, ymax + margin)

fig.subplots_adjust(wspace=0.05)
out = OUTPUT_DIR / f"fig_network_graph_spring_{PATIENT}_{PHASE}.pdf"
fig.savefig(out, bbox_inches="tight", dpi=200)
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# SEC 2.4, FIG 4c: Brain connectome (nilearn glass brain, axial view)
# ======================================================================
print("Brain connectome (glass brain)...")
from nilearn.plotting import plot_connectome as nilearn_plot_connectome

coords_mni = prepare_spatial_coordinates(metadata, scale="mm", center=False, to_mni=True)

fig, axes = plt.subplots(1, 6, figsize=(30, 6))

for k, band in enumerate(BRAIN_BANDS_NAMES):
    mat = matrices[band].copy()
    np.fill_diagonal(mat, 0)

    # Show top 15% of edges for clarity
    triu = mat[np.triu_indices_from(mat, k=1)]
    threshold = np.percentile(triu, 85)

    strengths = mat.sum(axis=1)
    s_norm = (strengths - strengths.min()) / (strengths.max() - strengths.min() + 1e-12)
    node_sizes = 15 + 80 * s_norm

    ax = axes[k]
    nilearn_plot_connectome(
        adjacency_matrix=mat,
        node_coords=coords_mni,
        edge_threshold=threshold,
        edge_cmap="Greys",
        edge_vmin=0, edge_vmax=mat.max(),
        node_color=strengths,
        node_size=node_sizes,
        display_mode="z",
        colorbar=False,
        title=None,
        axes=ax,
    )
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=16, fontweight="bold")

out = OUTPUT_DIR / f"fig_brain_connectome_{PATIENT}_{PHASE}.pdf"
fig.savefig(out, bbox_inches="tight", dpi=200)
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# SEC 2.4, FIG 5: Cross-patient metric comparison (heatmap)
# ======================================================================
print("Cross-patient metrics (heatmap)...")
cross_records = []
for patient in PATIENTS:
    p_data = load_msc_data_for_patient(patient, [PHASE])
    for band in BRAIN_BANDS_NAMES:
        mat = p_data.get(band, {}).get(PHASE)
        if mat is None:
            continue
        metrics = compute_network_metrics(mat)
        metrics["patient"] = patient
        metrics["band"] = band
        cross_records.append(metrics)

df_cross = pd.DataFrame(cross_records)

fig, axes = plt.subplots(2, 4, figsize=(20, 8), constrained_layout=True)
axes_flat = axes.ravel()

for m_idx, (metric, label) in enumerate(zip(METRIC_NAMES, METRIC_LABELS)):
    ax = axes_flat[m_idx]
    grid = np.full((len(PATIENTS), n_bands), np.nan)
    for i, patient in enumerate(PATIENTS):
        for j, band in enumerate(BRAIN_BANDS_NAMES):
            row = df_cross[(df_cross["patient"] == patient) & (df_cross["band"] == band)]
            if len(row) > 0:
                grid[i, j] = row[metric].values[0]

    ax.imshow(grid, cmap="viridis", aspect="auto")
    ax.set_box_aspect(1)

    for i in range(len(PATIENTS)):
        for j in range(n_bands):
            if not np.isnan(grid[i, j]):
                val = grid[i, j]
                fmt = f"{val:.3f}" if abs(val) < 10 else f"{val:.1f}"
                ax.text(j, i, fmt, ha="center", va="center", fontsize=8,
                        color="white" if val < (np.nanmax(grid) + np.nanmin(grid)) / 2 else "black")

    ax.set_xticks(range(n_bands))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES], fontsize=9)
    ax.set_yticks(range(len(PATIENTS)))
    ax.set_yticklabels(PATIENTS, fontsize=9, rotation=90, va="center")
    ax.set_title(label, fontsize=11, fontweight="bold")

out = FIGURES_ROOT / "report_fc_section" / f"fig_network_metrics_cross_patient_{PHASE}.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ── Cross-patient line plots (4 key metrics) ─────────────────────────
print("Cross-patient metrics (line plots)...")
KEY_METRICS = [
    ("mean_weight", "Mean edge weight"),
    ("clustering_coeff", "Clustering coeff."),
    ("strength_assort", "Strength assort."),
    ("weight_concentration_5pct", "Weight conc. (5%)"),
]

fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
axes_flat = axes.ravel()

for m_idx, (metric, label) in enumerate(KEY_METRICS):
    ax = axes_flat[m_idx]
    x = np.arange(n_bands)
    for patient in PATIENTS:
        vals = []
        for band in BRAIN_BANDS_NAMES:
            row = df_cross[(df_cross["patient"] == patient) & (df_cross["band"] == band)]
            vals.append(row[metric].values[0] if len(row) > 0 else np.nan)
        ax.plot(x, vals, marker="o", label=patient, lw=1.5, markersize=5)
    ax.set_xticks(x)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES], fontsize=10)
    ax.set_ylabel(label, fontsize=11)
    ax.set_title(label, fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)

out = FIGURES_ROOT / "report_fc_section" / f"fig_metrics_cross_patient_{PHASE}.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# SEC 2.4, FIG 6: Laplacian eigenvalue spectrum
# ======================================================================
print("Eigenvalue spectrum...")
fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
axes_flat = axes.ravel()

for k, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes_flat[k]
    mat = matrices[band].copy()
    np.fill_diagonal(mat, 0)

    D = np.diag(mat.sum(axis=1))
    L = D - mat
    eigvals = np.sort(np.linalg.eigvalsh(L))

    ax.plot(range(1, len(eigvals) + 1), eigvals, color=BAND_COLORS[band], lw=1.5)
    ax.axhline(eigvals[1], color="red", ls="--", lw=1.0, alpha=0.7,
               label=f"$\\lambda_2$ = {eigvals[1]:.3f}")
    ax.axhline(eigvals[-1], color="grey", ls=":", lw=1.0, alpha=0.7,
               label=f"$\\lambda_{{max}}$ = {eigvals[-1]:.1f}")

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14, fontweight="bold")
    ax.set_xlabel("Index" if k >= 3 else "", fontsize=11)
    ax.set_ylabel("$\\lambda$" if k % 3 == 0 else "", fontsize=11)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

out = OUTPUT_DIR / f"fig_eigenvalue_spectrum_{PATIENT}.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out}")


# ======================================================================
# SEC 2.3, FIG: Surrogate comparison (raw vs attenuated)
# ======================================================================
print("Surrogate comparison...")
BAND_CMP = "beta"
A_dense = load_msc_matrix(PATIENT, PHASE, BAND_CMP,
                          cache_root=MSC_CACHE, sparsify="none", n_surrogates=0, nperseg=NPERSEG)
A_soft = load_msc_matrix(PATIENT, PHASE, BAND_CMP,
                         cache_root=MSC_CACHE, sparsify="soft", n_surrogates=100, nperseg=NPERSEG)

if A_soft is not None:
    triu_idx = np.triu_indices_from(A_dense, k=1)
    dense_vals = A_dense[triu_idx]
    soft_vals = A_soft[triu_idx]
    diff_vals = dense_vals - soft_vals
    r = np.corrcoef(dense_vals, soft_vals)[0, 1]
    vmax_shared = 0.6

    fig = plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[1, 0.8], hspace=0.35, wspace=0.25)

    ax1 = fig.add_subplot(gs[0, 0])
    im1 = ax1.imshow(A_dense, cmap=CMAP_MATRIX, vmin=0, vmax=vmax_shared,
                     aspect="equal", interpolation="none")
    ax1.set_title(r"(a) Raw $A(\mathcal{B})$", fontsize=12, fontweight="bold")
    ax1.set_xlabel("Channel"); ax1.set_ylabel("Channel")

    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.imshow(A_soft, cmap=CMAP_MATRIX, vmin=0, vmax=vmax_shared,
                     aspect="equal", interpolation="none")
    ax2.set_title(r"(b) Attenuated $\hat{A}(\mathcal{B})$", fontsize=12, fontweight="bold")
    ax2.set_xlabel("Channel"); ax2.set_yticks([])

    cbar1 = fig.colorbar(im2, ax=[ax1, ax2], shrink=0.8, pad=0.02)
    cbar1.set_label("MSC", fontsize=11)

    ax3 = fig.add_subplot(gs[0, 2])
    ax3.scatter(dense_vals, soft_vals, s=1, alpha=0.3, color="0.4", rasterized=True)
    ax3.plot([0, vmax_shared], [0, vmax_shared], "r--", lw=1, label="identity")
    ax3.set_xlabel(r"Raw $A_{ij}$", fontsize=11)
    ax3.set_ylabel(r"Attenuated $\hat{A}_{ij}$", fontsize=11)
    ax3.set_title("(c) Edge-wise comparison", fontsize=12, fontweight="bold")
    ax3.set_xlim(0, vmax_shared); ax3.set_ylim(0, vmax_shared)
    ax3.set_aspect("equal")
    ax3.legend(fontsize=9)
    ax3.text(0.05, 0.92, f"$r$ = {r:.6f}", transform=ax3.transAxes, fontsize=10,
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8))

    ax4 = fig.add_subplot(gs[1, 0:2])
    ax4.hist(diff_vals, bins=80, color="0.6", edgecolor="white", lw=0.3, log=True)
    ax4.axvline(0, color="k", ls="-", lw=0.8)
    ax4.axvline(diff_vals.mean(), color="red", ls="--", lw=1.5,
                label=f"mean = {diff_vals.mean():.2e}")
    ax4.set_xlabel(r"$A_{ij} - \hat{A}_{ij}$", fontsize=11)
    ax4.set_ylabel("Count (log)", fontsize=11)
    ax4.set_title(r"(d) Distribution of edge-wise differences", fontsize=12, fontweight="bold")
    ax4.legend(fontsize=9)

    ax5 = fig.add_subplot(gs[1, 2])
    ax5.axis("off")
    stats_text = (
        f"Edge-wise difference statistics\n"
        f"{'─' * 35}\n"
        f"Mean diff:      {diff_vals.mean():.2e}\n"
        f"Max diff:       {diff_vals.max():.2e}\n"
        f"Std diff:       {diff_vals.std():.2e}\n"
        f"{'─' * 35}\n"
        f"Pearson r:      {r:.6f}\n"
        f"Edges changed\n"
        f"  (>1% diff):   {np.mean(np.abs(diff_vals) > 0.01):.1%}\n"
        f"  (>5% diff):   {np.mean(np.abs(diff_vals) > 0.05):.1%}\n"
        f"{'─' * 35}\n"
        f"Conclusion: soft sparsification\n"
        f"has negligible effect on this\n"
        f"recording ({PATIENT}, {PHASE}, {BAND_CMP})"
    )
    ax5.text(0.1, 0.95, stats_text, transform=ax5.transAxes, fontsize=10,
             verticalalignment="top", fontfamily="monospace",
             bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.9))

    out = OUTPUT_DIR / f"fig_surrogate_comparison_{PATIENT}.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
else:
    print("  Skipped (no soft-sparsified matrix found)")


print(f"\nDone! Figures saved to {ROOT / 'data' / 'figures' / 'report_fc_section'}")
