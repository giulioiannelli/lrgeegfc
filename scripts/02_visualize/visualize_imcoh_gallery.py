#!/usr/bin/env python3
"""Comprehensive MSC vs ImCoh visualization gallery.

Generates six figure panels comparing MSC and ImCoh FC methods:

1. Adjacency matrix heatmaps (side-by-side MSC vs ImCoh)
2. Entropy curves (S(tau) and C(tau) overlaid)
3. Eigenvalue spectra (Laplacian eigenvalues)
4. Dendrogram comparison (side-by-side)
5. Same-probe ratio across all patients and bands (summary bar chart)
6. Community structure at n=10 (node membership comparison)

Outputs saved to data/imcoh_figures/gallery/.

Run: python scripts/py/visualize_imcoh_gallery.py [-v]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.cluster.hierarchy import dendrogram, fcluster, optimal_leaf_ordering
from scipy.spatial.distance import squareform

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PATIENTS_4PHASE,
    PHASE_LABELS,
)
from lrg_eegfc.config.paths import (
    FIGURES_ROOT,
    IMCOH_CACHE as _IMCOH_CACHE,
    IMCOH_LRG_CACHE as _IMCOH_LRG_CACHE,
    LRG_CACHE as _LRG_CACHE,
    MSC_CACHE as _MSC_CACHE,
    SEEG_DATAPATH,
)
from lrg_eegfc.workflow.lrg import LRGResult

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PATIENTS = list(PATIENTS_4PHASE)
BANDS = list(BRAIN_BANDS_NAMES)
BAND_TEX = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]
PHASES = list(PHASE_LABELS)

MSC_CACHE = _MSC_CACHE
IMCOH_CACHE = _IMCOH_CACHE
MSC_LRG_CACHE = _LRG_CACHE
IMCOH_LRG_CACHE = _IMCOH_LRG_CACHE
DATASET_ROOT = SEEG_DATAPATH
OUT = FIGURES_ROOT / "imcoh" / "gallery"

from lrg_eegfc.config.const import FS_OVERRIDES as FS_MAP  # canonical

# Plotting defaults
plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 150,
})

# Distinct colours for the two methods
CLR_MSC = "#2166AC"    # blue
CLR_IMCOH = "#B2182B"  # red


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def get_nperseg(pat: str) -> int:
    fs = FS_MAP.get(pat, 2048.0)
    return int(fs * 2.0)


def load_channel_labels(patient: str) -> list[str]:
    """Load channel labels (monopolar), cleaned."""
    for ext, parser in [("txt", _parse_txt), ("csv", _parse_csv)]:
        fpath = DATASET_ROOT / patient / f"channel_labels.{ext}"
        if fpath.exists():
            labels = parser(fpath)
            if labels:
                return labels
    raise FileNotFoundError(f"No channel_labels for {patient} in {DATASET_ROOT}")


def _parse_txt(path: Path) -> list[str]:
    labels = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                label = line.split(",")[0] if "," in line else line
                labels.append(label.replace(" ", ""))
    return labels


def _parse_csv(path: Path) -> list[str]:
    labels = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i == 0 and "label" in line.lower():
                continue
            line = line.strip().strip('"')
            if line:
                label = line.split(",")[0] if "," in line else line
                labels.append(label.replace(" ", ""))
    return labels


def probe_from_label(label: str) -> str:
    """Extract probe name from a monopolar sEEG label: 'A1' -> 'A'."""
    label = label.replace("\u00ec", "'")
    m = re.match(r"([A-Za-z]+'?)", label)
    return m.group(1) if m else label


def load_fc_matrix(cache_dir: Path, patient: str, band: str, phase: str,
                   method: str) -> np.ndarray | None:
    """Load an FC matrix (.npy) from cache."""
    nperseg = get_nperseg(patient)
    fname = f"{band}_{phase}_{method}_sparsify-none_nperseg-{nperseg}.npy"
    path = cache_dir / patient / fname
    if not path.exists():
        return None
    return np.load(path)


def load_lrg(cache_dir: Path, patient: str, band: str, phase: str,
             fc_tag: str) -> LRGResult | None:
    """Load an LRG result from cache."""
    path = cache_dir / patient / f"{band}_{phase}_lrg_{fc_tag}.npz"
    if not path.exists():
        return None
    data = np.load(path)
    return LRGResult(
        ultrametric_matrix=data["ultrametric_matrix"],
        linkage_matrix=data["linkage_matrix"],
        entropy_tau=data["entropy_tau"],
        entropy_1_minus_S=data["entropy_1_minus_S"],
        entropy_C=data["entropy_C"],
        optimal_threshold=float(data["optimal_threshold"]),
        patient=str(data["patient"]),
        phase=str(data["phase"]),
        band=str(data["band"]),
        fc_method=str(data["fc_method"]),
        n_nodes=int(data["n_nodes"]),
    )


def save_fig(fig, name: str, dpi: int = 300):
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"  Saved: {path}")


def compute_same_probe_ratio(fc_matrix: np.ndarray, labels: list[str]) -> float:
    """Fraction of total edge weight that falls within same-probe pairs.

    Returns the ratio:  sum(w_ij for same-probe pairs) / sum(w_ij for all pairs).
    """
    N = fc_matrix.shape[0]
    probes = [probe_from_label(l) for l in labels[:N]]
    same_probe_weight = 0.0
    total_weight = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            w = abs(fc_matrix[i, j])
            total_weight += w
            if probes[i] == probes[j]:
                same_probe_weight += w
    return same_probe_weight / total_weight if total_weight > 0 else 0.0


def compute_laplacian_eigenvalues(fc_matrix: np.ndarray) -> np.ndarray:
    """Compute sorted eigenvalues of the graph Laplacian from an adjacency matrix."""
    # Ensure non-negative weights
    A = np.abs(fc_matrix.copy())
    np.fill_diagonal(A, 0.0)
    D = np.diag(A.sum(axis=1))
    L = D - A
    eigvals = np.linalg.eigvalsh(L)
    return np.sort(eigvals)


# ---------------------------------------------------------------------------
# Figure 1: Adjacency matrix heatmaps
# ---------------------------------------------------------------------------
def fig1_adjacency_heatmaps(patient="Pat_02", phase="rest_pre", band="alpha",
                            verbose=False):
    """Side-by-side MSC vs ImCoh adjacency matrix heatmaps."""
    if verbose:
        print(f"[Fig 1] Adjacency heatmaps: {patient} {phase} {band}")

    msc = load_fc_matrix(MSC_CACHE, patient, band, phase, "msc")
    imcoh = load_fc_matrix(IMCOH_CACHE, patient, band, phase, "imcoh")
    if msc is None or imcoh is None:
        print(f"  SKIP: missing data for {patient} {phase} {band}")
        return

    labels = load_channel_labels(patient)
    N = msc.shape[0]
    probes = [probe_from_label(l) for l in labels[:N]]
    unique_probes = sorted(set(probes))
    # Sort by probe for visualization
    sort_idx = np.argsort([unique_probes.index(p) * 1000 + i for i, p in enumerate(probes)])
    msc_sorted = msc[np.ix_(sort_idx, sort_idx)]
    imcoh_sorted = imcoh[np.ix_(sort_idx, sort_idx)]

    # Probe boundary positions
    sorted_probes = [probes[i] for i in sort_idx]
    boundaries = []
    for i in range(1, N):
        if sorted_probes[i] != sorted_probes[i - 1]:
            boundaries.append(i)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    for ax, mat, title, vmax_val in zip(
        axes,
        [msc_sorted, imcoh_sorted],
        ["MSC (Magnitude-Squared Coherence)", "ImCoh (Imaginary Coherence)"],
        [1.0, None],
    ):
        if vmax_val is None:
            vmax_val = np.percentile(mat[mat > 0], 99) if np.any(mat > 0) else 0.3
        im = ax.imshow(mat, cmap="inferno", vmin=0, vmax=vmax_val,
                        aspect="equal", interpolation="none")
        # Draw probe boundaries
        for b in boundaries:
            ax.axhline(b - 0.5, color="white", lw=0.5, alpha=0.6)
            ax.axvline(b - 0.5, color="white", lw=0.5, alpha=0.6)
        ax.set_title(title, fontsize=13)
        ax.set_xlabel("Channel (sorted by probe)")
        ax.set_ylabel("Channel (sorted by probe)")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Coherence")

    # Compute same-probe ratios for annotation
    spr_msc = compute_same_probe_ratio(msc, labels)
    spr_imcoh = compute_same_probe_ratio(imcoh, labels)
    axes[0].text(0.02, 0.98, f"Same-probe ratio: {spr_msc:.3f}",
                 transform=axes[0].transAxes, va="top", fontsize=10,
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
    axes[1].text(0.02, 0.98, f"Same-probe ratio: {spr_imcoh:.3f}",
                 transform=axes[1].transAxes, va="top", fontsize=10,
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

    fig.suptitle(
        f"FC Matrix Comparison  --  {patient}  {phase}  "
        f"{BRAIN_BAND_TEX_DICT[band]}",
        fontsize=15, y=1.02,
    )
    fig.tight_layout()
    save_fig(fig, "fig1_adjacency_heatmaps.png")


# ---------------------------------------------------------------------------
# Figure 2: Entropy curves
# ---------------------------------------------------------------------------
def fig2_entropy_curves(patient="Pat_02", phase="rest_pre", band="alpha",
                        verbose=False):
    """Overlay S(tau) and C(tau) from LRG for MSC vs ImCoh."""
    if verbose:
        print(f"[Fig 2] Entropy curves: {patient} {phase} {band}")

    msc_lrg = load_lrg(MSC_LRG_CACHE, patient, band, phase, "msc")
    imcoh_lrg = load_lrg(IMCOH_LRG_CACHE, patient, band, phase, "imcoh-abs")
    if msc_lrg is None or imcoh_lrg is None:
        print(f"  SKIP: missing LRG data for {patient} {phase} {band}")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: 1-S(tau)
    ax = axes[0]
    ax.plot(msc_lrg.entropy_tau, msc_lrg.entropy_1_minus_S,
            color=CLR_MSC, lw=2, label="MSC")
    ax.plot(imcoh_lrg.entropy_tau, imcoh_lrg.entropy_1_minus_S,
            color=CLR_IMCOH, lw=2, label="ImCoh")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\tau$ (Scale Parameter)")
    ax.set_ylabel(r"$1 - S(\tau)$ (Normalized Entropy)")
    ax.set_title(r"$1 - S(\tau)$  Entropy Curve")
    ax.set_ylim(-0.02, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel B: C(tau)
    ax = axes[1]
    ax.plot(msc_lrg.entropy_tau[:-1], msc_lrg.entropy_C,
            color=CLR_MSC, lw=2, label="MSC")
    ax.plot(imcoh_lrg.entropy_tau[:-1], imcoh_lrg.entropy_C,
            color=CLR_IMCOH, lw=2, label="ImCoh")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\tau$ (Scale Parameter)")
    ax.set_ylabel(r"$C(\tau)$ (Spectral Complexity)")
    ax.set_title(r"$C(\tau)$  Complexity Curve")
    ax.set_ylim(-0.02, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Mark optimal thresholds
    for ax_i in axes:
        ax_i.axvline(msc_lrg.optimal_threshold, color=CLR_MSC,
                     ls="--", lw=1, alpha=0.6)
        ax_i.axvline(imcoh_lrg.optimal_threshold, color=CLR_IMCOH,
                     ls="--", lw=1, alpha=0.6)

    fig.suptitle(
        f"LRG Entropy Curves  --  {patient}  {phase}  "
        f"{BRAIN_BAND_TEX_DICT[band]}",
        fontsize=15, y=1.02,
    )
    fig.tight_layout()
    save_fig(fig, "fig2_entropy_curves.png")


# ---------------------------------------------------------------------------
# Figure 3: Eigenvalue spectra
# ---------------------------------------------------------------------------
def fig3_eigenvalue_spectra(patient="Pat_02", phase="rest_pre", band="alpha",
                            verbose=False):
    """Compare Laplacian eigenvalue spectra of MSC vs ImCoh."""
    if verbose:
        print(f"[Fig 3] Eigenvalue spectra: {patient} {phase} {band}")

    msc = load_fc_matrix(MSC_CACHE, patient, band, phase, "msc")
    imcoh = load_fc_matrix(IMCOH_CACHE, patient, band, phase, "imcoh")
    if msc is None or imcoh is None:
        print(f"  SKIP: missing data for {patient} {phase} {band}")
        return

    eig_msc = compute_laplacian_eigenvalues(msc)
    eig_imcoh = compute_laplacian_eigenvalues(imcoh)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: Full spectrum
    ax = axes[0]
    ax.plot(range(len(eig_msc)), eig_msc, "o-", color=CLR_MSC,
            ms=3, lw=1, label="MSC", alpha=0.8)
    ax.plot(range(len(eig_imcoh)), eig_imcoh, "s-", color=CLR_IMCOH,
            ms=3, lw=1, label="ImCoh", alpha=0.8)
    ax.set_xlabel("Eigenvalue Index")
    ax.set_ylabel(r"$\lambda_k$")
    ax.set_title("Full Laplacian Spectrum")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel B: First 30 eigenvalues (zoom)
    n_zoom = min(30, len(eig_msc), len(eig_imcoh))
    ax = axes[1]
    ax.plot(range(n_zoom), eig_msc[:n_zoom], "o-", color=CLR_MSC,
            ms=4, lw=1.5, label="MSC", alpha=0.8)
    ax.plot(range(n_zoom), eig_imcoh[:n_zoom], "s-", color=CLR_IMCOH,
            ms=4, lw=1.5, label="ImCoh", alpha=0.8)
    ax.set_xlabel("Eigenvalue Index")
    ax.set_ylabel(r"$\lambda_k$")
    ax.set_title(f"First {n_zoom} Eigenvalues (zoom)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"Laplacian Eigenvalue Spectrum  --  {patient}  {phase}  "
        f"{BRAIN_BAND_TEX_DICT[band]}",
        fontsize=15, y=1.02,
    )
    fig.tight_layout()
    save_fig(fig, "fig3_eigenvalue_spectra.png")


# ---------------------------------------------------------------------------
# Figure 4: Dendrogram comparison
# ---------------------------------------------------------------------------
def fig4_dendrogram_comparison(patient="Pat_02", phase="rest_pre", band="alpha",
                               verbose=False):
    """Side-by-side MSC vs ImCoh dendrograms."""
    if verbose:
        print(f"[Fig 4] Dendrogram comparison: {patient} {phase} {band}")

    msc_lrg = load_lrg(MSC_LRG_CACHE, patient, band, phase, "msc")
    imcoh_lrg = load_lrg(IMCOH_LRG_CACHE, patient, band, phase, "imcoh-abs")
    if msc_lrg is None or imcoh_lrg is None:
        print(f"  SKIP: missing LRG data for {patient} {phase} {band}")
        return

    labels = load_channel_labels(patient)

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))

    for ax, lrg, fc_label, color in zip(
        axes,
        [msc_lrg, imcoh_lrg],
        ["MSC", "ImCoh"],
        [CLR_MSC, CLR_IMCOH],
    ):
        linkage = lrg.linkage_matrix.copy()
        condensed = lrg.ultrametric_matrix
        # Optimal leaf ordering
        try:
            linkage = optimal_leaf_ordering(linkage, condensed)
        except Exception:
            pass

        n_nodes = lrg.n_nodes
        use_labels = labels[:n_nodes] if len(labels) >= n_nodes else None

        dendrogram(
            linkage,
            ax=ax,
            orientation="top",
            labels=use_labels,
            no_labels=(n_nodes > 50),
            color_threshold=lrg.optimal_threshold,
            above_threshold_color="grey",
            leaf_font_size=5,
        )

        # Log-scale y axis with proper limits
        merge_heights = linkage[:, 2]
        tmin = merge_heights[merge_heights > 0].min() * 0.5
        tmax = merge_heights.max() * 2.0
        ax.set_yscale("log")
        ax.set_ylim(tmin, tmax)

        ax.axhline(lrg.optimal_threshold, color=color, ls="--", lw=2,
                    label=f"Optimal threshold = {lrg.optimal_threshold:.4f}")
        ax.set_title(f"{fc_label}", fontsize=14)
        ax.set_ylabel("Ultrametric Distance")
        ax.set_xlabel("Channel")
        ax.legend(fontsize=9)

    fig.suptitle(
        f"LRG Dendrogram Comparison  --  {patient}  {phase}  "
        f"{BRAIN_BAND_TEX_DICT[band]}",
        fontsize=15, y=1.02,
    )
    fig.tight_layout()
    save_fig(fig, "fig4_dendrogram_comparison.png")


# ---------------------------------------------------------------------------
# Figure 5: Same-probe ratio across all patients and bands
# ---------------------------------------------------------------------------
def fig5_same_probe_ratio(phase="rest_pre", verbose=False):
    """Bar chart: same-probe weight ratio for MSC vs ImCoh, all patients x bands."""
    if verbose:
        print(f"[Fig 5] Same-probe ratio summary (phase={phase})")

    n_patients = len(PATIENTS)
    n_bands = len(BANDS)
    msc_ratios = np.full((n_patients, n_bands), np.nan)
    imcoh_ratios = np.full((n_patients, n_bands), np.nan)

    for pi, pat in enumerate(PATIENTS):
        labels = load_channel_labels(pat)
        for bi, band in enumerate(BANDS):
            msc = load_fc_matrix(MSC_CACHE, pat, band, phase, "msc")
            imcoh = load_fc_matrix(IMCOH_CACHE, pat, band, phase, "imcoh")
            if msc is not None:
                msc_ratios[pi, bi] = compute_same_probe_ratio(msc, labels)
                if verbose:
                    print(f"    {pat} {band} MSC:   {msc_ratios[pi, bi]:.4f}")
            if imcoh is not None:
                imcoh_ratios[pi, bi] = compute_same_probe_ratio(imcoh, labels)
                if verbose:
                    print(f"    {pat} {band} ImCoh: {imcoh_ratios[pi, bi]:.4f}")

    # Create grouped bar chart: each patient group has 6 pairs of bars
    fig, ax = plt.subplots(figsize=(18, 6))

    x_positions = []
    msc_vals = []
    imcoh_vals = []
    tick_labels = []
    tick_positions = []

    group_width = n_bands * 2.5 + 2  # spacing between patient groups
    bar_width = 0.45

    for pi, pat in enumerate(PATIENTS):
        group_start = pi * group_width
        pat_label = pat.replace("Pat_0", "P")
        tick_positions.append(group_start + (n_bands - 0.5) * 1.1)

        for bi, band in enumerate(BANDS):
            x = group_start + bi * 1.1
            x_positions.append(x)
            msc_vals.append(msc_ratios[pi, bi])
            imcoh_vals.append(imcoh_ratios[pi, bi])
            if pi == 0:  # only for first patient to use as legend
                tick_labels.append(BAND_TEX[bi])

    x_arr = np.array(x_positions)
    msc_arr = np.array(msc_vals)
    imcoh_arr = np.array(imcoh_vals)

    bars_msc = ax.bar(x_arr - bar_width / 2, msc_arr, bar_width,
                       color=CLR_MSC, alpha=0.85, label="MSC", edgecolor="white",
                       linewidth=0.5)
    bars_imcoh = ax.bar(x_arr + bar_width / 2, imcoh_arr, bar_width,
                         color=CLR_IMCOH, alpha=0.85, label="ImCoh",
                         edgecolor="white", linewidth=0.5)

    # Band labels below each group
    for pi, pat in enumerate(PATIENTS):
        group_start = pi * group_width
        for bi in range(n_bands):
            x = group_start + bi * 1.1
            ax.text(x, -0.015, BAND_TEX[bi], ha="center", va="top",
                    fontsize=8, transform=ax.get_xaxis_transform())

    # Patient labels
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(PATIENTS, fontsize=11)
    ax.tick_params(axis="x", length=0, pad=20)

    ax.set_ylabel("Same-Probe Weight Ratio")
    ax.set_title(
        f"Same-Probe Edge Weight Ratio  --  MSC vs ImCoh  ({phase})",
        fontsize=14,
    )
    ax.legend(loc="upper right", fontsize=11)
    ax.grid(axis="y", alpha=0.3)

    # Mark Pat_03 as outlier
    pat03_idx = PATIENTS.index("Pat_03")
    pat03_start = pat03_idx * group_width
    pat03_end = pat03_start + (n_bands - 1) * 1.1 + 1
    ax.axvspan(pat03_start - 0.8, pat03_end, alpha=0.08, color="orange",
               label="Pat_03 (outlier)")
    ax.text(pat03_start + (n_bands - 1) * 1.1 / 2, ax.get_ylim()[1] * 0.95,
            "outlier\n(1024 Hz)", ha="center", fontsize=8, color="darkorange",
            fontstyle="italic")

    fig.tight_layout()
    save_fig(fig, "fig5_same_probe_ratio.png")


# ---------------------------------------------------------------------------
# Figure 6: Community structure at n=10
# ---------------------------------------------------------------------------
def fig6_community_structure(patient="Pat_02", phase="rest_pre", band="alpha",
                             n_communities=10, verbose=False):
    """Show community membership for MSC vs ImCoh at a fixed number of communities."""
    if verbose:
        print(f"[Fig 6] Community structure n={n_communities}: "
              f"{patient} {phase} {band}")

    msc_lrg = load_lrg(MSC_LRG_CACHE, patient, band, phase, "msc")
    imcoh_lrg = load_lrg(IMCOH_LRG_CACHE, patient, band, phase, "imcoh-abs")
    if msc_lrg is None or imcoh_lrg is None:
        print(f"  SKIP: missing LRG data for {patient} {phase} {band}")
        return

    labels = load_channel_labels(patient)
    N = msc_lrg.n_nodes

    # Cut dendrograms at n communities
    msc_comm = fcluster(msc_lrg.linkage_matrix, t=n_communities,
                        criterion="maxclust")
    imcoh_comm = fcluster(imcoh_lrg.linkage_matrix, t=n_communities,
                          criterion="maxclust")

    # Build probe information
    probes = [probe_from_label(l) for l in labels[:N]]
    unique_probes = sorted(set(probes))

    # Sort nodes by probe then by index within probe
    sort_idx = np.argsort([unique_probes.index(p) * 1000 + i
                           for i, p in enumerate(probes)])

    fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)

    cmap = matplotlib.colormaps.get_cmap("tab20").resampled(n_communities)

    for ax, comm, title in zip(
        axes,
        [msc_comm, imcoh_comm],
        [f"MSC Communities (n={n_communities})",
         f"ImCoh Communities (n={n_communities})"],
    ):
        sorted_comm = comm[sort_idx]
        sorted_labels = [labels[i] for i in sort_idx]
        sorted_probes = [probes[i] for i in sort_idx]

        # Colour each node by community
        colors = [cmap(c - 1) for c in sorted_comm]
        ax.bar(range(N), np.ones(N), color=colors, width=1.0, edgecolor="none")

        # Probe boundaries
        for i in range(1, N):
            if sorted_probes[i] != sorted_probes[i - 1]:
                ax.axvline(i - 0.5, color="black", lw=0.8, alpha=0.6)

        # Probe labels at midpoints
        prev = 0
        for i in range(1, N):
            if sorted_probes[i] != sorted_probes[i - 1]:
                mid = (prev + i) / 2
                ax.text(mid, 0.5, sorted_probes[prev], ha="center", va="center",
                        fontsize=7, fontweight="bold", rotation=90)
                prev = i
        mid = (prev + N) / 2
        ax.text(mid, 0.5, sorted_probes[prev], ha="center", va="center",
                fontsize=7, fontweight="bold", rotation=90)

        ax.set_ylabel(title)
        ax.set_ylim(0, 1)
        ax.set_yticks([])

    axes[1].set_xlabel("Channel (sorted by probe)")

    # Compute overlap: how many nodes agree on community assignment
    # Use Adjusted Rand Index
    from sklearn.metrics import adjusted_rand_score
    ari = adjusted_rand_score(msc_comm, imcoh_comm)

    # Add colourbar / legend for communities
    sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=0.5, vmax=n_communities + 0.5))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation="vertical", fraction=0.02, pad=0.02)
    cbar.set_label("Community ID")
    cbar.set_ticks(range(1, n_communities + 1))

    fig.suptitle(
        f"Community Structure (n={n_communities})  --  {patient}  {phase}  "
        f"{BRAIN_BAND_TEX_DICT[band]}\n"
        f"Adjusted Rand Index = {ari:.3f}",
        fontsize=14, y=1.04,
    )
    fig.subplots_adjust(right=0.92)
    save_fig(fig, "fig6_community_structure.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate MSC vs ImCoh visualization gallery."
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--patient", default="Pat_02",
                        help="Representative patient for single-case figures")
    parser.add_argument("--phase", default="rest_pre",
                        help="Phase for representative figures")
    parser.add_argument("--band", default="alpha",
                        help="Band for representative figures")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {OUT}")
    print()

    # 1. Adjacency heatmaps
    fig1_adjacency_heatmaps(args.patient, args.phase, args.band, args.verbose)

    # 2. Entropy curves
    fig2_entropy_curves(args.patient, args.phase, args.band, args.verbose)

    # 3. Eigenvalue spectra
    fig3_eigenvalue_spectra(args.patient, args.phase, args.band, args.verbose)

    # 4. Dendrogram comparison
    fig4_dendrogram_comparison(args.patient, args.phase, args.band, args.verbose)

    # 5. Same-probe ratio (all patients x bands)
    fig5_same_probe_ratio(args.phase, args.verbose)

    # 6. Community structure
    fig6_community_structure(args.patient, args.phase, args.band, verbose=args.verbose)

    print()
    print("Gallery complete.")


if __name__ == "__main__":
    main()
