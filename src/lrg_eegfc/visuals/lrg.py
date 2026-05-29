"""LRG (Laplacian Renormalization Group) visualization functions.

This module provides visualization functions for LRG analysis results including
entropy curves, dendrograms, ultrametric distances, and comprehensive panels.
"""

from pathlib import Path
from typing import List, Optional

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.io
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.cluster.hierarchy import dendrogram, fcluster, optimal_leaf_ordering
from scipy.spatial.distance import squareform

from lrg_eegfc.config.paths import CORR_CACHE, LRG_CACHE, MSC_CACHE, SEEG_DATAPATH, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.fc import load_fc_matrix as _load_fc_matrix

__all__ = [
    "compute_partition_stability_index",
    "plot_lrg_entropy_curves",
    "plot_lrg_dendrogram",
    "plot_lrg_dendrogram_shaft_colored",
    "plot_ultrametric_heatmap",
]
# plot_lrg_full_panel moved to lrg_panels.py on 2026-05-29 (Phase 4-B split 4/7).


def _load_channel_labels(patient: str, dataset_root: Path) -> List[str]:
    """Helper function to load channel labels.

    Tries multiple sources in order:
    1. channel_labels.txt (preferred - plain text)
    2. channel_labels.csv
    3. channel_labels.mat

    Labels are simplified by removing reference suffix (e.g., ",G2").
    """
    # Try plain text file first (most reliable)
    txt_file = dataset_root / patient / "channel_labels.txt"
    if txt_file.exists():
        try:
            with open(txt_file, "r") as f:
                labels = []
                for line in f:
                    line = line.strip()
                    if line:
                        # Remove reference suffix (e.g., ",G2") and clean up
                        if "," in line:
                            label = line.split(",")[0]
                        else:
                            label = line
                        # Remove spaces for compact dendrogram labels
                        label = label.replace(" ", "")
                        labels.append(label)
                if labels:
                    return labels
        except Exception:
            pass

    # Try CSV file
    csv_file = dataset_root / patient / "channel_labels.csv"
    if csv_file.exists():
        try:
            with open(csv_file, "r") as f:
                labels = []
                for i, line in enumerate(f):
                    if i == 0 and "label" in line.lower():
                        continue  # Skip header
                    line = line.strip().strip('"')
                    if line:
                        if "," in line:
                            label = line.split(",")[0]
                        else:
                            label = line
                        label = label.replace(" ", "")
                        labels.append(label)
                if labels:
                    return labels
        except Exception:
            pass

    # Fall back to .mat file
    try:
        label_file = dataset_root / patient / "channel_labels.mat"
        if label_file.exists():
            label_data = scipy.io.loadmat(label_file)
            # Try different possible field names
            if "channel_labels" in label_data:
                labels_raw = label_data["channel_labels"].flatten()
            elif "ChannelNames" in label_data:
                labels_raw = label_data["ChannelNames"].flatten()
            else:
                # Use first non-metadata field
                for key in label_data:
                    if not key.startswith("__"):
                        labels_raw = label_data[key].flatten()
                        break
            return [str(label[0]) if hasattr(label, '__getitem__') else str(label) for label in labels_raw]
    except Exception:
        pass
    return []


def _find_psi_optimal_partition(psi_values: np.ndarray, n_communities: np.ndarray) -> int:
    """Find optimal number of communities from PSI using peak detection.

    Strategy:
    1. Find all local maxima (peaks) in PSI
    2. If 2nd peak > max_peak / 2, use 2nd peak (hierarchical structure)
    3. Otherwise use the global maximum

    Parameters
    ----------
    psi_values : np.ndarray
        PSI values
    n_communities : np.ndarray
        Corresponding number of communities

    Returns
    -------
    int
        Optimal number of communities
    """
    if len(psi_values) < 3:
        return int(n_communities[np.argmax(psi_values)]) if len(psi_values) > 0 else 2

    # Find local maxima (peaks)
    # A peak is a point higher than both neighbors
    peaks_idx = []
    for i in range(1, len(psi_values) - 1):
        if psi_values[i] > psi_values[i - 1] and psi_values[i] > psi_values[i + 1]:
            peaks_idx.append(i)

    # Also check endpoints
    if psi_values[0] > psi_values[1]:
        peaks_idx.insert(0, 0)
    if psi_values[-1] > psi_values[-2]:
        peaks_idx.append(len(psi_values) - 1)

    if len(peaks_idx) == 0:
        # No peaks found, use global max
        return int(n_communities[np.argmax(psi_values)])

    # Sort peaks by PSI value (descending)
    peaks_idx = sorted(peaks_idx, key=lambda i: psi_values[i], reverse=True)

    max_psi = psi_values[peaks_idx[0]]

    # If there's a second peak and it's > max/2, use it (indicates hierarchical structure)
    if len(peaks_idx) > 1:
        second_psi = psi_values[peaks_idx[1]]
        if second_psi > max_psi / 2:
            # Sort the top 2 peaks by n_communities (prefer larger partition)
            top_two = sorted(peaks_idx[:2], key=lambda i: n_communities[i], reverse=True)
            return int(n_communities[top_two[0]])

    # Use global maximum
    return int(n_communities[peaks_idx[0]])


def compute_partition_stability_index(linkage_matrix: np.ndarray):
    """Compute Partition Stability Index from hierarchical clustering.

    Implements Eq. (2) from the multiscale community-detection paper:
        Ψ(n; τ) = N_norm * (log10 Δ_n - log10 Δ_{n+1})
    where Δ_n are dendrogram branch thresholds ordered from the root (single
    cluster) downward and
        N_norm = [log10 Δ_1 - log10 Δ_{n_max}]^{-1}
    normalizes by the total dendrogram length in log-space.

    Parameters
    ----------
    linkage_matrix : np.ndarray
        Linkage matrix from hierarchical clustering (N-1, 4).

    Returns
    -------
    psi_values : np.ndarray
        Partition Stability Index values using the paper's normalization.
    n_communities : np.ndarray
        Number of communities corresponding to each PSI value.
    """
    # Merge heights are non-decreasing; reverse so Δ_n starts at the root (largest gap)
    deltas = linkage_matrix[:, 2][::-1]

    # Guard against zeros/negatives before taking logs
    eps = 1e-12
    deltas = np.clip(deltas, eps, None)
    log_deltas = np.log10(deltas)

    if len(log_deltas) < 2:
        return np.array([]), np.array([])

    # Normalization constant over total dendrogram length (Eq. 2)
    denom = log_deltas[0] - log_deltas[-1]
    norm = 1.0 / denom if denom > 0 else 0.0

    # Ψ for consecutive gaps; n_communities starts at 2 (first split) up to N-1
    psi_values = norm * (log_deltas[:-1] - log_deltas[1:])
    n_communities = np.arange(2, 2 + len(psi_values))

    return psi_values, n_communities


def plot_lrg_entropy_curves(
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    cache_root: Path = LRG_CACHE,
    output_path: Optional[Path] = None,
    figsize: tuple = (10, 6),
) -> Path:
    """Plot LRG entropy curves (1-S and C) vs tau.

    This visualization shows how entropy metrics evolve across scales,
    revealing hierarchical transitions in the network.

    Parameters
    ----------
    patient : str
        Patient ID
    phase : str
        Experimental phase
    band : str
        Frequency band
    fc_method : str
        FC method ("corr" or "msc")
    cache_root : Path, optional
        Root directory for LRG cache
    output_path : Optional[Path], optional
        Path to save figure (auto-generated if None)
    figsize : tuple, optional
        Figure size

    Returns
    -------
    Path
        Path to the saved figure
    """
    # Load LRG result
    result = load_lrg_result(patient, phase, band, fc_method, cache_root)

    if result is None:
        from lrg_eegfc.workflow.lrg import get_lrg_cache_path
        cache_path = get_lrg_cache_path(patient, phase, band, fc_method, cache_root)
        raise FileNotFoundError(f"LRG result not found: {cache_path}")

    # Extract entropy data
    tau = result.entropy_tau
    S_norm = result.entropy_1_minus_S  # 1-S (normalized entropy)
    C = result.entropy_C  # Spectral complexity

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xscale("log")

    # Plot both curves
    ax.plot(tau, S_norm, label=r"$1-S$ (Normalized Entropy)", color="blue", lw=2)
    ax.plot(tau[:-1], C, label=r"$C$ (Spectral Complexity)", color="red", lw=2)

    # Formatting
    ax.set_xlabel(r"$\tau$ (Scale Parameter)", fontsize=12)
    ax.set_ylabel("Entropy Metrics", fontsize=12)
    ax.set_title(f"LRG Entropy - {patient} {phase} {band} ({fc_method})", fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)  # Both metrics normalized to [0, 1]

    # Generate output path if not provided
    if output_path is None:
        output_dir = FIGURES_ROOT / "lrg" / patient
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{band}_{phase}_lrg_{fc_method}_entropy.png"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save figure
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return output_path


def plot_lrg_dendrogram(
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    cache_root: Path = LRG_CACHE,
    output_path: Optional[Path] = None,
    orientation: str = "top",
    figsize: tuple = (12, 8),
    dataset_root: Path = SEEG_DATAPATH,
    optimal_leaf_order: bool = True,
    show_labels: bool = True,
    max_labels: int = 50,
) -> Path:
    """Plot hierarchical dendrogram with optimal threshold line.

    Parameters
    ----------
    patient : str
        Patient ID
    phase : str
        Experimental phase
    band : str
        Frequency band
    fc_method : str
        FC method ("corr" or "msc")
    cache_root : Path, optional
        Root directory for LRG cache
    output_path : Optional[Path], optional
        Path to save figure (auto-generated if None)
    orientation : str, optional
        Dendrogram orientation ("top", "right", "bottom", "left")
    figsize : tuple, optional
        Figure size
    dataset_root : Path, optional
        Root directory for patient data
    optimal_leaf_order : bool, optional
        Use optimal leaf ordering for better visualization
    show_labels : bool, optional
        Show channel labels
    max_labels : int, optional
        Maximum number of labels to show

    Returns
    -------
    Path
        Path to the saved figure
    """
    # Load LRG result
    result = load_lrg_result(patient, phase, band, fc_method, cache_root)

    if result is None:
        from lrg_eegfc.workflow.lrg import get_lrg_cache_path
        cache_path = get_lrg_cache_path(patient, phase, band, fc_method, cache_root)
        raise FileNotFoundError(f"LRG result not found: {cache_path}")

    linkage = result.linkage_matrix
    ultrametric_condensed = result.ultrametric_matrix
    optimal_th = result.optimal_threshold

    # Optional: Optimal leaf ordering for better visualization
    if optimal_leaf_order:
        try:
            linkage = optimal_leaf_ordering(linkage, ultrametric_condensed)
        except Exception:
            pass  # Fall back to default ordering if optimization fails

    # Load channel labels
    labels = _load_channel_labels(patient, dataset_root)
    if not labels:
        labels = [str(i) for i in range(result.n_nodes)]

    # Create dendrogram
    fig, ax = plt.subplots(figsize=figsize)

    dendro = dendrogram(
        linkage,
        ax=ax,
        orientation=orientation,
        labels=labels if (show_labels and len(labels) <= max_labels) else None,
        no_labels=(not show_labels or len(labels) > max_labels),
        color_threshold=optimal_th,
        above_threshold_color="black",  # Clusters above threshold are black
        leaf_font_size=8 if len(labels) <= max_labels else 5,
    )

    # Compute proper axis limits from merge heights (log-space padding)
    merge_heights = linkage[:, 2]
    tmin = merge_heights[merge_heights > 0].min() * 0.5
    tmax = merge_heights.max() * 2.0

    # Add optimal threshold line and set axis properties
    if orientation in ("top", "bottom"):
        ax.axhline(
            optimal_th,
            color="blue",
            linestyle="--",
            lw=2,
            label=f"Optimal Threshold = {optimal_th:.3f}",
        )
        ax.set_ylabel("Ultrametric Distance", fontsize=12)
        ax.set_xlabel("Channel Index", fontsize=12)
        ax.set_yscale("log")
        ax.set_ylim(tmin, tmax)
    elif orientation in ("right", "left"):
        ax.axvline(
            optimal_th,
            color="blue",
            linestyle="--",
            lw=2,
            label=f"Optimal Threshold = {optimal_th:.3f}",
        )
        ax.set_xlabel("Ultrametric Distance", fontsize=12)
        ax.set_ylabel("Channel Index", fontsize=12)
        ax.set_xscale("log")
        ax.set_xlim(tmin, tmax)

    ax.set_title(
        f"LRG Dendrogram - {patient} {phase} {band} ({fc_method})", fontsize=14
    )
    ax.legend(fontsize=10)

    # Generate output path if not provided
    if output_path is None:
        output_dir = FIGURES_ROOT / "lrg" / patient
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{band}_{phase}_lrg_{fc_method}_dendrogram.png"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save figure
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return output_path


def _psi_threshold_from_linkage(linkage_matrix, psi_optimal_n, n_nodes,
                                 fallback_threshold):
    """Compute the dendrogram cut height for a given n* (mirrors full_panel)."""
    merge_heights = linkage_matrix[:, 2]
    if 2 <= psi_optimal_n <= n_nodes:
        cut_idx = n_nodes - psi_optimal_n
        if 0 < cut_idx < len(merge_heights):
            return float((merge_heights[cut_idx - 1] + merge_heights[cut_idx]) / 2)
        if cut_idx == 0:
            return float(merge_heights[0] / 2)
    return float(fallback_threshold)


def plot_lrg_dendrogram_shaft_colored(
    ax,
    lrg_result,
    channel_labels: List[str],
    *,
    show_psi_cut: bool = True,
    show_xlabels: bool = True,
    leaf_font_size: float = 4.0,
    title: Optional[str] = None,
    branch_palette: str = "tab10",
    shaft_palette: str = "tab20",
    above_threshold_color: str = "0.55",
    optimal_leaf_order_flag: bool = False,
    ylim: Optional[tuple] = None,
    max_n_communities: Optional[int] = None,
    min_n_communities: Optional[int] = None,
    yscale: str = "log",
):
    """Render an LRG dendrogram on *ax* with two encodings:

    * Leaf tick-labels are coloured by electrode shaft identity
      (``shaft_palette``).
    * Branches at-and-below the Ψ-optimal cut are coloured by community
      using ``branch_palette``; branches above the cut are drawn in
      ``above_threshold_color``.

    The Y-axis is log-scaled in raw merge-height units; callers may impose
    a shared limit per column via ``ylim`` so that figures are visually
    comparable across fc_methods.

    Parameters
    ----------
    ax : matplotlib Axes
    lrg_result : LRGResult
        Result with ``linkage_matrix`` and ``optimal_threshold``.
    channel_labels : list of str
        Labels in the giant-component order used by the linkage matrix.
        Length must equal ``lrg_result.n_nodes``.
    show_psi_cut : bool
        Draw a horizontal dashed line at the Ψ-optimal threshold.
    show_xlabels : bool
        Whether to render leaf tick labels.
    leaf_font_size : float
        Font size of leaf tick labels.
    title : str, optional
        Panel title.
    branch_palette : str
        Matplotlib colormap name used to colour communities.
    shaft_palette : str
        Matplotlib colormap name used to colour leaves by shaft.
    above_threshold_color : str
        Colour for branches above the Ψ cut.
    optimal_leaf_order_flag : bool
        Apply ``scipy.cluster.hierarchy.optimal_leaf_ordering``.
    ylim : (low, high), optional
        Y-axis limits (raw merge heights). If None, derived from the
        linkage matrix with log padding.

    Returns
    -------
    dict with keys: ``psi_n``, ``psi_threshold``, ``leaves``, ``dendro``.
    """
    from scipy.cluster import hierarchy
    from lrg_eegfc.utils.probe import extract_probe_labels

    Z = lrg_result.linkage_matrix
    n_nodes = lrg_result.n_nodes

    if len(channel_labels) != n_nodes:
        if len(channel_labels) < n_nodes:
            raise ValueError(
                f"channel_labels has {len(channel_labels)} entries but "
                f"linkage matrix expects {n_nodes} leaves."
            )
        channel_labels = list(channel_labels[:n_nodes])
    else:
        channel_labels = list(channel_labels)

    if optimal_leaf_order_flag:
        try:
            Z = optimal_leaf_ordering(Z, lrg_result.ultrametric_matrix)
        except Exception:
            pass

    psi_values, psi_n_communities = compute_partition_stability_index(Z)
    if psi_values.size > 0:
        mask = np.ones_like(psi_n_communities, dtype=bool)
        if max_n_communities is not None:
            mask &= psi_n_communities <= max_n_communities
        if min_n_communities is not None:
            mask &= psi_n_communities >= min_n_communities
        if mask.any():
            psi_n = _find_psi_optimal_partition(
                psi_values[mask], psi_n_communities[mask]
            )
        else:
            psi_n = _find_psi_optimal_partition(psi_values, psi_n_communities)
    else:
        psi_n = max(2, int(np.unique(channel_labels).size // 2))
    psi_threshold = _psi_threshold_from_linkage(
        Z, psi_n, n_nodes, lrg_result.optimal_threshold,
    )

    # Branch palette (cycled across communities at the cut)
    n_colors = max(min(psi_n + 2, 20), 3)
    palette = [plt.matplotlib.colors.to_hex(c)
               for c in plt.get_cmap(branch_palette)(np.linspace(0, 1, n_colors))]
    hierarchy.set_link_color_palette(palette)

    dendro = dendrogram(
        Z,
        ax=ax,
        labels=channel_labels,
        no_labels=not show_xlabels,
        color_threshold=psi_threshold,
        above_threshold_color=above_threshold_color,
        leaf_font_size=leaf_font_size,
    )

    # Recolour leaf tick labels by shaft
    if show_xlabels:
        probe_ids_full = extract_probe_labels(channel_labels)
        unique_probes = sorted(set(probe_ids_full))
        shaft_cmap = plt.get_cmap(shaft_palette, max(len(unique_probes), 1))
        shaft_color = {p: shaft_cmap(i) for i, p in enumerate(unique_probes)}
        leaf_order = dendro["leaves"]
        ordered_probes = [probe_ids_full[i] for i in leaf_order]
        for tick_label, probe in zip(ax.get_xticklabels(), ordered_probes):
            tick_label.set_color(shaft_color[probe])

    # Y axis: log scale on raw merge heights — same recipe as
    # gen_phase_reorg_figures.py / wp1_figures.py (the canonical "good"
    # dendrogram style). Tight padding: 0.8× smallest, 1.05× largest.
    merge_heights = Z[:, 2]
    pos_heights = merge_heights[merge_heights > 0]
    if ylim is None:
        if pos_heights.size:
            if yscale == "log":
                ax.set_ylim(pos_heights[0] * 0.8, merge_heights[-1] * 1.05)
            else:
                span = merge_heights[-1] - pos_heights[0]
                pad = 0.03 * span
                ax.set_ylim(max(0.0, pos_heights[0] - pad),
                            merge_heights[-1] + pad)
    else:
        ax.set_ylim(*ylim)
    ax.set_yscale(yscale)

    if show_psi_cut:
        ax.axhline(
            psi_threshold,
            color="black",
            linestyle="--",
            linewidth=1.0,
            alpha=0.7,
        )

    if title is not None:
        ax.set_title(title, fontsize=10)

    # Restore default link palette so we don't leak state
    hierarchy.set_link_color_palette(None)

    return {
        "psi_n": psi_n,
        "psi_threshold": psi_threshold,
        "leaves": dendro["leaves"],
        "dendro": dendro,
    }


def plot_ultrametric_heatmap(
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    cache_root: Path = LRG_CACHE,
    output_path: Optional[Path] = None,
    figsize: tuple = (10, 10),
    dataset_root: Path = SEEG_DATAPATH,
    cmap: str = "viridis",
) -> Path:
    """Plot ultrametric distance matrix as heatmap.

    Parameters
    ----------
    patient : str
        Patient ID
    phase : str
        Experimental phase
    band : str
        Frequency band
    fc_method : str
        FC method ("corr" or "msc")
    cache_root : Path, optional
        Root directory for LRG cache
    output_path : Optional[Path], optional
        Path to save figure (auto-generated if None)
    figsize : tuple, optional
        Figure size
    dataset_root : Path, optional
        Root directory for patient data
    cmap : str, optional
        Colormap name

    Returns
    -------
    Path
        Path to the saved figure
    """
    # Load LRG result
    result = load_lrg_result(patient, phase, band, fc_method, cache_root)

    if result is None:
        from lrg_eegfc.workflow.lrg import get_lrg_cache_path
        cache_path = get_lrg_cache_path(patient, phase, band, fc_method, cache_root)
        raise FileNotFoundError(f"LRG result not found: {cache_path}")

    ultrametric_condensed = result.ultrametric_matrix

    # Convert condensed to square form
    ultrametric_square = squareform(ultrametric_condensed)

    # Load channel labels
    labels = _load_channel_labels(patient, dataset_root)
    if not labels:
        labels = [str(i) for i in range(ultrametric_square.shape[0])]

    # Plot symmetric heatmap
    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(
        ultrametric_square, cmap=cmap, interpolation="none", aspect="auto"
    )

    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label("Ultrametric Distance", fontsize=12)

    # Add labels if not too many
    if len(labels) <= 50:
        ax.set_xticks(range(len(labels)))
        ax.set_yticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=8)
        ax.set_yticklabels(labels, fontsize=8)
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    ax.set_title(
        f"Ultrametric Distance - {patient} {phase} {band} ({fc_method})",
        fontsize=14,
    )

    # Generate output path if not provided
    if output_path is None:
        output_dir = FIGURES_ROOT / "lrg" / patient
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{band}_{phase}_lrg_{fc_method}_ultrametric.png"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save figure
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return output_path
