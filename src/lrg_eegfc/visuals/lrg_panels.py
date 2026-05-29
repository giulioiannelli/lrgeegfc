"""Composite LRG full-panel figure.

Split out of `lrg.py` on 2026-05-29 (Phase 4-B split 4/7). Holds the
5-panel composite `plot_lrg_full_panel`. Individual-panel helpers
(entropy curves, dendrogram, ultrametric heatmap) stay in `lrg.py`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster

from lrg_eegfc.config.paths import CORR_CACHE, LRG_CACHE, SEEG_DATAPATH, FIGURES_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix as _load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result

from .lrg import (
    _find_psi_optimal_partition,
    _load_channel_labels,
    compute_partition_stability_index,
)


__all__ = ["plot_lrg_full_panel"]


def plot_lrg_full_panel(
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    cache_root: Path = LRG_CACHE,
    dataset_root: Path = SEEG_DATAPATH,
    output_path: Optional[Path] = None,
    figsize: tuple = (20, 12),
    verbose: bool = False,
    n_communities: Optional[int] = None,
    fc_matrix: Optional[np.ndarray] = None,
) -> Path:
    """Create comprehensive LRG analysis visualization - CORRECTED VERSION.

    Generates a 5-panel layout following FIGMNTGN notebooks exactly:
    - Panel (a): FC Matrix Heatmap
    - Panel (b): Entropy and Specific Heat (NO N multiplication)
    - Panel (c): Dendrogram (correct tmin/tmax/log scale from FIGMNTGN03)
    - Panel (d): Network with partition colors
    - Panel (e): PSI Plot (Partition Stability Index)

    Parameters
    ----------
    patient : str
        Patient identifier
    phase : str
        Recording phase
    band : str
        Frequency band
    fc_method : str
        FC method: "corr" or "msc"
    cache_root : Path
        Root directory for LRG cache
    dataset_root : Path
        Root directory for patient data
    output_path : Path, optional
        Output file path
    figsize : tuple
        Figure size
    verbose : bool
        Print progress
    n_communities : int, optional
        Fixed number of communities to use for partitioning.
        If provided, overrides PSI-based selection for cross-phase comparison.
    fc_matrix : np.ndarray, optional
        Pre-loaded FC matrix to use for panel (a) heatmap and network.
        If provided, skips auto-loading from cache. Useful for CReMa-validated
        or other custom matrices.

    Returns
    -------
    Path
        Path to saved figure
    """
    # Load LRG result
    lrg_result = load_lrg_result(patient, phase, band, fc_method, cache_root)

    if lrg_result is None:
        raise FileNotFoundError(
            f"LRG result not found: {patient} {phase} {band} ({fc_method})\n"
            f"Run: python src/compute_lrg_analysis.py --patient {patient} --fc-method {fc_method}"
        )

    if verbose:
        print(f"Loaded LRG result: {lrg_result.n_nodes} nodes")

    # Load FC matrix for visualization
    _METHOD_LABELS = {
        "corr": "Correlation",
        "msc": r"$\mathrm{MSC}$",
        "imcoh": r"$\mathrm{ImCoh}$",
        "imcoh_abs": r"$|\mathrm{ImCoh}|$",
        "imcoh_sq": r"$|\mathrm{ImCoh}|^2$",
    }
    percolation_threshold = None  # Will be set for correlation method
    if fc_matrix is not None:
        # Use pre-loaded FC matrix (e.g. CReMa-validated)
        method_label = _METHOD_LABELS.get(fc_method, fc_method.upper())
        if verbose:
            print(f"Using pre-loaded FC matrix ({fc_matrix.shape[0]} nodes)")
    else:
        fc_matrix = _load_fc_matrix(patient, phase, band, fc_method)
        method_label = _METHOD_LABELS.get(fc_method, fc_method.upper())

        # Load percolation threshold for correlation method
        if fc_method == "corr":
            meta_path = CORR_CACHE / patient / f"{band}_{phase}_corr_cleaned_meta.npz"
            if meta_path.exists():
                meta_data = np.load(meta_path, allow_pickle=True)
                percolation_threshold = float(meta_data["threshold"])
                if verbose:
                    print(f"Loaded percolation threshold: {percolation_threshold:.6f}")
            elif verbose:
                print("Warning: No percolation metadata found, will use all edges for layout")

    if fc_matrix is None:
        raise FileNotFoundError(f"FC matrix not found for {patient} {phase} {band} ({fc_method})")

    # Extract data from LRG result
    linkage_matrix = lrg_result.linkage_matrix
    ultrametric_matrix = lrg_result.ultrametric_matrix
    optimal_threshold = lrg_result.optimal_threshold
    entropy_tau = lrg_result.entropy_tau
    entropy_1_minus_S = lrg_result.entropy_1_minus_S
    entropy_C = lrg_result.entropy_C
    n_nodes = lrg_result.n_nodes

    # If LRG used only the giant component (n_nodes < matrix size),
    # extract the same subgraph so panels (a) and (d) match the dendrogram.
    _gc_nodes = None  # original indices of giant component nodes
    if n_nodes < fc_matrix.shape[0]:
        from lrgsglib.utils import get_giant_component as _get_gc
        _G_full = nx.from_numpy_array(fc_matrix)
        _G_giant = _get_gc(_G_full)
        _gc_nodes = sorted(_G_giant.nodes())
        fc_matrix = fc_matrix[np.ix_(_gc_nodes, _gc_nodes)]
        if verbose:
            print(f"Extracted giant component: {n_nodes}/{_G_full.number_of_nodes()} nodes")

    # Compute PSI
    psi_values, psi_n_communities = compute_partition_stability_index(linkage_matrix)

    # Determine target number of communities
    if n_communities is not None:
        # Use fixed number of communities (for cross-phase comparison)
        target_n = n_communities
        if verbose:
            print(f"Using fixed n_communities={target_n}")
    else:
        # Get optimal partition using PSI peak detection
        target_n = _find_psi_optimal_partition(psi_values, psi_n_communities)

    # Store for PSI plot annotation
    psi_optimal_n = target_n

    # Find the threshold corresponding to this number of communities
    # The threshold is between the (n-1)th and n-th merge heights
    merge_heights = linkage_matrix[:, 2]
    if target_n >= 2 and target_n <= n_nodes:
        # For n communities, we cut just above the (n_nodes - n)th merge
        cut_idx = n_nodes - target_n
        if cut_idx < len(merge_heights) and cut_idx > 0:
            # Threshold between this merge and the next
            psi_threshold = (merge_heights[cut_idx - 1] + merge_heights[cut_idx]) / 2
        elif cut_idx == 0:
            psi_threshold = merge_heights[0] / 2
        else:
            psi_threshold = optimal_threshold
    else:
        psi_threshold = optimal_threshold

    # Use computed threshold for partitioning
    optimal_clusters = fcluster(linkage_matrix, t=psi_threshold, criterion="distance")
    n_clusters = len(np.unique(optimal_clusters))

    if verbose:
        print(f"Optimal threshold: {optimal_threshold:.6f}")
        print(f"Number of clusters: {n_clusters}")

    # Load channel labels using the proper helper function
    labels_list = _load_channel_labels(patient, dataset_root)
    if labels_list and len(labels_list) >= n_nodes:
        if _gc_nodes is not None:
            # Map contiguous indices 0..n_nodes-1 to original giant component channels
            channel_labels = {i: labels_list[orig] for i, orig in enumerate(_gc_nodes)}
        else:
            channel_labels = {i: labels_list[i] for i in range(n_nodes)}
        if verbose:
            print(f"Loaded {len(labels_list)} channel labels from file")
    else:
        if verbose:
            print(f"Using default channel labels (found {len(labels_list) if labels_list else 0})")
        channel_labels = {i: f"Ch{i}" for i in range(n_nodes)}

    # Create figure with custom grid layout (matching FIGMNTGN01)
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(4, 6, figure=fig, hspace=0.3, wspace=0.4)

    # Define subplot areas
    ax_matrix = fig.add_subplot(gs[0:2, 0:2])  # Top-left: Matrix
    ax_entropy = fig.add_subplot(gs[0:2, 2:4])  # Top-mid: Entropy/Heat
    ax_dendro = fig.add_subplot(gs[0:4, 4:6])  # Right: Dendrogram (full height)
    ax_network = fig.add_subplot(gs[2:4, 0:3])  # Bottom-left: Network
    ax_psi = fig.add_subplot(gs[2:4, 3:4])  # Bottom-mid-right: PSI

    # No fig.suptitle on publication figures; identifying metadata
    # (patient, phase, band, method) lives in the file name.

    # -------------------------------------------------------------------------
    # Panel (a): FC Matrix Heatmap
    # -------------------------------------------------------------------------
    im_matrix = ax_matrix.imshow(
        fc_matrix,
        cmap="viridis",
        vmin=0,
        vmax=1,
        origin="upper",
        aspect="equal",
    )

    ax_matrix.set_title(f"(a) {method_label} Matrix", fontsize=12, fontweight="bold")
    ax_matrix.set_xlabel("Channel", fontsize=10)
    ax_matrix.set_ylabel("Channel", fontsize=10)

    # Colorbar
    cbar_matrix = plt.colorbar(im_matrix, ax=ax_matrix, fraction=0.046, pad=0.04)
    cbar_matrix.set_label(method_label, fontsize=10)

    # -------------------------------------------------------------------------
    # Panel (b): Entropy and Specific Heat (CORRECTED - NO N multiplication)
    # -------------------------------------------------------------------------
    ax_entropy_twin = ax_entropy.twinx()

    # Plot specific heat (CORRECTED: use entropy_C directly, NO multiplication)
    line_heat = ax_entropy.plot(
        entropy_tau[1:], entropy_C, "-", label="C (Specific Heat)", color="blue", linewidth=2
    )
    ax_entropy.set_ylabel("C (Specific Heat)", color="blue", fontsize=11)
    ax_entropy.tick_params(axis="y", labelcolor="blue")

    # Plot entropy (1 - S normalized)
    line_entropy = ax_entropy_twin.plot(
        entropy_tau, entropy_1_minus_S, "-", label="1-S (Entropy)", color="red", linewidth=2
    )
    ax_entropy_twin.set_ylabel("1-S (Normalized Entropy)", color="red", fontsize=11)
    ax_entropy_twin.tick_params(axis="y", labelcolor="red")

    # Set x-axis (log scale for tau)
    ax_entropy.set_xscale("log")
    ax_entropy.set_xlabel(r"$\tau$ (Diffusion Time)", fontsize=11)
    ax_entropy.set_title("(b) Thermodynamic Observables", fontsize=12, fontweight="bold")

    # Combined legend
    lines_heat, labels_heat = ax_entropy.get_legend_handles_labels()
    lines_entropy, labels_entropy = ax_entropy_twin.get_legend_handles_labels()
    ax_entropy_twin.legend(lines_heat + lines_entropy, labels_heat + labels_entropy, loc="best", fontsize=9).set_zorder(
        200
    )

    ax_entropy.grid(alpha=0.3)

    # -------------------------------------------------------------------------
    # Panel (c): Dendrogram (CORRECTED settings from FIGMNTGN03)
    # -------------------------------------------------------------------------
    # Generate color palette
    n_colors = min(n_clusters + 3, 20)
    palette = plt.cm.tab20(np.linspace(0, 1, n_colors))
    palette = [plt.matplotlib.colors.to_hex(c) for c in palette]

    from scipy.cluster import hierarchy

    hierarchy.set_link_color_palette(palette)

    # Create dendrogram
    node_list = list(range(n_nodes))
    labels_for_dendro = [channel_labels.get(n, f"Ch{n}") for n in node_list]

    dendro = dendrogram(
        linkage_matrix,
        ax=ax_dendro,
        color_threshold=psi_threshold,  # Use PSI-based threshold
        labels=labels_for_dendro,
        above_threshold_color="k",
        leaf_font_size=5,
        orientation="right",
    )

    # Set log scale and limits (log-space padding for readability)
    merge_heights = linkage_matrix[:, 2]
    tmin = merge_heights[merge_heights > 0].min() * 0.5
    tmax = merge_heights.max() * 2.0
    ax_dendro.set_xscale("log")
    ax_dendro.axvline(
        psi_threshold, color="b", linestyle="--", linewidth=2,
        label=f"PSI cut (n={n_clusters})"
    )
    ax_dendro.set_xlim(tmin, tmax)
    ax_dendro.set_xlabel(r"$\mathcal{D}/\mathcal{D}_{\max}$", fontsize=11)
    ax_dendro.set_title("(c) Hierarchical Tree", fontsize=12, fontweight="bold")
    ax_dendro.legend(fontsize=9)

    # -------------------------------------------------------------------------
    # Panel (d): Network with Partition Colors
    # Uses power-law edge scaling and thresholded backbone for layout
    # -------------------------------------------------------------------------
    # Extract node colors from dendrogram
    leaf_label_colors = {lbl: col for lbl, col in zip(dendro["ivl"], dendro["leaves_color_list"])}

    # Map colors to nodes
    node_colors = [leaf_label_colors.get(channel_labels.get(n, f"Ch{n}"), "gray") for n in node_list]

    # Build full graph from FC matrix
    G = nx.from_numpy_array(fc_matrix)

    # Edge width parameters (same as Figure 1)
    EDGE_POWER = 2.0  # Power law exponent
    MAX_WIDTH = 4.0  # Max edge width

    # Spectral layout from graph Laplacian eigenvectors
    pos = nx.spectral_layout(G)

    # Get edge weights for full graph
    edges = list(G.edges(data=True))
    weights = np.array([e[2]['weight'] for e in edges])

    # Power law scaling for width (same as Figure 1)
    scaled = np.power(weights, EDGE_POWER)
    widths = MAX_WIDTH * scaled

    # Color edges using same colormap as matrix (viridis)
    cmap_edges = plt.cm.get_cmap('viridis')
    alphas = np.power(weights, EDGE_POWER)
    edge_colors = [(*cmap_edges(w)[:3], a) for w, a in zip(weights, alphas)]

    # Draw network with power-law scaled edges
    nx.draw_networkx_nodes(
        G, pos, ax=ax_network,
        node_size=80,
        node_color=node_colors,
        alpha=0.9,
        edgecolors='white',
        linewidths=0.5,
    )

    nx.draw_networkx_edges(
        G, pos, ax=ax_network,
        width=widths,
        edge_color=edge_colors,
    )

    ax_network.set_title(
        f"(d) Network (PSI Partition: {n_clusters} communities)", fontsize=12, fontweight="bold"
    )
    ax_network.axis('off')

    # -------------------------------------------------------------------------
    # Panel (e): PSI Plot
    # -------------------------------------------------------------------------
    ax_psi.plot(psi_n_communities, psi_values, "-o", color="green", linewidth=2, markersize=4)

    # Mark all PSI peaks for reference
    if len(psi_values) > 2:
        # Find local maxima
        for i in range(1, len(psi_values) - 1):
            if psi_values[i] > psi_values[i - 1] and psi_values[i] > psi_values[i + 1]:
                ax_psi.plot(psi_n_communities[i], psi_values[i], 'ro', markersize=8, alpha=0.5)

    # Mark the selected partition
    partition_label = f"Selected: n={psi_optimal_n}"
    if n_communities is not None:
        partition_label += " (fixed)"
    if len(psi_values) > 0:
        ax_psi.axvline(
            psi_optimal_n,
            ls="--",
            c="blue",
            linewidth=2,
            label=partition_label,
        )
        # Mark the PSI value at selected partition
        if psi_optimal_n in psi_n_communities:
            idx = np.where(psi_n_communities == psi_optimal_n)[0][0]
            ax_psi.plot(psi_optimal_n, psi_values[idx], 'b*', markersize=15, zorder=10)

    ax_psi.set_xlabel(r"$n$ (Number of Communities)", fontsize=11)
    ax_psi.set_ylabel(r"$\Psi(n, \tau)$ (PSI)", fontsize=11)
    ax_psi.set_title("(e) Partition Stability Index", fontsize=12, fontweight="bold")
    ax_psi.legend(fontsize=9, loc="upper right")
    ax_psi.grid(alpha=0.3)
    if len(psi_n_communities) > 0:
        ax_psi.set_xlim(1, min(30, max(psi_n_communities)))

    # -------------------------------------------------------------------------
    # Save figure
    # -------------------------------------------------------------------------
    if output_path is None:
        output_dir = FIGURES_ROOT / "lrg" / patient
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{band}_{phase}_lrg_{fc_method}_full.png"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    if verbose:
        print(f"Saved LRG analysis: {output_path}")

    return output_path
