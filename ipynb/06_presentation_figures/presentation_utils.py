"""Utility functions for presentation figures.

This module contains specialized visualization and analysis functions
for the presentation notebooks.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
import networkx as nx
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, optimal_leaf_ordering
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score

# Professional color schemes
CMAP_MATRIX = "magma"  # For FC matrices - perceptually uniform, professional
CMAP_NETWORK = "viridis"  # For network edge coloring
CMAP_ULTRAMETRIC = "cividis"  # For ultrametric distances - colorblind friendly
CMAP_DIFF = "RdBu_r"  # For difference matrices

# Cluster color palette (tab20 is good for categorical data)
CLUSTER_PALETTE = plt.cm.tab20


def get_cluster_colors(labels: np.ndarray, palette=CLUSTER_PALETTE) -> List[str]:
    """Get consistent colors for cluster labels.

    Parameters
    ----------
    labels : np.ndarray
        Cluster labels (1-indexed from fcluster)
    palette : colormap
        Matplotlib colormap to use

    Returns
    -------
    List[str]
        List of hex color strings for each node
    """
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    colors = [palette(i / max(n_clusters, 1)) for i in range(n_clusters)]
    color_map = {label: plt.matplotlib.colors.to_hex(colors[i])
                 for i, label in enumerate(unique_labels)}
    return [color_map[label] for label in labels]


def compute_rho_matrix(L_normalized: np.ndarray, tau: float) -> np.ndarray:
    """Compute density matrix rho(tau) = exp(-tau * L).

    The density matrix is the Laplacian propagator at diffusion time tau.
    It represents how information spreads across the network.

    Parameters
    ----------
    L_normalized : np.ndarray
        Normalized Laplacian matrix
    tau : float
        Diffusion time parameter

    Returns
    -------
    np.ndarray
        Density matrix rho(tau)
    """
    from scipy.linalg import expm
    return expm(-tau * L_normalized)


def compute_normalized_laplacian(adjacency: np.ndarray) -> np.ndarray:
    """Compute normalized Laplacian from adjacency matrix.

    L_normalized = I - D^{-1/2} A D^{-1/2}

    Parameters
    ----------
    adjacency : np.ndarray
        Adjacency matrix (weighted)

    Returns
    -------
    np.ndarray
        Normalized Laplacian matrix
    """
    # Degree matrix
    degrees = np.sum(adjacency, axis=1)
    # Handle zero degrees
    degrees[degrees == 0] = 1.0
    D_inv_sqrt = np.diag(1.0 / np.sqrt(degrees))

    # Normalized Laplacian: L = I - D^{-1/2} A D^{-1/2}
    I = np.eye(adjacency.shape[0])
    L = I - D_inv_sqrt @ adjacency @ D_inv_sqrt
    return L


def create_rho_animation(
    adjacency: np.ndarray,
    output_path: Path,
    n_frames: int = 100,
    tau_range: Tuple[float, float] = (-3, 3),
    fps: int = 15,
    figsize: Tuple[int, int] = (8, 8),
    cmap: str = "inferno",
    title_template: str = r"$\rho(\tau)$ - Laplacian Propagator at $\tau$ = {tau:.3f}",
) -> Path:
    """Create animated GIF showing the filling of rho(tau) matrix.

    The colormap is normalized dynamically:
    - vmax = max value at tau_start (near identity, diagonal ~1)
    - vmin = min value at tau_end (uniform ~1/N)

    This shows the dynamic evolution as information spreads across the network.

    Parameters
    ----------
    adjacency : np.ndarray
        Adjacency matrix
    output_path : Path
        Output path for GIF
    n_frames : int
        Number of frames in animation
    tau_range : Tuple[float, float]
        Range of log10(tau) values
    fps : int
        Frames per second
    figsize : Tuple[int, int]
        Figure size
    cmap : str
        Colormap
    title_template : str
        Title template with {tau} placeholder

    Returns
    -------
    Path
        Path to saved GIF
    """
    # Compute normalized Laplacian
    L = compute_normalized_laplacian(adjacency)

    # Generate tau values (log-spaced)
    tau_values = np.logspace(tau_range[0], tau_range[1], n_frames)

    # Pre-compute all rho matrices for smooth animation
    print(f"Pre-computing {n_frames} rho matrices...")
    rho_matrices = []
    for tau in tau_values:
        rho = compute_rho_matrix(L, tau)
        rho_matrices.append(rho)

    # Dynamic colormap normalization:
    # vmax = max at beginning (tau→0, near identity)
    # vmin = min at end (tau→∞, uniform 1/N)
    vmax = rho_matrices[0].max()  # Start: diagonal ~1
    vmin = rho_matrices[-1].min()  # End: uniform ~1/N

    print(f"Colormap range: vmin={vmin:.4f}, vmax={vmax:.4f}")

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Initial plot with dynamic normalization
    im = ax.imshow(rho_matrices[0], cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal')
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$\rho_{ij}(\tau)$", fontsize=12)
    title = ax.set_title(title_template.format(tau=tau_values[0]), fontsize=14)
    ax.set_xlabel("Node", fontsize=11)
    ax.set_ylabel("Node", fontsize=11)

    def update(frame):
        im.set_array(rho_matrices[frame])
        title.set_text(title_template.format(tau=tau_values[frame]))
        return [im, title]

    print(f"Creating animation with {n_frames} frames at {fps} fps...")
    anim = animation.FuncAnimation(
        fig, update, frames=n_frames, interval=1000/fps, blit=True
    )

    # Save as GIF
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    writer = animation.PillowWriter(fps=fps)
    anim.save(str(output_path), writer=writer, dpi=100)
    plt.close(fig)

    print(f"Saved animation to {output_path}")
    return output_path


def compute_multiscale_ari(
    linkage1: np.ndarray,
    linkage2: np.ndarray,
    n_nodes: int,
    n_thresholds: int = 50,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute ARI between two dendrograms at multiple threshold levels.

    At each normalized threshold (0 to 1), cut both dendrograms and
    compute the Adjusted Rand Index between the resulting partitions.

    Parameters
    ----------
    linkage1 : np.ndarray
        First linkage matrix
    linkage2 : np.ndarray
        Second linkage matrix
    n_nodes : int
        Number of nodes
    n_thresholds : int
        Number of threshold levels to evaluate

    Returns
    -------
    thresholds : np.ndarray
        Normalized threshold values (0 to 1)
    ari_values : np.ndarray
        ARI values at each threshold
    """
    # Get max heights from both dendrograms
    max_height1 = linkage1[:, 2].max()
    max_height2 = linkage2[:, 2].max()

    # Use normalized thresholds (fraction of max height)
    thresholds = np.linspace(0.01, 1.0, n_thresholds)
    ari_values = []

    for t in thresholds:
        # Cut both dendrograms at the same normalized height
        threshold1 = t * max_height1
        threshold2 = t * max_height2

        labels1 = fcluster(linkage1, t=threshold1, criterion='distance')
        labels2 = fcluster(linkage2, t=threshold2, criterion='distance')

        # Compute ARI
        ari = adjusted_rand_score(labels1, labels2)
        ari_values.append(ari)

    return thresholds, np.array(ari_values)


def compute_multiscale_structural_distance(
    linkage1: np.ndarray,
    linkage2: np.ndarray,
    n_nodes: int,
    n_thresholds: int = 50,
) -> float:
    """Compute multiscale structural distance between two hierarchical trees.

    The distance is computed as the integral of (1 - ARI) across all scales:
    D = integral(1 - ARI(h)) dh

    This captures the total divergence between structures across all scales.
    A distance of 0 means identical structures at all scales.
    A distance of 1 means completely different at all scales.

    Parameters
    ----------
    linkage1 : np.ndarray
        First linkage matrix
    linkage2 : np.ndarray
        Second linkage matrix
    n_nodes : int
        Number of nodes
    n_thresholds : int
        Number of threshold levels

    Returns
    -------
    float
        Multiscale structural distance (0 to 1)
    """
    thresholds, ari_values = compute_multiscale_ari(
        linkage1, linkage2, n_nodes, n_thresholds
    )

    # Integrate (1 - ARI) using trapezoidal rule
    divergence = 1 - ari_values
    distance = np.trapz(divergence, thresholds)

    # Normalize by the integration range (0.01 to 1.0)
    distance = distance / (thresholds[-1] - thresholds[0])

    return distance


def compute_phase_distance_matrix_multiscale(
    lrg_results: Dict[str, "LRGResult"],
    phases: List[str],
    n_thresholds: int = 50,
) -> np.ndarray:
    """Compute pairwise multiscale structural distances between phases.

    Parameters
    ----------
    lrg_results : Dict[str, LRGResult]
        Dictionary of LRG results keyed by phase name
    phases : List[str]
        List of phase names in order
    n_thresholds : int
        Number of threshold levels

    Returns
    -------
    np.ndarray
        Distance matrix (n_phases x n_phases)
    """
    n_phases = len(phases)
    dist_matrix = np.zeros((n_phases, n_phases))

    for i, phase_i in enumerate(phases):
        for j, phase_j in enumerate(phases):
            if i < j:
                result_i = lrg_results[phase_i]
                result_j = lrg_results[phase_j]

                dist = compute_multiscale_structural_distance(
                    result_i.linkage_matrix,
                    result_j.linkage_matrix,
                    result_i.n_nodes,
                    n_thresholds,
                )
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist  # Symmetric

    return dist_matrix


def plot_network_with_partition(
    ax: plt.Axes,
    adjacency: np.ndarray,
    cluster_labels: np.ndarray,
    pos: Optional[Dict] = None,
    layout: str = "spring",
    k: float = 0.1,
    seed: int = 42,
    node_size: int = 100,
    edge_alpha: float = 0.3,
    edge_width_scale: float = 1.0,
    title: str = "",
) -> Dict:
    """Plot network with nodes colored by partition.

    Parameters
    ----------
    ax : plt.Axes
        Matplotlib axes
    adjacency : np.ndarray
        Adjacency matrix
    cluster_labels : np.ndarray
        Cluster labels for each node
    pos : Dict, optional
        Pre-computed node positions
    layout : str
        Layout algorithm: "spring" or "spectral"
    k : float
        Spring layout parameter (lower = tighter clusters)
    seed : int
        Random seed for layout
    node_size : int
        Node size
    edge_alpha : float
        Edge transparency
    edge_width_scale : float
        Scale factor for edge widths
    title : str
        Plot title

    Returns
    -------
    Dict
        Node positions (for reuse)
    """
    G = nx.from_numpy_array(adjacency)

    # Compute layout if not provided
    if pos is None:
        if layout == "spring":
            pos = nx.spring_layout(G, seed=seed, k=k, iterations=50)
        elif layout == "spectral":
            try:
                pos = nx.spectral_layout(G)
            except:
                pos = nx.spring_layout(G, seed=seed, k=k, iterations=50)
        else:
            pos = nx.spring_layout(G, seed=seed, k=k, iterations=50)

    # Get node colors from cluster labels
    node_colors = get_cluster_colors(cluster_labels)

    # Edge widths proportional to weight
    widths = [G[u][v]["weight"] * edge_width_scale for u, v in G.edges()]

    # Draw network
    nx.draw(
        G,
        pos=pos,
        ax=ax,
        node_color=node_colors,
        node_size=node_size,
        width=widths,
        edge_color="gray",
        alpha=edge_alpha,
        with_labels=False,
    )

    ax.set_title(title, fontsize=12, fontweight="bold")

    return pos


def plot_dendrogram_with_partition(
    ax: plt.Axes,
    linkage: np.ndarray,
    threshold: float,
    orientation: str = "top",
    labels: Optional[List[str]] = None,
    title: str = "",
    show_threshold_line: bool = True,
    log_scale: bool = True,
) -> Dict:
    """Plot dendrogram with optimal threshold line.

    Parameters
    ----------
    ax : plt.Axes
        Matplotlib axes
    linkage : np.ndarray
        Linkage matrix
    threshold : float
        Optimal threshold for coloring
    orientation : str
        Dendrogram orientation
    labels : List[str], optional
        Node labels
    title : str
        Plot title
    show_threshold_line : bool
        Whether to show threshold line
    log_scale : bool
        Whether to use log scale on distance axis

    Returns
    -------
    Dict
        Dendrogram result dictionary
    """
    from scipy.cluster import hierarchy

    # Color palette
    n_colors = 20
    palette = [plt.matplotlib.colors.to_hex(CLUSTER_PALETTE(i / n_colors))
               for i in range(n_colors)]
    hierarchy.set_link_color_palette(palette)

    # Create dendrogram
    dendro = dendrogram(
        linkage,
        ax=ax,
        orientation=orientation,
        color_threshold=threshold,
        above_threshold_color="black",
        labels=labels,
        leaf_font_size=6 if labels and len(labels) <= 50 else 4,
        no_labels=(labels is None or len(labels) > 50),
    )

    # Add threshold line
    if show_threshold_line:
        if orientation in ["top", "bottom"]:
            ax.axhline(threshold, color="blue", linestyle="--", lw=2,
                      label=f"Threshold = {threshold:.3f}")
            if log_scale:
                ax.set_yscale("log")
            ax.set_ylabel("Ultrametric Distance", fontsize=10)
        else:
            ax.axvline(threshold, color="blue", linestyle="--", lw=2,
                      label=f"Threshold = {threshold:.3f}")
            if log_scale:
                ax.set_xscale("log")
            ax.set_xlabel("Ultrametric Distance", fontsize=10)

    ax.set_title(title, fontsize=12, fontweight="bold")

    return dendro


def plot_triangular_distance_matrix(
    ax: plt.Axes,
    distance_matrix: np.ndarray,
    labels: List[str],
    title: str = "",
    cmap: str = "YlOrRd",
    annotate: bool = True,
    vmin: float = 0.0,
    vmax: Optional[float] = None,
) -> None:
    """Plot lower triangular distance matrix.

    Parameters
    ----------
    ax : plt.Axes
        Matplotlib axes
    distance_matrix : np.ndarray
        Square distance matrix
    labels : List[str]
        Labels for rows/columns
    title : str
        Plot title
    cmap : str
        Colormap
    annotate : bool
        Whether to annotate cells with values
    vmin : float
        Minimum value for colormap
    vmax : float, optional
        Maximum value for colormap
    """
    n = len(labels)

    # Create mask for upper triangle
    mask = np.triu(np.ones_like(distance_matrix, dtype=bool))

    # Masked array for plotting
    masked_data = np.ma.masked_array(distance_matrix, mask=mask)

    if vmax is None:
        vmax = np.nanmax(distance_matrix[~mask])

    # Plot
    im = ax.imshow(masked_data, cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Distance", fontsize=10)

    # Set ticks
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=10)
    ax.set_yticklabels(labels, fontsize=10)

    # Annotate cells
    if annotate:
        for i in range(n):
            for j in range(n):
                if i > j:  # Lower triangle only
                    ax.text(j, i, f"{distance_matrix[i, j]:.2f}",
                           ha="center", va="center", fontsize=9,
                           color="white" if distance_matrix[i, j] > vmax/2 else "black")

    ax.set_title(title, fontsize=12, fontweight="bold")


def create_four_phase_comparison(
    lrg_results: Dict[str, "LRGResult"],
    fc_matrices: Dict[str, np.ndarray],
    phases: List[str],
    band: str,
    output_path: Path,
    layout: str = "spring",
    k: float = 0.1,
    figsize: Tuple[int, int] = (20, 12),
) -> Path:
    """Create 4-column phase comparison figure.

    Top row: Networks with partition colors
    Bottom row: Dendrograms with partition

    Parameters
    ----------
    lrg_results : Dict[str, LRGResult]
        LRG results keyed by phase
    fc_matrices : Dict[str, np.ndarray]
        FC matrices keyed by phase
    phases : List[str]
        Phase names in order
    band : str
        Frequency band name
    output_path : Path
        Output path
    layout : str
        Network layout algorithm
    k : float
        Spring layout parameter
    figsize : Tuple[int, int]
        Figure size

    Returns
    -------
    Path
        Path to saved figure
    """
    n_phases = len(phases)
    fig, axes = plt.subplots(2, n_phases, figsize=figsize)

    # Compute shared layout based on first phase
    first_phase = phases[0]
    G_first = nx.from_numpy_array(fc_matrices[first_phase])
    if layout == "spring":
        shared_pos = nx.spring_layout(G_first, seed=42, k=k, iterations=50)
    else:
        try:
            shared_pos = nx.spectral_layout(G_first)
        except:
            shared_pos = nx.spring_layout(G_first, seed=42, k=k, iterations=50)

    for i, phase in enumerate(phases):
        result = lrg_results[phase]
        fc_matrix = fc_matrices[phase]

        # Get partition at optimal threshold
        cluster_labels = fcluster(
            result.linkage_matrix,
            t=result.optimal_threshold,
            criterion='distance'
        )
        n_clusters = len(np.unique(cluster_labels))

        # Top row: Network
        plot_network_with_partition(
            axes[0, i],
            fc_matrix,
            cluster_labels,
            pos=shared_pos,
            layout=layout,
            k=k,
            title=f"{phase}\n({n_clusters} communities)",
        )

        # Bottom row: Dendrogram
        plot_dendrogram_with_partition(
            axes[1, i],
            result.linkage_matrix,
            result.optimal_threshold,
            orientation="top",
            title="",
            log_scale=True,
        )

    fig.suptitle(f"Cross-Phase Comparison - {band.upper()}",
                 fontsize=16, fontweight="bold", y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return output_path
