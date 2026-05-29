"""3D Spatial Brain Network Visualization -- core Plotly / Matplotlib plots.

Holds the core electrode-network 3D plots (Plotly + Matplotlib):
`build_edge_traces`, `plot_spatial_network_3d`,
`plot_spatial_network_3d_mpl`, `plot_spatial_clusters_comparison`.

After the 2026-05-29 split (Phase 4-B split 6/7) this module also
re-exports the symbols that moved to `spatial_coords.py`
(`load_spatial_metadata`, `prepare_spatial_coordinates`) and to
`spatial_nilearn.py` (`view_brain_connectome`, `plot_brain_connectome`,
`plot_brain_connectome_at_n`), so callers importing from
`lrg_eegfc.visuals.spatial` continue to work unchanged.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.config.paths import CORR_CACHE, LRG_CACHE, MSC_CACHE, SEEG_DATAPATH, FIGURES_ROOT

# Re-export moved-out symbols so existing callers don't break.
from .spatial_coords import (  # noqa: F401
    load_spatial_metadata,
    prepare_spatial_coordinates,
)
from .spatial_nilearn import (  # noqa: F401
    view_brain_connectome,
    plot_brain_connectome,
    plot_brain_connectome_at_n,
)


__all__ = [
    # Coordinate loaders (re-exported from spatial_coords)
    "load_spatial_metadata",
    "prepare_spatial_coordinates",
    # Core 3D plots
    "build_edge_traces",
    "plot_spatial_network_3d",
    "plot_spatial_network_3d_mpl",
    "plot_spatial_clusters_comparison",
    # nilearn glass-brain wrappers (re-exported from spatial_nilearn)
    "view_brain_connectome",
    "plot_brain_connectome",
    "plot_brain_connectome_at_n",
]


def build_edge_traces(
    coords: np.ndarray,
    adjacency: np.ndarray,
    threshold: float = 0.0,
    max_edges: Optional[int] = None,
    edge_width: float = 1.0,
    edge_color: str = "rgba(100, 100, 100, 0.3)",
    weight_opacity: bool = True,
) -> go.Scatter3d:
    """Build Plotly edge traces for network visualization.

    Uses a single trace with None separators for efficient rendering of
    many edges.

    Parameters
    ----------
    coords : np.ndarray
        Node coordinates of shape (n_nodes, 3).
    adjacency : np.ndarray
        Adjacency/weight matrix of shape (n_nodes, n_nodes).
    threshold : float, optional
        Minimum edge weight to include. Default is 0.0.
    max_edges : int, optional
        Maximum number of edges to display (strongest edges kept).
        If None, display all edges above threshold.
    edge_width : float, optional
        Base line width for edges. Default is 1.0.
    edge_color : str, optional
        Edge color in rgba format. Default is semi-transparent gray.
    weight_opacity : bool, optional
        If True, scale edge opacity by weight. Default is True.

    Returns
    -------
    trace : go.Scatter3d
        Plotly trace containing all edge lines.
    """
    n_nodes = coords.shape[0]

    # Get upper triangle indices (avoid duplicates for symmetric matrix)
    i_upper, j_upper = np.triu_indices(n_nodes, k=1)
    weights = adjacency[i_upper, j_upper]

    # Filter by threshold
    mask = weights >= threshold
    i_edges = i_upper[mask]
    j_edges = j_upper[mask]
    edge_weights = weights[mask]

    # Limit number of edges if specified
    if max_edges is not None and len(edge_weights) > max_edges:
        top_indices = np.argsort(edge_weights)[-max_edges:]
        i_edges = i_edges[top_indices]
        j_edges = j_edges[top_indices]
        edge_weights = edge_weights[top_indices]

    # Build coordinate arrays with None separators
    x_edges: List[Optional[float]] = []
    y_edges: List[Optional[float]] = []
    z_edges: List[Optional[float]] = []

    for i, j in zip(i_edges, j_edges):
        x_edges.extend([coords[i, 0], coords[j, 0], None])
        y_edges.extend([coords[i, 1], coords[j, 1], None])
        z_edges.extend([coords[i, 2], coords[j, 2], None])

    trace = go.Scatter3d(
        x=x_edges,
        y=y_edges,
        z=z_edges,
        mode="lines",
        line=dict(color=edge_color, width=edge_width),
        hoverinfo="skip",
        name="Edges",
    )

    return trace


def _get_cluster_colors(
    labels: np.ndarray,
    colorscale: str = "Plotly",
) -> Tuple[List[str], Dict[int, str]]:
    """Generate colors for cluster labels.

    Parameters
    ----------
    labels : np.ndarray
        Cluster labels for each node.
    colorscale : str
        Plotly colorscale name.

    Returns
    -------
    colors : list
        Color string for each node.
    color_map : dict
        Mapping from cluster ID to color.
    """
    import plotly.colors as pc

    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)

    # Use qualitative colorscale for clusters
    if n_clusters <= 10:
        palette = pc.qualitative.Plotly
    elif n_clusters <= 20:
        palette = pc.qualitative.Alphabet
    else:
        palette = pc.qualitative.Light24 + pc.qualitative.Dark24

    color_map = {lab: palette[i % len(palette)] for i, lab in enumerate(unique_labels)}
    colors = [color_map[lab] for lab in labels]

    return colors, color_map


def plot_spatial_network_3d(
    patient: str,
    phase: str,
    band: str,
    fc_method: str = "imcoh_abs",
    dataset_root: Path = SEEG_DATAPATH,
    lrg_cache_root: Path = LRG_CACHE,
    fc_cache_root: Optional[Path] = None,
    edge_threshold: float = 0.3,
    max_edges: int = 500,
    node_size: int = 8,
    show_edges: bool = True,
    colorby: str = "cluster",
    title: Optional[str] = None,
    width: int = 900,
    height: int = 700,
) -> go.Figure:
    """Create interactive 3D network visualization with Plotly.

    Parameters
    ----------
    patient : str
        Patient identifier (e.g., "Pat_02").
    phase : str
        Recording phase (e.g., "rest_pre", "rest_post").
    band : str
        Frequency band (e.g., "beta", "alpha").
    fc_method : str, optional
        Functional connectivity method. Default ``"imcoh_abs"`` (current era).
        Other options: ``"msc"``, ``"corr"``, ``"imcoh"``, ``"imcoh_sq"``.
    dataset_root : Path, optional
        Root directory for patient data.
    lrg_cache_root : Path, optional
        Root directory for LRG cache.
    fc_cache_root : Path, optional
        Root directory for FC cache. If None, uses default based on fc_method.
    edge_threshold : float, optional
        Minimum edge weight to display. Default is 0.3.
    max_edges : int, optional
        Maximum number of edges to display. Default is 500.
    node_size : int, optional
        Node marker size. Default is 8.
    show_edges : bool, optional
        Whether to display network edges. Default is True.
    colorby : str, optional
        Node coloring mode: "cluster" (LRG clusters) or "atlas" (Desikan-Killiany).
        Default is "cluster".
    title : str, optional
        Plot title. Auto-generated if None.
    width : int, optional
        Figure width in pixels.
    height : int, optional
        Figure height in pixels.

    Returns
    -------
    fig : go.Figure
        Interactive Plotly figure.

    Raises
    ------
    FileNotFoundError
        If required data files are not found.
    """
    from lrg_eegfc.workflow.lrg import load_lrg_result

    # Load metadata with coordinates (using label normalization)
    metadata = load_spatial_metadata(patient, dataset_root)

    # Load LRG result for cluster labels
    lrg_result = load_lrg_result(patient, phase, band, fc_method, lrg_cache_root)
    if lrg_result is None:
        raise FileNotFoundError(
            f"LRG result not found for {patient} {phase} {band} ({fc_method})"
        )

    # Get cluster labels
    cluster_labels = fcluster(
        lrg_result.linkage_matrix,
        t=lrg_result.optimal_threshold,
        criterion="distance",
    )

    # Prepare coordinates (use MNI transform for proper brain positioning)
    coords = prepare_spatial_coordinates(metadata, scale="mm", center=False, to_mni=True)

    # Prepare node colors and hover info
    if colorby == "cluster":
        colors, color_map = _get_cluster_colors(cluster_labels)
        n_clusters = len(np.unique(cluster_labels))
    else:
        # Default to cluster coloring
        colors, color_map = _get_cluster_colors(cluster_labels)
        n_clusters = len(np.unique(cluster_labels))

    # Build hover text
    labels = metadata["label"].tolist() if "label" in metadata.columns else [f"Ch{i}" for i in range(len(coords))]
    hover_texts = []
    for i, label in enumerate(labels):
        text = f"<b>{label}</b><br>Cluster: {cluster_labels[i]}"
        if "Desikan-Killany" in metadata.columns:
            atlas = metadata["Desikan-Killany"].iloc[i]
            if pd.notna(atlas):
                # Parse atlas string (first region)
                region = str(atlas).split(",")[0].strip().strip('"')
                text += f"<br>Region: {region}"
        hover_texts.append(text)

    # Create figure
    fig = go.Figure()

    # Add edges if requested
    if show_edges:
        # Load FC matrix for edges
        from lrg_eegfc.workflow.fc import load_fc_matrix as _load_fc
        fc_kw = {}
        if fc_cache_root is not None:
            fc_kw["cache_root"] = fc_cache_root
        fc_matrix = _load_fc(patient, phase, band, fc_method, **fc_kw)

        if fc_matrix is not None:
            edge_trace = build_edge_traces(
                coords,
                fc_matrix,
                threshold=edge_threshold,
                max_edges=max_edges,
                edge_width=1.0,
                edge_color="rgba(150, 150, 150, 0.4)",
            )
            fig.add_trace(edge_trace)

    # Add nodes
    node_trace = go.Scatter3d(
        x=coords[:, 0],
        y=coords[:, 1],
        z=coords[:, 2],
        mode="markers",
        marker=dict(
            size=node_size,
            color=colors,
            line=dict(width=0.5, color="white"),
            opacity=0.9,
        ),
        text=hover_texts,
        hoverinfo="text",
        name="Electrodes",
    )
    fig.add_trace(node_trace)

    # Layout
    if title is None:
        title = f"{patient} {phase} {band} ({fc_method}) - {n_clusters} clusters"

    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(size=16)),
        width=width,
        height=height,
        scene=dict(
            xaxis=dict(title="X (mm)", showgrid=True, gridcolor="lightgray"),
            yaxis=dict(title="Y (mm)", showgrid=True, gridcolor="lightgray"),
            zaxis=dict(title="Z (mm)", showgrid=True, gridcolor="lightgray"),
            aspectmode="data",
        ),
        showlegend=False,
        margin=dict(l=0, r=0, t=40, b=0),
    )

    return fig


def plot_spatial_network_3d_mpl(
    patient: str,
    phase: str,
    band: str,
    fc_method: str = "imcoh_abs",
    dataset_root: Path = SEEG_DATAPATH,
    lrg_cache_root: Path = LRG_CACHE,
    fc_cache_root: Optional[Path] = None,
    edge_threshold: float = 0.3,
    max_edges: int = 200,
    node_size: int = 40,
    show_edges: bool = True,
    colorby: str = "cluster",
    title: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 8),
    output_path: Optional[Path] = None,
    elev: float = 20,
    azim: float = 45,
) -> Path:
    """Create static 3D network visualization with Matplotlib.

    Parameters
    ----------
    patient : str
        Patient identifier.
    phase : str
        Recording phase.
    band : str
        Frequency band.
    fc_method : str, optional
        FC method. Default ``"imcoh_abs"`` (current era). Other options:
        ``"msc"``, ``"corr"``, ``"imcoh"``, ``"imcoh_sq"``.
    dataset_root : Path, optional
        Root directory for patient data.
    lrg_cache_root : Path, optional
        Root directory for LRG cache.
    fc_cache_root : Path, optional
        Root directory for FC cache.
    edge_threshold : float, optional
        Minimum edge weight to display.
    max_edges : int, optional
        Maximum number of edges.
    node_size : int, optional
        Node marker size.
    show_edges : bool, optional
        Whether to show edges.
    colorby : str, optional
        Node coloring mode.
    title : str, optional
        Plot title.
    figsize : tuple, optional
        Figure size in inches.
    output_path : Path, optional
        Output file path. Auto-generated if None.
    elev : float, optional
        Elevation angle for 3D view.
    azim : float, optional
        Azimuth angle for 3D view.

    Returns
    -------
    output_path : Path
        Path to saved figure.
    """
    from matplotlib import cm

    from lrg_eegfc.workflow.lrg import load_lrg_result

    # Load data (using label normalization)
    metadata = load_spatial_metadata(patient, dataset_root)

    lrg_result = load_lrg_result(patient, phase, band, fc_method, lrg_cache_root)
    if lrg_result is None:
        raise FileNotFoundError(f"LRG result not found for {patient} {phase} {band}")

    cluster_labels = fcluster(
        lrg_result.linkage_matrix,
        t=lrg_result.optimal_threshold,
        criterion="distance",
    )

    coords = prepare_spatial_coordinates(metadata, scale="mm", center=False, to_mni=True)

    # Create figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    # Plot edges
    if show_edges:
        from lrg_eegfc.workflow.fc import load_fc_matrix as _load_fc
        fc_kw = {}
        if fc_cache_root is not None:
            fc_kw["cache_root"] = fc_cache_root
        fc_matrix = _load_fc(patient, phase, band, fc_method, **fc_kw)

        if fc_matrix is not None:
            n_nodes = coords.shape[0]
            i_upper, j_upper = np.triu_indices(n_nodes, k=1)
            weights = fc_matrix[i_upper, j_upper]
            mask = weights >= edge_threshold

            i_edges = i_upper[mask]
            j_edges = j_upper[mask]
            edge_weights = weights[mask]

            if max_edges is not None and len(edge_weights) > max_edges:
                top_idx = np.argsort(edge_weights)[-max_edges:]
                i_edges = i_edges[top_idx]
                j_edges = j_edges[top_idx]
                edge_weights = edge_weights[top_idx]

            for i, j, w in zip(i_edges, j_edges, edge_weights):
                ax.plot3D(
                    [coords[i, 0], coords[j, 0]],
                    [coords[i, 1], coords[j, 1]],
                    [coords[i, 2], coords[j, 2]],
                    color="gray",
                    alpha=0.2 + 0.3 * w,
                    linewidth=0.5,
                )

    # Plot nodes by cluster
    unique_clusters = np.unique(cluster_labels)
    n_clusters = len(unique_clusters)
    cmap = cm.get_cmap("tab20" if n_clusters <= 20 else "nipy_spectral")

    for idx, cluster_id in enumerate(unique_clusters):
        mask = cluster_labels == cluster_id
        color = cmap(idx / max(n_clusters - 1, 1))
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            coords[mask, 2],
            s=node_size,
            c=[color],
            label=f"Cluster {cluster_id}",
            alpha=0.85,
            edgecolors="white",
            linewidths=0.5,
        )

    # Labels and title
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_zlabel("Z (mm)")

    if title is None:
        title = f"{patient} {phase} {band} ({fc_method}) - {n_clusters} clusters"
    ax.set_title(title, fontsize=12)

    ax.view_init(elev=elev, azim=azim)

    if n_clusters <= 10:
        ax.legend(loc="upper left", fontsize=8, framealpha=0.8)

    # Save
    if output_path is None:
        output_dir = FIGURES_ROOT / "spatial" / patient
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{band}_{phase}_{fc_method}_spatial3d.png"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return output_path


def plot_spatial_clusters_comparison(
    patient: str,
    phase_a: str,
    phase_b: str,
    band: str,
    fc_method: str = "imcoh_abs",
    dataset_root: Path = SEEG_DATAPATH,
    lrg_cache_root: Path = LRG_CACHE,
    node_size: int = 8,
    title: Optional[str] = None,
    width: int = 1400,
    height: int = 600,
) -> go.Figure:
    """Create side-by-side 3D comparison of clusters between two phases.

    Parameters
    ----------
    patient : str
        Patient identifier.
    phase_a : str
        First phase (e.g., "rest_pre").
    phase_b : str
        Second phase (e.g., "rest_post").
    band : str
        Frequency band.
    fc_method : str, optional
        FC method. Default ``"imcoh_abs"`` (current era).
    dataset_root : Path, optional
        Root directory for patient data.
    lrg_cache_root : Path, optional
        Root directory for LRG cache.
    node_size : int, optional
        Node marker size.
    title : str, optional
        Overall figure title.
    width : int, optional
        Figure width in pixels.
    height : int, optional
        Figure height in pixels.

    Returns
    -------
    fig : go.Figure
        Interactive Plotly figure with two subplots.
    """
    from plotly.subplots import make_subplots

    from lrg_eegfc.workflow.lrg import load_lrg_result

    # Load metadata (using label normalization)
    metadata = load_spatial_metadata(patient, dataset_root)

    # Load LRG results
    result_a = load_lrg_result(patient, phase_a, band, fc_method, lrg_cache_root)
    result_b = load_lrg_result(patient, phase_b, band, fc_method, lrg_cache_root)

    if result_a is None or result_b is None:
        missing = []
        if result_a is None:
            missing.append(phase_a)
        if result_b is None:
            missing.append(phase_b)
        raise FileNotFoundError(f"LRG result not found for phases: {missing}")

    # Get cluster labels
    labels_a = fcluster(
        result_a.linkage_matrix, t=result_a.optimal_threshold, criterion="distance"
    )
    labels_b = fcluster(
        result_b.linkage_matrix, t=result_b.optimal_threshold, criterion="distance"
    )

    coords = prepare_spatial_coordinates(metadata, scale="mm", center=False, to_mni=True)

    # Get colors (use same palette for comparability)
    colors_a, _ = _get_cluster_colors(labels_a)
    colors_b, _ = _get_cluster_colors(labels_b)

    # Channel labels for hover
    ch_labels = (
        metadata["label"].tolist()
        if "label" in metadata.columns
        else [f"Ch{i}" for i in range(len(coords))]
    )

    hover_a = [f"<b>{ch}</b><br>Cluster: {c}" for ch, c in zip(ch_labels, labels_a)]
    hover_b = [f"<b>{ch}</b><br>Cluster: {c}" for ch, c in zip(ch_labels, labels_b)]

    # Create subplots
    fig = make_subplots(
        rows=1,
        cols=2,
        specs=[[{"type": "scene"}, {"type": "scene"}]],
        subplot_titles=[
            f"{phase_a} ({len(np.unique(labels_a))} clusters)",
            f"{phase_b} ({len(np.unique(labels_b))} clusters)",
        ],
        horizontal_spacing=0.02,
    )

    # Add traces
    fig.add_trace(
        go.Scatter3d(
            x=coords[:, 0],
            y=coords[:, 1],
            z=coords[:, 2],
            mode="markers",
            marker=dict(
                size=node_size,
                color=colors_a,
                line=dict(width=0.5, color="white"),
                opacity=0.9,
            ),
            text=hover_a,
            hoverinfo="text",
            name=phase_a,
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter3d(
            x=coords[:, 0],
            y=coords[:, 1],
            z=coords[:, 2],
            mode="markers",
            marker=dict(
                size=node_size,
                color=colors_b,
                line=dict(width=0.5, color="white"),
                opacity=0.9,
            ),
            text=hover_b,
            hoverinfo="text",
            name=phase_b,
        ),
        row=1,
        col=2,
    )

    # Shared camera settings
    camera = dict(eye=dict(x=1.5, y=1.5, z=1.0))
    scene_layout = dict(
        xaxis=dict(title="X (mm)", showgrid=True),
        yaxis=dict(title="Y (mm)", showgrid=True),
        zaxis=dict(title="Z (mm)", showgrid=True),
        aspectmode="data",
        camera=camera,
    )

    if title is None:
        title = f"{patient} {band} ({fc_method}): {phase_a} vs {phase_b}"

    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(size=16)),
        width=width,
        height=height,
        scene=scene_layout,
        scene2=scene_layout,
        showlegend=False,
        margin=dict(l=0, r=0, t=60, b=0),
    )

    return fig
