"""nilearn glass-brain wrappers for SEEG electrode networks.

Split out of `spatial.py` on 2026-05-29 (Phase 4-B split 6/7). Holds
`view_brain_connectome` (interactive nilearn `view_connectome`),
`plot_brain_connectome` (static glass-brain), and the
n-communities helper `plot_brain_connectome_at_n`.

Re-exported from `spatial.py` for backwards compatibility.
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.config.paths import CORR_CACHE, LRG_CACHE, MSC_CACHE, SEEG_DATAPATH

from .spatial_coords import load_spatial_metadata, prepare_spatial_coordinates


__all__ = [
    "view_brain_connectome",
    "plot_brain_connectome",
    "plot_brain_connectome_at_n",
]


def view_brain_connectome(
    patient: str,
    phase: str,
    band: str,
    fc_method: str = "imcoh_abs",
    dataset_root: Path = SEEG_DATAPATH,
    fc_cache_root: Optional[Path] = None,
    lrg_cache_root: Optional[Path] = None,
    edge_threshold: Optional[float] = None,
    edge_cmap: str = "bwr",
    symmetric_cmap: bool = True,
    linewidth: float = 6.0,
    node_color: Optional[Union[str, List, np.ndarray]] = "auto",
    node_size: float = 6.0,
    colorbar: bool = True,
    title: Optional[str] = None,
):
    """Interactive 3D brain connectome visualization using nilearn.

    Displays electrodes embedded in a 3D glass brain with connections
    colored by edge weights.

    Parameters
    ----------
    patient : str
        Patient identifier (e.g., "Pat_02").
    phase : str
        Recording phase (e.g., "rest_pre", "rest_post").
    band : str
        Frequency band (e.g., "beta", "alpha").
    fc_method : str, optional
        FC method. Default ``"imcoh_abs"`` (current era). Other options:
        ``"msc"``, ``"corr"``, ``"imcoh"``, ``"imcoh_sq"``.
    dataset_root : Path, optional
        Root directory for patient data.
    fc_cache_root : Path, optional
        Root directory for FC cache.
    lrg_cache_root : Path, optional
        Root directory for LRG cache (for cluster colors).
    edge_threshold : float or str, optional
        If None, show all edges. If float, minimum edge weight.
        If string like "25%", show top percentile.
    edge_cmap : str, optional
        Colormap for edges. Default is "bwr" (blue-white-red).
        Other options: "viridis", "plasma", "coolwarm", "RdBu_r".
    symmetric_cmap : bool, optional
        If True, colormap is symmetric around 0. Default is True.
    linewidth : float, optional
        Edge line width. Default is 6.0.
    node_color : str, list, or array, optional
        Node colors. "auto" uses cluster colors if LRG available,
        otherwise uniform. Can also pass explicit colors.
    node_size : float, optional
        Node marker size. Default is 6.0.
    colorbar : bool, optional
        Show colorbar for edge weights. Default is True.
    title : str, optional
        Plot title. Auto-generated if None.

    Returns
    -------
    view : nilearn.plotting.html_connectome.ConnectomeView
        Interactive HTML view. Can be displayed in Jupyter or saved.

    Examples
    --------
    >>> view = view_brain_connectome("Pat_02", "rest_pre", "beta")
    >>> view  # Display in Jupyter
    >>> view.save_as_html("connectome.html")  # Save to file
    """
    from nilearn.plotting import view_connectome

    # Load metadata and coordinates (transform to MNI-like space for brain viz)
    metadata = load_spatial_metadata(patient, dataset_root)
    coords = prepare_spatial_coordinates(metadata, scale="mm", center=False, to_mni=True)

    # Load FC matrix
    from lrg_eegfc.workflow.fc import load_fc_matrix as _load_fc
    fc_kw = {}
    if fc_cache_root is not None:
        fc_kw["cache_root"] = fc_cache_root
    fc_matrix = _load_fc(patient, phase, band, fc_method, **fc_kw)

    if fc_matrix is None:
        raise FileNotFoundError(
            f"FC matrix not found for {patient} {phase} {band} ({fc_method})"
        )

    # Handle node colors
    if node_color == "auto" and lrg_cache_root is not None:
        try:
            from lrg_eegfc.workflow.lrg import load_lrg_result
            lrg_result = load_lrg_result(
                patient, phase, band, fc_method, lrg_cache_root
            )
            if lrg_result is not None:
                cluster_labels = fcluster(
                    lrg_result.linkage_matrix,
                    t=lrg_result.optimal_threshold,
                    criterion="distance",
                )
                # Convert cluster labels to colors
                from matplotlib import cm
                n_clusters = len(np.unique(cluster_labels))
                cmap = cm.get_cmap("tab20" if n_clusters <= 20 else "nipy_spectral")
                unique = np.unique(cluster_labels)
                color_map = {c: cmap(i / max(len(unique) - 1, 1)) for i, c in enumerate(unique)}
                node_color = [color_map[c] for c in cluster_labels]
        except Exception:
            node_color = "auto"

    # Generate title
    if title is None:
        title = f"{patient} {phase} {band} ({fc_method})"

    # Create interactive view
    view = view_connectome(
        adjacency_matrix=fc_matrix,
        node_coords=coords,
        edge_threshold=edge_threshold,
        edge_cmap=edge_cmap,
        symmetric_cmap=symmetric_cmap,
        linewidth=linewidth,
        node_color=node_color,
        node_size=node_size,
        colorbar=colorbar,
        title=title,
    )

    return view


def plot_brain_connectome(
    patient: str,
    phase: str,
    band: str,
    fc_method: str = "imcoh_abs",
    dataset_root: Path = SEEG_DATAPATH,
    fc_cache_root: Optional[Path] = None,
    lrg_cache_root: Optional[Path] = None,
    edge_threshold: Optional[float] = None,
    edge_cmap: str = "bwr",
    node_color: Optional[Union[str, List, np.ndarray]] = "auto",
    node_size: float = 50.0,
    display_mode: str = "lzry",
    colorbar: bool = True,
    title: Optional[str] = None,
    output_path: Optional[Path] = None,
    figsize: Tuple[int, int] = (12, 4),
) -> Path:
    """Static brain connectome visualization using nilearn.

    Creates a multi-view glass brain plot with network connections.

    Parameters
    ----------
    patient : str
        Patient identifier.
    phase : str
        Recording phase.
    band : str
        Frequency band.
    fc_method : str, optional
        FC method. Default ``"imcoh_abs"`` (current era).
    dataset_root : Path, optional
        Root directory for patient data.
    fc_cache_root : Path, optional
        Root directory for FC cache.
    lrg_cache_root : Path, optional
        Root directory for LRG cache (for cluster colors).
    edge_threshold : float, optional
        Minimum edge weight to display. None shows all edges.
    edge_cmap : str, optional
        Colormap for edges.
    node_color : str, list, or array, optional
        Node colors. "auto" uses cluster colors if available.
    node_size : float, optional
        Node marker size.
    display_mode : str, optional
        Viewing planes: 'l'=left, 'r'=right, 'z'=axial, 'y'=coronal,
        'x'=sagittal. Default "lzry" shows 4 views.
    colorbar : bool, optional
        Show colorbar.
    title : str, optional
        Plot title.
    output_path : Path, optional
        Output file path. Auto-generated if None.
    figsize : tuple, optional
        Figure size in inches.

    Returns
    -------
    output_path : Path
        Path to saved figure.
    """
    from nilearn.plotting import plot_connectome

    # Load data (transform to MNI-like space for brain viz)
    metadata = load_spatial_metadata(patient, dataset_root)
    coords = prepare_spatial_coordinates(metadata, scale="mm", center=False, to_mni=True)

    # Load FC matrix
    from lrg_eegfc.workflow.fc import load_fc_matrix as _load_fc
    fc_kw = {}
    if fc_cache_root is not None:
        fc_kw["cache_root"] = fc_cache_root
    fc_matrix = _load_fc(patient, phase, band, fc_method, **fc_kw)

    if fc_matrix is None:
        raise FileNotFoundError(
            f"FC matrix not found for {patient} {phase} {band} ({fc_method})"
        )

    # Handle node colors
    if node_color == "auto" and lrg_cache_root is not None:
        try:
            from lrg_eegfc.workflow.lrg import load_lrg_result
            lrg_result = load_lrg_result(
                patient, phase, band, fc_method, lrg_cache_root
            )
            if lrg_result is not None:
                cluster_labels = fcluster(
                    lrg_result.linkage_matrix,
                    t=lrg_result.optimal_threshold,
                    criterion="distance",
                )
                from matplotlib import cm
                n_clusters = len(np.unique(cluster_labels))
                cmap = cm.get_cmap("tab20" if n_clusters <= 20 else "nipy_spectral")
                unique = np.unique(cluster_labels)
                color_map = {c: cmap(i / max(len(unique) - 1, 1)) for i, c in enumerate(unique)}
                node_color = [color_map[c] for c in cluster_labels]
        except Exception:
            node_color = "auto"

    # Generate title
    if title is None:
        title = f"{patient} {phase} {band} ({fc_method})"

    # Create figure
    fig = plt.figure(figsize=figsize)

    # Plot connectome
    plot_connectome(
        adjacency_matrix=fc_matrix,
        node_coords=coords,
        edge_threshold=edge_threshold,
        edge_cmap=edge_cmap,
        node_color=node_color,
        node_size=node_size,
        display_mode=display_mode,
        colorbar=colorbar,
        title=title,
        figure=fig,
    )

    # Save figure
    if output_path is None:
        output_dir = FIGURES_ROOT / "spatial" / patient
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{band}_{phase}_{fc_method}_brain_connectome.png"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return output_path


def plot_brain_connectome_at_n(
    axes,
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    n_clusters: int,
    *,
    dataset_root: Path = SEEG_DATAPATH,
    edge_top_pct: float = 10.0,
    display_mode: str = "z",
    palette: str = "tab10",
    node_size_range: Tuple[float, float] = (20.0, 100.0),
    title: Optional[str] = None,
    annotate_enrichment: bool = False,
    zoom_to_nodes: bool = True,
    zoom_pad_frac: float = 0.10,
):
    """Glass-brain projection coloured by community at a fixed *n_clusters* cut.

    Wraps :func:`nilearn.plotting.plot_connectome` so it can be embedded in a
    pre-existing matplotlib axes (typically a subplot in a grid).  The cut is
    obtained by ``scipy.cluster.hierarchy.fcluster(linkage, n_clusters,
    criterion='maxclust')`` rather than the LRG optimal threshold, which lets
    us directly compare the same scale across MSC and |ImCoh|.

    Parameters
    ----------
    axes : matplotlib.axes.Axes
        Target axes (single-view ``display_mode`` only).
    patient, phase, band, fc_method : str
        Standard FC routing.
    n_clusters : int
        Number of communities to cut the dendrogram at.
    dataset_root : Path
        Patient data root (for MNI metadata).
    edge_top_pct : float
        Keep only the top ``edge_top_pct`` % strongest edges (by weight).
    display_mode : str
        Single-view nilearn mode ("z" axial recommended for sEEG).
    palette : str
        Matplotlib colormap used for community colours.
    node_size_range : (low, high)
        Min/max marker size; nodes are scaled by row-sum strength.
    title : str, optional
        Panel title.
    annotate_enrichment : bool
        If True, write the same-shaft enrichment ratio on the axes.

    Returns
    -------
    dict with keys: ``community_labels``, ``enrichment``, ``n_nodes``,
    ``coords`` (subset used).
    """
    from nilearn.plotting import plot_connectome
    from lrg_eegfc.workflow.fc import load_fc_matrix as _load_fc
    from lrg_eegfc.workflow.lrg import load_lrg_result as _load_lrg
    from lrg_eegfc.utils.probe import (
        extract_probe_labels,
        compute_community_probe_enrichment,
    )

    fc_matrix = _load_fc(patient, phase, band, fc_method)
    if fc_matrix is None:
        raise FileNotFoundError(
            f"FC matrix missing: {patient} {phase} {band} ({fc_method})"
        )

    lrg = _load_lrg(patient, phase, band, fc_method)
    if lrg is None:
        raise FileNotFoundError(
            f"LRG result missing: {patient} {phase} {band} ({fc_method})"
        )

    metadata = load_spatial_metadata(patient, dataset_root)
    coords = prepare_spatial_coordinates(metadata, scale="mm", center=False, to_mni=True)
    channel_labels = list(metadata["label"])

    # If LRG used only the giant component, restrict everything to those nodes
    if lrg.n_nodes < fc_matrix.shape[0]:
        import networkx as nx
        from lrgsglib.utils import get_giant_component as _get_gc
        Gfull = nx.from_numpy_array(fc_matrix)
        Gg = _get_gc(Gfull)
        gc_nodes = sorted(Gg.nodes())
        fc_matrix = fc_matrix[np.ix_(gc_nodes, gc_nodes)]
        coords = coords[gc_nodes]
        channel_labels = [channel_labels[i] for i in gc_nodes]

    # Community partition at fixed n
    community_labels = fcluster(lrg.linkage_matrix, t=n_clusters, criterion="maxclust")

    # Node colours: rank communities by descending size and assign palette
    # in that order. This makes the largest community use palette[0] at
    # every cut, so colours are stable across n_clusters values within a
    # column.
    cmap = plt.get_cmap(palette)
    unique, counts = np.unique(community_labels, return_counts=True)
    rank_order = np.argsort(-counts, kind="stable")
    ranked_communities = [int(unique[i]) for i in rank_order]
    color_map = {c: cmap(i % cmap.N) for i, c in enumerate(ranked_communities)}
    node_colors = [color_map[int(c)] for c in community_labels]

    # Node sizes from row strength
    A = np.abs(fc_matrix).copy()
    np.fill_diagonal(A, 0.0)
    strengths = A.sum(axis=1)
    s_min, s_max = float(strengths.min()), float(strengths.max())
    if s_max > s_min:
        s_norm = (strengths - s_min) / (s_max - s_min)
    else:
        s_norm = np.zeros_like(strengths)
    node_sizes = node_size_range[0] + s_norm * (node_size_range[1] - node_size_range[0])

    # Top-k edge threshold
    triu = np.triu_indices_from(A, k=1)
    weights = A[triu]
    if weights.size == 0:
        edge_threshold_val = 0.0
    else:
        edge_threshold_val = float(np.percentile(weights, 100.0 - edge_top_pct))

    disp = plot_connectome(
        adjacency_matrix=A,
        node_coords=coords,
        edge_threshold=edge_threshold_val,
        edge_cmap="Greys",
        node_color=node_colors,
        node_size=node_sizes,
        display_mode=display_mode,
        colorbar=False,
        axes=axes,
        annotate=False,
    )

    if zoom_to_nodes:
        proj_ij = {"x": (1, 2), "y": (0, 2), "z": (0, 1)}.get(display_mode)
        if proj_ij is not None:
            pos2d = np.asarray(coords)[:, list(proj_ij)]
            xmin, ymin = pos2d.min(axis=0)
            xmax, ymax = pos2d.max(axis=0)
            xspan = max(xmax - xmin, 1.0)
            yspan = max(ymax - ymin, 1.0)
            pad = zoom_pad_frac * max(xspan, yspan)
            try:
                overlay_ax = list(disp.axes.values())[0].ax
                overlay_ax.set_xlim(xmin - pad, xmax + pad)
                overlay_ax.set_ylim(ymin - pad, ymax + pad)
            except (AttributeError, IndexError):
                pass

    # Same-shaft enrichment (returned for the companion .md, not drawn)
    probes = extract_probe_labels(channel_labels)
    enrichment = compute_community_probe_enrichment(community_labels, probes)

    if title is not None:
        axes.set_title(title, fontsize=9)

    if annotate_enrichment:
        axes.text(
            0.02, 0.02,
            f"shaft enrichment = {enrichment:.2f}×",
            transform=axes.transAxes,
            fontsize=7,
            color="black",
            ha="left", va="bottom",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1.5),
        )

    return {
        "community_labels": community_labels,
        "enrichment": enrichment,
        "n_nodes": len(community_labels),
        "coords": coords,
    }
