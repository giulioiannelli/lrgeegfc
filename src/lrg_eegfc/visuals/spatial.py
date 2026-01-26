"""3D Spatial Brain Network Visualization.

This module provides visualization functions for displaying SEEG networks
with nodes positioned at their real 3D anatomical coordinates from electrode
implant data.

Functions support both interactive Plotly visualizations and static
Matplotlib fallbacks for publication-quality figures.
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

__all__ = [
    "load_spatial_metadata",
    "prepare_spatial_coordinates",
    "build_edge_traces",
    "plot_spatial_network_3d",
    "plot_spatial_network_3d_mpl",
    "plot_spatial_clusters_comparison",
    "view_brain_connectome",
    "plot_brain_connectome",
]


def _normalize_label(label: str) -> str:
    """Normalize electrode label for matching.

    Handles differences like 'A1' vs 'A 1' vs 'A 1,G2' vs "G'1" vs "G' 1".
    """
    import re

    # Remove quotes and strip whitespace
    s = str(label).strip().strip('"').strip("'")
    # Remove suffix like ',G2'
    if "," in s:
        s = s.split(",")[0]
    # Remove ALL spaces (A 1 -> A1, G' 1 -> G'1)
    s = s.replace(" ", "")
    # Normalize prime/accent characters
    s = s.replace("ì", "'")  # Replace accented i with prime
    return s.strip()


def load_spatial_metadata(
    patient: str,
    dataset_root: Path,
) -> pd.DataFrame:
    """Load channel metadata with spatial coordinates.

    Handles label mismatches between channel_labels.csv and Implant_pat_*.csv
    by normalizing labels before merging.

    Parameters
    ----------
    patient : str
        Patient identifier (e.g., "Pat_02").
    dataset_root : Path
        Root directory containing patient folders.

    Returns
    -------
    metadata : pd.DataFrame
        DataFrame with columns: label, x, y, z, and optionally atlas info.

    Raises
    ------
    FileNotFoundError
        If required files are not found.
    """
    patient_path = Path(dataset_root) / patient

    # Find implant file
    try:
        patnum = int(patient.split("_")[-1])
    except ValueError:
        patnum = 0
    implant_csv = patient_path / f"Implant_pat_{patnum:02d}.csv"

    channel_labels_csv = patient_path / "channel_labels.csv"

    if not implant_csv.exists():
        raise FileNotFoundError(f"Implant file not found: {implant_csv}")
    if not channel_labels_csv.exists():
        raise FileNotFoundError(f"Channel labels not found: {channel_labels_csv}")

    # Load files
    # Implant file may have quoted Desikan-Killany field with internal commas
    # Only read the columns we need
    implant = pd.read_csv(implant_csv, usecols=["label", "x", "y", "z", "Desikan-Killany"])

    # Channel labels file may or may not have a header
    # Also may have extra columns (e.g., "F 1,G2" format)
    labels = pd.read_csv(channel_labels_csv)
    if "label" not in labels.columns:
        # No header - read again, use first column as label
        labels = pd.read_csv(channel_labels_csv, header=None)
        labels = labels.rename(columns={labels.columns[0]: "label"})
        # Keep only the label column
        labels = labels[["label"]]

    # Normalize labels for matching
    implant["_norm_label"] = implant["label"].apply(_normalize_label)
    labels["_norm_label"] = labels["label"].apply(_normalize_label)

    # Merge on normalized labels
    metadata = labels.merge(
        implant[["_norm_label", "x", "y", "z", "Desikan-Killany"]],
        on="_norm_label",
        how="left",
    )

    # Clean up coordinate columns (handle comma decimals)
    for col in ("x", "y", "z"):
        if col in metadata.columns:
            metadata[col] = (
                metadata[col]
                .astype(str)
                .str.replace(",", ".", regex=False)
                .replace("nan", np.nan)
                .astype(float)
            )

    # Drop helper column
    metadata = metadata.drop(columns=["_norm_label"])

    # Check how many matched
    n_matched = metadata["x"].notna().sum()
    n_total = len(metadata)
    if n_matched == 0:
        raise ValueError(
            f"No coordinates matched for {patient}. "
            f"Check label formats in implant and channel_labels files."
        )
    if n_matched < n_total:
        import warnings
        warnings.warn(
            f"Only {n_matched}/{n_total} channels matched coordinates for {patient}."
        )

    return metadata


def prepare_spatial_coordinates(
    metadata: pd.DataFrame,
    scale: str = "mm",
    center: bool = True,
    to_mni: bool = False,
) -> np.ndarray:
    """Extract and transform electrode coordinates from metadata.

    Parameters
    ----------
    metadata : pd.DataFrame
        Channel metadata with 'x', 'y', 'z' columns (in micrometers).
    scale : str, optional
        Output scale: "mm" (divide by 1000) or "um" (raw micrometers).
        Default is "mm".
    center : bool, optional
        If True, center coordinates around origin. Default is True.
        Ignored if to_mni=True.
    to_mni : bool, optional
        If True, apply estimated transform to MNI-like space for brain
        visualization. This shifts coordinates to match typical brain
        anatomy based on atlas labels. Default is False.

    Returns
    -------
    coords : np.ndarray
        Array of shape (n_channels, 3) with coordinates.

    Raises
    ------
    ValueError
        If required coordinate columns are missing or contain NaN values.
    """
    required_cols = ["x", "y", "z"]
    missing = [c for c in required_cols if c not in metadata.columns]
    if missing:
        raise ValueError(f"Missing coordinate columns: {missing}")

    coords = metadata[required_cols].to_numpy(dtype=float)

    if np.isnan(coords).any():
        n_nan = np.isnan(coords).any(axis=1).sum()
        raise ValueError(
            f"Found {n_nan} channels with missing coordinate values. "
            "Ensure implant metadata is complete."
        )

    if scale == "mm":
        coords = coords / 1000.0
    elif scale != "um":
        raise ValueError(f"Unknown scale '{scale}'. Use 'mm' or 'um'.")

    if to_mni:
        coords = _estimate_mni_transform(coords, metadata)
    elif center:
        coords = coords - coords.mean(axis=0)

    return coords


def _estimate_mni_transform(coords: np.ndarray, metadata: pd.DataFrame) -> np.ndarray:
    """Estimate transform from native space to MNI-like coordinates.

    Uses Desikan-Killiany atlas labels and their known MNI coordinates to
    estimate a translation that aligns native coordinates to MNI space.

    This is an approximation - for accurate MNI coordinates, use proper
    registration with patient MRI data.

    Parameters
    ----------
    coords : np.ndarray
        Coordinates in native space (mm).
    metadata : pd.DataFrame
        Metadata with Desikan-Killany atlas labels.

    Returns
    -------
    mni_coords : np.ndarray
        Coordinates in approximate MNI space.
    """
    # Known MNI coordinates for Desikan-Killiany regions (approximate centroids)
    # Values from standard MNI atlases
    MNI_REGION_COORDS = {
        # Left hemisphere cortical regions
        "ctx-lh-fusiform": (-35, -50, -18),
        "ctx-lh-parahippocampal": (-25, -25, -20),
        "ctx-lh-middletemporal": (-55, -25, -10),
        "ctx-lh-inferiortemporal": (-50, -30, -25),
        "ctx-lh-superiortemporal": (-55, -15, 0),
        "ctx-lh-transversetemporal": (-45, -20, 10),
        "ctx-lh-temporalpole": (-35, 15, -30),
        "ctx-lh-bankssts": (-55, -45, 5),
        "ctx-lh-entorhinal": (-25, -10, -30),
        "ctx-lh-insula": (-38, 0, 5),
        "ctx-lh-lateralorbitofrontal": (-30, 35, -15),
        "ctx-lh-medialorbitofrontal": (-8, 45, -15),
        "ctx-lh-parsorbitalis": (-40, 40, -10),
        "ctx-lh-parstriangularis": (-48, 30, 5),
        "ctx-lh-parsopercularis": (-50, 15, 10),
        "ctx-lh-rostralmiddlefrontal": (-35, 45, 20),
        "ctx-lh-caudalmiddlefrontal": (-40, 15, 45),
        "ctx-lh-superiorfrontal": (-15, 35, 45),
        "ctx-lh-frontalpole": (-10, 65, -5),
        "ctx-lh-precentral": (-40, -10, 55),
        "ctx-lh-postcentral": (-45, -25, 55),
        "ctx-lh-paracentral": (-10, -30, 65),
        "ctx-lh-supramarginal": (-55, -40, 35),
        "ctx-lh-inferiorparietal": (-45, -60, 40),
        "ctx-lh-superiorparietal": (-25, -60, 55),
        "ctx-lh-precuneus": (-10, -60, 40),
        "ctx-lh-cuneus": (-10, -85, 20),
        "ctx-lh-lingual": (-15, -70, -5),
        "ctx-lh-pericalcarine": (-10, -85, 5),
        "ctx-lh-lateraloccipital": (-40, -80, 10),
        "ctx-lh-rostralanteriorcingulate": (-8, 35, 10),
        "ctx-lh-caudalanteriorcingulate": (-8, 15, 30),
        "ctx-lh-posteriorcingulate": (-8, -40, 30),
        "ctx-lh-isthmuscingulate": (-10, -45, 10),
        "ctx-lh-hippocampus": (-28, -20, -15),
        "ctx-lh-amygdala": (-25, -5, -20),
        # Right hemisphere (mirror of left)
        "ctx-rh-fusiform": (35, -50, -18),
        "ctx-rh-parahippocampal": (25, -25, -20),
        "ctx-rh-middletemporal": (55, -25, -10),
        "ctx-rh-inferiortemporal": (50, -30, -25),
        "ctx-rh-superiortemporal": (55, -15, 0),
        "ctx-rh-transversetemporal": (45, -20, 10),
        "ctx-rh-temporalpole": (35, 15, -30),
        "ctx-rh-bankssts": (55, -45, 5),
        "ctx-rh-entorhinal": (25, -10, -30),
        "ctx-rh-insula": (38, 0, 5),
        "ctx-rh-lateralorbitofrontal": (30, 35, -15),
        "ctx-rh-medialorbitofrontal": (8, 45, -15),
        "ctx-rh-parsorbitalis": (40, 40, -10),
        "ctx-rh-parstriangularis": (48, 30, 5),
        "ctx-rh-parsopercularis": (50, 15, 10),
        "ctx-rh-rostralmiddlefrontal": (35, 45, 20),
        "ctx-rh-caudalmiddlefrontal": (40, 15, 45),
        "ctx-rh-superiorfrontal": (15, 35, 45),
        "ctx-rh-frontalpole": (10, 65, -5),
        "ctx-rh-precentral": (40, -10, 55),
        "ctx-rh-postcentral": (45, -25, 55),
        "ctx-rh-paracentral": (10, -30, 65),
        "ctx-rh-supramarginal": (55, -40, 35),
        "ctx-rh-inferiorparietal": (45, -60, 40),
        "ctx-rh-superiorparietal": (25, -60, 55),
        "ctx-rh-precuneus": (10, -60, 40),
        "ctx-rh-rostralanteriorcingulate": (8, 35, 10),
        "ctx-rh-caudalanteriorcingulate": (8, 15, 30),
        "ctx-rh-posteriorcingulate": (8, -40, 30),
        "ctx-rh-isthmuscingulate": (10, -45, 10),
        "ctx-rh-hippocampus": (28, -20, -15),
        "ctx-rh-amygdala": (25, -5, -20),
        # Subcortical
        "Left-Hippocampus": (-28, -20, -15),
        "Right-Hippocampus": (28, -20, -15),
        "Left-Amygdala": (-25, -5, -20),
        "Right-Amygdala": (25, -5, -20),
        "Left-Thalamus": (-12, -18, 8),
        "Right-Thalamus": (12, -18, 8),
        "Left-Caudate": (-12, 12, 10),
        "Right-Caudate": (12, 12, 10),
        "Left-Putamen": (-25, 5, 2),
        "Right-Putamen": (25, 5, 2),
        "Left-Accumbens-area": (-10, 10, -8),
        "Right-Accumbens-area": (10, 10, -8),
        "Left-Pallidum": (-18, 0, 0),
        "Right-Pallidum": (18, 0, 0),
        "Cerebellum-Cortex": (0, -55, -35),
        "Brain-Stem": (0, -30, -30),
    }

    def get_primary_region(dk_str):
        """Extract primary region from Desikan-Killiany string."""
        if pd.isna(dk_str):
            return None
        parts = str(dk_str).strip().split(",")
        if parts:
            return parts[0].strip().strip('"').strip()
        return None

    # Detect coordinate system orientation using atlas labels
    # In MNI: inferior structures (cerebellum, temporal pole) have negative Z
    #         superior structures (superior frontal) have positive Z
    # Some patients have inverted Z (more negative = more inferior)
    # Others have standard Z (more positive = more superior)

    z_inverted = False
    if "Desikan-Killany" in metadata.columns:
        # Find inferior regions (cerebellum, temporal pole, fusiform)
        inferior_z = []
        superior_z = []
        for idx, row in metadata.iterrows():
            region = get_primary_region(row.get("Desikan-Killany"))
            if region:
                region_lower = region.lower()
                if any(r in region_lower for r in ["cerebellum", "temporalpole", "fusiform", "entorhinal"]):
                    inferior_z.append(coords[idx, 2])
                elif any(r in region_lower for r in ["superiorfrontal", "precentral", "postcentral", "superiorparietal"]):
                    superior_z.append(coords[idx, 2])

        if inferior_z and superior_z:
            # If inferior regions have MORE negative Z than superior, it's inverted
            z_inverted = np.mean(inferior_z) < np.mean(superior_z)
        elif inferior_z:
            # If only inferior regions exist and Z is very negative, likely inverted
            z_inverted = np.mean(inferior_z) < -20
        elif superior_z:
            # If only superior regions exist and Z is very positive, likely standard
            z_inverted = np.mean(superior_z) < 0

    # Work with potentially flipped coordinates
    working_coords = coords.copy()

    # Collect matched region coordinates for calibration
    native_coords_matched = []
    mni_coords_matched = []

    if "Desikan-Killany" in metadata.columns:
        for idx, row in metadata.iterrows():
            region = get_primary_region(row.get("Desikan-Killany"))
            if region and region in MNI_REGION_COORDS:
                native_coords_matched.append(working_coords[idx])
                mni_coords_matched.append(MNI_REGION_COORDS[region])

    # Calculate translation offset
    if len(native_coords_matched) >= 3:
        native_arr = np.array(native_coords_matched)
        mni_arr = np.array(mni_coords_matched)
        # Use median to reduce effect of outliers
        translation = np.median(mni_arr - native_arr, axis=0)
    else:
        # Fallback: use hemisphere-based heuristic
        if "Desikan-Killany" in metadata.columns:
            atlas = metadata["Desikan-Killany"].dropna().astype(str)
            lh_count = sum("lh" in a or "Left" in a for a in atlas)
            rh_count = sum("rh" in a or "Right" in a for a in atlas)
            is_bilateral = lh_count > 10 and rh_count > 10
            is_predominantly_left = lh_count > rh_count
        else:
            is_bilateral = False
            is_predominantly_left = working_coords[:, 0].mean() < 0

        centroid = working_coords.mean(axis=0)

        if is_bilateral:
            target_centroid = np.array([0.0, -20.0, 10.0])
        elif is_predominantly_left:
            target_centroid = np.array([-35.0, -20.0, 10.0])
        else:
            target_centroid = np.array([35.0, -20.0, 10.0])

        translation = target_centroid - centroid

    # Apply translation
    mni_coords = working_coords + translation

    # Clamp to reasonable MNI bounds
    mni_coords[:, 0] = np.clip(mni_coords[:, 0], -75, 75)  # X: left-right
    mni_coords[:, 1] = np.clip(mni_coords[:, 1], -110, 75)  # Y: posterior-anterior
    mni_coords[:, 2] = np.clip(mni_coords[:, 2], -60, 85)  # Z: inferior-superior

    return mni_coords


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
    fc_method: str = "msc",
    dataset_root: Path = Path("data/stereoeeg_patients"),
    lrg_cache_root: Path = Path("data/lrg_cache"),
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
        Recording phase (e.g., "rsPre", "rsPost").
    band : str
        Frequency band (e.g., "beta", "alpha").
    fc_method : str, optional
        Functional connectivity method: "msc" or "corr". Default is "msc".
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
        if fc_method == "msc":
            from lrg_eegfc.workflow.msc import load_msc_matrix

            fc_root = fc_cache_root or Path("data/msc_cache")
            fc_matrix = load_msc_matrix(patient, phase, band, cache_root=fc_root)
        else:
            from lrg_eegfc.workflow.corr import load_corr_matrix

            fc_root = fc_cache_root or Path("data/corr_cache")
            fc_matrix = load_corr_matrix(
                patient, phase, band, cache_root=fc_root, filter_type="abs"
            )

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
    fc_method: str = "msc",
    dataset_root: Path = Path("data/stereoeeg_patients"),
    lrg_cache_root: Path = Path("data/lrg_cache"),
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
        FC method: "msc" or "corr". Default is "msc".
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
        if fc_method == "msc":
            from lrg_eegfc.workflow.msc import load_msc_matrix

            fc_root = fc_cache_root or Path("data/msc_cache")
            fc_matrix = load_msc_matrix(patient, phase, band, cache_root=fc_root)
        else:
            from lrg_eegfc.workflow.corr import load_corr_matrix

            fc_root = fc_cache_root or Path("data/corr_cache")
            fc_matrix = load_corr_matrix(
                patient, phase, band, cache_root=fc_root, filter_type="abs"
            )

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
        output_dir = Path("data/figures/spatial") / patient
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
    fc_method: str = "msc",
    dataset_root: Path = Path("data/stereoeeg_patients"),
    lrg_cache_root: Path = Path("data/lrg_cache"),
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
        First phase (e.g., "rsPre").
    phase_b : str
        Second phase (e.g., "rsPost").
    band : str
        Frequency band.
    fc_method : str, optional
        FC method. Default is "msc".
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


def view_brain_connectome(
    patient: str,
    phase: str,
    band: str,
    fc_method: str = "msc",
    dataset_root: Path = Path("data/stereoeeg_patients"),
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
        Recording phase (e.g., "rsPre", "rsPost").
    band : str
        Frequency band (e.g., "beta", "alpha").
    fc_method : str, optional
        FC method: "msc" or "corr". Default is "msc".
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
    >>> view = view_brain_connectome("Pat_02", "rsPre", "beta")
    >>> view  # Display in Jupyter
    >>> view.save_as_html("connectome.html")  # Save to file
    """
    from nilearn.plotting import view_connectome

    # Load metadata and coordinates (transform to MNI-like space for brain viz)
    metadata = load_spatial_metadata(patient, dataset_root)
    coords = prepare_spatial_coordinates(metadata, scale="mm", center=False, to_mni=True)

    # Load FC matrix
    if fc_method == "msc":
        from lrg_eegfc.workflow.msc import load_msc_matrix
        fc_root = fc_cache_root or Path("data/msc_cache")
        fc_matrix = load_msc_matrix(patient, phase, band, cache_root=fc_root)
    else:
        from lrg_eegfc.workflow.corr import load_corr_matrix
        fc_root = fc_cache_root or Path("data/corr_cache")
        fc_matrix = load_corr_matrix(
            patient, phase, band, cache_root=fc_root, filter_type="abs"
        )

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
    fc_method: str = "msc",
    dataset_root: Path = Path("data/stereoeeg_patients"),
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
        FC method. Default is "msc".
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
    if fc_method == "msc":
        from lrg_eegfc.workflow.msc import load_msc_matrix
        fc_root = fc_cache_root or Path("data/msc_cache")
        fc_matrix = load_msc_matrix(patient, phase, band, cache_root=fc_root)
    else:
        from lrg_eegfc.workflow.corr import load_corr_matrix
        fc_root = fc_cache_root or Path("data/corr_cache")
        fc_matrix = load_corr_matrix(
            patient, phase, band, cache_root=fc_root, filter_type="abs"
        )

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
        output_dir = Path("data/figures/spatial") / patient
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{band}_{phase}_{fc_method}_brain_connectome.png"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return output_path
