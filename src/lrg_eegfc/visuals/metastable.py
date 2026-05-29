"""Metastable nodes visualization using Sankey diagrams.

This module implements cluster evolution tracking and Sankey diagram generation
to visualize how network nodes transition between communities across different
hierarchical scales (tau values).

Two rendering backends share the same data-prep pipeline:

  * ``make_metastable_sankey_plotly`` — interactive HTML (Plotly)
  * ``make_metastable_sankey_mpl``    — static PDF (matplotlib alluvial
                                         with cubic-bezier ribbons)

Both consume ``SankeyFlowData`` produced by ``compute_sankey_flows``.

Based on FIGMNTGN04.ipynb Sankey diagram analysis; promoted to template
library 2026-05-26.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import plotly.graph_objects as go
import plotly.colors as colors
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform

__all__ = [
    "compute_clustering_across_tau",
    "analyze_node_trajectories",
    "track_cluster_changes",
    "compute_sankey_flows",
    "SankeyFlowData",
    "create_sankey_diagram",
    "make_metastable_sankey_plotly",
    "make_metastable_sankey_mpl",
    "DEFAULT_TAU_NCLUST_PAIRS",
]


DEFAULT_TAU_NCLUST_PAIRS = (
    (0.1, 15),
    (0.5, 8),
    (1.0, 5),
    (2.0, 3),
    (5.0, 2),
)


def compute_clustering_across_tau(
    linkage_matrix: np.ndarray,
    tau_values: np.ndarray,
    Gcc,
    method: str = "ward",
    scaling_factor: float = 0.99,
    n_clusters_list: list = None,
) -> Tuple[Dict, Dict]:
    """Compute clustering for different tau values.

    Parameters
    ----------
    linkage_matrix : np.ndarray
        Hierarchical linkage matrix (unused, kept for backwards compatibility)
    tau_values : np.ndarray
        Array of tau (diffusion time) values to test
    Gcc : networkx.Graph
        Network graph
    method : str
        Linkage method (default: "ward")
    scaling_factor : float
        For optimal threshold computation (used if n_clusters_list is None)
    n_clusters_list : list, optional
        Target number of clusters for each tau value. If provided, uses
        maxclust criterion for controlled hierarchical merging visualization.
        Should be decreasing (many clusters at small tau → few at large tau).

    Returns
    -------
    partitions : dict
        Mapping tau -> cluster labels array
    n_clusters : dict
        Mapping tau -> number of clusters
    """
    from lrgsglib.core import compute_laplacian_properties, compute_normalized_linkage, compute_optimal_threshold

    partitions = {}
    n_clusters = {}

    for i, tau in enumerate(tau_values):
        spectrum, L, rho, Trho, tau_used = compute_laplacian_properties(Gcc, tau=tau)
        dists = squareform(Trho)
        linkage, label_list, _ = compute_normalized_linkage(dists, Gcc, method=method)

        if n_clusters_list is not None:
            # Use fixed number of clusters for controlled visualization
            target_n = n_clusters_list[i]
            clusters = fcluster(linkage, target_n, criterion="maxclust")
        else:
            # Original behavior: compute optimal threshold at each tau
            clTh, *_ = compute_optimal_threshold(linkage, scaling_factor=scaling_factor)
            clusters = fcluster(linkage, clTh, criterion="distance")

        partitions[tau] = clusters
        n_clusters[tau] = len(np.unique(clusters))

    return partitions, n_clusters


def analyze_node_trajectories(
    partitions: Dict,
    tau_values: np.ndarray,
    node_labels: Optional[List[str]] = None,
) -> Tuple[Dict, Dict, Dict]:
    """Analyze node cluster membership trajectories.

    Parameters
    ----------
    partitions : dict
        Mapping tau -> cluster labels
    tau_values : ndarray
        Array of tau values
    node_labels : list, optional
        Node names/labels

    Returns
    -------
    node_trajectories : dict
        Mapping node_idx -> list of cluster assignments
    cluster_compositions : dict
        Mapping (tau, cluster) -> list of node indices
    trajectory_patterns : dict
        Unique trajectory patterns and their nodes
    """
    if node_labels is None:
        node_labels = [f"Node_{i}" for i in range(len(partitions[tau_values[0]]))]

    n_nodes = len(partitions[tau_values[0]])

    # Track each node's trajectory
    node_trajectories = {}
    for node_idx in range(n_nodes):
        trajectory = [partitions[tau][node_idx] for tau in tau_values]
        node_trajectories[node_idx] = trajectory

    # Track cluster compositions
    cluster_compositions = {}
    for tau in tau_values:
        clusters = partitions[tau]
        for cluster_id in np.unique(clusters):
            node_indices = np.where(clusters == cluster_id)[0]
            cluster_compositions[(tau, cluster_id)] = node_indices.tolist()

    # Find unique trajectory patterns
    trajectory_patterns = {}
    for node_idx, trajectory in node_trajectories.items():
        trajectory_tuple = tuple(trajectory)
        if trajectory_tuple not in trajectory_patterns:
            trajectory_patterns[trajectory_tuple] = []
        trajectory_patterns[trajectory_tuple].append(node_idx)

    return node_trajectories, cluster_compositions, trajectory_patterns


def track_cluster_changes(
    partitions: Dict,
    tau_values: np.ndarray,
    node_labels: Optional[List[str]] = None,
) -> Tuple[Dict, Dict, Dict]:
    """Track nodes that change clusters between consecutive tau steps.

    Parameters
    ----------
    partitions : dict
        Mapping tau -> cluster labels
    tau_values : ndarray
        Array of tau values
    node_labels : list, optional
        Node names

    Returns
    -------
    cluster_changes : dict
        Mapping (tau1, tau2) -> dict of node changes
    swapping_nodes : dict
        Nodes and their complete change history
    stability_analysis : dict
        Summary of node stability across tau values
    """
    if node_labels is None:
        node_labels = [f"Node_{i}" for i in range(len(partitions[tau_values[0]]))]

    cluster_changes = {}
    swapping_nodes = {node: [] for node in node_labels}

    # Track changes between consecutive tau steps
    for i in range(len(tau_values) - 1):
        tau1, tau2 = tau_values[i], tau_values[i + 1]
        clusters1, clusters2 = partitions[tau1], partitions[tau2]

        changes_at_step = {
            "tau_from": tau1,
            "tau_to": tau2,
            "nodes_changed": [],
            "nodes_stable": [],
            "change_details": {},
        }

        for node_idx, (c1, c2) in enumerate(zip(clusters1, clusters2)):
            node_name = node_labels[node_idx]

            if c1 != c2:
                change_info = {
                    "node": node_name,
                    "from_cluster": c1,
                    "to_cluster": c2,
                    "node_idx": node_idx,
                }
                changes_at_step["nodes_changed"].append(change_info)
                changes_at_step["change_details"][node_name] = change_info

                swapping_nodes[node_name].append(
                    {"tau_from": tau1, "tau_to": tau2, "from_cluster": c1, "to_cluster": c2}
                )
            else:
                changes_at_step["nodes_stable"].append(node_name)

        cluster_changes[(tau1, tau2)] = changes_at_step

    # Analyze stability
    stability_analysis = {}
    for node_name in node_labels:
        num_changes = len(swapping_nodes[node_name])
        stability_analysis[node_name] = {
            "total_changes": num_changes,
            "stability_score": 1 - (num_changes / (len(tau_values) - 1)),
            "change_history": swapping_nodes[node_name],
        }

    return cluster_changes, swapping_nodes, stability_analysis


@dataclass
class SankeyFlowData:
    """Shared data prep for both Sankey backends.

    Attributes
    ----------
    tau_values : np.ndarray
        The τ values, in column order (left-to-right).
    cluster_ids_per_tau : list[np.ndarray]
        ``cluster_ids_per_tau[i]`` is ``np.unique(partitions[tau_values[i]])``
        — the ordered cluster IDs that appear at column ``i``.
    sizes : list[np.ndarray]
        ``sizes[i][k]`` is the number of nodes in
        cluster ``cluster_ids_per_tau[i][k]`` at column ``i``.
    labels : list[str]
        Flat list of node labels, one per (τ, cluster) bar.
    hover_texts : list[str]
        Flat list of hover strings, one per (τ, cluster).
    sources, targets, values : list[int]
        Flat link lists (Plotly format).  Indices reference the flat
        ``labels`` list above.
    node_colors_rgb : list[tuple[float, float, float]]
        Per-(τ, cluster) RGB colours, normalised floats in [0, 1].
        Stable across renders.
    link_alpha : float
        Default alpha for ribbon colours.
    """

    tau_values: np.ndarray
    cluster_ids_per_tau: List[np.ndarray]
    sizes: List[np.ndarray]
    labels: List[str]
    hover_texts: List[str]
    sources: List[int]
    targets: List[int]
    values: List[int]
    node_colors_rgb: List[Tuple[float, float, float]]
    link_alpha: float = 0.6
    tau_to_offset: Dict[float, int] = field(default_factory=dict)
    cluster_compositions: Dict[Tuple[float, int], List[int]] = field(
        default_factory=dict
    )

    @property
    def n_columns(self) -> int:
        return len(self.tau_values)

    @property
    def node_count(self) -> int:
        return len(self.labels)


def _qualitative_palette_rgb() -> List[Tuple[float, float, float]]:
    """Plotly-Set3 + Pastel + Dark2 expanded to normalised RGB floats."""
    palette = (
        colors.qualitative.Set3
        + colors.qualitative.Pastel
        + colors.qualitative.Dark2
    )
    out: List[Tuple[float, float, float]] = []
    for c in palette:
        if c.startswith("rgb"):
            inner = c[c.index("(") + 1 : c.index(")")]
            r, g, b = (float(x) / 255.0 for x in inner.split(","))
        elif c.startswith("#"):
            hx = c.lstrip("#")
            r, g, b = (int(hx[j : j + 2], 16) / 255.0 for j in (0, 2, 4))
        else:
            r, g, b = (0.5, 0.5, 0.5)
        out.append((r, g, b))
    return out


def compute_sankey_flows(
    partitions: Dict,
    tau_values: np.ndarray,
    node_labels: Optional[Sequence[str]] = None,
    *,
    label_format: str = "τ={tau:.2f}\nC{cluster}",
    hover_max_nodes: int = 15,
    link_alpha: float = 0.6,
) -> SankeyFlowData:
    """Build backend-agnostic Sankey flow data.

    Shared by ``make_metastable_sankey_plotly`` and
    ``make_metastable_sankey_mpl``.

    Parameters
    ----------
    partitions : dict
        Mapping ``tau -> cluster labels`` (length N per τ).
    tau_values : array-like
        Ordered τ values; defines column order left-to-right.
    node_labels : sequence of str, optional
        Per-node names used in hover text.  Defaults to ``Node_i``.
    label_format : str
        Format string with ``{tau}`` and ``{cluster}`` placeholders for
        per-bar labels.  Default ``"τ={tau:.2f}\\nC{cluster}"``.
    hover_max_nodes : int
        Truncate hover lists longer than this with "+N more".
    link_alpha : float
        Default ribbon opacity (0–1).  Lower = airier; 0.6 matches the
        Plotly default.

    Returns
    -------
    SankeyFlowData
    """
    tau_values = np.asarray(tau_values)
    n_nodes = len(partitions[tau_values[0]])
    if node_labels is None:
        node_labels = [f"Node_{i}" for i in range(n_nodes)]
    node_labels = list(node_labels)

    _, cluster_compositions, _ = analyze_node_trajectories(
        partitions, tau_values, node_labels
    )

    palette = _qualitative_palette_rgb()

    cluster_ids_per_tau: List[np.ndarray] = []
    sizes: List[np.ndarray] = []
    labels: List[str] = []
    hover_texts: List[str] = []
    node_colors_rgb: List[Tuple[float, float, float]] = []
    tau_to_offset: Dict[float, int] = {}
    flat_index_of: Dict[Tuple[float, int], int] = {}
    node_count = 0

    for tau in tau_values:
        tau_key = float(tau)
        tau_to_offset[tau_key] = node_count
        clusters = partitions[tau]
        unique = np.unique(clusters)
        cluster_ids_per_tau.append(unique.copy())
        col_sizes = np.zeros(len(unique), dtype=int)
        for k, cl in enumerate(unique):
            members = cluster_compositions[(tau, cl)]
            col_sizes[k] = len(members)
            labels.append(label_format.format(tau=tau_key, cluster=cl))
            if len(members) <= hover_max_nodes:
                tail = "<br>".join(node_labels[i] for i in members)
            else:
                head = "<br>".join(node_labels[i] for i in members[:hover_max_nodes])
                tail = f"{head}<br>... and {len(members) - hover_max_nodes} more"
            hover_texts.append(
                f"τ = {tau_key:.3f}<br>Cluster {cl}<br>"
                f"{len(members)} nodes:<br>{tail}"
            )
            node_colors_rgb.append(palette[(int(cl) - 1) % len(palette)])
            flat_index_of[(tau_key, int(cl))] = node_count + k
        sizes.append(col_sizes)
        node_count += len(unique)

    sources: List[int] = []
    targets: List[int] = []
    values: List[int] = []
    for i in range(len(tau_values) - 1):
        tau1, tau2 = float(tau_values[i]), float(tau_values[i + 1])
        c1_arr, c2_arr = partitions[tau_values[i]], partitions[tau_values[i + 1]]
        trans: Dict[Tuple[int, int], int] = {}
        for c1, c2 in zip(c1_arr, c2_arr):
            key = (int(c1), int(c2))
            trans[key] = trans.get(key, 0) + 1
        for (c1, c2), count in trans.items():
            sources.append(flat_index_of[(tau1, c1)])
            targets.append(flat_index_of[(tau2, c2)])
            values.append(count)

    return SankeyFlowData(
        tau_values=tau_values,
        cluster_ids_per_tau=cluster_ids_per_tau,
        sizes=sizes,
        labels=labels,
        hover_texts=hover_texts,
        sources=sources,
        targets=targets,
        values=values,
        node_colors_rgb=node_colors_rgb,
        link_alpha=link_alpha,
        tau_to_offset=tau_to_offset,
        cluster_compositions={
            (float(t), int(c)): list(v) for (t, c), v in cluster_compositions.items()
        },
    )


def _rgb_to_rgba_string(rgb: Tuple[float, float, float], alpha: float) -> str:
    """``(r, g, b)`` floats in [0, 1] → ``"rgba(R,G,B,A)"`` Plotly string."""
    r, g, b = (int(round(c * 255)) for c in rgb)
    return f"rgba({r},{g},{b},{alpha:g})"


def make_metastable_sankey_plotly(
    flow: SankeyFlowData,
    *,
    title: str = "Cluster Evolution Across τ Values",
    width: int = 1400,
    height: int = 700,
    node_pad: int = 15,
    node_thickness: int = 20,
    font_size: int = 10,
) -> go.Figure:
    """Render a ``SankeyFlowData`` as a Plotly interactive Sankey.

    Use ``compute_sankey_flows`` to build the input.
    """
    node_colors = [
        _rgb_to_rgba_string(rgb, 1.0) for rgb in flow.node_colors_rgb
    ]
    link_colors = [
        _rgb_to_rgba_string(flow.node_colors_rgb[t], flow.link_alpha)
        for t in flow.targets
    ]
    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    pad=node_pad,
                    thickness=node_thickness,
                    line=dict(color="black", width=0.5),
                    label=flow.labels,
                    color=node_colors,
                    hovertemplate="%{customdata}<extra></extra>",
                    customdata=flow.hover_texts,
                ),
                link=dict(
                    source=flow.sources,
                    target=flow.targets,
                    value=flow.values,
                    color=link_colors,
                ),
            )
        ]
    )
    fig.update_layout(
        title_text=title, font_size=font_size, width=width, height=height
    )
    return fig


def create_sankey_diagram(
    partitions: Dict,
    tau_values: np.ndarray,
    node_labels: Optional[List[str]] = None,
    title: str = "Cluster Evolution Across τ Values",
    width: int = 1400,
    height: int = 700,
) -> go.Figure:
    """Back-compat shim around ``compute_sankey_flows`` +
    ``make_metastable_sankey_plotly``.  Kept for the legacy notebook
    and CLI call sites; new code should call the two-step path.
    """
    flow = compute_sankey_flows(partitions, np.asarray(tau_values), node_labels)
    return make_metastable_sankey_plotly(
        flow, title=title, width=width, height=height
    )


# ---------------------------------------------------------------------------
# matplotlib alluvial backend
# ---------------------------------------------------------------------------

def _ribbon_path(
    x0: float,
    y0_top: float,
    y0_bot: float,
    x1: float,
    y1_top: float,
    y1_bot: float,
):
    """Cubic-bezier alluvial ribbon between two cluster bars.

    Returns a ``matplotlib.path.Path`` that traces:
      A=(x0, y0_top) → smooth → D=(x1, y1_top) → straight → C=(x1, y1_bot)
      → smooth → B=(x0, y0_bot) → close.

    Control points are placed at the column midpoint, vertically aligned
    with the source/target ends — this gives a symmetric S-curve that
    matches Plotly's default ribbon shape.
    """
    from matplotlib.path import Path

    xm = 0.5 * (x0 + x1)
    verts = [
        (x0, y0_top),
        (xm, y0_top),
        (xm, y1_top),
        (x1, y1_top),
        (x1, y1_bot),
        (xm, y1_bot),
        (xm, y0_bot),
        (x0, y0_bot),
        (x0, y0_top),
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    return Path(verts, codes)


def make_metastable_sankey_mpl(
    flow: SankeyFlowData,
    *,
    ax=None,
    figsize: Tuple[float, float] = (12.0, 6.0),
    node_width: float = 0.06,
    node_gap_frac: float = 0.012,
    x_step: float = 1.0,
    ribbon_alpha: Optional[float] = None,
    node_edge_color: str = "black",
    node_edge_width: float = 0.5,
    show_bar_labels: bool = True,
    bar_label_fontsize: float = 7.0,
    show_column_labels: bool = True,
    column_label_format: str = r"$\tau$ = {tau:g}" + "\n({n} clusters)",
    column_label_fontsize: float = 10.0,
    column_label_pad: float = 0.05,
    title: Optional[str] = None,
    hide_axes: bool = True,
):
    """Render a ``SankeyFlowData`` as a static matplotlib alluvial.

    Polygon ribbons replace Plotly's path; cubic-bezier control points
    sit at the column midpoint for a symmetric S-curve.  Output is
    pure vector — never rasterise the ribbons.

    Parameters
    ----------
    flow : SankeyFlowData
        Output of ``compute_sankey_flows``.
    ax : matplotlib.axes.Axes, optional
        Axes to draw into.  If None, a new ``(fig, ax)`` is created at
        ``figsize`` and returned.
    figsize : (float, float)
        Used only when ``ax is None``.
    node_width : float
        Width of each cluster bar in x-axis units.  Default 0.06 of
        ``x_step``.  Increase for chunkier bars.
    node_gap_frac : float
        Vertical gap between adjacent cluster bars within a column,
        as a fraction of total node count.  Default 0.012 (≈1 row at
        N=80).
    x_step : float
        Horizontal spacing between columns.  Default 1.0.
    ribbon_alpha : float, optional
        Override ribbon opacity.  Default = ``flow.link_alpha``.
    show_bar_labels : bool
        Draw ``C{cluster} ({n})`` text on each bar.
    show_column_labels : bool
        Draw the per-column header (default ``"τ = X (k clusters)"``).
    column_label_format : str
        Format string with ``{tau}`` and ``{n}`` placeholders.
    title : str, optional
        Axis-level title.  No ``fig.suptitle``.
    hide_axes : bool
        Strip spines and ticks (alluvial figures rarely need them).

    Returns
    -------
    (fig, ax) : tuple
        ``fig`` is the parent figure of ``ax``.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, PathPatch

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    if ribbon_alpha is None:
        ribbon_alpha = flow.link_alpha

    n_cols = flow.n_columns
    n_total = int(flow.sizes[0].sum())
    gap = node_gap_frac * n_total

    # Cluster bar geometry — per (col, k) → (y_top, y_bottom) in node-count units.
    bar_bounds: List[List[Tuple[float, float]]] = []
    for col in range(n_cols):
        sizes = flow.sizes[col]
        k = len(sizes)
        total_gap = gap * max(0, k - 1)
        # Bars stack from y=0 at top down to y=-(N + total_gap) at bottom.
        col_bounds: List[Tuple[float, float]] = []
        y = 0.0
        for s in sizes:
            top = y
            bottom = y - s
            col_bounds.append((top, bottom))
            y = bottom - gap
        bar_bounds.append(col_bounds)

    # X coordinates for each column.
    x_centers = np.arange(n_cols, dtype=float) * x_step

    # Track per-bar remaining (top, bottom) consumption as we lay down ribbons.
    # Sources consume from the top down on the *right* side of the bar;
    # targets consume from the top down on the *left* side of the bar.
    src_cursor: List[List[float]] = [
        [bar_bounds[c][k][0] for k in range(len(flow.sizes[c]))]
        for c in range(n_cols)
    ]
    tgt_cursor: List[List[float]] = [
        [bar_bounds[c][k][0] for k in range(len(flow.sizes[c]))]
        for c in range(n_cols)
    ]

    # Map flat node index → (col, k_within_col).
    col_k_of: Dict[int, Tuple[int, int]] = {}
    flat = 0
    for col in range(n_cols):
        for k in range(len(flow.sizes[col])):
            col_k_of[flat] = (col, k)
            flat += 1

    # Order links by source column then by source bar order, then by target
    # bar order — gives a clean stacked-ribbon look without crossings within
    # a single source bar.
    link_idx = sorted(
        range(len(flow.values)),
        key=lambda ii: (
            col_k_of[flow.sources[ii]][0],
            col_k_of[flow.sources[ii]][1],
            col_k_of[flow.targets[ii]][1],
        ),
    )

    for ii in link_idx:
        src_flat = flow.sources[ii]
        tgt_flat = flow.targets[ii]
        v = flow.values[ii]
        col_s, k_s = col_k_of[src_flat]
        col_t, k_t = col_k_of[tgt_flat]

        y0_top = src_cursor[col_s][k_s]
        y0_bot = y0_top - v
        src_cursor[col_s][k_s] = y0_bot

        y1_top = tgt_cursor[col_t][k_t]
        y1_bot = y1_top - v
        tgt_cursor[col_t][k_t] = y1_bot

        x0 = x_centers[col_s] + 0.5 * node_width
        x1 = x_centers[col_t] - 0.5 * node_width

        rgb = flow.node_colors_rgb[tgt_flat]
        path = _ribbon_path(x0, y0_top, y0_bot, x1, y1_top, y1_bot)
        ax.add_patch(
            PathPatch(
                path,
                facecolor=rgb + (ribbon_alpha,),
                edgecolor="none",
                linewidth=0,
                zorder=1,
            )
        )

    # Draw cluster bars on top of ribbons.
    for col in range(n_cols):
        sizes = flow.sizes[col]
        for k, (top, bottom) in enumerate(bar_bounds[col]):
            flat = sum(len(flow.sizes[c]) for c in range(col)) + k
            rgb = flow.node_colors_rgb[flat]
            rect = Rectangle(
                (x_centers[col] - 0.5 * node_width, bottom),
                node_width,
                top - bottom,
                facecolor=rgb + (1.0,),
                edgecolor=node_edge_color,
                linewidth=node_edge_width,
                zorder=3,
            )
            ax.add_patch(rect)
            if show_bar_labels:
                cluster_id = int(flow.cluster_ids_per_tau[col][k])
                ax.text(
                    x_centers[col] + 0.5 * node_width + 0.01 * x_step,
                    0.5 * (top + bottom),
                    f"C{cluster_id} ({int(sizes[k])})",
                    ha="left",
                    va="center",
                    fontsize=bar_label_fontsize,
                    zorder=4,
                )

    if show_column_labels:
        # Column labels above the highest bar (y=0).
        y_label = column_label_pad * n_total
        for col in range(n_cols):
            tau = float(flow.tau_values[col])
            n_cl = len(flow.sizes[col])
            ax.text(
                x_centers[col],
                y_label,
                column_label_format.format(tau=tau, n=n_cl),
                ha="center",
                va="bottom",
                fontsize=column_label_fontsize,
            )

    # Padding around content.
    pad_x = 0.5 * node_width + 0.15 * x_step
    pad_y_top = (column_label_pad + 0.05) * n_total
    # Lowest bar bottom across columns.
    min_bottom = min(b[-1][1] for b in bar_bounds)
    pad_y_bot = 0.05 * n_total
    ax.set_xlim(x_centers[0] - pad_x, x_centers[-1] + pad_x)
    ax.set_ylim(min_bottom - pad_y_bot, pad_y_top)
    ax.set_aspect("auto")

    if hide_axes:
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
    if title is not None:
        ax.set_title(title, fontsize=11, loc="left", pad=8)

    return fig, ax
