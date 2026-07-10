"""Hierarchy-bundled chord plots (graph-tool curvy edges through LRG).

Split out of `network_templates.py` on 2026-05-29 (Phase 4-B split 7/7).
Holds the constants `HIERARCHY_*` + `DENDROGRAM_*`, the depth-2 chord
layout builder, the graph-tool renderer, and the matplotlib overlays
that combine into `plot_chord_with_dendrogram` and
`plot_chord_with_highlight`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, List, Literal, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy.cluster.hierarchy import fcluster


__all__ = [
    "HIERARCHY_DEFAULT_K_CLUSTERS",
    "HIERARCHY_DEFAULT_BETA",
    "HIERARCHY_DEFAULT_RENDER_PX",
    "HIERARCHY_DEFAULT_FIT_VIEW",
    "HIERARCHY_R_OUTER",
    "HIERARCHY_R_CLUSTER",
    "HIERARCHY_BG_RGB",
    "HIERARCHY_BG_ALPHA",
    "HIERARCHY_BG_PEN",
    "DENDROGRAM_DEFAULT_COLOR",
    "DENDROGRAM_DEFAULT_LW",
    "DENDROGRAM_R_OUTER",
    "DENDROGRAM_R_INNER",
    "build_chord_depth2_layout",
    "render_hierarchy_chord",
    "draw_circular_dendrogram_overlay",
    "draw_radial_leaf_labels",
    "plot_chord_with_dendrogram",
    "plot_chord_with_highlight",
]


# --- Hierarchy-bundled chord (graph-tool curvy edges through LRG) ---------
# Used by :func:`plot_chord_with_dendrogram` /
# :func:`plot_chord_with_highlight` and by the templates
# ``circular_dendrogram_network`` and ``chord_highlighted_edges``.
HIERARCHY_DEFAULT_K_CLUSTERS: int = 7
HIERARCHY_DEFAULT_BETA: float = 0.92        # bezier bundle tightness
HIERARCHY_DEFAULT_RENDER_PX: int = 2200     # PNG raster the chord renders to
HIERARCHY_DEFAULT_FIT_VIEW: float = 0.92    # graph-tool fit_view (data → canvas)
HIERARCHY_R_OUTER: float = 1.0              # leaves on the outer ring
HIERARCHY_R_CLUSTER: float = 0.55           # depth-2 cluster centroids

# Faint-grey "null backbone" defaults — every non-highlighted pair.
HIERARCHY_BG_RGB: Tuple[float, float, float] = (0.62, 0.62, 0.62)
HIERARCHY_BG_ALPHA: float = 0.085
HIERARCHY_BG_PEN: float = 0.42

# Circular-dendrogram overlay defaults (matplotlib polylines on top of
# the chord PNG embed).
DENDROGRAM_DEFAULT_COLOR: Tuple[float, float, float, float] = (0.12, 0.12, 0.12, 0.70)
DENDROGRAM_DEFAULT_LW: float = 0.75
DENDROGRAM_R_OUTER: float = 0.92            # just inside the leaf markers
DENDROGRAM_R_INNER: float = 0.10            # inner extent (root ≈ origin)
DENDROGRAM_N_ARC_POINTS: int = 36


def build_chord_depth2_layout(
    Z: NDArray,
    N: int,
    *,
    k_clusters: int = HIERARCHY_DEFAULT_K_CLUSTERS,
    r_outer: float = HIERARCHY_R_OUTER,
    r_cluster: float = HIERARCHY_R_CLUSTER,
):
    """Depth-2 hierarchy tree + equiangular leaf ring for the chord.

    Parameters
    ----------
    Z : scipy linkage matrix.
    N : number of leaves (≤ ``Z.shape[0] + 1``).
    k_clusters : depth-2 cut of ``Z`` (``maxclust``); leaves grouped
        under cluster T-junctions.
    r_outer : radius of the leaf ring (in chord-local space).
    r_cluster : radius of the cluster centroids.

    Returns
    -------
    dict with keys:

        ``t``                : graph-tool ``Graph`` (directed tree) with
                               ``N + n_cl + 1`` vertices (leaves 0..N-1,
                               clusters N..N+n_cl-1, root N+n_cl).
        ``root_vid``         : root vertex id in ``t``.
        ``tpos``             : ``t.new_vertex_property('vector<double>')``
                               with positions filled.
        ``theta_of_leaf``    : ``(N,)`` array, angular position of each
                               leaf (radians).
        ``dendro_order``     : leaf ids in ``leaves_list(Z)`` order.
        ``ordered_labels``   : the depth-2 cluster ids in the same order
                               they first appear walking ``dendro_order``.
        ``cluster_leaves``   : dict{cluster_id → list of leaf ids}.
    """
    import graph_tool.all as gt
    import math

    from scipy.cluster.hierarchy import fcluster, leaves_list

    k = int(min(max(2, k_clusters), N - 1))
    labels = fcluster(Z, t=k, criterion="maxclust")[:N]

    dendro_order = [int(x) for x in leaves_list(Z) if int(x) < N]
    seen: Dict[int, int] = {}
    cluster_leaves: Dict[int, List[int]] = {}
    for leaf in dendro_order:
        lab = int(labels[leaf])
        if lab not in seen:
            seen[lab] = len(seen)
            cluster_leaves[lab] = []
        cluster_leaves[lab].append(leaf)
    ordered_labels = sorted(cluster_leaves.keys(), key=lambda x: seen[x])

    n_cl = len(ordered_labels)
    t = gt.Graph(directed=True)
    t.add_vertex(N + n_cl + 1)
    root_vid = N + n_cl
    for j, lab in enumerate(ordered_labels):
        cluster_vid = N + j
        t.add_edge(t.vertex(root_vid), t.vertex(cluster_vid))
        for leaf in cluster_leaves[lab]:
            t.add_edge(t.vertex(cluster_vid), t.vertex(leaf))

    n_leaves_in_order = len(dendro_order)
    theta_of_leaf = np.zeros(N)
    for pos_in_order, leaf_id in enumerate(dendro_order):
        # Subtract π/2 so leaf 0 sits at the top (12 o'clock).
        theta_of_leaf[leaf_id] = (
            2 * math.pi * pos_in_order / n_leaves_in_order - math.pi / 2
        )

    tpos = t.new_vertex_property("vector<double>")
    for leaf_id in dendro_order:
        th = float(theta_of_leaf[leaf_id])
        tpos[t.vertex(leaf_id)] = [
            r_outer * math.cos(th), r_outer * math.sin(th),
        ]
    for j, lab in enumerate(ordered_labels):
        sum_cos = sum(math.cos(theta_of_leaf[l]) for l in cluster_leaves[lab])
        sum_sin = sum(math.sin(theta_of_leaf[l]) for l in cluster_leaves[lab])
        avg_th = math.atan2(sum_sin, sum_cos)
        tpos[t.vertex(N + j)] = [
            r_cluster * math.cos(avg_th), r_cluster * math.sin(avg_th),
        ]
    tpos[t.vertex(root_vid)] = [0.0, 0.0]

    return dict(
        t=t,
        root_vid=root_vid,
        tpos=tpos,
        theta_of_leaf=theta_of_leaf,
        dendro_order=dendro_order,
        ordered_labels=ordered_labels,
        cluster_leaves=cluster_leaves,
    )


def render_hierarchy_chord(
    out_path: Path,
    *,
    A: NDArray,
    probes: Sequence[str],
    Z: NDArray,
    edge_indices: NDArray,
    edge_rgba: NDArray,
    edge_pen: NDArray,
    edge_order: NDArray,
    vertex_rgba: Optional[NDArray] = None,
    vertex_size: float = 24.0,
    vertex_outline_rgba: Tuple[float, float, float, float] = (0.10, 0.10, 0.10, 0.95),
    k_clusters: int = HIERARCHY_DEFAULT_K_CLUSTERS,
    beta: float = HIERARCHY_DEFAULT_BETA,
    output_size: int = HIERARCHY_DEFAULT_RENDER_PX,
    fit_view: float = HIERARCHY_DEFAULT_FIT_VIEW,
) -> Dict:
    """Render a hierarchy-bundled chord network to a PNG via graph-tool.

    Pre-built edge arrays in upper-triangle pair order, so the caller
    fully owns the edge-style policy (default-style, highlight overlay,
    sign-aware, etc.).  ``vertex_rgba=None`` → :func:`shaft_colors`.

    Parameters
    ----------
    out_path : output PNG path (graph-tool writes here).
    A : FC matrix, ``(N, N)``.  Only ``A.shape[0]`` is used.
    probes : per-contact shaft labels.
    Z : scipy linkage matrix (the LRG dendrogram).
    edge_indices : ``(E, 2)`` array of ``(i, j)`` leaf-id pairs to draw.
    edge_rgba    : ``(E, 4)`` RGBA per drawn edge.
    edge_pen     : ``(E,)`` line widths per drawn edge.
    edge_order   : ``(E,)`` draw order (higher → on top).
    vertex_rgba  : ``(N, 4)`` RGBA per leaf, or ``None`` → shaft palette.
    vertex_size  : graph-tool vertex marker size.
    k_clusters   : depth-2 cut of ``Z`` for the chord backbone.
    beta         : bezier bundle tightness (0 chord, 1 fully bundled).
    output_size  : square PNG side in pixels.
    fit_view     : graph-tool ``fit_view`` (data bbox → canvas fraction).

    Returns
    -------
    dict with keys ``Z``, ``N``, ``theta_of_leaf``, ``R_outer``,
    ``output_size``, ``fit_view`` — feed straight into the overlay
    drawers below.
    """
    import graph_tool.all as gt

    N = int(min(A.shape[0], len(probes)))
    edge_indices = np.asarray(edge_indices, dtype=int)
    edge_rgba = np.asarray(edge_rgba, dtype=float)
    edge_pen = np.asarray(edge_pen, dtype=float)
    edge_order = np.asarray(edge_order, dtype=float)
    assert edge_rgba.shape == (edge_indices.shape[0], 4), \
        f"edge_rgba shape mismatch: {edge_rgba.shape} vs expected ({edge_indices.shape[0]}, 4)"

    layout = build_chord_depth2_layout(Z, N, k_clusters=k_clusters)
    t = layout["t"]
    tpos = layout["tpos"]
    theta_of_leaf = layout["theta_of_leaf"]

    g = gt.Graph(directed=False)
    g.add_vertex(N)

    ecolor = g.new_edge_property("vector<double>")
    pen = g.new_edge_property("double")
    eord = g.new_edge_property("double")

    for k_e, (i, j) in enumerate(edge_indices):
        if i >= N or j >= N or i == j:
            continue
        e = g.add_edge(int(i), int(j))
        ecolor[e] = [float(c) for c in edge_rgba[k_e]]
        pen[e] = float(edge_pen[k_e])
        eord[e] = float(edge_order[k_e])

    pos = g.new_vertex_property("vector<double>")
    for v in g.vertices():
        tv = t.vertex(int(v))
        pos[v] = [float(tpos[tv][0]), float(tpos[tv][1])]

    cts = gt.get_hierarchy_control_points(g, t, tpos, beta=beta)

    if vertex_rgba is None:
        shaft_clr = shaft_colors(probes[:N])
        vertex_rgba = np.array(
            [[c[0], c[1], c[2], c[3] if len(c) == 4 else 1.0] for c in shaft_clr]
        )
    vertex_rgba = np.asarray(vertex_rgba, dtype=float)

    vcolor = g.new_vertex_property("vector<double>")
    vsize = g.new_vertex_property("double")
    for v in g.vertices():
        vid = int(v)
        vcolor[v] = [float(c) for c in vertex_rgba[vid]]
        vsize[v] = float(vertex_size)

    gt.graph_draw(
        g, pos=pos,
        output=str(out_path),
        output_size=(int(output_size), int(output_size)),
        vertex_fill_color=vcolor,
        vertex_color=list(vertex_outline_rgba),
        vertex_size=vsize,
        edge_color=ecolor,
        edge_pen_width=pen,
        eorder=eord,
        edge_gradient=[],
        edge_control_points=cts,
        bg_color=[1.0, 1.0, 1.0, 1.0],
        fit_view=float(fit_view),
    )

    return dict(
        Z=np.asarray(Z),
        N=N,
        theta_of_leaf=theta_of_leaf,
        R_outer=HIERARCHY_R_OUTER,
        output_size=int(output_size),
        fit_view=float(fit_view),
    )


def _chord_px_transform(render_meta: Dict):
    """Build the ``(x, y) → pixel`` transform for overlay drawers.

    graph-tool's ``fit_view=fv`` scales the data bbox (here ≈ [-R, +R]²)
    to occupy ``fv * H`` pixels centered on the canvas; image coords
    have origin top-left and y growing DOWN — hence the y-flip.
    """
    H = float(render_meta["output_size"])
    fv = float(render_meta["fit_view"])
    R = float(render_meta["R_outer"])
    scale_px = fv * H / (2 * R)
    cx, cy = H / 2.0, H / 2.0

    def to_px(x: float, y: float) -> Tuple[float, float]:
        return cx + x * scale_px, cy - y * scale_px

    return to_px, (cx, cy), scale_px


def draw_circular_dendrogram_overlay(
    ax: plt.Axes,
    *,
    render_meta: Dict,
    color: Tuple[float, float, float, float] = DENDROGRAM_DEFAULT_COLOR,
    lw: float = DENDROGRAM_DEFAULT_LW,
    r_outer: float = DENDROGRAM_R_OUTER,
    r_inner: float = DENDROGRAM_R_INNER,
    n_arc_points: int = DENDROGRAM_N_ARC_POINTS,
    zorder: int = 7,
) -> None:
    """Draw the full circular dendrogram from ``render_meta['Z']``
    as matplotlib polylines on top of an imshow'd chord PNG.

    Each merge in ``Z`` contributes (a) two radial lines from each
    child's centroid inward to the merge radius, plus (b) one CCW arc
    joining the two child centroid angles at the merge radius.  The arc
    direction follows the children's leaves-list index ordering — no
    shortest-arc heuristic — so subtrees that span more than half the
    ring still connect through the correct side.

    Internal nodes sit at ``r = R_outer · (1 − h / h_max)``, with ``h``
    the merge height and ``h_max`` the deepest merge in ``Z``.
    """
    import math

    from scipy.cluster.hierarchy import leaves_list

    Z = render_meta["Z"]
    N = int(render_meta["N"])
    theta_of_leaf = render_meta["theta_of_leaf"]
    to_px, _, _ = _chord_px_transform(render_meta)

    h_max = float(Z[-1, 2])
    n_total = 2 * N - 1
    node_r = np.zeros(n_total)
    node_theta = np.zeros(n_total)
    node_idx_min = np.zeros(n_total, dtype=int)
    node_idx_max = np.zeros(n_total, dtype=int)

    leaf_order = [int(x) for x in leaves_list(Z) if int(x) < N]
    pos_in_order = {leaf_id: idx for idx, leaf_id in enumerate(leaf_order)}
    n_leaves = len(leaf_order)

    for i in range(N):
        node_r[i] = r_outer
        node_theta[i] = float(theta_of_leaf[i])
        slot = pos_in_order.get(i, i)
        node_idx_min[i] = slot
        node_idx_max[i] = slot

    def _angle_from_slot(slot_idx: float) -> float:
        return 2 * math.pi * slot_idx / n_leaves - math.pi / 2

    n_rows = int(Z.shape[0])
    for r in range(n_rows):
        node_id = N + r
        h = float(Z[r, 2])
        c1, c2 = int(Z[r, 0]), int(Z[r, 1])
        if c1 >= n_total or c2 >= n_total:
            continue
        node_r[node_id] = r_outer - (h / max(h_max, 1e-12)) * (r_outer - r_inner)
        node_idx_min[node_id] = min(node_idx_min[c1], node_idx_min[c2])
        node_idx_max[node_id] = max(node_idx_max[c1], node_idx_max[c2])
        node_theta[node_id] = _angle_from_slot(
            0.5 * (node_idx_min[node_id] + node_idx_max[node_id])
        )

    for r in range(n_rows):
        node_id = N + r
        c1, c2 = int(Z[r, 0]), int(Z[r, 1])
        if c1 >= n_total or c2 >= n_total:
            continue
        if node_idx_min[c1] <= node_idx_min[c2]:
            c_left, c_right = c1, c2
        else:
            c_left, c_right = c2, c1
        merge_r = node_r[node_id]
        t_left = node_theta[c_left]
        t_right = node_theta[c_right]

        for c in (c_left, c_right):
            cr = node_r[c]
            th = node_theta[c]
            x1, y1 = cr * math.cos(th), cr * math.sin(th)
            x2, y2 = merge_r * math.cos(th), merge_r * math.sin(th)
            px1, py1 = to_px(x1, y1)
            px2, py2 = to_px(x2, y2)
            ax.plot([px1, px2], [py1, py2],
                    color=color, lw=lw, zorder=zorder,
                    solid_capstyle="round")

        dth = (t_right - t_left) % (2 * math.pi)
        arc_ths = t_left + np.linspace(0.0, dth, n_arc_points)
        arc_x = merge_r * np.cos(arc_ths)
        arc_y = merge_r * np.sin(arc_ths)
        arc_px = [to_px(x, y) for x, y in zip(arc_x, arc_y)]
        ax.plot([p[0] for p in arc_px], [p[1] for p in arc_px],
                color=color, lw=lw, zorder=zorder,
                solid_capstyle="round")


def draw_radial_leaf_labels(
    ax: plt.Axes,
    *,
    render_meta: Dict,
    labels: Sequence[str],
    offset_frac: float = 0.015,
    fontsize: float = 4.8,
    color: str = "0.18",
    zorder: int = 8,
) -> None:
    """Draw rotated contact labels around the outer leaf ring of a
    chord PNG that's already imshow'd in ``ax``.

    Standard chord-diagram idiom: outward-pointing on the right half,
    flipped to point inward on the left half so all text stays
    upright-readable.
    """
    import math

    N = int(render_meta["N"])
    theta_of_leaf = render_meta["theta_of_leaf"]
    H = float(render_meta["output_size"])
    to_px, (cx, cy), _ = _chord_px_transform(render_meta)
    offset_px = H * float(offset_frac)

    for i in range(N):
        if i >= len(labels):
            break
        lab = str(labels[i])
        if not lab:
            continue
        lx, ly = math.cos(theta_of_leaf[i]), math.sin(theta_of_leaf[i])
        lpx, lpy = to_px(lx, ly)
        ux, uy = (lpx - cx), (lpy - cy)
        norm = math.hypot(ux, uy)
        if norm < 1e-6:
            continue
        ux /= norm; uy /= norm
        tx, ty = lpx + ux * offset_px, lpy + uy * offset_px
        angle_deg = math.degrees(math.atan2(-(lpy - cy), lpx - cx))
        if -90 < angle_deg <= 90:
            rotation, ha = angle_deg, "left"
        else:
            rotation, ha = angle_deg + 180, "right"
        ax.text(tx, ty, lab,
                rotation=rotation, rotation_mode="anchor",
                ha=ha, va="center",
                fontsize=fontsize, color=color,
                zorder=zorder, clip_on=False)


def _chord_default_edge_arrays(
    A: NDArray,
    probes: Sequence[str],
    *,
    gamma: float = 2.0,
    width_range: Tuple[float, float] = (0.0, 3.5),
    alpha_range: Tuple[float, float] = (0.0, 0.85),
) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """Build the default edge arrays for ``plot_chord_with_dendrogram``.

    Same-probe pairs painted in the shaft's tab20 colour, cross-probe
    in mid-gray.  Width and alpha scale as ``t = (|w|/|w|_max)^gamma``
    so the strong cross-probe edges read through the dendrogram
    overlay without overwhelming it.

    Defaults reflect the same diagnosis as ``matrix_plus_network``:
    γ=2 with zero-floor width/alpha makes the bottom ~half of edges
    genuinely invisible, leaving the strong-edge skeleton (top ~30%)
    as visible context against the dendrogram polyline overlay.  The
    earlier defaults (γ=0.55, non-zero floors) produced a dense dark
    spaghetti through the chord centre that drowned the dendrogram;
    γ=3 over-shot in the other direction and lost the FC context
    entirely.
    """
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = np.abs(A[r, c])
    w_max = float(w.max()) if w.size and w.max() > 0 else 1.0
    t = (w / w_max) ** float(gamma)
    widths = width_range[0] + t * (width_range[1] - width_range[0])
    alphas = alpha_range[0] + t * (alpha_range[1] - alpha_range[0])

    shaft_map = _shaft_color_map(probes[:N])
    rgba = np.zeros((r.size, 4), dtype=float)
    for k in range(r.size):
        i, j = int(r[k]), int(c[k])
        if probes[i] == probes[j]:
            rgb = shaft_map[probes[i]][:3]
        else:
            rgb = CROSS_PROBE_RGB
        rgba[k] = [rgb[0], rgb[1], rgb[2], float(alphas[k])]

    idx = np.stack([r, c], axis=1)
    return idx, rgba, widths, w  # w doubles as draw-order ⇒ strong on top


def plot_chord_with_dendrogram(
    ax: plt.Axes,
    A: NDArray,
    probes: Sequence[str],
    lrg,
    *,
    k_clusters: int = HIERARCHY_DEFAULT_K_CLUSTERS,
    beta: float = HIERARCHY_DEFAULT_BETA,
    render_px: int = HIERARCHY_DEFAULT_RENDER_PX,
    fit_view: float = HIERARCHY_DEFAULT_FIT_VIEW,
    show_dendrogram: bool = True,
    dendrogram_kwargs: Optional[Dict] = None,
    labels: Optional[Sequence[str]] = None,
    label_kwargs: Optional[Dict] = None,
    edge_gamma: float = 2.0,
    edge_width_range: Tuple[float, float] = (0.0, 3.5),
    edge_alpha_range: Tuple[float, float] = (0.0, 0.85),
    vertex_size: float = 24.0,
    tmp_dir: Optional[Path] = None,
) -> Dict:
    """Composite: chord PNG + (optional) circular-dendrogram overlay +
    (optional) radial leaf labels.

    Template 1 in this family.  Read the matrix once with the canonical
    same-probe / cross-probe rule, then overlay the full LRG dendrogram
    (every merge in ``lrg.linkage_matrix``) as polylines so the
    communication-hierarchy structure reads through the curvy edges.

    Returns the ``render_meta`` dict (also written to ``ax`` as the
    embedded PNG).
    """
    import tempfile

    Z = np.asarray(lrg.linkage_matrix)
    idx, rgba, pen, eord = _chord_default_edge_arrays(
        A, probes,
        gamma=edge_gamma,
        width_range=edge_width_range,
        alpha_range=edge_alpha_range,
    )

    cleanup = None
    if tmp_dir is None:
        cleanup = tempfile.TemporaryDirectory(prefix="chord_dendro_")
        tmp_dir = Path(cleanup.name)
    out_png = Path(tmp_dir) / "chord.png"
    try:
        render_meta = render_hierarchy_chord(
            out_png,
            A=A, probes=probes, Z=Z,
            edge_indices=idx, edge_rgba=rgba,
            edge_pen=pen, edge_order=eord,
            vertex_size=vertex_size,
            k_clusters=k_clusters, beta=beta,
            output_size=render_px, fit_view=fit_view,
        )
        from matplotlib.image import imread
        ax.imshow(imread(out_png), interpolation="lanczos")
    finally:
        if cleanup is not None:
            cleanup.cleanup()
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    if show_dendrogram:
        draw_circular_dendrogram_overlay(
            ax, render_meta=render_meta, **(dendrogram_kwargs or {}),
        )
    if labels is not None:
        draw_radial_leaf_labels(
            ax, render_meta=render_meta, labels=labels,
            **(label_kwargs or {}),
        )
    return render_meta


def plot_chord_with_highlight(
    ax: plt.Axes,
    A: NDArray,
    probes: Sequence[str],
    lrg,
    *,
    highlights: Sequence[Dict],
    k_clusters: int = HIERARCHY_DEFAULT_K_CLUSTERS,
    beta: float = HIERARCHY_DEFAULT_BETA,
    render_px: int = HIERARCHY_DEFAULT_RENDER_PX,
    fit_view: float = HIERARCHY_DEFAULT_FIT_VIEW,
    background_rgb: Tuple[float, float, float] = HIERARCHY_BG_RGB,
    background_alpha: float = HIERARCHY_BG_ALPHA,
    background_pen: float = HIERARCHY_BG_PEN,
    cross_probe_only: bool = False,
    show_dendrogram: bool = False,
    dendrogram_kwargs: Optional[Dict] = None,
    labels: Optional[Sequence[str]] = None,
    label_kwargs: Optional[Dict] = None,
    vertex_size: float = 24.0,
    vertex_rgba: Optional[NDArray] = None,
    tmp_dir: Optional[Path] = None,
) -> Dict:
    """Composite: chord with HIGHLIGHTED edge subsets on top of a faint
    grey "null backbone" of every other pair.

    Template 2 in this family.  Generalises the trace/anti pattern from
    ``preprint_07_beta_rho_split_figure_test2.py`` panel (c):

    Parameters
    ----------
    highlights : list of dicts, one per highlight class.  Each dict:

        ``pairs``        : array of ``(i, j)`` tuples or 1D index array
                           into the upper-triangle pair order (the
                           canonical numpy ordering
                           ``np.triu_indices(N, k=1)``).  REQUIRED.
        ``rgb``          : highlight colour (3-tuple).  REQUIRED.
        ``alpha_range``  : (min, max) alpha; alpha scales with the
                           highlight's magnitudes if provided, else
                           uses ``(min + max)/2``.  Default ``(0.45, 0.95)``.
        ``width_range``  : (min, max) line width.  Default ``(1.4, 6.4)``.
        ``magnitudes``   : optional per-pair magnitude (positive) used
                           to scale alpha and width within the
                           highlight class.  ``None`` → all members
                           drawn at width_range[1] and
                           ``(alpha_range[0] + alpha_range[1])/2``.
        ``gamma``        : optional γ in ``t = (m/m_max)**gamma``.
                           Default 0.55 (mild compression).
        ``draw_order_offset`` : float added on top of the per-pair
                           magnitude so multiple classes can be
                           layered.  Default 1000 × (1-based class index).

    cross_probe_only : if True, the BACKGROUND draws only cross-probe
        pairs (same-probe filtered out).  Highlight pairs are drawn
        verbatim either way.  Default False.

    show_dendrogram : default False — the highlight focus IS the
        figure; turn on if you want both at once (Template 1+2 hybrid).
    """
    import tempfile

    N = int(min(A.shape[0], len(probes)))
    Z = np.asarray(lrg.linkage_matrix)
    r, c = np.triu_indices(N, k=1)
    n_pairs = r.size

    # ----- background: every pair drawn as faint grey, low pen, low draw order.
    keep_bg = np.ones(n_pairs, dtype=bool)
    if cross_probe_only:
        for k in range(n_pairs):
            if probes[int(r[k])] == probes[int(c[k])]:
                keep_bg[k] = False

    bg_idx = np.stack([r[keep_bg], c[keep_bg]], axis=1)
    bg_rgba = np.tile(
        [background_rgb[0], background_rgb[1], background_rgb[2],
         float(background_alpha)],
        (bg_idx.shape[0], 1),
    )
    bg_pen = np.full(bg_idx.shape[0], float(background_pen))
    bg_eord = np.zeros(bg_idx.shape[0])

    # ----- highlights: stacked on top, one class at a time.
    hl_idx_list = [bg_idx]
    hl_rgba_list = [bg_rgba]
    hl_pen_list = [bg_pen]
    hl_eord_list = [bg_eord]

    def _resolve_pairs(pairs, N) -> NDArray:
        arr = np.asarray(pairs)
        if arr.ndim == 1:
            # interpret as flat indices into triu_indices(N, k=1)
            return np.stack([r[arr], c[arr]], axis=1)
        if arr.ndim == 2 and arr.shape[1] == 2:
            return arr.astype(int)
        raise ValueError(
            f"highlight['pairs'] must be 1D (flat triu indices) or "
            f"(M, 2) (i,j) pairs, got shape {arr.shape}"
        )

    for k_class, h in enumerate(highlights, start=1):
        pairs = _resolve_pairs(h["pairs"], N)
        rgb = h["rgb"]
        width_range = tuple(h.get("width_range", (1.4, 6.4)))
        alpha_range = tuple(h.get("alpha_range", (0.45, 0.95)))
        gamma = float(h.get("gamma", 0.55))
        draw_offset = float(h.get("draw_order_offset", 1000.0 * k_class))
        mags = h.get("magnitudes", None)

        M = pairs.shape[0]
        if mags is None or M == 0:
            t_vec = np.full(M, 0.5)
        else:
            mags_arr = np.asarray(mags, dtype=float)
            assert mags_arr.shape == (M,), \
                f"magnitudes shape {mags_arr.shape} ≠ pairs count ({M},)"
            m_max = float(mags_arr.max()) if mags_arr.size and mags_arr.max() > 0 else 1.0
            t_vec = (mags_arr / m_max) ** gamma

        widths = width_range[0] + t_vec * (width_range[1] - width_range[0])
        alphas = alpha_range[0] + t_vec * (alpha_range[1] - alpha_range[0])
        rgba = np.tile([rgb[0], rgb[1], rgb[2], 0.0], (M, 1))
        rgba[:, 3] = alphas
        if mags is None:
            eord = np.full(M, draw_offset)
        else:
            eord = draw_offset + np.asarray(mags, dtype=float)

        hl_idx_list.append(pairs)
        hl_rgba_list.append(rgba)
        hl_pen_list.append(widths)
        hl_eord_list.append(eord)

    idx = np.concatenate(hl_idx_list, axis=0)
    rgba = np.concatenate(hl_rgba_list, axis=0)
    pen = np.concatenate(hl_pen_list, axis=0)
    eord = np.concatenate(hl_eord_list, axis=0)

    cleanup = None
    if tmp_dir is None:
        cleanup = tempfile.TemporaryDirectory(prefix="chord_highlight_")
        tmp_dir = Path(cleanup.name)
    out_png = Path(tmp_dir) / "chord.png"
    try:
        render_meta = render_hierarchy_chord(
            out_png,
            A=A, probes=probes, Z=Z,
            edge_indices=idx, edge_rgba=rgba,
            edge_pen=pen, edge_order=eord,
            vertex_size=vertex_size,
            vertex_rgba=vertex_rgba,
            k_clusters=k_clusters, beta=beta,
            output_size=render_px, fit_view=fit_view,
        )
        from matplotlib.image import imread
        ax.imshow(imread(out_png), interpolation="lanczos")
    finally:
        if cleanup is not None:
            cleanup.cleanup()
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    if show_dendrogram:
        draw_circular_dendrogram_overlay(
            ax, render_meta=render_meta, **(dendrogram_kwargs or {}),
        )
    if labels is not None:
        draw_radial_leaf_labels(
            ax, render_meta=render_meta, labels=labels,
            **(label_kwargs or {}),
        )
    return render_meta
