"""Edge styling, node decoration, and graph-tool panel rendering helpers.

Promoted from ``scripts/01_compute/figures_for_notes/_shared.py`` on
2026-05-29 (per the library-naming meta-rule locked 2026-05-28).

This module is the canonical home for FC-method-aware edge scaling
(``scale_edge_weights``), matplotlib edge drawing with same-probe vs.
cross-probe highlighting (``draw_network_edges``), and graph-tool panel
renderers driven either by nested-SBM inference (``render_sbm_panel``)
or by an LRG dendrogram (``render_lrg_panel``).
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from lrg_eegfc.config.const import (
    EDGE_ALPHA_RANGE,
    EDGE_RANK_ALPHA,
    EDGE_RANK_GAMMA,
    EDGE_RANK_WIDTH,
    EDGE_WIDTH_RANGE,
    LAYOUT_SEED_DEFAULT,
)

from .network_layouts import compute_percolation_threshold


__all__ = [
    "EDGE_GAMMA",
    "scale_edge_weights",
    "draw_network_edges",
    "render_sbm_panel",
    "render_lrg_panel",
]


#: Gamma exponent per FC method for edge width/alpha scaling.
#: Higher gamma -> only top-weight edges are visually prominent.
EDGE_GAMMA = {
    "msc": 2.0,
    "imcoh": 6.0,
    "imcoh_abs": 6.0,
    "imcoh_sq": 6.0,
    "corr": 2.0,
}


def scale_edge_weights(
    weights: NDArray,
    fc_method: str = "imcoh_abs",
) -> tuple[NDArray, NDArray]:
    """Scale edge weights to visual width and alpha using method-aware config.

    Uses a gamma transform t = (w / w_max) ** gamma so that weak edges
    get small widths and high-weight edges stand out. Works well for
    near-uniform distributions (ImCoh) at gamma >= 6 and for heavy-tailed
    ones (MSC) at gamma ~ 2.
    """
    w_min, w_max = EDGE_WIDTH_RANGE.get(fc_method, (0.15, 4.0))
    a_min, a_max = EDGE_ALPHA_RANGE.get(fc_method, (0.03, 0.9))
    gamma = EDGE_GAMMA.get(fc_method, 2.0)

    w_norm = weights / (weights.max() + 1e-30)
    t = w_norm ** gamma

    widths = w_min + t * (w_max - w_min)
    alphas = a_min + t * (a_max - a_min)
    return widths, alphas


def _probe_color_map(probes: Sequence[str]) -> dict:
    """Return {probe_id: rgba} using the same tab20 palette as the nodes."""
    uniq = sorted(set(probes))
    cmap = plt.get_cmap("tab20", len(uniq))
    return {p: cmap(i) for i, p in enumerate(uniq)}


def draw_network_edges(
    ax,
    pos_arr: NDArray,
    mat: NDArray,
    fc_method: str = "imcoh_abs",
    probe_labels: Optional[Sequence[str]] = None,
    highlight_same_probe: bool = True,
    base_color=(0.0, 0.0, 0.0),
    gamma: Optional[float] = None,
    wmin: Optional[float] = None,
    wmax: Optional[float] = None,
    amin: Optional[float] = None,
    amax: Optional[float] = None,
    scaling: str = "rank",
    value_range: Optional[tuple[float, float]] = None,
):
    """Edge drawing per method (gamma, width and alpha from
    :mod:`config.const`, overridable per call).

    *scaling*:
      - ``"rank"`` (default): ``t = ((rank+1)/N)**gamma`` -- scale-invariant
        across patients but visual density depends on edge count.
      - ``"value"``: per-patient min-max normalise weights to [0,1], then
        ``t = w_norm**gamma`` -- same width/alpha mapping applies to the
        same *relative* weight level across patients.

    If *highlight_same_probe* is ``True`` and *probe_labels* is given,
    same-probe edges are drawn in the shaft's tab20 colour at full alpha
    (so they stand out); cross-probe edges keep *base_color* with the
    value/rank-based alpha fade.
    """
    from matplotlib.collections import LineCollection

    _g = EDGE_RANK_GAMMA.get(fc_method, 2.0)
    _wmn, _wmx = EDGE_RANK_WIDTH.get(fc_method, (0.15, 4.0))
    _amn, _amx = EDGE_RANK_ALPHA.get(fc_method, (0.06, 0.9))
    gamma = _g if gamma is None else gamma
    wmin = _wmn if wmin is None else wmin
    wmax = _wmx if wmax is None else wmax
    amin = _amn if amin is None else amin
    amax = _amx if amax is None else amax

    N = mat.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = mat[r, c]
    active = np.where(w > 0)[0]
    if active.size == 0:
        return
    w_act = w[active]
    ii, jj = r[active], c[active]

    if highlight_same_probe and probe_labels is not None:
        same_all = np.array(
            [probe_labels[ii[k]] == probe_labels[jj[k]]
             for k in range(len(active))],
            dtype=bool,
        )
        probe_cmap = _probe_color_map(probe_labels)
    else:
        same_all = np.zeros(len(active), dtype=bool)
        probe_cmap = {}

    cross_mask = ~same_all
    t = np.zeros_like(w_act, dtype=float)

    def _value_t(w_sub: np.ndarray) -> np.ndarray:
        if w_sub.size == 0:
            return np.zeros(0)
        lo, hi = (float(np.percentile(w_sub, 1)),
                  float(np.percentile(w_sub, 99)))
        if hi <= lo:
            return np.zeros_like(w_sub)
        wn = np.clip((w_sub - lo) / (hi - lo), 0.0, 1.0)
        return wn ** gamma

    def _rank_t(w_sub: np.ndarray) -> np.ndarray:
        if w_sub.size == 0:
            return np.zeros(0)
        ranks = np.argsort(np.argsort(w_sub)).astype(float)
        return ((ranks + 1) / len(ranks)) ** gamma

    if scaling == "value":
        if value_range is not None:
            lo, hi = float(value_range[0]), float(value_range[1])
            if hi > lo:
                wn = np.clip((w_act - lo) / (hi - lo), 0.0, 1.0)
                t = wn ** gamma
        else:
            t[cross_mask] = _value_t(w_act[cross_mask])
            t[~cross_mask] = _value_t(w_act[~cross_mask])
    else:
        t[cross_mask] = _rank_t(w_act[cross_mask])
        t[~cross_mask] = _rank_t(w_act[~cross_mask])

    widths = wmin + t * (wmax - wmin)
    alphas = amin + t * (amax - amin)

    base_rgb = np.array(base_color, dtype=float)

    # Bucket (lw, alpha) per RGB so ~7k edges collapse into a handful of
    # LineCollections -- keeps vector PDFs small.
    N_LW_BINS = 12
    N_AL_BINS = 12
    ALPHA_INVISIBLE = 0.01

    def _emit_bucketed(segs_arr, widths_arr, alphas_arr, rgb_arr, zorder):
        if len(segs_arr) == 0:
            return
        keep = alphas_arr >= ALPHA_INVISIBLE
        if not keep.any():
            return
        segs_arr = segs_arr[keep]
        widths_arr = widths_arr[keep]
        alphas_arr = alphas_arr[keep]
        rgb_arr = rgb_arr[keep]
        w_lo, w_hi = float(widths_arr.min()), float(widths_arr.max())
        a_lo, a_hi = float(alphas_arr.min()), float(alphas_arr.max())
        if w_hi > w_lo:
            w_idx = np.clip(((widths_arr - w_lo) / (w_hi - w_lo)
                              * (N_LW_BINS - 1)).round().astype(int),
                            0, N_LW_BINS - 1)
            w_bins = w_lo + (w_hi - w_lo) * (np.arange(N_LW_BINS) / (N_LW_BINS - 1))
        else:
            w_idx = np.zeros(len(widths_arr), dtype=int)
            w_bins = np.array([w_lo])
        if a_hi > a_lo:
            a_idx = np.clip(((alphas_arr - a_lo) / (a_hi - a_lo)
                              * (N_AL_BINS - 1)).round().astype(int),
                            0, N_AL_BINS - 1)
            a_bins = a_lo + (a_hi - a_lo) * (np.arange(N_AL_BINS) / (N_AL_BINS - 1))
        else:
            a_idx = np.zeros(len(alphas_arr), dtype=int)
            a_bins = np.array([a_lo])
        rgb_q = np.round(rgb_arr * 255).astype(np.uint8)
        keys = (rgb_q[:, 0].astype(np.int64) << 32
                | rgb_q[:, 1].astype(np.int64) << 24
                | rgb_q[:, 2].astype(np.int64) << 16
                | (w_idx.astype(np.int64) & 0xFF) << 8
                | (a_idx.astype(np.int64) & 0xFF))
        unique_keys, inverse = np.unique(keys, return_inverse=True)
        for b, key in enumerate(unique_keys):
            mask = inverse == b
            r_ = (key >> 32) & 0xFF
            g_ = (key >> 24) & 0xFF
            bl = (key >> 16) & 0xFF
            wi = (key >> 8) & 0xFF
            ai = key & 0xFF
            color = (float(r_) / 255.0, float(g_) / 255.0, float(bl) / 255.0,
                     float(a_bins[ai]))
            ax.add_collection(LineCollection(
                segs_arr[mask],
                colors=[color],
                linewidths=float(w_bins[wi]),
                zorder=zorder, rasterized=False,
            ))

    for layer_mask, layer_zorder, force_alpha in (
        (cross_mask, 1, None),
        (~cross_mask, 2, 1.0),
    ):
        if not layer_mask.any():
            continue
        sub_idx = np.where(layer_mask)[0]
        sub_t = t[sub_idx]
        order = sub_idx[np.argsort(sub_t)]
        segs = np.stack([pos_arr[ii[order]], pos_arr[jj[order]]], axis=1)
        widths_o = widths[order]
        alphas_o = (np.full(len(order), force_alpha)
                    if force_alpha is not None else alphas[order])
        rgb = np.tile(base_rgb, (len(order), 1))
        if layer_zorder == 2:
            for k_pos, k in enumerate(order):
                col = probe_cmap.get(probe_labels[ii[k]])
                if col is not None:
                    rgb[k_pos] = col[:3]
        _emit_bucketed(segs, widths_o, alphas_o, rgb, layer_zorder)


def render_sbm_panel(
    output_path: Path,
    A: NDArray,
    probe_labels: Optional[list[str]] = None,
    shaft_colors: Optional[list] = None,
    *,
    seed: int = LAYOUT_SEED_DEFAULT,
    threshold: Optional[float] = None,
    threshold_scale: float = 1.0,
    highlight_same_probe: bool = True,
    output_size: tuple = (600, 600),
    vertex_size: float = 10.0,
    edge_width_range: tuple = (0.3, 3.5),
    title: Optional[str] = None,
) -> None:
    """Render ONE nested-SBM panel directly via graph-tool's native
    ``state.draw()`` to a PDF file.

    Edges below the percolation threshold (min-weight edge of the maximum
    spanning tree, scaled by ``threshold_scale``) are hidden. SBM
    inference runs on the FULL graph; thresholding is for display only.

    Compose multiple panels into a grid with pdfjam after rendering::

        pdfjam --nup 3x2 panel_*.pdf -o fig_E2_sbm.pdf
    """
    import graph_tool.all as gt
    import gc

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gt.seed_rng(seed)
    N = A.shape[0]

    g = gt.Graph(directed=False)
    g.add_vertex(N)
    weight = g.new_edge_property("double")
    ri, ci = np.triu_indices(N, k=1)
    vals = A[ri, ci]
    mask = vals > 0
    elist = np.column_stack([ri[mask], ci[mask], vals[mask]])
    g.add_edge_list(elist, eprops=[weight])

    state = gt.minimize_nested_blockmodel_dl(
        g, state_args=dict(recs=[weight], rec_types=["real-exponential"])
    )

    theta = (threshold if threshold is not None
             else compute_percolation_threshold(A)) * threshold_scale
    efilter = g.new_edge_property("bool")
    for e in g.edges():
        efilter[e] = float(weight[e]) >= theta
    g.set_edge_filter(efilter)

    cross_color = [0.1, 0.1, 0.1, 0.5]
    ecolor = g.new_edge_property("vector<double>")
    draw_order = g.new_edge_property("double")
    w_max = float(max((float(weight[e]) for e in g.edges()), default=1.0))
    same_offset = w_max + 1.0
    for e in g.edges():
        u = int(e.source()); v = int(e.target())
        same = (highlight_same_probe and probe_labels is not None
                and probe_labels[u] == probe_labels[v])
        if same and shaft_colors is not None:
            c = shaft_colors[u]
            ecolor[e] = ([c[0], c[1], c[2], 0.95]
                         if len(c) == 3 else list(c))
            draw_order[e] = same_offset + float(weight[e])
        elif same:
            ecolor[e] = [0.85, 0.15, 0.15, 0.95]
            draw_order[e] = same_offset + float(weight[e])
        else:
            ecolor[e] = cross_color
            draw_order[e] = float(weight[e])

    vcolor = g.new_vertex_property("vector<double>")
    for v in g.vertices():
        if shaft_colors is not None:
            c = shaft_colors[int(v)]
            vcolor[v] = [c[0], c[1], c[2], 1.0] if len(c) == 3 else list(c)
        else:
            vcolor[v] = [0.5, 0.5, 0.5, 1.0]

    pen = gt.prop_to_size(
        weight, mi=edge_width_range[0], ma=edge_width_range[1],
        power=1, log=True,
    )

    state.draw(
        output=str(output_path),
        output_size=output_size,
        vertex_fill_color=vcolor,
        vertex_color=[0, 0, 0, 1],
        vertex_size=vertex_size,
        edge_color=ecolor,
        edge_pen_width=pen,
        eorder=draw_order,
        edge_gradient=[],
    )

    del state, g, weight, efilter, ecolor, vcolor, pen
    gc.collect()


def render_lrg_panel(
    output_path: Path,
    A: NDArray,
    lrg_result,
    probe_labels: Optional[list[str]] = None,
    shaft_colors: Optional[list] = None,
    *,
    seed: int = LAYOUT_SEED_DEFAULT,
    threshold: Optional[float] = None,
    threshold_scale: float = 1.0,
    highlight_same_probe: bool = True,
    output_size: tuple = (600, 600),
    vertex_size: float = 10.0,
    edge_width_range: tuple = (0.03, 2.5),
    edge_alpha_range: tuple = (0.10, 0.95),
    edge_width_power: float = 4.0,
    beta: float = 0.8,
    k_clusters: int = 30,
) -> None:
    """Render ONE panel via graph-tool using the LRG dendrogram as the
    hierarchy tree (no SBM fit).

    Nodes are placed on the outer ring in dendrogram order via
    ``radial_tree_layout`` on the linkage-derived tree -- same-cluster
    nodes land angularly adjacent by construction. Edges curve through
    the LRG hierarchy via ``get_hierarchy_control_points``.

    ``lrg_result`` must be an ``LRGResult`` (from
    :func:`lrg_eegfc.workflow.lrg.load_lrg_result`).
    """
    import graph_tool.all as gt
    import gc
    from scipy.cluster.hierarchy import fcluster, leaves_list

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gt.seed_rng(seed)

    Z = np.asarray(lrg_result.linkage_matrix)
    M = int(Z.shape[0]) + 1
    if M != A.shape[0]:
        A = A[:M, :M]
        if probe_labels is not None:
            probe_labels = list(probe_labels)[:M]
        if shaft_colors is not None:
            shaft_colors = list(shaft_colors)[:M]
    N = M

    g = gt.Graph(directed=False)
    g.add_vertex(N)
    weight = g.new_edge_property("double")
    ri, ci = np.triu_indices(N, k=1)
    vals = A[ri, ci]
    mask = vals > 0
    elist = np.column_stack([ri[mask], ci[mask], vals[mask]])
    g.add_edge_list(elist, eprops=[weight])

    k = int(min(max(2, k_clusters), N - 1))
    labels = fcluster(Z, t=k, criterion="maxclust")
    dendro_order = leaves_list(Z)
    seen = {}
    cluster_leaves: dict = {}
    for leaf in dendro_order:
        lab = int(labels[leaf])
        if lab not in seen:
            seen[lab] = len(seen)
            cluster_leaves[lab] = []
        cluster_leaves[lab].append(int(leaf))
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
    root = t.vertex(root_vid)

    tpos = gt.radial_tree_layout(t, root)
    cts = gt.get_hierarchy_control_points(g, t, tpos, beta=beta)

    pos = g.new_vertex_property("vector<double>")
    for v in g.vertices():
        tv = t.vertex(int(v))
        pos[v] = [float(tpos[tv][0]), float(tpos[tv][1])]

    if highlight_same_probe and probe_labels is not None:
        A_cross = A.copy()
        for i in range(N):
            for j in range(i + 1, N):
                if probe_labels[i] == probe_labels[j]:
                    A_cross[i, j] = 0.0
                    A_cross[j, i] = 0.0
    else:
        A_cross = A
    theta = (threshold if threshold is not None
             else compute_percolation_threshold(A_cross)) * threshold_scale
    efilter = g.new_edge_property("bool")
    for e in g.edges():
        u = int(e.source()); v = int(e.target())
        same = (highlight_same_probe and probe_labels is not None
                and probe_labels[u] == probe_labels[v])
        if same:
            efilter[e] = True
        else:
            efilter[e] = float(weight[e]) >= theta
    g.set_edge_filter(efilter)

    edge_list = list(g.edges())
    same_flags = np.array([
        (highlight_same_probe and probe_labels is not None
         and probe_labels[int(e.source())] == probe_labels[int(e.target())])
        for e in edge_list
    ], dtype=bool)
    w_vals = np.array([float(weight[e]) for e in edge_list])
    if w_vals.size:
        w_hi = float(np.percentile(w_vals, 95))
        t_vals = (np.clip(w_vals / w_hi, 0.0, 1.0) ** edge_width_power
                  if w_hi > 0 else np.zeros_like(w_vals))
    else:
        t_vals = np.zeros(0)
    a_mi, a_ma = float(edge_alpha_range[0]), float(edge_alpha_range[1])
    alphas = a_mi + t_vals * (a_ma - a_mi)
    w_mi, w_ma = float(edge_width_range[0]), float(edge_width_range[1])
    pen_vals = w_mi + t_vals * (w_ma - w_mi)

    ecolor = g.new_edge_property("vector<double>")
    draw_order = g.new_edge_property("double")
    w_max = float(w_vals.max()) if w_vals.size else 1.0
    same_offset = w_max + 1.0
    for idx, e in enumerate(edge_list):
        u = int(e.source()); v = int(e.target())
        a = float(alphas[idx])
        if same_flags[idx] and shaft_colors is not None:
            c = shaft_colors[u]
            ecolor[e] = [c[0], c[1], c[2], a] if len(c) == 3 else list(c[:3]) + [a]
            draw_order[e] = same_offset + float(weight[e])
        elif same_flags[idx]:
            ecolor[e] = [0.85, 0.15, 0.15, a]
            draw_order[e] = same_offset + float(weight[e])
        else:
            ecolor[e] = [0.1, 0.1, 0.1, a]
            draw_order[e] = float(weight[e])

    vcolor = g.new_vertex_property("vector<double>")
    for v in g.vertices():
        if shaft_colors is not None:
            c = shaft_colors[int(v)]
            vcolor[v] = [c[0], c[1], c[2], 1.0] if len(c) == 3 else list(c)
        else:
            vcolor[v] = [0.5, 0.5, 0.5, 1.0]

    pen = g.new_edge_property("double")
    for idx, e in enumerate(edge_list):
        pen[e] = float(pen_vals[idx])

    gt.graph_draw(
        g, pos=pos,
        output=str(output_path),
        output_size=output_size,
        vertex_fill_color=vcolor,
        vertex_color=[0, 0, 0, 1],
        vertex_size=vertex_size,
        edge_color=ecolor,
        edge_pen_width=pen,
        eorder=draw_order,
        edge_gradient=[],
        edge_control_points=cts,
    )

    del g, weight, efilter, ecolor, vcolor, pen, t, tpos, cts, pos
    gc.collect()
