"""Network-drawing templates — layouts, edge rendering, node decoration.

Library home for the figure class introduced under
``.agents/guides/05_plotting/network_templates/``.  The shipped
templates (``single_network``, ``row_per_phase``, mosaics,
``mosaic_layout_compare``, ``lrg_seeded``) are thin scripts that
compose helpers from this module.

Three orthogonal layers, each independently swappable:

1. **Layouts** — pure ``(A, ...) → (N, 2)`` functions.  Eleven variants
   spanning NetworkX (``spring``, ``kk``, ``spectral``,
   ``laplacian_pca``, ``community_grouped``, ``backbone_guided``,
   ``circular_by_shaft``, ``lrg_kk``) and graph-tool (``sfdp``,
   ``arf``, ``lrg_sfdp``).
2. **Edge renderer** — :func:`draw_gamma_edges` with γ-power on width
   and alpha.  ``coloring="probe"`` (default for magnitude FC):
   same-probe red / cross-probe gray.  ``coloring="signed"`` for raw
   signed FC: diverging ``RdBu_r`` over signed weights.
3. **Node decoration** — :func:`shaft_colors`, :func:`community_colors`,
   :func:`draw_nodes`.

Style invariants (project-wide):

- Full vector PDFs — never call ``set_rasterized(True)``.
- No ``fig.suptitle``; no watermark by default.
- Same-probe / cross-probe edge contrast for any magnitude FC.

Three-layer composition:

>>> from lrg_eegfc.workflow.fc import load_fc_matrix
>>> from lrg_eegfc.visuals.network_templates import (
...     layout_kk, draw_gamma_edges, draw_nodes, shaft_colors,
...     load_probe_labels,
... )
>>> A = load_fc_matrix("Pat_05", "rest_pre", "beta", "imcoh_abs")
>>> probes = load_probe_labels("Pat_05")
>>> pos = layout_kk(A)
>>> fig, ax = plt.subplots(figsize=(5, 5))
>>> draw_gamma_edges(ax, pos, A, coloring="probe", probes=probes)
>>> draw_nodes(ax, pos, c=shaft_colors(probes))
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, List, Literal, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
from numpy.typing import NDArray
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.config.const import (
    BRAIN_BAND_TEX_DICT,
    BRAIN_BANDS_NAMES,
    PHASE_LABELS,
)
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.fc.msc.sparsify import disparity_filter
from lrg_eegfc.utils.io.patient import load_patient_metadata
from lrg_eegfc.utils.probe import extract_probe_labels
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result

# Layout + chord helpers moved out on 2026-05-29 (Phase 4-B split 7/7).
# Imported at module top so functions in this file can use them as default
# argument values (``DEFAULT_GALLERY``, etc.) at definition time. A full
# re-export block at the bottom of this file extends backwards compatibility
# to the rest of the symbol surface.
from .network_layouts import (
    compute_layout,
    DEFAULT_GALLERY,
    is_lrg_layout,
    LAYOUT_REGISTRY,
    layout_arf,
    layout_backbone_guided,
    layout_circular_by_shaft,
    layout_community_grouped,
    layout_kk,
    layout_laplacian_pca,
    layout_lrg_kk,
    layout_lrg_sfdp,
    layout_sfdp,
    layout_spectral,
    layout_spring,
)
from .network_chord import (
    DENDROGRAM_DEFAULT_COLOR,
    DENDROGRAM_DEFAULT_LW,
    DENDROGRAM_R_INNER,
    DENDROGRAM_R_OUTER,
    HIERARCHY_BG_ALPHA,
    HIERARCHY_BG_PEN,
    HIERARCHY_BG_RGB,
    HIERARCHY_DEFAULT_BETA,
    HIERARCHY_DEFAULT_FIT_VIEW,
    HIERARCHY_DEFAULT_K_CLUSTERS,
    HIERARCHY_DEFAULT_RENDER_PX,
    HIERARCHY_R_CLUSTER,
    HIERARCHY_R_OUTER,
    build_chord_depth2_layout,
    draw_circular_dendrogram_overlay,
    draw_radial_leaf_labels,
    plot_chord_with_dendrogram,
    plot_chord_with_highlight,
    render_hierarchy_chord,
)


__all__ = [
    # Constants
    "SAME_PROBE_RGB",
    "CROSS_PROBE_RGB",
    "DEFAULT_GAMMA",
    "DEFAULT_WIDTH_RANGE",
    "DEFAULT_ALPHA_RANGE",
    "SIGNED_CMAP",
    "DEFAULT_NODE_SIZE",
    "DEFAULT_LAYOUT",
    "DEFAULT_COLORING",
    "MAGNITUDE_FC_METHODS",
    "SIGNED_FC_METHODS",
    # Adaptors / loaders
    "load_probe_labels",
    "nx_to_gt",
    "matrix_to_gt",
    # Layouts
    "layout_spring",
    "layout_kk",
    "layout_spectral",
    "layout_laplacian_pca",
    "layout_community_grouped",
    "layout_backbone_guided",
    "layout_circular_by_shaft",
    "layout_lrg_kk",
    "layout_sfdp",
    "layout_arf",
    "layout_lrg_sfdp",
    "LAYOUT_REGISTRY",
    "DEFAULT_GALLERY",
    "is_lrg_layout",
    "compute_layout",
    # Preprocessing
    "rank_transform",
    "power_transform",
    "disparity_backbone",
    "top_k_backbone",
    # Edge renderer
    "draw_gamma_edges",
    # Node decoration
    "shaft_colors",
    "community_colors",
    "draw_nodes",
    # Top-level templates
    "plot_fc_network",
    "plot_fc_network_row",
    "plot_fc_network_grid",
    "plot_fc_network_lrg",
    "plot_layout_gallery",
    # Hierarchy-bundled chord (graph-tool curvy edges through LRG)
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
    # Matrix + network (cmap-linked)
    "spring_auto_k",
    "plot_fc_matrix_and_network",
    "plot_fc_matrix_and_network_rows",
]


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SAME_PROBE_RGB: Tuple[float, float, float] = (0.8, 0.2, 0.2)
CROSS_PROBE_RGB: Tuple[float, float, float] = (0.35, 0.35, 0.35)

DEFAULT_GAMMA: float = 6.0
DEFAULT_WIDTH_RANGE: Tuple[float, float] = (0.15, 4.0)
DEFAULT_ALPHA_RANGE: Tuple[float, float] = (0.03, 0.9)

SIGNED_CMAP: str = "RdBu_r"

DEFAULT_NODE_SIZE: int = 18
DEFAULT_LAYOUT: str = "spring"
DEFAULT_COLORING: str = "probe"

# Magnitude FC methods → ``coloring="probe"`` default; signed methods →
# ``coloring="signed"`` (diverging cmap, no probe red/gray).
MAGNITUDE_FC_METHODS: Tuple[str, ...] = ("imcoh_abs", "imcoh_sq", "msc")
SIGNED_FC_METHODS: Tuple[str, ...] = ("corr", "imcoh")

ColoringMode = Literal["probe", "signed", "weight", "shaft", "cmap"]
NodeColorMode = Literal["shaft", "community", "uniform"]


def _shaft_color_map(probes: Sequence[str]):
    """probe id → tab20 RGBA, matching :func:`shaft_colors`."""
    uniq = sorted(set(probes))
    cmap = plt.get_cmap("tab20", max(len(uniq), 3))
    return {p: cmap(i % 20) for i, p in enumerate(uniq)}


# ---------------------------------------------------------------------------
# Loaders / adaptors
# ---------------------------------------------------------------------------

def load_probe_labels(patient: str, *, dataset_root: Path = SEEG_DATAPATH) -> List[str]:
    """Return probe identifiers (one per channel) for ``patient``.

    Reads ``channel_labels.csv`` via :func:`load_patient_metadata`,
    applies any per-patient channel drops, and runs
    :func:`extract_probe_labels` to convert raw labels (``"A 1,G2"``,
    ``"G' 3"``, …) to probe ids (``"A"``, ``"G'"``).
    """
    meta = load_patient_metadata(patient, dataset_root)
    if meta is None:
        raise RuntimeError(
            f"No channel/implant metadata for {patient} under {dataset_root}"
        )
    return extract_probe_labels(meta["label"].tolist())


def matrix_to_gt(A: NDArray):
    """Build a ``graph_tool.Graph`` from an adjacency matrix.

    Returns ``(g, edge_weight_property)``.  Only the upper triangle is
    iterated (graphs are undirected).  Edges with zero or negative
    weight are dropped — graph-tool layouts are written for non-negative
    edge weights.  For signed FC, threshold or magnitude-transform
    before calling.
    """
    import graph_tool.all as gt
    N = A.shape[0]
    g = gt.Graph(directed=False)
    g.add_vertex(N)
    ew = g.new_edge_property("double")
    r, c = np.triu_indices(N, k=1)
    mask = A[r, c] > 0
    if mask.any():
        g.add_edge_list(np.column_stack([r[mask], c[mask]]))
        ew.a = A[r[mask], c[mask]]
    return g, ew


def nx_to_gt(G: nx.Graph):
    """Convert a ``networkx.Graph`` to a ``graph_tool.Graph``.

    Edge weights are read from the ``"weight"`` attribute (default
    ``1.0``).  Returns ``(g, ew)`` like :func:`matrix_to_gt`.
    """
    N = G.number_of_nodes()
    A = nx.to_numpy_array(G, nodelist=range(N), weight="weight")
    return matrix_to_gt(A)


# ---------------------------------------------------------------------------
# Preprocessing — operate on A; return modified A.  Drawing always uses
# the *original* A so panels stay comparable across preprocessing
# variants.
# ---------------------------------------------------------------------------

def rank_transform(A: NDArray, alpha: float = 8.0) -> NDArray:
    """Replace upper-triangle weights with ``(rank / E) ** alpha``.

    Uniform-by-rank in ``(0, 1]``, then power-amplified.  Useful to
    redistribute mass for spring/KK on near-uniform weights (raw
    ImCoh).  Lower triangle mirrored, diagonal zero.
    """
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    active = w > 0
    if not active.any():
        return A.copy()
    w_pos = w[active]
    ranks = np.argsort(np.argsort(w_pos)).astype(float)
    ranks_norm = (ranks + 1) / len(ranks)
    new_w = np.zeros_like(w)
    new_w[active] = ranks_norm ** alpha
    B = np.zeros_like(A)
    B[r, c] = new_w
    B = B + B.T
    return B


def power_transform(A: NDArray, alpha: float = 2.0) -> NDArray:
    """Return ``(A / A.max()) ** alpha``.  Sharpens dynamic range
    without reordering.  Diagonal zeroed.
    """
    B = np.abs(A.copy())
    np.fill_diagonal(B, 0.0)
    m = B.max()
    if m <= 0:
        return B
    return (B / m) ** alpha


def disparity_backbone(A: NDArray, alpha: float = 0.05) -> NDArray:
    """Disparity-filter backbone (Serrano et al. 2009).

    Wraps :func:`lrg_eegfc.utils.fc.msc.sparsify.disparity_filter`.
    Keeps every locally-significant edge regardless of absolute
    magnitude.
    """
    return disparity_filter(np.abs(A), alpha=alpha)


def top_k_backbone(A: NDArray, keep_frac: float = 0.05) -> NDArray:
    """Keep the top ``keep_frac`` of upper-triangle edges by weight.

    Quick-and-dirty backbone; less principled than disparity filter.
    """
    B = np.abs(A.copy())
    np.fill_diagonal(B, 0.0)
    r, c = np.triu_indices(B.shape[0], k=1)
    w = B[r, c]
    if w.size == 0:
        return B
    n_keep = max(1, int(np.ceil(keep_frac * w.size)))
    threshold = np.partition(w, -n_keep)[-n_keep]
    M = np.zeros_like(B)
    keep = w >= threshold
    M[r[keep], c[keep]] = w[keep]
    M = M + M.T
    return M


# ---------------------------------------------------------------------------
# Edge renderer
# ---------------------------------------------------------------------------

def draw_gamma_edges(
    ax: plt.Axes,
    pos: NDArray,
    A: NDArray,
    *,
    coloring: ColoringMode = "probe",
    probes: Optional[Sequence[str]] = None,
    gamma: float = DEFAULT_GAMMA,
    width_range: Tuple[float, float] = DEFAULT_WIDTH_RANGE,
    alpha_range: Tuple[float, float] = DEFAULT_ALPHA_RANGE,
    t_mode: Literal["magnitude", "rank"] = "magnitude",
    same_probe_rgb: Tuple[float, float, float] = SAME_PROBE_RGB,
    cross_probe_rgb: Tuple[float, float, float] = CROSS_PROBE_RGB,
    same_probe_alpha: float = 0.9,
    signed_cmap: str = SIGNED_CMAP,
    signed_vrange: Optional[Tuple[float, float]] = None,
    cmap: Optional[object] = None,
    norm: Optional[object] = None,
    zorder: int = 1,
    min_alpha: float = 0.0,
) -> Optional["LineCollection"]:
    """Render edges with γ-power-scaled width and alpha.

    The recipe (canonical across the project):

    - ``t = (|w| / |w|_max) ** gamma``
    - ``width = width_range[0] + t · Δwidth``
    - ``alpha = alpha_range[0] + t · Δalpha``
    - sort ascending → strongest edges drawn last (on top of weak)

    Coloring modes:

    - ``"probe"`` (default for magnitude FC) — same-probe pairs
      red, cross-probe pairs gray.  Requires ``probes=``.  Reveals
      same-probe sEEG bias (probe-bias guide).
    - ``"shaft"`` — same-probe pairs painted in **the shaft's tab20
      colour** (matching the node colour produced by
      :func:`shaft_colors`); cross-probe pairs gray.  Requires
      ``probes=``.  Lets the reader see at a glance which shaft has
      dense intra-shaft coupling without losing the cross-probe
      backdrop.
    - ``"signed"`` — diverging ``RdBu_r`` colormap on signed weight,
      symmetric around 0.  ``signed_vrange`` overrides the default
      symmetric ``(-|w|_max, |w|_max)``.
    - ``"weight"`` — sequential gray scale by ``|w|``; no probe
      info.  Used for cohort comparisons where probe geometry differs.
    - ``"cmap"`` — arbitrary ``cmap + norm`` on ``|w|``.  The canonical
      "matrix-attached" coloring: pass the same ``cmap`` and ``norm``
      used by the partner ``imshow`` so the matrix and network share
      the colour scale.  Requires ``cmap=`` (string or
      ``matplotlib.colors.Colormap``); ``norm=`` falls back to
      :class:`matplotlib.colors.Normalize(vmin=0, vmax=|w|_max)`.

    Returns the ``LineCollection`` (so callers can attach a colorbar
    when ``coloring="signed"`` or ``"cmap"``).  ``None`` if no edges to draw.

    NEVER passes ``rasterized=True`` — vector PDFs handle ~6.5k thin
    lines fine and the project rule is "no rasterization".
    """
    A = np.asarray(A)
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w_signed = A[r, c]
    w_abs = np.abs(w_signed)
    active = np.where(w_abs > 0)[0]
    if active.size == 0:
        return None

    w_act_abs = w_abs[active]
    w_act_signed = w_signed[active]
    if t_mode == "rank":
        # Replace |w| with its normalized rank in (0, 1] before γ.
        # Result: t-distribution is uniform pre-γ, so γ acts on a
        # well-conditioned input regardless of how compressed the
        # raw weight distribution is.
        ranks = np.argsort(np.argsort(w_act_abs)).astype(float)
        t_pre = (ranks + 1) / len(ranks)
    elif t_mode == "magnitude":
        t_pre = w_act_abs / w_act_abs.max()
    else:
        raise ValueError(f"Unknown t_mode {t_mode!r}")
    t = t_pre ** gamma
    widths = width_range[0] + t * (width_range[1] - width_range[0])
    alphas = alpha_range[0] + t * (alpha_range[1] - alpha_range[0])
    # Global min_alpha filter applies to modes whose edge alpha tracks t
    # across the whole edge set (probe / signed / weight).  Shaft mode
    # overrides same-probe alpha to a fixed value so the global filter
    # would incorrectly drop same-probe edges based on their raw t-driven
    # alpha; shaft mode handles its own per-pass filtering below.
    if min_alpha > 0.0 and coloring != "shaft":
        keep = alphas >= min_alpha
        active = active[keep]
        w_act_abs = w_act_abs[keep]
        w_act_signed = w_act_signed[keep]
        t = t[keep]
        widths = widths[keep]
        alphas = alphas[keep]
        if active.size == 0:
            return None
    order = np.argsort(w_act_abs)

    segments = []
    linew: List[float] = []
    colors: List[Tuple[float, float, float, float]] = []

    if coloring == "probe":
        if probes is None:
            raise ValueError("coloring='probe' requires probes=")
        for oi in order:
            idx = active[oi]
            i, j = int(r[idx]), int(c[idx])
            segments.append([pos[i], pos[j]])
            linew.append(float(widths[oi]))
            rgb = same_probe_rgb if probes[i] == probes[j] else cross_probe_rgb
            colors.append((*rgb, float(alphas[oi])))

    elif coloring == "shaft":
        if probes is None:
            raise ValueError("coloring='shaft' requires probes=")
        shaft_map = _shaft_color_map(probes)
        # Two-pass bucketed rendering — keeps every edge in vector form
        # while collapsing the 6.5k strokes into ~50 LineCollections.
        # Each bucket emits ONE multi-subpath PDF stroke instead of N
        # individual strokes (10-20× file-size savings).
        #
        # Bucketing snaps each edge's (alpha, width) to a fixed
        # 32×16 grid; edges sharing a bucket are drawn with the
        # bucket's representative scalar alpha + scalar width.  The
        # visual difference is invisible to the eye (alphas snap to
        # ~3% bins, widths to ~6% bins) but the PDF byte cost
        # drops by an order of magnitude.
        #
        # Same-probe edges: width quantized, alpha FIXED — drawn ON TOP.
        # Cross-probe edges: both quantized; filtered by ``min_alpha``
        # (drops ghosts at the alpha-floor).
        N_ALPHA_BINS = 32
        N_WIDTH_BINS = 16
        a_lo, a_hi = alpha_range
        w_lo, w_hi = width_range

        def _bin_idx(values, n_bins, vlo, vhi):
            denom = (vhi - vlo) if (vhi - vlo) > 1e-12 else 1.0
            idx = np.floor((np.asarray(values) - vlo) / denom * n_bins).astype(int)
            return np.clip(idx, 0, n_bins - 1)

        def _bin_value(idx, n_bins, vlo, vhi):
            return vlo + (idx + 0.5) / n_bins * (vhi - vlo)

        # Split into same/cross, applying min_alpha filter to cross only.
        same_idx_list = []
        cross_idx_list = []
        for oi in order:
            idx = active[oi]
            i, j = int(r[idx]), int(c[idx])
            if probes[i] == probes[j]:
                same_idx_list.append((i, j, oi))
            else:
                if min_alpha > 0.0 and float(alphas[oi]) < min_alpha:
                    continue
                cross_idx_list.append((i, j, oi))

        # ---- cross-probe pass: bucket on (alpha_bin, width_bin) ----
        if cross_idx_list:
            ois = np.array([t[2] for t in cross_idx_list])
            a_bins = _bin_idx(alphas[ois], N_ALPHA_BINS, a_lo, a_hi)
            w_bins = _bin_idx(widths[ois], N_WIDTH_BINS, w_lo, w_hi)
            keys = a_bins * N_WIDTH_BINS + w_bins
            for key in np.unique(keys):
                mask = keys == key
                a_b = int(key // N_WIDTH_BINS)
                w_b = int(key % N_WIDTH_BINS)
                a_rep = float(_bin_value(a_b, N_ALPHA_BINS, a_lo, a_hi))
                w_rep = float(_bin_value(w_b, N_WIDTH_BINS, w_lo, w_hi))
                segs = [
                    [pos[cross_idx_list[m][0]], pos[cross_idx_list[m][1]]]
                    for m in np.where(mask)[0]
                ]
                ax.add_collection(LineCollection(
                    segs, colors=[(*cross_probe_rgb, a_rep)],
                    linewidths=w_rep, zorder=zorder,
                ))

        # ---- same-probe pass: bucket on (shaft_id, width_bin) ----
        if same_idx_list:
            ois = np.array([t[2] for t in same_idx_list])
            ii = np.array([t[0] for t in same_idx_list])
            shaft_ids = np.array([probes[k] for k in ii])
            w_bins = _bin_idx(widths[ois], N_WIDTH_BINS, w_lo, w_hi)
            for shaft_id in np.unique(shaft_ids):
                shaft_mask = shaft_ids == shaft_id
                rgba = shaft_map[shaft_id]
                shaft_w_bins = w_bins[shaft_mask]
                shaft_indices = np.where(shaft_mask)[0]
                for w_b in np.unique(shaft_w_bins):
                    sub_mask = shaft_w_bins == w_b
                    w_rep = float(_bin_value(int(w_b), N_WIDTH_BINS, w_lo, w_hi))
                    segs = [
                        [pos[same_idx_list[shaft_indices[m]][0]],
                         pos[same_idx_list[shaft_indices[m]][1]]]
                        for m in np.where(sub_mask)[0]
                    ]
                    ax.add_collection(LineCollection(
                        segs, colors=[(rgba[0], rgba[1], rgba[2], float(same_probe_alpha))],
                        linewidths=w_rep, zorder=zorder + 1,
                    ))

        return None  # multiple bucketed collections; no single ref

    elif coloring == "signed":
        cmap = plt.get_cmap(signed_cmap)
        if signed_vrange is None:
            vmax = float(w_act_abs.max())
            vmin = -vmax
        else:
            vmin, vmax = signed_vrange
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
        for oi in order:
            idx = active[oi]
            i, j = int(r[idx]), int(c[idx])
            segments.append([pos[i], pos[j]])
            linew.append(float(widths[oi]))
            rgba = cmap(norm(w_act_signed[oi]))
            colors.append((rgba[0], rgba[1], rgba[2], float(alphas[oi])))

    elif coloring == "weight":
        cmap_w = plt.get_cmap("Greys")
        for oi in order:
            idx = active[oi]
            i, j = int(r[idx]), int(c[idx])
            segments.append([pos[i], pos[j]])
            linew.append(float(widths[oi]))
            rgba = cmap_w(t[oi] * 0.7 + 0.3)  # avoid pure white
            colors.append((rgba[0], rgba[1], rgba[2], float(alphas[oi])))

    elif coloring == "cmap":
        if cmap is None:
            raise ValueError(
                "coloring='cmap' requires cmap= (string or "
                "matplotlib.colors.Colormap)"
            )
        from matplotlib.colors import Colormap, Normalize
        cmap_obj = cmap if isinstance(cmap, Colormap) else plt.get_cmap(cmap)
        if norm is None:
            norm = Normalize(vmin=0.0, vmax=float(w_act_abs.max()))
        for oi in order:
            idx = active[oi]
            i, j = int(r[idx]), int(c[idx])
            segments.append([pos[i], pos[j]])
            linew.append(float(widths[oi]))
            rgba = cmap_obj(norm(w_act_abs[oi]))
            colors.append((rgba[0], rgba[1], rgba[2], float(alphas[oi])))

    else:
        raise ValueError(f"Unknown coloring {coloring!r}")

    lc = LineCollection(segments, colors=colors, linewidths=linew, zorder=zorder)
    ax.add_collection(lc)
    return lc


# ---------------------------------------------------------------------------
# Node decoration
# ---------------------------------------------------------------------------

def shaft_colors(probes: Sequence[str]) -> List[Tuple[float, float, float, float]]:
    """tab20 colour per electrode shaft.  Repeats if > 20 unique probes."""
    uniq = sorted(set(probes))
    cmap = plt.get_cmap("tab20", max(len(uniq), 3))
    return [cmap(uniq.index(p) % 20) for p in probes]


def community_colors(labels: Sequence[int]) -> List[Tuple[float, float, float, float]]:
    """tab20 colour per LRG community id (1-indexed)."""
    arr = np.asarray(labels, dtype=int)
    K = max(int(arr.max()), 2)
    cmap = plt.get_cmap("tab20", max(K, 3))
    return [cmap((int(l) - 1) % 20) for l in arr]


def draw_nodes(
    ax: plt.Axes,
    pos: NDArray,
    *,
    c,
    s: int = DEFAULT_NODE_SIZE,
    edgecolors: str = "white",
    linewidths: float = 0.3,
    zorder: int = 5,
) -> None:
    ax.scatter(
        pos[:, 0], pos[:, 1],
        c=c, s=s,
        edgecolors=edgecolors, linewidths=linewidths,
        zorder=zorder,
    )


def _finalize_axes(ax: plt.Axes, pos: NDArray, *, margin: float = 0.08) -> None:
    """Standard limits + ``axis('off')`` for network panels."""
    ax.axis("off")
    if pos.size == 0:
        return
    xmin, xmax = pos[:, 0].min(), pos[:, 0].max()
    ymin, ymax = pos[:, 1].min(), pos[:, 1].max()
    span = max(xmax - xmin, ymax - ymin, 1e-9)
    ax.set_xlim(xmin - margin * span, xmax + margin * span)
    ax.set_ylim(ymin - margin * span, ymax + margin * span)


# ---------------------------------------------------------------------------
# Top-level templates
# ---------------------------------------------------------------------------

def _coloring_default_for(fc_method: str) -> ColoringMode:
    if fc_method in SIGNED_FC_METHODS:
        return "signed"
    return "probe"


def _node_color_array(
    *,
    node_color: NodeColorMode,
    probes: Sequence[str],
    community_labels: Optional[Sequence[int]] = None,
):
    if node_color == "shaft":
        return shaft_colors(probes)
    if node_color == "community":
        if community_labels is None:
            raise ValueError("node_color='community' requires community_labels=")
        return community_colors(community_labels)
    if node_color == "uniform":
        return ["lightblue"] * len(probes)
    raise ValueError(f"Unknown node_color {node_color!r}")


def plot_fc_network(
    A: NDArray,
    probes: Sequence[str],
    *,
    ax: Optional[plt.Axes] = None,
    layout: str = DEFAULT_LAYOUT,
    layout_kwargs: Optional[dict] = None,
    coloring: Optional[ColoringMode] = None,
    fc_method: Optional[str] = None,
    community_labels: Optional[Sequence[int]] = None,
    ultrametric_matrix: Optional[NDArray] = None,
    node_color: NodeColorMode = "shaft",
    node_size: int = DEFAULT_NODE_SIZE,
    gamma: float = DEFAULT_GAMMA,
    width_range: Tuple[float, float] = DEFAULT_WIDTH_RANGE,
    alpha_range: Tuple[float, float] = DEFAULT_ALPHA_RANGE,
    t_mode: Literal["magnitude", "rank"] = "magnitude",
    min_alpha: float = 0.0,
    title: Optional[str] = None,
    title_fontsize: int = 11,
) -> Tuple[plt.Axes, Optional["LineCollection"]]:
    """Single-panel network drawing — the atomic template.

    Composes the three layers: ``compute_layout`` →
    :func:`draw_gamma_edges` → :func:`draw_nodes`.  Returns
    ``(ax, edge_collection)``; the edge collection is non-None only
    when ``coloring="signed"`` so callers can attach a colorbar.

    ``coloring`` defaults to magnitude (``"probe"``) if ``fc_method``
    is in :data:`MAGNITUDE_FC_METHODS`, else ``"signed"``.  Pass
    explicit ``coloring=`` to override.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))

    if coloring is None:
        coloring = _coloring_default_for(fc_method) if fc_method else "probe"

    layout_kwargs = layout_kwargs or {}
    pos = compute_layout(
        layout, A,
        probes=probes,
        community_labels=community_labels,
        ultrametric_matrix=ultrametric_matrix,
        **layout_kwargs,
    )

    lc = draw_gamma_edges(
        ax, pos, A,
        coloring=coloring,
        probes=probes if coloring in ("probe", "shaft") else None,
        gamma=gamma, width_range=width_range, alpha_range=alpha_range,
        t_mode=t_mode,
        min_alpha=min_alpha,
    )

    nc = _node_color_array(
        node_color=node_color, probes=probes,
        community_labels=community_labels,
    )
    draw_nodes(ax, pos, c=nc, s=node_size)

    _finalize_axes(ax, pos)

    if title:
        ax.set_title(title, fontsize=title_fontsize, fontweight="bold")

    return ax, lc


def _load_inputs(
    patient: str,
    band: str,
    phase: str,
    fc_method: str,
    *,
    dataset_root: Path = SEEG_DATAPATH,
) -> Tuple[NDArray, List[str]]:
    A = load_fc_matrix(patient, phase, band, fc_method)
    if A is None:
        raise FileNotFoundError(
            f"No {fc_method} FC matrix for {patient} / {band} / {phase}"
        )
    probes = load_probe_labels(patient, dataset_root=dataset_root)
    N = min(len(probes), A.shape[0])
    A = A[:N, :N]
    probes = list(probes[:N])
    return np.asarray(A), probes


def plot_fc_network_row(
    patient: str,
    band: str,
    *,
    fc_method: str = "imcoh_abs",
    phases: Sequence[str] = PHASE_LABELS,
    bands: Optional[Sequence[str]] = None,
    layout: str = DEFAULT_LAYOUT,
    layout_kwargs: Optional[dict] = None,
    coloring: Optional[ColoringMode] = None,
    node_color: NodeColorMode = "shaft",
    node_size: int = 14,
    panel_size: float = 3.5,
    show_titles: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
) -> plt.Figure:
    """Row of network panels — one row × {phases | bands}.

    If ``bands`` is given, iterate over bands at fixed ``phase=phases[0]``
    (single-band case is meaningless).  Otherwise iterate over
    ``phases`` at fixed ``band``.
    """
    iterate_bands = bands is not None
    panels = list(bands) if iterate_bands else list(phases)
    n = len(panels)

    if figsize is None:
        figsize = (panel_size * n, panel_size)
    fig, axes = plt.subplots(1, n, figsize=figsize)
    if n == 1:
        axes = [axes]

    coloring_default = coloring or _coloring_default_for(fc_method)

    for k, label in enumerate(panels):
        ax = axes[k]
        if iterate_bands:
            cur_band, cur_phase = label, phases[0]
        else:
            cur_band, cur_phase = band, label
        try:
            A, probes = _load_inputs(patient, cur_band, cur_phase, fc_method)
        except FileNotFoundError as exc:
            ax.text(0.5, 0.5, f"missing\n{exc}", ha="center", va="center",
                    transform=ax.transAxes, fontsize=8, color="red")
            ax.axis("off")
            continue
        plot_fc_network(
            A, probes,
            ax=ax, layout=layout, layout_kwargs=layout_kwargs,
            coloring=coloring_default,
            fc_method=fc_method,
            node_color=node_color, node_size=node_size,
        )
        if show_titles:
            if iterate_bands:
                ttl = BRAIN_BAND_TEX_DICT.get(cur_band, cur_band)
            else:
                ttl = cur_phase.replace("_", "$_{") + "}$" if "_" in cur_phase else cur_phase
            ax.set_title(ttl, fontsize=11, fontweight="bold")

    fig.tight_layout()
    return fig


def plot_fc_network_grid(
    patient: str,
    *,
    fc_method: str = "imcoh_abs",
    phases: Sequence[str] = PHASE_LABELS,
    bands: Sequence[str] = tuple(BRAIN_BANDS_NAMES),
    layout: str = DEFAULT_LAYOUT,
    layout_kwargs: Optional[dict] = None,
    coloring: Optional[ColoringMode] = None,
    node_color: NodeColorMode = "shaft",
    node_size: int = 12,
    panel_size: float = 2.4,
    wspace: float = 0.05,
    hspace: float = 0.1,
    show_titles: bool = True,
) -> plt.Figure:
    """Mosaic of network panels — phases × bands for one patient."""
    n_rows = len(phases)
    n_cols = len(bands)
    figsize = (panel_size * n_cols, panel_size * n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_rows == 1:
        axes = np.array([axes])
    if n_cols == 1:
        axes = axes.reshape(-1, 1)

    coloring_default = coloring or _coloring_default_for(fc_method)

    for i, phase in enumerate(phases):
        for j, band in enumerate(bands):
            ax = axes[i, j]
            try:
                A, probes = _load_inputs(patient, band, phase, fc_method)
            except FileNotFoundError as exc:
                ax.text(0.5, 0.5, f"missing", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8, color="red")
                ax.axis("off")
                continue
            plot_fc_network(
                A, probes,
                ax=ax, layout=layout, layout_kwargs=layout_kwargs,
                coloring=coloring_default,
                fc_method=fc_method,
                node_color=node_color, node_size=node_size,
            )
            if show_titles and i == 0:
                ax.set_title(BRAIN_BAND_TEX_DICT.get(band, band),
                             fontsize=11, fontweight="bold")
            if show_titles and j == 0:
                ax.text(-0.06, 0.5, phase.replace("_", "_"),
                        rotation=90, transform=ax.transAxes,
                        ha="right", va="center",
                        fontsize=10, fontweight="bold")

    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    return fig


def plot_layout_gallery(
    patient: str,
    band: str,
    phase: str,
    *,
    fc_method: str = "imcoh_abs",
    layouts: Sequence[str] = DEFAULT_GALLERY,
    n_communities: int = 10,
    coloring: Optional[ColoringMode] = None,
    node_color: NodeColorMode = "shaft",
    node_size: int = 16,
    panel_size: float = 3.0,
    n_cols: int = 4,
    show_titles: bool = True,
) -> plt.Figure:
    """One panel per layout — the "test all layouts" mosaic.

    Skips LRG layouts gracefully if no LRG cache is present for this
    (patient, band, phase, fc_method).
    """
    A, probes = _load_inputs(patient, band, phase, fc_method)

    # Try to fetch LRG community labels for the LRG-seeded layouts.
    community_labels: Optional[NDArray] = None
    try:
        lrg = load_lrg_result(patient, phase, band, fc_method)
        if lrg is not None and lrg.linkage_matrix is not None:
            N = A.shape[0]
            community_labels = fcluster(
                lrg.linkage_matrix, t=n_communities, criterion="maxclust",
            )[:N]
            ultrametric_matrix = getattr(lrg, "ultrametric_matrix", None)
        else:
            ultrametric_matrix = None
    except Exception:
        community_labels = None
        ultrametric_matrix = None

    layouts = list(layouts)
    if community_labels is None:
        layouts = [name for name in layouts if not is_lrg_layout(name)]

    n = len(layouts)
    n_rows = (n + n_cols - 1) // n_cols
    figsize = (panel_size * n_cols, panel_size * n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = np.atleast_2d(axes).flatten()

    coloring_default = coloring or _coloring_default_for(fc_method)

    for k, name in enumerate(layouts):
        ax = axes[k]
        try:
            pos = compute_layout(
                name, A, probes=probes,
                community_labels=community_labels,
                ultrametric_matrix=ultrametric_matrix,
            )
        except Exception as exc:
            ax.text(0.5, 0.5, f"{name}\nFAILED:\n{exc}",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=9, color="red")
            ax.axis("off")
            continue
        draw_gamma_edges(
            ax, pos, A,
            coloring=coloring_default,
            probes=probes if coloring_default in ("probe", "shaft") else None,
        )
        nc = _node_color_array(
            node_color=node_color, probes=probes,
            community_labels=community_labels,
        )
        draw_nodes(ax, pos, c=nc, s=node_size)
        _finalize_axes(ax, pos)
        if show_titles:
            ax.set_title(name, fontsize=11, fontweight="bold")

    for k in range(n, n_rows * n_cols):
        axes[k].axis("off")

    fig.tight_layout()
    return fig


def plot_fc_network_lrg(
    patient: str,
    band: str,
    phase: str,
    *,
    fc_method: str = "imcoh_abs",
    n_communities: int = 10,
    backend: Literal["nx", "gt"] = "gt",
    separation: float = 5.0,
    sfdp_gamma: float = 0.1,
    coloring: Optional[ColoringMode] = None,
    node_color: NodeColorMode = "community",
    node_size: int = 24,
    figsize: Tuple[float, float] = (6, 6),
) -> plt.Figure:
    """LRG-seeded network drawing — the cluster figure.

    ``backend="gt"`` → :func:`layout_lrg_sfdp` (graph-tool, fast,
    scales to thousands of nodes).
    ``backend="nx"`` → :func:`layout_lrg_kk` (NetworkX KK, no extra
    dep, slower at N≈115).

    Default ``node_color="community"`` because the figure's whole point
    is to show LRG modules.  Switch to ``"shaft"`` to verify the
    modules aren't probe-trivial.
    """
    A, probes = _load_inputs(patient, band, phase, fc_method)
    lrg = load_lrg_result(patient, phase, band, fc_method)
    if lrg is None or lrg.linkage_matrix is None:
        raise FileNotFoundError(
            f"No LRG result cached for {patient} / {band} / {phase} / {fc_method}"
        )
    N = A.shape[0]
    community_labels = fcluster(
        lrg.linkage_matrix, t=n_communities, criterion="maxclust",
    )[:N]

    fig, ax = plt.subplots(figsize=figsize)
    layout_name = "lrg_sfdp" if backend == "gt" else "lrg_kk"
    layout_kwargs = {"gamma": sfdp_gamma} if backend == "gt" else {}
    ultrametric_matrix = getattr(lrg, "ultrametric_matrix", None)
    plot_fc_network(
        A, probes,
        ax=ax, layout=layout_name,
        layout_kwargs=layout_kwargs,
        coloring=coloring,
        fc_method=fc_method,
        community_labels=community_labels,
        ultrametric_matrix=ultrametric_matrix,
        node_color=node_color,
        node_size=node_size,
    )
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Matrix + network (cmap-linked) — Pat_02/fig1_msc_alpha-style figure
# ---------------------------------------------------------------------------
#
# Composes :func:`lrg_eegfc.visuals.fc_templates.plot_fc_adjacency` on
# the left with a network drawing on the right whose edges are coloured
# by the SAME ``cmap + norm`` as the imshow.  Reader's eye reads both
# panels in the same colour language: the dark/bright FC cells map
# 1-to-1 onto the dark/bright edges.
#
# Modular: layout, edge γ-power recipe, coloring mode, backend
# (NetworkX / graph-tool), and tick-label mode are independent knobs.

def spring_auto_k(
    N: int,
    *,
    scale: Literal["msc_era", "kbase_sqrtN", "fixed", "inv_sqrtN"] = "msc_era",
    k_base: float = 5.0,
    k_fixed: float = 0.1,
) -> float:
    """Spring ``k`` chooser for dense weighted FC graphs.

    Four rules:

    - ``msc_era`` (default) — **``k = 0.1`` regardless of ``N``**.
      Matches the legacy ``plot_msc_and_network`` recipe that produced
      the structured MSC presentation figures
      (``data/outputs/figures/presentation_figures/Pat_02/fig1_msc_*``).
      Paired with ``iterations ≈ 50`` (NOT 600) the strong heavy-
      tailed edges pull tight clusters together before the spring
      relaxation equilibrates everything to a ball.  This is the
      load-bearing trick for matrix+network figures on
      heavy-tailed FC.
    - ``kbase_sqrtN`` — ``k = k_base / √N`` (FR canonical).  Reads
      gracefully across N=80–130 with **fully-converged** iterations,
      but on heavy-tailed dense FC it equilibrates the strong-edge
      structure into a ball.  Use only when you want a converged
      layout; combine with ``iterations ≥ 600``.
    - ``inv_sqrtN`` — ``k = 1 / √N``.  NetworkX's documented default.
    - ``fixed`` — ``k = k_fixed``.  Use when sweeping ``k`` for a
      diagnostic.

    Returns the raw ``k`` value to pass to ``layout_spring(k=...)``.
    """
    if scale == "msc_era":
        return 0.1
    if scale == "kbase_sqrtN":
        return float(k_base) / float(np.sqrt(N))
    if scale == "inv_sqrtN":
        return 1.0 / float(np.sqrt(N))
    if scale == "fixed":
        return float(k_fixed)
    raise ValueError(f"Unknown spring k scale: {scale!r}")


def _fc_method_cmap_default(fc_method: Optional[str]) -> str:
    """Default colormap for the matrix + network panel.

    Magnitude FC (``imcoh_abs``, ``imcoh_sq``, ``msc``) → ``magma``;
    signed FC (``corr``, ``imcoh``) → ``RdBu_r`` (diverging).
    """
    if fc_method in SIGNED_FC_METHODS:
        return SIGNED_CMAP
    return "magma"


def plot_fc_matrix_and_network(
    A: NDArray,
    probes: Sequence[str],
    *,
    ax_matrix: plt.Axes,
    ax_network: plt.Axes,
    fc_method: Optional[str] = None,
    band: Optional[str] = None,
    cmap: Optional[str] = None,
    vmin: Optional[float] = 0.0,
    vmax: Optional[float] = None,
    log_scale: bool = False,
    tick_labels: Literal["generic", "index", "chnames"] = "generic",
    channel_labels: Optional[Sequence[str]] = None,
    matrix_colorbar: bool = True,
    layout: str = "spring",
    layout_kwargs: Optional[Dict] = None,
    pos_override: Optional[NDArray] = None,
    coloring: ColoringMode = "cmap",
    node_color: NodeColorMode = "shaft",
    community_labels: Optional[NDArray] = None,
    ultrametric_matrix: Optional[NDArray] = None,
    node_size: int = DEFAULT_NODE_SIZE,
    network_colorbar: bool = False,
    gamma: float = 2.0,
    width_range: Tuple[float, float] = (0.0, 4.0),
    alpha_range: Tuple[float, float] = (0.0, 1.0),
    min_alpha: float = 0.0,
    t_mode: Literal["magnitude", "rank"] = "magnitude",
    finalize_network: bool = True,
) -> Dict:
    """One row of ``(matrix, network)`` with cmap-linked edges.

    The network's edges are coloured by the same ``cmap + norm`` as the
    matrix imshow when ``coloring="cmap"`` (the default and the figure's
    whole point).

    Parameters
    ----------
    A : ``(N, N)`` FC matrix.
    probes : per-contact shaft labels.
    ax_matrix, ax_network : pre-allocated axes the caller owns.
    fc_method, band : forwarded to ``plot_fc_adjacency`` for the
        canonical colourbar label.
    cmap : colormap shared by the matrix and the network edges.  When
        ``None`` → :func:`_fc_method_cmap_default`.
    vmin / vmax / log_scale : matrix ``imshow`` controls.  The network
        edges use the *same* numeric range so the two panels share a
        colour scale.
    layout : any name in :data:`LAYOUT_REGISTRY` (``spring``, ``kk``,
        ``spectral``, ``laplacian_pca``, ``lrg_kk``, ``lrg_sfdp``, …).
    layout_kwargs : forwarded to the layout function (e.g.
        ``{"k": 0.4, "iterations": 1000}`` for spring).
    coloring : one of ``probe / shaft / signed / weight / cmap``.
        Default ``cmap`` — paint edges by ``cmap(norm(|w|))``, locked
        to the matrix.
    node_color : ``shaft / community / uniform``.  Default ``shaft``.
    network_colorbar : if True, add a colorbar to ``ax_network`` (only
        useful when ``coloring="cmap"`` and the matrix colorbar is
        suppressed via ``matrix_colorbar=False``).
    finalize_network : if True (default), call ``ax_network.axis('off')``
        and zoom to the layout bbox with a uniform margin.

    Returns
    -------
    dict with keys ``pos``, ``im_matrix``, ``cmap``, ``norm``,
    ``edge_collection``, ``n_edges``, ``n_nodes``.
    """
    from matplotlib.colors import LogNorm, Normalize

    from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency

    A = np.asarray(A, dtype=float)
    N = A.shape[0]

    if cmap is None:
        cmap = _fc_method_cmap_default(fc_method)

    A_for_imshow = A.copy()
    np.fill_diagonal(A_for_imshow, 0.0)
    if vmax is None:
        vmax = float(np.nanmax(A_for_imshow)) if np.any(np.isfinite(A_for_imshow)) else 1.0

    _, _, im = plot_fc_adjacency(
        A,
        ax=ax_matrix,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        log_scale=log_scale,
        tick_labels=tick_labels,
        channel_labels=channel_labels,
        colorbar=matrix_colorbar,
        fc_method=fc_method,
        band=band,
    )

    # Lock the network's edge norm to the matrix's norm.
    if log_scale:
        positives = A_for_imshow[A_for_imshow > 0]
        log_vmin = float(positives.min()) if positives.size else 1e-6
        if vmin is not None and vmin > 0:
            log_vmin = max(log_vmin, float(vmin))
        norm = LogNorm(vmin=log_vmin, vmax=float(vmax))
    else:
        norm = Normalize(vmin=float(vmin or 0.0), vmax=float(vmax))

    if pos_override is not None:
        pos = np.asarray(pos_override)
    else:
        pos = compute_layout(
            layout, A_for_imshow,
            probes=probes,
            community_labels=community_labels,
            ultrametric_matrix=ultrametric_matrix,
            **(layout_kwargs or {}),
        )

    lc = draw_gamma_edges(
        ax_network, pos, A_for_imshow,
        coloring=coloring,
        probes=probes,
        cmap=cmap,
        norm=norm,
        gamma=gamma,
        width_range=width_range,
        alpha_range=alpha_range,
        t_mode=t_mode,
        min_alpha=min_alpha,
    )

    nc = _node_color_array(
        node_color=node_color,
        probes=probes,
        community_labels=community_labels,
    )
    draw_nodes(ax_network, pos, c=nc, s=node_size)

    if finalize_network:
        _finalize_axes(ax_network, pos)

    if network_colorbar and coloring == "cmap":
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax_network)
        cax = divider.append_axes("right", size="4%", pad=0.05)
        cb = plt.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=6)

    n_edges = int((np.triu(A_for_imshow, k=1) > 0).sum())
    return dict(
        pos=pos,
        im_matrix=im,
        cmap=cmap,
        norm=norm,
        edge_collection=lc,
        n_edges=n_edges,
        n_nodes=N,
    )


def plot_fc_matrix_and_network_rows(
    patient: str,
    band: str,
    *,
    phases: Sequence[str] = PHASE_LABELS,
    fc_method: str = "imcoh_abs",
    cmap: Optional[str] = None,
    vmin: Optional[float] = 0.0,
    vmax: Optional[float] = None,
    log_scale: bool = False,
    tick_labels: Literal["generic", "index", "chnames"] = "generic",
    layout: str = "spring",
    layout_kwargs: Optional[Dict] = None,
    shared_layout: bool = True,
    coloring: ColoringMode = "cmap",
    node_color: NodeColorMode = "shaft",
    node_size: int = DEFAULT_NODE_SIZE,
    n_communities: Optional[int] = None,
    gamma: float = 2.0,
    width_range: Tuple[float, float] = (0.0, 4.0),
    alpha_range: Tuple[float, float] = (0.0, 1.0),
    min_alpha: float = 0.0,
    figsize_per_row: Tuple[float, float] = (10.5, 4.4),
    shared_scale: bool = True,
) -> Tuple[plt.Figure, NDArray]:
    """Row-per-phase mosaic: one ``(matrix, network)`` row per phase.

    Mirrors the layout of the reference figure
    ``presentation_figures/Pat_02/fig1_msc_alpha.pdf`` (two phases) but
    generalised to any subset of ``PHASE_LABELS``.

    Defaults reproduce the MSC-presentation recipe that made the
    structure visible in the legacy figure:

    - ``layout="spring"`` with ``k=0.1`` / ``iterations=100`` (set by
      the calling script via ``layout_kwargs``).
    - ``gamma=2``  →  ``widths = max_width · t²``,
                       ``alphas = alpha_max · t²``  with
                       ``t = |w|/|w|_max``.
    - ``width_range=(0, 4)``, ``alpha_range=(0, 1)`` — full range
      starting from zero so weak edges become invisible and only the
      structural strong-edge skeleton survives visually.
    - ``shared_layout=True`` — layout computed once on the cross-phase
      AVERAGE FC, then every phase row uses ``shared_pos``.  Node
      positions stay identical across rows so side-by-side comparison
      is meaningful; only edges (widths + colours) change between
      phases.

    ``shared_scale=True`` (default) → cohort-wide ``vmax`` across the
    selected phases so the colorbar range is directly comparable.

    Returns ``(fig, axes)`` with ``axes`` an ``(n_phases, 2)`` array
    (column 0 = matrix, column 1 = network).
    """
    A_per_phase = {
        ph: load_fc_matrix(patient=patient, phase=ph,
                            band=band, fc_method=fc_method)
        for ph in phases
    }
    probes = load_probe_labels(patient)

    community_labels = None
    ultrametric_matrix = None
    if is_lrg_layout(layout) or node_color == "community":
        anchor_phase = "rest_pre" if "rest_pre" in phases else phases[0]
        lrg = load_lrg_result(patient, anchor_phase, band, fc_method)
        if lrg is None or lrg.linkage_matrix is None:
            raise FileNotFoundError(
                f"No LRG cache for {patient}/{band}/{anchor_phase}/"
                f"{fc_method}; needed for layout={layout!r} or "
                f"node_color='community'."
            )
        N_lrg = next(iter(A_per_phase.values())).shape[0]
        community_labels = fcluster(
            lrg.linkage_matrix,
            t=n_communities or 10,
            criterion="maxclust",
        )[:N_lrg]
        ultrametric_matrix = getattr(lrg, "ultrametric_matrix", None)

    if shared_scale and vmax is None:
        vmax_candidates = []
        for A in A_per_phase.values():
            Az = A.copy()
            np.fill_diagonal(Az, 0.0)
            vmax_candidates.append(float(np.nanmax(Az)))
        vmax = max(vmax_candidates)

    # ---- shared layout from the cross-phase average FC ----
    # The load-bearing trick of the MSC-era presentation figure:
    # compute the spring layout ONCE on the mean of all selected
    # phases, then re-use those positions for every row.  Without this
    # each row's spring relaxation lands in a different local minimum
    # and same-node comparison across rows is meaningless.
    pos_shared = None
    if shared_layout and len(phases) > 1:
        A_avg = np.mean(
            [A_per_phase[p] for p in phases], axis=0,
        )
        np.fill_diagonal(A_avg, 0.0)
        pos_shared = compute_layout(
            layout, A_avg,
            probes=probes,
            community_labels=community_labels,
            ultrametric_matrix=ultrametric_matrix,
            **(layout_kwargs or {}),
        )

    n_rows = len(phases)
    fig, axes = plt.subplots(
        n_rows, 2,
        figsize=(figsize_per_row[0], figsize_per_row[1] * n_rows),
        gridspec_kw={"width_ratios": [1.0, 1.1]},
    )
    if n_rows == 1:
        axes = np.atleast_2d(axes)

    for row, phase in enumerate(phases):
        ax_mat = axes[row, 0]
        ax_net = axes[row, 1]
        plot_fc_matrix_and_network(
            A_per_phase[phase], probes,
            ax_matrix=ax_mat, ax_network=ax_net,
            fc_method=fc_method, band=band,
            cmap=cmap, vmin=vmin,
            vmax=vmax if shared_scale else None,
            log_scale=log_scale,
            tick_labels=tick_labels,
            layout=layout, layout_kwargs=layout_kwargs,
            pos_override=pos_shared,
            coloring=coloring, node_color=node_color,
            community_labels=community_labels,
            ultrametric_matrix=ultrametric_matrix,
            node_size=node_size, gamma=gamma,
            width_range=width_range, alpha_range=alpha_range,
            min_alpha=min_alpha,
        )

    return fig, axes
