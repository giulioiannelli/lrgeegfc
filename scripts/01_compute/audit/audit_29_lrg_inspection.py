#!/usr/bin/env python3
"""Audit 29 — LRG inspection (visual ground-truth before any measure).

Plot, per (patient, band), the LRG primitives at τ = 1/λ_max for the
four phases side-by-side. Each phase uses its OWN natural dendrogram
leaf order (no forced cross-phase reordering). Cluster the dendrogram
at a fixed target k (default 20) — Ψ-based cuts were retired because
in our continuous-spectrum case all merge-height gaps are comparable
and Ψ collapses to either trivial (k=2) or all-singletons.

Layout (Family A, per (patient, band)):

    Row 1: ρ̂(τ_max) heatmap                       (linear scale)
    Row 2: 𝒟(τ_max) = (1 - δ_ij)/ρ̂(τ_max) heatmap (raw, NOT cophenetic)
    Row 3: dendrogram (own natural order, fixed-k cut line)
    Row 4: codebar — cluster id at fixed-k cut, tight under dendrogram

Family B: per-patient overview, 6 bands × 4 phases (dendrogram + codebar).

Cluster ids are re-ranked by size (1 = biggest cluster, 2 = next, …)
so colours mean the same thing across phases.

Usage
-----
    python audit_29_lrg_inspection.py                  # smoke test (Pat_06 alpha)
    python audit_29_lrg_inspection.py --patient Pat_06
    python audit_29_lrg_inspection.py --cohort
    python audit_29_lrg_inspection.py --cohort --k 30
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm, to_hex
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.ticker import LogFormatterMathtext, LogLocator
from scipy.cluster.hierarchy import dendrogram, fcluster

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BAND_TEX_DICT,
    BRAIN_BANDS_NAMES,
    PATIENTS_4PHASE,
    PHASE_SHORT_LABELS,
)
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, SEEG_DATAPATH
from lrg_eegfc.visuals.fc_templates import probe_groups
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrgsglib.plotlib.colorbars import imshow_colorbar_caxdivider

PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
FC_METHOD = "imcoh_abs"
OUT_ROOT = ROOT / "data" / "audit" / "lrg_inspection"

# Width-ratio padding so the row colorbar steals from the OUTER margin,
# not from the rightmost imshow (mirrors `fc_templates._ROW_CB_FRAC`).
_ROW_CB_FRAC = 0.06


def _style_log_colorbar(cb) -> None:
    """Decade major ticks (decade-only labels) + smaller minor ticks
    that are explicitly UNLABELED.

    Matplotlib's defaults on a narrow cax sometimes render minor
    tick labels (``2×10ⁿ``, ``3×10ⁿ``, …) that overflow the cax
    width and appear truncated. Force ``labelOnlyBase=True`` on the
    major formatter and a ``NullFormatter`` on the minor axis so
    only the integer decades carry labels. Also set tick lengths so
    major > minor.
    """
    ax = cb.ax
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.yaxis.set_major_formatter(
        LogFormatterMathtext(base=10.0, labelOnlyBase=True)
    )
    ax.yaxis.set_minor_locator(
        LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=200)
    )
    ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.tick_params(which="major", length=3.0, width=0.7, labelsize=6)
    ax.tick_params(which="minor", length=1.2, width=0.5)


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------


def compute_rho(eigenvalues: np.ndarray, eigenvectors: np.ndarray) -> np.ndarray:
    """ρ̂(τ) at τ = 1/λ_max from cached spectrum."""
    lam = np.asarray(eigenvalues, dtype=float)
    lam_max = float(lam.max())
    if lam_max <= 0:
        raise ValueError("λ_max ≤ 0 — degenerate Laplacian spectrum")
    tau = 1.0 / lam_max
    weights = np.exp(-tau * lam)
    Z = float(weights.sum())
    rho = (eigenvectors * weights[None, :]) @ eigenvectors.T
    rho /= Z
    return 0.5 * (rho + rho.T)


def compute_d_tau(rho: np.ndarray) -> np.ndarray:
    """𝒟(τ) = (1 - δ_ij)/ρ̂_ij(τ). Off-diagonal raw communication distance."""
    with np.errstate(divide="ignore", invalid="ignore"):
        D = 1.0 / rho
    np.fill_diagonal(D, 0.0)
    return D


def fixed_k_cut_height(linkage_matrix: np.ndarray, k_target: int) -> float:
    """Cophenetic height that yields exactly ``k_target`` flat clusters.

    With ``N`` leaves and ``N − 1`` merges sorted ascending, ``k_target``
    clusters are produced by cutting strictly above merge index
    ``N − k_target − 1`` and strictly below merge index ``N − k_target``.
    Returns the midpoint between those two heights for visual placement.
    """
    heights = np.sort(np.asarray(linkage_matrix[:, 2], dtype=float))
    N = linkage_matrix.shape[0] + 1
    k = max(1, min(int(k_target), N))
    if k >= N:
        return float(heights[0]) * 0.5
    if k == 1:
        return float(heights[-1]) * 1.05
    return 0.5 * (float(heights[N - k - 1]) + float(heights[N - k]))


def rank_clusters_by_size(cluster_ids: np.ndarray) -> np.ndarray:
    """Renumber clusters so id 1 = largest, 2 = second-largest, etc.

    Stabilises palette assignment across phases: when the colormap
    indexes by ``rank - 1``, the biggest cluster is always palette[0],
    the next palette[1], and so on — so a viewer comparing four
    phases sees consistent colour semantics.
    """
    unique, counts = np.unique(cluster_ids, return_counts=True)
    # Descending by size, deterministic tie-break by original id.
    order = np.lexsort((unique, -counts))
    remap = {int(old): rank for rank, old in enumerate(unique[order], start=1)}
    return np.asarray([remap[int(c)] for c in cluster_ids], dtype=int)


# ---------------------------------------------------------------------------
# Cluster coloring (matched between dendrogram links and codebar leaves)
# ---------------------------------------------------------------------------


def _cluster_palette(n_clusters: int) -> List[Tuple[float, float, float, float]]:
    """Distinct colour per cluster up to ``n_clusters`` (no cycling).

    - ``n ≤ 10``  → ``tab10`` (most discriminable for small palettes).
    - ``n ≤ 20``  → ``tab20``.
    - ``n  > 20`` → evenly-spaced hues from ``hsv``, so palette length
      always matches ``n``.
    """
    n = max(int(n_clusters), 1)
    if n <= 10:
        cmap = plt.get_cmap("tab10")
        return [cmap(i % 10) for i in range(n)]
    if n <= 20:
        cmap = plt.get_cmap("tab20")
        return [cmap(i % 20) for i in range(n)]
    cmap = plt.get_cmap("hsv", n)
    return [cmap(i) for i in range(n)]


def _first_merge_heights(Z: np.ndarray, n_leaves: int) -> np.ndarray:
    """For each leaf, the height at which it joins its first parent merge.

    Returns an array of length ``n_leaves``. Used to overdraw leaf
    stems with the leaf's own cluster colour (scipy paints stems with
    the parent merge's colour, which is grey for cross-cluster merges).
    """
    out = np.full(n_leaves, np.nan, dtype=float)
    for i in range(Z.shape[0]):
        for child in (int(Z[i, 0]), int(Z[i, 1])):
            if child < n_leaves and np.isnan(out[child]):
                out[child] = float(Z[i, 2])
    return out


def lca_node_matrix(Z: np.ndarray, n_leaves: int) -> np.ndarray:
    """Lowest-common-ancestor node id for every leaf pair.

    ``LCA[i, j]`` ∈ ``[0, 2N − 2]``: leaf id when ``i == j``, internal
    node id (``≥ N``) otherwise. Computed in O(N²) by walking each
    merge once and recording the cross-product of its two child
    leaf-sets.
    """
    leaves_below: List[List[int]] = [None] * (2 * n_leaves - 1)
    for i in range(n_leaves):
        leaves_below[i] = [i]
    LCA = np.full((n_leaves, n_leaves), -1, dtype=np.int32)
    for i in range(Z.shape[0]):
        nid = n_leaves + i
        l = int(Z[i, 0])
        r = int(Z[i, 1])
        leaves_below[nid] = leaves_below[l] + leaves_below[r]
        left = np.asarray(leaves_below[l], dtype=np.int32)
        right = np.asarray(leaves_below[r], dtype=np.int32)
        LCA[np.ix_(left, right)] = nid
        LCA[np.ix_(right, left)] = nid
    np.fill_diagonal(LCA, np.arange(n_leaves, dtype=np.int32))
    return LCA


def overdraw_leaf_stems(
    ax: plt.Axes,
    leaf_order: Sequence[int],
    cluster_ids: np.ndarray,
    palette: Sequence,
    first_heights: np.ndarray,
    *,
    lw: float = 1.0,
) -> None:
    """Overdraw each leaf's vertical stem in its own cluster colour.

    scipy's ``dendrogram`` draws every U-shape (left-stem,
    horizontal, right-stem) in a single colour from
    ``link_color_func``; when a merge crosses cluster boundaries we
    paint that whole U grey, which leaves the LEAVES grey too.
    Overdrawing the stems on top restores per-leaf cluster colour.
    """
    for j, leaf_id in enumerate(leaf_order):
        x = 5.0 + 10.0 * j  # scipy's leaf x-coordinate convention
        h = first_heights[int(leaf_id)]
        if not np.isfinite(h):
            continue
        cid = int(cluster_ids[int(leaf_id)])
        col = to_hex(palette[(cid - 1) % len(palette)])
        ax.plot([x, x], [0.0, h], color=col, lw=lw, zorder=10, solid_capstyle="butt")


def _node_colors_by_cluster(
    Z: np.ndarray,
    n_leaves: int,
    cluster_ids: np.ndarray,
    palette: List,
    gray: str = "lightgray",
) -> Dict[int, str]:
    """Map every node id (leaves: 0..N-1; internals: N..2N-2) to a
    hex color: cluster's color when all leaves below share one cluster
    id, else gray. Used as ``link_color_func`` for scipy dendrogram and
    for the codebar.
    """
    leaves_below: List[List[int]] = [None] * (2 * n_leaves - 1)
    for i in range(n_leaves):
        leaves_below[i] = [i]
    for i in range(Z.shape[0]):
        nid = n_leaves + i
        left, right = int(Z[i, 0]), int(Z[i, 1])
        leaves_below[nid] = leaves_below[left] + leaves_below[right]

    out: Dict[int, str] = {}
    for nid in range(2 * n_leaves - 1):
        cids = {int(cluster_ids[l]) for l in leaves_below[nid]}
        if len(cids) == 1:
            cid = cids.pop()
            out[nid] = to_hex(palette[(cid - 1) % len(palette)])
        else:
            out[nid] = gray
    return out


# ---------------------------------------------------------------------------
# Channel labels (probe groups) for axis ticks
# ---------------------------------------------------------------------------


def load_channel_labels(patient: str) -> List[str]:
    csv = SEEG_DATAPATH / patient / "channel_labels.csv"
    return pd.read_csv(csv)["label"].astype(str).tolist()


def apply_probe_ticks(
    ax: plt.Axes,
    channel_labels: Sequence[str],
    *,
    is_first_col: bool,
    show_ylabel: bool = True,
) -> None:
    """Ticks at probe-family midpoints; first-col only gets ylabel."""
    groups = probe_groups(list(channel_labels))
    major = [(s + e) / 2 for _, s, e in groups]
    labels = [p for p, _, _ in groups]
    boundaries = [s - 0.5 for _, s, _ in groups[1:]]

    ax.set_xticks(major)
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_xticks(boundaries, minor=True)
    ax.set_xlabel("probe", fontsize=7)

    if is_first_col:
        ax.set_yticks(major)
        ax.set_yticklabels(labels, fontsize=6)
        ax.set_yticks(boundaries, minor=True)
        if show_ylabel:
            ax.set_ylabel("probe", fontsize=7)
    else:
        ax.tick_params(axis="y", which="both", left=False, labelleft=False)

    ax.tick_params(which="major", length=0)
    ax.tick_params(which="minor", length=2.5, color="0.5", width=0.6)


# ---------------------------------------------------------------------------
# Family A — per (patient, band) deep view
# ---------------------------------------------------------------------------


def _draw_codebar(
    ax: plt.Axes,
    cluster_ids: np.ndarray,
    leaf_order: Sequence[int],
    palette: List,
) -> None:
    arr = np.zeros((1, len(leaf_order), 4))
    for k, l in enumerate(leaf_order):
        cid = int(cluster_ids[l])
        arr[0, k] = palette[(cid - 1) % len(palette)]
    ax.imshow(arr, aspect="auto", interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def render_family_a(patient: str, band: str, out_path: Path, *, k_target: int = 20) -> None:
    results = {}
    for phase in PHASES:
        r = load_lrg_result(patient, phase, band, FC_METHOD, IMCOH_LRG_CACHE)
        if r is None:
            print(f"[skip] missing LRG cache: {patient} {band} {phase}")
            return
        results[phase] = r

    Ns = {p: results[p].n_nodes for p in PHASES}
    if len(set(Ns.values())) != 1:
        print(f"[skip] {patient} {band}: phase N mismatch {Ns}")
        return
    N = Ns[PHASES[0]]

    # Compute primitives per phase
    rho = {p: compute_rho(results[p].eigenvalues, results[p].eigenvectors) for p in PHASES}
    D = {p: compute_d_tau(rho[p]) for p in PHASES}

    # Per-phase natural leaf order + fixed-k partition
    cuts: Dict[str, float] = {}
    cluster_ids: Dict[str, np.ndarray] = {}
    node_colors: Dict[str, Dict[int, str]] = {}
    palettes: Dict[str, List] = {}
    for p in PHASES:
        Z = results[p].linkage_matrix
        cuts[p] = fixed_k_cut_height(Z, k_target)
        ids_raw = fcluster(Z, k_target, criterion="maxclust")
        ids = rank_clusters_by_size(ids_raw)
        cluster_ids[p] = ids
        n_clusters = int(np.unique(ids).size)
        palettes[p] = _cluster_palette(n_clusters)
        node_colors[p] = _node_colors_by_cluster(Z, N, ids, palettes[p])

    # Try to load channel labels for probe-group ticks; fall back gracefully.
    try:
        chnames = load_channel_labels(patient)
        # Pat_10 task rows [53,54,55] dropped at load → channel_labels.csv may
        # have a different length than N. If lengths disagree, fall back.
        use_chnames = len(chnames) == N
    except Exception:
        use_chnames = False
        chnames = []

    # ------------------------------------------------------------------
    # Figure layout via nested gridspecs:
    #   outer: 3 rows  (ρ̂ | 𝒟 | dendro+codebar block)
    #   row 2 inner: 2 sub-rows (dendrogram | codebar) with hspace=0
    # Last col widened by `_ROW_CB_FRAC` so the colorbar steals from
    # the outer margin, not from the rightmost imshow.
    # ------------------------------------------------------------------
    width_ratios = [1.0] * 3 + [1.0 + _ROW_CB_FRAC]
    panel_w = 2.4
    fig_w = panel_w * (3 + 1.0 + _ROW_CB_FRAC) + 0.6
    outer_h = [1.0, 1.0, 0.85]
    fig_h = panel_w * sum(outer_h) + 0.6

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs_outer = GridSpec(
        3, 1, figure=fig,
        height_ratios=outer_h,
        hspace=0.18,
        left=0.06,
        right=0.965,
        top=0.965,
        bottom=0.04,
    )

    def _row_axes(outer_cell, *, n_cols=4):
        gs = GridSpecFromSubplotSpec(
            1, n_cols, subplot_spec=outer_cell,
            width_ratios=width_ratios, wspace=0.06,
        )
        return [fig.add_subplot(gs[0, j]) for j in range(n_cols)]

    axes_rho = _row_axes(gs_outer[0])
    axes_D = _row_axes(gs_outer[1])

    # Dendrogram + codebar (zero hspace)
    gs_tree = GridSpecFromSubplotSpec(
        2, 4, subplot_spec=gs_outer[2],
        width_ratios=width_ratios,
        height_ratios=[0.92, 0.08],
        hspace=0.0,
        wspace=0.06,
    )
    axes_dendro = [fig.add_subplot(gs_tree[0, j]) for j in range(4)]
    axes_codebar = [fig.add_subplot(gs_tree[1, j]) for j in range(4)]

    # --- Row 1: ρ̂ (log scale, diagonal masked) --------------------------
    rho_disp = {}
    for p in PHASES:
        d = rho[p].copy()
        np.fill_diagonal(d, np.nan)
        d = np.where(d > 0, d, np.nan)  # LogNorm needs strictly positive
        rho_disp[p] = d
    flat_rho = np.concatenate([rho_disp[p].ravel() for p in PHASES])
    finite = flat_rho[np.isfinite(flat_rho)]
    # Clip vmin/vmax to robust percentiles so the noise floor doesn't
    # spread the dynamic range across many empty decades.
    rho_min = float(np.percentile(finite, 1)) if finite.size else 1e-12
    rho_max = float(np.percentile(finite, 99.5)) if finite.size else 1.0
    rho_norm = LogNorm(vmin=rho_min, vmax=rho_max)
    last_im = None
    for j, phase in enumerate(PHASES):
        ax = axes_rho[j]
        im = ax.imshow(
            rho_disp[phase],
            norm=rho_norm,
            cmap="magma",
            aspect="equal",
            interpolation="nearest",
        )
        last_im = im
        ax.set_title(PHASE_SHORT_LABELS.get(phase, phase), fontsize=10)
        if use_chnames:
            apply_probe_ticks(ax, chnames, is_first_col=(j == 0), show_ylabel=False)
        else:
            ax.set_xticks([])
            ax.set_yticks([])
        if j == 0:
            ax.set_ylabel(r"$\hat{\rho}(\tau_{\max})$  (log)", fontsize=10)
    _, _, cb = imshow_colorbar_caxdivider(last_im, axes_rho[-1], size="4%", pad=0.06)
    _style_log_colorbar(cb)

    # --- Row 2: 𝒟(τ) = (1-δ)/ρ̂  (log scale, diagonal masked) ------------
    D_disp = {}
    for p in PHASES:
        d = D[p].copy()
        d = np.where(d > 0, d, np.nan)  # diag was 0, mask it
        D_disp[p] = d
    flat_D = np.concatenate([D_disp[p].ravel() for p in PHASES])
    finite = flat_D[np.isfinite(flat_D)]
    D_min = float(np.percentile(finite, 0.5)) if finite.size else 1e-12
    D_max = float(np.percentile(finite, 99.5)) if finite.size else 1.0
    D_norm = LogNorm(vmin=D_min, vmax=D_max)
    last_im = None
    for j, phase in enumerate(PHASES):
        ax = axes_D[j]
        im = ax.imshow(
            D_disp[phase],
            norm=D_norm,
            cmap="viridis",
            aspect="equal",
            interpolation="nearest",
        )
        last_im = im
        if use_chnames:
            apply_probe_ticks(ax, chnames, is_first_col=(j == 0), show_ylabel=False)
        else:
            ax.set_xticks([])
            ax.set_yticks([])
        if j == 0:
            ax.set_ylabel(r"$\mathcal{D}(\tau_{\max})$  (log)", fontsize=10)
    _, _, cb = imshow_colorbar_caxdivider(last_im, axes_D[-1], size="4%", pad=0.06)
    _style_log_colorbar(cb)

    # --- Row 3a: Dendrograms (own natural leaf order, colored by cluster) -
    leaf_orders: Dict[str, List[int]] = {}
    for j, phase in enumerate(PHASES):
        ax = axes_dendro[j]
        Z = results[phase].linkage_matrix
        colors = node_colors[phase]

        def link_color_fn(node_id, _colors=colors):
            return _colors.get(int(node_id), "lightgray")

        info = dendrogram(
            Z,
            ax=ax,
            no_labels=True,
            color_threshold=0,
            above_threshold_color="lightgray",
            link_color_func=link_color_fn,
        )
        leaf_orders[phase] = list(info["leaves"])
        # Overdraw leaf stems in each leaf's own cluster colour.
        first_h = _first_merge_heights(Z, N)
        overdraw_leaf_stems(
            ax, info["leaves"], cluster_ids[phase], palettes[phase], first_h,
        )
        ax.axhline(cuts[phase], color="red", linestyle="--", lw=0.7, alpha=0.8)
        heights = Z[:, 2]
        ax.set_ylim(heights[0] * 0.8, heights[-1] * 1.05)
        ax.set_xticks([])
        ax.tick_params(axis="y", labelsize=6)
        if j == 0:
            ax.set_ylabel("dendrogram", fontsize=9)
        else:
            ax.tick_params(axis="y", labelleft=False)

    # --- Row 3b: Codebar — directly under dendrogram (hspace=0) -----------
    for j, phase in enumerate(PHASES):
        ax = axes_codebar[j]
        _draw_codebar(ax, cluster_ids[phase], leaf_orders[phase], palettes[phase])
        if j == 0:
            ax.set_ylabel(
                rf"$k={k_target}$",
                fontsize=8, rotation=0, ha="right", va="center", labelpad=10,
            )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Family C — per (patient, band) dendrogram-structure projection
# ---------------------------------------------------------------------------


def render_family_c(patient: str, band: str, out_path: Path, *, k_target: int = 20) -> None:
    """Color each ``(i, j)`` cell by the cluster colour of its LCA in
    the dendrogram. Within-cluster pairs get that cluster's colour;
    cross-cluster pairs get the LCA's grey. This is a pure
    structural projection — no FC magnitudes — to test whether the
    dendrogram hierarchy maps interpretably onto the matrix.
    """
    results = {}
    for phase in PHASES:
        r = load_lrg_result(patient, phase, band, FC_METHOD, IMCOH_LRG_CACHE)
        if r is None:
            print(f"[skip] missing LRG cache: {patient} {band} {phase}")
            return
        results[phase] = r

    Ns = {p: results[p].n_nodes for p in PHASES}
    if len(set(Ns.values())) != 1:
        print(f"[skip] {patient} {band}: phase N mismatch {Ns}")
        return
    N = Ns[PHASES[0]]

    # Channel labels for probe-group ticks
    try:
        chnames = load_channel_labels(patient)
        use_chnames = len(chnames) == N
    except Exception:
        use_chnames = False
        chnames = []

    # Per-phase: cluster ids, palette, node colors, LCA matrix → RGBA image
    images: Dict[str, np.ndarray] = {}
    for p in PHASES:
        Z = results[p].linkage_matrix
        ids_raw = fcluster(Z, k_target, criterion="maxclust")
        ids = rank_clusters_by_size(ids_raw)
        n_clusters = int(np.unique(ids).size)
        palette = _cluster_palette(n_clusters)
        node_colors = _node_colors_by_cluster(Z, N, ids, palette)
        # Vectorised RGBA lookup
        n_nodes = 2 * N - 1
        rgba = np.zeros((n_nodes, 4), dtype=float)
        for nid in range(n_nodes):
            rgba[nid] = mcolors.to_rgba(node_colors[nid])
        LCA = lca_node_matrix(Z, N)
        images[p] = rgba[LCA]

    panel_w = 2.6
    fig_w = panel_w * 4 + 0.4
    fig_h = panel_w + 0.6
    fig, axes = plt.subplots(
        1, 4,
        figsize=(fig_w, fig_h),
        gridspec_kw={"wspace": 0.06},
    )
    for j, phase in enumerate(PHASES):
        ax = axes[j]
        ax.imshow(images[phase], aspect="equal", interpolation="nearest")
        ax.set_title(PHASE_SHORT_LABELS.get(phase, phase), fontsize=10)
        if use_chnames:
            apply_probe_ticks(ax, chnames, is_first_col=(j == 0), show_ylabel=True)
        else:
            ax.set_xticks([])
            ax.set_yticks([])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Family B — per-patient overview (6 bands × 4 phases dendrograms + codebars)
# ---------------------------------------------------------------------------


def render_family_b(patient: str, out_path: Path, *, k_target: int = 20) -> None:
    bands = list(BRAIN_BANDS_NAMES)
    n_bands = len(bands)

    panel_w = 2.0
    fig_w = panel_w * 4 + 0.8
    fig_h = panel_w * 1.0 * n_bands + 0.5

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs_outer = GridSpec(
        n_bands, 1, figure=fig,
        hspace=0.18, left=0.06, right=0.985, top=0.97, bottom=0.02,
    )

    for bi, band in enumerate(bands):
        results = {}
        skip = False
        for phase in PHASES:
            r = load_lrg_result(patient, phase, band, FC_METHOD, IMCOH_LRG_CACHE)
            if r is None:
                skip = True
                break
            results[phase] = r

        gs_band = GridSpecFromSubplotSpec(
            2, 4, subplot_spec=gs_outer[bi],
            height_ratios=[0.92, 0.08],
            hspace=0.0, wspace=0.06,
        )
        axes_d = [fig.add_subplot(gs_band[0, j]) for j in range(4)]
        axes_c = [fig.add_subplot(gs_band[1, j]) for j in range(4)]

        if skip:
            for ax in axes_d + axes_c:
                ax.axis("off")
            axes_d[0].text(
                0.5, 0.5, f"missing cache: {band}", transform=axes_d[0].transAxes
            )
            continue

        N = results[PHASES[0]].n_nodes

        for j, phase in enumerate(PHASES):
            ax_d = axes_d[j]
            ax_c = axes_c[j]
            Z = results[phase].linkage_matrix
            cut = fixed_k_cut_height(Z, k_target)
            ids_raw = fcluster(Z, k_target, criterion="maxclust")
            ids = rank_clusters_by_size(ids_raw)
            n_clusters = int(np.unique(ids).size)
            palette = _cluster_palette(n_clusters)
            colors = _node_colors_by_cluster(Z, N, ids, palette)

            def link_color_fn(node_id, _colors=colors):
                return _colors.get(int(node_id), "lightgray")

            info = dendrogram(
                Z,
                ax=ax_d,
                no_labels=True,
                color_threshold=0,
                above_threshold_color="lightgray",
                link_color_func=link_color_fn,
            )
            first_h = _first_merge_heights(Z, N)
            overdraw_leaf_stems(ax_d, info["leaves"], ids, palette, first_h)
            ax_d.axhline(cut, color="red", linestyle="--", lw=0.7, alpha=0.8)
            heights = Z[:, 2]
            ax_d.set_ylim(heights[0] * 0.8, heights[-1] * 1.05)
            ax_d.set_xticks([])
            ax_d.tick_params(axis="y", labelsize=6)
            if bi == 0:
                ax_d.set_title(PHASE_SHORT_LABELS.get(phase, phase), fontsize=10)
            if j == 0:
                tex = BRAIN_BAND_TEX_DICT.get(band, band)
                ax_d.set_ylabel(tex, fontsize=10)
            else:
                ax_d.tick_params(axis="y", labelleft=False)

            _draw_codebar(ax_c, ids, list(info["leaves"]), palette)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------


def run_patient(patient: str, *, only_band: Optional[str] = None, k_target: int = 20) -> None:
    out_dir = OUT_ROOT / patient
    out_dir.mkdir(parents=True, exist_ok=True)

    bands = [only_band] if only_band else list(BRAIN_BANDS_NAMES)
    for band in bands:
        out_a = out_dir / f"{band}.pdf"
        out_c = out_dir / f"{band}_dendro_proj.pdf"
        print(f"[A] {patient} {band} (k={k_target}) → {out_a.relative_to(ROOT)}")
        render_family_a(patient, band, out_a, k_target=k_target)
        print(f"[C] {patient} {band} (k={k_target}) → {out_c.relative_to(ROOT)}")
        render_family_c(patient, band, out_c, k_target=k_target)

    if only_band is None:
        out_b = out_dir / "dendrograms_overview.pdf"
        print(f"[B] {patient} overview (k={k_target}) → {out_b.relative_to(ROOT)}")
        render_family_b(patient, out_b, k_target=k_target)


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--patient", default=None)
    p.add_argument("--band", default=None)
    p.add_argument("--cohort", action="store_true")
    p.add_argument("--k", type=int, default=20, help="Target number of flat clusters (maxclust).")
    args = p.parse_args(argv)

    if args.cohort:
        for pat in PATIENTS_4PHASE:
            run_patient(pat, k_target=args.k)
    elif args.patient:
        run_patient(args.patient, only_band=args.band, k_target=args.k)
    else:
        run_patient("Pat_06", only_band="alpha", k_target=args.k)

    return 0


if __name__ == "__main__":
    sys.exit(main())
