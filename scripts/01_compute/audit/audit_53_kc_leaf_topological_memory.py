#!/usr/bin/env python3
"""Audit 53 — KC leaf-topological memory (λ = 0).

Per-leaf, per-pair visualization of KC topology persistence on each
phase's dendrogram. Surfaces *which leaves remember their partners*
across `task_test → rest_post` directly on the tree, and exposes the
band asymmetry between trace bands and ergodic bands at the per-pair
level.

For each pair `(i, j)`:
    m_φ(i, j) = number of edges from the root of `Z_φ` to MRCA(i, j)
    T(i, j)   = ( |m_tt − m_post| ≤ δ ) ∧ ( |m_pre − m_tt| > δ )

Per leaf:
    f_i = | { j : T(i, j) = 1 } | / (n − 1)

Default focal leaf: argmax f_i (smallest-id tiebreak). Default top-M
partners ranked by m_tt(i*, j) descending (closest pairs first → short
paths → readable visual).

Layout: 2 rows × 3 columns. Top row = dendrograms at pre / tt / post
with the focal leaf and its top-M partners highlighted by colored
verticals on the path from leaf to MRCA, plus the MRCA's horizontal U
connector. Bottom row = per-leaf f_i histogram strip aligned with the
dendrogram leaf order; focal bar saturated, partner bars desaturated,
others gray.

Scope report:
    .agents/guides/task-persistence-investigation/2026-05-07_kc-leaf-topological-memory.md

Outputs:
    data/audit/kc_leaf_topological_memory/figures/
        {patient}__{band}__delta-{δ}__topM-{M}.pdf
    data/audit/kc_leaf_topological_memory/per_leaf_scores.csv
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.tree_distance import (
    _build_parent_height_depth,
    kc_vectors,
)
from lrg_eegfc.workflow.lrg import load_lrg_result


PHASES = ("rest_pre", "task_test", "rest_post")
PHASE_LABELS = {
    "rest_pre":  "rest pre",
    "task_test": "task",
    "rest_post": "rest post",
}

ACCENT       = "#d62728"   # default single-focal accent
ACCENT_PART  = "#fa8a87"
GRAY_LINK    = "#bababa"
GRAY_LEAF    = "#c8c8c8"
GRAY_BAR     = "#dadada"

# Multi-focal palette. Each entry = (focal_color, partner_color). Up to
# 6 focals supported; beyond that the palette wraps and visual collisions
# are unavoidable.
FOCAL_PALETTE: list[tuple[str, str]] = [
    ("#d62728", "#fa8a87"),   # red
    ("#1f77b4", "#9bc1de"),   # blue
    ("#2ca02c", "#9dd49d"),   # green
    ("#ff7f0e", "#ffc28a"),   # orange
    ("#9467bd", "#c1a7d6"),   # purple
    ("#17becf", "#9be0ea"),   # cyan
]

OUT_DIR = ROOT / "data/audit/kc_leaf_topological_memory"


# ---------------------------------------------------------------------------
# Pair-index helpers (scipy squareform upper-triangle row-major)
# ---------------------------------------------------------------------------

def pair_index(n: int, i: int, j: int) -> int:
    """Condensed-vector index for unordered pair (i, j) under scipy's
    upper-triangle row-major iteration: pairs are (0,1), (0,2), ...,
    (0,n-1), (1,2), ..., (n-2,n-1)."""
    if i > j:
        i, j = j, i
    return n * i - (i * (i + 1)) // 2 + (j - i - 1)


def subtree_sizes(Z: np.ndarray) -> np.ndarray:
    """Subtree size for every cluster id in `Z`. Length 2n-1; leaves have
    size 1, internal cluster `n + k` has size = sum of children's sizes."""
    n = Z.shape[0] + 1
    sizes = np.ones(2 * n - 1, dtype=np.int64)
    for k in range(n - 1):
        a, b = int(Z[k, 0]), int(Z[k, 1])
        sizes[n + k] = sizes[a] + sizes[b]
    return sizes


def pair_mrca_subtree_sizes(Z: np.ndarray) -> np.ndarray:
    """Condensed-form vector (length n*(n-1)/2) whose entry for pair (i,j)
    in scipy upper-triangle row-major order is the subtree size of
    MRCA(i, j) in `Z`. Mirrors the iteration in `kc_vectors` so the
    indexing stays scipy-condensed-aligned.
    """
    from lrg_eegfc.utils.metrics.tree_distance import (
        _ancestor_chain,
        _build_parent_height_depth,
    )
    n = Z.shape[0] + 1
    parent, _, _ = _build_parent_height_depth(Z)
    sizes = subtree_sizes(Z)
    chains = [_ancestor_chain(i, parent) for i in range(n)]
    chain_len = np.array([len(c) for c in chains], dtype=np.int64)

    npairs = n * (n - 1) // 2
    out = np.zeros(npairs, dtype=np.int64)
    idx = 0
    for i in range(n):
        ci = chains[i]
        for j in range(i + 1, n):
            cj = chains[j]
            limit = min(chain_len[i], chain_len[j])
            mrca = ci[0]
            for k in range(limit):
                if ci[k] == cj[k]:
                    mrca = ci[k]
                else:
                    break
            out[idx] = sizes[mrca]
            idx += 1
    return out


def expand_per_leaf_score(predicate_vec: np.ndarray, n: int) -> np.ndarray:
    """f_i = |{j ≠ i : T(i,j)}| / (n-1) from a condensed predicate vector."""
    counts = np.zeros(n, dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if predicate_vec[idx]:
                counts[i] += 1
                counts[j] += 1
            idx += 1
    return counts.astype(float) / max(n - 1, 1)


# ---------------------------------------------------------------------------
# Dendrogram coordinate model
# ---------------------------------------------------------------------------

def compute_leaf_x_positions(Z: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """Return (x_pos, leaf_order). x_pos has length 2n-1; cluster id
    indexes both leaves (0..n-1) and internal nodes (n..2n-2). Internal
    nodes are positioned at the midpoint of their two children's x.
    Leaves are at scipy's default 10-unit spacing starting at 5.
    """
    n = Z.shape[0] + 1
    R = dendrogram(Z, no_plot=True)
    leaf_order = list(R["leaves"])

    x_pos = np.zeros(2 * n - 1, dtype=float)
    for pos, leaf in enumerate(leaf_order):
        x_pos[leaf] = 5.0 + 10.0 * pos
    for k in range(n - 1):
        a, b = int(Z[k, 0]), int(Z[k, 1])
        x_pos[n + k] = 0.5 * (x_pos[a] + x_pos[b])
    return x_pos, leaf_order


def find_mrca(parent: np.ndarray, leaf_i: int, leaf_j: int) -> int:
    """Lowest common ancestor of two leaves in the linkage tree."""
    seen_i = set()
    c = leaf_i
    while c != -1:
        seen_i.add(c)
        c = int(parent[c])
    c = leaf_j
    while c != -1:
        if c in seen_i:
            return c
        c = int(parent[c])
    return -1


def path_segments_to_mrca(
    leaf: int,
    mrca: int,
    parent: np.ndarray,
    height: np.ndarray,
    x_pos: np.ndarray,
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Vertical segments from `leaf` up to (but not including) `mrca`.

    Each segment is at x = x_pos[c] for cluster c on the chain, going
    from y = height[c] to y = height[parent[c]]. The leaf-to-first-merge
    segment has y ∈ [0, h_first_parent]; subsequent internal-to-internal
    segments have y ∈ [h_c, h_parent]. Stops at the child of mrca.
    """
    segs = []
    c = leaf
    while c != mrca:
        p = int(parent[c])
        if p == -1:
            break
        segs.append((
            (x_pos[c], height[c]),
            (x_pos[c], height[p]),
        ))
        if p == mrca:
            break
        c = p
    return segs


def pick_top_focals(
    f: np.ndarray,
    T: np.ndarray,
    m_tt: np.ndarray,
    n: int,
    K: int,
    top_m: int,
    *,
    max_overlap_frac: float = 0.30,
) -> list[dict]:
    """Greedy disjoint focal-leaf selection.

    For each accepted focal `i_k`, the union {focal} ∪ top-M partners is
    added to the `used` set. A new candidate `i` is rejected if its own
    {i} ∪ top-M partner-pattern overlaps `used` by more than
    `max_overlap_frac`. Stops when K focals accepted or no more
    candidates have f > 0.
    """
    accepted: list[dict] = []
    used: set[int] = set()
    order = np.argsort(-f, kind="stable")
    for cand in order.tolist():
        if len(accepted) >= K:
            break
        if f[cand] <= 0.0:
            break
        if cand in used:
            continue
        # Top-M partners of cand
        partners_pred = []
        for j in range(n):
            if j == cand:
                continue
            k = pair_index(n, cand, j)
            if T[k]:
                partners_pred.append(j)
        partners_pred.sort(key=lambda j: -int(m_tt[pair_index(n, cand, j)]))
        top_M = partners_pred[: top_m]
        pattern = set(top_M) | {cand}
        if used:
            ov = len(pattern & used) / max(len(pattern), 1)
            if ov > max_overlap_frac:
                continue
        accepted.append({
            "focal": cand,
            "f_focal": float(f[cand]),
            "top_M": list(top_M),
            "n_partners_pred": len(partners_pred),
        })
        used |= pattern
    return accepted


def mrca_horizontal(
    Z: np.ndarray,
    mrca: int,
    n: int,
    x_pos: np.ndarray,
    height: np.ndarray,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Horizontal connector at the MRCA's U: from x of one child to x of
    the other child, both at h_{MRCA}."""
    k = mrca - n
    a, b = int(Z[k, 0]), int(Z[k, 1])
    return ((x_pos[a], height[mrca]), (x_pos[b], height[mrca]))


# ---------------------------------------------------------------------------
# Per-focal page (single-focal Pat_07-style layout)
# ---------------------------------------------------------------------------

def _render_focal_page(
    fc: dict,
    focals: list,
    rank: int,
    f: np.ndarray,
    n: int,
    coord_by_phase: dict,
    m_tt: np.ndarray,
    *,
    patient: str,
    band: str,
    band_tex: str,
    top_m: int,
    delta: int,
    max_mrca_frac: float,
    s_max: int,
) -> "plt.Figure":
    """Render one focal's pattern as the single-focal Pat_07-style page:
    top row = 3 dendrograms with this focal's path bundle highlighted,
    bottom row = 3 histograms of f_i with this focal+partners colored.

    The other focals on the cell are summarised in a small grey legend
    strip at the bottom of the page so the reader can navigate to the
    other pages knowing how this pattern fits into the full top-K.
    """
    focal = fc["focal"]
    top_M = fc["top_M"]
    color = fc["color"]
    color_partner = fc["color_partner"]

    fig = plt.figure(figsize=(11.5, 6.4))
    gs = fig.add_gridspec(
        2, 3,
        height_ratios=[1.8, 0.9],
        left=0.05, right=0.99, top=0.90, bottom=0.10,
        hspace=0.28, wspace=0.06,
    )

    den_axes: list = []
    his_axes: list = []
    for c, ph in enumerate(PHASES):
        sharey_d = den_axes[0] if den_axes else None
        ax_d = fig.add_subplot(gs[0, c], sharey=sharey_d)
        den_axes.append(ax_d)
        sharey_h = his_axes[0] if his_axes else None
        ax_h = fig.add_subplot(gs[1, c], sharex=ax_d, sharey=sharey_h)
        his_axes.append(ax_h)

    # ---- Dendrograms ----
    for c, ph in enumerate(PHASES):
        ax = den_axes[c]
        coord = coord_by_phase[ph]
        Z = coord["Z"]
        parent = coord["parent"]
        height = coord["height"]
        x_pos  = coord["x_pos"]

        dendrogram(
            Z, ax=ax, no_labels=True,
            color_threshold=0,
            above_threshold_color=GRAY_LINK,
            link_color_func=lambda _id: GRAY_LINK,
        )
        merge_heights = sorted(set(Z[:, 2]))
        tmin = 0.0
        tmax = merge_heights[-1] * 1.05
        ax.set_ylim(tmin, tmax)
        ax.set_xticks([])
        ax.tick_params(axis="y", labelsize=7)
        for sn in ("top", "right", "bottom"):
            ax.spines[sn].set_visible(False)
        ax.spines["left"].set_linewidth(0.5)
        ax.set_title(PHASE_LABELS[ph], fontsize=11, pad=4)
        if c == 0:
            ax.set_ylabel("merge height", fontsize=9)
        else:
            ax.tick_params(axis="y", labelleft=False)
            ax.spines["left"].set_visible(False)

        # Path tracing for this single focal
        if top_M:
            rgba_path = to_rgba(color, 0.55)
            rgba_mrca = to_rgba(color, 0.80)
            segs: list = []
            colors: list = []
            for j in top_M:
                mrca = find_mrca(parent, focal, j)
                if mrca < n:
                    continue
                segs_i = path_segments_to_mrca(focal, mrca, parent, height, x_pos)
                segs_j = path_segments_to_mrca(j,     mrca, parent, height, x_pos)
                mrca_h = mrca_horizontal(Z, mrca, n, x_pos, height)
                for s in segs_i:
                    segs.append(s); colors.append(rgba_path)
                for s in segs_j:
                    segs.append(s); colors.append(rgba_path)
                segs.append(mrca_h); colors.append(rgba_mrca)
            if segs:
                ax.add_collection(
                    LineCollection(segs, colors=colors, linewidths=1.6, zorder=4)
                )

        # Leaf bars: focal saturated, partners desaturated, rest gray.
        for lf in range(n):
            if lf == focal:
                ax.vlines(x_pos[lf], tmin, height[int(parent[lf])],
                          color=color, lw=1.9, zorder=5)
            elif lf in top_M:
                ax.vlines(x_pos[lf], tmin, height[int(parent[lf])],
                          color=color_partner, lw=1.4, zorder=4)

        # Caret marker above the focal leaf (this focal only).
        ax.scatter([x_pos[focal]], [tmax * 0.985], marker="v",
                   color=color, s=46, zorder=6,
                   edgecolors="white", linewidths=0.6)

    # ---- Histograms (per phase, leaf order matching that phase's dendrogram) ----
    bar_h = np.array([f[lf] for lf in range(n)])
    ymax = max(0.05, float(bar_h.max()) * 1.15)
    for c, ph in enumerate(PHASES):
        ax = his_axes[c]
        x_pos = coord_by_phase[ph]["x_pos"]
        bar_x = np.array([x_pos[lf] for lf in range(n)])

        bar_colors = []
        for lf in range(n):
            if lf == focal:
                bar_colors.append(color)
            elif lf in top_M:
                bar_colors.append(color_partner)
            else:
                bar_colors.append(GRAY_BAR)

        ax.bar(bar_x, bar_h, width=8.6, color=bar_colors,
               edgecolor="none", zorder=2)
        ax.set_ylim(0.0, ymax)
        ax.set_xticks([])
        ax.set_xlim(0.0, 10.0 * n)
        ax.tick_params(axis="y", labelsize=7)
        for sn in ("top", "right", "bottom"):
            ax.spines[sn].set_visible(False)
        ax.spines["left"].set_linewidth(0.5)
        for y in (0.05, 0.10, 0.15, 0.25, 0.5, 0.75):
            if y < ymax:
                ax.axhline(y, color="#dddddd", lw=0.4, zorder=1)
        if c == 0:
            ax.set_ylabel(r"$f_i$", fontsize=9)
        else:
            ax.tick_params(axis="y", labelleft=False)
            ax.spines["left"].set_visible(False)

    # ---- Title + page footer ----
    title = (
        rf"{patient}  {band_tex}   pattern {rank+1}/{len(focals)}   "
        rf"focal $i^\star={focal}$,  $f^\star={fc['f_focal']:.2f}$,  "
        rf"$|P|={fc['n_partners_pred']}$,  top-{top_m} drawn"
    )
    fig.text(0.5, 0.955, title, ha="center", va="center", fontsize=11)

    other_focals = [
        rf"#{i+1}: $i={ofc['focal']}\,(f={ofc['f_focal']:.2f})$"
        for i, ofc in enumerate(focals) if i != rank
    ]
    if other_focals:
        legend = "other patterns this cell — " + "  ·  ".join(other_focals)
    else:
        legend = "this is the only pattern surviving the predicate at this cell"
    fig.text(0.5, 0.030, legend, ha="center", va="center",
             fontsize=8, color="#666")
    fig.text(
        0.5, 0.005,
        rf"$\delta={delta}$,  $s_{{\max}}={max_mrca_frac:.2f}\,n={s_max}$,  "
        rf"$n={n}$ leaves",
        ha="center", va="center", fontsize=7, color="#888",
    )

    return fig


# ---------------------------------------------------------------------------
# Per-cell render
# ---------------------------------------------------------------------------

def render_cell(
    patient: str,
    band: str,
    *,
    delta: int = 1,
    top_m: int = 8,
    min_m: int = 2,
    max_mrca_frac: float = 0.50,
    top_n_focals: int = 4,
    focal_overlap_frac: float = 0.30,
    min_partners: int = 5,
    focal_override: Optional[int] = None,
) -> dict:
    Z_by_phase = {}
    for ph in PHASES:
        res = load_lrg_result(patient, ph, band, fc_method="imcoh_abs")
        if res is None or res.linkage_matrix is None:
            raise FileNotFoundError(f"missing LRG result {patient} {ph} {band}")
        Z_by_phase[ph] = np.asarray(res.linkage_matrix)

    n = int(Z_by_phase["rest_pre"].shape[0]) + 1
    for ph in PHASES:
        n_ph = int(Z_by_phase[ph].shape[0]) + 1
        if n_ph != n:
            raise ValueError(f"{patient} {band}: leaf count differs across phases ({n} vs {n_ph})")

    # KC m-vectors (topology only) + per-pair MRCA subtree sizes
    m_pre  = kc_vectors(Z_by_phase["rest_pre"])[0]
    m_tt   = kc_vectors(Z_by_phase["task_test"])[0]
    m_post = kc_vectors(Z_by_phase["rest_post"])[0]
    size_tt   = pair_mrca_subtree_sizes(Z_by_phase["task_test"])
    size_post = pair_mrca_subtree_sizes(Z_by_phase["rest_post"])

    # Predicate. The size cap (max_mrca_frac · n) requires the MRCA
    # subtree in BOTH tt and post to be fine-grained — without it the
    # predicate is dominated by near-root anchor pairs (one leaf in a
    # tiny subtree absorbed at one merge near the root produces O(n)
    # m-matched but visually-empty "paths span the whole tree" hits;
    # see Pat_07 θ leaf-38 diagnostic, scope-report Caveats §).
    s_max = int(np.ceil(max_mrca_frac * n))
    matched_tt_post = np.abs(m_tt - m_post) <= delta
    differs_pre_tt  = np.abs(m_pre - m_tt) > delta
    nontrivial_m    = (m_tt >= min_m) & (m_post >= min_m)
    fine_grained    = (size_tt <= s_max) & (size_post <= s_max)
    T = matched_tt_post & differs_pre_tt & nontrivial_m & fine_grained

    # Per-leaf score
    f = expand_per_leaf_score(T, n)

    # Focal selection: explicit override → 1 focal; otherwise greedy
    # disjoint top-K. Each focal record carries (focal, f_focal,
    # top_M, n_partners_pred).
    if focal_override is not None:
        cand = int(focal_override)
        partners_pred = []
        for j in range(n):
            if j == cand:
                continue
            k = pair_index(n, cand, j)
            if T[k]:
                partners_pred.append(j)
        partners_pred.sort(key=lambda j: -int(m_tt[pair_index(n, cand, j)]))
        focals = [{
            "focal": cand,
            "f_focal": float(f[cand]),
            "top_M": partners_pred[: top_m],
            "n_partners_pred": len(partners_pred),
        }]
    else:
        focals = pick_top_focals(
            f, T, m_tt, n,
            K=top_n_focals, top_m=top_m,
            max_overlap_frac=focal_overlap_frac,
        )
        focals = [fc for fc in focals if fc["n_partners_pred"] >= min_partners]
    if not focals:
        # Predicate empty → record sentinel focal at argmax-of-f for the
        # caret marker, but with empty top_M so no paths are drawn.
        empty_focal = int(np.argmax(f))
        focals = [{
            "focal": empty_focal,
            "f_focal": float(f[empty_focal]),
            "top_M": [],
            "n_partners_pred": 0,
        }]

    # Per-focal palette
    for k, fc in enumerate(focals):
        sat, des = FOCAL_PALETTE[k % len(FOCAL_PALETTE)]
        fc["color"] = sat
        fc["color_partner"] = des

    # Resolve a single (focal_idx, role) per leaf for bar coloring.
    # Priority: lower focal_idx wins (highest f focal first); within
    # one focal, focal-leaf > partner > none.
    leaf_assign: list[Optional[tuple[int, str]]] = [None] * n
    for k, fc in enumerate(focals):
        if leaf_assign[fc["focal"]] is None:
            leaf_assign[fc["focal"]] = (k, "focal")
    for k, fc in enumerate(focals):
        for j in fc["top_M"]:
            if leaf_assign[j] is None:
                leaf_assign[j] = (k, "partner")

    # Per-phase coordinate model + path tracing
    coord_by_phase = {}
    for ph in PHASES:
        Z = Z_by_phase[ph]
        parent, height, _ = _build_parent_height_depth(Z)
        x_pos, leaf_order = compute_leaf_x_positions(Z)
        coord_by_phase[ph] = {
            "Z": Z,
            "parent": parent,
            "height": height,
            "x_pos": x_pos,
            "leaf_order": leaf_order,
        }

    # ---- Figure: multi-page PDF, one page per focal, each page is the
    # single-focal Pat_07-style layout (2 rows × 3 cols: dendrograms +
    # f_i histograms aligned to each phase's leaf order).
    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)

    smax_tag = f"{int(round(max_mrca_frac * 100)):02d}"
    n_focals_tag = f"K{len(focals)}"
    out_path = OUT_DIR / "figures" / (
        f"{patient}__{band}__delta-{delta}__topM-{top_m}__smax-{smax_tag}__{n_focals_tag}.pdf"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(out_path) as pdf:
        for k, fc in enumerate(focals):
            fig = _render_focal_page(
                fc, focals, k, f, n,
                coord_by_phase, m_tt,
                patient=patient, band=band, band_tex=band_tex,
                top_m=top_m, delta=delta, max_mrca_frac=max_mrca_frac,
                s_max=s_max,
            )
            pdf.savefig(fig)
            plt.close(fig)

    return {
        "patient": patient,
        "band": band,
        "n": n,
        "delta": delta,
        "top_m": top_m,
        "n_focals": len(focals),
        "focals": focals,
        "f_mean": float(f.mean()),
        "f_max": float(f.max()),
        "out_path": out_path,
        "f_per_leaf": f,
    }


# ---------------------------------------------------------------------------
# CLI / batch
# ---------------------------------------------------------------------------

DEFAULT_CELLS = [
    ("Pat_07", "beta"),     # multiscale trace flagship
    ("Pat_07", "theta"),    # ergodic counterexample
]

COHORT = (
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
)
TRACE_BANDS = ("beta",)
ERGODIC_BANDS = ("theta", "high_gamma")


def main(argv: list[str] | None = None) -> int:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--patient", default=None,
                   help="restrict to a single patient")
    p.add_argument("--band", default=None,
                   help="restrict to a single band")
    p.add_argument("--cohort", action="store_true",
                   help="sweep all 10 patients × {beta, theta, high_gamma}")
    p.add_argument("--delta", type=int, default=1,
                   help="integer tolerance on m-match (default 1: allow "
                        "off-by-one MRCA-depth match in tt vs post; pre "
                        "must differ by > δ)")
    p.add_argument("--top-m", type=int, default=8,
                   help="top-M partners per focal (default 8)")
    p.add_argument("--top-n-focals", type=int, default=4,
                   help="top-K disjoint focals per cell (default 4)")
    p.add_argument("--focal-overlap", type=float, default=0.30,
                   help="max fraction overlap between successive focals' "
                        "patterns (default 0.30)")
    p.add_argument("--min-m", type=int, default=2,
                   help="drop pairs whose MRCA depth from root in tt is "
                        "below this (default 2 = exclude near-root pairs)")
    p.add_argument("--max-mrca-frac", type=float, default=0.50,
                   help="MRCA-subtree fraction cap in tt and post "
                        "(default 0.50: keep mid-scale and finer matches)")
    p.add_argument("--min-partners", type=int, default=5,
                   help="drop pages with fewer partners than this (default 5: "
                        "ensures each rendered pattern reads as a bundle)")
    p.add_argument("--focal", type=int, default=None,
                   help="override focal selection with a single leaf id")
    args = p.parse_args(argv)

    if args.cohort:
        cells = [(pat, band)
                 for pat in COHORT
                 for band in TRACE_BANDS + ERGODIC_BANDS]
    elif args.patient is not None and args.band is not None:
        cells = [(args.patient, args.band)]
    else:
        cells = DEFAULT_CELLS

    rows = []
    summary = []
    for pat, band in cells:
        try:
            info = render_cell(
                pat, band,
                delta=args.delta, top_m=args.top_m, min_m=args.min_m,
                max_mrca_frac=args.max_mrca_frac,
                top_n_focals=args.top_n_focals,
                focal_overlap_frac=args.focal_overlap,
                min_partners=args.min_partners,
                focal_override=args.focal,
            )
        except Exception as exc:
            print(f"  {pat:7s} {band:11s}  ERROR: {exc}")
            continue
        focal_str = " ".join(
            f"{fc['focal']:3d}({fc['f_focal']:.2f},|P|={fc['n_partners_pred']})"
            for fc in info["focals"]
        )
        print(
            f"  {pat:7s} {band:11s}  "
            f"n={info['n']:3d}  K={info['n_focals']}  focals=[{focal_str}]  "
            f"f_mean={info['f_mean']:.3f}  -> {info['out_path'].name}"
        )

        # CSV: one row per focal
        for k, fc in enumerate(info["focals"]):
            summary.append({
                "patient": pat, "band": band,
                "focal_rank": k,
                "focal": fc["focal"],
                "f_focal": fc["f_focal"],
                "n_partners_pred": fc["n_partners_pred"],
                "top_M_leaves": ",".join(map(str, fc["top_M"])),
                "delta": args.delta,
                "top_m": args.top_m,
                "max_mrca_frac": args.max_mrca_frac,
            })
        # CSV: one row per leaf x cell
        for lf, val in enumerate(info["f_per_leaf"]):
            rows.append({
                "patient": pat, "band": band, "leaf": lf,
                "f_i": float(val),
                "delta": args.delta,
            })

    if rows:
        out_csv = OUT_DIR / "per_leaf_scores.csv"
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(out_csv, index=False)
        print(f"\nwrote {out_csv.relative_to(ROOT)}")
    if summary:
        out_summary = OUT_DIR / "focal_summary.csv"
        pd.DataFrame(summary).to_csv(out_summary, index=False)
        print(f"wrote {out_summary.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
