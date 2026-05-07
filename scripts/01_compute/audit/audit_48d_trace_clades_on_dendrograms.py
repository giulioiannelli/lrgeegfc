#!/usr/bin/env python3
"""Audit 48d -- visualise trace / reset / persist / rearrange clades ON the
4-phase dendrograms.

Four modes (matching the canonical trace/anchor/reset/emergent taxonomy
plus the user-named 'persist' = canonical 'anchor'):

  --mode trace        anchor=task_test
                      keep C with bestJ(C, rs_post)>=J_high AND bestJ(C, rs_pre)<=J_low
                      visual: scatter -> cluster -> persist
                      (TRACE: emergent in task AND persists in post)

  --mode reset        anchor=rs_pre
                      keep C with bestJ(C, rs_post)>=J_high AND bestJ(C, task_test)<=J_low
                      visual: cluster -> scatter -> cluster
                      (RESET: in pre, scattered by task, returns in post)

  --mode persist      anchor=rs_pre  (= canonical ANCHOR)
                      keep C with bestJ(C, task_test)>=J_high AND bestJ(C, rs_post)>=J_high
                      visual: cluster -> cluster -> cluster
                      (PERSIST: module never moves through phases)

  --mode rearrange    anchor=rs_post
                      keep C with bestJ(C, rs_pre)<=J_low AND bestJ(C, task_test)<=J_low
                      visual: scatter -> scatter -> cluster
                      (EMERGENT-IN-POST: emergent only in rs_post, not in task)

CRITICAL: bestJ(C, T) = max over ALL internal-node leafsets of T (NOT
just same-cut-level clusters). Using same-cut-level Jaccard would
mis-classify persist as reset when C's leaves stay tightly grouped
inside a LARGER task_test cluster.

The four modes are MUTUALLY EXCLUSIVE at the predicate level: trace
requires absence in pre; reset requires absence in task; persist
requires presence in both pre AND task AND post; rearrange requires
absence in both pre AND task. No clade can satisfy two predicates.

Each mode produces one PDF per (patient, band) cell with at least one
matching clade.

Output dirs:
  data/reports/section_5_lrg_trace/14_rf_clade_persistence/figures/
    trace_clades_on_dendrograms/
    reset_clades_on_dendrograms/
    persist_clades_on_dendrograms/
    rearrange_clades_on_dendrograms/

By default runs only Pat_06 β as a pilot; --all to run the cohort.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.tree import tree_internal_nodes
from lrg_eegfc.workflow.lrg import load_lrg_result

PATIENTS = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
PHASE_TITLES = {
    "rest_pre": r"rs$_{\rm pre}$",
    "task_learn": "task learn",
    "task_test": "task test",
    "rest_post": r"rs$_{\rm post}$",
}
TRACE_COLORS = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e"]
DEFAULT_COLOR = "#bbbbbb"

OUT_BASE = ROOT / "data" / "reports" / "section_5_lrg_trace" / "14_rf_clade_persistence" / "figures"
MODE_DIRS = {
    "trace": OUT_BASE / "trace_clades_on_dendrograms",
    "reset": OUT_BASE / "reset_clades_on_dendrograms",
    "persist": OUT_BASE / "persist_clades_on_dendrograms",
    "rearrange": OUT_BASE / "rearrange_clades_on_dendrograms",
}
MODE_HEADLINE = {
    "trace": r"TRACE: scatter $\to$ cluster $\to$ persist",
    "reset": r"RESET: cluster $\to$ scatter $\to$ cluster",
    "persist": r"PERSIST (= anchor): cluster across all phases",
    "rearrange": r"REARRANGE: scatter $\to$ scatter $\to$ cluster",
}


# ---------------------------------------------------------------------------
# Trace-clade enumeration
# ---------------------------------------------------------------------------

def jaccard(A: frozenset, B: frozenset) -> float:
    if not A or not B:
        return 0.0
    inter = len(A & B)
    uni = len(A | B)
    return inter / uni if uni else 0.0


def _enumerate_clades_at(Z, k):
    L = fcluster(Z, t=k, criterion="maxclust")
    return [frozenset(np.where(L == lab)[0].tolist()) for lab in np.unique(L)]


def _select_disjoint(cand, top_n):
    cand.sort(reverse=True)
    selected = []
    used = set()
    for tup in cand:
        C = tup[-1]
        if len(C & used) <= 0.2 * len(C):
            selected.append(tup)
            used |= C
            if len(selected) >= top_n:
                break
    return selected


def all_subtree_leafsets(Z) -> list[frozenset]:
    """All internal-node leafsets of Z (one per internal node, N-1 total)."""
    return [n["leaves"] for n in tree_internal_nodes(Z)]


def best_jaccard_against(C: frozenset, leafsets: list[frozenset]) -> float:
    """Max Jaccard between C and any leafset in the list. The 'is C
    represented anywhere in this tree at any scale' question."""
    return max((jaccard(C, S) for S in leafsets), default=0.0)


def find_clades(mode, Z_pre, Z_tt, Z_post,
                k_grid=(10, 14, 18, 22, 26, 32, 40),
                J_high=0.65, J_low=0.30,
                size_min=6, size_max=35, top_n=3):
    """Unified mode dispatcher with all-internal-node Jaccard checks.

    Predicates (bestJ = max Jaccard over ALL internal-node leafsets of T):

      trace      anchor=task_test    bestJ(C, rs_post) >= J_high
                                     bestJ(C, rs_pre)  <= J_low

      reset      anchor=rs_pre       bestJ(C, rs_post) >= J_high
                                     bestJ(C, task_test) <= J_low

      persist    anchor=rs_pre       bestJ(C, task_test) >= J_high
                                     bestJ(C, rs_post)   >= J_high

      rearrange  anchor=rs_post      bestJ(C, rs_pre)    <= J_low
                                     bestJ(C, task_test) <= J_low

    Score by mode:
      trace:     J_post - J_pre        (max separation)
      reset:     J_post - J_tt
      persist:   min(J_tt, J_post)     (weakest leg = bottleneck)
      rearrange: 1 - max(J_pre, J_tt)  (max isolation)
    """
    pre_subs = all_subtree_leafsets(Z_pre)
    tt_subs = all_subtree_leafsets(Z_tt)
    post_subs = all_subtree_leafsets(Z_post)

    if mode in ("trace",):
        anchor_Z = Z_tt
    elif mode in ("reset", "persist"):
        anchor_Z = Z_pre
    elif mode == "rearrange":
        anchor_Z = Z_post
    else:
        raise ValueError(f"unknown mode: {mode}")

    cand = []
    seen_cset = set()
    for k in k_grid:
        for C in _enumerate_clades_at(anchor_Z, k):
            if not (size_min <= len(C) <= size_max):
                continue
            if C in seen_cset:
                continue
            seen_cset.add(C)
            bj_pre = best_jaccard_against(C, pre_subs)
            bj_tt = best_jaccard_against(C, tt_subs)
            bj_post = best_jaccard_against(C, post_subs)

            if mode == "trace":
                if bj_post >= J_high and bj_pre <= J_low:
                    cand.append((bj_post - bj_pre, bj_post, bj_pre, k, C))
            elif mode == "reset":
                if bj_post >= J_high and bj_tt <= J_low:
                    cand.append((bj_post - bj_tt, bj_post, bj_tt, k, C))
            elif mode == "persist":
                if bj_tt >= J_high and bj_post >= J_high:
                    cand.append((min(bj_tt, bj_post), bj_tt, bj_post, k, C))
            elif mode == "rearrange":
                if bj_pre <= J_low and bj_tt <= J_low:
                    isolation = 1.0 - max(bj_pre, bj_tt)
                    cand.append((isolation, bj_pre, bj_tt, k, C))
    return _select_disjoint(cand, top_n)


# ---------------------------------------------------------------------------
# Dendrogram coloring -- leaves BELOW each trace clade colored in same color
# ---------------------------------------------------------------------------

def descendants_cache(Z: np.ndarray) -> dict[int, frozenset[int]]:
    """Map node_id -> frozenset of leaf indices for every internal node id."""
    n = Z.shape[0] + 1
    cache: dict[int, frozenset[int]] = {}
    for i in range(Z.shape[0]):
        l, r = int(Z[i, 0]), int(Z[i, 1])
        ll = frozenset([l]) if l < n else cache[l]
        rr = frozenset([r]) if r < n else cache[r]
        cache[n + i] = ll | rr
    return cache


def make_link_color_func(Z: np.ndarray, trace_leafsets: list[frozenset[int]],
                          colors: list[str], default: str = DEFAULT_COLOR):
    n = Z.shape[0] + 1
    desc = descendants_cache(Z)

    def color(node_id: int) -> str:
        if node_id < n:
            return default
        leaves = desc[node_id]
        for cset, c in zip(trace_leafsets, colors):
            if leaves.issubset(cset):
                return c
        return default

    return color


def leaf_color_array(n: int, trace_leafsets: list[frozenset[int]],
                     colors: list[str]) -> list[str]:
    out = [DEFAULT_COLOR] * n
    for cset, c in zip(trace_leafsets, colors):
        for leaf in cset:
            out[leaf] = c
    return out


def leaf_first_merge_heights(Z: np.ndarray) -> np.ndarray:
    """For each leaf, return the height of the merge that first incorporates it.

    Each leaf appears as a child in exactly one row of Z (after that, only
    the internal node it merged into is referenced). The merge order in Z
    is ascending by height, so the first appearance is the smallest-height
    merge involving the leaf.
    """
    n = Z.shape[0] + 1
    heights = np.full(n, np.nan, dtype=float)
    for i in range(Z.shape[0]):
        l, r = int(Z[i, 0]), int(Z[i, 1])
        h = float(Z[i, 2])
        if l < n and np.isnan(heights[l]):
            heights[l] = h
        if r < n and np.isnan(heights[r]):
            heights[r] = h
    return heights


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def plot_one_cell(patient: str, band: str, mode: str = "trace"):
    Z_dict = {}
    for phase in PHASES:
        res = load_lrg_result(patient, phase, band, "imcoh_abs")
        if res is None:
            print(f"  [{patient} {band} {mode}] {phase} cache miss, skip")
            return None
        Z_dict[phase] = res.linkage_matrix
    n = Z_dict["task_test"].shape[0] + 1

    selected = find_clades(mode, Z_dict["rest_pre"], Z_dict["task_test"],
                            Z_dict["rest_post"])
    if not selected:
        print(f"  [{patient} {band} {mode}] no clades at strict gates -- relaxing")
        selected = find_clades(mode, Z_dict["rest_pre"], Z_dict["task_test"],
                                Z_dict["rest_post"], J_high=0.55, J_low=0.40)
    if not selected:
        print(f"  [{patient} {band} {mode}] no clades found")
        return None
    trace_leafsets = [c for _, _, _, _, c in selected]
    colors = TRACE_COLORS[:len(trace_leafsets)]

    # Plot 4 dendrograms in 1 row
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.4),
                              gridspec_kw={"wspace": 0.05})

    leaf_colors = leaf_color_array(n, trace_leafsets, colors)

    for j, phase in enumerate(PHASES):
        ax = axes[j]
        Z = Z_dict[phase]
        lc_func = make_link_color_func(Z, trace_leafsets, colors)
        d = dendrogram(
            Z, ax=ax, link_color_func=lc_func,
            no_labels=True, above_threshold_color=DEFAULT_COLOR,
            count_sort="ascending",
        )
        leaves_in_order = d["leaves"]
        # x position of each leaf (scipy default unit = 10, centered at +5)
        leaf_x = {int(leaf): (idx * 10 + 5) for idx, leaf in enumerate(leaves_in_order)}
        merge_h = leaf_first_merge_heights(Z)

        # Overlay colored vertical lines for EVERY trace leaf, regardless of
        # whether its first-merge sibling is in the same clade. This makes
        # scattered leaves visible on the dendrogram itself (not only on the
        # barcode strip).
        for cset, c in zip(trace_leafsets, colors):
            for leaf in cset:
                x = leaf_x[int(leaf)]
                h = merge_h[int(leaf)]
                if not np.isnan(h):
                    ax.plot([x, x], [0, h], color=c, lw=1.4, zorder=4,
                            solid_capstyle="butt")

        # Barcode strip at the bottom indicating leaf cluster membership
        ymin, ymax = ax.get_ylim()
        strip_y = ymin - 0.04 * (ymax - ymin)
        strip_h = 0.025 * (ymax - ymin)
        for x_idx, leaf in enumerate(leaves_in_order):
            c = leaf_colors[int(leaf)]
            if c != DEFAULT_COLOR:
                ax.add_patch(plt.Rectangle(
                    ((x_idx + 0.5) * 10 - 4, strip_y - strip_h),
                    8, strip_h, facecolor=c, edgecolor="none",
                    clip_on=False, zorder=5,
                ))
        ax.set_ylim(strip_y - strip_h * 1.2, ymax)
        ax.set_xticks([])
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(PHASE_TITLES[phase], fontsize=11)
        if j == 0:
            ax.set_ylabel("merge height")

    if mode == "trace":
        label_fmt = "clade #{i} (|L|={n}, J$_{{post}}$={a:.2f}, J$_{{pre}}$={b:.2f}, k={k})"
    elif mode == "reset":
        label_fmt = "clade #{i} (|L|={n}, J$_{{post}}$={a:.2f}, J$_{{tt}}$={b:.2f}, k={k})"
    elif mode == "persist":
        label_fmt = "clade #{i} (|L|={n}, J$_{{tt}}$={a:.2f}, J$_{{post}}$={b:.2f}, k={k})"
    else:  # rearrange
        label_fmt = "clade #{i} (|L|={n}, J$_{{pre}}$={a:.2f}, J$_{{tt}}$={b:.2f}, k={k})"
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=c,
                       label=label_fmt.format(i=i+1, n=len(cs), a=a_val,
                                                b=b_val, k=kk))
        for i, (c, (_, a_val, b_val, kk, cs)) in enumerate(zip(TRACE_COLORS, selected))
    ]
    fig.text(0.5, 0.97, MODE_HEADLINE[mode], ha="center", va="top",
             fontsize=10, fontweight="bold")
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.02),
               ncol=len(handles), frameon=False, fontsize=9)

    out_dir = MODE_DIRS[mode]
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{patient}_{band}.pdf"
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  [{patient} {band} {mode}] {len(selected)} clades -> {out_path.name}")
    return out_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true",
                    help="Run all 10 patients × 6 bands (default: pilot Pat_06 β only)")
    ap.add_argument("--patient", default="Pat_06")
    ap.add_argument("--band", default="beta")
    ap.add_argument("--mode", choices=("trace", "reset", "persist", "rearrange", "all"),
                    default="trace")
    args = ap.parse_args()

    modes = ("trace", "reset", "persist", "rearrange") if args.mode == "all" else (args.mode,)
    counts = {m: 0 for m in modes}
    for mode in modes:
        if args.all:
            for pat in PATIENTS:
                for band in BRAIN_BANDS_NAMES:
                    out = plot_one_cell(pat, band, mode=mode)
                    if out is not None:
                        counts[mode] += 1
        else:
            out = plot_one_cell(args.patient, args.band, mode=mode)
            if out is not None:
                counts[mode] += 1
    print("\n=== summary ===")
    for m, c in counts.items():
        n_total = len(PATIENTS) * len(BRAIN_BANDS_NAMES) if args.all else 1
        print(f"  {m}: {c} / {n_total} cells with at least 1 matching clade")


if __name__ == "__main__":
    main()
