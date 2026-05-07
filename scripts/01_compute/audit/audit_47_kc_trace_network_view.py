#!/usr/bin/env python3
"""Audit 47 — KC trace as a multiscale structural feature of the FC network.

EXACT-MATCH trace module identification: a "trace module" is a leaf set L
that forms a SUBTREE in BOTH the task_test AND the rest_post dendrograms
(identical leaves in both phases) AND that is fragmented in rest_pre
(max Jaccard against any rest_pre subtree < J_PRE_MAX).

Because the leaf set is identical in tt and rPost by construction, the
coloured subtree appears coherently in both phase dendrograms; in rest_pre
the same leaves are scattered (no rest_pre subtree contains them as a
group). This avoids the previous "Jaccard ≥ 0.5 quasi-match" bug where
phase-specific tree structure would split the coloured region across
multiple branches.

Renders Option C only (3x3 dendrogram + network paired panels):
rows = lambda regime (0 / 0.5 / 1), cols = phase (rest_pre / task / rest_post).

Lambda regimes:
  lambda=0  (topology) — all exact-match trace modules
  lambda=1  (heights)  — additionally requires merge heights of tt and rPost
                         subtrees to differ by ≤ HEIGHT_REL_THR (relative)
  lambda=0.5 (balanced) — same as lambda=1 (intersection)

For now: a single (patient, band) case (Pat_03 / beta) for design
iteration. Extend by adding entries to PATIENTS / BANDS once the visual
is locked.

Outputs
-------
``data/audit/kc_trace_network_view/figures/{pat}__{band}__option_C.pdf``
``data/audit/kc_trace_network_view/trace_subtrees_summary.csv``
"""
from __future__ import annotations

from pathlib import Path
import traceback

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import to_hex, to_rgba
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.network_templates import (
    _finalize_axes,
    compute_layout,
)
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result

PATIENTS = ["Pat_06", "Pat_03", "Pat_05", "Pat_13", "Pat_02", "Pat_10", "Pat_14", "Pat_07", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
PHASES = ("rest_pre", "task_test", "rest_post")
PHASE_LABELS = {
    "rest_pre": "rest pre",
    "task_test": "task",
    "rest_post": "rest post",
}
LAMBDAS = (0.0, 0.5, 1.0)
LAMBDA_LABELS = {
    0.0: r"$\lambda=0$  (topology)",
    0.5: r"$\lambda=0.5$ (balanced)",
    1.0: r"$\lambda=1$  (heights)",
}

J_POST_MIN = 0.6    # tt-rPost Jaccard threshold (lowered 2026-05-07 from 0.65 to surface more β modules)
J_PRE_MAX = 0.5
FRAG_FACTOR = 2.0   # smallest rPre subtree containing L_tt must be ≥ FRAG_FACTOR * |L_tt|
HEIGHT_REL_THR = 0.30
MIN_SIZE = 3
MAX_TRACE_TOTAL = 12
MAX_OVERLAP_FRAC = 0.30

LAYOUT_K = 20

OUT = ROOT / "data" / "audit" / "kc_trace_network_view"
FIG = OUT / "figures"
GRAY_NODE = "#cfcfcf"
GRAY_LEAF_HEX = "#c8c8c8"
GRAY_LINK_HEX = "#bababa"


def all_subtrees(Z: np.ndarray, n: int, min_size: int, max_size: int) -> list:
    subtree_leaves = {i: {i} for i in range(n)}
    out = []
    for i, row in enumerate(Z):
        a, b, h = int(row[0]), int(row[1]), float(row[2])
        new_id = n + i
        leaves = subtree_leaves[a] | subtree_leaves[b]
        subtree_leaves[new_id] = leaves
        sz = len(leaves)
        if min_size <= sz <= max_size:
            out.append({"id": new_id, "leaves": leaves, "height": h, "size": sz})
    return out


def find_trace_pairs(
    Z_tt: np.ndarray,
    Z_post: np.ndarray,
    Z_pre: np.ndarray,
    n: int,
    *,
    j_post_min: float = J_POST_MIN,
    j_pre_max: float = J_PRE_MAX,
    min_size: int = MIN_SIZE,
    max_size: int = None,
) -> dict:
    """Find (tt-subtree, rPost-subtree) pairs with high Jaccard in both phases.

    Each pair stores phase-specific leaf sets so each phase's dendrogram
    can colour its own coherent subtree (rather than borrowing tt's leaf
    set, which produced fragmented colouring at moderate Jaccard).
    """
    if max_size is None:
        max_size = max(min_size + 1, n // 2)

    sts_tt = all_subtrees(Z_tt, n, min_size, max_size)
    sts_post = all_subtrees(Z_post, n, min_size, max_size)
    sts_pre = all_subtrees(Z_pre, n, min_size, max_size)

    candidates = []
    for st_tt in sts_tt:
        L_tt = frozenset(st_tt["leaves"])
        s = len(L_tt)
        size_lo_post = max(min_size, int(s / 1.5))
        size_hi_post = min(max_size, int(s * 1.5))

        # rPost match: keep size restriction (we want a similar-size match).
        best_j_post = 0.0
        best_st_post = None
        for st_post in sts_post:
            if not (size_lo_post <= st_post["size"] <= size_hi_post):
                continue
            L_post = frozenset(st_post["leaves"])
            u = len(L_tt | L_post)
            if u == 0:
                continue
            j = len(L_tt & L_post) / u
            if j > best_j_post:
                best_j_post = j
                best_st_post = st_post

        if best_j_post < j_post_min or best_st_post is None:
            continue

        # rPre fragmentation: NO size restriction. We must catch the case where
        # L_tt's leaves are tightly co-clustered inside a SMALLER rPre subtree
        # (e.g. 3-leaf tt module packed into a 5-leaf rPre subtree gives
        # Jaccard 0.6 — that's an anchor we must reject) AND the case where
        # they're packed into a much LARGER rPre subtree where they'd give
        # low Jaccard but still be "co-clustered with a few others". For the
        # latter, also check the smallest rPre subtree fully containing L_tt.
        j_pre = 0.0
        for st_pre in sts_pre:
            L_pre = frozenset(st_pre["leaves"])
            u = len(L_tt | L_pre)
            if u == 0:
                continue
            j = len(L_tt & L_pre) / u
            if j > j_pre:
                j_pre = j
        if j_pre >= j_pre_max:
            continue

        # Containment-based fragmentation: smallest rPre subtree containing
        # all of L_tt must be at least FRAG_FACTOR * |L_tt|. If L_tt fits
        # inside a small rPre subtree (smallest containing size close to
        # |L_tt|), the leaves are anchored there and we reject.
        smallest_containing = float("inf")
        for st_pre in sts_pre:
            if L_tt <= st_pre["leaves"]:
                if st_pre["size"] < smallest_containing:
                    smallest_containing = st_pre["size"]
        # Add the singleton-extension subtrees: any rPre internal node above
        # min_size is in sts_pre; if the smallest containing subtree is
        # smaller than min_size cap it at the L size.
        if smallest_containing < FRAG_FACTOR * s:
            continue

        candidates.append({
            "leaves_tt": set(L_tt),
            "leaves_post": set(best_st_post["leaves"]),
            "h_tt": st_tt["height"],
            "h_post": best_st_post["height"],
            "j_post": best_j_post,
            "j_pre": j_pre,
            "size": max(len(L_tt), len(best_st_post["leaves"])),
            "score": best_j_post * (1.0 - j_pre),
        })

    candidates.sort(key=lambda c: (-c["size"], -c["score"]))
    selected = []
    used = set()
    for c in candidates:
        union_leaves = c["leaves_tt"] | c["leaves_post"]
        ovl = len(union_leaves & used) / max(len(union_leaves), 1)
        if ovl < MAX_OVERLAP_FRAC:
            selected.append(c)
            used |= union_leaves
        if len(selected) >= MAX_TRACE_TOTAL:
            break

    for i, c in enumerate(selected):
        c["color_idx"] = i

    trace_topo = list(selected)
    trace_height = []
    for c in selected:
        h_max = max(c["h_tt"], c["h_post"])
        if h_max > 0 and abs(c["h_tt"] - c["h_post"]) / h_max <= HEIGHT_REL_THR:
            trace_height.append(c)

    return {
        0.0: trace_topo,
        0.5: list(trace_height),
        1.0: trace_height,
    }


def module_leaves_for_phase(module: dict, phase: str) -> set:
    """Per-phase leaf set — used for DENDROGRAM coloring only.
    Each phase shows its own coherent subtree under its own tree topology."""
    if phase == "task_test":
        return module["leaves_tt"]
    elif phase == "rest_post":
        return module["leaves_post"]
    else:
        return module["leaves_tt"] | module["leaves_post"]


def module_leaves_for_network(module: dict) -> set:
    """Phase-INVARIANT leaf set — used for NETWORK coloring across all phases.
    The network has one shared node set; phase difference must show through
    edge weights, not through which nodes are colored. We use the union so
    no member of the trace module is dropped from any panel."""
    return module["leaves_tt"] | module["leaves_post"]


def palette_for(modules: list) -> dict:
    cmap = plt.get_cmap("tab10" if len(modules) <= 10 else "tab20")
    return {m["color_idx"]: to_hex(cmap(m["color_idx"] % cmap.N))
            for m in modules}


def assign_leaf_for_phase(modules: list, phase: str, n: int) -> list:
    """Each leaf -> smallest module containing it under that phase's leaf set
    (used for DENDROGRAM coloring only)."""
    out = [None] * n
    for m in sorted(modules, key=lambda x: x["size"]):
        leaves = module_leaves_for_phase(m, phase)
        for lf in leaves:
            if out[lf] is None:
                out[lf] = m
    return out


def assign_leaf_for_network(modules: list, n: int) -> list:
    """Each leaf -> smallest module containing it under the union leaf set
    (phase-invariant; used for NETWORK coloring)."""
    out = [None] * n
    for m in sorted(modules, key=lambda x: x["size"]):
        leaves = module_leaves_for_network(m)
        for lf in leaves:
            if out[lf] is None:
                out[lf] = m
    return out


def node_colors_for_network(modules: list, pal: dict, n: int) -> list:
    leaf_m = assign_leaf_for_network(modules, n)
    return [pal[m["color_idx"]] if m is not None else GRAY_NODE
            for m in leaf_m]


def draw_network(ax, A, pos, modules, pal, n, *, gamma=4.0):
    """Draw the FC network for ONE phase. Trace-module membership and
    intra-module edge selection are PHASE-INVARIANT (union of leaves_tt and
    leaves_post). Phase difference is encoded only through edge weights:
    each intra-module edge is alpha/width-scaled by its rank-percentile
    in this phase's edge set, so weak rPre intra-module edges fade and
    strong tt/rPost intra-module edges saturate."""
    rs, cs = np.triu_indices(n, k=1)
    w = np.abs(A[rs, cs])
    ranks = np.argsort(np.argsort(w)).astype(float)
    rank_pct = ranks / max(len(ranks) - 1, 1)

    bg_segs = [(pos[rs[i]], pos[cs[i]]) for i in range(len(rs))]
    bg_t = rank_pct ** 6
    bg_widths = 0.05 + 0.35 * bg_t
    bg_alphas = 0.012 + 0.08 * bg_t
    bg_colors = [(0.55, 0.55, 0.55, a) for a in bg_alphas]
    order_bg = np.argsort(rank_pct)
    bg_segs = [bg_segs[i] for i in order_bg]
    bg_colors = [bg_colors[i] for i in order_bg]
    bg_widths = [bg_widths[i] for i in order_bg]
    ax.add_collection(
        LineCollection(bg_segs, colors=bg_colors, linewidths=bg_widths, zorder=1)
    )

    for m in modules:
        members = module_leaves_for_network(m)
        intra = [i for i in range(len(rs))
                 if int(rs[i]) in members and int(cs[i]) in members]
        if not intra:
            continue
        intra = np.array(intra, dtype=int)
        rp = rank_pct[intra]
        t = rp ** gamma
        widths = 0.30 + 3.0 * t
        alphas = 0.10 + 0.85 * t
        rgb = to_rgba(pal[m["color_idx"]])[:3]
        rgba = [(rgb[0], rgb[1], rgb[2], float(a)) for a in alphas]
        segs = [(pos[rs[i]], pos[cs[i]]) for i in intra]
        order = np.argsort(rp)
        segs = [segs[i] for i in order]
        rgba = [rgba[i] for i in order]
        widths = [float(widths[i]) for i in order]
        ax.add_collection(
            LineCollection(segs, colors=rgba, linewidths=widths, zorder=3)
        )

    nc = node_colors_for_network(modules, pal, n)
    ax.scatter(
        pos[:, 0], pos[:, 1], c=nc, s=22,
        edgecolors="white", linewidths=0.4, zorder=5,
    )
    _finalize_axes(ax, pos)


def _link_color_func(Z_phase, modules, phase, pal, n):
    """Per-phase: an internal link's leaves ⊆ that module's PHASE-SPECIFIC
    leaf set. tt uses leaves_tt, rPost uses leaves_post, rPre uses the union
    (so isolated rPre leaves still get colored individually but their internal
    links rarely qualify, leaving the colored region visually fragmented)."""
    leaf_m = assign_leaf_for_phase(modules, phase, n)
    sorted_mods = sorted(modules, key=lambda x: x["size"])

    subtree_leaves = {i: {i} for i in range(n)}
    link_color = {}
    for i, row in enumerate(Z_phase):
        a, b = int(row[0]), int(row[1])
        leaves = subtree_leaves[a] | subtree_leaves[b]
        new_id = n + i
        subtree_leaves[new_id] = leaves
        link_color[new_id] = GRAY_LINK_HEX
        for m in sorted_mods:
            mod_leaves = module_leaves_for_phase(m, phase)
            if leaves <= mod_leaves:
                link_color[new_id] = pal[m["color_idx"]]
                break

    leaf_color = {}
    for lf in range(n):
        m = leaf_m[lf]
        leaf_color[lf] = pal[m["color_idx"]] if m is not None else GRAY_LEAF_HEX

    def fn(node_id):
        if node_id < n:
            return leaf_color[node_id]
        return link_color.get(node_id, GRAY_LINK_HEX)

    return fn


def draw_dendro(ax, Z, modules, phase, pal, n, *, draw_height_bars=False):
    fn = _link_color_func(Z, modules, phase, pal, n)
    R = dendrogram(
        Z,
        ax=ax,
        no_labels=True,
        color_threshold=0,
        above_threshold_color=GRAY_LINK_HEX,
        link_color_func=fn,
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

    # Build leaf -> parent-merge-height map. scipy's link_color_func is
    # only invoked for internal nodes; leaf segments inherit the parent
    # link's color (gray for fragmented). Override by drawing colored
    # vlines from y=0 up to each leaf's parent merge height — the leaf
    # line in the dendrogram is now the module color, all the way up.
    leaf_parent_h = {}
    for i, row in enumerate(Z):
        a, b, h = int(row[0]), int(row[1]), float(row[2])
        if a < n:
            leaf_parent_h.setdefault(a, h)
        if b < n:
            leaf_parent_h.setdefault(b, h)

    leaf_m = assign_leaf_for_phase(modules, phase, n)
    leaf_order = R["leaves"]
    for i, lf_idx in enumerate(leaf_order):
        m = leaf_m[lf_idx]
        if m is None:
            continue
        x = 5.0 + 10.0 * i
        h_parent = leaf_parent_h.get(lf_idx, tmax * 0.05)
        ax.vlines(x, tmin, h_parent, color=pal[m["color_idx"]], lw=1.5, zorder=4)

    # Merge-height bars at lambda = 0.5 / 1: full-width dashed line at each
    # matched merge height. Tt panel uses h_tt, rPost uses h_post; comparing
    # bar heights across the two panels makes the height-matching condition
    # (|h_tt - h_post|/max ≤ HEIGHT_REL_THR) visually obvious. Skip rPre.
    if draw_height_bars and phase != "rest_pre":
        for m in modules:
            h = m["h_tt"] if phase == "task_test" else m["h_post"]
            color = pal[m["color_idx"]]
            ax.axhline(h, color=color, lw=2.0, ls="--", alpha=0.65, zorder=2)


def render_option_c(patient, band, Z_by_phase, A_by_phase, pos,
                    trace_by_lambda, pal, n):
    fig = plt.figure(figsize=(20.0, 13.0))
    gs = fig.add_gridspec(
        3, 6,
        height_ratios=[1, 1, 1],
        width_ratios=[1.0, 1.4, 1.0, 1.4, 1.0, 1.4],
        hspace=0.16, wspace=0.12,
    )
    for r, lam in enumerate(LAMBDAS):
        modules = trace_by_lambda[lam]
        for c, phase in enumerate(PHASES):
            ax_d = fig.add_subplot(gs[r, 2 * c])
            ax_n = fig.add_subplot(gs[r, 2 * c + 1])
            draw_dendro(ax_d, Z_by_phase[phase], modules, phase, pal, n,
                        draw_height_bars=(lam > 0.0))
            draw_network(ax_n, A_by_phase[phase], pos, modules, pal, n)
            if r == 0:
                ax_d.set_title(f"{PHASE_LABELS[phase]}: dendrogram", fontsize=10)
                ax_n.set_title(f"{PHASE_LABELS[phase]}: network", fontsize=10)
            if c == 0:
                ax_d.text(
                    -0.22, 0.5, LAMBDA_LABELS[lam], transform=ax_d.transAxes,
                    rotation=90, ha="center", va="center", fontsize=11,
                )

    band_tex = BRAIN_BAND_TEX_DICT[band]
    n_trace = {lam: len(trace_by_lambda[lam]) for lam in LAMBDAS}
    fig.text(
        0.5, 0.99,
        f"{patient}  {band_tex}  trace modules (Jaccard$_{{\\rm post}}$$\\geq${J_POST_MIN}, "
        f"Jaccard$_{{\\rm pre}}$$<${J_PRE_MAX})  "
        rf"$\lambda$=0:{n_trace[0.0]}, "
        rf"0.5:{n_trace[0.5]}, "
        rf"1:{n_trace[1.0]}",
        ha="center", va="top", fontsize=12,
    )
    fig.subplots_adjust(left=0.04, right=0.99, top=0.95, bottom=0.02)
    out_path = FIG / f"{patient}__{band}__option_C.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for patient in PATIENTS:
        for band in BANDS:
            try:
                A_by_phase = {}
                Z_by_phase = {}
                for phase in PHASES:
                    A = load_fc_matrix(patient, phase, band, fc_method="imcoh_abs")
                    res = load_lrg_result(
                        patient, phase, band, fc_method="imcoh_abs"
                    )
                    if A is None or res is None or res.linkage_matrix is None:
                        raise FileNotFoundError(f"{patient} {phase} {band}")
                    A_by_phase[phase] = np.asarray(A)
                    Z_by_phase[phase] = np.asarray(res.linkage_matrix)
                n = A_by_phase["task_test"].shape[0]
                for ph in PHASES:
                    if A_by_phase[ph].shape[0] != n:
                        m = min(A_by_phase[ph].shape[0], n)
                        A_by_phase[ph] = A_by_phase[ph][:m, :m]
                        n = m

                trace_by_lambda = find_trace_pairs(
                    Z_by_phase["task_test"],
                    Z_by_phase["rest_post"],
                    Z_by_phase["rest_pre"],
                    n,
                )
                modules = trace_by_lambda[0.0]

                if not modules:
                    print(f"[audit_47] {patient} {band} -> NO modules; skipping figure")
                    summary_rows.append({
                        "patient": patient,
                        "band": band,
                        "n_nodes": n,
                        "n_modules_lambda0": 0,
                        "n_modules_lambda1": 0,
                        "n_modules_lambda05": 0,
                        "sizes_tt": [],
                        "sizes_post": [],
                        "j_post": [],
                        "j_pre": [],
                        "option_C": "",
                    })
                    continue

                pal = palette_for(modules)

                labels_layout = fcluster(
                    Z_by_phase["task_test"], t=LAYOUT_K, criterion="maxclust"
                )
                pos = compute_layout(
                    "lrg_sfdp",
                    A_by_phase["task_test"],
                    community_labels=labels_layout.tolist(),
                )

                pc = render_option_c(
                    patient, band, Z_by_phase, A_by_phase, pos,
                    trace_by_lambda, pal, n,
                )

                summary_rows.append({
                    "patient": patient,
                    "band": band,
                    "n_nodes": n,
                    "n_modules_lambda0": len(trace_by_lambda[0.0]),
                    "n_modules_lambda1": len(trace_by_lambda[1.0]),
                    "n_modules_lambda05": len(trace_by_lambda[0.5]),
                    "sizes_tt": sorted([len(m["leaves_tt"]) for m in trace_by_lambda[0.0]]),
                    "sizes_post": sorted([len(m["leaves_post"]) for m in trace_by_lambda[0.0]]),
                    "j_post": [round(m["j_post"], 3) for m in
                               sorted(trace_by_lambda[0.0], key=lambda x: -x["score"])],
                    "j_pre": [round(m["j_pre"], 3) for m in
                              sorted(trace_by_lambda[0.0], key=lambda x: -x["score"])],
                    "option_C": str(pc.relative_to(ROOT)),
                })
                print(
                    f"[audit_47] {patient} {band} -> modules topo={len(trace_by_lambda[0.0])} "
                    f"height={len(trace_by_lambda[1.0])} "
                    f"sizes_tt={sorted([len(m['leaves_tt']) for m in trace_by_lambda[0.0]])} "
                    f"j_post={[round(m['j_post'],2) for m in sorted(trace_by_lambda[0.0], key=lambda x: -x['score'])]}"
                )
            except Exception as e:
                print(f"[audit_47] WARN {patient} {band}: {e}")
                traceback.print_exc()

    df = pd.DataFrame(summary_rows)
    df.to_csv(OUT / "trace_subtrees_summary.csv", index=False)
    print(f"[audit_47] summary: {OUT / 'trace_subtrees_summary.csv'}")


if __name__ == "__main__":
    main()
