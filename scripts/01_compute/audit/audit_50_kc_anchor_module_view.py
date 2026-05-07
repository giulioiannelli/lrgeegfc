#!/usr/bin/env python3
"""Audit 50 — KC anchor modules: coherent across all three phases.

Fourth class in the cross-phase taxonomy (alongside trace / reset /
rearrangement). User term: "persist". An anchor module is a leaf
set that has a matching coherent subtree in EVERY phase — rest_pre,
task, AND rest_post. Anchors are the **null** of the trace story:
modules that the task does NOT touch.

  Gate 1: rPre↔tt   J(L_pre, L_tt) >= J_ANCHOR = 0.6
  Gate 2: rPre↔rPost J(L_pre, L_post) >= J_ANCHOR = 0.6
  Gate 3: tt↔rPost  J(L_tt, L_post) >= J_ANCHOR = 0.6
  Gate 4: size match all three within /1.5..*1.5 of each other
  Gate 5: size floor |L_pre| >= MIN_SIZE = 3

Maximizer pick: among (u_tt, u_post) candidates per u_pre, choose the
pair maximizing geometric mean of the three pairwise Jaccards.

Renders Option C (3x3 dendrogram + network paired panels):
rows = lambda regime (0 / 0.5 / 1), cols = phase (rest_pre / task / rest_post).

All three dendrograms color the matched subtree coherently. Networks
look similar across phases (modules don't reorganize). Merge-height
bars at lambda = 0.5 / 1 in all three phases.

Empty-cell skip + linear y-axis (audit_47 conventions).

Outputs
-------
``data/audit/kc_anchor_module_view/figures/{pat}__{band}__option_C.pdf``
``data/audit/kc_anchor_module_view/anchor_subtrees_summary.csv``

Scope: ``.agents/guides/task-persistence-investigation/2026-05-07_kc-anchor-modules.md``
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

PATIENTS = ["Pat_06", "Pat_03", "Pat_05", "Pat_13", "Pat_02", "Pat_10", "Pat_14", "Pat_07", "Pat_15", "Pat_08"]
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

J_ANCHOR = 0.6
HEIGHT_REL_THR = 0.30
MIN_SIZE = 3
SIZE_RATIO = 1.5
MAX_ANCHOR_TOTAL = 12
MAX_OVERLAP_FRAC = 0.30

LAYOUT_K = 20

OUT = ROOT / "data" / "audit" / "kc_anchor_module_view"
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


def find_anchor_triples(
    Z_pre: np.ndarray,
    Z_tt: np.ndarray,
    Z_post: np.ndarray,
    n: int,
    *,
    j_anchor: float = J_ANCHOR,
    min_size: int = MIN_SIZE,
    max_size: int = None,
) -> dict:
    """Find triple matches across rPre / tt / rPost at Jaccard >= j_anchor."""
    if max_size is None:
        max_size = max(min_size + 1, n // 2)

    sts_pre = all_subtrees(Z_pre, n, min_size, max_size)
    sts_tt = all_subtrees(Z_tt, n, min_size, max_size)
    sts_post = all_subtrees(Z_post, n, min_size, max_size)

    candidates = []
    for st_pre in sts_pre:
        L_pre = frozenset(st_pre["leaves"])
        s = len(L_pre)
        s_lo = max(min_size, int(s / SIZE_RATIO))
        s_hi = min(max_size, int(s * SIZE_RATIO))

        best_score = 0.0
        best_st_tt = None
        best_st_post = None
        best_j_pre_tt = best_j_pre_post = best_j_tt_post = 0.0

        for st_tt in sts_tt:
            if not (s_lo <= st_tt["size"] <= s_hi):
                continue
            L_tt = frozenset(st_tt["leaves"])
            u_pt = len(L_pre | L_tt)
            if u_pt == 0:
                continue
            j_pre_tt = len(L_pre & L_tt) / u_pt
            if j_pre_tt < j_anchor:
                continue

            for st_post in sts_post:
                if not (s_lo <= st_post["size"] <= s_hi):
                    continue
                L_post = frozenset(st_post["leaves"])
                u_pp = len(L_pre | L_post)
                if u_pp == 0:
                    continue
                j_pre_post = len(L_pre & L_post) / u_pp
                if j_pre_post < j_anchor:
                    continue
                u_tp = len(L_tt | L_post)
                if u_tp == 0:
                    continue
                j_tt_post = len(L_tt & L_post) / u_tp
                if j_tt_post < j_anchor:
                    continue

                score = (j_pre_tt * j_pre_post * j_tt_post) ** (1.0 / 3.0)
                if score > best_score:
                    best_score = score
                    best_st_tt = st_tt
                    best_st_post = st_post
                    best_j_pre_tt = j_pre_tt
                    best_j_pre_post = j_pre_post
                    best_j_tt_post = j_tt_post

        if best_st_tt is None:
            continue

        candidates.append({
            "leaves_pre": set(L_pre),
            "leaves_tt": set(best_st_tt["leaves"]),
            "leaves_post": set(best_st_post["leaves"]),
            "h_pre": st_pre["height"],
            "h_tt": best_st_tt["height"],
            "h_post": best_st_post["height"],
            "j_pre_tt": best_j_pre_tt,
            "j_pre_post": best_j_pre_post,
            "j_tt_post": best_j_tt_post,
            "size": s,
            "score": best_score,
        })

    candidates.sort(key=lambda c: (-c["size"], -c["score"]))
    selected = []
    used = set()
    for c in candidates:
        union_leaves = c["leaves_pre"] | c["leaves_tt"] | c["leaves_post"]
        ovl = len(union_leaves & used) / max(len(union_leaves), 1)
        if ovl < MAX_OVERLAP_FRAC:
            selected.append(c)
            used |= union_leaves
        if len(selected) >= MAX_ANCHOR_TOTAL:
            break

    for i, c in enumerate(selected):
        c["color_idx"] = i

    anchor_topo = list(selected)
    anchor_height = []
    for c in selected:
        h_max = max(c["h_pre"], c["h_tt"], c["h_post"])
        h_min = min(c["h_pre"], c["h_tt"], c["h_post"])
        if h_max > 0 and (h_max - h_min) / h_max <= HEIGHT_REL_THR:
            anchor_height.append(c)

    return {
        0.0: anchor_topo,
        0.5: list(anchor_height),
        1.0: anchor_height,
    }


def module_leaves_for_phase(module: dict, phase: str) -> set:
    if phase == "rest_pre":
        return module["leaves_pre"]
    elif phase == "task_test":
        return module["leaves_tt"]
    else:  # rest_post
        return module["leaves_post"]


def module_leaves_for_network(module: dict) -> set:
    return module["leaves_pre"] | module["leaves_tt"] | module["leaves_post"]


def palette_for(modules: list) -> dict:
    cmap = plt.get_cmap("tab10" if len(modules) <= 10 else "tab20")
    return {m["color_idx"]: to_hex(cmap(m["color_idx"] % cmap.N))
            for m in modules}


def assign_leaf_for_phase(modules: list, phase: str, n: int) -> list:
    out = [None] * n
    for m in sorted(modules, key=lambda x: x["size"]):
        leaves = module_leaves_for_phase(m, phase)
        for lf in leaves:
            if out[lf] is None:
                out[lf] = m
    return out


def assign_leaf_for_network(modules: list, n: int) -> list:
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

    # Anchor: all three phases have a matched height. Draw the bar in
    # all three dendrograms (rPre at h_pre, tt at h_tt, rPost at h_post).
    if draw_height_bars:
        for m in modules:
            h = (m["h_pre"] if phase == "rest_pre"
                 else m["h_tt"] if phase == "task_test"
                 else m["h_post"])
            color = pal[m["color_idx"]]
            ax.axhline(h, color=color, lw=2.0, ls="--", alpha=0.65, zorder=2)


def render_option_c(patient, band, Z_by_phase, A_by_phase, pos,
                    anchor_by_lambda, pal, n):
    fig = plt.figure(figsize=(20.0, 13.0))
    gs = fig.add_gridspec(
        3, 6,
        height_ratios=[1, 1, 1],
        width_ratios=[1.0, 1.4, 1.0, 1.4, 1.0, 1.4],
        hspace=0.16, wspace=0.12,
    )
    for r, lam in enumerate(LAMBDAS):
        modules = anchor_by_lambda[lam]
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
    n_anchor = {lam: len(anchor_by_lambda[lam]) for lam in LAMBDAS}
    fig.text(
        0.5, 0.99,
        f"{patient}  {band_tex}  anchor (persist) modules "
        f"(Jaccard$_{{\\rm pairwise}}$$\\geq${J_ANCHOR})  "
        rf"$\lambda$=0:{n_anchor[0.0]}, "
        rf"0.5:{n_anchor[0.5]}, "
        rf"1:{n_anchor[1.0]}",
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
                n = A_by_phase["rest_pre"].shape[0]
                for ph in PHASES:
                    if A_by_phase[ph].shape[0] != n:
                        m = min(A_by_phase[ph].shape[0], n)
                        A_by_phase[ph] = A_by_phase[ph][:m, :m]
                        n = m

                anchor_by_lambda = find_anchor_triples(
                    Z_by_phase["rest_pre"],
                    Z_by_phase["task_test"],
                    Z_by_phase["rest_post"],
                    n,
                )
                modules = anchor_by_lambda[0.0]

                if not modules:
                    print(f"[audit_50] {patient} {band} -> NO modules; skipping figure")
                    summary_rows.append({
                        "patient": patient,
                        "band": band,
                        "n_nodes": n,
                        "n_modules_lambda0": 0,
                        "n_modules_lambda1": 0,
                        "n_modules_lambda05": 0,
                        "sizes_pre": [],
                        "j_pre_tt": [],
                        "j_pre_post": [],
                        "j_tt_post": [],
                        "option_C": "",
                    })
                    continue

                pal = palette_for(modules)

                # Layout seeded with rPre partition (anchor baseline).
                labels_layout = fcluster(
                    Z_by_phase["rest_pre"], t=LAYOUT_K, criterion="maxclust"
                )
                pos = compute_layout(
                    "lrg_sfdp",
                    A_by_phase["rest_pre"],
                    community_labels=labels_layout.tolist(),
                )

                pc = render_option_c(
                    patient, band, Z_by_phase, A_by_phase, pos,
                    anchor_by_lambda, pal, n,
                )

                summary_rows.append({
                    "patient": patient,
                    "band": band,
                    "n_nodes": n,
                    "n_modules_lambda0": len(anchor_by_lambda[0.0]),
                    "n_modules_lambda1": len(anchor_by_lambda[1.0]),
                    "n_modules_lambda05": len(anchor_by_lambda[0.5]),
                    "sizes_pre": sorted([len(m["leaves_pre"]) for m in anchor_by_lambda[0.0]]),
                    "j_pre_tt": [round(m["j_pre_tt"], 3) for m in
                                 sorted(anchor_by_lambda[0.0], key=lambda x: -x["score"])],
                    "j_pre_post": [round(m["j_pre_post"], 3) for m in
                                   sorted(anchor_by_lambda[0.0], key=lambda x: -x["score"])],
                    "j_tt_post": [round(m["j_tt_post"], 3) for m in
                                  sorted(anchor_by_lambda[0.0], key=lambda x: -x["score"])],
                    "option_C": str(pc.relative_to(ROOT)),
                })
                print(
                    f"[audit_50] {patient} {band} -> anchor modules topo={len(anchor_by_lambda[0.0])} "
                    f"height={len(anchor_by_lambda[1.0])} "
                    f"sizes_pre={sorted([len(m['leaves_pre']) for m in anchor_by_lambda[0.0]])}"
                )
            except Exception as e:
                print(f"[audit_50] WARN {patient} {band}: {e}")
                traceback.print_exc()

    df = pd.DataFrame(summary_rows)
    df.to_csv(OUT / "anchor_subtrees_summary.csv", index=False)
    print(f"[audit_50] summary: {OUT / 'anchor_subtrees_summary.csv'}")


if __name__ == "__main__":
    main()
