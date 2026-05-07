#!/usr/bin/env python3
"""Audit 49 — KC rearrangement modules: rPost-only emergent subtrees.

Counter-example to trace. A rearrangement module is a leaf set L_post
that is coherent ONLY in rest_post — fragmented in BOTH rest_pre AND
task_test. Anchored at rPost subtrees:

  Gate 1: rPre fragment    forall pre subtree T: J(L_post, T) < J_DISRUPT = 0.5
  Gate 2: rPre containment |C_pre(L_post)| >= FRAG_FACTOR * |L_post|
  Gate 3: tt  fragment     forall tt  subtree T: J(L_post, T) < J_DISRUPT = 0.5
  Gate 4: tt  containment  |C_tt(L_post)|  >= FRAG_FACTOR * |L_post|
  Gate 5: size floor       |L_post| >= MIN_SIZE = 3

No matched pair in rPre or tt — both fragment. KC lambda regimes do
not apply (no height-comparison target). Option C 3x3 grid: rows
are size strata (size >= 3 / >= 5 / >= 8) instead of lambda regimes.
The visualization is the cohort control on the trace claim: a small
rearrangement count supports "rPost coherence is task-locked"; a
large count weakens it.

Empty-cell skip: if no modules, no figure is rendered.
Linear y-axis on dendrograms.

Outputs
-------
``data/audit/kc_rearrangement_module_view/figures/{pat}__{band}__option_C.pdf``
``data/audit/kc_rearrangement_module_view/rearrangement_subtrees_summary.csv``

Scope: ``.agents/guides/task-persistence-investigation/2026-05-07_kc-rearrangement-modules.md``

NOTE: audit_47 + audit_48 + audit_49 share substantial scaffolding.
Library promotion to ``lrg_eegfc.visuals.kc_module_view`` is the planned
follow-up refactor.
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

J_DISRUPT = 0.5
FRAG_FACTOR = 2.0
MIN_SIZE = 3
MAX_REARRANGE_TOTAL = 12
MAX_OVERLAP_FRAC = 0.30

LAYOUT_K = 20

OUT = ROOT / "data" / "audit" / "kc_rearrangement_module_view"
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


def find_rearrangement_modules(
    Z_pre: np.ndarray,
    Z_tt: np.ndarray,
    Z_post: np.ndarray,
    n: int,
    *,
    j_disrupt: float = J_DISRUPT,
    min_size: int = MIN_SIZE,
    max_size: int = None,
) -> list:
    """Find rPost-only coherent subtrees (rPre + tt both fragment).

    Anchored at rPost. Returns a list of dicts with leaves_post and
    h_post; no matched pair in either rPre or tt.
    """
    if max_size is None:
        max_size = max(min_size + 1, n // 2)

    sts_pre = all_subtrees(Z_pre, n, min_size, max_size)
    sts_tt = all_subtrees(Z_tt, n, min_size, max_size)
    sts_post = all_subtrees(Z_post, n, min_size, max_size)

    candidates = []
    for st_post in sts_post:
        L_post = frozenset(st_post["leaves"])
        s = len(L_post)

        # Gate 1: rPre fragmentation, leaf-set side.
        j_pre = 0.0
        for st_pre in sts_pre:
            L_pre = frozenset(st_pre["leaves"])
            u = len(L_post | L_pre)
            if u == 0:
                continue
            j = len(L_post & L_pre) / u
            if j > j_pre:
                j_pre = j
        if j_pre >= j_disrupt:
            continue

        # Gate 2: rPre containment fragmentation.
        smallest_containing_pre = float("inf")
        for st_pre in sts_pre:
            if L_post <= st_pre["leaves"]:
                if st_pre["size"] < smallest_containing_pre:
                    smallest_containing_pre = st_pre["size"]
        if smallest_containing_pre < FRAG_FACTOR * s:
            continue

        # Gate 3: tt fragmentation, leaf-set side.
        j_tt = 0.0
        for st_tt in sts_tt:
            L_tt = frozenset(st_tt["leaves"])
            u = len(L_post | L_tt)
            if u == 0:
                continue
            j = len(L_post & L_tt) / u
            if j > j_tt:
                j_tt = j
        if j_tt >= j_disrupt:
            continue

        # Gate 4: tt containment fragmentation.
        smallest_containing_tt = float("inf")
        for st_tt in sts_tt:
            if L_post <= st_tt["leaves"]:
                if st_tt["size"] < smallest_containing_tt:
                    smallest_containing_tt = st_tt["size"]
        if smallest_containing_tt < FRAG_FACTOR * s:
            continue

        candidates.append({
            "leaves_post": set(L_post),
            "h_post": st_post["height"],
            "j_pre": j_pre,
            "j_tt": j_tt,
            "size": s,
            "score": (1.0 - j_pre) * (1.0 - j_tt),
        })

    candidates.sort(key=lambda c: (-c["size"], -c["score"]))
    selected = []
    used = set()
    for c in candidates:
        ovl = len(c["leaves_post"] & used) / max(len(c["leaves_post"]), 1)
        if ovl < MAX_OVERLAP_FRAC:
            selected.append(c)
            used |= c["leaves_post"]
        if len(selected) >= MAX_REARRANGE_TOTAL:
            break

    for i, c in enumerate(selected):
        c["color_idx"] = i

    return selected


def module_leaves(module: dict) -> set:
    return module["leaves_post"]


def palette_for(modules: list) -> dict:
    cmap = plt.get_cmap("tab10" if len(modules) <= 10 else "tab20")
    return {m["color_idx"]: to_hex(cmap(m["color_idx"] % cmap.N))
            for m in modules}


def assign_leaf(modules: list, n: int) -> list:
    out = [None] * n
    for m in sorted(modules, key=lambda x: x["size"]):
        for lf in module_leaves(m):
            if out[lf] is None:
                out[lf] = m
    return out


def node_colors(modules: list, pal: dict, n: int) -> list:
    leaf_m = assign_leaf(modules, n)
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
        members = module_leaves(m)
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

    nc = node_colors(modules, pal, n)
    ax.scatter(
        pos[:, 0], pos[:, 1], c=nc, s=22,
        edgecolors="white", linewidths=0.4, zorder=5,
    )
    _finalize_axes(ax, pos)


def _link_color_func(Z_phase, modules, phase, pal, n):
    """Phase-aware coloring. rPost shows coherent colored subtree;
    rPre and tt show only scattered leaf bars (parent links remain gray
    because the leaves are fragmented across multiple internal subtrees)."""
    leaf_m = assign_leaf(modules, n)
    sorted_mods = sorted(modules, key=lambda x: x["size"])

    subtree_leaves = {i: {i} for i in range(n)}
    link_color = {}
    for i, row in enumerate(Z_phase):
        a, b = int(row[0]), int(row[1])
        leaves = subtree_leaves[a] | subtree_leaves[b]
        new_id = n + i
        subtree_leaves[new_id] = leaves
        link_color[new_id] = GRAY_LINK_HEX
        if phase == "rest_post":
            for m in sorted_mods:
                if leaves <= m["leaves_post"]:
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

    leaf_m = assign_leaf(modules, n)
    leaf_order = R["leaves"]
    for i, lf_idx in enumerate(leaf_order):
        m = leaf_m[lf_idx]
        if m is None:
            continue
        x = 5.0 + 10.0 * i
        h_parent = leaf_parent_h.get(lf_idx, tmax * 0.05)
        ax.vlines(x, tmin, h_parent, color=pal[m["color_idx"]], lw=1.5, zorder=4)

    # Merge-height bar at h_post in rPost dendrogram only.
    if draw_height_bars and phase == "rest_post":
        for m in modules:
            color = pal[m["color_idx"]]
            ax.axhline(m["h_post"], color=color, lw=2.0, ls="--", alpha=0.65, zorder=2)


# Row stratification for the 3x3 option-C grid. Rearrangement has no
# matched-pair lambda regime, so rows are size strata instead:
#   row 0: all modules (size >= 3)
#   row 1: medium+ modules (size >= 5)
#   row 2: large modules (size >= 8)
ROW_LABELS = (
    r"size $\geq$ 3 (all)",
    r"size $\geq$ 5",
    r"size $\geq$ 8",
)
ROW_MIN_SIZES = (3, 5, 8)


def render_option_c(patient, band, Z_by_phase, A_by_phase, pos,
                    modules, pal, n):
    fig = plt.figure(figsize=(20.0, 13.0))
    gs = fig.add_gridspec(
        3, 6,
        height_ratios=[1, 1, 1],
        width_ratios=[1.0, 1.4, 1.0, 1.4, 1.0, 1.4],
        hspace=0.16, wspace=0.12,
    )
    for r, min_size in enumerate(ROW_MIN_SIZES):
        row_modules = [m for m in modules if m["size"] >= min_size]
        for c, phase in enumerate(PHASES):
            ax_d = fig.add_subplot(gs[r, 2 * c])
            ax_n = fig.add_subplot(gs[r, 2 * c + 1])
            draw_dendro(ax_d, Z_by_phase[phase], row_modules, phase, pal, n,
                        draw_height_bars=True)
            draw_network(ax_n, A_by_phase[phase], pos, row_modules, pal, n)
            if r == 0:
                ax_d.set_title(f"{PHASE_LABELS[phase]}: dendrogram", fontsize=10)
                ax_n.set_title(f"{PHASE_LABELS[phase]}: network", fontsize=10)
            if c == 0:
                ax_d.text(
                    -0.22, 0.5, ROW_LABELS[r], transform=ax_d.transAxes,
                    rotation=90, ha="center", va="center", fontsize=11,
                )

    band_tex = BRAIN_BAND_TEX_DICT[band]
    n_per_row = [len([m for m in modules if m["size"] >= s]) for s in ROW_MIN_SIZES]
    fig.text(
        0.5, 0.99,
        f"{patient}  {band_tex}  rearrangement modules "
        f"(Jaccard$_{{\\rm pre}}$, Jaccard$_{{\\rm tt}}$ both $<${J_DISRUPT})  "
        f"n={n_per_row[0]}, $\\geq$5:{n_per_row[1]}, $\\geq$8:{n_per_row[2]}",
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
                n = A_by_phase["rest_post"].shape[0]
                for ph in PHASES:
                    if A_by_phase[ph].shape[0] != n:
                        m = min(A_by_phase[ph].shape[0], n)
                        A_by_phase[ph] = A_by_phase[ph][:m, :m]
                        n = m

                modules = find_rearrangement_modules(
                    Z_by_phase["rest_pre"],
                    Z_by_phase["task_test"],
                    Z_by_phase["rest_post"],
                    n,
                )

                if not modules:
                    print(f"[audit_49] {patient} {band} -> NO modules; skipping figure")
                    summary_rows.append({
                        "patient": patient,
                        "band": band,
                        "n_nodes": n,
                        "n_modules": 0,
                        "sizes_post": [],
                        "j_pre": [],
                        "j_tt": [],
                        "figure": "",
                    })
                    continue

                pal = palette_for(modules)

                # Layout seeded with rPost partition (the anchor phase for rearrangement).
                labels_layout = fcluster(
                    Z_by_phase["rest_post"], t=LAYOUT_K, criterion="maxclust"
                )
                pos = compute_layout(
                    "lrg_sfdp",
                    A_by_phase["rest_post"],
                    community_labels=labels_layout.tolist(),
                )

                pc = render_option_c(
                    patient, band, Z_by_phase, A_by_phase, pos,
                    modules, pal, n,
                )

                summary_rows.append({
                    "patient": patient,
                    "band": band,
                    "n_nodes": n,
                    "n_modules": len(modules),
                    "sizes_post": sorted([m["size"] for m in modules]),
                    "j_pre": [round(m["j_pre"], 3) for m in
                              sorted(modules, key=lambda x: -x["score"])],
                    "j_tt": [round(m["j_tt"], 3) for m in
                             sorted(modules, key=lambda x: -x["score"])],
                    "figure": str(pc.relative_to(ROOT)),
                })
                print(
                    f"[audit_49] {patient} {band} -> rearrangement modules={len(modules)} "
                    f"sizes={sorted([m['size'] for m in modules])} "
                    f"j_pre={[round(m['j_pre'],2) for m in sorted(modules, key=lambda x: -x['score'])]} "
                    f"j_tt={[round(m['j_tt'],2) for m in sorted(modules, key=lambda x: -x['score'])]}"
                )
            except Exception as e:
                print(f"[audit_49] WARN {patient} {band}: {e}")
                traceback.print_exc()

    df = pd.DataFrame(summary_rows)
    df.to_csv(OUT / "rearrangement_subtrees_summary.csv", index=False)
    print(f"[audit_49] summary: {OUT / 'rearrangement_subtrees_summary.csv'}")


if __name__ == "__main__":
    main()
