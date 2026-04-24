#!/usr/bin/env python3
"""Diagnostic: visualize π-identified task-born-retained modules.

Pick (patient, band) with the strongest S_CBR in that band; identify the
top-K rest_post internal nodes scored by π(η) · (match_tt − match_rpre);
plot the four phase dendrograms side by side with matching subtrees
highlighted in shared colors across phases.

The key diagnostic: in rest_post and task_test the highlighted module
should appear as a single coherent subtree; in rest_pre the same leaves
should fragment across multiple disjoint sub-subtrees — showing the
module did not pre-exist. If the matched tt-subtree sits at a very
different h_rel than the rpost-subtree, the Jaccard-only match is
incomplete and we would need a height-agreement factor (read: this
figure decides it).
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.utils.metrics import jaccard_leafsets
from lrg_eegfc.workflow.lrg import load_lrg_result


PATIENT = "Pat_06"
BAND = "delta"
K_MODULES = 4
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHASE_TEX = {"rest_pre": "rest$_\\mathrm{pre}$",
             "task_learn": "task$_\\mathrm{learn}$",
             "task_test": "task$_\\mathrm{test}$",
             "rest_post": "rest$_\\mathrm{post}$"}
MODULE_COLORS = ["#e63946", "#1d3557", "#2a9d8f", "#f4a261",
                 "#8338ec", "#6a994e"]
BACKGROUND_GRAY = "#bfbfbf"


# ──────────────────────────── tree utilities ────────────────────────────

def nodes_with_prominence(Z: np.ndarray, h_floor_frac: float = 1e-6) -> list[dict]:
    """Internal nodes indexed by Z row, each with id, h, size, leaves, π.

    π(η) = log(h_parent) − log(h_self). Root has parent_row = -1 and π = nan.
    h_self is floored at ``h_floor_frac · D_max`` to avoid log(0) for
    near-identical leaf pairs.
    """
    from scipy.cluster.hierarchy import to_tree
    n_leaves = Z.shape[0] + 1
    n_int = Z.shape[0]
    root = to_tree(Z, rd=False)
    nodes: list[dict | None] = [None] * n_int

    def walk(nd) -> frozenset[int]:
        if nd.is_leaf():
            return frozenset([nd.id])
        left = walk(nd.get_left())
        right = walk(nd.get_right())
        leaves = left | right
        row = nd.id - n_leaves
        nodes[row] = {
            "row": row, "id": nd.id, "h": float(nd.dist),
            "size": len(leaves), "leaves": leaves,
        }
        return leaves
    walk(root)

    dmax = float(Z[-1, 2])
    h_floor = h_floor_frac * dmax

    # parent lookup
    for j in range(n_int):
        for col in (0, 1):
            child = int(Z[j, col])
            if child >= n_leaves:
                nodes[child - n_leaves]["parent_row"] = j
    nodes[n_int - 1]["parent_row"] = -1

    for nd in nodes:
        p = nd["parent_row"]
        if p < 0:
            nd["pi"] = np.nan
            nd["h_parent"] = np.nan
        else:
            nd["h_parent"] = float(Z[p, 2])
            nd["pi"] = np.log(nd["h_parent"]) - np.log(max(nd["h"], h_floor))
        nd["h_rel"] = nd["h"] / dmax if dmax > 0 else 0.0
    return nodes


def best_match(leafset: frozenset[int], ref_nodes: list[dict]) -> tuple[int, float]:
    """Return (row_of_best_match, J). Row = −1 if nothing overlaps."""
    best_j, best_row = 0.0, -1
    for rn in ref_nodes:
        j = jaccard_leafsets(leafset, rn["leaves"])
        if j > best_j:
            best_j, best_row = j, rn["row"]
    return best_row, best_j


# ──────────────────────────── module selection ────────────────────────────

def select_top_modules(nodes_rpost, nodes_tt, nodes_rpre, k: int,
                        n_leaves: int) -> list[dict]:
    """Rank rpost internal nodes by π · size · (1 − size/N) · (match_tt − match_rpre).

    ``size · (1 − size/N)`` is the scale-balanced weight (peaks at N/2), so
    leaf-pairs and root-adjacent nodes are suppressed. ``π`` is the
    log-branch-length (topological stability). Product ≈ "mid-size,
    long-branch" weight — the intuitive notion of a module.
    """
    ranked = []
    N = float(n_leaves)
    for η in nodes_rpost:
        if not np.isfinite(η["pi"]):
            continue
        s = float(η["size"])
        bal = s * (1.0 - s / N)
        if bal <= 0:
            continue
        r_tt, j_tt = best_match(η["leaves"], nodes_tt)
        r_rpre, j_rpre = best_match(η["leaves"], nodes_rpre)
        delta = j_tt - j_rpre
        score = η["pi"] * bal * delta
        if score <= 0:
            continue
        ranked.append({**η, "tt_row": r_tt, "j_tt": j_tt,
                       "rpre_row": r_rpre, "j_rpre": j_rpre,
                       "delta": delta, "bal": bal, "score": score})
    ranked.sort(key=lambda d: d["score"], reverse=True)
    return ranked[:k]


# ──────────────────────────── plotting helpers ────────────────────────────

def descendants(nodes: list[dict], root_row: int) -> set[int]:
    """All internal-node row-ids in the subtree rooted at ``root_row``."""
    if root_row < 0:
        return set()
    n_leaves = len(nodes) + 1
    out = set()
    stack = [root_row]
    # reverse lookup by leafset: we walk Z via children ids.
    Zrows = {nd["row"]: nd for nd in nodes}
    # We need the Z itself to descend. Re-derive from node leafsets: for
    # any internal row r, its children internal-rows are the two children
    # Z[r, 0] and Z[r, 1] if they are >= n_leaves.
    # But we don't have Z here; store it globally below.
    raise NotImplementedError  # replaced by inline version in phase_style


def phase_style(Z: np.ndarray, highlight_rows: dict[int, str]
                 ) -> dict:
    """Return per-row color map {internal_node_id : color} for link_color_func.

    ``highlight_rows``: dict row → color. All descendants of each listed row
    inherit its color.
    """
    n_leaves = Z.shape[0] + 1
    color_by_id = {}

    def paint(row, color):
        # BFS over children in Z.
        stack = [row]
        while stack:
            r = stack.pop()
            color_by_id[n_leaves + r] = color
            for col in (0, 1):
                child = int(Z[r, col])
                if child >= n_leaves:
                    stack.append(child - n_leaves)

    for row, color in highlight_rows.items():
        paint(row, color)
    return color_by_id


def highlight_subset_rows(Z: np.ndarray, target_leaves: frozenset[int]
                           ) -> list[int]:
    """Row-ids of internal nodes whose leafset is a subset of ``target_leaves``.

    Used for rest_pre where the selected module has NO coherent counterpart:
    we colour the maximal sub-subtrees that fit entirely inside the module's
    leaves — they will appear as multiple disjoint patches on rpre.
    """
    n_leaves = Z.shape[0] + 1
    # Build leafset per row bottom-up from Z.
    rows_leaves: list[frozenset[int]] = [frozenset()] * Z.shape[0]
    for r in range(Z.shape[0]):
        out = set()
        for col in (0, 1):
            child = int(Z[r, col])
            if child < n_leaves:
                out.add(child)
            else:
                out |= rows_leaves[child - n_leaves]
        rows_leaves[r] = frozenset(out)
    # maximal subsets: rows whose leaves ⊆ target_leaves AND whose parent's
    # leafset is NOT ⊆ target_leaves.
    parent_row = [-1] * Z.shape[0]
    for r in range(Z.shape[0]):
        for col in (0, 1):
            child = int(Z[r, col])
            if child >= n_leaves:
                parent_row[child - n_leaves] = r
    out_rows = []
    for r in range(Z.shape[0]):
        if rows_leaves[r].issubset(target_leaves) and len(rows_leaves[r]) >= 2:
            p = parent_row[r]
            if p < 0 or not rows_leaves[p].issubset(target_leaves):
                out_rows.append(r)
    return out_rows


def plot_phase(ax, Z: np.ndarray, modules: list[dict], which: str,
                title: str) -> None:
    """Plot one phase's dendrogram with per-module highlighting.

    ``which`` = 'rpost' | 'tt' | 'tl' | 'rpre'.
    """
    # Build highlight_rows per module.
    highlight_rows: dict[int, str] = {}
    for mi, m in enumerate(modules):
        color = MODULE_COLORS[mi % len(MODULE_COLORS)]
        if which == "rpost":
            highlight_rows[m["row"]] = color
        elif which == "tt":
            if m["j_tt"] > 0.2 and m["tt_row"] >= 0:
                highlight_rows[m["tt_row"]] = color
        elif which == "tl":
            # find best tl match for this module's leaves on the fly
            r, j = best_match(m["leaves"], modules[0]["_tl_nodes"]) if "_tl_nodes" in m else (-1, 0.0)
            if j > 0.2 and r >= 0:
                highlight_rows[r] = color
        elif which == "rpre":
            # In rpre, highlight maximal subset-subtrees of the module's leaves.
            for r in highlight_subset_rows(Z, m["leaves"]):
                highlight_rows[r] = color

    color_by_id = phase_style(Z, highlight_rows)

    def link_color(k):
        return color_by_id.get(k, BACKGROUND_GRAY)

    dendrogram(
        Z, ax=ax,
        no_labels=True,
        color_threshold=0,
        above_threshold_color=BACKGROUND_GRAY,
        link_color_func=link_color,
    )
    ax.set_yscale("log")
    dmax = float(Z[-1, 2])
    # Bound y: min merge / D_max to 1 (log scale).
    hmin = max(Z[Z[:, 2] > 0, 2].min() if (Z[:, 2] > 0).any() else 1e-3, 1e-3)
    ax.set_ylim(hmin * 0.8, dmax * 1.05)
    ax.set_title(title, fontsize=11)
    ax.set_xticks([])
    ax.tick_params(axis="y", labelsize=9)


# ──────────────────────────── main ────────────────────────────

def _load_Z(phase: str) -> np.ndarray:
    r = load_lrg_result(PATIENT, phase, BAND, "imcoh_abs", IMCOH_LRG_CACHE)
    assert r is not None, f"missing LRG for {PATIENT} {phase} {BAND}"
    return np.asarray(r.linkage_matrix)


def main() -> None:
    Zs = {ph: _load_Z(ph) for ph in PHASES}
    nodes = {ph: nodes_with_prominence(Zs[ph]) for ph in PHASES}

    n_leaves = Zs["rest_post"].shape[0] + 1
    modules = select_top_modules(
        nodes["rest_post"], nodes["task_test"], nodes["rest_pre"],
        K_MODULES, n_leaves,
    )
    # Attach tl-nodes list lazily via the first module (used by plot_phase).
    for m in modules:
        m["_tl_nodes"] = nodes["task_learn"]

    print(f"\nTop {K_MODULES} task-born-retained modules ({PATIENT}, {BAND}):")
    print(f"  row  size   π   h_rel(rpost)   match_tt   match_rpre   score")
    for mi, m in enumerate(modules):
        mcolor = MODULE_COLORS[mi]
        print(f"  {m['row']:3d}  {m['size']:3d}  {m['pi']:.3f}   "
              f"{m['h_rel']:.3f}          {m['j_tt']:.3f}      "
              f"{m['j_rpre']:.3f}       {m['score']:.3f}   ← {mcolor}")

    # Report the tt-matched node's h_rel for each module — this is the key
    # diagnostic on whether height-agreement is needed.
    print("\nHeight comparison (rpost vs best-matching tt subtree):")
    print(f"  mod  h_rel(rpost)   h_rel(tt)   Δh_rel   J(leaves)")
    for mi, m in enumerate(modules):
        if m["tt_row"] >= 0:
            h_tt = nodes["task_test"][m["tt_row"]]["h_rel"]
        else:
            h_tt = np.nan
        dh = abs(m["h_rel"] - h_tt)
        print(f"  {mi}    {m['h_rel']:.3f}          {h_tt:.3f}       "
              f"{dh:.3f}    {m['j_tt']:.3f}")

    # ─── plot ───
    fig, axes = plt.subplots(1, 4, figsize=(18.0, 5.0), dpi=160, sharey=True)
    phase_to_ax = dict(zip(PHASES, axes))

    plot_phase(phase_to_ax["rest_pre"], Zs["rest_pre"], modules, "rpre",
                PHASE_TEX["rest_pre"])
    plot_phase(phase_to_ax["task_learn"], Zs["task_learn"], modules, "tl",
                PHASE_TEX["task_learn"])
    plot_phase(phase_to_ax["task_test"], Zs["task_test"], modules, "tt",
                PHASE_TEX["task_test"])
    plot_phase(phase_to_ax["rest_post"], Zs["rest_post"], modules, "rpost",
                PHASE_TEX["rest_post"])

    for ax in axes:
        ax.set_xlabel("channels (dendrogram order)", fontsize=10)
    axes[0].set_ylabel(r"merge height $h$ (log)", fontsize=10)

    # Title + legend
    fig.suptitle(
        f"{PATIENT}, band = {BRAIN_BAND_TEX_DICT[BAND]} — top-{K_MODULES} "
        r"task-born-retained modules in rest$_\mathrm{post}$ "
        r"(scored by $\pi\cdot(J_{tt}-J_{rpre})$). "
        "Same colour = same leaf set across phases. "
        r"In rest$_\mathrm{pre}$ these leaves fragment.",
        fontsize=11, y=1.02,
    )
    legend_handles = []
    for mi, m in enumerate(modules):
        from matplotlib.patches import Patch
        legend_handles.append(Patch(
            facecolor=MODULE_COLORS[mi],
            label=(f"M{mi + 1}: size={m['size']}, "
                   f"$J_{{tt}}$={m['j_tt']:.2f}, "
                   f"$J_{{rpre}}$={m['j_rpre']:.2f}, "
                   f"$\\pi$={m['pi']:.2f}"),
        ))
    fig.legend(handles=legend_handles, loc="lower center",
                bbox_to_anchor=(0.5, -0.03), ncol=len(modules), fontsize=9,
                frameon=False)

    out_dir = REPORTS_ROOT / "imcoh_vi" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    base = out_dir / f"diag_pi_modules_{PATIENT}_{BAND}"
    for ext in ("pdf", "png"):
        fig.savefig(f"{base}.{ext}", bbox_inches="tight")
        print(f"saved {base}.{ext}")
    plt.close(fig)


if __name__ == "__main__":
    main()
