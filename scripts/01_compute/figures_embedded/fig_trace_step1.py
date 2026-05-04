#!/usr/bin/env python3
"""Step 1 — sanity check for the task-trace claim, on ONE patient × ONE band.

For every internal node ``v`` of the rest_post dendrogram, compute

    trace_score(v) = J*(leaves(v), T_task_test) − J*(leaves(v), T_rest_pre)

Pick the ``v*`` with the highest score. Plot the four phase dendrograms
(rest_pre / task_learn / task_test / rest_post) side-by-side with the
leaves of ``v*`` coloured the same in all four. Goal: see whether
those leaves are scattered in rest_pre and clustered in the task and
rest_post panels.

If we can't find one clean visible example for one (patient, band),
the cohort claim is hopeless. This is the precondition for anything
else.

Default: Pat_06, alpha. CLI overrides ``--patient`` and ``--band``.

Output: ``data/reports/imcoh_mrl/figures/trace_step1_<patient>_<band>.pdf``
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.utils.metrics.tree import (
    tree_internal_nodes,
    jaccard_leafsets,
)
from lrg_eegfc.workflow.lrg import load_lrg_result


PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
HIGHLIGHT = "#d62728"   # task-trace leaves
NEUTRAL = "#888888"     # all other leaves


def load_Z(patient: str, phase: str, band: str) -> np.ndarray:
    r = load_lrg_result(patient, phase, band, fc_method="imcoh_abs",
                        cache_root=IMCOH_LRG_CACHE)
    if r is None:
        raise FileNotFoundError(f"No LRG cache for {patient} {phase} {band}")
    return np.asarray(r.linkage_matrix)


def best_match(S: frozenset[int], V_other: list[dict]) -> tuple[float, dict | None]:
    """Return (J*, best_node) for S vs the internal nodes of another tree."""
    best_j = 0.0
    best_u = None
    for u in V_other:
        j = jaccard_leafsets(S, u["leaves"])
        if j > best_j:
            best_j = j
            best_u = u
    return best_j, best_u


def find_best_trace_node(Z_post: np.ndarray, Z_test: np.ndarray,
                         Z_pre: np.ndarray, k_min: int = 5,
                         score_floor: float = 0.5,
                         ) -> tuple[dict, float, float, dict, dict, list[dict]]:
    """For every internal v ∈ T_post (size ≥ k_min), compute
    trace_score = J*(leaves(v), T_test) − J*(leaves(v), T_pre).

    Return the candidate v that:
      1. maximises ``size`` among all v with ``score >= score_floor``,
         i.e. the LARGEST visible trace module;
      2. or, if no v clears ``score_floor``, the global argmax score.

    Also return the top-10 candidates (by score) for diagnostics.
    """
    V_post = [v for v in tree_internal_nodes(Z_post)
              if k_min <= v["size"] <= Z_post.shape[0]]
    V_test = tree_internal_nodes(Z_test)
    V_pre  = tree_internal_nodes(Z_pre)
    rows = []
    for v in V_post:
        j_test, u_test = best_match(v["leaves"], V_test)
        j_pre,  u_pre  = best_match(v["leaves"], V_pre)
        rows.append({
            "v": v, "size": v["size"], "h_rel": v["h_rel"],
            "j_test": j_test, "j_pre": j_pre,
            "score": j_test - j_pre,
            "u_test": u_test, "u_pre": u_pre,
        })
    if not rows:
        raise RuntimeError("No eligible internal nodes in T_post.")

    above = [r for r in rows if r["score"] >= score_floor]
    if above:
        chosen = max(above, key=lambda r: (r["size"], r["score"]))
    else:
        chosen = max(rows, key=lambda r: r["score"])

    top = sorted(rows, key=lambda r: r["score"], reverse=True)[:10]
    return (chosen["v"], chosen["j_test"], chosen["j_pre"],
            chosen["u_test"], chosen["u_pre"], top)


def best_match_in_phase(S: frozenset[int], Z: np.ndarray
                        ) -> tuple[float, dict | None]:
    return best_match(S, tree_internal_nodes(Z))


def plot_one_dendrogram(ax: plt.Axes, Z: np.ndarray, S: frozenset[int],
                        match_node: dict | None, match_J: float,
                        title: str) -> None:
    """Plot one phase dendrogram with leaves of S highlighted and the
    best-matching subtree marked by a horizontal bar at its merge height.
    """
    n_leaves = Z.shape[0] + 1
    # Plot in log-scale; suppress automatic colouring.
    dres = dendrogram(
        Z,
        ax=ax,
        no_labels=True,
        color_threshold=0,
        above_threshold_color=NEUTRAL,
        leaf_rotation=0,
    )
    leaves_in_display_order = dres["leaves"]
    x_of_leaf = {leaf: 5 + 10 * i
                 for i, leaf in enumerate(leaves_in_display_order)}

    merge_heights = sorted(Z[:, 2].tolist())
    h_min = merge_heights[0]
    h_max = merge_heights[-1]
    y_lo = h_min * 0.8
    y_hi = h_max * 1.05

    # Vertical highlight bars at the x-position of each leaf in S, going
    # from y_lo up to a fixed fraction of the plot. Visible against grey
    # dendrogram lines.
    bar_top = h_min  # rises up to the lowest merge — distinct from dendro
    for leaf in S:
        x = x_of_leaf.get(leaf)
        if x is None:
            continue
        ax.vlines(x, y_lo, bar_top, color=HIGHLIGHT, lw=1.4,
                  alpha=0.95, zorder=5)

    # Mark the matched subtree only when the match is meaningful (J ≥ 0.5).
    # At low J the best-matching subtree is large by construction (a big
    # enclosing node is the best fit for scattered leaves), so the box
    # is visual noise. Below J=0.5 we let the scattered red leaf bars
    # speak for themselves.
    if match_node is not None and match_J >= 0.5:
        match_leaves = match_node["leaves"]
        xs = [x_of_leaf[l] for l in match_leaves if l in x_of_leaf]
        if xs:
            x_lo, x_hi = min(xs) - 4, max(xs) + 4
            h = match_node["h"]
            ax.add_patch(plt.Rectangle(
                (x_lo, y_lo), x_hi - x_lo, h - y_lo,
                facecolor=HIGHLIGHT, alpha=0.10, linewidth=0, zorder=2,
            ))
            ax.hlines(h, x_lo, x_hi, color=HIGHLIGHT, lw=2.4, zorder=6)

    # Annotation: J value + module size in the top-left corner of the panel.
    ax.text(0.02, 0.98, f"J = {match_J:.2f}",
            transform=ax.transAxes, fontsize=10, color=HIGHLIGHT,
            va="top", ha="left", fontweight="bold")

    if h_min > 0:
        ax.set_yscale("log")
    ax.set_ylim(y_lo, y_hi)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(f"contacts (n={n_leaves}, leaf order)", fontsize=9)
    ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", default="Pat_06")
    ap.add_argument("--band", default="alpha")
    ap.add_argument("--k-min", type=int, default=5)
    ap.add_argument("--score-floor", type=float, default=0.5)
    args = ap.parse_args()

    patient = args.patient
    band = args.band
    k_min = int(args.k_min)
    score_floor = float(args.score_floor)

    Zs = {phase: load_Z(patient, phase, band) for phase in PHASES}

    # 1. Find the LARGEST rest_post subtree with trace_score >= score_floor.
    v_star, J_test_v, J_pre_v, u_test, u_pre, top10 = find_best_trace_node(
        Z_post=Zs["rest_post"],
        Z_test=Zs["task_test"],
        Z_pre=Zs["rest_pre"],
        k_min=k_min,
        score_floor=score_floor,
    )
    S = v_star["leaves"]
    score = J_test_v - J_pre_v
    n_leaves = Zs["rest_pre"].shape[0] + 1

    print(f"[step1] {patient} / {band}")
    print(f"[step1] picked rest_post subtree:  size={v_star['size']}/"
          f"{n_leaves}  h_rel={v_star['h_rel']:.3f}")
    print(f"[step1] J*(post-subtree, T_task_test) = {J_test_v:.3f}")
    print(f"[step1] J*(post-subtree, T_rest_pre)  = {J_pre_v:.3f}")
    print(f"[step1] trace score = {score:+.3f}")
    print(f"[step1] leaves: {sorted(S)}")
    print()
    print("[step1] top-10 candidates by raw score (post → task_test − rest_pre):")
    print(f"[step1]   {'size':>4s}  {'h_rel':>5s}  {'j_test':>6s}  "
          f"{'j_pre':>6s}  {'score':>6s}")
    for r in top10:
        print(f"[step1]   {r['size']:4d}  {r['h_rel']:5.3f}  "
              f"{r['j_test']:6.3f}  {r['j_pre']:6.3f}  {r['score']:+.3f}")

    # Also compute J* of S against task_learn (same logic, for figure label).
    J_learn_v, u_learn = best_match_in_phase(S, Zs["task_learn"])

    # 2. Plot the 4 dendrograms with S highlighted.
    fig, axes = plt.subplots(1, 4, figsize=(18.0, 4.6), dpi=150,
                             sharey=False)
    titles = [
        f"rest_pre  (J = {J_pre_v:.2f})",
        f"task_learn  (J = {J_learn_v:.2f})",
        f"task_test  (J = {J_test_v:.2f})",
        f"rest_post  (J = 1.00, native)",
    ]
    matched_nodes = [u_pre, u_learn, u_test, v_star]
    matched_Js    = [J_pre_v, J_learn_v, J_test_v, 1.0]
    for ax, phase, title, match_node, match_J in zip(
        axes, PHASES, titles, matched_nodes, matched_Js
    ):
        plot_one_dendrogram(ax, Zs[phase], S, match_node, match_J, title)

    axes[0].set_ylabel(r"merge height $h$ (log)", fontsize=10)

    # No suptitle (per project rule). Sidecar prose carries the context.
    fig.tight_layout()
    out_dir = REPORTS_ROOT / "imcoh_mrl" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"trace_step1_{patient}_{band}.pdf"
    # Rasterise the dendrogram line collections inside the PDF for size.
    for ax in axes:
        for coll in ax.collections:
            coll.set_rasterized(True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"[step1] saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
