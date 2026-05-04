#!/usr/bin/env python3
"""Step 4 — multi-patient dendrogram reprojection for one band.

For the chosen band (default ``alpha``), pick the top-K patients by
trace score (from ``trace_sweep.csv``). For each patient, re-find the
best rest_post subtree (same logic as Step 1) and plot the four phase
dendrograms in a row, with the trace leaves coloured.

K rows × 4 columns. PDF only, dendrogram line-collections rasterised.

Reads:  ``data/reports/imcoh_mrl/trace_sweep.csv``
Writes: ``data/reports/imcoh_mrl/figures/trace_multipatient_<band>.pdf``
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
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
HIGHLIGHT = "#d62728"
NEUTRAL = "#888888"


def load_Z(patient: str, phase: str, band: str) -> np.ndarray:
    r = load_lrg_result(patient, phase, band, fc_method="imcoh_abs",
                        cache_root=IMCOH_LRG_CACHE)
    return np.asarray(r.linkage_matrix)


def best_match(S: frozenset[int], V_other: list[dict]
               ) -> tuple[float, dict | None]:
    best_j = 0.0
    best_u = None
    for u in V_other:
        j = jaccard_leafsets(S, u["leaves"])
        if j > best_j:
            best_j = j
            best_u = u
    return best_j, best_u


def plot_one_dendrogram(ax: plt.Axes, Z: np.ndarray, S: frozenset[int],
                        match_node: dict | None, match_J: float,
                        title: str, ylabel: str | None = None) -> None:
    n_leaves = Z.shape[0] + 1
    dres = dendrogram(
        Z, ax=ax, no_labels=True, color_threshold=0,
        above_threshold_color=NEUTRAL, leaf_rotation=0,
    )
    leaves_in_display = dres["leaves"]
    x_of_leaf = {leaf: 5 + 10 * i for i, leaf in enumerate(leaves_in_display)}
    merge_heights = sorted(Z[:, 2].tolist())
    h_min = merge_heights[0]
    h_max = merge_heights[-1]
    y_lo = h_min * 0.8
    y_hi = h_max * 1.05
    bar_top = h_min  # rise to first merge

    for leaf in S:
        x = x_of_leaf.get(leaf)
        if x is None:
            continue
        ax.vlines(x, y_lo, bar_top, color=HIGHLIGHT, lw=1.4,
                  alpha=0.95, zorder=5)

    # Only draw the matched-subtree box when the match is meaningful
    # (J ≥ 0.5). At low J the matched subtree is large by construction
    # (a big enclosing node is the best fit for scattered leaves), so
    # the box would be a visual noise rather than a marker.
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

    ax.text(0.02, 0.97, f"J = {match_J:.2f}",
            transform=ax.transAxes, fontsize=9, color=HIGHLIGHT,
            va="top", ha="left", fontweight="bold")

    if h_min > 0:
        ax.set_yscale("log")
    ax.set_ylim(y_lo, y_hi)
    ax.set_title(title, fontsize=10)
    ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9)
    else:
        ax.set_yticklabels([])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default="alpha")
    ap.add_argument("--top-k", type=int, default=4)
    args = ap.parse_args()
    band = args.band
    top_k = int(args.top_k)

    df = pd.read_csv(REPORTS_ROOT / "imcoh_mrl" / "trace_sweep.csv")
    sub = df[df["band"] == band].sort_values("score", ascending=False).head(top_k)
    if sub.empty:
        raise SystemExit(f"No sweep rows for band={band}")

    chosen = []
    for _, r in sub.iterrows():
        pat = str(r["patient"])
        S = frozenset(int(x) for x in json.loads(r["leaves"]))
        chosen.append({
            "patient": pat, "leaves": S, "size": int(r["size"]),
            "h_rel": float(r["h_rel"]), "score": float(r["score"]),
        })
        print(f"[fig4] {pat}  size={int(r['size'])}  "
              f"h_rel={r['h_rel']:.3f}  score={r['score']:+.3f}")

    n_rows = len(chosen)
    fig, axes = plt.subplots(n_rows, 4, figsize=(18.0, 3.4 * n_rows),
                             dpi=150, squeeze=False)

    for i, c in enumerate(chosen):
        pat = c["patient"]
        S = c["leaves"]
        Zs = {phase: load_Z(pat, phase, band) for phase in PHASES}

        # Compute J* and matched nodes per phase.
        Js = {}
        Us = {}
        for phase in PHASES:
            j, u = best_match(S, tree_internal_nodes(Zs[phase]))
            Js[phase] = j
            Us[phase] = u

        for col, phase in enumerate(PHASES):
            ax = axes[i, col]
            ylabel = (f"{pat}\nsize={c['size']}, h_rel={c['h_rel']:.2f}\n"
                      r"merge $h$ (log)") if col == 0 else None
            ttl = phase.replace("_", "\\_") if False else phase
            ax.set_title(ttl, fontsize=10)
            plot_one_dendrogram(ax, Zs[phase], S,
                                Us[phase], Js[phase],
                                title=ttl, ylabel=ylabel)

    fig.tight_layout()
    out = REPORTS_ROOT / "imcoh_mrl" / "figures" / f"trace_multipatient_{band}.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    for ax in axes.flatten():
        for coll in ax.collections:
            coll.set_rasterized(True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"[fig4] saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
