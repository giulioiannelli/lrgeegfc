#!/usr/bin/env python3
"""Audit 38b — Visual spotlight of task-anchored TRACE subtrees.

For each (patient, band) cell, pick the top-K TRACE-classified anchor subtrees
from audit_38's anchor_classification.csv (ranked by a_trace * size, so the
most confident large anchors come first). For each top anchor S, plot the
THREE main phase dendrograms (rest_pre, task_test, rest_post) side-by-side
with the leaves of S highlighted in distinct colors. Visual demonstration:
in TT and RPost the leaves cluster in a tight subtree (TRACE); in RPre they
are scattered (the corner-distance soft classification's "low COH_rpre"
captures this).

This is the figure the manuscript needs to make the trace claim concrete:
"these N leaves go from scattered in rest_pre to tight modules in task_test
and rest_post."

Output
------
``data/audit/per_patient_hierarchy_task_anchored/spotlight/{Pat_XX}/{band}_top{K}.pdf``
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import (  # noqa: E402
    BAND_ORDER, PHASE_TEX, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE,
    BACKGROUND_GRAY, _load_Z,
)

ANCHOR_CSV = ROOT / "data" / "audit" / "per_patient_hierarchy_task_anchored" / "anchor_classification.csv"
OUT_DIR = ROOT / "data" / "audit" / "per_patient_hierarchy_task_anchored" / "spotlight"
SPOTLIGHT_PHASES = ("rest_pre", "task_test", "rest_post")
TOP_K = 3
SPOTLIGHT_COLORS = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd"]


def parse_leaves(s: str) -> set[int]:
    return {int(x) for x in str(s).split("|") if x}


def link_color_func(Z: np.ndarray, leafset: set[int],
                     color: str) -> dict[int, str]:
    """Color internal links red iff every leaf below is in leafset."""
    n_leaves = Z.shape[0] + 1
    rows_leaves: list[set[int]] = [set() for _ in range(Z.shape[0])]
    for r in range(Z.shape[0]):
        out = set()
        for col in (0, 1):
            child = int(Z[r, col])
            if child < n_leaves:
                out.add(child)
            else:
                out |= rows_leaves[child - n_leaves]
        rows_leaves[r] = out
    out_map = {}
    for r in range(Z.shape[0]):
        nid = n_leaves + r
        leaves = rows_leaves[r]
        if leaves <= leafset:           # subset → all in leafset
            out_map[nid] = color
        else:
            out_map[nid] = BACKGROUND_GRAY
    return out_map


def plot_three_phase(Zs, leafset, color, title_prefix, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.4))
    for ax, ph in zip(axes, SPOTLIGHT_PHASES):
        Z = Zs[ph]
        if Z is None:
            ax.text(0.5, 0.5, "missing", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_axis_off()
            continue
        link_colors = link_color_func(Z, leafset, color)
        n_leaves = Z.shape[0] + 1
        leaf_colors = [color if i in leafset else BACKGROUND_GRAY for i in range(n_leaves)]
        # dendrogram doesn't take per-leaf colors directly; build via xticks after
        dendrogram(
            Z, ax=ax, no_labels=True, color_threshold=0.0,
            link_color_func=lambda nid: link_colors.get(nid, BACKGROUND_GRAY),
            above_threshold_color=BACKGROUND_GRAY,
        )
        # Color leaf x-tick markers
        leaf_order = []
        for t in ax.get_xticklabels():
            try:
                leaf_order.append(int(t.get_text()))
            except ValueError:
                pass
        # Add a colored bar at the bottom showing which leaves are in the spotlight
        ymin, ymax = ax.get_ylim()
        bar_y = ymin - 0.05 * (ymax - ymin)
        for i, l in enumerate(leaf_order):
            if l in leafset:
                ax.scatter(5 + 10 * i, bar_y, marker="s", s=14,
                           color=color, clip_on=False, zorder=10)
        ax.set_title(f"{PHASE_TEX.get(ph, ph)}", fontsize=11)
        ax.set_xticks([])
    fig.suptitle(title_prefix, fontsize=10, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(ANCHOR_CSV)
    df["dominant"] = df[["a_trace", "a_persist", "a_reset", "a_rearrange"]].idxmax(axis=1).str.replace("a_", "")
    df_trace = df[df["dominant"] == "trace"].copy()
    df_trace["score"] = df_trace["a_trace"] * df_trace["size"].astype(float)

    n_emitted = 0
    for pat in PATIENTS_4PHASE:
        for band in BAND_ORDER:
            sub = df_trace[(df_trace["patient"] == pat) & (df_trace["band"] == band)]
            if sub.empty:
                continue
            Zs = {ph: _load_Z(pat, ph, band) for ph in SPOTLIGHT_PHASES}
            top = sub.nlargest(TOP_K, "score")
            for rank, (_, r) in enumerate(top.iterrows(), start=1):
                leafset = parse_leaves(r["leaves"])
                color = SPOTLIGHT_COLORS[rank - 1]
                title = (rf"{pat}  $\mathbf{{{BRAIN_BAND_TEX_DICT[band].strip('$')}}}$"
                         rf"  --  TRACE anchor #{rank} | size = {len(leafset)} | "
                         rf"$a_{{trace}} = {r['a_trace']:.2f}$ | "
                         rf"$\mathrm{{COH}}_{{\mathrm{{RPre}}}} = {r['COH_rpre']:.2f}$, "
                         rf"$\mathrm{{COH}}_{{\mathrm{{RPost}}}} = {r['COH_rpost']:.2f}$")
                out_path = OUT_DIR / pat / f"{band}_anchor{rank}.pdf"
                plot_three_phase(Zs, leafset, color, title, out_path)
                n_emitted += 1
    print(f"[audit_38b] emitted {n_emitted} spotlight figures at {OUT_DIR}")


if __name__ == "__main__":
    main()
