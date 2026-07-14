#!/usr/bin/env python3
r"""talk_fig_measures_ladder -- detect != discriminate, the full measure ladder (S4/S6 talk).

The talk-clean, dark-slide version of new_results_sec1/fig_controls_ladder: a dot-matrix of
every trace read-out (rows) x every band (cols), each cell gated on the SAME mst@0.20 backbone
against the SAME matched-strength null with the SAME representativeness gate (cohort p<0.05 AND
leave-one-out AND >=3/10 patients above their own null; _common.representativeness_gate).

FILL = representative clear (band colour) · RING = cohort-sig but not representative ·
open grey = null. Dot AREA = #patients above their own null. Reading (T_test / trace):
  raw FC edges -> beta only (the carrier; raw sees it, so NOT "blind to beta")
  node strength / clustering / resistance -> nothing (the simple descriptors are blind)
  graph geodesic -> low_gamma only (a NON-carrier band -- wrong answer)
  cophenetic (ours) -> the ONLY read-out representative in BOTH carriers (alpha, beta)
  [Grassmann row appended once the whole-graph MS numbers land -- flagged as WHOLE-GRAPH]

Dark-slide mode: all ink white, transparent background, colored data untouched.
Reads : data/sparsified_arc/controls_ladder_apples/per_cell.csv  (functional == "T_test")
Writes: data/outputs/figures/talk/fig_measures_ladder.{png,pdf}   (transparent, Canva-ready)
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[1] /
                       "01_compute" / "figures_embedded" / "new_results_sec1"))
import _common as C  # noqa: E402

C.use_lrg_style()

# dark-slide mode: all ink white, colored data untouched
INK = "white"
plt.rcParams.update({
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
    "axes.titlecolor": INK, "xtick.color": INK, "ytick.color": INK,
})
SEP = "0.6"          # separator + null-ring grey, lightened for dark
F = "T_test"
OUT = C.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_measures_ladder"


def _cell_dark(ax, j, y, band, res):
    """plot_gate_cell, but null/degenerate greys lightened for a dark background."""
    lvl, n = res["level"], res["n_above"]
    col = C.band_color(band)
    if lvl == "degenerate":
        ax.scatter(j, y, marker="_", s=190, color=SEP, lw=2.2, zorder=3)
    elif lvl == "rep":
        ax.scatter(j, y, s=C.gate_area(n), facecolor=col, edgecolor="white", lw=0.8, zorder=4)
    elif lvl == "fragile":
        ax.scatter(j, y, s=C.gate_area(n), facecolor="none", edgecolor=col, lw=2.1, zorder=3)
    else:
        ax.scatter(j, y, s=C.gate_area(n), facecolor="none", edgecolor=SEP, lw=1.2, zorder=2)


def main():
    df = pd.read_csv(C.LADDER / "per_cell.csv")
    descs, bands = C.LADDER_DESCS, list(C.BANDS)
    n_desc = len(descs)

    fig, ax = plt.subplots(figsize=(9.6, 6.0))
    # faint separator setting off the nested-hierarchy (cophenetic) row at the bottom
    ax.axhline(len(C.LADDER_HIER) - 0.5, color=SEP, lw=1.0, ls=(0, (4, 3)), zorder=0)

    for i, (desc, _) in enumerate(descs):
        y = n_desc - 1 - i                       # first descriptor on top
        for j, band in enumerate(bands):
            _cell_dark(ax, j, y, band, C.representativeness_gate(df, desc, band, F))

    ax.set_xlim(-0.7, len(bands) - 0.3)
    ax.set_ylim(-0.7, n_desc - 0.3)
    ax.set_xticks(range(len(bands)))
    ax.set_xticklabels([C.BTeX[b] for b in bands], fontsize=17)
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    ax.set_yticks([n_desc - 1 - i for i in range(n_desc)])
    ax.set_yticklabels([lab for _, lab in descs], fontsize=13)
    # bold the "ours" cophenetic row label
    for t, (desc, _) in zip(ax.get_yticklabels(), descs[::-1]):
        if desc in C.LADDER_HIER:
            t.set_fontweight("bold")
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)

    # legend: level (fill/ring/null) + size = #patients above null
    size_h = [plt.scatter([], [], s=C.gate_area(k), facecolor=SEP, edgecolor="white",
                          lw=0.6, label=f"{k}/10") for k in (2, 5, 8)]
    lvl_h = [
        Line2D([], [], marker="o", ls="none", ms=11, mfc=SEP, mec="white",
               label="representative"),
        Line2D([], [], marker="o", ls="none", ms=11, mfc="none", mec=SEP, mew=2.1,
               label="sig. not representative"),
        Line2D([], [], marker="o", ls="none", ms=11, mfc="none", mec=SEP,
               label="null"),
    ]
    leg1 = fig.legend(handles=size_h, title="patients above own null",
                      loc="lower left", bbox_to_anchor=(0.10, -0.02), ncol=3,
                      frameon=False, fontsize=10, title_fontsize=10, handletextpad=0.2,
                      columnspacing=1.1)
    fig.add_artist(leg1)
    fig.legend(handles=lvl_h, loc="lower right", bbox_to_anchor=(0.99, -0.02),
               ncol=3, frameon=False, fontsize=10, handletextpad=0.3, columnspacing=1.1)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}.{ext}", transparent=True, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"wrote {OUT}.png + .pdf")
    for desc, lab in descs:
        cells = {b: C.representativeness_gate(df, desc, b, F)["level"] for b in bands}
        reps = [b for b, lv in cells.items() if lv == "rep"]
        print(f"  {lab:22s} representative: {reps or '—'}")


if __name__ == "__main__":
    main()
