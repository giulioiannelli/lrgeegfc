#!/usr/bin/env python3
r"""talk_fig_epi_relational_marker -- dark-slide wrap of fig:epi_a (slide 17).

The cohort seizure-zone marker: strength-residual seed-affinity AUC per band, one dot per
patient (open = single fine scale, filled = multiscale; green = beats own matched-strength
fake-SOZ null), cohort-median bar per band. Reuses the paper figure's draw() verbatim so the
numbers can never drift; only the ink is flipped white for a dark Canva slide.

Reads : data/sparsified_arc/epi_arc_mst020/per_cell.csv  (via fig_epi_a_relational_marker.draw)
Writes: data/outputs/figures/talk/fig_epi_relational_marker.{png,pdf}  (transparent, dark)
"""
from __future__ import annotations
import importlib.util
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FE = Path(__file__).resolve().parents[1] / "01_compute" / "figures_embedded"
spec = importlib.util.spec_from_file_location("fig_epi_a", FE / "fig_epi_a_relational_marker.py")
epi_a = importlib.util.module_from_spec(spec)
sys.modules["fig_epi_a"] = epi_a
spec.loader.exec_module(epi_a)

INK = "white"
plt.rcParams.update({
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
    "axes.titlecolor": INK, "xtick.color": INK, "ytick.color": INK,
})
OUT = epi_a.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_epi_relational_marker"


def main():
    fig = plt.figure(figsize=(9.2, 5.4))
    epi_a.draw(fig)
    # the single-scale open dots use facecolor="white" (invisible on dark) -> re-tint to slate
    for ax in fig.axes:
        for pc in ax.collections:
            fc = pc.get_facecolor()
            if len(fc) and tuple(round(v, 2) for v in fc[0][:3]) == (1.0, 1.0, 1.0):
                pc.set_facecolor("#2a2f36")     # dark slate so the open ring reads on dark
    # match the legend swatches to the dark plot: single-scale -> slate, cohort-median -> visible
    for lg in fig.legends:
        for h in lg.legend_handles:
            fcw = h.get_markerfacecolor() if hasattr(h, "get_markerfacecolor") else None
            if fcw in ("white", "#ffffff", (1.0, 1.0, 1.0, 1.0)):
                h.set_markerfacecolor("#2a2f36")
            if fcw == "#444":
                h.set_markerfacecolor("#c9ced6")
                h.set_color("#c9ced6")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}.{ext}", transparent=True, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"wrote {OUT}.png + .pdf")


if __name__ == "__main__":
    main()
