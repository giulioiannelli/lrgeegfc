#!/usr/bin/env python3
r"""talk_fig_epi_two_populations -- dark-slide wrap of fig:epi_c (slide 17).

The cohort quantitative payoff: a calibrated leave-one-patient-out detector. Panel a = LOPO
detector AUC per patient -- 8 patients are the co-diffusing community (propagator strong, up
to 0.98), 2 are right-hemisphere hubs where the community read fails (Pat_15 0.57, Pat_10 0.49)
and node strength rescues them. Panel b = affinity -> per-patient switch (Pat_15 0.21 -> 0.93);
cohort 0.72 -> 0.75, 8/10 -> 9/10, Pat_10 unrecoverable. Reuses the paper draw() verbatim; only
the ink is flipped white for a dark Canva slide.

Reads : data/audit/epi_propagator_detector/*.csv  (via fig_epi_c_two_populations.draw)
Writes: data/outputs/figures/talk/fig_epi_two_populations.{png,pdf}  (transparent, dark)
"""
from __future__ import annotations
import importlib.util
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FE = Path(__file__).resolve().parents[1] / "01_compute" / "figures_embedded"
spec = importlib.util.spec_from_file_location("fig_epi_c", FE / "fig_epi_c_two_populations.py")
epi_c = importlib.util.module_from_spec(spec)
sys.modules["fig_epi_c"] = epi_c
spec.loader.exec_module(epi_c)

INK = "white"
plt.rcParams.update({
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
    "axes.titlecolor": INK, "xtick.color": INK, "ytick.color": INK,
})
OUT = epi_c.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_epi_two_populations"


def main():
    fig = plt.figure(figsize=(12.4, 5.6))
    epi_c.draw(fig)
    # any light annotation box (facecolor white / near-white) -> dark slate so text reads white on it
    for ax in fig.axes:
        for txt in ax.texts:
            bb = txt.get_bbox_patch()
            if bb is not None:                        # boxed annotation: dark slate box + white ink
                bb.set_facecolor("#2a2f36")
                bb.set_edgecolor("#5a616b")
                txt.set_color("white")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}.{ext}", transparent=True, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"wrote {OUT}.png + .pdf")


if __name__ == "__main__":
    main()
