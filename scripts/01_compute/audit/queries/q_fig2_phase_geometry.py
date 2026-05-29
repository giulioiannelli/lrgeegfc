#!/usr/bin/env python3
"""q_bundle figure fig2 -- phase geometry (paired strip + IQR box).

Split out of q_bundle_figures.py on 2026-05-29 (Phase 4-B split 3/7).
"""
from __future__ import annotations

import sys
from pathlib import Path

# Folder-local imports
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd  # noqa: F401
from matplotlib.patches import Ellipse, Patch, Rectangle, Arc  # noqa: F401
from matplotlib.lines import Line2D  # noqa: F401
from matplotlib.colors import Normalize  # noqa: F401
from scipy.stats import pearsonr, spearmanr  # noqa: F401

from _q_bundle_shared import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, NODE_SIZE_DEFAULT, N_COHORT,
    PAIR_COLORS, PAIR_LABELS, PAIR_ORDER, PATIENTS_4PHASE,
    PHASE_SHORT, ROOT,
    compute_network_layout, draw_network_edges, _probe_color_map,
    extract_probe_labels, imshow_colorbar_caxdivider, load_bundle_data,
    load_channel_labels, load_fc_matrix, lookup_pair, resolve_substrate,
)


def main() -> None:
    SUBSTRATE, SUFFIX, SRC, OUT = resolve_substrate()
    FC_METHOD_FIG5 = SUBSTRATE
    print(f"=== q_fig2 for substrate: {SUBSTRATE} ===")
    distances_long, td, contrast = load_bundle_data(SRC)
    # In the original q_bundle_figures.py, ``sub_S`` was defined once at
    # module top after the data load (line 133) and reused across figs.
    sub_S = distances_long[distances_long.distance == "S"]

    # ----------------------------------------------------------------------
    # Fig 2 — Phase geometry (paired strip + IQR box, sorted phase pairs)
    # ----------------------------------------------------------------------
    fig, axes = plt.subplots(1, 6, figsize=(15.6, 3.6), squeeze=False)

    for c, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[0, c]
        pair_data = {pr: lookup_pair(sub_S, band, pr[0], pr[1], "S")
                     for pr in PAIR_ORDER}
        common = sorted(set.intersection(*[set(d.index) for d in pair_data.values()]))
        # Per-band sort: small → large by cohort median
        band_order = sorted(
            PAIR_ORDER,
            key=lambda pr: float(np.median(pair_data[pr].loc[common].values)),
        )
        xs = np.arange(len(band_order))

        # Per-patient connecting lines (in band-sorted order)
        for p in common:
            ys = [pair_data[pr][p] for pr in band_order]
            ax.plot(xs, ys, color="#bbb", lw=0.45, alpha=0.55, zorder=1)

        # Per-pair: IQR box + median bar + strip dots
        for i, pr in enumerate(band_order):
            vals = pair_data[pr].loc[common].values
            col = PAIR_COLORS[pr]
            q1, med, q3 = np.percentile(vals, [25, 50, 75])

            # IQR box (translucent fill)
            ax.add_patch(Rectangle(
                (i - 0.20, q1), 0.40, q3 - q1,
                facecolor=col, alpha=0.20, edgecolor=col,
                linewidth=0.7, zorder=2,
            ))
            # Median line
            ax.plot([i - 0.24, i + 0.24], [med, med],
                    color=col, lw=2.2, zorder=4)

            # Strip dots (deterministic jitter)
            rs = np.random.RandomState(42 + i + c * 13)
            jitter = rs.uniform(-0.11, 0.11, size=len(vals))
            is_p03 = np.asarray(common) == "Pat_03"
            ax.scatter(np.full_like(vals, i)[~is_p03] + jitter[~is_p03],
                       vals[~is_p03], c="white",
                       edgecolor=col, linewidth=1.0, s=26,
                       alpha=0.95, zorder=3)
            if is_p03.any():
                ax.scatter(np.full_like(vals, i)[is_p03] + jitter[is_p03],
                           vals[is_p03], marker="^", c="white",
                           edgecolor=col, linewidth=1.0, s=34,
                           alpha=0.95, zorder=3)

        ax.set_xticks(xs)
        ax.set_xticklabels([PAIR_LABELS[pr] for pr in band_order],
                           fontsize=7.2, rotation=22, ha="right")
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
        ax.tick_params(axis="y", labelsize=7)
        if c == 0:
            ax.set_ylabel(rf"$d_S$  (n={N_COHORT} patients)", fontsize=9)
        ax.grid(axis="y", color="#eee", lw=0.5, zorder=0)
        ax.set_axisbelow(True)

    handles_f2 = [Patch(facecolor=PAIR_COLORS[pr], alpha=0.55,
                        edgecolor=PAIR_COLORS[pr],
                        label=PAIR_LABELS[pr])
                  for pr in PAIR_ORDER]
    handles_f2 += [
        Line2D([0], [0], color="#bbb", lw=0.7, label="patient (paired)"),
        Line2D([0], [0], color="k", lw=2.0, label="cohort median"),
    ]
    fig.legend(handles=handles_f2, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=len(handles_f2),
               frameon=False, fontsize=8.2)
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    fig.savefig(OUT / f"fig2_phase_geometry{SUFFIX}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT / f'fig2_phase_geometry{SUFFIX}.pdf'}")


if __name__ == "__main__":
    main()
