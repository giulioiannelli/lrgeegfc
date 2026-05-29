#!/usr/bin/env python3
"""q_bundle figure figS2 -- T_d swarm per band.

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
    print(f"=== q_figS2 for substrate: {SUBSTRATE} ===")
    distances_long, td, contrast = load_bundle_data(SRC)


    # ----------------------------------------------------------------------
    # Fig S2 — T_d swarm per band, INDEPENDENT y-axes
    # ----------------------------------------------------------------------
    td_swarm = pd.read_csv(SRC / "Td_per_patient_per_band.csv")

    fig, axes = plt.subplots(1, 6, figsize=(15.4, 3.6), squeeze=False)
    fig.subplots_adjust(top=0.88, bottom=0.18, left=0.05, right=0.985,
                        wspace=0.30)

    for c_idx, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[0, c_idx]
        sub = td_swarm[td_swarm.band == band]
        pat_arr = sub.patient.values
        is_p03 = pat_arr == "Pat_03"
        for x_pos, dist in enumerate(["S", "P", "F"]):
            vals = sub[dist].values
            rs = np.random.RandomState(11 + c_idx * 7 + x_pos * 13)
            jitter = rs.uniform(-0.13, 0.13, size=len(vals))
            # In-pool patients (not Pat_03)
            in_pool_mask = ~is_p03
            pool_vals = vals[in_pool_mask]
            pool_x = np.full(in_pool_mask.sum(), x_pos) + jitter[in_pool_mask]
            ax.scatter(
                pool_x, pool_vals,
                # T_d > 0 = trace (green); T_d < 0 = anti (red).
                c=["#2ca02c" if v > 0 else "#d62728" for v in pool_vals],
                edgecolor="k", linewidth=0.4, s=42, zorder=3, alpha=0.95,
            )
            # Pat_03 — orange triangle
            if is_p03.any():
                p03_vals = vals[is_p03]
                p03_x = np.full(is_p03.sum(), x_pos) + jitter[is_p03]
                ax.scatter(p03_x, p03_vals, marker="^", c="#ff7f0e",
                           edgecolor="k", linewidth=0.4, s=58, zorder=4)
            # Cohort median bar
            med = float(np.median(vals))
            ax.plot([x_pos - 0.27, x_pos + 0.27], [med, med],
                    color="k", lw=2.2, zorder=5)
            # n_trace / n_tot annotation, anchored just inside the top (T_d > 0).
            n_trace = int((vals > 0).sum())
            n_tot = len(vals)
            ax.annotate(f"{n_trace}/{n_tot}",
                        xy=(x_pos, 1.0), xycoords=("data", "axes fraction"),
                        xytext=(0, -2), textcoords="offset points",
                        ha="center", va="top",
                        fontsize=7.2, fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.18",
                                  facecolor="white", edgecolor="#bbb",
                                  alpha=0.92))
        ax.axhline(0, color="#999", lw=0.6, ls="--", zorder=1)
        ax.set_xlim(-0.55, 2.55)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels([r"$d_S$", r"$d_P$", r"$d_F$"], fontsize=9)
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10, pad=6)
        if c_idx == 0:
            ax.set_ylabel(r"$T_d$  (positive $=$ trace)",
                          fontsize=9, labelpad=2)
        ax.tick_params(labelsize=7)
        ax.grid(axis="y", color="#eee", lw=0.5, zorder=0)
        ax.set_axisbelow(True)

    handles_s2 = [
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="#2ca02c", markeredgecolor="k",
               markersize=8, label=r"patient with $T_d>0$ (trace)"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="#d62728", markeredgecolor="k",
               markersize=8, label=r"patient with $T_d<0$"),
        Line2D([0], [0], marker="^", color="none",
               markerfacecolor="#ff7f0e", markeredgecolor="k",
               markersize=9, label="Pat_03 (1024 Hz)"),
        Line2D([0], [0], color="k", lw=2.2, label="cohort median"),
        Line2D([0], [0], color="#999", lw=0.8, ls="--", label="zero line"),
    ]
    fig.legend(handles=handles_s2, loc="lower center",
               bbox_to_anchor=(0.5, 0.005), ncol=len(handles_s2),
               frameon=False, fontsize=8.0)

    fig.savefig(OUT / f"figS2_Td_swarm_all_distances{SUFFIX}.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT / f'figS2_Td_swarm_all_distances{SUFFIX}.pdf'}")


if __name__ == "__main__":
    main()
