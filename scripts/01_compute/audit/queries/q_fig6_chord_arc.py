#!/usr/bin/env python3
"""q_bundle figure fig6 -- 4-phase geometry as chord/arc diagrams.

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
from matplotlib.colors import Colormap, LinearSegmentedColormap, Normalize  # noqa: F401
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
    print(f"=== q_fig6 for substrate: {SUBSTRATE} ===")
    distances_long, td, contrast = load_bundle_data(SRC)



    # ----------------------------------------------------------------------
    # Fig 6 — 4-phase geometry as chord/arc diagrams (no redundant matrix)
    # ----------------------------------------------------------------------
    # A 4-phase distance "matrix" has 4×4 = 16 cells but only 6 unique pair
    # distances (4 choose 2). The diagonal is trivial (=0) and the lower
    # triangle is symmetric to the upper. Replacing the matrix with a
    # chord diagram per (band, distance) cell:
    #   - 4 phase markers on a baseline preserve the temporal order
    #     RPre → TL → TT → RPost
    #   - 6 arcs above connect each pair
    #   - arc COLOUR encodes the cohort-median distance value
    #   - arc HEIGHT encodes the temporal gap (adjacent / skip-1 / skip-2)
    geom = pd.read_csv(SRC / "cohort_geometry_4phase_summary.csv")

    PHASES_ORDER = ["rest_pre", "task_learn", "task_test", "rest_post"]
    PHASE_X = {ph: i for i, ph in enumerate(PHASES_ORDER)}
    PHASE_LABEL_SHORT = {"rest_pre": "RPre", "task_learn": "TL",
                         "task_test": "TT", "rest_post": "RPost"}
    DISTANCE_KEYS = ["S", "P", "F"]
    DISTANCE_LABEL = {
        "S": r"$d_S$  Spearman",
        "P": r"$d_P$  Pearson",
        "F": r"$d_F$  Frobenius",
    }

    fig = plt.figure(figsize=(15.4, 5.8))
    gs = fig.add_gridspec(
        3, 6,
        width_ratios=[1, 1, 1, 1, 1, 1],
        height_ratios=[1, 1, 1],
        wspace=0.32, hspace=0.32,
        top=0.91, bottom=0.10, left=0.05, right=0.985,
    )
    REF_LINES = [0.0, 0.25, 0.5, 0.75, 1.0]

    # Custom 2-colour colormap: vivid green (close) → neutral gray midpoint
    # → vivid red (far). NO yellow waypoint, so the perceived semantic is
    # binary "good/bad" with a clean neutral middle.
    cmap_chord = LinearSegmentedColormap.from_list(
        "close_gray_far",
        [(0.00, "#1a9850"),   # vivid green = close
         (0.50, "#9e9e9e"),   # neutral gray midpoint (no yellow)
         (1.00, "#d73027")],  # vivid red   = far
        N=256,
    )


    # scipy.optimize.minimize was imported at module-top in the original
    # q_bundle_figures.py; recreate the alias inside main() for the split.
    from scipy.optimize import minimize as _spo_minimize

    def _solve_1d_spring_layout(panel_distances, n_phases=4):
        """Weighted 1D MDS: place `n_phases` phase points on a line so the
        pairwise distances |p_i - p_j| match `panel_distances` as closely
        as possible, with **spring stiffness w_ij = 1 / d_ij**² (closer
        pairs pull MUCH harder than far ones — amplifies the difference
        between near and distant pairs). Monotonic order p_0 < p_1 < ...
        < p_{n-1} is enforced via a positive-delta parameterisation.

        panel_distances: list of (i, j, d_ij) with i < j.
        Returns: array of n_phases positions, p_0 = 0.
        """
        def _stress(deltas):
            steps = np.abs(deltas)
            p = np.concatenate([[0.0], np.cumsum(steps)])
            s = 0.0
            for i, j, d in panel_distances:
                w = 1.0 / max(d, 1e-3) ** 2     # 1/d^2 stiffness
                s += w * (abs(p[i] - p[j]) - d) ** 2
            return s

        mean_d = np.mean([d for _, _, d in panel_distances])
        x0 = np.full(n_phases - 1, mean_d)
        result = _spo_minimize(
            _stress, x0, method="Nelder-Mead",
            options={"xatol": 1e-7, "fatol": 1e-9, "maxiter": 8000},
        )
        deltas = np.abs(result.x)
        return np.concatenate([[0.0], np.cumsum(deltas)])


    from mpl_toolkits.axes_grid1 import make_axes_locatable as _make_axes_locatable

    for r_idx, dist in enumerate(DISTANCE_KEYS):
        sub = geom[geom.distance == dist]

        for c_idx, band in enumerate(BRAIN_BANDS_NAMES):
            ax = fig.add_subplot(gs[r_idx, c_idx])
            panel = sub[sub.band == band]

            panel_distances = []
            for _, row in panel.iterrows():
                i = PHASE_X[row.phase_A]
                j = PHASE_X[row.phase_B]
                if i > j:
                    i, j = j, i
                panel_distances.append((i, j, float(row["median"])))

            # Per-PANEL colour normalisation (relative within panel —
            # distances are not directly comparable across cells).
            d_values = np.array([d for _, _, d in panel_distances])
            vmin = float(d_values.min())
            vmax = float(d_values.max())
            if vmax <= vmin:
                vmax = vmin + 1e-6
            norm_d = Normalize(vmin=vmin, vmax=vmax)

            # 1-D spring layout (stiffness 1/d^2) — rescaled to [0, 1].
            raw_pos = _solve_1d_spring_layout(panel_distances, n_phases=4)
            if raw_pos.max() > 0:
                positions = raw_pos / raw_pos.max()
            else:
                positions = np.linspace(0, 1, 4)

            # Reference dashed grid at 0, .25, .5, .75, 1
            for ref_x in REF_LINES:
                ax.axvline(ref_x, color="#cccccc", lw=0.45, ls=(0, (2, 2)),
                           alpha=0.7, zorder=0)

            # Arcs
            for i, j, d in panel_distances:
                x_i, x_j = positions[i], positions[j]
                x_lo, x_hi = sorted([x_i, x_j])
                arc_w = x_hi - x_lo
                arc_h = max(0.55 * arc_w, 0.05)
                colour = cmap_chord(float(norm_d(d)))
                arc = Arc(
                    ((x_lo + x_hi) / 2, 0),
                    width=arc_w, height=2 * arc_h,
                    theta1=0, theta2=180,
                    edgecolor=colour, linewidth=2.2,
                    zorder=3,
                )
                ax.add_patch(arc)
                ax.text(
                    (x_lo + x_hi) / 2, arc_h + 0.015,
                    f"{d:.2f}",
                    ha="center", va="bottom",
                    fontsize=5.5, color=colour, fontweight="bold",
                    zorder=4,
                )

            # Phase markers + labels at their spring-equilibrated positions
            for px_i, ph in enumerate(PHASES_ORDER):
                ax.scatter([positions[px_i]], [0], s=38, c="white",
                           edgecolor="k", linewidth=0.9, zorder=5)
            if r_idx == len(DISTANCE_KEYS) - 1:
                for px_i, ph in enumerate(PHASES_ORDER):
                    ax.text(positions[px_i], -0.10,
                            PHASE_LABEL_SHORT[ph],
                            ha="center", va="top", fontsize=7)

            # Per-row label (rotated) on the leftmost panel only
            if c_idx == 0:
                ax.text(
                    -0.32, 0.35, DISTANCE_LABEL[dist],
                    rotation=90, ha="right", va="center",
                    fontsize=10, fontweight="bold",
                    transform=ax.transData,
                )

            if r_idx == 0:
                ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10, pad=4)

            ax.set_xlim(-0.10, 1.10)
            ax.set_ylim(-0.18, 0.78)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")
            for spine in ax.spines.values():
                spine.set_visible(False)

            # Per-PANEL colorbar — small, just shows vmin/vmax for that
            # cell so each panel's colour gradient is interpretable on its
            # own scale.
            divider = _make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.04)
            sm = plt.cm.ScalarMappable(cmap=cmap_chord, norm=norm_d)
            sm.set_array([])
            cbar = fig.colorbar(sm, cax=cax)
            cbar.set_ticks([vmin, vmax])
            cbar.set_ticklabels([f"{vmin:.2f}", f"{vmax:.2f}"])
            cbar.minorticks_off()
            cbar.ax.tick_params(labelsize=5.5, pad=0.5, length=1.5)

    # Footer legend explaining the encoding
    handles_f6 = [
        Line2D([0], [0], color="#1a9850", lw=2.6,
               label="arc COLOUR = within-panel distance (green=close, red=far)"),
        Line2D([0], [0], color="#888", lw=2.6,
               label=r"phase POSITIONS = 1-D spring layout (stiffness = $1/d^{2}$)"),
        Line2D([0], [0], color="#bbb", lw=0.8, ls="--",
               label="reference grid at 0, .25, .5, .75, 1"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor="k",
               markersize=8, label="phase node"),
    ]
    fig.legend(handles=handles_f6, loc="lower center",
               bbox_to_anchor=(0.5, 0.005), ncol=len(handles_f6),
               frameon=False, fontsize=8.0)

    fig.savefig(OUT / f"fig6_4phase_geometry{SUFFIX}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT / f'fig6_4phase_geometry{SUFFIX}.pdf'}")


if __name__ == "__main__":
    main()
