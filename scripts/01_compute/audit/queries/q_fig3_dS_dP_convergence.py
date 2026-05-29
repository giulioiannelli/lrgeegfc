#!/usr/bin/env python3
"""q_bundle figure fig3 -- d_S x d_P convergence.

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
    print(f"=== q_fig3 for substrate: {SUBSTRATE} ===")
    distances_long, td, contrast = load_bundle_data(SRC)


    # ----------------------------------------------------------------------
    # Fig 3 — d_S × d_P convergence (numbered markers, trace-zone shading)
    # ----------------------------------------------------------------------
    fig, axes = plt.subplots(1, 6, figsize=(15.6, 3.0), squeeze=False)

    for c, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[0, c]
        sub = td[td.band == band].dropna(subset=["S", "P"])
        if sub.empty:
            ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11); continue
        xs_v = sub.S.values
        ys_v = sub.P.values
        pats = sub.patient.values

        pad_pct = 0.10
        xmin, xmax = float(xs_v.min()), float(xs_v.max())
        ymin, ymax = float(ys_v.min()), float(ys_v.max())
        lo = min(xmin, ymin) - pad_pct * max(abs(xmin), abs(xmax),
                                             abs(ymin), abs(ymax), 0.05)
        hi = max(xmax, ymax) + pad_pct * max(abs(xmin), abs(xmax),
                                             abs(ymin), abs(ymax), 0.05)
        # Symmetrize a bit to keep zero centered-ish
        extent = max(abs(lo), abs(hi))
        lo, hi = -extent, extent

        # Trace zone (upper-right quadrant, both T_d > 0)
        ax.add_patch(Rectangle((0, 0), hi, hi,
                               facecolor="#2ca02c", alpha=0.07, zorder=0))
        # Crosshair at origin
        ax.axhline(0, color="#999", lw=0.6, ls="--", zorder=1)
        ax.axvline(0, color="#999", lw=0.6, ls="--", zorder=1)
        # Identity line
        ax.plot([lo, hi], [lo, hi], color="#444", lw=0.7, ls=":",
                alpha=0.7, zorder=1)

        # Cohort mean (star)
        mu_x, mu_y = float(xs_v.mean()), float(ys_v.mean())
        ax.scatter([mu_x], [mu_y], marker="*", s=240, c="#1f77b4",
                   edgecolor="k", linewidth=0.9, zorder=5, alpha=0.95)

        # Numbered patient circles
        for x, y, p in zip(xs_v, ys_v, pats):
            nid = p.replace("Pat_", "")
            in_trace = (x > 0) and (y > 0)
            is_p03 = (p == "Pat_03")
            if is_p03:
                edge = "#ff7f0e"
            elif in_trace:
                edge = "#1b5e20"
            else:
                edge = "#666"
            ax.scatter([x], [y], marker="o", s=320, c="white",
                       edgecolor=edge, linewidth=1.1, zorder=3)
            ax.text(x, y, nid, ha="center", va="center",
                    fontsize=7.2, fontweight="bold", color=edge,
                    zorder=4)

        n_trace = int(((xs_v > 0) & (ys_v > 0)).sum())
        n_total = len(xs_v)
        ax.text(0.05, 0.95, f"trace: {n_trace}/{n_total}",
                transform=ax.transAxes, fontsize=8.5, va="top",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="#999", alpha=0.92))

        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
        ax.set_xlabel(r"$T_d^{(d_S)}$  rank-only", fontsize=8.5)
        if c == 0:
            ax.set_ylabel(r"$T_d^{(d_P)}$  mag-weighted", fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(color="#eee", lw=0.4, zorder=0)
        ax.set_axisbelow(True)

    handles_f3 = [
        Line2D([0], [0], marker="*", color="none",
               markerfacecolor="#1f77b4", markeredgecolor="k",
               markersize=14, label="cohort mean"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor="#1b5e20",
               markersize=12, label="patient (in trace zone)"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor="#666",
               markersize=12, label="patient (off zone)"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor="#ff7f0e",
               markersize=12, label="Pat_03 (1024 Hz)"),
        Patch(facecolor="#2ca02c", alpha=0.13,
              label=r"trace zone ($T_d^{(d_S)}>0$ AND $T_d^{(d_P)}>0$)"),
    ]
    fig.legend(handles=handles_f3, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=len(handles_f3),
               frameon=False, fontsize=8.0)
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    fig.savefig(OUT / f"fig3_dS_dP_convergence{SUFFIX}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT / f'fig3_dS_dP_convergence{SUFFIX}.pdf'}")


if __name__ == "__main__":
    main()
