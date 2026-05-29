#!/usr/bin/env python3
"""q_bundle figure fig4 -- structural vs drift T_d^(d_S) vs T_d^(d_F).

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
    print(f"=== q_fig4 for substrate: {SUBSTRATE} ===")
    distances_long, td, contrast = load_bundle_data(SRC)



    # ----------------------------------------------------------------------
    # Fig 4 — Structural vs drift (T_d^(d_S) vs T_d^(d_F))
    # ----------------------------------------------------------------------
    fig, axes = plt.subplots(1, 6, figsize=(15.6, 3.0), squeeze=False)

    for c, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[0, c]
        sub = td[td.band == band].dropna(subset=["S", "F"])
        if sub.empty:
            ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11); continue
        xs_v = sub.S.values
        ys_v = sub.F.values
        pats = sub.patient.values

        # Symmetric square axes
        extent = max(np.abs(xs_v).max(), np.abs(ys_v).max())
        extent = extent * 1.15  # padding
        lo, hi = -extent, extent

        # Trace zone (upper-right quadrant, both T_d > 0)
        ax.add_patch(Rectangle((0, 0), hi, hi,
                               facecolor="#2ca02c", alpha=0.07, zorder=0))
        # Crosshair at origin
        ax.axhline(0, color="#999", lw=0.6, ls="--", zorder=1)
        ax.axvline(0, color="#999", lw=0.6, ls="--", zorder=1)
        # Identity line — dotted
        ax.plot([lo, hi], [lo, hi], color="#444", lw=0.8, ls=":",
                alpha=0.75, zorder=1)

        # Cohort 1σ covariance ellipse — shows the elongation along
        # identity == strong rank/amplitude correlation
        mu_x, mu_y = float(xs_v.mean()), float(ys_v.mean())
        cov = np.cov(xs_v, ys_v)
        vals, vecs = np.linalg.eigh(cov)
        order = np.argsort(vals)[::-1]
        vals = vals[order]
        vecs = vecs[:, order]
        angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
        width, height = 2.0 * np.sqrt(np.maximum(vals, 0.0))
        if width > 0 and height > 0:
            ell = Ellipse(
                (mu_x, mu_y), width, height, angle=angle,
                facecolor="#888", alpha=0.18, edgecolor="#444",
                linewidth=0.7, zorder=2,
            )
            ax.add_patch(ell)

        # Cohort mean star (atop ellipse)
        ax.scatter([mu_x], [mu_y], marker="*", s=240, c="#1f77b4",
                   edgecolor="k", linewidth=0.9, zorder=5, alpha=0.95)

        # Numbered patient circles (last 2 digits)
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

        # Inset: per-band Spearman ρ + sign agreement (the load-bearing
        # numbers for the "this is not drift" call)
        rho = float(contrast.loc[band, "spearman_rho_S_F"])
        sign_str = str(contrast.loc[band, "sign_agree_S_F"])
        n_trace = int(((xs_v > 0) & (ys_v > 0)).sum())
        n_total = len(xs_v)
        ax.text(0.05, 0.95,
                (f"trace: {n_trace}/{n_total}\n"
                 rf"$\rho_S$ = {rho:+.2f}" + "\n"
                 f"sign: {sign_str}"),
                transform=ax.transAxes, fontsize=7.5, va="top",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="#999", alpha=0.92))

        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
        ax.set_xlabel(r"$T_d^{(d_S)}$  rank-only", fontsize=8.5)
        if c == 0:
            ax.set_ylabel(r"$T_d^{(d_F)}$  amplitude-only", fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(color="#eee", lw=0.4, zorder=0)
        ax.set_axisbelow(True)

    handles_f4 = [
        Line2D([0], [0], marker="*", color="none",
               markerfacecolor="#1f77b4", markeredgecolor="k",
               markersize=14, label="cohort mean"),
        Patch(facecolor="#888", alpha=0.18, edgecolor="#444",
              label=r"cohort 1$\sigma$ ellipse (elongation $\sim$ correlation)"),
        Line2D([0], [0], color="#444", lw=0.9, ls=":",
               label="identity (perfect rank/amp agreement)"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor="#1b5e20",
               markersize=12, label="patient (in trace zone)"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor="#666",
               markersize=12, label="patient (off zone)"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor="#ff7f0e",
               markersize=12, label="Pat_03 (1024 Hz)"),
    ]
    fig.legend(handles=handles_f4, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=3,
               frameon=False, fontsize=8.0)
    fig.tight_layout(rect=[0, 0.10, 1, 1])
    fig.savefig(OUT / f"fig4_structural_vs_drift{SUFFIX}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT / f'fig4_structural_vs_drift{SUFFIX}.pdf'}")


if __name__ == "__main__":
    main()
