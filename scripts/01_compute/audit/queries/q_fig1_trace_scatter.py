#!/usr/bin/env python3
"""q_bundle figure fig1 -- trace scatter on d_S.

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
    print(f"=== q_fig1 for substrate: {SUBSTRATE} ===")
    distances_long, td, contrast = load_bundle_data(SRC)

    # ----------------------------------------------------------------------
    # Fig 1 — Trace scatter on d_S
    # ----------------------------------------------------------------------
    fig, axes = plt.subplots(1, 6, figsize=(15.6, 2.95), squeeze=False)
    sub_S = distances_long[distances_long.distance == "S"]

    for c, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[0, c]
        x_s = lookup_pair(sub_S, band, "rest_pre",  "task_test", "S")
        y_s = lookup_pair(sub_S, band, "task_test", "rest_post", "S")
        common = sorted(set(x_s.index) & set(y_s.index))
        if not common:
            ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
            continue
        x = x_s.loc[common].values
        y = y_s.loc[common].values

        span = max(x.max(), y.max()) - min(x.min(), y.min())
        lo = min(x.min(), y.min()) - 0.05 * span
        hi = max(x.max(), y.max()) + 0.05 * span

        # Trace zone (below identity) — light red shade
        poly_x = [lo, hi, hi]
        poly_y = [lo, lo, hi]
        ax.fill(poly_x, poly_y, color="#d62728", alpha=0.06, zorder=0)

        # Identity line
        ax.plot([lo, hi], [lo, hi], color="k", lw=0.7, ls="--",
                alpha=0.55, zorder=1)

        # Cohort 1σ covariance ellipse
        mu_x, mu_y = float(x.mean()), float(y.mean())
        cov = np.cov(x, y)
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

        # Per-patient circles (Pat_03 distinct)
        pat_codes = np.asarray(common)
        is_p03 = pat_codes == "Pat_03"
        ax.scatter(x[~is_p03], y[~is_p03], c="white",
                   edgecolor="#1f77b4", linewidth=1.1, s=34, zorder=3)
        if is_p03.any():
            ax.scatter(x[is_p03], y[is_p03], marker="^", c="white",
                       edgecolor="#ff7f0e", linewidth=1.2, s=46, zorder=4)

        # Cohort mean star (atop everything)
        ax.scatter([mu_x], [mu_y], marker="*", s=180, c="#2ca02c",
                   edgecolor="k", linewidth=0.9, zorder=5)

        # n_trace annotation
        n_trace = int((y < x).sum())
        n_total = len(common)
        ax.text(0.05, 0.95, f"trace: {n_trace}/{n_total}",
                transform=ax.transAxes, fontsize=8.5, va="top",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="#999", alpha=0.92))

        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
        ax.set_xlabel(r"$d_S$(RPre, TT)", fontsize=8.5)
        if c == 0:
            ax.set_ylabel(r"$d_S$(TT, RPost)", fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(color="#eee", lw=0.4, zorder=0)
        ax.set_axisbelow(True)

    handles_f1 = [
        Line2D([0], [0], marker="*", color="none",
               markerfacecolor="#2ca02c", markeredgecolor="k",
               markersize=12, label="cohort mean"),
        Patch(facecolor="#888", alpha=0.18, edgecolor="#444",
              label=r"cohort 1$\sigma$ ellipse"),
        Patch(facecolor="#d62728", alpha=0.10,
              label="trace zone (below identity)"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor="#1f77b4",
               markersize=8, label=f"patient (n={N_COHORT - 1})"),
        Line2D([0], [0], marker="^", color="none",
               markerfacecolor="white", markeredgecolor="#ff7f0e",
               markersize=9, label="Pat_03 (1024 Hz)"),
    ]
    fig.legend(handles=handles_f1, loc="lower center",
               bbox_to_anchor=(0.5, -0.03), ncol=len(handles_f1),
               frameon=False, fontsize=8.2)
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(OUT / f"fig1_trace_scatter_dS{SUFFIX}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT / f'fig1_trace_scatter_dS{SUFFIX}.pdf'}")


if __name__ == "__main__":
    main()
