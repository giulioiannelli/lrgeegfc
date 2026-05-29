#!/usr/bin/env python3
"""q_bundle figure figS3 -- significance counts stoplight heatmap.

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
    print(f"=== q_figS3 for substrate: {SUBSTRATE} ===")
    distances_long, td, contrast = load_bundle_data(SRC)



    # ----------------------------------------------------------------------
    # Fig S3 — Significance counts (stoplight heatmap, instant read)
    # Skipped for substrates without audit_25 cohort_summary.csv.
    # ----------------------------------------------------------------------
    _cohort_sum_csv = SRC / "cohort_summary.csv"
    if not _cohort_sum_csv.exists():
        print(f"skipping figS3 (no {_cohort_sum_csv.name} for {SUBSTRATE})")
    else:
        cohort_sum = pd.read_csv(_cohort_sum_csv)

        PAIR_S3 = [
            ("rest_pre", "task_test"),
            ("rest_pre", "rest_post"),
            ("task_test", "rest_post"),
        ]
        PAIR_S3_LABEL = {
            ("rest_pre", "task_test"):  r"RPre$\rightarrow$TT",
            ("rest_pre", "rest_post"):  r"RPre$\rightarrow$RPost",
            ("task_test", "rest_post"): r"TT$\rightarrow$RPost",
        }
        DIST_ORDER_S3 = ["S", "P", "F"]

        n_rows = len(PAIR_S3) * len(DIST_ORDER_S3)   # 9
        n_cols = len(BRAIN_BANDS_NAMES)              # 6

        matrix = np.full((n_rows, n_cols), np.nan)
        row_labels = []
        for pi, pair in enumerate(PAIR_S3):
            for di, d in enumerate(DIST_ORDER_S3):
                row_labels.append(rf"$d_{d}$  ·  {PAIR_S3_LABEL[pair]}")
                for bi, band in enumerate(BRAIN_BANDS_NAMES):
                    sel = cohort_sum[
                        (cohort_sum.distance == d) & (cohort_sum.band == band)
                        & (cohort_sum.phase_A == pair[0])
                        & (cohort_sum.phase_B == pair[1])
                    ]
                    if not sel.empty:
                        matrix[pi * 3 + di, bi] = float(sel.iloc[0].n_plus)

        # Stoplight cmap with hard thresholds: <5 red, 5-7 amber, >=8 green
        cmap_st = LinearSegmentedColormap.from_list(
            "stoplight",
            [(0.0, "#c0392b"),
             (0.49, "#c0392b"),
             (0.50, "#f1c40f"),
             (0.79, "#f1c40f"),
             (0.80, "#27ae60"),
             (1.00, "#27ae60")],
            N=256,
        )
        norm_st = Normalize(vmin=0, vmax=10)

        fig, ax = plt.subplots(figsize=(11.5, 5.6))
        fig.subplots_adjust(top=0.92, bottom=0.18, left=0.22, right=0.98)
        im = ax.imshow(matrix, cmap=cmap_st, norm=norm_st, aspect="auto",
                       interpolation="nearest")

        # Annotate each cell with the count
        for i in range(n_rows):
            for j in range(n_cols):
                val = int(matrix[i, j])
                col_text = "white" if val < 5 or val >= 8 else "#222"
                ax.text(j, i, f"{val}/10",
                        ha="center", va="center",
                        fontsize=10, fontweight="bold", color=col_text)

        # Black separator lines between phase-pair groups
        for k in [3, 6]:
            ax.axhline(k - 0.5, color="k", lw=1.4)

        ax.set_xticks(range(n_cols))
        ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                           fontsize=10)
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(row_labels, fontsize=8.5)
        ax.tick_params(axis="x", labelsize=10, pad=4)
        ax.tick_params(axis="y", labelsize=8.5, pad=2)
        ax.set_title(r"Patients with $Z > 2$ vs the within-RPre split-half null"
                     r"  (out of 10)  —  cohort screening only, NOT the trace test",
                     fontsize=10.5, pad=8)

        # Stoplight legend at bottom
        handles_s3 = [
            Patch(facecolor="#27ae60",
                  label=r"$n_+ \geq 8/10$  passes cohort screen"),
            Patch(facecolor="#f1c40f",
                  label=r"$5 \leq n_+ \leq 7/10$  marginal"),
            Patch(facecolor="#c0392b",
                  label=r"$n_+ \leq 4/10$  fails"),
        ]
        fig.legend(handles=handles_s3, loc="lower center",
                   bbox_to_anchor=(0.5, 0.015), ncol=len(handles_s3),
                   frameon=False, fontsize=9)

        fig.savefig(OUT / f"figS3_significance_counts{SUFFIX}.pdf", bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {OUT / f'figS3_significance_counts{SUFFIX}.pdf'}")


if __name__ == "__main__":
    main()
