#!/usr/bin/env python3
"""q_bundle figure figS1 -- Z inflation diagnostic.

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
    print(f"=== q_figS1 for substrate: {SUBSTRATE} ===")
    distances_long, td, contrast = load_bundle_data(SRC)



    # ----------------------------------------------------------------------
    # Fig S1 — Z inflation diagnostic (only when null artefacts exist;
    # audit_25 within-RPre split-half null is currently only run on
    # `imcoh_abs`, so figS1 is skipped for `imcoh_sq`).
    # ----------------------------------------------------------------------
    _per_patient_csv = SRC / "per_patient_with_scale.csv"
    if not _per_patient_csv.exists():
        print(f"skipping figS1 (no {_per_patient_csv.name} for {SUBSTRATE})")
    else:
        per_patient = pd.read_csv(_per_patient_csv)

        PAIR_S1 = [
            ("rest_pre", "task_test"),
            ("rest_pre", "rest_post"),
            ("task_test", "rest_post"),
        ]
        PAIR_S1_COLOUR = {
            ("rest_pre", "task_test"):  "#d62728",   # red
            ("rest_pre", "rest_post"):  "#1f77b4",   # blue
            ("task_test", "rest_post"): "#2ca02c",   # green
        }
        PAIR_S1_LABEL = {
            ("rest_pre", "task_test"):  r"$d_{\mathrm{obs}}$  RPre$\rightarrow$TT",
            ("rest_pre", "rest_post"):  r"$d_{\mathrm{obs}}$  RPre$\rightarrow$RPost",
            ("task_test", "rest_post"): r"$d_{\mathrm{obs}}$  TT$\rightarrow$RPost",
        }

        fig, axes = plt.subplots(3, 6, figsize=(15.6, 7.0), squeeze=False)
        fig.subplots_adjust(top=0.93, bottom=0.10, left=0.06, right=0.985,
                            wspace=0.18, hspace=0.32)

        for r_idx, dist in enumerate(DISTANCE_KEYS):
            for c_idx, band in enumerate(BRAIN_BANDS_NAMES):
                ax = axes[r_idx, c_idx]
                sub = per_patient[(per_patient.distance == dist)
                                  & (per_patient.band == band)]
                for p_idx, p in enumerate(PATIENTS_4PHASE):
                    pdata = sub[sub.patient == p]
                    if pdata.empty:
                        continue
                    null_q1 = float(pdata.iloc[0].null_q1)
                    null_q3 = float(pdata.iloc[0].null_q3)
                    null_median = float(pdata.iloc[0].null_median)
                    ax.add_patch(Rectangle(
                        (null_q1, p_idx - 0.32), null_q3 - null_q1, 0.64,
                        facecolor="#cccccc", edgecolor="#888",
                        linewidth=0.4, alpha=0.75, zorder=1,
                    ))
                    ax.plot([null_median, null_median],
                            [p_idx - 0.38, p_idx + 0.38],
                            color="k", lw=0.7, zorder=2)
                    for _, row in pdata.iterrows():
                        pair = (row.phase_A, row.phase_B)
                        ax.scatter([row.d_obs], [p_idx], s=22,
                                   c=PAIR_S1_COLOUR.get(pair, "#444"),
                                   edgecolor="k", linewidth=0.35,
                                   zorder=3)
                ax.set_yticks(np.arange(len(PATIENTS_4PHASE)))
                if c_idx == 0:
                    ax.set_yticklabels([p.replace("Pat_", "")
                                        for p in PATIENTS_4PHASE],
                                       fontsize=6)
                else:
                    ax.set_yticklabels([])
                ax.set_ylim(-0.6, len(PATIENTS_4PHASE) - 0.4)
                ax.invert_yaxis()
                ax.tick_params(labelsize=6.5)
                if r_idx == 0:
                    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
                ax.set_xlabel(rf"$d_{dist}$", fontsize=8, labelpad=2)
                if c_idx == 0:
                    ax.text(-0.40, 0.5, DISTANCE_LABEL[dist],
                            rotation=90, ha="right", va="center",
                            fontsize=10, fontweight="bold",
                            transform=ax.transAxes)
                ax.grid(axis="x", color="#eee", lw=0.4, zorder=0)
                ax.set_axisbelow(True)

        handles_s1 = [
            Patch(facecolor="#cccccc", edgecolor="#888",
                  label=r"within-RPre null  Q1$-$Q3 (50 splits)"),
            Line2D([0], [0], color="k", lw=0.9, label="null median"),
            Line2D([0], [0], marker="o", color="none",
                   markerfacecolor=PAIR_S1_COLOUR[("rest_pre", "task_test")],
                   markeredgecolor="k", markersize=7,
                   label=PAIR_S1_LABEL[("rest_pre", "task_test")]),
            Line2D([0], [0], marker="o", color="none",
                   markerfacecolor=PAIR_S1_COLOUR[("rest_pre", "rest_post")],
                   markeredgecolor="k", markersize=7,
                   label=PAIR_S1_LABEL[("rest_pre", "rest_post")]),
            Line2D([0], [0], marker="o", color="none",
                   markerfacecolor=PAIR_S1_COLOUR[("task_test", "rest_post")],
                   markeredgecolor="k", markersize=7,
                   label=PAIR_S1_LABEL[("task_test", "rest_post")]),
        ]
        fig.legend(handles=handles_s1, loc="lower center",
                   bbox_to_anchor=(0.5, 0.015), ncol=len(handles_s1),
                   frameon=False, fontsize=8.0)

        fig.savefig(OUT / f"figS1_z_inflation{SUFFIX}.pdf", bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {OUT / f'figS1_z_inflation{SUFFIX}.pdf'}")


if __name__ == "__main__":
    main()
