#!/usr/bin/env python3
"""Section 3 Figure H — Dendrograms MSC vs |ImCoh| (3 patient/band combos).

Produces a 2x3 grid:
  rows  = fc_method  (MSC top, |ImCoh|=imcoh_abs bottom)
  cols  = (patient, band): (Pat_02, beta), (Pat_05, theta), (Pat_08, alpha)
  phase = rest_pre

Each panel: average-linkage dendrogram with leaf tick-labels coloured by
electrode shaft (tab20) and branches at-and-below the Ψ-optimal cut coloured
by community (tab10). Horizontal dashed line at the Ψ cut. Y-axis on log
scale, shared limits per column.

Output:
  data/outputs/figures/section3/fig_H/fig_H1_dendrograms_MSC_vs_ImCoh_3patients.pdf
  (mirrored under section3/for_writing_agent/)
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt

from _shared import apply_pub_style, save_fig, load_channel_labels
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.visuals.lrg import plot_lrg_dendrogram_shaft_colored


SECTION3_ROOT = FIGURES_ROOT / "section3"
FIG_H_DIR = SECTION3_ROOT / "fig_H"
WRITING_DIR = SECTION3_ROOT / "for_writing_agent"

FC_METHODS = [("msc", r"$\mathrm{MSC}$"),
              ("imcoh_abs", r"$|\mathrm{ImCoh}|$")]

# Per-row y-axis scale. Easy to flip back to "log" for both rows.
YSCALES = {"msc": "log", "imcoh_abs": "linear"}

TRIPLETS = [
    ("Pat_02", "beta", "rest_pre"),
    ("Pat_05", "theta", "rest_pre"),
    ("Pat_08", "alpha", "rest_pre"),
]


def _gather_ylim(lrg_results):
    """Per-column shared (tmin, tmax) using the canonical 0.8×/1.05× padding
    on the sorted merge heights — matches gen_phase_reorg_figures.py."""
    heights = []
    for r in lrg_results:
        if r is None:
            continue
        h = r.linkage_matrix[:, 2]
        h = h[h > 0]
        if h.size:
            heights.append(h)
    if not heights:
        return None
    arr = np.concatenate(heights)
    arr.sort()
    return float(arr[0] * 0.8), float(arr[-1] * 1.05)


def main(verbose: bool = False) -> Path:
    apply_pub_style()
    FIG_H_DIR.mkdir(parents=True, exist_ok=True)
    WRITING_DIR.mkdir(parents=True, exist_ok=True)

    n_rows, n_cols = len(FC_METHODS), len(TRIPLETS)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(4.6 * n_cols, 3.6 * n_rows),
        squeeze=False,
    )

    # Each panel uses its OWN tight (tmin, tmax) — sharing across rows
    # forces ImCoh's tmin to MSC's much smaller smallest-merge value, which
    # leaves a huge dead band below the lowest ImCoh branch. Per-column
    # sharing within a single fc_method would also distort, so we let the
    # helper compute its own limit per panel.

    for col, (pat, band, phase) in enumerate(TRIPLETS):
        labels = load_channel_labels(pat)
        for row, (fc_method, label) in enumerate(FC_METHODS):
            ax = axes[row, col]
            lrg = load_lrg_result(pat, phase, band, fc_method)
            if lrg is None:
                ax.set_axis_off()
                ax.text(
                    0.5, 0.5,
                    f"missing: {pat}\n{band} {phase}\n{fc_method}",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color="red",
                )
                continue

            info = plot_lrg_dendrogram_shaft_colored(
                ax, lrg, labels[: lrg.n_nodes],
                show_xlabels=True,
                leaf_font_size=3.5,
                ylim=None,
                branch_palette="tab20",
                min_n_communities=20,
                max_n_communities=max(20, lrg.n_nodes // 5),
                yscale=YSCALES.get(fc_method, "log"),
            )

            band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
            note = " (1024 Hz)" if pat == "Pat_03" else ""
            ax.set_title(
                f"{pat} — {band_tex} {phase} — {label}{note}\n"
                rf"$n^*={info['psi_n']}$",
                fontsize=10,
            )
            if row == n_rows - 1:
                ax.set_xlabel("channels (coloured by shaft)", fontsize=9)
            if col == 0:
                yscale = YSCALES.get(fc_method, "log")
                ax.set_ylabel(rf"$\Delta$ ({yscale})")
            else:
                ax.set_ylabel("")
            if verbose:
                print(f"  [{pat} {band} {phase} {fc_method}] "
                      f"n*={info['psi_n']}, "
                      f"psi_thr={info['psi_threshold']:.3g}")

    fig.tight_layout()

    out = FIG_H_DIR / "fig_H1_dendrograms_MSC_vs_ImCoh_3patients.pdf"
    save_fig(fig, out)
    shutil.copy(out, WRITING_DIR / out.name)
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(verbose=args.verbose)
