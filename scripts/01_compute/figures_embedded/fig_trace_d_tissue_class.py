#!/usr/bin/env python3
r"""fig:trace_d — which tissue carries the trace, band by band (Results §1, para d).

CORE MESSAGE: labelling every node pair by the tissue it links and reading the cophenetic
trace per band gives a full carrier map, not a single number. Two dissociations stand out.
(1) The beta trace is a GRAY-MATTER, NON-DISEASED phenomenon: it concentrates in
gray<->gray coupling (rho_sym=0.235, clears the strength null) and in non-SOZ coupling
(0.168, clears), but NOT in white<->white or SOZ-internal coupling. (2) alpha is the
mirror image: it is carried by white matter (0.170, clears) and, above all, by
SOZ-internal coupling (0.410, the single strongest cell, clears) — alpha recruits the
diseased core that beta leaves out. delta (gray) and high_gamma (white) each clear one
within-tissue cell; theta carries nothing anywhere.

Cell = cohort rho_sym; star = clears the matched-strength null (Wilcoxon p<0.05).
Left panel stratifies by gray/white matter; right panel by the seizure-onset zone.

Reads : data/audit/epi_wm_stratified_rhosym/pairclass_rhosym_cohort.csv
Writes: data/reports/results_section1/fig_trace_d_tissue_class.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

SRC = ROOT / "data/audit/epi_wm_stratified_rhosym/pairclass_rhosym_cohort.csv"
OUT = ROOT / "data/reports/results_section1/fig_trace_d_tissue_class.pdf"

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
# (title, [(config, column-label), ...])
PANELS = [
    ("gray / white matter", "wm",
     [("gray_gray", "gray\n↔ gray"), ("cross", "gray\n↔ white"), ("wm_wm", "white\n↔ white")]),
    ("seizure-onset zone", "epi",
     [("nonepi_nonepi", "non-SOZ\n↔ non-SOZ"), ("cross", "SOZ\n↔ rest"), ("epi_epi", "SOZ\n↔ SOZ")]),
]
CMAP = plt.get_cmap("cividis")
VMIN, VMAX = -0.06, 0.42


def _text_color(rgba):
    r, g, b = rgba[:3]
    return "white" if (0.299 * r + 0.587 * g + 0.114 * b) < 0.55 else "0.1"


def main():
    d = pd.read_csv(SRC)
    norm = Normalize(VMIN, VMAX)

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.6))

    for ax, (title, strat, cols) in zip(axes, PANELS):
        sub = d[(d.stratify == strat)]
        M = np.full((len(BANDS), len(cols)), np.nan)
        clears = np.zeros_like(M, dtype=bool)
        for i, b in enumerate(BANDS):
            for j, (cfg, _) in enumerate(cols):
                r = sub[(sub.band == b) & (sub.config == cfg)]
                if len(r):
                    M[i, j] = float(r.iloc[0].cohort_obs_rho_sym)
                    clears[i, j] = float(r.iloc[0].ms_wilcoxon_p) < 0.05

        ax.imshow(M, cmap=CMAP, norm=norm, aspect="auto")
        for i in range(len(BANDS)):
            for j in range(len(cols)):
                if np.isnan(M[i, j]):
                    continue
                tc = _text_color(CMAP(norm(M[i, j])))
                ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center",
                        color=tc, fontsize=13,
                        fontweight="bold" if clears[i, j] else "normal")
                if clears[i, j]:
                    ax.scatter(j + 0.34, i - 0.32, s=120, marker="*",
                               facecolor="white", edgecolor="black", linewidth=0.7,
                               zorder=6, clip_on=False)

        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([c[1] for c in cols], fontsize=12)
        ax.set_yticks(range(len(BANDS)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=18)
        ax.set_title(title, fontsize=14, pad=10)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)

    axes[0].text(-0.30, 1.06, r"$\mathbf{d}$", transform=axes[0].transAxes,
                 fontsize=17, va="bottom", ha="left")

    # shared colorbar (never imshow_colorbar_caxdivider across panels) + star note
    cax = fig.add_axes([0.30, 0.06, 0.42, 0.032])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=CMAP), cax=cax,
                      orientation="horizontal")
    cb.set_label(r"held trace  $\rho^{\mathrm{coph}}$  (cohort)", fontsize=12)
    cb.ax.tick_params(labelsize=10)
    fig.text(0.5, 0.155, r"$\bigstar$  clears matched-strength null",
             ha="center", va="center", fontsize=11, color="0.25")

    fig.subplots_adjust(left=0.10, right=0.97, top=0.90, bottom=0.24, wspace=0.28)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:trace_d — band x tissue-class carrier map\n")
    for _, r in d[d.ms_wilcoxon_p < 0.05].iterrows():
        print(f"  clears: {r.stratify:4s} {r.config:14s} {r.band:11s} "
              f"rho={r.cohort_obs_rho_sym:+.3f} ms_p={r.ms_wilcoxon_p:.3f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
