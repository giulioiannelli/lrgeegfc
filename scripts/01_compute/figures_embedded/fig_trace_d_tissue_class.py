#!/usr/bin/env python3
r"""fig:trace_d — which tissue carries the beta trace (Results §1, para d).

CORE MESSAGE: labelling every node pair by the tissue it links and reading the beta
cophenetic trace, the held trace is a GRAY-MATTER, NON-DISEASED phenomenon. It concentrates
in gray<->gray coupling (rho_sym=0.24, clears the strength null) and in mixed gray<->white
coupling (0.19, clears), but falls off in white<->white coupling (0.10, does not clear) and
carries no more than a random set of pairs inside the seizure-onset zone (SOZ<->SOZ 0.08,
does not clear). The physiology of the engaged cortex holds the trace, not the pathology that
placed the electrodes.

Two styles (same source CSV; the compact beta forest is the default and the paper panel):
    python fig_trace_d_tissue_class.py                 # -> fig_trace_d_tissue_class.pdf  (beta forest)
    python fig_trace_d_tissue_class.py --style matrix  # -> fig_trace_d_tissue_matrix.pdf (6-band x class, Supp.)
The matrix variant is the full band x tissue-class carrier map (all six bands, both
stratifications) and belongs in the Supplement; the forest is the beta-only paper panel.

Reads : data/audit/epi_wm_stratified_rhosym/pairclass_rhosym_cohort.csv
Writes: data/preprint/figures/results_section1/fig_trace_d_tissue_{class,matrix}.pdf
"""
from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

ROOT = setup_script_env()
use_lrg_style()

SRC = ROOT / "data/audit/epi_wm_stratified_rhosym/pairclass_rhosym_cohort.csv"
OUT = {
    "forest": ROOT / "data/preprint/figures/results_section1/fig_trace_d_tissue_class.pdf",
    "matrix": ROOT / "data/preprint/figures/results_section1/fig_trace_d_tissue_matrix.pdf",
}

C_CARRY, C_NULL = "#2a9d5c", "#9a9a9a"     # (kept for the matrix variant)

# General band x tissue-class carrier map: bands are color-coded markers, tissue
# classes are the y-rows (top -> bottom).
BANDS_D = ["beta", "alpha", "low_gamma", "high_gamma", "delta", "theta"]
# Band colours come from the canonical library palette (band_color); never
# hardcode them here — the whole codebase restyles from lrg_eegfc.visuals.styles.
TISSUE_ROWS = [
    (r"gray $\leftrightarrow$ gray", "wm", "gray_gray"),
    (r"gray $\leftrightarrow$ white", "wm", "cross"),
    (r"white $\leftrightarrow$ white", "wm", "wm_wm"),
    ("non-SOZ $\\leftrightarrow$\nnon-SOZ", "epi", "nonepi_nonepi"),
    (r"SOZ $\leftrightarrow$ SOZ", "epi", "epi_epi"),
]


def build_panel(ax, draw_letter=True):
    """General band x tissue-class carrier map onto ``ax`` (importable by the mosaic).

    Each tissue class is a row; each band a colored marker at its cohort rho_sym.
    Filled marker = clears the matched-strength null; open = does not.
    ``draw_letter=False`` lets a compound assembler own the panel letter instead.
    """
    d = pd.read_csv(SRC)
    n = len(TISSUE_ROWS)
    boff = dict(zip(BANDS_D, np.linspace(0.30, -0.30, len(BANDS_D))))  # beta at row top

    ax.axvline(0.0, color="0.55", lw=1.1, ls="--", zorder=1)
    for k, (label, strat, cfg) in enumerate(TISSUE_ROWS):
        y0 = n - 1 - k                                          # first row on top
        if k % 2:
            ax.axhspan(y0 - 0.5, y0 + 0.5, color="0.5", alpha=0.05, zorder=0)
        sub = d[(d.stratify == strat) & (d.config == cfg)]
        for band in BANDS_D:
            r = sub[sub.band == band]
            if not len(r):
                continue
            rho = float(r.iloc[0].cohort_obs_rho_sym)
            clears = float(r.iloc[0].ms_wilcoxon_p) < 0.05
            col, y = band_color(band), y0 + boff[band]
            if clears:
                ax.scatter(rho, y, s=150, color=col, edgecolor="black",
                           linewidth=1.0, zorder=5)
            else:
                ax.scatter(rho, y, s=64, facecolor="white", edgecolor=col,
                           linewidth=1.6, zorder=4)

    # flag the standout: alpha concentrates in the seizure-onset zone
    a_soz = d[(d.stratify == "epi") & (d.config == "epi_epi") & (d.band == "alpha")]
    if len(a_soz):
        arho = float(a_soz.iloc[0].cohort_obs_rho_sym)
        a_col = band_color("alpha", shade=0.72)   # darken the bright yellow-green for text
        ax.annotate(r"$\alpha$ holds SOZ", xy=(arho, 0 + boff["alpha"]),
                    xytext=(arho - 0.14, 0.60), fontsize=10.5, color=a_col,
                    ha="center", va="bottom", fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color="black", lw=1.3, shrinkB=6))

    ax.set_yticks(range(n))
    ax.set_yticklabels([lab for lab, _, _ in reversed(TISSUE_ROWS)], fontsize=12.5)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlim(-0.10, 0.50)
    ax.set_xlabel(r"held trace  $\rho^{\mathrm{coph}}$   (cohort)", fontsize=14)
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if draw_letter:
        ax.text(-0.02, 1.03, r"$\mathbf{b}$", transform=ax.transAxes, fontsize=17,
                va="bottom", ha="left")


def _build_forest():
    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    build_panel(ax)
    band_h = [Line2D([0], [0], marker="o", ls="", mfc=band_color(b), mec="black",
                     ms=10, label=BRAIN_BAND_TEX_DICT.get(b, b)) for b in BANDS_D]
    style_h = [
        Line2D([0], [0], marker="o", ls="", mfc="0.4", mec="black", ms=11,
               label="clears matched-strength null"),
        Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.4", mew=1.6, ms=10,
               label="does not clear"),
    ]
    fig.legend(handles=band_h + style_h, loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=4, frameon=False, fontsize=10, handletextpad=0.4, columnspacing=1.2)
    fig.subplots_adjust(left=0.195, right=0.96, top=0.90, bottom=0.26)
    return fig, None


def draw(target):
    """Render the beta forest panel into a Figure or SubFigure ``target``.

    Mirrors ``_build_forest`` but uses gridspec margins (SubFigure has no
    ``subplots_adjust``) so the panel is embeddable in a mosaic.
    """
    gs = target.add_gridspec(1, 1, left=0.26, right=0.96, top=0.90, bottom=0.26)
    ax = target.add_subplot(gs[0, 0])
    build_panel(ax)
    band_h = [Line2D([0], [0], marker="o", ls="", mfc=band_color(b), mec="black",
                     ms=10, label=BRAIN_BAND_TEX_DICT.get(b, b)) for b in BANDS_D]
    style_h = [
        Line2D([0], [0], marker="o", ls="", mfc="0.4", mec="black", ms=11,
               label="clears matched-strength null"),
        Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.4", mew=1.6, ms=10,
               label="does not clear"),
    ]
    target.legend(handles=band_h + style_h, loc="lower center", bbox_to_anchor=(0.5, -0.03),
                  ncol=4, frameon=False, fontsize=10, handletextpad=0.4, columnspacing=1.2)


# --------------------------------------------------------------------------- matrix (Supp.)
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
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


def _build_matrix():
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
    cax = fig.add_axes([0.30, 0.06, 0.42, 0.032])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=CMAP), cax=cax,
                      orientation="horizontal")
    cb.set_label(r"held trace  $\rho^{\mathrm{coph}}$  (cohort)", fontsize=12)
    cb.ax.tick_params(labelsize=10)
    fig.text(0.5, 0.155, r"$\bigstar$  clears matched-strength null",
             ha="center", va="center", fontsize=11, color="0.25")
    fig.subplots_adjust(left=0.10, right=0.97, top=0.90, bottom=0.24, wspace=0.28)
    return fig


def main(style="forest"):
    if style == "matrix":
        fig = _build_matrix()
    else:
        fig, _ = _build_forest()
    OUT[style].parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT[style], bbox_inches="tight", transparent=True)
    plt.close(fig)

    d = pd.read_csv(SRC)
    print(f"fig:trace_d — tissue carrier ({style})\n")
    for _, r in d[d.ms_wilcoxon_p < 0.05].sort_values(
            "cohort_obs_rho_sym", ascending=False).iterrows():
        print(f"  clears: {r.band:11s} {r.stratify:4s}/{r.config:14s} "
              f"rho={r.cohort_obs_rho_sym:+.3f} ms_p={r.ms_wilcoxon_p:.3f}")
    print(f"\nwrote {OUT[style]}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--style", choices=["forest", "matrix"], default="forest",
                    help="forest (default, beta paper panel) or matrix (6-band, Supp.)")
    args = ap.parse_args()
    main(args.style)
