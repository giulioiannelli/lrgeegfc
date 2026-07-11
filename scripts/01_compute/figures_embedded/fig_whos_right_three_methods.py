#!/usr/bin/env python3
r"""fig:who's-right — three connectivity read-outs give three different per-band verdicts.

Talk slide 14 (the SETUP beat). The same task->rest connectivity, read three ways, disagrees
about which bands carry a trace. Each row is one method's OWN matched-strength per-band verdict,
shown BEFORE the strict nulls (drift + localization, slides 15-17) resolve the contest:

  raw edges (pairwise)        ~ undifferentiated: every band sits near threshold -> "no structure"
  Grassmann (spectral)        -> delta, beta, low-gamma
  cophenetic (multiscale)     -> alpha, beta

beta is the ONLY band lit by all three (the eventual survivor); otherwise the methods dissociate.
The figure asserts nothing about who is right -- that is the nulls' job on the next two slides.

Disc size = -log10(p_matched-strength); FILLED (band colour) = clears p<0.05, HOLLOW = does not.
Everything is surfaced from locked matched-strength cohort summaries (n=10); nothing recomputed:
  raw        data/audit/pairwise_descriptor_ladder/cohort_summary.csv  (raw_fc, T_test) gate_p
  Grassmann  data/audit/grassmann_cluster_extent/cohort_summary.csv     cluster_p_cluster_mass
  cophenetic data/audit/rho_sym_gate/cohort_summary.csv                 gate_p_sym

Writes: data/outputs/figures/talk/fig_whos_right_three_methods.pdf
"""
from __future__ import annotations

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

use_lrg_style()

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
ALPHA = 0.05
PFLOOR = 1e-3                      # -log10 ceiling 3.0 (empirical-null floor, R>=200)

# (row label, csv (rel to data/audit), p-column, optional row filter) — top row drawn first
METHODS = [
    ("raw edges\n(pairwise)",
     "pairwise_descriptor_ladder/cohort_summary.csv", "gate_p",
     lambda df: df[(df.descriptor == "raw_fc") & (df.functional == "T_test")]),
    ("Grassmann\n(spectral subspace)",
     "grassmann_cluster_extent/cohort_summary.csv", "cluster_p_cluster_mass", None),
    ("cophenetic\n(multiscale hierarchy)",
     "rho_sym_gate/cohort_summary.csv", "gate_p_sym", None),
]


def load_p():
    out = []
    for label, rel, pcol, filt in METHODS:
        df = pd.read_csv(ROOT / "data" / "audit" / rel)
        if filt is not None:
            df = filt(df)
        df = df.set_index("band")
        out.append((label, {b: float(df.loc[b, pcol]) for b in BANDS}))
    return out


def neglog(p):
    return -np.log10(max(float(p), PFLOOR))


def size_of(nl):
    return 60 + 340 * min(nl, 3.0) / 3.0      # marker area in pts^2


def main():
    data = load_p()
    nrow = len(data)

    fig, ax = plt.subplots(figsize=(7.8, 3.7))

    # subtle beta-column cue (the band the eye should hold onto — lit in every row)
    bi = BANDS.index("beta")
    ax.axvspan(bi - 0.46, bi + 0.46, color="#f2c14e", alpha=0.12, zorder=0)
    ax.annotate("lit by all three", xy=(bi, nrow - 0.5), xytext=(bi, nrow - 0.28),
                ha="center", va="bottom", fontsize=7.5, color="#8a6d1a")

    # faint per-row baselines
    for r in range(nrow):
        ax.axhline(nrow - 1 - r, color="0.90", lw=0.8, zorder=0.5)

    for r, (label, pmap) in enumerate(data):
        y = nrow - 1 - r                                    # top row = first method
        for c, b in enumerate(BANDS):
            p = pmap[b]
            nl = neglog(p)
            col = band_color(b)
            if p < ALPHA:                                    # clears the matched-strength null
                ax.scatter(c, y, s=size_of(nl), c=col, edgecolors="#20232a",
                           linewidths=1.1, zorder=4)
            else:                                            # does not clear
                ax.scatter(c, y, s=size_of(nl), facecolors="white", edgecolors=col,
                           linewidths=1.5, alpha=0.95, zorder=3)

    # axes / labels
    ax.set_xlim(-0.65, len(BANDS) - 0.35)
    ax.set_ylim(-0.6, nrow - 0.35)
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS], fontsize=14)
    for tick, b in zip(ax.get_xticklabels(), BANDS):
        tick.set_color(band_color(b))
        tick.set_fontweight("bold")
    ax.set_yticks(range(nrow))
    ax.set_yticklabels([lab for lab, _ in data][::-1], fontsize=9.5)
    # emphasise our method (bottom row)
    ax.get_yticklabels()[0].set_fontweight("bold")
    ax.tick_params(length=0)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.set_axisbelow(True)

    # legend: size = significance, fill = clears null
    ref = [(0.05, "p = 0.05"), (0.01, "p = 0.01"), (0.001, r"p $\leq$ 0.001")]
    size_h = [Line2D([], [], marker="o", ls="", markerfacecolor="0.55",
                     markeredgecolor="0.25",
                     markersize=np.sqrt(size_of(neglog(p))) * 0.9, label=lab)
              for p, lab in ref]
    fill_h = [
        Line2D([], [], marker="o", ls="", markerfacecolor="#6b7078",
               markeredgecolor="#20232a", markersize=9,
               label="clears matched-strength null (p < 0.05)"),
        Line2D([], [], marker="o", ls="", markerfacecolor="white",
               markeredgecolor="#6b7078", markersize=9, label="does not clear"),
    ]
    leg1 = fig.legend(handles=size_h, loc="lower center", bbox_to_anchor=(0.30, -0.11),
                      ncol=3, frameon=False, fontsize=8, handletextpad=0.2,
                      columnspacing=1.1, title="disc size = $-\\log_{10}p$",
                      title_fontsize=8)
    fig.add_artist(leg1)
    fig.legend(handles=fill_h, loc="lower center", bbox_to_anchor=(0.76, -0.105),
               ncol=1, frameon=False, fontsize=8, handletextpad=0.4)

    fig.tight_layout(rect=(0, 0.04, 1, 1))
    out = ROOT / "data" / "outputs" / "figures" / "talk" / "fig_whos_right_three_methods.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)

    # stdout audit table (numbers stay out of the figure; project rule)
    print("=== who's-right per-band matched-strength verdicts (n=10) ===")
    print(f"{'band':<11}" + "".join(f"{lab.splitlines()[0]:>13}" for lab, _ in data))
    for b in BANDS:
        cells = []
        for _, pmap in data:
            p = pmap[b]
            cells.append(f"{p:.3f}{'*' if p < ALPHA else ' '}")
        print(f"{b:<11}" + "".join(f"{c:>13}" for c in cells))
    print(f"\n[fig] -> {out}")


if __name__ == "__main__":
    main()
