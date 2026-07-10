#!/usr/bin/env python3
r"""fig:enc_dynamics — the §2 "dynamics" compound (Fig 3): attractor over forest.

COMPOUND, built entirely in matplotlib (no LaTeX minipage). Two stacked sub-figures:

  a  the consolidation trace as a live windowed state-space trajectory (two exemplars +
     the verified beta inference-specific cohort strip)         -- fig_arc_attractor_3d
  b  the encoding-vs-inference decomposition forest             -- fig_arc_a_...forest

Both sub-panels are drawn by the SINGLE-figure scripts' own helpers, imported here and
composed onto two matplotlib SubFigures; nothing is re-derived. Edit the science in the
single-figure scripts; edit only the mosaic here.

Reads : (delegated) data/audit/windowed_attractor/*, data/audit/consolidation_arc_rhosym/*
Writes: data/preprint/figures/results_section2/fig_arc_enc_dynamics.pdf
"""
from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
import fig_arc_a_decomposition_forest as forest  # noqa: E402
import fig_arc_attractor_3d as att  # noqa: E402

OUT = ROOT / "data/preprint/figures/results_section2/fig_arc_enc_dynamics.pdf"


def draw_attractor(sf):
    """Sub-figure a: two 3-D portraits + time-gradient bar + beta inference cohort strip."""
    nb = pd.read_csv(att.NULL)
    b = nb[nb.band == "beta"]
    verd = pd.DataFrame({"patient": b.patient.values, "obs": b["T_infspec_pe_obs"].values,
                         "surr": b["T_infspec_pe_surr_p50"].values,
                         "p": b["T_infspec_pe_p"].values})

    # row 1 is a dedicated spacer band: the time colorbar drops into it, clear of the
    # 3-D plots above and the cohort strip below.
    gs = sf.add_gridspec(3, 2, height_ratios=[2.55, 0.42, 1.05], hspace=0.05,
                         wspace=0.02, left=0.02, right=0.99, top=1.0, bottom=0.09)
    axL = sf.add_subplot(gs[0, 0], projection="3d")
    axR = sf.add_subplot(gs[0, 1], projection="3d")
    att.portrait(axL, att.HERO_TRACE, r"Pat_06 $\cdot$ tracer")
    att.portrait(axR, att.HERO_RESET, r"Pat_10 $\cdot$ resetter")
    axc = sf.add_subplot(gs[2, :])
    n_tr, n = att.cohort_strip(axc, verd)

    # time gradient bar (VECTOR pcolormesh) in the spacer band, phase transitions marked
    cax = sf.add_axes([0.30, 0.368, 0.40, 0.019])
    g = np.linspace(0, 1, 256)[None, :]
    cax.pcolormesh(np.linspace(0, 1, 257), np.array([0, 1]), g, cmap=att.TCMAP,
                   rasterized=False, shading="flat")
    bl = []
    for pat in (att.HERO_TRACE, att.HERO_RESET):
        _, ph, _ = att.embed_ei(pat)
        cnt = np.array([np.sum(ph == p) for p in att.PH], float)
        bl.append(np.cumsum(cnt) / cnt.sum())
    cum = np.mean(bl, axis=0)
    mids = 0.5 * (np.concatenate([[0.0], cum[:-1]]) + cum)
    cax.vlines(cum[:-1], 0, 1, colors="white", lw=1.6, zorder=3)
    cax.set_yticks([]); cax.set_xlim(0, 1); cax.set_xticks(mids)
    cax.set_xticklabels([att.PHLAB[p] for p in att.PH], fontsize=11)
    cax.tick_params(length=0)
    cax.set_title("window time  (one long recording)", fontsize=11, pad=4)
    for s in cax.spines.values():
        s.set_visible(False)

    handles = [Line2D([0], [0], marker="o", ls="", mfc=att.C_TRACE, mec="white", ms=10,
                      label=r"clears $\beta$ inference null"),
               Line2D([0], [0], marker="o", ls="", mfc=att.C_RESET, mec="white", ms=10,
                      label="within null"),
               Line2D([0], [0], marker="o", ls="", mfc="none", mec="0.25", ms=10,
                      label="phase-average (ghost basin)")]
    sf.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005),
              ncol=3, frameon=False, fontsize=9.5, handletextpad=0.5, columnspacing=1.8)
    sf.text(0.008, 0.985, r"$\mathbf{a}$", fontsize=20, va="top", ha="left",
            fontweight="bold")
    return n_tr, n


def draw_forest(sf):
    """Sub-figure b (right of a): encoding echo (top) over inference-specific (bottom),
    stacked and compact so the whole figure reads landscape. Violin style matches the
    §1 band figure (gate-tinted green when passing, per-band markers)."""
    cohort = pd.read_csv(forest.BASE / "arc_cohort_verdict.csv")
    per_pat = pd.read_csv(forest.BASE / "arc_per_patient.csv")
    per_null = pd.read_csv(forest.BASE / "arc_null_per_patient.csv")

    enc = cohort[cohort.functional == "T_learn"].sort_values("med_obs", ascending=False)
    order = enc.band.tolist()
    ypos = {b: len(order) - 1 - i for i, b in enumerate(order)}

    gs = sf.add_gridspec(2, 1, left=0.20, right=0.965, top=0.955, bottom=0.135,
                         hspace=0.34)
    for row, (func, val_col, p_col, xlab) in enumerate(forest.PANELS):
        ax = sf.add_subplot(gs[row, 0])
        forest.draw_panel(ax, func, val_col, p_col, order, ypos, cohort, per_pat,
                          per_null, dot_s=68, med_s=148, p_fs=9.5)
        ax.set_xlabel(xlab, fontsize=11.5)
        ax.tick_params(axis="x", labelsize=10.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_yticks(list(ypos.values()))
        ax.set_yticklabels([forest.BRAIN_BAND_TEX_DICT[b] for b in order], fontsize=15)
        ax.set_ylim(-0.8, len(order) - 0.2)
        ax.tick_params(axis="y", length=0)

    handles = [
        Line2D([0], [0], marker="o", ls="", mfc=forest.C_TRACE, mec="white", ms=9,
               label="clears own null (trace)"),
        Line2D([0], [0], marker="o", ls="", mfc=forest.C_ANTI, mec="white", ms=9,
               label=r"reset ($T<0$)"),
        Line2D([0], [0], marker="o", ls="", mfc=forest.C_UNDET, mec="white", ms=9,
               label="undetermined"),
        Line2D([0], [0], marker="o", ls="-", color=forest.GATE_PASS, mfc=forest.GATE_PASS,
               mec="white", ms=10, lw=2.6, label="cohort median + gate $p$"),
        Patch(facecolor=forest.GATE_PASS, edgecolor=forest.GATE_PASS, alpha=0.20,
              label="passing-band distribution"),
    ]
    sf.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005),
              ncol=2, frameon=False, fontsize=9, handletextpad=0.5, columnspacing=1.2)
    sf.text(0.015, 0.99, r"$\mathbf{b}$", fontsize=20, va="top", ha="left",
            fontweight="bold")


def main():
    fig = plt.figure(figsize=(17.6, 9.7))
    sf_a, sf_b = fig.subfigures(1, 2, width_ratios=[1.72, 1.0], wspace=0.01)
    n_tr, n = draw_attractor(sf_a)
    draw_forest(sf_b)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True, pad_inches=0.05)
    plt.close(fig)
    print("fig:enc_dynamics — compound (attractor | decomposition forest, side by side)")
    print(f"  attractor cohort: {n_tr}/{n} clear the beta inference null")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
