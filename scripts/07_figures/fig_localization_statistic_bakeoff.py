#!/usr/bin/env python3
"""Verdict figure — localization cohort-statistic bake-off on mst@0.20.

Row 1: per-target heatmaps (systems x bands), colour = -log10(best-scale
sign-consistency whole-grid BH q); a white ring marks LOO-robust survivors
(q<0.05 AND LOO-worst<0.05). low_gamma lights up; beta/alpha/delta stay dark.
Row 2: the artifact — the four cohort statistics' q for a DEAD pre-registered
cell (encoding beta->OFC) vs a LIVE cell (trace low_gamma->PFC): Wilcoxon is
DOA in both (the scripts-16/17 bug), the median statistics separate them.

Reads data/sparsified_arc/localization_bakeoff_mst020/all_cells.csv (script 23).
PDF only, project style. No suptitle.
"""
from __future__ import annotations
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.visuals.styles import use_lrg_style, band_color
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT

use_lrg_style()

SRC = ROOT / "data/sparsified_arc/localization_bakeoff_mst020/all_cells.csv"
OUT = ROOT / "data/outputs/figures/localization_bakeoff"
BANDS = ["delta", "alpha", "beta", "low_gamma"]
SYSTEMS = ["OFC", "cingulate", "insula", "MTL", "PFC", "sensorimotor",
           "lateral_temporal", "parietal", "occipital"]
TARGETS = ["trace", "encoding", "inference"]
QSIG = 0.05


def load():
    df = df0 = pd.read_csv(SRC)
    df = df[df.grouping == "system"].copy()
    # whole-grid BH per (band) over all systems x scales x targets
    df["sign_qg"] = np.nan
    for band, idx in df.groupby("band").groups.items():
        s = df.loc[idx]; m = s.signconsist_p.notna()
        df.loc[s.index[m], "sign_qg"] = bh_fdr(s.loc[m, "signconsist_p"].values)
    return df


def best_cells(df):
    """Best-scale (min sign p) per (target, system, band): qg + loo."""
    g = (df.sort_values("signconsist_p")
         .groupby(["target", "unit", "band"], as_index=False).first())
    return g.set_index(["target", "unit", "band"])


def panel_heatmaps(fig, gs, cells):
    axes = []
    vmax = 1.7  # -log10(0.02)
    for ti, t in enumerate(TARGETS):
        ax = fig.add_subplot(gs[0, ti])
        M = np.full((len(SYSTEMS), len(BANDS)), np.nan)
        surv = np.zeros_like(M, bool)
        for i, u in enumerate(SYSTEMS):
            for j, b in enumerate(BANDS):
                key = (t, u, b)
                if key in cells.index:
                    r = cells.loc[key]
                    q = float(r.sign_qg)
                    M[i, j] = -np.log10(max(q, 1e-3))
                    surv[i, j] = (q < QSIG) and (float(r.loo_worst_signconsist) < QSIG)
        im = ax.imshow(M, aspect="auto", cmap="magma", vmin=0, vmax=vmax,
                       origin="upper")
        for i in range(len(SYSTEMS)):
            for j in range(len(BANDS)):
                if surv[i, j]:
                    ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                           ec="white", lw=2.2))
        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS], rotation=0)
        for j, b in enumerate(BANDS):
            ax.get_xticklabels()[j].set_color(band_color(b))
        if ti == 0:
            ax.set_yticks(range(len(SYSTEMS))); ax.set_yticklabels(SYSTEMS)
        else:
            ax.set_yticks([])
        ax.set_title(t, fontsize=11)
        axes.append((ax, im))
    # shared colorbar (row-spanning -> explicit add_axes, NOT caxdivider)
    cax = fig.add_axes([0.92, 0.56, 0.014, 0.30])
    cb = fig.colorbar(axes[-1][1], cax=cax)
    cb.set_label(r"$-\log_{10}\,q$ (sign-consistency, whole-grid BH)", fontsize=9)
    cb.ax.axhline(-np.log10(QSIG), color="white", lw=1.4)
    return axes


def panel_bakeoff(fig, gs, df):
    """Four-statistic q for a DEAD vs LIVE cell — shows Wilcoxon is the broken one."""
    ax = fig.add_subplot(gs[1, :])
    exemplars = [("encoding", "OFC", "beta", "encoding  β→OFC  (pre-registered, DEAD)"),
                 ("trace", "PFC", "low_gamma", "trace  low-γ→PFC  (recovered, LIVE)")]
    stats = [("pooled_p", "pooled"), ("wilcoxon_p", "Wilcoxon\n(scripts 16/17)"),
             ("stouffer_p", "Stouffer"), ("signconsist_p", "sign-consistency")]
    x = np.arange(len(stats)); w = 0.36
    for e, (t, u, b, lab) in enumerate(exemplars):
        s = df[(df.target == t) & (df.unit == u) & (df.band == b)]
        best = s.loc[s.signconsist_p.idxmin()]
        # per-cell qg for each statistic at that cell's scale (whole-grid recompute)
        qs = []
        for pc, _ in stats:
            qcol = pc.replace("_p", "_qgc")
            dd = df[df.band == b].copy()
            dd[qcol] = bh_fdr(dd[pc].fillna(1).values)
            row = dd[(dd.target == t) & (dd.unit == u) & (dd.scale_idx == best.scale_idx)]
            qs.append(float(row[qcol].iloc[0]) if len(row) else np.nan)
        ax.bar(x + (e - 0.5) * w, [-np.log10(max(q, 1e-3)) for q in qs], w,
               color=(band_color(b)), alpha=0.55 + 0.3 * e, label=lab,
               edgecolor="k", lw=0.5)
    ax.axhline(-np.log10(QSIG), color="0.3", ls="--", lw=1.0)
    ax.text(len(stats) - 0.4, -np.log10(QSIG) + 0.04, "q=0.05", fontsize=8, color="0.3")
    ax.set_xticks(x); ax.set_xticklabels([s[1] for s in stats], fontsize=9)
    ax.set_ylabel(r"$-\log_{10}\,q$", fontsize=9)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.set_title("cohort statistic decides the verdict — Wilcoxon fails BOTH; "
                 "median statistics separate dead (β) from live (low-γ)", fontsize=9.5)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = load()
    cells = best_cells(df)
    fig = plt.figure(figsize=(11, 8.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1.0], hspace=0.42, wspace=0.10,
                          left=0.11, right=0.90, top=0.93, bottom=0.09)
    panel_heatmaps(fig, gs, cells)
    panel_bakeoff(fig, gs, df)
    out = OUT / "fig_localization_statistic_bakeoff.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    print(f"[fig] -> {out}")


if __name__ == "__main__":
    main()
