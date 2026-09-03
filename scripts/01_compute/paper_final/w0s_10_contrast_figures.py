#!/usr/bin/env python3
"""Figure for lane W0-S part B -- why the contrast's scale slope is not evidence.

One panel-of-record, three columns, in the order the argument is made.

Left -- **a cohort-median slope has no degrees of freedom.** For each object and
band, the observed Spearman between the cohort-median margin profile and log s,
drawn on top of the null band that the *same statistic* occupies when a held-out
matched-strength realization is promoted to the observed slot. READING RULE: a
marker inside the grey band is a number the statistic produces on data with no
task information in it.

Middle -- **the effective number of independent scales, against its own noise
floor.** Observed n_eff per object and band, against the band the same
quantity occupies on held-out draws of that same object. READING RULE: a bar
below its own floor has bought nothing; a noisier object decorrelates its scales
for free, so only this comparison separates information from noise.

Right -- **what matched-strength misses.** The direction-agnostic scale
dependence slope_abs for the real data, its matched-strength null, and the
ordered sham -- a no-task five-phase arc carved from one rest recording with
block order preserved. READING RULE: where the sham sits above the
matched-strength null, matched-strength is too generous for a claim about the
SHAPE of a curve along the axis, however sound it remains for the HEIGHT of a
margin.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import band_color, use_lrg_style

ROOT = setup_script_env()

BASE = Path(os.environ.get(
    "W0S_OUT_CONTRAST",
    ROOT / "data" / "paper_final" / "lane_s_scale" / "contrast"))
FIGS = Path(os.environ.get(
    "W0S_OUT_ANALYSIS",
    ROOT / "data" / "paper_final" / "lane_s_scale")) / "figures"
#: Objects carried in the figure, in argument order. The two raw functionals the
#: contrast is built from are shown beside it so "the contrast decorrelates more
#: than its parts" is visible rather than asserted.
OBJECTS = ("T_test", "T_learn", "T_infspec", "C_learn_minus_test")
LABEL = {"T_test": r"$T_{\mathrm{test}}$", "T_learn": r"$T_{\mathrm{learn}}$",
         "T_infspec": r"$T_{\mathrm{infspec}}$",
         "C_learn_minus_test": r"$T_{\mathrm{learn}}-T_{\mathrm{test}}$"}
BANDS = ("delta", "theta", "alpha", "beta", "low_gamma", "high_gamma")
SHAM_BANDS = ("delta", "theta", "alpha", "beta")


def _tex(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


def main() -> None:
    use_lrg_style()
    FIGS.mkdir(parents=True, exist_ok=True)
    ms = pd.read_csv(BASE / "cohort_median_slope_null.csv")
    ss = pd.read_csv(BASE / "contrast_scale_structure.csv")

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
    x = np.arange(len(OBJECTS) * len(BANDS), dtype=float)
    x = x.reshape(len(OBJECTS), len(BANDS))
    x = x + np.arange(len(OBJECTS))[:, None] * 1.2          # gap between objects

    # ---- left: cohort-median slope against its own null band -------------- #
    ax = axes[0]
    for i, o in enumerate(OBJECTS):
        for j, b in enumerate(BANDS):
            r = ms[(ms.object == o) & (ms.band == b)]
            if r.empty:
                continue
            r = r.iloc[0]
            lo, hi = -r.rho_median_null_abs_p95, r.rho_median_null_abs_p95
            ax.add_patch(plt.Rectangle((x[i, j] - .42, lo), .84, hi - lo,
                                       facecolor="0.82", edgecolor="none",
                                       zorder=0))
            ax.plot([x[i, j] - .42, x[i, j] + .42],
                    [r.rho_median_null_absmed] * 2, color="0.55", lw=.8, zorder=1)
            ax.plot([x[i, j] - .42, x[i, j] + .42],
                    [-r.rho_median_null_absmed] * 2, color="0.55", lw=.8, zorder=1)
            ax.plot(x[i, j], r.rho_cohort_median, "o", ms=5.5,
                    color=band_color(b), mec="k", mew=.5, zorder=3)
    ax.axhline(0, color="k", lw=.6)
    ax.set_ylabel(r"$\rho$(cohort-median margin, $\log s$)")
    ax.set_ylim(-1.05, 1.05)

    # ---- middle: n_eff against its own noise floor ------------------------ #
    ax = axes[1]
    for i, o in enumerate(OBJECTS):
        for j, b in enumerate(BANDS):
            r = ss[(ss.object == o) & (ss.band == b)]
            if r.empty:
                continue
            r = r.iloc[0]
            ax.add_patch(plt.Rectangle((x[i, j] - .42, r.n_eff_null_med), .84,
                                       r.n_eff_null_p95 - r.n_eff_null_med,
                                       facecolor="0.82", edgecolor="none",
                                       zorder=0))
            ax.plot([x[i, j] - .42, x[i, j] + .42],
                    [r.n_eff_null_med] * 2, color="0.55", lw=.8, zorder=1)
            ax.bar(x[i, j], r.n_eff_pr, width=.7, color=band_color(b),
                   edgecolor="k", lw=.4, zorder=2)
    ax.axhline(1.0, color="k", lw=.6, ls=":")
    ax.set_ylabel(r"effective independent scales $n_{\mathrm{eff}}$")

    # ---- right: real vs matched-strength vs ordered sham ------------------ #
    ax = axes[2]
    xs = np.arange(len(OBJECTS) * len(SHAM_BANDS), dtype=float)
    xs = xs.reshape(len(OBJECTS), len(SHAM_BANDS))
    xs = xs + np.arange(len(OBJECTS))[:, None] * 1.2
    for i, o in enumerate(OBJECTS):
        for j, b in enumerate(SHAM_BANDS):
            r = ss[(ss.object == o) & (ss.band == b)]
            if r.empty:
                continue
            r = r.iloc[0]
            ax.bar(xs[i, j], r.slope_abs, width=.7, color=band_color(b),
                   edgecolor="k", lw=.4, zorder=1)
            ax.plot([xs[i, j] - .42, xs[i, j] + .42],
                    [r.slope_abs_null_med] * 2, color="k", lw=1.4, zorder=3)
            for src, mk in (("rest_pre", "^"), ("rest_post", "v")):
                k = f"sham_{src}_sham_slope_abs"
                if k in r.index and np.isfinite(r[k]):
                    ax.plot(xs[i, j], r[k], mk, ms=4.5, color="w", mec="k",
                            mew=.9, zorder=4)
    ax.set_ylabel(r"scale dependence $\overline{|\rho_k(\mathrm{margin},\log s)|}$")

    for ax, xx, bnds in ((axes[0], x, BANDS), (axes[1], x, BANDS),
                         (axes[2], xs, SHAM_BANDS)):
        ax.set_xticks(xx.ravel())
        ax.set_xticklabels([_tex(b) for b in bnds] * len(OBJECTS), fontsize=7)
        ax.set_xlim(xx.min() - .9, xx.max() + .9)
        for i, o in enumerate(OBJECTS):
            ax.text(xx[i].mean(), -0.155, LABEL[o], transform=
                    ax.get_xaxis_transform(), ha="center", va="top", fontsize=8)
        ax.tick_params(axis="x", length=0)

    h = [plt.Line2D([], [], marker="s", ls="", color="0.82", ms=9,
                    label="held-out matched-strength null (median to 95th pct)"),
         plt.Line2D([], [], color="k", lw=1.4,
                    label="matched-strength null median"),
         plt.Line2D([], [], marker="^", ls="", color="w", mec="k", mew=.9, ms=6,
                    label="ordered sham, rest_pre"),
         plt.Line2D([], [], marker="v", ls="", color="w", mec="k", mew=.9, ms=6,
                    label="ordered sham, rest_post")]
    fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, -0.10),
               ncol=len(h), frameon=False, fontsize=8)
    fig.tight_layout()
    out = FIGS / "fig_contrast_scale_slope_vs_null.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[w0s-fig] wrote {out}", flush=True)


if __name__ == "__main__":
    main()
