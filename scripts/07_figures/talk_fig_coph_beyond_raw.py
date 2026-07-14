#!/usr/bin/env python3
r"""talk_fig_coph_beyond_raw -- higher-order beta, complementary not redundant (S6).

Is the cophenetic hierarchy just a repackaging of the pairwise edges? No -- neither view
subsumes the other:
  coph|raw = cophenetic structure BEYOND raw FC (partial, matched-strength null)
  raw|coph = raw-FC structure BEYOND cophenetic
Left  : beta only, both residuals vs scale s -- coph|raw GROWS toward the coarse scale
        (+0.15 -> +0.24): the higher-order beta structure is a COARSE-scale phenomenon;
        raw|coph stays positive too (the edges keep their own persistent structure).
Right : all bands at the coarse scale s=30 -- coph|raw and raw|coph both positive across
        the board = two complementary lenses; beta's coph|raw is the tallest.

Reads : data/sparsified_arc/coph_beyond_raw_s{01.0,05.6,30.0}/cohort.csv
Writes: data/outputs/figures/talk/fig_coph_beyond_raw.png  (transparent, Canva-ready)
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] /
                       "01_compute" / "figures_embedded" / "new_results_sec1"))
import _common as C  # noqa: E402

C.use_lrg_style()

# dark-slide mode: all ink (text/axes/spines/ticks) white, colored data untouched
INK = "white"
plt.rcParams.update({
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
    "axes.titlecolor": INK, "xtick.color": INK, "ytick.color": INK,
})

SCALES = [(1.0, "coph_beyond_raw_s01.0"), (5.6, "coph_beyond_raw_s05.6"),
          (30.0, "coph_beyond_raw_s30.0")]
BANDS = ["delta", "alpha", "beta", "low_gamma"]     # the 4 bands this deliverable ran
C_COPH, C_RAW = "#4fc274", "#8fb0d8"      # brightened for a dark background
CAP = "0.8"                                # muted light grey (secondary text)
GRID = "0.65"                              # reference lines on dark
OUT = C.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_coph_beyond_raw.png"


def load():
    rows = []
    for s, d in SCALES:
        f = C.SA / d / "cohort.csv"
        df = pd.read_csv(f); df["s"] = s
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def main():
    df = load()
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.4, 5.4),
                                   gridspec_kw=dict(width_ratios=[1.0, 1.25], wspace=0.28))

    # -- left: beta only, coph|raw & raw|coph vs scale --
    b = df[df.band == "beta"]
    for meas, col, lab in [("coph|raw", C_COPH, r"coph $\mid$ raw  (higher-order)"),
                           ("raw|coph", C_RAW, r"raw $\mid$ coph")]:
        g = b[b.measure == meas].sort_values("s")
        axL.plot(g.s, g.obs_med, "-o", color=col, lw=2.4, ms=8, label=lab, zorder=3,
                 mfc=col, mec="white", mew=0.8)
        for _, r in g.iterrows():
            if r["q"] < 0.05:
                axL.scatter(r.s, r.obs_med, s=190, facecolor="none", edgecolor=col,
                            lw=1.6, zorder=2)
    axL.set_xscale("log"); axL.axhline(0, color=GRID, lw=0.8)
    axL.set_xlabel(C.XLAB_S, fontsize=12)
    axL.set_ylabel("residual trace effect  (beyond the other read-out)", fontsize=11)
    axL.set_title(r"$\beta$: cophenetic's extra grows coarse", fontsize=13,
                  color=C.band_color("beta"))
    axL.annotate("", xy=(29, 0.238), xytext=(6.2, 0.150),
                 arrowprops=dict(arrowstyle="->", color=C_COPH, lw=1.6))
    axL.text(13.5, 0.215, "higher-order,\ncoarse-scale", color=C_COPH, fontsize=10,
             ha="center", va="center")
    axL.legend(loc="lower left", frameon=False, fontsize=10.5)
    for sp in ("top", "right"):
        axL.spines[sp].set_visible(False)

    # -- right: all bands at coarse s=30 --
    c = df[np.isclose(df.s, 30.0)]
    x = np.arange(len(BANDS)); w = 0.36
    for k, (meas, col, off) in enumerate([("coph|raw", C_COPH, -w / 2),
                                          ("raw|coph", C_RAW, +w / 2)]):
        vals, qs = [], []
        for band in BANDS:
            r = c[(c.band == band) & (c.measure == meas)]
            vals.append(float(r.obs_med.iloc[0]) if len(r) else np.nan)
            qs.append(float(r["q"].iloc[0]) if len(r) else 1.0)
        bars = axR.bar(x + off, vals, w, color=col, edgecolor="white", lw=0.6,
                       label=meas, zorder=2, alpha=0.92)
        for xi, v, q in zip(x + off, vals, qs):
            if q < 0.05:
                axR.text(xi, v + 0.006, "*", ha="center", va="bottom", color=col, fontsize=14)
    axR.axhline(0, color=GRID, lw=0.8)
    axR.set_xticks(x); axR.set_xticklabels([C.BTeX.get(b, b) for b in BANDS], fontsize=14)
    axR.set_ylabel("residual trace effect", fontsize=11)
    axR.set_title(r"coarse scale ($s=30$): each sees what the other misses",
                  fontsize=13, color=CAP)
    axR.legend(loc="upper right", frameon=False, fontsize=10.5,
               title=r"* clears null ($q<.05$)", title_fontsize=9)
    for sp in ("top", "right"):
        axR.spines[sp].set_visible(False)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("beta coph|raw vs scale:", b[b.measure == "coph|raw"].sort_values("s")
          [["s", "obs_med", "q"]].to_dict("records"))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
