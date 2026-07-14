#!/usr/bin/env python3
r"""talk_fig_pedestal_vs_peak -- detect != discriminate (S4).

Whole-task trace (T_test), cohort effect per band, two read-outs:
  raw FC edges  -> a BROAD positive PEDESTAL (all 6 bands positive; clears 4/6, and its
                   theta effect is a top-3 -- it registers movement everywhere, so it
                   cannot NAME the carrier bands);
  cophenetic    -> a SHARP alpha/beta PEAK (theta goes negative) -- it discriminates.
theta is the internal control: raw's 2nd/3rd-largest effect, cophenetic's only negative.
Filled marker = clears the matched-strength gate (p<0.05); open = does not.

Headline raw = DENSE FC (simple FC as normally computed). The same-graph control (raw read
on the identical mst@0.20 backbone) is printed for the caption: it stays non-selective too,
so the selectivity is the hierarchy read-out, not the sparsification.

Reads : data/sparsified_arc/controls_ladder/cohort_gate.csv         (raw = DENSE, headline)
        data/sparsified_arc/controls_ladder_apples/cohort_gate.csv  (raw on backbone, control)
Writes: data/outputs/figures/talk/fig_pedestal_vs_peak.png  (transparent, Canva-ready)
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

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
DENSE = C.SA / "controls_ladder" / "cohort_gate.csv"
APPLES = C.SA / "controls_ladder_apples" / "cohort_gate.csv"
C_RAW, C_COPH = "#8fb0d8", "#4fc274"      # brightened for a dark background
CAP = "0.8"                                # muted light grey (secondary text)
GRID = "0.65"                              # reference lines on dark
OUT = C.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_pedestal_vs_peak.png"


def _row(df, desc, band):
    """Return (effect above matched-strength null, gate_p) for one cell."""
    g = df[(df.descriptor == desc) & (df.band == band) & (df.functional == "T_test")]
    return float(g.med_obs.iloc[0]) - float(g.med_surr.iloc[0]), float(g.gate_p.iloc[0])


def main():
    dd = pd.read_csv(DENSE)
    ap = pd.read_csv(APPLES)
    x = np.arange(len(BANDS))
    fig, ax = plt.subplots(figsize=(9.8, 5.6))
    ax.axhline(0.0, color=GRID, lw=0.9, zorder=1)

    for i, band in enumerate(BANDS):
        ro, rp = _row(dd, "raw_fc", band)
        co, cp = _row(dd, "coph_meso", band)
        xr, xc = i - 0.15, i + 0.15
        ax.plot([xr, xr], [0, ro], color=C_RAW, lw=2.0, zorder=2, alpha=0.9)
        ax.plot([xc, xc], [0, co], color=C_COPH, lw=2.0, zorder=2, alpha=0.9)
        ax.scatter(xr, ro, s=130, facecolor=C_RAW if rp < 0.05 else "none",
                   edgecolor=C_RAW, lw=1.8, zorder=4)
        ax.scatter(xc, co, s=130, facecolor=C_COPH if cp < 0.05 else "none",
                   edgecolor=C_COPH, lw=1.8, zorder=4)

    # theta = internal control: box it
    ti = BANDS.index("theta")
    ax.axvspan(ti - 0.42, ti + 0.42, color="0.85", alpha=0.16, zorder=0)
    ax.text(ti, 0.085, "internal\ncontrol", ha="center", va="center",
            fontsize=9.5, color=CAP)

    ax.set_xticks(x)
    ax.set_xticklabels([C.BTeX.get(b, b) for b in BANDS], fontsize=15)
    ax.set_ylabel(r"cohort trace effect  (obs $-$ matched-strength null)", fontsize=12)
    ax.set_xlim(-0.6, len(BANDS) - 0.4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x", length=0)

    from matplotlib.lines import Line2D
    leg = [Line2D([0], [0], color=C_RAW, lw=2.4, marker="o", ms=9, mfc=C_RAW,
                  label="raw FC edges  —  broad pedestal (fires 4/6)"),
           Line2D([0], [0], color=C_COPH, lw=2.4, marker="o", ms=9, mfc=C_COPH,
                  label=r"cophenetic  —  sharp $\alpha/\beta$ peak (2/6)"),
           Line2D([0], [0], color=CAP, lw=0, marker="o", ms=9, mfc="none", mec=CAP,
                  label="open = misses matched-strength gate")]
    ax.legend(handles=leg, loc="upper left", frameon=False, fontsize=10.5,
              bbox_to_anchor=(0.0, 1.0))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight", dpi=200)
    plt.close(fig)

    print("DENSE (headline)  raw_fc / coph_meso  T_test:")
    for b in BANDS:
        ro, rp = _row(dd, "raw_fc", b); co, cp = _row(dd, "coph_meso", b)
        print(f"  {b:11s} raw {ro:+.3f} p={rp:.3f}{'*' if rp<.05 else ' '}   "
              f"coph {co:+.3f} p={cp:.3f}{'*' if cp<.05 else ' '}")
    print("APPLES same-graph control  raw_fc T_test (raw on backbone):")
    for b in BANDS:
        ro, rp = _row(ap, "raw_fc", b)
        print(f"  {b:11s} raw|backbone {ro:+.3f} p={rp:.3f}{'*' if rp<.05 else ' '}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
