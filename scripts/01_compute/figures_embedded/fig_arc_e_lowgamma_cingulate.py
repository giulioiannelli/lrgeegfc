#!/usr/bin/env python3
r"""fig:arc_e — a focal low-gamma memory trace the whole-brain view misses (§2, R2.5).

CORE MESSAGE (the method showcase): low-gamma carries NO whole-brain cohort trace, yet
reading one region at a time -- an ABSOLUTE within-region concordance, no baseline
subtracted -- surfaces a genuine, focal ENCODING (memory) trace on cingulate edges and on
no other system. rho_sym = +0.25 against a surrogate median near zero (BH q = 0.035 across
systems, positive in all 8 patients sampled there), rising to +0.39 with the SOZ excluded
-- a LARGER effect on FEWER contacts, so not a product of epileptic tissue. It is specific
to memory: on the same cingulate edges the inference component is null and the standard
task change does not clear. The clearest case of the hierarchy resolving a single
structure that whole-brain averaging conceals.

Reads : data/audit/inference_localization_rhosym/within_system_trace_rhosym_{include,exclude}.csv
Writes: data/preprint/figures/results_section2/fig_arc_e_lowgamma_cingulate.pdf
"""
from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import band_color, use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
try:
    from fig_trace_b_ofc_localization import pooled_coords_systems  # type: ignore
    from nilearn.plotting import plot_glass_brain
    _BRAIN_OK = True
except Exception as exc:  # pragma: no cover
    print(f"[warn] glass-brain machinery unavailable ({exc}); shipping panel a only")
    _BRAIN_OK = False

BASE = ROOT / "data/audit/inference_localization_rhosym"
OUT = ROOT / "data/preprint/figures/results_section2/fig_arc_e_lowgamma_cingulate.pdf"

BAND, TARGET = "low_gamma", "encoding"
WINNER = "cingulate"
DROP = {"other", "non_anatomical"}
# highlight uses the low-gamma BAND colour (blue), matching the band chip in the compound
# (Fig 3 palette identity: slow=red -> fast=blue), not the encoding-green.
C_WIN, C_NS = band_color("low_gamma", 0.80), "#9a9a9a"
Q_SIG = 0.05


def load(inc: str) -> pd.DataFrame:
    df = pd.read_csv(BASE / f"within_system_trace_rhosym_{inc}.csv")
    d = df[(df.band == BAND) & (df.target == TARGET) &
           (~df.system.isin(DROP))].copy()
    return d.set_index("system")


def draw_lollipop(ax, inc, exc):
    order = inc.sort_values("median_obs_rho").index.tolist()  # winner at top
    y0 = np.arange(len(order))
    ax.axvline(0, color="0.55", lw=1.0, zorder=1)
    for y, s in zip(y0, order):
        rec = inc.loc[s]
        sig = rec.bh_q < Q_SIG
        col = C_WIN if (s == WINNER and sig) else (C_WIN if sig else C_NS)
        ax.plot([0, rec.median_obs_rho], [y, y], color=col, lw=2.2, zorder=2,
                solid_capstyle="round")
        # SOZ-excluded replicate (pale shadow): shows the effect grows on fewer contacts
        if s in exc.index:
            ax.scatter([exc.loc[s, "median_obs_rho"]], [y + 0.17], s=18, marker="o",
                       facecolors="none", edgecolors=col, linewidths=0.9, alpha=0.6,
                       zorder=3)
        ax.scatter([rec.median_obs_rho], [y],
                   s=165 if s == WINNER else 95, marker="o",
                   facecolors=col if sig else "white", edgecolors=col,
                   linewidths=1.6 if sig else 1.4, zorder=4)
        if sig:
            ax.annotate(r"$\ast$", (rec.median_obs_rho, y), textcoords="offset points",
                        xytext=(9, 3), fontsize=13, color=col, zorder=5)
    # winner annotation: q only (patient count + SOZ-excl growth removed -- they crowd the
    # branch above; both are carried by the caption / prose)
    wy = order.index(WINNER)
    wr = inc.loc[WINNER]
    ax.annotate(rf"$q={wr.bh_q:.3f}$", (wr.median_obs_rho, wy),
                textcoords="offset points", xytext=(0, 15), fontsize=9.5,
                color=C_WIN, ha="center")

    ax.set_yticks(y0)
    ax.set_yticklabels(order)
    for lab in ax.get_yticklabels():
        if lab.get_text() == WINNER:
            lab.set_fontweight("bold")
            lab.set_color(C_WIN)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlim(-0.12, 0.47)
    ax.set_xlabel(r"within-region memory trace  $\rho_{\mathrm{sym}}$"
                  r"   (absolute, no baseline subtracted)", fontsize=11.5)
    ax.tick_params(axis="y", length=0)
    ax.spines[["top", "right"]].set_visible(False)
    # the "averaging hides it" note: no whole-brain low-gamma cohort trace
    ax.text(0.015, -0.62, "whole-brain low-$\\gamma$: no cohort trace",
            fontsize=9.2, color="0.4", style="italic", ha="left")


def draw_brain(fig, rect, coords, systems, sig_win):
    disp = plot_glass_brain(None, display_mode="lzr", axes=rect, figure=fig)
    m_win = np.isin(systems, [WINNER])
    m_rest = ~m_win
    if m_rest.any():
        disp.add_markers(coords[m_rest], marker_color=C_NS, marker_size=6, alpha=0.18)
    if sig_win and m_win.any():
        disp.add_markers(coords[m_win], marker_color=C_WIN, marker_size=34, alpha=0.95,
                         edgecolors="black", linewidths=1.2)
    return disp


def main():
    inc, exc = load("include"), load("exclude")
    sig_win = bool(inc.loc[WINNER, "bh_q"] < Q_SIG)

    fig = plt.figure(figsize=(11.6, 5.5))
    axa = fig.add_axes([0.14, 0.18, 0.40, 0.74])
    draw_lollipop(axa, inc, exc)
    # tile letter supplied by the LaTeX mosaic in results_sec_2.tex

    if _BRAIN_OK:
        try:
            coords, systems = pooled_coords_systems()
            draw_brain(fig, (0.57, 0.10, 0.41, 0.82), coords, systems, sig_win)
            # tile letter supplied by the LaTeX mosaic
            fig.text(0.775, 0.13, r"memory (encoding) $\rightarrow$ cingulate, low-$\gamma$",
                     ha="center", fontsize=11.5, color=C_WIN, fontweight="bold")
        except Exception as exc_b:  # pragma: no cover
            print(f"[warn] brain panel failed ({exc_b}); shipping panel a only")

    handles = [
        Line2D([], [], color=C_WIN, marker="o", ls="", mfc=C_WIN, mec="black",
               mew=1.2, ms=11, label=r"cingulate memory trace, BH $q<0.05$ ($\ast$)"),
        Line2D([], [], color=C_NS, marker="o", ls="", mfc="white", mec=C_NS,
               mew=1.4, ms=9, label=r"other system, n.s. ($q\geq0.05$)"),
        Line2D([], [], color=C_WIN, marker="o", ls="", mfc="none", mec=C_WIN,
               mew=0.9, ms=6, alpha=0.6, label="SOZ-excluded replicate"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=3, frameon=False, fontsize=9.5, handletextpad=0.5,
               columnspacing=1.5)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight")
    plt.close(fig)

    print(f"fig:arc_e — low-gamma {TARGET} within-region trace")
    for s in inc.sort_values("median_obs_rho", ascending=False).index:
        print(f"  {s:16s} rho_inc={inc.loc[s,'median_obs_rho']:+.3f} "
              f"rho_exc={exc.loc[s,'median_obs_rho']:+.3f} "
              f"q={inc.loc[s,'bh_q']:.3f} n_pos={int(inc.loc[s,'n_pos'])}/"
              f"{int(inc.loc[s,'K_implanted'])}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
