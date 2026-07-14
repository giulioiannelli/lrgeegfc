#!/usr/bin/env python3
r"""fig_spatial_reach -- the physical length the diffusion CONNECTS, as a function of the
diffusion scale, from the finest functional cluster to the whole implant span.

One panel. ell(s) = max{ ||x_i - x_j|| : K_ij(s) >= 1/N } in mm -- the largest physical
distance between two contacts that COMMUNICATE (heat from one has reached its equilibrium
share at the other), K = e^{-tau L}, s = tau*lambda_max. A physical length, NOT a heat-
weighted average. It rises from the finest connected functional cluster and reaches the
implant span when the two farthest contacts equilibrate.

Anchored on the two true physical scales -- contact pitch (electrode spacing, ~3.5 mm) and
implant span (two farthest contacts, ~9-11 cm) -- and marked at s=1 = tau_min = 1/lambda_max,
the raw<->communication boundary (for s<1, e^{-tau L} ~ I - tau L reads single raw edges only).

Reads (CSV only): spatial_reach_mst020/cohort.csv, per_cell.csv.
Writes: data/preprint/figures/new_results_sec1/fig_spatial_reach.{pdf,png}
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import _common as C


def main():
    C.use_lrg_style()
    coh = pd.read_csv(C.REACH / "cohort.csv")
    pc = pd.read_csv(C.REACH / "per_cell.csv")

    pitch = float(pc["pitch_mm"].median())        # electrode spacing (a single edge)
    span = float(pc["span_mm"].median())           # two farthest contacts = whole space
    smin = float(coh.s.min())

    fig, ax = plt.subplots(figsize=(7.8, 5.3))

    # --- two TRUE physical endpoints ---
    ax.axhline(pitch, color="0.30", lw=1.1, ls=(0, (1, 1.6)), zorder=1)
    ax.axhline(span, color="0.30", lw=1.1, ls=(0, (6, 2, 1, 2)), zorder=1)
    ax.text(smin, span * 1.05, rf"implant span $\approx{span:.0f}$ mm  "
            "(the two farthest contacts $=$ whole space)", fontsize=8.4,
            color="0.30", ha="left", va="bottom")
    ax.text(smin, pitch * 0.82, rf"contact pitch $\approx{pitch:.1f}$ mm  "
            "(electrode spacing $=$ a single edge)", fontsize=8.4,
            color="0.30", ha="left", va="top")

    # --- raw (s<1) regime shading + tau_min = 1/lambda_max boundary ---
    ax.axvspan(smin * 0.9, 1.0, color="0.5", alpha=0.06, lw=0, zorder=0)
    ax.axvline(1.0, color="0.35", lw=1.1, ls=(0, (4, 2)), zorder=1)
    ax.text(smin * 1.12, span * 0.62, "raw / single-edge\nregime  ($s<1$):\n"
            r"$e^{-\tau L}\approx I-\tau L$", fontsize=8.0, color="0.40",
            ha="left", va="center")
    ax.text(1.16, span * 0.40, r"$s=1$: $\tau_{\min}=1/\lambda_{\max}$"
            "\n(raw $\\leftrightarrow$ communication)", fontsize=8.0, color="0.30",
            ha="left", va="center")

    # --- connected-length curves ell(s) [mm]; only where >=5 patients have CONNECTED
    # (below that the cohort median is 1-2 patients / not meaningful; ell=0 dead-zone
    # patients are NaN, never averaged in as "0 mm") ---
    MIN_CONN = 5
    for band in C.BANDS:
        sub = coh[(coh.band == band) & (coh.n_conn >= MIN_CONN)].sort_values("s")
        col = C.band_color(band)
        ax.fill_between(sub.s, sub.ell_lo, sub.ell_hi, color=col, alpha=0.07, lw=0)
        ax.plot(sub.s, sub.ell_med, color=col, lw=2.1, label=C.BTeX[band], zorder=3)

    ax.annotate("saturates at the span once the two\nfarthest contacts equilibrate "
                "($s\\gtrsim2$):\nphysical reach is exhausted, the rest is functional",
                xy=(45, span), xytext=(6.5, span * 0.34), fontsize=8.3, color="0.15",
                ha="left", va="top",
                arrowprops=dict(arrowstyle="-|>", color="0.15", lw=1.0,
                                connectionstyle="arc3,rad=0.25"))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(smin * 0.9, 190)
    ax.set_ylim(pitch * 0.7, span * 1.4)
    ax.set_yticks([pitch, 10, 30, span])
    ax.set_yticklabels([f"{pitch:.1f}", "10", "30", f"{span:.0f}"])
    ax.set_xlabel(C.XLAB_S, fontsize=13)
    ax.set_ylabel(r"connected physical length   $\ell(s)$   [mm]", fontsize=12)
    ax.set_title("Physical length the diffusion connects vs scale: finest functional "
                 "cluster $\\to$ implant span\n(mst@0.20, rest$_\\mathrm{post}$; cohort "
                 "median $\\pm$ IQR)", fontsize=10.0, fontweight="bold", loc="left", pad=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(loc="lower right", ncol=2, frameon=False, fontsize=9.5,
              handlelength=1.4, columnspacing=1.2, title="band", title_fontsize=9.5)

    C.FIGDIR.mkdir(parents=True, exist_ok=True)
    out = C.FIGDIR / "fig_spatial_reach.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=200, bbox_inches="tight")   # PNG on request
    plt.close(fig)
    print(f"[fig] {out} (+ .png)")


if __name__ == "__main__":
    main()
