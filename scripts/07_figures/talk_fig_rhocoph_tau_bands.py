#!/usr/bin/env python3
r"""talk_fig_rhocoph_tau_bands -- cohort rho^coph vs scale tau, all bands (the "lasting
trace" band-comparison tool).

One curve per band: cohort-median rho_sym^coph(s) across the 16-scale mst@0.20 diffusion
sweep (s = tau*lambda_max), with the matched-strength surrogate floor and a filled marker
at every scale that clears the cohort gate (p<0.05). Reading it as a band comparator:
  beta  -> held HIGH at every scale (scale-invariant; 16/16),
  alpha -> a mesoscale bump (12/16),
  theta / low-gamma -> ride the surrogate floor (null; the internal controls),
  delta / high-gamma -> only fine-scale, decaying (scale-artifacts).

Reads : data/sparsified_arc/ms_mst020/cohort_gate.csv (band x 16 scales; obs_med, surr_med, gate_p).
Writes: data/outputs/figures/talk/fig_rhocoph_tau_bands.png  (transparent, Canva-ready)
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

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
CARRIER = {"alpha", "beta"}
CAP = "0.55"
OUT = C.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_rhocoph_tau_bands.png"


def main():
    df = pd.read_csv(C.MS / "cohort_gate.csv")
    fig, ax = plt.subplots(figsize=(9.6, 6.0))

    surr_lo, surr_hi = [], []
    for band in BANDS:
        g = df[df.band == band].sort_values("s")
        s, obs, gp = g.s.values, g.obs_med.values, g.gate_p.values
        col = C.band_color(band)
        carrier = band in CARRIER
        lw = 2.6 if carrier else 1.6
        alpha = 1.0 if carrier else 0.75
        ax.plot(s, obs, "-", color=col, lw=lw, alpha=alpha, zorder=3 if carrier else 2)
        sig = gp < 0.05
        ax.scatter(s[sig], obs[sig], s=46 if carrier else 30, color=col,
                   edgecolor="white", lw=0.7, zorder=4)
        ax.scatter(s[~sig], obs[~sig], s=22, facecolor="none", edgecolor=col,
                   lw=1.0, alpha=0.7, zorder=3)
        surr_lo.append(g.surr_med.values); surr_hi.append(g.surr_med.values)

    # surrogate floor: envelope of the per-band matched-strength medians
    s0 = df[df.band == "beta"].sort_values("s").s.values
    slo = np.min(np.vstack(surr_lo), axis=0); shi = np.max(np.vstack(surr_hi), axis=0)
    ax.fill_between(s0, slo, shi, color="0.75", alpha=0.35, zorder=1, lw=0)
    ax.plot(s0, np.median(np.vstack(surr_hi), axis=0), color="0.6", lw=1.0, ls=":", zorder=1)
    ax.text(s0[-1], np.median(np.vstack(surr_hi), axis=0)[-1], "  matched-strength null",
            color=CAP, fontsize=9, va="center", ha="left")

    ax.set_xscale("log")
    ax.axhline(0.0, color="0.7", lw=0.8, zorder=0)
    ax.set_xlabel(C.XLAB_S, fontsize=13)
    ax.set_ylabel(C.YLAB_RHO + r"   (cohort median)", fontsize=13)
    ax.set_xlim(s0[0] * 0.9, s0[-1] * 1.15)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # per-band annotations (place near each curve's right end)
    def _lab(band, text, dy=0.0, ha="left"):
        g = df[df.band == band].sort_values("s")
        ax.text(g.s.values[-1] * 1.02, g.obs_med.values[-1] + dy, text,
                color=C.band_color(band), fontsize=12, va="center", ha=ha,
                fontweight="bold" if band in CARRIER else "normal")
    _lab("beta", r"$\beta$  held at every scale", dy=+0.015)
    _lab("alpha", r"$\alpha$  mesoscale", dy=-0.02)
    # theta null callout near its mid-scale
    gt = df[df.band == "theta"].sort_values("s")
    ax.annotate(r"$\theta$  null (control)", xy=(gt.s.values[7], gt.obs_med.values[7]),
                xytext=(gt.s.values[7], gt.obs_med.values[7] - 0.10), color=C.band_color("theta"),
                fontsize=11, ha="center",
                arrowprops=dict(arrowstyle="-", color=C.band_color("theta"), lw=1.0))

    # legend (band swatches, ordered slow->fast)
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color=C.band_color(b), lw=3,
                      label=C.BTeX.get(b, b)) for b in BANDS]
    ax.legend(handles=handles, loc="upper center", ncol=6, frameon=False,
              fontsize=11, bbox_to_anchor=(0.5, 1.08), columnspacing=1.3, handlelength=1.4)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight", dpi=200)
    plt.close(fig)
    for band in BANDS:
        g = df[df.band == band]
        print(f"  {band:11s} sig {int((g.gate_p < 0.05).sum())}/16   "
              f"obs_med range [{g.obs_med.min():+.2f}, {g.obs_med.max():+.2f}]")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
