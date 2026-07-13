#!/usr/bin/env python3
"""fig_tau_role -- the "role of tau" money figure for the sparsified §1 rebuild.

Two stacked stories, one figure:
  (top)  sparsification switches the multiscale ladder ON. Left: specific heat C(tau)
         of one exemplar network, dense (single collapse peak) vs mst@0.20 (a peak
         ladder). Right: interior-peak count vs backbone density across sparsifiers
         -- the ladder turns on below density ~= 0.2.
  (bottom) the trace is tau-resolved. Cohort rho_sym^coph(s) per band with the
         matched-strength null envelope and the scales that clear the cohort gate,
         so beta reads scale-broad, alpha mesoscale, the rest null.

Reads (CSV only):
  sparsification_recovery/multiscale.csv   (interior-peak count vs density)
  ms_mst020/{cohort_gate.csv, per_patient_scale.csv}
and recomputes one exemplar C(tau) fresh from FC via the shared backbone layer.
"""
from __future__ import annotations
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import _common as C
from lrg_eegfc.utils.fc.heat_multiscale import entropy_specific_heat, specific_heat_peaks

EXEMPLAR = ("Pat_02", "beta")     # clean multi-peak ladder on the flagship band
DENSE, MST = "0.35", C.band_color("beta")   # dense grey vs mst@0.20 accent


def ctau_panel(ax):
    """C(tau) dense vs mst@0.20 for one exemplar -- single peak -> ladder."""
    pat, band = EXEMPLAR
    W = C.load_phase(pat, "task_test", band)
    for lab, frac, col, lw in [("dense", 1.0, DENSE, 1.8), ("mst@0.20", C.FRAC, MST, 2.4)]:
        ev, _ = C.eig_backbone(W, frac)
        res = entropy_specific_heat(ev)
        pk = specific_heat_peaks(res)
        ax.plot(res["s"], res["C"], "-", color=col, lw=lw,
                label=f"{lab}  ({pk['n_peaks']} peak{'s' if pk['n_peaks'] != 1 else ''})")
        ax.plot(pk["s_peaks"], np.interp(pk["s_peaks"], res["s"], res["C"]),
                "v", color=col, ms=7, mec="none")
    ax.axvline(1.0, color="0.8", lw=0.7, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel(C.XLAB_S)
    ax.set_ylabel(r"specific heat $C(\tau)$")
    ax.set_title(rf"sparsification turns on multiscale  ({C.BTeX[band]}, {pat.replace('Pat_','P')})",
                 fontsize=11, fontweight="bold", loc="left")
    ax.legend(frameon=False, fontsize=9, loc="upper right")


def onset_panel(ax):
    """Interior-peak count vs backbone density across sparsifiers (cohort mean)."""
    ms = pd.read_csv(C.RECOV / "multiscale.csv")
    g = ms.groupby("method").agg(density=("density", "mean"),
                                 n_peaks=("n_peaks", "mean")).sort_values("density")
    ax.axvspan(0.0, 0.20, color=MST, alpha=0.08, lw=0)
    ax.plot(g.density, g.n_peaks, "-o", color="0.25", lw=1.6, ms=5, zorder=3)
    # flag the chosen backbone
    if "mst0.20" in g.index:
        r = g.loc["mst0.20"]
        ax.plot(r.density, r.n_peaks, "o", ms=11, mec=MST, mfc="none", mew=2.2, zorder=4)
        ax.annotate("mst@0.20", (r.density, r.n_peaks), textcoords="offset points",
                    xytext=(6, 8), color=MST, fontsize=9, fontweight="bold")
    ax.axhline(1.0, color="0.7", lw=0.7, ls="--")
    ax.axvline(0.20, color=MST, lw=0.8, ls=":")
    ax.set_xlabel("backbone density")
    ax.set_ylabel(r"interior $C(\tau)$ peaks")
    ax.set_title("multiscale ladder turns on below density $\\approx 0.2$",
                 fontsize=11, fontweight="bold", loc="left")


def _classify(fire_s, clears_report):
    """Honest shape label keyed to the actual verdict gate.

    The claim is set by the cohort gate at the reporting scale s~=5.6 (beta 0.001,
    alpha 0.007 clear; delta/high_gamma/theta/low_gamma fail). A CLAIMED band is
    then sub-labelled by breadth (scale-broad vs mesoscale); a band that does not
    clear is patchy / coarse-only (collapse tail) / null -- never over-claimed.
    """
    n = len(fire_s)
    coarse_only = n > 0 and np.min(fire_s) >= 40.0
    if clears_report:
        return "scale-broad" if n >= 14 else "mesoscale"
    if coarse_only:
        return "coarse-only"
    return "patchy" if n >= 3 else "null"


def band_panel(ax, band, gate, pp):
    col = C.band_color(band)
    g = gate[gate.band == band].sort_values("s")
    p = pp[pp.band == band]
    env = p.groupby("s").agg(p50=("surr_p50", "median"), p95=("surr_p95", "median"))
    obs = p.groupby("s").agg(q1=("obs_rho", lambda x: x.quantile(0.25)),
                             q3=("obs_rho", lambda x: x.quantile(0.75)))
    ax.fill_between(env.index, env.p50, env.p95, color="0.6", alpha=0.30, lw=0,
                    label="null $p_{50}$–$p_{95}$")
    ax.fill_between(obs.index, obs.q1, obs.q3, color=col, alpha=0.16, lw=0)
    ax.plot(g.s, g.obs_med, "-o", color=col, lw=2.3, ms=4, label="obs median")
    sig = g[g.gate_p < 0.05]
    n_clear = len(sig)
    if not sig.empty:
        yb = g.obs_med.max() * 1.10 + 0.01
        ax.plot(sig.s, np.full(len(sig), yb), "v", color=col, ms=6, mec="none")
    ax.axhline(0, color="0.75", lw=0.6)
    ax.axvline(1.0, color="0.85", lw=0.7, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel(C.XLAB_S)
    ax.set_ylabel(C.YLAB_RHO)
    gr = g.iloc[(g.s - C.S_REPORT).abs().argmin()]
    ax.set_title(rf"{C.BTeX[band]}  [{_classify(sig.s.values, gr.gate_p < 0.05)}]  {n_clear}/16  "
                 rf"$p_{{\min}}$={g.gate_p.min():.3f}",
                 color=col, fontsize=11, fontweight="bold", loc="left")


def main():
    C.use_lrg_style()
    gate = pd.read_csv(C.MS / "cohort_gate.csv")
    pp = pd.read_csv(C.MS / "per_patient_scale.csv")

    fig = plt.figure(figsize=(15, 12))
    outer = GridSpec(2, 1, height_ratios=[1.0, 1.9], hspace=0.30, figure=fig)
    top = outer[0].subgridspec(1, 2, wspace=0.24)
    ctau_panel(fig.add_subplot(top[0, 0]))
    onset_panel(fig.add_subplot(top[0, 1]))
    grid = outer[1].subgridspec(2, 3, wspace=0.30, hspace=0.42)
    band_axes = [fig.add_subplot(grid[k // 3, k % 3]) for k in range(len(C.BANDS))]
    for ax, band in zip(band_axes, C.BANDS):
        band_panel(ax, band, gate, pp)
    # single shared legend for the band grid, figure-level
    h, l = band_axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=len(h), frameon=False, fontsize=10)

    C.FIGDIR.mkdir(parents=True, exist_ok=True)
    out = C.FIGDIR / "fig_tau_role.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {out}")


if __name__ == "__main__":
    main()
