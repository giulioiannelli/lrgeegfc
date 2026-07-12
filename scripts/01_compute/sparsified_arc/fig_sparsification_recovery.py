#!/usr/bin/env python3
"""Figure: rho_sym(s) trace vs diffusion scale, across sparsifiers, per band.

The recovery money plot. One panel per band; three sparsifiers overlaid
(dense=old scheme, mst0.20=minimal-modification density backbone, perc=parameter-
free percolation). Cohort median +/- IQR over 10 patients. Shows which bands carry
a trace that is INVARIANT to the sparsification choice (beta: flat-positive
everywhere) vs scale-tuned (alpha: a peak that needs the strong-edge backbone) vs
collapse artifacts (delta: coarse-only on the dense graph). Reads per-scale, never
a scalar (feedback_report_tau_dependence_no_scalar_collapse).
"""
from __future__ import annotations
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

OUT = ROOT / "data" / "sparsified_arc" / "sparsification_recovery"
FIG = ROOT / "data" / "sparsified_arc" / "figures" / "sparsification_recovery"
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TeX = {"delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
            "beta": r"$\beta$", "low_gamma": r"$\gamma_{\mathrm{low}}$",
            "high_gamma": r"$\gamma_{\mathrm{high}}$"}
METHODS = [("mst1.00", "dense (old)", "0.15", "-"),
           ("mst0.20", "mst@0.20", "#c1121f", "-"),
           ("perc", "percolation", "#0353a4", "-")]


def med_iqr(x, col="rho_sym"):
    g = x.groupby("s")[col]
    return g.median(), g.quantile(0.25), g.quantile(0.75)


def main():
    use_lrg_style()
    df = pd.read_csv(OUT / "rho_sym_sweep.csv")
    FIG.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), sharex=True)
    for ax, band in zip(axes.ravel(), BANDS):
        for meth, lab, col, ls in METHODS:
            x = df[(df.method == meth) & (df.band == band)]
            if x.empty:
                continue
            m, lo, hi = med_iqr(x)
            s = m.index.values
            ax.fill_between(s, lo.values, hi.values, color=col, alpha=0.10, lw=0)
            ax.plot(s, m.values, ls, color=col, lw=2.2, label=lab)
        ax.axhline(0, color="0.6", lw=0.7)
        ax.axvline(1.0, color="0.8", lw=0.7, ls=":")
        ax.set_xscale("log")
        ax.set_title(BAND_TeX[band], color=band_color(band), fontweight="bold",
                     fontsize=15, loc="left")
        ax.set_xlabel(r"diffusion scale  $s=\tau\lambda_{\max}$")
        ax.set_ylabel(r"$\rho_{\mathrm{sym}}^{\mathrm{coph}}$ (trace)")
        ax.set_ylim(-0.12, 0.34)
    axes.ravel()[0].legend(frameon=False, fontsize=10, loc="upper right")
    fig.tight_layout()
    out = FIG / "rho_sym_vs_scale_by_sparsifier.pdf"
    fig.savefig(out, transparent=True); plt.close(fig)
    print(f"[fig] {out}")


if __name__ == "__main__":
    main()
