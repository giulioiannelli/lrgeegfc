#!/usr/bin/env python3
"""D1 figures -- entropy Shat(tau) & specific heat C(tau) on the percolation backbone.

  1. entropy_summary.pdf -- cohort-median Shat(s) and C(s) per band (one clean
     panel per band) + the fraction of cells with a multiscale (>1) C-peak.
  2. entropy_curves_{band}.pdf -- control grid, one panel per (patient x phase):
     Shat(s) (0->1, left axis) and C(s) (right axis) vs s=tau*lambda_max; s=1
     (=tau_min) marked; C-peaks dotted. Shows S spans the full 0->1 range.
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

ENT = ROOT / "data" / "sparsified_arc" / "entropy"
FIGS = ROOT / "data" / "sparsified_arc" / "figures" / "entropy"
PHASES = ("A", "B", "task_test", "rest_post")


def _load(pat, phase, band):
    f = ENT / band / f"{pat}_{phase}.npz"
    return np.load(f) if f.exists() else None


def _interp_median(band, key, s_common):
    vals = []
    for pat in COHORT:
        for ph in PHASES:
            z = _load(pat, ph, band)
            if z is None:
                continue
            vals.append(np.interp(s_common, z["s"], z[key], left=np.nan, right=np.nan))
    return np.nanmedian(np.array(vals), axis=0) if vals else np.full_like(s_common, np.nan)


def summary_fig(df):
    use_lrg_style()
    s_common = np.logspace(0, np.log10(200), 200)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    for ax, band in zip(axes.ravel(), BANDS):
        col = band_color(band)
        S = _interp_median(band, "S", s_common)
        C = _interp_median(band, "C", s_common)
        ax.plot(s_common, S, "-", color=col, lw=2.2, label=r"$\hat S(\tau)$")
        ax.set_ylim(-0.02, 1.05); ax.set_xscale("log")
        ax.axvline(1.0, color="0.6", lw=0.8, ls=":")
        ax2 = ax.twinx()
        ax2.plot(s_common, C, "--", color="0.4", lw=1.6)
        ax2.set_ylabel(r"$C(\tau)$ (grey)", color="0.4", fontsize=9)
        ax2.tick_params(axis="y", labelcolor="0.4")
        frac_multi = (df[df.band == band].n_peaks > 1).mean() * 100
        ax.set_title(f"{band}   —   {frac_multi:.0f}% cells multiscale (>1 C-peak)",
                     fontweight="bold", fontsize=10, color=col)
        ax.set_xlabel(r"scale $s=\tau\,\lambda_{\max}$")
        ax.set_ylabel(r"$\hat S(\tau)=S/\ln N$")
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    out = FIGS / "entropy_summary.pdf"
    fig.savefig(out, transparent=True); plt.close(fig); print(f"[fig] {out}")


def per_band_grids():
    use_lrg_style()
    for band in BANDS:
        fig, axes = plt.subplots(len(COHORT), len(PHASES),
                                 figsize=(4 * len(PHASES), 2.1 * len(COHORT)),
                                 sharex=True)
        col = band_color(band)
        for i, pat in enumerate(COHORT):
            for j, ph in enumerate(PHASES):
                ax = axes[i, j]
                z = _load(pat, ph, band)
                if z is None:
                    ax.set_axis_off(); continue
                ax.plot(z["s"], z["S"], "-", color=col, lw=1.5)
                ax.set_ylim(-0.02, 1.05); ax.set_xscale("log")
                ax.axvline(1.0, color="0.6", lw=0.7, ls=":")
                ax2 = ax.twinx()
                ax2.plot(z["s"], z["C"], "--", color="0.45", lw=1.0)
                ax2.set_yticks([])
                for sp in np.atleast_1d(z["s_peaks"]):
                    ax2.axvline(float(sp), color="#c0392b", lw=0.6, alpha=0.6)
                if i == 0:
                    ax.set_title(ph, fontsize=9)
                if j == 0:
                    ax.set_ylabel(pat.replace("Pat_", "P"), fontsize=8)
        fig.supxlabel(r"scale $s=\tau\lambda_{\max}$   —   $\hat S$ colour (0→1), $C$ grey, peaks red, $s{=}1$ dotted",
                      fontsize=9)
        fig.tight_layout()
        out = FIGS / f"entropy_curves_{band}.pdf"
        fig.savefig(out, transparent=True); plt.close(fig); print(f"[fig] {out}")


def main():
    df = pd.read_csv(ENT / "summary.csv")
    summary_fig(df)
    per_band_grids()
    print(f"[D1-figs] -> {FIGS}")


if __name__ == "__main__":
    main()
