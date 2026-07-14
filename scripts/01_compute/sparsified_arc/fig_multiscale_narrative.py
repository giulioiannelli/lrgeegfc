#!/usr/bin/env python3
"""Two minimal talk figures for "the multiscale nature of inference".

FIG 1  scale-signature panel: cohort median concordance vs diffusion scale s, for
       trace (T_test) / encoding (T_learn) / inference-specific (T_infspec_pe), bands
       beta/alpha/delta. Filled marker = cohort gate p<.05 at that scale. Vertical line
       at s=1 (the local / tau_min scale). Shows: encoding fires already at the local
       scale; inference is silent locally and EMERGES at the mesoscale.
FIG 2  method ladder (beta): matched-strength gate p for each descriptor x component.
       Shows the cophenetic hierarchy catches beta encoding where raw/clustering/geodesic/
       resistance/strength do not — the multiscale read enriches the single-scale views.

Sources (cached, no recompute):
  data/sparsified_arc/enc_inf_arc_mst020/cohort_gate_vs_s.csv
  data/sparsified_arc/controls_ladder/cohort_gate.csv
PDF only, vector, transparent.
"""
from __future__ import annotations
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

ENC = ROOT / "data/sparsified_arc/enc_inf_arc_mst020/cohort_gate_vs_s.csv"
LAD = ROOT / "data/sparsified_arc/controls_ladder/cohort_gate.csv"
OUT = ROOT / "data/sparsified_arc/figures"


def fig_scale_signatures():
    d = pd.read_csv(ENC)
    panels = [("T_test", "trace"), ("T_learn", "encoding"), ("T_infspec_pe", "inference-specific")]
    bands = ["beta", "alpha", "delta"]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.1), sharey=True)
    for ax, (func, title) in zip(axes, panels):
        for b in bands:
            x = d[(d.functional == func) & (d.band == b)].sort_values("s")
            if x.empty:
                continue
            s, y, p = x.s.values, x.obs_median.values, x.gate_p.values
            col = band_color(b)
            ax.plot(s, y, "-", color=col, lw=1.6, alpha=0.9, zorder=2)
            sig = p < 0.05
            ax.scatter(s[sig], y[sig], s=34, color=col, zorder=3, label=b)
            ax.scatter(s[~sig], y[~sig], s=26, facecolor="white", edgecolor=col,
                       lw=1.3, zorder=3)
        ax.axhline(0, color="0.7", lw=0.7)
        ax.axvline(1.0, color="0.55", lw=0.9, ls=":")
        ax.set_xscale("log")
        ax.set_title(title, fontsize=11, fontweight="bold", loc="left")
        ax.set_xlabel(r"diffusion scale  $s=\tau\lambda_{\max}$")
    axes[0].set_ylabel(r"cohort median  $\rho_{\rm sym}$")
    axes[0].annotate("local\n$\\tau_{\\min}$", (1.0, axes[0].get_ylim()[1]),
                     fontsize=7.5, color="0.4", ha="center", va="top")
    h = [plt.Line2D([], [], color=band_color(b), lw=2, marker="o", ls="-") for b in bands]
    fig.legend(h, bands, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    p = OUT / "fig_scale_signatures.pdf"
    fig.savefig(p, transparent=True); plt.close(fig)
    return p


def fig_method_ladder():
    d = pd.read_csv(LAD)
    order = ["coph_meso", "coph_taumin", "geodesic", "clustering", "raw_fc",
             "resistance", "strength"]
    labels = {"coph_meso": "cophenetic (meso)", "coph_taumin": r"cophenetic ($\tau_{\min}$)",
              "geodesic": "geodesic", "clustering": "clustering", "raw_fc": "raw pairwise FC",
              "resistance": "spectral resistance", "strength": "node strength"}
    cols = [("T_test", "trace"), ("T_learn", "encoding"), ("T_infspec_pe", "inference")]
    b = d[d.band == "beta"]
    M = np.full((len(order), len(cols)), np.nan)
    P = np.full_like(M, np.nan)
    for i, desc in enumerate(order):
        for j, (func, _) in enumerate(cols):
            r = b[(b.descriptor == desc) & (b.functional == func)]
            if not r.empty:
                P[i, j] = float(r.gate_p.iloc[0])
                M[i, j] = -np.log10(max(float(r.gate_p.iloc[0]), 1e-3))
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    im = ax.imshow(M, cmap="YlGnBu", vmin=0, vmax=3, aspect="auto")
    for i in range(len(order)):
        for j in range(len(cols)):
            if np.isnan(P[i, j]):
                continue
            txt = f"{P[i,j]:.3f}"
            sig = P[i, j] < 0.05
            ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                    color="white" if M[i, j] > 1.6 else "0.15",
                    fontweight="bold" if sig else "normal")
            if sig:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor="crimson", lw=2.0, zorder=4))
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([c[1] for c in cols])
    ax.set_yticks(range(len(order))); ax.set_yticklabels([labels[o] for o in order])
    ax.set_title(r"$\beta$: matched-strength gate $p$ per read-out", fontsize=10.5,
                 fontweight="bold", loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(r"$-\log_{10}p$", fontsize=9)
    ax.axhline(1.5, color="0.3", lw=1.0)   # separate cophenetic (multiscale) from the rest
    fig.tight_layout()
    p = OUT / "fig_method_ladder.pdf"
    fig.savefig(p, transparent=True); plt.close(fig)
    return p


def main():
    use_lrg_style()
    OUT.mkdir(parents=True, exist_ok=True)
    p1 = fig_scale_signatures()
    p2 = fig_method_ladder()
    print(f"[fig1 scale-signatures] -> {p1}")
    print(f"[fig2 method-ladder]    -> {p2}")


if __name__ == "__main__":
    main()
