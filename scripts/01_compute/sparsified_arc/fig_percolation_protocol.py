#!/usr/bin/env python3
"""D0 figures -- percolation protocol validation, for visual control.

Two outputs under data/sparsified_arc/percolation/figs/:
  1. percolation_stability.pdf -- the one-glance protocol verdict:
     (A) edge-fraction at theta* per band (moderate cycle-rich band shaded);
     (B) same per patient (exposes the two structured outliers Pat_14/Pat_15);
     (C) cycle rank vs edge-fraction -- every cell well above the tree line;
     (D) example P_inf/E_inf vs theta curves (moderate / near-tree / near-dense).
  2. percolation_curves_{band}.pdf -- exhaustive control grid, one panel per
     (patient x phase): P_inf (solid) & E_inf (dashed) vs theta/w_max, theta* marked.
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

PERC = ROOT / "data" / "sparsified_arc" / "percolation"
FIGS = ROOT / "data" / "sparsified_arc" / "figures" / "percolation"
PHASES = ("A", "B", "task_test", "rest_post")
MOD_LO, MOD_HI = 5.0, 35.0                     # moderate cycle-rich band (% edges)


def _load(pat, phase, band):
    f = PERC / band / f"{pat}_{phase}.npz"
    return np.load(f) if f.exists() else None


def stability_fig(df):
    use_lrg_style()
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    ef = df.edge_frac_star * 100

    # (A) per band
    ax = axes[0, 0]
    ax.axhspan(MOD_LO, MOD_HI, color="0.85", zorder=0, label="moderate cycle-rich")
    for i, b in enumerate(BANDS):
        y = ef[df.band == b].values
        x = np.full(y.size, i) + np.random.default_rng(i).uniform(-0.15, 0.15, y.size)
        ax.scatter(x, y, s=14, color=band_color(b), alpha=0.7, edgecolor="none")
        ax.plot([i - 0.28, i + 0.28], [np.median(y)] * 2, color="k", lw=2)
    ax.set_xticks(range(len(BANDS))); ax.set_xticklabels(BANDS, rotation=30, ha="right")
    ax.set_ylabel(r"edge fraction at $\theta^*$  [%]")
    ax.set_title("(A) surviving-edge fraction per band", fontweight="bold", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    # (B) per patient
    ax = axes[0, 1]
    ax.axhspan(MOD_LO, MOD_HI, color="0.85", zorder=0)
    for i, p in enumerate(COHORT):
        y = ef[df.patient == p].values
        x = np.full(y.size, i) + np.random.default_rng(i).uniform(-0.15, 0.15, y.size)
        c = "#c0392b" if p in ("Pat_14", "Pat_15") else "0.35"
        ax.scatter(x, y, s=14, color=c, alpha=0.7, edgecolor="none")
        ax.plot([i - 0.28, i + 0.28], [np.median(y)] * 2, color="k", lw=2)
    ax.set_xticks(range(len(COHORT)))
    ax.set_xticklabels([p.replace("Pat_", "") for p in COHORT])
    ax.set_xlabel("patient"); ax.set_ylabel(r"edge fraction at $\theta^*$  [%]")
    ax.set_title("(B) per patient — 14 near-tree, 15 near-dense (structured)",
                 fontweight="bold", fontsize=10)

    # (C) cycle rank vs edge fraction
    ax = axes[1, 0]
    for b in BANDS:
        x = ef[df.band == b].values
        y = df.cycle_rank[df.band == b].values
        ax.scatter(x, y, s=14, color=band_color(b), alpha=0.7, edgecolor="none", label=b)
    ax.axhline(1, color="#c0392b", lw=1.2, ls="--")
    ax.text(ax.get_xlim()[1], 1.3, "tree (cycle rank = 0)", color="#c0392b",
            ha="right", va="bottom", fontsize=8)
    ax.set_yscale("log"); ax.set_xlabel(r"edge fraction at $\theta^*$  [%]")
    ax.set_ylabel("cycle rank  (E − N + 1)")
    ax.set_title("(C) always cycle-rich, never a tree", fontweight="bold", fontsize=10)
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="lower right")

    # (D) example percolation curves
    ax = axes[1, 1]
    examples = [("Pat_05", "rest_post", "beta", "moderate (Pat_05 β)", "#2e7d32"),
                ("Pat_14", "rest_post", "beta", "near-tree (Pat_14 β)", "#c0392b"),
                ("Pat_15", "rest_post", "alpha", "near-dense (Pat_15 α)", "#1565c0")]
    for pat, ph, b, lab, col in examples:
        z = _load(pat, ph, b)
        if z is None:
            continue
        x = z["theta"] / z["w_max"]
        ax.plot(x, z["p_inf"], "-", color=col, lw=2, label=f"P∞  {lab}")
        ax.plot(x, z["e_inf"], ":", color=col, lw=1.5)
        ax.axvline(z["theta_star"] / z["w_max"], color=col, lw=0.8, alpha=0.5)
    ax.set_xlabel(r"threshold  $\theta / w_{\max}$")
    ax.set_ylabel(r"$P_\infty$ (solid) ,  $E_\infty$ (dotted)")
    ax.set_title("(D) giant-component erosion; θ* = last P∞=1",
                 fontweight="bold", fontsize=10)
    ax.legend(frameon=False, fontsize=7, loc="upper right")

    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    out = FIGS / "percolation_stability.pdf"
    fig.savefig(out, transparent=True); plt.close(fig)
    print(f"[fig] {out}")


def per_band_grids():
    use_lrg_style()
    for b in BANDS:
        fig, axes = plt.subplots(len(COHORT), len(PHASES),
                                 figsize=(4 * len(PHASES), 2.1 * len(COHORT)),
                                 sharex=True, sharey=True)
        for i, pat in enumerate(COHORT):
            for j, ph in enumerate(PHASES):
                ax = axes[i, j]
                z = _load(pat, ph, b)
                if z is None:
                    ax.set_axis_off(); continue
                x = z["theta"] / z["w_max"]
                ax.plot(x, z["p_inf"], "-", color=band_color(b), lw=1.6)
                ax.plot(x, z["e_inf"], ":", color="0.35", lw=1.2)
                ax.axvline(z["theta_star"] / z["w_max"], color="#c0392b", lw=0.9)
                ax.text(0.96, 0.9, f"{float(z['edge_frac_star'])*100:.0f}%",
                        transform=ax.transAxes, ha="right", va="top", fontsize=7,
                        color="#c0392b")
                ax.set_ylim(-0.02, 1.05)
                if i == 0:
                    ax.set_title(ph, fontsize=9)
                if j == 0:
                    ax.set_ylabel(pat.replace("Pat_", "P"), fontsize=8)
        fig.suptitle("")  # publication rule: no suptitle; band encoded in filename
        fig.supxlabel(r"threshold $\theta/w_{\max}$   —   $P_\infty$ solid, $E_\infty$ dotted, $\theta^*$ red",
                      fontsize=9)
        fig.tight_layout()
        out = FIGS / f"percolation_curves_{b}.pdf"
        fig.savefig(out, transparent=True); plt.close(fig)
        print(f"[fig] {out}")


def main():
    df = pd.read_csv(PERC / "summary.csv")
    stability_fig(df)
    per_band_grids()
    print(f"[D0-figs] -> {FIGS}")


if __name__ == "__main__":
    main()
