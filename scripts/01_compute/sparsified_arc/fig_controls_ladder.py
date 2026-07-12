#!/usr/bin/env python3
"""Controls ladder figure: descriptor x band gate significance, per functional.

2x2 heatmaps (T_test, T_learn, T_infspec, T_infspec|e). Rows = descriptors in
structural order (raw pairwise -> node scalars -> topological -> spectral -> LRG
cophenetic single/multi), cols = bands, cell = -log10(cohort gate_p) capped at 3;
p<0.05 cells ringed. Shows band-selectivity (alpha/beta, others silent) is unique to
the LRG cophenetic rungs, while raw-FC fires broadly and spectral/topological fail.
Reads data/sparsified_arc/controls_ladder/cohort_gate.csv.
"""
from __future__ import annotations
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.visuals.styles import use_lrg_style

OUT = ROOT / "data" / "sparsified_arc" / "controls_ladder"
FIG = ROOT / "data" / "sparsified_arc" / "figures" / "controls_ladder"
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BL = [r"$\delta$", r"$\theta$", r"$\alpha$", r"$\beta$", r"$\gamma_l$", r"$\gamma_h$"]
DESCS = ["raw_fc", "strength", "clustering", "geodesic", "resistance",
         "coph_taumin", "coph_meso"]
DL = ["raw FC", "strength", "clustering", "geodesic", "resistance (spectral)",
      "coph @ $\\tau_{\\min}$", "coph @ meso"]
FUNCS = [("T_test", "whole-task trace"), ("T_learn", "encoding"),
         ("T_infspec", "inference-specific"), ("T_infspec_pe", "inference | encoding")]


def main():
    use_lrg_style()
    g = pd.read_csv(OUT / "cohort_gate.csv")
    FIG.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    for ax, (f, lab) in zip(axes.ravel(), FUNCS):
        M = np.full((len(DESCS), len(BANDS)), np.nan)
        for i, d in enumerate(DESCS):
            for j, b in enumerate(BANDS):
                r = g[(g.functional == f) & (g.descriptor == d) & (g.band == b)]
                if not r.empty:
                    M[i, j] = -np.log10(max(float(r.gate_p.iloc[0]), 1e-4))
        im = ax.imshow(M, aspect="auto", cmap="magma", vmin=0, vmax=3.0)
        for i in range(len(DESCS)):
            for j in range(len(BANDS)):
                if np.isfinite(M[i, j]) and M[i, j] >= -np.log10(0.05):
                    ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                           edgecolor="cyan", lw=2.2))
        ax.set_xticks(range(len(BANDS))); ax.set_xticklabels(BL, fontsize=12)
        ax.set_yticks(range(len(DESCS))); ax.set_yticklabels(DL, fontsize=10)
        ax.set_title(f"{f} — {lab}", fontweight="bold", fontsize=12, loc="left")
        ax.axhline(4.5, color="white", lw=1.5)  # separate simple vs LRG cophenetic
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(r"$-\log_{10}$ gate $p$", fontsize=9)
    fig.tight_layout()
    out = FIG / "controls_ladder_gate.pdf"
    fig.savefig(out, transparent=True); plt.close(fig); print(f"[fig] {out}")


if __name__ == "__main__":
    main()
