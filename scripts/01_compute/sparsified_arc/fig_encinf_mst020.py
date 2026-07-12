#!/usr/bin/env python3
"""Encoding-vs-inference on mst@0.20: per-scale matched-strength gate, 4 functionals.

Reads data/sparsified_arc/enc_inf_arc_mst020/cohort_gate_vs_s.csv. 2x2 panels, one
per functional (T_test whole-task, T_learn encoding echo, T_infspec inference-specific,
T_infspec_pe inference-specific|encoding). Each: -log10(cohort gate_p) vs diffusion
scale s, one band-coloured line per band; dashed line = p=0.05. Reads per-scale, never
scale-max (the scale-max forking-path inflates every band).
"""
from __future__ import annotations
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

OUT = ROOT / "data" / "sparsified_arc" / "enc_inf_arc_mst020"
FIG = ROOT / "data" / "sparsified_arc" / "figures" / "enc_inf_arc_mst020"
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BTeX = {"delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$", "beta": r"$\beta$",
        "low_gamma": r"$\gamma_{\rm low}$", "high_gamma": r"$\gamma_{\rm high}$"}
FUNCS = [("T_test", "whole-task trace"), ("T_learn", "encoding echo"),
         ("T_infspec", "inference-specific"), ("T_infspec_pe", "inference-specific | encoding")]


def main():
    use_lrg_style()
    g = pd.read_csv(OUT / "cohort_gate_vs_s.csv")
    FIG.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    for ax, (f, lab) in zip(axes.ravel(), FUNCS):
        for b in BANDS:
            x = g[(g.functional == f) & (g.band == b)].sort_values("s")
            if x.empty:
                continue
            y = -np.log10(np.clip(x.gate_p.values, 1e-4, 1.0))
            ax.plot(x.s.values, y, "-o", color=band_color(b), lw=2, ms=3.5, label=BTeX[b])
        ax.axhline(-np.log10(0.05), color="0.5", ls="--", lw=1)
        ax.set_xscale("log")
        ax.set_title(f"{f} — {lab}", fontweight="bold", fontsize=12, loc="left")
        ax.set_xlabel(r"$s=\tau\lambda_{\max}$")
        ax.set_ylabel(r"$-\log_{10}$ gate $p$")
        ax.set_ylim(0, 3.4)
    axes.ravel()[0].legend(frameon=False, fontsize=9, ncol=2, loc="upper right")
    fig.tight_layout()
    out = FIG / "encinf_gate_vs_scale.pdf"
    fig.savefig(out, transparent=True); plt.close(fig); print(f"[fig] {out}")


if __name__ == "__main__":
    main()
