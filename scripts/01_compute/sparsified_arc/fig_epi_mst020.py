#!/usr/bin/env python3
"""Epilepsy (SOZ) marker on mst@0.20: AUC(s) tau-sweep + multiscale-vs-single summary.

LEFT: cohort-median ROC-AUC vs diffusion scale s (per band, interpolated to a common
log-s grid), dotted = single-scale baseline s=10 (heat_t5). RIGHT: multiscale AUC per
band (bar) with n/10 beating the strength-matched fake-SOZ null annotated; the LRG
propagator is a SOZ detector strongest in delta / low_gamma (and beta). Reads
data/sparsified_arc/epi_arc_mst020/.
"""
from __future__ import annotations
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

OUT = ROOT / "data" / "sparsified_arc" / "epi_arc_mst020"
FIG = ROOT / "data" / "sparsified_arc" / "figures" / "epi_arc_mst020"
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BTeX = {"delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$", "beta": r"$\beta$",
        "low_gamma": r"$\gamma_{\rm low}$", "high_gamma": r"$\gamma_{\rm high}$"}
GS = np.logspace(0, np.log10(60), 40)


def main():
    use_lrg_style()
    df = pd.read_csv(OUT / "per_cell.csv")
    FIG.mkdir(parents=True, exist_ok=True)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14, 5.6))

    for b in BANDS:
        curves = []
        for pat in df[df.band == b].patient:
            f = OUT / b / f"{pat}.npz"
            if not f.exists():
                continue
            d = np.load(f)
            curves.append(np.interp(GS, d["s"], d["auc"]))
        if not curves:
            continue
        C = np.array(curves)
        med = np.nanmedian(C, 0)
        axL.plot(GS, med, "-", color=band_color(b), lw=2.2, label=BTeX[b])
    axL.axhline(0.5, color="0.6", lw=0.7)
    axL.axvline(10, color="0.75", lw=0.8, ls=":")
    axL.set_xscale("log")
    axL.set_xlabel(r"diffusion scale $s=\tau\lambda_{\max}$")
    axL.set_ylabel("SOZ ROC-AUC (cohort median)")
    axL.set_title("propagator SOZ marker vs scale", fontweight="bold", fontsize=12, loc="left")
    axL.legend(frameon=False, fontsize=9, ncol=2)

    x = np.arange(len(BANDS))
    med = [df[df.band == b].auc_multi.median() for b in BANDS]
    beat = [int((df[df.band == b].p_null_multi < 0.05).sum()) for b in BANDS]
    cols = [band_color(b) for b in BANDS]
    axR.bar(x, med, color=cols, alpha=0.85)
    for xi, (m, nb) in enumerate(zip(med, beat)):
        axR.text(xi, m + 0.01, f"{nb}/10", ha="center", fontsize=9, fontweight="bold")
    axR.axhline(0.5, color="0.6", lw=0.7)
    axR.set_xticks(x); axR.set_xticklabels([BTeX[b] for b in BANDS])
    axR.set_ylabel("multiscale SOZ ROC-AUC (median)")
    axR.set_ylim(0.4, 0.9)
    axR.set_title("multiscale marker + n/10 beat matched fake-SOZ null",
                  fontweight="bold", fontsize=11, loc="left")
    fig.tight_layout()
    out = FIG / "epi_marker_mst020.pdf"
    fig.savefig(out, transparent=True); plt.close(fig); print(f"[fig] {out}")


if __name__ == "__main__":
    main()
