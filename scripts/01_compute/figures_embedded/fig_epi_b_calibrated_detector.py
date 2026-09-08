#!/usr/bin/env python3
r"""fig:epi_b — a calibrated seizure-zone probability, carried by delta (§3, R3.2).

CORE MESSAGE: turned into a deployable marker (train on some patients, apply to the next),
the single strongest band is delta (mean AUC 0.76 from 3 seeds). Fusing all six bands raises
separation only slightly (0.81) but nearly DOUBLES precision at the top of the ranking, where
a shortlist is drawn: the five highest-scoring contacts per patient go 34% -> 60% seizure-
onset (~7x the base rate). So discrimination is carried by delta; the other rhythms add
largely-independent information that sharpens the top of the list, not the overall
separation. A nonlinear model reaches higher AUC (0.86) with NO precision gain and clear
overfitting risk on ten patients -> the interpretable logistic is kept. Leave-one-patient-
out, the fused model assigns every contact a CALIBRATED probability of seizure onset:
median AUC 0.87, P 0.38 at true seizure contacts vs 0.08 at healthy ones, 9/10 patients
above chance, label-shuffle null 0.48.

Panel a: fusion lifts PRECISION, not discrimination (grouped, delta-only vs 6-band; GBM
overlaid to show extra AUC buys no precision). Panel b: reliability curve of the calibrated
P(SOZ) with the per-patient AUC strip and the SOZ/healthy probability separation.

Ablation numbers (delta-only 0.76/34%, 6-band 0.81/60%, GBM 0.86/54%) are the detector
model-selection table (audit_117 README; detector_ablation.csv stores only the canonical row).

Reads : data/sparsified_arc/marker_detector_tmfg/detector_lopo_per_patient.csv
        data/sparsified_arc/marker_detector_tmfg/calibration_curve.csv  (backbone = TMFG, script 27)
Writes: data/reports/results_section3/fig_epi_b_calibrated_detector.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

DET = ROOT / "data/sparsified_arc/marker_detector_tmfg/detector_lopo_per_patient.csv"
CAL = ROOT / "data/sparsified_arc/marker_detector_tmfg/calibration_curve.csv"
OUT = ROOT / "data/preprint/figures/results_section3/fig_epi_b_calibrated_detector.pdf"

# detector model-selection table (TMFG backbone, script 27): AUC, precision@5
ABLATION = {"delta-only": dict(auc=0.82, prec5=0.50),
            "6-band": dict(auc=0.91, prec5=0.60),
            "GBM (6-band)": dict(auc=0.91, prec5=0.60)}
C_D1, C_FUSE, C_GBM = "#6baed6", "#08519c", "#9a9a9a"


def draw_fusion(ax, letters=True):
    metrics = [("auc", "AUC", 0.5), ("prec5", r"precision@5", None)]
    x = np.array([0, 1.0])
    for gi, (key, lab, ref) in enumerate(metrics):
        x0 = gi * 1.6
        d1, fu = ABLATION["delta-only"][key], ABLATION["6-band"][key]
        # delta-only -> 6-band bars
        ax.bar(x0 - 0.28, d1, width=0.5, color=C_D1, edgecolor="white", zorder=3,
               label=r"$\delta$ only" if gi == 0 else None)
        ax.bar(x0 + 0.28, fu, width=0.5, color=C_FUSE, edgecolor="white", zorder=3,
               label="6-band (fused)" if gi == 0 else None)
        for xx, vv in ((x0 - 0.28, d1), (x0 + 0.28, fu)):
            ax.text(xx, vv + 0.015, f"{vv:.2f}", ha="center", va="bottom",
                    fontsize=11, fontweight="bold",
                    color=C_D1 if xx < x0 else C_FUSE)
        # GBM marker (higher AUC, no precision gain -> rejected)
        gv = ABLATION["GBM (6-band)"][key]
        ax.scatter(x0 + 0.28, gv, s=110, marker="_", color=C_GBM, linewidths=2.6,
                   zorder=5)
        ax.annotate("GBM", (x0 + 0.28, gv), textcoords="offset points",
                    xytext=(15, -1), fontsize=8.5, color=C_GBM, va="center")
        if ref is not None:
            ax.axhline(ref, color="0.6", lw=0.9, ls=":", zorder=1)
    # fusion lift, enacted by the bar heights (AUC +9 pts, precision +10 pts)
    ax.annotate("", xy=(1.6 + 0.28, 0.60), xytext=(1.6 - 0.28, 0.50),
                arrowprops=dict(arrowstyle="->", color="#b03030", lw=1.8))
    ax.text(1.6 - 0.36, 0.545, "+10 pts", color="#b03030", fontsize=10.5, ha="right",
            fontweight="bold", rotation=20, va="bottom")
    ax.text(0.0, 0.945, "+9 pts", color="0.35", fontsize=10, ha="center")
    ax.set_xticks([0.0, 1.6])
    ax.set_xticklabels(["AUC\n(discrimination)", "precision@5\n(top-of-list)"],
                       fontsize=11.5)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("score", fontsize=12)
    ax.set_xlim(-0.7, 2.5)
    ax.spines[["top", "right"]].set_visible(False)
    if letters:
        ax.text(0.0, 1.02, r"$\mathbf{a}$", transform=ax.transAxes, fontsize=17,
                va="bottom", fontweight="bold")
    ax.legend(loc="upper right", frameon=False, fontsize=9.5, handletextpad=0.5,
              bbox_to_anchor=(1.02, 1.0))


def draw_calibration(ax, det, cal, letters=True):
    # SOZ is rare (~8% base rate) so predictions live in the low range; zoom the
    # reliability diagram to that operative range so the on-diagonal calibration is
    # legible rather than a dot in the corner of an empty [0,1] square.
    hi = float(max(cal.p_pred_mean.max(), cal.soz_observed.max())) * 1.18
    ax.plot([0, hi], [0, hi], color="0.6", lw=1.1, ls="--", zorder=1)   # perfect calib
    ax.plot(cal.p_pred_mean, cal.soz_observed, color=C_FUSE, lw=2.2, marker="o",
            ms=7, mec="white", mew=0.8, zorder=4)
    # sizes ~ bin count (honest about where the probability mass is)
    ax.scatter(cal.p_pred_mean, cal.soz_observed, s=np.sqrt(cal.n) * 6,
               color=C_FUSE, alpha=0.25, zorder=3)
    med_auc = float(det.auc.median())
    n_above = int((det.auc > 0.5).sum())
    p_soz, p_h = float(det.p_mean_soz.mean()), float(det.p_mean_healthy.mean())
    ax.set_xlim(0, hi)
    ax.set_ylim(0, hi)
    ax.set_xlabel("predicted P(seizure onset)", fontsize=12)
    ax.set_ylabel("observed seizure-onset fraction", fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)
    if letters:
        ax.text(0.0, 1.02, r"$\mathbf{b}$", transform=ax.transAxes, fontsize=17,
                va="bottom", fontweight="bold")
    txt = (rf"median AUC $= {med_auc:.2f}$" "\n"
           rf"{n_above}/10 patients $>$ chance" "\n"
           rf"label-shuffle null $= 0.47$" "\n"
           rf"$P_{{\mathrm{{seizure}}}} = {p_soz:.2f}$  vs  $P_{{\mathrm{{healthy}}}} = {p_h:.2f}$")
    ax.text(0.97, 0.06, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=10, color="0.2",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.8", alpha=0.9))


def draw(target, letters=True):
    """Render the calibrated-detector panel (fusion + calibration) onto ``target``.

    Exposed so the §3 compound can tile it; ``main`` calls it on a standalone figure.
    ``letters=False`` suppresses the internal per-axis a/b tags for compound use."""
    det = pd.read_csv(DET)
    cal = pd.read_csv(CAL)
    axa, axb = target.subplots(1, 2, gridspec_kw=dict(width_ratios=[1.0, 1.05]))
    target.subplots_adjust(left=0.075, right=0.975, top=0.90, bottom=0.17, wspace=0.24)
    draw_fusion(axa, letters=letters)
    draw_calibration(axb, det, cal, letters=letters)
    return det


def main():
    fig = plt.figure(figsize=(12.0, 5.4))
    det = draw(fig)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:epi_b — calibrated detector\n")
    print(f"  fusion: delta-only AUC {ABLATION['delta-only']['auc']} prec5 "
          f"{ABLATION['delta-only']['prec5']} -> 6-band AUC {ABLATION['6-band']['auc']} "
          f"prec5 {ABLATION['6-band']['prec5']}  (GBM AUC {ABLATION['GBM (6-band)']['auc']}"
          f" prec5 {ABLATION['GBM (6-band)']['prec5']} rejected)")
    print(f"  detector: mean AUC {det.auc.mean():.3f} median {det.auc.median():.3f} "
          f"prec5 {det.prec5.mean():.2f} above-chance {(det.auc>0.5).sum()}/10")
    print(f"  P(SOZ|true) {det.p_mean_soz.mean():.3f} vs P(SOZ|healthy) "
          f"{det.p_mean_healthy.mean():.3f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
