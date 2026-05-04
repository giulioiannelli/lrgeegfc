#!/usr/bin/env python3
"""TAM persistence figure — three panels.

Panel A: per-band cohort dotplot of mean Δ = J_post − J_pre across TAMs.
         One dot per (patient, band).
Panel B: scatter of all TAMs in (J_pre, J_post) space, coloured by band.
         Diagonal = symmetric presence; above-diagonal = trace-direction;
         high (J_pre, J_post) ≈ (0.8, 0.8) corner = anatomy-shared modules.
Panel C: histogram of TAM persistence Δ distribution per band.

Reads:  ``trace_tam_summary.csv``, ``trace_tam_per_module.csv``
Writes: ``data/reports/imcoh_mrl/figures/trace_anchor.pdf``
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT


IN_DIR = REPORTS_ROOT / "imcoh_mrl"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT

BAND_COLORS = {
    "delta":      "#1f77b4",
    "theta":      "#ff7f0e",
    "alpha":      "#2ca02c",
    "beta":       "#d62728",
    "low_gamma":  "#9467bd",
    "high_gamma": "#8c564b",
}


def main() -> None:
    sm = pd.read_csv(IN_DIR / "trace_tam_summary.csv")
    pm = pd.read_csv(IN_DIR / "trace_tam_per_module.csv")

    fig = plt.figure(figsize=(15.5, 5.0), dpi=150)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.1, 1.4])
    ax_A = fig.add_subplot(gs[0, 0])
    ax_B = fig.add_subplot(gs[0, 1])
    ax_C = fig.add_subplot(gs[0, 2])

    # ───── Panel A — per-band cohort dotplot of mean Δ ─────
    rng = np.random.default_rng(0)
    for j, band in enumerate(BANDS):
        sub = sm[sm["band"] == band]
        if sub.empty:
            continue
        vals = sub["mean_delta"].to_numpy()
        xs = j + rng.uniform(-0.15, 0.15, size=len(vals))
        cols = [BAND_COLORS[band]] * len(vals)
        ax_A.scatter(xs, vals, s=44, color=cols, alpha=0.8,
                     edgecolor="black", linewidth=0.3, zorder=4)
        med = float(np.median(vals))
        q1, q3 = float(np.quantile(vals, 0.25)), float(np.quantile(vals, 0.75))
        ax_A.vlines(j, q1, q3, color="#404040", lw=2.4, zorder=5, alpha=0.85)
        ax_A.hlines(med, j - 0.22, j + 0.22, color="#404040", lw=2.6, zorder=6)
        n_pos = int((vals > 0).sum())
        ax_A.text(j, 1.03, f"{n_pos}/{len(vals)}",
                  ha="center", va="bottom", fontsize=10, color="#404040",
                  fontweight="bold",
                  transform=ax_A.get_xaxis_transform())
    ax_A.set_xticks(range(len(BANDS)))
    ax_A.set_xticklabels([TEX[b] for b in BANDS], fontsize=11)
    ax_A.set_xlim(-0.5, len(BANDS) - 0.5)
    ax_A.axhline(0, color="black", lw=0.6, zorder=1)
    ax_A.set_ylim(-0.18, 0.18)
    ax_A.set_ylabel(r"per-patient mean $\Delta = J_{\mathrm{post}} - J_{\mathrm{pre}}$",
                    fontsize=10)
    ax_A.set_xlabel("band", fontsize=10)
    ax_A.set_title(r"A — Cohort mean TAM persistence $\Delta$ (one dot/patient)",
                   fontsize=10.5)

    # ───── Panel B — (J_pre, J_post) scatter of all TAMs ─────
    for band in BANDS:
        sub = pm[pm["band"] == band]
        ax_B.scatter(sub["J_pre"], sub["J_post"],
                     s=8, color=BAND_COLORS[band], alpha=0.4,
                     edgecolor="none", zorder=3, label=TEX[band],
                     rasterized=True)
    # Trace box: J_pre < J_low, J_post >= J_high
    ax_B.add_patch(Rectangle((0, 0.7), 0.5, 0.3,
                             facecolor="none", edgecolor="#d62728",
                             linewidth=1.6, linestyle="--", zorder=4))
    ax_B.text(0.05, 0.98, "trace zone\n(J_pre<0.5,\n J_post≥0.7)",
              fontsize=8, color="#d62728", va="top")
    ax_B.plot([0, 1], [0, 1], color="#888", lw=0.8, ls="--", zorder=2)
    ax_B.set_xlim(0, 1); ax_B.set_ylim(0, 1)
    ax_B.set_xlabel(r"$J_{\mathrm{pre}}$ (TAM matched in rest_pre)", fontsize=10)
    ax_B.set_ylabel(r"$J_{\mathrm{post}}$ (TAM matched in rest_post)", fontsize=10)
    ax_B.legend(loc="lower right", fontsize=8, frameon=False, ncol=2,
                markerscale=1.4, handlelength=1.0, columnspacing=0.6,
                handletextpad=0.4)
    ax_B.set_title(rf"B — All {len(pm)} TAMs in $(J_\mathrm{{pre}}, J_\mathrm{{post}})$ space",
                   fontsize=10.5)
    ax_B.grid(alpha=0.25)

    # ───── Panel C — Δ histograms per band ─────
    bins = np.linspace(-0.6, 0.6, 41)
    for band in BANDS:
        sub = pm[pm["band"] == band]
        ax_C.hist(sub["delta"], bins=bins, histtype="step", lw=1.6,
                  color=BAND_COLORS[band], label=TEX[band], density=True)
    ax_C.axvline(0, color="black", lw=0.8, zorder=1)
    ax_C.set_xlabel(r"$\Delta = J_{\mathrm{post}} - J_{\mathrm{pre}}$ (per TAM)",
                    fontsize=10)
    ax_C.set_ylabel("density (TAMs)", fontsize=10)
    ax_C.set_xlim(-0.6, 0.6)
    ax_C.legend(loc="upper right", fontsize=8, frameon=False, ncol=2)
    ax_C.set_title(r"C — Per-TAM $\Delta$ distributions per band",
                   fontsize=10.5)

    fig.tight_layout()
    out = OUT_DIR / "trace_anchor.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
