#!/usr/bin/env python3
"""Cohort figures for the continuous-trace matrix.

Two figures:

(1) Cross-band dotplot (single panel, 6 bands × 10 patients):
    `cohort_all_bands.pdf` — H3-style band-by-band reading.

(2) Per-band 3-panel cohort figure (one PDF per band):
    `cohort_<band>.pdf` —
      Panel A: per-patient ρ dotplot for this band.
      Panel B: pooled-cohort hex-bin of normalised
               (Δ_task, Δ_rest) across all 10 patients.
      Panel C: per-patient stacked sign-agreement bar
               (σ>0 / σ=0 / σ<0).

Reads:  data/reports/imcoh_continuous_trace/per_cell_summary.csv
        data/reports/imcoh_continuous_trace/per_pair/{patient}_{band}.npz
Writes: data/reports/imcoh_continuous_trace/figures/cohort_*.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT


IN_DIR = REPORTS_ROOT / "imcoh_continuous_trace"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT


def render_all_bands(summary: pd.DataFrame) -> None:
    """Single panel — 6 bands, 10 patients per band, ρ on y-axis."""
    fig, ax = plt.subplots(figsize=(9.5, 5.0), dpi=150)
    rng = np.random.default_rng(0)
    for j, band in enumerate(BANDS):
        sub = summary[summary["band"] == band]
        rhos = sub["rho"].to_numpy()
        if len(rhos) == 0:
            continue
        xs = j + rng.uniform(-0.16, 0.16, size=len(rhos))
        cols = ["#2ca02c" if r > 0 else "#d62728" for r in rhos]
        ax.scatter(xs, rhos, s=38, color=cols, alpha=0.85,
                   edgecolor="black", linewidth=0.3, zorder=3)
        med = float(np.median(rhos))
        q1, q3 = float(np.quantile(rhos, 0.25)), float(np.quantile(rhos, 0.75))
        ax.vlines(j, q1, q3, color="#404040", lw=2.2, zorder=4, alpha=0.85)
        ax.hlines(med, j - 0.22, j + 0.22, color="#404040", lw=2.5, zorder=5)
        n_pos = int((rhos > 0).sum())
        ax.text(j, 1.02, f"{n_pos}/{len(rhos)}",
                ha="center", va="bottom", fontsize=10,
                color="#2ca02c", fontweight="bold",
                transform=ax.get_xaxis_transform())

    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([TEX[b] for b in BANDS], fontsize=12)
    ax.set_xlim(-0.5, len(BANDS) - 0.5)
    ax.axhline(0, color="black", lw=0.7)
    ax.set_ylim(-0.25, 1.0)
    ax.set_ylabel(r"$\rho = \mathrm{Spearman}(\Delta_{\mathrm{task}}, "
                  r"\Delta_{\mathrm{rest}})$", fontsize=11)
    ax.set_xlabel("band", fontsize=11)
    handles = [
        plt.Line2D([0], [0], marker='o', color="#2ca02c", lw=0,
                   markersize=7, markeredgecolor="black",
                   markeredgewidth=0.3, label=r"$\rho > 0$"),
        plt.Line2D([0], [0], marker='o', color="#d62728", lw=0,
                   markersize=7, markeredgecolor="black",
                   markeredgewidth=0.3, label=r"$\rho < 0$"),
        plt.Line2D([0], [0], color="#404040", lw=2.2,
                   label="cohort median + IQR"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=9)
    ax.set_title(r"Per-patient continuous-trace $\rho$ per band  "
                 r"(numbers above: patients with $\rho>0$)", fontsize=10.5)
    fig.tight_layout()
    out = OUT_DIR / "cohort_all_bands.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}")
    plt.close(fig)


def render_band_cohort(band: str, summary: pd.DataFrame) -> None:
    sub = summary[summary["band"] == band].copy()
    if sub.empty:
        return

    # Pool per-pair vectors with per-patient normalisation.
    pooled_task = []
    pooled_rest = []
    bar_records = []
    for _, r in sub.iterrows():
        pat = str(r["patient"])
        pair_path = IN_DIR / "per_pair" / f"{pat}_{band}.npz"
        if not pair_path.exists():
            continue
        pp = np.load(pair_path)
        dT = pp["dD_task"]; dR = pp["dD_rest"]; sig = pp["sigma"]
        dmax = max(float(r["dmax_pre"]), 1e-9)
        pooled_task.append(dT / dmax)
        pooled_rest.append(dR / dmax)
        bar_records.append({
            "patient": pat,
            "pos": float((sig > 0).mean()),
            "zero": float((sig == 0).mean()),
            "neg": float((sig < 0).mean()),
            "rho": float(r["rho"]),
        })
    pooled_task = np.concatenate(pooled_task) if pooled_task else np.array([])
    pooled_rest = np.concatenate(pooled_rest) if pooled_rest else np.array([])

    fig = plt.figure(figsize=(15.0, 4.6), dpi=150)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.2, 1.4])
    ax_A = fig.add_subplot(gs[0, 0])
    ax_B = fig.add_subplot(gs[0, 1])
    ax_C = fig.add_subplot(gs[0, 2])

    # --- A. per-patient ρ dotplot ---
    rng = np.random.default_rng(0)
    rhos = sub["rho"].to_numpy()
    xs = rng.uniform(-0.18, 0.18, size=len(rhos))
    cols = ["#2ca02c" if r > 0 else "#d62728" for r in rhos]
    ax_A.scatter(xs, rhos, s=44, color=cols, alpha=0.9,
                 edgecolor="black", linewidth=0.3, zorder=3)
    med = float(np.median(rhos))
    q1, q3 = float(np.quantile(rhos, 0.25)), float(np.quantile(rhos, 0.75))
    ax_A.vlines(0, q1, q3, color="#404040", lw=2.4, zorder=4)
    ax_A.hlines(med, -0.25, 0.25, color="#404040", lw=2.6, zorder=5)
    ax_A.set_xlim(-0.5, 0.5)
    ax_A.axhline(0, color="black", lw=0.7)
    ax_A.set_xticks([])
    ax_A.set_ylabel(r"$\rho$ per patient", fontsize=10)
    n_pos = int((rhos > 0).sum())
    ax_A.set_title(rf"A — $\rho$ cohort, {n_pos}/{len(rhos)} positive, "
                   rf"median $={med:+.3f}$", fontsize=10.5)

    # --- B. pooled hex-bin ---
    if len(pooled_task) > 0:
        lim = max(np.abs(pooled_task).max(), np.abs(pooled_rest).max())
        ax_B.axhline(0, color="#888", lw=0.6, zorder=1)
        ax_B.axvline(0, color="#888", lw=0.6, zorder=1)
        ax_B.plot([-lim, lim], [-lim, lim], color="#444", lw=0.8,
                  ls="--", zorder=1)
        hb = ax_B.hexbin(pooled_task, pooled_rest, gridsize=60,
                         cmap="viridis", bins="log", mincnt=1,
                         rasterized=True)
        ax_B.set_xlim(-lim, lim); ax_B.set_ylim(-lim, lim)
        from scipy.stats import spearmanr
        rho_pool = float(spearmanr(pooled_task, pooled_rest).statistic)
        ax_B.set_aspect("equal")
        ax_B.set_xlabel(r"$\Delta_{\mathrm{task}} / d_{\max}^{\mathrm{pre}}$",
                        fontsize=10)
        ax_B.set_ylabel(r"$\Delta_{\mathrm{rest}} / d_{\max}^{\mathrm{pre}}$",
                        fontsize=10)
        ax_B.set_title(rf"B — pooled {len(pooled_task)} pairs across {len(sub)} "
                       rf"patients, $\rho_{{\mathrm{{pool}}}} = {rho_pool:+.3f}$",
                       fontsize=10.5)
        fig.colorbar(hb, ax=ax_B, fraction=0.046, pad=0.04, shrink=0.85,
                     label="log density")

    # --- C. per-patient stacked σ bars ---
    bar_df = pd.DataFrame(bar_records).sort_values("rho", ascending=False)
    pats = bar_df["patient"].tolist()
    pos_b = bar_df["pos"].to_numpy()
    zero_b = bar_df["zero"].to_numpy()
    neg_b = bar_df["neg"].to_numpy()
    xs = np.arange(len(pats))
    ax_C.bar(xs, pos_b, color="#f4c430", label="σ>0 (same)",
             edgecolor="black", linewidth=0.3)
    ax_C.bar(xs, zero_b, bottom=pos_b, color="#cccccc", label="σ=0",
             edgecolor="black", linewidth=0.3)
    ax_C.bar(xs, neg_b, bottom=pos_b + zero_b, color="#7e3a8a",
             label="σ<0 (opp)", edgecolor="black", linewidth=0.3)
    ax_C.axhline(0.5, color="black", lw=0.6, ls="--")
    ax_C.set_xticks(xs)
    ax_C.set_xticklabels(pats, rotation=45, fontsize=8, ha="right")
    ax_C.set_ylim(0, 1)
    ax_C.set_ylabel("fraction of pairs", fontsize=10)
    ax_C.set_title(r"C — Per-patient sign-agreement (sorted by $\rho$)",
                   fontsize=10.5)
    ax_C.legend(loc="upper right", frameon=False, fontsize=8)

    fig.tight_layout()
    out = OUT_DIR / f"cohort_{band}.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", default=None,
                    help="Comma-sep band names (default all 6 + cross-band).")
    args = ap.parse_args()
    bands = args.bands.split(",") if args.bands else BANDS

    summary = pd.read_csv(IN_DIR / "per_cell_summary.csv")
    render_all_bands(summary)
    for band in bands:
        render_band_cohort(band, summary)


if __name__ == "__main__":
    main()
