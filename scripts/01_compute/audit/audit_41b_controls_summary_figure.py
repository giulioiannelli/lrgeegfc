#!/usr/bin/env python3
"""Audit 41b -- summary heatmap + paired-line figure for measure 11 controls.

Two figures:
1. Heatmap of all 21 (probe, variant, band) cells with -log10(p) color
   and Bonferroni / BH-surviving annotations.
2. Per-cell paired-line plot showing each patient's real T_d connected
   to their null T_d, for the 5 Bonferroni-surviving cells.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

OUT_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "lrg_controls"
RPT_FIG = ROOT / "data" / "reports" / "section_5_lrg_trace" / "11_lrg_global_probe_controls" / "figures"
SUMMARY_CSV = ROOT / "data" / "audit" / "lrg_global_probe_controls" / "cohort_controls_summary.csv"

DR_NULL = ROOT / "data" / "audit" / "lrg_global_probe_controls" / "drank_null.csv"
KC_NULL = ROOT / "data" / "audit" / "lrg_global_probe_controls" / "kc_null.csv"
G_NULL = ROOT / "data" / "audit" / "lrg_global_probe_controls" / "grassmann_null.csv"

DR_REAL = ROOT / "data" / "reports" / "section_5_lrg_trace" / "02_d_rank_triangle" / "tables" / "Td_per_patient_per_band.csv"
KC_REAL = ROOT / "data" / "reports" / "section_5_lrg_trace" / "03_kc_lambda_triangle" / "tables" / "Td_per_patient_per_band_lambda.csv"
G_REAL = ROOT / "data" / "reports" / "section_5_lrg_trace" / "04_grassmann_triangle" / "tables" / "Td_per_patient_per_band_k.csv"

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {"delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
            "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$"}


def make_heatmap():
    df = pd.read_csv(SUMMARY_CSV)
    # Pivot: rows = (probe, variant), cols = bands, values = -log10(p) and pass flags
    df["row_label"] = df["probe"] + " " + df["variant"]
    # Order: D-rank first (d_S, d_P, d_F), then KC (λ=0, λ=0.5, λ=1), then Grassmann
    row_order = [
        "drank d_S", "drank d_P", "drank d_F",
        "kc lambda=0.0", "kc lambda=0.5", "kc lambda=1.0",
        f"grassmann k=13",
    ]
    # MTC inline: Bonferroni m=42 (7 probes × 6 bands) and BH-FDR q≤0.05.
    pvec = df["wilcoxon_p_real_lt_null"].to_numpy()
    M = len(pvec)
    bonf_thr = 0.05 / M
    df["passes_bonf"] = pvec < bonf_thr
    # BH-FDR: sort ascending, find max k with p_k <= k * 0.05 / M
    sort_idx = np.argsort(pvec)
    p_sorted = pvec[sort_idx]
    bh_pass = np.zeros(M, dtype=bool)
    threshold_k = -1
    for k in range(M):
        if p_sorted[k] <= (k + 1) * 0.05 / M:
            threshold_k = k
    if threshold_k >= 0:
        bh_pass[sort_idx[: threshold_k + 1]] = True
    df["passes_q05"] = bh_pass
    print(f"[heatmap] M={M} cells | Bonferroni threshold={bonf_thr:.5f} | "
          f"#bonf={int(df.passes_bonf.sum())} | #BH-q05={int(df.passes_q05.sum())}")

    pmat = np.full((len(row_order), len(BANDS)), np.nan)
    nmat = np.full((len(row_order), len(BANDS)), np.nan)
    bonf = np.zeros((len(row_order), len(BANDS)), dtype=bool)
    bh = np.zeros((len(row_order), len(BANDS)), dtype=bool)
    for i, rl in enumerate(row_order):
        for j, b in enumerate(BANDS):
            sub = df[(df.row_label == rl) & (df.band == b)]
            if not sub.empty:
                row = sub.iloc[0]
                p = row["wilcoxon_p_real_lt_null"]
                pmat[i, j] = -np.log10(p) if p > 0 else 4.0
                nmat[i, j] = row["n_real_below_null"]
                bonf[i, j] = bool(row["passes_bonf"])
                bh[i, j] = bool(row["passes_q05"])

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    im = ax.imshow(pmat, cmap="YlOrRd", vmin=0, vmax=3.5, aspect="auto")
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=12)
    pretty = {
        "drank d_S": r"D-rank $d_S$",
        "drank d_P": r"D-rank $d_P$",
        "drank d_F": r"D-rank $d_F$",
        "kc lambda=0.0": r"KC $\lambda=0$ (topology)",
        "kc lambda=0.5": r"KC $\lambda=0.5$ (balanced)",
        "kc lambda=1.0": r"KC $\lambda=1$ (heights)",
        "grassmann k=13": r"Grassmann $k=13$",
    }
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels([pretty[r] for r in row_order], fontsize=10)
    ax.set_xlabel("band", fontsize=11)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$-\log_{10}(p)$  (Wilcoxon real $<$ null)", fontsize=10)
    # Annotate each cell with n/N and significance stars
    for i in range(len(row_order)):
        for j in range(len(BANDS)):
            n = int(nmat[i, j]) if not np.isnan(nmat[i, j]) else 0
            txt = f"{n}/10"
            if bonf[i, j]:
                txt += "\n***"
            elif bh[i, j]:
                txt += "\n*"
            color = "white" if pmat[i, j] > 1.6 else "black"
            ax.text(j, i, txt, ha="center", va="center", fontsize=9, color=color)
    # Horizontal separators between probe families
    for sep in [2.5, 5.5]:
        ax.axhline(sep, color="black", lw=0.6)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "controls_heatmap.pdf")
    fig.savefig(RPT_FIG / "controls_heatmap.pdf")
    plt.close(fig)
    print(f"Wrote heatmap to {RPT_FIG / 'controls_heatmap.pdf'}")


def make_paired_lines():
    """Per-patient real vs null connected lines, for the 5 Bonferroni-surviving cells."""
    real_dr = pd.read_csv(DR_REAL)
    null_dr = pd.read_csv(DR_NULL)
    real_kc = pd.read_csv(KC_REAL)
    null_kc = pd.read_csv(KC_NULL)
    real_g = pd.read_csv(G_REAL)
    null_g = pd.read_csv(G_NULL)

    cells = [
        ("KC β λ=0 (topology)",  "kc",  "beta",      "T_KC", 0.0,         "null_T_lam0.0"),
        ("KC β λ=0.5 (balanced)", "kc",  "beta",      "T_KC", 0.5,         "null_T_lam0.5"),
        ("KC β λ=1 (heights)",   "kc",  "beta",      "T_KC", 1.0,         "null_T_lam1.0"),
        ("D-rank β d_F (magnitude)", "drank", "beta", "T_F",  None,        "null_T_F"),
        ("Grassmann low_γ k=13",  "grassmann", "low_gamma", "T_E1", None,  "null_T_E1"),
    ]
    fig, axes = plt.subplots(1, 5, figsize=(15, 4.2), sharey=False)
    PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
                "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
    for ax, (label, probe, band, real_col, lam, null_col) in zip(axes, cells):
        if probe == "drank":
            real_sub = real_dr[real_dr.band == band].set_index("patient")
            null_sub = null_dr[null_dr.band == band].set_index("patient")
        elif probe == "kc":
            real_sub = real_kc[(real_kc.band == band) & (real_kc.lam == lam)].set_index("patient")
            null_sub = null_kc[null_kc.band == band].set_index("patient")
        else:  # grassmann
            real_sub = real_g[(real_g.band == band) & (real_g.k == 13)].set_index("patient")
            null_sub = null_g[null_g.band == band].set_index("patient")
        for pat in PATIENTS:
            if pat in real_sub.index and pat in null_sub.index:
                r = float(real_sub.loc[pat, real_col])
                n = float(null_sub.loc[pat, null_col])
                color = "#1a4f73" if r < n else "#d62728"
                ax.plot([0, 1], [n, r], color=color, alpha=0.7, lw=1.2)
                ax.scatter([0, 1], [n, r], color=color, s=15, zorder=3)
                ax.annotate(pat[-2:], (1.02, r), fontsize=7, alpha=0.8, va="center")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["null\n(within-rest)", "real\n(task triangle)"], fontsize=9)
        ax.set_xlim(-0.15, 1.4)
        ax.axhline(0, color="0.5", lw=0.4, ls=":")
        ax.set_title(label, fontsize=10)
        if ax is axes[0]:
            ax.set_ylabel(r"$T_d$  (negative = trace direction)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "controls_paired_lines.pdf")
    fig.savefig(RPT_FIG / "controls_paired_lines.pdf")
    plt.close(fig)
    print(f"Wrote paired-lines to {RPT_FIG / 'controls_paired_lines.pdf'}")


if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RPT_FIG.mkdir(parents=True, exist_ok=True)
    make_heatmap()
    make_paired_lines()
