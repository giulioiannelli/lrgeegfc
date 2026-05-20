#!/usr/bin/env python3
"""Audit 54b — Figures for epi cross-phase rigidity (Direction A).

Reads the three CSVs from `audit_54_epi_rigidity_compute.py` and emits
three publication-style figures into `data/audit/epi_rigidity/figures/`:

  fig_01_cohort_class_heatmap.pdf
      3 panels (one per lambda); rows = {anchor, trace, reset, rearrange};
      cols = bands; cell value = number of patients with z_class > 1.96.
      Stars where BH q < 0.05.

  fig_02_Td_E_box_and_dot.pdf
      3 x 6 grid (lambda x band). Boxplot of per-patient Td_E with dots.
      Horizontal red dashed line at 0; trace zone (Td_E < 0) lightly shaded.
      Stars where BH q < 0.05 on Td_E vs 0 cohort Wilcoxon.

  fig_03_pairmask_vs_induced.pdf
      Side-by-side at alpha/lambda=0 (the audit_48 fig_05 cell):
      left panel = pair-mask Td (audit_48), right panel = induced Td (audit_54).
      Dots per patient labelled. Operationalisation sensitivity figure.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT


COMPUTE_DIR = ROOT / "data" / "audit" / "epi_rigidity"
FIG_DIR = COMPUTE_DIR / "figures"
AUDIT_48_DIR = ROOT / "data" / "audit" / "epileptic_n10_revisit"

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
KC_LAMBDAS = [0.0, 0.5, 1.0]
CLASS_NAMES = ["anchor", "trace", "reset", "rearrange"]


def _band_label(b: str) -> str:
    return BRAIN_BAND_TEX_DICT.get(b, b)


def fig_01_cohort_class_heatmap(per_pat: pd.DataFrame,
                                  cohort: pd.DataFrame,
                                  variant: str = "all") -> Path:
    p = per_pat[per_pat.variant == variant]
    c = cohort[cohort.variant == variant]

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6),
                              constrained_layout=True)
    n_total = int(p.patient.nunique())

    cmap = plt.get_cmap("YlGnBu")
    vmax = float(n_total)

    for col, lam in enumerate(KC_LAMBDAS):
        ax = axes[col]
        mat = np.zeros((len(CLASS_NAMES), len(BANDS)), dtype=float)
        for i, cls in enumerate(CLASS_NAMES):
            for j, band in enumerate(BANDS):
                sub = p[(p.lam == lam) & (p.band == band)]
                z = sub[f"z_class_{cls}"].to_numpy()
                mat[i, j] = float(np.sum(z > 1.96))
        im = ax.imshow(mat, vmin=0, vmax=vmax, cmap=cmap, aspect="auto")
        # Annotate counts
        for i in range(len(CLASS_NAMES)):
            for j in range(len(BANDS)):
                ax.text(j, i, f"{int(mat[i, j])}",
                         ha="center", va="center",
                         color=("white" if mat[i, j] > vmax * 0.6 else "black"),
                         fontsize=8)
        # Stars where BH q < 0.05
        for i, cls in enumerate(CLASS_NAMES):
            for j, band in enumerate(BANDS):
                row = c[(c.lam == lam) & (c.band == band) & (c.cls == cls)]
                if not row.empty:
                    q = float(row.q_bh.iloc[0])
                    if np.isfinite(q) and q < 0.05:
                        ax.text(j + 0.32, i - 0.32, "*",
                                 color="red", fontsize=12, fontweight="bold",
                                 ha="center", va="center")

        ax.set_xticks(range(len(BANDS)),
                       [_band_label(b) for b in BANDS])
        ax.set_yticks(range(len(CLASS_NAMES)), CLASS_NAMES)
        ax.set_title(rf"$\lambda = {lam:.1f}$", fontsize=10)

    # Figure-level colorbar (right side)
    cb = fig.colorbar(im, ax=axes, location="right", shrink=0.8, pad=0.02)
    cb.set_label(rf"# patients with $z_{{class}} > 1.96$ (n = {n_total})")

    out = FIG_DIR / "fig_01_cohort_class_heatmap.pdf"
    fig.savefig(out, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return out


def fig_02_Td_E_box_and_dot(per_pat: pd.DataFrame,
                              trace: pd.DataFrame,
                              variant: str = "all") -> Path:
    p = per_pat[per_pat.variant == variant]
    c = trace[trace.variant == variant]

    fig, axes = plt.subplots(len(KC_LAMBDAS), len(BANDS),
                              figsize=(13, 6.5),
                              constrained_layout=True,
                              squeeze=False)

    rng = np.random.default_rng(0)
    for r, lam in enumerate(KC_LAMBDAS):
        for col, band in enumerate(BANDS):
            ax = axes[r, col]
            sub = p[(p.lam == lam) & (p.band == band)]
            vals = sub.Td_E.dropna().to_numpy()
            ax.axhspan(-3.5, 0.0, color="#a6cee3", alpha=0.18, zorder=0)
            ax.axhline(0.0, color="red", linestyle="--", linewidth=0.8,
                        zorder=1)
            if vals.size > 0:
                ax.boxplot([vals], positions=[0], widths=0.55,
                            patch_artist=True,
                            boxprops=dict(facecolor="#dddddd", edgecolor="black"),
                            medianprops=dict(color="black"),
                            zorder=2)
                jitter = rng.uniform(-0.15, 0.15, size=vals.size)
                ax.scatter(jitter, vals, s=28, color="#1f78b4",
                            edgecolors="white", linewidths=0.6, zorder=3)
            # Stars
            row = c[(c.lam == lam) & (c.band == band)]
            if not row.empty:
                q = float(row.q_Td_vs_zero.iloc[0])
                if np.isfinite(q) and q < 0.05:
                    ax.text(0, ax.get_ylim()[1] * 0.92, "*",
                             ha="center", color="red", fontsize=14,
                             fontweight="bold")
            ax.set_xticks([])
            ax.set_xlim(-0.6, 0.6)
            ax.set_ylim(-3.5, 3.5)
            if r == 0:
                ax.set_title(_band_label(band), fontsize=10)
            if col == 0:
                ax.set_ylabel(rf"$T_d^E$ ($\lambda = {lam:.1f}$)", fontsize=9)
            else:
                ax.set_yticklabels([])

    out = FIG_DIR / "fig_02_Td_E_box_and_dot.pdf"
    fig.savefig(out, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return out


def fig_03_pairmask_vs_induced(per_pat: pd.DataFrame) -> Path:
    """Operationalisation sensitivity: pair-mask (audit_48) vs induced
    (audit_54) Td values per patient, at the alpha/lambda=0 cell.
    """
    csv48 = AUDIT_48_DIR / "M3_kc_epi_subtree.csv"
    if not csv48.exists():
        return Path("(skipped: audit_48 CSV missing)")
    df48 = pd.read_csv(csv48)
    rows48 = []
    for (pat, band, lam), g in df48.groupby(["patient", "band", "lam"]):
        d_a = g[(g.phase_a == "rest_pre") & (g.phase_b == "task_test")].d_kc_epi
        d_b = g[(g.phase_a == "task_test") & (g.phase_b == "rest_post")].d_kc_epi
        if len(d_a) == 1 and len(d_b) == 1:
            rows48.append(dict(patient=pat, band=band, lam=float(lam),
                                Td=float(d_b.iloc[0]) - float(d_a.iloc[0])))
    df48_T = pd.DataFrame(rows48)

    target_band = "alpha"
    target_lam = 0.0

    sub54 = per_pat[(per_pat.variant == "all") & (per_pat.band == target_band) &
                     (per_pat.lam == target_lam)]
    sub48 = df48_T[(df48_T.band == target_band) & (df48_T.lam == target_lam)]

    common = sorted(set(sub54.patient) & set(sub48.patient))
    if not common:
        return Path("(skipped: no overlapping patients)")

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0),
                              constrained_layout=True, sharey=True)

    for ax, df, title in [
        (axes[0], sub48.set_index("patient").loc[common],
         "audit_48 pair-mask"),
        (axes[1], sub54.set_index("patient").loc[common],
         "audit_54 induced subtree"),
    ]:
        if "Td" in df.columns:
            vals = df["Td"].to_numpy()
        else:
            vals = df["Td_E"].to_numpy()
        ax.axhspan(-2.5, 0.0, color="#a6cee3", alpha=0.18, zorder=0)
        ax.axhline(0.0, color="red", linestyle="--", linewidth=0.8)
        x = np.arange(len(common))
        colours = ["#33a02c" if v < 0 else "#e31a1c" for v in vals]
        ax.scatter(x, vals, s=70, c=colours, edgecolors="black",
                    linewidths=0.6, zorder=3)
        for xi, pat, v in zip(x, common, vals):
            ax.text(xi, v, "  " + pat.replace("Pat_", ""),
                     fontsize=8, va="center")
        n_neg = int((vals < 0).sum())
        ax.set_xticks(x, [pat.replace("Pat_", "") for pat in common],
                       rotation=0, fontsize=8)
        ax.set_xlabel("patient")
        ax.set_title(rf"{title} — $\alpha$, $\lambda=0$, n_neg={n_neg}/{len(vals)}",
                      fontsize=10)
    axes[0].set_ylabel(r"$T_d^{\,\mathrm{epi}} = d_{KC}(tt, post) - d_{KC}(pre, tt)$")

    out = FIG_DIR / "fig_03_pairmask_vs_induced.pdf"
    fig.savefig(out, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    per_pat = pd.read_csv(COMPUTE_DIR / "M_class_per_patient.csv")
    cohort = pd.read_csv(COMPUTE_DIR / "M_class_cohort.csv")
    trace = pd.read_csv(COMPUTE_DIR / "M_trace_score_cohort.csv")

    out1 = fig_01_cohort_class_heatmap(per_pat, cohort, variant="all")
    out2 = fig_02_Td_E_box_and_dot(per_pat, trace, variant="all")
    out3 = fig_03_pairmask_vs_induced(per_pat)

    print(f"[audit_54b] wrote: {out1}")
    print(f"[audit_54b] wrote: {out2}")
    print(f"[audit_54b] wrote: {out3}")


if __name__ == "__main__":
    main()
