#!/usr/bin/env python3
"""§5.7 Figure 2 — LRG-CTM headline under matched-strength null.

Re-bases the headline on the audit_round3 `ctm_headline.pdf` boxplot
layout (per-band three boxes + per-patient connectors + asterisks on
the controlled box) but swaps the control battery from the
within-baseline triple (split / drift / x-probe) to the
matched-strength comparison:

    ρ_split        — observed Spearman cophenetic ρ on D(τ) (audit_63)
    ρ_surr         — matched-strength surrogate per-patient p50
    ρ_split (epi-X) — α only, audit_68

Per-patient connectors join each patient's observed ρ_split to its
own-surrogate median (and to its epi-X observation where the box
applies). Cohort-paired one-sided Wilcoxon "obs > surr" p-value is
printed above each band's observed box with the same asterisk
convention as the headline.

Inputs:
    data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
    data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
    data/audit/alpha_epi_exclusion/per_patient.csv
    data/audit/alpha_epi_exclusion/cohort_summary.csv

Output:
    data/reports/section_5_matched_strength_refinement/figures/
        fig_2_ctm_headline_matched_strength.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from lrg_eegfc.visuals.styles import use_lrg_style

from _shared_ms import (
    ALL_BANDS, ALPHA_EPIX_DIR, BAND_LABELS, CLR_EPIX, FIG_DIR, LRG_CTM_DIR,
)


CLR_OBS_BOX = "#1f77b4"
CLR_SURR_BOX = "#888888"
CLR_EPIX_BOX = CLR_EPIX
CLR_PRO = "#1a7c3e"
CLR_ANTI = "#c0392b"


def asterisks(p: float) -> str:
    if not np.isfinite(p):
        return ""
    if p <= 0.001:
        return "***"
    if p <= 0.01:
        return "**"
    if p <= 0.05:
        return "*"
    return ""


def main() -> Path:
    use_lrg_style()

    per_pat = pd.read_csv(LRG_CTM_DIR / "per_patient_per_band.csv")
    cohort = pd.read_csv(LRG_CTM_DIR / "cohort_summary.csv")
    epiX_pp = pd.read_csv(ALPHA_EPIX_DIR / "per_patient.csv")
    epiX_cohort = pd.read_csv(ALPHA_EPIX_DIR / "cohort_summary.csv")

    bands = [b for b in ALL_BANDS if b in set(cohort.band)]
    x_band = np.arange(len(bands))
    width = 0.26

    fig, ax = plt.subplots(figsize=(13.5, 5.4))

    # ---- boxes ------------------------------------------------------------
    box_specs = [
        ("obs",   -width, CLR_OBS_BOX),
        ("surr",   0.0,   CLR_SURR_BOX),
        ("epiX", +width, CLR_EPIX_BOX),
    ]

    for spec, off, color in box_specs:
        data = []
        for b in bands:
            sub = per_pat[per_pat.band == b]
            if spec == "obs":
                vals = sub.obs_rho.values
            elif spec == "surr":
                vals = sub.surr_p50.values
            elif spec == "epiX":
                # α only — empty list for other bands so the slot is still
                # there but nothing draws.
                if b == "alpha":
                    vals = epiX_pp.obs_rho.values
                else:
                    vals = np.array([])
            data.append(vals)
        positions = x_band + off
        # Drop empties for boxplot but keep alignment via separate calls.
        for xi, vals in zip(positions, data):
            if vals.size == 0:
                continue
            ax.boxplot(
                vals, positions=[xi], widths=width * 0.85,
                patch_artist=True, showmeans=False, showfliers=False,
                boxprops=dict(facecolor=color, alpha=0.50, edgecolor="#444"),
                medianprops=dict(color="#111", lw=1.4),
                whiskerprops=dict(color="#666"),
                capprops=dict(color="#666"),
            )

    # ---- per-patient connectors --------------------------------------------
    pats = sorted(per_pat.patient.unique())
    for xi, b in enumerate(bands):
        sub = per_pat[per_pat.band == b].set_index("patient")
        for p in pats:
            if p not in sub.index:
                continue
            obs_val = float(sub.loc[p, "obs_rho"])
            surr_val = float(sub.loc[p, "surr_p50"])
            xs = [xi - width, xi]
            ys = [obs_val, surr_val]
            if b == "alpha":
                epi_row = epiX_pp[epiX_pp.patient == p]
                if not epi_row.empty:
                    xs.append(xi + width)
                    ys.append(float(epi_row.obs_rho.iloc[0]))
            line_c = CLR_PRO if obs_val > 0 else CLR_ANTI
            ax.plot(xs, ys, color=line_c, alpha=0.45, lw=0.8, zorder=3)
            ax.scatter(xs, ys, s=10, color=line_c, alpha=0.65,
                       edgecolor="white", linewidth=0.3, zorder=4)

    # ---- significance text above the observed box --------------------------
    for xi, b in enumerate(bands):
        c_row = cohort[cohort.band == b].iloc[0]
        p = float(c_row.paired_wilcoxon_p)
        ratio = (abs(c_row.obs_median_rho)
                 / max(abs(c_row.surr_median_rho_median), 1e-12))
        vals = per_pat[per_pat.band == b].obs_rho.dropna().values
        if not vals.size:
            continue
        q3 = float(np.nanpercentile(vals, 75))
        top = max(q3, float(np.nanmax(vals)))
        anchor = top + 0.04
        stars = asterisks(p)
        col = "#c0392b" if p <= 0.05 else "#444"
        txt = rf"$p={p:.3f}${stars}" "\n" rf"ratio $={ratio:.1f}\times$"
        weight = "bold" if stars else "normal"
        ax.text(xi - width, anchor, txt, ha="center", va="bottom",
                fontsize=8, color=col, fontweight=weight)

        if b == "alpha":
            ex = epiX_cohort.iloc[0]
            ex_p = float(ex.paired_wilcoxon_p)
            ex_ratio = (abs(ex.obs_median_rho)
                        / max(abs(ex.surr_median_rho_median), 1e-12))
            ex_vals = epiX_pp.obs_rho.dropna().values
            ex_top = max(float(np.nanpercentile(ex_vals, 75)),
                          float(np.nanmax(ex_vals))) + 0.04
            ex_stars = asterisks(ex_p)
            ex_col = "#c0392b" if ex_p <= 0.05 else "#444"
            ax.text(xi + width, ex_top,
                    rf"$p={ex_p:.3f}${ex_stars}" "\n"
                    rf"ratio $={ex_ratio:.1f}\times$",
                    ha="center", va="bottom",
                    fontsize=8, color=ex_col,
                    fontweight="bold" if ex_stars else "normal")

    ax.axhline(0, color="0.3", lw=0.7, ls="--")
    ax.set_xticks(x_band)
    ax.set_xticklabels([BAND_LABELS[b] for b in bands], fontsize=12)
    ax.set_ylabel(r"$\rho_{\mathrm{split}}^{\mathrm{LRG}}$  on $D(\tau)$",
                  fontsize=11)

    # y-limits from p1/p99 over all rendered values to give headroom for
    # the asterisks.
    all_v = []
    for b in bands:
        all_v.append(per_pat[per_pat.band == b].obs_rho.values)
        all_v.append(per_pat[per_pat.band == b].surr_p50.values)
    all_v.append(epiX_pp.obs_rho.values)
    all_v = np.concatenate(all_v)
    qmin = float(np.nanpercentile(all_v, 1))
    qmax = float(np.nanpercentile(all_v, 99))
    ax.set_ylim(qmin - 0.18, qmax + 0.28)
    ax.spines[["top", "right"]].set_visible(False)

    handles = [
        mpatches.Patch(facecolor=CLR_OBS_BOX, alpha=0.50, edgecolor="#444",
                        label=r"$\rho_{\mathrm{split}}$  (observed)"),
        mpatches.Patch(facecolor=CLR_SURR_BOX, alpha=0.50, edgecolor="#444",
                        label=r"$\rho_{\mathrm{surr}}$  (matched-strength p50)"),
        mpatches.Patch(facecolor=CLR_EPIX_BOX, alpha=0.50, edgecolor="#444",
                        label=r"$\rho_{\mathrm{split}}$  (epi-excluded, $\alpha$)"),
        mlines.Line2D([], [], color=CLR_PRO, marker="o", lw=0.8, alpha=0.7,
                       markersize=4, label=r"patient: pro ($\rho_{\mathrm{split}}>0$)"),
        mlines.Line2D([], [], color=CLR_ANTI, marker="o", lw=0.8, alpha=0.7,
                       markersize=4, label=r"patient: anti ($\rho_{\mathrm{split}}\leq 0$)"),
    ]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.03),
               ncol=5, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 1))

    out = FIG_DIR / "fig_2_ctm_headline_matched_strength.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
