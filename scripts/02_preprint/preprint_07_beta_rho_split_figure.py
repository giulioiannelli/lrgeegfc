#!/usr/bin/env python3
"""Per-pair persistence (`ρ_split^coph`) headline preprint figure (per band).

Three-panel composite, redesigned 2026-05-22 to put the per-patient
z-score above the matched-strength noise floor as the primary visual
element. The cohort decile tilt and the drift/split/xprobe paired
slope corroborate.

Panel (a)  per-patient z above own matched-strength surrogate
           (R = 200, 4-cycle ±δ rewiring). Patients sorted by z
           descending; grey vertical band marks the universal
           noise-floor zone |z| ≤ 2; stems coloured green for
           z > +2 (trace), red for z < −2 (anti), grey otherwise.
           At a glance β shows 7/10 patients past z = +2 with the
           remaining three sitting inside / just outside the band —
           the band-discriminator quantity is the cohort being
           uniformly *above* the band, never below it.

Panel (b)  cohort-pooled binned regression of within-patient ranks,
           with the ten per-patient curves shown underneath as
           thin grey lines. Cohort line = thick blue + bootstrap
           CI; shaded band = within-patient shuffle null. No
           quadrant tinting, no corner labels, no callout box.

Panel (c)  per-patient ρ_drift → ρ_split → ρ_xprobe paired slope
           graph (replaces the 3-boxplot null triangle). Drift→split
           rise visualises C2 (trace above rest-only drift); split↔
           xprobe parallel visualises C4 (cross-probe non-degradation).

Usage
-----
    python preprint_07_beta_rho_split_figure.py             # default: beta
    python preprint_07_beta_rho_split_figure.py --band alpha
    python preprint_07_beta_rho_split_figure.py --band low_gamma

Inputs
------
data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
data/audit/ctm_triangle/Td_per_patient_per_band.csv
data/audit/ctm_triangle/cohort_summary.csv
data/audit/ctm_triangle/c4_wilcoxon_cohort.csv
data/reports/imcoh_continuous_trace/per_pair_split/<Pat>_<band>.npz

Output (PDF only, no PNG sibling)
---------------------------------
data/preprint/figures/<band>/per_pair_trace/fig_<band>_rho_split.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()


COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]

LRG_CTM_DIR = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"
CTM_TRI_DIR = ROOT / "data" / "audit" / "ctm_triangle"
PAIR_SPLIT_DIR = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"

CLR_OBS = "#1f3d6e"
CLR_PRO = "#1a7c3e"
CLR_ANTI = "#c0392b"

# Slope-graph palette (panel c)
CLR_RHO_SPLIT = "#1f3d6e"
CLR_RHO_DRIFT = "#a48b22"
CLR_RHO_XPROBE = "#5d3a8f"


def short(pat: str) -> str:
    return pat.replace("Pat_", "P")


def load_pair_data(pat: str, band: str) -> dict:
    npz_path = PAIR_SPLIT_DIR / f"{pat}_{band}.npz"
    d = np.load(npz_path)
    return dict(
        dD_task=np.asarray(d["dD_task"]),
        dD_rest=np.asarray(d["dD_rest"]),
    )


def panel_a_zstrip(ax: plt.Axes, per_pat: pd.DataFrame, cohort_row: pd.Series,
                   band: str, band_tex: str) -> None:
    """Per-patient z above own matched-strength surrogate, sorted by z.

    Each patient is one horizontal stem from x = 0 to the observed z.
    A grey vertical band at |z| ≤ 2 marks the universal noise-floor
    zone (because z normalises by the patient's own surrogate std,
    the band is shared across patients). Stem + dot are coloured by
    sign of excess: green for z > +2 (trace), red for z < −2 (anti),
    grey otherwise.
    """
    sub = per_pat[per_pat.band == band].set_index("patient").loc[COHORT]
    z = sub["obs_z"].values
    order = np.argsort(z)[::-1]
    patients_sorted = [COHORT[i] for i in order]
    z_sorted = z[order]

    n = len(z_sorted)
    y = np.arange(n)[::-1]

    ax.axvspan(-2, 2, color="0.88", alpha=0.55, zorder=0)
    ax.axvline(0, color="0.35", lw=0.8, zorder=1)
    ax.axvline(2, color="0.55", lw=0.7, ls=":", zorder=1)
    ax.axvline(-2, color="0.55", lw=0.7, ls=":", zorder=1)

    for yi, zi in zip(y, z_sorted):
        if zi > 2:
            clr = CLR_PRO
        elif zi < -2:
            clr = CLR_ANTI
        else:
            clr = "0.45"
        ax.plot([0, zi], [yi, yi], color=clr, lw=1.7, alpha=0.78,
                zorder=2, solid_capstyle="round")
        ax.scatter([zi], [yi], s=95, color=clr, edgecolor="white",
                   linewidth=1.1, zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels([short(p) for p in patients_sorted], fontsize=9)
    ax.set_ylim(-1.0, n - 0.35)

    z_max = max(float(z_sorted.max()) + 1.0, 3.5)
    z_min = min(float(z_sorted.min()) - 1.0, -3.0)
    ax.set_xlim(z_min, z_max)
    ax.set_xlabel(r"$z$ above own matched-strength surrogate",
                  fontsize=10)

    ax.text(0, -0.75, r"noise floor ($|z| \leq 2$)",
            ha="center", va="top",
            fontsize=8, color="0.30", style="italic", zorder=4)

    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(rf"(a) {band_tex} per-patient signal above noise floor",
                 fontsize=10.5, loc="left", pad=4)

    n_above_p95 = str(cohort_row["n_above_surrogate"])
    n_above_2 = int((z_sorted > 2).sum())
    p_val = float(cohort_row["paired_wilcoxon_p"])
    annot = (
        rf"{n_above_2}/10 at $z > 2$"
        rf"$\quad\mid\quad$"
        rf"{n_above_p95} above own surrogate p95"
        rf"$\quad\mid\quad$"
        rf"Wilcoxon $p = {p_val:.3f}$"
    )
    ax.text(0.5, -0.20, annot, transform=ax.transAxes,
            ha="center", va="top", fontsize=8.6, color="0.20")


def panel_b_tilt(ax: plt.Axes, cohort_row: pd.Series, band: str,
                 band_tex: str) -> None:
    """Cohort + per-patient rank tilt — polished, no clutter."""
    rng = np.random.default_rng(20260520)
    n_bins = 10
    bin_edges = np.linspace(0.0, 1.0, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    rank_data = []
    for pat in COHORT:
        pair = load_pair_data(pat, band)
        n = len(pair["dD_task"])
        if n < 2:
            continue
        t_rank = (rankdata(pair["dD_task"]) - 1) / (n - 1)
        r_rank = (rankdata(pair["dD_rest"]) - 1) / (n - 1)
        bin_id = np.minimum(np.digitize(t_rank, bin_edges[1:-1]), n_bins - 1)
        rank_data.append(dict(r_rank=r_rank, bin_id=bin_id, n=n))

    P = len(rank_data)
    per_pat_bin = np.full((P, n_bins), np.nan)
    for i, d in enumerate(rank_data):
        for k in range(n_bins):
            mask = d["bin_id"] == k
            if mask.any():
                per_pat_bin[i, k] = d["r_rank"][mask].mean()
    cohort_bin = np.nanmean(per_pat_bin, axis=0)

    B = 2000
    boot = np.empty((B, n_bins))
    for b in range(B):
        idx = rng.integers(0, P, P)
        boot[b] = np.nanmean(per_pat_bin[idx], axis=0)
    ci_lo = np.percentile(boot, 2.5, axis=0)
    ci_hi = np.percentile(boot, 97.5, axis=0)

    n_perm = 500
    null_curve = np.empty((n_perm, n_bins))
    for p in range(n_perm):
        per_pat_bin_null = np.full((P, n_bins), np.nan)
        for i, d in enumerate(rank_data):
            r_shuf = rng.permutation(d["r_rank"])
            for k in range(n_bins):
                mask = d["bin_id"] == k
                if mask.any():
                    per_pat_bin_null[i, k] = r_shuf[mask].mean()
        null_curve[p] = np.nanmean(per_pat_bin_null, axis=0)
    null_p5 = np.percentile(null_curve, 2.5, axis=0)
    null_p95 = np.percentile(null_curve, 97.5, axis=0)

    for i in range(P):
        ax.plot(bin_centers, per_pat_bin[i], color="0.60", lw=0.9,
                alpha=0.55, zorder=2)

    ax.fill_between(bin_centers, null_p5, null_p95, color="0.55",
                    alpha=0.28, linewidth=0, zorder=1.5)
    ax.axhline(0.5, color="0.40", lw=0.7, ls="--", alpha=0.55, zorder=1)

    ax.fill_between(bin_centers, ci_lo, ci_hi, color=CLR_OBS, alpha=0.30,
                    linewidth=0, zorder=3)
    ax.plot(bin_centers, cohort_bin, color=CLR_OBS, lw=2.4, marker="o",
            markersize=7.5, markerfacecolor=CLR_OBS,
            markeredgecolor="white", markeredgewidth=1.0, zorder=4)

    obs_half_span = max(
        float(np.nanmax(np.abs(cohort_bin - 0.5))),
        float(np.nanmax(np.abs(ci_lo - 0.5))),
        float(np.nanmax(np.abs(ci_hi - 0.5))),
        float(np.nanmax(np.abs(per_pat_bin - 0.5))),
    )
    y_half = max(0.12, obs_half_span + 0.020)
    ax.set_ylim(0.5 - y_half, 0.5 + y_half)
    ax.set_xlim(0.0, 1.0)

    ax.set_xticks([0.05, 0.5, 0.95])
    ax.set_xticklabels(["0.05", "0.50", "0.95"])
    ax.set_xlabel(r"within-patient rank of $\Delta_{\mathrm{task}}(i,j)$",
                  fontsize=10)
    ax.set_ylabel(r"cohort mean rank of $\Delta_{\mathrm{rest}}(i,j)$",
                  fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(rf"(b) cohort pair-rank tilt at {band_tex}",
                 fontsize=10.5, loc="left", pad=4)

    null_mid = float((null_p5[0] + null_p95[0]) / 2)
    ax.text(0.04, null_mid, "null", ha="left", va="center",
            fontsize=8, color="0.25", style="italic", zorder=2.5)
    ax.text(0.96, cohort_bin[-1] + 0.005,
            "cohort", ha="right", va="bottom",
            fontsize=8.6, color=CLR_OBS, fontweight="bold", zorder=5)

    rho_cohort = float(cohort_row["obs_median_rho"])
    delta_lr = float(cohort_bin[-1] - cohort_bin[0])
    n_pairs_total = int(sum(d["n"] for d in rank_data))
    annot = (
        rf"cohort $\rho = {rho_cohort:+.3f}$"
        rf"$\quad\mid\quad$"
        rf"$\Delta_{{10-1}} = {delta_lr:+.3f}$"
        rf"$\quad\mid\quad$"
        rf"$N = {n_pairs_total}$ pairs"
    )
    ax.text(0.5, -0.20, annot, transform=ax.transAxes,
            ha="center", va="top", fontsize=8.6, color="0.20")


def panel_c_slope(ax: plt.Axes, td_per_pat: pd.DataFrame,
                  td_cohort_row: pd.Series, c4_row: pd.Series,
                  band: str, band_tex: str) -> None:
    """Per-patient drift / split / xprobe paired slope graph.

    Three x-positions: ρ_drift, ρ_split^coph, ρ_xprobe. Each patient
    is a thin grey 3-point polyline; cohort medians are large
    coloured dots joined by a thick black line. The drift→split rise
    visualises C2; the split↔xprobe parallel visualises C4.
    """
    sub = td_per_pat[td_per_pat.band == band].set_index("patient").loc[COHORT]
    rho_drift = sub["rho_null_drift"].values
    rho_split = sub["rho_split"].values
    rho_xp = sub["rho_split_cross_probe"].values

    xs = np.array([0.0, 1.0, 2.0])

    for d_val, s_val, xp_val in zip(rho_drift, rho_split, rho_xp):
        ax.plot(xs, [d_val, s_val, xp_val], color="0.55", lw=0.9,
                alpha=0.50, solid_capstyle="round", zorder=2)
        ax.scatter(xs, [d_val, s_val, xp_val], s=22, color="0.45",
                   alpha=0.70, zorder=3, edgecolor="white", linewidth=0.4)

    med = [float(np.median(rho_drift)),
           float(np.median(rho_split)),
           float(np.median(rho_xp))]
    cohort_colors = [CLR_RHO_DRIFT, CLR_RHO_SPLIT, CLR_RHO_XPROBE]
    ax.plot(xs, med, color="0.12", lw=2.0, zorder=4)
    for x, m, c in zip(xs, med, cohort_colors):
        ax.scatter([x], [m], s=180, color=c, edgecolor="0.12",
                   linewidth=1.5, zorder=5)

    ax.axhline(0, color="0.55", lw=0.7, ls="--", zorder=1)

    ax.set_xticks(xs)
    ax.set_xticklabels([
        r"$\rho_{\mathrm{drift}}$",
        r"$\rho_{\mathrm{split}}^{\mathrm{coph}}$",
        r"$\rho_{\mathrm{xprobe}}$",
    ], fontsize=10)
    ax.set_xlim(-0.35, 2.35)
    ax.set_ylabel(r"Spearman $\rho$", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(rf"(c) {band_tex} drift / split / xprobe",
                 fontsize=10.5, loc="left", pad=4)

    p1 = float(td_cohort_row["wilcoxon_split_gt_0_p"])
    p2 = float(td_cohort_row["wilcoxon_split_gt_drift_p"])
    c4_paired_p = float(c4_row["paired_wilcoxon_p_split_gt_xprobe"])
    annot = (
        rf"C1 $p = {p1:.3f}$"
        rf"$\quad\mid\quad$"
        rf"C2 $p = {p2:.3f}$"
        rf"$\quad\mid\quad$"
        rf"C4 $p = {c4_paired_p:.3f}$"
    )
    ax.text(0.5, -0.20, annot, transform=ax.transAxes,
            ha="center", va="top", fontsize=8.6, color="0.20")


def main(band: str = "beta") -> Path:
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")

    per_pat = pd.read_csv(LRG_CTM_DIR / "per_patient_per_band.csv")
    cohort = pd.read_csv(LRG_CTM_DIR / "cohort_summary.csv")
    cohort_row = cohort[cohort.band == band].iloc[0]
    td_per_pat = pd.read_csv(CTM_TRI_DIR / "Td_per_patient_per_band.csv")
    td_cohort = pd.read_csv(CTM_TRI_DIR / "cohort_summary.csv")
    td_cohort_row = td_cohort[td_cohort.band == band].iloc[0]
    c4_cohort = pd.read_csv(CTM_TRI_DIR / "c4_wilcoxon_cohort.csv")
    c4_row = c4_cohort[c4_cohort.band == band].iloc[0]

    out_dir = ROOT / "data" / "preprint" / "figures" / band / "per_pair_trace"
    out_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(16.0, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.75], wspace=0.30)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    panel_a_zstrip(ax_a, per_pat, cohort_row, band, band_tex)
    panel_b_tilt(ax_b, cohort_row, band, band_tex)
    panel_c_slope(ax_c, td_per_pat, td_cohort_row, c4_row, band, band_tex)

    fig.tight_layout(rect=(0, 0.04, 1, 0.97))

    out = out_dir / f"fig_{band}_rho_split.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--band", default="beta",
                        choices=["delta", "theta", "alpha", "beta",
                                 "low_gamma", "high_gamma"])
    args = parser.parse_args()
    main(args.band)
