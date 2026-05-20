#!/usr/bin/env python3
"""Per-pair persistence (`ρ_split^coph`) headline preprint figure (per band).

Two-panel composite for any band's per-pair paragraph of the preprint:

Panel (a)  per-patient observed ρ_split^coph vs own matched-strength
           surrogate p50, paired connectors. Cohort median + matched-
           strength gap ratio annotated.

Panel (b)  per-pair Δ_task(i,j) vs Δ_rest(i,j) scatter for the cohort-
           median exemplar patient (from selection.csv if present,
           otherwise auto-picked as the patient whose obs_rho is closest
           to the cohort median).

Usage
-----
    python preprint_07_beta_rho_split_figure.py             # default: beta
    python preprint_07_beta_rho_split_figure.py --band alpha
    python preprint_07_beta_rho_split_figure.py --band low_gamma

Inputs
------
data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
data/audit/ctm_triangle/cohort_summary.csv
data/audit/ctm_per_pair_scatter/tables/selection.csv  (β/α/θ only)
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
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.probe import extract_probe_labels
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
SEL_CSV = ROOT / "data" / "audit" / "ctm_per_pair_scatter" / "tables" / "selection.csv"
SEEG_RAW = ROOT / "data" / "raw" / "stereoeeg_patients"

CLR_OBS = "#1f3d6e"
CLR_SURR = "#7d7d7d"
CLR_PRO = "#1a7c3e"
CLR_ANTI = "#c0392b"
CLR_SP = "#bababa"
CLR_XP = "#1f3d6e"

# Null-triangle palette (panel c)
CLR_RHO_SPLIT = "#1f3d6e"
CLR_RHO_DRIFT = "#a48b22"
CLR_RHO_XPROBE = "#5d3a8f"


def short(pat: str) -> str:
    return pat.replace("Pat_", "P")


def load_pair_data(pat: str, band: str) -> dict:
    npz_path = PAIR_SPLIT_DIR / f"{pat}_{band}.npz"
    d = np.load(npz_path)
    iu_i = d["iu_i"].astype(int)
    iu_j = d["iu_j"].astype(int)
    n_contacts = int(max(iu_i.max(), iu_j.max())) + 1
    ch = pd.read_csv(SEEG_RAW / pat / "channel_labels.csv")
    probes = np.array(extract_probe_labels(ch["label"].tolist()))[:n_contacts]
    same_probe = probes[iu_i] == probes[iu_j]
    return dict(
        dD_task=np.asarray(d["dD_task"]),
        dD_rest=np.asarray(d["dD_rest"]),
        same_probe=same_probe,
    )


def panel_a_strip(ax: plt.Axes, per_pat: pd.DataFrame, cohort: pd.Series,
                  band: str, band_tex: str) -> None:
    """Per-patient observed vs own matched-strength surrogate strip plot.

    For each patient: horizontal box-and-whisker over the R=200
    strength-preserving surrogates (p5–p25–p50–p75–p95) overlaid with the
    filled observed `ρ_split^coph` marker.
    """
    sub = per_pat[per_pat.band == band].set_index("patient").loc[COHORT]
    n = len(sub)
    y = np.arange(n)[::-1]

    obs = sub["obs_rho"].values
    p5 = sub["surr_p5"].values
    p25 = sub["surr_p25"].values
    p50 = sub["surr_p50"].values
    p75 = sub["surr_p75"].values
    p95 = sub["surr_p95"].values

    box_h = 0.46
    surr_fill = "#dcdcdc"
    for yi, _p5, _p25, _p50, _p75, _p95 in zip(y, p5, p25, p50, p75, p95):
        # p5–p95 whisker
        ax.plot([_p5, _p95], [yi, yi], color=CLR_SURR, lw=1.1,
                solid_capstyle="round", zorder=2)
        # p25–p75 box
        ax.add_patch(mpatches.Rectangle(
            (_p25, yi - box_h / 2.0), _p75 - _p25, box_h,
            facecolor=surr_fill, edgecolor=CLR_SURR, linewidth=1.0,
            zorder=3))
        # p50 tick
        ax.plot([_p50, _p50], [yi - box_h / 2.0, yi + box_h / 2.0],
                color=CLR_SURR, lw=1.6, zorder=4)

    pro_mask = obs > 0
    ax.scatter(obs[pro_mask], y[pro_mask], s=72, marker="o",
               color=CLR_PRO, edgecolor="white", linewidth=0.8, zorder=6)
    ax.scatter(obs[~pro_mask], y[~pro_mask], s=72, marker="o",
               color=CLR_ANTI, edgecolor="white", linewidth=0.8, zorder=6)

    obs_med = float(cohort["obs_median_rho"])
    surr_med = float(cohort["surr_median_rho_median"])
    ratio = abs(obs_med) / max(abs(surr_med), 1e-12)
    p_val = float(cohort["paired_wilcoxon_p"])
    n_above = str(cohort["n_above_surrogate"])

    ax.axvline(0, color="0.55", lw=0.7, ls="--", zorder=1)
    ax.axvline(obs_med, color=CLR_OBS, lw=1.1, ls=":", zorder=1)
    ax.axvline(surr_med, color=CLR_SURR, lw=1.0, ls=":", zorder=1)

    ax.set_yticks(y)
    ax.set_yticklabels([short(p) for p in COHORT], fontsize=8.5)
    ax.set_xlabel(r"$\rho_{\mathrm{split}}^{\mathrm{coph}}$",
                  fontsize=11)
    ax.set_title(rf"(a) per-patient {band_tex} persistence on $D_{{\mathrm{{coph}}}}$",
                 fontsize=10.5, loc="left", pad=4)

    xmin = min(obs.min(), p5.min(), -0.05) - 0.04
    xmax = max(obs.max(), p95.max(),  0.05) + 0.10
    ax.set_xlim(xmin, xmax)
    ax.spines[["top", "right"]].set_visible(False)

    obs_trace = int((obs > 0).sum())
    txt = (
        rf"cohort median (obs) $= {obs_med:+.3f}$"
        "\n"
        rf"cohort median (surr) $= {surr_med:+.3f}$"
        "\n"
        rf"ratio $= {ratio:.1f}\times$"
        "\n"
        rf"$n_{{\rho>0}}={obs_trace}/10$, "
        rf"$n_{{>\mathrm{{surr}}}} = {n_above}$"
        "\n"
        rf"Wilcoxon $p = {p_val:.3f}$"
    )
    ax.text(0.98, 0.02, txt, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8.5,
            bbox=dict(facecolor="white", edgecolor="#bbb",
                      boxstyle="round,pad=0.32"))


def pick_exemplar_patient(band: str, per_pat: pd.DataFrame, cohort_row: pd.Series,
                          sel: pd.DataFrame) -> str:
    """Return the cohort-median exemplar patient for a given band.

    Prefer selection.csv if it has a row for this band; otherwise pick
    the patient whose obs_rho is closest to the cohort median.
    """
    sel_band = sel[sel.band == band]
    if not sel_band.empty:
        return sel_band.iloc[0].patient
    sub = per_pat[per_pat.band == band].set_index("patient").loc[COHORT]
    cohort_med = float(cohort_row["obs_median_rho"])
    idx = (sub["obs_rho"] - cohort_med).abs().idxmin()
    return idx


def panel_b_scatter(ax: plt.Axes, sel: pd.DataFrame, per_pat: pd.DataFrame,
                    cohort_row: pd.Series, band: str, band_tex: str) -> None:
    """Per-pair Δ_task vs Δ_rest for the cohort-median exemplar patient."""
    pat = pick_exemplar_patient(band, per_pat, cohort_row, sel)
    pair = load_pair_data(pat, band)
    sp = pair["same_probe"]
    cp = ~sp
    rho_full, _ = spearmanr(pair["dD_task"], pair["dD_rest"])
    rho_x, _ = spearmanr(pair["dD_task"][cp], pair["dD_rest"][cp])

    lim = float(np.quantile(np.abs(
        np.concatenate([pair["dD_task"], pair["dD_rest"]])), 0.985)) * 1.06

    ax.scatter(pair["dD_task"][sp], pair["dD_rest"][sp],
               s=6, color=CLR_SP, alpha=0.50, edgecolors="none", zorder=1)
    ax.scatter(pair["dD_task"][cp], pair["dD_rest"][cp],
               s=7, color=CLR_XP, alpha=0.55, edgecolors="none", zorder=2)
    ax.axhline(0, color="#cccccc", lw=0.6, zorder=0)
    ax.axvline(0, color="#cccccc", lw=0.6, zorder=0)
    ax.plot([-lim, lim], [-lim, lim], color="#888888", lw=0.6, ls="--",
            zorder=0)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\Delta_{\mathrm{task}}(i,j) = "
                  r"D_{\mathrm{coph}}^{\mathrm{taskT}} - "
                  r"D_{\mathrm{coph}}^{\mathrm{rsPre,A}}$",
                  fontsize=10)
    ax.set_ylabel(r"$\Delta_{\mathrm{rest}}(i,j) = "
                  r"D_{\mathrm{coph}}^{\mathrm{rsPost}} - "
                  r"D_{\mathrm{coph}}^{\mathrm{rsPre,B}}$",
                  fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

    ax.set_title(rf"(b) per-pair geometry ({pat}, {band_tex}; "
                 rf"cohort-median exemplar)",
                 fontsize=10.5, loc="left", pad=4)

    n_total = pair["dD_task"].size
    n_same = int(sp.sum())
    n_cross = int(cp.sum())
    txt = (
        rf"$\rho_{{\mathrm{{split}}}} = {rho_full:+.3f}$"
        "\n"
        rf"$\rho_{{\mathrm{{xprobe}}}} = {rho_x:+.3f}$"
        "\n"
        rf"$N = {n_total}$ pairs "
        rf"({n_same} same / {n_cross} cross)"
    )
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, ha="left", va="top",
            fontsize=9,
            bbox=dict(facecolor="white", edgecolor="#bbb",
                      boxstyle="round,pad=0.32"))


def panel_c_null_triangle(ax: plt.Axes, td_per_pat: pd.DataFrame,
                          td_cohort_row: pd.Series, band: str,
                          band_tex: str) -> None:
    """Per-patient null triangle (`ρ_split^coph`, `ρ_drift`, `ρ_xprobe`)
    rendered as three boxplots, addressing C1 (split > 0), C2 (split
    > drift, paired) and C4 (cross-probe sign + count).
    """
    sub = td_per_pat[td_per_pat.band == band].set_index("patient").loc[COHORT]
    rho_split = sub["rho_split"].values
    rho_drift = sub["rho_null_drift"].values
    rho_xp = sub["rho_split_cross_probe"].values

    data = [rho_split, rho_drift, rho_xp]
    colors = [CLR_RHO_SPLIT, CLR_RHO_DRIFT, CLR_RHO_XPROBE]
    cats = [r"$\rho_{\mathrm{split}}^{\mathrm{coph}}$",
            r"$\rho_{\mathrm{drift}}$",
            r"$\rho_{\mathrm{xprobe}}$"]
    xs = np.array([0.0, 1.0, 2.0])

    bp = ax.boxplot(data, positions=xs, widths=0.55, patch_artist=True,
                    showfliers=False, zorder=2,
                    medianprops=dict(color="white", lw=2.0),
                    whiskerprops=dict(color="0.30", lw=1.2),
                    capprops=dict(color="0.30", lw=1.2),
                    boxprops=dict(lw=1.2))
    for patch, clr in zip(bp["boxes"], colors):
        patch.set_facecolor(clr)
        patch.set_edgecolor("0.20")
        patch.set_alpha(0.85)

    # Per-patient dots overlay (jittered)
    rng = np.random.default_rng(20260519)
    for x, vals, clr in zip(xs, data, colors):
        jx = rng.uniform(-0.13, 0.13, size=len(vals))
        ax.scatter(np.full_like(vals, x) + jx, vals,
                   s=22, color=clr, edgecolor="white", linewidth=0.5,
                   alpha=0.95, zorder=4)

    # Per-patient connectors (split → drift, split → xprobe)
    for s, d, xp in zip(rho_split, rho_drift, rho_xp):
        ax.plot([xs[0], xs[1]], [s, d], color="#cccccc", lw=0.5,
                alpha=0.6, zorder=1)
        ax.plot([xs[0], xs[2]], [s, xp], color="#cccccc", lw=0.5,
                alpha=0.6, zorder=1)

    ax.axhline(0, color="0.55", lw=0.7, ls="--", zorder=0)

    ax.set_xticks(xs)
    ax.set_xticklabels(cats, fontsize=10)
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylabel(r"Spearman $\rho$", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(rf"(c) {band_tex} null triangle (C1 / C2 / C4)",
                 fontsize=10.5, loc="left", pad=4)

    p1 = float(td_cohort_row["wilcoxon_split_gt_0_p"])
    p2 = float(td_cohort_row["wilcoxon_split_gt_drift_p"])
    n_xp = int(td_cohort_row["n_trace_xprobe_int"])
    rho_xp_med = float(td_cohort_row["rho_xprobe_median"])
    rho_sp_med = float(td_cohort_row["rho_split_median"])
    sign_match = (np.sign(rho_xp_med) == np.sign(rho_sp_med)) and (n_xp >= 6)
    txt = (
        rf"C1 $\rho_{{\mathrm{{split}}}}>0$: $p = {p1:.4f}$"
        "\n"
        rf"C2 $\rho_{{\mathrm{{split}}}}>\rho_{{\mathrm{{drift}}}}$: $p = {p2:.4f}$"
        "\n"
        rf"C4 $\rho_{{\mathrm{{xprobe}}}}={rho_xp_med:+.3f}$, "
        rf"$n={n_xp}/10$"
        "\n"
        rf"sign + $n$ check: {'pass' if sign_match else 'fail'}"
    )
    ax.text(0.97, 0.03, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.5,
            bbox=dict(facecolor="white", edgecolor="#bbb",
                      boxstyle="round,pad=0.32"))


def main(band: str = "beta") -> Path:
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")

    per_pat = pd.read_csv(LRG_CTM_DIR / "per_patient_per_band.csv")
    cohort = pd.read_csv(LRG_CTM_DIR / "cohort_summary.csv")
    cohort_row = cohort[cohort.band == band].iloc[0]
    sel = pd.read_csv(SEL_CSV)
    td_per_pat = pd.read_csv(CTM_TRI_DIR / "Td_per_patient_per_band.csv")
    td_cohort = pd.read_csv(CTM_TRI_DIR / "cohort_summary.csv")
    td_cohort_row = td_cohort[td_cohort.band == band].iloc[0]

    out_dir = ROOT / "data" / "preprint" / "figures" / band / "per_pair_trace"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Layout: 1×3. (a) per-patient strip + (b) per-pair scatter + (c) null
    # triangle boxplots. (c) is narrower since it carries only three categories.
    fig = plt.figure(figsize=(18.0, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 0.62], wspace=0.30)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    panel_a_strip(ax_a, per_pat, cohort_row, band, band_tex)
    panel_b_scatter(ax_b, sel, per_pat, cohort_row, band, band_tex)
    panel_c_null_triangle(ax_c, td_per_pat, td_cohort_row, band, band_tex)

    handles = [
        mpatches.Patch(facecolor="#dcdcdc", edgecolor=CLR_SURR,
                       label="matched-strength surrogate IQR (p25–p75)"),
        mlines.Line2D([], [], color=CLR_SURR, lw=1.6,
                      label="surrogate p50 / whisker (p5–p95)"),
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=8,
                      markerfacecolor=CLR_PRO, markeredgecolor="white",
                      label=r"observed: pro ($\rho_{\mathrm{split}}^{\mathrm{coph}} > 0$)"),
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=8,
                      markerfacecolor=CLR_ANTI, markeredgecolor="white",
                      label=r"observed: anti ($\rho_{\mathrm{split}}^{\mathrm{coph}} \leq 0$)"),
        mpatches.Patch(facecolor=CLR_SP, label="same-probe pairs"),
        mpatches.Patch(facecolor=CLR_XP, alpha=0.6, label="cross-probe pairs"),
        mpatches.Patch(facecolor=CLR_RHO_SPLIT, alpha=0.85,
                       label=r"$\rho_{\mathrm{split}}^{\mathrm{coph}}$ per-patient"),
        mpatches.Patch(facecolor=CLR_RHO_DRIFT, alpha=0.85,
                       label=r"$\rho_{\mathrm{drift}}$ (rest-only null)"),
        mpatches.Patch(facecolor=CLR_RHO_XPROBE, alpha=0.85,
                       label=r"$\rho_{\mathrm{xprobe}}$ (cross-probe only)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=9,
               frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.06, 1, 1))

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
