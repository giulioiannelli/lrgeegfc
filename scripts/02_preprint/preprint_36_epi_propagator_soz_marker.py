#!/usr/bin/env python3
"""preprint_36 — epileptic diffusion-community finding + SOZ marker figures.

Three PDFs for the epilepsy section (orthogonal to the cross-phase trace story):

1. fig_epi_propagator_controls.pdf — Part A. Per band, the epi diffusion community
   z against its controls C2 (EE−NN contrast), C3 (strength-matched random subset)
   and C5 (strength+spatial-matched random subset), all at the C3-headline τ. Bars
   coloured by cohort Wilcoxon significance. δ/β/low-γ (and α) survive C3+C5;
   high-γ is the clean internal null.

2. fig_epi_soz_marker.pdf — the marker. (a) LOPO cross-patient classifier AUC,
   relational vs strength-only (strength is at chance cross-patient). (b) within-
   patient masking recovery AUC, strength-residual relational vs strength.
   (c) calibration of P(SOZ) pooled over the airtight bands.

3. fig_epi_occult_discovery.pdf — discovery. For example patients, the top occult
   (unlabelled) candidates ranked by cross-band P(SOZ), point colour = distance to
   nearest labelled SOZ, with the known-SOZ P(SOZ) band for reference. Hypotheses.

Inputs (read-only, no recompute):
  data/audit/epi_propagator_blocks/{propagator_blocks_cohort,spatial_null_cohort}.csv
  data/audit/epi_marker_relational/{classifier_lopo,masked_recovery_cohort,
        classifier_calibration,occult_candidates}.csv
PDF only, full vector, transparent.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BLK = ROOT / "data/audit/epi_propagator_blocks"
REL = ROOT / "data/audit/epi_marker_relational"
OUT = ROOT / "data/preprint/figures/all_bands"
OUT.mkdir(parents=True, exist_ok=True)

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
AIRTIGHT = ["delta", "beta", "low_gamma"]
C_SIG = "#1b7837"      # significant (green)
C_NS = "#bdbdbd"       # non-significant (grey)
C_REL = "#2166ac"      # relational marker (blue)
C_STR = "#b2182b"      # strength baseline (red)


def _tex(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


# --------------------------------------------------------------------------
# Figure 1 — Part A controls
# --------------------------------------------------------------------------
def fig_controls():
    coh = pd.read_csv(BLK / "propagator_blocks_cohort.csv")
    spat = pd.read_csv(BLK / "spatial_null_cohort.csv")
    controls = [("z_EE_minus_NN", "C2  EE−NN"),
                ("z_rand_EE", "C3  strength-matched"),
                ("z_spat_EE", "C5  strength+spatial")]
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.1), sharey=True)
    for ax, (key, title) in zip(axes, controls):
        zs, sig = [], []
        for b in BANDS:
            cb = coh[coh.band == b]
            # headline τ = argmax of the decisive C3 z
            if cb.empty or cb["med_z_rand_EE"].isna().all():
                zs.append(np.nan); sig.append(False); continue
            ti = int(cb.loc[cb.med_z_rand_EE.idxmax(), "tau_idx"])
            if key == "z_spat_EE":
                row = spat[(spat.band == b) & (spat.tau_idx == ti)]
                z = float(row.med_z_spat_EE.iloc[0]) if not row.empty else np.nan
                p = float(row.wilcoxon_p_spat_EE.iloc[0]) if not row.empty else np.nan
            else:
                row = cb[cb.tau_idx == ti]
                z = float(row[f"med_{key}"].iloc[0])
                p = float(row[f"wilcoxon_p_{key.replace('z_', '')}"].iloc[0])
            zs.append(z); sig.append(np.isfinite(p) and p < 0.05)
        x = np.arange(len(BANDS))
        ax.bar(x, zs, color=[C_SIG if s else C_NS for s in sig],
               edgecolor="black", linewidth=0.5, width=0.74)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([_tex(b) for b in BANDS], rotation=45, ha="right")
        ax.text(0.03, 0.95, title, transform=ax.transAxes, va="top", fontsize=9)
    axes[0].set_ylabel(r"cohort median $z$ (epi block)")
    h = [plt.Rectangle((0, 0), 1, 1, fc=C_SIG, ec="black"),
         plt.Rectangle((0, 0), 1, 1, fc=C_NS, ec="black")]
    fig.legend(h, ["cohort Wilcoxon p < 0.05", "n.s."], loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    out = OUT / "fig_epi_propagator_controls.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  -> {out.name}")


# --------------------------------------------------------------------------
# Figure 2 — the marker (classifier + recovery + calibration)
# --------------------------------------------------------------------------
def fig_marker():
    clf = pd.read_csv(REL / "classifier_lopo.csv")
    rec = pd.read_csv(REL / "masked_recovery_cohort.csv")
    cal = pd.read_csv(REL / "classifier_calibration.csv")
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.2))
    x = np.arange(len(BANDS))
    w = 0.38

    # (a) LOPO classifier AUC: relational vs strength
    ax = axes[0]
    rel = [float(clf[clf.band == b].auc_relational.iloc[0]) if not clf[clf.band == b].empty else np.nan for b in BANDS]
    str_ = [float(clf[clf.band == b].auc_strength.iloc[0]) if not clf[clf.band == b].empty else np.nan for b in BANDS]
    ax.bar(x - w / 2, rel, w, color=C_REL, edgecolor="black", linewidth=0.4, label="relational")
    ax.bar(x + w / 2, str_, w, color=C_STR, edgecolor="black", linewidth=0.4, label="strength")
    ax.axhline(0.5, color="black", lw=0.8, ls="--")
    ax.set_ylim(0.3, 0.9)
    ax.set_ylabel("LOPO ROC-AUC")
    ax.text(0.03, 0.96, "a  cross-patient classifier", transform=ax.transAxes, va="top", fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels([_tex(b) for b in BANDS], rotation=45, ha="right")

    # (b) masking recovery: resid vs strength
    ax = axes[1]
    res = [float(rec[rec.band == b].med_auc_resid_all.iloc[0]) if not rec[rec.band == b].empty else np.nan for b in BANDS]
    rst = [float(rec[rec.band == b].med_auc_strength_all.iloc[0]) if not rec[rec.band == b].empty else np.nan for b in BANDS]
    ax.bar(x - w / 2, res, w, color=C_REL, edgecolor="black", linewidth=0.4, label="relational (strength-resid.)")
    ax.bar(x + w / 2, rst, w, color=C_STR, edgecolor="black", linewidth=0.4, label="strength")
    ax.axhline(0.5, color="black", lw=0.8, ls="--")
    ax.set_ylim(0.3, 0.9)
    ax.set_ylabel("masking recovery AUC")
    ax.text(0.03, 0.96, "b  within-patient recovery", transform=ax.transAxes, va="top", fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels([_tex(b) for b in BANDS], rotation=45, ha="right")

    # (c) calibration pooled over airtight bands
    ax = axes[2]
    ca = cal[cal.band.isin(AIRTIGHT)].copy()
    # bin by predicted prob across pooled airtight bins
    if not ca.empty:
        order = ca.sort_values("mean_pred")
        ax.plot(order.mean_pred, order.obs_freq, "o-", color=C_REL, ms=4, lw=1.2)
    ax.plot([0, 0.5], [0, 0.5], color="black", lw=0.8, ls="--")
    ax.set_xlabel("predicted P(SOZ)")
    ax.set_ylabel("observed SOZ frequency")
    ax.text(0.03, 0.96, "c  calibration (δ/β/low-γ)", transform=ax.transAxes, va="top", fontsize=9)

    handles_a = [plt.Rectangle((0, 0), 1, 1, fc=C_REL, ec="black"),
                 plt.Rectangle((0, 0), 1, 1, fc=C_STR, ec="black")]
    fig.legend(handles_a, ["relational diffusion marker", "node strength (hubness)"],
               loc="lower center", bbox_to_anchor=(0.5, -0.03), ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    out = OUT / "fig_epi_soz_marker.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  -> {out.name}")


# --------------------------------------------------------------------------
# Figure 3 — occult discovery
# --------------------------------------------------------------------------
def _short_region(r):
    r = str(r).replace("ctx-lh-", "L-").replace("ctx-rh-", "R-")
    return r[:18]


def fig_discovery():
    cand = pd.read_csv(REL / "occult_candidates.csv")
    pats = ["Pat_06", "Pat_08", "Pat_13"]
    pats = [p for p in pats if p in set(cand.patient)]
    fig, axes = plt.subplots(1, len(pats), figsize=(4.0 * len(pats), 3.4))
    if len(pats) == 1:
        axes = [axes]
    dmax = np.nanpercentile(cand.dist_nearest_soz_um.to_numpy() / 1000.0, 95)
    sc = None
    for ax, pat in zip(axes, pats):
        g = cand[cand.patient == pat]
        epi_p = g[g.is_epi].P_mean_airtight.dropna()
        top = g[g.is_top_candidate].sort_values("P_mean_airtight")
        if epi_p.size:
            ax.axvspan(np.nanpercentile(epi_p, 25), np.nanpercentile(epi_p, 75),
                       color="#fddbc7", alpha=0.6, zorder=0)
            ax.axvline(epi_p.median(), color=C_STR, lw=1.0, ls=":", zorder=1)
        y = np.arange(len(top))
        ax.hlines(y, 0, top.P_mean_airtight, color="#888888", lw=0.8, zorder=2)
        if {"P_mean_lo", "P_mean_hi"}.issubset(top.columns):
            ax.errorbar(top.P_mean_airtight, y,
                        xerr=[top.P_mean_airtight - top.P_mean_lo,
                              top.P_mean_hi - top.P_mean_airtight],
                        fmt="none", ecolor="#555555", elinewidth=0.8,
                        capsize=2, zorder=2.5)
        sc = ax.scatter(top.P_mean_airtight, y, c=top.dist_nearest_soz_um / 1000.0,
                        cmap="turbo", vmin=0, vmax=dmax, s=46, edgecolor="black",
                        linewidth=0.5, zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{_short_region(r)}" for r in top.region], fontsize=7)
        ax.set_xlim(0, max(0.4, float(top.P_mean_airtight.max()) * 1.15) if len(top) else 1)
        ax.set_xlabel("cross-band P(SOZ)")
        ax.text(0.03, 1.02, pat.replace("_", " "), transform=ax.transAxes, va="bottom", fontsize=9)
    if sc is not None:
        cb = fig.colorbar(sc, ax=axes, fraction=0.025, pad=0.02)
        cb.set_label("dist. to nearest SOZ (mm)")
    fig.subplots_adjust(wspace=0.62)
    out = OUT / "fig_epi_occult_discovery.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  -> {out.name}")


def main():
    print("[preprint_36] epi propagator SOZ-marker figures:")
    fig_controls()
    fig_marker()
    fig_discovery()
    print(f"[preprint_36] done -> {OUT}/")


if __name__ == "__main__":
    main()
