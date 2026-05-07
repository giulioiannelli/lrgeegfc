#!/usr/bin/env python3
"""Audit 50 -- §5.5 anatomy figure (3-band version w/ heterogeneity disclosure).

Updated `anatomy_2d.pdf` with three explicit requirements pinned to §5.5:

  1. α + β + γ_l enrichment bars + axial MNI scatter (audit_45's version
     dropped α).
  2. Anti-aligned patients Pat_07 + Pat_15 shown with a distinct marker
     (X) so the reader sees that they contribute few cortical trace-
     leaves, mostly Wm/Unk.
  3. Pat_02 trace-leaves at γ_l ctx-lh-fusiform are flagged with a
     thicker yellow halo so the reader can see Pat_02 carries 9/13 of
     the γ_l fusiform peak.

Overwrites both:
    data/outputs/figures/section_5_lrg_trace/headline/anatomy_2d.pdf
    data/reports/section_5_lrg_trace/headline/figures/anatomy_2d.pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import hypergeom

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
ANTI_ALIGNED = {"Pat_07", "Pat_15"}
ANAT_CSV = ROOT / "data" / "audit" / "lrg_localization_anatomy" / "per_trace_leaf.csv"
OUT_FIG = (ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace"
           / "headline" / "anatomy_2d.pdf")
OUT_RPT = (ROOT / "data" / "reports" / "section_5_lrg_trace" / "headline"
           / "figures" / "anatomy_2d.pdf")

BANDS = ("alpha", "beta", "low_gamma")
BAND_TEX = {"alpha": r"$\alpha$", "beta": r"$\beta$",
            "low_gamma": r"$\gamma_l$"}
BAND_COLOR = {"alpha": "#a3b35c", "beta": "#d62728", "low_gamma": "#e89c40"}


def load_cohort_contacts() -> pd.DataFrame:
    rows = []
    for pat in PATIENTS:
        impl = pd.read_csv(ROOT / "data" / "raw" / "stereoeeg_patients" / pat
                           / f"implant_pat_{pat[-2:]}.csv")
        dk_col = next(c for c in impl.columns if c.startswith("Desikan"))
        for _, r in impl.iterrows():
            region = str(r[dk_col]).split(",")[0].strip()
            rows.append({"patient": pat, "region": region,
                         "x": float(r["x"]), "y": float(r["y"]),
                         "z": float(r["z"])})
    return pd.DataFrame(rows)


def enrichment_for_band(df_leaf: pd.DataFrame, base: pd.DataFrame,
                        band: str) -> pd.DataFrame:
    """Per-region enrichment with hypergeometric p (same recipe as audit_45)."""
    sub = df_leaf[(df_leaf.band == band) & (~df_leaf.region.isin(["Wm", "Unk"]))]
    base_cort = base[~base.region.isin(["Wm", "Unk"])]
    K = len(sub)
    N_total = len(base_cort)
    base_rate = K / N_total
    base_per_reg = base_cort.groupby("region").size().rename("n_contacts").reset_index()
    per_reg = (sub.groupby("region")
               .agg(n_trace=("leaf_id", "count"),
                    n_pat_trace=("patient", "nunique"))
               .reset_index())
    m = per_reg.merge(base_per_reg, on="region", how="left")
    m["trace_rate"] = m.n_trace / m.n_contacts
    m["enrichment"] = m.trace_rate / base_rate
    m["p_hyper"] = m.apply(
        lambda r: float(hypergeom.sf(r.n_trace - 1, N_total, K, int(r.n_contacts))),
        axis=1)
    m = m[(m.n_contacts >= 5) & (m.n_pat_trace >= 2) & (m.enrichment > 1.0)]
    return m.sort_values("enrichment", ascending=True).tail(8), base_rate


def make_figure() -> None:
    df = pd.read_csv(ANAT_CSV)
    base = load_cohort_contacts()
    base_cort = base[~base.region.isin(["Wm", "Unk"])]

    enrich = {}
    base_rate = {}
    for b in BANDS:
        e, br = enrichment_for_band(df, base, b)
        enrich[b] = e
        base_rate[b] = br

    fig = plt.figure(figsize=(17.5, 10.0))
    # 3 columns for panels + 1 thin column for shared colorbar.
    gs = fig.add_gridspec(2, 4, width_ratios=[1.0, 1.0, 1.0, 0.04],
                          height_ratios=[1.0, 1.25],
                          hspace=0.35, wspace=0.32)

    # ------------------------------------------------------------------
    # Top row: enrichment bar charts (α / β / γ_l)
    # ------------------------------------------------------------------
    for col_idx, band in enumerate(BANDS):
        ax = fig.add_subplot(gs[0, col_idx])
        agg = enrich[band]
        if agg.empty:
            ax.text(0.5, 0.5, "no enriched regions", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(BAND_TEX[band] + " trace-leaves",
                         fontsize=11, color=BAND_COLOR[band], fontweight="bold")
            ax.spines[["top", "right"]].set_visible(False)
            continue
        colors = []
        alphas = []
        for _, r in agg.iterrows():
            reg = r.region
            if reg.startswith("ctx-lh"):
                col = "#2c5b96"
            elif reg.startswith("ctx-rh"):
                col = "#c83737"
            elif reg in ("Hip", "Amy") or "Hippocampus" in reg or "Amygdala" in reg:
                col = "#7a4ba8"
            else:
                col = "#7f7f7f"
            colors.append(col)
            p = r.p_hyper
            alphas.append(1.0 if p < 0.05 else 0.6 if p < 0.10 else 0.3)
        bars = ax.barh(range(len(agg)), agg.enrichment, color=colors,
                       edgecolor="0.2", linewidth=0.6, height=0.74)
        for b_obj, a in zip(bars, alphas):
            b_obj.set_alpha(a)
        ax.axvline(1.0, color="0.3", lw=1.0, ls="--", zorder=0)
        for i, (_, r) in enumerate(agg.iterrows()):
            stars = ("***" if r.p_hyper < 0.001 else
                     "**" if r.p_hyper < 0.01 else
                     "*" if r.p_hyper < 0.05 else "n.s.")
            ax.text(r.enrichment + 0.08, i,
                    f"{r.trace_rate*100:.0f}% rate "
                    f"({int(r.n_trace)}/{int(r.n_contacts)} contacts, "
                    f"{int(r.n_pat_trace)} pts) [p={r.p_hyper:.3f} {stars}]",
                    va="center", fontsize=7.5, color="0.15")
        ax.set_yticks(range(len(agg)))
        ax.set_yticklabels(agg.region.tolist(), fontsize=8.5)
        ax.set_xlabel("enrichment", fontsize=9)
        title = (rf"{BAND_TEX[band]} trace-leaf enrichment "
                 rf"(baseline {base_rate[band]*100:.1f}%)")
        ax.set_title(title, fontsize=10.5, color=BAND_COLOR[band],
                     fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlim(0, max(agg.enrichment.max() * 1.65, 4.0))

    # ------------------------------------------------------------------
    # Bottom row: axial scatter — three bands
    #     pro-trace patients = circle, anti-aligned (Pat_07/Pat_15) = X
    #     dot color = region enrichment
    #     Pat_02 γ_l fusiform leaves = additional yellow halo
    # ------------------------------------------------------------------
    last_sc = None
    for col_idx, band in enumerate(BANDS):
        ax = fig.add_subplot(gs[1, col_idx])
        sub = df[df.band == band].copy()
        et = enrich[band].set_index("region") if not enrich[band].empty else None
        sub["region_enrichment"] = (sub.region.map(et["enrichment"]) if et is not None
                                    else 0.0)
        sub["region_enrichment"] = sub.region_enrichment.fillna(0.0)
        sub["is_anti"] = sub.patient.isin(ANTI_ALIGNED)

        # Background: all cortical sEEG contacts (very faint)
        ax.scatter(base_cort.x / 1000, base_cort.y / 1000, color="0.92",
                   s=3, alpha=0.45, edgecolors="none", zorder=1)
        ax.axvline(0, color="0.5", lw=0.6, ls="--", zorder=1)

        # Pat_02 γ_l fusiform halo (only on γ_l panel)
        if band == "low_gamma":
            halo = sub[(sub.patient == "Pat_02") & (sub.region == "ctx-lh-fusiform")]
            if not halo.empty:
                ax.scatter(halo.x / 1000, halo.y / 1000, s=240, marker="o",
                           facecolor="none", edgecolor="#f0c14b", linewidths=2.6,
                           zorder=4)

        # Pro-trace dots — non-enriched grey
        non_enr = sub[(sub.region_enrichment <= 1.0) & (~sub.is_anti)]
        ax.scatter(non_enr.x / 1000, non_enr.y / 1000, color="0.55", s=30,
                   alpha=0.55, edgecolor="0.25", linewidths=0.4, zorder=3)
        # Pro-trace dots — enriched colored
        enr = sub[(sub.region_enrichment > 1.0) & (~sub.is_anti)].copy()
        sc = None
        if not enr.empty:
            sc = ax.scatter(enr.x / 1000, enr.y / 1000, c=enr.region_enrichment,
                            cmap="Reds", vmin=1.0, vmax=4.0, s=95,
                            alpha=0.95, edgecolor="black", linewidths=0.7,
                            zorder=5)
            last_sc = sc
        # Anti-aligned X markers (always drawn last)
        anti = sub[sub.is_anti]
        if not anti.empty:
            ax.scatter(anti.x / 1000, anti.y / 1000, color="#3a3a3a", s=70,
                       marker="x", linewidths=1.6, alpha=0.85, zorder=6)

        ax.text(0.04, 0.96, "L", transform=ax.transAxes, fontsize=13,
                color="0.3", fontweight="bold", ha="left", va="top")
        ax.text(0.96, 0.96, "R", transform=ax.transAxes, fontsize=13,
                color="0.3", fontweight="bold", ha="right", va="top")
        ax.text(0.5, 0.97, "anterior", transform=ax.transAxes,
                fontsize=8, color="0.4", ha="center", va="top")
        ax.text(0.5, 0.04, "posterior", transform=ax.transAxes,
                fontsize=8, color="0.4", ha="center", va="bottom")

        ax.set_xlabel(r"MNI $i$ (mm)", fontsize=9.5)
        if col_idx == 0:
            ax.set_ylabel(r"MNI $j$ (mm)", fontsize=9.5)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"axial view — {BAND_TEX[band]}",
                     fontsize=10, color=BAND_COLOR[band], fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(True, alpha=0.18, lw=0.4)
        ax.set_xlim(-90, 90)
        ax.set_ylim(-100, 80)

    # Shared colorbar to the right of the bottom row
    if last_sc is not None:
        cbar_ax = fig.add_subplot(gs[1, 3])
        cb = fig.colorbar(last_sc, cax=cbar_ax)
        cb.set_label("region enrichment", fontsize=9, labelpad=4)
        cb.ax.tick_params(labelsize=8)

    # Shared figure-level legend
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color="#2c5b96", label="LH cortex"),
        plt.Rectangle((0, 0), 1, 1, color="#c83737", label="RH cortex"),
        plt.Rectangle((0, 0), 1, 1, color="#7a4ba8", label="subcortical (Hip/Amy)"),
        plt.Line2D([0], [0], color="0.3", lw=1.0, ls="--",
                   label="cohort baseline (enrichment = 1)"),
        plt.Line2D([0], [0], marker="o", color="black", lw=0,
                   markerfacecolor="0.55", markersize=7,
                   label="pro-trace, non-enriched region"),
        plt.Line2D([0], [0], marker="x", color="#3a3a3a", lw=0, markersize=8,
                   markeredgewidth=1.6,
                   label=r"anti-aligned patient (Pat_07, Pat_15)"),
        plt.Line2D([0], [0], marker="o", color="#f0c14b", lw=0,
                   markerfacecolor="none", markersize=11, markeredgewidth=2.3,
                   label=r"Pat_02 leaves at $\gamma_l$ fusiform"),
    ]
    fig.legend(handles=legend_handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.015), ncol=4, frameon=False, fontsize=8.5)

    fig.subplots_adjust(left=0.115, right=0.985, top=0.95, bottom=0.085)
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    OUT_RPT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, bbox_inches="tight")
    fig.savefig(OUT_RPT, bbox_inches="tight")
    plt.close(fig)
    print(f"[50] wrote {OUT_FIG}")
    print(f"[50] wrote {OUT_RPT}")


if __name__ == "__main__":
    make_figure()
