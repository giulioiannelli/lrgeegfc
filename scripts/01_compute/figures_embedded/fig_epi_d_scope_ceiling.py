#!/usr/bin/env python3
r"""fig:epi_d — a mechanistic SOZ marker, calibrated for triage (§3, R3.4, the closer).

CORE MESSAGE (enact the ceiling, do not narrate it): the marker is a triage aid, not a
standalone localizer. Its top-five precision is ~60% -- roughly a sevenfold enrichment over
the per-patient base rate -- and precision falls back toward the base rate as the shortlist
lengthens (the ceiling is structural: target rarity x ranker quality). It applies to gray
and white matter alike: 41% of the seizure-onset contacts lie in white matter, which the
marker ranks without a gray-matter prior. Its highest-ranked UNMARKED contacts -- in
hippocampus, amygdala, and medial orbitofrontal cortex -- are candidate occult sites for
prospective testing (hypotheses, none outcome-validated); any that are genuine would raise
the precision reported here.

Panel a: precision@k vs the per-patient base-rate band (the triage number + its ceiling),
with the seizure-onset tissue composition (gray/white) beneath. Panel b: the top unmarked
candidates by P(SOZ), colored by anatomical class.

WM composition (53/129 = 41% of true SOZ are white-matter-labelled) is the detector README
figure (audit_117); tissue labels are not in the per-node CSV.

Reads : data/audit/epi_propagator_detector/detector_node_predictions.csv
        data/audit/epi_propagator_detector/occult_candidates_graymatter.csv
Writes: data/reports/results_section3/fig_epi_d_scope_ceiling.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

DET = ROOT / "data/audit/epi_propagator_detector/detector_node_predictions.csv"
OCC = ROOT / "data/audit/epi_propagator_detector/occult_candidates_graymatter.csv"
OUT = ROOT / "data/reports/results_section3/fig_epi_d_scope_ceiling.pdf"

KMAX = 15
N_OCC = 10
WM_N, SOZ_N = 53, 129                      # audit_117 README: 41% of SOZ are white-matter
C_PREC, C_BASE = "#08519c", "#9a9a9a"
C_WM, C_GM = "#b8860b", "#8a9199"
C_MTL, C_OFC, C_OTHER = "#2b8cbe", "#d1491c", "#b0b0b0"


def region_class(region: str) -> str:
    r = region.lower()
    if any(k in r for k in ("hip", "amygdal", "entorhinal", "parahippo")):
        return "MTL"
    if "orbitofrontal" in r:
        return "OFC"
    return "other"


def pretty_region(region: str) -> str:
    r = region.replace("ctx-lh-", "").replace("ctx-rh-", "")
    return {"Hip": "hippocampus", "medialorbitofrontal": "medial-OFC",
            "lateralorbitofrontal": "lateral-OFC", "superiortemporal": "sup.temporal",
            "parsopercularis": "pars operc.", "postcentral": "postcentral"}.get(r, r)


def precision_at_k(nd, kmax):
    ks = np.arange(1, kmax + 1)
    prec = np.zeros(kmax)
    for _, g in nd.groupby("patient"):
        gs = g.sort_values("p_soz", ascending=False).is_soz.to_numpy()
        for i, k in enumerate(ks):
            prec[i] += gs[:k].mean()
    return ks, prec / nd.patient.nunique()


def main():
    nd = pd.read_csv(DET)
    occ = pd.read_csv(OCC).nlargest(N_OCC, "p_soz").reset_index(drop=True)

    ks, prec = precision_at_k(nd, KMAX)
    base = nd.groupby("patient").is_soz.mean()
    b25, b50, b75 = base.quantile([0.25, 0.5, 0.75])
    p5 = prec[4]
    enrich_med = float((nd.groupby("patient").apply(
        lambda g: g.nlargest(5, "p_soz").is_soz.mean() / g.is_soz.mean())).median())

    fig = plt.figure(figsize=(12.6, 5.8))
    axa = fig.add_axes([0.075, 0.32, 0.40, 0.60])     # precision@k
    axw = fig.add_axes([0.075, 0.13, 0.40, 0.085])    # tissue composition bar
    axb = fig.add_axes([0.60, 0.13, 0.375, 0.79])     # occult candidates

    # -- panel a: precision@k + base-rate band --------------------------------
    axa.axhspan(b25, b75, color=C_BASE, alpha=0.30, zorder=0)
    axa.axhline(b50, color=C_BASE, lw=1.2, ls="--", zorder=1)
    axa.plot(ks, prec, color=C_PREC, lw=2.6, marker="o", ms=6, mec="white",
             mew=0.8, zorder=4)
    axa.scatter([5], [p5], s=200, marker="*", color=C_PREC, edgecolor="white",
                linewidths=0.8, zorder=6)
    axa.annotate(rf"precision@5 $= {p5:.2f}$" "\n" rf"$\approx{enrich_med:.0f}\times$ base rate",
                 (5, p5), textcoords="offset points", xytext=(16, 6), fontsize=11,
                 color=C_PREC, fontweight="bold")
    axa.text(KMAX, b50, "  base rate\n  (per-patient)", va="center", ha="left",
             fontsize=9, color="0.4")
    axa.set_xlim(0.5, KMAX + 2.2)
    axa.set_ylim(0, max(prec) * 1.18)
    axa.set_xlabel("shortlist length $k$ (contacts per patient)", fontsize=12)
    axa.set_ylabel("precision@$k$", fontsize=12)
    axa.spines[["top", "right"]].set_visible(False)
    axa.text(0.0, 1.03, r"$\mathbf{a}$", transform=axa.transAxes, fontsize=17,
             va="bottom", fontweight="bold")

    # -- tissue composition bar (gray vs white matter) ------------------------
    wm_frac = WM_N / SOZ_N
    axw.barh([0], [wm_frac], color=C_WM, edgecolor="white", height=0.8)
    axw.barh([0], [1 - wm_frac], left=[wm_frac], color=C_GM, edgecolor="white",
             height=0.8)
    axw.text(wm_frac / 2, 0, f"white matter\n{wm_frac*100:.0f}%", ha="center",
             va="center", fontsize=9, color="white", fontweight="bold")
    axw.text(wm_frac + (1 - wm_frac) / 2, 0, f"gray matter\n{(1-wm_frac)*100:.0f}%",
             ha="center", va="center", fontsize=9, color="white", fontweight="bold")
    axw.set_xlim(0, 1)
    axw.set_yticks([])
    axw.set_xticks([])
    for s in axw.spines.values():
        s.set_visible(False)
    axw.set_xlabel(f"seizure-onset contacts by tissue "
                   f"($n={SOZ_N}$; ranked without a gray-matter prior)", fontsize=9.5)

    # -- panel b: occult candidates (hypotheses) ------------------------------
    occ = occ.iloc[::-1].reset_index(drop=True)       # highest at top
    y = np.arange(len(occ))
    for yi, row in zip(y, occ.itertuples()):
        cls = region_class(row.region)
        col = {"MTL": C_MTL, "OFC": C_OFC, "other": C_OTHER}[cls]
        axb.barh(yi, row.p_soz, color=col, edgecolor="white", height=0.72, zorder=3)
        axb.text(row.p_soz - 0.01, yi, f"{row.p_soz:.2f}", ha="right", va="center",
                 fontsize=8.5, color="white", fontweight="bold")
    axb.set_yticks(y)
    axb.set_yticklabels([f"{r.patient.replace('Pat_','P')} · {pretty_region(r.region)}"
                         for r in occ.itertuples()], fontsize=9.5)
    for tick, row in zip(axb.get_yticklabels(), occ.itertuples()):
        cls = region_class(row.region)
        if cls != "other":
            tick.set_color({"MTL": C_MTL, "OFC": C_OFC}[cls])
            tick.set_fontweight("bold")
    axb.set_xlim(0, 1.05)
    axb.set_xlabel("P(seizure onset) at unmarked contacts", fontsize=11.5)
    axb.set_ylim(-0.6, len(occ) - 0.4)
    axb.tick_params(axis="y", length=0)
    axb.spines[["top", "right"]].set_visible(False)
    axb.text(0.0, 1.03, r"$\mathbf{b}$   candidate occult sites (hypotheses)",
             transform=axb.transAxes, fontsize=12.5, va="bottom", fontweight="bold")

    handles = [
        Line2D([0], [0], color=C_PREC, lw=2.6, marker="*", ms=13, mfc=C_PREC,
               label=r"precision@$k$ (triage)"),
        Patch(facecolor=C_BASE, alpha=0.30, label="per-patient base-rate band"),
        Patch(facecolor=C_MTL, label="mesial temporal (hippocampus/amygdala)"),
        Patch(facecolor=C_OFC, label="orbitofrontal"),
        Patch(facecolor=C_OTHER, label="other cortex"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.035),
               ncol=5, frameon=False, fontsize=9.2, handletextpad=0.5,
               columnspacing=1.3)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:epi_d — scope & ceiling\n")
    print(f"  precision@5 = {p5:.3f}  (~{enrich_med:.1f}x per-patient median base rate)")
    print(f"  base rate per-patient: q25={b25:.3f} median={b50:.3f} q75={b75:.3f}")
    print(f"  white-matter SOZ: {WM_N}/{SOZ_N} = {WM_N/SOZ_N*100:.0f}%")
    print("  top occult candidates:")
    for r in occ.iloc[::-1].itertuples():
        print(f"    {r.patient} {pretty_region(r.region):14s} P={r.p_soz:.3f} [{region_class(r.region)}]")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
