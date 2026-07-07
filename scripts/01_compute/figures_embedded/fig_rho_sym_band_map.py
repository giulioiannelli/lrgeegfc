#!/usr/bin/env python3
r"""Band map of the cophenetic trace: cohort verdict vs per-patient reality.

Brutally honest version. The cohort gate does NOT cleanly separate the bands:
at n=10 the signed-rank p is discrete, and low-gamma / high-gamma / delta sit at
p=0.080 -- one patient short of the p<0.05 line (they are 6-7/10 positive; beta
is itself only 7/10 and clears on MAGNITUDE, not count). So this gate alone
cannot license "trace only in alpha/beta"; it can only say alpha is most
consistent (8/10), theta is absent (3/10, anti), and beta/gamma/delta share a
marginal 6-7/10 band. The beta-flagship claim lives on magnitude + multi-probe
convergence (gate + Grassmann + OFC), shown elsewhere -- not here.

  A  COHORT VERDICT -- magnitude (median rho_sym) x consistency (rank-biserial
     of the gate). THREE tiers, not a binary: clears (p<.05), MARGINAL
     (.05<=p<.15, the amber zone -> low/high-gamma, delta), absent (theta).
     Marker area = trace share of cross-phase structure (taxonomy comp_trace):
     beta is dominant (26%) and clear; alpha thin (3%) but clear; gamma/delta
     marginal; theta anti.

  B  PER-PATIENT REALITY -- patient x band rho_sym, three cell states:
     solid box = clears own matched-strength null (TRACE), dotted box = clears
     the LOWER tail (significant ANTI / reset), diagonal hatch = below the noise
     threshold (not distinguishable from null). Shows the off-scale tracers are
     2-3 broad tracers (P06 all six, P05 the gamma pair), and the bottom four
     patients (P10/P14/P13/P15) are non-tracers or anti.

Reads
-----
- data/audit/rho_sym_gate/per_patient_per_band.csv       (obs_rho, surr_p50, obs_p_one_sided)
- data/audit/rho_sym_gate/cohort_summary.csv             (verdict_sym)
- data/audit/cross_phase_taxonomy_rhosym/cohort_summary.csv  (comp_trace_med = %struct)

Writes
------
- data/reports/rho_sym_band_map/fig_rho_sym_band_map.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from scipy.stats import rankdata, wilcoxon

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style
from lrgsglib.plotlib import imshow_colorbar_caxdivider

ROOT = setup_script_env()
use_lrg_style()

GATE = ROOT / "data/audit/rho_sym_gate/per_patient_per_band.csv"
COH = ROOT / "data/audit/rho_sym_gate/cohort_summary.csv"
TAX = ROOT / "data/audit/cross_phase_taxonomy_rhosym/cohort_summary.csv"
OUT = ROOT / "data/reports/rho_sym_band_map/fig_rho_sym_band_map.pdf"

TIER_C = {"clear": "#1f77b4", "marginal": "#e0902f", "fail": "#9e9e9e"}
TRACE_P, ANTI_P = 0.05, 0.95     # per-patient own-null thresholds (upper/lower tail)
RB_CLEAR = 0.636                 # rank-biserial <-> Wilcoxon p=0.05 at n=10 (W+=45)
RB_MARG = 0.45                   # <-> p~0.13 (W+=40): lower edge of the marginal zone


def rank_biserial(d):
    nz = d[d != 0]
    r = rankdata(np.abs(nz))
    return (r[nz > 0].sum() - r[nz < 0].sum()) / (r.sum())


def _tier(p):
    return "clear" if p < 0.05 else ("marginal" if p < 0.15 else "fail")


def _band_table(g, coh, tax):
    rows = []
    for b, sub in g.groupby("band"):
        o = sub.obs_rho.to_numpy()
        d = o - sub.surr_p50.to_numpy()
        p = float(wilcoxon(d, alternative="greater").pvalue)
        rows.append(dict(
            band=b, med=float(np.median(o)), rb=rank_biserial(d),
            share=float(tax.loc[b, "comp_trace_med"]),
            n_pos=int((d > 0).sum()), gate_p=p, tier=_tier(p),
        ))
    return pd.DataFrame(rows).sort_values("med", ascending=False).reset_index(drop=True)


# ------------------------------------------------------------------ panel A
def panel_A(ax, T):
    ax.axvspan(RB_CLEAR, 1.0, color=TIER_C["clear"], alpha=0.06, zorder=0)
    ax.axvspan(RB_MARG, RB_CLEAR, color=TIER_C["marginal"], alpha=0.08, zorder=0)
    ax.axvline(RB_CLEAR, color=TIER_C["clear"], ls="--", lw=1.0, alpha=0.6, zorder=1)
    ax.axhline(0, color="0.6", ls="--", lw=0.8, zorder=1)
    for _, r in T.iterrows():
        c = TIER_C[r.tier]
        ax.scatter([r.rb], [r.med], s=40 + r.share * 2500, marker="o", zorder=4,
                   facecolor=c if r.tier != "fail" else "white", edgecolor=c,
                   lw=1.8, alpha=0.95)
    lab = {  # (dx, dy, ha, va)
        "beta": (0.045, 0.012, "left", "bottom"),
        "alpha": (0.05, 0.004, "left", "center"),
        "low_gamma": (0.04, 0.011, "left", "bottom"),
        "high_gamma": (-0.045, 0.004, "right", "center"),
        "delta": (0.045, -0.006, "left", "top"),
        "theta": (0.05, 0.004, "left", "bottom"),
    }
    for _, r in T.iterrows():
        dx, dy, ha, va = lab[r.band]
        ax.annotate(BRAIN_BAND_TEX_DICT[r.band], (r.rb, r.med),
                    xytext=(r.rb + dx, r.med + dy), fontsize=14, color=TIER_C[r.tier],
                    ha=ha, va=va, zorder=6, fontweight="bold")
    ax.text(0.845, 0.215, "clears\n(p<.05)", fontsize=7.3,
            color=TIER_C["clear"], ha="center", va="top", style="italic")
    ax.text((RB_MARG + RB_CLEAR) / 2, 0.215, "marginal\n(p≈.08)", fontsize=7.3,
            color=TIER_C["marginal"], ha="center", va="top", style="italic")
    ax.set_xlabel(r"cohort consistency   (rank-biserial of the gate)")
    ax.set_ylabel(r"magnitude   median $\rho_{\mathrm{sym}}$")
    ax.set_xlim(-0.45, 0.92)
    ax.set_ylim(-0.075, 0.225)
    for s_pct in (0.03, 0.10, 0.26):
        ax.scatter([], [], s=40 + s_pct * 2500, facecolor="0.55", edgecolor="0.3",
                   label=f"{s_pct*100:.0f}%")
    leg = ax.legend(loc="lower right", frameon=False, labelspacing=1.5,
                    borderpad=1.0, handletextpad=1.2, fontsize=7.5,
                    title="trace share\nof structure", title_fontsize=7.5)
    leg._legend_box.align = "left"
    ax.text(0.015, 0.97, r"$\mathbf{A}$", transform=ax.transAxes, fontsize=13,
            va="top", ha="left")


# ------------------------------------------------------------------ panel B
def panel_B(ax, g, T):
    bands = list(T.band)
    rho = g.pivot(index="patient", columns="band", values="obs_rho")[bands]
    pv = g.pivot(index="patient", columns="band", values="obs_p_one_sided")[bands]
    trace, anti = pv < TRACE_P, pv > ANTI_P
    order = sorted(rho.index, key=lambda p: (int(trace.loc[p].sum()),
                                             float(rho.loc[p].mean())), reverse=True)
    rho, trace, anti = rho.loc[order], trace.loc[order], anti.loc[order]
    M = rho.to_numpy()
    im = ax.imshow(M, cmap="coolwarm", norm=TwoSlopeNorm(0.0, -0.3, 0.9), aspect="auto")
    nrow, ncol = M.shape
    for i in range(nrow):
        for j in range(ncol):
            v = M[i, j]
            ax.text(j, i, f"{v*100:+.0f}", ha="center", va="center", fontsize=6.3,
                    color="white" if abs(v) > 0.33 else "0.15", zorder=6)
            xy = (j - 0.5, i - 0.5)
            if trace.iloc[i, j]:                       # clears null: TRACE
                ax.add_patch(Rectangle(xy, 1, 1, fill=False, ec="black", lw=1.7, zorder=5))
            elif anti.iloc[i, j]:                      # clears lower tail: ANTI
                ax.add_patch(Rectangle(xy, 1, 1, fill=False, ec="#5b1a1a", ls=":",
                                       lw=1.5, zorder=5))
            else:                                      # below threshold: NOISE
                ax.add_patch(Rectangle(xy, 1, 1, fill=False, hatch="////",
                                       ec="0.55", lw=0.0, zorder=4))
    ax.set_xticks(range(ncol)); ax.set_yticks(range(nrow))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in bands], fontsize=13)
    ax.set_yticklabels([p.replace("Pat_", "P") for p in order], fontsize=8)
    for i, c in enumerate(trace.sum(axis=1).to_numpy()):
        ax.text(ncol - 0.35, i, f"{c}", ha="left", va="center", fontsize=8,
                fontweight="bold", clip_on=False)
    ax.set_xlabel(r"band   ($\rho_{\mathrm{sym}}\times100$;  right col = # bands traced)",
                  fontsize=8.5)
    ax.text(0.02, 0.985, r"$\mathbf{B}$", transform=ax.transAxes, fontsize=13,
            va="top", ha="left", color="0.1")
    imshow_colorbar_caxdivider(im, ax, size="4%", pad=0.5)
    return order


def main():
    g = pd.read_csv(GATE)
    coh = pd.read_csv(COH).set_index("band")
    tax = pd.read_csv(TAX).set_index("band")
    T = _band_table(g, coh, tax)

    print("cohort verdict — the gate does NOT cleanly separate the bands\n")
    hdr = f"{'band':10s} {'median':>7s} {'#pos':>5s} {'gate_p':>7s} {'rank-bis':>8s} {'%struct':>7s}  tier"
    print(hdr); print("-" * len(hdr))
    for _, r in T.iterrows():
        print(f"{r.band:10s} {r.med:+7.3f} {r.n_pos:4d}/10 {r.gate_p:7.3f} {r.rb:+8.3f} "
              f"{r.share*100:6.1f}%  {r.tier}")

    fig = plt.figure(figsize=(12.8, 5.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.22, 1.0], wspace=0.22)
    panel_A(fig.add_subplot(gs[0, 0]), T)
    order = panel_B(fig.add_subplot(gs[0, 1]), g, T)
    print("\npatient order by breadth:", [p.replace('Pat_', 'P') for p in order])

    tiers = [Line2D([0], [0], marker="o", color=TIER_C["clear"], ls="", ms=10,
                    mec=TIER_C["clear"], label="clears gate (p<.05)"),
             Line2D([0], [0], marker="o", color=TIER_C["marginal"], ls="", ms=10,
                    mec=TIER_C["marginal"], label="marginal (p≈.08)"),
             Line2D([0], [0], marker="o", color="white", mec=TIER_C["fail"], ls="",
                    ms=10, label="absent / anti")]
    cells = [Patch(fc="white", ec="black", lw=1.6, label="trace (clears null)"),
             Patch(fc="white", ec="#5b1a1a", ls=":", lw=1.4, label="anti (reset)"),
             Patch(fc="white", ec="0.55", hatch="////", label="noise (below threshold)")]
    fig.legend(handles=tiers, loc="lower center", bbox_to_anchor=(0.28, -0.07),
               ncol=3, frameon=False, fontsize=8, title="panel A — band tier",
               title_fontsize=8, columnspacing=1.3, handletextpad=0.4)
    fig.legend(handles=cells, loc="lower center", bbox_to_anchor=(0.79, -0.07),
               ncol=3, frameon=False, fontsize=8, title="panel B — cell",
               title_fontsize=8, columnspacing=1.3, handletextpad=0.4)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
