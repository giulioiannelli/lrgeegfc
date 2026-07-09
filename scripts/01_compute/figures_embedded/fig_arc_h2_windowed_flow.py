#!/usr/bin/env python3
r"""fig:arc_h2 — the trace forms in time: a windowed dynamical flow (§2, R2 dynamics).

CORE MESSAGE: dynamize the beta consolidation trace. Slide a 30 s window across the
whole recording, read the LRG cophenetic state per window, and project it onto the SAME
encoding/inference plane. The state FLOWS: baseline hovers at the origin, the task drives
it out along encoding then up the inference axis, and offline it either stays displaced
(a tracer) or relaxes back (a resetter). Two patients make the contrast: Pat_08 does the
task and its rest_post cloud stays lifted (retains the inferred order); Pat_10 makes the
same outbound journey but relaxes back toward baseline (releases it). Across the cohort,
the beta rest_post displacement tracks each patient's VERIFIED matched-strength inference
verdict (windowed vs phase-scale corr = 0.68): the six patients that clear their own null
sit displaced from baseline, the rest return.

This is a VISUALIZATION of a matched-strength-verified effect (audit_152), not a new test:
overlapping windows are autocorrelated, so no per-window significance is claimed. The
window state reduces to the cached phase FC in the full-window limit (audit_124).

Panels a,b: per-window beta state coloured by phase, with the phase-mean trajectory
(rest_pre -> task_learn -> task_test -> rest_post). Panel c: cohort beta persistence
displacement (rest_pre -> rest_post per patient), coloured by the verified own-null gate.

Reads : data/audit/windowed_consolidation_flow/{windowed_flow_per_window,windowed_flow_endpoints}.csv
        data/audit/consolidation_arc_rhosym/arc_null_per_patient.csv   (per-patient gate)
Writes: data/reports/results_section2/fig_arc_h2_windowed_flow.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

WIN = ROOT / "data/audit/windowed_consolidation_flow"
NULL = ROOT / "data/audit/consolidation_arc_rhosym/arc_null_per_patient.csv"
OUT = ROOT / "data/reports/results_section2/fig_arc_h2_windowed_flow.pdf"

VARIANT, BAND = "B", "beta"
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
PH_COL = {"rest_pre": "#9a9a9a", "task_learn": "#2c7fb8",
          "task_test": "#6a51a3", "rest_post": "#d1491c"}
PH_LAB = {"rest_pre": r"rest$_{\mathrm{pre}}$", "task_learn": r"task$_{\mathrm{learn}}$",
          "task_test": r"task$_{\mathrm{test}}$", "rest_post": r"rest$_{\mathrm{post}}$"}
C_TRACE, C_UNDET = "#3a9a4f", "#b9b9b9"


def _arrow(ax, p0, p1, color, lw, alpha=1.0, ms=12, zorder=6):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, lw=lw,
                                 color=color, alpha=alpha, shrinkA=1, shrinkB=2,
                                 zorder=zorder, capstyle="round"))


def draw_exemplar(ax, pw, ep, pat, tag):
    ax.axhspan(0.0, 0.42, color="#eaf3ea", zorder=0)
    ax.axhline(0.0, color="0.6", lw=0.9, zorder=1)
    ax.axvline(0.0, color="0.85", lw=0.8, zorder=1)
    d = pw[pw.patient == pat]
    for ph in PHASES:                                    # per-window cloud, by phase
        s = d[d.phase == ph]
        ax.scatter(s.x, s.y, s=26, color=PH_COL[ph], alpha=0.32, edgecolor="none",
                   zorder=3)
    means = {ph: (ep[(ep.patient == pat) & (ep.phase == ph)].x_mean.values[0],
                  ep[(ep.patient == pat) & (ep.phase == ph)].y_mean.values[0])
             for ph in PHASES}
    for a, b in zip(PHASES[:-1], PHASES[1:]):            # trajectory of phase means
        _arrow(ax, means[a], means[b], "0.25", 2.4, alpha=0.9)
    for ph in PHASES:
        ax.scatter(*means[ph], s=180, color=PH_COL[ph], edgecolor="white",
                   linewidth=1.6, zorder=8)
    # emphasize the offline landing
    ax.annotate(PH_LAB["rest_post"], means["rest_post"], textcoords="offset points",
                xytext=(8, 8), fontsize=11, color=PH_COL["rest_post"], fontweight="bold")
    ax.set_xlim(-0.17, 0.52)
    ax.set_ylim(-0.17, 0.42)
    ax.set_xlabel(r"encoding  $\rho_S(\Delta, e)$", fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(tag, fontsize=13, fontweight="bold", pad=6)


def main():
    pw = pd.read_csv(WIN / "windowed_flow_per_window.csv")
    ep = pd.read_csv(WIN / "windowed_flow_endpoints.csv")
    nb = pd.read_csv(NULL)
    pw = pw[(pw.variant == VARIANT) & (pw.band == BAND)]
    ep = ep[(ep.variant == VARIANT) & (ep.band == BAND)]
    gate = nb[nb.band == BAND].set_index("patient")["T_infspec_pe_p"]

    fig, axes = plt.subplots(1, 3, figsize=(15.2, 5.4),
                             gridspec_kw=dict(width_ratios=[1.0, 1.0, 1.12]))
    draw_exemplar(axes[0], pw, ep, "Pat_08",
                  r"Pat_08 — does the task, then $\bf{retains}$")
    draw_exemplar(axes[1], pw, ep, "Pat_10",
                  r"Pat_10 — same journey, then $\bf{relaxes\ back}$")
    axes[0].set_ylabel(r"inference-specific  $\rho_S(\Delta, f\!\mid\! e)$", fontsize=12)

    # -------- panel c: cohort beta persistence displacement, gated ------------
    axc = axes[2]
    axc.axhspan(0.0, 0.34, color="#eaf3ea", zorder=0)
    axc.axhline(0.0, color="0.6", lw=0.9, zorder=1)
    axc.axvline(0.0, color="0.85", lw=0.8, zorder=1)
    axc.text(0.02, 0.325, "inference-specific retained", fontsize=10,
             color="#3a7a44", va="top", ha="left", style="italic")
    n_tr = 0
    for pat in COHORT:
        pre = ep[(ep.patient == pat) & (ep.phase == "rest_pre")]
        post = ep[(ep.patient == pat) & (ep.phase == "rest_post")]
        if not len(pre) or not len(post):
            continue
        p0 = (pre.x_mean.values[0], pre.y_mean.values[0])
        p1 = (post.x_mean.values[0], post.y_mean.values[0])
        g = float(gate.get(pat, np.nan))
        col = C_TRACE if g < 0.05 else C_UNDET
        n_tr += int(g < 0.05)
        _arrow(axc, p0, p1, col, 1.9, alpha=0.7, zorder=5)
        axc.scatter(*p0, s=32, color="0.6", edgecolor="white", linewidth=0.5, zorder=4)
        axc.scatter(*p1, s=150, color=col, edgecolor="white", linewidth=1.2, zorder=7)
        if pat in ("Pat_08", "Pat_10"):                  # tie back to the exemplars only
            axc.annotate(pat.replace("Pat_", "P"), p1, textcoords="offset points",
                         xytext=(8, -2), fontsize=9.5, color=col, fontweight="bold")
    axc.set_xlim(-0.05, 0.36)
    axc.set_ylim(-0.05, 0.34)
    axc.set_xlabel(r"encoding  $\rho_S(\Delta, e)$", fontsize=12)
    axc.spines[["top", "right"]].set_visible(False)
    axc.set_title(rf"cohort $\beta$: rest$_{{\mathrm{{pre}}}}\!\to\!$rest$_{{\mathrm{{post}}}}$"
                  rf"  ({n_tr}/10 clear own null)", fontsize=12.5, fontweight="bold", pad=6)

    handles = [
        Line2D([0], [0], marker="o", ls="", mfc=PH_COL["rest_pre"], mec="white", ms=10,
               label=r"rest$_{\mathrm{pre}}$"),
        Line2D([0], [0], marker="o", ls="", mfc=PH_COL["task_learn"], mec="white", ms=10,
               label=r"task$_{\mathrm{learn}}$"),
        Line2D([0], [0], marker="o", ls="", mfc=PH_COL["task_test"], mec="white", ms=10,
               label=r"task$_{\mathrm{test}}$"),
        Line2D([0], [0], marker="o", ls="", mfc=PH_COL["rest_post"], mec="white", ms=10,
               label=r"rest$_{\mathrm{post}}$  (each dot = one 30 s window)"),
        Line2D([0], [0], marker="o", ls="", mfc=C_TRACE, mec="white", ms=10,
               label=r"c: patient clears own inference null"),
        Line2D([0], [0], marker="o", ls="", mfc=C_UNDET, mec="white", ms=10,
               label=r"c: within own null"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=3, frameon=False, fontsize=10.2, handletextpad=0.5, columnspacing=1.8)

    fig.subplots_adjust(left=0.055, right=0.99, top=0.9, bottom=0.19, wspace=0.20)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:arc_h2 — windowed dynamical flow (beta, variant B)\n")
    for pat in ("Pat_08", "Pat_10"):
        row = " -> ".join(
            f"{ph[:4]}({ep[(ep.patient==pat)&(ep.phase==ph)].x_mean.values[0]:+.2f},"
            f"{ep[(ep.patient==pat)&(ep.phase==ph)].y_mean.values[0]:+.2f})" for ph in PHASES)
        print(f"  {pat}: {row}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
