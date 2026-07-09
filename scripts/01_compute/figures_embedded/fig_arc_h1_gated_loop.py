#!/usr/bin/env python3
r"""fig:arc_h1 — online vs offline: every rhythm reasons, only beta keeps it (§2, R2.1/R2.2).

CORE MESSAGE (the disentangling in two moves):
  (a) ONLINE. During the task every rhythm makes the same journey through the
      encoding/inference plane -- out along the encoding axis (rest_pre -> task_learn),
      then UP the inference-specific axis (task_learn -> task_test). All six bands reach
      the inference peak: reasoning online is universal.
  (b) OFFLINE. What persists is not. At rest_post the loops relax back -- and only BETA
      lands still lifted on the inference axis (clears its matched-strength null,
      p = 0.010). Alpha / delta relax back ONTO the encoding axis (y ~ 0): they keep the
      seen premises (encoding echo clears, p = 0.014 / 0.042) but shed the inferred order.
      Low-gamma / theta / high-gamma keep neither.

Panel a shows the shared outbound journey (rest_pre -> task_learn -> task_test), full
scale, so the universal rise is visible. Panel b zooms to the rest_post LANDING, where
the band-specific result lives; landings sit EXACTLY on the verified rho_sym
(T_learn, T_infspec_pe) [audit_163 cross-check = machine precision]. Per-patient beta
dots coloured by each patient's own matched-strength inference null; the shaded band is
the inference-retained zone (y > 0).

Reads : data/audit/consolidation_arc_rhosym/arc_phase_trajectory.csv   (audit_163)
        data/audit/consolidation_arc_rhosym/arc_cohort_verdict.csv      (gate p)
        data/audit/consolidation_arc_rhosym/arc_null_per_patient.csv    (per-patient gate)
Writes: data/reports/results_section2/fig_arc_h1_gated_loop.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BASE = ROOT / "data/audit/consolidation_arc_rhosym"
OUT = ROOT / "data/reports/results_section2/fig_arc_h1_gated_loop.pdf"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
ORDER = ["beta", "alpha", "delta", "low_gamma", "theta", "high_gamma"]
OUTBOUND = ["rest_pre", "task_learn", "task_test"]

# qualitative band palette — NO greens (green is reserved for the inference gate)
C_BAND = {"beta": "#d1491c", "alpha": "#7b3294", "delta": "#1a6fbf",
          "low_gamma": "#0d8f8f", "theta": "#c9a227", "high_gamma": "#9a6a2f"}
C_TRACE, C_UNDET, C_ANTI = "#3a9a4f", "#b9b9b9", "#d1352b"
GATE_PASS, GATE_GREY = "#1f7a34", "#8a8a8a"
# explicit rest_post label offsets (points) to avoid collisions near the axis
LBL_OFF = {"beta": (11, 6), "alpha": (11, -15), "delta": (2, 11),
           "low_gamma": (2, -17), "theta": (-9, -17), "high_gamma": (-33, 3)}


def _verdict_color(p: float) -> str:
    if not np.isfinite(p):
        return C_UNDET
    return C_TRACE if p < 0.05 else (C_ANTI if p > 0.95 else C_UNDET)


def _arrow(ax, p0, p1, color, lw, alpha=1.0, ms=13, zorder=4):
    ax.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle="-|>", mutation_scale=ms, lw=lw, color=color,
        alpha=alpha, shrinkA=0, shrinkB=2, zorder=zorder, capstyle="round"))


def main():
    traj = pd.read_csv(BASE / "arc_phase_trajectory.csv")
    verd = pd.read_csv(BASE / "arc_cohort_verdict.csv")
    nullp = pd.read_csv(BASE / "arc_null_per_patient.csv")
    med = (traj.groupby(["band", "phase"])[["x", "y"]].median()
           .reset_index().set_index(["band", "phase"]))

    def gate_p(band, func):
        r = verd[(verd.band == band) & (verd.functional == func)]
        return float(r.wilcoxon_p.values[0]) if len(r) else np.nan

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.6, 6.7),
                                   gridspec_kw=dict(width_ratios=[1.0, 1.18]))

    # ================= panel a — ONLINE (shared outbound journey) =============
    axA.axhspan(0.0, 0.72, color="#f4f1e8", zorder=0)
    axA.axhline(0.0, color="0.6", lw=1.0, zorder=1)
    axA.axvline(0.0, color="0.85", lw=0.8, zorder=1)
    for band in ORDER:                                  # identity not needed here: all rise
        col = C_BAND[band]
        pts = {ph: (med.loc[(band, ph), "x"], med.loc[(band, ph), "y"]) for ph in OUTBOUND}
        _arrow(axA, pts["rest_pre"], pts["task_learn"], col, 1.6, alpha=0.45)
        _arrow(axA, pts["task_learn"], pts["task_test"], col, 2.3, alpha=0.9)
        axA.scatter(*pts["task_learn"], s=30, facecolor="white", edgecolor=col,
                    linewidth=1.5, zorder=5)
        axA.scatter(*pts["task_test"], s=88, color=col, edgecolor="white",
                    linewidth=1.2, zorder=6)
    axA.scatter(0, 0, s=60, color="0.45", zorder=7)
    axA.annotate(r"rest$_{\mathrm{pre}}$", (0, 0), textcoords="offset points",
                 xytext=(8, -3), fontsize=10, color="0.4")
    axA.annotate("all six rhythms\nrise to the\ninference peak", (0.33, 0.55),
                 xytext=(-0.24, 0.52), fontsize=11.5, color="0.25", style="italic",
                 ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color="0.55", lw=1.2,
                                 connectionstyle="arc3,rad=-0.2"))
    axA.annotate(r"task$_{\mathrm{learn}}$ (encoding pole)", (0.60, -0.13),
                 xytext=(0.60, -0.26), fontsize=10, color="0.4", ha="center", va="top",
                 arrowprops=dict(arrowstyle="->", color="0.6", lw=1.0))
    axA.set_xlim(-0.34, 0.82)
    axA.set_ylim(-0.32, 0.72)
    axA.set_xlabel(r"encoding alignment  $\rho_S(\Delta,\,e)$", fontsize=12.5)
    axA.set_ylabel(r"inference-specific alignment  $\rho_S(\Delta,\,f\!\mid\! e)$", fontsize=12.5)
    axA.spines[["top", "right"]].set_visible(False)
    axA.text(0.0, 1.02, r"$\mathbf{a}$   online — every rhythm reasons",
             transform=axA.transAxes, fontsize=14, va="bottom", fontweight="bold")

    # ================= panel b — OFFLINE (rest_post landing, zoomed) ==========
    yb_hi = 0.44
    axB.axhspan(0.0, yb_hi, color="#eaf3ea", zorder=0)
    axB.axhline(0.0, color="0.6", lw=1.0, zorder=1)
    axB.axvline(0.0, color="0.85", lw=0.8, zorder=1)
    axB.text(0.02, yb_hi - 0.01, "inference-specific retained", fontsize=10,
             color="#3a7a44", va="top", ha="left", style="italic")

    # per-patient beta dots (the inference tracer band), coloured by own null
    ppb = traj[(traj.band == "beta") & (traj.phase == "rest_post")].set_index("patient")
    nb = nullp[nullp.band == "beta"].set_index("patient")
    for pat in COHORT:
        if pat not in ppb.index:
            continue
        pc = _verdict_color(float(nb.loc[pat, "T_infspec_pe_p"])) if pat in nb.index else C_UNDET
        axB.scatter(ppb.loc[pat, "x"], ppb.loc[pat, "y"], s=52, color=pc,
                    edgecolor="white", linewidth=0.7, zorder=5, alpha=0.85)

    axB.scatter(0, 0, s=60, color="0.45", zorder=7)
    axB.annotate(r"rest$_{\mathrm{pre}}$", (0, 0), textcoords="offset points",
                 xytext=(-3, -15), fontsize=10, color="0.4", ha="center")
    for band in ORDER:
        x, y = med.loc[(band, "rest_post"), "x"], med.loc[(band, "rest_post"), "y"]
        col = C_BAND[band]
        inf_gated = gate_p(band, "T_infspec_pe") < 0.05
        enc_gated = gate_p(band, "T_learn") < 0.05
        _arrow(axB, (0, 0), (x, y), col, 2.0, alpha=0.5, zorder=4)
        # fill = inference gate (green / grey); edge = encoding gate (band ring / white)
        fill = GATE_PASS if inf_gated else GATE_GREY
        edge = col if enc_gated else "white"
        axB.scatter(x, y, s=(285 if inf_gated else 180), marker="D", color=fill,
                    edgecolor=edge, linewidth=(2.6 if enc_gated else 1.4), zorder=8)
        star = r"$^\star$" if inf_gated else ""
        axB.annotate(rf"{BRAIN_BAND_TEX_DICT[band]}{star}", (x, y),
                     textcoords="offset points", xytext=LBL_OFF[band], fontsize=14,
                     color=col, fontweight="bold")
    axB.set_xlim(-0.34, 0.60)
    axB.set_ylim(-0.20, yb_hi)
    axB.set_xlabel(r"encoding alignment  $=\ T_{\mathrm{learn}}$", fontsize=12.5)
    axB.set_ylabel(r"inference-specific  $=\ T_{\mathrm{inf}\cdot e}$", fontsize=12.5)
    axB.spines[["top", "right"]].set_visible(False)
    axB.text(0.0, 1.02, r"$\mathbf{b}$   offline — only $\beta$ keeps the inferred order",
             transform=axB.transAxes, fontsize=14, va="bottom", fontweight="bold")
    axB.text(0.985, 0.03,
             r"$\star$ clears matched-strength inference null" "\n"
             r"$\beta$: $p=0.010$   ·   encoding echo $\alpha/\beta/\delta$: $p=0.014/0.032/0.042$",
             transform=axB.transAxes, fontsize=9.6, va="bottom", ha="right", color="0.2",
             bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.85", alpha=0.92))

    handles = [
        Line2D([0], [0], marker="D", ls="", mfc=GATE_PASS, mec="white", ms=13,
               label=r"inference-gated (green fill) — $\beta$ only"),
        Line2D([0], [0], marker="D", ls="", mfc=GATE_GREY, mec=C_BAND["alpha"], mew=2.4,
               ms=11, label=r"encoding-gated (coloured ring) — $\alpha/\beta/\delta$"),
        Line2D([0], [0], marker="D", ls="", mfc=GATE_GREY, mec="white", ms=11,
               label="ungated"),
        Line2D([0], [0], marker="o", ls="", mfc=C_TRACE, mec="white", ms=9,
               label=r"$\beta$ patient clears own null"),
        Line2D([0], [0], marker="o", ls="", mfc=C_UNDET, mec="white", ms=9,
               label=r"$\beta$ patient within null"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.055),
               ncol=3, frameon=False, fontsize=10.0, handletextpad=0.5,
               columnspacing=1.7)

    fig.subplots_adjust(left=0.065, right=0.985, top=0.92, bottom=0.16, wspace=0.24)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:arc_h1 — online vs offline consolidation plane\n")
    for band in ORDER:
        print(f"  {band:11s} test_y={med.loc[(band,'task_test'),'y']:+.2f}  "
              f"post=({med.loc[(band,'rest_post'),'x']:+.3f},{med.loc[(band,'rest_post'),'y']:+.3f})  "
              f"enc_p={gate_p(band,'T_learn'):.3f} inf_p={gate_p(band,'T_infspec_pe'):.3f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
