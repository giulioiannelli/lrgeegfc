#!/usr/bin/env python3
r"""fig:arc_a — the encoding-vs-inference decomposition forest (Results §2, R2.1).

CORE MESSAGE (the flagship split): the persistent reorganization decomposes into an
ENCODING echo (what persists of the premises the patient was SHOWN, T_learn = corr(e,p))
and an INFERENCE-SPECIFIC component (what persists of the relations the patient had to
REASON out, encoding held constant, T_infspec.e = corr(f,p|e)). Both persist, in
different rhythms:
  - encoding echo persists in alpha, beta AND delta (each clears its own matched-strength
    null: p = 0.014 / 0.032 / 0.042);
  - the inference-specific component persists in BETA ALONE (p = 0.010; next-smallest
    band p = 0.25).
So beta keeps not only the seen premises but a component tied specifically to the
inferred order.

Each dot is one patient's rho_sym-based T, coloured by whether it clears that patient's
OWN strength-matched null. The cohort median is a lollipop coloured by the cohort gate
(green = clears matched-strength, p < 0.05). Same band order in both panels so the
collapse from broad (encoding) to beta-only (inference) is read left-to-right.

Reads : data/audit/consolidation_arc_rhosym/arc_cohort_verdict.csv
        data/audit/consolidation_arc_rhosym/arc_per_patient.csv
        data/audit/consolidation_arc_rhosym/arc_null_per_patient.csv
Writes: data/reports/results_section2/fig_arc_a_decomposition_forest.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BASE = ROOT / "data/audit/consolidation_arc_rhosym"
OUT = ROOT / "data/reports/results_section2/fig_arc_a_decomposition_forest.pdf"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
Y_OFF = dict(zip(COHORT, np.linspace(-0.30, 0.30, len(COHORT))))

C_TRACE, C_UNDET, C_ANTI = "#3a9a4f", "#b9b9b9", "#d1352b"
GATE_PASS, GATE_FAIL = "#1f7a34", "#565656"

# the two components: (functional column, per-patient value col, own-null p col, label)
PANELS = [
    ("T_learn", "T_learn", "T_learn_p",
     r"encoding echo   $T_{\mathrm{learn}}=\mathrm{corr}(e,\,p)$"),
    ("T_infspec_pe", "T_infspec_pe", "T_infspec_pe_p",
     r"inference-specific   $T_{\mathrm{inf}\cdot e}=\mathrm{corr}(f,\,p\mid e)$"),
]


def _verdict_color(p: float) -> str:
    if p < 0.05:
        return C_TRACE
    if p > 0.95:
        return C_ANTI
    return C_UNDET


def draw_panel(ax, func, val_col, p_col, order, ypos, cohort, per_pat, per_null):
    ax.axvline(0.0, color="0.55", lw=1.1, ls="--", zorder=1)
    cv = cohort[cohort.functional == func].set_index("band")
    for b in order:
        y0 = ypos[b]
        # per-patient dots, coloured by own-null clearance
        pv = per_pat[per_pat.band == b].set_index("patient")
        nl = per_null[per_null.band == b].set_index("patient")
        for pat in COHORT:
            if pat not in pv.index:
                continue
            ax.scatter(pv.loc[pat, val_col], y0 + Y_OFF[pat], s=115,
                       color=_verdict_color(float(nl.loc[pat, p_col])),
                       edgecolor="white", linewidth=0.9, zorder=4)
        # cohort median lollipop, coloured by the matched-strength gate
        med = float(cv.loc[b, "med_obs"])
        gate = float(cv.loc[b, "wilcoxon_p"])
        mcol = GATE_PASS if gate < 0.05 else GATE_FAIL
        ax.plot([0.0, med], [y0, y0], color=mcol, lw=3.2, zorder=5,
                solid_capstyle="round")
        ax.scatter(med, y0, s=175, color=mcol, edgecolor="white", linewidth=1.3,
                   zorder=6)
        ax.text(0.015, y0 + 0.40, f"$p={gate:.3f}$", color=mcol, fontsize=11.5,
                ha="left", va="center",
                fontweight="bold" if gate < 0.05 else "normal")


def main():
    cohort = pd.read_csv(BASE / "arc_cohort_verdict.csv")
    per_pat = pd.read_csv(BASE / "arc_per_patient.csv")
    per_null = pd.read_csv(BASE / "arc_null_per_patient.csv")

    # shared band order: by the encoding echo median (the "what persists" ordering)
    enc = cohort[cohort.functional == "T_learn"].sort_values("med_obs", ascending=False)
    order = enc.band.tolist()
    ypos = {b: len(order) - 1 - i for i, b in enumerate(order)}

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 6.3), sharey=True)
    for ax, (func, val_col, p_col, xlab) in zip(axes, PANELS):
        draw_panel(ax, func, val_col, p_col, order, ypos, cohort, per_pat, per_null)
        ax.set_xlabel(xlab, fontsize=13)
        ax.tick_params(axis="x", labelsize=12)
        ax.spines[["top", "right"]].set_visible(False)

    # neutral band labels: each panel's lollipop colour carries its own gate verdict,
    # so the shared axis must not privilege one panel's result over the other.
    axes[0].set_yticks(list(ypos.values()))
    axes[0].set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in order], fontsize=19)
    axes[0].set_ylim(-0.75, len(order) - 0.25)
    axes[0].tick_params(axis="y", length=0)

    fig.text(0.075, 0.965, r"$\mathbf{a}$", fontsize=17, va="top", fontweight="bold")
    fig.text(0.545, 0.965, r"$\mathbf{b}$", fontsize=17, va="top", fontweight="bold")

    handles = [
        Line2D([0], [0], marker="o", ls="", mfc=C_TRACE, mec="white", ms=11,
               label="patient clears own null"),
        Line2D([0], [0], marker="o", ls="", mfc=C_UNDET, mec="white", ms=11,
               label="within null"),
        Line2D([0], [0], marker="o", ls="", mfc=C_ANTI, mec="white", ms=11,
               label="anti (reset)"),
        Line2D([0], [0], marker="o", ls="-", color=GATE_PASS, mfc=GATE_PASS,
               mec="white", ms=13, lw=3.2, label="cohort median + gate $p$"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=4, frameon=False, fontsize=11, handletextpad=0.5,
               columnspacing=1.8)

    fig.subplots_adjust(left=0.11, right=0.975, top=0.94, bottom=0.17, wspace=0.08)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:arc_a — encoding-vs-inference decomposition forest\n")
    for func in ("T_learn", "T_infspec_pe"):
        print(f"  [{func}]")
        cv = cohort[cohort.functional == func].set_index("band")
        for b in order:
            print(f"    {b:11s} med={cv.loc[b,'med_obs']:+.3f} "
                  f"gate_p={cv.loc[b,'wilcoxon_p']:.3f} "
                  f"n_above={int(cv.loc[b,'n_above'])}/10")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
