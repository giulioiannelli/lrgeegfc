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
Writes: data/preprint/figures/results_section2/fig_arc_a_decomposition_forest.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import band_marker, use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BASE = ROOT / "data/audit/consolidation_arc_rhosym"
OUT = ROOT / "data/preprint/figures/results_section2/fig_arc_a_decomposition_forest.pdf"

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


def _verdict_color(v: float, p: float) -> str:
    """Green if the patient clears its own trace null; red if it reverts (reset,
    T < 0); grey if weak / undetermined."""
    if p < 0.05:
        return C_TRACE
    if v < 0:
        return C_ANTI
    return C_UNDET


def draw_panel(ax, func, val_col, p_col, order, ypos, cohort, per_pat, per_null,
               dot_s=88, med_s=185, p_fs=10.5):
    """Per-band violin forest (matches the §1 band figure): gate-tinted violin (green
    when the band clears its matched-strength null), per-band marker glyph for the
    patient dots and the cohort median, dot colour by verdict (reset -> red)."""
    ax.axvline(0.0, color="0.55", lw=1.1, ls="--", zorder=1)
    cv = cohort[cohort.functional == func].set_index("band")
    for b in order:
        y0 = ypos[b]
        gate = float(cv.loc[b, "wilcoxon_p"])
        gcol = GATE_PASS if gate < 0.05 else GATE_FAIL
        mk = band_marker(b)
        pv = per_pat[per_pat.band == b].set_index("patient")
        nl = per_null[per_null.band == b].set_index("patient")
        pats = [pat for pat in COHORT if pat in pv.index]
        vals = np.array([float(pv.loc[pat, val_col]) for pat in pats])
        vals = vals[np.isfinite(vals)]
        # gate-tinted violin: green fill + crisp outline when the band passes, grey else
        if vals.size >= 4:
            parts = ax.violinplot([vals], positions=[y0], vert=False, widths=0.86,
                                  showextrema=False, bw_method=0.55)
            for body in parts["bodies"]:
                body.set_facecolor(gcol); body.set_edgecolor("none")
                body.set_alpha(0.15); body.set_zorder(2)
                verts = body.get_paths()[0].vertices
                ax.plot(verts[:, 0], verts[:, 1], color=gcol, lw=1.0, alpha=0.6,
                        zorder=2.3, solid_joinstyle="round")
        # per-patient dots: per-band marker, verdict colour (reset -> red)
        for pat in pats:
            v = float(pv.loc[pat, val_col])
            if not np.isfinite(v):
                continue
            ax.scatter(v, y0 + 0.66 * Y_OFF[pat], s=dot_s, marker=mk,
                       color=_verdict_color(v, float(nl.loc[pat, p_col])),
                       edgecolor="white", linewidth=0.8, zorder=4)
        # cohort median: gate-coloured stem from 0 + per-band marker + gate p
        med = float(cv.loc[b, "med_obs"])
        ax.plot([0.0, med], [y0, y0], color=gcol, lw=3.0, zorder=5,
                solid_capstyle="round")
        ax.scatter(med, y0, s=med_s, marker=mk, color=gcol, edgecolor="white",
                   linewidth=1.3, zorder=6)
        ax.text(med, y0 + 0.45, f"$p={gate:.3f}$", color=gcol, fontsize=p_fs,
                ha="center", va="center",
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

    # panel letters are supplied by the LaTeX mosaic in results_sec_2.tex

    handles = [
        Line2D([0], [0], marker="o", ls="", mfc=C_TRACE, mec="white", ms=11,
               label="clears own null (trace)"),
        Line2D([0], [0], marker="o", ls="", mfc=C_ANTI, mec="white", ms=11,
               label=r"reset ($T<0$)"),
        Line2D([0], [0], marker="o", ls="", mfc=C_UNDET, mec="white", ms=11,
               label="undetermined"),
        Line2D([0], [0], marker="o", ls="-", color=GATE_PASS, mfc=GATE_PASS,
               mec="white", ms=12, lw=3.0, label="cohort median + gate $p$"),
        Patch(facecolor=GATE_PASS, edgecolor=GATE_PASS, alpha=0.20,
              label="passing-band distribution"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=5, frameon=False, fontsize=11, handletextpad=0.5,
               columnspacing=1.5)

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
