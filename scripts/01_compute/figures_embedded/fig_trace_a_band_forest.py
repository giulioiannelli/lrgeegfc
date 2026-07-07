#!/usr/bin/env python3
r"""fig:trace_a — per-band per-patient trace forest (Results §1, paragraph a).

CORE MESSAGE: only beta and alpha hold the cophenetic trace at the cohort level,
and they win on CONSISTENCY, not peak amplitude. Each dot is a patient's rho_sym,
coloured by whether it clears its OWN strength-matched null. The cohort median is a
lollipop coloured by the cohort gate (green = trace, p<0.05). low_gamma carries the
single strongest tracers (Pat_05 +0.86 > beta's max +0.54) but is offset by anti
patients and fails; theta sits negative.

NOTE on the gate: the Wilcoxon gate p uses ALL 10 patients' (obs - surrogate),
including the ones individually indistinguishable from null — dropping them would be
selection bias. Per-patient own-null bars are shown in the companion 6-panel forest.

Two styles (same data, same median + gate p; the forest is the default):
    python fig_trace_a_band_forest.py                 # -> fig_trace_a_band_forest.pdf
    python fig_trace_a_band_forest.py --style violin  # -> fig_trace_a_band_violin.pdf
The violin overlays the per-band distribution (KDE of the 10 patients) behind the
same dots + median lollipop.

Reads : data/audit/rho_sym_gate/per_patient_per_band.csv
        data/audit/rho_sym_gate/cohort_summary.csv
Writes: data/reports/results_section1/fig_trace_a_band_{forest,violin}.pdf
"""
from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

GATE = ROOT / "data/audit/rho_sym_gate/per_patient_per_band.csv"
SUMMARY = ROOT / "data/audit/rho_sym_gate/cohort_summary.csv"
OUT = {
    "forest": ROOT / "data/reports/results_section1/fig_trace_a_band_forest.pdf",
    "violin": ROOT / "data/reports/results_section1/fig_trace_a_band_violin.pdf",
}

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
Y_OFF = dict(zip(COHORT, np.linspace(-0.30, 0.30, len(COHORT))))

C_TRACE, C_UNDET, C_ANTI = "#3a9a4f", "#b9b9b9", "#d1352b"
GATE_PASS, GATE_FAIL = "#1f7a34", "#565656"
TIER_C = {"beta": GATE_PASS, "alpha": GATE_PASS, "low_gamma": GATE_FAIL,
          "high_gamma": GATE_FAIL, "delta": GATE_FAIL, "theta": GATE_FAIL}


def _verdict_color(p):
    if p < 0.05:
        return C_TRACE
    if p > 0.95:
        return C_ANTI
    return C_UNDET


def _add_violins(ax, g, s, order, ypos):
    """Per-band KDE violin (horizontal), tinted by the cohort gate verdict."""
    data = [g[g.band == b].obs_rho.to_numpy(float) for b in order]
    positions = [ypos[b] for b in order]
    parts = ax.violinplot(data, positions=positions, vert=False, widths=0.74,
                          showextrema=False, showmedians=False)
    for b, body in zip(order, parts["bodies"]):
        c = GATE_PASS if s.loc[b, "gate_p_sym"] < 0.05 else GATE_FAIL
        body.set_facecolor(c)
        body.set_edgecolor(c)
        body.set_alpha(0.20)
        body.set_zorder(2)


def _add_patient_dots(ax, g, order, ypos):
    for b in order:
        y0 = ypos[b]
        for _, r in g[g.band == b].iterrows():
            ax.scatter(r.obs_rho, y0 + Y_OFF[r.patient], s=125,
                       color=_verdict_color(r.obs_p_one_sided),
                       edgecolor="white", linewidth=0.9, zorder=4)


def _add_median_and_p(ax, s, order, ypos):
    for b in order:
        y0 = ypos[b]
        med = s.loc[b, "obs_median_rho_sym"]
        p = s.loc[b, "gate_p_sym"]
        mcol = GATE_PASS if p < 0.05 else GATE_FAIL
        ax.plot([0.0, med], [y0, y0], color=mcol, lw=3.2, zorder=5, solid_capstyle="round")
        ax.scatter(med, y0, s=185, color=mcol, edgecolor="white", linewidth=1.3, zorder=6)
        # p-values in a tidy column hugging the zero line, in the clear gap above each row
        ax.text(0.02, y0 + 0.42, f"$p={p:.3f}$", color=mcol, fontsize=12.5,
                ha="left", va="center", fontweight="bold" if p < 0.05 else "normal")


def main(style="forest"):
    g = pd.read_csv(GATE)
    s = pd.read_csv(SUMMARY).set_index("band")

    order = s.sort_values("obs_median_rho_sym", ascending=False).index.tolist()
    ypos = {b: len(order) - 1 - i for i, b in enumerate(order)}

    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    ax.axvline(0.0, color="0.55", lw=1.1, ls="--", zorder=1)

    if style == "violin":
        _add_violins(ax, g, s, order, ypos)
    _add_patient_dots(ax, g, order, ypos)
    _add_median_and_p(ax, s, order, ypos)

    ax.set_yticks(list(ypos.values()))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in order], fontsize=20)
    for tick, b in zip(ax.get_yticklabels(), order):
        tick.set_color(TIER_C[b])
    ax.set_ylim(-0.75, len(order) - 0.25)
    ax.set_xlim(-0.52, 0.97)
    ax.set_xlabel(r"held trace  $\rho^{\mathrm{coph}}$   "
                  r"(rest$_{\mathrm{post}}$ vs task$_{\mathrm{test}}$)", fontsize=16)
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", length=0)
    ax.text(0.0, 1.02, r"$\mathbf{a}$", transform=ax.transAxes, fontsize=17,
            va="bottom", ha="left")

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
    if style == "violin":
        handles.append(Patch(facecolor=GATE_PASS, edgecolor=GATE_PASS, alpha=0.20,
                             label="per-band distribution (KDE)"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=2 if style == "forest" else 3, frameon=False, fontsize=11.5,
               handletextpad=0.5, columnspacing=1.8)

    fig.subplots_adjust(left=0.12, right=0.96, top=0.95, bottom=0.17)
    OUT[style].parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT[style], bbox_inches="tight", transparent=True)
    plt.close(fig)

    print(f"fig:trace_a — band {style}\n")
    for b in order:
        gb = g[g.band == b]
        nt = int((gb.obs_p_one_sided < 0.05).sum())
        na = int((gb.obs_p_one_sided > 0.95).sum())
        print(f"{b:11s} med={s.loc[b,'obs_median_rho_sym']:+.3f} "
              f"gate_p={s.loc[b,'gate_p_sym']:.3f} trace={nt} anti={na}")
    print(f"\nwrote {OUT[style]}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--style", choices=["forest", "violin"], default="forest",
                    help="forest (default) or violin variant")
    args = ap.parse_args()
    main(args.style)
