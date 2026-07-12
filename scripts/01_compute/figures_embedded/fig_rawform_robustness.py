#!/usr/bin/env python3
"""Talk robustness figure — the cross-phase trace is a CONNECTIVITY property,
not a diffusion artifact.

Reads audit_176 (rho_sym matched-strength gate on 4 cophenetic distance arms) and
renders two PDFs:

  fig_robustness_propagator_vs_raw.pdf
      Per-band, per-patient matched-strength-corrected trace
      (rho_sym_obs - rho_sym_surr_p50) for the LRG propagator D=1/K(tau_min) vs a
      plain raw-connectivity distance D=-log A. Cohort median bar + MS-gate stars.
      Message: propagator and raw distance AGREE -> alpha+beta clear under both,
      low_gamma/nulls fail under both. The diffusion machinery is not load-bearing
      for the trace.

  fig_robustness_distance_matrix.pdf
      6 bands x 4 distances gate grid (gate p + CLEAR/fail). Isolates the reciprocal
      1/A as the one pathological form (fails alpha, spuriously clears low_gamma).
      The 'why we do not use 1/A' backup slide.

Source data: data/audit/rawform_robustness/{cohort_summary,per_patient_per_band}.csv
Build:       scripts/01_compute/figures_embedded/fig_rawform_robustness.py
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc_context
from matplotlib.patches import Rectangle

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.paths import DATA_ROOT, FIGURES_ROOT
from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

use_lrg_style()

SRC = DATA_ROOT / "audit" / "rawform_robustness"
OUT = FIGURES_ROOT / "robustness"
OUT.mkdir(parents=True, exist_ok=True)

# arm registry: csv-suffix -> (talk label, colour, marker)
C_PROP, C_RAW = "#1b6f8a", "#c1440e"            # diffusion teal / raw burnt-orange
ARMS_MAIN = [
    ("diff",      r"propagator  $1/K(\tau_{\min})$", C_PROP, "o"),
    ("raw_nlogA", r"raw  $-\log A$",                 C_RAW,  "D"),
]
BANDS = list(BRAIN_BANDS_NAMES)                  # slow -> fast


def _row_label(b: str) -> str:
    return BRAIN_BAND_TEX_DICT[b]


def _lab_colour(b: str) -> str:
    # darken the bright turbo midtones (alpha yellow, beta green) so labels stay legible
    return band_color(b, shade=0.72 if b in ("alpha", "beta", "low_gamma") else 0.9)


def load():
    pp = pd.read_csv(SRC / "per_patient_per_band.csv")
    gate = pd.read_csv(SRC / "cohort_summary.csv").set_index("band")
    return pp, gate


# --------------------------------------------------------------------------- #
def fig_main(pp, gate):
    """Per-band per-patient MS-corrected trace: propagator vs raw -logA."""
    rng = np.random.default_rng(0)
    off = 0.19
    with rc_context({"font.size": 12.5, "axes.labelsize": 14,
                     "xtick.labelsize": 15, "ytick.labelsize": 11.5,
                     "legend.fontsize": 12}):
        fig, ax = plt.subplots(figsize=(7.9, 4.0))
        ax.axhline(0.0, color="0.45", lw=1.0, ls=(0, (4, 3)), zorder=1)

        y_top = None
        for xi, b in enumerate(BANDS):
            sub = pp[pp.band == b]
            for (suf, _lab, col, mk), sgn in zip(ARMS_MAIN, (-1, +1)):
                x0 = xi + sgn * off
                d = (sub[f"obs_rho_{suf}"] - sub[f"surr_p50_{suf}"]).to_numpy()
                xj = x0 + rng.uniform(-0.055, 0.055, size=d.size)
                ax.scatter(xj, d, s=26, facecolor=col, edgecolor="white",
                           linewidth=0.5, marker=mk, alpha=0.9, zorder=3)
                med = float(np.median(d))
                ax.hlines(med, x0 - 0.10, x0 + 0.10, color=col, lw=3.0, zorder=4)
                # MS-gate star only where the cohort Wilcoxon clears .05
                if float(gate.loc[b, f"gate_p_{suf}"]) < 0.05:
                    ax.plot(x0, 0.99, marker="*", ms=15, color=col,
                            markeredgecolor="white", markeredgewidth=0.5,
                            transform=ax.get_xaxis_transform(), zorder=5, clip_on=False)

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([_row_label(b) for b in BANDS])
        for t, b in zip(ax.get_xticklabels(), BANDS):
            t.set_color(_lab_colour(b))
            t.set_fontweight("bold")
        ax.set_ylabel("cross-phase trace\n"
                      r"(MS-corrected $\rho_{\mathrm{sym}}$)", labelpad=6)
        ax.set_xlim(-0.55, len(BANDS) - 0.45)
        ax.margins(y=0.10)
        ax.spines[["top", "right"]].set_visible(False)

        # figure-level legend: two arms + star key
        from matplotlib.lines import Line2D
        handles = [Line2D([0], [0], marker=mk, color="none", markerfacecolor=col,
                          markeredgecolor="white", markersize=9, label=lab)
                   for (suf, lab, col, mk) in ARMS_MAIN]
        handles.append(Line2D([0], [0], marker="*", color="none", markerfacecolor="0.25",
                              markeredgecolor="white", markersize=13,
                              label=r"clears MS gate ($p{<}0.05$)"))
        fig.legend(handles=handles, loc="lower center",
                   bbox_to_anchor=(0.5, -0.01), ncol=3, frameon=False,
                   fontsize=11, handletextpad=0.3, columnspacing=1.1)
        fig.subplots_adjust(bottom=0.22, left=0.155, right=0.975, top=0.95)
        p = OUT / "fig_robustness_propagator_vs_raw.pdf"
        fig.savefig(p, transparent=True)
        plt.close(fig)
        return p


# --------------------------------------------------------------------------- #
def fig_matrix(gate):
    """6 bands x 4 distances gate grid; isolates the pathological 1/A."""
    cols = [("diff",      r"$1/K(\tau_{\min})$", "propagator"),
            ("raw_nlogA", r"$-\log A$",          "raw"),
            ("raw_1mA",   r"$1-A$",              "raw"),
            ("raw_recip", r"$1/A$",              "reciprocal")]
    gap = 0.34   # visual separation before the pathological 1/A column
    xpos = [0, 1, 2, 3 + gap]
    rows = BANDS[::-1]                              # slow at top -> fast at bottom
    C_CLEAR, C_FAIL = "#2e7d32", "0.78"
    with rc_context({"font.size": 12}):
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        for yi, b in enumerate(rows):
            ax.text(-0.95, yi, _row_label(b), ha="right", va="center",
                    fontsize=17, fontweight="bold", color=_lab_colour(b))
            for (suf, _tex, _kind), xp in zip(cols, xpos):
                p = float(gate.loc[b, f"gate_p_{suf}"])
                clear = p < 0.05
                fc = C_CLEAR if clear else C_FAIL
                ax.add_patch(Rectangle((xp - 0.44, yi - 0.44), 0.88, 0.88,
                                       facecolor=fc, edgecolor="white", lw=2,
                                       alpha=0.95 if clear else 0.55, zorder=2))
                ax.text(xp, yi + 0.10, f"{p:.3f}", ha="center", va="center",
                        color="white" if clear else "0.25",
                        fontsize=11.5, fontweight="bold" if clear else "normal", zorder=3)
                ax.text(xp, yi - 0.22, "clear" if clear else "n.s.", ha="center",
                        va="center", color="white" if clear else "0.4",
                        fontsize=8.5, zorder=3)
        # column headers
        for (suf, tex, kind), xp in zip(cols, xpos):
            ax.text(xp, len(rows) - 0.32, tex, ha="center", va="bottom", fontsize=15)
            ax.text(xp, len(rows) - 0.02, kind, ha="center", va="bottom",
                    fontsize=9.5, color="0.45", style="italic")
        # divider before the pathological reciprocal column
        xd = (xpos[2] + xpos[3]) / 2
        ax.axvline(xd, color="0.6", lw=1.1, ls=(0, (3, 3)), zorder=1)
        ax.set_xlim(-1.9, xpos[-1] + 0.6)
        ax.set_ylim(-0.7, len(rows) + 0.55)
        ax.axis("off")
        fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
        p = OUT / "fig_robustness_distance_matrix.pdf"
        fig.savefig(p, transparent=True)
        plt.close(fig)
        return p


def main():
    pp, gate = load()
    p1 = fig_main(pp, gate)
    p2 = fig_matrix(gate)
    print(f"[fig_rawform_robustness] wrote:\n  {p1}\n  {p2}")


if __name__ == "__main__":
    main()
