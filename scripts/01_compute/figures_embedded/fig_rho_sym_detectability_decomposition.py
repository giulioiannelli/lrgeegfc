#!/usr/bin/env python3
r"""Detectability decomposition of the cophenet trace — why "no-trace" is the wrong label.

The cohort gate asks "is the trace cohort-CONSISTENT?" and conflates two very
different reasons a band can fail: patients showing ANTI (reset) vs patients
below the detection floor (UNDETERMINED). This figure separates them. Each band's
10 patients are split by whether they clear their OWN matched-strength null:

  TRACE        obs p < .05  (rho_sym above own strength null)
  ANTI         obs p > .95  (below own null — significant reset)
  undetermined otherwise    (below the detection floor — no verdict)

Result: every band except theta leans TRACE among detectable patients (67-100%);
low-gamma / high-gamma / delta miss the cohort gate on UNDETERMINED patients, not
anti ones — so "no trace" overstates; the honest label is "subset trace, not
cohort-consistent". theta is the only genuinely balanced band (2 trace / 2 anti /
6 undetermined) — "absent/undetermined", not "anti". Overall 29 trace vs 9 anti
(p=8e-4) vs ~3+3 expected under the null: the per-patient traces are real and
directional. Caveat printed to stdout: per-band trace>anti is NOT individually
significant for gamma/delta (small n) — cohort-consistency (beta/alpha) stays the
stronger claim; "subset trace" != "cohort trace".

Reads  : data/audit/rho_sym_gate/per_patient_per_band.csv
Writes : data/reports/rho_sym_band_map/fig_rho_sym_detectability_decomposition.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import binomtest

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

GATE = ROOT / "data/audit/rho_sym_gate/per_patient_per_band.csv"
OUT = ROOT / "data/reports/rho_sym_band_map/fig_rho_sym_detectability_decomposition.pdf"

C_TRACE, C_UNDET, C_ANTI = "#3a9a4f", "#d0d0d0", "#d1352b"
# band tier -> label colour (clear / marginal-subset / absent)
TIER_C = {"clear": "#1a1a1a", "marginal": "#c07d1f", "absent": "#8a8a8a"}
TIER = {"beta": "clear", "alpha": "clear", "low_gamma": "marginal",
        "high_gamma": "marginal", "delta": "marginal", "theta": "absent"}


def _counts(g):
    rows = []
    for b, s in g.groupby("band"):
        nt = int((s.obs_p_one_sided < 0.05).sum())
        na = int((s.obs_p_one_sided > 0.95).sum())
        nu = len(s) - nt - na
        det = nt + na
        rows.append(dict(band=b, trace=nt, undet=nu, anti=na,
                         rate=(nt / det if det else np.nan),
                         p=(binomtest(nt, det, 0.5, alternative="greater").pvalue
                            if det else np.nan)))
    df = pd.DataFrame(rows).sort_values("trace", ascending=True).reset_index(drop=True)
    return df


def main():
    g = pd.read_csv(GATE)
    T = _counts(g)

    print("detectability decomposition (of 10 patients per band)\n")
    for _, r in T.sort_values("trace", ascending=False).iterrows():
        star = "" if np.isnan(r.p) else ("  (trace>anti p=%.3f)" % r.p)
        print(f"  {r.band:11s} trace={r.trace}  undet={r.undet}  anti={r.anti}  "
              f"trace|detect={r.rate:.0%}{star}")
    print("\nCAVEAT: per-band trace>anti is NOT individually significant for "
          "gamma/delta (small n); cohort-consistency (beta/alpha) stays the "
          "stronger claim. 'subset trace' != 'cohort trace'.")

    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    y = np.arange(len(T))
    ax.barh(y, T.trace, color=C_TRACE, edgecolor="white", lw=0.8, label="trace (clears own null)")
    ax.barh(y, T.undet, left=T.trace, color=C_UNDET, edgecolor="white", lw=0.8,
            label="undetermined (below floor)")
    ax.barh(y, T.anti, left=T.trace + T.undet, color=C_ANTI, edgecolor="white", lw=0.8,
            label="anti / reset (clears lower tail)")

    for i, r in T.iterrows():
        if r.trace:
            ax.text(r.trace / 2, i, f"{r.trace}", ha="center", va="center",
                    color="white", fontsize=9, fontweight="bold")
        if r.anti:
            ax.text(r.trace + r.undet + r.anti / 2, i, f"{r.anti}", ha="center",
                    va="center", color="white", fontsize=9, fontweight="bold")
        # trace|detectable rate at the right margin
        if not np.isnan(r.rate):
            ax.text(10.25, i, f"{r.rate:.0%} of\ndetectable", ha="left", va="center",
                    fontsize=6.8, color=C_TRACE if r.rate > 0.5 else "0.4")

    ax.set_yticks(y)
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in T.band], fontsize=14)
    for tick, b in zip(ax.get_yticklabels(), T.band):
        tick.set_color(TIER_C[TIER[b]])
    ax.set_xlim(0, 10)
    ax.set_xticks(range(0, 11, 2))
    ax.set_xlabel("patients (of 10)")
    ax.set_ylim(-0.6, len(T) - 0.4)
    # cohort-tier annotation on the far right
    ax.text(13.2, len(T) - 1, "cohort-consistent\n(clears gate)", fontsize=6.8,
            color=TIER_C["clear"], ha="center", va="center", style="italic")
    ax.text(13.2, 1.5, "subset trace\n(marginal gate)", fontsize=6.8,
            color=TIER_C["marginal"], ha="center", va="center", style="italic")
    ax.text(13.2, 0, "absent /\nundetermined", fontsize=6.8,
            color=TIER_C["absent"], ha="center", va="center", style="italic")

    handles = [Patch(fc=C_TRACE, label="trace (clears own null)"),
               Patch(fc=C_UNDET, label="undetermined (below floor)"),
               Patch(fc=C_ANTI, label="anti / reset")]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.30),
              ncol=3, frameon=False, fontsize=8, handletextpad=0.5, columnspacing=1.4)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
