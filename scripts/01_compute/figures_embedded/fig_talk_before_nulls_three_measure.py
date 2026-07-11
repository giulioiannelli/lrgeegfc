#!/usr/bin/env python3
r"""fig (talk slide 13): three reads, three answers — the observed band profile BEFORE any null.

Three horizontal-bar panels, one per read-out of the SAME connectivity, sharing
the band order (delta -> high_gamma, slow -> fast):

  * raw edges (pairwise |ImCoh| overlap)   -- fires in every band, near-flat:
                                              the edge level cannot separate a
                                              band that traces from one that does not.
  * Grassmann subspace (global modes)       -- observed subspace-rotation retention.
  * cophenetic hierarchy (multiscale)       -- observed rho_sym cross-phase retention.

Each panel is normalized to its own maximum so the SHAPE (which bands are loud) is
comparable; the native observed value is annotated on every bar. A small check marks
bands whose cohort is consistently in the trace direction BEFORE the matched-strength /
drift nulls (raw, coph: one-sided Wilcoxon of per-patient effect vs 0; Grassmann: the
observed contiguous cohort-significant run exceeds the null p95 run-length). beta is the
single band that is loudest under ALL THREE reads; the reads disagree on everything else.
That disagreement is the "who's right?" the null slide answers.

Reads : data/audit/raw_fc_matched_strength/{cohort_summary_all_bands,per_patient_per_band}.csv
        data/audit/rho_sym_gate/{cohort_summary,per_patient_per_band}.csv
        data/audit/grassmann_cluster_extent/cohort_summary.csv
Writes: data/outputs/figures/talk/fig_before_nulls_three_measure.pdf
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

ROOT = setup_script_env()
use_lrg_style()

RAW_C = ROOT / "data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv"
RAW_P = ROOT / "data/audit/raw_fc_matched_strength/per_patient_per_band.csv"
COPH_C = ROOT / "data/audit/rho_sym_gate/cohort_summary.csv"
COPH_P = ROOT / "data/audit/rho_sym_gate/per_patient_per_band.csv"
GRASS_C = ROOT / "data/audit/grassmann_cluster_extent/cohort_summary.csv"
OUT = ROOT / "data/outputs/figures/talk/fig_before_nulls_three_measure.pdf"

# slow -> fast (top of each panel = delta); reversed on the y-axis so delta sits on top
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
INK = "#2b2f36"


def _wilcox_pos(csv, band):
    df = pd.read_csv(csv)
    v = df[df.band == band]["obs_rho"].to_numpy()
    try:
        return wilcoxon(v, alternative="greater").pvalue
    except ValueError:
        return np.nan


def build():
    rawc = pd.read_csv(RAW_C).set_index("band")
    cophc = pd.read_csv(COPH_C).set_index("band")
    grassc = pd.read_csv(GRASS_C).set_index("band")

    panels = []
    # (title, subtitle, native-value dict, consistency-flag dict, value-format, has_sign)
    raw_val = {b: float(rawc.loc[b, "obs_median_rho"]) for b in BANDS}
    raw_ok = {b: _wilcox_pos(RAW_P, b) < 0.05 for b in BANDS}
    panels.append(("raw edges", "pairwise · fires everywhere", raw_val, raw_ok,
                   lambda x: f"{x:+.2f}", True))

    grass_val = {b: float(grassc.loc[b, "obs_cluster_mass_neglog10p"]) for b in BANDS}
    grass_ok = {b: float(grassc.loc[b, "obs_longest_run"])
                > float(grassc.loc[b, "null_p95_LR"]) for b in BANDS}
    panels.append(("Grassmann subspace", "global modes", grass_val, grass_ok,
                   lambda x: f"{x:.0f}", False))

    coph_val = {b: float(cophc.loc[b, "obs_median_rho_sym"]) for b in BANDS}
    coph_ok = {b: _wilcox_pos(COPH_P, b) < 0.05 for b in BANDS}
    panels.append(("cophenetic hierarchy", "multiscale · nested", coph_val, coph_ok,
                   lambda x: f"{x:+.2f}", True))
    return panels


def main():
    panels = build()
    y = np.arange(len(BANDS))[::-1]      # delta on top

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.3), sharey=True)

    for ax, (title, sub, val, ok, fmt, signed) in zip(axes, panels):
        vmax = max(abs(v) for v in val.values())
        for yi, b in zip(y, BANDS):
            raw = val[b]
            norm = raw / vmax
            col = band_color(b, 0.95)
            ax.barh(yi, norm, height=0.66, color=col,
                    edgecolor=INK if b == "beta" else "none",
                    linewidth=2.2 if b == "beta" else 0, zorder=3)
            # native value + consistency check
            xt = norm + (0.04 if norm >= 0 else -0.04)
            ha = "left" if norm >= 0 else "right"
            chk = r"  $\checkmark$" if ok[b] else ""
            ax.text(xt, yi, fmt(raw) + chk, va="center", ha=ha,
                    fontsize=10.5, color=INK if ok[b] else "#9aa0a7",
                    fontweight="bold" if b == "beta" else "normal", zorder=4)
        if signed:
            ax.axvline(0, color="0.6", lw=1.0, zorder=1)
        ax.set_xlim(-1.15 if signed else -0.02, 1.5)
        ax.set_title(title, fontsize=14, color=INK, pad=14)
        ax.text(0.5, 1.005, sub, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=9, color="#6a7079")
        ax.set_xticks([])
        for s in ("top", "right", "bottom"):
            ax.spines[s].set_visible(False)
        ax.spines["left"].set_visible(False)

    axes[0].set_yticks(y)
    axes[0].set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=15)
    for yi, b in zip(y, BANDS):
        axes[0].get_yticklabels()[list(y).index(yi)].set_color(band_color(b, 0.95))

    fig.text(0.5, -0.02,
             r"observed cohort effect, before the strength / drift nulls (bar normalized within read; value native)   ·   "
             r"$\checkmark$ = the read's own basic cohort test   ·   "
             r"$\beta$ is loudest under all three — the reads disagree on the rest",
             ha="center", va="center", fontsize=9.0, color="#6a7079")

    fig.subplots_adjust(left=0.09, right=0.98, top=0.86, bottom=0.10, wspace=0.28)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("=== before-nulls three-measure band profile ===")
    for title, sub, val, ok, fmt, signed in panels:
        print(f"\n{title} ({sub})")
        for b in BANDS:
            print(f"  {b:11s} {fmt(val[b]):>7s}  {'consistent' if ok[b] else 'not'}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
