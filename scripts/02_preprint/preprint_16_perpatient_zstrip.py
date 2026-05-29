#!/usr/bin/env python3
"""Per-patient β-band ρ_split^coph signal vs matched-strength surrogate
(standalone, multi-style).

Standalone version of panel (a) from
``preprint_07_beta_rho_split_figure_test2``.  Two render styles selectable
via ``--style``:

* ``stem`` (default, legacy):
    Each patient = one horizontal stem from x = 0 to the observed z-score
    above the patient's own matched-strength surrogate.  Grey vertical
    band marks the universal noise-floor zone ``|z| ≤ 2``; stems colored
    grass-green for ``z > 2`` (trace), brick-red for ``z < -2`` (anti),
    grey otherwise.  Sorted by z descending.  Reads at a glance.

* ``whisker``:
    Each patient = one horizontal box-and-whisker showing the surrogate
    Spearman-ρ distribution (box = p25-p75, whiskers = p5-p95, vertical
    tick = p50) with the observed ρ as a colored dot.  Coloured by
    whether the observed ρ beats the patient's own surrogate p95.
    Shows the actual ρ scales (ρ ∈ [-1, 1]) instead of normalised z;
    useful when reviewers ask "how far above noise, in real units?"

Usage
-----
    python preprint_16_perpatient_zstrip.py                       # beta, stem
    python preprint_16_perpatient_zstrip.py --style whisker
    python preprint_16_perpatient_zstrip.py --band alpha --style whisker

Inputs
------
data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv

Output (PDF only)
-----------------
data/preprint/figures/<band>/per_pair_trace/
    fig_<band>_perpatient_<style>.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()


COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]

LRG_CTM_DIR = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"

CLR_TRACE = "#5fa844"   # grass green
CLR_ANTI = "#c0392b"    # brick red
CLR_NEUTRAL = "0.55"

STYLES = ("stem", "whisker")


def _sort_descending_by_z(per_pat: pd.DataFrame, band: str) -> pd.DataFrame:
    sub = per_pat[per_pat.band == band].set_index("patient").loc[COHORT]
    order = np.argsort(sub["obs_z"].values)[::-1]
    return sub.iloc[order]


def render_stem(ax: plt.Axes, per_pat: pd.DataFrame, band: str,
                band_tex: str) -> None:
    """Horizontal-stem style — panel (a) of fig_beta_rho_split_test2.

    Each patient is one horizontal stem from x = 0 to the observed z;
    grey vertical band at ``|z| ≤ 2`` marks the universal noise-floor
    zone; coloured by sign of excess past the floor.
    """
    sub = _sort_descending_by_z(per_pat, band)
    z = sub["obs_z"].values
    patients_sorted = sub.index.tolist()

    if "obs_p_one_sided" in sub.columns:
        frac = 1.0 - sub["obs_p_one_sided"].values
    elif "obs_p_emp" in sub.columns:
        frac = 1.0 - sub["obs_p_emp"].values
    elif "frac_above_surr" in sub.columns:
        frac = sub["frac_above_surr"].values
    else:
        frac = np.full(len(sub), np.nan)

    n = len(z)
    y = np.arange(n)[::-1]

    ax.axvspan(-2, 2, color="0.92", alpha=0.55, zorder=0)
    ax.axvline(0, color="0.35", lw=0.8, zorder=1)

    for yi, zi in zip(y, z):
        if zi > 2:
            clr = CLR_TRACE
        elif zi < -2:
            clr = CLR_ANTI
        else:
            clr = CLR_NEUTRAL
        ax.plot([0, zi], [yi, yi], color=clr, lw=1.9, alpha=0.85,
                zorder=2, solid_capstyle="round")
        ax.scatter([zi], [yi], s=110, color=clr, edgecolor="white",
                   linewidth=1.2, zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels(patients_sorted, fontsize=9.5)
    ax.set_ylim(-0.6, n - 0.4)

    z_max = max(float(np.nanmax(z)) + 1.2, 3.5)
    z_min = min(float(np.nanmin(z)) - 1.0, -3.0)
    ax.set_xlim(z_min, z_max)
    ax.set_xlabel(r"$z$ above own matched-strength surrogate", fontsize=10)

    if np.isfinite(frac).any():
        x0 = z_max - 0.35 * (z_max - z_min) * 0.04
        bar_max_w = (z_max - z_min) * 0.08
        for yi, fr in zip(y, frac):
            if not np.isfinite(fr):
                continue
            fr_clip = float(np.clip(fr, 0.0, 1.0))
            bar_clr = CLR_TRACE if fr_clip >= 0.95 else CLR_NEUTRAL
            ax.add_patch(plt.Rectangle(
                (x0 - bar_max_w, yi - 0.18),
                bar_max_w * fr_clip, 0.36,
                facecolor=bar_clr, edgecolor="none",
                alpha=0.75, zorder=2,
            ))
            ax.add_patch(plt.Rectangle(
                (x0 - bar_max_w, yi - 0.18),
                bar_max_w, 0.36,
                facecolor="none", edgecolor="0.65",
                lw=0.5, zorder=2.5,
            ))

    ax.spines[["top", "right"]].set_visible(False)


def render_whisker(ax: plt.Axes, per_pat: pd.DataFrame, band: str,
                   band_tex: str) -> None:
    """Box-and-whisker style — surrogate ρ distribution per patient with
    observed ρ as a coloured dot.

    Box = surrogate p25-p75, whiskers = surrogate p5-p95, vertical tick
    inside the box = surrogate p50.  Observed ρ rendered as a dot on
    top.  Coloured green if observed beats own p95, red if below own p5,
    grey otherwise.  Sorted by ``obs_z`` descending — same order as the
    stem style so the two figures are visually aligned.
    """
    sub = _sort_descending_by_z(per_pat, band)
    patients_sorted = sub.index.tolist()
    obs_rho = sub["obs_rho"].values
    surr_p5 = sub["surr_p5"].values
    surr_p25 = sub["surr_p25"].values
    surr_p50 = sub["surr_p50"].values
    surr_p75 = sub["surr_p75"].values
    surr_p95 = sub["surr_p95"].values
    z = sub["obs_z"].values

    n = len(obs_rho)
    y = np.arange(n)[::-1]
    box_h = 0.34
    whisker_h = 0.18

    x_min = float(min(np.nanmin(surr_p5), np.nanmin(obs_rho)) - 0.06)
    x_max = float(max(np.nanmax(surr_p95), np.nanmax(obs_rho)) + 0.04)

    ax.axvline(0, color="0.35", lw=0.8, zorder=1)

    z_abs_max = float(max(np.nanmax(np.abs(z)), 1.0))

    for yi, p5, p25, p50, p75, p95, ro, zi in zip(
        y, surr_p5, surr_p25, surr_p50, surr_p75, surr_p95, obs_rho, z,
    ):
        # Whiskers: line from p5 to p25 and p75 to p95
        ax.plot([p5, p25], [yi, yi], color="0.40", lw=1.1, zorder=2,
                solid_capstyle="round")
        ax.plot([p75, p95], [yi, yi], color="0.40", lw=1.1, zorder=2,
                solid_capstyle="round")
        for cap_x in (p5, p95):
            ax.plot([cap_x, cap_x], [yi - whisker_h, yi + whisker_h],
                    color="0.40", lw=0.9, zorder=2)
        # Box: p25 to p75
        ax.add_patch(plt.Rectangle(
            (p25, yi - box_h / 2), p75 - p25, box_h,
            facecolor="0.82", edgecolor="0.30", lw=0.9, zorder=2.5,
        ))
        # Median tick inside the box
        ax.plot([p50, p50], [yi - box_h / 2, yi + box_h / 2],
                color="0.10", lw=1.4, zorder=3)
        # Sign-driven palette + null-boundary anchor for the excess line
        if ro > p95:
            clr = CLR_TRACE
            line_x0 = p95
            has_line = True
        elif ro < p5:
            clr = CLR_ANTI
            line_x0 = p5
            has_line = True
        else:
            clr = CLR_NEUTRAL
            line_x0 = None
            has_line = False
        # Excess line — geometric encoding of z = (ρ - μ_surr) / σ_surr:
        # line length already shows the absolute excess in ρ units; line
        # thickness reinforces |z| so a z=8 patient gets a fat coloured
        # arm and a z=2 patient gets a thin one.  Soft halo behind it for
        # |z| > 2 so the trace patients pop further from the page.
        if has_line:
            line_lw = float(np.clip(1.4 + 0.55 * abs(zi), 1.4, 5.5))
            if abs(zi) > 2:
                ax.plot([line_x0, ro], [yi, yi], color=clr,
                        lw=line_lw + 4.0, alpha=0.18,
                        zorder=3.6, solid_capstyle="round")
            ax.plot([line_x0, ro], [yi, yi], color=clr,
                    lw=line_lw, alpha=0.92,
                    zorder=4, solid_capstyle="round")
        # Observed-ρ dot — radius scales with |z| (capped at z_abs_max)
        # so the strongest patients project further from the page.
        z_frac = min(abs(zi) / z_abs_max, 1.0)
        dot_s = float(60 + 220 * z_frac)
        ax.scatter([ro], [yi], s=dot_s, color=clr, edgecolor="white",
                   linewidth=1.4, zorder=5)

    ax.set_yticks(y)
    ax.set_yticklabels(patients_sorted, fontsize=9.5)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlim(x_min, x_max)
    ax.set_xlabel(r"$\rho_{\mathrm{split}}^{\mathrm{coph}}$  "
                  r"(box = surrogate p25-p75, whiskers = p5-p95, "
                  r"dot = observed, line + size $\propto |z|$)",
                  fontsize=9.5)

    ax.spines[["top", "right"]].set_visible(False)


def main(band: str = "beta", style: str = "stem") -> Path:
    if style not in STYLES:
        raise ValueError(f"style must be in {STYLES}, got {style!r}")
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")

    per_pat = pd.read_csv(LRG_CTM_DIR / "per_patient_per_band.csv")

    out_dir = ROOT / "data" / "preprint" / "figures" / band / "per_pair_trace"
    out_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(7.5, 5.0))
    gs = fig.add_gridspec(
        1, 1, left=0.13, right=0.98, top=0.96, bottom=0.12,
    )
    ax = fig.add_subplot(gs[0, 0])

    if style == "stem":
        render_stem(ax, per_pat, band, band_tex)
    elif style == "whisker":
        render_whisker(ax, per_pat, band, band_tex)

    out = out_dir / f"fig_{band}_perpatient_{style}.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--band", default="beta",
                        choices=["delta", "theta", "alpha", "beta",
                                 "low_gamma", "high_gamma"])
    parser.add_argument("--style", default="stem",
                        choices=list(STYLES),
                        help="Render style: 'stem' (legacy z-stems) or "
                             "'whisker' (boxplot of surrogate ρ + "
                             "observed ρ dot).")
    args = parser.parse_args()
    main(args.band, args.style)
