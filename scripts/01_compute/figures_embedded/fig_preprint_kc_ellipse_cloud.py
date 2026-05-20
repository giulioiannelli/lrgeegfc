"""KC cohort fingerprint — band ellipse-cloud in (-T_KC^0, -T_KC^1) plane.

Augments the existing kc_decomposition idiom with cohort variance:
each band gets a 1σ ellipse around its 10 per-patient (real_T_KC^0,
real_T_KC^1) pairs, sign-flipped so trace direction = upper-right
quadrant. Band-level median dot kept (size ∝ −log10 p_min from the
cohort summary).

β should sit in the trace quadrant with most of its ellipse on the
trace side; θ / γ_l ellipses overlap the origin.

Data sources:
  data/audit/section5_v2_kc_controls/per_patient_table.csv
  data/audit/section5_v2_round3_redo/tables/kc_decomposition_points.csv
Output:
  data/reports/preprint/cohort_measures/kc_ellipse_cloud.pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()


ROOT = Path(__file__).resolve().parents[3]
PER_PAT = ROOT / "data/audit/section5_v2_kc_controls/per_patient_table.csv"
BAND_LVL = ROOT / "data/audit/section5_v2_round3_redo/tables/kc_decomposition_points.csv"
OUT_DIR = ROOT / "data/reports/preprint/cohort_measures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_PDF = OUT_DIR / "kc_ellipse_cloud.pdf"

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_COLORS = {
    "delta": "#7f7f7f",
    "theta": "#9e9e9e",
    "alpha": "#56a0d3",
    "beta": "#c0392b",
    "low_gamma": "#bababa",
    "high_gamma": "#cccccc",
}
HIGHLIGHT = {"beta", "alpha"}


def _ellipse_from_xy(x: np.ndarray, y: np.ndarray, color: str, alpha: float):
    if x.size < 2:
        return None
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals = vals[order]; vecs = vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width, height = 2.0 * np.sqrt(np.maximum(vals, 0.0))
    if width <= 0 or height <= 0:
        return None
    return Ellipse(
        (x.mean(), y.mean()), width, height, angle=angle,
        facecolor=color, alpha=alpha, edgecolor=color, lw=1.0, zorder=2,
    )


def render(per_pat: pd.DataFrame, band_lvl: pd.DataFrame, out_pdf: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 5.6))

    # Compute symmetric extent
    sub0 = per_pat[per_pat["lam"] == 0.0]
    sub1 = per_pat[per_pat["lam"] == 1.0]
    pivot0 = sub0.set_index(["band", "patient"])["real_T_KC"]
    pivot1 = sub1.set_index(["band", "patient"])["real_T_KC"]

    all_x = -pivot0.values  # sign-flip
    all_y = -pivot1.values
    extent = max(np.nanmax(np.abs(all_x)), np.nanmax(np.abs(all_y))) * 1.15

    # Trace-quadrant background shading
    ax.axhspan(0, extent, xmin=0.5, xmax=1.0, color="#1a7c3e", alpha=0.10, zorder=0)
    ax.axhline(0, color="0.4", lw=0.7, ls="--", zorder=1)
    ax.axvline(0, color="0.4", lw=0.7, ls="--", zorder=1)

    # Per-band ellipse + dots + median marker
    for band in BAND_ORDER:
        if band not in pivot0.index.get_level_values(0):
            continue
        # Per-patient pairs
        x_pat = -pivot0.loc[band].values
        y_pat = -pivot1.loc[band].values
        finite = np.isfinite(x_pat) & np.isfinite(y_pat)
        x_pat = x_pat[finite]; y_pat = y_pat[finite]
        if x_pat.size == 0:
            continue

        is_hl = band in HIGHLIGHT
        color = BAND_COLORS[band]
        ell_alpha = 0.20 if is_hl else 0.12
        dot_alpha = 0.55 if is_hl else 0.30

        ell = _ellipse_from_xy(x_pat, y_pat, color, ell_alpha)
        if ell is not None:
            ax.add_patch(ell)
        ax.scatter(x_pat, y_pat, s=18, color=color, alpha=dot_alpha,
                   edgecolor="white", linewidth=0.4, zorder=3)

        # Band-level median marker (size ∝ -log10 p_min)
        row = band_lvl[band_lvl["band"] == band]
        if len(row):
            r = row.iloc[0]
            T0 = float(r["T_topology"])  # already sign-flipped in source
            T1 = float(r["T_heights"])
            p_min = float(r["p_min"])
            size = float(np.clip(80 + 240.0 * (-np.log10(max(p_min, 1e-6)) / 3.0),
                                 80, 460))
            edge = "#7d1a0e" if (p_min <= 0.05) else "#444"
            face = color if is_hl else "white"
            ax.scatter([T0], [T1], s=size, facecolor=face,
                       edgecolor=edge, linewidth=1.6, zorder=5)
            # Annotation
            dx = 0.030 * extent
            dy = 0.030 * extent
            label = BRAIN_BAND_TEX_DICT[band]
            weight = "bold" if is_hl else "normal"
            fs = 13 if is_hl else 10
            ax.annotate(label, (T0, T1), xytext=(T0 + dx, T1 + dy),
                        fontsize=fs, fontweight=weight, color=color)

    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_xlabel(r"$-T_{KC}(\lambda=0)$  —  topology, trace direction →")
    ax.set_ylabel(r"$-T_{KC}(\lambda=1)$  —  heights,  trace direction ↑")
    ax.text(0.96, 0.96, "TRACE", ha="right", va="top",
            fontsize=11, fontweight="bold", color="#1a7c3e",
            transform=ax.transAxes)
    ax.text(0.04, 0.04, "anti-trace", ha="left", va="bottom",
            fontsize=9, fontstyle="italic", color="#c0392b",
            transform=ax.transAxes)
    ax.set_aspect("equal")
    ax.spines[["top", "right"]].set_visible(False)

    handles = [
        Patch(facecolor="#1a7c3e", alpha=0.15, edgecolor="none",
              label="trace quadrant"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor="#7d1a0e",
               markersize=9, markeredgewidth=1.6,
               label=r"band median ($p \leq 0.05$)"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor="#444",
               markersize=9, markeredgewidth=1.6,
               label=r"band median ($p > 0.05$)"),
        Patch(facecolor="#c0392b", alpha=0.20, edgecolor="#c0392b",
              label=r"cohort 1$\sigma$ ellipse (per band)"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="#c0392b", markeredgecolor="white",
               markersize=6, markeredgewidth=0.4, alpha=0.55,
               label="per-patient T_KC"),
    ]
    ax.legend(handles=handles, loc="lower right",
              bbox_to_anchor=(1.0, 0.06), frameon=False, ncol=1,
              handletextpad=0.6, labelspacing=0.5)

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


def main() -> None:
    per_pat = pd.read_csv(PER_PAT)
    band_lvl = pd.read_csv(BAND_LVL)
    render(per_pat, band_lvl, OUT_PDF)


if __name__ == "__main__":
    main()
