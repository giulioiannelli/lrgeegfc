"""KC cohort fingerprint — sign-count 2D plane.

Replaces the magnitude-driven ellipse-cloud view with a scale-free
visualization of the actual cohort claim: per band, plot
(n_below_null at λ=0, n_below_null at λ=1) ∈ [0, 10]² as a single dot.

Trace zone = upper-right [≥7, ≥7] corner shaded green. Diagonal y=x
reference. Each band labeled. β at (10, 10); other bands cluster near
the chance corner (5, 5).

Data source: data/audit/section5_v2_kc_controls/per_patient_table.csv
Output:      data/reports/preprint/cohort_measures/kc_sign_count_plane.pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()


ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "data/audit/section5_v2_kc_controls/per_patient_table.csv"
OUT_DIR = ROOT / "data/reports/preprint/cohort_measures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_PDF = OUT_DIR / "kc_sign_count_plane.pdf"

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
HIGHLIGHT = {"beta", "alpha"}
BAND_COLORS = {
    "delta": "#888888", "theta": "#888888", "alpha": "#56a0d3",
    "beta": "#c0392b", "low_gamma": "#888888", "high_gamma": "#888888",
}


def render(df: pd.DataFrame, out_pdf: Path) -> None:
    counts = []
    for band in BAND_ORDER:
        n0 = int(df[(df["band"] == band) & (df["lam"] == 0.0)]["real_below_null"].sum())
        n1 = int(df[(df["band"] == band) & (df["lam"] == 1.0)]["real_below_null"].sum())
        counts.append((band, n0, n1))

    fig, ax = plt.subplots(figsize=(5.6, 5.6))

    # Trace zone (≥7, ≥7) — light green wedge (top-right)
    ax.fill([6.5, 10.5, 10.5, 6.5], [6.5, 6.5, 10.5, 10.5],
            color="#1a7c3e", alpha=0.16, zorder=0)
    # Anti-trace zone (≤3, ≤3) — same-weight red wedge (bottom-left)
    ax.fill([-0.5, 3.5, 3.5, -0.5], [-0.5, -0.5, 3.5, 3.5],
            color="#c0392b", alpha=0.16, zorder=0)
    # Chance lines (5/10 = binomial midpoint per axis)
    ax.axhline(5, color="0.7", lw=0.5, ls=":", zorder=1)
    ax.axvline(5, color="0.7", lw=0.5, ls=":", zorder=1)
    # Reference diagonal y = x (topology agreement = heights agreement);
    # NOT a null reference — the null/chance point is (5, 5).
    ax.plot([0, 10], [0, 10], color="0.85", lw=0.5, ls=":", zorder=1)
    # Explicit chance marker at (5, 5)
    ax.scatter([5], [5], marker="x", s=80, color="0.55",
               linewidth=1.4, zorder=2)
    ax.annotate("chance\n(5, 5)", (5, 5), xytext=(5.5, 4.5),
                fontsize=8, color="0.45", ha="left", va="top")

    # Per-band markers
    for band, n0, n1 in counts:
        is_hl = band in HIGHLIGHT
        color = BAND_COLORS[band]
        ax.scatter([n0], [n1], s=420 if is_hl else 220,
                   facecolor=color if is_hl else "white",
                   edgecolor=color, linewidth=2.0 if is_hl else 1.4,
                   zorder=4, alpha=0.95)
        # Label offset
        dx = 0.35 if n0 < 9 else -0.35
        dy = 0.30 if n1 < 9 else -0.30
        ha = "left" if dx > 0 else "right"
        va = "bottom" if dy > 0 else "top"
        ax.annotate(BRAIN_BAND_TEX_DICT[band], (n0, n1),
                    xytext=(n0 + dx, n1 + dy), fontsize=14 if is_hl else 11,
                    fontweight="bold" if is_hl else "normal",
                    ha=ha, va=va, color=color)

    # Corner annotations — symmetric on both ends of the diagonal
    ax.text(0.97, 0.97, "TRACE\n(10/10)",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=11, fontweight="bold", color="#1a7c3e")
    ax.text(0.03, 0.03, "ANTI-TRACE\n(0/10)",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=11, fontweight="bold", color="#c0392b")
    # Mirrored direction arrows from chance to both corners.
    ax.annotate("", xy=(9.5, 9.5), xytext=(5.5, 5.5),
                arrowprops=dict(arrowstyle="->", color="#1a7c3e",
                                lw=1.2, alpha=0.7))
    ax.annotate("", xy=(0.5, 0.5), xytext=(4.5, 4.5),
                arrowprops=dict(arrowstyle="->", color="#c0392b",
                                lw=1.2, alpha=0.7))
    ax.text(7.5, 7.0, "stronger trace →", rotation=45, color="#1a7c3e",
            fontsize=8, ha="center", va="bottom", fontstyle="italic")
    ax.text(2.5, 3.0, "← stronger anti-trace", rotation=45, color="#c0392b",
            fontsize=8, ha="center", va="top", fontstyle="italic")

    ax.set_xlim(-0.5, 10.5)
    ax.set_ylim(-0.5, 10.5)
    ax.set_xticks(range(0, 11))
    ax.set_yticks(range(0, 11))
    ax.set_xlabel(r"$n_{\mathrm{below\,null}}^{\lambda=0}$  —  topology")
    ax.set_ylabel(r"$n_{\mathrm{below\,null}}^{\lambda=1}$  —  heights")
    ax.set_aspect("equal")
    ax.grid(color="#eee", lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    handles = [
        Patch(facecolor="#1a7c3e", alpha=0.16, label="trace zone (≥7/10 below null on both)"),
        Patch(facecolor="#c0392b", alpha=0.16, label="anti-trace zone (≤3/10 below null on both)"),
    ]
    ax.legend(handles=handles, loc="lower right",
              bbox_to_anchor=(1.0, 0.06), frameon=False)

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


def main() -> None:
    render(pd.read_csv(SRC), OUT_PDF)


if __name__ == "__main__":
    main()
