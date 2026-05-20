"""Grassmann cohort fingerprint — n_trace(k)/10 overlay across all bands.

The "cohort answer" view that complements the existing
grassmann_principal_angles 2×3 small-multiples figure: for each band,
plot the per-k cohort-agreement count (number of patients with
T_E1 < 0 = trace direction) as a line. β + α highlighted, others faint.

Reader sees:
  - which bands carry sustained cohort agreement above k > 12
  - β / α as the load-bearing bands (sustained 6–8/10)
  - θ / γ_l hovering near the n=5 chance line

Data sources:
  data/audit/section5_v2_round3_redo/tables/grassmann_principal_angles_per_patient.csv
Output:
  data/reports/preprint/cohort_measures/grassmann_cohort_answer.pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()


ROOT = Path(__file__).resolve().parents[3]
PER_PAT = ROOT / "data/audit/section5_v2_round3_redo/tables/grassmann_principal_angles_per_patient.csv"
OUT_DIR = ROOT / "data/reports/preprint/cohort_measures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_PDF = OUT_DIR / "grassmann_cohort_answer.pdf"

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_COLORS = {
    "delta":     "#cccccc",
    "theta":     "#999999",
    "alpha":     "#56a0d3",
    "beta":      "#c0392b",
    "low_gamma": "#bababa",
    "high_gamma": "#888888",
}
HIGHLIGHT = {"beta", "alpha"}


def cohort_n_trace(df: pd.DataFrame) -> pd.DataFrame:
    """Count patients with T_E1 < 0 per (band, k) — de-dupe across modes."""
    sub = df.drop_duplicates(subset=["patient", "band", "k"])
    agg = (
        sub.assign(is_trace=lambda d: (d["T_E1"] < 0).astype(int))
        .groupby(["band", "k"], as_index=False)
        .agg(n_trace=("is_trace", "sum"), n_pat=("patient", "nunique"))
    )
    return agg


def render(agg: pd.DataFrame, out_pdf: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 4.4))

    k_min = int(agg["k"].min())
    k_max = int(agg["k"].max())

    # Reference shading: cold zone (k < 12), majority threshold n=7
    ax.axvspan(k_min - 0.5, 12.5, color="#eeeeee", alpha=0.5, zorder=0)
    ax.axhline(7, color="#1a7c3e", lw=0.7, ls="--", zorder=1)
    ax.axhline(5, color="0.6", lw=0.6, ls=":", zorder=1)

    # Plot non-highlight bands first (so highlights overlay)
    for band in BAND_ORDER:
        sub = agg[agg["band"] == band].sort_values("k")
        if not len(sub) or band in HIGHLIGHT:
            continue
        ax.plot(sub["k"], sub["n_trace"], color=BAND_COLORS[band],
                lw=1.0, alpha=0.7, zorder=2,
                label=BRAIN_BAND_TEX_DICT[band])

    for band in [b for b in BAND_ORDER if b in HIGHLIGHT]:
        sub = agg[agg["band"] == band].sort_values("k")
        if not len(sub):
            continue
        ax.plot(sub["k"], sub["n_trace"], color=BAND_COLORS[band],
                lw=2.2, alpha=0.95, zorder=3,
                label=BRAIN_BAND_TEX_DICT[band])

    # Annotations
    ax.text(12.5, 9.6, r"$k > 12$ trace zone",
            ha="left", va="top", fontsize=9, color="#444",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="0.85", alpha=0.9))
    ax.text(k_max - 1, 7.18, r"$n = 7/10$ majority",
            ha="right", va="bottom", fontsize=8, color="#1a7c3e")
    ax.text(k_max - 1, 5.18, r"$n = 5/10$ chance",
            ha="right", va="bottom", fontsize=8, color="0.5")

    ax.set_xlim(k_min - 0.5, k_max + 0.5)
    ax.set_ylim(-0.4, 10.4)
    ax.set_yticks([0, 2, 5, 7, 10])
    ax.set_xlabel(r"spectral cutoff $k$")
    ax.set_ylabel(r"$n_{\mathrm{trace}}\,/\,10$")
    ax.spines[["top", "right"]].set_visible(False)

    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.0),
              ncol=3, frameon=False, columnspacing=1.2,
              handlelength=1.4, handletextpad=0.5)

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


def main() -> None:
    df = pd.read_csv(PER_PAT)
    agg = cohort_n_trace(df)
    render(agg, OUT_PDF)


if __name__ == "__main__":
    main()
