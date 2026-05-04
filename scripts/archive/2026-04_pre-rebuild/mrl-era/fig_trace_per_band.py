#!/usr/bin/env python3
"""Step 3 (revised) — per-band cohort dotplot with permutation-null overlay.

For each band: a column of dots (one per patient) plotted at the
observed trace score, coloured by empirical p-value from the
size-matched permutation null in ``trace_null.csv``:

  - green       p < 0.01     (strong cohort-wide trace)
  - amber       p < 0.05     (significant)
  - grey        n.s.

The grey horizontal band shows the null distribution's mean ± p95
range pooled across patients within that band. Numbers above each
column report cohort tallies (#p<0.05 / N).

Reads:  ``data/reports/imcoh_mrl/trace_sweep.csv``,
        ``data/reports/imcoh_mrl/trace_null.csv``
Writes: ``data/reports/imcoh_mrl/figures/trace_per_band.pdf``
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT


IN_DIR = REPORTS_ROOT / "imcoh_mrl"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT

C_STRONG = "#2ca02c"     # p < 0.01
C_SIG    = "#ff9f1c"     # p < 0.05
C_NS     = "#b0b0b0"     # n.s.
C_NULL   = "#cccccc"     # null shading


def main() -> None:
    sweep = pd.read_csv(IN_DIR / "trace_sweep.csv")
    null  = pd.read_csv(IN_DIR / "trace_null.csv")

    # Merge null stats onto sweep rows.
    df = sweep.merge(null[["patient", "band", "null_mean", "null_p95",
                           "p_empirical", "z_observed"]],
                     on=["patient", "band"], how="left")

    fig, ax = plt.subplots(figsize=(10.0, 5.2), dpi=150)

    rng = np.random.default_rng(0)
    for j, band in enumerate(BANDS):
        sub = df[df["band"] == band]
        if sub.empty:
            continue

        # Null shading: mean ± p95 across patients in this band.
        null_lo = float(sub["null_mean"].mean()) - 0  # null mean ≈ 0
        null_hi = float(sub["null_p95"].mean())
        ax.fill_betweenx([null_lo, null_hi], j - 0.42, j + 0.42,
                         color=C_NULL, alpha=0.35, zorder=1)

        scores = sub["score"].to_numpy()
        sizes = sub["size"].to_numpy()
        ps = sub["p_empirical"].to_numpy()
        xs = j + rng.uniform(-0.15, 0.15, size=len(scores))

        cols = []
        for p in ps:
            if p < 0.01: cols.append(C_STRONG)
            elif p < 0.05: cols.append(C_SIG)
            else: cols.append(C_NS)

        ax.scatter(xs, scores,
                   s=18 + 1.6 * np.maximum(sizes, 1),
                   color=cols, alpha=0.9, edgecolor="black",
                   linewidth=0.3, zorder=4)

        # Cohort median (for reference)
        med = float(np.median(scores))
        ax.hlines(med, j - 0.22, j + 0.22, color="#404040", lw=2.2, zorder=5)

        # Tally on top
        n_05 = int((ps < 0.05).sum())
        n_01 = int((ps < 0.01).sum())
        ax.text(j, 1.04, f"{n_01}/{len(ps)} (p<.01)",
                ha="center", va="bottom", fontsize=9, color=C_STRONG,
                transform=ax.get_xaxis_transform(), fontweight="bold")
        ax.text(j, 1.10, f"{n_05}/{len(ps)} (p<.05)",
                ha="center", va="bottom", fontsize=9, color=C_SIG,
                transform=ax.get_xaxis_transform())

    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([TEX[b] for b in BANDS], fontsize=12)
    ax.set_xlim(-0.5, len(BANDS) - 0.5)
    ax.axhline(0, color="black", lw=0.6, ls="-", zorder=1)
    ax.set_ylim(-0.05, 1.08)
    ax.set_ylabel("trace score = J*(post→test) − J*(post→pre)", fontsize=11)
    ax.set_xlabel("band", fontsize=11)

    # Marker-size legend: convert scatter ``s`` (area in pts²) → markersize
    # (diameter in pts) for plt.Line2D, i.e. markersize = sqrt(s).
    def _ms_for_size(size: int) -> float:
        return float(np.sqrt(18 + 1.6 * size))

    handles = [
        plt.Line2D([0], [0], marker='o', color=C_STRONG, lw=0,
                   markersize=8, markeredgecolor="black",
                   markeredgewidth=0.3, label="p < 0.01"),
        plt.Line2D([0], [0], marker='o', color=C_SIG, lw=0,
                   markersize=8, markeredgecolor="black",
                   markeredgewidth=0.3, label="p < 0.05"),
        plt.Line2D([0], [0], marker='o', color=C_NS, lw=0,
                   markersize=8, markeredgecolor="black",
                   markeredgewidth=0.3, label="n.s."),
        plt.Rectangle((0, 0), 1, 1, facecolor=C_NULL, alpha=0.5,
                      label="null mean → p95 (size-matched random subsets)"),
        plt.Line2D([0], [0], color="#404040", lw=2.2, label="cohort median"),
        plt.Line2D([0], [0], lw=0, label=" "),  # spacer
        plt.Line2D([0], [0], lw=0, label="dot size ∝ #contacts in module:"),
        plt.Line2D([0], [0], marker='o', color="#404040", lw=0,
                   markersize=_ms_for_size(5), markeredgecolor="black",
                   markeredgewidth=0.3, label="5 contacts"),
        plt.Line2D([0], [0], marker='o', color="#404040", lw=0,
                   markersize=_ms_for_size(15), markeredgecolor="black",
                   markeredgewidth=0.3, label="15 contacts"),
        plt.Line2D([0], [0], marker='o', color="#404040", lw=0,
                   markersize=_ms_for_size(30), markeredgecolor="black",
                   markeredgewidth=0.3, label="30 contacts"),
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
              frameon=False, fontsize=9)

    ax.set_title("Trace score per patient × band, against size-matched null  "
                 "(dot size ∝ module size; tally above = significant patients)",
                 fontsize=10.5)

    fig.tight_layout()
    out = OUT_DIR / "trace_per_band.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
