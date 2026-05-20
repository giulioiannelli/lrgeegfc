"""KC cohort fingerprint — signed-delta caterpillar (sign-prominent, scale-tamed).

Per (patient, band, λ), compute the signed log-compressed delta:

  Δ̃ = sign(real - null) * log10(1 + |real - null|)

The log compression suppresses single-patient outliers (some real_T_KC
reach ±20 while the cohort sign is unanimous); the sign prominence
keeps the cohort's directional message intact.

Two side-by-side panels (λ=0 topology, λ=1 heights). Each panel: 6
horizontal strips (one per band), 10 dots per strip (one per patient),
sorted by Δ̃ within the strip. Red = trace direction (Δ̃ < 0), grey
otherwise. β strip: 10 red dots all on the left side of zero. γ_l / θ:
dots straddle zero.

Data source: data/audit/section5_v2_kc_controls/per_patient_table.csv
Output:      data/reports/preprint/cohort_measures/kc_caterpillar.pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()


ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "data/audit/section5_v2_kc_controls/per_patient_table.csv"
OUT_DIR = ROOT / "data/reports/preprint/cohort_measures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_PDF = OUT_DIR / "kc_caterpillar.pdf"

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
HIGHLIGHT = {"beta", "alpha"}


def signed_log(x: np.ndarray) -> np.ndarray:
    return np.sign(x) * np.log10(1.0 + np.abs(x))


def render(df: pd.DataFrame, out_pdf: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), sharey=True)
    panel_titles = [
        r"$\lambda = 0$  —  topology",
        r"$\lambda = 1$  —  heights",
    ]
    lams = [0.0, 1.0]

    all_vals = []
    for lam in lams:
        sub = df[df["lam"] == lam]
        delta = sub["real_T_KC"].values - sub["null_T_KC"].values
        all_vals.append(signed_log(delta))
    extent = float(np.max(np.abs(np.concatenate(all_vals)))) * 1.05

    band_y = {b: i for i, b in enumerate(BAND_ORDER)}

    for ax, lam, title in zip(axes, lams, panel_titles):
        sub = df[df["lam"] == lam]
        ax.axvline(0, color="0.3", lw=0.7, ls="--", zorder=2)
        # Trace zone (left of zero) — light red
        ax.axvspan(-extent, 0, color="#c0392b", alpha=0.05, zorder=0)

        for band in BAND_ORDER:
            sb = sub[sub["band"] == band].copy()
            if not len(sb):
                continue
            d = sb["real_T_KC"].values - sb["null_T_KC"].values
            d_log = signed_log(d)
            order = np.argsort(d_log)
            d_log = d_log[order]
            y = band_y[band]
            colors = ["#a82a1a" if v < 0 else "#888" for v in d_log]
            jitter = np.linspace(-0.18, 0.18, len(d_log))
            is_hl = band in HIGHLIGHT
            # Optional connecting line to indicate strip
            ax.plot([d_log.min(), d_log.max()], [y, y],
                    color="0.85", lw=0.6, zorder=1)
            ax.scatter(d_log, np.full_like(d_log, y, dtype=float) + jitter,
                       c=colors, s=60 if is_hl else 38,
                       edgecolor="white" if is_hl else "0.6",
                       linewidth=0.6 if is_hl else 0.3,
                       zorder=4, alpha=0.95)
            # Cohort median tick
            med = float(np.median(d_log))
            ax.plot([med, med], [y - 0.30, y + 0.30],
                    color="black", lw=1.4, zorder=5)
            # Annotation: n_below_null
            n_below = int((d < 0).sum())
            ax.text(extent * 0.98, y, f"  {n_below}/10",
                    ha="left", va="center",
                    fontsize=9 if not is_hl else 11,
                    fontweight="bold" if is_hl else "normal",
                    color="#a82a1a" if n_below >= 7 else "#444")

        ax.set_xlim(-extent, extent * 1.18)
        ax.set_ylim(len(BAND_ORDER) - 0.5, -0.5)
        ax.set_yticks(range(len(BAND_ORDER)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER],
                           fontsize=12)
        ax.set_xlabel(r"$\mathrm{sgn}(r-n)\,\cdot\,\log_{10}(1+|r-n|)$  "
                      r"  $\quad r = T_{KC}^{\mathrm{real}}, n = T_{KC}^{\mathrm{null}}$",
                      fontsize=8)
        ax.set_title(title)
        ax.tick_params(axis="y", which="both", length=0)
        ax.spines[["top", "right"]].set_visible(False)

    # Annotation legend
    ax_l = axes[0]
    ax_l.text(-extent * 0.95, -0.45, "← trace", color="#a82a1a",
              fontweight="bold", ha="left", va="bottom", fontsize=9)
    ax_l.text(extent * 0.95, -0.45, "anti-trace →", color="#888",
              ha="right", va="bottom", fontsize=9)

    handles = [
        Line2D([], [], marker="o", linestyle="None", markersize=7,
               markerfacecolor="#a82a1a", markeredgecolor="white",
               label=r"patient: real < null (trace)"),
        Line2D([], [], marker="o", linestyle="None", markersize=7,
               markerfacecolor="#888", markeredgecolor="0.6",
               label=r"patient: real ≥ null"),
        Line2D([], [], marker="|", linestyle="None", markersize=14,
               markeredgecolor="black", markeredgewidth=1.4,
               label="band median"),
    ]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.06), ncol=3, frameon=False)

    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


def main() -> None:
    render(pd.read_csv(SRC), OUT_PDF)


if __name__ == "__main__":
    main()
