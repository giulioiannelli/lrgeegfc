#!/usr/bin/env python3
"""Per-band Δρ(k) trajectories across patients — multiscale persistence.

The scientific question: at each dendrogram cut k, does task-induced
block persistence survive across patients in a consistent way? Bands
where trajectories cluster above zero across a coherent k range carry
a multiscale memory signature; bands where trajectories diverge or
straddle zero do not.

Layout: 6 band subplots (2×3 grid). In each:
  * 9 thin patient curves of Δρ(k) = ρ_task(k) − ρ_inert(k), k ∈ [2, 49]
  * Median curve + IQR shaded band
  * Zero line
  * Green shading in k-ranges where all 9 patients have Δρ > 0
    (strict unanimity), orange where ≥8/9, grey where <8/9

Reads:  data/reports/imcoh_vi/h2d_persistence_raw.csv
Writes: data/reports/imcoh_vi/figures/band_persistence_trajectories.{pdf,png}
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT


IN_DIR = REPORTS_ROOT / "imcoh_vi"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT

LINE_COLOR   = "#5b7db3"
MEDIAN_COLOR = "#1a1a1a"
BAND_COLOR   = "#c8d6ea"
UNANIM_FULL  = "#b4e5b4"   # all 9 patients positive
UNANIM_NEAR  = "#ffd9a8"   # 8/9 positive
ZERO_COLOR   = "#888888"
THETA_COLOR  = "#c44e4e"


def main() -> None:
    raw = pd.read_csv(IN_DIR / "h2d_persistence_raw.csv")
    # Pivot to (patient, band) x k matrices
    ks = sorted(raw["k"].unique())
    patients = sorted(raw["patient"].unique())

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.0), dpi=160, sharex=True, sharey=True)
    axes = axes.ravel()

    for ib, band in enumerate(BANDS):
        ax = axes[ib]
        bdf = raw[raw["band"] == band]
        if bdf.empty:
            ax.set_title(f"{TEX[band]} (no data)", fontsize=11)
            continue
        mat = np.full((len(patients), len(ks)), np.nan)
        for ip, pat in enumerate(patients):
            sub = bdf[bdf["patient"] == pat].set_index("k")
            for jk, k in enumerate(ks):
                if k in sub.index:
                    mat[ip, jk] = sub.at[k, "delta_rho"]

        # Background shading by unanimity level at each k
        for jk, k in enumerate(ks):
            col = mat[:, jk]
            finite = col[np.isfinite(col)]
            if finite.size == 0:
                continue
            pos = int((finite > 0).sum())
            tot = finite.size
            if pos == tot:
                shade = UNANIM_FULL
            elif pos >= tot - 1:
                shade = UNANIM_NEAR
            else:
                continue
            ax.axvspan(k - 0.5, k + 0.5, color=shade, alpha=0.55, zorder=0,
                        edgecolor="none")

        # Individual patient trajectories
        for ip, pat in enumerate(patients):
            y = mat[ip]
            ok = np.isfinite(y)
            if ok.sum() >= 3:
                ax.plot(np.array(ks)[ok], y[ok], color=LINE_COLOR,
                        alpha=0.4, lw=0.9, zorder=2)

        # Median ± IQR
        med = np.nanmedian(mat, axis=0)
        q1 = np.nanpercentile(mat, 25, axis=0)
        q3 = np.nanpercentile(mat, 75, axis=0)
        ax.fill_between(ks, q1, q3, color=BAND_COLOR, alpha=0.7, zorder=3)
        med_color = THETA_COLOR if band == "theta" else MEDIAN_COLOR
        ax.plot(ks, med, color=med_color, lw=2.0, zorder=4,
                label="median" if ib == 0 else None)

        ax.axhline(0, color=ZERO_COLOR, lw=0.8, linestyle="--", zorder=1)

        # Annotate longest strict-unanimity run
        unan = np.array([
            np.isfinite(mat[:, jk]).any()
            and (mat[np.isfinite(mat[:, jk]), jk] > 0).all()
            for jk in range(len(ks))
        ])
        longest = _longest_run(unan)
        if longest[2] > 0:
            kmin, kmax, ln = ks[longest[0]], ks[longest[1]], longest[2]
            ax.text(0.98, 0.04,
                    f"longest unan run: k={kmin}–{kmax}  (L={ln})",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=8.5, color="#333",
                    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))

        ax.set_title(TEX[band], fontsize=13)
        if ib % 3 == 0:
            ax.set_ylabel(r"$\Delta\rho(k) = \rho_{task} - \rho_{inert}$",
                           fontsize=10)
        if ib >= 3:
            ax.set_xlabel("k (cut scale)", fontsize=10)
        ax.set_xlim(min(ks), max(ks))

    # Global legend
    from matplotlib.patches import Patch
    handles = [
        plt.Line2D([0], [0], color=LINE_COLOR, lw=0.9, alpha=0.8,
                   label="per-patient Δρ(k)"),
        plt.Line2D([0], [0], color=MEDIAN_COLOR, lw=2.0, label="median"),
        Patch(facecolor=BAND_COLOR, alpha=0.7, label="IQR"),
        Patch(facecolor=UNANIM_FULL, alpha=0.55, label="k with 9/9 patients Δρ > 0"),
        Patch(facecolor=UNANIM_NEAR, alpha=0.55, label="k with 8/9 patients Δρ > 0"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               bbox_to_anchor=(0.5, -0.02), fontsize=9.5)
    # No fig.suptitle on publication figures — context lives in the file name.
    fig.tight_layout()

    for ext in ("pdf", "png"):
        out = OUT_DIR / f"band_persistence_trajectories.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)


def _longest_run(mask: np.ndarray) -> tuple[int, int, int]:
    best = (0, -1, 0)
    cur_start = None
    for i, v in enumerate(mask):
        if v and cur_start is None:
            cur_start = i
        elif not v and cur_start is not None:
            length = i - cur_start
            if length > best[2]:
                best = (cur_start, i - 1, length)
            cur_start = None
    if cur_start is not None:
        length = len(mask) - cur_start
        if length > best[2]:
            best = (cur_start, len(mask) - 1, length)
    return best


if __name__ == "__main__":
    main()
