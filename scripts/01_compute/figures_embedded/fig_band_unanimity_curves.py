#!/usr/bin/env python3
"""Per-band unanimity curves — one line per band, k on x-axis.

The "all-in-one" band-comparison figure the current multi-panel
landscape buries. For each (band, k) we plot the fraction of patients
with Δρ(k) > 0; bands that carry the trace rise quickly above 0.7
across a coherent k range; ergodic bands stay near 0.5 and wander.

Two panels:
  (a) Fraction positive (H2d block-persistence)
  (b) Fraction with ρ_task > 2·ρ_inert (stronger "memory" criterion)

Reads:  data/reports/imcoh_vi/h2d_persistence_raw.csv
Writes: data/reports/imcoh_vi/figures/band_unanimity_curves.{pdf,png}
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

# Band colors — θ distinct (red), others in a perceptual spread
BAND_COLORS = {
    "delta":      "#3b6e9c",
    "theta":      "#c44e4e",   # ergodic candidate — red
    "alpha":      "#5ea85e",
    "beta":       "#b28ad1",
    "low_gamma":  "#e5a24b",
    "high_gamma": "#7d5a50",
}


def _frac_positive_curve(mat: np.ndarray) -> np.ndarray:
    """Per-column fraction of patients with value > 0, ignoring NaN."""
    finite = np.isfinite(mat)
    pos = finite & (mat > 0)
    n = finite.sum(axis=0)
    p = pos.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, p / n, np.nan)


def _frac_strong_curve(rho_task, rho_inert, thresh=2.0) -> np.ndarray:
    """Per-k fraction of patients where ρ_task > thresh · ρ_inert."""
    finite = np.isfinite(rho_task) & np.isfinite(rho_inert)
    strong = finite & (rho_task > thresh * rho_inert)
    n = finite.sum(axis=0)
    s = strong.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, s / n, np.nan)


def main() -> None:
    raw = pd.read_csv(IN_DIR / "h2d_persistence_raw.csv")
    ks = sorted(raw["k"].unique())
    patients = sorted(raw["patient"].unique())

    # Build per-band (patient × k) matrices of Δρ, ρ_task, ρ_inert
    per_band: dict[str, dict[str, np.ndarray]] = {}
    for band in BRAIN_BANDS_NAMES:
        bdf = raw[raw["band"] == band]
        d  = np.full((len(patients), len(ks)), np.nan)
        rt = np.full_like(d, np.nan)
        ri = np.full_like(d, np.nan)
        for ip, pat in enumerate(patients):
            sub = bdf[bdf["patient"] == pat].set_index("k")
            for jk, k in enumerate(ks):
                if k in sub.index:
                    d[ip, jk]  = sub.at[k, "delta_rho"]
                    rt[ip, jk] = sub.at[k, "rho_task"]
                    ri[ip, jk] = sub.at[k, "rho_inert"]
        per_band[band] = {"delta": d, "rho_task": rt, "rho_inert": ri}

    fig = plt.figure(figsize=(12.0, 9.0), dpi=160)
    gs = fig.add_gridspec(2, 2, width_ratios=[3.2, 1.0],
                          height_ratios=[1.0, 1.0], wspace=0.25, hspace=0.35)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])]
    ax_rank = fig.add_subplot(gs[:, 1])

    # ── Panel (a): fraction of patients with Δρ(k) > 0 ───────────────
    ax = axes[0]
    for band in BRAIN_BANDS_NAMES:
        curve = _frac_positive_curve(per_band[band]["delta"])
        ax.plot(ks, curve, color=BAND_COLORS[band], lw=2.2,
                label=BRAIN_BAND_TEX_DICT[band])
    ax.axhline(0.5, color="#888", lw=0.8, linestyle="--")
    ax.axhline(7/9, color="#000", lw=0.5, linestyle=":")
    ax.axhline(9/9, color="#000", lw=0.5, linestyle=":")
    ax.text(ks[-1] + 0.5, 7/9, "7/9", fontsize=8, va="center")
    ax.text(ks[-1] + 0.5, 9/9, "9/9", fontsize=8, va="center")
    ax.text(ks[-1] + 0.5, 0.5, "chance", fontsize=8, va="center", color="#666")
    ax.set_ylabel(r"fraction of patients with $\Delta\rho(k) > 0$", fontsize=11)
    ax.set_ylim(0.0, 1.05)
    ax.set_title(
        "(a) Per-k patient-unanimity of block persistence, per band "
        "(n=9). Bands that carry a multiscale trace sit high across "
        "a coherent k range; ergodic bands stay near chance.",
        loc="left", fontsize=10,
    )
    ax.legend(ncol=6, loc="lower right", frameon=False, fontsize=10)

    # ── Panel (b): stronger criterion — ρ_task > 2·ρ_inert ───────────
    ax = axes[1]
    for band in BRAIN_BANDS_NAMES:
        curve = _frac_strong_curve(per_band[band]["rho_task"],
                                   per_band[band]["rho_inert"])
        ax.plot(ks, curve, color=BAND_COLORS[band], lw=2.2,
                label=BRAIN_BAND_TEX_DICT[band])
    ax.axhline(0.5, color="#888", lw=0.8, linestyle="--")
    ax.axhline(7/9, color="#000", lw=0.5, linestyle=":")
    ax.set_xlabel("k (dendrogram cut scale)", fontsize=11)
    ax.set_ylabel(r"fraction with $\rho_{task} > 2\,\rho_{inert}$",
                   fontsize=11)
    ax.set_ylim(0.0, 1.05)
    ax.set_title(
        "(b) Stronger criterion: fraction of patients whose task-induced "
        "pairs persist at >2× the baseline rate. Separates trace-bands "
        "from ergodic bands more sharply.",
        loc="left", fontsize=10,
    )

    # ── Right panel: band ranking by trace strength ──────────────────
    strength = {}
    for band in BRAIN_BANDS_NAMES:
        c = _frac_strong_curve(per_band[band]["rho_task"],
                               per_band[band]["rho_inert"])
        strength[band] = float(np.nanmean(c))
    ranked = sorted(strength.items(), key=lambda kv: kv[1], reverse=True)
    y = np.arange(len(ranked))
    vals = [v for _, v in ranked]
    colors = [BAND_COLORS[b] for b, _ in ranked]
    labels = [BRAIN_BAND_TEX_DICT[b] for b, _ in ranked]
    ax_rank.barh(y, vals, color=colors, edgecolor="black", linewidth=0.6)
    ax_rank.set_yticks(y)
    ax_rank.set_yticklabels(labels, fontsize=12)
    ax_rank.invert_yaxis()
    ax_rank.axvline(0.5, color="#888", lw=0.8, linestyle="--")
    ax_rank.set_xlim(0.0, 1.0)
    ax_rank.set_xlabel("mean fraction across k", fontsize=10)
    ax_rank.set_title(
        "Band ranking\nby trace strength\n"
        r"(mean over k of panel b)",
        fontsize=10,
    )
    for yi, v in zip(y, vals):
        ax_rank.text(v + 0.02, yi, f"{v:.2f}", va="center", fontsize=9)

    # No fig.suptitle on publication figures — context lives in the file name.
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"band_unanimity_curves.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
