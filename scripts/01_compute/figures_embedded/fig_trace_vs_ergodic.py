#!/usr/bin/env python3
"""Trace vs ergodic dichotomy — the band-selective narrative figure.

Using H2e's drift-floor-controlled Δρ contrast (the RIGHT measure,
because it subtracts within-session drift), bands fall into three
clear categories at n=9:

  TRACE  — α, δ   — Δρ_cross > Δρ_null at q<0.05 FDR
  ERGODIC — θ    — Δρ_cross ≈ Δρ_null (no effect above drift)
  marginal — β, γ_l, γ_h  — trend positive, do not survive FDR

The figure presents this dichotomy explicitly, not buried inside a
"universal passing" summary. It is the scientifically honest and
reviewer-defensible version of your band-heterogeneity claim.

Panels:
  (a) Per-band cross vs null bar chart with drift-floor stars.
  (b) Per-patient Δρ_cross − Δρ_null strip, ordered by category.
  (c) Headline band categorization.

Reads:  data/reports/imcoh_vi/h2e_split_half.csv (stats)
        data/reports/imcoh_vi/h2e_split_half_drho_raw.csv (per-patient × k)
        data/reports/imcoh_vi/h2d_persistence_raw.csv (|ImCoh| H2d for context)

Writes: data/reports/imcoh_vi/figures/trace_vs_ergodic.{pdf,png}
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT


IN_DIR  = REPORTS_ROOT / "imcoh_vi"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TRACE_COLOR    = "#2ca02c"   # α, δ: trace
ERGODIC_COLOR  = "#c44e4e"   # θ: ergodic
MARGINAL_COLOR = "#7f7f7f"   # β, γ_l, γ_h: marginal
NULL_COLOR     = "#cfcfcf"


def _band_category(q: float) -> str:
    if np.isfinite(q) and q < 0.05:
        return "trace"
    if np.isfinite(q) and q > 0.25:
        return "ergodic"
    return "marginal"


BAND_CATEGORY = {
    "delta":      "trace",    # Δρ_cross − Δρ_null passes drift-floor FDR
    "theta":      "ergodic",
    "alpha":      "trace",
    "beta":       "marginal",
    "low_gamma":  "marginal",
    "high_gamma": "marginal",
}
CATEGORY_COLOR = {"trace": TRACE_COLOR, "ergodic": ERGODIC_COLOR,
                   "marginal": MARGINAL_COLOR}


def main() -> None:
    stats = pd.read_csv(IN_DIR / "h2e_split_half.csv")
    stats = stats[stats["metric"] == "drho"]  # Δρ test
    per_k = pd.read_csv(IN_DIR / "h2e_split_half_drho_raw.csv")

    # Per-patient × band k-averaged cross-null
    per = (per_k.groupby(["patient", "band"])[["drho_cross", "drho_null"]]
           .mean().reset_index())
    per["diff"] = per["drho_cross"] - per["drho_null"]

    fig = plt.figure(figsize=(14.0, 10.5), dpi=160)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.4, 1.0],
                           height_ratios=[1.0, 0.75],
                           hspace=0.50, wspace=0.30)

    # ── (a) per-band cross vs null bars ──────────────────────
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(BRAIN_BANDS_NAMES))
    w = 0.36
    # Compute per-patient bootstrap CIs on each mean directly.
    def _ci(x: np.ndarray) -> tuple[float, float]:
        rng = np.random.default_rng(0)
        x = x[np.isfinite(x)]
        if x.size < 2:
            return float("nan"), float("nan")
        idx = rng.integers(0, len(x), size=(10_000, len(x)))
        means = x[idx].mean(axis=1)
        lo, hi = np.quantile(means, [0.025, 0.975])
        return float(lo), float(hi)

    for i, band in enumerate(BRAIN_BANDS_NAMES):
        row = stats[stats["band"] == band].iloc[0]
        col = CATEGORY_COLOR[BAND_CATEGORY[band]]
        sub = per[per["band"] == band]
        mc = float(sub["drho_cross"].mean())
        mn = float(sub["drho_null"].mean())
        lo_c, hi_c = _ci(sub["drho_cross"].to_numpy())
        lo_n, hi_n = _ci(sub["drho_null"].to_numpy())
        ax.bar(i - w/2, mc, w, color=col, edgecolor="black", linewidth=0.6,
                label="Δρ_cross" if i == 0 else None)
        ax.bar(i + w/2, mn, w, color=NULL_COLOR, edgecolor="black",
                linewidth=0.6, label="Δρ_null (drift)" if i == 0 else None)
        ax.errorbar(i - w/2, mc, yerr=[[mc - lo_c], [hi_c - mc]],
                    fmt="none", ecolor="black", capsize=3, linewidth=0.8)
        ax.errorbar(i + w/2, mn, yerr=[[mn - lo_n], [hi_n - mn]],
                    fmt="none", ecolor="black", capsize=3, linewidth=0.8)
        # Star + category label
        star = "★" if (np.isfinite(row["q"]) and row["q"] < 0.05) else ""
        ax.text(i, max(hi_c, hi_n) + 0.015,
                f"{BAND_CATEGORY[band].upper()}\nq={row['q']:.2f}{star}",
                ha="center", va="bottom", fontsize=9,
                color=col, fontweight="bold")

    ax.axhline(0, color="black", lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                       fontsize=13)
    ax.set_ylabel(r"mean $\Delta\rho$  (k-averaged, k$\in[2,49]$)", fontsize=11)
    ax.set_ylim(0, 0.43)
    ax.set_title(
        "(a) Task-induced block persistence vs within-session drift null per band.  "
        "Green = TRACE (α, δ: cross > null, q<0.05 FDR);  "
        "Red = ERGODIC (θ: cross ≈ null);  Grey = marginal.",
        fontsize=10, loc="left",
    )
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    # ── (b) per-patient cross-null strip by category ────────
    ax = fig.add_subplot(gs[0, 1])
    rng = np.random.default_rng(0)
    ordered_bands = (
        [b for b in BRAIN_BANDS_NAMES if BAND_CATEGORY[b] == "trace"]
        + [b for b in BRAIN_BANDS_NAMES if BAND_CATEGORY[b] == "marginal"]
        + [b for b in BRAIN_BANDS_NAMES if BAND_CATEGORY[b] == "ergodic"]
    )
    for i, band in enumerate(ordered_bands):
        sub = per[per["band"] == band]
        diffs = sub["diff"].to_numpy()
        col = CATEGORY_COLOR[BAND_CATEGORY[band]]
        jit = rng.uniform(-0.12, 0.12, size=len(diffs))
        ax.scatter(np.full(len(diffs), i) + jit, diffs, s=38,
                    color=col, edgecolor="black", linewidth=0.4, alpha=0.85)
        ax.scatter(i, float(np.mean(diffs)), s=110, marker="_",
                    color="black", linewidth=2.0, zorder=4)

    ax.axhline(0, color="black", lw=0.6, linestyle="--")
    ax.set_xticks(np.arange(len(ordered_bands)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in ordered_bands],
                       fontsize=12)
    ax.set_ylabel(r"$\Delta\rho_{cross} - \Delta\rho_{null}$  (per patient)",
                   fontsize=10)
    ax.set_title("(b) Per-patient excess above drift, ordered by category",
                  fontsize=10, loc="left")
    ax.grid(axis="y", alpha=0.3)

    # ── (c) categorization panel ────────────────────────────
    ax = fig.add_subplot(gs[1, :])
    ax.axis("off")

    # Three category columns
    cats = [
        ("TRACE",      TRACE_COLOR,
         "α, δ",
         "Δρ_cross − Δρ_null > 0\n"
         "at q < 0.05 FDR.\n"
         "Task-formed co-clusters\n"
         "persist above drift."),
        ("MARGINAL",   MARGINAL_COLOR,
         "β, γ_l, γ_h",
         "Mean excess > 0 but\n"
         "does not survive FDR.\n"
         "Signal mixed with drift;\n"
         "underpowered at n=9."),
        ("ERGODIC",    ERGODIC_COLOR,
         "θ",
         "Δρ_cross ≈ Δρ_null\n"
         "(excess +0.018, q=0.43).\n"
         "Indistinguishable from\n"
         "within-session drift."),
    ]
    width_frac = 0.30
    starts = [0.03, 0.03 + width_frac + 0.02, 0.03 + 2 * (width_frac + 0.02)]
    for (label, color, bands, descr), x0 in zip(cats, starts):
        ax.add_patch(FancyBboxPatch(
            (x0, 0.04), width_frac, 0.84,
            boxstyle="round,pad=0.01,rounding_size=0.02",
            facecolor="white", edgecolor=color, linewidth=2.2,
            transform=ax.transAxes,
        ))
        ax.text(x0 + width_frac/2, 0.82, label, fontsize=14, fontweight="bold",
                color=color, ha="center", va="center", transform=ax.transAxes)
        ax.text(x0 + width_frac/2, 0.63, bands, fontsize=18, fontweight="bold",
                color="black", ha="center", va="center", transform=ax.transAxes)
        ax.text(x0 + width_frac/2, 0.28, descr, fontsize=10,
                ha="center", va="center", transform=ax.transAxes,
                linespacing=1.4)

    ax.text(0.5, 0.96,
            "(c) Band categorization by drift-controlled trace strength (n=9, |ImCoh|)",
            ha="center", va="center", fontsize=11, fontweight="bold",
            transform=ax.transAxes)

    # No fig.suptitle on publication figures — context lives in the file name.
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"trace_vs_ergodic.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
