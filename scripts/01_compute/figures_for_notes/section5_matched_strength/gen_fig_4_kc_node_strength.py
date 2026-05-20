#!/usr/bin/env python3
"""§5.7 Figure 4 — KC tree-distance decomposition with matched-strength null
rescaling.

Re-bases the figure on the audit_round3 `kc_decomposition.pdf` two-panel
layout:

  • Left — band points on the (−T_KC(λ=0), −T_KC(λ=1)) plane with the
    trace quadrant shaded green. Each band carries two dots connected
    by an arrow:
        ◦ a hollow grey dot at the matched-strength surrogate location
          (cohort-median surrogate −T_KC at each λ) — what node-strength
          preservation alone reproduces.
        ◦ a filled dot at the observed −T_KC location; the arrow shows
          the wiring-specific *residual*. Fill is rust if the cohort-
          paired Wilcoxon "obs < surr" is significant at any λ.
  • Right — per-band grouped bars `n_below_own_surrogate / 10` at λ=0
    and λ=1 with asterisks above cells that survive cohort-paired
    matched-strength.

The figure makes the §5.7 prose claim visible directly: the surrogate
ghost dot is "trace reproduced by strength scrambling" and the arrow
to the observed dot is the wiring-specific component on top.

Inputs:
    data/audit/kc_matched_strength_surrogate/cohort_summary.csv

Output:
    data/reports/section_5_matched_strength_refinement/figures/
        fig_4_kc_decomposition_matched_strength.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.lines as mlines
import matplotlib.pyplot as plt

from lrg_eegfc.visuals.styles import use_lrg_style

from _shared_ms import (
    ALL_BANDS, BAND_LABELS, FIG_DIR, KC_DIR,
)


LAMBDAS = [0.0, 1.0]


def asterisks(p: float) -> str:
    if not np.isfinite(p):
        return ""
    if p <= 0.001:
        return "***"
    if p <= 0.01:
        return "**"
    if p <= 0.05:
        return "*"
    return ""


def main() -> Path:
    use_lrg_style()

    cohort = pd.read_csv(KC_DIR / "cohort_summary.csv")
    bands = [b for b in ALL_BANDS if b in set(cohort.band)]

    # Build per-band record: observed and surrogate −T_KC, plus p-values.
    recs = []
    for b in bands:
        r0 = cohort[(cohort.band == b)
                    & np.isclose(cohort["lambda"], 0.0)].iloc[0]
        r1 = cohort[(cohort.band == b)
                    & np.isclose(cohort["lambda"], 1.0)].iloc[0]
        T0_obs = -float(r0.obs_median_T_KC)
        T1_obs = -float(r1.obs_median_T_KC)
        T0_surr = -float(r0.surr_median_T_KC_per_patient_median)
        T1_surr = -float(r1.surr_median_T_KC_per_patient_median)
        p0 = float(r0.paired_wilcoxon_p)
        p1 = float(r1.paired_wilcoxon_p)
        # n_patients_below_own_surrogate is stored as "k/10" string.
        n0 = int(str(r0.n_patients_below_own_surrogate).split("/")[0])
        n1 = int(str(r1.n_patients_below_own_surrogate).split("/")[0])
        recs.append({
            "band": b,
            "T0_obs": T0_obs, "T1_obs": T1_obs,
            "T0_surr": T0_surr, "T1_surr": T1_surr,
            "p0": p0, "p1": p1, "p_min": min(p0, p1),
            "n0_below": n0, "n1_below": n1,
        })
    R = pd.DataFrame(recs)

    fig = plt.figure(figsize=(13.2, 5.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.4, 1.0], wspace=0.30)
    ax = fig.add_subplot(gs[0, 0])

    xmax = max(R["T0_obs"].abs().max(), R["T0_surr"].abs().max(), 0.5) * 1.32
    ymax = max(R["T1_obs"].abs().max(), R["T1_surr"].abs().max(), 0.5) * 1.32

    # Trace quadrant background (NE quadrant in the −T plane)
    ax.axhspan(0, ymax, xmin=0.5, xmax=1.0, color="#1a7c3e", alpha=0.13)
    ax.axhline(0, color="0.4", lw=0.7, ls="--")
    ax.axvline(0, color="0.4", lw=0.7, ls="--")

    # Per-band annotation offsets in axes-fraction (matches reference).
    ann_offsets = {
        "delta":      (+0.025, -0.030),
        "theta":      (+0.025, +0.025),
        "alpha":      (+0.025, -0.040),
        "beta":       (+0.025, -0.035),
        "low_gamma":  (+0.025, +0.025),
    }

    for _, row in R.iterrows():
        is_sig = (row["p0"] <= 0.05) or (row["p1"] <= 0.05)
        face = "#c0392b" if is_sig else "white"
        edge = "#c0392b" if is_sig else "#444"
        size = 80 + 220.0 * (-np.log10(max(row["p_min"], 1e-3)) / 2.0)
        size = float(np.clip(size, 80, 360))

        # surrogate ghost — small hollow grey dot
        ax.scatter(row["T0_surr"], row["T1_surr"], s=64,
                   facecolor="#dddddd", edgecolor="#777777", lw=1.0,
                   zorder=2)
        # arrow surrogate → observed: the wiring-specific residual
        ax.annotate("",
                    xy=(row["T0_obs"], row["T1_obs"]),
                    xytext=(row["T0_surr"], row["T1_surr"]),
                    arrowprops=dict(
                        arrowstyle="->", color="#444",
                        lw=0.9, alpha=0.85,
                        shrinkA=4, shrinkB=4,
                    ),
                    zorder=2.5)
        # observed point
        ax.scatter(row["T0_obs"], row["T1_obs"], s=size,
                   facecolor=face, edgecolor=edge, lw=1.4, zorder=3)

        dx_frac, dy_frac = ann_offsets.get(row["band"], (0.05, 0.05))
        dx_data = dx_frac * 2 * xmax
        dy_data = dy_frac * 2 * ymax
        ha = "left" if dx_frac >= 0 else "right"
        va = "bottom" if dy_frac >= 0 else "top"
        weight = "bold" if is_sig else "normal"
        fs = 13 if is_sig else 11
        ax.annotate(BAND_LABELS[row["band"]],
                    (row["T0_obs"], row["T1_obs"]),
                    xytext=(row["T0_obs"] + dx_data,
                            row["T1_obs"] + dy_data),
                    fontsize=fs, fontweight=weight, ha=ha, va=va,
                    color="#222")

    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-ymax, ymax)
    ax.set_xlabel(r"$-T_{KC}(\lambda=0)$ — pure topology, "
                  r"trace direction $\rightarrow$", fontsize=10)
    ax.set_ylabel(r"$-T_{KC}(\lambda=1)$ — pure heights, "
                  r"$\uparrow$ trace direction", fontsize=10)

    ax.text(0.04, 0.96, "heights only", ha="left", va="top",
            fontsize=9, color="#444", fontstyle="italic",
            transform=ax.transAxes)
    ax.text(0.96, 0.04, "topology only", ha="right", va="bottom",
            fontsize=9, color="#444", fontstyle="italic",
            transform=ax.transAxes)
    ax.text(0.04, 0.04, "anti-trace", ha="left", va="bottom",
            fontsize=9, color="#c0392b", fontstyle="italic",
            transform=ax.transAxes)
    ax.spines[["top", "right"]].set_visible(False)

    handles = [
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=8,
                      markerfacecolor="#dddddd", markeredgecolor="#777",
                      label="matched-strength surrogate  "
                            r"$-T_{KC}^{\mathrm{p50}}$"),
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=9,
                      markerfacecolor="white", markeredgecolor="#444",
                      label=r"observed $-T_{KC}$  ($p > 0.05$)"),
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=9,
                      markerfacecolor="#c0392b", markeredgecolor="#c0392b",
                      label=r"observed $-T_{KC}$  ($p \leq 0.05$)"),
        mlines.Line2D([], [], color="#444", lw=1.0,
                      label="arrow $=$ wiring-specific residual"),
    ]
    ax.legend(handles=handles, loc="lower right",
              bbox_to_anchor=(1.0, 0.10),
              frameon=False, ncol=1,
              handletextpad=0.5, labelspacing=0.5,
              fontsize=8)

    # ---- right panel: n_below_own_surrogate bars --------------------------
    ax_r = fig.add_subplot(gs[0, 1])
    x = np.arange(len(R))
    w = 0.36
    band_list = list(R["band"])

    col_lam0 = ["#c0392b" if R.iloc[j]["p0"] <= 0.05 else "#4C72B0"
                 for j in range(len(R))]
    col_lam1 = ["#c0392b" if R.iloc[j]["p1"] <= 0.05 else "#DD8452"
                 for j in range(len(R))]
    edge_lam0 = ["#7d1a0e" if R.iloc[j]["p0"] <= 0.05 else "#2F4A78"
                  for j in range(len(R))]
    edge_lam1 = ["#7d1a0e" if R.iloc[j]["p1"] <= 0.05 else "#9C5B36"
                  for j in range(len(R))]
    bars0 = ax_r.bar(x - w / 2, R["n0_below"], width=w,
                     color=col_lam0, edgecolor=edge_lam0, lw=0.8,
                     label=r"$\lambda=0$ (topology / node-strength proxy)")
    bars1 = ax_r.bar(x + w / 2, R["n1_below"], width=w,
                     color=col_lam1, edgecolor=edge_lam1, lw=0.8,
                     label=r"$\lambda=1$ (heights / wiring-specific)")

    # Reference line at n = 5
    ax_r.axhline(5, color="0.4", lw=0.6, ls="--", zorder=1)
    ax_r.text(-0.55, 5, "n = 5", fontsize=8, color="0.4",
              ha="right", va="center")
    # Numerical labels above each bar
    for bar, v in zip(list(bars0) + list(bars1),
                      list(R["n0_below"]) + list(R["n1_below"])):
        ax_r.text(bar.get_x() + bar.get_width() / 2,
                  bar.get_height() + 0.18,
                  str(int(v)), ha="center", va="bottom", fontsize=8.5,
                  color="0.15")

    # Asterisks above cells with cohort-paired p ≤ 0.05.
    for j in range(len(R)):
        for off, p in [(-w / 2, R.iloc[j]["p0"]), (+w / 2, R.iloc[j]["p1"])]:
            stars = asterisks(p)
            if stars:
                ax_r.text(j + off, 9.6, stars, ha="center", va="bottom",
                          fontsize=12, fontweight="bold", color="#7d1a0e")

    ax_r.set_xticks(x)
    ax_r.set_xticklabels([BAND_LABELS[b] for b in R["band"]], fontsize=12)
    ax_r.set_ylim(0, 10.5)
    ax_r.set_ylabel(r"$n_{<\,\mathrm{own\;surrogate}}\,/\,10$",
                    fontsize=11)
    ax_r.legend(fontsize=8.5, frameon=False, loc="upper left",
                bbox_to_anchor=(0.0, 1.0))
    ax_r.spines[["top", "right"]].set_visible(False)
    ax_r.tick_params(axis="x", which="both", length=0)

    fig.tight_layout()
    out = FIG_DIR / "fig_4_kc_decomposition_matched_strength.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
