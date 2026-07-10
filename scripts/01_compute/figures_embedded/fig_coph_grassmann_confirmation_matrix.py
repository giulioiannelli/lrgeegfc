#!/usr/bin/env python3
r"""fig (talk R-4): coph x Grassmann confirmation matrix — the double dissociation.

A 6-band x 2-probe truth table. For each band, does the trace clear the
matched-strength null under EACH of the two strength-independent LRG readouts?

  * cophenetic  (local / tree)     -- ultrametric merge-height reorganization,
                                       the rho_sym cross-phase gate (audit_150).
  * Grassmann   (global / subspace)-- rotation of the dominant diffusion
                                       eigenmodes, cluster-extent gate (audit_66).

Cell = -log10(gate p) in a per-probe colour ramp when the band CLEARS that
probe's matched-strength gate (p < 0.05, marked with a uniform check); pale grey
when it fails. A right-margin verdict labels each row BOTH / local-only /
global-only / neither. beta is the ONLY band in BOTH -- the "why you should
believe it" panel: the two probes genuinely dissociate on every other band.

This REPLACES the stale `preprint_34_bands_coph_grassmann_dissociation_map.py`
COPHENETIC axis, which read the pre-rho_sym `rho_split`
(`matched_strength_surrogate_split_baseline/`). Here the cophenetic gate is the
current rho_sym gate (`rho_sym_gate/cohort_summary.csv`); the dissociation is
unchanged in shape and tighter in the numbers.

Reads : data/audit/rho_sym_gate/cohort_summary.csv          (cophenetic gate_p_sym)
        data/audit/grassmann_cluster_extent/cohort_summary.csv (Grassmann cluster-mass p + LOO)
Writes: data/outputs/figures/talk/fig_coph_grassmann_confirmation_matrix.pdf
"""
from __future__ import annotations

import argparse

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

use_lrg_style()

COPH_SRC = ROOT / "data/audit/rho_sym_gate/cohort_summary.csv"
GRASS_SRC = ROOT / "data/audit/grassmann_cluster_extent/cohort_summary.csv"
OUT = ROOT / "data/outputs/figures/talk/fig_coph_grassmann_confirmation_matrix.pdf"

ALPHA = 0.05
GX = 1.18                          # Grassmann column x-offset (gap separates the two lenses)
R_SURR = 200                       # both gates: R=200 -> empirical floor 1/(R+1)
PMIN = 1.0 / (R_SURR + 1)          # 0.00498 -> -log10 ceiling 2.303
NEGLOG_MAX = -np.log10(PMIN)

# verdict-tiered row order (BOTH -> local -> global -> neither), matches the talk mock
ROW_ORDER = ["beta", "alpha", "low_gamma", "delta", "high_gamma", "theta"]

# per-probe ramps (cool = local/tree, warm = global/subspace) + verdict hues,
# echoing the quadrant-map palette for cross-figure consistency
COPH_CMAP = LinearSegmentedColormap.from_list("coph", ["#dce7f0", "#2c5f8a"])
GRASS_CMAP = LinearSegmentedColormap.from_list("grass", ["#fbe6cf", "#e67e22"])
FAIL_FC = "#eceef0"
VERDICT_COL = {"BOTH": "#1f7a1f", "local only": "#2c5f8a",
               "global only": "#e67e22", "neither": "#8a9096"}
INK = "#2b2f36"


def _neglog(p: float) -> float:
    return -np.log10(max(float(p), PMIN))


def load_rows() -> pd.DataFrame:
    coph = pd.read_csv(COPH_SRC).set_index("band")
    grass = pd.read_csv(GRASS_SRC).set_index("band")
    rows = []
    for b in ROW_ORDER:
        cp = float(coph.loc[b, "gate_p_sym"])
        gp = float(grass.loc[b, "cluster_p_cluster_mass"])
        gloo = float(grass.loc[b, "cluster_p_mass_loo_max"])
        c_ok, g_ok = cp < ALPHA, gp < ALPHA
        verdict = ("BOTH" if c_ok and g_ok else "local only" if c_ok
                   else "global only" if g_ok else "neither")
        rows.append(dict(band=b, coph_p=cp, grass_p=gp, grass_loo=gloo,
                         coph_ok=c_ok, grass_ok=g_ok,
                         grass_fragile=g_ok and gloo >= ALPHA, verdict=verdict))
    return pd.DataFrame(rows)


def _check(ax, cx, cy, color):
    """small vector check-mark (top-right of a cell), font-independent."""
    xs = np.array([cx - 0.055, cx - 0.012, cx + 0.075])
    ys = np.array([cy + 0.010, cy - 0.035, cy + 0.070])
    ax.plot(xs, ys, color=color, lw=2.4, solid_capstyle="round",
            solid_joinstyle="round", zorder=6)


def _cell(ax, col, row_y, p, ok, cmap, fragile, annotate):
    x0 = col
    fc = cmap(0.30 + 0.70 * _neglog(p) / NEGLOG_MAX) if ok else FAIL_FC
    edge = "#c0392b" if fragile else "white"
    ax.add_patch(Rectangle((x0 + 0.04, row_y + 0.04), 0.92, 0.92,
                           facecolor=fc, edgecolor=edge,
                           linewidth=2.0 if fragile else 1.0, zorder=2))
    if ok:
        _check(ax, x0 + 0.80, row_y + 0.66, "white" if _neglog(p) > 1.6 else INK)
    if annotate:
        txt = f"{p:.3f}".lstrip("0")
        tcol = "white" if (ok and _neglog(p) > 1.4) else (INK if ok else "#9aa0a7")
        ax.text(x0 + 0.44, row_y + 0.40, txt, ha="center", va="center",
                fontsize=12.5, fontweight="bold" if ok else "normal",
                color=tcol, zorder=5)
        if fragile:
            ax.text(x0 + 0.44, row_y + 0.17, "LOO-fragile", ha="center",
                    va="center", fontsize=6.6, style="italic",
                    color="#c0392b", zorder=5)


def main(annotate: bool = True) -> None:
    R = load_rows()
    n = len(R)

    fig, ax = plt.subplots(figsize=(6.9, 6.8))
    ax.axis("off")

    # highlight the unique BOTH row (beta)
    for i, r in R.iterrows():
        if r.verdict == "BOTH":
            y = n - 1 - i
            ax.add_patch(FancyBboxPatch((-1.62, y - 0.02), 4.72, 1.04,
                         boxstyle="round,pad=0.02,rounding_size=0.10",
                         facecolor="#fff6da", edgecolor="#e5c454", lw=1.1,
                         zorder=1))

    for i, r in R.iterrows():
        y = n - 1 - i
        _cell(ax, 0, y, r.coph_p, r.coph_ok, COPH_CMAP, False, annotate)
        _cell(ax, GX, y, r.grass_p, r.grass_ok, GRASS_CMAP, r.grass_fragile, annotate)
        # band row label (canonical palette)
        star = r"$\,\bigstar$" if r.verdict == "BOTH" else ""
        ax.text(-0.18, y + 0.5, BRAIN_BAND_TEX_DICT[r.band] + star,
                ha="right", va="center", fontsize=17,
                color=band_color(r.band, 0.95), fontweight="bold", zorder=4)
        # verdict tag (right margin)
        vc = VERDICT_COL[r.verdict]
        ax.text(GX + 1.06, y + 0.5, r.verdict, ha="left", va="center",
                fontsize=11.5, color=vc,
                fontweight="bold" if r.verdict == "BOTH" else "normal", zorder=4)

    # column headers
    ax.text(0.5, n + 0.16, "cophenetic", ha="center", va="bottom",
            fontsize=11.5, color="#2c5f8a", fontweight="bold")
    ax.text(0.5, n + 0.05, "local · tree", ha="center", va="top",
            fontsize=8, color="#5a6068")
    ax.text(GX + 0.5, n + 0.16, "Grassmann", ha="center", va="bottom",
            fontsize=11.5, color="#e67e22", fontweight="bold")
    ax.text(GX + 0.5, n + 0.05, "global · subspace", ha="center", va="top",
            fontsize=8, color="#5a6068")

    # compact key
    fig.text(0.5, 0.028,
             r"$\checkmark$ clears matched-strength null (p < 0.05)   ·   "
             r"shade $\propto -\log_{10}$ gate p   ·   "
             r"only $\beta$ is confirmed by BOTH readouts",
             ha="center", va="center", fontsize=8.2, color="#6a7079")

    ax.set_xlim(-1.70, 3.55)
    ax.set_ylim(-0.15, n + 0.55)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight")
    plt.close(fig)

    # stdout audit trail
    print("\n=== coph (rho_sym) x Grassmann confirmation matrix ===")
    print(f"{'band':<11}{'coph_p':>8}{'coph':>6}{'grass_p':>9}{'grass':>7}"
          f"{'LOO':>7}   verdict")
    for _, r in R.iterrows():
        print(f"{r.band:<11}{r.coph_p:>8.3f}{'  v' if r.coph_ok else '  -':>6}"
              f"{r.grass_p:>9.3f}{'  v' if r.grass_ok else '  -':>7}"
              f"{r.grass_loo:>7.3f}   {r.verdict}"
              f"{'  (LOO-fragile)' if r.grass_fragile else ''}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-annotate", action="store_true",
                    help="hide the per-cell gate p-values (pure visual truth table)")
    args = ap.parse_args()
    main(annotate=not args.no_annotate)
