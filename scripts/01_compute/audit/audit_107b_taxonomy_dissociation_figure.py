#!/usr/bin/env python3
"""audit_107b — figures for the cross-phase taxonomy dissociation.

Reads the audit_105 decomposition + audit_107 localization CSVs and renders:
  fig 1 (dissociation): per band, signal x system heatmap of the demeaned
        per-system endpoint mean M_obs, with geometry-baseline q<0.05 cells
        outlined. Shows whether trace / reset / anchor live in DIFFERENT systems.
  fig 2 (robustness): non-shaft vs shaft-collapsed top-system q per channel.

PDF only, vector, use_lrg_style, no suptitle, figure-level legend. Diverging
color is outline-gated by significance so the (neutral) white centre is not
load-bearing (feedback_no_near_white_cmaps).

Usage: python audit_107b_taxonomy_dissociation_figure.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
OUT = ROOT / "data" / "audit" / "cross_phase_taxonomy"
FIG = OUT / "figures"

SIGNAL_ORDER = ["anchor", "trace", "persist", "reset", "reset_dom"]
SIGNAL_TEX = {"anchor": "anchor\n(rigid)", "trace": "trace\n(concordance)",
              "persist": r"persist ($\varphi_1^2$)", "reset": r"reset ($\varphi_2^2$)",
              "reset_dom": r"reset$-$trace"}
SYSTEM_ORDER = ["OFC", "MTL", "cingulate", "insula", "lateral_temporal", "PFC",
                "sensorimotor", "parietal", "occipital", "subcortical_other", "other"]


def heatmap_band(ax, df, band, gran="system"):
    sub = df[(df.band == band) & (df.granularity == gran)]
    systems = [s for s in SYSTEM_ORDER if s in set(sub.unit)]
    M = np.full((len(SIGNAL_ORDER), len(systems)), np.nan)
    Q = np.ones_like(M)
    for i, sg in enumerate(SIGNAL_ORDER):
        for j, sy in enumerate(systems):
            row = sub[(sub.signal == sg) & (sub.unit == sy)]
            if not row.empty:
                M[i, j] = float(row.M_obs.iloc[0])
                Q[i, j] = float(row.q_geom.iloc[0])
    vmax = np.nanmax(np.abs(M)) if np.isfinite(M).any() else 1.0
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for i in range(len(SIGNAL_ORDER)):
        for j in range(len(systems)):
            if Q[i, j] < 0.05:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor="black", lw=2.0))
            if Q[i, j] < 0.10 and not np.isnan(M[i, j]):
                ax.text(j, i, "*" if Q[i, j] < 0.05 else "·", ha="center",
                        va="center", fontsize=9,
                        color="black" if abs(M[i, j]) < 0.6 * vmax else "white")
    ax.set_xticks(range(len(systems)))
    ax.set_xticklabels(systems, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(SIGNAL_ORDER)))
    ax.set_yticklabels([SIGNAL_TEX[s] for s in SIGNAL_ORDER], fontsize=8)
    ax.set_title(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=11)
    return im


def main():
    use_lrg_style()
    inc = OUT / "localization_include.csv"
    if not inc.exists():
        print(f"missing {inc}; run audit_107 first")
        return
    df = pd.read_csv(inc)
    bands = [b for b in ["beta", "alpha", "low_gamma"] if b in set(df.band)]

    fig, axes = plt.subplots(1, len(bands), figsize=(4.2 * len(bands), 3.6))
    if len(bands) == 1:
        axes = [axes]
    im = None
    for ax, band in zip(axes, bands):
        im = heatmap_band(ax, df, band)
    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label(r"demeaned per-system endpoint mean $M_{\mathrm{obs}}$")
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "dissociation_heatmap.pdf", transparent=True,
                bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG/'dissociation_heatmap.pdf'}")

    # robustness: non-shaft vs shaft-collapsed (beta), per channel best q
    sc = OUT / "localization_shaftcollapsed_include.csv"
    if sc.exists():
        dsc = pd.read_csv(sc)
        print("\nbeta per-channel best system (non-shaft -> shaft-collapsed):")
        for sg in SIGNAL_ORDER:
            a = df[(df.band == "beta") & (df.granularity == "system") & (df.signal == sg)]
            b = dsc[(dsc.band == "beta") & (dsc.granularity == "system") & (dsc.signal == sg)]
            if a.empty:
                continue
            ar = a.loc[a.p_geom.idxmin()]
            br = b.loc[b.p_geom.idxmin()] if not b.empty else None
            bs = f"{br.unit}(q={br.q_geom:.3f})" if br is not None else "-"
            print(f"  {sg:10s}: {ar.unit}(q={ar.q_geom:.3f})  ->  {bs}")


if __name__ == "__main__":
    main()
