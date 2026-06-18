#!/usr/bin/env python3
"""Epi-marker classifier verdict — LOPO discrimination + whole-shaft masking.

Two panels, all six bands.

(a) **LOPO cross-patient discrimination** (PR-AUC; chance = prevalence, dashed).
    Grouped bars per band for the trace model vs the transferable baselines
    (node strength, label-free electrode geometry) and the full model. Reading:
    *the trace bar sitting at the prevalence line with electrode geometry well
    above it is the whole result — epi recoverability is geometric, not a
    cross-phase-trace marker.*

(b) **Whole-shaft masking recovery** (ROC; chance = 0.5, dashed). Grouped bars
    per band for trace / strength / geometry / the label tautology. Reading:
    *the label-tautology bar collapsing below 0.5 is the probe-shortcut control
    firing; geometry stays high (label-free); trace hugs chance.*

No in-axes numeric text (project rule); the numbers live in
`data/audit/epi_marker/clf_*` and `README_classifier.md`.

Input  : data/audit/epi_marker/{clf_lopo_per_band.csv, clf_masking_per_band.csv}
Output : data/preprint/figures/all_bands/fig_epi_marker_verdict.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

# consistent, saturated colours (no near-white interiors — feedback rule)
C = {
    "trace": "#d1495b",        # crimson — the candidate marker
    "strength": "#2c5f8a",     # blue — strength baseline
    "probe_geom": "#3a7d44",   # green — label-free electrode geometry
    "all": "#6a4c93",          # purple — full model
    "probe_label": "#8a8a8a",  # grey — leaky tautology
}
LBL = {
    "trace": "trace (z ρ, z coh, z mode-load)",
    "strength": "node strength",
    "probe_geom": "electrode geometry (label-free)",
    "all": "all features",
    "probe_label": "shaft tautology (label-leaky)",
}
PANEL_A = ["trace", "strength", "probe_geom", "all"]      # LOPO PR-AUC
PANEL_B = ["trace", "strength", "probe_geom", "probe_label"]  # masked ROC


def _grouped_bars(ax, df_lookup, models, value_key, color):
    x = np.arange(len(BANDS))
    w = 0.8 / len(models)
    for mi, m in enumerate(models):
        vals = [df_lookup(b, m, value_key) for b in BANDS]
        ax.bar(x + (mi - (len(models) - 1) / 2) * w, vals, width=w,
               color=color[m], edgecolor="none", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS])
    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(x=0.01)


def main() -> None:
    base = ROOT / "data" / "audit" / "epi_marker"
    lopo = pd.read_csv(base / "clf_lopo_per_band.csv")
    mask = pd.read_csv(base / "clf_masking_per_band.csv")
    if lopo.empty or mask.empty:
        raise SystemExit("[preprint_33] empty inputs; run audit_80 first")

    prevalence = float(lopo["prevalence"].iloc[0])

    def lopo_pr(b, m, _):
        r = lopo[lopo.band == b]
        return float(r[f"{m}_pr_auc"].iloc[0]) if not r.empty else np.nan

    def mask_roc(b, m, _):
        r = mask[(mask.band == b) & (mask.feature_set == m)]
        return float(r["masked_roc_auc"].iloc[0]) if not r.empty else np.nan

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(13.2, 4.8))

    # (a) LOPO PR-AUC
    _grouped_bars(axa, lopo_pr, PANEL_A, "pr_auc", C)
    axa.axhline(prevalence, color="0.35", lw=1.1, ls="--", zorder=4)
    axa.set_ylabel("LOPO PR-AUC")
    axa.set_ylim(0, max(0.25, lopo[[f"{m}_pr_auc" for m in PANEL_A]]
                        .to_numpy().max() * 1.15))
    axa.text(-0.45, prevalence, "chance", color="0.35", va="bottom",
             ha="left", fontsize=9)
    axa.set_title("(a) cross-patient discrimination", fontsize=11, loc="left")

    # (b) masked ROC
    _grouped_bars(axb, mask_roc, PANEL_B, "masked_roc_auc", C)
    axb.axhline(0.5, color="0.35", lw=1.1, ls="--", zorder=4)
    axb.set_ylabel("whole-shaft-masked ROC-AUC")
    axb.set_ylim(0.3, 0.75)
    axb.text(-0.45, 0.5, "chance", color="0.35", va="bottom", ha="left",
             fontsize=9)
    axb.set_title("(b) whole-shaft masking recovery", fontsize=11, loc="left")

    # figure-level legend (union of both panels' models, stable order)
    order = ["trace", "strength", "probe_geom", "all", "probe_label"]
    handles = [Patch(facecolor=C[m], edgecolor="none", label=LBL[m])
               for m in order]
    handles.append(Line2D([], [], color="0.35", lw=1.1, ls="--",
                          label="chance (prevalence / 0.5)"))
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.10))
    fig.tight_layout(rect=(0, 0.04, 1, 1))

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_epi_marker_verdict.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[preprint_33] -> {out}")


if __name__ == "__main__":
    main()
