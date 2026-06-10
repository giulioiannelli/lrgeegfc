#!/usr/bin/env python3
"""preprint_37 — label-free cross-patient SOZ localization: narrows, cannot select.

One PDF (fig_epi_community_localization.pdf), the honest Direction-B result:

(a) SOZ precision per band — prevalence (chance) vs the ORACLE best label-free
    community (within-patient narrowing, ~4–5× lift) vs the CROSS-PATIENT SELECTED
    top community (label-free selection, ≈ chance). Narrowing works; selection fails.
(b) Recall of the oracle best community per band — the single label-free community
    captures ~50–60% of the SOZ (the narrowing payoff), with the SOZ-community #1
    selection rate vs chance (1/K) annotated as the failure of automatic selection.

Inputs (read-only): data/audit/epi_community_localization/{community_features,
        localization_per_patient}.csv. PDF only, full vector, transparent.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

SRC = ROOT / "data/audit/epi_community_localization"
OUT = ROOT / "data/preprint/figures/all_bands"
OUT.mkdir(parents=True, exist_ok=True)

BANDS = ["delta", "alpha", "beta", "low_gamma"]
C_CHANCE = "#bdbdbd"   # prevalence / chance (grey)
C_ORACLE = "#1b7837"   # label-free narrowing (green)
C_SELECT = "#b2182b"   # label-free cross-patient selection (red)


def _tex(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


def main():
    comm = pd.read_csv(SRC / "community_features.csv")
    loc = pd.read_csv(SRC / "localization_per_patient.csv")

    prev, oracle_p, oracle_r, sel_p, soz1, kmed = {}, {}, {}, {}, {}, {}
    for b in BANDS:
        cb = comm[comm.band == b]
        op, orc, pv = [], [], []
        for pat, g in cb.groupby("patient"):
            if g.n_soz.max() < 3:
                continue
            best = g.loc[g.n_soz.idxmax()]
            op.append(best.soz_precision)
            orc.append(best.n_soz / best.n_soz_total)
            pv.append(best.prevalence)
        oracle_p[b] = np.median(op) if op else np.nan
        oracle_r[b] = np.median(orc) if orc else np.nan
        prev[b] = np.median(pv) if pv else np.nan
        lb = loc[loc.band == b]
        sel_p[b] = float(lb.top_precision.median()) if not lb.empty else np.nan
        soz1[b] = float(lb.top_is_soz_comm.mean()) if not lb.empty else np.nan
        kmed[b] = float(lb.K.median()) if not lb.empty else 8

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.3))
    x = np.arange(len(BANDS))
    w = 0.27

    ax = axes[0]
    ax.bar(x - w, [prev[b] for b in BANDS], w, color=C_CHANCE, edgecolor="black",
           linewidth=0.4, label="prevalence (chance)")
    ax.bar(x, [oracle_p[b] for b in BANDS], w, color=C_ORACLE, edgecolor="black",
           linewidth=0.4, label="oracle community (label-free narrowing)")
    ax.bar(x + w, [sel_p[b] for b in BANDS], w, color=C_SELECT, edgecolor="black",
           linewidth=0.4, label="selected community (cross-patient, label-free)")
    ax.set_ylabel("SOZ precision of community")
    ax.set_xticks(x); ax.set_xticklabels([_tex(b) for b in BANDS])
    ax.text(0.03, 0.96, "a  narrowing works, selection fails", transform=ax.transAxes,
            va="top", fontsize=9)

    ax = axes[1]
    ax.bar(x, [oracle_r[b] for b in BANDS], 0.5, color=C_ORACLE, edgecolor="black",
           linewidth=0.4)
    ax.set_ylabel("SOZ recall of best label-free community")
    ax.set_ylim(0, 1)
    ax.set_xticks(x); ax.set_xticklabels([_tex(b) for b in BANDS])
    for xi, b in zip(x, BANDS):
        ax.text(xi, oracle_r[b] + 0.02, f"#1:{soz1[b]:.0%}\nvs {1/kmed[b]:.0%}",
                ha="center", va="bottom", fontsize=6.5, color=C_SELECT)
    ax.text(0.03, 0.96, "b  one community holds the majority of SOZ",
            transform=ax.transAxes, va="top", fontsize=9)

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, -0.04), ncol=3,
               frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    out = OUT / "fig_epi_community_localization.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[preprint_37] -> {out.name}")


if __name__ == "__main__":
    main()
