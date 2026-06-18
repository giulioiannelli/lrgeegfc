#!/usr/bin/env python3
"""preprint_38 — the deployment marker: given a few labelled SOZ, predict the rest.

One PDF (fig_epi_seed_prediction.pdf), three honest panels:

(a) CONCRETE prediction (Pat_08, δ, 3 seeds): non-seed contacts ranked by the
    strength-residual propagator-affinity score; filled = true remaining SOZ (hit),
    open = false positive; dashed line = top-m cutoff (m = #remaining SOZ). Shows the
    actual prediction, not just an AUC — and that the recovered SOZ sit far from the
    seeds (distance annotated), so it is not trivial proximity.
(b) TRIAGE curve: cohort-median precision@top-m vs seed count k (2/3/5/8) per band,
    against the prevalence floor (chance) and the node-strength baseline. The honest
    "how many seeds do you need" read — δ/low-γ/β enrich 3–5×, α ≈ chance.
(c) TWO-POPULATION split (δ, k=3): per-patient propagator precision@top-m vs the
    node-strength precision. Above the diagonal = propagator wins (the SOZ is a
    diffusion community); below = strength wins (the SOZ are hubs — Pat_15/10/07).

Inputs (read-only) data/audit/epi_seed_prediction/:
    seed_prediction_showcase_k3.csv, seed_prediction_per_patient_k3.csv
  + data/audit/epi_marker_relational/seed_curve_cohort.csv
PDF only, full vector, transparent.
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

SP = ROOT / "data/audit/epi_seed_prediction"
SC = ROOT / "data/audit/epi_marker_relational"
OUT = ROOT / "data/preprint/figures/all_bands"
OUT.mkdir(parents=True, exist_ok=True)

BANDS = ["delta", "beta", "low_gamma", "alpha"]
BAND_C = {"delta": "#1b7837", "beta": "#2166ac", "low_gamma": "#762a83", "alpha": "#b35806"}
C_HIT = "#1b7837"     # true SOZ recovered
C_MISS = "#bdbdbd"    # false positive
C_STR = "#b2182b"     # strength baseline
C_PROP = "#2166ac"    # propagator


def _tex(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


def panel_concrete(ax, pat="Pat_08", band="delta"):
    df = pd.read_csv(SP / "seed_prediction_showcase_k3.csv")
    d = df[(df.patient == pat) & (df.band == band)].sort_values("rank")
    if d.empty:
        ax.set_visible(False)
        return None
    m = int(d["m"].iloc[0])
    y = np.arange(len(d))[::-1]                 # rank 1 at top
    col = [C_HIT if t else C_MISS for t in d.is_SOZ]
    ax.hlines(y, 0, d.score, color=col, lw=2.2, zorder=1)
    ax.scatter(d.score, y, c=col, s=46, zorder=2, edgecolor="black", linewidth=0.4)
    for yi, (_, r) in zip(y, d.iterrows()):
        tag = f"{r.label.strip()}  ({r.dist_to_seed_mm:.0f} mm)" if r.is_SOZ else r.label.strip()
        ax.text(d.score.max() * 1.04, yi, tag, va="center", fontsize=6.2,
                color=(C_HIT if r.is_SOZ else "0.45"))
    cut = len(d) - m - 0.5
    ax.axhline(cut, ls="--", color="black", lw=0.8)
    ax.text(d.score.max() * 0.5, cut + 0.4, f"top-{m} cutoff", fontsize=6.5, va="bottom")
    ax.set_yticks([])
    ax.set_xlabel("strength-residual propagator affinity to seeds")
    ax.set_xlim(0, d.score.max() * 1.5)
    ax.text(0.02, 1.02, f"a  {pat}, {_tex(band)}: 3 seeds → ranked SOZ prediction",
            transform=ax.transAxes, va="bottom", fontsize=8.5)
    return d


def panel_triage(ax):
    coh = pd.read_csv(SC / "seed_curve_cohort.csv")
    for b in BANDS:
        d = coh[coh.band == b].sort_values("k_seeds")
        if d.empty:
            continue
        ax.plot(d.k_seeds, d.med_precision_at_m, "-o", color=BAND_C[b], ms=4,
                lw=1.8, label=_tex(b))
    # prevalence floor + strength baseline (band-averaged for a single reference line)
    prev = coh.groupby("k_seeds").med_prevalence.median()
    strv = coh.groupby("k_seeds").med_precision_strength.median()
    ax.plot(prev.index, prev.values, ":", color="black", lw=1.0, label="prevalence (chance)")
    ax.plot(strv.index, strv.values, "--", color=C_STR, lw=1.0, label="strength baseline")
    ax.set_xlabel("number of labelled seed SOZ (k)")
    ax.set_ylabel("precision @ top-m")
    ax.set_xticks([2, 3, 5, 8])
    ax.text(0.02, 1.02, "b  triage value vs seed count", transform=ax.transAxes,
            va="bottom", fontsize=8.5)


def panel_two_pop(ax, band="delta"):
    pp = pd.read_csv(SP / "seed_prediction_per_patient_k3.csv")
    d = pp[pp.band == band].copy()
    lim = max(d.precision_at_m.max(), d.precision_strength.max()) * 1.15 + 1e-3
    ax.plot([0, lim], [0, lim], "-", color="0.7", lw=0.9, zorder=0)
    for _, r in d.iterrows():
        win = r.precision_at_m >= r.precision_strength
        ax.scatter(r.precision_strength, r.precision_at_m, s=42,
                   c=(C_PROP if win else C_STR), edgecolor="black", linewidth=0.4, zorder=2)
        ax.annotate(r.patient.replace("Pat_", ""),
                    (r.precision_strength, r.precision_at_m),
                    textcoords="offset points", xytext=(3, 3), fontsize=6)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("strength precision @ top-m")
    ax.set_ylabel("propagator precision @ top-m")
    ax.text(0.02, 1.02, f"c  two populations ({_tex(band)})", transform=ax.transAxes,
            va="bottom", fontsize=8.5)
    ax.text(0.04, 0.94, "propagator wins\n(SOZ = community)", transform=ax.transAxes,
            ha="left", va="top", fontsize=6.3, color=C_PROP)
    ax.text(0.96, 0.06, "strength wins\n(SOZ = hubs)", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6.3, color=C_STR)


def main():
    fig = plt.figure(figsize=(11.4, 3.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.0], wspace=0.32)
    panel_concrete(fig.add_subplot(gs[0, 0]))
    axb = fig.add_subplot(gs[0, 1])
    panel_triage(axb)
    panel_two_pop(fig.add_subplot(gs[0, 2]))

    h, l = axb.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, -0.06),
               ncol=len(h), frameon=False, fontsize=7.5)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    out = OUT / "fig_epi_seed_prediction.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[preprint_38] -> {out.name}")


if __name__ == "__main__":
    main()
