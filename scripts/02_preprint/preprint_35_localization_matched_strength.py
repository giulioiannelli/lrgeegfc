#!/usr/bin/env python3
"""preprint_35 — matched-strength-verified trace localization figures.

Two PDFs visualizing the verified verdict in
`data/audit/localization_atlas/README.md`:

1. fig_localization_matched_strength.pdf — survivor map. For the 3 trace bands
   (beta, alpha, low_gamma), horizontal bars of -log10(matched-strength p) per
   anatomical SYSTEM, reference line at p=0.05. Bars coloured by trace sign
   (M_obs > 0 = trace direction). Shows beta -> orbitofrontal, alpha -> frontal
   operculum, low_gamma weak.

2. fig_localization_beta_two_null.pdf — why the strength control matters (beta).
   (a) per system, the LABEL-SHUFFLE significance (audit_81) vs the MATCHED-STRENGTH
   significance (audit_83): hippocampus/MTL clear the label shuffle but NOT matched
   strength; OFC clears both. (b) strength-residualization (audit_84): raw vs
   strength-residualized concentration for the key units -- the survivors (and Hip)
   are not strength artifacts (residualization barely changes them).

Inputs: data/audit/localization_atlas/{matched_strength_include.csv,
        atlas_cophenet_include.csv, strength_residual_include.csv}.
Read-only on CSVs; no recompute. PDF only, full vector, transparent.
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

ATLAS = ROOT / "data/audit/localization_atlas"
OUT = ROOT / "data/preprint/figures/all_bands"
OUT.mkdir(parents=True, exist_ok=True)

TRACE_BANDS = ["beta", "alpha", "low_gamma"]
KMIN = 3
P_REF = -np.log10(0.05)

# saturated, non-near-white palette (feedback_no_near_white_cmaps)
C_TRACE = "#1b7837"      # trace-direction, significant (green)
C_ANTI = "#762a83"       # anti-direction, significant (purple)
C_NS = "#bdbdbd"         # non-significant (grey)


def neglog(p):
    return -np.log10(np.clip(np.asarray(p, float), 1e-6, 1.0))


def load_systems(csv, pcol):
    df = pd.read_csv(ATLAS / csv)
    df = df[(df.granularity == "system") & (df.K_implanted >= KMIN)
            & (df.unit != "other")].copy()
    df["nlp"] = neglog(df[pcol])
    return df


# ---------------------------------------------------------------------------
# Figure 1 — matched-strength survivor map across trace bands
# ---------------------------------------------------------------------------
def fig_survivor_map():
    ms = load_systems("matched_strength_include.csv", "matched_strength_p")
    # consistent system order: by beta -log10 p (descending)
    order = (ms[ms.band == "beta"].sort_values("nlp", ascending=True)
             ["unit"].tolist())
    extra = [u for u in ms.unit.unique() if u not in order]
    order = extra + order

    fig, axes = plt.subplots(1, len(TRACE_BANDS),
                             figsize=(3.1 * len(TRACE_BANDS), 4.4), sharey=True)
    for ax, band in zip(axes, TRACE_BANDS):
        sub = ms[ms.band == band].set_index("unit").reindex(order)
        y = np.arange(len(order))
        vals = sub["nlp"].to_numpy()
        sig = sub["matched_strength_p"].to_numpy() < 0.05
        pos = sub["M_obs"].to_numpy() > 0
        colors = [C_TRACE if (s and p) else C_ANTI if (s and not p) else C_NS
                  for s, p in zip(sig, pos)]
        # low-support units (K<4) are hatched + outlined so a fragile "survivor"
        # (e.g. beta occipital, K=3, single-patient driven) does not read as solid
        ksub = sub["K_implanted"].to_numpy()
        hatches = ["//" if (np.isfinite(k) and k < 4) else "" for k in ksub]
        bars = ax.barh(y, np.nan_to_num(vals), color=colors,
                       edgecolor="0.25", linewidth=0.5)
        for bar, ht in zip(bars, hatches):
            if ht:
                bar.set_hatch(ht)
        ax.axvline(P_REF, ls="--", lw=0.9, color="0.35")
        ax.set_yticks(y)
        ax.set_yticklabels(order)
        ax.set_xlabel(r"$-\log_{10} p_{\mathrm{MS}}$")
        btex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
        ax.text(0.97, 0.03, btex, transform=ax.transAxes, ha="right", va="bottom")
        ax.set_xlim(0, max(2.5, np.nanmax(ms["nlp"]) * 1.05))
    handles = [plt.Rectangle((0, 0), 1, 1, color=c)
               for c in (C_TRACE, C_ANTI, C_NS)]
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor="white",
                                 edgecolor="0.25", hatch="//"))
    fig.legend(handles, ["trace-direction, $p_{\\mathrm{MS}}<0.05$",
                         "anti-direction, $p_{\\mathrm{MS}}<0.05$",
                         "not significant", r"low support ($K<4$)"],
               loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    out = OUT / "fig_localization_matched_strength.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# Figure 2 — beta: two-null comparison + strength-residualization
# ---------------------------------------------------------------------------
def fig_beta_two_null():
    """beta: label-shuffle vs matched-strength significance per system.

    The crux figure: MTL (and others) clear the label-shuffle null but collapse
    under the strength-preserving matched-strength null; OFC clears both. This is
    why the strength control is decisive (the hippocampal/MTL concentration is
    real vs a naive shuffle but sub-threshold once node strength is held fixed).
    """
    ls = load_systems("atlas_cophenet_include.csv", "perm_p")
    ms = load_systems("matched_strength_include.csv", "matched_strength_p")
    lsb = ls[ls.band == "beta"].set_index("unit")
    msb = ms[ms.band == "beta"].set_index("unit")
    units = [u for u in msb.index if u in lsb.index]
    units = sorted(units, key=lambda u: msb.loc[u, "nlp"])      # by MS significance

    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    y = np.arange(len(units))
    h = 0.38
    lsv = [lsb.loc[u, "nlp"] for u in units]
    msv = [msb.loc[u, "nlp"] for u in units]
    ax.barh(y + h / 2, lsv, height=h, color="#9ecae1", edgecolor="0.3",
            linewidth=0.4, label="label-shuffle null (sampling)")
    ax.barh(y - h / 2, msv, height=h, color="#016450", edgecolor="0.3",
            linewidth=0.4, label="matched-strength null")
    ax.axvline(P_REF, ls="--", lw=0.9, color="0.35")
    ax.set_yticks(y)
    ax.set_yticklabels(units)
    ax.set_xlabel(r"$-\log_{10} p$")
    ax.set_xlim(0, max(3.1, max(lsv + msv) * 1.05))
    btex = BRAIN_BAND_TEX_DICT.get("beta", r"$\beta$")
    ax.text(0.97, 0.03, btex, transform=ax.transAxes, ha="right", va="bottom")
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    fig.tight_layout()
    out = OUT / "fig_localization_beta_two_null.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  wrote {out}")


def main():
    fig_survivor_map()
    fig_beta_two_null()
    print(f"Outputs -> {OUT}/")


if __name__ == "__main__":
    main()
