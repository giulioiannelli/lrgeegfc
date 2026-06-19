#!/usr/bin/env python3
"""audit_110b — the inference-mark localization ON THE BRAIN + per-system panel.

The figure for audit_110: within β, the offline trace splits by CONTENT —
the standard trace and the ENCODING component (memorised pairs) over-accumulate
in ORBITOFRONTAL CORTEX, but the INFERENCE-specific component (the figured-out
relations, controlling encoding) over-accumulates in the CINGULATE. A double
dissociation: the cingulate carries INFERENCE in β but ENCODING in α/low-γ.

FRAMING (LOCKED, feedback_localization_is_overexpression_not_container): every
panel shows the trace as DISTRIBUTED across the implant footprint (grey markers /
all 9 systems) with the hotspot as an ACCUMULATION glowing on top — never an
exclusive container. Colour is the SYSTEM effect (audit_76/108/109 lesson:
single-contact cophenetic is too noisy to map per-contact), so a contact glows by
how strongly its anatomical SYSTEM accumulates the signal.

Reuses audit_109's glass-brain machinery verbatim (pooled MNI coords, draw_row,
colorbar/legend) and audit_110's observed concordance (per-patient-normalised for
the system panel). NO recompute of the trace, NO tree-cutting.

Two PDFs (data/audit/inference_localization/figures/):
  brain_inference_dissociation_beta.pdf  3 rows (standard / encoding / inference)
                                         × 4-view glass brain, β: OFC → OFC →
                                         cingulate migration.
  systems_inference_dissociation.pdf     per-system normalised effect, β, the
                                         distributed field + the two hotspots;
                                         + cingulate role-flip across α/β/low-γ.

Usage:
    python audit_110b_inference_localization_figure.py
"""
from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
# reuse the glass-brain machinery + the observed concordance (import, no fork)
from audit_109_taxonomy_brain_figure import (  # type: ignore
    pooled_coords_systems, draw_row, row_label, add_colorbar,
    ROW_W, ROW_LEFT,
)
from audit_110_inference_mark_localization import obs_targets, COHORT  # type: ignore
from audit_83_localization_matched_strength import unit_means_from_s, DROP  # type: ignore

OUT = ROOT / "data" / "audit" / "inference_localization"
FIG = OUT / "figures"

# β rows read from the R=1000 CSV if present (sharper p), else R=200.
ROW_TARGETS = [("standard", "standard trace"),
               ("encoding", "encoding"),
               ("inference_pe", "inference (| encoding)")]
PANEL_TARGETS = ["standard", "encoding", "inference_pe"]
PANEL_COLORS = {"standard": "#7a7a7a", "encoding": "#1f6fb2",
                "inference_pe": "#d1491c"}
PANEL_LABEL = {"standard": "standard trace", "encoding": "encoding",
               "inference_pe": "inference ( | encoding)"}
# system display order: paralimbic/limbic core first, then the rest
SYS_ORDER = ["OFC", "cingulate", "MTL", "insula", "lateral_temporal", "PFC",
             "parietal", "sensorimotor", "occipital", "subcortical_other"]
BANDS_FLIP = ["alpha", "beta", "low_gamma"]
ROW_BOT3 = [0.70, 0.40, 0.10]
ROW_H3 = 0.27


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #
def load_inf_csv(band):
    """Per-system matched-strength rows for a band; β prefers the R=1000 CSV."""
    r1000 = OUT / "inference_mark_R1000_include.csv"
    base = OUT / "inference_mark_include.csv"
    if band == "beta" and r1000.exists():
        df = pd.read_csv(r1000)
        if not df[(df.band == "beta")].empty:
            return df, "R=1000"
    return pd.read_csv(base), "R=200"


def inf_system_effect(df, band, target, q_thresh=0.05):
    """system -> min(-log10 p, 3); + BH-significant positive set (for the ring)."""
    s = df[(df.band == band) & (df.target == target)
           & (df.granularity == "system")]
    val, sig = {}, set()
    for r in s.itertuples():
        p = max(float(r.matched_strength_p), 1e-4)
        val[r.unit] = min(-np.log10(p), 3.0)
        if float(r.bh_q_system) < q_thresh and (r.M_obs - r.surr_median) > 0:
            sig.add(r.unit)
    return val, sig


def per_patient_norm(band):
    """dict target -> {system -> list of per-patient RMS-normalised unit-means}.
    Observed only (no surrogate); each patient's 9-system vector is scaled to unit
    RMS so patients are comparable — the distributed field + hotspot, honestly."""
    out = {t: {} for t in PANEL_TARGETS}
    for pat in COHORT:
        try:
            t_obs, iu_i, iu_j = obs_targets(pat, band)
        except FileNotFoundError:
            continue
        uv = load_channel_regions(pat)["system"].to_numpy().astype(str)
        for t in PANEL_TARGETS:
            om = unit_means_from_s(t_obs[t], iu_i, iu_j, uv, None, DROP["system"])
            if not om:
                continue
            arr = np.array(list(om.values()), float)
            rms = np.sqrt(np.mean(arr ** 2)) or 1.0
            for s_, v in om.items():
                out[t].setdefault(s_, []).append(v / rms)
    return out


# --------------------------------------------------------------------------- #
# figure 1: brain dissociation (β, 3 rows)
# --------------------------------------------------------------------------- #
def fig_brain(coords, systems):
    use_lrg_style()
    df, rtag = load_inf_csv("beta")
    fig = plt.figure(figsize=(13.5, 9.8))
    for (target, label), bottom in zip(ROW_TARGETS, ROW_BOT3):
        val, sig = inf_system_effect(df, "beta", target)
        draw_row(fig, (ROW_LEFT, bottom, ROW_W, ROW_H3), coords, systems, val, sig)
        row_label(fig, bottom, label)
        top = sorted(val.items(), key=lambda kv: -kv[1])[:2]
        print(f"  beta {target:13s} ({rtag}) top: "
              + ", ".join(f"{u}(-log10p={v:.2f})" for u, v in top)
              + f"  ringed: {sorted(sig) or '(none)'}")
    add_colorbar(fig, "system accumulation  ($-\\log_{10}p_{\\mathrm{MS}}$)")
    h = [Line2D([], [], marker="o", ls="", mfc="none", mec="black", mew=1.2,
                ms=10, label="system clears BH q<0.05 across the 9 systems"),
         Line2D([], [], marker="o", ls="", mfc="0.72", mec="none", ms=7,
                label="not a hotspot (p≥0.05) / white-matter contact")]
    fig.legend(handles=h, loc="lower center", frameon=False, fontsize=9, ncol=2,
               bbox_to_anchor=(0.5, 0.002))
    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / "brain_inference_dissociation_beta.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


# --------------------------------------------------------------------------- #
# figure 2: per-system panel (β field+hotspots) + cingulate band-flip
# --------------------------------------------------------------------------- #
def _star(df, band, target, unit):
    s = df[(df.band == band) & (df.target == target)
           & (df.granularity == "system") & (df.unit == unit)]
    if s.empty:
        return ""
    p = float(s.iloc[0]["matched_strength_p"])
    return "*" if p < 0.05 else ""


def fig_systems(beta_norm):
    use_lrg_style()
    dfb, _ = load_inf_csv("beta")
    fig = plt.figure(figsize=(12.5, 5.6))
    axL = fig.add_axes([0.07, 0.16, 0.58, 0.74])
    axR = fig.add_axes([0.73, 0.16, 0.24, 0.74])

    # --- LEFT: β, all systems, three targets (distributed field + hotspots) ---
    syss = [s for s in SYS_ORDER if s in beta_norm["inference_pe"]]
    y = np.arange(len(syss))[::-1]
    off = {"standard": +0.24, "encoding": 0.0, "inference_pe": -0.24}
    for t in PANEL_TARGETS:
        med = [np.median(beta_norm[t].get(s, [np.nan])) for s in syss]
        axL.plot([0, 0], [y.min() - 0.5, y.max() + 0.5], color="0.8", lw=0.8, zorder=0)
        for yi, s_, m in zip(y, syss, med):
            axL.plot([0, m], [yi + off[t], yi + off[t]], color=PANEL_COLORS[t],
                     lw=1.6, alpha=0.85, zorder=1)
            star = _star(dfb, "beta", t, s_)
            axL.plot(m, yi + off[t], "o", color=PANEL_COLORS[t], ms=7,
                     mec="black" if star else "none", mew=1.1, zorder=2)
            if star:
                axL.annotate("∗", (m, yi + off[t]), textcoords="offset points",
                             xytext=(6, -1), fontsize=12, color=PANEL_COLORS[t])
    axL.set_yticks(y)
    axL.set_yticklabels(syss)
    axL.set_xlabel("per-patient–normalised cohort-median accumulation "
                   "(distributed field → 0; hotspot → +)")
    axL.set_ylim(y.min() - 0.6, y.max() + 0.6)
    handles = [Line2D([], [], color=PANEL_COLORS[t], marker="o", ls="-", ms=7,
                      label=PANEL_LABEL[t]) for t in PANEL_TARGETS]
    handles.append(Line2D([], [], color="black", marker="o", ls="", mfc="none",
                          mew=1.1, ms=8, label="matched-strength p<0.05"))
    fig.legend(handles=handles, loc="lower center", frameon=False, ncol=4,
               fontsize=9, bbox_to_anchor=(0.5, -0.02))

    # --- RIGHT: cingulate role-flip across bands (encoding vs inference) ---
    xb = np.arange(len(BANDS_FLIP))
    for t, mk in (("encoding", "encoding"), ("inference_pe", "inference")):
        med = []
        for band in BANDS_FLIP:
            nb = per_patient_norm(band) if band != "beta" else beta_norm
            med.append(np.median(nb[t].get("cingulate", [np.nan])))
        axR.plot(xb, med, "-o", color=PANEL_COLORS[t], ms=7, lw=1.8, label=mk)
    axR.axhline(0, color="0.8", lw=0.8)
    axR.set_xticks(xb)
    axR.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS_FLIP])
    axR.set_title("cingulate", fontsize=11)
    axR.set_ylabel("normalised accumulation")

    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / "systems_inference_dissociation.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def main():
    coords, systems = pooled_coords_systems()
    print(f"\npooled {len(coords)} contacts "
          f"({int((systems != 'non_anatomical').sum())} anatomical)\n")
    print("== brain_inference_dissociation_beta ==")
    fig_brain(coords, systems)
    print("\n== systems_inference_dissociation ==")
    beta_norm = per_patient_norm("beta")
    fig_systems(beta_norm)


if __name__ == "__main__":
    main()
