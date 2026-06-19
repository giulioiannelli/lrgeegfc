#!/usr/bin/env python3
"""audit_109 — the cross-phase taxonomy ON THE BRAIN (glass-brain dissociation).

The anatomical trace figure (the catchy figure the taxonomy was really for).
Pools all 10 patients' contacts into approximate MNI space and lights up each
anatomical system by how strongly it carries a cross-phase channel, straight
from the verified geometry-baseline localization (audit_107). Shows the beta
DISSOCIATION: trace -> orbitofrontal, reset -> lateral temporal, anchor ->
insula — anatomically distinct systems.

Encoding (per contact):
  colour  = -log10(geometry-baseline p) for the contact's anatomical SYSTEM —
            the enrichment of channel signal on the EDGES INCIDENT to that system
            (audit_107 edge->endpoint incidence). NOT a per-node classification:
            the taxonomy is per-PAIR, and a system/node has three independent
            incidence scores (trace / reset / anchor), not one label. Bright =
            incident edges carry the channel; white-matter / unknown = grey.
  ring    = the system passed BH-FDR within (band, channel), q < 0.05.

This is a result-at-system-granularity map (the granularity audit_107 actually
verified), NOT a per-contact field — single-contact cophenetic values are too
noisy to read (audit_76/108 lesson). Colour/ring come straight from
localization_include.csv; the only thing computed here is the cohort MNI
position of each contact (approximate per-patient DK-anchored affine in
visuals.spatial_coords — adequate for a cohort glass brain, not mm-precise).
NO tree-cutting, NO recompute of the trace.

Two figures (data/audit/cross_phase_taxonomy/figures/):
  brain_dissociation_beta.pdf  3 rows (trace / reset / anchor) x 4-view glass
                               brain, beta.
  brain_trace_bands.pdf        trace channel across alpha / beta / low_gamma —
                               beta lights OFC; alpha/low_gamma trace is
                               non-significant => beta-specificity.

Honest scope: trace -> OFC is the verified/locked result; reset and anchor are
new and UNVERIFIED (geometry baseline only).

Usage:
    python audit_109_taxonomy_brain_figure.py
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from nilearn.plotting import plot_glass_brain

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.spatial_coords import (
    load_spatial_metadata,
    prepare_spatial_coordinates,
)
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
OUT = ROOT / "data" / "audit" / "cross_phase_taxonomy"
FIG = OUT / "figures"
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]

# (row label, localization-CSV signal); the beta dissociation
CHANNELS = [("trace", "trace"), ("reset", "reset"), ("anchor", "anchor")]
TRACE_BANDS = ["alpha", "beta", "low_gamma"]      # bands audit_107 localized
DISPLAY_MODE = "lyrz"
CMAP = matplotlib.colormaps["plasma"]
C_MIN = 1.301                                      # -log10 p >= 1.301 (p<0.05) glows
NORM = Normalize(vmin=C_MIN, vmax=3.0)             # -log10 p, capped at 3
ROW_H, ROW_W, ROW_LEFT = 0.30, 0.83, 0.05
ROW_BOT = [0.68, 0.35, 0.02]                       # top, middle, bottom row


# --------------------------------------------------------------------------- #
# cohort contacts in MNI space (band-independent)
# --------------------------------------------------------------------------- #
def pooled_coords_systems():
    coords_all, systems_all = [], []
    for pat in COHORT:
        md = load_spatial_metadata(pat, SEEG_DATAPATH)
        rdf = load_channel_regions(pat)
        if len(md) != len(rdf):
            print(f"  [skip] {pat}: len md={len(md)} rdf={len(rdf)}")
            continue
        xyz = md[["x", "y", "z"]].to_numpy(float)
        finite = np.all(np.isfinite(xyz), axis=1)
        if finite.sum() == 0:
            print(f"  [skip] {pat}: no finite coordinates")
            continue
        md_f = md.loc[finite].reset_index(drop=True)
        coords = prepare_spatial_coordinates(md_f, scale="mm", center=False,
                                             to_mni=True)
        coords_all.append(np.asarray(coords, float))
        systems_all.append(rdf["system"].to_numpy().astype(str)[finite])
    return np.vstack(coords_all), np.concatenate(systems_all)


# --------------------------------------------------------------------------- #
# per-system effect from audit_107 (geometry baseline)
# --------------------------------------------------------------------------- #
def load_loc():
    loc = pd.read_csv(OUT / "localization_include.csv")
    return loc[loc.granularity == "system"]


def system_effect(loc, band, signal, q_thresh=0.05):
    """value dict system -> min(-log10 p, 3); plus the BH-significant enriched set."""
    s = loc[(loc.band == band) & (loc.signal == signal)]
    val, sig = {}, set()
    for r in s.itertuples():
        p = max(float(r.p_geom), 1e-4)
        val[r.unit] = min(-np.log10(p), 3.0)
        if r.q_geom < q_thresh and (r.M_obs - r.null_median) > 0:
            sig.add(r.unit)
    return val, sig


# --------------------------------------------------------------------------- #
# one glass-brain row
# --------------------------------------------------------------------------- #
def draw_row(fig, rect, coords, systems, valdict, sig_units):
    disp = plot_glass_brain(None, display_mode=DISPLAY_MODE, axes=rect,
                            figure=fig)
    c = np.array([valdict.get(s, 0.0) for s in systems])
    glow = c >= C_MIN
    # everything not enriched (white-matter, unknown, non-significant system):
    # faint grey context so the implant footprint is visible
    if (~glow).any():
        disp.add_markers(coords[~glow], marker_color="0.72", marker_size=6,
                         alpha=0.28)
    # enriched systems glow on a saturated ramp (strongest = bright yellow)
    if glow.any():
        disp.add_markers(coords[glow], marker_color=CMAP(NORM(c[glow])),
                         marker_size=26, alpha=0.97)
    # BH-significant systems: black ring
    if sig_units:
        m = np.isin(systems, list(sig_units))
        if m.any():
            disp.add_markers(coords[m], marker_color="none", marker_size=66,
                             edgecolors="black", linewidths=1.2)
    return disp


def row_label(fig, bottom, text):
    fig.text(0.018, bottom + ROW_H / 2, text, rotation=90, va="center",
             ha="center", fontsize=13, fontweight="bold")


def add_colorbar(fig, label):
    cax = fig.add_axes([0.90, 0.36, 0.014, 0.30])
    cb = fig.colorbar(ScalarMappable(norm=NORM, cmap=CMAP), cax=cax)
    cb.set_label(label, fontsize=9)
    cb.set_ticks([1.301, 2, 3])
    cb.set_ticklabels(["1.3\n(p.05)", "2", "≥3"])


def add_legend(fig):
    h = [Line2D([], [], marker="o", ls="", mfc="none", mec="black", mew=1.2,
                ms=10, label="system passes BH q<0.05 (geometry baseline)"),
         Line2D([], [], marker="o", ls="", mfc="0.72", mec="none", ms=7,
                label="non-enriched (p≥0.05) / white-matter contact")]
    fig.legend(handles=h, loc="lower center", frameon=False, fontsize=9, ncol=2,
               bbox_to_anchor=(0.5, 0.002))


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def fig_dissociation(coords, systems, loc):
    use_lrg_style()
    fig = plt.figure(figsize=(13.5, 9.8))
    for (name, signal), bottom in zip(CHANNELS, ROW_BOT):
        val, sig = system_effect(loc, "beta", signal)
        draw_row(fig, (ROW_LEFT, bottom, ROW_W, ROW_H), coords, systems, val, sig)
        row_label(fig, bottom, name)
        qd = {u: float(loc[(loc.band == "beta") & (loc.signal == signal)
                           & (loc.unit == u)].q_geom.iloc[0]) for u in sig}
        print(f"  beta {name:7s} ringed: "
              f"{', '.join(f'{u}(q={q:.02g})' for u, q in sorted(qd.items())) or '(none)'}")
    add_colorbar(fig, "incident-edge enrichment  ($-\\log_{10}p$)")
    add_legend(fig)
    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / "brain_dissociation_beta.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def fig_trace_bands(coords, systems, loc):
    use_lrg_style()
    fig = plt.figure(figsize=(13.5, 9.8))
    for band, bottom in zip(TRACE_BANDS, ROW_BOT):
        val, sig = system_effect(loc, band, "trace")
        draw_row(fig, (ROW_LEFT, bottom, ROW_W, ROW_H), coords, systems, val, sig)
        row_label(fig, bottom, BRAIN_BAND_TEX_DICT.get(band, band))
        print(f"  trace {band:10s} ringed: {', '.join(sorted(sig)) or '(none)'}")
    add_colorbar(fig, "trace-edge enrichment  ($-\\log_{10}p$)")
    add_legend(fig)
    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / "brain_trace_bands.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def main():
    loc = load_loc()
    coords, systems = pooled_coords_systems()
    print(f"\npooled {len(coords)} contacts across {len(COHORT)} patients "
          f"({int((systems != 'non_anatomical').sum())} anatomical)")
    print("\n== brain_dissociation_beta ==")
    fig_dissociation(coords, systems, loc)
    print("\n== brain_trace_bands ==")
    fig_trace_bands(coords, systems, loc)


if __name__ == "__main__":
    main()
