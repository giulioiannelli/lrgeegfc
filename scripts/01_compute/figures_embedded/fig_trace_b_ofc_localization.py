#!/usr/bin/env python3
r"""fig:trace_b — where the beta trace lives, on the brain (Results §1, para b).

CORE MESSAGE: pooling all ten implants into one glass brain and lighting each contact by
how its anatomical system carries the beta cophenetic trace (relative to each patient's
own whole-brain baseline), the trace has a clear address. It CONCENTRATES ventrally, in
orbitofrontal cortex (OFC, green, ringed) — and is DEPLETED dorsally and posteriorly, in
prefrontal and sensorimotor cortex (orange, ringed). "Enriched" = the system holds more of
the held trace than the whole-brain baseline; "depleted" = less. Rings mark the three
systems whose direction survives the full robustness battery (R=1000 + leave-one-out +
shaft-collapse). OFC is the ventral-frontal green cluster that emerges from the stack.

HONEST NUANCE (tested, not hidden): occipital, cingulate and lateral-temporal also clear
the R=1000 enrichment tail — but when each is put through the IDENTICAL leave-one-out +
shaft battery (audit_162), their enrichment collapses on a single-patient drop (worst-drop
BH q = 0.54 occipital, 0.17 cingulate, 0.999 lateral-temporal), whereas OFC stays stable
(worst-drop q = 0.06, always the top system). So they are single-patient-driven, not
robustness-untested; occipital additionally rests on only three sparse nodes. Only OFC /
PFC / sensorimotor survive as cohort-robust and are ringed.

Colour = system direction (green enriched / orange depleted, BH q<0.05; grey = n.s. or
white matter). Saturated + ring = cohort-robust (R=1000 + LOO + shaft, audit_155/162);
pale = clears the gate but single-patient-driven under LOO.

Reads : data/audit/localization_atlas_rhosym/matched_strength_rhosym_include.csv
        data/audit/localization_atlas_rhosym/beta_all_systems_robustness_R1000.csv
Writes: data/preprint/figures/results_section1/fig_trace_b_ofc_localization.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from nilearn.plotting import plot_glass_brain

from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.spatial_coords import (
    load_spatial_metadata,
    prepare_spatial_coordinates,
)
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

SRC = ROOT / "data/audit/localization_atlas_rhosym/matched_strength_rhosym_include.csv"
OUT = ROOT / "data/preprint/figures/results_section1/fig_trace_b_ofc_localization.pdf"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
DROP = {"other", "non_anatomical"}
ROBUST = {"OFC", "PFC", "sensorimotor"}       # survive R=1000 + LOO + shaft (audit_155)
DISPLAY_MODE = "lzr"                            # left-sagittal, axial, right-sagittal
R = 200
C_ENR, C_DEP, C_NS = "#2a9d5c", "#d1791f", "#9a9a9a"
C_ENR_PALE, C_DEP_PALE = "#88c9a2", "#e6bd8f"   # not-robustness-tested (secondary)


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
            continue
        md_f = md.loc[finite].reset_index(drop=True)
        coords = prepare_spatial_coordinates(md_f, scale="mm", center=False, to_mni=True)
        coords_all.append(np.asarray(coords, float))
        systems_all.append(rdf["system"].to_numpy().astype(str)[finite])
    return np.vstack(coords_all), np.concatenate(systems_all)


def system_direction():
    """system -> ('enriched'|'depleted'|'ns'); BH-FDR over the 9 a-priori systems."""
    d = pd.read_csv(SRC)
    d = d[(d.band == "beta") & (d.granularity == "system") & (~d.unit.isin(DROP))].copy()
    enr = d.M_obs > d.surr_median
    p_dir = np.where(enr, d.matched_strength_p, d.matched_strength_p_lower)
    p_dir = np.clip(p_dir, 1.0 / (R + 1), 1.0)
    q = np.asarray(bh_fdr(p_dir.tolist()))
    out = {}
    for unit, is_enr, qq in zip(d.unit, enr.values, q):
        out[unit] = ("enriched" if is_enr else "depleted") if qq < 0.05 else "ns"
    return out


def legend_handles(compact=False):
    """The 5-tier system-direction legend, built fresh so the standalone AND the flat
    compound can both place it. ``compact`` shortens the labels for the tight fig2 tile."""
    if compact:
        labels = ["OFC (enriched)", "PFC / SM (depleted)", "enriched, 1-patient",
                  "depleted, 1-patient", "n.s. / white matter"]
    else:
        labels = ["OFC — enriched, robust", "PFC / sensorimotor — depleted, robust",
                  "enriched, single-patient-driven", "depleted, single-patient-driven",
                  "n.s. / white matter"]
    spec = [(C_ENR, 12), (C_DEP, 12), (C_ENR_PALE, 9), (C_DEP_PALE, 9), (C_NS, 7)]
    return [Line2D([], [], marker="o", ls="", mfc=c, mec="white", ms=m, label=lab)
            for (c, m), lab in zip(spec, labels)]


def draw_brain(target, rect=(0.02, 0.14, 0.96, 0.84)):
    """Pooled glass brain + system-direction markers into ``target`` at ``rect`` (figure
    coords). ``target`` MUST be a real Figure: nilearn resolves ``figure=`` to the root and
    lays the brain in ROOT coordinates, so a SubFigure's frame is ignored (brains would
    stack). The flat compound therefore passes the root Figure and a controlled rect."""
    direction = system_direction()
    coords, systems = pooled_coords_systems()

    robust = np.isin(systems, list(ROBUST))
    d_of = np.array([direction.get(s, "ns") for s in systems])
    # four foreground tiers: robust emerge (saturated+large); non-robust sig are pale
    m_enr_R = (d_of == "enriched") & robust
    m_dep_R = (d_of == "depleted") & robust
    m_enr_p = (d_of == "enriched") & ~robust
    m_dep_p = (d_of == "depleted") & ~robust
    m_bg = d_of == "ns"

    disp = plot_glass_brain(None, display_mode=DISPLAY_MODE, axes=rect, figure=target)
    # draw order: faint footprint -> pale non-robust -> saturated robust on top
    if m_bg.any():
        disp.add_markers(coords[m_bg], marker_color=C_NS, marker_size=6, alpha=0.22)
    if m_enr_p.any():
        disp.add_markers(coords[m_enr_p], marker_color=C_ENR_PALE, marker_size=16, alpha=0.75)
    if m_dep_p.any():
        disp.add_markers(coords[m_dep_p], marker_color=C_DEP_PALE, marker_size=16, alpha=0.75)
    if m_dep_R.any():
        disp.add_markers(coords[m_dep_R], marker_color=C_DEP, marker_size=44, alpha=0.98)
    if m_enr_R.any():
        disp.add_markers(coords[m_enr_R], marker_color=C_ENR, marker_size=44, alpha=0.98)

    print(f"fig:trace_b — beta trace on the brain "
          f"({len(coords)} contacts, {int((systems!='non_anatomical').sum())} anatomical)\n")
    for s in sorted(set(systems)):
        if direction.get(s) in ("enriched", "depleted"):
            n = int((systems == s).sum())
            print(f"  {s:16s} {direction[s]:9s} n={n:3d} {'RING' if s in ROBUST else ''}")
    print("\nTESTED (audit_162): occipital/cingulate/lateral_temporal clear the R=1000 tail "
          "but collapse under LOO (worst-drop q 0.54/0.17/0.999 = single-patient-driven); "
          "OFC stable (worst-drop 0.06). Only OFC/PFC/sensorimotor ringed as cohort-robust.")
    return disp


def draw(target):
    draw_brain(target)
    target.text(0.012, 0.96, r"$\mathbf{a}$", fontsize=17, va="top", ha="left",
                fontweight="bold")
    target.legend(handles=legend_handles(), loc="lower center", bbox_to_anchor=(0.5, -0.03),
                  ncol=3, frameon=False, fontsize=9.4, handletextpad=0.4,
                  columnspacing=1.2)


def main():
    fig = plt.figure(figsize=(10.2, 3.9))
    draw(fig)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
