#!/usr/bin/env python3
"""β anatomy figure: distributed cortical localization at both LRG probes.

Two-row composite (six panels total).

Row 0 — D_coph per-pair multiscale probe (top-decile of |Δρ_split^coph|).
  (a) Glass-brain LR-sagittal + axial + coronal views with cohort sEEG
      contacts. Markers in any of the 7 locked D_coph regions
      (passes A1 hypergeometric AND A3 matched-strength, excluding the
      catch-all `Unk` bucket) are coloured per region; the remaining
      contacts are drawn faint grey for context.
  (b) 3D MNI scatter — significant-region contacts only, marker area
      proportional to enrichment ratio.
  (c) Horizontal bar chart of the 7 regions ranked by enrichment ratio,
      annotated with hypergeometric `q_BH` and matched-strength
      `p_empirical`.

Row 1 — Grassmann subspace probe aggregated over S(β).
  S(β) = { k : p_k(β) < 0.05 } = 40 cluster-extent-significant cells
  (span k=[21, 90], non-contiguous), the support of the cluster-mass
  statistic T_G^*. Per-node participation is the unweighted average
  ||U_k^T e_i||² across S(β); the retired contiguous-significant window
  K*(β)=[27,55] is not used (locked methodology, ANATOMY_LEDGER.md
  2026-05-19 pm).
  (d) Glass-brain view, markers in the 7 locked Grassmann regions
      (Hippocampus + 6 cortical, A3 only — A1 is structurally sparse on
      the per-node-participation endpoint of this probe).
  (e) 3D MNI scatter of the same Grassmann regions.
  (f) Horizontal bar chart annotated with matched-strength `p_empirical`.

Outputs (PDF only)
------------------
data/preprint/figures/beta/anatomy/fig_beta_anatomy_brain.pdf

Inputs
------
data/audit/anatomy_beta_cophenet/cohort_summary.csv
data/audit/anatomy_beta_grassmann_clusterext/cohort_summary.csv
data/raw/stereoeeg_patients/Pat_NN/implant_pat_NN.csv
data/raw/stereoeeg_patients/Pat_NN/channel_labels.csv
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.spatial import (
    load_spatial_metadata,
    prepare_spatial_coordinates,
)
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()


COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]

ANAT_COPH_CSV = ROOT / "data" / "audit" / "anatomy_beta_cophenet" / "cohort_summary.csv"
# Use the cluster-extent set (all p<0.05 k cells from grassmann_cluster_extent
# /per_k_obs_p.csv) rather than the longest-contiguous-run window. The C3
# Grassmann gate (CONTROLS.md) is the disjunctive cluster-extent permutation —
# anatomy is scoped to the same set so the localization reads the actual
# significant k-mass, not just the LR sub-window.
ANAT_GRASS_CSV = (ROOT / "data" / "audit"
                  / "anatomy_beta_grassmann_clusterext" / "cohort_summary.csv")

# Locked region selections (CSV column `region`).
# Filter excludes the catch-all "Unk" / non-DK buckets.
COPH_REGIONS = [
    "ctx-rh-rostralanteriorcingulate",
    "ctx-lh-isthmuscingulate",
    "ctx-rh-insula",
    "ctx-lh-parahippocampal",
    "ctx-lh-superiorfrontal",
    "ctx-rh-postcentral",
    "ctx-lh-entorhinal",
]
GRASS_REGIONS = [
    "Hip",
    "ctx-rh-rostralmiddlefrontal",
    "ctx-lh-lateralorbitofrontal",
    "ctx-rh-medialorbitofrontal",
    "ctx-lh-middletemporal",
    "ctx-lh-superiortemporal",
    "ctx-lh-insula",
]

# 14 unique categorical hues — every region gets its OWN colour across both
# probes. The two palettes never share a hex value.
PALETTE_COPH = [
    "#1f3d6e",  # navy           — R rostral ant cingulate
    "#0d5c1c",  # forest green   — L isthmus cingulate
    "#dc7633",  # burnt orange   — R insula
    "#7d3c98",  # royal purple   — L parahippocampal
    "#a04000",  # brown          — L superior frontal
    "#566573",  # slate          — R postcentral
    "#2e86c1",  # sky blue       — L entorhinal
]
PALETTE_GRASS = [
    "#cb4335",  # brick red      — Hippocampus
    "#f1c40f",  # golden yellow  — R rostral middle frontal
    "#e84393",  # magenta        — L lateral OFC
    "#16a085",  # jade           — R medial OFC
    "#5b2c6f",  # plum           — L middle temporal
    "#239b56",  # emerald        — L superior temporal
    "#85929e",  # cool gray      — L insula
]

# Pretty labels for the region rows
REGION_PRETTY = {
    "ctx-lh-isthmuscingulate": "L isthmus cingulate",
    "ctx-rh-rostralanteriorcingulate": "R rostral ant cingulate",
    "ctx-lh-parahippocampal": "L parahippocampal",
    "ctx-lh-entorhinal": "L entorhinal",
    "ctx-rh-insula": "R insula",
    "ctx-rh-postcentral": "R postcentral",
    "ctx-lh-superiorfrontal": "L superior frontal",
    "Hip": "Hippocampus",
    "ctx-lh-insula": "L insula",
    "ctx-lh-middletemporal": "L middle temporal",
    "ctx-lh-superiortemporal": "L superior temporal",
    "ctx-lh-lateralorbitofrontal": "L lateral OFC",
    "ctx-rh-medialorbitofrontal": "R medial OFC",
    "ctx-rh-rostralmiddlefrontal": "R rostral mid frontal",
}

NONSIG_COLOR = "#c8c8c8"


# ---------------------------------------------------------------------------
# Cohort contact pool
# ---------------------------------------------------------------------------
def _primary_region(dk_str) -> str | None:
    """Mirror of spatial._estimate_mni_transform's region parser."""
    if pd.isna(dk_str):
        return None
    parts = str(dk_str).strip().split(",")
    if not parts:
        return None
    return parts[0].strip().strip('"').strip()


def collect_cohort_contacts(cohort: list[str]) -> pd.DataFrame:
    rows = []
    for pat in cohort:
        meta = load_spatial_metadata(pat, SEEG_DATAPATH)
        meta = meta.dropna(subset=["x", "y", "z"]).reset_index(drop=True)
        coords = prepare_spatial_coordinates(meta, scale="mm",
                                             center=False, to_mni=True)
        regions = meta["Desikan-Killany"].apply(_primary_region).values
        for i, region in enumerate(regions):
            rows.append(dict(
                patient=pat,
                idx=i,
                x=coords[i, 0],
                y=coords[i, 1],
                z=coords[i, 2],
                region=region,
            ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Panels — glass brain
# ---------------------------------------------------------------------------
def panel_glass_brain(ax, contacts: pd.DataFrame, regions: list[str],
                      palette: list[str], display_mode: str = "lzry",
                      title: str | None = None,
                      marker_size_sig: float = 28.0,
                      marker_size_bg: float = 6.0) -> dict:
    from nilearn.plotting import plot_glass_brain

    disp = plot_glass_brain(
        None,
        figure=ax.figure,
        axes=ax,
        display_mode=display_mode,
        annotate=False,
        plot_abs=False,
        colorbar=False,
    )

    sig_mask = contacts["region"].isin(regions).values
    bg_coords = contacts.loc[~sig_mask, ["x", "y", "z"]].values
    if bg_coords.size > 0:
        disp.add_markers(
            marker_coords=bg_coords,
            marker_color=NONSIG_COLOR,
            marker_size=marker_size_bg,
            alpha=0.6,
        )

    sig_counts: dict[str, int] = {}
    for region, color in zip(regions, palette):
        sel = contacts.loc[contacts["region"] == region,
                           ["x", "y", "z"]].values
        sig_counts[region] = int(len(sel))
        if sel.size == 0:
            continue
        disp.add_markers(
            marker_coords=sel,
            marker_color=color,
            marker_size=marker_size_sig,
            alpha=0.95,
        )

    if title is not None:
        ax.set_title(title, loc="left", pad=4)
    return dict(display=disp, sig_counts=sig_counts)


# ---------------------------------------------------------------------------
# Panels — 3D scatter
# ---------------------------------------------------------------------------
def panel_3d_scatter(ax, contacts: pd.DataFrame, regions: list[str],
                     palette: list[str], enrichment_by_region: dict[str, float],
                     title: str | None = None,
                     marker_size_sig_range: tuple[float, float] = (35.0, 110.0),
                     marker_size_bg: float = 4.0) -> None:
    sig_mask = contacts["region"].isin(regions).values
    bg = contacts.loc[~sig_mask, ["x", "y", "z"]].values
    if bg.size > 0:
        ax.scatter(bg[:, 0], bg[:, 1], bg[:, 2],
                   s=marker_size_bg, c=NONSIG_COLOR,
                   alpha=0.30, depthshade=False, zorder=1)

    if enrichment_by_region:
        e_vals = list(enrichment_by_region.values())
        e_lo, e_hi = float(min(e_vals)), float(max(e_vals))
    else:
        e_lo, e_hi = 1.0, 2.0
    s_lo, s_hi = marker_size_sig_range

    def _size(e):
        if e_hi == e_lo:
            return 0.5 * (s_lo + s_hi)
        t = (e - e_lo) / (e_hi - e_lo)
        return s_lo + t * (s_hi - s_lo)

    for region, color in zip(regions, palette):
        sel = contacts.loc[contacts["region"] == region,
                           ["x", "y", "z"]].values
        if sel.size == 0:
            continue
        e = enrichment_by_region.get(region, 1.0)
        ax.scatter(sel[:, 0], sel[:, 1], sel[:, 2],
                   s=_size(e), c=color, alpha=0.95,
                   edgecolor="white", linewidth=0.6,
                   depthshade=True, zorder=3)

    # Tight bounds on the actual contact cloud (not a cubic shell, so the
    # axes box fills with data instead of empty cube margins).
    all_xyz = contacts[["x", "y", "z"]].values
    lo = all_xyz.min(axis=0)
    hi = all_xyz.max(axis=0)
    pad = 0.04 * (hi - lo)
    ax.set_xlim(lo[0] - pad[0], hi[0] + pad[0])
    ax.set_ylim(lo[1] - pad[1], hi[1] + pad[1])
    ax.set_zlim(lo[2] - pad[2], hi[2] + pad[2])
    span = hi - lo
    try:
        ax.set_box_aspect(span / span.max())
    except (AttributeError, TypeError):
        pass
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xlabel("x", fontsize=8, labelpad=-12)
    ax.set_ylabel("y", fontsize=8, labelpad=-12)
    ax.set_zlabel("z", fontsize=8, labelpad=-12)
    ax.view_init(elev=18, azim=-65)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_edgecolor("#dddddd")
        axis.line.set_color("#dddddd")
    if title is not None:
        ax.set_title(title, loc="left", pad=2, fontsize=10)


# ---------------------------------------------------------------------------
# Panels — region bar chart
# ---------------------------------------------------------------------------
def panel_region_bars(ax, df: pd.DataFrame, regions: list[str],
                      palette: list[str], counts: dict[str, int],
                      value_col: str = "enrichment",
                      annotate_cols: tuple[str, ...] = ("q_bh", "p_empirical"),
                      title: str | None = None,
                      vmax: float | None = None) -> None:
    sub = df.set_index("region").loc[regions].copy()
    order = sub[value_col].sort_values(ascending=False).index.tolist()
    sub = sub.loc[order]
    palette_ordered = [palette[regions.index(r)] for r in order]

    y = np.arange(len(sub))[::-1]
    vals = sub[value_col].values.astype(float)
    bars = ax.barh(y, vals, color=palette_ordered, edgecolor="0.10",
                   linewidth=0.8, alpha=0.95, zorder=3)

    if vmax is None:
        vmax = float(vals.max()) * 1.55

    for yi, r, v, clr in zip(y, sub.index, vals, palette_ordered):
        text = REGION_PRETTY.get(r, r)
        # Luminance-aware label colour so light bars get dark text and
        # dark bars get white text.
        rgb = matplotlib.colors.to_rgb(clr)
        lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
        tcolor = "white" if lum < 0.55 else "#1c1c1c"
        ax.text(0.04, yi, text, ha="left", va="center",
                fontsize=9, color=tcolor, fontweight="bold", zorder=5)

    for yi, r, v in zip(y, sub.index, vals):
        parts = []
        if "q_bh" in annotate_cols:
            q = float(sub.loc[r, "q_bh"])
            if q < 1.0:
                parts.append(rf"$q_{{\mathrm{{BH}}}}\!=\!{q:.0e}$" if q < 1e-3
                             else rf"$q_{{\mathrm{{BH}}}}\!=\!{q:.3f}$")
        if "p_empirical" in annotate_cols:
            p = float(sub.loc[r, "p_empirical"])
            parts.append(rf"$p_{{\mathrm{{emp}}}}\!=\!{p:.3f}$")
        parts.append(rf"$n\!=\!{counts.get(r, 0)}$")
        ax.text(v + vmax * 0.015, yi, "   ".join(parts),
                ha="left", va="center", fontsize=8, color="0.15", zorder=5)

    ax.axvline(1.0, color="0.55", lw=0.7, ls="--", zorder=2)
    ax.set_yticks([])
    ax.set_xlim(0.0, vmax)
    ax.set_xlabel("enrichment ratio")
    if title is not None:
        ax.set_title(title, loc="left", pad=4)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(left=False)


# ---------------------------------------------------------------------------
# Figure assembly
# ---------------------------------------------------------------------------
def build_figure() -> Path:
    out_dir = ROOT / "data" / "preprint" / "figures" / "beta" / "anatomy"
    out_dir.mkdir(parents=True, exist_ok=True)

    coph = pd.read_csv(ANAT_COPH_CSV)
    grass = pd.read_csv(ANAT_GRASS_CSV)
    contacts = collect_cohort_contacts(COHORT)

    fig = plt.figure(figsize=(13.0, 15.0))
    gs = GridSpec(4, 2, figure=fig,
                  height_ratios=[0.55, 1.0, 0.55, 1.0],
                  width_ratios=[1.0, 1.25],
                  hspace=0.18, wspace=0.06,
                  left=0.02, right=0.99,
                  top=0.975, bottom=0.045)

    ax_brain_dc = fig.add_subplot(gs[0, :])
    ax_3d_dc = fig.add_subplot(gs[1, 0], projection="3d")
    ax_bar_dc = fig.add_subplot(gs[1, 1])
    ax_brain_gr = fig.add_subplot(gs[2, :])
    ax_3d_gr = fig.add_subplot(gs[3, 0], projection="3d")
    ax_bar_gr = fig.add_subplot(gs[3, 1])

    coph_enrichment = (coph.set_index("region")
                       .loc[COPH_REGIONS, "enrichment"].to_dict())
    grass_enrichment = (grass.set_index("region")
                        .loc[GRASS_REGIONS, "enrichment"].to_dict())

    dc_info = panel_glass_brain(
        ax_brain_dc, contacts, COPH_REGIONS, PALETTE_COPH,
        title=r"(a) $D_{\mathrm{coph}}$ glass-brain (LR sagittal + axial + coronal)")
    panel_3d_scatter(
        ax_3d_dc, contacts, COPH_REGIONS, PALETTE_COPH, coph_enrichment,
        title=r"(b) 3D MNI scatter — marker $\propto$ enrichment")
    panel_region_bars(
        ax_bar_dc, coph, COPH_REGIONS, PALETTE_COPH, dc_info["sig_counts"],
        value_col="enrichment",
        annotate_cols=("q_bh", "p_empirical"),
        title=r"(c) $D_{\mathrm{coph}}$ region enrichment "
              r"(hypergeometric $\times$ matched-strength)",
        vmax=4.2,
    )

    gr_info = panel_glass_brain(
        ax_brain_gr, contacts, GRASS_REGIONS, PALETTE_GRASS,
        title=r"(d) Grassmann glass-brain — cluster-extent set "
              r"(40 cells, $p_k\!<\!0.05$)")
    panel_3d_scatter(
        ax_3d_gr, contacts, GRASS_REGIONS, PALETTE_GRASS, grass_enrichment,
        title=r"(e) 3D MNI scatter — Hippocampus + 6 cortical")
    panel_region_bars(
        ax_bar_gr, grass, GRASS_REGIONS, PALETTE_GRASS, gr_info["sig_counts"],
        value_col="enrichment",
        annotate_cols=("p_empirical",),
        title=r"(f) Grassmann region enrichment "
              r"(matched-strength only)",
        vmax=4.2,
    )

    leg_handles = [
        mlines.Line2D([], [], marker="o", linestyle="None",
                      markersize=10, markerfacecolor=NONSIG_COLOR,
                      markeredgecolor="white",
                      label="cohort contacts (non-trace regions)"),
        mpatches.Patch(facecolor="white", edgecolor="0.20",
                       label=r"medial-temporal overlap: parahippocampal "
                             r"($D_{\mathrm{coph}}$) $\leftrightarrow$ "
                             r"Hippocampus (Grassmann)"),
        mpatches.Patch(facecolor="white", edgecolor="0.20",
                       label=r"insula overlap: R "
                             r"($D_{\mathrm{coph}}$) $\leftrightarrow$ "
                             r"L (Grassmann) — distinct DK labels"),
    ]
    fig.legend(handles=leg_handles, loc="lower center", ncol=3,
               frameon=False, bbox_to_anchor=(0.5, 0.01),
               fontsize=9)

    out = out_dir / "fig_beta_anatomy_brain.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    build_figure()
