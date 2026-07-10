#!/usr/bin/env python3
r"""fig:trace_g (hero) — the held trace crossing the gate: theta -> alpha -> beta (3D).

Three identical pooled-cohort brains inside a translucent pial shell, ONE per band,
left-to-right in ascending gate order: theta (null, gate p=0.78) -> alpha (thin,
p=0.024) -> beta (full, p=0.032). Every sphere is one sEEG contact, sized by its own
trace magnitude |rho_sym|: GREEN = carrier (its cophenetic task->rest move beats its OWN
strength-matched null; every carrier drawn identically), RED = reset/anti (the move
reverses), faint grey = the sea (no held trace). BLUE arcs are the strongest concordant
carrier<->carrier co-movements per patient -- the top FRAC_CC (7%) of that patient's
significant, concordant (dD_task*dD_rest > 0) carrier-carrier pairs by concordance-
product magnitude, so the arc count scales with the number of carriers (~carrier^2): few
in theta, many in beta. Read left to right: the green carrier field + its blue backbone
MATERIALISE as the band crosses the matched-strength gate, while theta stays a balanced
green/red speckle (carriers cancelled by reversers, mean trace strength ~0). Anatomical
stubs (FRONT/BACK/TOP/BOTTOM) orient a reader who cannot rotate the view.

WHY A BAND CONTRAST. Raw per-pair concordance does NOT separate the bands (theta 53% vs
beta 61% concordant pairs); the selectivity lives ONLY in the per-node matched-strength
comparison. Node colour/size and the backbone are all driven by rho_sym-vs-null
(per_node.csv polarity); theta renders flat because its nodes genuinely fail their own
strength-matched null, not because of styling.

ILLUSTRATIVE ARCS, NODE-LEVEL CLAIM. Per-pair arcs are noisy; the load-bearing evidence
is the per-node carrier/reset field (and the cohort gate + beta->OFC, Fig. 2), never a
single arc. MNI is a heuristic affine (schematic).

Reads : data/audit/per_node_trace_decomposition_rhosym/per_node.csv  (all 6 bands)
        data/reports/imcoh_continuous_trace/per_pair_split/{Pat}_{band}.npz
        fsaverage5 pial mesh (nilearn, cached)
Writes: data/preprint/figures/results_section1/fig_trace_g_pairglow_3d.{pdf,png,html}
"""
from __future__ import annotations

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.visuals.brain3d import pial_mesh, spheres_mesh, bezier_arcs
from lrg_eegfc.visuals.spatial_coords import (
    load_spatial_metadata, prepare_spatial_coordinates,
)

PN = ROOT / "data/audit/per_node_trace_decomposition_rhosym/per_node.csv"
PP = ROOT / "data/reports/imcoh_continuous_trace/per_pair_split"
OUT = ROOT / "data/preprint/figures/results_section1/fig_trace_g_pairglow_3d"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["theta", "alpha", "beta"]                       # ascending gate order
BAND_TEX = {"theta": "θ (4–8 Hz)", "alpha": "α (8–13 Hz)", "beta": "β (13–30 Hz)"}
BAND_GATE = {"theta": "gate p = 0.78 · null",
             "alpha": "gate p = 0.024 · thin",
             "beta": "gate p = 0.032 · full"}
FRAC_CC = 0.07                                           # fraction of avail. carrier<->carrier arcs (more lines)
CAP_CC = 46
C_TRACE, C_RESET, C_SEA, C_LINK = "#12b06a", "#e0342b", "#9aa0a6", "#2e6be6"  # green / red / grey / BLUE arcs
C_PIAL = "#868d97"                                        # darker pial shell (was near-white)
_NODE_LIGHT = dict(ambient=0.82, diffuse=0.40, specular=0.06, roughness=0.5)  # soft luminous, no halo
DOT_SCALE = 1.6                                           # global node-size multiplier (bigger contacts)
LINK_W_GLOW, LINK_W_CORE = 5.6, 2.4                       # blue co-movement arc widths (glow / core)
NODE_EDGE = "#111111"                                     # thin dark rim so contacts read against the shell

# MNI/RAS poles (config convention: +y FRONT, +z TOP). Only the sagittal-plane axes
# (anterior-posterior, superior-inferior) are labelled — left-right is the depth axis in
# this near-lateral view and would clip across panels.
AXIS_POLES = [(1, +1, "FRONT", True), (1, -1, "BACK", True),
              (2, +1, "TOP", True), (2, -1, "BOTTOM", False)]

pn_all = pd.read_csv(PN)
RHO95 = float(np.percentile(np.abs(pn_all[pn_all.band.isin(BANDS)].rho_sym), 95))


def _marker_px(rho, pol):
    """Screen-space marker diameter (px) per contact ~ its own |rho_sym|; carriers largest."""
    f = np.clip(np.abs(rho) / RHO95, 0, 1) ** 1.3
    base = np.where(pol == "carrier", 6.0 + 7.0 * f,
                    np.where(pol == "anti", 5.5 + 4.5 * f, 3.4))
    return base * DOT_SCALE


def build_panel(band):
    """Pool all ten implants for ``band``; return node arrays + green tracer edges."""
    sub = pn_all[pn_all.band == band]
    coords_l, pol_l, rho_l, edges = [], [], [], []
    off = 0
    for pat in COHORT:
        s = sub[sub.patient == pat].sort_values("node_idx")
        md = load_spatial_metadata(pat, SEEG_DATAPATH)
        if len(md) != len(s):
            continue
        xyz = md[["x", "y", "z"]].to_numpy(float)
        fin = np.all(np.isfinite(xyz), axis=1)
        coords = np.asarray(prepare_spatial_coordinates(
            md.loc[fin].reset_index(drop=True), scale="mm", center=False, to_mni=True), float)
        pol = s.polarity.to_numpy().astype(str)
        rho = s.rho_sym.to_numpy(float)
        local = np.full(len(s), -1, int); local[fin] = np.arange(fin.sum())
        coords_l.append(coords); pol_l.append(pol[fin]); rho_l.append(rho[fin])

        d = np.load(PP / f"{pat}_{band}.npz")
        i, j = d["iu_i"], d["iu_j"]
        prod = d["dD_task"] * d["dD_rest"]
        car = pol == "carrier"
        keep = (d["sigma"] > 0) & fin[i] & fin[j] & car[i] & car[j] & (prod > 0)
        ix = np.where(keep)[0]
        if ix.size:                              # fraction of avail. carrier-carrier pairs (~carrier^2)
            k = min(CAP_CC, max(3, int(round(FRAC_CC * ix.size))))
            top = ix[np.argsort(prod[ix])[-k:]]
            for a, b in zip(i[top], j[top]):
                edges.append((off + local[a], off + local[b]))
        off += fin.sum()
    return dict(coords=np.vstack(coords_l), pol=np.concatenate(pol_l),
                rho=np.concatenate(rho_l), edges=edges)


def panel_traces(P):
    """Ghost sea, blue backbone, red resets, green carriers — nodes as thin-edged markers."""
    c, pol, rho = P["coords"], P["pol"], P["rho"]
    px = _marker_px(rho, pol)
    tr = []
    xs, ys, zs = bezier_arcs(c, P["edges"], CENTER)                     # blue links, UNDER the nodes
    tr.append(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", hoverinfo="skip",
                           showlegend=False, opacity=0.09,
                           line=dict(color=C_LINK, width=LINK_W_GLOW)))  # link glow
    tr.append(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", hoverinfo="skip",
                           showlegend=False, opacity=0.78,
                           line=dict(color=C_LINK, width=LINK_W_CORE)))  # link core

    def _nodes(mask, color, opacity, edge_w, name):
        return go.Scatter3d(x=c[mask, 0], y=c[mask, 1], z=c[mask, 2], mode="markers",
                            hoverinfo="skip", showlegend=False, opacity=opacity, name=name,
                            marker=dict(size=px[mask], color=color,
                                        line=dict(color=NODE_EDGE, width=edge_w)))
    tr.append(_nodes(pol == "neutral", C_SEA, 0.32, 0.0, "sea"))        # faint sea, no rim
    tr.append(_nodes(pol == "anti", C_RESET, 0.97, 0.8, "reset"))       # thin dark rim
    tr.append(_nodes(pol == "carrier", C_TRACE, 1.0, 0.8, "carrier"))
    return tr


def orientation_traces(lo, hi, poles=AXIS_POLES):
    """Short outward stubs + anatomical labels so a static view is orientable."""
    c = 0.5 * (lo + hi)
    gap, ext = 7.0, 17.0
    sx, sy, sz = [], [], []
    big_pos, big_txt, sml_pos, sml_txt = [], [], [], []
    for ax, sgn, txt, big in poles:
        edge = hi[ax] if sgn > 0 else lo[ax]
        p0 = c.copy(); p0[ax] = edge + sgn * gap
        p1 = c.copy(); p1[ax] = edge + sgn * (gap + ext)
        sx += [p0[0], p1[0], None]; sy += [p0[1], p1[1], None]; sz += [p0[2], p1[2], None]
        lp = c.copy(); lp[ax] = edge + sgn * (gap + ext + 8)
        (big_pos if big else sml_pos).append(lp)
        (big_txt if big else sml_txt).append(txt)
    out = [go.Scatter3d(x=sx, y=sy, z=sz, mode="lines", hoverinfo="skip",
                        showlegend=False, line=dict(color="#7b828b", width=3))]
    for pos, txt, size, col in [(big_pos, big_txt, 12, "#2a2e34"),
                                (sml_pos, sml_txt, 9.5, "#9096a0")]:
        if pos:
            pos = np.array(pos)
            out.append(go.Scatter3d(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2], mode="text",
                                    text=txt, hoverinfo="skip", showlegend=False,
                                    textfont=dict(size=size, color=col, family="Arial Black")))
    return out


CENTER = np.array([0.0, -18.0, 12.0])                    # ~brain centroid; arcs bow away


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)

    fig = make_subplots(rows=1, cols=len(BANDS), horizontal_spacing=0.028,
                        specs=[[{"type": "scene"} for _ in BANDS]])
    counts = {}
    bbox = None
    for k, band in enumerate(BANDS, start=1):
        P = build_panel(band)
        if bbox is None:
            bbox = (P["coords"].min(0) - 2, P["coords"].max(0) + 2)
        counts[band] = (int((P["pol"] == "carrier").sum()),
                        int((P["pol"] == "anti").sum()),
                        int((P["pol"] == "neutral").sum()), len(P["edges"]))
        fig.add_trace(pial_mesh(color=C_PIAL, opacity=0.11), row=1, col=k)
        for t in panel_traces(P):
            fig.add_trace(t, row=1, col=k)
        for t in orientation_traces(*bbox):
            fig.add_trace(t, row=1, col=k)

    # one shared legend (dummy traces on the first scene)
    leg = [("carrier — holds trace", C_TRACE, "m"),
           ("reset — reverses", C_RESET, "m"),
           ("sea — no trace", C_SEA, "m"),
           ("strongest co-movements", C_LINK, "l")]
    for name, color, kind in leg:
        common = dict(x=[None], y=[None], z=[None], name=name, showlegend=True)
        fig.add_trace(go.Scatter3d(mode="markers", marker=dict(size=16, color=color, line=dict(color=NODE_EDGE, width=1.0)), **common)
                      if kind == "m" else
                      go.Scatter3d(mode="lines", line=dict(color=color, width=9), **common),
                      row=1, col=1)

    cam = dict(eye=dict(x=0.80, y=-0.76, z=0.44), up=dict(x=0, y=0, z=1))   # zoomed to fill each band
    fig.update_scenes(xaxis=dict(visible=False), yaxis=dict(visible=False),
                      zaxis=dict(visible=False), aspectmode="data", camera=cam,
                      bgcolor="rgba(0,0,0,0)")

    xc = [(2 * k + 1) / (2 * len(BANDS)) for k in range(len(BANDS))]  # domain centres
    for band, x in zip(BANDS, xc):
        fig.add_annotation(xref="paper", yref="paper", x=x, y=0.995, xanchor="center",
                           yanchor="top", showarrow=False, text=f"<b>{BAND_TEX[band]}</b>",
                           font=dict(size=18, color="#222"))
        fig.add_annotation(xref="paper", yref="paper", x=x, y=0.95, xanchor="center",
                           yanchor="top", showarrow=False, text=BAND_GATE[band],
                           font=dict(size=12.5, color="#777"))
    fig.add_annotation(xref="paper", yref="paper", x=0.004, y=0.995, xanchor="left",
                       yanchor="top", showarrow=False, text="<b>a</b>",
                       font=dict(size=24, color="#111"))
    fig.add_annotation(xref="paper", yref="paper", x=0.5, y=0.055, xanchor="center",
                       yanchor="bottom", showarrow=False, align="center",
                       font=dict(size=11, color="#555"),
                       text="matched-strength gate crossed left→right · sphere size ∝ each contact's own "
                            "|ρ<sub>sym</sub>| · green = holds the trace, red = reverses · blue = strongest held "
                            "co-movements (top 7% per patient) · pooled 10 implants")

    fig.update_layout(showlegend=True,
                      legend=dict(orientation="h", x=0.5, xanchor="center", y=0.015,
                                  yanchor="bottom", font=dict(size=12),
                                  itemsizing="constant", bgcolor="rgba(0,0,0,0)"),
                      paper_bgcolor="rgba(0,0,0,0)", margin=dict(l=0, r=0, t=6, b=64),
                      width=2040, height=740)

    fig.write_html(str(OUT) + ".html", include_plotlyjs="cdn")
    fig.write_image(str(OUT) + ".png", scale=2)
    try:
        fig.write_image(str(OUT) + ".pdf", scale=3)
        pdf_ok = True
    except Exception as e:                       # kaleido PDF can be finicky; PNG is the fallback
        pdf_ok = False
        print(f"[warn] PDF export failed ({e}); PNG written, wrap it if needed")

    for band in BANDS:
        nc, na, ns, ne = counts[band]
        print(f"{band:6s}: carriers(green)={nc:3d}  resets(red)={na:3d}  sea={ns:4d}  arcs={ne:3d}")
    print(f"wrote {OUT}.{{html,png{',pdf' if pdf_ok else ''}}}")


def main_vertical():
    """Vertical variant: the same three brains stacked in a 3-row x 1-col column,
    theta on TOP -> alpha MIDDLE -> beta BOTTOM (ascending gate order, top->bottom).
    Writes fig_trace_g_pairglow_3d_vertical.{html,png,pdf} for the left column of a
    compound figure. All helpers, constants and camera are shared with ``main()``."""
    OUT_V = ROOT / "data/preprint/figures/results_section1/fig_trace_g_pairglow_3d_vertical"
    OUT_V.parent.mkdir(parents=True, exist_ok=True)

    fig = make_subplots(rows=len(BANDS), cols=1, vertical_spacing=0.02,
                        specs=[[{"type": "scene"}] for _ in BANDS])
    counts = {}
    bbox = None
    for k, band in enumerate(BANDS, start=1):
        P = build_panel(band)
        if bbox is None:
            bbox = (P["coords"].min(0) - 2, P["coords"].max(0) + 2)
        counts[band] = (int((P["pol"] == "carrier").sum()),
                        int((P["pol"] == "anti").sum()),
                        int((P["pol"] == "neutral").sum()), len(P["edges"]))
        fig.add_trace(pial_mesh(color=C_PIAL, opacity=0.11), row=k, col=1)
        for t in panel_traces(P):
            fig.add_trace(t, row=k, col=1)
        # stacked mini-brains: keep only the horizontal FRONT/BACK poles (TOP/BOTTOM
        # would collide with the per-brain titles); FRONT/BACK orient the sagittal view.
        for t in orientation_traces(*bbox, poles=[(1, +1, "FRONT", True), (1, -1, "BACK", True)]):
            fig.add_trace(t, row=k, col=1)

    # one shared legend (dummy traces on the first scene)
    leg = [("carrier — holds trace", C_TRACE, "m"),
           ("reset — reverses", C_RESET, "m"),
           ("sea — no trace", C_SEA, "m"),
           ("strongest co-movements", C_LINK, "l")]
    for name, color, kind in leg:
        common = dict(x=[None], y=[None], z=[None], name=name, showlegend=True)
        fig.add_trace(go.Scatter3d(mode="markers", marker=dict(size=16, color=color, line=dict(color=NODE_EDGE, width=1.0)), **common)
                      if kind == "m" else
                      go.Scatter3d(mode="lines", line=dict(color=color, width=9), **common),
                      row=1, col=1)

    # The brain footprint at this oblique angle is LANDSCAPE (~1.85:1, BACK→FRONT horizontal).
    # Zoomed in so it fills the full-width band below; the band aspect (~1.8) is matched to the
    # brain so it fills both width and height with only a hair of vertical margin (no bottom crop).
    cam = dict(eye=dict(x=1.00, y=-0.95, z=0.55), up=dict(x=0, y=0, z=1))
    fig.update_scenes(xaxis=dict(visible=False), yaxis=dict(visible=False),
                      zaxis=dict(visible=False), aspectmode="data", camera=cam,
                      bgcolor="rgba(0,0,0,0)")

    # FULL-WIDTH stacked bands (no side strip): each landscape brain grows to fill the whole panel
    # width, killing the left/right white. Three equal ~0.29-tall bands with small gaps, and a
    # floor at 0.07 so the bottom (beta) brain clears the legend beneath it. Titles go in-corner.
    # (band, scene-id, x-domain, y-domain)
    DOMS = [("theta", "scene",  [0.02, 0.98], [0.680, 0.980]),
            ("alpha", "scene2", [0.02, 0.98], [0.370, 0.670]),
            ("beta",  "scene3", [0.02, 0.98], [0.060, 0.360])]
    fig.update_layout(**{sid: dict(domain=dict(x=xd, y=yd)) for _, sid, xd, yd in DOMS})

    # band title + gate subtitle in each band's TOP-LEFT corner (over the sparse occipital/back
    # wedge — the brain's superior edge slopes away there, leaving white for the label).
    # x indented to 0.055 so the compound's bold panel letter "a" keeps the very top-left corner
    # (at x~0.008) without colliding with the theta title; alpha/beta share the same indent.
    for band, _sid, _xd, yd in DOMS:
        fig.add_annotation(xref="paper", yref="paper", x=0.055, y=yd[1] - 0.002, xanchor="left",
                           yanchor="top", showarrow=False, text=f"<b>{BAND_TEX[band]}</b>",
                           font=dict(size=23, color="#222"))
        fig.add_annotation(xref="paper", yref="paper", x=0.055, y=yd[1] - 0.040, xanchor="left",
                           yanchor="top", showarrow=False, text=BAND_GATE[band],
                           font=dict(size=13, color="#777"))
    # (panel letter is drawn by the compound assembler; bottom caption removed)

    fig.update_layout(showlegend=True,
                      legend=dict(orientation="h", x=0.5, xanchor="center", y=0.004,
                                  yanchor="bottom", font=dict(size=25),
                                  itemsizing="constant", bgcolor="rgba(0,0,0,0)"),
                      paper_bgcolor="rgba(0,0,0,0)", margin=dict(l=0, r=0, t=2, b=10),
                      width=760, height=2050)

    fig.write_html(str(OUT_V) + ".html", include_plotlyjs="cdn")
    fig.write_image(str(OUT_V) + ".png", scale=2)
    try:
        fig.write_image(str(OUT_V) + ".pdf", scale=3)
        pdf_ok = True
    except Exception as e:                       # kaleido PDF can be finicky; PNG is the fallback
        pdf_ok = False
        print(f"[warn] PDF export failed ({e}); PNG written, wrap it if needed")

    for band in BANDS:
        nc, na, ns, ne = counts[band]
        print(f"{band:6s}: carriers(green)={nc:3d}  resets(red)={na:3d}  sea={ns:4d}  arcs={ne:3d}")
    print(f"wrote {OUT_V}.{{html,png{',pdf' if pdf_ok else ''}}}")


if __name__ == "__main__":
    if "vertical" in sys.argv:
        main_vertical()
    else:
        main()
