#!/usr/bin/env python3
r"""fig:trace_B_3d — the beta cophenetic trace as a FLOWY 3D brain (talk 'wow' asset).

Two rotatable 3D brains (beta | alpha), all ten implants pooled inside a translucent
pial shell. Contacts sized by their own null-aware rho_sym trace magnitude; carrier =
amber, anti = blue, neutral = faint (the sea). OFC contacts get a green ring. The hero
is the CURVED co-movement bundle: top-5% concordant carrier<->carrier cophenetic
co-movements per patient, drawn as bezier arcs bowing off the cortex, GREEN + glowing
when OFC-incident, faint amber otherwise.

Beta = a green OFC bundle. Alpha = the same carrier backbone with no OFC bundle
(control). Exports interactive HTML + a static PNG snapshot (kaleido).

HONEST: MNI is a heuristic affine (schematic); L/R contact imbalance (766/403) is
implant coverage, not lateralised trace; pairs are within-patient (block-diagonal);
the statistical claim is per-node polarity vs own surrogate + the locked OFC
enrichment (audit_155, q=0.010-0.015), not this pooled display.

Reads : data/audit/per_node_trace_decomposition_rhosym/per_node.csv
        data/reports/imcoh_continuous_trace/per_pair_split/{Pat}_{band}.npz
        fsaverage5 pial mesh (nilearn, cached)
Writes: data/preprint/figures/_drafts/fig_beta_ofc_pairglow_3d_DRAFT.{html,png}
"""
from __future__ import annotations

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from nilearn import datasets, surface

from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.visuals.spatial_coords import (
    load_spatial_metadata, prepare_spatial_coordinates,
)
PN = ROOT / "data/audit/per_node_trace_decomposition_rhosym/per_node.csv"
PP = ROOT / "data/reports/imcoh_continuous_trace/per_pair_split"
OUT = ROOT / "data/preprint/figures/_drafts/fig_beta_ofc_pairglow_3d_DRAFT"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
EDGE_TOP_PCT = 5.0
C_CAR, C_ANTI, C_NEU, C_OFC, C_OTH = "#e0742e", "#3b7fbf", "#9aa0a6", "#12b06a", "#d8a06a"

pn_all = pd.read_csv(PN)
RHO95 = float(np.percentile(np.abs(pn_all[pn_all.band.isin(["alpha", "beta"])].rho_sym), 95))
CENTER = np.array([0.0, -18.0, 12.0])          # ~brain centroid, arcs bow away from it


def brain_mesh():
    fs = datasets.fetch_surf_fsaverage("fsaverage5")
    vl, fl = surface.load_surf_mesh(fs["pial_left"])
    vr, fr = surface.load_surf_mesh(fs["pial_right"])
    V = np.vstack([vl, vr]); F = np.vstack([fl, fr + len(vl)])
    return go.Mesh3d(x=V[:, 0], y=V[:, 1], z=V[:, 2], i=F[:, 0], j=F[:, 1], k=F[:, 2],
                     color="#d2d6dc", opacity=0.05, flatshading=False, hoverinfo="skip",
                     lighting=dict(ambient=0.74, diffuse=0.20, specular=0.0),
                     showscale=False, showlegend=False, name="brain")


def _radius(rho, pol):
    """node sphere radius in mm (data space -> scales with zoom, floor = always visible)."""
    f = np.clip(np.abs(rho) / RHO95, 0, 1) ** 1.3
    r = np.where(pol == "carrier", 1.9 + 2.6 * f,
                 np.where(pol == "anti", 1.6 + 1.6 * f, 1.5))
    return r


def _unit_sphere(nlat=7, nlon=11):
    th = np.linspace(0, np.pi, nlat); ph = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    V = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1).reshape(-1, 3)
    F = []
    for a in range(nlat - 1):
        for b in range(nlon):
            b2 = (b + 1) % nlon
            v00, v01, v10, v11 = a * nlon + b, a * nlon + b2, (a + 1) * nlon + b, (a + 1) * nlon + b2
            F += [[v00, v10, v11], [v00, v11, v01]]
    return V, np.array(F)


_SPHV, _SPHF = _unit_sphere()


def spheres_mesh(centers, radii, color, opacity=1.0, name="", lo=None):
    """One Mesh3d of many mm-space spheres (constant screen behaviour: grows with zoom)."""
    n, S = len(centers), len(_SPHV)
    V = (centers[:, None, :] + radii[:, None, None] * _SPHV[None, :, :]).reshape(-1, 3)
    F = (_SPHF[None, :, :] + (np.arange(n) * S)[:, None, None]).reshape(-1, 3)
    return go.Mesh3d(x=V[:, 0], y=V[:, 1], z=V[:, 2], i=F[:, 0], j=F[:, 1], k=F[:, 2],
                     color=color, opacity=opacity, flatshading=False, hoverinfo="skip",
                     lighting=lo or dict(ambient=0.55, diffuse=0.7, specular=0.18, roughness=0.5),
                     showscale=False, showlegend=False, name=name)


def bezier_arcs(coords, edges, n=18, lift=0.32):
    """None-separated bezier polylines bowing outward from CENTER."""
    xs, ys, zs = [], [], []
    t = np.linspace(0, 1, n)
    for a, b in edges:
        pa, pb = coords[a], coords[b]
        mid = 0.5 * (pa + pb)
        out = mid - CENTER
        nrm = np.linalg.norm(out) or 1.0
        ctrl = mid + (lift * np.linalg.norm(pb - pa) + 6.0) * out / nrm
        curve = (np.outer((1 - t) ** 2, pa) + np.outer(2 * (1 - t) * t, ctrl)
                 + np.outer(t ** 2, pb))
        xs += [*curve[:, 0], None]; ys += [*curve[:, 1], None]; zs += [*curve[:, 2], None]
    return xs, ys, zs


def build_panel(band):
    sub = pn_all[pn_all.band == band]
    coords_l, pol_l, ofc_l, rho_l = [], [], [], []
    e_ofc, e_oth = [], []
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
        sysn = s.system.to_numpy().astype(str); pol = s.polarity.to_numpy().astype(str)
        rho = s.rho_sym.to_numpy(float); is_ofc = sysn == "OFC"
        local = np.full(len(s), -1, int); local[fin] = np.arange(fin.sum())
        coords_l.append(coords); pol_l.append(pol[fin]); ofc_l.append(is_ofc[fin]); rho_l.append(rho[fin])
        d = np.load(PP / f"{pat}_{band}.npz")
        i, j = d["iu_i"], d["iu_j"]; prod = d["dD_task"] * d["dD_rest"]
        # top-X% concordant movers (NOT carrier-restricted: matches the OFC-incidence
        # enrichment the claim is about; beta x1.46 vs alpha x1.12 over base ~0.10)
        keep = (d["sigma"] > 0) & fin[i] & fin[j]
        if keep.sum():
            keep &= prod >= np.percentile(prod[keep], 100 - EDGE_TOP_PCT)
        for a, b in zip(i[keep], j[keep]):
            (e_ofc if (is_ofc[a] or is_ofc[b]) else e_oth).append((off + local[a], off + local[b]))
        off += fin.sum()
    return dict(coords=np.vstack(coords_l), pol=np.concatenate(pol_l),
                ofc=np.concatenate(ofc_l), rho=np.concatenate(rho_l), e_ofc=e_ofc, e_oth=e_oth)


def panel_traces(P):
    c, pol, ofc, rho = P["coords"], P["pol"], P["ofc"], P["rho"]
    tr = []
    # faint non-OFC mover field (context; subsampled so it stays a whisper)
    xs, ys, zs = bezier_arcs(c, P["e_oth"][::5])
    tr.append(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", hoverinfo="skip", showlegend=False,
                           line=dict(color=C_OTH, width=1.0), opacity=0.10, name="other mover"))
    # OFC-incident bundle: glow halo + bright core
    xs, ys, zs = bezier_arcs(c, P["e_ofc"])
    tr.append(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", hoverinfo="skip", showlegend=False,
                           line=dict(color=C_OFC, width=10), opacity=0.09, name="OFC glow"))
    tr.append(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", hoverinfo="skip", showlegend=False,
                           line=dict(color=C_OFC, width=3.0), opacity=0.85, name="OFC co-movement"))
    # nodes as mm-space spheres (scale with zoom; visible floor so none vanish)
    r = _radius(rho, pol)
    m = pol == "neutral"
    tr.append(spheres_mesh(c[m], r[m], C_NEU, opacity=0.45, name="neutral"))
    m = pol == "anti"
    tr.append(spheres_mesh(c[m], r[m], C_ANTI, opacity=0.9, name="anti"))
    # OFC green halo (translucent, larger) behind fills
    tr.append(spheres_mesh(c[ofc], r[ofc] + 1.4, C_OFC, opacity=0.28, name="OFC halo",
                           lo=dict(ambient=0.9, diffuse=0.1, specular=0.0)))
    m = pol == "carrier"
    tr.append(spheres_mesh(c[m], r[m], C_CAR, opacity=1.0, name="carrier"))
    return tr


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    mesh = brain_mesh()
    Pb, Pa = build_panel("beta"), build_panel("alpha")

    fig = make_subplots(rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]],
                        subplot_titles=("β (13–30 Hz) · movers reach OFC (×1.5)",
                                        "α (8–13 Hz) · movers at OFC baseline (×1.1)"),
                        horizontal_spacing=0.01)
    for col, P in [(1, Pb), (2, Pa)]:
        fig.add_trace(mesh, row=1, col=col)
        for t in panel_traces(P):
            fig.add_trace(t, row=1, col=col)

    # ---- legend (dummy traces -> swatches) + one-line encoding note ----
    leg = [("contact — holds the trace (carrier)", C_CAR, "m"),
           ("contact — reverses (anti)", C_ANTI, "m"),
           ("contact — no trace (the sea)", C_NEU, "m"),
           ("OFC contact (green halo)", C_OFC, "m"),
           ("co-movement reaching OFC", C_OFC, "l"),
           ("other strong co-movement", C_OTH, "l")]
    for name, color, kind in leg:
        common = dict(x=[None], y=[None], z=[None], name=name, showlegend=True)
        fig.add_trace(go.Scatter3d(mode="markers", marker=dict(size=11, color=color), **common)
                      if kind == "m" else
                      go.Scatter3d(mode="lines", line=dict(color=color, width=6), **common),
                      row=1, col=1)
    fig.add_annotation(xref="paper", yref="paper", x=0.5, y=0.015, xanchor="center",
                       yanchor="bottom", showarrow=False, align="center", width=1500,
                       font=dict(size=11, color="#444"),
                       text="sphere size ∝ each contact's own trace magnitude |ρ<sub>sym</sub>| "
                            "(brain-mm, grows as you zoom) · links = top-5% strongest concordant "
                            "cophenetic co-movements per patient (of ~6800 pairs, within-patient), "
                            "green = touching OFC")

    cam = dict(eye=dict(x=1.02, y=-1.18, z=0.42), up=dict(x=0, y=0, z=1))
    scene = dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
                 aspectmode="data", camera=cam, bgcolor="rgba(0,0,0,0)")
    fig.update_layout(scene=scene, scene2=scene, showlegend=True,
                      legend=dict(orientation="h", x=0.5, xanchor="center", y=0.075,
                                  yanchor="bottom", font=dict(size=11.5),
                                  itemsizing="constant", bgcolor="rgba(255,255,255,0.55)"),
                      paper_bgcolor="white", margin=dict(l=0, r=0, t=34, b=64),
                      width=1560, height=650)

    fig.write_html(str(OUT) + ".html", include_plotlyjs="cdn")
    fig.write_image(str(OUT) + ".png", scale=2)
    print(f"beta  OFC/other arcs = {len(Pb['e_ofc'])}/{len(Pb['e_oth'])}, carriers {(Pb['pol']=='carrier').sum()}")
    print(f"alpha OFC/other arcs = {len(Pa['e_ofc'])}/{len(Pa['e_oth'])}, carriers {(Pa['pol']=='carrier').sum()}")
    print(f"wrote {OUT}.html and {OUT}.png")


if __name__ == "__main__":
    main()
