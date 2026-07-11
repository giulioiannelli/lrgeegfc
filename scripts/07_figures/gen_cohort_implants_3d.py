#!/usr/bin/env python3
r"""Cohort implant overlay — 3D pooled brain, contacts as mm-space spheres, one colour per patient.

Slide 03 (sEEG dataset) hero. Replaces the flat nilearn glass-brain overlay
(``cohort_implants_overlay.pdf``) with the results-section 3D style: all ten implants pooled
inside one translucent pial shell, every contact drawn as a real mm-space SPHERE
(``visuals.brain3d.spheres_mesh`` — grows with zoom, never a flat dot), tinted by PATIENT.
The point the slide makes — sparse coverage, only where the clinic placed electrodes,
heterogeneous across patients — is exactly what a per-patient-coloured point cloud shows.

Rendered as THREE angular views of the same brain (left hemisphere · superior · right
hemisphere) so a static slide still conveys the 3D spread. Interactive HTML keeps full rotation.

HONEST: MNI is a heuristic affine (schematic, approximate); purely electrode geometry — no
connectivity, no measure. n = 10.

Reads : per-patient spatial metadata (visuals.spatial_coords.load_spatial_metadata)
Writes: data/outputs/figures/implant_in_brain/cohort_implants_overlay_3d.{html,png}
Usage : python scripts/07_figures/gen_cohort_implants_3d.py [--qa /path/qa.png]
"""
from __future__ import annotations

import sys

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from lrg_eegfc.config.paths import SEEG_DATAPATH, FIGURES_ROOT
from lrg_eegfc.visuals.brain3d import pial_mesh, spheres_mesh
from lrg_eegfc.visuals.spatial_coords import (
    load_spatial_metadata, prepare_spatial_coordinates,
)

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]

# one vivid, well-separated colour per patient (Trubetskoy distinct-20, no near-white)
PAL = {
    "Pat_02": "#e6194B", "Pat_03": "#3cb44b", "Pat_05": "#4363d8",
    "Pat_06": "#f58231", "Pat_07": "#911eb4", "Pat_08": "#279fbf",
    "Pat_10": "#f032e6", "Pat_13": "#9A6324", "Pat_14": "#469990",
    "Pat_15": "#808000",
}
R_CONTACT = 1.7                 # sphere radius (mm)
SHELL_COL, SHELL_OP = "#8b929c", 0.09

# three angular views of the SAME brain (MNI: +x RIGHT, +y FRONT, +z TOP)
VIEWS = [
    ("left hemisphere",  dict(eye=dict(x=-1.16, y=-0.15, z=0.20), up=dict(x=0, y=0, z=1))),
    ("superior",         dict(eye=dict(x=0.0,  y=-0.03, z=1.20), up=dict(x=0, y=1, z=0))),
    ("right hemisphere", dict(eye=dict(x=1.16, y=-0.15, z=0.20), up=dict(x=0, y=0, z=1))),
]
# overlapping scene x-domains pull the three brains toward centre (their cubic scenes
# pad the non-cubic brain with whitespace; overlapping the empty margins closes the gap)
DOMAINS = [dict(x=[0.00, 0.46], y=[0, 1]),
           dict(x=[0.26, 0.74], y=[0, 1]),
           dict(x=[0.54, 1.00], y=[0, 1])]
TITLE_X = [0.23, 0.50, 0.77]                      # match shifted brain centres


def patient_coords(pat):
    """Finite per-contact MNI(mm) coordinates for one patient."""
    md = load_spatial_metadata(pat, SEEG_DATAPATH)
    xyz = md[["x", "y", "z"]].to_numpy(float)
    fin = np.all(np.isfinite(xyz), axis=1)
    coords = np.asarray(prepare_spatial_coordinates(
        md.loc[fin].reset_index(drop=True), scale="mm", center=False, to_mni=True), float)
    return coords


def legend_traces(counts):
    """Dummy Scatter3d markers -> one legend swatch per patient (Mesh3d can't legend)."""
    out = []
    for pat in COHORT:
        out.append(go.Scatter3d(
            x=[None], y=[None], z=[None], mode="markers",
            marker=dict(size=9, color=PAL[pat]), showlegend=True,
            name=f"{pat.replace('Pat_', 'P')} ({counts[pat]})"))
    return out


def main():
    data = {p: patient_coords(p) for p in COHORT}
    counts = {p: len(c) for p, c in data.items()}
    total = sum(counts.values())

    fig = make_subplots(
        rows=1, cols=3, specs=[[{"type": "scene"}] * 3],
        subplot_titles=[v[0] for v in VIEWS], horizontal_spacing=0.0)

    for col in (1, 2, 3):
        fig.add_trace(pial_mesh(color=SHELL_COL, opacity=SHELL_OP), row=1, col=col)
        for pat in COHORT:
            c = data[pat]
            fig.add_trace(
                spheres_mesh(c, np.full(len(c), R_CONTACT), PAL[pat],
                             opacity=0.96, name=pat),
                row=1, col=col)

    for t in legend_traces(counts):           # legend swatches ride on scene 1
        fig.add_trace(t, row=1, col=1)

    base = dict(xaxis=dict(visible=False), yaxis=dict(visible=False),
                zaxis=dict(visible=False), aspectmode="data",
                bgcolor="rgba(0,0,0,0)")
    fig.update_layout(
        scene={**base, "camera": VIEWS[0][1], "domain": DOMAINS[0]},
        scene2={**base, "camera": VIEWS[1][1], "domain": DOMAINS[1]},
        scene3={**base, "camera": VIEWS[2][1], "domain": DOMAINS[2]},
        paper_bgcolor="rgba(0,0,0,0)", margin=dict(l=0, r=0, t=40, b=96),
        showlegend=True,
        legend=dict(orientation="h", x=0.5, xanchor="center", y=0.03,
                    yanchor="bottom", font=dict(size=16.5, color="#1b1e23"),
                    itemsizing="constant", bgcolor="rgba(0,0,0,0)",
                    tracegroupgap=0),
        width=1740, height=680)
    fig.add_annotation(
        xref="paper", yref="paper", x=0.5, y=-0.02, xanchor="center", yanchor="top",
        showarrow=False, font=dict(size=14, color="#6f757c"),
        text=f"n = 10 patients · {total} contacts · one colour per patient · "
             "approximate MNI space")
    for a, tx in zip(fig.layout.annotations[:3], TITLE_X):   # titles: darker, larger, re-centred
        a.font = dict(size=21, color="#1b1e23", family="Arial Black")
        a.x = tx

    outdir = FIGURES_ROOT / "implant_in_brain"
    outdir.mkdir(parents=True, exist_ok=True)
    html = outdir / "cohort_implants_overlay_3d.html"
    png = outdir / "cohort_implants_overlay_3d.png"
    fig.write_html(str(html), include_plotlyjs="cdn")
    fig.write_image(str(png), scale=2)
    print(f"contacts per patient: " + ", ".join(f"{p.replace('Pat_','P')}={counts[p]}" for p in COHORT))
    print(f"total {total} contacts, {len(COHORT)} patients")
    print(f"wrote {html}")
    print(f"wrote {png}")

    if "--qa" in sys.argv:                     # white-matte composite for review
        from PIL import Image
        qa = sys.argv[sys.argv.index("--qa") + 1]
        im = Image.open(png).convert("RGBA")
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        Image.alpha_composite(bg, im).convert("RGB").save(qa)
        print(f"wrote QA {qa}")


if __name__ == "__main__":
    main()
