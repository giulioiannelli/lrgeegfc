#!/usr/bin/env python3
r"""fig:trace-laterality-brain — WHY the cohort varies: the β trace is left/bilateral.

╔══════════════════════════════════════════════════════════════════════════════════╗
║ WHAT THIS ANSWERS: not "where is the trace" but "why do some patients carry the    ║
║ β→OFC trace and others reset?". The whole cohort's contacts are pooled into one    ║
║ MNI brain; each contact is tinted by its PATIENT's β trace net (carrier−anti rate, ║
║ the stable ρ=.95 per-patient trait, audit_148). Warm clouds = carriers, cool       ║
║ clouds = resetters. The reader SEES the answer: the two cool clouds (Pat_10,       ║
║ Pat_15) are the ones implanted purely on the RIGHT.                                ║
║                                                                                    ║
║ HONESTY (do not overclaim "right hemisphere is inert"):                            ║
║  • Pooled left net +0.34 vs right +0.05 — but the low right is DRIVEN by the two   ║
║    right-only reset patients. Right CONTACTS in bilateral patients carry fine      ║
║    (+0.36 excl. Pat_10/15). So the cohort-variance driver is coverage/laterality   ║
║    at the PATIENT level, not an inert right hemisphere.                            ║
║  • There is ALSO a milder within-patient left-lean in the two bilateral patients   ║
║    (Pat_05 L−R +0.45, Pat_06 +0.21) — a real within-subject effect on top.         ║
║  • CONTROL: α is bilaterally SYMMETRIC (L +0.145 ≈ R +0.144). Running α proves      ║
║    β's left-lateralisation is not a generic implant-sampling artifact.             ║
║  • Descriptive spatial aggregate, n=10, approximate MNI. NOT a new test — it        ║
║    visualises the audit_148 trait; β→OFC is verified elsewhere (MS+BH, q=0.010).   ║
╚══════════════════════════════════════════════════════════════════════════════════╝

Companion note (numbers + reading key): fig_trace_laterality_brain_README.md.

Reads : data/audit/per_node_trace_decomposition_rhosym/per_node.csv
        per-patient spatial metadata (visuals.spatial_coords.load_spatial_metadata)
Writes: data/preprint/figures/_drafts/fig_trace_laterality_brain_{band}_DRAFT.html
"""
from __future__ import annotations

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.visuals.brain3d import pial_mesh, orientation_traces
from lrg_eegfc.visuals.spatial_coords import (
    load_spatial_metadata, prepare_spatial_coordinates,
)

PN = ROOT / "data/audit/per_node_trace_decomposition_rhosym/per_node.csv"
OUTDIR = ROOT / "data/preprint/figures/_drafts"
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["beta", "alpha", "low_gamma"]   # beta = headline; alpha = symmetric control
BAND_SYM = {"delta": "δ", "theta": "θ", "alpha": "α", "beta": "β",
            "low_gamma": "γ_low", "high_gamma": "γ_high"}

SHELL_COL = "#8b929c"
SHELL_OP = 0.10

# diverging net colour: reset (blue) -> neutral grey -> carrier (orange). No white.
CSCALE = [[0.0, "#1f4a86"], [0.28, "#4f8fd0"], [0.5, "#c9ccd1"],
          [0.72, "#e8934a"], [1.0, "#bd511a"]]


# ---- pool cohort contacts, tint each by its patient's trace net -------------
def load_cohort(band):
    """Per-contact coords (MNI), patient-net tint, |t_rel| size, patient idx; and
    per-patient (centroid, net, side, n) summaries."""
    pn = pd.read_csv(PN)
    sub = pn[pn.band == band]
    C, NET, TREL, PI = [], [], [], []
    pats = []
    for pi, pat in enumerate(COHORT):
        s = sub[sub.patient == pat].sort_values("node_idx")
        if not len(s):
            continue
        md = load_spatial_metadata(pat, SEEG_DATAPATH)
        if len(md) != len(s):
            continue
        fin = np.all(np.isfinite(md[["x", "y", "z"]].to_numpy(float)), axis=1)
        coords = np.asarray(prepare_spatial_coordinates(
            md.loc[fin].reset_index(drop=True), scale="mm", center=False,
            to_mni=True), float)
        pol = s.polarity.to_numpy().astype(str)[fin]
        trel = pd.to_numeric(s.t_rel, errors="coerce").to_numpy(float)[fin]
        net = float((pol == "carrier").mean() - (pol == "anti").mean())
        C.append(coords); PI.append(np.full(len(coords), pi))
        NET.append(np.full(len(coords), net)); TREL.append(trel)
        side = "R" if (coords[:, 0] > 0).mean() > 0.5 else "L"
        pats.append(dict(pat=pat, cen=coords.mean(0), net=net, n=len(coords),
                         side=side, fracR=float((coords[:, 0] > 0).mean())))
    return (np.vstack(C), np.concatenate(NET), np.concatenate(TREL),
            np.concatenate(PI).astype(int), pd.DataFrame(pats))


# ---- render -----------------------------------------------------------------
def contact_cloud(coords, net, trel, clim, band):
    a = np.abs(trel); p95 = np.nanpercentile(a, 95) or 1.0
    size = 2.8 + 4.8 * np.clip(np.nan_to_num(a) / p95, 0, 1)
    sym = BAND_SYM.get(band, band)
    return go.Scatter3d(
        x=coords[:, 0], y=coords[:, 1], z=coords[:, 2], mode="markers",
        marker=dict(size=size, color=net, colorscale=CSCALE, cmin=-clim, cmax=clim,
                    opacity=0.88, line=dict(width=0),
                    colorbar=dict(title=f"patient {sym} trace net<br>(+ carrier / − reset)",
                                  thickness=14, len=0.55, x=0.94)),
        hoverinfo="skip", showlegend=False, name="contacts")


def patient_labels(pats, clim):
    """Small labelled markers at each implant centroid (identity only; the net is
    already carried by the contact-cloud hue, so these stay small and unobtrusive)."""
    cen = np.vstack(pats.cen.to_numpy())
    return go.Scatter3d(
        x=cen[:, 0], y=cen[:, 1], z=cen[:, 2], mode="markers+text",
        marker=dict(symbol="diamond", size=5, color=pats.net.to_numpy(),
                    colorscale=CSCALE, cmin=-clim, cmax=clim,
                    line=dict(width=1.2, color="#15181c"), opacity=1.0),
        text=[p.replace("Pat_", "P") for p in pats.pat],
        textposition="top center",
        textfont=dict(size=11, color="#1b1e23", family="Arial Black"),
        hoverinfo="skip", showlegend=False, name="patients")


def make_fig(coords, net, trel, pats, band):
    clim = float(np.abs(pats.net).max()) or 1e-3
    lo, hi = coords.min(0) - 12, coords.max(0) + 12
    fig = go.Figure([pial_mesh(color=SHELL_COL, opacity=SHELL_OP),
                     contact_cloud(coords, net, trel, clim, band),
                     patient_labels(pats, clim)])
    for t in orientation_traces(lo, hi):
        fig.add_trace(t)
    fig.update_layout(
        scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False),
                   zaxis=dict(visible=False), aspectmode="data",
                   camera=dict(eye=dict(x=0.35, y=-1.95, z=0.85),
                               up=dict(x=0, y=0, z=1))),  # near-frontal: L/R split reads
        paper_bgcolor="rgba(0,0,0,0)", margin=dict(l=0, r=0, t=0, b=0),
        showlegend=False)
    return fig, clim


def report(coords, pats, band):
    side = np.where(coords[:, 0] < 0, "L", "R")
    print(f"\n=== {band} ===  contacts pooled: {len(coords)}  patients: {len(pats)}")
    print(pats.assign(net=pats.net.round(3), fracR=pats.fracR.round(2))
          [["pat", "side", "fracR", "n", "net"]]
          .sort_values("net", ascending=False).to_string(index=False))
    for s in ("L", "R"):
        sub = pats[pats.side == s]
        print(f"  side {s}: patients={list(sub.pat)}  mean net={sub.net.mean():+.3f}")


def build(band):
    coords, net, trel, pi, pats = load_cohort(band)
    report(coords, pats, band)
    fig, clim = make_fig(coords, net, trel, pats, band)
    out = OUTDIR / f"fig_trace_laterality_brain_{band}_DRAFT.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(out), include_plotlyjs="cdn")
    print(f"wrote {out.name}  clim={clim:.3f}")


if __name__ == "__main__":
    avail = set(pd.read_csv(PN, usecols=["band"]).band.unique())
    for b in BANDS:
        if b in avail:
            build(b)
