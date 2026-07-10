#!/usr/bin/env python3
r"""fig:trace-heat-brain — the cohort trace as a 3-D heat field on the brain.

╔══════════════════════════════════════════════════════════════════════════════════╗
║ READ FIRST (do NOT forget): this is a SINGLE-CONTACT view. In it, β looks SPREAD   ║
║ OUT / posterior — NOT concentrated in OFC. That is NOT a contradiction of the      ║
║ "β → OFC" headline; it is a different, less-relevant measurement. β → OFC lives in ║
║ (a) co-moving region PAIRS (pairglow, paper Fig 1: OFC-incident co-movements ×1.70;║
║ α ×1.02, high-γ ×0.79) and (b) the GATED per-system localizer (paper Fig 2, OFC    ║
║ q=0.010). This single-contact heat misses OFC because raw contact-level carrier    ║
║ rate is dominated by a few patients' occipital contacts (OFC 0.37 ≈ cohort 0.32);  ║
║ OFC only emerges after matched-strength gating + per-patient demeaning, which the  ║
║ single-patient posterior spikes fail. So it is raw-vs-gated AND node-vs-pair.      ║
║ Do NOT relabel this figure "β → OFC"; its honest message is that β is delocalised  ║
║ at the single-contact level while α/high-γ are more orbital. Full note + numbers:  ║
║ data/preprint/figures/_drafts/fig_trace_heat_brain_README.md                       ║
╚══════════════════════════════════════════════════════════════════════════════════╝

Every contact from all ten patients is registered into a common MNI-like space and
its per-node trace call is deposited onto a 3-D grid, Gaussian-smoothed, and rendered
as a glowing volume inside a translucent fsaverage shell. Instead of ten sparse brains
or one bar per patient, the WHOLE cohort becomes one continuous field: warm where the
band's structure is RETAINED into offline rest (carrier), cool where it REVERTS (anti),
transparent where the cohort has no coverage.

Per-node signed weight = the gated polarity CALL (audit_152 rho_sym decomposition,
per_node.csv):
    w = +1 carrier (trace kept)   -1 anti (reset)   0 neutral / below own null
Neutral nodes contribute nothing, so the field is the local net of nodes that actually
cleared their per-node null -- not a faint global sign bias. (The continuous
p_minus - p_plus is available but leans ~65% positive from near-zero neutrals, which
washes out the anatomy.)

The rendered field is a spatial ENRICHMENT: local net-trace rate minus the cohort
baseline, so the colour asks "does the trace concentrate HERE more than the cohort
average?" -- not just "is it positive" (carriers outnumber antis ~4:1 everywhere, so a
raw mean is warm all over and hides the anatomy).
    m(v) = (smoothed[carrier - anti] / smoothed[coverage]) - global_net_rate
Warm = the trace concentrates here more than the cohort average; cool = reset zone.

Honesty guardrails (descriptive spatial aggregate, NOT a new inferential claim; the
beta->OFC localization is verified elsewhere under matched-strength + BH):
  - the heat is masked to voxels the cohort actually samples; cold => possibly
    UNCOVERED, not trace-free. A faint coverage envelope shows the sampled territory.
  - MNI registration is the approximate Desikan-Killiany transform (visuals.spatial_coords);
    good enough to pool, not stereotactic.
  - overlapping smoothing kernels are not independent; no per-voxel significance claimed.

Runs for every band with a per-node decomposition (per_node.csv: alpha, beta, low_gamma;
delta/theta/high_gamma have no node-level decomposition and would need a separate compute).

Reads : data/audit/per_node_trace_decomposition_rhosym/per_node.csv
        per-patient spatial metadata (visuals.spatial_coords.load_spatial_metadata)
Writes: data/preprint/figures/_drafts/fig_trace_heat_brain_{band}_DRAFT.html
"""
from __future__ import annotations

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.ndimage import gaussian_filter

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.visuals.brain3d import pial_mesh, orientation_traces
from lrg_eegfc.visuals.spatial_coords import (
    load_spatial_metadata, prepare_spatial_coordinates,
)

PN = ROOT / "data/audit/per_node_trace_decomposition_rhosym/per_node.csv"
OUTDIR = ROOT / "data/preprint/figures/_drafts"
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]

H_MM = 3.5              # voxel edge (mm)
SIGMA_MM = 9.0         # smoothing kernel sd (mm)
PAD_MM = 8.0           # grid padding beyond pial bbox
COV_FLOOR = 0.06       # fraction of peak coverage below which a voxel is "uncovered"

SHELL_COL = "#9aa1aa"  # pial shell colour (darker than the near-invisible default)
SHELL_OP = 0.14        # pial shell opacity (readable silhouette, still see-through)

# diverging: anti (blue) -> transparent zero -> carrier (orange), no near-white ink
CSCALE = [[0.0, "#22508c"], [0.28, "#4f8fd0"], [0.5, "#eeeeec"],
          [0.72, "#e8934a"], [1.0, "#bf541a"]]
OPAC = [[0.0, 0.9], [0.40, 0.03], [0.5, 0.0], [0.60, 0.03], [1.0, 0.9]]

# ---- pool cohort nodes into common MNI space --------------------------------
def load_cohort_nodes(band):
    """coords (M,3) MNI mm, w (M,) signed weight, pat (M,) patient index."""
    pn = pd.read_csv(PN)
    sub = pn[pn.band == band]
    C, W, P, S = [], [], [], []
    kept = []
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
        pol = s.polarity.to_numpy().astype(str)
        w = np.where(pol == "carrier", 1.0, np.where(pol == "anti", -1.0, 0.0))[fin]
        C.append(coords); W.append(w); P.append(np.full(len(w), pi))
        S.append(s.system.to_numpy().astype(str)[fin])
        kept.append(pat)
    return (np.vstack(C), np.concatenate(W), np.concatenate(P).astype(int),
            np.concatenate(S), kept)


def brain_bbox():
    m = pial_mesh()
    v = np.stack([np.asarray(m.x), np.asarray(m.y), np.asarray(m.z)], 1)
    return v.min(0) - PAD_MM, v.max(0) + PAD_MM


def deposit(coords, weights, lo, dims):
    """accumulate weights into the nearest voxel of an (nx,ny,nz) grid."""
    idx = np.floor((coords - lo) / H_MM).astype(int)
    idx = np.clip(idx, 0, np.array(dims) - 1)
    g = np.zeros(dims, float)
    np.add.at(g, (idx[:, 0], idx[:, 1], idx[:, 2]), weights)
    return g


# ---- build the enrichment field ---------------------------------------------
def build_fields(coords, w):
    lo, hi = brain_bbox()
    dims = tuple(int(np.ceil((hi[a] - lo[a]) / H_MM)) + 1 for a in range(3))
    sig = SIGMA_MM / H_MM

    num = gaussian_filter(deposit(coords, w, lo, dims), sig)            # signed mass
    den = gaussian_filter(deposit(coords, np.ones(len(w)), lo, dims), sig)  # coverage
    covered = den > COV_FLOOR * den.max()
    global_net = float(w.sum() / len(w))                               # cohort baseline
    local = num / np.where(den > 0, den, 1)
    m = np.where(covered, local - global_net, 0.0)                     # enrichment

    gx = lo[0] + H_MM * np.arange(dims[0])
    gy = lo[1] + H_MM * np.arange(dims[1])
    gz = lo[2] + H_MM * np.arange(dims[2])
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
    return dict(X=X, Y=Y, Z=Z, m=m, den=den, covered=covered, lo=lo,
                dims=dims, global_net=global_net)


def sample_field_at_nodes(F, coords):
    """value of the enrichment field at each contact's voxel (for anatomy check)."""
    idx = np.clip(np.floor((coords - F["lo"]) / H_MM).astype(int), 0,
                  np.array(F["dims"]) - 1)
    return F["m"][idx[:, 0], idx[:, 1], idx[:, 2]]


# ---- render -----------------------------------------------------------------
def volume_trace(F, value, clim):
    return go.Volume(
        x=F["X"].ravel(), y=F["Y"].ravel(), z=F["Z"].ravel(),
        value=value.ravel(), isomin=-clim, isomax=clim,
        opacityscale=OPAC, colorscale=CSCALE, cmin=-clim, cmax=clim,
        surface_count=13, showscale=True, caps=dict(x_show=False, y_show=False,
                                                     z_show=False),
        colorbar=dict(title="trace  (+ kept / − reset)", thickness=14, len=0.5,
                      x=0.94), name="heat", hoverinfo="skip")


def coverage_envelope(F):
    d = F["den"]
    lvl = COV_FLOOR * d.max()
    return go.Isosurface(
        x=F["X"].ravel(), y=F["Y"].ravel(), z=F["Z"].ravel(), value=d.ravel(),
        isomin=lvl, isomax=lvl, surface_count=1, showscale=False, opacity=0.035,
        colorscale=[[0, "#8a9098"], [1, "#8a9098"]],
        caps=dict(x_show=False, y_show=False, z_show=False),
        name="coverage", hoverinfo="skip")


def scene(fig):
    fig.update_layout(
        scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False),
                   zaxis=dict(visible=False), aspectmode="data",
                   camera=dict(eye=dict(x=1.55, y=-1.35, z=0.75))),
        paper_bgcolor="rgba(0,0,0,0)", margin=dict(l=0, r=0, t=0, b=0),
        showlegend=False)


def make_fig(F):
    clim = float(np.percentile(np.abs(F["m"][F["covered"]]), 98)) or 1e-3
    fig = go.Figure([pial_mesh(color=SHELL_COL, opacity=SHELL_OP),
                     coverage_envelope(F), volume_trace(F, F["m"], clim)])
    hi = np.asarray(F["lo"]) + H_MM * (np.asarray(F["dims"]) - 1)
    for t in orientation_traces(F["lo"], hi):
        fig.add_trace(t)
    scene(fig)
    return fig, clim


def render(F, band):
    fig, clim = make_fig(F)
    out = OUTDIR / f"fig_trace_heat_brain_{band}_DRAFT"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(out) + ".html", include_plotlyjs="cdn")
    print(f"wrote {out.name}.html  clim={clim:.3f}")


def report(F, coords, w, sys, band, kept):
    m, cov = F["m"], F["covered"]
    peak = np.unravel_index(np.argmax(np.where(cov, m, -np.inf)), m.shape)
    pxyz = np.array([F["lo"][a] + H_MM * peak[a] for a in range(3)])
    print(f"\n=== {band} ===  patients pooled: {len(kept)}  nodes: {len(w)}")
    print(f"carrier: {(w > 0).sum()}  anti: {(w < 0).sum()}  "
          f"cohort baseline net-rate: {F['global_net']:+.3f}")
    print(f"covered voxels: {cov.sum()}/{cov.size}   peak enrichment m={m[peak]:+.3f} "
          f"at MNI~({pxyz[0]:+.0f},{pxyz[1]:+.0f},{pxyz[2]:+.0f})")
    # does the field reproduce the carrier-by-region table? (anatomy sanity)
    fv = sample_field_at_nodes(F, coords)
    df = pd.DataFrame({"system": sys, "enrich": fv})
    g = df.groupby("system").agg(enrich=("enrich", "mean"), n=("enrich", "size"))
    g = g[g.n >= 15].sort_values("enrich", ascending=False)
    print("field enrichment by system (mean over contacts, n>=15):")
    for r in g.itertuples():
        print(f"  {r.Index:18s} {r.enrich:+.3f}  (n={r.n})")


def build(band):
    coords, w, pat, sys, kept = load_cohort_nodes(band)
    F = build_fields(coords, w)
    report(F, coords, w, sys, band, kept)
    render(F, band)


if __name__ == "__main__":
    avail = set(pd.read_csv(PN, usecols=["band"]).band.unique())
    bands = [b for b in BRAIN_BANDS_NAMES if b in avail]
    print(f"bands with per-node decomposition: {bands}")
    for b in bands:
        build(b)
