#!/usr/bin/env python3
r"""fig:epi_seed_spread (hero, panel a of §3) — the SOZ found as a spreading process.

CORE MESSAGE (R3.1, made visual): the seizure-onset zone is not a contact-by-contact
label; it is a co-diffusing group. Put heat on THREE known seizure contacts on ONE
electrode (the seeds) and let it diffuse on the |ImCoh| connectivity graph via the
heat-kernel propagator e^{-tau L} at the coarse scale tau_5 = 10/lambda_max. The heat
does not stay local: it flows to the OTHER, unseeded seizure contacts on a DIFFERENT
electrode and largely spares healthy tissue. Because the |ImCoh| substrate is zero-lag-
immune (Nolte 2004), this spread cannot be volume conduction or same-shaft leakage.

Exemplar Pat_08 (delta, rest_post), a left-temporal implant. The three onset contacts on
electrode T are the seeds; the six onset contacts on electrode W are the off-shaft ground
truth the spread must rediscover (off-shaft AUC 0.99 here; cohort delta 0.72, 8/10 in
panel b). A DETECTION read, four glyphs:
  - GOLD diamond      = seed (a known onset contact used to place heat);
  - node COLOUR+SIZE  = strength-residualised seed affinity RANK (the heat reaching a
                        contact, hub/strength component regressed out -> the score);
  - CRIMSON ring      = true seizure onset, unseeded (the ground truth) -- warm inside = a
                        HIT (rediscovered), cool inside would be a miss;
  - STEEL open square  = a false alarm (healthy contact the score ranks in the shortlist);
  - warm ARC          = seed-group -> onset diffusion flow, width/opacity ∝ heat-kernel
                        affinity the target receives (a real measure, not decoration).
Read together: crimson rings sitting on warm nodes = the onset zone recovered; the few
steel squares = the honest false-alarm cost (triage, ~60% top-list precision, panel c).

ILLUSTRATIVE, ONE PATIENT. Mechanism on one exemplar; the cohort claim is the quantitative
panels (relational AUC, calibrated detector) with their matched-strength / label-shuffle
nulls. MNI is a heuristic affine (schematic); Pat_08 is a left-only implant so the left
pial shell is drawn.

Reads : |ImCoh| FC (imcoh_abs, Pat_08 delta rest_post); implant SOZ labels + coords.
Writes: data/preprint/figures/results_section3/fig_epi_seed_spread_3d.{html,png,pdf}
"""
from __future__ import annotations

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

import numpy as np
import plotly.graph_objects as go
from scipy.linalg import eigh
from scipy.stats import mannwhitneyu, rankdata

from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.visuals.brain3d import pial_mesh, bezier_arcs, orientation_traces
from lrg_eegfc.visuals.spatial_coords import (
    load_spatial_metadata, prepare_spatial_coordinates,
)
from lrg_eegfc.workflow.fc import load_fc_matrix

OUT = ROOT / "data/preprint/figures/results_section3/fig_epi_seed_spread_3d"

PATIENT, BAND, PHASE = "Pat_08", "delta", "rest_post"
HEMI = "left"                           # Pat_08 is a left-only implant
TAU_FRAC = 10.0                         # tau_5 = 10 / lambda_max (coarse marker scale)
N_SEED = 3                              # k=3 known onset seeds (matches the detector)
TOP_ARC = 7                             # top-K shortlist (false-alarm markers + precision)
ARC_MODE = "topk"                       # arcs drawn: "topk" (paper: strongest flows, label-blind),
#                                         "onset" (seed->every onset, a star), or "community" (the
#                                         onset↔onset co-diffusion web — each links to its partners)
ARC_KNN = 2                             # community mode: strongest heat partners per contact

C_SEED = "#f2c200"                      # gold: known seed seizure contacts
C_TRUTH = "#0a0c0f"                     # near-black ring: true (unseeded) onset — ground truth
C_FALARM = "#586170"                    # steel: false alarm (healthy in the shortlist)
C_FLOW = "#e8681f"                      # warm: seed -> onset diffusion flow
C_ARCGLOW = "#d43b2b"                   # arc glow (crimson)
C_PIAL = "#9aa0aa"
NODE_EDGE = "#20242a"
# heat colourscale for the affinity rank — cool slate (low) -> orange -> red (high);
# no near-white end (project rule).
HEAT = [[0.00, "#39424c"], [0.32, "#4f7f88"], [0.58, "#e2a63a"],
        [0.82, "#e2652a"], [1.00, "#cf2f24"]]
CENTER = np.array([-40.0, -18.0, 10.0])  # left-hemisphere centroid, arcs bow outward from it


def _clean(W):
    W = np.asarray(W, float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, None)
    return 0.5 * (W + W.T)


def compute():
    """Heat-kernel seed affinity for the exemplar; return node arrays + provenance."""
    W = _clean(load_fc_matrix(PATIENT, PHASE, BAND, "imcoh_abs"))
    N = W.shape[0]
    m = build_epi_masks(PATIENT)
    soz = np.where(m.epi_mask)[0]
    probes = np.asarray(m.probes, dtype=object)

    md = load_spatial_metadata(PATIENT, SEEG_DATAPATH)
    xyz = md[["x", "y", "z"]].to_numpy(float)
    fin = np.all(np.isfinite(xyz), axis=1)
    coords = np.asarray(prepare_spatial_coordinates(
        md.loc[fin].reset_index(drop=True), scale="mm", center=False, to_mni=True), float)

    # combinatorial Laplacian heat kernel at the coarse marker scale tau_5
    d = W.sum(1)
    L = np.diag(d) - W
    ev, U = eigh(L)
    tau = TAU_FRAC / ev[-1]
    K = (U * np.exp(-tau * ev)) @ U.T

    # SEED the seizure electrode carrying the FEWEST onset contacts (a small, known focus),
    # capped at k=3; the onset contacts on every OTHER electrode are the off-shaft ground
    # truth the spread must rediscover. This keeps the seed set small and the discovery
    # off-shaft (proximity removed by construction, and by the zero-lag-immune substrate).
    soz_shafts, counts = np.unique(probes[soz], return_counts=True)
    seed_shaft = soz_shafts[np.argmin(counts)]
    seeds = soz[probes[soz] == seed_shaft][:N_SEED]
    held_soz = soz[probes[soz] != seed_shaft]              # distant / off-shaft SOZ
    healthy = np.setdiff1d(np.arange(N), soz)

    aff = K[:, seeds].mean(1)                              # raw heat reaching each contact
    aff_r = aff - np.polyval(np.polyfit(d, aff, 1), d)     # strength-residual affinity (score)

    # off-shaft discovery AUC (held SOZ vs healthy), for the caption / verification
    u, _ = mannwhitneyu(aff_r[held_soz], aff_r[healthy], alternative="greater")
    auc = u / (len(held_soz) * len(healthy))

    # SHORTLIST = the top-K non-seed contacts by score (the detection read); its healthy members
    # are the false alarms and its onset fraction is the precision — always computed for the caption.
    nonseed = np.setdiff1d(np.arange(N), seeds)
    shortlist = nonseed[np.argsort(aff_r[nonseed])[::-1][:TOP_ARC]]
    false_alarm = np.intersect1d(shortlist, healthy)       # healthy the heat actually reaches
    prec = np.isin(shortlist, held_soz).mean()
    # ARCS = the diffusion links drawn.
    #  "topk"      : strongest flows from seeds, label-blind (paper default; honest false alarms).
    #  "onset"     : one flow seed->EACH unseeded onset contact (a star from the seeds).
    #  "community" : the co-diffusing COMMUNITY web — each seizure contact linked to its ARC_KNN
    #                strongest heat-kernel partners AMONG the seizure set, so onset contacts tie to
    #                EACH OTHER (not spokes from one node). width ∝ K[i,j]. This is the actual claim.
    if ARC_MODE == "community":
        # co-diffusion community backbone: a MAXIMUM spanning tree over the seizure set on the
        # heat-kernel affinity (== MST on 1/K) -> every onset threaded into ONE connected community
        # with no hub-star and minimal crossings; then add each contact's single strongest extra
        # partner for a little web richness. Edge weight = K[i,j].
        from scipy.sparse.csgraph import minimum_spanning_tree
        S = np.concatenate([seeds, held_soz]).astype(int)
        KS = K[np.ix_(S, S)].copy()
        D = 1.0 / (KS + 1e-12)
        np.fill_diagonal(D, 0.0)
        mst = minimum_spanning_tree(D).toarray()
        edges = {}
        for a in range(len(S)):
            for b in range(len(S)):
                if mst[a, b] > 0:
                    edges[(min(S[a], S[b]), max(S[a], S[b]))] = float(K[S[a], S[b]])
        for a in range(len(S)):                          # + 1 strongest extra partner each (web)
            order_b = np.argsort(-KS[a])
            for b in order_b:
                if b != a:
                    edges.setdefault((min(S[a], S[b]), max(S[a], S[b])), float(K[S[a], S[b]]))
                    break
        arcs = [(int(a), int(b), w) for (a, b), w in edges.items()]
    else:
        arc_dest = held_soz if ARC_MODE == "onset" else shortlist
        arcs = [(int(seeds[np.argmax(K[t, seeds])]), int(t), float(K[t, seeds].max()))
                for t in arc_dest]

    keep = fin                                             # only plot finite-coord contacts
    idx = np.where(keep)[0]
    remap = {int(g): i for i, g in enumerate(idx)}
    return dict(
        coords=coords, aff_r=aff_r[keep], N=int(keep.sum()),
        is_seed=np.isin(idx, seeds),
        is_soz_ring=np.isin(idx, held_soz),                # ALL unseeded onset (hits + misses)
        is_falarm=np.isin(idx, false_alarm),
        arcs=[(remap[s], remap[t], w) for s, t, w in arcs if keep[s] and keep[t]],
        n_seed=len(seeds), n_held=len(held_soz), seed_shaft=str(seed_shaft),
        held_shafts=sorted(set(probes[held_soz].tolist())),
        auc=float(auc), prec=float(prec), n_fa=int(len(false_alarm)),
        n_arc_onset=int(np.isin(shortlist, held_soz).sum()), n_soz=len(soz),
    )


def _rank(a):
    """Percentile rank in [0,1] — discrimination is rank-based (AUC) and the raw residual
    magnitudes are tiny, so rank gives the honest, legible dynamic range."""
    return (rankdata(a) - 1) / max(len(a) - 1, 1)


def _arc_trace(coords, src, tgt, width, opacity, color):
    xs, ys, zs = bezier_arcs(coords, [(src, tgt)], CENTER, lift=0.36)
    return go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", hoverinfo="skip",
                        showlegend=False, opacity=opacity,
                        line=dict(color=color, width=width))


def build_traces(P):
    """diffusion arcs (flow ∝ affinity) + scored nodes (rank heat) + truth rings + FP squares + seeds."""
    c, a = P["coords"], P["aff_r"]
    seed, ring, fa = P["is_seed"], P["is_soz_ring"], P["is_falarm"]
    scored = ~seed
    rk = np.zeros(len(a))
    rk[scored] = _rank(a[scored])
    px = 8.0 + 13.0 * rk ** 0.85                            # node screen size ∝ score rank

    tr = []
    # arcs: the top-affinity diffusion flows (hits AND false alarms), each from its strongest
    # seed; width & opacity ∝ the heat-kernel flow K[s,t] the target receives (a real measure)
    ws = np.array([w for _, _, w in P["arcs"]], float)
    lo, hi = ws.min(), ws.max()
    span = (hi - lo) or 1.0
    for (src, tgt, w) in P["arcs"]:
        f = (w - lo) / span                                # 0..1 relative flow
        tr.append(_arc_trace(c, src, tgt, 5.0 + 2.5 * f, 0.12, C_ARCGLOW))          # glow
        tr.append(_arc_trace(c, src, tgt, 2.0 + 3.0 * f, 0.5 + 0.42 * f, C_FLOW))    # core

    # scored contacts, coloured + sized by seed-affinity RANK (one carries the colourbar)
    tr.append(go.Scatter3d(
        x=c[scored, 0], y=c[scored, 1], z=c[scored, 2], mode="markers",
        hoverinfo="skip", showlegend=False,
        marker=dict(size=px[scored], color=rk[scored], colorscale=HEAT,
                    cmin=0.0, cmax=1.0, line=dict(color=NODE_EDGE, width=0.6),
                    showscale=True,
                    colorbar=dict(title=dict(
                        text="seed-affinity rank<br>(strength-residual)", side="right",
                        font=dict(size=13)),
                        thickness=15, len=0.5, x=0.99, y=0.5,
                        tickvals=[0, 0.5, 1.0], ticktext=["low", "mid", "high"],
                        tickfont=dict(size=11)))))

    # false alarms: compact steel "x" on healthy contacts inside the flagged shortlist —
    # a small glyph ON the node (no big outline) so it never clutters the dense cluster.
    tr.append(go.Scatter3d(
        x=c[fa, 0], y=c[fa, 1], z=c[fa, 2], mode="markers",
        hoverinfo="skip", showlegend=False,
        marker=dict(size=px[fa] * 0.52, color=C_FALARM, symbol="x",
                    line=dict(color="#2a2f36", width=1.0))))

    # ALL unseeded seizure onset (ground truth) gets a bold near-black ring hugging the node,
    # recovered OR not — so failures are visible: a warm-filled ring = a hit, a cooler ring
    # that no strong flow reaches = a near-miss. Near-black ≠ any heat colour and reads on
    # white paper.
    tr.append(go.Scatter3d(
        x=c[ring, 0], y=c[ring, 1], z=c[ring, 2], mode="markers",
        hoverinfo="skip", showlegend=False,
        marker=dict(size=px[ring] + 4.0, color="rgba(0,0,0,0)", symbol="circle-open",
                    line=dict(color=C_TRUTH, width=4.5))))

    # gold seed diamonds on top
    tr.append(go.Scatter3d(
        x=c[seed, 0], y=c[seed, 1], z=c[seed, 2], mode="markers",
        hoverinfo="skip", showlegend=False,
        marker=dict(size=12, color=C_SEED, symbol="diamond",
                    line=dict(color="#3a2e00", width=1.4))))
    return tr


def _legend_traces():
    L = [("seed — known onset (×3)", C_SEED, "diamond"),
         ("seizure onset, unseeded — fill = recovery", C_TRUTH, "ring"),
         ("false alarm — healthy the heat reaches", C_FALARM, "x"),
         ("strongest diffusion flow (∝ heat)", C_FLOW, "line")]
    out = []
    for name, col, kind in L:
        if kind == "line":
            out.append(go.Scatter3d(x=[None], y=[None], z=[None], mode="lines",
                                    name=name, showlegend=True,
                                    line=dict(color=col, width=8)))
        elif kind == "ring":
            out.append(go.Scatter3d(x=[None], y=[None], z=[None], mode="markers",
                                    name=name, showlegend=True,
                                    marker=dict(size=15, color="rgba(0,0,0,0)",
                                                symbol="circle-open",
                                                line=dict(color=col, width=4.5))))
        elif kind == "x":
            out.append(go.Scatter3d(x=[None], y=[None], z=[None], mode="markers",
                                    name=name, showlegend=True,
                                    marker=dict(size=11, color=col, symbol="x",
                                                line=dict(color="#eef1f4", width=1.2))))
        else:
            out.append(go.Scatter3d(x=[None], y=[None], z=[None], mode="markers",
                                    name=name, showlegend=True,
                                    marker=dict(size=13, color=col, symbol="diamond",
                                                line=dict(color="#3a2e00", width=1.2))))
    return out


def make_figure(P):
    fig = go.Figure()
    pm = pial_mesh(color=C_PIAL, opacity=0.11, hemi=HEMI)
    fig.add_trace(pm)
    for t in build_traces(P):
        fig.add_trace(t)
    # anatomical labels anchor to the electrode bounding box (the zoom crops the pial, so
    # pial-anchored A-P labels would fly off-frame) and point outward = brain orientation.
    # Top (axial) view: keep the in-plane poles FRONT/BACK (anterior-posterior) and LEFT/RIGHT
    # (lateral-medial); drop TOP/BOTTOM, which point at / away from the camera.
    elo = P["coords"].min(0) - 4
    ehi = P["coords"].max(0) + 4
    for t in orientation_traces(elo, ehi):
        if t.mode != "text":
            continue                                       # drop the stub lines (declutter)
        keep = [k for k, s in enumerate(t.text) if s not in ("TOP", "BOTTOM")]
        if not keep:
            continue
        t.x = tuple(t.x[k] for k in keep)
        t.y = tuple(t.y[k] for k in keep)
        t.z = tuple(t.z[k] for k in keep)
        t.text = tuple(t.text[k] for k in keep)
        fig.add_trace(t)
    for t in _legend_traces():
        fig.add_trace(t)

    # top (axial) oblique: the SOZ depth electrodes run medial->lateral (along x), so a lateral
    # view looks down their length and their contacts stack. Viewing from ABOVE (large +z, small
    # x) fans each electrode's contacts out horizontally and separates the anterior seed shaft
    # from the posterior onset shaft vertically (FRONT up). up = +y so anterior is toward the top.
    cam = dict(eye=dict(x=-0.28, y=-0.42, z=0.92), center=dict(x=-0.13, y=-0.05, z=-0.02),
               up=dict(x=0, y=1, z=0))
    fig.update_scenes(xaxis=dict(visible=False), yaxis=dict(visible=False),
                      zaxis=dict(visible=False), aspectmode="data", camera=cam,
                      bgcolor="rgba(0,0,0,0)")
    fig.update_layout(
        showlegend=True,
        legend=dict(orientation="h", x=0.5, xanchor="center", y=0.005, yanchor="bottom",
                    font=dict(size=13), itemsizing="constant", bgcolor="rgba(0,0,0,0)"),
        paper_bgcolor="rgba(0,0,0,0)", margin=dict(l=0, r=0, t=34, b=52),
        width=1180, height=820)
    return fig


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    P = compute()
    fig = make_figure(P)

    fig.write_html(str(OUT) + ".html", include_plotlyjs="cdn")
    fig.write_image(str(OUT) + ".png", scale=2)
    pdf_ok = True
    try:
        fig.write_image(str(OUT) + ".pdf", scale=3)
    except Exception as e:
        pdf_ok = False
        print(f"[warn] PDF export failed ({e}); PNG written")

    print("fig:epi_seed_spread — heat spreads from 3 seeds to the distant SOZ\n")
    print(f"  {PATIENT} {BAND} {PHASE}: N={P['N']} contacts, {P['n_soz']} seizure-onset")
    print(f"  seeds: {P['n_seed']} on electrode {P['seed_shaft']}")
    print(f"  unseeded onset marked (all, hits+misses): {P['n_held']} on {P['held_shafts']}")
    print(f"  off-shaft discovery AUC (held SOZ vs healthy, strength-residual) = {P['auc']:.3f}")
    print(f"  top-{TOP_ARC} flows drawn: {P['n_arc_onset']} onset + {P['n_fa']} false alarms "
          f"(precision {P['prec']:.2f})")
    print(f"\nwrote {OUT}.{{html,png{',pdf' if pdf_ok else ''}}}")


if __name__ == "__main__":
    main()
