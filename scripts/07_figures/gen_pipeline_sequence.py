#!/usr/bin/env python3
r"""Talk slide 11 — the network-analysis pipeline as SEVEN separate clean panels.

One script, seven standalone high-res transparent PNGs (Canva import; + white-matte QA
PNGs). Composed left-to-right in Canva with arrows, so each panel is a self-contained
icon with no cross-panel decoration. Every panel uses ONE (patient, band) -- default
Pat_05 / beta -- so the sequence tells a single coherent story:

    1  pipeline_seq_1_brain.png            glass-brain projection of the implant (electrodes)
    2  pipeline_seq_2_timeseries.png       a handful of raw sEEG contacts (rest_pre window)
    3  pipeline_seq_3_imcoh_matrix.png      the |ImCoh| functional-connectivity matrix (dense)
    4  pipeline_seq_4_backbone_network.png  the mst@0.20 sparse backbone (spring layout)
    5  pipeline_seq_5_specific_heat.png     C(tau) & Shat(tau) on the backbone, 4 tau marked,
                                            propagator K(tau)=exp(-tau L) inset at those 4 tau
    6  pipeline_seq_6_dendrogram.png        the cophenetic (LRG) dendrogram of the backbone
    7  pipeline_seq_7_rho.png               the trace estimator symbol rho (the output)

The scientific point of the sequence (2026-07 sparsification arc): the |ImCoh| graph is
FULLY CONNECTED and therefore single-scale -- the diffusion propagator is degenerate at
the fine scale (audit_174). Sparsifying to a connected, cycle-rich backbone (mst@0.20 =
maximum spanning tree union the strongest ~20% of edges) is what makes the propagator
telescope: C(tau) grows a mesoscale ladder (audit_175, Villegas 2025). Panels 4-5 make
that step visible; the propagator inset shows K(tau) spreading from localized (tau_min) to
uniform as tau is swept.

Nothing here is recomputed that a cache already holds for the raw FC; the backbone,
spectrum, specific heat and tree are computed fresh from the cached FC (cheap: N~100-130,
a few seconds total). House rules: use_lrg_style(), transparent PNG for Canva + white-matte
QA PNG, no suptitle, no rasterisation of vector artists, plt.close per panel.

Usage:
    /home/giulio/Documents/miniconda3/envs/lapbrain/bin/python \
        scripts/07_figures/gen_pipeline_sequence.py --patient Pat_05 --band beta
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette

from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.config.paths import SEEG_DATAPATH, FIGURES_ROOT
from lrg_eegfc.config.const import FS_OVERRIDES, DEFAULT_SAMPLE_RATE
from lrg_eegfc.utils.io import load_timeseries
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency
from lrg_eegfc.visuals.network_templates import load_probe_labels
from lrg_eegfc.visuals.spatial_coords import (
    load_spatial_metadata,
    prepare_spatial_coordinates,
)
from lrg_eegfc.utils.probe import extract_probe_labels
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction, backbone_density
from lrg_eegfc.utils.fc.heat_multiscale import (
    laplacian_eig,
    entropy_specific_heat,
    specific_heat_peaks,
    linkage_at_scale,
)

OUTDIR = FIGURES_ROOT / "talk"
QA = Path(
    "/tmp/claude-1000/-home-giulio-Documents-research-neural-networks-lrgeegfc/"
    "20c35621-890b-48ca-8413-660e8b128ed1/scratchpad"
)

BACKBONE_FRAC = 0.20                 # mst@0.20 = MST union strongest 20% of edges

# cohort (panel 1 establishing shot) — one vivid colour per patient (Trubetskoy distinct-20)
COHORT_10 = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
             "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PAT_PAL = {
    "Pat_02": "#e6194B", "Pat_03": "#3cb44b", "Pat_05": "#4363d8", "Pat_06": "#f58231",
    "Pat_07": "#911eb4", "Pat_08": "#279fbf", "Pat_10": "#f032e6", "Pat_13": "#9A6324",
    "Pat_14": "#469990", "Pat_15": "#808000",
}

# palette
DARK = "#1f2a37"          # raw sEEG traces / dendrogram ink / electrodes
BRACKET_INK = "#2b2f36"   # tree brackets
C_COL = "#c1440e"         # specific-heat curve (warm)
S_COL = "#2b6cb0"         # entropy curve (cool)
MARK_CMAP = plt.get_cmap("plasma")   # the 4 tau markers / inset frames


def _save(fig, name: str) -> Path:
    """High-res TRANSPARENT PNG for Canva import + a white-matte QA PNG, then close."""
    OUTDIR.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    png = OUTDIR / f"{name}.png"
    fig.savefig(png, bbox_inches="tight", transparent=True, pad_inches=0.03, dpi=300)
    fig.savefig(QA / f"{name}.png", bbox_inches="tight", facecolor="white", dpi=150)
    plt.close(fig)
    print(f"      wrote {png}", flush=True)
    return png


# ============================ Panel 1 — brain electrodes ============================
def _patient_mni(pat: str) -> np.ndarray:
    """Finite per-contact MNI(mm) coordinates for one patient."""
    md = load_spatial_metadata(pat, SEEG_DATAPATH)
    xyz = md[["x", "y", "z"]].to_numpy(float)
    fin = np.all(np.isfinite(xyz), axis=1)
    return np.asarray(prepare_spatial_coordinates(
        md.loc[fin].reset_index(drop=True), scale="mm", center=False, to_mni=True), float)


def panel1_brain(patient: str, *, brain: str = "cohort") -> Path:
    """3D pial brain, contacts as mm-space spheres, ONE (left-hemisphere) view.

    ``brain='cohort'`` → all 10 implants inside one translucent shell, one colour per
    patient (the establishing "here is our sEEG coverage" shot, same style as the slide-03
    cohort overlay). ``brain='patient'`` → only the example patient's contacts (dark).
    Rendered via plotly + kaleido to a transparent PNG (matches the 3D results style, not
    the flat nilearn glass brain).
    """
    import plotly.graph_objects as go
    from lrg_eegfc.visuals.brain3d import pial_mesh, spheres_mesh
    from PIL import Image

    pats = COHORT_10 if brain == "cohort" else [patient]
    # matte lighting — no specular highlight (kills the shiny "reference" glare on each bead)
    matte = dict(ambient=0.82, diffuse=0.40, specular=0.0, roughness=1.0, fresnel=0.0)
    fig = go.Figure()
    fig.add_trace(pial_mesh(color="#8b929c", opacity=0.09))
    ntot = 0
    for p in pats:
        c = _patient_mni(p)
        ntot += len(c)
        col = PAT_PAL.get(p, DARK) if brain == "cohort" else DARK
        fig.add_trace(spheres_mesh(c, np.full(len(c), 1.8), col, opacity=1.0, name=p,
                                   lighting=matte))
    # zoomed out so the frontal/temporal poles are not clipped at the frame edge
    view = dict(eye=dict(x=-1.55, y=-0.20, z=0.26), up=dict(x=0, y=0, z=1))   # left hemisphere
    fig.update_layout(
        scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False),
                   zaxis=dict(visible=False), aspectmode="data",
                   bgcolor="rgba(0,0,0,0)", camera=view),
        paper_bgcolor="rgba(0,0,0,0)", margin=dict(l=8, r=8, t=8, b=8),
        showlegend=False, width=900, height=820)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    png = OUTDIR / "pipeline_seq_1_brain.png"
    fig.write_image(str(png), scale=2)
    im = Image.open(png).convert("RGBA")                 # white-matte QA composite
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
    Image.alpha_composite(bg, im).convert("RGB").save(QA / "pipeline_seq_1_brain.png")
    print(f"      3D brain ({brain}, {len(pats)} implant(s), {ntot} contacts, left view)  "
          f"wrote {png}", flush=True)
    return png


# ============================ Panel 2 — raw sEEG ============================
def panel2_timeseries(patient: str, *, n_contacts: int, win_s: float, start_s: float) -> Path:
    """A few stacked raw sEEG contacts over a short rest_pre window (common gain)."""
    fs = FS_OVERRIDES.get(patient, DEFAULT_SAMPLE_RATE)
    ts = np.asarray(load_timeseries(patient, "rest_pre", SEEG_DATAPATH), dtype=float)
    if ts.shape[0] < ts.shape[1]:                    # -> (n_samples, n_channels)
        ts = ts.T
    n_samp, n_chan = ts.shape
    n = int(round(win_s * fs))
    s0 = int(round(min(start_s, max(0.0, n_samp / fs - win_s)) * fs))
    seg = ts[s0:s0 + n]
    seg = seg - seg.mean(0, keepdims=True)
    sel = np.linspace(0, n_chan - 1, n_contacts).round().astype(int)
    spacing = 6.0 * float(np.median(seg[:, sel].std(0))) or 1.0
    t = np.arange(n) / fs

    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    for k, ci in enumerate(sel):
        y0 = (n_contacts - 1 - k) * spacing
        ax.plot(t, seg[:, ci] + y0, lw=0.6, color=DARK, solid_capstyle="round")
    ax.set_yticks([])
    ax.set_xlim(float(t[0]), float(t[-1]))
    ax.set_xlabel("time (s)")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    print(f"      window {s0 / fs:.1f}-{s0 / fs + win_s:.1f}s @ {int(fs)} Hz  "
          f"contacts={list(map(int, sel))}", flush=True)
    return _save(fig, "pipeline_seq_2_timeseries")


# ======================= Panel 3 — |ImCoh| FC matrix =======================
def _load_fc_or_die(patient: str, phase: str, band: str) -> np.ndarray:
    try:
        A = load_fc_matrix(patient, phase, band, "imcoh_abs")
    except FileNotFoundError as e:
        raise SystemExit(f"!! missing FC cache for {patient}/{phase}/{band} (imcoh_abs): {e}")
    if A is None:
        raise SystemExit(f"!! missing FC cache for {patient}/{phase}/{band} (imcoh_abs)")
    return np.asarray(A, dtype=float)


def panel3_fc(A: np.ndarray) -> Path:
    """The |ImCoh| functional-connectivity matrix, log colour, generic FC_{ij} colorbar."""
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    plot_fc_adjacency(
        A, ax=ax, fc_method="imcoh_abs", band=None,
        colorbar_label=r"$\mathrm{FC}_{ij}$", tick_labels="generic", log_scale=True,
    )
    return _save(fig, "pipeline_seq_3_imcoh_matrix")


# ==================== Panel 4 — mst@0.20 backbone network ====================
def _shaft_colors(probes: list[str]) -> list:
    """One tab20 colour per sEEG shaft (probe), for the node markers."""
    shafts = extract_probe_labels(probes)
    uniq = sorted(set(shafts))
    cmap = plt.get_cmap("tab20", max(len(uniq), 3))
    cd = {s: cmap(i % 20) for i, s in enumerate(uniq)}
    return [cd[s] for s in shafts]


def panel4_backbone(patient: str, B: np.ndarray, Z: np.ndarray, dens: float,
                    *, beta: float = 0.88, edge_cmap: str = "magma_r",
                    ring_order: str = "shaft") -> Path:
    """The mst@0.20 backbone as a circular hierarchical-edge-bundled chord (graph-tool).

    Leaves sit evenly on a ring; the backbone edges are Holten-bundled THROUGH the LRG
    communication hierarchy ``Z`` (``get_hierarchy_control_points``, beta), so coupling
    within a community routes as tight cool arcs and cross-community coupling sweeps across
    the disk — the "cool graph-tool chord" look. ONLY the circular layout + bundled arcs +
    shaft-coloured node beads; NO dendrogram overlay, colorbars, or labels. Rendered with
    cairo to a transparent PNG.

    ``ring_order``:
      - ``"shaft"``  — beads sorted by sEEG probe, so each shaft is a contiguous coloured
        arc segment (anatomical reading). Output ``..._shaft.png``.
      - ``"tree"``   — beads in the LRG hierarchy's own leaf order (``leaves_list(Z)``), i.e.
        the order the drawing algorithm itself produces; communities become contiguous and
        the bundling tightens (the classic HEB look). Beads stay shaft-coloured so you can
        see how each probe scatters across communities. Output ``..._auto.png``.
    """
    import math
    import cairo
    import graph_tool.all as gt
    import graph_tool.draw as gtd
    from scipy.cluster.hierarchy import leaves_list

    N = B.shape[0]
    probes = list(load_probe_labels(patient)[:N])
    shafts = list(extract_probe_labels(probes))
    RAD = 100.0

    if ring_order == "shaft":
        # same-probe contacts contiguous (one coloured arc per shaft); ties by contact index
        order = sorted(range(N), key=lambda i: (shafts[i], i))
        suffix, tag = "_shaft", "shaft-sorted"
    else:
        # the drawing algorithm's own order = the LRG hierarchy leaf order (communities
        # contiguous, bundles tightest). Beads still shaft-coloured.
        order = list(map(int, leaves_list(Z)))
        suffix, tag = "_auto", "tree/auto-ordered"
    posmap = {int(l): i for i, l in enumerate(order)}
    members = {i: [i] for i in range(N)}
    for k in range(N - 1):
        members[N + k] = members[int(Z[k, 0])] + members[int(Z[k, 1])]
    depth = {2 * N - 2: 0}
    stack = [2 * N - 2]
    while stack:
        nd = stack.pop()
        if nd >= N:
            for ch in (int(Z[nd - N, 0]), int(Z[nd - N, 1])):
                depth[ch] = depth[nd] + 1
                stack.append(ch)
    maxd = max(depth.values())

    # circular-mean angle: internal nodes sit at the centroid of their leaves even when the
    # tree cluster is scattered around the (shaft-sorted) ring -> bundling stays sensible.
    leaf_ang = {i: 2 * math.pi * posmap[i] / N for i in range(N)}

    def angof(nd):
        a = [leaf_ang[l] for l in members[nd]]
        return math.atan2(sum(math.sin(x) for x in a), sum(math.cos(x) for x in a))

    t = gt.Graph(directed=True)
    t.add_vertex(N + (N - 1))
    for k in range(N - 1):
        t.add_edge(t.vertex(N + k), t.vertex(int(Z[k, 0])))
        t.add_edge(t.vertex(N + k), t.vertex(int(Z[k, 1])))
    tpos = t.new_vertex_property("vector<double>")
    for v in range(2 * N - 1):
        rad = RAD if v < N else RAD * depth[v] / maxd
        a = angof(v)
        tpos[t.vertex(v)] = [rad * math.cos(a), rad * math.sin(a)]

    # backbone edges, bundled through the hierarchy
    r, c = np.triu_indices(N, 1)
    w = B[r, c]
    keep = np.where(w > 0)[0]
    g = gt.Graph(directed=False)
    g.add_vertex(N)
    for kk in keep:
        g.add_edge(g.vertex(int(r[kk])), g.vertex(int(c[kk])))
    cts = gt.get_hierarchy_control_points(g, t, tpos, beta=beta, is_tree=False)
    vpos = g.own_property(tpos)

    ei, ej, wk = r[keep], c[keep], w[keep]
    rank = (np.argsort(np.argsort(wk)) + 1) / max(len(keep), 1)
    tw = rank ** 2.0                                   # emphasise strong edges (gentler)
    # guarantee EVERY node shows at least one clear arc: mark each node's strongest incident
    # edge (the spanning backbone reaches all nodes, so no node is truly disconnected — only
    # its weak arcs were near-invisible). Those get a visibility floor and draw on top.
    strongest = np.zeros(len(keep), bool)
    best = {}
    for idx in range(len(keep)):
        for nd in (int(ei[idx]), int(ej[idx])):
            if nd not in best or wk[idx] > best[nd][0]:
                best[nd] = (wk[idx], idx)
    for _, idx in best.values():
        strongest[idx] = True

    mag = plt.colormaps[edge_cmap]
    rng = np.random.default_rng(0)
    zrand = rng.permutation(len(keep))                 # interleave strong/weak in the stack
    ecol = g.new_edge_property("vector<double>")
    epw = g.new_edge_property("double")
    eord = g.new_edge_property("double")
    for idx, e in enumerate(g.edges()):
        col = mag(0.16 + 0.55 * tw[idx])
        alpha = 0.14 + 0.62 * tw[idx]
        width = 0.45 + 4.6 * tw[idx]
        order_z = float(zrand[idx])
        if strongest[idx]:                             # one visible arc per node
            alpha = max(alpha, 0.60)
            width = max(width, 1.7)
            order_z += len(keep)                       # draw on top of the faint web
        ecol[e] = [col[0], col[1], col[2], float(alpha)]
        epw[e] = float(width)
        eord[e] = order_z

    leafcol = _shaft_colors(probes)

    W = 1100
    surf = cairo.ImageSurface(cairo.FORMAT_ARGB32, W, W)   # transparent background
    cr = cairo.Context(surf)
    half = RAD * 1.16
    pad = W * 0.04
    scale = (W - 2 * pad) / (2 * half)
    cr.save()
    cr.translate(W / 2.0, W / 2.0)
    cr.scale(scale, -scale)
    epw_u = g.new_edge_property("double")
    for e in g.edges():
        epw_u[e] = float(epw[e]) / scale
    gtd.cairo_draw(g, vpos, cr, edge_control_points=cts, edge_color=ecol,
                   edge_pen_width=epw_u, eorder=eord,
                   vertex_fill_color=[0, 0, 0, 0], vertex_color=[0, 0, 0, 0],
                   vertex_pen_width=0.0, vertex_size=1.0)
    cr.restore()

    spacing = 2 * math.pi * RAD / N * scale
    R_node = spacing * 0.62
    for l in range(N):
        a = angof(l)
        dx = W / 2.0 + RAD * math.cos(a) * scale
        dy = W / 2.0 - RAD * math.sin(a) * scale
        rgb = leafcol[l]
        cr.new_sub_path()
        cr.arc(dx, dy, R_node, 0.0, 2 * math.pi)
        cr.set_source_rgb(float(rgb[0]), float(rgb[1]), float(rgb[2]))
        cr.fill_preserve()
        cr.set_source_rgba(0.12, 0.12, 0.14, 0.85)
        cr.set_line_width(max(0.4, R_node * 0.10))
        cr.stroke()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    png = OUTDIR / f"pipeline_seq_4_backbone_network{suffix}.png"
    surf.write_to_png(str(png))
    from PIL import Image                               # white-matte QA composite
    im = Image.open(png).convert("RGBA")
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
    Image.alpha_composite(bg, im).convert("RGB").save(QA / png.name)
    print(f"      HEB chord ({tag}): {N} nodes, {len(keep)} bundled edges (density {dens:.2f})  "
          f"wrote {png}", flush=True)
    return png


# =============== Panel 5 — specific heat + entropy + propagator inset ===============
def _propagator(ev: np.ndarray, V: np.ndarray, tau: float) -> np.ndarray:
    """Heat-kernel propagator K = exp(-tau L) = V diag(exp(-tau ev)) V^T (>= 0 exactly)."""
    return (V * np.exp(-tau * ev)) @ V.T


def panel5_specific_heat(ev: np.ndarray, V: np.ndarray, order: list[int]) -> tuple[Path, dict]:
    """C(tau) & Shat(tau) on the backbone; 4 timescales marked; propagator inset at each.

    The propagator insets share the (tau_min-tree) leaf order ``order`` so the eye can
    watch the same block structure emerge and then wash out as tau grows: at s=1 (tau_min)
    K is diagonal-dominant (localized); at large s it homogenizes toward the uniform ground
    state -- diffusion has "spread" across the graph. Each inset is per-matrix normalized
    (structure, not absolute scale) and framed in its marker colour.
    """
    res = entropy_specific_heat(ev)
    s, C, S = res["s"], res["C"], res["S"]
    lam_max, lam2 = res["lambda_max"], res["lambda_2"]
    s_hi_win = lam_max / lam2 if lam2 > 1e-30 else s.max()
    Cn = C / max(float(C.max()), 1e-12)

    pk = specific_heat_peaks(res)

    # 4 marker scales anchored to the C(tau) structure: s=1 (tau_min, fine) -> the
    # specific-heat peak(s) (the mesoscale ladder the sparsification unlocks) -> the valley
    # between them, so the four propagator insets are FOUR distinct spread states (sparse ->
    # block -> more fill -> near-uniform), all inside the resolution window. Falls back to a
    # geometric spread when <2 interior peaks are present.
    coll = np.where(S <= 0.1)[0]
    s_end = min(float(s[coll[0]]) if coll.size else float(s.max()), float(s.max()))
    peaks = np.sort(pk["s_peaks"]) if pk["s_peaks"].size else np.array([])
    peaks = peaks[(peaks > 1.3) & (peaks < s_end)]
    if peaks.size >= 2:
        p1, p2 = float(peaks[0]), float(peaks[-1])
        marks = [1.0, p1, float(np.sqrt(p1 * p2)), p2]     # fine · peak1 · valley · peak2
    elif peaks.size == 1:
        p1 = float(peaks[0])
        marks = [1.0, float(np.sqrt(p1)), p1, min(p1 * 3.0, s_end)]
    else:
        marks = list(np.geomspace(1.0, max(s_end, 4.0), 4))
    s_marks = np.array(sorted(dict.fromkeys(marks))[:4])
    mark_cols = [MARK_CMAP(0.08 + 0.80 * i / 3) for i in range(len(s_marks))]

    fig, ax = plt.subplots(figsize=(7.8, 3.3))         # wide landscape: curve left, insets right
    fig.subplots_adjust(left=0.075, right=0.60, top=0.93, bottom=0.17)

    ax.plot(s, Cn, color=C_COL, lw=2.0, label=r"$C(\tau)$")
    ax.plot(s, S, color=S_COL, lw=1.6, ls="--", label=r"$\hat{S}(\tau)$")
    ax.axvspan(1.0, s_hi_win, color="0.85", alpha=0.30, zorder=0)   # resolution window
    for sm, col in zip(s_marks, mark_cols):
        ax.axvline(sm, color=col, lw=1.4, alpha=0.9, zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(float(s.min()), float(s.max()))
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel(r"$s=\tau\,\lambda_{\max}$")
    ax.set_ylabel("norm.")
    ax.legend(fontsize=8.5, frameon=False, loc="center right")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # propagator insets in a 2x2 block to the right, colour-matched to the markers. Show the
    # OFF-diagonal communicability (diagonal zeroed) on a SHARED power-normalised scale
    # so the eye reads the spread: sparse near-neighbour reach at s=1 (only backbone edges
    # lit) -> multi-step fill-in at the mesoscale -> near-uniform at the collapse scale.
    from matplotlib.colors import PowerNorm
    Ks = []
    for sm in s_marks:
        K = _propagator(ev, V, sm / lam_max)[np.ix_(order, order)]
        np.fill_diagonal(K, 0.0)
        Ks.append(K)
    vmax = max(float(np.percentile(K[K > 0], 99.5)) if np.any(K > 0) else 1.0 for K in Ks)
    norm = PowerNorm(gamma=0.5, vmin=0.0, vmax=vmax)
    # 2x2 grid of propagator insets on the RIGHT (wide landscape). Each inset framed in its
    # mark colour = the matching vertical line on C(tau); NO s= titles -- the frame colour is
    # the key, the imshow colour carries the value. Reading order TL,TR,BL,BR = ascending s.
    iw = 0.170
    ih = iw * 7.8 / 3.3                                 # square boxes given the figure aspect
    x0, xgap, ytop, ygap = 0.635, 0.015, 0.55, 0.03
    xs = [x0, x0 + iw + xgap]
    ys = [ytop, ytop - ih - ygap]
    grid = [(xs[0], ys[0]), (xs[1], ys[0]), (xs[0], ys[1]), (xs[1], ys[1])]
    for (gx, gy), col, K in zip(grid, mark_cols, Ks):
        iax = fig.add_axes([gx, gy, iw, ih])
        iax.imshow(K, cmap="magma", norm=norm, interpolation="nearest")
        iax.set_xticks([]); iax.set_yticks([])
        for sp in iax.spines.values():
            sp.set_color(col); sp.set_linewidth(2.4)

    info = dict(n_peaks=int(pk["n_peaks"]), s_marks=[float(x) for x in s_marks],
                s_hi_win=float(s_hi_win))
    print(f"      C(tau): {pk['n_peaks']} peak(s); s_marks={info['s_marks']}", flush=True)
    return _save(fig, "pipeline_seq_5_specific_heat"), info


# ===================== Panel 6 — cophenetic dendrogram (backbone) =====================
# distinct, non-near-white partition palette (matches the fc band feel; not tab10 red/gray)
PART_PALETTE = ["#c1440e", "#2b6cb0", "#2e8b57", "#8e44ad",
                "#d19a00", "#0e7c7b", "#b5179e", "#5a6b1f"]


def panel6_dendrogram(Z: np.ndarray, cut_k: int = 12) -> tuple[Path, dict]:
    """Normalised UPGMA tree of the backbone with ONE hand-selected illustrative cut.

    Purely illustrative pipeline schematic: heights are normalised by the max merge height
    (dimensionless ultrametric distance in (0, 1]) and a SINGLE horizontal cut is placed by
    hand at a height that carves the tree into ``cut_k`` visible partitions, labelled
    psi_n(tau) (the LRG's natural partition scale of the diffusion hierarchy). Clades below
    the cut are coloured, the trunk above is grey. The real cross-phase comparison uses the
    full cophenetic distances, not this partition — the cut is a teaching device.
    """
    h = Z[:, 2]
    hmax = float(h.max())
    Zn = Z.copy()
    Zn[:, 2] = Zn[:, 2] / hmax                    # normalise: ultrametric distance in (0, 1]
    hs = np.sort(Zn[:, 2])                         # ascending merge heights (N-1 of them)
    n = Zn.shape[0] + 1                            # leaves
    # cut between the (n-cut_k)-th and (n-cut_k+1)-th smallest heights -> exactly cut_k clusters
    j = n - cut_k
    cut = float(np.sqrt(hs[j - 1] * hs[j])) if 1 <= j < hs.size else float(hs[-1] * 0.9)

    set_link_color_palette(list(PART_PALETTE))
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    dendrogram(Zn, ax=ax, no_labels=True, color_threshold=cut,
               above_threshold_color="#000000")               # trunk above the cut = black
    ax.axhline(cut, ls=(0, (5, 3)), lw=1.4, color=BRACKET_INK, zorder=5)
    ax.set_yscale("log")
    ax.set_ylim(float(Zn[:, 2].min()) * 0.7, 1.30)            # full tree (top merge = 1)
    ax.set_ylabel("ultrametric distance (norm.)")
    ax.set_xticks([])
    for sp in ("top", "right", "bottom"):
        ax.spines[sp].set_visible(False)
    set_link_color_palette(None)
    print(f"      illustrative cut -> {cut_k} partitions (norm height {cut:.2e})", flush=True)
    return _save(fig, "pipeline_seq_6_dendrogram"), dict(cut_norm=cut, cut_k=cut_k)


# ============ Panel 6 (v2) — dendrogram TENSOR: trees across tau, stacked ============
def panel6_dendrogram_tensor(ev: np.ndarray, V: np.ndarray,
                             *, n_cards: int = 4, cut_k: int = 12) -> Path:
    """v2 of the tree panel: ``n_cards`` UPGMA dendrograms at increasing diffusion time τ,
    each drawn 2D but laid on oblique-projected 'slices' receding up-and-right into the page
    — a tensor stack (messy-but-nice).

    Scales span the FULL diffusion range, from the finest (``Shat~0.98``, many communities)
    to the collapse (``Shat~0.02``, all nodes in ONE cluster = the ground state). As τ grows
    the propagator homogenises and communities MERGE, so the hand cut coarsens
    (``cut_k → 1``): the front (lower-left, opaque, finest τ) slice carves many partitions and
    each deeper slice fewer, the last (coarsest, largest-τ) slice being the single-cluster
    collapse the LRG renormalises to — the textbook multiscale coarsening (more diffusion =
    more merging = fewer clusters). Deeper slices are colour-matched to the panel-5 propagator
    insets (``MARK_CMAP``) and faded with depth; each uses that scale's native leaf order and
    its own normalised log-height. Purely illustrative; cards labelled by τ (not s).

    NB β is delocalised (the *natural* community count is ~2 at every resolved scale, then a
    singleton comb before collapse); the stepped ``cut_k`` is an illustrative teaching cut,
    not a balanced-community claim.
    """
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Polygon
    from matplotlib.colors import to_rgb

    nc = int(n_cards)
    lam_max = float(ev[-1])
    # scales span the FULL diffusion range: finest (Shat~0.98, many communities) down to the
    # collapse (Shat~0.02, all nodes in ONE cluster = the ground state). As tau grows the
    # propagator homogenises, communities MERGE, and the tree coarsens -> the last (coarsest,
    # largest-tau) slice is the single-cluster collapse the LRG renormalises to.
    res = entropy_specific_heat(ev)
    Sarr, sarr = res["S"], res["s"]
    idxc = np.where(Sarr <= 0.02)[0]
    s_collapse = float(sarr[idxc[0]]) if idxc.size else float(sarr[-1])
    # head cards at ascending order-of-magnitude tau (10^-1, 10^0, 10^1, ...); the last card is
    # the Shat->0 collapse (single cluster), shown as tau -> infinity. The tau numbers are
    # order-of-magnitude scale markers only, not precise values.
    exps = np.arange(nc - 1) - 1                             # [-1, 0, 1, ...]
    s_head = np.clip(10.0 ** exps.astype(float) * lam_max, 1.0, None)   # tau~10^e, but s>=1
    s_vals = np.append(s_head, s_collapse)
    tau_vals = s_vals / lam_max
    card_labels = [rf"$\tau\!\sim\!10^{{{int(e)}}}$" for e in exps] + [r"$\tau\!\to\!\infty$"]
    # partitions coarsen with tau: many at fine tau -> 1 (all nodes together) at the collapse
    cut_ks = np.round(np.geomspace(cut_k, 1, nc)).astype(int)
    cut_ks[-1] = 1

    # oblique parallel projection: card-local (u∈[0,1] leaf span, v∈[0,1] log-height),
    # depth d = card index. screen = u*E_u + v*E_v + d*E_d. The depth axis is laid MORE
    # horizontal (larger angle off the face-on/orthogonal view) so the stack reads as a wide,
    # tilted deck of slices rather than an upright pile.
    E_u = np.array([1.00, -0.10])         # cards tilt back a touch (tabletop feel)
    E_v = np.array([0.00, 1.000])         # height straight up
    E_d = np.array([0.62, 0.34])          # deeper cards recede mostly RIGHTWARD (~45° laid back)
    cols_scale = [MARK_CMAP(0.08 + 0.80 * i / max(nc - 1, 1)) for i in range(nc)]

    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    allpts = []
    ncl = []
    for i in range(nc - 1, -1, -1):       # draw back (coarsest) first, front (finest) last (on top)
        s = float(s_vals[i]); ck = int(cut_ks[i])
        Z = linkage_at_scale(ev, V, s)
        n = Z.shape[0] + 1
        hs = np.sort(Z[:, 2])
        j = n - ck
        cut = float(np.sqrt(hs[j - 1] * hs[j])) if 1 <= j < hs.size else float(hs[-1] * 0.9)
        set_link_color_palette(list(PART_PALETTE))
        dn = dendrogram(Z, no_plot=True, color_threshold=cut, above_threshold_color="#000000")
        set_link_color_palette(None)
        ic = np.asarray(dn["icoord"], float)
        dc = np.asarray(dn["dcoord"], float)
        lcolors = dn["color_list"]
        pos = dc[dc > 0]
        if pos.size == 0:
            continue
        collapse = ck <= 1                                    # Shat~0 slice: one cluster
        ncl.append(1 if collapse else len(set(c for c in lcolors if c != "k")))
        u = ic / ic.max()                                     # leaf span -> [0,1]
        lo, hi = np.log10(pos.min()), np.log10(dc.max())      # log-height -> [0,1] (leaves at 0)
        v = np.where(dc > 0, (np.log10(np.where(dc > 0, dc, 1.0)) - lo) / (hi - lo), 0.0)
        v = np.clip(v, 0.0, 1.0)

        d = float(i)
        depth = d / max(nc - 1, 1)
        alpha = 1.0 - 0.50 * depth          # front pops, deeper slices recede
        lw = 1.5 - 0.6 * depth
        z = 10 * (nc - i)

        cu = np.array([0.0, 1.0, 1.0, 0.0]); cv = np.array([0.0, 0.0, 1.0, 1.0])
        cp = cu[:, None] * E_u + cv[:, None] * E_v + d * E_d          # slab corners
        allpts.append(cp)
        sc = to_rgb(cols_scale[i])
        ax.add_patch(Polygon(cp, closed=True,
                             facecolor=(*sc, 0.05 * alpha + 0.015),
                             edgecolor=(*sc, 0.45 * alpha), lw=1.1, zorder=z))

        segs = []
        for k in range(ic.shape[0]):
            segs.append(u[k][:, None] * E_u + v[k][:, None] * E_v + d * E_d)   # (4,2) ⊓ path
        allpts.append(np.concatenate(segs, 0))
        # fine/intermediate slices keep the informative partition palette + black trunk; the
        # collapse slice (Shat~0) is a single hue = all nodes in one cluster (the ground state).
        seg_c = (to_rgb(cols_scale[i]) if collapse
                 else [c if c != "k" else "#000000" for c in lcolors])
        ax.add_collection(LineCollection(segs, colors=seg_c, linewidths=lw, alpha=alpha,
                                         zorder=z + 1, capstyle="round", joinstyle="round"))

        tl = 1.0 * E_v + d * E_d                                       # card top-left corner
        ax.text(tl[0] - 0.02, tl[1] + 0.03, card_labels[i], color=cols_scale[i], fontsize=11,
                ha="right", va="bottom", zorder=z + 2, fontweight="bold")

    # shared height-axis cue along the front card's (vertical) left edge: the dendrogram
    # merge-height symbol only (D-script, NOT D^coph — avoids confounding the tree height
    # with the cophenetic-distance matrix used in the cross-phase comparison).
    ax.annotate("", xy=(-0.05, 1.0), xytext=(-0.05, 0.0),
                arrowprops=dict(arrowstyle="->", color=DARK, lw=1.3))
    ax.text(-0.135, 0.5, r"$\mathcal{D}$", rotation=90,
            va="center", ha="center", fontsize=14, color=DARK)

    P = np.concatenate([p.reshape(-1, 2) for p in allpts], 0)
    xmin, ymin = P.min(0); xmax, ymax = P.max(0)
    mx, my = 0.05 * (xmax - xmin), 0.05 * (ymax - ymin)
    ax.set_xlim(xmin - mx - 0.16, xmax + mx * 3)                       # left room for the axis cue
    ax.set_ylim(ymin - my, ymax + my * 3)                             # top room for the τ labels
    ax.set_aspect("equal")
    ax.axis("off")
    print(f"      tensor: {nc} trees, τ={[round(float(x),2) for x in tau_vals]}, "
          f"partitions(front→back)={list(reversed(ncl))}", flush=True)
    return _save(fig, "pipeline_seq_6_dendrogram_tensor")


# ===================== Panel 7 — the trace estimator symbol =====================
def panel7_rho() -> Path:
    """The cross-phase cophenetic correlation symbol -- the pipeline's output."""
    fig, ax = plt.subplots(figsize=(2.6, 2.0))
    ax.text(0.5, 0.5, r"$\rho^{\mathrm{coph}}$", ha="center", va="center",
            fontsize=66, color=DARK)
    ax.axis("off")
    return _save(fig, "pipeline_seq_7_rho")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Seven-panel network-analysis pipeline sequence for talk slide 11."
    )
    ap.add_argument("--patient", default="Pat_05")
    ap.add_argument("--band", default="beta")
    ap.add_argument("--n-contacts", type=int, default=6, help="panel 2: stacked sEEG contacts")
    ap.add_argument("--win-s", type=float, default=2.0, help="panel 2: window length (s)")
    ap.add_argument("--start-s", type=float, default=60.0, help="panel 2: window start (s)")
    ap.add_argument("--brain", choices=["cohort", "patient"], default="cohort",
                    help="panel 1: all-implant cohort overlay (default) or the example patient only")
    ap.add_argument("--cut-k", type=int, default=12, help="panel 6: illustrative number of partitions at the psi_n cut")
    ap.add_argument("--tensor-cards", type=int, default=4, help="panel 6 v2: number of stacked dendrograms (3-4)")
    a = ap.parse_args()

    t0 = time.perf_counter()
    use_lrg_style()
    print(f"[pipeline-seq] patient={a.patient} band={a.band}  ->  {OUTDIR}", flush=True)

    print(f"[1/7] brain electrodes (3D, {a.brain})", flush=True)
    panel1_brain(a.patient, brain=a.brain)

    print("[2/7] raw sEEG timeseries", flush=True)
    panel2_timeseries(a.patient, n_contacts=a.n_contacts, win_s=a.win_s, start_s=a.start_s)

    print("[3/7] |ImCoh| FC matrix (dense)", flush=True)
    A = _load_fc_or_die(a.patient, "rest_pre", a.band)
    panel3_fc(A)

    print(f"[4/6] mst@{BACKBONE_FRAC:.2f} backbone chord network (HEB) — 2 ring orders", flush=True)
    B = mst_union_top_fraction(A, BACKBONE_FRAC)
    dens = backbone_density(B)
    # shared spectrum + tau_min tree (chord bundling hierarchy + panels 5 & 6 leaf order)
    ev, V = laplacian_eig(B)
    Z0 = linkage_at_scale(ev, V, 1.0)                 # s=1 => tau=1/lambda_max
    order = list(dendrogram(Z0, no_plot=True)["leaves"])
    panel4_backbone(a.patient, B, Z0, dens, ring_order="shaft")   # v1: anatomical
    panel4_backbone(a.patient, B, Z0, dens, ring_order="tree")    # v2: algorithm-ordered
    print(f"      backbone density = {dens:.3f}", flush=True)

    print("[5/6] specific heat C(tau) & entropy + propagator inset", flush=True)
    panel5_specific_heat(ev, V, order)

    print("[6/6] cophenetic dendrogram (backbone) + illustrative psi_n cut", flush=True)
    panel6_dendrogram(Z0, cut_k=a.cut_k)
    print("[6/6 v2] dendrogram TENSOR — trees stacked across tau, 3D", flush=True)
    panel6_dendrogram_tensor(ev, V, n_cards=a.tensor_cards, cut_k=a.cut_k)

    # NB: the trace-estimator symbol (rho^coph) lives on the dedicated rho-measure slide,
    # NOT this pipeline slide — panel7_rho() is retained in the module but not emitted here.
    print(f"[done] {time.perf_counter() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
