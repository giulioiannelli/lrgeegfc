#!/usr/bin/env python3
r"""Pipeline hero — network + brain + hierarchical tree in ONE 3D plotly scene.

Slide "network-analysis pipeline": the whole machine in a single image.
  BOTTOM  a translucent pial brain with the implant contacts as mm-space spheres
          and the |ImCoh| functional network drawn as arcs between them (the graph);
  ROOTS   faint threads rising from every contact up to a leaf (diffusion distils
          the graph into a hierarchy);
  TOP      the cophenetic dendrogram growing UPWARD out of the head — the multiscale
          hierarchy the diffusion propagator reads off the network.

The glue between the three parts is COLOUR: a handful of mesoscale communities
(a distance cut of the cophenetic tree) tint the contact spheres, their roots, and
their clade in the tree the same hue; everything else is neutral grey. So the eye
follows a coloured brain module up its roots into a coloured branch of the tree.

Faithful to the real pipeline: contacts, the cached |ImCoh|_abs FC, and the cached
LRG cophenetic linkage (leaf i == channel i == FC row i — verified n_nodes == N).
One patient / band / phase; purely illustrative (no measure, no null). Approx MNI.

Reads : FC (workflow.fc.load_fc_matrix), LRG linkage (workflow.lrg.load_lrg_result),
        spatial metadata (visuals.spatial_coords), channel labels.
Writes: data/outputs/figures/talk/pipeline_brain_tree_3d.{html,png}
Usage : python scripts/07_figures/gen_pipeline_brain_tree_3d.py
          [--patient Pat_05] [--band beta] [--phase rest_pre]
          [--n-comms 6] [--cut-frac 0.68] [--min-size 4] [--n-edges 150]
          [--roots all|comms|none] [--qa out.png]
"""
from __future__ import annotations

import argparse
from collections import defaultdict

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import plotly.graph_objects as go
from scipy.cluster.hierarchy import fcluster, to_tree

from lrg_eegfc.config.paths import SEEG_DATAPATH, FIGURES_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.visuals.brain3d import pial_mesh, spheres_mesh, bezier_arcs
from lrg_eegfc.visuals.spatial_coords import (
    load_spatial_metadata, prepare_spatial_coordinates,
)

# vivid, well-separated community hues (grey = background / no mesoscale community)
PAL = ["#e6194B", "#4363d8", "#3cb44b", "#f58231", "#911eb4",
       "#00a3a3", "#f032e6", "#9A6324"]
GREY_NODE, GREY_TREE = "#a7adb6", "#c4c9d0"
SHELL_COL, SHELL_OP = "#8b929c", 0.06
EDGE_COL = "#5b6b82"
R_CONTACT = 1.7


def _hex_rgba(hx, a):
    hx = hx.lstrip("#")
    r, g, b = (int(hx[i:i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{a})"


def load_aligned(pat, band, phase, fc_method):
    """FC, cophenetic linkage, and per-channel MNI coords, all on channel index i."""
    A = load_fc_matrix(pat, phase, band, fc_method)
    np.fill_diagonal(A, 0.0)
    res = load_lrg_result(pat, phase, band, fc_method)
    if res is None:
        raise SystemExit(f"no cached LRG for {pat} {band} {phase} {fc_method}")
    Z = res.linkage_matrix
    n = A.shape[0]
    if res.n_nodes != n:
        raise SystemExit(f"leaf/channel misalignment: n_nodes={res.n_nodes} N={n}")
    md = load_spatial_metadata(pat, SEEG_DATAPATH)
    xyz = md[["x", "y", "z"]].to_numpy(float)
    fin = np.all(np.isfinite(xyz), axis=1)
    coords = np.full((n, 3), np.nan)
    coords[fin] = np.asarray(prepare_spatial_coordinates(
        md.loc[fin].reset_index(drop=True), scale="mm", center=False, to_mni=True),
        float)
    return A, Z, coords, fin


def communities(Z, n, cut_frac, n_comms, min_size):
    """Top-`n_comms` largest clades from a distance cut -> {leaf: colour}."""
    hmax = Z[:, 2].max()
    cl = fcluster(Z, cut_frac * hmax, criterion="distance")
    sizes = np.bincount(cl)
    order = [c for c in np.argsort(sizes)[::-1] if c > 0 and sizes[c] >= min_size]
    top = order[:n_comms]
    cmap = {c: PAL[i % len(PAL)] for i, c in enumerate(top)}
    node_col = np.array([cmap.get(cl[i], GREY_NODE) for i in range(n)], dtype=object)
    return cl, set(top), cmap, node_col


def build_tree_geometry(Z, cl, topset, cmap, coords, fin, xmap, hmap, Z_BASE, YT):
    """Leaf ranks + coloured dendrogram U-links in the x-z plane at y=YT.

    Leaves are ordered by a valid child-swap: at every merge the child whose
    contacts sit further LEFT in the brain is drawn first, so each clade rises
    over its own brain region and the roots run vertically instead of crossing.
    """
    tree = to_tree(Z, rd=False)

    mxc = {}

    def mean_x(nd):
        if nd.id in mxc:
            return mxc[nd.id]
        xs = [coords[l, 0] for l in nd.pre_order(lambda x: x.id) if fin[l]]
        mxc[nd.id] = float(np.mean(xs)) if xs else 0.0
        return mxc[nd.id]

    def ordered_leaves(nd):
        if nd.is_leaf():
            return [nd.id]
        a, b = (nd.left, nd.right) if mean_x(nd.left) <= mean_x(nd.right) \
            else (nd.right, nd.left)
        return ordered_leaves(a) + ordered_leaves(b)

    leaf_order = ordered_leaves(tree)
    rank = {lid: i for i, lid in enumerate(leaf_order)}

    xcache = {}

    def node_x(nd):
        if nd.id in xcache:
            return xcache[nd.id]
        v = xmap(rank[nd.id]) if nd.is_leaf() else 0.5 * (node_x(nd.left) + node_x(nd.right))
        xcache[nd.id] = v
        return v

    def subtree_col(nd):
        cs = {cl[l] for l in nd.pre_order(lambda x: x.id)}
        if len(cs) == 1 and next(iter(cs)) in topset:
            return cmap[next(iter(cs))]
        return GREY_TREE

    segs = defaultdict(lambda: ([], [], []))       # colour -> (xs, ys, zs)

    def draw(nd):
        if nd.is_leaf():
            return
        xL, xR = node_x(nd.left), node_x(nd.right)
        zN = hmap(nd.dist)
        zL = hmap(nd.left.dist) if not nd.left.is_leaf() else Z_BASE
        zR = hmap(nd.right.dist) if not nd.right.is_leaf() else Z_BASE
        col = subtree_col(nd)
        xs, ys, zs = segs[col]
        xs += [xL, xL, None, xR, xR, None, xL, xR, None]
        zs += [zL, zN, None, zR, zN, None, zN, zN, None]
        ys += [YT, YT, None] * 3
        draw(nd.left)
        draw(nd.right)

    draw(tree)
    return rank, segs


def bezier_roots(coords, fin, rank, cl, node_col, cmap, topset, xmap,
                 YT, Z_BASE, brain_top, which):
    """Threads from each contact up to its leaf, BUNDLED per community.

    Each community's roots pass through a shared gather point (its brain
    centroid, lifted toward the tree) so they read as one tidy coloured stream
    from a brain module into its clade — not a crossing hairball. Grey
    background contacts are dropped unless ``which='all'``.
    """
    gather = {}
    for c in topset:
        idx = [i for i in range(len(coords)) if fin[i] and cl[i] == c]
        cen = coords[idx].mean(0)
        gather[c] = np.array([cen[0], YT, brain_top + 0.60 * (Z_BASE - brain_top)])

    groups = defaultdict(lambda: ([], [], []))
    t = np.linspace(0, 1, 18)
    for i in range(len(coords)):
        if not fin[i]:
            continue
        col = node_col[i]
        if col == GREY_NODE:
            if which != "all":
                continue
            ctrl = np.array([coords[i, 0], YT, brain_top + 6.0])       # no bundling
        else:
            ctrl = gather[cl[i]]
        p = coords[i]
        q = np.array([xmap(rank[i]), YT, Z_BASE])
        curve = (np.outer((1 - t) ** 2, p) + np.outer(2 * (1 - t) * t, ctrl)
                 + np.outer(t ** 2, q))
        xs, ys, zs = groups[col]
        xs += [*curve[:, 0], None]
        ys += [*curve[:, 1], None]
        zs += [*curve[:, 2], None]
    return groups


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", default="Pat_05")
    ap.add_argument("--band", default="beta")
    ap.add_argument("--phase", default="rest_pre")
    ap.add_argument("--fc-method", default="imcoh_abs")
    ap.add_argument("--n-comms", type=int, default=6)
    ap.add_argument("--cut-frac", type=float, default=0.68)
    ap.add_argument("--min-size", type=int, default=4)
    ap.add_argument("--n-edges", type=int, default=120)
    ap.add_argument("--roots", choices=["all", "comms", "none"], default="comms")
    # default = the "idea rising, seen from BELOW" view: camera lateral (-x),
    # slightly behind (+y, which mirrors leaf-x so orange sits left), and BELOW
    # the tree (-z) looking UP — so we see the hierarchy's underside floating
    # up-left while the whole brain reads as a lateral silhouette, unobscured.
    ap.add_argument("--eye", default="-1.50,1.00,-0.45",
                    help="camera eye 'x,y,z' (MNI: +x R, +y FRONT, +z TOP)")
    ap.add_argument("--up", default="0,0,1")
    ap.add_argument("--center", default="0,0,0.28")
    ap.add_argument("--qa")
    args = ap.parse_args()
    eye = dict(zip("xyz", (float(v) for v in args.eye.split(","))))
    up = dict(zip("xyz", (float(v) for v in args.up.split(","))))
    cen = dict(zip("xyz", (float(v) for v in args.center.split(","))))

    A, Z, coords, fin = load_aligned(args.patient, args.band, args.phase, args.fc_method)
    n = A.shape[0]
    cl, topset, cmap, node_col = communities(
        Z, n, args.cut_frac, args.n_comms, args.min_size)

    fc = coords[fin]
    lo, hi = fc.min(0), fc.max(0)
    center = fc.mean(0)
    XLO, XHI = lo[0] - 16, hi[0] + 16
    YT = center[1]
    brain_top = hi[2]
    Z_BASE = brain_top + 20.0
    TREE_H = 88.0
    hmin, hmax = Z[:, 2].min(), Z[:, 2].max()
    h0, h1 = hmin * 0.8, hmax * 1.05

    def xmap(rk):
        return XLO + (rk / (n - 1)) * (XHI - XLO)

    def hmap(h):
        zf = (np.log(max(h, h0)) - np.log(h0)) / (np.log(h1) - np.log(h0))
        return Z_BASE + float(np.clip(zf, 0, 1)) * TREE_H

    rank, segs = build_tree_geometry(
        Z, cl, topset, cmap, coords, fin, xmap, hmap, Z_BASE, YT)

    fig = go.Figure()

    # 1) brain shell
    fig.add_trace(pial_mesh(color=SHELL_COL, opacity=SHELL_OP))

    # 2) FC network edges (top-|w| arcs bowing outward inside the shell)
    Wt = np.triu(A, 1)
    iu = np.triu_indices(n, 1)
    w = Wt[iu]
    keep = [k for k in np.argsort(w)[::-1] if fin[iu[0][k]] and fin[iu[1][k]]]
    edges = [(int(iu[0][k]), int(iu[1][k])) for k in keep[:args.n_edges]]
    xs, ys, zs = bezier_arcs(np.nan_to_num(coords), edges, center, n=16, lift=0.10)
    fig.add_trace(go.Scatter3d(
        x=xs, y=ys, z=zs, mode="lines", hoverinfo="skip", showlegend=False,
        line=dict(color=_hex_rgba(EDGE_COL, 0.36), width=1.5)))

    # 3) roots: contact -> leaf (bundled per community)
    if args.roots != "none":
        roots = bezier_roots(coords, fin, rank, cl, node_col, cmap, topset,
                             xmap, YT, Z_BASE, brain_top, args.roots)
        for col, (rx, ry, rz) in roots.items():
            a = 0.12 if col == GREY_NODE else 0.50
            fig.add_trace(go.Scatter3d(
                x=rx, y=ry, z=rz, mode="lines", hoverinfo="skip", showlegend=False,
                line=dict(color=_hex_rgba(col, a), width=1.3)))

    # 4) contact spheres, one Mesh3d per colour group
    for col in [GREY_NODE, *[cmap[c] for c in topset]]:
        idx = [i for i in range(n) if fin[i] and node_col[i] == col]
        if not idx:
            continue
        c = coords[idx]
        fig.add_trace(spheres_mesh(c, np.full(len(c), R_CONTACT), col,
                                   opacity=0.97))

    # 5) dendrogram U-links (grey first, coloured clades on top)
    for col in sorted(segs, key=lambda c: c != GREY_TREE):
        sx, sy, sz = segs[col]
        fig.add_trace(go.Scatter3d(
            x=sx, y=sy, z=sz, mode="lines", hoverinfo="skip", showlegend=False,
            line=dict(color=col, width=(3.4 if col != GREY_TREE else 2.4))))

    base = dict(xaxis=dict(visible=False), yaxis=dict(visible=False),
                zaxis=dict(visible=False), aspectmode="data",
                bgcolor="rgba(0,0,0,0)")
    fig.update_layout(
        scene={**base, "camera": dict(eye=eye, up=up, center=cen)},
        paper_bgcolor="rgba(0,0,0,0)", showlegend=False,
        margin=dict(l=0, r=0, t=6, b=6), width=1040, height=1040)

    outdir = FIGURES_ROOT / "talk"
    outdir.mkdir(parents=True, exist_ok=True)
    stem = "pipeline_brain_tree_3d"
    html, png = outdir / f"{stem}.html", outdir / f"{stem}.png"
    fig.write_html(str(html), include_plotlyjs="cdn")
    fig.write_image(str(png), scale=2)
    ncol = sum(1 for i in range(n) if fin[i] and node_col[i] != GREY_NODE)
    print(f"{args.patient} {args.band} {args.phase}: {n} contacts "
          f"({int(fin.sum())} placed), {len(edges)} edges, "
          f"{len(topset)} communities ({ncol} coloured contacts)")
    print("wrote", html)
    print("wrote", png)

    if args.qa:
        from PIL import Image
        im = Image.open(png).convert("RGBA")
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        Image.alpha_composite(bg, im).convert("RGB").save(args.qa)
        print("wrote QA", args.qa)


if __name__ == "__main__":
    main()
