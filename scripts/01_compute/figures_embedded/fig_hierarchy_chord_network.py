#!/usr/bin/env python3
r"""fig:hierarchy-chord-network — graph-tool-style radial bundled network + branch-origin
hierarchy overlay (§2, single patient).

The "hierarchical community detection" look (graph-tool demo): nodes on a circle in LRG leaf
order, edges routed as wavy Holten-bundled Beziers through the community tree, coloured by
CONNECTION WEIGHT (|ImCoh|) on a magma_r heat scale (strong = dark, weak = pale) on a white
ground. Overlaid OUTSIDE the ring: the branch-origin dendrogram of the SAME LRG tree, every
merge coloured by the phase that WROTE its cophenetic height —
    amber = encoding (premises shown)   teal = inference (reasoning)   grey = baseline/reverted
    thickness / vividness  ∝  persistence into rest_post
Leaf nodes are tinted to MATCH: each contact takes the origin colour aggregated over the
persistent merges on its root-path (which stage's held structure it sits under).

Reuses the bundling path (get_hierarchy_control_points + cairo_draw) and the attribution
logic from fig_branch_origin_tree (DRY). LRG tree recomputed with the audit-65 estimator so
network, overlay and the standalone branch-origin tree are byte-identical.

Honest register: the network is descriptive (|ImCoh| weights); the encoding/inference
attribution is an illustrative magnitude split of two cophenetic increments (see
fig_branch_origin_tree), NOT a gated test — the cohort claim is the matched-strength gate.

Parametrized: --patients --bands --phase --k-comm --topk --beta --node-color.
Writes: data/preprint/figures/_drafts/fig_hierarchy_chord_network_{pat}_{band}_{phase}_DRAFT.pdf
"""
from __future__ import annotations

import argparse
import re
import sys

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import math
import cairo
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import graph_tool.all as gt
import graph_tool.draw as gtd
from scipy.cluster.hierarchy import fcluster, leaves_list

from lrg_eegfc.visuals.network_templates import _load_inputs

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
from audit_63_split_baseline_surrogate import load_phase_fc          # type: ignore
from audit_65_kc_matched_strength_surrogate import lrg_linkage       # type: ignore
from fig_branch_origin_tree import (                                 # type: ignore  (DRY)
    attribute, link_color, coph_square, leaf_origin_weights,
    AMBER, TEAL, GREY, PHASES,
)

OUTDIR = ROOT / "data/preprint/figures/_drafts"
BG = (1.0, 1.0, 1.0)                            # white ground
HEB_RADIUS = 100.0                             # ring radius (cairo_draw is scale-sensitive)
LABEL_GAP = 6.0                                # gap ring between node beads and dendrogram
DENDRO_PAD = 78.0                              # outward depth of the overlaid dendrogram
# edge-weight heat: sample only a LIGHT window of the (reversed) colormap so strong edges
# land on a mid magenta/purple, never near-black -> the bundle reads lighter/airier.
EDGE_C0, EDGE_C1 = 0.14, 0.60                   # cmap sample window (0 = palest end)
EDGE_A0, EDGE_A1 = 0.06, 0.78                   # per-edge alpha ramp weak -> strong
NODE_R_FRAC = 0.60                              # bead radius as a fraction of node spacing
STRIP = 0.12                                    # extra height (fraction of W) for colorbars
# encoding<->inference node palette: amber (encoding) -> pale neutral (tie) -> teal
# (inference). Routes THROUGH the pale midpoint so it never passes the muddy amber/teal
# blend; ends match the dendrogram's AMBER/TEAL so nodes read against the branches.
EI_PALE = np.array([0.82, 0.81, 0.78])
BAND_SYM = {"delta": "δ", "theta": "θ", "alpha": "α", "beta": "β",
            "low_gamma": "γ_low", "high_gamma": "γ_high"}


def split_label(lab):
    """'A10' -> ('A', '10'): non-digit shaft prefix + numeric index (drawn as a superscript)."""
    m = re.match(r"^(\D*)(\d*)$", lab)
    return (m.group(1), m.group(2)) if m else (lab, "")


def linkage_to_gt_tree(Z, N):
    t = gt.Graph(directed=True)
    t.add_vertex(N + (N - 1))
    for k in range(N - 1):
        t.add_edge(t.vertex(N + k), t.vertex(int(Z[k, 0])))
        t.add_edge(t.vertex(N + k), t.vertex(int(Z[k, 1])))
    return t, N + (N - 1) - 1


def heb_layout(t, Z, N, radius=HEB_RADIUS):
    """Leaves EVENLY on the outer circle (dendrogram order); internal nodes pulled toward
    the centre by tree depth -> the graph-tool community circle for any unbalanced tree."""
    order = list(leaves_list(Z))
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
    tpos = t.new_vertex_property("vector<double>")
    for v in range(2 * N - 1):
        rad = radius if v < N else radius * depth[v] / maxd
        ps = [posmap[l] for l in members[v]]
        th = 2 * math.pi * (0.5 * (min(ps) + max(ps))) / N
        tpos[t.vertex(v)] = [rad * math.cos(th), rad * math.sin(th)]
    return tpos, order, posmap, members


def reference_scales(ref_pat, ref_band):
    """Relative-persistence colour anchors from a reference (patient, band):
    (pscale_rel, tscale_rel) = (p75/medH, tot75/medH). Colouring ANOTHER figure against
    these puts it on an ABSOLUTE (cross-figure) scale — so a weak / no-trace patient renders
    grey instead of self-normalising to look as vivid as a strong one. medH divides out the
    per-patient cophenetic scale so the anchor is comparable across patients/bands."""
    Zr = lrg_linkage(load_phase_fc(ref_pat, "rest_post", ref_band))
    Nr = Zr.shape[0] + 1
    cophr = {ph: coph_square(load_phase_fc(ref_pat, ph, ref_band))[0] for ph in PHASES}
    efpr, _ = attribute(Zr, cophr)
    _, _, totr = leaf_origin_weights(Zr, Nr, efpr)
    medH = np.median(Zr[:, 2])
    p75 = np.percentile([abs(p) for _, _, p in efpr.values()], 75)
    return p75 / medH, (np.percentile(totr, 75) / medH)


def leaf_origin_colors(Z, N, efp, members, scale="rank", tscale=None):
    """Per-leaf tint on a CONTINUOUS encoding<->inference axis (amber -> pale -> teal),
    NOT a hard teal/amber cut. For each leaf, aggregate encoding (|e|) and inference (|f|)
    magnitudes over the persistent merges on its root-path (weighted by |persistence| p),
    and form the inference fraction inf/(inf+enc). The raw fraction sits in a tight band
    around 0.5 (every leaf shares the same high-persistence spine merges near the root),
    so a hard cut painted ~93% of nodes an identical teal. We instead spread it:

      scale="rank"   : rank-normalise the fraction across leaves -> full amber..teal spread.
                       A RELATIVE (within-patient) ordering: which contacts lean more toward
                       inference vs encoding than the others. ~half read each way BY
                       CONSTRUCTION (chosen for legibility of the gradient).
      scale="signed" : keep the true 0.5 tie at the pale midpoint and robustly contrast-
                       stretch the deviation (90th-pct scaling). HONEST: colour reflects the
                       absolute lean, so a mostly-inference cohort reads mostly-teal.

    Persistence sets vividness: weakly-persisted leaves stay pale. Illustrative tint matched
    to the branch-origin tree — not a gated per-node statistic."""
    enc, inf, tot = leaf_origin_weights(Z, N, efp)
    inf_frac = inf / (inf + enc + 1e-12)                    # 0 = encoding, 1 = inference
    if scale == "signed":
        d = inf_frac - 0.5
        h = np.clip(0.5 + 0.5 * d / (np.percentile(np.abs(d), 90) + 1e-12), 0, 1)
    else:                                                   # rank
        h = np.argsort(np.argsort(inf_frac)) / max(1, N - 1)
    ei = LinearSegmentedColormap.from_list("ei", [AMBER, EI_PALE, TEAL])
    if tscale is None:
        tscale = np.percentile(tot, 75) or 1e-9
    cols = []
    for l in range(N):
        hue = np.array(ei(float(h[l]))[:3])
        pw = np.clip(tot[l] / tscale, 0, 1)
        cols.append((1 - pw ** 0.7) * EI_PALE + pw ** 0.7 * hue)   # pale when weakly held
    return cols


def contact_labels(probes, N):
    """sEEG contact names (shaft letter + running index within shaft), FC-aligned."""
    counts, labels = {}, []
    for p in probes[:N]:
        counts[p] = counts.get(p, 0) + 1
        labels.append(f"{p}{counts[p]}")
    return labels


def build(pat, band, phase, k_comm, topk, beta, gamma, size, node_color, edge_cmap,
          origin_scale, ref=None):
    A, probes = _load_inputs(pat, band, phase, "imcoh_abs")
    A = np.asarray(A, float); np.fill_diagonal(A, 0.0)
    N = A.shape[0]
    labels = contact_labels(probes, N)
    Z = lrg_linkage(load_phase_fc(pat, phase, band))

    # --- branch-origin attribution across the four phases (same Z skeleton) ---
    coph = {ph: coph_square(load_phase_fc(pat, ph, band))[0] for ph in PHASES}
    efp, members_attr = attribute(Z, coph)
    pscale = np.percentile([abs(p) for _, _, p in efp.values()], 75) or 1e-9
    tscale_ref = None                                          # per-figure node scale
    if ref:                                                    # absolute reference scale
        rp, rb = (s.strip() for s in ref.split(","))
        ps_rel, ts_rel = reference_scales(rp, rb)
        medH = np.median(Z[:, 2])
        pscale = ps_rel * medH                                 # weak trace -> grey + thin
        tscale_ref = ts_rel * medH                             # weak trace -> pale nodes

    # --- edges: keep top-`topk` (0 = all) ---
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    pos = np.where(w > 0)[0]
    order_e = pos[np.argsort(w[pos])[::-1]]
    keep = order_e[:topk] if topk and topk < len(order_e) else order_e
    rank = (np.argsort(np.argsort(w[keep])) + 1) / len(keep)
    tw = rank ** gamma
    mag = plt.colormaps[edge_cmap]

    g = gt.Graph(directed=False); g.add_vertex(N)
    meta = []
    for kk in range(len(keep)):
        i, j = int(r[keep[kk]]), int(c[keep[kk]])
        g.add_edge(g.vertex(i), g.vertex(j)); meta.append(kk)

    # --- even-circle HEB layout + Holten bundling ---
    t, root = linkage_to_gt_tree(Z, N)
    tpos, order, posmap, members = heb_layout(t, Z, N)
    cts = gt.get_hierarchy_control_points(g, t, tpos, beta=beta, is_tree=False)
    vpos = g.own_property(tpos)

    # z-order is RANDOMISED (seeded), not weight-sorted: otherwise every heavy/dark edge
    # is painted last and the whole bundle reads as "dark on top, faint below". A random
    # draw order interleaves strong and weak so no single weight class dominates the stack.
    rng = np.random.default_rng(0)
    zrand = rng.permutation(len(keep))
    ecol = g.new_edge_property("vector<double>")
    epw = g.new_edge_property("double")
    eord = g.new_edge_property("double")
    for idx, (e, kk) in enumerate(zip(g.edges(), meta)):
        col = mag(EDGE_C0 + (EDGE_C1 - EDGE_C0) * tw[kk])
        ecol[e] = [col[0], col[1], col[2], float(EDGE_A0 + (EDGE_A1 - EDGE_A0) * tw[kk])]
        epw[e] = float(0.35 + 6.5 * tw[kk])
        eord[e] = float(zrand[idx])

    # --- node colours: origin (match hierarchy) or community ---
    if node_color == "community":
        comm = fcluster(Z, k_comm, criterion="maxclust")
        pal = plt.colormaps["turbo"]
        seen = []
        for l in order:
            if comm[int(l)] not in seen:
                seen.append(comm[int(l)])
        cmap_c = {cid: pal(0.05 + 0.9 * i / max(1, len(seen) - 1))
                  for i, cid in enumerate(seen)}
        leafcol = [np.array(cmap_c[comm[v]][:3]) for v in range(N)]
    else:
        leafcol = leaf_origin_colors(Z, N, efp, members_attr, scale=origin_scale,
                                     tscale=tscale_ref)
    # graph-tool's own vertices are drawn invisible; the visible nodes are the labelled
    # beads we paint by hand (device coords) so the letter+index sits cleanly inside.

    tree_base = HEB_RADIUS + LABEL_GAP            # dendrogram grows outward from here
    r_max = tree_base + DENDRO_PAD

    # span-weighted radii: each merge consumes radial room ∝ its angular span, so the
    # wide (coarse, near-full-circle) spine splits get generously separated while the
    # tiny leaf-joins pack near the ring. (A raw-height / uniform-rank map crams the
    # coarse merges on top of each other.)
    span = np.array([(max(posmap[l] for l in members[N + k])
                      - min(posmap[l] for l in members[N + k]) + 1) / N
                     for k in range(N - 1)])
    horder = np.argsort(Z[:, 2])                  # merges low -> high (fine -> coarse)
    cum = np.zeros(N - 1); acc = 0.0
    for k in horder:
        acc += span[k]; cum[k] = acc
    cum /= cum.max()
    ranknorm = np.empty(N - 1); ranknorm[horder] = np.arange(N - 1) / (N - 2)
    merge_radius = tree_base + DENDRO_PAD * (0.72 * cum + 0.28 * ranknorm)

    def radof(nd):
        return tree_base if nd < N else float(merge_radius[nd - N])

    def angof(nd):
        ps = [posmap[l] for l in members[nd]]
        return 2 * math.pi * (0.5 * (min(ps) + max(ps))) / N

    def render(cr, W):
        cr.set_source_rgb(*BG); cr.paint()
        half = r_max * 1.05
        pad = W * 0.03
        scale = (W - 2 * pad) / (2 * half)
        cr.save()
        cr.translate(W / 2.0, W / 2.0)
        cr.scale(scale, -scale)
        # network edges + nodes (nodes sit on the ring at HEB_RADIUS)
        epw_u = g.new_edge_property("double")
        for e in g.edges():
            epw_u[e] = float(epw[e]) / scale
        gtd.cairo_draw(g, vpos, cr, edge_control_points=cts, edge_color=ecol,
                       edge_pen_width=epw_u, eorder=eord,
                       vertex_fill_color=[0, 0, 0, 0], vertex_color=[0, 0, 0, 0],
                       vertex_pen_width=0.0, vertex_size=1.0)
        # --- branch-origin dendrogram, drawn OUTSIDE the ring ---
        segs = []
        for k in range(N - 1):
            nd = N + k
            col, _lw, pw = link_color(*efp[nd], pscale)
            segs.append((pw, nd, int(Z[k, 0]), int(Z[k, 1]), col))
        for pw, nd, c1, c2, col in sorted(segs, key=lambda s: s[0]):
            rm = radof(nd)
            cr.set_source_rgba(float(col[0]), float(col[1]), float(col[2]), 1.0)
            cr.set_line_width((0.7 + 3.2 * pw) / scale)
            for ch in (c1, c2):
                a, rc = angof(ch), radof(ch)
                cr.move_to(rc * math.cos(a), rc * math.sin(a))
                cr.line_to(rm * math.cos(a), rm * math.sin(a))
                cr.stroke()
            a1, a2 = angof(c1), angof(c2)
            steps = max(2, int(abs(a2 - a1) / (2 * math.pi) * 240))
            aa = np.linspace(a1, a2, steps)
            cr.move_to(rm * math.cos(aa[0]), rm * math.sin(aa[0]))
            for a in aa[1:]:
                cr.line_to(rm * math.cos(a), rm * math.sin(a))
            cr.stroke()
        cr.restore()

        # --- node beads with the contact label INSIDE: big shaft letter + small numeric
        #     index as a superscript, so the bead stays compact. Drawn in device coords
        #     (upright, unmirrored); text ink flips light/dark with the bead luminance. ---
        cx = W / 2.0
        spacing = 2 * math.pi * HEB_RADIUS / N * scale        # device px between neighbours
        R_node = spacing * NODE_R_FRAC
        f_big = R_node * 1.34
        f_sm = f_big * 0.56
        cr.select_font_face("sans-serif", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        for l in range(N):
            a = angof(l)
            dx = cx + HEB_RADIUS * math.cos(a) * scale
            dy = cx - HEB_RADIUS * math.sin(a) * scale
            rgb = leafcol[l]
            lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
            ink = (0.99, 0.99, 0.99) if lum < 0.45 else (0.10, 0.10, 0.13)
            cr.save()
            cr.translate(dx, dy)
            cr.new_sub_path()                     # else arc() draws a line from the
            cr.arc(0.0, 0.0, R_node, 0.0, 2 * math.pi)   # previous bead -> stray chords
            cr.set_source_rgb(float(rgb[0]), float(rgb[1]), float(rgb[2]))
            cr.fill_preserve()
            cr.set_source_rgba(0.12, 0.12, 0.14, 0.85)
            cr.set_line_width(max(0.4, R_node * 0.09))
            cr.stroke()
            letter, num = split_label(labels[l])
            cr.set_source_rgb(*ink)
            cr.set_font_size(f_big); tl = cr.text_extents(letter)
            cr.set_font_size(f_sm); tn = cr.text_extents(num) if num else None
            total = tl.x_advance + (tn.x_advance if num else 0.0)
            x0 = -total * 0.5
            cr.set_font_size(f_big)
            cr.move_to(x0, f_big * 0.34)                       # baseline -> centred cap
            cr.show_text(letter)
            if num:
                cr.set_font_size(f_sm)
                cr.move_to(x0 + tl.x_advance, f_big * 0.34 - f_big * 0.40)   # raised index
                cr.show_text(num)
            cr.restore()

        # --- two colorbars in the bottom strip: edge weight + node origin axis ---
        barY = W * 1.01
        hb = W * 0.026
        wb = W * 0.30
        ftitle = max(7.0, W * 0.0155)
        flab = ftitle * 0.9
        cr.select_font_face("sans-serif", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

        def cbar(x, cmap, lo, hi, title, left, right):
            grad = cairo.LinearGradient(x, barY, x + wb, barY)
            for i in range(25):
                t = i / 24.0
                c = cmap(lo + (hi - lo) * t)
                grad.add_color_stop_rgb(t, c[0], c[1], c[2])
            cr.rectangle(x, barY, wb, hb); cr.set_source(grad); cr.fill()
            cr.set_source_rgba(0.20, 0.20, 0.22, 0.85)
            cr.set_line_width(max(0.5, W * 0.0008))
            cr.rectangle(x, barY, wb, hb); cr.stroke()
            cr.set_source_rgb(0.12, 0.12, 0.15)
            cr.set_font_size(ftitle)
            te = cr.text_extents(title)
            cr.move_to(x + wb / 2 - te.width / 2, barY - ftitle * 0.55); cr.show_text(title)
            cr.set_font_size(flab)
            cr.move_to(x, barY + hb + flab * 1.2); cr.show_text(left)
            te = cr.text_extents(right)
            cr.move_to(x + wb - te.width, barY + hb + flab * 1.2); cr.show_text(right)

        ei2 = LinearSegmentedColormap.from_list("ei", [AMBER, EI_PALE, TEAL])
        cbar(W * 0.13, plt.colormaps[edge_cmap], EDGE_C0, EDGE_C1,
             r"connection weight  |ImCoh|", "weak", "strong")
        cbar(W * 0.57, ei2, 0.0, 1.0, "branch origin", "encoding", "inference")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    tag = f"_absref-{ref.split(',')[0].strip()}" if ref else ""
    stem = f"fig_hierarchy_chord_network_{pat}_{band}_{phase}{tag}_DRAFT"
    surf = cairo.PDFSurface(str(OUTDIR / f"{stem}.pdf"), float(size), float(size) * (1 + STRIP))
    render(cairo.Context(surf), float(size)); surf.finish()
    qa = "/tmp/claude-1000/-home-giulio-Documents-research-neural-networks-lrgeegfc/c2197e94-c496-4e4e-97cc-89e43564e7da/scratchpad"
    isz = int(size * 2.4)
    isurf = cairo.ImageSurface(cairo.FORMAT_ARGB32, isz, int(isz * (1 + STRIP)))
    render(cairo.Context(isurf), float(isz))
    isurf.write_to_png(f"{qa}/chordnet_{pat}_{band}_{phase}{tag}.png")
    frac_inf = np.mean([abs(f) > abs(e) for e, f, _ in efp.values()])
    print(f"wrote {stem}.pdf   N={N} edges={len(keep)} inf-merges={frac_inf:.2f} "
          f"nodes={node_color}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--patients", nargs="+", default=["Pat_06"])
    ap.add_argument("--bands", nargs="+", default=["beta"])
    ap.add_argument("--phase", default="rest_post")
    ap.add_argument("--k-comm", type=int, default=9)
    ap.add_argument("--topk", type=int, default=0, help="0 = ALL positive-weight edges")
    ap.add_argument("--beta", type=float, default=0.88)
    ap.add_argument("--gamma", type=float, default=4.0)
    ap.add_argument("--node-color", choices=["origin", "community"], default="origin",
                    help="origin = tint to match the branch-origin hierarchy (default)")
    ap.add_argument("--origin-scale", choices=["rank", "signed"], default="signed",
                    help="origin tint spread: signed = absolute lean, tie at pale (default, "
                         "honest — a no-inference patient goes amber); rank = relative "
                         "amber..teal spread (prettier but fakes a spread for null bands)")
    ap.add_argument("--edge-cmap", default="magma_r",
                    help="edge-weight colormap (sampled in a light window). e.g. rocket_r, "
                         "mako_r, YlGnBu, BuPu, PuRd, cividis_r")
    ap.add_argument("--size", type=int, default=1300)
    ap.add_argument("--ref", default=None,
                    help="absolute colour anchor 'PATIENT,BAND' (e.g. Pat_06,beta): colour "
                         "this figure on the reference's persistence scale so a weak / "
                         "no-trace patient renders grey instead of self-normalising")
    a = ap.parse_args()
    for pat in a.patients:
        for band in a.bands:
            build(pat, band, a.phase, a.k_comm, a.topk, a.beta, a.gamma, a.size,
                  a.node_color, a.edge_cmap, a.origin_scale, a.ref)


if __name__ == "__main__":
    main()
