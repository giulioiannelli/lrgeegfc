#!/usr/bin/env python3
"""diag_lrg_radial_curvy — graph-tool radial layout + LRG-routed Bezier edges.

The "iconic" graph-tool visualisation: nodes placed on a radial tree
derived from the LRG dendrogram, edges drawn as Bezier curves whose
control points pass through the internal tree nodes (same-cluster
edges arc inside the cluster's wedge; cross-cluster edges sweep
through the tree interior toward the LCA).

Renders directly via Cairo PDF backend (no matplotlib).  Same-probe
edges are shaft-coloured at fixed alpha; cross-probe edges are gray
with t-scaled alpha.  Edge thickness / alpha use the same rank-based
γ recipe locked in for the matplotlib templates.

Run:
    conda activate lapbrain
    python scripts/01_compute/diagnostics/diag_lrg_radial_curvy.py \
        --patient Pat_05 --band beta --phase rest_pre
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import cairo
import numpy as np
import graph_tool.all as gt
import graph_tool.draw as gtd

from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.visuals.network_templates import (
    _load_inputs,
    _shaft_color_map,
    matrix_to_gt,
)
from lrg_eegfc.workflow.lrg import load_lrg_result


def _linkage_to_gt_tree(linkage_matrix: np.ndarray, N: int):
    """Build a directed graph-tool tree from a scipy linkage matrix.

    Vertices 0..N-1 are leaves (the original nodes; share IDs with the
    main graph).  Vertices N..2N-2 are internal merge nodes in scipy's
    convention.  Vertex 2N-2 is the root.

    Edge convention: parent → child, per graph-tool's
    ``get_hierarchy_control_points`` docstring (edges point from root
    to leaves; leaves have out-degree 0; root is t.vertex(2N-2)).
    """
    n_internal = N - 1
    total = N + n_internal
    t = gt.Graph(directed=True)
    t.add_vertex(total)
    for k in range(n_internal):
        c1 = int(linkage_matrix[k, 0])
        c2 = int(linkage_matrix[k, 1])
        parent = N + k
        t.add_edge(t.vertex(parent), t.vertex(c1))
        t.add_edge(t.vertex(parent), t.vertex(c2))
    return t, total - 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument("--band", default="beta")
    parser.add_argument("--phase", default="rest_pre")
    parser.add_argument("--fc-method", default="imcoh_abs")
    parser.add_argument("--gamma", type=float, default=20.0)
    parser.add_argument("--width-min", type=float, default=0.5)
    parser.add_argument("--width-max", type=float, default=8.0)
    parser.add_argument("--alpha-min", type=float, default=0.10)
    parser.add_argument("--alpha-max", type=float, default=0.95)
    parser.add_argument("--same-probe-alpha", type=float, default=0.9)
    parser.add_argument(
        "--vertex-size", type=float, default=12.0,
        help="Vertex diameter in graph-tool's vertex-size units.  "
             "Below ~11 falls under graph-tool's sub-pixel threshold "
             "with our scaled CTM and renders as nothing.  Default 12.",
    )
    parser.add_argument(
        "--vertex-pen-width", type=float, default=0.8,
        help="Black outline width around each vertex.  Default 0.8.",
    )
    parser.add_argument("--output-size", type=int, default=900)
    parser.add_argument("--bias", type=float, default=0.85,
                        help="Bezier control-point bias; closer to 1 = "
                             "tighter curve toward LCA, 0.5 = mid-routing.")
    parser.add_argument(
        "--cross-min-alpha", type=float, default=0.13,
        help="Drop cross-probe edges below this rendered alpha "
             "(saves PDF weight; same-probe edges always kept).  "
             "Default 0.13 → drops bottom ~80%% of cross-probe at γ=20 "
             "rank.  Set 0.0 to draw every edge.",
    )
    parser.add_argument(
        "--straight-threshold", type=float, default=0.005,
        help="Edges with t (post-γ rank) below this are drawn as STRAIGHT "
             "LINES (no Bezier).  Bucketed by quantized (alpha, width) so "
             "thousands of fog edges collapse to a handful of multi-subpath "
             "strokes — keeps every edge in vector form at a fraction of "
             "the file weight.  Default 0.005; set 0 to disable, set 1.0 "
             "to make ALL edges straight.",
    )
    parser.add_argument(
        "--format", default="pdf", choices=["pdf", "svg"],
        help="Output vector format.  Default pdf; svg can be smaller for "
             "heavily curved Bezier paths.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "network_templates" / "_diagnostic"),
    )
    args = parser.parse_args()

    A, probes = _load_inputs(args.patient, args.band, args.phase, args.fc_method)
    lrg = load_lrg_result(args.patient, args.phase, args.band, args.fc_method)
    if lrg is None or lrg.linkage_matrix is None:
        raise SystemExit(
            f"No LRG cache for {args.patient}/{args.band}/{args.phase}/"
            f"{args.fc_method}."
        )
    N = A.shape[0]

    # --- compute t/width/alpha for every potential edge -----------------
    r, c = np.triu_indices(N, k=1)
    w_full = A[r, c]
    pos_w_idx = np.where(w_full > 0)[0]
    w_pos = w_full[pos_w_idx]
    ranks = np.argsort(np.argsort(w_pos)).astype(float)
    t_pre = (ranks + 1) / len(ranks)
    t_arr_pos = t_pre ** args.gamma
    t_arr = np.zeros_like(w_full)
    t_arr[pos_w_idx] = t_arr_pos
    widths_arr = args.width_min + t_arr * (args.width_max - args.width_min)
    alphas_arr = args.alpha_min + t_arr * (args.alpha_max - args.alpha_min)

    # --- main graph: keep all same-probe + cross-probe above threshold --
    # Filtering BEFORE graph construction is the file-size lever — Cairo
    # PDF emits ~450 bytes per Bezier curve, so dropping invisible
    # cross-probe edges saves megabytes.  Same-probe edges (the
    # biological signal) are never filtered.
    g = gt.Graph(directed=False)
    g.add_vertex(N)
    edge_meta: list[tuple[int, int, int]] = []  # (i, j, k_idx)
    n_same = 0
    n_cross_kept = 0
    n_cross_dropped = 0
    for k_idx, (ii, jj) in enumerate(zip(r, c)):
        if w_full[k_idx] <= 0:
            continue
        ii = int(ii); jj = int(jj)
        if probes[ii] == probes[jj]:
            n_same += 1
            g.add_edge(g.vertex(ii), g.vertex(jj))
            edge_meta.append((ii, jj, k_idx))
        else:
            if float(alphas_arr[k_idx]) < args.cross_min_alpha:
                n_cross_dropped += 1
                continue
            n_cross_kept += 1
            g.add_edge(g.vertex(ii), g.vertex(jj))
            edge_meta.append((ii, jj, k_idx))
    print(f"  same-probe: {n_same} (all kept), cross-probe: "
          f"{n_cross_kept} kept / {n_cross_dropped} dropped "
          f"(threshold α≥{args.cross_min_alpha})")

    # --- LRG dendrogram tree --------------------------------------------
    t, root = _linkage_to_gt_tree(lrg.linkage_matrix, N)

    # --- radial layout on the tree --------------------------------------
    tpos = gt.radial_tree_layout(t, t.vertex(root))
    # is_tree=False triggers general path-finding (works on our scipy-
    # linkage-derived binary tree where the tree-only path fails).
    cts = gt.get_hierarchy_control_points(g, t, tpos, beta=args.bias, is_tree=False)

    pos = g.own_property(tpos)

    # --- per-edge styling on the filtered graph -------------------------
    shaft_map = _shaft_color_map(probes)
    edge_color = g.new_edge_property("vector<double>")
    edge_pen_width = g.new_edge_property("double")

    # graph-tool edges iterate in insertion order; edge_meta[i] matches
    # the i-th edge of g.
    for e, (i, j, k_idx) in zip(g.edges(), edge_meta):
        if probes[i] == probes[j]:
            rgb = shaft_map[probes[i]]
            edge_color[e] = [rgb[0], rgb[1], rgb[2], float(args.same_probe_alpha)]
            edge_pen_width[e] = float(widths_arr[k_idx])
        else:
            edge_color[e] = [0.35, 0.35, 0.35, float(alphas_arr[k_idx])]
            edge_pen_width[e] = float(widths_arr[k_idx])

    vcolor = g.new_vertex_property("vector<double>")
    for v_idx in range(N):
        rgb = shaft_map[probes[v_idx]]
        vcolor[g.vertex(v_idx)] = [rgb[0], rgb[1], rgb[2], 1.0]

    # --- order edges so strong/coloured edges draw on top --------------
    eorder = g.new_edge_property("int")
    for e, (i, j, k_idx) in zip(g.edges(), edge_meta):
        if probes[i] == probes[j]:
            eorder[e] = 1_000_000  # all same-probe on top
        else:
            eorder[e] = int(t_arr[k_idx] * 1000)

    # --- compute leaf radii so we can overlay reference rings ----------
    # ``pos`` is the leaf-restricted view of ``tpos``; root is at the
    # tree's root-vertex position.  All leaves' Euclidean distances to
    # the root give the radial spread.
    root_xy = np.array(tpos[t.vertex(root)])
    leaf_xy = np.array([list(pos[g.vertex(v)]) for v in range(N)])
    leaf_radii = np.linalg.norm(leaf_xy - root_xy, axis=1)
    r_min = float(leaf_radii.min())
    r_mean = float(leaf_radii.mean())
    r_max = float(leaf_radii.max())
    print(f"  leaf-radius: min={r_min:.3f} mean={r_mean:.3f} max={r_max:.3f}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    suffix_filter = ""
    if args.cross_min_alpha > 0:
        suffix_filter += f"_xalpha{args.cross_min_alpha:g}"
    if args.straight_threshold > 0:
        suffix_filter += f"_str{args.straight_threshold:g}"
    out = out_dir / (
        f"lrg_radial_curvy_{args.patient}_{args.band}_{args.phase}_"
        f"{args.fc_method}_g{args.gamma:g}{suffix_filter}.{args.format}"
    )

    # Render directly to a Cairo PDF surface so we can paint reference
    # circles into the same context.  Avoids matplotlib's mpl-backend
    # render-bake bug with set_xlim/ylim after gt.graph_draw.
    W = H = float(args.output_size)
    surface = cairo.PDFSurface(str(out), W, H)
    cr = cairo.Context(surface)

    # White background.
    cr.set_source_rgb(1.0, 1.0, 1.0)
    cr.paint()

    # Layout extent: pick a square bbox centered on root_xy with side
    # 2 * (r_max * 1.18) so the outermost reference ring + a small
    # margin fit inside.
    half_extent = r_max * 1.18
    canvas_pad = W * 0.05
    canvas_inner = W - 2 * canvas_pad
    scale = canvas_inner / (2 * half_extent)

    # Compose the transform: canvas center → scaled, y-flipped → root_xy at origin.
    cr.save()
    cr.translate(W / 2.0, H / 2.0)
    cr.scale(scale, -scale)
    cr.translate(-root_xy[0], -root_xy[1])

    # Reference rings — drawn FIRST (under the graph).  Line width set
    # in *data* coords; divide by scale to get the desired pixel width.
    cr.set_source_rgba(0.4, 0.4, 0.4, 1.0)
    cr.set_line_width(0.6 / scale)
    cr.set_dash([5.0 / scale, 5.0 / scale])
    for r_ref in (r_min, r_mean, r_max):
        cr.arc(root_xy[0], root_xy[1], r_ref, 0.0, 2.0 * math.pi)
        cr.stroke()
    cr.set_dash([])

    # Graph drawing — uses the current CTM, so positions in `pos` will
    # be applied through our transform.  Vertex size is in *device*
    # (pixel) units after the CTM.
    # ``vertex_size`` here is in data (user-CTM) coords so it scales with
    # the canvas.  ``edge_pen_width`` values, on the other hand, look
    # device-pixel-sized after our scaled CTM unless we divide by scale.
    edge_pen_width_user = g.new_edge_property("double")
    for e in g.edges():
        edge_pen_width_user[e] = float(edge_pen_width[e]) / scale

    # ---- Same-probe-only filter; let graph-tool draw the curves -------
    # Bypassing the manual Bezier emission (it was producing wrong
    # curves due to a frame-mapping issue).  Use gtd.cairo_draw
    # directly on a GraphView filtered to keep only intra-probe edges
    # (~500 of 6903) — that's the structurally meaningful subset and
    # the file weight is naturally small.
    same_probe_filter = g.new_edge_property("bool")
    n_same = 0
    for e, (i, j, _k) in zip(g.edges(), edge_meta):
        keep = probes[i] == probes[j]
        same_probe_filter[e] = keep
        if keep:
            n_same += 1
    print(f"  → drawing {n_same} same-probe edges (cross-probe omitted for now)")

    g_view = gt.GraphView(g, efilt=same_probe_filter)
    gtd.cairo_draw(
        g_view, pos, cr,
        edge_control_points=cts,
        edge_color=edge_color,
        edge_pen_width=edge_pen_width_user,
        eorder=eorder,
        vertex_fill_color=vcolor,
        vertex_color=[0.0, 0.0, 0.0, 1.0],
        vertex_pen_width=args.vertex_pen_width,
        vertex_size=args.vertex_size,
    )
    cr.restore()

    # Ring labels in device (pixel) coordinates so font size is sane.
    cr.select_font_face("sans-serif")
    cr.set_font_size(8.0)
    cr.set_source_rgba(0.4, 0.4, 0.4, 1.0)
    for r_ref, label in ((r_min, "min"), (r_mean, "mean"), (r_max, "max")):
        # Map (root_xy + r_ref, root_xy[1]) through CTM:
        # x_dev = W/2 + r_ref * scale, y_dev = H/2 (no flip needed for label y)
        x_dev = W / 2.0 + r_ref * scale + 4.0
        y_dev = H / 2.0
        cr.move_to(x_dev, y_dev)
        cr.show_text(f"{label} ({r_ref:.2f})")

    surface.finish()
    size_kb = out.stat().st_size / 1024
    print(f"\nwrote {out}")
    print(f"size: {size_kb:.0f} KB")


if __name__ == "__main__":
    main()
