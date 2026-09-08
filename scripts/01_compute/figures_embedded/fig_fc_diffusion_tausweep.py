#!/usr/bin/env python
"""Brain-FC diffusion-τ sweep — row-1 treatment on the mst@0.20 SPARSIFIED graph.

The complete |ImCoh| graph is degenerate under LRG (D=1/K ≡ D=1/A, audit_174), so
this figure runs on the **mst@0.20 backbone** (maximum spanning tree ∪ strongest 20%
of edges) — the connected, spanning, cycle-rich regime LRG is built for. THREE views
of that ONE object, laid out 1-on-top-2-below to leave room for a second slide image:
  (1) TOP  — the fixed backbone connectome as a graph on the MDS communication bubble.
             Edge WIDTH/OPACITY ∝ the |ImCoh| weight (set ONCE, τ-independent);
             edge COLOUR = ρ(τ) flux through each fixed link (same cmap+clim as the
             matrix). Nodes recolour by their current community.
  (2) LEFT — the LRG density operator ρ(τ)=e^{-τL̂}/Tr[e^{-τL̂}] as a matrix, reordered
             into this τ's leaf order, tinted by the community blocks (nested
             communities reprojected on the propagator).
  (3) RIGHT— the τ-morphing UPGMA tree with a FIXED communication-distance cut line at
             h0=1.012·n. The y-axis is FIXED (log) so the cut is a stable reference:
             the tree deflates DOWN through it as τ grows, communities merge many→one,
             and the nodes recolour in step.

τ sweeps the canonical fine scale PAST the characteristic scale into the degenerate
regime so the hierarchy is seen to dissolve:
    τ_min = 1/λ_max      (fine, the canonical LRG τ)
    τ*    = specific-heat C(τ) peak (marked on the τ bar)
    τ_end = first τ where the tree is essentially FLAT (relative merge-height spread
            < 3%): ρ(τ)→uniform 1/N, every pair equidistant, tree collapses to h≈n.

Node colour is GENEALOGICAL: each community takes the hue of the largest fine-τ atom
it contains, so a cluster keeps its colour while it keeps its members and, on a merge,
the union inherits the bigger atom's colour (small joins big) — no flicker.
Combinatorial L̂ throughout.

Outputs (data/outputs/figures/lrg_diffusion_zoom/):
    fig_fc_diffusion_tausweep.mp4             #f6f6f8 τ-sweep animation (talk asset)
    fig_fc_diffusion_tausweep_snapshots.pdf   transparent 3-stage contact sheet

Run from the repo root (inside lapbrain), or set LRGEEGFC_DATA_ROOT:
    python scripts/01_compute/figures_embedded/fig_fc_diffusion_tausweep.py
    python .../fig_fc_diffusion_tausweep.py --snap    # snapshots only (fast)
"""
from __future__ import annotations

import os
import sys
import time
import shutil
import tempfile
import subprocess
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc_context
from matplotlib.collections import LineCollection
from matplotlib.ticker import MaxNLocator, NullLocator
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import squareform
from scipy.interpolate import BSpline
import imageio_ffmpeg
plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_lrg_diffusion_zoom_demo import (          # noqa: E402
    laplacian_eig, tau_at, EDGE_CMAP, REAL)
from fig_lrg_multiscale_combined_demo import _dist_tau  # noqa: E402
from fig_lrg_dendrogram_cut_demo import (          # noqa: E402
    cophenetic_linkage, dendrogram_segments, curved_paths,
    relabel_by_leaforder, pal, GREY)
from lrg_eegfc.utils.metrics.functional_tree_distance import tau_star_from_C  # noqa: E402
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction  # noqa: E402
from lrg_eegfc.config.paths import FIGURES_ROOT      # noqa: E402
from lrg_eegfc.workflow.fc import load_fc_matrix      # noqa: E402
from lrg_eegfc.visuals.styles import use_lrg_style    # noqa: E402

NFRAMES = 240                   # long, smooth τ sweep (slower playback, no skipped frames)
FPS = 20                        # nice framerate; slowness comes from NFRAMES, not a low fps
HOLD = 36                       # frames held on the (spanning) end
DPI = 120                       # animation-frame DPI (talk overlay; ↓ = smaller/smoother file)
FRAC = 0.20                     # mst@0.20 backbone = MST ∪ strongest 20% of edges
H0_FACTOR = 1.012               # comm-distance cut = H0_FACTOR·n (n = fully-mixed D)
SPAN_MIN_DEC = 0.5              # stop the sweep while the tree still spans ≥ this many log-decades
VMAX_SMOOTH = 5                 # moving-average window on the per-frame colour scale (kills flicker)
OV_ALPHA = 0.40                 # matrix community-tint alpha (heat shows through)
ELW_MIN = 0.12                  # floor on edge linewidth so weak |ImCoh| links stay visible
BG = "#f6f6f8"                  # opaque-fallback background (matches a light slide)
GRID = "#9aa3af"                # spine / track grey
CUT_C = "#8b949e"               # h0 cut-line colour
TRANSP = {"figure.facecolor": "none", "axes.facecolor": "none",
          "savefig.facecolor": "none"}                   # transparent frames (real alpha)
SYS_FFMPEG = "/usr/bin/ffmpeg"  # full ffmpeg (VP9 alpha); bundled imageio one can't do alpha


def entropy_C_from_eigs(w, steps=600, t1=-2, t2=5):
    """C(τ) on the combinatorial spectrum w — mirrors the exact formula of
    lrgsglib.utils.lrg.infocomm.entropy (kept eigenvalue-based so C(τ) sits on
    the SAME L̂ as the propagator we draw). Returns (τ grid, C = -dS/dlogτ)."""
    w = np.where(np.abs(w) > 1e-12, w, 0.0)
    N = len(w)
    t = np.logspace(t1, t2, steps)
    S = np.empty(steps)
    with np.errstate(divide="ignore", invalid="ignore"):
        for i, tau in enumerate(t):
            r = np.exp(-tau * w)
            rho = r / np.nansum(r)
            xlogx = np.where(rho > 0, rho * np.log(rho), 0.0)
            S[i] = -np.nansum(xlogx) / np.log(N)               # normalised, 1→0
    C = np.log(N) * np.diff(1.0 - S) / np.diff(np.log(t))     # specific heat ≥0
    return t, C


def tau_span_floor(w, V, tau_star, span_min_dec=SPAN_MIN_DEC, reach=40.0, probes=80):
    """First τ ≥ τ* where the UPGMA tree's LOG-height spread drops below
    `span_min_dec` decades, i.e. log10(Zh.max) − log10(Zh.min) < span_min_dec.
    We STOP the sweep here instead of running to full degeneracy: as τ→∞ every
    D_ij → n so the merge heights converge to one value and the tree flattens to a
    razor-thin band (no y-scaling can make a zero-spread tree fill the panel). Ending
    while the tree still spans ≥ span_min_dec decades keeps the dendrogram full-height
    to the last frame, and h0 (=1.012·n) stays INSIDE [Zh.min, Zh.max] there (the tree
    top is still above n), so the cut line never leaves the panel."""
    grid = tau_star * np.logspace(0.0, np.log10(reach), probes)
    for t in grid:
        Zh = linkage(squareform(_dist_tau(w, V, t), checks=False), "average")[:, 2]
        top, bot = float(Zh.max()), float(Zh.min())
        if top <= 0 or bot <= 0 or (np.log10(top) - np.log10(bot)) < span_min_dec:
            return float(t)
    return float(grid[-1])


HEB_BETA = 0.86                 # Holten bundling strength (0 = straight chords, 1 = tight trunks)
HEB_SAMP = 40                   # points sampled per bundled edge curve


def _tree_topology(Z, n):
    """Parents, per-node leaf members, and depth for a scipy linkage Z (n leaves)."""
    parent = np.full(2 * n - 1, -1, int)
    members = [None] * (2 * n - 1)
    for i in range(n):
        members[i] = [i]
    for k in range(n - 1):
        a, b, p = int(Z[k, 0]), int(Z[k, 1]), n + k
        parent[a] = parent[b] = p
        members[p] = members[a] + members[b]
    depth = np.zeros(2 * n - 1, int)
    stack = [(2 * n - 2, 0)]
    while stack:
        nd, dd = stack.pop()
        depth[nd] = dd
        if nd >= n:
            stack.append((int(Z[nd - n, 0]), dd + 1))
            stack.append((int(Z[nd - n, 1]), dd + 1))
    return parent, members, depth


def _lca_path(i, j, parent):
    """Tree-node path i → LCA(i,j) → j (list of node ids), the bundling route."""
    up_i = [i]
    while parent[up_i[-1]] != -1:
        up_i.append(parent[up_i[-1]])
    idx_i = {nd: k for k, nd in enumerate(up_i)}
    up_j = [j]
    while up_j[-1] not in idx_i:
        up_j.append(parent[up_j[-1]])
    lca = up_j[-1]
    return up_i[:idx_i[lca] + 1] + up_j[:-1][::-1]


def _bundle_curve(Q, beta, n_samp):
    """Holten-bundle a control polyline Q (m×2) toward its straight chord, then draw
    a clamped cubic B-spline that USES Q as CONTROL points (not interpolation targets).
    A B-spline curve stays inside the convex hull of its control points — and every Q
    lies within the unit disk (leaves at r=1, internal nodes at r≤1) — so the arc can
    NEVER bulge outside the ring. (An *interpolating* spline overshoots at the sharp
    leaf turns and pokes outside the circle — the bug this fixes.) Clamped knots make
    the curve pass exactly through the two ring endpoints Q[0], Q[-1]."""
    m = len(Q)
    if m >= 3:
        t = np.linspace(0.0, 1.0, m)[:, None]
        Q = beta * Q + (1.0 - beta) * (Q[0] + t * (Q[-1] - Q[0]))
    keep = np.concatenate([[True], (np.abs(np.diff(Q, axis=0)).sum(1) > 1e-9)])
    Q = Q[keep]
    m = len(Q)
    u = np.linspace(0.0, 1.0, n_samp)
    if m < 2:
        return np.repeat(Q, n_samp, axis=0)
    if m == 2:
        return Q[0] + u[:, None] * (Q[1] - Q[0])
    k = min(3, m - 1)                                     # spline degree
    knots = np.concatenate((np.zeros(k), np.linspace(0.0, 1.0, m - k + 1), np.ones(k)))
    sx = BSpline(knots, Q[:, 0], k)                       # Q = control points ⇒ hull-bounded
    sy = BSpline(knots, Q[:, 1], k)
    return np.column_stack([sx(u), sy(u)])


def heb_bundled_layout(Z, n, ei, ej, beta=HEB_BETA, n_samp=HEB_SAMP):
    """Circular hierarchical-edge-bundling layout (Holten 2006) on the LRG tree Z.
    Leaves sit evenly on a ring in the tree's OWN leaf order (communities contiguous,
    trunks tightest); internal nodes sit at radius ∝ depth and the circular-mean angle
    of their leaves, so every backbone edge routes i → common-ancestor trunk → j as a
    smooth bundled arc. Reproduces `gen_pipeline_sequence.panel4_backbone`'s look in
    pure numpy so the arcs are recolourable per τ (LineCollection). Geometry is FIXED
    (τ-independent structure); only the edge FLOW colour evolves.

    Returns (P_leaf n×2 on the unit ring, list of n_samp×2 bundled arcs aligned to ei/ej)."""
    order = [int(x) for x in leaves_list(Z)]
    posmap = np.empty(n, int)
    posmap[order] = np.arange(n)
    parent, members, depth = _tree_topology(Z, n)
    maxd = max(int(depth.max()), 1)
    leaf_ang = 2.0 * np.pi * posmap / n
    node_pos = np.zeros((2 * n - 1, 2))
    for v in range(2 * n - 1):
        a = leaf_ang[members[v]]
        ang = np.arctan2(np.sin(a).sum(), np.cos(a).sum())    # circular mean
        rad = 1.0 if v < n else depth[v] / maxd               # leaves outer, root at centre
        node_pos[v] = (rad * np.cos(ang), rad * np.sin(ang))
    segs = [_bundle_curve(node_pos[_lca_path(int(i), int(j), parent)], beta, n_samp)
            for i, j in zip(ei, ej)]
    return node_pos[:n], segs


def prepare():
    A = np.asarray(load_fc_matrix(REAL["patient"], REAL["phase"], REAL["band"],
                                  "imcoh_abs"), float).copy()
    np.fill_diagonal(A, 0.0)
    A = np.clip(A, 0.0, None)
    n = A.shape[0]
    # SPARSIFY: the complete graph is LRG-degenerate (audit_174); run on the
    # mst@0.20 backbone — connected, spanning, cycle-rich — the regime LRG needs.
    A = mst_union_top_fraction(A, FRAC)
    dens = float((np.triu(A, 1) > 0).sum()) / (n * (n - 1) / 2.0)
    w, V = laplacian_eig(A)
    lam_max = float(w.max())
    tau_min = 1.0 / lam_max
    tC, C = entropy_C_from_eigs(w)
    tau_star, i_star, interior = tau_star_from_C(tC, C)
    tau_end = tau_span_floor(w, V, tau_star)              # stop while the tree still spans
    s_star = float(np.log(tau_star / tau_min) / np.log(tau_end / tau_min))
    print(f"  {REAL['patient']} β rest_pre: n={n}  mst@{FRAC:.2f} backbone "
          f"density={dens:.3f}  λ_max={lam_max:.3f}", flush=True)
    print(f"  τ_min = 1/λ_max = {tau_min:.4g}   ({tau_min * lam_max:.2f}/λ_max)",
          flush=True)
    print(f"  τ*    = {tau_star:.4g}   ({tau_star * lam_max:.2f}/λ_max)   "
          f"[C-peak interior={interior}]  → bar fraction {s_star:.2f}", flush=True)
    print(f"  τ_end = {tau_end:.4g}   ({tau_end * lam_max:.2f}/λ_max)   "
          f"[tree flat / degenerate]", flush=True)

    h0 = H0_FACTOR * n                                    # fixed communication-distance cut
    print(f"  h0 = {H0_FACTOR:.3f}·n = {h0:.1f}   (fixed comm-distance threshold; "
          f"clusters merge as diffusion shrinks D below it)", flush=True)
    # GENEALOGY colouring. Atoms = the finest-τ communities (h0 cut on the τ=1/λ_max
    # tree), fixed ONCE and size-ranked (rank 0 = the biggest atom, the global core).
    # A community at ANY τ is painted by the LARGEST atom it contains (min atom-rank),
    # so a cluster keeps its colour while it keeps its members and, on a merge, the
    # union inherits the bigger atom's colour — small joins big, with no flicker.
    Z_ref = cophenetic_linkage(w, V)                      # UPGMA at τ = 1/λ_max
    atom = fcluster(Z_ref, t=h0, criterion="distance")    # finest-τ communities = atoms
    labs, sizes = np.unique(atom, return_counts=True)
    order = np.argsort(-sizes, kind="stable")             # biggest atom first → rank 0
    atom_rank = {int(labs[o]): r for r, o in enumerate(order)}
    node_atom_rank = np.array([atom_rank[int(a)] for a in atom], int)   # per node, 0 = biggest
    atom_color = np.array([pal(r) for r in range(len(labs))])          # rank → tab20 hue (cycled)
    print(f"  {len(labs)} fine-τ atoms   (biggest = {int(sizes.max())} nodes → the "
          f"hue everything coalesces to)", flush=True)

    ei, ej = np.nonzero(np.triu(A, 1) > 0)                # BACKBONE edges only (sparse)
    P, net_segs = heb_bundled_layout(Z_ref, n, ei, ej)    # circular HEB: ring + bundled arcs
    print(f"  circular HEB layout: {n} beads on the ring, {len(ei)} edges bundled "
          f"through the τ_min hierarchy (β={HEB_BETA})", flush=True)
    # FIXED edge weights = the backbone |ImCoh| connectome (the INPUT graph). Width and
    # opacity ∝ weight, computed ONCE — they never change with τ. The propagator ρ(τ)
    # (the diffusion OUTPUT) lives in the matrix + as edge COLOUR; never as edge width.
    wn = A[ei, ej] / (float(A[ei, ej].max()) + 1e-12)
    elw = np.maximum(ELW_MIN, 3.1 * wn ** 1.8)          # ∝ weight, floored so weak links stay visible
    ealpha = 0.10 + 0.66 * wn ** 1.25                     # opacity ∝ weight (sparse ⇒ higher floor)
    taus = tau_at(tau_min, tau_end, np.linspace(0.0, 1.0, NFRAMES))
    taus_static = tau_at(tau_min, tau_end, np.array([0.0, s_star, 1.0]))
    # SMOOTHED colour scale: the raw per-frame p99(ρ) jitters frame-to-frame and makes
    # the flux colours (edges + matrix) flicker ("a scatti"). Precompute vmax(τ) and
    # moving-average it in log space so the colour ramp evolves smoothly.
    vmax_raw = np.array([max(float(np.percentile(
        density_matrix(w, V, t)[ei, ej], 99.0)), 1e-12) for t in taus])
    kw = max(1, VMAX_SMOOTH | 1)                          # odd window
    logv = np.log(vmax_raw)
    logv = np.convolve(np.pad(logv, kw // 2, mode="edge"), np.ones(kw) / kw, "valid")
    vmax_sched = np.exp(logv[:NFRAMES])
    return dict(w=w, V=V, n=n, P=P, ei=ei, ej=ej, net_segs=net_segs, h0=h0,
                node_atom_rank=node_atom_rank, atom_color=atom_color,
                elw=elw, ealpha=ealpha, vmax_sched=vmax_sched,
                taus=taus, taus_static=taus_static, s_star=s_star,
                lam_max=lam_max, tau_min=tau_min, tau_star=tau_star, tau_end=tau_end)


# ----------------------------------------------------------------------------
# per-τ state (one hierarchy drives edges + matrix + blocks + tree together)
# ----------------------------------------------------------------------------
def density_matrix(w, V, tau):
    """ρ(τ) = e^{-τL̂} / Tr[e^{-τL̂}] — the Villegas LRG density operator (the
    heat-kernel propagator NORMALISED by its trace Z(τ)=Σ e^{-τλ}). Off-diagonal
    kept (information flow); entrywise ≥ 0. The trace is a scalar, so ρ(τ) and the
    bare kernel e^{-τL̂} are identical under the per-frame contrast scaling."""
    ev = np.exp(-tau * w)
    K = (V * ev[None, :]) @ V.T
    rho = K / ev.sum()
    np.fill_diagonal(rho, 0.0)
    return np.clip(rho, 0.0, None)


def community_colors(comm, node_atom_rank, atom_color):
    """Genealogy colouring — per-node RGBA. Each current-τ community is painted by
    the LARGEST fine-τ atom it contains (smallest atom-rank ⇒ biggest atom; rank 0 =
    the global core). The colour is a deterministic function of the member SET, so a
    cluster that keeps the same members keeps its colour across every τ, and on a
    merge the union inherits the bigger atom's colour (small joins big). No mode-of-
    bins flicker, no reference-giant palette collapse."""
    ncol = np.zeros((len(comm), 4))
    for c in np.unique(comm):
        m = comm == c
        ncol[m] = atom_color[int(node_atom_rank[m].min())]
    return ncol


def frame_state(d, tau, vmax=None):
    """Everything that EVOLVES with τ: the ρ(τ) matrix reordered into THIS τ's leaf
    order, the per-edge ρ(τ) flow values, the community tint overlay, per-node
    community colours, and the morphing dendrogram segments. Edge SHAPE + WIDTH are
    FIXED (the |ImCoh| structure set once in prepare()); only the edge COLOUR evolves
    — = ρ(τ) flux through each fixed link (same cmap + clim as the matrix panel). The
    dendrogram y-axis is ADAPTIVE per frame (returns `ylim`) so the tree always fits.
    `vmax` (the colour clim) is passed from the smoothed schedule to avoid flicker;
    if None it falls back to this frame's raw p99(ρ)."""
    n = d["n"]
    rho = density_matrix(d["w"], d["V"], tau)
    if vmax is None:
        vmax = max(float(np.percentile(rho[d["ei"], d["ej"]], 99.0)), 1e-12)  # matrix clim
    # morphing tree at this τ, scipy's proper order → no crossing segments
    Z = linkage(squareform(_dist_tau(d["w"], d["V"], tau), checks=False), "average")
    leaves_f, segs, dh_f, drep_f, _ = dendrogram_segments(Z)
    # ADAPTIVE (log) y-range that keeps THIS τ's min & max merge heights within view
    # (the fixed axis clipped the fine tree and left the flat tree at half-panel). Feet
    # clamp to a floor just below the lowest split (short stems). The h0 cut line is
    # kept in frame; it now moves with the rescaling axis instead of being anchored.
    Zh = Z[:, 2]
    hmin, hmax = float(Zh.min()), float(Zh.max())
    floor = hmin / 1.18
    ylim = (floor, hmax * 1.18)                           # tight: tree fills the panel; h0 sits inside
    for s in segs:
        s[:, 1] = np.maximum(s[:, 1], floor)
    # communities: a FIXED communication-distance threshold h0. Diffusion shrinks
    # every D_ij = 1/K_ij(τ) as τ grows, so clusters merge on their own — MANY small
    # communities at fine τ → ONE at degeneracy — continuously, with NO imposed
    # k-schedule (a height-gap fails: one core module + a rim ⇒ gap/Ψ give one blob).
    comm = relabel_by_leaforder(fcluster(Z, t=d["h0"], criterion="distance"), leaves_f)
    ncol = community_colors(comm, d["node_atom_rank"], d["atom_color"])           # genealogy hues
    below = dh_f < d["h0"]
    tcol = np.tile(GREY, (len(segs), 1))
    tcol[below] = ncol[leaves_f[drep_f[below]]]           # colour merges under the cut
    # matrix reordered into THIS τ's leaves + community tint (both track the tree)
    rho_ord = rho[np.ix_(leaves_f, leaves_f)]
    lab_ord = comm[leaves_f]
    same = lab_ord[:, None] == lab_ord[None, :]
    ov = np.zeros((n, n, 4))
    ov[..., :3] = ncol[leaves_f][:, :3][:, None, :]
    ov[..., 3] = np.where(same, OV_ALPHA, 0.0)
    return dict(vmax=vmax, rho_ord=rho_ord, overlay=ov,
                edge_val=rho[d["ei"], d["ej"]],
                ncol=ncol, segs=segs, tcol=tcol, ylim=ylim)


def _net_axes(ax, P, pad=0.1):
    ax.set_xlim(P[:, 0].min() - pad, P[:, 0].max() + pad)
    ax.set_ylim(P[:, 1].min() - pad, P[:, 1].max() + pad)
    ax.set_aspect("equal"); ax.axis("off")


def edge_flow_colors(edge_val, vmax, ealpha, cmap):
    """Per-edge RGBA for the network. The propagator FLOW ρ(τ)_ij is mapped through
    the SAME cmap + clim [0, vmax] as the matrix panel, so an edge's colour equals its
    cell in the ρ(τ) matrix — 'what flux passes through this link at time τ'. Shape,
    width and opacity stay the FIXED |ImCoh| structure (from prepare()); only the
    COLOUR carries the time-varying flow."""
    c = cmap(np.clip(edge_val / vmax, 0.0, 1.0))
    c[:, 3] = ealpha
    return c


def _style_tree_axis(ax, d):
    """LOG dendrogram axis (static parts) + the h0 cut line. The y-limits are set by
    the CALLER per frame (adaptive fit); the cut line is drawn in data coords so it
    tracks the rescaling axis."""
    ax.set_yscale("log")
    ax.set_xlim(0, 10 * d["n"])
    ax.set_xticks([])
    ax.yaxis.set_major_locator(NullLocator())             # NO tick text — even sub-decade
    ax.yaxis.set_minor_locator(NullLocator())             # (log auto-labels minor ticks < 1 decade)
    ax.axhline(d["h0"], color=CUT_C, lw=1.1, ls=(0, (5, 4)), zorder=1)  # comm-distance cut
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)


# ----------------------------------------------------------------------------
# animation (TRANSPARENT frames → WebM-alpha + #f6f6f8 MP4 · no text · single ROW)
# ----------------------------------------------------------------------------
# figure 13.0 × 4.5 in; the three panels are square (width_in == height_in) in one row.
NET = [0.012, 0.066, 0.300, 0.867]      # circular HEB network (square, 3.90")
MAT = [0.335, 0.066, 0.300, 0.867]      # ρ(τ) propagator matrix (square)
DEN = [0.658, 0.066, 0.300, 0.867]      # τ-morphing UPGMA tree (square)
BAR = [0.968, 0.066, 0.020, 0.867]      # slim τ progress bar at the right edge


def _encode(frames_dir: Path, w: int, h: int, out_webm: Path, out_mp4: Path):
    """PNG sequence → (1) TRANSPARENT WebM (VP9 yuva420p, alpha_mode=1 — needs the full
    system ffmpeg; the bundled imageio build silently drops alpha) and (2) an opaque
    #f6f6f8 MP4 fallback (H.264, alpha flattened onto the slide colour). Matplotlib's
    Animation.save() hard-composites to white, so real alpha REQUIRES this manual route."""
    ff = SYS_FFMPEG if os.path.exists(SYS_FFMPEG) else imageio_ffmpeg.get_ffmpeg_exe()
    src = ["-framerate", str(FPS), "-i", str(frames_dir / "f%05d.png")]
    if os.path.exists(SYS_FFMPEG):                        # transparent WebM (alpha)
        subprocess.run([ff, "-hide_banner", "-loglevel", "error", "-y", *src,
                        "-c:v", "libvpx-vp9", "-pix_fmt", "yuva420p", "-b:v", "0",
                        "-crf", "30", "-auto-alt-ref", "0", "-row-mt", "1",
                        "-metadata:s:v:0", "alpha_mode=1", str(out_webm)], check=True)
        print(f"done: {out_webm.name}  ({out_webm.stat().st_size / 1e6:.1f} MB, "
              f"transparent VP9)", flush=True)
    else:
        print("  ⚠ no system ffmpeg — skipping transparent WebM (bundled build "
              "can't encode alpha)", flush=True)
    # opaque #f6f6f8 fallback: flatten alpha onto the slide colour, even dims for yuv420p
    subprocess.run([ff, "-hide_banner", "-loglevel", "error", "-y", *src,
                    "-filter_complex",
                    f"color=c=0xf6f6f8:s={w}x{h}:r={FPS}[bg];"
                    "[bg][0:v]overlay=shortest=1,format=yuv420p",
                    "-c:v", "libx264", "-crf", "20", "-preset", "slow",
                    "-movflags", "+faststart", str(out_mp4)], check=True)
    print(f"done: {out_mp4.name}  ({out_mp4.stat().st_size / 1e6:.1f} MB, opaque "
          f"#f6f6f8)", flush=True)


def make_animation(d, out_webm: Path, out_mp4: Path):
    n = d["n"]
    ecmap = plt.get_cmap("turbo")
    tmp = Path(tempfile.mkdtemp(prefix="tausweep_"))
    with rc_context(TRANSP):
        fig = plt.figure(figsize=(13.0, 4.5), dpi=DPI)
        ax_net = fig.add_axes(NET)
        ax_bar = fig.add_axes(BAR)
        ax_mat = fig.add_axes(MAT)
        ax_den = fig.add_axes(DEN)
        st_init = frame_state(d, d["taus"][0], vmax=d["vmax_sched"][0])
        ecol0 = edge_flow_colors(st_init["edge_val"], st_init["vmax"], d["ealpha"], ecmap)
        pc = LineCollection(d["net_segs"], colors=ecol0, linewidths=d["elw"],
                            capstyle="round", rasterized=True)   # rasterise the heavy graph
        ax_net.add_collection(pc)
        sc = ax_net.scatter(d["P"][:, 0], d["P"][:, 1], s=46,
                            c=st_init["ncol"], edgecolors="#0b0e17",
                            linewidths=0.6, zorder=3, rasterized=True)
        _net_axes(ax_net, d["P"])
        im = ax_mat.imshow(np.zeros((n, n)), cmap=ecmap, vmin=0, vmax=1,
                           interpolation="nearest")
        ov_im = ax_mat.imshow(np.zeros((n, n, 4)), interpolation="nearest", zorder=3)
        ax_mat.set_xticks([]); ax_mat.set_yticks([])
        for s in ax_mat.spines.values():
            s.set_color(GRID)
        dlc = LineCollection(st_init["segs"], linewidths=1.2)
        ax_den.add_collection(dlc)
        _style_tree_axis(ax_den, d)                       # log axis + h0 cut line
        ax_den.set_ylim(*st_init["ylim"])
        ax_bar.set_xlim(0, 1); ax_bar.set_ylim(0, 1); ax_bar.axis("off")
        ax_bar.plot([0.5, 0.5], [0.16, 0.84], color=GRID, lw=3, solid_capstyle="round")
        y_star = 0.16 + 0.68 * d["s_star"]                # τ* mark on the bar
        ax_bar.plot([0.3, 0.7], [y_star, y_star], color=CUT_C, lw=1.1, alpha=0.85, zorder=2)
        marker, = ax_bar.plot([0.5], [0.16], "o", color="#e08d00", ms=11, zorder=3)
        t0 = time.time()

        total = NFRAMES + HOLD
        for frame in range(total):
            f = min(frame, NFRAMES - 1)
            s = f / (NFRAMES - 1)
            tau = d["taus"][f]
            st = frame_state(d, tau, vmax=d["vmax_sched"][f])   # smoothed clim (no flicker)
            pc.set_color(edge_flow_colors(st["edge_val"], st["vmax"], d["ealpha"], ecmap))
            sc.set_facecolor(st["ncol"])                  # nodes recolour by community
            im.set_data(st["rho_ord"]); im.set_clim(0, st["vmax"])
            ov_im.set_data(st["overlay"])
            dlc.set_segments(st["segs"]); dlc.set_color(st["tcol"])
            ax_den.set_ylim(*st["ylim"])                  # adaptive: tree always fits
            marker.set_data([0.5], [0.16 + 0.68 * s])
            fig.savefig(tmp / f"f{frame:05d}.png", transparent=True, dpi=DPI)
            if frame % 24 == 0 or frame == total - 1:
                print(f"  [frame {frame + 1:3d}/{total}] s={s:4.2f} τ={tau:.4g} "
                      f"elapsed={time.time() - t0:5.1f}s", flush=True)
        w_px, h_px = (fig.canvas.get_width_height())
        plt.close(fig)
    out_webm.parent.mkdir(parents=True, exist_ok=True)
    _encode(tmp, w_px, h_px, out_webm, out_mp4)
    shutil.rmtree(tmp, ignore_errors=True)


# ----------------------------------------------------------------------------
# snapshots (transparent PDF) — 3 rows × 3 τ stages
# ----------------------------------------------------------------------------
def make_snapshots(d, out_pdf: Path):
    ecmap = plt.get_cmap("turbo")
    stages = ("fine  τ=1/λ_max", "τ*  (C-peak)", "degenerate  (flat tree)")
    rc = {"figure.facecolor": "none", "savefig.facecolor": "none",
          "axes.facecolor": "none", "text.color": "#1f2937"}
    with rc_context(rc):
        fig, axes = plt.subplots(3, 3, figsize=(10.0, 10.4))
        fig.subplots_adjust(left=0.065, right=0.985, top=0.945, bottom=0.02,
                            wspace=0.13, hspace=0.14)
        for ci, tau in enumerate(d["taus_static"]):
            st = frame_state(d, tau)
            axn = axes[0, ci]
            eflow = edge_flow_colors(st["edge_val"], st["vmax"], d["ealpha"], ecmap)
            axn.add_collection(LineCollection(d["net_segs"], colors=eflow,
                                              linewidths=d["elw"], capstyle="round"))
            axn.scatter(d["P"][:, 0], d["P"][:, 1], s=30,
                        c=st["ncol"], edgecolors=(0, 0, 0, 0.35), linewidths=0.4,
                        zorder=3)
            _net_axes(axn, d["P"])
            axn.set_title(f"{stages[ci]}", fontsize=11, color="#1f2937")
            axm = axes[1, ci]
            axm.imshow(st["rho_ord"], cmap=ecmap, vmin=0, vmax=st["vmax"],
                       interpolation="nearest")
            axm.imshow(st["overlay"], interpolation="nearest", zorder=3)
            axm.set_xticks([]); axm.set_yticks([])
            for sp in axm.spines.values():
                sp.set_color(GRID)
            axd = axes[2, ci]
            axd.add_collection(LineCollection(st["segs"], colors=st["tcol"], linewidths=1.0))
            _style_tree_axis(axd, d)                       # log axis + h0 cut line
            axd.set_ylim(*st["ylim"])                      # adaptive per-frame fit
        for r, lab in enumerate(("structure · ρ(τ) flow · communities",
                                 "ρ(τ) propagator · community blocks",
                                 "τ-morphing UPGMA tree · fixed h0 cut")):
            axes[r, 0].set_ylabel(lab, fontsize=10, color="#1f2937", labelpad=8)
            axes[r, 0].axis("on"); axes[r, 0].set_xticks([]); axes[r, 0].set_yticks([])
            if r != 2:
                for sp in axes[r, 0].spines.values():
                    sp.set_visible(False)
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf, transparent=True)
        plt.close(fig)
    print(f"done: {out_pdf.name}", flush=True)


def main():
    use_lrg_style()
    out_dir = Path(FIGURES_ROOT) / "lrg_diffusion_zoom"
    d = prepare()
    make_snapshots(d, out_dir / "fig_fc_diffusion_tausweep_snapshots.pdf")
    if "--snap" not in sys.argv:
        make_animation(d, out_dir / "fig_fc_diffusion_tausweep.webm",
                       out_dir / "fig_fc_diffusion_tausweep.mp4")
    print(f"\nall assets in {out_dir}", flush=True)


if __name__ == "__main__":
    main()
