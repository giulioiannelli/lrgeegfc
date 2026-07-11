#!/usr/bin/env python
"""Brain-FC diffusion-τ sweep — row-1 treatment applied to the complete FC graph.

New standalone figure. THREE distinct views of the same object, NOT the same
quantity thrice: (1) the fixed |ImCoh| connectome as a graph — edge width/opacity
∝ the actual FC weight, never changing with τ; (2) the LRG density operator
ρ(τ)=e^{-τL̂}/Tr[e^{-τL̂}] (the diffusion OUTPUT) as a matrix, swept over τ; (3) a
τ-morphing UPGMA tree. The edges are the INPUT weights; ρ(τ) belongs ONLY in the
matrix — do not restate the propagator on the edges. Nodes recolour by the current
community, so the τ-evolution shows as communities coalescing on the fixed graph.

τ sweeps from the canonical fine scale PAST the characteristic scale into the
degenerate regime, so the hierarchy is seen to dissolve:
    τ_min = 1/λ_max      (fine, the canonical LRG τ)
    τ*    = specific-heat C(τ) peak (marked on the τ bar, functional_tree_distance)
    τ_end = first τ where the tree is essentially FLAT (relative merge-height
            spread < 10%): ρ(τ) → uniform 1/N, every pair equidistant, tree collapses.

Everything τ-dependent is tied to ONE per-τ hierarchy so the panels move together:
  · CURVED Bézier arcs (curved_paths) on the MDS communication bubble
    (lrg_communication_layout); straight chords look wrong on a complete graph. Edge
    width/opacity ∝ the FIXED |ImCoh| weight (computed once); edges never restate ρ(τ).
  · the ρ(τ) matrix uses row-1's EXACT [0, p99] normalisation, so it VISIBLY warms
    and restructures fine→coarse (a fixed [p5,p99] stretch renormalises every frame
    and reads as "stuck at fixed τ").
  · the dendrogram MORPHS: UPGMA(D(τ)) recomputed each frame, drawn in scipy's OWN
    proper leaf order for that τ (no crossing horizontal segments). Drawn in LINEAR
    scale with y-limits that TRACK the current tree (bottom = ε below the finest
    split → every split shown, no long leaf stems; top = the root) so it always
    fills the panel and is seen to compress, then flatten, as τ grows.
  · communities come from a FIXED communication-distance threshold h0 = 1.012·n
    (n = the fully-mixed distance, since D_ij → n as τ → ∞). Diffusion shrinks every
    D_ij(τ) as τ grows, so clusters merge on their own — MANY at fine τ → ONE at
    degeneracy — continuously, with NO imposed k-schedule. (A height-gap fails here:
    this tree is one core module + a rim, so gap/Ψ/target-k all give one giant blob,
    and a literal biggest-gap shatters to singletons; a fixed absolute cut lets the
    diffusion do the coalescing.) The matrix is REORDERED into that τ's leaf order so
    blocks stay contiguous, tinted at α=0.40 over ρ(τ) (heat shows through). Node colour
    is GENEALOGICAL: each community takes the hue of the largest fine-τ atom it
    contains, so a cluster keeps its colour as long as it keeps its members and,
    when clusters merge, the union inherits the bigger one's colour (small joins
    big) — no flicker.
Combinatorial L̂ throughout.

Outputs (data/outputs/figures/lrg_diffusion_zoom/):
    fig_fc_diffusion_tausweep.mp4             white τ-sweep animation (talk asset)
    fig_fc_diffusion_tausweep_snapshots.pdf   transparent 3-stage contact sheet

Run from the repo root (inside lapbrain), or set LRGEEGFC_DATA_ROOT:
    python scripts/01_compute/figures_embedded/fig_fc_diffusion_tausweep.py
    python .../fig_fc_diffusion_tausweep.py --snap    # snapshots only (fast)
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc_context
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.ticker import MaxNLocator
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import imageio_ffmpeg
plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_lrg_diffusion_zoom_demo import (          # noqa: E402
    laplacian_eig, tau_at, EDGE_CMAP, REAL)
from fig_lrg_multiscale_combined_demo import (     # noqa: E402
    lrg_communication_layout, _dist_tau)
from fig_lrg_dendrogram_cut_demo import (          # noqa: E402
    cophenetic_linkage, dendrogram_segments, curved_paths,
    relabel_by_leaforder, pal, GREY)
from lrg_eegfc.utils.metrics.functional_tree_distance import tau_star_from_C  # noqa: E402
from lrg_eegfc.config.paths import FIGURES_ROOT      # noqa: E402
from lrg_eegfc.workflow.fc import load_fc_matrix      # noqa: E402
from lrg_eegfc.visuals.styles import use_lrg_style    # noqa: E402

NFRAMES = 96
FPS = 12
HOLD = 14
DPI = 150                       # animation-frame DPI (talk overlay; ↓ = smaller file)
H0_FACTOR = 1.012               # fixed comm-distance cut = H0_FACTOR·n (n = fully-mixed D)
OV_ALPHA = 0.40                 # matrix community-tint alpha (heat shows through)
ELW_MIN = 0.12                  # floor on edge linewidth so weak |ImCoh| links stay visible
DARK = {"figure.facecolor": "#0b0e17", "axes.facecolor": "#0b0e17",
        "savefig.facecolor": "#0b0e17", "text.color": "#e6edf3",
        "axes.edgecolor": "#9aa3af", "axes.labelcolor": "#8b949e"}
WHITE = {"figure.facecolor": "#f6f6f8", "axes.facecolor": "#f6f6f8",
         "savefig.facecolor": "#f6f6f8"}


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


def tau_degenerate(w, V, tau_star, reach=40.0, probes=26, flat=0.10):
    """First τ ≥ τ* where the UPGMA tree is essentially FLAT — relative merge-
    height spread (Zh.max − Zh.min)/Zh.max < `flat`. Lets the sweep run PAST the
    characteristic scale into the degenerate regime, where ρ(τ) → uniform 1/N and
    every pair is equidistant so the whole hierarchy collapses to one height.
    Falls back to reach·τ* if the tree never flattens inside the probe window."""
    grid = tau_star * np.logspace(0.0, np.log10(reach), probes)
    for t in grid:
        Zh = linkage(squareform(_dist_tau(w, V, t), checks=False), "average")[:, 2]
        top = float(Zh.max())
        if top <= 0 or (top - float(Zh.min())) / top < flat:
            return float(t)
    return float(grid[-1])


def prepare():
    A = np.asarray(load_fc_matrix(REAL["patient"], REAL["phase"], REAL["band"],
                                  "imcoh_abs"), float).copy()
    np.fill_diagonal(A, 0.0)
    A = np.clip(A, 0.0, None)
    n = A.shape[0]
    w, V = laplacian_eig(A)
    lam_max = float(w.max())
    tau_min = 1.0 / lam_max
    tC, C = entropy_C_from_eigs(w)
    tau_star, i_star, interior = tau_star_from_C(tC, C)
    tau_end = tau_degenerate(w, V, tau_star, flat=0.03)   # sweep PAST τ* to flat tree
    s_star = float(np.log(tau_star / tau_min) / np.log(tau_end / tau_min))
    print(f"  Pat_05 β rest_pre: n={n}  λ_max={lam_max:.3f}", flush=True)
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

    P = lrg_communication_layout(w, V)                    # the bubble (not hairball)
    ei, ej = np.triu_indices(n, 1)                        # ALL edges (complete)
    paths = curved_paths(P, ei, ej, rad=0.26)             # curved Bézier arcs (bowed)
    # FIXED edge weights = the |ImCoh| connectome itself (the INPUT graph). Width and
    # opacity ∝ weight, computed ONCE — they never change with τ. The propagator ρ(τ)
    # (the diffusion OUTPUT) lives only in the matrix; edges must not restate it.
    wn = A[ei, ej] / (float(A[ei, ej].max()) + 1e-12)
    elw = np.maximum(ELW_MIN, 3.1 * wn ** 1.8)          # ∝ weight, floored so weak links stay visible
    ealpha = 0.04 + 0.66 * wn ** 1.25                     # opacity ∝ weight
    taus = tau_at(tau_min, tau_end, np.linspace(0.0, 1.0, NFRAMES))
    taus_static = tau_at(tau_min, tau_end, np.array([0.0, s_star, 1.0]))
    return dict(w=w, V=V, n=n, P=P, ei=ei, ej=ej, paths=paths, h0=h0,
                node_atom_rank=node_atom_rank, atom_color=atom_color,
                elw=elw, ealpha=ealpha,
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


def frame_state(d, tau):
    """Everything that EVOLVES with τ: the ρ(τ) matrix reordered into THIS τ's leaf
    order, the per-edge ρ(τ) flow values, the community tint overlay, per-node
    community colours, and the morphing dendrogram. Edge SHAPE + WIDTH are FIXED
    (the |ImCoh| structure set once in prepare()); only the edge COLOUR evolves —
    = ρ(τ) flux through each fixed link (same cmap + clim as the matrix panel)."""
    n = d["n"]
    rho = density_matrix(d["w"], d["V"], tau)
    vmax = max(float(np.percentile(rho[d["ei"], d["ej"]], 99.0)), 1e-12)  # matrix clim
    # morphing tree at this τ, scipy's proper order → no crossing segments
    Z = linkage(squareform(_dist_tau(d["w"], d["V"], tau), checks=False), "average")
    leaves_f, segs, dh_f, drep_f, _ = dendrogram_segments(Z)
    Zh = Z[:, 2]
    # LINEAR y-range that TRACKS this τ's tree so it always fills the panel; leaf
    # feet clamped just below the lowest split (ε) → every split shown, no long
    # stems. As τ coarsens the heights compress and the axis compresses with them
    # (the tree stays full-height instead of shrinking), collapsing flat at τ_end.
    top = float(Zh.max()); span = top - float(Zh.min())
    pad = max(span * 0.06, top * 1e-3)                    # nonzero even when flat
    ymin = float(Zh.min()) - pad
    for s in segs:
        s[:, 1] = np.maximum(s[:, 1], ymin)
    ylim = (ymin, top + max(span * 0.03, top * 5e-3))
    # communities: a FIXED communication-distance threshold h0. Diffusion shrinks
    # every D_ij = 1/K_ij(τ) as τ grows, so clusters merge on their own — MANY small
    # communities at fine τ → ONE at degeneracy — continuously, with NO imposed
    # k-schedule. (Reading a height-gap fails: this tree is one core module + a rim,
    # so gap/Ψ/target-k all give one giant blob; a literal biggest-gap shatters to
    # singletons. A fixed absolute cut instead lets the diffusion do the coalescing.)
    comm = relabel_by_leaforder(fcluster(Z, t=d["h0"], criterion="distance"), leaves_f)
    h_cut = d["h0"]
    ncol = community_colors(comm, d["node_atom_rank"], d["atom_color"])           # genealogy hues (largest atom wins)
    below = dh_f < h_cut
    tcol = np.tile(GREY, (len(segs), 1))
    tcol[below] = ncol[leaves_f[drep_f[below]]]
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


# ----------------------------------------------------------------------------
# animation (white MP4 · no text)
# ----------------------------------------------------------------------------
def make_animation(d, out_vid: Path):
    n = d["n"]
    ecmap = plt.get_cmap("turbo")
    with rc_context(WHITE):
        fig = plt.figure(figsize=(12.4, 4.52), dpi=90)
        gs = fig.add_gridspec(1, 4, width_ratios=[1.15, 0.92, 1.0, 0.16],
                              left=0.02, right=0.985, top=0.97, bottom=0.06,
                              wspace=0.14)
        ax_net, ax_mat, ax_den, ax_bar = (fig.add_subplot(gs[0, c]) for c in range(4))
        st_init = frame_state(d, d["taus"][0])            # initial flow colours
        ecol0 = edge_flow_colors(st_init["edge_val"], st_init["vmax"], d["ealpha"], ecmap)
        pc = PathCollection(d["paths"], facecolors="none", edgecolors=ecol0,
                            linewidths=d["elw"])          # FIXED shape+width; colour = ρ(τ) flow
        ax_net.add_collection(pc)
        sc = ax_net.scatter(d["P"][:, 0], d["P"][:, 1], s=82,
                            c=np.tile(GREY, (n, 1)), edgecolors="#0b0e17",
                            linewidths=0.6, zorder=3)
        _net_axes(ax_net, d["P"])
        im = ax_mat.imshow(np.zeros((n, n)), cmap=ecmap, vmin=0, vmax=1,
                           interpolation="nearest", animated=True)
        ov_im = ax_mat.imshow(np.zeros((n, n, 4)), interpolation="nearest",
                              zorder=3, animated=True)
        ax_mat.set_xticks([]); ax_mat.set_yticks([])
        for s in ax_mat.spines.values():
            s.set_color("#9aa3af")
        st0 = frame_state(d, d["taus"][0])                # valid initial geometry
        dlc = LineCollection(st0["segs"], linewidths=1.2)
        ax_den.add_collection(dlc)
        ax_den.set_xlim(0, 10 * n); ax_den.set_ylim(*st0["ylim"])   # LINEAR, tracks τ
        ax_den.set_xticks([])
        ax_den.set_yticks([])
        for s in ("top", "right"):
            ax_den.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax_den.spines[s].set_color("#9aa3af")
        ax_bar.set_xlim(0, 1); ax_bar.set_ylim(0, 1); ax_bar.axis("off")
        ax_bar.plot([0.5, 0.5], [0.16, 0.84], color="#9aa3af", lw=3,
                    solid_capstyle="round")
        y_star = 0.16 + 0.68 * d["s_star"]                # τ* mark on the bar
        ax_bar.plot([0.34, 0.66], [y_star, y_star], color="#8b949e", lw=1.1,
                    alpha=0.75, zorder=2)
        marker, = ax_bar.plot([0.5], [0.16], "o", color="#e08d00", ms=10, zorder=3)
        t0 = time.time()

        def update(frame):
            f = min(frame, NFRAMES - 1)
            s = f / (NFRAMES - 1)
            tau = d["taus"][f]
            st = frame_state(d, tau)
            pc.set_edgecolors(edge_flow_colors(st["edge_val"], st["vmax"],
                                               d["ealpha"], ecmap))   # edge colour = ρ(τ) flow
            sc.set_facecolor(st["ncol"])                  # nodes recolour by community
            im.set_data(st["rho_ord"]); im.set_clim(0, st["vmax"])
            ov_im.set_data(st["overlay"])
            dlc.set_segments(st["segs"]); dlc.set_color(st["tcol"])
            ax_den.set_ylim(*st["ylim"])
            marker.set_data([0.5], [0.16 + 0.68 * s])
            if frame % 16 == 0 or frame == NFRAMES + HOLD - 1:
                print(f"  [render {frame + 1:3d}/{NFRAMES + HOLD}] "
                      f"s={s:4.2f} τ={tau:.4g} elapsed={time.time() - t0:5.1f}s",
                      flush=True)
            return [pc, sc, im, ov_im, dlc, marker]

        anim = FuncAnimation(fig, update, frames=NFRAMES + HOLD, blit=False)
        out_vid.parent.mkdir(parents=True, exist_ok=True)
        print(f"writing MP4 -> {out_vid}", flush=True)
        writer = FFMpegWriter(fps=FPS, codec="libx264",
                              extra_args=["-vf",
                                          "pad=ceil(iw/2)*2:ceil(ih/2)*2:0:0:"
                                          "color=0xf6f6f8",
                                          "-pix_fmt", "yuv420p", "-crf", "18",
                                          "-preset", "medium"])
        anim.save(out_vid, writer=writer, dpi=DPI)
        plt.close(fig)
    print(f"done: {out_vid.name}  ({out_vid.stat().st_size / 1e6:.1f} MB)", flush=True)


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
            axn.add_collection(PathCollection(d["paths"], facecolors="none",
                                              edgecolors=eflow, linewidths=d["elw"]))
            axn.scatter(d["P"][:, 0], d["P"][:, 1], s=58,
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
                sp.set_color("#9aa3af")
            axd = axes[2, ci]
            axd.add_collection(LineCollection(st["segs"], colors=st["tcol"], linewidths=1.0))
            axd.set_xlim(0, 10 * d["n"]); axd.set_ylim(*st["ylim"])   # LINEAR, tracks τ
            axd.set_xticks([])
            axd.yaxis.set_major_locator(MaxNLocator(nbins=3))
            axd.tick_params(labelsize=7, colors="#1f2937")
            for sp in ("top", "right"):
                axd.spines[sp].set_visible(False)
            for sp in ("left", "bottom"):
                axd.spines[sp].set_color("#9aa3af")
        for r, lab in enumerate(("structure · ρ(τ) flow · communities",
                                 "ρ(τ) propagator · community blocks",
                                 "τ-morphing UPGMA tree")):
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
        make_animation(d, out_dir / "fig_fc_diffusion_tausweep.mp4")
    print(f"\nall assets in {out_dir}", flush=True)


if __name__ == "__main__":
    main()
