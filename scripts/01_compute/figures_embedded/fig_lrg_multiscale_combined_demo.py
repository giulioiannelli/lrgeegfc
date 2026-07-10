#!/usr/bin/env python
"""LRG multiscale combined demo — one latent hierarchy, two readouts.

Third figure for the 20-min talk. It fuses the two sibling demos into a single
side-by-side story about *what "multiscale" means* in each regime:

    Row 1 — SPARSE hierarchical modular network (structure is TOPOLOGICAL).
        Information flow as diffusion time τ slides: the heat-kernel propagator
        K(τ)=e^{-τ L̂} re-projected on the links, next to the K(τ) matrix.
        Nested communities resolve level-by-level as τ grows — the correct
        picture when the hierarchy lives in the *edge set*.

    Row 2 — REAL brain FC (imcoh_abs), a COMPLETE weighted connectome
        (structure lives in edge-weight HETEROGENEITY, not topology). The
        multiscale we actually USE: at the canonical τ=1/λ_max the propagator
        already encodes the whole hierarchy, so we do not sweep τ — we slide the
        dendrogram CUT height h. Three linked panels:
          · the network in an LRG communication-distance MDS layout (PCoA on
            1/K(τ_max)), curved edges, nodes coloured by community at h — a
            non-trivial single "bubble" whose lobes are the communication groups;
          · the connectivity matrix (leaf-ordered) with the community blocks
            projected onto it — coloured blocks appear / merge as the cut fills;
          · the cophenetic dendrogram with the growing cut.
        One shared palette: a community is the SAME colour in all three panels.

Both rows sweep fine → coarse together (small τ ↔ many communities; large τ ↔
one community), tied by a single vertical scale bar.

Reuses every primitive from the two siblings (build_hmn, laplacian_eig,
propagator, tau_bounds/at, cophenetic_linkage, dendrogram_segments, curved_paths,
cut_height, _labels_at, pal, _optimize_gif); only the assembly, the LRG
diffusion-map layout and the matrix-block overlay are new. Combinatorial
Laplacian only (project rule).

Outputs (data/outputs/figures/lrg_diffusion_zoom/):
    lrg_multiscale_combined.gif             animated combined sweep (talk asset)
    lrg_multiscale_combined_snapshots.pdf   static 3-stage contact sheet (vector)

Run (inside lapbrain):
    python scripts/01_compute/figures_embedded/fig_lrg_multiscale_combined_demo.py
    python .../fig_lrg_multiscale_combined_demo.py --check   # QA: 3 PNG frames
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import rc_context
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
from scipy.linalg import eigh
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
import imageio_ffmpeg
plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_lrg_diffusion_zoom_demo import (              # noqa: E402  (sibling reuse)
    build_hmn, build_fc_heterogeneous, laplacian_eig, propagator,
    tau_bounds, tau_at, REAL, SEED, NODE_CMAP, EDGE_CMAP, mid_id)
from fig_lrg_dendrogram_cut_demo import (              # noqa: E402  (sibling reuse)
    cophenetic_linkage, curved_paths, dendrogram_segments,
    cut_height, _labels_at, pal, GREY)
from lrg_eegfc.config.paths import FIGURES_ROOT        # noqa: E402
from lrg_eegfc.workflow.fc import load_fc_matrix        # noqa: E402
from lrg_eegfc.visuals.styles import use_lrg_style      # noqa: E402

NFRAMES = 96                             # sweep frames (slow, smooth)
FPS = 12
HOLD = 16                                # extra held frames on the coarse end
DPI = 150                                # true-colour MP4 → SHARP text at this res
MAT_CMAP = "magma"                       # matrix backgrounds (dark-friendly)


# ----------------------------------------------------------------------------
# NEW primitive: LRG diffusion-map layout — the non-trivial "one bubble"
# ----------------------------------------------------------------------------
def lrg_communication_layout(w, V, tau=None):
    """Classical MDS (PCoA) on the √-compressed LRG communication distance
    D_ij = 1/K_ij(τ) at the canonical τ = 1/λ_max — an LRG-native 2-D embedding.

    The propagator communication distance is the SAME quantity whose UPGMA tree
    we cut, so the layout and the dendrogram are consistent: nodes close in the
    plane merge early. MDS lays the complete weighted connectome out as a single
    organic "bubble" whose lobes are the communication-distance communities. A
    force layout would collapse the complete graph to a hairball; the raw
    heat-kernel eigenmode (Fiedler) embedding localizes on a few low-strength
    nodes and degenerates to a 1-D spike — the √ tames the 1/K dynamic range so
    the bulk spreads into a genuine 2-D cloud."""
    if tau is None:
        tau = 1.0 / w.max()
    K = propagator(w, V, tau)
    with np.errstate(divide="ignore"):
        D = 1.0 / K
    D[~np.isfinite(D)] = np.nanmax(D[np.isfinite(D)])
    np.fill_diagonal(D, 0.0)
    delta = np.sqrt(0.5 * (D + D.T))                  # √-compressed, symmetric
    n = D.shape[0]
    J = np.eye(n) - 1.0 / n
    B = -0.5 * J @ (delta ** 2) @ J                   # double-centred Gram
    ev, evec = eigh(B)
    idx = np.argsort(ev)[::-1][:2]                    # top-2 principal coordinates
    P = evec[:, idx] * np.sqrt(np.clip(ev[idx], 0.0, None))
    P -= P.mean(0)
    P /= (np.abs(P).max() + 1e-12)
    return P


def _community_matrix_rgb(labels, leaves, Mbright):
    """Direct community colouring of the leaf-ordered connectivity matrix:
    within-community cells take the community hue (value ∝ |FC| via Mbright),
    between-community cells stay a faint grey. Returns an (n, n, 3) RGB image —
    as the cut merges communities the coloured blocks grow in the *pixels*
    (no overlaid squares). Same idea as the diffusion snapshots' direct heatmap,
    but hue-coded by the current partition instead of by K value."""
    lab = labels[leaves]                              # community per matrix slot
    same = lab[:, None] == lab[None, :]
    ci = pal(lab)[:, :3]                              # (n, 3) hue per row
    b = Mbright[..., None]
    within = ci[:, None, :] * b                       # row-community hue × strength
    between = b * np.array([0.30, 0.32, 0.38])        # faint grey structure
    return np.where(same[..., None], within, between)


def _dist_tau(w, V, tau):
    """Symmetric communication distance D_ij = 1/K_ij(τ), finite off-diagonal."""
    K = propagator(w, V, tau)
    with np.errstate(divide="ignore"):
        D = 1.0 / K
    D[~np.isfinite(D)] = np.nanmax(D[np.isfinite(D)])
    np.fill_diagonal(D, 0.0)
    return 0.5 * (D + D.T)


def _tree_geometry_fixed_order(Z, xpos, n, floor, leaf_group):
    """Dendrogram polylines drawn with a FIXED leaf x-order so the tree morphs
    smoothly as τ (hence Z) changes frame-to-frame. Returns (segments, per-merge
    group id): the id is leaf_group's value if every leaf under the merge shares
    it, else -1 (→ grey). Lets row 1 show 'a different tree for every τ'."""
    cx = list(xpos) + [0.0] * (n - 1)
    ch = [floor] * n + [0.0] * (n - 1)
    grp = list(np.asarray(leaf_group, int)) + [-2] * (n - 1)
    segs, gids = [], []
    for m in range(Z.shape[0]):
        a, b = int(Z[m, 0]), int(Z[m, 1])
        H = max(float(Z[m, 2]), floor)
        xa, xb = cx[a], cx[b]
        ha, hb = max(ch[a], floor), max(ch[b], floor)
        segs.append(np.array([[xa, ha], [xa, H], [xb, H], [xb, hb]]))
        g = grp[a] if (grp[a] == grp[b] and grp[a] >= 0) else -1
        gids.append(g)
        cx[n + m] = 0.5 * (xa + xb)
        ch[n + m] = H
        grp[n + m] = g
    return segs, np.array(gids)


# ----------------------------------------------------------------------------
# per-graph preparation
# ----------------------------------------------------------------------------
def _process_hmn(A):
    """Row 1: sparse HMN, diffusion-τ sweep. Network (own links) + K(τ) matrix +
    a dendrogram RECOMPUTED at every τ (one tree per diffusion time)."""
    n = A.shape[0]
    w, V = laplacian_eig(A)
    t_min, t_max = tau_bounds(w)
    taus = tau_at(t_min, t_max, np.linspace(0.0, 1.0, NFRAMES))
    pos = nx.spring_layout(nx.from_numpy_array(A), weight="weight",
                           seed=SEED, k=1.6 / np.sqrt(n), iterations=300)
    P = np.array([pos[i] for i in range(n)])
    ei, ej = np.where(np.triu(A, 1) > 0)              # sparse own links
    iu, ju = np.triu_indices(n, 1)
    # τ-dependent dendrogram: FIXED leaf x-order (reference τ) + global height
    # range across the whole sweep, so the morphing tree stays framed & legible.
    ref = float(np.sqrt(t_min * t_max))
    order = leaves_list(linkage(squareform(_dist_tau(w, V, ref), checks=False),
                                "average"))
    xpos = np.zeros(n); xpos[order] = np.arange(n) * 10.0 + 5.0
    allh = np.concatenate([
        linkage(squareform(_dist_tau(w, V, t), checks=False), "average")[:, 2]
        for t in tau_at(t_min, t_max, np.linspace(0, 1, 14))])
    dfloor = float(allh[allh > 0].min()) * 0.5
    dtop = float(allh.max())
    return dict(n=n, w=w, V=V, P=P, ei=ei, ej=ej, iu=iu, ju=ju, order=order,
                segs=np.stack([P[ei], P[ej]], axis=1),
                node_rgba=NODE_CMAP(mid_id), node_size=min(34.0, 2400.0 / n),
                taus=taus, xpos=xpos, dfloor=dfloor, dtop=dtop,
                taus_static=tau_at(t_min, t_max, np.array([0.14, 0.5, 0.9])))


def _process_real(A, k_ref):
    """Row 2: dense brain FC, dendrogram-cut sweep. LRG diffusion-map bubble
    layout, curved strong edges, leaf-ordered matrix for the block overlay."""
    n = A.shape[0]
    w, V = laplacian_eig(A)
    Z = cophenetic_linkage(w, V)
    leaves, dseg, dh, drep, _ = dendrogram_segments(Z)
    floor = dh.min() * 0.5
    for s in dseg:                                    # clamp leaf feet for log-y
        s[:, 1] = np.maximum(s[:, 1], floor)
    P = lrg_communication_layout(w, V)                # NON-TRIVIAL BUBBLE
    ai, aj = np.where(np.triu(A, 1) > 0)              # ALL links (complete graph)
    ew = A[ai, aj]
    ealpha = (ew / ew.max()) ** 1.6                   # per-edge base opacity ∝ weight
    M = A[np.ix_(leaves, leaves)]                     # leaf-ordered connectivity
    mvmax = float(np.percentile(M, 97.0))
    Mbright = np.clip(M / max(mvmax, 1e-12), 0, 1) ** 0.6   # cell brightness ∝ |FC|
    k_start = n                                       # START FROM THE SMALLEST MERGE
    ks = np.round(np.linspace(k_start, 1, NFRAMES)).astype(int)
    ks_static = np.array([max(2, round(n / 3)), max(2, round(n / 14)), 1])
    return dict(n=n, Z=Z, Zh=Z[:, 2], leaves=leaves, dseg=dseg, dh=dh,
                drep=drep, floor=floor, P=P, ei=ai, ej=aj, ealpha=ealpha,
                paths=curved_paths(P, ai, aj, rad=0.08), M=M, mvmax=mvmax,
                Mbright=Mbright, node_size=min(40.0, 2800.0 / n),
                ks=ks, ks_static=ks_static, k_start=int(k_start))


def prepare():
    A_hmn = build_hmn(np.random.default_rng(SEED))
    try:
        A_real = np.asarray(load_fc_matrix(
            REAL["patient"], REAL["phase"], REAL["band"], "imcoh_abs"), float).copy()
        np.fill_diagonal(A_real, 0.0)
        A_real = np.clip(A_real, 0.0, None)
        print(f"  loaded real FC {REAL['patient']} {REAL['band']} "
              f"{REAL['phase']}: {A_real.shape[0]} nodes", flush=True)
    except Exception as exc:                                    # pragma: no cover
        print(f"  (real FC load failed: {exc}; synthetic fallback)", flush=True)
        A_real = build_fc_heterogeneous(np.random.default_rng(SEED))
    return {"hmn": _process_hmn(A_hmn),
            "real": _process_real(A_real, REAL["k"])}


# ----------------------------------------------------------------------------
# figure assembly (shared by the animation and the QA frame dump)
# ----------------------------------------------------------------------------
DARK = {"figure.facecolor": "#0b0e17", "axes.facecolor": "#0b0e17",
        "savefig.facecolor": "#0b0e17", "text.color": "#e6edf3",
        "axes.edgecolor": "#30363d", "axes.labelcolor": "#8b949e"}
LIGHT = {"figure.facecolor": "none", "axes.facecolor": "none",
         "savefig.facecolor": "none", "text.color": "#111827",
         "axes.edgecolor": "#6b7280", "axes.labelcolor": "#374151"}
CUT_C = "#e5006e"                        # dendrogram cut / rank bar (pink)
TAU_C = "#e08d00"                        # diffusion-time bar (amber)
BAR_Y0, BAR_Y1 = 0.16, 0.84             # vertical progress-bar track

# theme colours consumed by build_figure / helpers (flipped by _set_theme)
THEME = dict(fg="#e6edf3", muted="#8b949e", label="#c9d1d9", spine="#30363d",
             node_edge="#0b0e17", edge_cmap=EDGE_CMAP, transparent=False)


def _set_theme(light):
    """Dark (on-screen / dark slide) vs light-transparent (drop on a light
    slide: transparent bg, black text, dark-on-light network edges)."""
    if light:
        THEME.update(fg="#111827", muted="#374151", label="#1f2937",
                     spine="#6b7280", node_edge="#111827",
                     edge_cmap=plt.get_cmap("cividis_r"), transparent=True)
    else:
        THEME.update(fg="#e6edf3", muted="#8b949e", label="#c9d1d9",
                     spine="#30363d", node_edge="#0b0e17",
                     edge_cmap=EDGE_CMAP, transparent=False)


def build_figure(data):
    """Create the 2×4 figure and return (fig, update). Columns: network │
    connectivity matrix │ dendrogram │ progress bar. The rows differ in what the
    dendrogram *is*: row 1 (sparse HMN) sweeps τ and RECOMPUTES the tree at every
    diffusion time — a family of trees, one per τ; row 2 (dense FC) shows the
    SINGLE cophenetic tree at τ=1/λ_max and sweeps the cut. That contrast is the
    corrected point: τ-family of trees (topological) vs one hierarchy at τ_max."""
    hmn, real = data["hmn"], data["real"]
    fig = plt.figure(figsize=(13.6, 7.6), dpi=90)
    gs = fig.add_gridspec(2, 4, width_ratios=[1.10, 0.92, 1.0, 0.26],
                          left=0.05, right=0.985, top=0.845, bottom=0.055,
                          wspace=0.15, hspace=0.20)
    ax_hn, ax_hm, ax_hd, ax_tb = (fig.add_subplot(gs[0, c]) for c in range(4))
    ax_rn, ax_rm, ax_rd, ax_hb = (fig.add_subplot(gs[1, c]) for c in range(4))

    # ---- Row 1: HMN network (propagator on links) + K(τ) matrix + τ-tree ----
    hn_lc = LineCollection(hmn["segs"], linewidths=0.0, capstyle="round")
    ax_hn.add_collection(hn_lc)
    ax_hn.scatter(hmn["P"][:, 0], hmn["P"][:, 1], s=hmn["node_size"],
                  c=hmn["node_rgba"], edgecolors=THEME["node_edge"],
                  linewidths=0.4, zorder=3)
    _frame_off(ax_hn, hmn["P"])
    hm_im = ax_hm.imshow(np.zeros((hmn["n"], hmn["n"])), cmap=MAT_CMAP,
                         vmin=0, vmax=1, interpolation="nearest", animated=True)
    _matrix_axes(ax_hm)
    hd_lc = LineCollection([np.array([[0.0, hmn["dfloor"]], [0.0, hmn["dfloor"]]])],
                           linewidths=1.0)
    ax_hd.add_collection(hd_lc)
    ax_hd.set_yscale("log"); ax_hd.set_xlim(0, 10 * hmn["n"])
    ax_hd.set_ylim(hmn["dfloor"] * 0.9, hmn["dtop"] * 1.3)
    _dendro_axes(ax_hd)
    ax_hd.text(0.035, 0.955, "a new tree for every τ", transform=ax_hd.transAxes,
               ha="left", va="top", color=TAU_C, fontsize=8.5)
    tau_marker = _vbar(ax_tb, TAU_C, "τ", fine_coarse=True)

    # ---- Row 2: FC network (LRG bubble) + block matrix + fixed tree + cut ----
    rn_pc = PathCollection(real["paths"], facecolors="none",
                           edgecolors=GREY * [1, 1, 1, 0.14],
                           linewidths=0.2 + 0.9 * real["ealpha"])
    ax_rn.add_collection(rn_pc)
    rn_sc = ax_rn.scatter(real["P"][:, 0], real["P"][:, 1], s=real["node_size"],
                          c=np.tile(GREY, (real["n"], 1)), edgecolors=THEME["node_edge"],
                          linewidths=0.4, zorder=3)
    _frame_off(ax_rn, real["P"])
    rm_im = ax_rm.imshow(np.zeros((real["n"], real["n"], 3)),
                         interpolation="nearest", animated=True)
    _matrix_axes(ax_rm, n=real["n"])
    rd_lc = LineCollection(real["dseg"], linewidths=1.1)
    ax_rd.add_collection(rd_lc)
    rd_cut = ax_rd.axhline(real["floor"], color=CUT_C, lw=1.6, ls="--")
    ax_rd.set_yscale("log"); ax_rd.set_xlim(0, 10 * real["n"])
    ax_rd.set_ylim(real["floor"] * 0.9, real["dh"].max() * 1.3)
    _dendro_axes(ax_rd)
    ax_rd.text(0.035, 0.955, "one tree at τ_max · cut ↑", transform=ax_rd.transAxes,
               ha="left", va="top", color=CUT_C, fontsize=8.5)
    rd_kcount = ax_rd.text(0.96, 0.955, "", transform=ax_rd.transAxes, ha="right",
                           va="top", color=THEME["fg"], fontsize=9.5)
    rank_marker = _vbar(ax_hb, CUT_C, "h", fine_coarse=False)

    # ---- figure-level headers / row labels / thesis line ----
    fig.text(0.5, 0.965, "one latent hierarchy — two readouts", ha="center",
             color=THEME["fg"], fontsize=13, fontweight="bold")
    fig.text(0.5, 0.928,
             "sparse · topological → sweep diffusion time τ        vs.        "
             "dense · weighted → one tree at τ_max, sweep the cut h",
             ha="center", color=THEME["muted"], fontsize=9.5)
    fig.text(0.205, 0.882, "network", ha="center", color=THEME["muted"], fontsize=9)
    fig.text(0.485, 0.882, "connectivity matrix", ha="center",
             color=THEME["muted"], fontsize=9)
    fig.text(0.760, 0.882, "dendrogram", ha="center", color=THEME["muted"], fontsize=9)
    fig.text(0.018, 0.650, "sparse HMN  ·  diffusion τ", rotation=90,
             va="center", ha="center", color=THEME["label"], fontsize=9.5)
    fig.text(0.018, 0.235, "brain FC  ·  dendrogram cut", rotation=90,
             va="center", ha="center", color=THEME["label"], fontsize=9.5)

    art = dict(hn_lc=hn_lc, hm_im=hm_im, hd_lc=hd_lc, rn_pc=rn_pc, rn_sc=rn_sc,
               rm_im=rm_im, rd_lc=rd_lc, rd_cut=rd_cut, rd_kcount=rd_kcount,
               tau_marker=tau_marker, rank_marker=rank_marker)

    def update(frame):
        f = min(frame, NFRAMES - 1)                   # clamp → hold on coarse end
        s = f / (NFRAMES - 1)
        changed = []
        # Row 1 — HMN diffusion: network + K(τ) matrix + τ-dependent tree
        tau = hmn["taus"][f]
        K = propagator(hmn["w"], hmn["V"], tau)
        vmax = max(np.percentile(K[hmn["iu"], hmn["ju"]], 99.0), 1e-9)
        val = np.clip(K[hmn["ei"], hmn["ej"]] / vmax, 0, 1)
        rgba = THEME["edge_cmap"](val); rgba[:, 3] = val ** 0.75
        art["hn_lc"].set_color(rgba)
        art["hn_lc"].set_linewidth(0.15 + 2.6 * val)
        art["hm_im"].set_data(K[np.ix_(hmn["order"], hmn["order"])])
        art["hm_im"].set_clim(0, vmax)
        D = _dist_tau(hmn["w"], hmn["V"], tau)         # recompute the tree at THIS τ
        Ztau = linkage(squareform(D, checks=False), "average")
        hts = Ztau[:, 2]
        top, floor_f = float(hts.max()), float(hts.min()) * 0.4
        dsegs, gids = _tree_geometry_fixed_order(Ztau, hmn["xpos"], hmn["n"],
                                                 floor_f, mid_id)
        dcol = np.tile(GREY, (len(dsegs), 1))
        mono = gids >= 0
        if mono.any():
            dcol[mono] = NODE_CMAP(gids[mono])
        art["hd_lc"].set_segments(dsegs)
        art["hd_lc"].set_color(dcol)
        ax_hd.set_ylim(floor_f * 0.85, top * 1.15)     # rescale so the tree always
        changed += [art["hn_lc"], art["hm_im"], art["hd_lc"]]  #  fills the y-axis
        # Row 2 — FC dendrogram cut: fixed τ_max tree, moving cut
        k = int(real["ks"][f])
        labels = _labels_at(real, k)
        h = cut_height(real, k)
        art["rn_sc"].set_facecolor(pal(labels))
        # ALL links drawn: faint weighted grey wash + within-community glow
        same = labels[real["ei"]] == labels[real["ej"]]
        ecol = np.empty((real["ei"].size, 4))
        ecol[:] = GREY
        ecol[:, 3] = real["ealpha"] * 0.16
        ecol[same] = pal(labels[real["ei"][same]])
        ecol[same, 3] = np.clip(real["ealpha"][same] * 0.85, 0.04, 0.9)
        art["rn_pc"].set_edgecolor(ecol)
        # matrix: DIRECT community colouring (blocks grow in the pixels)
        art["rm_im"].set_data(_community_matrix_rgb(labels, real["leaves"],
                                                    real["Mbright"]))
        below = real["dh"] < h
        dcol2 = np.tile(GREY, (len(real["dseg"]), 1))
        dcol2[below] = pal(labels[real["leaves"][real["drep"][below]]])
        art["rd_lc"].set_color(dcol2)
        art["rd_cut"].set_ydata([h, h])
        art["rd_kcount"].set_text(f"{labels.max() + 1} communities")
        changed += [art["rn_sc"], art["rn_pc"], art["rm_im"],
                    art["rd_lc"], art["rd_cut"], art["rd_kcount"]]
        # two independent bars: τ (diffusion, row 1) and h (cut by rank, row 2)
        art["tau_marker"].set_data([0.5], [BAR_Y0 + (BAR_Y1 - BAR_Y0) * s])
        s_rank = (real["n"] - k) / (real["n"] - 1)     # merge rank fraction
        art["rank_marker"].set_data([0.5], [BAR_Y0 + (BAR_Y1 - BAR_Y0) * s_rank])
        changed += [art["tau_marker"], art["rank_marker"]]
        if frame % 16 == 0 or frame == NFRAMES + HOLD - 1:
            print(f"  [render {frame + 1:3d}/{NFRAMES + HOLD}] s={s:4.2f} "
                  f"k={k}", flush=True)
        return changed

    return fig, update


def _frame_off(ax, P, pad=0.12):
    ax.set_xlim(P[:, 0].min() - pad, P[:, 0].max() + pad)
    ax.set_ylim(P[:, 1].min() - pad, P[:, 1].max() + pad)
    ax.set_aspect("equal"); ax.axis("off")


def _matrix_axes(ax, n=None):
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color(THEME["spine"])


def _dendro_axes(ax):
    ax.set_xticks([]); ax.tick_params(axis="y", labelsize=7, colors=THEME["muted"])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(THEME["spine"])


def _vbar(ax, color, label, fine_coarse=False):
    """Thin vertical progress bar in its own narrow column; returns the marker."""
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.plot([0.5, 0.5], [BAR_Y0, BAR_Y1], color=THEME["spine"], lw=3,
            solid_capstyle="round", zorder=1)
    ax.text(0.5, 0.955, label, ha="center", va="bottom", color=color,
            fontsize=13, fontweight="bold")
    if fine_coarse:
        ax.text(0.5, BAR_Y1 + 0.01, "coarse", ha="center", va="bottom",
                color=THEME["muted"], fontsize=6.5, rotation=90)
        ax.text(0.5, BAR_Y0 - 0.01, "fine", ha="center", va="top",
                color=THEME["muted"], fontsize=6.5, rotation=90)
    m, = ax.plot([0.5], [BAR_Y0], "o", color=color, ms=11, zorder=3)
    return m


# ----------------------------------------------------------------------------
# animation + QA frame dump
# ----------------------------------------------------------------------------
def make_animation(data, out_mp4: Path, out_gif: Path | None = None):
    """Primary deliverable: a true-colour MP4 (SHARP text — no palette
    quantization). Optional GIF fallback (per-frame adaptive 256-colour)."""
    total = NFRAMES + HOLD
    out_mp4.parent.mkdir(parents=True, exist_ok=True)
    with rc_context(DARK):
        fig, update = build_figure(data)
        t0 = time.time()
        anim = FuncAnimation(fig, update, frames=total, blit=False)
        print(f"writing MP4 -> {out_mp4}", flush=True)
        writer = FFMpegWriter(fps=FPS, codec="libx264",
                              extra_args=["-pix_fmt", "yuv420p", "-crf", "18",
                                          "-preset", "medium"])
        anim.save(out_mp4, writer=writer, dpi=DPI)
        print(f"  MP4 render wall-clock {time.time() - t0:.1f}s "
              f"({out_mp4.stat().st_size / 1e6:.1f} MB)", flush=True)
        if out_gif is not None:
            print(f"writing GIF -> {out_gif}", flush=True)
            anim.save(out_gif, writer=PillowWriter(fps=FPS), dpi=110)
        plt.close(fig)
    if out_gif is not None:
        _sharpen_gif(out_gif)
        print(f"done: {out_gif.name}  ({out_gif.stat().st_size / 1e6:.1f} MB)",
              flush=True)
    print(f"done: {out_mp4.name}  ({total} frames @ {FPS} fps, dpi {DPI})",
          flush=True)


def _sharpen_gif(path: Path):
    """Re-encode with a per-frame adaptive 256-colour palette — far sharper text
    than a shared low-colour palette. Temp-file safe."""
    try:
        from PIL import Image
        im = Image.open(path)
        rgb = []
        try:
            while True:
                rgb.append(im.copy().convert("RGB"))
                im.seek(im.tell() + 1)
        except EOFError:
            pass
        pframes = []
        for fr in rgb:
            pf = fr.quantize(colors=256, method=Image.MEDIANCUT, dither=Image.NONE)
            pf.info.pop("transparency", None)
            pframes.append(pf)
        tmp = path.with_suffix(".sharp.gif")
        pframes[0].save(tmp, save_all=True, append_images=pframes[1:],
                        duration=int(1000 / FPS), loop=0, optimize=True)
        tmp.replace(path)
    except Exception as exc:                         # pragma: no cover
        print(f"  (gif sharpen skipped: {exc})", flush=True)


def _transparent_gif(frames, out_path: Path, fps):
    """Assemble RGBA PIL frames into a GIF with a TRANSPARENT background
    (1-bit alpha: a reserved palette index marks fully-see-through pixels)."""
    from PIL import Image
    pframes = []
    for im in frames:
        alpha = im.getchannel("A")
        p = im.convert("RGB").quantize(colors=255, method=Image.MEDIANCUT,
                                       dither=Image.NONE)          # leave idx 255
        mask = alpha.point(lambda a: 255 if a < 128 else 0)        # low α → clear
        p.paste(255, mask)
        p.info["transparency"] = 255
        pframes.append(p)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pframes[0].save(out_path, save_all=True, append_images=pframes[1:],
                    duration=int(1000 / fps), loop=0, disposal=2,
                    transparency=255, optimize=False)


def render_transparent_gif(data, out_path: Path, dpi=110):
    """Light-theme, TRANSPARENT-background animated GIF for a light slide
    (MP4 cannot carry transparency, so the animated transparent asset is a GIF).
    Black text, dark-on-light network edges."""
    import io
    from PIL import Image
    _set_theme(light=True)
    total, t0 = NFRAMES + HOLD, time.time()
    with rc_context(LIGHT):
        fig, update = build_figure(data)
        frames = []
        for f in range(total):
            update(f)
            buf = io.BytesIO()
            fig.savefig(buf, format="png", transparent=True, dpi=dpi)
            buf.seek(0)
            frames.append(Image.open(buf).convert("RGBA"))
            if f % 16 == 0 or f == total - 1:
                print(f"  [transp {f + 1:3d}/{total}] {time.time() - t0:5.1f}s",
                      flush=True)
        plt.close(fig)
    _transparent_gif(frames, out_path, FPS)
    print(f"done: {out_path.name}  ({out_path.stat().st_size / 1e6:.1f} MB, "
          f"transparent bg, black text)", flush=True)


def save_frames(data, frames, out_dir: Path, light=False):
    """QA: render a few frames to PNG at the MP4 dpi (scratch — not a
    deliverable). light=True previews the transparent theme on a white matte."""
    _set_theme(light)
    out_dir.mkdir(parents=True, exist_ok=True)
    with rc_context(LIGHT if light else DARK):
        fig, update = build_figure(data)
        for fr in frames:
            update(fr)
            t0 = time.time()
            p = out_dir / f"combined_frame_{fr:03d}.png"
            fig.savefig(p, dpi=DPI, facecolor="white" if light else "#0b0e17")
            print(f"  wrote {p.name}  ({time.time() - t0:.2f}s/frame)", flush=True)
        plt.close(fig)


# ----------------------------------------------------------------------------
# static contact sheet (vector PDF, light theme) — 6 rows × 3 scale stages
# ----------------------------------------------------------------------------
def make_static(data, out_path: Path, transparent=False):
    hmn, real = data["hmn"], data["real"]
    ecmap = plt.get_cmap("turbo")
    stages = ("fine", "intermediate", "coarse")
    rc = ({"figure.facecolor": "none", "savefig.facecolor": "none"} if transparent
          else {"figure.facecolor": "white", "savefig.facecolor": "white"})
    with rc_context(rc):
        fig, axes = plt.subplots(6, 3, figsize=(9.2, 17.2))
        fig.subplots_adjust(left=0.10, right=0.985, top=0.97, bottom=0.012,
                            wspace=0.10, hspace=0.17)
        # rows 0-2: HMN network + K(τ) matrix + τ-dependent dendrogram over 3 τ
        for ci, tau in enumerate(hmn["taus_static"]):
            K = propagator(hmn["w"], hmn["V"], tau)
            vmax = max(np.percentile(K[hmn["iu"], hmn["ju"]], 99.0), 1e-9)
            val = np.clip(K[hmn["ei"], hmn["ej"]] / vmax, 0, 1)
            axn = axes[0, ci]
            rgba = ecmap(val); rgba[:, 3] = val ** 0.8
            axn.add_collection(LineCollection(hmn["segs"], colors=rgba,
                                              linewidths=0.1 + 2.2 * val,
                                              capstyle="round"))
            axn.scatter(hmn["P"][:, 0], hmn["P"][:, 1], s=0.5 * hmn["node_size"],
                        c=hmn["node_rgba"], edgecolors="white", linewidths=0.3, zorder=3)
            _frame_off(axn, hmn["P"])
            axes[1, ci].imshow(K[np.ix_(hmn["order"], hmn["order"])], cmap=ecmap,
                               vmin=0, vmax=vmax, interpolation="nearest")
            axes[1, ci].set_xticks([]); axes[1, ci].set_yticks([])
            # τ-dependent dendrogram (recomputed at this τ, y-axis rescaled to fill)
            D = _dist_tau(hmn["w"], hmn["V"], tau)
            Ztau = linkage(squareform(D, checks=False), "average")
            hts = Ztau[:, 2]; floor_f = float(hts.min()) * 0.4
            dsegs, gids = _tree_geometry_fixed_order(Ztau, hmn["xpos"], hmn["n"],
                                                     floor_f, mid_id)
            dcol = np.tile([0.5, 0.5, 0.55, 1.0], (len(dsegs), 1))
            mono = gids >= 0
            if mono.any():
                dcol[mono] = NODE_CMAP(gids[mono])
            axd1 = axes[2, ci]
            axd1.add_collection(LineCollection(dsegs, colors=dcol, linewidths=0.9))
            axd1.set_yscale("log"); axd1.set_xlim(0, 10 * hmn["n"])
            axd1.set_ylim(floor_f * 0.85, float(hts.max()) * 1.15)
            axd1.set_xticks([]); axd1.tick_params(labelsize=7)
            for sp in ("top", "right"):
                axd1.spines[sp].set_visible(False)
            axes[0, ci].set_title(f"{stages[ci]} (τ)", fontsize=10)
        _blank_ylabel(axes[0, 0], "HMN · network", on=True)
        _blank_ylabel(axes[1, 0], "HMN · K(τ) matrix", on=True)
        _blank_ylabel(axes[2, 0], "HMN · dendrogram at τ (morphs)", on=True)
        # rows 3-5: FC network + block matrix + fixed dendrogram over 3 cuts
        for ci, k in enumerate(real["ks_static"]):
            labels = _labels_at(real, k)
            same = labels[real["ei"]] == labels[real["ej"]]
            ecol = np.tile([0.7, 0.7, 0.72, 0.10], (real["ei"].size, 1))
            ecol[same] = pal(labels[real["ei"][same]]) * [1, 1, 1, 0.6]
            axn = axes[3, ci]
            axn.add_collection(PathCollection(real["paths"], facecolors="none",
                                              edgecolors=ecol, linewidths=0.5))
            axn.scatter(real["P"][:, 0], real["P"][:, 1], s=0.5 * real["node_size"],
                        c=pal(labels), edgecolors="white", linewidths=0.3, zorder=3)
            _frame_off(axn, real["P"])
            axn.set_title(f"{k} comm.", fontsize=10)
            axm = axes[4, ci]
            axm.imshow(_community_matrix_rgb(labels, real["leaves"], real["Mbright"]),
                       interpolation="nearest")
            axm.set_xticks([]); axm.set_yticks([])
            axd = axes[5, ci]
            below = real["dh"] < cut_height(real, k)
            dcol = np.tile([0.4, 0.44, 0.5, 1.0], (len(real["dseg"]), 1))
            dcol[below] = pal(labels[real["leaves"][real["drep"][below]]])
            axd.add_collection(LineCollection(real["dseg"], colors=dcol, linewidths=0.9))
            axd.axhline(cut_height(real, k), color="#c1121f", lw=1.0, ls="--")
            axd.set_yscale("log"); axd.set_xlim(0, 10 * real["n"])
            axd.set_ylim(real["floor"] * 0.9, real["dh"].max() * 1.3)
            axd.set_xticks([]); axd.tick_params(labelsize=7)
            for sp in ("top", "right"):
                axd.spines[sp].set_visible(False)
        _blank_ylabel(axes[3, 0], "brain FC · LRG-layout network", on=True)
        _blank_ylabel(axes[4, 0], "brain FC · block matrix", on=True)
        _blank_ylabel(axes[5, 0], "brain FC · one tree at τ_max · cut", on=True)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, transparent=transparent)
        plt.close(fig)
    print(f"done: {out_path.name}", flush=True)


def _blank_ylabel(ax, label, on=False):
    if on:
        ax.axis("on"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    ax.set_ylabel(label, fontsize=9.5, rotation=90, labelpad=10)


def main():
    use_lrg_style()
    out_dir = Path(FIGURES_ROOT) / "lrg_diffusion_zoom"
    print("preparing graphs + eigendecompositions ...", flush=True)
    data = prepare()
    for tag, d in (("HMN", data["hmn"]), ("FC ", data["real"])):
        print(f"  {tag}: N={d['n']}", flush=True)
    if "--check" in sys.argv:
        scratch = Path("/tmp/claude-1000/-home-giulio-Documents-research-neural-"
                       "networks-lrgeegfc/093cb33a-9501-421a-b1c9-5c82211b57c0/"
                       "scratchpad")
        light = "--light" in sys.argv or "--transparent" in sys.argv
        save_frames(data, [1, NFRAMES // 4, NFRAMES // 2,
                           3 * NFRAMES // 4, NFRAMES - 1], scratch, light=light)
        return
    # dark, sharp H.264 video — for on-screen / dark slides
    _set_theme(light=False)
    make_animation(data, out_dir / "lrg_multiscale_combined.mp4")
    make_static(data, out_dir / "lrg_multiscale_combined_snapshots.pdf")
    # light, TRANSPARENT background + black text — to drop on a light slide
    make_static(data, out_dir / "lrg_multiscale_combined_transparent.pdf",
                transparent=True)
    render_transparent_gif(data, out_dir / "lrg_multiscale_combined_transparent.gif")
    print(f"\nall assets in {out_dir}", flush=True)


if __name__ == "__main__":
    main()
