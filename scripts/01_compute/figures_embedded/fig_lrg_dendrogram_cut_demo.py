#!/usr/bin/env python
"""LRG dendrogram-cut demo — the multiscale we actually USE.

Companion to fig_lrg_diffusion_zoom_demo.py. That figure sweeps diffusion time
τ; THIS one makes the point that we do NOT need a τ-sweep. At the canonical
τ = 1/λ_max the propagator K(τ) already encodes the whole hierarchy: the
communication distance D_ij = 1/K_ij(τ) → UPGMA dendrogram, and the *multiscale*
structure is read by sliding the dendrogram CUT HEIGHT h from all-leaves-apart
(each node its own community) up to all-merged (one community).

For each graph — a synthetic HMN (known hierarchy) and a REAL brain FC matrix
(imcoh_abs, a complete weighted connectome) — the animation shows, side by side:
    left  : the network (graph-tool sfdp layout, curved edges), nodes coloured
            by the community they belong to at the current cut height h;
    right : the cophenetic dendrogram with a horizontal cut that grows, its
            branches coloured to match the network communities.
As h grows the colours coalesce until every node is one colour — the multiscale
nature of a graph that "lives at one scale", encoded in edge-weight heterogeneity.

Outputs (data/outputs/figures/lrg_diffusion_zoom/):
    lrg_dendrogram_cut.gif            animated cut-height sweep (talk asset)
    lrg_dendrogram_cut_snapshots.pdf  static 3-cut contact sheet (vector)

Run (inside lapbrain):
    python scripts/01_compute/figures_embedded/fig_lrg_dendrogram_cut_demo.py
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc_context
from matplotlib.path import Path as MPath
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_lrg_diffusion_zoom_demo import (              # noqa: E402  (sibling reuse)
    build_hmn, laplacian_eig, propagator, REAL, SEED, _optimize_gif)
from lrg_eegfc.config.paths import FIGURES_ROOT       # noqa: E402
from lrg_eegfc.workflow.fc import load_fc_matrix       # noqa: E402
from lrg_eegfc.visuals.styles import use_lrg_style     # noqa: E402

GRAPHS = ("Hierarchical modular network (synthetic)",
          f"Brain FC · {REAL['band']} · {REAL['patient']} {REAL['phase']}")
NFRAMES = 54
FPS = 16
TAB20 = np.array([plt.get_cmap("tab20")(i) for i in range(20)])
GREY = np.array([0.42, 0.46, 0.52, 1.0])


def pal(idx):
    """Categorical colour(s) from tab20 (cycled)."""
    return TAB20[np.asarray(idx) % 20]


# ----------------------------------------------------------------------------
# cophenetic dendrogram at the canonical τ = 1/λ_max
# ----------------------------------------------------------------------------
def cophenetic_linkage(w, V):
    """D_ij = 1/K_ij at τ = 1/λ_max (canonical LRG τ) → UPGMA linkage."""
    K = propagator(w, V, 1.0 / w.max())          # off-diagonal >= 0
    with np.errstate(divide="ignore"):
        D = 1.0 / K
    D[~np.isfinite(D)] = np.nanmax(D[np.isfinite(D)])
    np.fill_diagonal(D, 0.0)
    D = 0.5 * (D + D.T)
    return linkage(squareform(D, checks=False), method="average")


def relabel_by_leaforder(raw, leaves):
    """Relabel community ids by first appearance along the leaf order so the
    palette is assigned left→right — matches scipy's dendrogram colouring."""
    remap, out = {}, np.zeros_like(raw)
    for node in leaves:
        remap.setdefault(raw[node], len(remap))
        out[node] = remap[raw[node]]
    return out


# ----------------------------------------------------------------------------
# graph-tool layout + curved matplotlib edges
# ----------------------------------------------------------------------------
def sfdp_positions(A, groups):
    """graph-tool sfdp layout (weighted, community-grouped) → (n, 2) array."""
    import graph_tool.all as gt
    n = A.shape[0]
    g = gt.Graph(directed=False)
    g.add_vertex(n)
    ew = g.new_edge_property("double")
    ei, ej = np.where(np.triu(A, 1) > 0)
    for a, b in zip(ei.tolist(), ej.tolist()):
        e = g.add_edge(a, b)
        ew[e] = float(A[a, b])
    grp = g.new_vertex_property("int")
    grp.a = np.asarray(groups, int)
    pos = gt.sfdp_layout(g, eweight=ew, groups=grp, C=0.5)
    P = np.array([list(pos[v]) for v in g.vertices()], float)
    P -= P.mean(0)
    P /= np.abs(P).max() + 1e-9
    return P


def _knn_backbone(A, m=6):
    """Symmetric top-m-per-node backbone — a sparse skeleton of a dense graph
    so sfdp can separate communities instead of collapsing to a hairball."""
    n = A.shape[0]
    B = np.zeros_like(A)
    for i in range(n):
        idx = np.argsort(A[i])[::-1][:m]
        B[i, idx] = A[i, idx]
    return np.maximum(B, B.T)


def curved_paths(P, ei, ej, rad=0.14):
    """One quadratic-Bézier Path per edge (gentle curve for a 'nice' look)."""
    paths = []
    for a, b in zip(ei, ej):
        p0, p1 = P[a], P[b]
        d = p1 - p0
        ctrl = 0.5 * (p0 + p1) + rad * np.array([-d[1], d[0]])
        paths.append(MPath([p0, ctrl, p1],
                           [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3]))
    return paths


def dendrogram_segments(Z):
    """scipy dendrogram geometry → per-link polyline (4 pts), merge height,
    and a representative leaf position (leftmost leaf under the link)."""
    dd = dendrogram(Z, no_plot=True)
    leaves = np.array(dd["leaves"])              # node index at each x-position
    segs, heights, rep = [], [], []
    for xs, ys in zip(dd["icoord"], dd["dcoord"]):
        segs.append(np.column_stack([xs, ys]))   # (x0,y0)->(x0,yt)->(x1,yt)->(x1,y1)
        heights.append(ys[1])                    # merge height (top of the ⊓)
        rep.append(int((min(xs) - 5.0) / 10.0))  # leftmost leaf position
    return leaves, segs, np.array(heights), np.array(rep, int), dd


# ----------------------------------------------------------------------------
# assembly
# ----------------------------------------------------------------------------
def _process_graph(A, is_real, k_ref):
    n = A.shape[0]
    w, V = laplacian_eig(A)
    Z = cophenetic_linkage(w, V)
    leaves, dseg, dh, drep, _ = dendrogram_segments(Z)
    floor = dh.min() * 0.5
    for s in dseg:                               # clamp leaf feet for log y-axis
        s[:, 1] = np.maximum(s[:, 1], floor)
    # reference communities → spatial groups for the layout (fixed)
    groups = relabel_by_leaforder(fcluster(Z, t=k_ref, criterion="maxclust"), leaves)
    # layout graph: the sparse HMN as-is; a kNN backbone for the dense connectome
    # (sfdp on a complete graph collapses to a hairball)
    A_layout = A if not is_real else _knn_backbone(A, m=6)
    P = sfdp_positions(A_layout, groups)
    # drawn edges: all for the sparse HMN, strong subset for the dense connectome
    ai, aj = np.where(np.triu(A, 1) > 0)
    if is_real:
        keep = A[ai, aj] >= np.quantile(A[ai, aj], 0.85)   # top ~15% (clarity)
        ai, aj = ai[keep], aj[keep]
    # sweep by TARGET community count so the count decreases smoothly — heights
    # are top-heavy for a dense connectome, so a height-linear sweep stalls
    k_start = min(30, n - 1)
    ks = np.round(np.linspace(k_start, 1, NFRAMES)).astype(int)
    ks_static = np.array([k_start, max(2, round(k_start / 6)), 1])
    return dict(n=n, Z=Z, Zh=Z[:, 2], leaves=leaves, dseg=dseg, dh=dh, drep=drep,
                P=P, floor=floor, ei=ai, ej=aj, paths=curved_paths(P, ai, aj),
                node_size=min(34.0, 2400.0 / n), ks=ks, ks_static=ks_static)


def cut_height(d, k):
    """Cophenetic height of the cut that yields exactly k clusters (log-mid of
    the two bracketing merges) — where the horizontal cut line is drawn."""
    n, Zh = d["n"], d["Zh"]
    if k >= n:
        return d["floor"]
    if k <= 1:
        return d["dh"].max() * 1.1
    return float(np.sqrt(Zh[n - k - 1] * Zh[n - k]))


def prepare():
    A_hmn = build_hmn(np.random.default_rng(SEED))
    A_real = np.asarray(load_fc_matrix(
        REAL["patient"], REAL["phase"], REAL["band"], "imcoh_abs"), float).copy()
    np.fill_diagonal(A_real, 0.0)
    A_real = np.clip(A_real, 0.0, None)
    print(f"  real FC {REAL['patient']} {REAL['band']} {REAL['phase']}: "
          f"{A_real.shape[0]} nodes", flush=True)
    return {GRAPHS[0]: _process_graph(A_hmn, False, 4),
            GRAPHS[1]: _process_graph(A_real, True, REAL["k"])}


def _labels_at(d, k):
    """Relabelled community id per node at target community count k."""
    return relabel_by_leaforder(
        fcluster(d["Z"], t=k, criterion="maxclust"), d["leaves"])


# ----------------------------------------------------------------------------
# animation
# ----------------------------------------------------------------------------
def make_animation(data, out_path: Path):
    dark = {"figure.facecolor": "#0b0e17", "axes.facecolor": "#0b0e17",
            "savefig.facecolor": "#0b0e17", "text.color": "#e6edf3",
            "axes.edgecolor": "#30363d", "axes.labelcolor": "#8b949e"}
    with rc_context(dark):
        fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.6), dpi=90,
                                 gridspec_kw=dict(width_ratios=[1.35, 1.0]))
        fig.subplots_adjust(left=0.02, right=0.97, top=0.85, bottom=0.10,
                            wspace=0.10, hspace=0.22)
        art = {}
        for r, name in enumerate(GRAPHS):
            d = data[name]
            ax_net, ax_den = axes[r, 0], axes[r, 1]
            # network: static curved edges + recolourable nodes
            pc = PathCollection(d["paths"], facecolors="none",
                                edgecolors=GREY * [1, 1, 1, 0.16], linewidths=0.5)
            ax_net.add_collection(pc)
            sc = ax_net.scatter(d["P"][:, 0], d["P"][:, 1], s=d["node_size"],
                                c=np.tile(GREY, (d["n"], 1)),
                                edgecolors="#0b0e17", linewidths=0.4, zorder=3)
            pad = 0.12
            ax_net.set_xlim(d["P"][:, 0].min() - pad, d["P"][:, 0].max() + pad)
            ax_net.set_ylim(d["P"][:, 1].min() - pad, d["P"][:, 1].max() + pad)
            ax_net.set_aspect("equal"); ax_net.axis("off")
            ax_net.set_title(name, color="#e6edf3", fontsize=10.5, pad=4,
                             loc="left")
            # dendrogram: recolourable links + moving cut line
            dlc = LineCollection(d["dseg"], linewidths=1.1)
            ax_den.add_collection(dlc)
            cut = ax_den.axhline(d["floor"], color="#f72585", lw=1.6, ls="--")
            ax_den.set_yscale("log")
            ax_den.set_xlim(0, 10 * d["n"])
            ax_den.set_ylim(d["floor"] * 0.9, d["dh"].max() * 1.3)
            ax_den.set_xticks([])
            ax_den.tick_params(axis="y", labelsize=7, colors="#8b949e")
            for s in ax_den.spines.values():
                s.set_color("#30363d")
            ax_den.spines["top"].set_visible(False)
            ax_den.spines["right"].set_visible(False)
            kcount = ax_den.text(0.97, 0.95, "", transform=ax_den.transAxes,
                                 ha="right", va="top", color="#e6edf3", fontsize=10)
            art[name] = dict(pc=pc, sc=sc, dlc=dlc, cut=cut, kcount=kcount, ecol=None)
        fig.text(0.30, 0.925, "network — coloured by community at cut h",
                 ha="center", color="#8b949e", fontsize=9)
        fig.text(0.78, 0.925, "cophenetic dendrogram  ·  cut h ↑",
                 ha="center", color="#8b949e", fontsize=9)
        htxt = fig.text(0.5, 0.955, "", ha="center", color="#e6edf3",
                        fontsize=12, fontweight="bold")
        t0 = time.time()

        def update(frame):
            changed = [htxt]
            for name in GRAPHS:
                d, a = data[name], art[name]
                labels = _labels_at(d, d["ks"][frame])
                h = cut_height(d, d["ks"][frame])
                a["sc"].set_facecolor(pal(labels))
                # edges: within-community glow in community colour, else faint grey
                same = labels[d["ei"]] == labels[d["ej"]]
                ecol = np.tile(GREY * [1, 1, 1, 0.10], (d["ei"].size, 1))
                ecol[same] = pal(labels[d["ei"][same]]) * [1, 1, 1, 0.55]
                a["pc"].set_edgecolor(ecol)
                # dendrogram links: below cut → community colour, above → grey
                below = d["dh"] < h
                dcol = np.tile(GREY, (len(d["dseg"]), 1))
                reps = d["leaves"][d["drep"][below]]
                dcol[below] = pal(labels[reps])
                a["dlc"].set_color(dcol)
                a["cut"].set_ydata([h, h])
                a["kcount"].set_text(f"{labels.max() + 1} communities")
                changed += [a["sc"], a["pc"], a["dlc"], a["cut"], a["kcount"]]
            htxt.set_text(f"dendrogram cut  h  —  scale {frame / (NFRAMES - 1):4.0%}")
            if frame % 12 == 0 or frame == NFRAMES - 1:
                print(f"  [render {frame + 1:2d}/{NFRAMES}] "
                      f"elapsed={time.time() - t0:5.1f}s", flush=True)
            return changed

        anim = FuncAnimation(fig, update, frames=NFRAMES, blit=False)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"writing GIF -> {out_path}", flush=True)
        anim.save(out_path, writer=PillowWriter(fps=FPS))
        plt.close(fig)
    _optimize_gif(out_path)
    print(f"done: {out_path.name}  ({out_path.stat().st_size / 1e6:.1f} MB)",
          flush=True)


# ----------------------------------------------------------------------------
# static contact sheet (vector PDF, light theme)
# ----------------------------------------------------------------------------
def make_static(data, out_path: Path):
    with rc_context({"figure.facecolor": "white", "savefig.facecolor": "white"}):
        fig, axes = plt.subplots(2, 4, figsize=(11.5, 6.2),
                                 gridspec_kw=dict(width_ratios=[1, 1, 1, 1.1]))
        fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.03,
                            wspace=0.12, hspace=0.18)
        titles = ("fine cut", "intermediate", "coarse cut")
        for gi, name in enumerate(GRAPHS):
            d = data[name]
            for ci, k in enumerate(d["ks_static"]):
                labels = _labels_at(d, k)
                ax = axes[gi, ci]
                same = labels[d["ei"]] == labels[d["ej"]]
                ecol = np.tile([0.7, 0.7, 0.72, 0.12], (d["ei"].size, 1))
                ecol[same] = pal(labels[d["ei"][same]]) * [1, 1, 1, 0.6]
                ax.add_collection(PathCollection(d["paths"], facecolors="none",
                                                 edgecolors=ecol, linewidths=0.5))
                ax.scatter(d["P"][:, 0], d["P"][:, 1], s=0.5 * d["node_size"],
                           c=pal(labels), edgecolors="white", linewidths=0.3, zorder=3)
                pad = 0.12
                ax.set_xlim(d["P"][:, 0].min() - pad, d["P"][:, 0].max() + pad)
                ax.set_ylim(d["P"][:, 1].min() - pad, d["P"][:, 1].max() + pad)
                ax.set_aspect("equal"); ax.axis("off")
                if gi == 0:
                    ax.set_title(f"{titles[ci]} ({labels.max() + 1} comm.)",
                                 fontsize=10)
                else:
                    ax.set_title(f"{labels.max() + 1} comm.", fontsize=9)
            axd = axes[gi, 3]
            dlc = LineCollection(d["dseg"], colors="0.35", linewidths=0.9)
            axd.add_collection(dlc)
            for k in d["ks_static"]:
                axd.axhline(cut_height(d, k), color="#c1121f", lw=0.9, ls="--")
            axd.set_yscale("log")
            axd.set_xlim(0, 10 * d["n"])
            axd.set_ylim(d["floor"] * 0.9, d["dh"].max() * 1.3)
            axd.set_xticks([]); axd.tick_params(labelsize=7)
            for s in ("top", "right"):
                axd.spines[s].set_visible(False)
            axd.set_ylabel(name, fontsize=9)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, transparent=False)
        plt.close(fig)
    print(f"done: {out_path.name}", flush=True)


def main():
    use_lrg_style()
    out_dir = Path(FIGURES_ROOT) / "lrg_diffusion_zoom"
    print("preparing graphs + cophenetic dendrograms ...", flush=True)
    data = prepare()
    make_static(data, out_dir / "lrg_dendrogram_cut_snapshots.pdf")
    make_animation(data, out_dir / "lrg_dendrogram_cut.gif")
    print(f"\nall assets in {out_dir}", flush=True)


if __name__ == "__main__":
    main()
