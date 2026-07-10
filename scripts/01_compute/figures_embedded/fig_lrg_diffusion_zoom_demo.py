#!/usr/bin/env python
"""LRG diffusion-zoom demo — information flow as diffusion time τ slides.

Synthetic illustration for the 20-min talk (I-3 -> I-4 -> M-3 pivot):
the heat-kernel propagator K(τ) = exp(-τ L̂) re-projected onto the links,
swept over diffusion time τ, on two graphs that share the SAME latent
3-level hierarchy:

    (A) a hierarchical modular network (HMN) — sparse, block-structured,
        synthetic, with a KNOWN ground-truth hierarchy;
    (B) a REAL brain functional-connectivity matrix — imcoh_abs, one
        patient / band / phase — a COMPLETE weighted graph (~100% dense, no
        sparse blocks), the correlator-operator / Villegas "outlier" regime.

The point of the figure: as τ grows the propagator resolves nested
communities level by level. Panel (A) VALIDATES that LRG recovers a known
hierarchy; panel (B) APPLIES the identical operator to real brain FC and
reveals multiscale communities inside a dense weighted correlator — so LRG
is not confined to sparse modular graphs. The real matrix is reordered by
its own LRG dendrogram (UPGMA on the communication distance at τ=1/λ_max)
so the discovered blocks read on the diagonal. Combinatorial Laplacian only
(project rule).

Outputs (to data/outputs/figures/lrg_diffusion_zoom/):
    lrg_diffusion_zoom.gif             animated τ-sweep (talk asset)
    lrg_diffusion_zoom_snapshots.pdf   static 3-τ contact sheet (vector)

Run (inside the lapbrain env):
    python scripts/01_compute/figures_embedded/fig_lrg_diffusion_zoom_demo.py
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc_context
from matplotlib.collections import LineCollection
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.linalg import eigh
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import squareform

from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.visuals.styles import use_lrg_style

# ----------------------------------------------------------------------------
# geometry of the shared latent hierarchy
# ----------------------------------------------------------------------------
N = 64                    # 64 nodes = 4 mid-blocks x 4 fine-blocks x 4 nodes
FINE = 4                  # nodes per finest block  -> 16 fine blocks
MID = 16                  # nodes per mid block      ->  4 mid blocks
SEED = 7

fine_id = np.arange(N) // FINE      # 0..15
mid_id = np.arange(N) // MID        # 0..3


def _pair_level(i: int, j: int) -> int:
    """0 = same fine block, 1 = same mid block, 2 = only same graph (coarse)."""
    if fine_id[i] == fine_id[j]:
        return 0
    if mid_id[i] == mid_id[j]:
        return 1
    return 2


# base coupling strength by hierarchical level (decreasing with distance)
BASE_W = np.array([1.00, 0.35, 0.12])


# ----------------------------------------------------------------------------
# graph builders
# ----------------------------------------------------------------------------
def build_hmn(rng: np.random.Generator) -> np.ndarray:
    """Sparse hierarchical modular network: dense within fine blocks, rarer
    across, with a guaranteed connected backbone."""
    p_level = np.array([0.90, 0.16, 0.025])  # connection prob per level (sparse)
    A = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            lv = _pair_level(i, j)
            if rng.random() < p_level[lv]:
                w = BASE_W[lv] * (0.7 + 0.6 * rng.random())
                A[i, j] = A[j, i] = w
    # guarantee connectivity: bridge any disconnected components weakly
    G = nx.from_numpy_array(A)
    comps = list(nx.connected_components(G))
    for a, b in zip(comps[:-1], comps[1:]):
        i, j = min(a), min(b)
        A[i, j] = A[j, i] = BASE_W[2] * 0.5
    return A


def build_fc_heterogeneous(rng: np.random.Generator) -> np.ndarray:
    """Complete graph, every pair coupled, heterogeneous node strengths.
    The 3-level hierarchy survives only as a *gradient* buried in a dense
    weighted matrix — the brain-FC regime."""
    s = rng.lognormal(mean=0.0, sigma=0.5, size=N)   # heterogeneous strengths
    A = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            lv = _pair_level(i, j)
            w = BASE_W[lv] * np.sqrt(s[i] * s[j]) * (0.6 + 0.8 * rng.random())
            A[i, j] = A[j, i] = w
    return A


# ----------------------------------------------------------------------------
# LRG primitive: combinatorial Laplacian + heat-kernel propagator
# ----------------------------------------------------------------------------
def laplacian_eig(A: np.ndarray):
    """Combinatorial L̂ = D̂ - Â, eigendecomposed ONCE (symmetric PSD)."""
    L = np.diag(A.sum(1)) - A
    w, V = eigh(L)
    w = np.clip(w, 0.0, None)
    return w, V


def propagator(w: np.ndarray, V: np.ndarray, tau: float) -> np.ndarray:
    """K(τ) = exp(-τ L̂) = V diag(e^{-τλ}) Vᵀ  (entrywise >= 0: -L̂ is Metzler)."""
    K = (V * np.exp(-tau * w)[None, :]) @ V.T
    np.fill_diagonal(K, 0.0)            # keep the OFF-diagonal information flow
    return np.clip(K, 0.0, None)


def tau_bounds(w: np.ndarray):
    """τ range tied to the spectrum: t_min sees nearest neighbours only,
    t_max sits at the coarse-community (Fiedler) scale — stopping BEFORE
    full uniform mixing so the coarsest communities stay visible."""
    wpos = w[w > 1e-9]
    lam2, lam_max = wpos.min(), w.max()
    return 0.15 / lam_max, 1.6 / lam2


def tau_at(t_min: float, t_max: float, s) -> np.ndarray:
    """Log-interpolate τ for zoom fraction s in [0, 1] (0 = fine, 1 = coarse)."""
    return t_min * (t_max / t_min) ** np.asarray(s, float)


def dendro_order_and_communities(w: np.ndarray, V: np.ndarray, k: int):
    """LRG dendrogram (UPGMA on communication distance D=1/K at τ=1/λ_max):
    return the leaf order (to reveal blocks on the diagonal) and k community
    labels relabelled by first appearance along that order."""
    K = propagator(w, V, 1.0 / w.max())
    with np.errstate(divide="ignore"):
        D = 1.0 / K
    D[~np.isfinite(D)] = np.nanmax(D[np.isfinite(D)])
    np.fill_diagonal(D, 0.0)
    D = 0.5 * (D + D.T)
    Z = linkage(squareform(D, checks=False), method="average")
    order = leaves_list(Z)
    raw = fcluster(Z, t=k, criterion="maxclust")
    remap, labels = {}, np.zeros(len(raw), int)
    for idx in order:
        remap.setdefault(raw[idx], len(remap))
        labels[idx] = remap[raw[idx]]
    return order, labels


def _cat_cmap(ncat: int) -> mpl.colors.ListedColormap:
    base = plt.get_cmap("tab10")
    return mpl.colors.ListedColormap([base(i % 10) for i in range(max(ncat, 1))])


def _circular_layout(order: np.ndarray) -> np.ndarray:
    """Nodes evenly on a circle in dendrogram-leaf order — communities become
    contiguous coloured arcs and edges become chords: a clean, balanced layout
    for a DENSE connectome, where force-directed collapses to a hairball and a
    ring-of-blobs is lopsided when one module dominates."""
    n = order.size
    ang = np.empty(n)
    ang[order] = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return np.c_[np.cos(ang), np.sin(ang)]


# ----------------------------------------------------------------------------
# assembly — one synthetic HMN (known hierarchy) + one REAL brain FC matrix
# ----------------------------------------------------------------------------
REAL = dict(patient="Pat_05", band="beta", phase="rest_pre", k=6)

GRAPHS = ("Hierarchical modular network (synthetic)",
          f"Brain FC · {REAL['band']} · {REAL['patient']} {REAL['phase']}")


def _process_graph(A: np.ndarray, is_real: bool, k) -> dict:
    """Per-graph prep: Laplacian eig, spring layout, node communities +
    matrix reorder, own edge set. For the real graph, communities and the
    diagonal-block ordering come from its LRG dendrogram; for the synthetic
    HMN they are the known ground-truth mid-blocks in natural order."""
    n = A.shape[0]
    w, V = laplacian_eig(A)
    t_min, t_max = tau_bounds(w)
    ei, ej = np.where(np.triu(A, 1) > 0)      # OWN links: sparse HMN / dense FC
    if is_real:
        # dense complete graph → circular layout in dendrogram order
        # (spring → hairball); communities read as contiguous coloured arcs
        order, labels = dendro_order_and_communities(w, V, k)
        P = _circular_layout(order)
        ncmap = _cat_cmap(int(labels.max()) + 1)
    else:
        # sparse structured graph → force-directed reads well
        order, labels, ncmap = np.arange(n), mid_id, NODE_CMAP
        pos = nx.spring_layout(nx.from_numpy_array(A), weight="weight",
                               seed=SEED, k=1.6 / np.sqrt(n), iterations=300)
        P = np.array([pos[i] for i in range(n)])
    iu, ju = np.triu_indices(n, 1)
    return dict(
        A=A, w=w, V=V, n=n, P=P, iu=iu, ju=ju,
        ei=ei, ej=ej, n_edges=int(ei.size),
        order=order, node_rgba=ncmap(labels),
        node_size=min(30.0, 2200.0 / n),
        segs=np.stack([P[ei], P[ej]], axis=1),
        taus=tau_at(t_min, t_max, np.linspace(0.0, 1.0, NFRAMES)),
        taus_static=tau_at(t_min, t_max, np.array([0.18, 0.52, 0.86])))


def prepare():
    rng = np.random.default_rng(SEED)
    A_hmn = build_hmn(rng)
    try:
        A_real = np.asarray(load_fc_matrix(
            REAL["patient"], REAL["phase"], REAL["band"], "imcoh_abs"), float).copy()
        np.fill_diagonal(A_real, 0.0)
        A_real = np.clip(A_real, 0.0, None)
        print(f"  loaded real FC {REAL['patient']} {REAL['band']} "
              f"{REAL['phase']}: {A_real.shape[0]} nodes", flush=True)
    except Exception as exc:                                    # pragma: no cover
        print(f"  (real FC load failed: {exc}; synthetic fallback)", flush=True)
        A_real = build_fc_heterogeneous(rng)
    return {GRAPHS[0]: _process_graph(A_hmn, False, None),
            GRAPHS[1]: _process_graph(A_real, True, REAL["k"])}


# node colours anchor the ground-truth mid-blocks (4 groups)
NODE_CMAP = mpl.colors.ListedColormap(["#4cc9f0", "#f72585", "#ffd166", "#8ac926"])
EDGE_CMAP = plt.get_cmap("magma")
NFRAMES = 60
FPS = 18


# ----------------------------------------------------------------------------
# animation
# ----------------------------------------------------------------------------
def make_animation(data, out_path: Path):
    dark = {
        "figure.facecolor": "#0b0e17", "axes.facecolor": "#0b0e17",
        "savefig.facecolor": "#0b0e17", "text.color": "#e6edf3",
        "axes.edgecolor": "#30363d", "axes.labelcolor": "#e6edf3",
    }
    with rc_context(dark):
        fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.4), dpi=90)
        fig.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.10,
                            wspace=0.06, hspace=0.20)
        artists = {}
        for r, name in enumerate(GRAPHS):
            d = data[name]
            ax_net, ax_mat = axes[r, 0], axes[r, 1]
            # --- network (propagator re-projected on the graph's OWN links) ---
            Pg = d["P"]
            lc = LineCollection(d["segs"], linewidths=0.0, capstyle="round")
            ax_net.add_collection(lc)
            ax_net.scatter(Pg[:, 0], Pg[:, 1], s=d["node_size"], c=d["node_rgba"],
                           edgecolors="#0b0e17", linewidths=0.4, zorder=3)
            pad = 0.08
            ax_net.set_xlim(Pg[:, 0].min() - pad, Pg[:, 0].max() + pad)
            ax_net.set_ylim(Pg[:, 1].min() - pad, Pg[:, 1].max() + pad)
            ax_net.set_aspect("equal"); ax_net.axis("off")
            # --- matrix K(τ) (reordered by the graph's own community order) ---
            im = ax_mat.imshow(np.zeros((d["n"], d["n"])), cmap=EDGE_CMAP, vmin=0,
                               vmax=1, interpolation="nearest", animated=True)
            ax_mat.set_xticks([]); ax_mat.set_yticks([])
            for s in ax_mat.spines.values():
                s.set_color("#30363d")
            ax_net.set_title(name, color="#e6edf3", fontsize=11, pad=4,
                             loc="left")
            artists[name] = dict(lc=lc, im=im)
        # column headers (once, figure-level — no collision with the τ title)
        fig.text(0.25, 0.915, "propagator on the links", ha="center",
                 color="#8b949e", fontsize=9)
        fig.text(0.747, 0.915, "K(τ) = e^{-τ L̂}", ha="center",
                 color="#8b949e", fontsize=9)
        # τ progress bar + label spanning the bottom
        cax = fig.add_axes([0.28, 0.035, 0.44, 0.016])
        cax.set_xlim(0, 1); cax.set_ylim(0, 1); cax.axis("off")
        cax.axhline(0.5, color="#30363d", lw=2)
        marker, = cax.plot([0], [0.5], "o", color="#f72585", ms=8, zorder=4)
        cax.text(0.0, 1.9, "fine · local", color="#8b949e", fontsize=8, ha="left")
        cax.text(1.0, 1.9, "coarse · global", color="#8b949e", fontsize=8,
                 ha="right")
        tau_txt = fig.text(0.5, 0.965, "", ha="center", color="#e6edf3",
                           fontsize=12, fontweight="bold")

        t0 = time.time()

        def update(frame):
            s = frame / (NFRAMES - 1)
            changed = [marker, tau_txt]
            for name in GRAPHS:
                d, art = data[name], artists[name]
                tau = d["taus"][frame]
                K = propagator(d["w"], d["V"], tau)
                vmax = max(np.percentile(K[d["iu"], d["ju"]], 99.0), 1e-9)
                val = np.clip(K[d["ei"], d["ej"]] / vmax, 0, 1)   # flow on own edges
                rgba = EDGE_CMAP(val)
                rgba[:, 3] = val ** 0.75              # fade weak links
                art["lc"].set_color(rgba)
                art["lc"].set_linewidth(0.15 + 2.6 * val)
                art["im"].set_data(K[np.ix_(d["order"], d["order"])])
                art["im"].set_clim(0, vmax)
                changed += [art["lc"], art["im"]]
            marker.set_data([s], [0.5])
            tau_txt.set_text("diffusion time  τ  —  scale " f"{s:4.0%}")
            if frame % 12 == 0 or frame == NFRAMES - 1:
                el = time.time() - t0
                print(f"  [render {frame + 1:2d}/{NFRAMES}] s={s:4.2f} "
                      f"elapsed={el:5.1f}s", flush=True)
            return changed

        anim = FuncAnimation(fig, update, frames=NFRAMES, blit=False)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"writing GIF -> {out_path}", flush=True)
        anim.save(out_path, writer=PillowWriter(fps=FPS))
        plt.close(fig)
    _optimize_gif(out_path)
    size_mb = out_path.stat().st_size / 1e6
    print(f"done: {out_path.name}  ({size_mb:.1f} MB, {NFRAMES} frames @ {FPS} fps)",
          flush=True)


def _optimize_gif(path: Path):
    """Re-encode with a single shared 96-colour palette (no dither) to roughly
    halve the file size. Temp-file safe: a failure never corrupts the good GIF."""
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
        pal = rgb[0].quantize(colors=96, method=Image.MEDIANCUT)
        pframes = []
        for f in rgb:
            pf = f.quantize(palette=pal, dither=Image.NONE)
            pf.info.pop("transparency", None)        # avoid P-mode save bug
            pframes.append(pf)
        tmp = path.with_suffix(".opt.gif")
        pframes[0].save(tmp, save_all=True, append_images=pframes[1:],
                        duration=int(1000 / FPS), loop=0, optimize=True)
        if tmp.stat().st_size < path.stat().st_size:
            tmp.replace(path)
        else:
            tmp.unlink()
    except Exception as exc:                         # pragma: no cover
        print(f"  (gif optimize skipped: {exc})", flush=True)


# ----------------------------------------------------------------------------
# static contact sheet (vector PDF, light theme) — the M-3 figure
# ----------------------------------------------------------------------------
def make_static(data, out_path: Path):
    ecmap = plt.get_cmap("turbo")
    labels = ("fine (τ small)", "intermediate", "coarse (τ large)")
    with rc_context({"figure.facecolor": "white", "savefig.facecolor": "white"}):
        fig, axes = plt.subplots(4, 3, figsize=(8.6, 11.0))
        fig.subplots_adjust(left=0.09, right=0.98, top=0.95, bottom=0.03,
                            wspace=0.08, hspace=0.14)
        row_of = {(0, "net"): 0, (0, "mat"): 1, (1, "net"): 2, (1, "mat"): 3}
        for gi, name in enumerate(GRAPHS):
            d = data[name]
            Pg = d["P"]
            for ci, tau in enumerate(d["taus_static"]):
                K = propagator(d["w"], d["V"], tau)
                vmax = max(np.percentile(K[d["iu"], d["ju"]], 99.0), 1e-9)
                val = np.clip(K[d["ei"], d["ej"]] / vmax, 0, 1)   # flow on own edges
                # network
                axn = axes[row_of[(gi, "net")], ci]
                rgba = ecmap(val); rgba[:, 3] = val ** 0.8
                lc = LineCollection(d["segs"], colors=rgba,
                                    linewidths=0.1 + 2.2 * val, capstyle="round")
                axn.add_collection(lc)
                axn.scatter(Pg[:, 0], Pg[:, 1], s=0.55 * d["node_size"],
                            c=d["node_rgba"], edgecolors="white", linewidths=0.3,
                            zorder=3)
                pad = 0.08
                axn.set_xlim(Pg[:, 0].min() - pad, Pg[:, 0].max() + pad)
                axn.set_ylim(Pg[:, 1].min() - pad, Pg[:, 1].max() + pad)
                axn.set_aspect("equal"); axn.axis("off")
                # matrix (reordered by the graph's own community order)
                axm = axes[row_of[(gi, "mat")], ci]
                axm.imshow(K[np.ix_(d["order"], d["order"])], cmap=ecmap,
                           vmin=0, vmax=vmax, interpolation="nearest")
                axm.set_xticks([]); axm.set_yticks([])
                if row_of[(gi, "net")] == 0:
                    axn.set_title(labels[ci], fontsize=11)
            axes[row_of[(gi, "net")], 0].set_ylabel(
                name, fontsize=10.5, rotation=90, labelpad=12)
            axes[row_of[(gi, "net")], 0].axis("on")
            axes[row_of[(gi, "net")], 0].set_xticks([])
            axes[row_of[(gi, "net")], 0].set_yticks([])
            for sp in axes[row_of[(gi, "net")], 0].spines.values():
                sp.set_visible(False)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, transparent=False)
        plt.close(fig)
    print(f"done: {out_path.name}", flush=True)


def main():
    use_lrg_style()
    out_dir = Path(FIGURES_ROOT) / "lrg_diffusion_zoom"
    print("preparing graphs + eigendecompositions ...", flush=True)
    data = prepare()
    for name in GRAPHS:
        d = data[name]
        w = d["w"]
        n_pairs = d["n"] * (d["n"] - 1) // 2
        print(f"  {name}: N={d['n']}  edges={d['n_edges']}/{n_pairs} "
              f"({100 * d['n_edges'] / n_pairs:.0f}% density)  "
              f"λ2={w[w>1e-9].min():.4f}  λmax={w.max():.2f}", flush=True)
    make_static(data, out_dir / "lrg_diffusion_zoom_snapshots.pdf")
    make_animation(data, out_dir / "lrg_diffusion_zoom.gif")
    print(f"\nall assets in {out_dir}", flush=True)


if __name__ == "__main__":
    main()
