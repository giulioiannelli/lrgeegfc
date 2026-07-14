#!/usr/bin/env python3
r"""talk_fig_scale_coarse_grains_space -- slide 14b: what the scale s DOES, on the brain.

THREE renormalization regimes as the diffusion scale s = tau*lambda_max grows:
  MICRO  (~O(n) clusters, avg ~2-3 nodes)  -- each contact ~its own population
  MESO   (~n/10 clusters, avg ~10 nodes, BALANCED groups)
  MACRO  (~2 groups)                        -- the implant as two integrated blocks
Communities = spectral clustering on the diffusion propagator rho = e^{-tau L} at scale
s (an INDEPENDENT partition of each operator -- different tau, different operator,
different partition; never one tree cut at three heights).

Graph: the DENSE |ImCoh| connectivity (NOT the mst@0.20 backbone). This is a deliberate
choice for a CONCEPT slide: only the dense graph yields BALANCED micro/meso/macro
(cv ~0.1-0.5) -- the sparse backbone is core-dominated (cv ~1.4-3, one ~50-node core
that never splits into balanced 10-node blocks) and cannot show a clean 3-regime
coarse-graining. Communities here are still reasonably local (within-community spatial
spread ~0.5 of the overall extent). Colour = size-rank at each scale (biggest -> same
hue across panels). Drawn ON the 3D brain: coloured contacts + the strongest intra-
community edges, arcs bowed inward so the contacts read on top.

Real data (Pat_05, largest bilateral implant; beta; rest_pre). Illustrative (scales /
levels chosen for clarity, not a statistic). No hulls, no dendrogram, NO text on the
asset -- three tight-cropped brains only; the slide carries the s-arrow and labels.

Writes: data/outputs/figures/talk/fig_scale_coarse_grains_space.png     (3-brain strip)
        data/outputs/figures/talk/_scale_brain_{micro,meso,macro}.png   (panels)
"""
from __future__ import annotations
import sys
from collections import defaultdict

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import plotly.graph_objects as go
from sklearn.cluster import SpectralClustering

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_D = ROOT / "scripts" / "01_compute" / "figures_embedded"
sys.path.insert(0, str(_D))
from fig_lrg_diffusion_zoom_demo import laplacian_eig            # noqa: E402

from lrg_eegfc.config.paths import SEEG_DATAPATH, FIGURES_ROOT   # noqa: E402
from lrg_eegfc.workflow.fc import load_fc_matrix                 # noqa: E402
from lrg_eegfc.visuals.brain3d import pial_mesh, spheres_mesh, bezier_arcs  # noqa: E402
from lrg_eegfc.visuals.spatial_coords import (                   # noqa: E402
    load_spatial_metadata, prepare_spatial_coordinates)
from lrg_eegfc.visuals.styles import use_lrg_style               # noqa: E402

use_lrg_style()

PAT, BAND, PHASE = "Pat_05", "beta", "rest_pre"
# three regimes (scale s = tau*lambda_max, k communities): micro / meso / macro
TAGS = ("micro", "meso", "macro")
SCALES = [(0.4, 50), (1.2, 12), (6.0, 2)]
N_EDGES = 420                # strongest |ImCoh| edges considered for the drawn network
SHELL, SHELL_OP = "#8c9299", 0.11    # grey glass shell
R_CONTACT = 2.4
GAP = 12                     # transparent px between the three tight-cropped brains
# arcs bow INWARD (into the volume, behind the contacts) so the nodes read on top
ARC_LIFT, ARC_BULGE = -0.25, -2.0
# anterosuperior-right 3/4 oblique -> reads as a 3D volume (was near-frontal)
CAM = dict(eye=dict(x=0.95, y=-1.5, z=0.8), up=dict(x=0, y=0, z=1),
           center=dict(x=0, y=-0.05, z=-0.05))
OUTDIR = FIGURES_ROOT / "talk"
OUT = OUTDIR / "fig_scale_coarse_grains_space.png"


def _rgba(c, a):
    r, g, b = (int(255 * x) for x in c[:3])
    return f"rgba({r},{g},{b},{a})"


def make_palette(nc):
    """nc distinct hues, ordered for a DARK background: bright warm/cyan/green colours
    lead (they colour the biggest communities), dark blues pushed to the end where only
    tiny micro-dots reach them."""
    from matplotlib.colors import to_rgb
    lead = ["#ff7f0e", "#17becf", "#2ca02c", "#e6194b", "#f1c40f", "#e377c2",
            "#00e0b0", "#ff6ec7", "#98df8a", "#ffbb78", "#c5b0d5", "#9edae5"]
    cols = [to_rgb(c) for c in lead]
    for name in ("tab20b", "tab20c", "tab20"):        # fill for the many micro-dots
        cm = plt.get_cmap(name)
        cols += [cm(i)[:3] for i in range(cm.N)]
    return np.array([cols[i % len(cols)] for i in range(nc)])


def load():
    """Dense |ImCoh| graph: Laplacian eig + strongest-N edges + mm coords."""
    A = np.asarray(load_fc_matrix(PAT, PHASE, BAND, "imcoh_abs"), float).copy()
    np.fill_diagonal(A, 0.0)
    A = np.clip(A, 0.0, None)
    n = A.shape[0]
    w, V = laplacian_eig(A)
    md = load_spatial_metadata(PAT, SEEG_DATAPATH)
    xyz = md[["x", "y", "z"]].to_numpy(float)
    fin = np.all(np.isfinite(xyz), axis=1)
    coords = np.full((n, 3), np.nan)
    coords[fin] = np.asarray(prepare_spatial_coordinates(
        md.loc[fin].reset_index(drop=True), scale="mm", center=False, to_mni=True), float)
    iu = np.triu_indices(n, 1)
    keep = [k for k in np.argsort(A[iu])[::-1] if fin[iu[0][k]] and fin[iu[1][k]]][:N_EDGES]
    edges = [(int(iu[0][k]), int(iu[1][k])) for k in keep]
    return A, w, V, coords, fin, n, edges


def partition(w, V, s, k):
    """k LOCAL communities: spectral clustering on the diffusion propagator at scale s.
    An independent partition of THIS operator (different tau -> different partition)."""
    tau = s / w.max()
    rho = (V * np.exp(-tau * w)) @ V.T
    rho = np.clip((rho + rho.T) / 2, 0.0, None)
    return SpectralClustering(n_clusters=k, affinity="precomputed",
                              assign_labels="discretize", random_state=0).fit_predict(rho)


def color_by_size(comm, palette):
    """Colour each community by its size-rank AT THIS SCALE: biggest -> palette[0].
    The growing core stays the same hue across scales (it is always the biggest),
    while the other communities read as distinct colours -- works for independent
    per-scale partitions where fine atoms do not nest into coarse ones."""
    labs, sizes = np.unique(comm, return_counts=True)
    order = np.argsort(-sizes, kind="stable")                 # biggest first
    rank = {int(labs[o]): r for r, o in enumerate(order)}
    return [palette[rank[int(c)]] for c in comm]


def render_panel(coords, fin, edges, comm, ncol, out_png):
    center = coords[fin].mean(0)
    fig = go.Figure()
    fig.add_trace(pial_mesh(color=SHELL, opacity=SHELL_OP))

    # strongest |ImCoh| edges, intra-community, in the cluster colour;
    # arcs bow inward so they sit behind the contacts
    groups = defaultdict(list)
    for i, j in edges:
        if fin[i] and fin[j] and comm[i] == comm[j]:
            groups[tuple(ncol[i])].append((i, j))
    for col, es in groups.items():
        xs, ys, zs = bezier_arcs(np.nan_to_num(coords), es, center, n=16,
                                 lift=ARC_LIFT, bulge=ARC_BULGE)
        fig.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", hoverinfo="skip",
                                   showlegend=False,
                                   line=dict(color=_rgba(col, 0.55), width=2.2)))
    by_col = defaultdict(list)
    for i in range(len(coords)):
        if fin[i]:
            by_col[tuple(ncol[i])].append(i)
    for col, idx in by_col.items():
        fig.add_trace(spheres_mesh(coords[idx], np.full(len(idx), R_CONTACT),
                                   _rgba(col, 1.0), opacity=1.0))
    fig.update_layout(
        scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False),
                   zaxis=dict(visible=False), aspectmode="data",
                   bgcolor="rgba(0,0,0,0)", camera=CAM),
        paper_bgcolor="rgba(0,0,0,0)", showlegend=False,
        margin=dict(l=0, r=0, t=0, b=0), width=760, height=760)
    fig.write_image(str(out_png), scale=2)


def _crop_alpha(img, pad=6):
    """Tight-crop an RGBA float image to its non-transparent bounding box."""
    a = img[..., 3] > 0.02
    ys, xs = np.where(a)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + 1 + pad, img.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + 1 + pad, img.shape[1])
    return img[y0:y1, x0:x1]


def compose(pngs, out):
    """Tight-crop each brain, centre-pad to a common box, hstack with a hair gap.
    No text, transparent background -- the slide carries the s-arrow and labels."""
    imgs = [_crop_alpha(plt.imread(p)) for p in pngs]
    H = max(im.shape[0] for im in imgs)
    W = max(im.shape[1] for im in imgs)
    gapcol = np.zeros((H, GAP, 4))
    parts = []
    for i, im in enumerate(imgs):
        canvas = np.zeros((H, W, 4))
        y0, x0 = (H - im.shape[0]) // 2, (W - im.shape[1]) // 2
        canvas[y0:y0 + im.shape[0], x0:x0 + im.shape[1]] = im
        parts.append(canvas)
        if i < len(imgs) - 1:
            parts.append(gapcol)
    plt.imsave(out, np.clip(np.concatenate(parts, axis=1), 0, 1))


def main():
    A, w, V, coords, fin, n, edges = load()
    palette = make_palette(max(k for _, k in SCALES))

    OUTDIR.mkdir(parents=True, exist_ok=True)
    pngs = []
    print(f"{PAT} {BAND} {PHASE}: n={n} placed={int(fin.sum())} edges={len(edges)}")
    for tag, (s, k) in zip(TAGS, SCALES):
        comm = partition(w, V, s, k)
        ncol = color_by_size(comm, palette)
        cnt = np.bincount(comm)
        sizes = sorted(cnt[cnt > 0], reverse=True)
        png = OUTDIR / f"_scale_brain_{tag}.png"
        render_panel(coords, fin, edges, comm, ncol, png)
        pngs.append(png)
        print(f"  {tag:5s} s={s:5.1f} k={k:2d}  n_comm={len(sizes)}  avg={n/len(sizes):.1f}  "
              f"sizes={[int(x) for x in sizes][:14]}")

    compose(pngs, OUT)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
