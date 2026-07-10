"""Reusable plotly primitives for 3D brain / electrode renderings.

Pure-geometry helpers shared by the pooled-implant "pairglow" brains (the talk
asset and the Results §1 hero). General names — no manuscript-local tokens —
so any 3D-brain figure can import them::

    from lrg_eegfc.visuals.brain3d import pial_mesh, spheres_mesh, bezier_arcs

Heavy deps (plotly, nilearn) live here and NOT in ``visuals/__init__`` so that
importing the package does not pull them in — only 3D-figure scripts do.
"""
from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np
import plotly.graph_objects as go


def unit_sphere(nlat: int = 7, nlon: int = 11) -> Tuple[np.ndarray, np.ndarray]:
    """Vertices / triangular faces of a low-poly unit sphere (for instancing)."""
    th = np.linspace(0, np.pi, nlat)
    ph = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    V = np.stack(
        [np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1
    ).reshape(-1, 3)
    F = []
    for a in range(nlat - 1):
        for b in range(nlon):
            b2 = (b + 1) % nlon
            v00, v01 = a * nlon + b, a * nlon + b2
            v10, v11 = (a + 1) * nlon + b, (a + 1) * nlon + b2
            F += [[v00, v10, v11], [v00, v11, v01]]
    return V, np.array(F)


_SPHV, _SPHF = unit_sphere()


def spheres_mesh(centers, radii, color, opacity: float = 1.0, name: str = "",
                 lighting: Optional[dict] = None) -> go.Mesh3d:
    """One ``Mesh3d`` of many mm-space spheres.

    Spheres live in data (brain-mm) space, so they grow when the viewer zooms
    in and never shrink away — the right behaviour for contacts inside a shell.
    """
    centers = np.asarray(centers, float)
    radii = np.asarray(radii, float)
    n, S = len(centers), len(_SPHV)
    V = (centers[:, None, :] + radii[:, None, None] * _SPHV[None, :, :]).reshape(-1, 3)
    F = (_SPHF[None, :, :] + (np.arange(n) * S)[:, None, None]).reshape(-1, 3)
    return go.Mesh3d(
        x=V[:, 0], y=V[:, 1], z=V[:, 2], i=F[:, 0], j=F[:, 1], k=F[:, 2],
        color=color, opacity=opacity, flatshading=False, hoverinfo="skip",
        lighting=lighting or dict(ambient=0.55, diffuse=0.7, specular=0.18,
                                  roughness=0.5),
        showscale=False, showlegend=False, name=name)


def bezier_arcs(coords, edges, center, n: int = 18, lift: float = 0.32):
    """None-separated quadratic-bezier polylines bowing outward from ``center``.

    Returns flat ``(xs, ys, zs)`` lists (``None`` between arcs) ready for a
    single ``Scatter3d(mode='lines')`` trace.
    """
    coords = np.asarray(coords, float)
    center = np.asarray(center, float)
    xs, ys, zs = [], [], []
    t = np.linspace(0, 1, n)
    for a, b in edges:
        pa, pb = coords[a], coords[b]
        mid = 0.5 * (pa + pb)
        out = mid - center
        nrm = np.linalg.norm(out) or 1.0
        ctrl = mid + (lift * np.linalg.norm(pb - pa) + 6.0) * out / nrm
        curve = (np.outer((1 - t) ** 2, pa) + np.outer(2 * (1 - t) * t, ctrl)
                 + np.outer(t ** 2, pb))
        xs += [*curve[:, 0], None]
        ys += [*curve[:, 1], None]
        zs += [*curve[:, 2], None]
    return xs, ys, zs


def pial_mesh(surf: str = "fsaverage5", color: str = "#d2d6dc",
              opacity: float = 0.05, lighting: Optional[dict] = None) -> go.Mesh3d:
    """Translucent bilateral pial shell (nilearn fsaverage) as a ``Mesh3d``."""
    from nilearn import datasets, surface
    fs = datasets.fetch_surf_fsaverage(surf)
    vl, fl = surface.load_surf_mesh(fs["pial_left"])
    vr, fr = surface.load_surf_mesh(fs["pial_right"])
    V = np.vstack([vl, vr])
    F = np.vstack([fl, fr + len(vl)])
    return go.Mesh3d(
        x=V[:, 0], y=V[:, 1], z=V[:, 2], i=F[:, 0], j=F[:, 1], k=F[:, 2],
        color=color, opacity=opacity, flatshading=False, hoverinfo="skip",
        lighting=lighting or dict(ambient=0.74, diffuse=0.20, specular=0.0),
        showscale=False, showlegend=False, name="pial")


# MNI/RAS anatomical poles: +x RIGHT  -x LEFT  +y FRONT  -y BACK  +z TOP  -z BOTTOM
# (verified via OFC y+38 / occipital y-60 / sensorimotor z+33). ``big`` flags the
# four in-plane cardinal labels rendered larger/darker than TOP/BOTTOM.
_AXIS_POLES = [(0, +1, "RIGHT", True), (0, -1, "LEFT", True),
               (1, +1, "FRONT", True), (1, -1, "BACK", True),
               (2, +1, "TOP", False), (2, -1, "BOTTOM", False)]


def orientation_traces(lo, hi, gap: float = 8.0, ext: float = 24.0):
    """Outward axis stubs + anatomical labels around a bbox ``[lo, hi]`` (MNI/RAS).

    Returns a list of ``Scatter3d`` traces (one line trace + two text traces) so a
    rotatable 3D brain is never disorienting: FRONT/BACK/LEFT/RIGHT/TOP/BOTTOM.
    Shared by every pooled-brain figure (trace heat, laterality, ...).
    """
    lo = np.asarray(lo, float)
    hi = np.asarray(hi, float)
    c = 0.5 * (lo + hi)
    sx, sy, sz = [], [], []
    pos, txt, big = [], [], []
    for ax, sgn, t, isbig in _AXIS_POLES:
        edge = hi[ax] if sgn > 0 else lo[ax]
        p0 = c.copy(); p0[ax] = edge + sgn * gap
        p1 = c.copy(); p1[ax] = edge + sgn * (gap + ext)
        sx += [p0[0], p1[0], None]; sy += [p0[1], p1[1], None]; sz += [p0[2], p1[2], None]
        lp = c.copy(); lp[ax] = edge + sgn * (gap + ext + 7)
        pos.append(lp); txt.append(t); big.append(isbig)
    traces = [go.Scatter3d(x=sx, y=sy, z=sz, mode="lines", hoverinfo="skip",
                           showlegend=False, line=dict(color="#6f767f", width=4))]
    pos = np.array(pos); big = np.array(big)
    for mask, size, col in [(big, 17, "#20242a"), (~big, 12, "#8a9098")]:
        i = np.where(mask)[0]
        traces.append(go.Scatter3d(
            x=pos[i, 0], y=pos[i, 1], z=pos[i, 2], mode="text",
            text=[txt[k] for k in i], hoverinfo="skip", showlegend=False,
            textfont=dict(size=size, color=col, family="Arial Black")))
    return traces
