#!/usr/bin/env python3
r"""talk_fig_epi_seed_spread_3d -- dark-slide wrap of the SOZ seed-spread brain (slide 17 HERO).

The 3D co-diffusing-community read, made for a dark Canva slide. Reuses the paper module's
compute() + make_figure() verbatim (Pat_08 delta rest_post: heat placed on 3 known onset
contacts on ONE electrode diffuses on the zero-lag-immune |ImCoh| graph and rediscovers the
UNSEEDED onset contacts on a DIFFERENT electrode; off-shaft AUC ~0.99) so no number can drift.
Only the presentation is flipped for dark: the near-black ground-truth ring -> near-white, and
every plotly font (legend, colorbar, orientation labels) -> white.

Reads : |ImCoh| FC (imcoh_abs, Pat_08 delta rest_post) + implant SOZ labels/coords.
Writes: data/outputs/figures/talk/fig_epi_seed_spread_3d.{html,png}  (transparent, dark)
"""
from __future__ import annotations
import importlib.util
import sys
from pathlib import Path

FE = Path(__file__).resolve().parents[1] / "01_compute" / "figures_embedded"
spec = importlib.util.spec_from_file_location("fig_epi_seed3d", FE / "fig_epi_seed_spread_3d.py")
m = importlib.util.module_from_spec(spec)
sys.modules["fig_epi_seed3d"] = m
spec.loader.exec_module(m)

# --- exemplar: the co-diffusing onset community (Pat_14 delta rest_post, 10 unseeded onset) ---
m.PATIENT, m.BAND, m.PHASE, m.HEMI = "Pat_14", "delta", "rest_post", "left"
m.CENTER = __import__("numpy").array([-38.0, -20.0, 8.0])   # left-hemi centroid, arcs bow out
m.ARC_MODE = "community"    # onset↔onset heat web (each contact tied to its strongest partners),
m.ARC_KNN = 2               #   NOT a star of spokes from the seeds
m.TOP_ARC = 7               # a few ✕ = healthy contacts that also couple in (false positives)

# --- dark-slide overrides (applied BEFORE build so legend + node ring both pick them up) ---
m.C_TRUTH = "#2fe36e"        # correctly-identified onset: drawn as a SOLID bright-green disc
#                              (see _solid_green_onset — plotly renders solid fills reliably but
#                              open-circle strokes almost invisibly). Green = onset the heat found,
#                              size still ∝ affinity; distinct from gold seeds and steel ✕ false alarms.
m.NODE_EDGE = "#11141a"      # keep a dark node edge (nodes are heat-coloured)
m.C_PIAL = "#767c86"         # dim pial shell so it recedes behind the data on a dark slide

OUT = m.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_epi_seed_spread_3d"


def _whiten(fig):
    """Flip text elements white for a dark background -- SURGICALLY.

    Only touch (a) the single real colorbar (the one marker with showscale=True) and
    (b) text traces. A global ``update_layout(font=...)`` activates the template's default
    coloraxis and paints a phantom plasma colorbar, so set fonts element-by-element instead.
    """
    fig.update_layout(legend=dict(font=dict(color="white")),
                      title=dict(font=dict(color="white")))
    for tr in fig.data:
        mk = getattr(tr, "marker", None)
        if mk is not None and getattr(mk, "showscale", None) is True:
            cb = mk.colorbar
            if cb.title is not None:
                cb.title.font = dict(color="white", size=13)
            cb.tickfont = dict(color="white", size=11)
        if getattr(tr, "mode", None) and "text" in tr.mode:
            tr.textfont = dict(color="white", size=13)
    return fig


def _zoom_and_declutter(fig, P, margin=20.0, msize=1.7):
    """Zoom the scene to the ELECTRODE cluster (not the whole brain) so the contacts read big,
    enlarge every node, and strip in-figure prose (legend off, colorbar title removed) -- the
    words live on the slide, per the talk brief."""
    import numpy as np
    c = P["coords"]
    lo, hi = c.min(0) - margin, c.max(0) + margin
    fig.update_scenes(
        xaxis=dict(visible=False, range=[lo[0], hi[0]]),
        yaxis=dict(visible=False, range=[lo[1], hi[1]]),
        zaxis=dict(visible=False, range=[lo[2], hi[2]]),
        aspectmode="data")
    for tr in fig.data:
        mk = getattr(tr, "marker", None)
        if mk is None or getattr(tr, "type", None) != "scatter3d":
            continue
        s = mk.size
        if s is None:
            continue
        mk.size = tuple(float(v) * msize for v in np.atleast_1d(s)) if np.ndim(s) else float(s) * msize
        if getattr(mk, "showscale", None) is True and mk.colorbar is not None:
            mk.colorbar.title = None                     # drop the multi-word colorbar title
            mk.colorbar.len = 0.42
            mk.colorbar.x = 0.97
    fig.update_layout(showlegend=False, margin=dict(l=0, r=0, t=8, b=8))
    return fig


def _fix_seed_count(fig, n_seed):
    """The paper legend hardcodes '(x3)'; rewrite it to the actual seed count for this patient."""
    for tr in fig.data:
        nm = getattr(tr, "name", None)
        if nm and nm.startswith("seed "):
            tr.name = f"seed — known onset (×{n_seed})"
    return fig


def _solid_green_onset(fig, fill="#2fe36e", edge="#0c1a10", grow=1.20):
    """Mark the correctly-identified onset with SOLID bright-green discs — not open rings.

    Plotly's 3D renderer draws solid marker fills reliably but ``circle-open`` strokes almost
    invisibly (line-width is effectively ignored on export), so the paper's ground-truth RING reads
    as an ultrathin line on the slide. Here we convert that real ring trace (transparent-fill
    ``circle-open``, not the hidden legend stub) into a bold filled green node with a dark edge,
    grown ~20% so it clearly stands off the heat field: unmistakable = "the heat found this onset".
    Size still ∝ seed-affinity (inherited from the trace), so a weakly-reached onset stays smaller."""
    import numpy as np
    for tr in fig.data:
        mk = getattr(tr, "marker", None)
        if mk is None or getattr(mk, "symbol", None) != "circle-open" or mk.line is None:
            continue
        xs = getattr(tr, "x", None)
        if xs is None or not any(v is not None for v in np.atleast_1d(xs)):
            continue                                            # skip the legend stub
        mk.symbol = "circle"                                     # SOLID fill (reliably rendered)
        mk.color = fill
        mk.line.color = edge
        mk.line.width = 1.6
        if mk.size is not None:
            mk.size = tuple(float(v) * grow for v in np.atleast_1d(mk.size))
    return fig


def _snap_inside_pial(P, hemi="left", inward=0.06):
    """Cosmetic-only: pull any contact that lands OUTSIDE the schematic template pial back to
    just inside its surface, so no pin appears to float outside the brain on the slide.

    The MNI mapping here is a heuristic affine (schematic), so a contact can project a few mm
    past the smoothed template pial. Coordinates are used ONLY for drawing — the diffusion runs
    on the FC graph, not on positions — so this moves no number. Each out-of-hull contact is sent
    to its nearest pial vertex, nudged ``inward`` of the way toward the pial centroid.
    """
    import numpy as np
    from scipy.spatial import Delaunay
    from fig_epi_seed3d import pial_mesh  # same module the figure draws from
    pm = pial_mesh(color="#888", opacity=0.1, hemi=hemi)
    V = np.column_stack([np.asarray(pm.x, float), np.asarray(pm.y, float), np.asarray(pm.z, float)])
    ctr = V.mean(0)
    c = P["coords"].copy()
    outside = Delaunay(V).find_simplex(c) < 0
    for i in np.where(outside)[0]:
        vn = V[np.linalg.norm(V - c[i], axis=1).argmin()]      # nearest pial vertex
        c[i] = vn * (1.0 - inward) + ctr * inward              # just inside the surface
    if outside.any():
        print(f"  snapped {int(outside.sum())} out-of-pial contact(s) just inside the shell "
              f"(cosmetic; coords are display-only)")
    P = dict(P); P["coords"] = c
    return P


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    P = _snap_inside_pial(m.compute(), hemi=m.HEMI)
    fig = _solid_green_onset(_zoom_and_declutter(_whiten(m.make_figure(P)), P, msize=1.35))
    fig.write_html(str(OUT) + ".html", include_plotlyjs="cdn")
    fig.write_image(str(OUT) + ".png", scale=2)
    print("talk fig:epi_seed_spread (dark, zoomed) — heat spreads from the seeds to the distant SOZ")
    print(f"  {m.PATIENT} {m.BAND} {m.PHASE}: off-shaft discovery AUC = {P['auc']:.3f}, "
          f"top-flow precision = {P['prec']:.2f}")
    print(f"wrote {OUT}.png + .html")


if __name__ == "__main__":
    main()
