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

# --- dark-slide overrides (applied BEFORE build so legend + node ring both pick them up) ---
m.C_TRUTH = "#f4f6f9"        # ground-truth onset ring: near-white so it reads on a dark slide
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


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    P = m.compute()
    fig = _whiten(m.make_figure(P))
    fig.write_html(str(OUT) + ".html", include_plotlyjs="cdn")
    fig.write_image(str(OUT) + ".png", scale=2)
    print("talk fig:epi_seed_spread (dark) — heat spreads from 3 seeds to the distant SOZ")
    print(f"  {m.PATIENT} {m.BAND} {m.PHASE}: off-shaft discovery AUC = {P['auc']:.3f}, "
          f"top-flow precision = {P['prec']:.2f}")
    print(f"wrote {OUT}.png + .html")


if __name__ == "__main__":
    main()
