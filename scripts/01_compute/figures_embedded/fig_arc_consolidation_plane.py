#!/usr/bin/env python3
r"""fig:arc_consolidation_plane — the whole-cohort consolidation arc as a light-trail flow (§2).

CORE MESSAGE (the inference-specific arc, all 10 patients at once): every patient is shown
the premises (a leap out along the ENCODING axis), then reasons (a drift up along the
INFERENCE axis); at offline rest the cohort SPLITS — the six who clear the strength-matched
beta inference null HOLD their displaced state up along inference (warm, stay bright), the
rest RETURN to baseline (cool, fade out). "Held, not replayed", drawn as glowing trails.

    x = encoding axis  = state . (task_learn - rest_pre)         (leap when shown premises)
    y = inference axis = state . (task_test - task_learn)_|_x    (drift while reasoning)
    rest_pre = origin (shared baseline); rest_post above y=0 = inference retained.

2-D COHORT companion to the 3-D single-exemplar attractor (fig_arc_attractor_3d): same
encoding/inference basis, but all ten four-phase arcs, static, talk/publication-ready.

HONESTY (same register as the attractor): the arc geometry is the ILLUSTRATIVE projection
of the phase LRG states; the load-bearing statistic is the matched-strength inference gate
T_infspec_pe (beta p=0.010, 6/10), which is what the trail COLOUR encodes and what the
decomposition forest (Fig 3b) reports. No new claim rests on the geometry. Dark background is
a deliberate aesthetic (talk hero); rendered opaque (not the default transparent).

Reads : data/audit/consolidation_arc_rhosym/arc_phase_trajectory.csv  (x,y per phase)
        data/audit/consolidation_arc_rhosym/arc_null_per_patient.csv   (T_infspec_pe_p)
Writes: data/preprint/figures/_drafts/fig_arc_consolidation_plane_{band}_DRAFT.pdf
"""
from __future__ import annotations

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

from lrg_eegfc.visuals.styles import use_lrg_style

ARC = ROOT / "data/audit/consolidation_arc_rhosym"
OUTDIR = ROOT / "data/preprint/figures/_drafts"
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
GATE = 0.05

BG = "#0c111c"                                   # deep navy background
INK = "#c9d2e0"                                  # light text/axes on dark
FAINT = "#3a4560"                                # faint guide lines
# trail colourmaps: dark(rest) -> bright(reasoning/rest); warm=held, cool=returns
CM_TRACER = LinearSegmentedColormap.from_list("tr", ["#5a2f0c", "#c9781f", "#ffca5e"])
CM_RESET = LinearSegmentedColormap.from_list("re", ["#123a46", "#2f7f93", "#8fe0ec"])
GLOW_TR, GLOW_RE = to_rgb("#ff9e2c"), to_rgb("#59c6d8")
PHASE_LABEL = {"rest_pre": "rest", "task_learn": "shown\npremises",
               "task_test": "reasoning", "rest_post": "offline rest"}


def load(band):
    tr = pd.read_csv(ARC / "arc_phase_trajectory.csv")
    nn = pd.read_csv(ARC / "arc_null_per_patient.csv")
    tr = tr[tr.band == band]
    pval = (nn[nn.band == band].set_index("patient").T_infspec_pe_p).to_dict()
    return {pat: dict(xy=g.set_index("phase").loc[PHASES][["x", "y"]].to_numpy(float),
                      p=float(pval[pat]))
            for pat, g in tr.groupby("patient")}


def smooth(xy, n_per=70, tension=0.72):
    """Cardinal (Hermite) spline through the four phase points; local + interpolating."""
    P = np.asarray(xy, float)
    P = np.vstack([P[0] + (P[0] - P[1]), P, P[-1] + (P[-1] - P[-2])])
    m = (1 - tension) * (P[2:] - P[:-2]) / 2.0
    t = np.linspace(0, 1, n_per)
    h = (2 * t**3 - 3 * t**2 + 1, t**3 - 2 * t**2 + t,
         -2 * t**3 + 3 * t**2, t**3 - t**2)
    out = [np.outer(h[0], P[i]) + np.outer(h[1], m[i - 1])
           + np.outer(h[2], P[i + 1]) + np.outer(h[3], m[i])
           for i in range(1, len(P) - 2)]
    return np.vstack(out)


def glow_arc(ax, xy, cmap, glow_rgb, held, z=3):
    """One glowing trail: soft wide halo + a colour-gradient core whose alpha rises over
    time. HELD trails stay bright to the end; RETURNED trails are recessive and fade, so
    the amber 'held' story reads on top of a quiet cool substrate."""
    p = smooth(xy)
    seg = np.stack([p[:-1], p[1:]], axis=1)
    n = len(seg)
    halo = [(7.0, 0.05), (3.6, 0.09)] if held else [(5.0, 0.03), (2.8, 0.05)]
    for lw, a in halo:
        ax.add_collection(LineCollection(seg, colors=[(*glow_rgb, a)] * n,
                                         linewidths=lw, capstyle="round", zorder=z))
    t = np.linspace(0, 1, n)
    cols = cmap(t)
    if held:
        cols[:, 3] = np.clip(np.linspace(0.12, 0.98, n), 0, 1)
        lw_core = 1.9
    else:                                                          # recede + fade
        fade = np.clip(1 - np.maximum(0, t - 0.5) / 0.5 * 0.8, 0.18, 1)
        cols[:, 3] = np.clip(np.linspace(0.07, 0.55, n) * fade, 0, 1)
        lw_core = 1.3
    ax.add_collection(LineCollection(seg, colors=cols, linewidths=lw_core,
                                     capstyle="round", zorder=z + 1))
    return p


def endpoint(ax, xy, glow_rgb, held, z=8):
    r = xy[-1]
    ax.scatter([r[0]], [r[1]], s=300 if held else 150,
               color=[(*glow_rgb, 0.18 if held else 0.10)], zorder=z, edgecolors="none")
    ax.scatter([r[0]], [r[1]], s=110 if held else 46,
               color=[(*glow_rgb, 1.0 if held else 0.55)],
               edgecolors=BG, linewidths=1.2, zorder=z + 1)


def axis_arrow(ax, lo, hi):
    """Faint encoding-> / inference-> reference arrows from the origin."""
    ax.add_patch(FancyArrowPatch((0, 0), (hi[0] * 0.96, 0), arrowstyle="-|>",
                 mutation_scale=13, color=FAINT, lw=1.1, zorder=1))
    ax.add_patch(FancyArrowPatch((0, 0), (0, hi[1] * 0.96), arrowstyle="-|>",
                 mutation_scale=13, color=FAINT, lw=1.1, zorder=1))
    ax.text(hi[0] * 0.96, -0.03 * (hi[1] - lo[1]), "encoding →", color=INK,
            fontsize=11, ha="right", va="top", alpha=0.85)
    ax.text(0.02 * (hi[0] - lo[0]), hi[1] * 0.96, "inference →", color=INK,
            fontsize=11, ha="left", va="top", rotation=90, alpha=0.85)


def build(band):
    use_lrg_style()
    data = load(band)
    fig, ax = plt.subplots(figsize=(7.4, 6.7))
    fig.patch.set_facecolor(BG); ax.set_facecolor(BG)

    ax.axhline(0, color=FAINT, lw=0.8, alpha=0.6, zorder=0)

    for d in data.values():                                       # trails, held last (on top)
        held = d["p"] < GATE
        glow_arc(ax, d["xy"], CM_TRACER if held else CM_RESET,
                 GLOW_TR if held else GLOW_RE, held, z=4 if held else 3)
    for d in data.values():
        endpoint(ax, d["xy"], GLOW_TR if d["p"] < GATE else GLOW_RE, d["p"] < GATE)

    # cohort-median arc: bright neutral ribbon with phase labels + arrowhead
    med = np.stack([np.median([data[p]["xy"][i] for p in data], axis=0)
                    for i in range(4)])
    mp = smooth(med)
    mseg = np.stack([mp[:-1], mp[1:]], axis=1)
    ax.add_collection(LineCollection(mseg, colors=[(*to_rgb("#f2efe6"), 0.10)] * len(mseg),
                                     linewidths=6.5, capstyle="round", zorder=6))
    ax.add_collection(LineCollection(mseg, colors=[(*to_rgb("#f7f3e8"), 0.92)] * len(mseg),
                                     linewidths=2.6, capstyle="round", zorder=7))
    ax.add_patch(FancyArrowPatch(mp[-7], mp[-1], arrowstyle="-|>", mutation_scale=17,
                                 color="#f7f3e8", lw=0, zorder=8))
    off = {"rest_pre": (-0.03, -0.05, "right", "top"),
           "task_learn": (0.0, -0.10, "center", "top"),
           "task_test": (0.0, 0.09, "center", "bottom"),
           "rest_post": (-0.10, 0.05, "right", "bottom")}
    for i, ph in enumerate(PHASES):
        x, y = med[i]; dx, dy, ha, va = off[ph]
        ax.scatter([x], [y], s=26, color="#f7f3e8", zorder=8)
        ax.annotate(PHASE_LABEL[ph], (x, y), (x + dx, y + dy), ha=ha, va=va,
                    fontsize=10.5, color="#f2efe6", fontweight="bold",
                    linespacing=0.9, zorder=9)

    xy_all = np.vstack([d["xy"] for d in data.values()])
    lo, hi = xy_all.min(0), xy_all.max(0)
    axis_arrow(ax, lo, hi)
    ax.set_xlim(lo[0] - 0.12, hi[0] + 0.14)
    ax.set_ylim(lo[1] - 0.12, hi[1] + 0.12)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    n_tr = sum(d["p"] < GATE for d in data.values())
    handles = [Line2D([], [], color="#ffca5e", lw=3, label=f"holds inference  ({n_tr}/10)"),
               Line2D([], [], color="#8fe0ec", lw=3, label=f"returns to baseline  ({10 - n_tr}/10)")]
    leg = ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=10,
                    labelcolor=INK, handlelength=1.7)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / f"fig_arc_consolidation_plane_{band}_DRAFT.pdf"
    fig.savefig(out, bbox_inches="tight", facecolor=BG, pad_inches=0.08)
    qa = "/tmp/claude-1000/-home-giulio-Documents-research-neural-networks-lrgeegfc/c2197e94-c496-4e4e-97cc-89e43564e7da/scratchpad"
    fig.savefig(f"{qa}/arcplane_{band}.png", dpi=150, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"wrote {out.name}  tracers={n_tr}/10")


if __name__ == "__main__":
    for b in ["beta", "alpha"]:
        build(b)
