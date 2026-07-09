#!/usr/bin/env python3
r"""fig:arc_ribbons_3d — the consolidation disentangling as six ribbons in space (§2, R2).

ONE 3D object, no flat panel. Each frequency band is a smooth ribbon flowing through
space; the parameter along the ribbon is the consolidation phase
(rest_pre -> task_learn -> task_test -> rest_post). The three spatial axes are:

    x  (left -> right)  = phase progression (the journey)
    y  (into the page)  = encoding      rho_S(Delta, e)
    z  (up)             = inference-specific   partial rho_S(Delta, f | e)

Reading, entirely from the SHAPE (no on-figure prose):
  * all six ribbons leave the origin together and swing up to a common inference ridge
    during the task  -> ENCODING / online engagement is broadband;
  * on the offline return (task_test -> rest_post) five ribbons DISSOLVE back toward the
    baseline floor (they fade to a ghost: no verified trace), while the beta ribbon stays
    vivid and lands on a small lifted node  -> only beta carries an offline trace.

Opacity on the return leg encodes the VERIFIED result, not raw height: a band's return
leg stays solid iff it clears its matched-strength rho_sym inference gate (beta, p=.010);
every other band's return leg fades. This is honest by construction -- the retained beta
coordinate is small in magnitude (~0.09) but is the only one that survives the null, so we
encode survival by vividness rather than by an exaggerated height.

Endpoints are the verified rho_sym waypoints (audit_163; rest_post == audit_152 exactly).
Splines interpolate only BETWEEN the four measured phase waypoints (marked as spheres);
the four spheres are the data.

Reads : data/audit/consolidation_arc_rhosym/arc_phase_trajectory.csv   (patient,band,phase,x,y)
        data/audit/consolidation_arc_rhosym/arc_cohort_verdict.csv     (inference gate per band)
Writes: data/reports/results_section2/fig_arc_ribbons_3d.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

ARC = ROOT / "data/audit/consolidation_arc_rhosym/arc_phase_trajectory.csv"
VERD = ROOT / "data/audit/consolidation_arc_rhosym/arc_cohort_verdict.csv"
OUT = ROOT / "data/reports/results_section2/fig_arc_ribbons_3d.pdf"

PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
TPH = {"rest_pre": 0.0, "task_learn": 1.0, "task_test": 2.0, "rest_post": 3.0}
BANDS = ["delta", "theta", "alpha", "low_gamma", "high_gamma", "beta"]  # beta last = on top
HERO = "beta"
COL = {
    "delta": "#2c4bd0", "theta": "#159a9a", "alpha": "#7b5cc4",
    "low_gamma": "#c04ab0", "high_gamma": "#8a8f98", "beta": "#ec1c24",
}
GLYPH = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "low_gamma": r"$\gamma_{L}$", "high_gamma": r"$\gamma_{H}$", "beta": r"$\beta$",
}
GATE_A = 0.05                       # inference gate threshold (matched-strength rho_sym)
NS = 140                            # spline samples along a ribbon


def _spline(tk, vk, ts):
    """Natural cubic through the 4 phase waypoints (tk=[0,1,2,3])."""
    from scipy.interpolate import CubicSpline
    return CubicSpline(tk, vk, bc_type="natural")(ts)


def _curve(sub):
    """(t, encoding, inference) dense curve from one band's 4 waypoints."""
    tk = np.array([TPH[p] for p in PHASES])
    xk = np.array([sub[sub.phase == p].x.iloc[0] for p in PHASES])   # encoding
    yk = np.array([sub[sub.phase == p].y.iloc[0] for p in PHASES])   # inference
    ts = np.linspace(0.0, 3.0, NS)
    return ts, _spline(tk, xk, ts), _spline(tk, yk, ts), (tk, xk, yk)


def _alpha_profile(ts, retained):
    """Per-vertex opacity: fade-in on the outbound leg; on the return leg stay solid
    if the band clears its gate, else dissolve toward a ghost."""
    a = np.where(ts <= 2.0, 0.45 + 0.55 * (ts / 2.0), 1.0)           # outbound fade-in
    ret = (ts - 2.0)                                                 # 0..1 on return leg
    if not retained:
        a = np.where(ts > 2.0, 1.0 - 0.80 * ret, a)                 # dissolve to 0.20
    return np.clip(a, 0.0, 1.0)


def _ribbon(ax, ts, xs, ys, zs, color, retained, bloom=True):
    """Glowing ribbon: soft bloom underlay + alpha-graded crisp core."""
    pts = np.column_stack([xs, ys, zs]).reshape(-1, 1, 3)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    a = _alpha_profile(ts, retained)
    aseg = 0.5 * (a[:-1] + a[1:])
    rgba = np.tile(np.array(to_rgba(color)), (len(segs), 1))
    rgba[:, 3] = aseg
    if bloom:                                                       # neon underlay
        m_out = ts <= 2.02
        ax.plot(xs[m_out], ys[m_out], zs[m_out], color=color, lw=8.5,
                alpha=0.06 if not retained else 0.09, solid_capstyle="round", zorder=2)
        if retained:                                               # hero glows throughout
            ax.plot(xs, ys, zs, color=color, lw=5.0, alpha=0.13,
                    solid_capstyle="round", zorder=2)
    lc = Line3DCollection(segs, colors=rgba, linewidths=3.2 if retained else 2.4)
    lc.set_capstyle("round")
    lc.set_zorder(5)
    ax.add_collection3d(lc)


def main():
    df = pd.read_csv(ARC)
    verd = pd.read_csv(VERD)
    gate = {b: float(verd[(verd.band == b) & (verd.functional == "T_infspec_pe")]
                      .wilcoxon_p.iloc[0]) for b in BANDS}
    med = (df.groupby(["band", "phase"], as_index=False)[["x", "y"]].median())

    fig = plt.figure(figsize=(12.6, 8.4))
    ax = fig.add_subplot(111, projection="3d")

    # ---- baseline floor at z=0 (anchors "height = inference") -------------------
    ymin, ymax = -0.06, 0.78
    floor = [[(-0.05, ymin, 0.0), (3.35, ymin, 0.0), (3.35, ymax, 0.0), (-0.05, ymax, 0.0)]]
    pc = Poly3DCollection(floor, facecolor="#c9ccd1", alpha=0.10, edgecolor="none")
    pc.set_zorder(0)
    ax.add_collection3d(pc)

    post_nodes = []
    for band in BANDS:
        sub = med[med.band == band]
        ts, xs_enc, ys_inf, (tk, xk, yk) = _curve(sub)
        # spatial mapping: X=phase, Y=encoding, Z=inference
        X, Y, Z = ts, xs_enc, ys_inf
        retained = gate[band] < GATE_A
        _ribbon(ax, ts, X, Y, Z, COL[band], retained)
        # phase waypoint spheres (the data); grow small->large along the journey
        a_wp = _alpha_profile(tk, retained)
        for i, p in enumerate(PHASES):
            ax.scatter(tk[i], xk[i], yk[i], s=34 + 44 * i,
                       color=COL[band], edgecolor="white", linewidth=1.1,
                       alpha=max(0.28, a_wp[i]), depthshade=False,
                       zorder=7 if retained else 6)
        post_nodes.append((band, xk[-1], yk[-1], retained))

    # ---- hero payoff: halo + dropline at beta's offline landing -----------------
    hb = next(n for n in post_nodes if n[0] == HERO)
    _, hx, hy, _ = hb
    ax.scatter([3.0], [hx], [hy], s=520, color=COL[HERO], alpha=0.16,
               edgecolor="none", depthshade=False, zorder=6)
    ax.plot([3.0, 3.0], [hx, hx], [hy, 0.0], color=COL[HERO], lw=1.1, ls=(0, (1, 2)),
            alpha=0.45, zorder=4)

    # ---- tip glyphs at the offline (rest_post) end ------------------------------
    for band, px, py, retained in post_nodes:
        ax.text(3.14, px, py + (0.035 if band == HERO else 0.0), GLYPH[band],
                color=COL[band], fontsize=17 if retained else 13,
                fontweight="bold" if retained else "normal",
                alpha=1.0 if retained else 0.6, ha="left", va="center", zorder=9)

    # ---- axes: minimal, clean, no box clutter -----------------------------------
    ax.set_xlim(-0.05, 3.4)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(-0.28, 0.66)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels([r"rest$_\mathrm{pre}$", r"task$_\mathrm{learn}$",
                        r"task$_\mathrm{test}$", r"rest$_\mathrm{post}$"], fontsize=11)
    ax.set_yticks([0.0, 0.3, 0.6])
    ax.set_zticks([-0.2, 0.0, 0.2, 0.4, 0.6])
    ax.tick_params(labelsize=9, pad=0.5)
    ax.set_ylabel("encoding", fontsize=12, labelpad=6)
    ax.set_zlabel("inference-specific", fontsize=12, labelpad=4)
    # transparent panes, faint grid
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((1, 1, 1, 0))
        axis.line.set_color((0.7, 0.7, 0.7, 0.5))
        axis._axinfo["grid"].update(color=(0.85, 0.85, 0.85, 0.5), linewidth=0.5)
    ax.view_init(elev=20, azim=-60)
    ax.set_box_aspect((2.35, 1.25, 1.05), zoom=1.06)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True, pad_inches=0.05)
    plt.close(fig)

    print("fig:arc_ribbons_3d — six-band consolidation ribbons\n")
    print(f"{'band':<12}{'post(enc,inf)':>18}{'inf gate p':>12}  {'return leg':>10}")
    for band in BANDS:
        sub = med[med.band == band]
        px = sub[sub.phase == 'rest_post'].x.iloc[0]
        py = sub[sub.phase == 'rest_post'].y.iloc[0]
        leg = "SOLID" if gate[band] < GATE_A else "fade"
        print(f"{band:<12}{f'({px:+.2f},{py:+.2f})':>18}{gate[band]:>12.3f}  {leg:>10}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
