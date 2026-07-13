#!/usr/bin/env python3
r"""fig:arc_attractor — the consolidation trace as a live state-space trajectory (§2, R2).

Concatenate the four phases into one long recording, slide a 30 s window across it (all
rhythms resolved), read the LRG cophenetic state per window, and draw the WHOLE
time-ordered trajectory through a 3-D state space. The phase you are in EMERGES from the
motion alone -- nothing is coloured by phase, only by time:

    encoding axis  x = state . (task_learn - rest_pre)         (the leap the task makes)
    inference axis y = state . (task_test  - task_learn)_|x    (the drift within the task)
    residual axis  z = leading residual direction

    rest_pre   : the state OSCILLATES around baseline (a tight ball at the origin)
    task_learn : a JUMP out along encoding
    task_test  : a slighter DRIFT along inference
    rest_post  : TRACER stays near the task attractor ; RESETTER returns to baseline

The phase-average positions we characterised before (the basins) are over-laid as faint
ghosts, so the live trajectory carries the dynamics and the averages just annotate it. The
1:1:1 cube keeps the encoding leap large and the inference drift honestly small.

Because both rests share the recording's time-order, a slow drift would carry EVERY
rest_post outward; resetters returning to baseline is the built-in drift control. The
windowed per-patient separation is a noisy proxy for the verified trace (cohort corr with
the gate ~ -0.37), so the two trajectories ILLUSTRATE two clean, verified exemplars; the
cohort strip is the statistic of record -- matched-strength rho_sym inference-specific
trace per patient (audit_152/163), 6/10 clearing. Overlapping windows are autocorrelated;
no per-window significance is claimed.

Reads : data/audit/windowed_attractor/{pat}_beta_pca.npz              (audit_164)
        data/audit/consolidation_arc_rhosym/arc_null_per_patient.csv  (verified gate)
Writes: data/preprint/figures/results_section2/fig_arc_attractor_3d.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.colors import Normalize, to_rgba
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.ndimage import gaussian_filter1d
from scipy.stats import gaussian_kde
from contourpy import contour_generator

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

ATT = ROOT / "data/audit/windowed_attractor"
NULL = ROOT / "data/audit/consolidation_arc_rhosym/arc_null_per_patient.csv"
OUT = ROOT / "data/preprint/figures/results_section2/fig_arc_attractor_3d.pdf"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PH = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHLAB = {"rest_pre": r"rest$_\mathrm{pre}$", "task_learn": r"task$_\mathrm{learn}$",
         "task_test": r"task$_\mathrm{test}$", "rest_post": r"rest$_\mathrm{post}$"}
TASK = ("task_learn", "task_test")
HERO_TRACE, HERO_RESET = "Pat_06", "Pat_10"
GATE_A = 0.05
C_TRACE, C_RESET = "#2f9e44", "#adb5bd"
TCMAP = cm.plasma                                    # window-time colormap


def embed_ei(pat, band="beta"):
    """Encoding x, inference y, residual z; centred on the rest_pre baseline."""
    d = np.load(ATT / f"{pat}_{band}_pca.npz", allow_pickle=True)
    sc = d["scores"].astype(float)
    ph = d["phase"]
    tc = d["t_center"].astype(float)
    pre = sc[ph == "rest_pre"].mean(0)
    learn = sc[ph == "task_learn"].mean(0)
    test = sc[ph == "task_test"].mean(0)
    a1 = learn - pre
    a1 /= np.linalg.norm(a1) + 1e-9
    f = test - learn
    f = f - (f @ a1) * a1
    a2 = f / (np.linalg.norm(f) + 1e-9)
    S0 = sc - pre
    R = S0 - np.outer(S0 @ a1, a1) - np.outer(S0 @ a2, a2)
    Rc = R - R.mean(0)
    _, _, Vt = np.linalg.svd(Rc, full_matrices=False)
    return np.c_[S0 @ a1, S0 @ a2, Rc @ Vt[0]], ph, tc


def time_order(ph, tc):
    """Global window order: phases in sequence, by t_center within each phase."""
    idx = []
    for p in PH:
        m = np.where(ph == p)[0]
        idx.extend(m[np.argsort(tc[m])])
    return np.asarray(idx)


def _smooth_phase(Tp, sigma=2.2):
    """Gaussian-kernel smoothing within a phase: turns window jitter into a genuinely
    smooth flow while preserving the phase's mean position (jumps between phases stay
    sharp because smoothing is per-phase, edges held)."""
    if len(Tp) < 3:
        return Tp
    return np.stack([gaussian_filter1d(Tp[:, j], sigma, mode="nearest")
                     for j in range(3)], axis=1)


def _spread_surface(pts, n_theta=44, n_phi=24, level_q=0.24, tpad=1.8):
    """A star-shaped KDE isosurface around the cloud centroid: for each direction on
    the sphere, the radius at which the kernel density drops to a low quantile of its
    at-point values. Density-driven, so it undulates with the actual points (wiggly)
    and hugs their spread instead of imposing a smooth ellipsoid."""
    c = pts.mean(0)
    val, vec = np.linalg.eigh(np.cov((pts - c).T))
    val = np.clip(val, 1e-9, None)
    rmax = tpad * float(np.sqrt(val.max()))
    th = np.linspace(0, 2 * np.pi, n_theta)
    phg = np.linspace(0, np.pi, n_phi)
    TH, PHG = np.meshgrid(th, phg)
    ux, uy, uz = np.sin(PHG) * np.cos(TH), np.sin(PHG) * np.sin(TH), np.cos(PHG)
    U = np.stack([ux.ravel(), uy.ravel(), uz.ravel()], 1)
    try:
        kde = gaussian_kde(pts.T)
        level = float(np.quantile(kde(pts.T), level_q))
        ts = np.linspace(0.0, rmax, 28)
        grid = c[None, None, :] + ts[None, :, None] * U[:, None, :]
        D = kde(grid.reshape(-1, 3).T).reshape(U.shape[0], ts.size)
        cnt = np.cumprod(D >= level, axis=1).astype(bool).sum(1)
        r = np.where(cnt > 0, ts[np.clip(cnt - 1, 0, ts.size - 1)],
                     ts[np.argmax(D, axis=1)])
    except Exception:                                        # singular / tiny cloud
        minv = vec @ np.diag(1.0 / val) @ vec.T
        r = 1.6 / np.sqrt(np.einsum("md,dk,mk->m", U, minv, U))
    r = r.reshape(TH.shape)
    r = gaussian_filter1d(r, 1.0, axis=1, mode="wrap")       # smooth the wiggle
    r = gaussian_filter1d(r, 0.8, axis=0, mode="nearest")
    r[0, :], r[-1, :] = r[0, :].mean(), r[-1, :].mean()      # clean poles
    return c[0] + r * ux, c[1] + r * uy, c[2] + r * uz


def _ghost_basin(ax, pts, color, alpha=0.16, zorder=2):
    """Wiggly, lit, translucent spread surface (density isosurface; see _spread_surface)."""
    X, Y, Z = _spread_surface(pts)
    ax.plot_surface(X, Y, Z, color=color, alpha=alpha, linewidth=0.25,
                    edgecolor=(1, 1, 1, 0.07), antialiased=True, shade=True,
                    zorder=zorder)


def _draw_on_wall(ax, seg, plane, off, col, lw, al):
    """Embed one 2-D contour segment onto a fixed wall (x=, y=, or z=off)."""
    if plane == "z":
        ax.plot(seg[:, 0], seg[:, 1], off, color=col, lw=lw, alpha=al, zorder=0)
    elif plane == "y":
        ax.plot(seg[:, 0], np.full(len(seg), off), seg[:, 1], color=col, lw=lw,
                alpha=al, zorder=0)
    else:
        ax.plot(np.full(len(seg), off), seg[:, 0], seg[:, 1], color=col, lw=lw,
                alpha=al, zorder=0)


def _plane_contours(ax, XY, plane, off, col, lims2d, n_levels=5):
    """Concentric 2-D KDE level curves of a projection, drawn on one wall. The phase
    colour deepens (thicker, more opaque) toward higher density, so each shadow reads
    as a nested topographic density map rather than a single ellipse."""
    if XY.shape[0] < 5:
        return
    try:
        kde = gaussian_kde(XY.T)
    except Exception:
        return
    (a0, a1), (b0, b1) = lims2d
    ga, gb = np.linspace(a0, a1, 64), np.linspace(b0, b1, 64)
    GA, GB = np.meshgrid(ga, gb)
    D = kde(np.stack([GA.ravel(), GB.ravel()])).reshape(GA.shape)
    dmax = float(D.max())
    if dmax <= 0:
        return
    cg = contour_generator(ga, gb, D, line_type="Separate")
    for k, lv in enumerate(np.linspace(0.12 * dmax, 0.9 * dmax, n_levels)):
        frac = (k + 1) / n_levels                        # inner (denser) -> stronger
        lw, al = 0.5 + 1.2 * frac, 0.13 + 0.42 * frac
        for seg in cg.lines(lv):
            _draw_on_wall(ax, seg, plane, off, col, lw, al)


def _project_walls(ax, P, T, ph, pcol, lo, hi):
    """Shadow projections on the three walls: concentric KDE density level-curves per
    phase on each face + the trajectory's floor shadow -- position and spread become
    readable from all three sides as nested contour maps."""
    x0, y1, z0 = lo[0], hi[1], lo[2]
    xl, yl, zl = (lo[0], hi[0]), (lo[1], hi[1]), (lo[2], hi[2])
    ax.plot(T[:, 0], T[:, 1], z0, color="0.5", lw=1.3, alpha=0.22, zorder=0)
    for p in PH:
        col, pw = pcol[p], P[ph == p]
        _plane_contours(ax, pw[:, [0, 1]], "z", z0, col, (xl, yl))   # floor  (x,y)
        _plane_contours(ax, pw[:, [0, 2]], "y", y1, col, (xl, zl))   # back   (x,z)
        _plane_contours(ax, pw[:, [1, 2]], "x", x0, col, (yl, zl))   # side   (y,z)


def portrait(ax, pat, tag):
    """Load a patient's dense embedding and render its state-space portrait."""
    P, ph, tc = embed_ei(pat)
    render_portrait(ax, P, ph, tc, tag)


def render_portrait(ax, P, ph, tc, tag, subtitle=None):
    """Render one state-space portrait from a precomputed embedding ``(P, ph, tc)``.

    Split out from :func:`portrait` so callers with their own embedding (e.g. the
    mst@0.20 windowed-trajectory gallery) reuse the identical blob / wall-contour /
    glowing-flow rendering instead of forking it."""
    ax.computed_zorder = False                          # our explicit zorders decide layering
    # smoothed, time-ordered path: within-phase Gaussian smoothing turns window jitter
    # into a smooth flow; phase-to-phase jumps stay sharp (per-phase smoothing).
    Tsm = []
    for p in PH:
        m = np.where(ph == p)[0]
        m = m[np.argsort(tc[m])]
        Tsm.append(_smooth_phase(P[m]))
    T = np.concatenate(Tsm)
    tn = np.linspace(0.0, 1.0, len(T))

    # each phase cloud takes the colormap's central colour at that phase's median time
    counts = np.array([np.sum(ph == p) for p in PH], float)
    ends = np.cumsum(counts) / counts.sum()
    starts = np.concatenate([[0.0], ends[:-1]])
    pcol = {p: TCMAP(0.5 * (starts[i] + ends[i])) for i, p in enumerate(PH)}

    # tight 1:1:1 limits FIRST -- the walls (projections) live at these planes
    lo, hi = np.percentile(P, [2, 98], axis=0)
    span = hi - lo
    lo, hi = lo - 0.12 * span, hi + 0.12 * span
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])

    _project_walls(ax, P, T, ph, pcol, lo, hi)          # side projections, zorder 0
    for p in PH:                                         # wiggly spread surface + cloud
        m = ph == p
        _ghost_basin(ax, P[m], pcol[p], alpha=0.13, zorder=2)
        ax.scatter(P[m, 0], P[m, 1], P[m, 2], color=[pcol[p]], s=17, alpha=0.62,
                   edgecolor="white", linewidth=0.25, depthshade=False, zorder=5)

    pts = T.reshape(-1, 1, 3)                            # live flow, coloured by time
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    cols = TCMAP(0.5 * (tn[:-1] + tn[1:]))
    glow = Line3DCollection(segs, colors=cols, linewidths=8.5, alpha=0.13)
    glow.set_capstyle("round"); glow.set_zorder(6)
    core = Line3DCollection(segs, colors=cols, linewidths=3.1)
    core.set_capstyle("round"); core.set_zorder(7)
    ax.add_collection3d(glow); ax.add_collection3d(core)

    for p in PH:                                         # phase-average markers -- always on top
        c = P[ph == p].mean(0)
        ax.scatter(*c, s=175, facecolor="white", edgecolor=pcol[p], linewidth=2.9,
                   depthshade=False, zorder=20)

    ax.text2D(0.5, 0.99, tag, transform=ax.transAxes, ha="center", va="top",
              fontsize=14, fontweight="bold")
    if subtitle:
        ax.text2D(0.5, 0.945, subtitle, transform=ax.transAxes, ha="center", va="top",
                  fontsize=8.5, color="0.35")
    ax.set_xlabel("encoding", fontsize=13, labelpad=-2)
    ax.set_ylabel("inference", fontsize=13, labelpad=-2)
    ax.set_zlabel("residual", fontsize=11.5, labelpad=-8)
    ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
    ax.tick_params(length=0)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((0.984, 0.988, 0.996, 0.65))
        axis.line.set_color((0.62, 0.62, 0.66, 0.7))
        axis._axinfo["grid"].update(color=(0.80, 0.80, 0.85, 0.55), linewidth=0.6)
    ax.view_init(elev=18, azim=-58)
    ax.set_box_aspect((1, 1, 1), zoom=1.02)


def cohort_strip(ax, verd):
    v = verd.sort_values("obs").reset_index(drop=True)
    ax.axvline(0.0, color="0.7", lw=0.9, zorder=1)
    ys = np.arange(len(v))
    for y, (pat, obs, surr, p) in zip(ys, v[["patient", "obs", "surr", "p"]]
                                      .itertuples(index=False, name=None)):
        col = C_TRACE if p < GATE_A else C_RESET
        big = pat in (HERO_TRACE, HERO_RESET)
        ax.plot([surr, obs], [y, y], color=col, lw=1.4, alpha=0.45, zorder=2)
        ax.scatter(surr, y, s=42, color="0.55", marker="|", zorder=3)
        ax.scatter(obs, y, s=235 if big else 130, color=col, zorder=5,
                   edgecolor="#111" if big else "white", linewidth=2.0 if big else 1.0)
        ax.text(-0.075, y, pat.replace("Pat_", "P"), ha="right", va="center",
                fontsize=9 if big else 8.3, color="0.15" if big else "0.5",
                fontweight="bold" if big else "normal")
        if big:
            ax.text(obs, y + 0.42, "tracer" if pat == HERO_TRACE else "resetter",
                    ha="center", va="bottom", fontsize=8.5, color=col, style="italic")
    ax.set_xlim(-0.095, 0.42)
    ax.set_ylim(-0.7, len(v) - 0.1)
    ax.set_yticks([])
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4])
    ax.tick_params(labelsize=9)
    ax.set_xlabel(r"$\beta$ inference-specific trace  $\rho_\mathrm{sym}$ "
                  r"(matched-strength; $\vert$ = own null)", fontsize=10.5)
    ax.spines[["top", "right", "left"]].set_visible(False)
    return int((v.p < GATE_A).sum()), len(v)


def main():
    nb = pd.read_csv(NULL)
    b = nb[nb.band == "beta"]
    verd = pd.DataFrame({"patient": b.patient.values, "obs": b["T_infspec_pe_obs"].values,
                         "surr": b["T_infspec_pe_surr_p50"].values, "p": b["T_infspec_pe_p"].values})

    fig = plt.figure(figsize=(13.2, 9.6))
    gs = GridSpec(2, 2, height_ratios=[2.55, 1.0], hspace=0.0, wspace=0.02,
                  left=0.02, right=0.99, top=1.02, bottom=0.09)
    axL = fig.add_subplot(gs[0, 0], projection="3d")
    axR = fig.add_subplot(gs[0, 1], projection="3d")
    portrait(axL, HERO_TRACE, r"Pat_06 $\cdot$ tracer")
    portrait(axR, HERO_RESET, r"Pat_10 $\cdot$ resetter")
    axc = fig.add_subplot(gs[1, :])
    n_tr, n = cohort_strip(axc, verd)

    # time gradient bar (VECTOR pcolormesh) with the phase transitions marked on it
    cax = fig.add_axes([0.32, 0.345, 0.36, 0.015])
    g = np.linspace(0, 1, 256)[None, :]
    cax.pcolormesh(np.linspace(0, 1, 257), np.array([0, 1]), g, cmap=TCMAP,
                   rasterized=False, shading="flat")
    bl = []                                             # mean phase boundaries, two exemplars
    for pat in (HERO_TRACE, HERO_RESET):
        _, ph, _ = embed_ei(pat)
        cnt = np.array([np.sum(ph == p) for p in PH], float)
        bl.append(np.cumsum(cnt) / cnt.sum())
    cum = np.mean(bl, axis=0)
    mids = 0.5 * (np.concatenate([[0.0], cum[:-1]]) + cum)
    cax.vlines(cum[:-1], 0, 1, colors="white", lw=1.6, zorder=3)   # phase transitions
    cax.set_yticks([]); cax.set_xlim(0, 1)
    cax.set_xticks(mids)
    cax.set_xticklabels([PHLAB[p] for p in PH], fontsize=8.6)
    cax.tick_params(length=0)
    cax.set_title("window time  (one long recording)", fontsize=8.8, pad=3)
    for s in cax.spines.values():
        s.set_visible(False)

    handles = [Line2D([0], [0], marker="o", ls="", mfc=C_TRACE, mec="white", ms=10,
                      label=r"clears $\beta$ inference null"),
               Line2D([0], [0], marker="o", ls="", mfc=C_RESET, mec="white", ms=10,
                      label="within null"),
               Line2D([0], [0], marker="o", ls="", mfc="none", mec="0.25", ms=10,
                      label="phase-average (ghost basin)")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.01),
               ncol=3, frameon=False, fontsize=9.5, handletextpad=0.5, columnspacing=1.8)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True, pad_inches=0.04)
    plt.close(fig)

    print("fig:arc_attractor — live windowed state-space trajectory\n")
    print(f"cohort: {n_tr}/{n} clear the verified beta inference null (matched-strength)")
    for pat in (HERO_TRACE, HERO_RESET):
        P, ph, tc = embed_ei(pat)
        cen = {p: P[ph == p].mean(0) for p in PH}
        row = "  ".join(f"{p[:4]}({cen[p][0]:+.0f},{cen[p][1]:+.0f})" for p in PH)
        print(f"  {pat}: {row}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
