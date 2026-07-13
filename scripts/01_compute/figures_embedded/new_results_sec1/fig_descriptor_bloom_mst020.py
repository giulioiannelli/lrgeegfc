#!/usr/bin/env python3
r"""fig_descriptor_bloom_mst020 -- the whole descriptor ladder as patient-coherency blooms.

HEAD: one radial bloom per graph descriptor, all on the **mst@0.20** backbone at the
reporting scale **s=5.6** under the SINGLE matched-strength null (drift is discarded --
task is directional, so the drift null collapses onto the alternative). Six band sectors
around a hub; each of the ten patients is a petal grown from a T=0 waterline (outward =
trace kept, inward = reverted). A band's wedge glows only when it passes the
REPRESENTATIVE gate (cohort p<0.05 AND leave-one-out AND >=3/10 patients above their own
null; _common.representativeness_gate) -- so a band that merely lights the naive cohort
test but is single-patient-driven does NOT bloom. Scan the six panels against the
cophenetic one (bottom-right):

  raw FC edges  ->  beta blooms (a CONVERGENCE band -- raw backbone edges carry it too),
                    alpha only flickers (fragile ring, not representative).
  strength / clustering / resistance  ->  no representative bloom in any band.
  graph geodesic  ->  a single bloom in low_gamma, a NON-carrier band.
  cophenetic    ->  the glow COLLAPSES onto alpha+beta and nothing else. Same patients,
                    same null, same scale -- only the read-out changes. alpha is
                    cophenetic-EXCLUSIVE; the hierarchy is the one read-out holding BOTH
                    carriers with no false positives.

Descriptor pipeline: controls_ladder_apples (ALL descriptors computed on the mst@0.20
backbone -- apples-to-apples). Petal length is the per-patient T_test; the cohort verdict
is the wedge glow (representative vs matched-strength). The contribution is the
SELECTIVITY + scale-shape of the cophenetic read-out, NOT that it detects a trace the
others miss for beta (they don't) -- it is exclusive for alpha and clean for both.

Reads : data/sparsified_arc/controls_ladder_apples/per_cell.csv (frac=0.2, s_meso=5.6,
        matched-strength R=200; functional T_test).
Writes: data/preprint/figures/new_results_sec1/fig_descriptor_bloom_mst020.pdf
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path as MPath

import _common as C
from lrg_eegfc.config.const import BRAIN_BANDS_NAMES

C.use_lrg_style()

# the full ladder as blooms (ladder order; cophenetic last = the payoff panel)
PANELS = [
    ("raw_fc",     "raw FC edges",           r"$\rho^{\mathrm{raw}}$"),
    ("strength",   "node strength",          r"$\rho^{\mathrm{str}}$"),
    ("clustering", "clustering coef.",       r"$\rho^{\mathrm{clust}}$"),
    ("geodesic",   "graph geodesic",         r"$\rho^{\mathrm{geo}}$"),
    ("resistance", "resistance (spectral)",  r"$\rho^{\mathrm{res}}$"),
    ("coph_meso",  r"cophenetic  ($s=5.6$)", r"$\rho^{\mathrm{coph}}_{\mathrm{sym}}$"),
]
FUNC = "T_test"
R_HUB, R_OUT = 0.30, 1.00
C_ANTI = np.array([0.49, 0.52, 0.57])            # cool grey for reverted petals
BETA_I = BRAIN_BANDS_NAMES.index("beta")


# ---- pure geometry (radial bloom primitives) ---------------------------------
def r_of(T, lo, hi):
    return R_HUB + (R_OUT - R_HUB) * (np.clip(T, lo, hi) - lo) / (hi - lo)


def sigfrac(p):
    return float(np.clip((-np.log10(max(p, 1e-6)) - 0.5) / 1.5, 0.0, 1.0))


def petal_outline(theta_c, half_w, r0, r1, bulge=0.55):
    verts = [
        (theta_c - half_w, r0),
        (theta_c - half_w * bulge, r0 + (r1 - r0) * 0.45),
        (theta_c - half_w * 0.12, r0 + (r1 - r0) * 0.86),
        (theta_c, r1),
        (theta_c + half_w * 0.12, r0 + (r1 - r0) * 0.86),
        (theta_c + half_w * bulge, r0 + (r1 - r0) * 0.45),
        (theta_c + half_w, r0),
    ]
    xy = [(r * np.cos(t), r * np.sin(t)) for t, r in verts]
    xy.append(xy[0])
    codes = [MPath.MOVETO] + [MPath.CURVE4] * 6 + [MPath.LINETO]
    return MPath(xy, codes)


def wedge_gradient(ax, cdeg, color, amax, gamma=1.5, N=40):
    edges = np.linspace(R_HUB, R_OUT, N + 1)
    for a in range(N):
        t = (a + 0.5) / N
        ax.add_patch(Wedge((0, 0), edges[a + 1], cdeg - 30, cdeg + 30,
                           width=edges[a + 1] - edges[a], facecolor=color,
                           alpha=amax * t ** gamma, ec="none", zorder=0.2))


def ring(ax, r, **kw):
    a = np.linspace(0, 2 * np.pi, 400)
    ax.plot(r * np.cos(a), r * np.sin(a), **kw)


def annulus(ax, r0, r1, **kw):
    a = np.linspace(0, 2 * np.pi, 400)
    xo, yo = r1 * np.cos(a), r1 * np.sin(a)
    xi, yi = r0 * np.cos(a[::-1]), r0 * np.sin(a[::-1])
    ax.fill(np.concatenate([xo, xi]), np.concatenate([yo, yi]), **kw)


def clip_cap(ax, th, color, outward=True):
    r = (R_OUT + 0.02) if outward else (R_HUB - 0.02)
    s = 0.022 if outward else -0.022
    ang = np.array([th - 0.013, th, th + 0.013])
    rad = np.array([r, r + s, r])
    ax.plot(rad * np.cos(ang), rad * np.sin(ang), color=color, lw=1.1,
            zorder=5, solid_capstyle="round")


# ---- data --------------------------------------------------------------------
def load_descriptor(pcfull, descriptor):
    """per-(patient,band) obs/null/p + per-band representative verdict for one descriptor."""
    pc = pcfull[pcfull.descriptor == descriptor]
    df = pc[["patient", "band"]].copy()
    df["obs"] = pc[f"{FUNC}_obs"].astype(float)
    df["surr"] = pc[f"{FUNC}_surr_p50"].astype(float)
    df["pp"] = pc[f"{FUNC}_p"].astype(float)
    verdict = {b: C.representativeness_gate(pcfull, descriptor, b, FUNC)
               for b in BRAIN_BANDS_NAMES}
    return df, verdict


def _glow(res):
    """Wedge-glow strength keyed to the representative gate (not raw cohort p)."""
    p = res.get("cohort_p", 1.0)
    if res["level"] == "rep":
        return 0.14 + 0.52 * sigfrac(p)
    if res["level"] == "fragile":
        return 0.09
    return 0.05


# ---- one bloom panel ---------------------------------------------------------
def draw_bloom(ax, df, verdict, ring_q, lo, hi):
    ax.set_aspect("equal")
    ax.axis("off")
    r0 = r_of(0.0, lo, hi)

    for i, band in enumerate(BRAIN_BANDS_NAMES):
        c = 90 + (BETA_I - i) * 60
        wedge_gradient(ax, c, C.band_color(band), _glow(verdict[band]))

    for bnd in (0, 60, 120, 180, 240, 300):
        a = np.deg2rad(bnd)
        ax.plot([R_HUB * np.cos(a), (R_OUT + 0.015) * np.cos(a)],
                [R_HUB * np.sin(a), (R_OUT + 0.015) * np.sin(a)],
                color="#b9bdc4", lw=0.7, alpha=0.6, zorder=0.45)

    ticks = [t for t in np.round(np.arange(-0.4, 0.61, 0.2), 1)
             if lo + 0.03 < t < hi - 0.005 and abs(t) > 1e-9]
    for Tg in ticks:
        ring(ax, r_of(Tg, lo, hi), color="#d0d3d9", lw=0.5, ls=(0, (1, 6)), zorder=0.3)
    nz = df["surr"].to_numpy(float)
    annulus(ax, r_of(np.percentile(nz, 10), lo, hi), r_of(np.percentile(nz, 90), lo, hi),
            color="#8a9096", alpha=0.22, lw=0, zorder=0.35)
    ring(ax, r0, color="#5c6166", lw=1.05, zorder=0.5)                 # waterline T=0
    for Tg in [0.0] + ticks:
        lab = "0" if Tg == 0.0 else f"{Tg:+.1f}"
        rr = r_of(Tg, lo, hi)
        ax.text(-rr, 0.0, lab, fontsize=6.5, color="#8f959d", ha="center", va="center",
                zorder=6, rotation=90,
                bbox=dict(boxstyle="round,pad=0.05", fc="white", ec="none", alpha=0.75))
    ax.text(-(R_OUT + 0.17), 0.0, ring_q, fontsize=11, color="#3f454c",
            ha="center", va="center", rotation=90, zorder=6)

    # rim band labels + cohort p (bold only where REPRESENTATIVE)
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        theta0 = np.deg2rad(90 + (BETA_I - i) * 60)
        res = verdict[band]
        pv = res.get("cohort_p", np.nan)
        rep = res["level"] == "rep"
        ax.text((R_OUT + 0.13) * np.cos(theta0), (R_OUT + 0.13) * np.sin(theta0),
                C.BTeX[band], fontsize=17 if rep else 13, ha="center", va="center",
                zorder=6, color="#141414" if rep else "#565b62",
                fontweight="bold" if rep else "normal")
        ptxt = (rf"$p=\mathbf{{{pv:.3f}}}$" if rep else f"$p={pv:.2f}$")
        ax.text((R_OUT + 0.25) * np.cos(theta0), (R_OUT + 0.25) * np.sin(theta0),
                ptxt, fontsize=9 if rep else 7.5, ha="center", va="center", zorder=6,
                color="#111111" if rep else "#9aa0a7",
                fontweight="bold" if rep else "normal")

    ax.add_patch(plt.Circle((0, 0), R_HUB * 0.72, fc="#f6f7f9", ec="#d5d8dd",
                            lw=0.7, zorder=3))
    ax.text(0, 0.030, "task", fontsize=8.0, ha="center", va="center", color="#3f454c",
            zorder=4, fontweight="bold")
    ax.text(0, -0.050, r"$\rightarrow$ rest", fontsize=6.8, ha="center", va="center",
            color="#5a6068", zorder=4)

    # petals (one per patient)
    petal_span = np.deg2rad(60) * 0.82
    offs = np.linspace(-petal_span / 2, petal_span / 2, len(C.COHORT))
    half_w = (petal_span / len(C.COHORT)) * 0.46
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        theta0 = np.deg2rad(90 + (BETA_I - i) * 60)
        hue = np.asarray(to_rgb(C.band_color(band)))
        res = verdict[band]
        rep = res["level"] == "rep"
        pv = res.get("cohort_p", 1.0)
        em = (0.42 + 0.58 * sigfrac(pv)) if rep else 0.40
        sub = df[df.band == band].sort_values("obs")
        obs = sub["obs"].to_numpy(float)
        pp = sub["pp"].to_numpy(float)
        med = float(np.median(obs)) if len(obs) else 0.0
        aa = np.linspace(theta0 - petal_span / 2, theta0 + petal_span / 2, 60)
        rr = r_of(med, lo, hi)
        ax.plot(rr * np.cos(aa), rr * np.sin(aa), color=hue if rep else "#8f949b",
                lw=3.2 if rep else 1.6, alpha=0.98 if rep else 0.5, zorder=1.7,
                solid_capstyle="round")
        for k in range(len(obs)):
            pos = obs[k] > 0
            g = float(np.clip((0.5 - pp[k]) / 0.5, 0.0, 1.0))
            base = hue if pos else C_ANTI
            th = theta0 + offs[k]
            r1 = r_of(obs[k], lo, hi)
            if pos and g > 0.4:
                pg = petal_outline(th, half_w * 1.85, r0, r1)
                ax.add_patch(PathPatch(pg, fc=base, ec="none",
                                       alpha=(0.05 + 0.22 * g) * em, zorder=2.0))
            p = petal_outline(th, half_w, r0, r1)
            fa = ((0.34 + 0.5 * g) if pos else (0.28 + 0.12 * g)) * em
            ax.add_patch(PathPatch(p, fc=(*base, fa), ec=(*base, 0.85 * em),
                                   lw=0.5 + 0.5 * g, zorder=2.5 + 0.5 * g,
                                   joinstyle="round"))
            if obs[k] > hi or obs[k] < lo:
                clip_cap(ax, th, base, outward=obs[k] > 0)

    lim = R_OUT + 0.36
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)


def main():
    pcfull = pd.read_csv(C.LADDER / "per_cell.csv")
    data = {d: load_descriptor(pcfull, d) for d, _, _ in PANELS}
    allobs = np.concatenate([data[d][0]["obs"].to_numpy(float) for d, _, _ in PANELS])
    hi = float(max(np.percentile(allobs, 97), 0.30)) * 1.08
    lo = float(min(np.percentile(allobs, 3), -0.05)) * 1.06

    fig, axes = plt.subplots(2, 3, figsize=(18.5, 11.6))
    axes = axes.ravel()
    letters = "abcdef"
    for ax, (descriptor, title, ring_q), letter in zip(axes, PANELS, letters):
        df, verdict = data[descriptor]
        draw_bloom(ax, df, verdict, ring_q, lo, hi)
        ax.set_title(title, fontsize=14, fontweight="bold", pad=4)
        ax.text(0.01, 0.99, rf"$\mathbf{{{letter}}}$", transform=ax.transAxes,
                fontsize=16, va="top", ha="left", fontweight="bold")

    fig.text(0.5, 0.055, "petals = patients   ·   length = per-patient trace, inward = "
             "reverted   ·   wedge glow = REPRESENTATIVE cohort significance "
             r"(gate $+$ LOO $+$ consistency)   ·   mst@0.20, $s=5.6$",
             ha="center", va="center", fontsize=9.5, color="#8f959d")

    fig.subplots_adjust(left=0.03, right=0.98, top=0.94, bottom=0.09, hspace=0.20, wspace=0.06)
    C.FIGDIR.mkdir(parents=True, exist_ok=True)
    out = C.FIGDIR / "fig_descriptor_bloom_mst020.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)

    print("fig_descriptor_bloom_mst020 -- full ladder blooms (T_test, mst@0.20 s=5.6):")
    for descriptor, title, _ in PANELS:
        _, verdict = data[descriptor]
        rep = [b for b in BRAIN_BANDS_NAMES if verdict[b]["level"] == "rep"]
        frag = [b for b in BRAIN_BANDS_NAMES if verdict[b]["level"] == "fragile"]
        print(f"  {title:22s} representative: {', '.join(rep) or '(none)':22s}"
              f" fragile: {', '.join(frag) or '(none)'}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
