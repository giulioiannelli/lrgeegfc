#!/usr/bin/env python3
r"""talk_fig_descriptor_bloom -- the contrast measures as patient-coherency blooms (S6 talk).

The dark-slide, three-panel cut matching the Grassmann-contrast figure (Detect != discriminate):
raw FC edges (pairwise baseline) · Grassmann chordal (the field-standard spectral-subspace
measure) · cophenetic (multiscale, ours). Each panel shows the PER-PATIENT spread the contrast
dot-matrix cannot.

CAVEAT (honest): raw + cophenetic are the mst@0.20 backbone trace at scale s=5.6 (rho_sym units);
Grassmann is a WHOLE-GRAPH signed subspace-distance trace T_G, aggregated over the k-grid -- a
different construction in different units, so the Grassmann panel carries its OWN radial scale and
its wedge glow is the matched-strength CLUSTER-EXTENT-over-k verdict (beta + low_gamma clear,
delta LOO-fragile, theta/alpha/high_gamma null), NOT a per-k mean (high_gamma has a high mean but
no coherent significant cluster). Petal length is each patient's mean-over-k T_G.

Each panel: six band sectors around a task->rest hub; ten patients as petals grown from a T=0
waterline (outward = trace kept, inward = reverted). A wedge glows only when the band clears its
cohort gate. Reading, left to right:
  raw FC edges  -> beta blooms (a CONVERGENCE band); alpha only flickers (fragile).
  Grassmann     -> beta + low_gamma clear (the spectral read agrees on beta, adds low_gamma,
                   MISSES alpha) -- contrastive to ours.
  cophenetic    -> the glow COLLAPSES onto alpha+beta; alpha is cophenetic-EXCLUSIVE.

Dark-slide mode: ink lightened for a dark background, transparent PNG, colored data untouched.
Reads : data/sparsified_arc/controls_ladder_apples/per_cell.csv (raw, cophenetic; s_meso=5.6),
        data/audit/epi_stratified/grassmann_per_patient_per_k.csv (config=full; per-patient T_G),
        data/audit/grassmann_cluster_extent/cohort_summary.csv     (cluster-extent verdict + LOO).
Writes: data/outputs/figures/talk/fig_descriptor_bloom_contrast.{png,pdf}
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path as MPath

sys.path.insert(0, str(Path(__file__).resolve().parents[1] /
                       "01_compute" / "figures_embedded" / "new_results_sec1"))
import _common as C  # noqa: E402
from lrg_eegfc.config.const import BRAIN_BANDS_NAMES  # noqa: E402

C.use_lrg_style()

# dark-slide mode: all ink white; hardcoded greys lightened; colored data untouched
INK = "white"
plt.rcParams.update({
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
    "axes.titlecolor": INK, "xtick.color": INK, "ytick.color": INK,
})

# dark-friendly palette (originals were dark-on-white; flipped for a dark slide)
SPOKE   = "#7f858d"      # radial band separators
TICKRING = "#5f656d"     # dotted T-tick rings
NULLBAND = "#aeb4bc"     # matched-strength null annulus (soft light band)
WATER    = "#cfd4db"     # T=0 waterline
TICKLAB  = "#c2c7cf"     # T-axis numerals (no white pill on dark)
RINGQ    = "#e8eaee"     # rho-symbol axis label
RIM_REP, RIM_NON   = "white", "#aab0b8"     # band rim label
PVAL_REP, PVAL_NON = "white", "#9aa0a7"     # cohort-p under the rim label
NONREP   = "#9297a0"     # non-representative petals + median arc
FOOT     = "#9298a1"     # footer caption

# the three contrast read-outs: raw (baseline) · Grassmann (field standard) · cophenetic (ours)
PANELS = [
    ("raw_fc",     "raw FC edges\n(pairwise · baseline)",  r"$\rho^{\mathrm{raw}}$"),
    ("grassmann",  "Grassmann chordal\n(spectral · whole graph)", r"$T_{\mathrm{G}}$"),
    ("coph_meso",  r"cophenetic  ($s=5.6$)" + "\n(multiscale · ours)",
                   r"$\rho^{\mathrm{coph}}_{\mathrm{sym}}$"),
]
FUNC = "T_test"
R_HUB, R_OUT = 0.30, 1.00
C_ANTI = np.array([0.55, 0.58, 0.63])            # cool grey for reverted petals (lightened)
BETA_I = BRAIN_BANDS_NAMES.index("beta")
OUT = C.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_descriptor_bloom_contrast"


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


def load_grassmann():
    """per-(patient,band) whole-graph Grassmann trace T_G (mean over the k-grid) + the
    matched-strength CLUSTER-EXTENT cohort verdict. T_G is a signed subspace-distance trace
    (not a rho), so this panel gets its OWN radial scale; the wedge glow is the cluster-extent
    -over-k verdict (beta + low_gamma clear, delta LOO-fragile, else null), NOT a per-k mean."""
    per = pd.read_csv(C.ROOT / "data" / "audit" / "epi_stratified" /
                      "grassmann_per_patient_per_k.csv")
    per = per[(per.config == "full") & (per.substrate == "grassmann")]
    rows = []
    for (pat, band), g in per.groupby(["patient", "band"]):
        frac_above = float((g.obs_stat.values > g.surr_p50.values).mean())
        rows.append(dict(patient=pat, band=band, obs=float(g.obs_stat.mean()),
                         surr=float(g.surr_p50.mean()),
                         pp=float(np.clip(1.0 - frac_above, 0.0, 1.0))))
    df = pd.DataFrame(rows)
    coh = pd.read_csv(C.ROOT / "data" / "audit" / "grassmann_cluster_extent" /
                      "cohort_summary.csv")
    verdict = {}
    for band in BRAIN_BANDS_NAMES:
        r = coh[coh.band == band].iloc[0]
        strong = str(r.verdict_cluster_extent) == "strong"
        loo_ok = float(r.cluster_p_mass_loo_max) < 0.05
        lvl = "rep" if (strong and loo_ok) else ("fragile" if strong else "null")
        sub = df[df.band == band]
        verdict[band] = dict(level=lvl, cohort_p=float(r.cluster_p_cluster_mass),
                             n_above=int((sub.obs.values > sub.surr.values).sum()),
                             loo_worst=float(r.cluster_p_mass_loo_max))
    return df, verdict


def _glow(res):
    """Wedge-glow strength keyed to the representative gate (not raw cohort p)."""
    p = res.get("cohort_p", 1.0)
    if res["level"] == "rep":
        return 0.16 + 0.56 * sigfrac(p)          # slightly hotter for a dark slide
    if res["level"] == "fragile":
        return 0.10
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
                color=SPOKE, lw=0.7, alpha=0.6, zorder=0.45)

    ticks = [t for t in np.round(np.arange(-0.4, 0.61, 0.2), 1)
             if lo + 0.03 < t < hi - 0.005 and abs(t) > 1e-9]
    for Tg in ticks:
        ring(ax, r_of(Tg, lo, hi), color=TICKRING, lw=0.5, ls=(0, (1, 6)), zorder=0.3)
    nz = df["surr"].to_numpy(float)
    annulus(ax, r_of(np.percentile(nz, 10), lo, hi), r_of(np.percentile(nz, 90), lo, hi),
            color=NULLBAND, alpha=0.20, lw=0, zorder=0.35)
    ring(ax, r0, color=WATER, lw=1.05, zorder=0.5)                    # waterline T=0
    for Tg in [0.0] + ticks:
        lab = "0" if Tg == 0.0 else f"{Tg:+.1f}"
        rr = r_of(Tg, lo, hi)
        ax.text(-rr, 0.0, lab, fontsize=9.5, color=TICKLAB, ha="center", va="center",
                zorder=6, rotation=90)
    ax.text(-(R_OUT + 0.42), 0.0, ring_q, fontsize=19, color=RINGQ,
            ha="center", va="center", rotation=90, zorder=6)

    # rim band labels + cohort p: CLEARS (rep bold / fragile) -> band colour; null -> white
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        theta0 = np.deg2rad(90 + (BETA_I - i) * 60)
        res = verdict[band]
        pv = res.get("cohort_p", np.nan)
        rep = res["level"] == "rep"
        clears = res["level"] in ("rep", "fragile")     # band colour; only null -> white
        col = C.band_color(band)
        ax.text((R_OUT + 0.17) * np.cos(theta0), (R_OUT + 0.17) * np.sin(theta0),
                C.BTeX[band], fontsize=27 if rep else 21, ha="center", va="center",
                zorder=6, color=col if clears else "white",
                fontweight="bold" if rep else "normal")
        ptxt = (rf"$p=\mathbf{{{pv:.3f}}}$" if clears else f"$p={pv:.2f}$")
        ax.text((R_OUT + 0.39) * np.cos(theta0), (R_OUT + 0.39) * np.sin(theta0),
                ptxt, fontsize=14 if rep else 12, ha="center", va="center", zorder=6,
                color="white", fontweight="bold" if rep else "normal")

    ax.add_patch(plt.Circle((0, 0), R_HUB * 0.74, fc="none", ec=WATER, lw=1.1, zorder=3))
    ax.text(0, 0.045, "task", fontsize=12, ha="center", va="center", color="white",
            zorder=4, fontweight="bold")
    ax.text(0, -0.070, r"$\rightarrow$ rest", fontsize=10, ha="center", va="center",
            color=TICKLAB, zorder=4)

    # petals (one per patient)
    petal_span = np.deg2rad(60) * 0.82
    offs = np.linspace(-petal_span / 2, petal_span / 2, len(C.COHORT))
    half_w = (petal_span / len(C.COHORT)) * 0.46
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        theta0 = np.deg2rad(90 + (BETA_I - i) * 60)
        hue = np.asarray(to_rgb(C.band_color(band)))
        res = verdict[band]
        rep = res["level"] == "rep"
        clears = res["level"] in ("rep", "fragile")
        pv = res.get("cohort_p", 1.0)
        em = (0.42 + 0.58 * sigfrac(pv)) if rep else (0.50 if clears else 0.40)
        sub = df[df.band == band].sort_values("obs")
        obs = sub["obs"].to_numpy(float)
        pp = sub["pp"].to_numpy(float)
        med = float(np.median(obs)) if len(obs) else 0.0
        aa = np.linspace(theta0 - petal_span / 2, theta0 + petal_span / 2, 60)
        rr = r_of(med, lo, hi)
        ax.plot(rr * np.cos(aa), rr * np.sin(aa), color=hue if clears else "white",
                lw=3.4 if rep else (2.3 if clears else 1.5),
                alpha=0.98 if clears else 0.55, zorder=1.7, solid_capstyle="round")
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

    lim = R_OUT + 0.60
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)


def main():
    pcfull = pd.read_csv(C.LADDER / "per_cell.csv")
    data = {d: (load_grassmann() if d == "grassmann" else load_descriptor(pcfull, d))
            for d, _, _ in PANELS}

    def _rng(obs):
        hi = float(max(np.percentile(obs, 97), 0.30)) * 1.08
        lo = float(min(np.percentile(obs, 3), -0.05)) * 1.06
        return lo, hi
    # raw + cophenetic share a rho-unit scale; Grassmann (whole-graph T_G) gets its own
    rho_obs = np.concatenate([data["raw_fc"][0]["obs"].to_numpy(float),
                              data["coph_meso"][0]["obs"].to_numpy(float)])
    ranges = {"raw_fc": _rng(rho_obs), "coph_meso": _rng(rho_obs),
              "grassmann": _rng(data["grassmann"][0]["obs"].to_numpy(float))}

    fig, axes = plt.subplots(1, 3, figsize=(18.6, 7.2))
    axes = np.atleast_1d(axes).ravel()
    for ax, (descriptor, title, ring_q) in zip(axes, PANELS):
        df, verdict = data[descriptor]
        lo, hi = ranges[descriptor]
        draw_bloom(ax, df, verdict, ring_q, lo, hi)
        ax.set_title(title, fontsize=20, fontweight="bold", pad=8, color="white")

    fig.text(0.5, 0.045, "petals = patients   ·   length = per-patient trace, inward = "
             "reverted   ·   wedge glow = cohort significance   ·   raw & cophenetic: "
             r"mst@0.20 at $s=5.6$   ·   Grassmann: whole graph, mean over $k$",
             ha="center", va="center", fontsize=12, color=FOOT)

    fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.10, wspace=0.14)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}.{ext}", transparent=True, dpi=200)
    plt.close(fig)

    print("talk_fig_descriptor_bloom -- contrast subset blooms (T_test, mst@0.20 s=5.6):")
    for descriptor, title, _ in PANELS:
        _, verdict = data[descriptor]
        rep = [b for b in BRAIN_BANDS_NAMES if verdict[b]["level"] == "rep"]
        frag = [b for b in BRAIN_BANDS_NAMES if verdict[b]["level"] == "fragile"]
        lab = title.split("\n")[0]
        print(f"  {lab:20s} representative: {', '.join(rep) or '(none)':16s}"
              f" fragile: {', '.join(frag) or '(none)'}")
    print(f"wrote {OUT}.png + .pdf")


if __name__ == "__main__":
    main()
