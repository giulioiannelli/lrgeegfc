#!/usr/bin/env python3
r"""fig:N2 trace blooms — task-trace scalars as a radial bloom, generalized.

Six band sectors around a hub. Each petal = one patient. Petal length = that
patient's task-trace correlation for the chosen functional, grown from a T=0
waterline (outward = trace kept, inward = reverted). A band's wedge glows in
proportion to its cohort matched-strength significance (-log10 Wilcoxon p), so the
band(s) that actually clear the null bloom out of the ring; the rest hug it.

Functionals (Spearman rank traces across region-pairs; positive = TRACE, audit_152):
    e = D_learn - D_pre    g = D_test - D_pre    f = D_test - D_learn    p = D_post - D_pre
    T_test       = rho_s(g, p)             overall reinstatement
    T_learn      = rho_s(e, p)             encoding retained (alpha/beta, broadband)
    T_infspec    = rho_s(f, p)             inference retained (uncontrolled)
    T_infspec_pe = partial rho_s(f, p | e) inference retained, encoding out -> BETA-ONLY
The radial value is this trace correlation T (built FROM rho_sym cophenetic
distances); it is NOT rho_sym itself. Cohort p read from the locked verdict file.

One extra, non-arc functional -- the MAIN cross-phase gate itself:
    rho_sym_gate = per-patient rho_sym (task_test vs rest_post cophenetic trace),
                   the R-2 "which bands hold" measure. Here the radial value IS
                   rho_sym, and the cohort p is gate_p_sym (Wilcoxon vs
                   matched-strength). Reads data/audit/rho_sym_gate/ instead of
                   the consolidation-arc files. beta / alpha bloom (p .032 / .024);
                   gamma_low / gamma_high / delta hug the ring (p .080); theta absent.

Two builders, each runs over the requested functionals:
  build_bloom     -> headline (petal colour = band)
  build_coverage  -> diagnostic (petal colour = dominant gray-matter territory;
                     EXPLORATORY, no subsetting, no new null)

Reads : data/audit/consolidation_arc_rhosym/arc_null_per_patient.csv (per-patient obs/surr/p)
        data/audit/consolidation_arc_rhosym/arc_cohort_verdict.csv    (cohort Wilcoxon p)
        data/audit/rho_sym_gate/per_patient_per_band.csv              (gate: per-patient rho_sym)
        data/audit/rho_sym_gate/cohort_summary.csv                    (gate: cohort gate_p_sym)
        data/audit/per_node_trace_decomposition_rhosym/per_node.csv   (coverage)
Writes: data/preprint/figures/_drafts/fig_reasoning_bloom{,_encoding,_reinstatement,_inference_raw,_gate}{,_coverage}_DRAFT.pdf
"""
from __future__ import annotations

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path as MPath
from matplotlib.lines import Line2D

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()

SRC = ROOT / "data/audit/consolidation_arc_rhosym/arc_null_per_patient.csv"
VERDICT = ROOT / "data/audit/consolidation_arc_rhosym/arc_cohort_verdict.csv"
GATE_SRC = ROOT / "data/audit/rho_sym_gate/per_patient_per_band.csv"
GATE_VERDICT = ROOT / "data/audit/rho_sym_gate/cohort_summary.csv"
PN_COV = ROOT / "data/audit/per_node_trace_decomposition_rhosym/per_node.csv"
OUTDIR = ROOT / "data/preprint/figures/_drafts"

# dual-null (matched-strength AND drift) sources ---------------------------------
DRIFT_LADDER = ROOT / "data/audit/drift_null_ladder/per_patient_per_band.csv"
PAIR_LADDER = ROOT / "data/audit/pairwise_descriptor_ladder/per_patient_per_band.csv"
PAIR_COHORT = ROOT / "data/audit/pairwise_descriptor_ladder/cohort_summary.csv"
GARC = ROOT / "data/audit/grassmann_inference_arc/per_patient.csv"
GMS_PERK = ROOT / "data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv"
GCLUSTER = ROOT / "data/audit/grassmann_cluster_extent/cohort_summary.csv"
GDRIFT = ROOT / "data/audit/grassmann_drift_null/cohort_summary.csv"
DRIFT_COL = "#b0322d"  # drift waterline (must be overcome IN ADDITION to MS)

# functional specs: which trace, how to label it, output stem -------------------
FUNCS = {
    "rho_sym_gate": dict(
        hub="task", ring_q=r"$\rho_{\mathrm{sym}}$",
        axis="task structure kept\n(task -> rest)",
        stem="fig_reasoning_bloom_gate"),
    "T_infspec_pe": dict(
        hub="inference", ring_q=r"$\rho_s(f,\,p\,|\,e)$",
        axis="inference-specific\nstructure kept",
        stem="fig_reasoning_bloom"),
    "T_learn": dict(
        hub="encoding", ring_q=r"$\rho_s(e,\,p)$",
        axis="encoding\nstructure kept",
        stem="fig_reasoning_bloom_encoding"),
    "T_test": dict(
        hub="task", ring_q=r"$\rho_s(g,\,p)$",
        axis="overall reinstatement\n(task -> rest)",
        stem="fig_reasoning_bloom_reinstatement"),
    "T_infspec": dict(
        hub="inference", ring_q=r"$\rho_s(f,\,p)$",
        axis="inference structure kept\n(encoding NOT removed)",
        stem="fig_reasoning_bloom_inference_raw"),
}

R_HUB, R_OUT = 0.30, 1.00
C_ANTI = np.array([0.49, 0.52, 0.57])       # cool grey for reverted petals
T_LO, T_HI = -0.42, 0.46                     # set per functional in build()

# coverage (diagnostic) ---------------------------------------------------------
SUPER2MACRO = {"PFC": "frontal", "limbic": "limbic", "lateral_temporal": "temporal",
               "parietal": "posterior", "sensorimotor": "posterior",
               "occipital": "posterior"}
MACRO_ORDER = ["frontal", "limbic", "temporal", "posterior"]
MACRO_COL = {"frontal": "#d1495b", "limbic": "#8a5fbf", "temporal": "#2a9d8f",
             "posterior": "#3d6fb4"}
MACRO_LAB = {"frontal": "frontal (PFC)", "limbic": "limbic (OFC/cing/MTL)",
             "temporal": "temporal / insula", "posterior": "parietal / SM / occ."}


# ---- primitives ---------------------------------------------------------------
def r_of(T):
    t = np.clip(T, T_LO, T_HI)
    return R_HUB + (R_OUT - R_HUB) * (t - T_LO) / (T_HI - T_LO)


def band_hue(i):
    return np.array(plt.get_cmap("turbo")(0.06 + 0.88 * i / 5.0))[:3]


def sigfrac(p):
    return float(np.clip((-np.log10(max(p, 1e-6)) - 0.5) / 1.5, 0.0, 1.0))


def petal_outline(theta_c, half_w, r0, r1, bulge=0.55):
    """slim triangular petal: base at r0 (+-half_w), tip at r1 (r1<r0 => inward)."""
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
    """radial-gradient annular wedge: ~transparent at hub, amax at rim."""
    edges = np.linspace(R_HUB, R_OUT, N + 1)
    for a in range(N):
        t = (a + 0.5) / N
        ax.add_patch(Wedge((0, 0), edges[a + 1], cdeg - 30, cdeg + 30,
                           width=edges[a + 1] - edges[a], facecolor=color,
                           alpha=amax * t ** gamma, ec="none", zorder=0.2))


def clip_cap(ax, th, color, outward=True):
    """small chevron just past the rim (or hub) marking an off-scale petal."""
    r = (R_OUT + 0.02) if outward else (R_HUB - 0.02)
    s = 0.022 if outward else -0.022
    ang = np.array([th - 0.013, th, th + 0.013])
    rad = np.array([r, r + s, r])
    ax.plot(rad * np.cos(ang), rad * np.sin(ang), color=color, lw=1.1,
            zorder=5, solid_capstyle="round")


def ring(ax, r, **kw):
    a = np.linspace(0, 2 * np.pi, 400)
    ax.plot(r * np.cos(a), r * np.sin(a), **kw)


def annulus(ax, r0, r1, **kw):
    a = np.linspace(0, 2 * np.pi, 400)
    xo, yo = r1 * np.cos(a), r1 * np.sin(a)
    xi, yi = r0 * np.cos(a[::-1]), r0 * np.sin(a[::-1])
    ax.fill(np.concatenate([xo, xi]), np.concatenate([yo, yi]), **kw)


def load_func(func):
    if func == "rho_sym_gate":                       # main cross-phase gate (audit_150)
        df = pd.read_csv(GATE_SRC)
        out = df[["patient", "band"]].copy()
        out["obs"] = df["obs_rho"].astype(float)     # per-patient rho_sym itself
        out["surr"] = df["surr_p50"].astype(float)   # matched-strength null median
        out["pp"] = df["obs_p_one_sided"].astype(float)
        return out
    df = pd.read_csv(SRC)
    out = df[["patient", "band"]].copy()
    out["obs"] = df[f"{func}_obs"].astype(float)
    out["surr"] = df[f"{func}_surr_p50"].astype(float)
    out["pp"] = df[f"{func}_p"].astype(float)
    return out


def cohort_p(func):
    if func == "rho_sym_gate":                       # gate_p_sym = Wilcoxon vs matched-strength
        v = pd.read_csv(GATE_VERDICT)
        return {r.band: float(r.gate_p_sym) for r in v.itertuples()}
    v = pd.read_csv(VERDICT)
    vv = v[v.functional == func]
    return {r.band: float(r.wilcoxon_p) for r in vv.itertuples()}


def scale_ticks():
    cand = np.round(np.arange(-0.8, 0.81, 0.2), 1)
    ticks = [float(t) for t in cand if T_LO + 0.03 < t < T_HI - 0.005 and t != 0.0]
    return ticks


# ---- shared scaffold (wash, rings, separators, ticks, hub) --------------------
def draw_scaffold(ax, df, cohp, spec):
    beta_i = BRAIN_BANDS_NAMES.index("beta")
    r0 = r_of(0.0)

    # band-gradient + cohort-significance wash (turbo palette kept faintly visible
    # everywhere for frequency-continuity; extra brightness = significance)
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        c = 90 + (beta_i - i) * 60
        wedge_gradient(ax, c, band_hue(i), 0.12 + 0.52 * sigfrac(cohp[band]))

    # radial sector separators (band boundaries)
    for bnd in (0, 60, 120, 180, 240, 300):
        a = np.deg2rad(bnd)
        ax.plot([R_HUB * np.cos(a), (R_OUT + 0.015) * np.cos(a)],
                [R_HUB * np.sin(a), (R_OUT + 0.015) * np.sin(a)],
                color="#b9bdc4", lw=0.7, alpha=0.6, zorder=0.45)

    # dashed scale rings + numeric value at ONE clean spoke (180 deg, alpha|theta gap)
    for Tg in scale_ticks():
        ring(ax, r_of(Tg), color="#d0d3d9", lw=0.5, ls=(0, (1, 6)), zorder=0.3)
    surr = df["surr"].to_numpy(float)
    nz0, nz1 = np.percentile(surr, 10), np.percentile(surr, 90)
    annulus(ax, r_of(nz0), r_of(nz1), color="#8a9096", alpha=0.22, lw=0, zorder=0.35)
    ring(ax, r0, color="#5c6166", lw=1.05, zorder=0.5)          # waterline T=0

    for Tg in [0.0] + scale_ticks():
        lab = "0" if Tg == 0.0 else f"{Tg:+.1f}"
        rr = r_of(Tg)
        ax.text(-rr, 0.0, lab, fontsize=7.0, color="#8f959d", ha="center",
                va="center", zorder=6, rotation=90,
                bbox=dict(boxstyle="round,pad=0.05", fc="white", ec="none",
                          alpha=0.75))
    # what the radial value IS (quantity), at the outer end of that spoke
    ax.text(-(R_OUT + 0.15), 0.0, spec["ring_q"], fontsize=11, color="#3f454c",
            ha="center", va="center", rotation=90, zorder=6)
    ax.text(-(R_OUT + 0.15), 0.30, "trace " + r"$\rightarrow$", fontsize=6.8,
            color="#9298a0", ha="center", va="center", rotation=90, zorder=6)

    # rim band labels + cohort p (band that clears is bold/dark)
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        theta0 = np.deg2rad(90 + (beta_i - i) * 60)
        sig = cohp[band] < 0.05
        ax.text((R_OUT + 0.13) * np.cos(theta0), (R_OUT + 0.13) * np.sin(theta0),
                BRAIN_BAND_TEX_DICT[band], fontsize=17 if sig else 14,
                ha="center", va="center", zorder=6,
                color="#141414" if sig else "#565b62",
                fontweight="bold" if sig else "normal")
        pv = cohp[band]
        ptxt = (rf"$p=\mathbf{{{pv:.3f}}}$" if sig else f"$p={pv:.2f}$")
        ax.text((R_OUT + 0.255) * np.cos(theta0), (R_OUT + 0.255) * np.sin(theta0),
                ptxt, fontsize=9.5 if sig else 8, ha="center", va="center",
                zorder=6, color="#111111" if sig else "#9aa0a7",
                fontweight="bold" if sig else "normal")

    # hub
    ax.add_patch(plt.Circle((0, 0), R_HUB * 0.72, fc="#f6f7f9", ec="#d5d8dd",
                            lw=0.7, zorder=3))
    ax.text(0, 0.030, spec["hub"], fontsize=8.6, ha="center", va="center",
            color="#3f454c", zorder=4, fontweight="bold")
    ax.text(0, -0.050, r"$\rightarrow$ offline rest", fontsize=7.4, ha="center",
            va="center", color="#5a6068", zorder=4)
    return beta_i, r0


def _set_range(df, cohp):
    """radial range scaled to the SIGNIFICANT-band signal: the outer axis is set
    by the bands that clear the cohort null, so null-band single-patient spikes
    clip to the rim (with an off-scale chevron) instead of stretching the scale."""
    global T_LO, T_HI
    o = df["obs"].to_numpy(float)
    sig = [b for b in BRAIN_BANDS_NAMES if cohp[b] < 0.05]
    # scale to the significant-band signal; if no band clears, robust p90 fallback
    hi = df[df.band.isin(sig)]["obs"].max() if sig else np.percentile(o, 90)
    T_HI = float(max(hi, 0.30)) * 1.10
    T_LO = float(min(np.percentile(o, 5), -0.05)) * 1.06


# ---- headline bloom (petal colour = band) -------------------------------------
def build_bloom(func):
    spec = FUNCS[func]
    df = load_func(func)
    cohp = cohort_p(func)
    _set_range(df, cohp)

    fig, ax = plt.subplots(figsize=(7.8, 7.8))
    ax.set_aspect("equal")
    ax.axis("off")
    beta_i, r0 = draw_scaffold(ax, df, cohp, spec)

    petal_span = np.deg2rad(60) * 0.82
    offs = np.linspace(-petal_span / 2, petal_span / 2, 10)
    half_w = (petal_span / 10) * 0.46
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        theta0 = np.deg2rad(90 + (beta_i - i) * 60)
        hue = band_hue(i)
        em = 0.42 + 0.58 * sigfrac(cohp[band])
        sub = df[df.band == band].sort_values("obs")
        obs = sub["obs"].to_numpy(float)
        pp = sub["pp"].to_numpy(float)
        # median arc BELOW the petals (zorder < petals)
        med = float(np.median(obs))
        aa = np.linspace(theta0 - petal_span / 2, theta0 + petal_span / 2, 60)
        rr = r_of(med)
        sig = cohp[band] < 0.05
        ax.plot(rr * np.cos(aa), rr * np.sin(aa), color=hue if sig else "#8f949b",
                lw=3.4 if sig else 1.7, alpha=0.98 if sig else 0.5, zorder=1.7,
                solid_capstyle="round")
        for k in range(len(obs)):
            pos = obs[k] > 0
            g = float(np.clip((0.5 - pp[k]) / 0.5, 0.0, 1.0))
            base = hue if pos else C_ANTI
            th = theta0 + offs[k]
            r1 = r_of(obs[k])
            if pos and g > 0.4:
                pg = petal_outline(th, half_w * 1.85, r0, r1)
                ax.add_patch(PathPatch(pg, fc=base, ec="none",
                                       alpha=(0.05 + 0.22 * g) * em, zorder=2.0))
            p = petal_outline(th, half_w, r0, r1)
            fa = ((0.34 + 0.5 * g) if pos else (0.28 + 0.12 * g)) * em
            ax.add_patch(PathPatch(p, fc=(*base, fa), ec=(*base, 0.85 * em),
                                   lw=0.5 + 0.5 * g, zorder=2.5 + 0.5 * g,
                                   joinstyle="round"))
            if obs[k] > T_HI or obs[k] < T_LO:
                clip_cap(ax, th, base, outward=obs[k] > 0)

    lim = R_OUT + 0.36
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    out = OUTDIR / f"{spec['stem']}_DRAFT.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    _report(func, df, cohp)
    print(f"wrote {out}")


# ---- coverage-annotated diagnostic (petal colour = dominant gray territory) ----
def dominant_macro():
    pn = pd.read_csv(PN_COV)
    g0 = pn[(pn.band == "beta") & (pn.tissue == "gray")]
    macro = {}
    for p, g in g0.groupby("patient"):
        cnt = g.supersystem.map(SUPER2MACRO).dropna().value_counts()
        macro[p] = cnt.index[0] if len(cnt) else "limbic"
    return macro


def build_coverage(func):
    spec = FUNCS[func]
    df = load_func(func)
    cohp = cohort_p(func)
    _set_range(df, cohp)
    macro = dominant_macro()
    cohort = sorted(df.patient.unique(),
                    key=lambda p: (MACRO_ORDER.index(macro[p]), p))

    fig, ax = plt.subplots(figsize=(7.9, 8.3))
    ax.set_aspect("equal")
    ax.axis("off")
    beta_i, r0 = draw_scaffold(ax, df, cohp, spec)

    petal_span = np.deg2rad(60) * 0.82
    offs = np.linspace(-petal_span / 2, petal_span / 2, len(cohort))
    half_w = (petal_span / len(cohort)) * 0.46
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        theta0 = np.deg2rad(90 + (beta_i - i) * 60)
        sub = df[df.band == band].set_index("patient")
        med = float(np.median(sub.loc[cohort, "obs"]))
        aa = np.linspace(theta0 - petal_span / 2, theta0 + petal_span / 2, 50)
        rr = r_of(med)
        sig = cohp[band] < 0.05
        ax.plot(rr * np.cos(aa), rr * np.sin(aa),
                color="#2c3036" if sig else "#9aa0a7", lw=2.8 if sig else 1.3,
                alpha=0.9 if sig else 0.5, zorder=1.7, solid_capstyle="round")
        for slot, p in enumerate(cohort):
            obs = float(sub.loc[p, "obs"])
            g = float(np.clip((0.5 - float(sub.loc[p, "pp"])) / 0.5, 0.0, 1.0))
            col = MACRO_COL[macro[p]]
            th = theta0 + offs[slot]
            r1 = r_of(obs)
            pos = obs > 0
            if pos and g > 0.4:
                pg = petal_outline(th, half_w * 1.85, r0, r1)
                ax.add_patch(PathPatch(pg, fc=col, ec="none",
                                       alpha=0.06 + 0.20 * g, zorder=2.0))
            pt = petal_outline(th, half_w, r0, r1)
            fa = (0.32 + 0.55 * g) if pos else (0.20 + 0.15 * g)
            ax.add_patch(PathPatch(pt, fc=col, alpha=fa, ec="white", lw=0.4,
                                   zorder=2.5 + 0.5 * g))
            if obs > T_HI or obs < T_LO:
                clip_cap(ax, th, col, outward=obs > 0)

    present = [m for m in MACRO_ORDER if m in set(macro.values())]
    handles = [Line2D([], [], marker="^", ls="", mfc=MACRO_COL[m], mec="white",
                      ms=12, label=MACRO_LAB[m]) for m in present]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.045),
               ncol=len(present), frameon=False, fontsize=8.8,
               title="petal colour = dominant gray-matter territory (implant coverage)",
               title_fontsize=8.0)
    fig.text(0.5, 0.018, "slot fixed across bands  ·  length = trace, inward = "
             "reverted  ·  wash = cohort significance  ·  EXPLORATORY",
             ha="center", va="center", fontsize=7.0, color="#9298a0")

    lim = R_OUT + 0.36
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    out = OUTDIR / f"{spec['stem']}_coverage_DRAFT.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def _report(func, df, cohp):
    print(f"\n=== {func} ===  T range [{T_LO:+.2f},{T_HI:+.2f}]")
    print("band        median   p<.05/10   wilcox_p")
    for band in BRAIN_BANDS_NAMES:
        g = df[df.band == band]
        med = float(g["obs"].median())
        nsig = int((g["pp"] < 0.05).sum())
        star = "  <== BLOOM" if cohp[band] < 0.05 else ""
        print(f"  {band:11s} {med:+.3f}    {nsig}/10       {cohp[band]:.4f}{star}")


# ---- DUAL-NULL bloom (must clear matched-strength AND drift) -------------------
# One picture per measure. Per-patient triangles as before; the matched-strength
# null is the grey annulus; the DRIFT null is a red dashed waterline per wedge.
# A wedge glows only if the band clears BOTH gates -> effective p = max(p_MS, p_drift)
# (weakest link). Raw clears MS in beta but NOT drift -> drowns; coph/Grassmann
# clear both only at beta -> beta is the sole band lit in both multiscale probes.
DUAL = {
    "raw": dict(hub="raw edges", ring_q=r"$\rho^{\mathrm{raw}}$",
                stem="fig_dual_bloom_raw"),
    "coph": dict(hub="cophenetic", ring_q=r"$\rho^{\mathrm{coph}}$",
                 stem="fig_dual_bloom_coph"),
    "grassmann": dict(hub="Grassmann", ring_q=r"$T_G$",
                      stem="fig_dual_bloom_grassmann"),
}


def _drift_gate(col_real, col_drift):
    """cohort drift gate (paired one-sided Wilcoxon real>drift) + drift median."""
    d = pd.read_csv(DRIFT_LADDER)
    dp, dm = {}, {}
    for b in BRAIN_BANDS_NAMES:
        g = d[d.band == b]
        rr = g[col_real].to_numpy(float)
        dd = g[col_drift].to_numpy(float)
        try:
            dp[b] = float(wilcoxon(rr - dd, alternative="greater").pvalue)
        except ValueError:
            dp[b] = 1.0
        dm[b] = float(np.median(dd))
    return dp, dm


def load_dual(measure):
    spec = DUAL[measure]
    if measure == "coph":
        pp = pd.read_csv(GATE_SRC)
        df = pp[["patient", "band"]].copy()
        df["obs"] = pp["obs_rho"].astype(float)
        df["surr"] = pp["surr_p50"].astype(float)
        df["pp"] = pp["obs_p_one_sided"].astype(float)
        cs = pd.read_csv(GATE_VERDICT)
        ms_p = {r.band: float(r.gate_p_sym) for r in cs.itertuples()}
        drift_p, drift_med = _drift_gate("coph_real_Ttest", "coph_drift_Ttest")
    elif measure == "raw":
        lad = pd.read_csv(PAIR_LADDER)
        lad = lad[lad.descriptor == "raw_fc"]
        df = lad[["patient", "band"]].copy()
        df["obs"] = lad["T_test_obs"].astype(float)
        df["surr"] = lad["T_test_surr_p50"].astype(float)
        df["pp"] = lad["T_test_p"].astype(float)
        cs = pd.read_csv(PAIR_COHORT)
        cs = cs[(cs.descriptor == "raw_fc") & (cs.functional == "T_test")]
        ms_p = {r.band: float(r.gate_p) for r in cs.itertuples()}
        drift_p, drift_med = _drift_gate("raw_real_Ttest", "raw_drift_Ttest")
    else:  # grassmann
        arc = pd.read_csv(GARC)
        arc = arc[arc.functional == "onl"]
        df = arc[["patient", "band"]].copy()
        df["obs"] = arc["median_obs_over_k"].astype(float)
        df["surr"] = (arc["median_obs_over_k"]
                      - arc["median_obs_minus_surr_over_k"]).astype(float)
        ppmap = pd.read_csv(GMS_PERK).groupby(
            ["patient", "band"])["obs_p_one_sided_upper"].median()
        df["pp"] = [float(ppmap.get((r.patient, r.band), 0.5))
                    for r in df.itertuples()]
        ms_p = {r.band: float(r.cluster_p_cluster_mass)
                for r in pd.read_csv(GCLUSTER).itertuples()}
        gd = pd.read_csv(GDRIFT)
        drift_p = {r.band: float(r.drift_gate_p) for r in gd.itertuples()}
        drift_med = {r.band: float(r.drift_median) for r in gd.itertuples()}
    eff = {b: max(ms_p[b], drift_p[b]) for b in BRAIN_BANDS_NAMES}
    return df, ms_p, drift_p, drift_med, eff, spec


def build_dual_bloom(measure):
    global T_LO, T_HI
    df, ms_p, drift_p, drift_med, eff, spec = load_dual(measure)
    _set_range(df, eff)                                   # scale to clears-both bands
    dvals = list(drift_med.values())
    T_LO = min(T_LO, min(dvals) * 1.10)
    T_HI = max(T_HI, max(dvals) * 1.10)

    fig, ax = plt.subplots(figsize=(7.8, 7.8))
    ax.set_aspect("equal")
    ax.axis("off")
    beta_i, r0 = draw_scaffold(ax, df, eff, spec)         # wedge glow = clears-BOTH

    petal_span = np.deg2rad(60) * 0.82
    offs = np.linspace(-petal_span / 2, petal_span / 2, 10)
    half_w = (petal_span / 10) * 0.46
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        theta0 = np.deg2rad(90 + (beta_i - i) * 60)
        hue = band_hue(i)
        both = eff[band] < 0.05
        em = 0.42 + 0.58 * sigfrac(eff[band])
        aa = np.linspace(theta0 - petal_span / 2, theta0 + petal_span / 2, 60)
        # DRIFT waterline (second bar the median must clear)
        rd = r_of(drift_med[band])
        ax.plot(rd * np.cos(aa), rd * np.sin(aa), color=DRIFT_COL, lw=1.7,
                ls=(0, (4, 2)), alpha=0.9, zorder=1.55, solid_capstyle="round")
        # median obs arc: band-hued+bold only if it clears BOTH nulls
        sub = df[df.band == band].sort_values("obs")
        obs = sub["obs"].to_numpy(float)
        pp = sub["pp"].to_numpy(float)
        med = float(np.median(obs))
        rr = r_of(med)
        ax.plot(rr * np.cos(aa), rr * np.sin(aa), color=hue if both else "#8f949b",
                lw=3.6 if both else 1.7, alpha=0.98 if both else 0.5, zorder=1.7,
                solid_capstyle="round")
        for k in range(len(obs)):
            pos = obs[k] > 0
            g = float(np.clip((0.5 - pp[k]) / 0.5, 0.0, 1.0))
            base = hue if pos else C_ANTI
            th = theta0 + offs[k]
            r1 = r_of(obs[k])
            if pos and g > 0.4:
                pg = petal_outline(th, half_w * 1.85, r0, r1)
                ax.add_patch(PathPatch(pg, fc=base, ec="none",
                                       alpha=(0.05 + 0.22 * g) * em, zorder=2.0))
            p = petal_outline(th, half_w, r0, r1)
            fa = ((0.34 + 0.5 * g) if pos else (0.28 + 0.12 * g)) * em
            ax.add_patch(PathPatch(p, fc=(*base, fa), ec=(*base, 0.85 * em),
                                   lw=0.5 + 0.5 * g, zorder=2.5 + 0.5 * g,
                                   joinstyle="round"))
            if obs[k] > T_HI or obs[k] < T_LO:
                clip_cap(ax, th, base, outward=obs[k] > 0)

    fig.text(0.5, 0.05, "petals = patients   ·   grey band = matched-strength null   "
             "·   red dashed = drift null   ·   a wedge glows only if it clears BOTH",
             ha="center", va="center", fontsize=7.4, color="#9298a0")
    lim = R_OUT + 0.36
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    out = OUTDIR / f"{spec['stem']}_DRAFT.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"\n=== dual {measure} ===  T[{T_LO:+.2f},{T_HI:+.2f}]")
    print("band        MS_p   drift_p   eff_p   clears")
    for b in BRAIN_BANDS_NAMES:
        print(f"  {b:11s}{ms_p[b]:6.3f}{drift_p[b]:9.3f}{eff[b]:8.3f}"
              f"   {'BOTH' if eff[b] < 0.05 else '-'}")
    print(f"wrote {out}")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--func", nargs="*", default=list(FUNCS),
                    choices=list(FUNCS), help="which functional(s) to build")
    ap.add_argument("--no-coverage", action="store_true",
                    help="build only the headline bloom, skip the coverage diagnostic")
    ap.add_argument("--dual", nargs="*", choices=list(DUAL), default=None,
                    help="build dual-null (MS + drift) blooms for these measures")
    args = ap.parse_args()
    if args.dual is not None:
        for m in (args.dual or list(DUAL)):
            build_dual_bloom(m)
    else:
        for f in args.func:
            build_bloom(f)
            if not args.no_coverage:
                build_coverage(f)
