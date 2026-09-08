#!/usr/bin/env python3
r"""talk_fig_trace_multiscale — the slide-13 ("A lasting, multiscale trace") result figure set.

Two talk-clean figures, both from the settled mst@0.20 matched-strength cache
(``data/sparsified_arc/ms_mst020/``). Redesigned 2026-07-13 (user): the multiscale story
is carried by the (band x scale) gate MAP, and the per-patient spread is folded INTO the
scale sweep as a split violin, so both figures now express the tau-dependence honestly.

  Fig 1  trace_scale_map.png
         the (band x scale) trace map: -log10(cohort matched-strength gate p) over all 16
         diffusion scales, one row per band, with the p=0.05 contour drawn. This is the
         honest per-scale gate (NO scalar collapse): beta lights up scale-BROAD, alpha a
         mesoscale band, theta / low_gamma stay dark (NO COHORT CONSENSUS -- not trace-free;
         individual patients still trace, low_gamma strongly), delta / high_gamma flare only in
         the coarse-collapse tail. Answers WHERE the cohort-consistent trace lives.

  Fig 2  trace_band_scale_violin.png
         per-band SPLIT violin across the tau-sweep. For each band the LEFT half is the
         10-patient distribution at that band's LEAST-traced scale (argmin cohort-median
         Delta) and the RIGHT half is the distribution at its MOST-traced scale (argmax),
         plotted as Delta = rho_sym - matched-strength null so 0 == the null. Only the two
         medians survive, each carrying its scale label. This is the tau-dependence ENVELOPE
         (worst..best), NOT a best-scale verdict: the verdict is the gate breadth N/16 (the
         tag under each band) and Fig 1. The one band whose BOTH halves clear the null is
         beta = scale-invariant; alpha (and the patient-specific delta/high_gamma) clear
         only at their best scale; theta/low_gamma never reach cohort consensus -- heterogeneous
         (strong but inconsistent per-patient traces), NOT trace-free.

  Fig 3  trace_band_scale_curves.png
         per-band rho_sym(s) CURVES over the tau-sweep, one per band, with the same full/empty
         gate grammar as the encoding/inference scale-payoff figure (fig_encinf_scale_payoff):
         FILLED marker = the cohort clears the matched-strength gate at that scale (p<0.05),
         OPEN marker = ns. The MAGNITUDE companion to Fig 1 -- the gate map shows -log10(p), this
         shows the actual cross-phase similarity and its scale profile, so band-selectivity reads
         at a glance: beta rides high and all-filled (16/16, scale-invariant), alpha high with a
         mesoscale filled run (12/16), delta/high_gamma partial, theta/low_gamma low and all-open
         (0/16, NO COHORT CONSENSUS -- not trace-free; the per-patient spread is Fig 2). Grey floor
         = matched-strength null envelope.

Register: cohort claim = the matched-strength Wilcoxon gate READ PER SCALE (Fig 1 / Fig 3 dots /
the N/16 tags), NOT the violin/curve height, NOT a best-tau scalar. beta is DELOCALIZED (no anatomy
here).

Writes: data/outputs/figures/talk/trace_scale_map.{png,pdf}
        data/outputs/figures/talk/trace_band_scale_violin.{png,pdf}
        data/outputs/figures/talk/trace_band_scale_curves.{png,pdf}
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patheffects as pe
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

_STROKE = [pe.withStroke(linewidth=2.4, foreground="white")]      # legibility on any band colour

sys.path.insert(0, str(Path(__file__).resolve().parents[1] /
                       "01_compute" / "figures_embedded" / "new_results_sec1"))
import _common as C  # noqa: E402

C.use_lrg_style()

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
CARRIER = ("beta", "alpha")                       # cohort carriers (bold)
NOCONSENSUS = ("theta", "low_gamma")              # no cohort consensus: theta hugs the floor,
#                                                   low_gamma is heterogeneous (strong but inconsistent)
OUTDIR = C.ROOT / "data" / "outputs" / "figures" / "talk"
GREY = "0.5"


def _load_gate():
    return pd.read_csv(C.MS / "cohort_gate.csv")


def _sig_count(g):
    return int((g.gate_p < 0.05).sum())


def _fmt_s(s):
    """Compact scale label."""
    return f"$s{{=}}{s:.1f}$" if s < 10 else f"$s{{=}}{s:.0f}$"


# ─────────────────────────── Fig 1 — the (band × scale) trace map ───────────────────────────
def fig_scale_map():
    """-log10(matched-strength gate p) over band × scale, p=0.05 contour drawn.

    The honest per-scale gate: no collapse to a single tau. beta scale-broad, alpha
    mesoscale, theta/low_gamma dark, delta/high_gamma coarse-tail only.
    """
    gate = _load_gate()
    s_vals = np.sort(gate.s.unique())
    M = np.full((len(BANDS), len(s_vals)), np.nan)
    for i, b in enumerate(BANDS):
        gb = gate[gate.band == b]
        for _, r in gb.iterrows():
            j = int(np.argmin(np.abs(s_vals - r.s)))
            M[i, j] = -np.log10(max(r.gate_p, 1e-4))

    fig, ax = plt.subplots(figsize=(11.0, 5.0))
    x = np.log10(s_vals)
    dx = x[-1] - x[0]
    extent = [x[0], x[-1], len(BANDS) - 0.5, -0.5]
    P05 = -np.log10(0.05)
    im = ax.imshow(M, aspect="auto", cmap="magma", norm=Normalize(0, 3), extent=extent,
                   interpolation="nearest", zorder=1)

    # desaturate the sub-threshold (p >= 0.05) cells so the lit region = the trace zone
    veil = np.zeros((*M.shape, 4))
    veil[..., :3] = 0.92
    veil[..., 3] = np.where(M < P05, 0.68, 0.0)
    ax.imshow(veil, aspect="auto", extent=extent, interpolation="nearest", zorder=2)

    # scale landmarks
    ax.axvline(0.0, color="0.25", lw=1.0, ls=":", alpha=0.8, zorder=3)               # s = tau_min
    ax.axvline(np.log10(C.S_REPORT), color="#1f9fd0", lw=1.6, ls="--", alpha=0.95, zorder=3)

    ax.set_yticks(range(len(BANDS)))
    ax.set_yticklabels([C.BTeX.get(b, b) for b in BANDS], fontsize=17)
    for t, b in zip(ax.get_yticklabels(), BANDS):
        if b in CARRIER:
            t.set_fontweight("bold")
    ax.set_xlabel(r"$\log_{10}\,s=\log_{10}(\tau\lambda_{\max})$      fine $\rightarrow$ coarse",
                  fontsize=13)
    ax.tick_params(axis="x", labelsize=11)

    cb = fig.colorbar(im, ax=ax, pad=0.03, fraction=0.043)
    cb.set_label(r"$-\log_{10}\,p_{\mathrm{matched\text{-}strength}}$", fontsize=11)
    cb.ax.axhline(P05, color="white", lw=1.8)
    cb.ax.text(0.5, P05, "$p{=}0.05$", transform=cb.ax.get_yaxis_transform(),
               ha="center", va="center", fontsize=8.5, color="0.15", clip_on=False,
               bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))

    _save(fig, "trace_scale_map")
    for b in BANDS:
        print(f"  {b:11s} {_sig_count(gate[gate.band == b])}/16 scales sig")


# ───────────────────────── Fig 2 — per-band split violin over the τ-sweep ─────────────────────────
def _trim(vals, q=(10.0, 90.0)):
    """Central [q_lo, q_hi] percentile window — kills single-outlier violin necks (n=10)."""
    return float(np.percentile(vals, q[0])), float(np.percentile(vals, q[1]))


def _kde_half(vals, width):
    """KDE half-profile clipped to the central-percentile data window (no fake outlier
    necks, no KDE overshoot). Returns (yg, dens)."""
    vals = np.asarray(vals, float)
    vals = vals[~np.isnan(vals)]
    if len(vals) < 2 or np.ptp(vals) < 1e-9:
        return None, None
    lo, hi = _trim(vals)
    if hi - lo < 1e-9:
        return None, None
    yg = np.linspace(lo, hi, 200)
    dens = gaussian_kde(vals, bw_method=0.35)(yg)
    return yg, dens / dens.max() * width


def _half(ax, vals, xc, side, color, width, *, filled, carrier, sig):
    """One half of the split violin. filled=True -> shaded 'violin'; False -> outline 'line'.
    The median is drawn SOLID + filled dot when that scale clears the matched-strength gate
    (sig=True) and DASHED + open ring when it does not (sig=False) — so 'above the null line'
    can never be mistaken for 'traces'. Returns (median, density_at_median)."""
    vals = np.asarray(vals, float); vals = vals[~np.isnan(vals)]
    med = float(np.median(vals))
    sign = -1.0 if side == "left" else 1.0
    yg, dens = _kde_half(vals, width)
    dm = width * 0.55 if yg is None else float(np.interp(med, yg, dens))
    if yg is not None:
        if filled:
            ax.fill_betweenx(yg, xc, xc + sign * dens, color=color,
                             alpha=0.55 if carrier else 0.22, lw=0, zorder=3)
            ax.plot(xc + sign * dens, yg, color=color, lw=1.4 if carrier else 0.9,
                    alpha=0.85, zorder=3.4)
        else:                                              # 'line' half = outline only
            ax.plot(xc + sign * dens, yg, color=color, lw=2.6 if carrier else 1.4,
                    alpha=0.95 if carrier else 0.55, zorder=3.6, solid_capstyle="round")
    lw = 4.0 if carrier else 2.5
    if sig:                                                # clears the gate -> solid + filled dot
        ax.plot([xc, xc + sign * dm], [med, med], color=color, lw=lw, zorder=5,
                solid_capstyle="round", path_effects=_STROKE)
        ax.plot([xc + sign * dm], [med], "o", ms=15 if carrier else 12, mfc=color,
                mec="white", mew=1.6, zorder=6)
    else:                                                  # fails the gate -> dashed + open ring
        ax.plot([xc, xc + sign * dm], [med, med], color=color, lw=lw * 0.55,
                ls=(0, (2.2, 1.7)), alpha=0.7, zorder=5)
        ax.plot([xc + sign * dm], [med], "o", ms=13 if carrier else 10.5, mfc="white",
                mec=color, mew=2.6, zorder=6)
    return med, dm


def _rank_scales(gate, band):
    """(s, p) of the STRONGEST-trace scale (min gate p) and WEAKEST-trace scale (max gate p)."""
    g = gate[gate.band == band].sort_values("s").reset_index(drop=True)
    i_hi = int(np.argmin(g.gate_p.values))                # strongest trace
    i_lo = int(np.argmax(g.gate_p.values))                # weakest trace
    return (float(g.s[i_hi]), float(g.gate_p[i_hi])), (float(g.s[i_lo]), float(g.gate_p[i_lo]))


def fig_band_scale_violin():
    """Split violin per band over the tau-sweep, halves RANKED BY THE GATE (not by Delta, which
    is upward-biased and drifts to the coarse-collapse tail). LEFT filled 'violin' = the band's
    WEAKEST-trace scale (max gate p); RIGHT outline 'line' = its STRONGEST-trace scale (min gate
    p). y = Delta = obs_rho - matched-strength null so 0 == null. Each median is SOLID+dot if it
    clears the gate at that scale, DASHED+ring if not: beta = two solid (scale-invariant),
    theta/low_gamma = two rings (NO COHORT CONSENSUS -- traces present per patient but they don't
    cohere; NOT trace-free), alpha/delta/high_gamma = one of each. The ring means the cohort didn't
    hold, NEVER 'no trace'. The verdict is the gate, shown BOTH as the dot/ring and the N/16."""
    pp = pd.read_csv(C.MS / "per_patient_scale.csv")
    gate = _load_gate()

    def dvals(band, s):
        return (pp[(pp.band == band) & np.isclose(pp.s, s)].eval("obs_rho - surr_p50")).values

    # gate-ranked strongest/weakest scale per band + trimmed ylim
    sel = {}
    los, his = [], []
    for band in BANDS:
        (s_hi, p_hi), (s_lo, p_lo) = _rank_scales(gate, band)
        sel[band] = dict(s_hi=s_hi, p_hi=p_hi, s_lo=s_lo, p_lo=p_lo)
        for s in (s_lo, s_hi):
            lo, hi = _trim(dvals(band, s))
            los.append(lo); his.append(hi)
    ylo, yhi = min(los) - 0.05, max(his) + 0.17           # headroom for the reading key
    ytag = ylo + 0.015

    fig, ax = plt.subplots(figsize=(11.4, 6.0))
    ax.axhline(0.0, color="0.45", lw=1.4, ls="--", zorder=1)
    ax.text(-0.66, 0.006, "matched-strength null", color="0.45", fontsize=10,
            va="bottom", ha="left", style="italic")

    width = 0.42
    for k, band in enumerate(BANDS):
        col = C.band_color(band)
        carrier = band in CARRIER
        S = sel[band]
        sig_hi, sig_lo = S["p_hi"] < 0.05, S["p_lo"] < 0.05

        m_lo, dm_lo = _half(ax, dvals(band, S["s_lo"]), k, "left", col, width,
                            filled=True, carrier=carrier, sig=sig_lo)     # weakest = filled
        m_hi, dm_hi = _half(ax, dvals(band, S["s_hi"]), k, "right", col, width,
                            filled=False, carrier=carrier, sig=sig_hi)    # strongest = outline

        # scale label riding INSIDE each half, just above the median (dark grey = legible on any hue)
        fs = 10.5 if carrier else 9
        lab_c = "0.12" if carrier else "0.4"
        ax.text(k - dm_lo * 0.5, m_lo + 0.016, _fmt_s(S["s_lo"]), va="bottom", ha="center",
                fontsize=fs, color=lab_c, fontweight="bold" if carrier else "normal",
                path_effects=_STROKE, zorder=6)
        ax.text(k + dm_hi * 0.5, m_hi + 0.016, _fmt_s(S["s_hi"]), va="bottom", ha="center",
                fontsize=fs, color=lab_c, fontweight="bold" if carrier else "normal",
                path_effects=_STROKE, zorder=6)

        # gate breadth = the cohort verdict (below each band)
        n16 = _sig_count(gate[gate.band == band])
        ax.text(k, ytag, f"{n16}/16", va="bottom", ha="center",
                fontsize=12 if carrier else 10, color=col if carrier else "0.5",
                fontweight="bold" if carrier else "normal", zorder=6, path_effects=_STROKE)

    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([C.BTeX.get(b, b) for b in BANDS], fontsize=18)
    for t, b in zip(ax.get_xticklabels(), BANDS):
        if b in CARRIER:
            t.set_fontweight("bold")
    ax.set_xlim(-0.72, len(BANDS) - 0.28)
    ax.set_ylim(ylo, yhi)
    ax.set_ylabel(C.YLAB_RHO + r"$\,-\,$null   (per patient)", fontsize=13)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # compact reading key (top-left, above the violins)
    ax.text(-0.66, yhi,
            r"$\blacktriangleleft$ filled $=$ weakest  ·  "
            r"outline $\blacktriangleright$ $=$ strongest scale  ·  violin $=$ all 10 patients"
            "\n" r"$\bullet$ cohort clears matched-strength  ·  $\circ$ no cohort consensus  ·  "
            r"$N/16$ $=$ scales the cohort holds",
            va="top", ha="left", fontsize=9.5, color="0.4")

    _save(fig, "trace_band_scale_violin")
    for band in BANDS:
        S = sel[band]
        print(f"  {band:11s} strongest s={S['s_hi']:6.1f} p={S['p_hi']:.3f}"
              f"{'●' if S['p_hi']<0.05 else '○'}   "
              f"weakest s={S['s_lo']:6.1f} p={S['p_lo']:.3f}{'●' if S['p_lo']<0.05 else '○'}")


# ─────────────── Fig 3 — per-band ρ_sym(s) curves over the τ-sweep (filled/open gate) ───────────────
def fig_band_scale_curves():
    r"""Per-band cophenetic trace ρ_sym as a function of diffusion scale s = τ·λ_max, one curve per
    band, **filled marker = the cohort clears the matched-strength gate at that scale (p<0.05), open
    marker = ns** — the same full/empty grammar as the encoding/inference scale-payoff figure
    (``fig_encinf_scale_payoff``). This is the MAGNITUDE companion to Fig 1: the gate map shows
    −log₁₀(p), this shows the actual cross-phase similarity and its scale profile, so band-selectivity
    reads at a glance — β rides high and all-filled (16/16, scale-invariant), α high with a filled
    mesoscale run (12/16), δ / γ_high partial, θ / γ_low low and all-open (0/16, no cohort consensus —
    NOT trace-free; the per-patient spread is Fig 2). y = raw ρ_sym (whole-graph reinstatement, whose
    ρ_sym approaches the split-half reliability ceiling); grey floor = matched-strength null envelope.
    Register: the ●/○ IS the per-scale gate (never a curve-height read, never a patient count)."""
    gate = _load_gate()
    s_all = np.sort(gate.s.unique())
    # matched-strength null floor: envelope of the per-band surrogate medians (grey = null level)
    surr_env = np.array([gate[np.isclose(gate.s, s)].surr_med.max() for s in s_all])

    fig, ax = plt.subplots(figsize=(11.3, 5.7))
    ax.axhline(0.0, color="0.82", lw=0.8, zorder=0)
    ax.fill_between(s_all, 0, np.clip(surr_env, 0, None), color="0.72", alpha=0.35, lw=0, zorder=1)
    ax.plot(s_all, surr_env, color="0.6", lw=1.0, zorder=1)
    ax.axvline(C.S_REPORT, color="0.78", lw=1.1, ls=(0, (2, 3)), zorder=0)                    # s_report

    # bands drawn weakest→strongest so the carriers sit on top; legend ranked by breadth N/16
    order = sorted(BANDS, key=lambda b: _sig_count(gate[gate.band == b]))
    handles = []
    for band in order:
        gb = gate[gate.band == band].sort_values("s")
        ss, obs, gp = gb.s.values, gb.obs_med.values, gb.gate_p.values
        col = C.band_color(band)
        carrier = band in CARRIER
        sig = gp < 0.05
        n16 = int(sig.sum())
        ax.plot(ss, obs, "-", color=col, lw=3.2 if carrier else 1.9,
                alpha=0.97 if carrier else 0.78, zorder=(5 if carrier else 3),
                solid_capstyle="round")
        ax.scatter(ss[sig], obs[sig], s=78 if carrier else 54, color=col, edgecolor="white",
                   lw=1.6, zorder=(7 if carrier else 6))                       # clears gate → filled
        ax.scatter(ss[~sig], obs[~sig], s=54 if carrier else 40, facecolor="white", edgecolor=col,
                   lw=1.9, alpha=0.95, zorder=(6 if carrier else 4))           # ns → open ring
        handles.append(Line2D([0], [0], color=col, lw=3.0 if carrier else 2.0, marker="o",
                              mfc=col, mec="white", ms=9 if carrier else 7.5,
                              label=rf"{C.BTeX.get(band, band)}  {n16}/16"))

    ax.set_xscale("log")
    ax.set_xlim(0.9, 210)
    ax.set_xticks([1, 10, 100]); ax.set_xticklabels(["1", "10", "100"])
    ax.set_xlabel(r"diffusion scale  $s=\tau\lambda_{\max}$   (fine $\rightarrow$ coarse)", fontsize=13)
    ax.set_ylabel(C.YLAB_RHO + r"   (cohort median)", fontsize=13)
    ax.tick_params(labelsize=11)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # bands legend OUTSIDE the axes (ranked by breadth, carriers first) — colour + N/16
    leg = fig.legend(handles=handles[::-1], loc="center left", bbox_to_anchor=(0.795, 0.60),
                     frameon=False, fontsize=14, handletextpad=0.6, labelspacing=0.55,
                     title="band   scales/16", title_fontsize=11)
    leg._legend_box.align = "left"
    # filled/open + null reading key
    key = [
        Line2D([0], [0], color="0.35", lw=0, marker="o", mfc="0.35", mec="white", ms=9,
               label=r"cohort clears matched-strength ($p<0.05$)"),
        Line2D([0], [0], color="0.35", lw=0, marker="o", mfc="white", mec="0.35", mew=1.9, ms=9,
               label="ns (no cohort consensus)"),
        Line2D([0], [0], color="0.66", lw=6, alpha=0.5, label="matched-strength null"),
    ]
    fig.legend(handles=key, loc="lower center", bbox_to_anchor=(0.45, -0.02), ncol=3,
               frameon=False, fontsize=10, handletextpad=0.5, columnspacing=1.5)
    fig.subplots_adjust(left=0.085, right=0.79, top=0.975, bottom=0.16)

    _save(fig, "trace_band_scale_curves")
    for band in BANDS:
        gb = gate[gate.band == band].sort_values("s")
        n16 = int((gb.gate_p < 0.05).sum())
        print(f"  {band:11s} {n16}/16  ρ_sym peak={gb.obs_med.max():.3f} "
              f"@s={gb.s.values[int(np.argmax(gb.obs_med.values))]:.1f}  "
              f"(s_report ρ_sym={gb[np.isclose(gb.s, C.S_REPORT)].obs_med.values[0]:.3f})")


def _save(fig, stem):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUTDIR / f"{stem}.{ext}", transparent=True, bbox_inches="tight",
                    pad_inches=0.05, dpi=300)
    plt.close(fig)
    print(f"wrote {OUTDIR / stem}.png  +  .pdf")


if __name__ == "__main__":
    fig_scale_map()
    fig_band_scale_violin()
    fig_band_scale_curves()
