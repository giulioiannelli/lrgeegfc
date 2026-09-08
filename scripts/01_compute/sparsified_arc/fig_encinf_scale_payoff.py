#!/usr/bin/env python3
r"""fig_encinf_scale_payoff -- the scientific payoff of the encoding/inference slide:
encoding persists SCALE-BROADLY, inference-specific persistence is MESOSCALE-EMERGENT.

HEAD: the inferred order (encoding partialled out) leaves NO significant trace at the finest
scale (beta s1 gate_p=.097, alpha s1 gate_p=.28) and becomes significant only as the network
coarse-grains -- the signature of an abstraction (a re-grouping that exists only when you zoom
out), in BOTH alpha and beta (7/16 scales each). Encoding, by contrast, persists at EVERY scale
in beta (16/16). Read per-scale, null-referenced; matched-strength is the only null.

Two panels (beta, alpha). Per diffusion scale s = tau*lambda_max, the cohort-median concordance
with what persists into rest -- encoding T_learn = 1/2[rho(e,p)+rho(e2,p2)] and inference-specific
T_infspec_pe = 1/2[pr(f,p|e)+pr(f,p2|e2)] -- over the matched-strength null (grey floor = surrogate
median). Filled marker = clears the gate at that scale (gate_p < 0.05); open marker = ns.

Y-AXIS is expressed as a FRACTION of the split-half reliability rel(s) = Spearman(D_A, D_B)
(the same reachable-ceiling anchor the measure-slide bars are cut at; rel_vs_s.csv). Every plotted
quantity -- observed curves AND the matched-strength null -- is divided by rel(s), so 1.0 = the
reliability ceiling and the null floor stays honest. gate_p is UNCHANGED (dividing obs and null by
the same per-scale positive rel leaves the per-patient ranking, hence significance, identical).
NB the curves peak at ~0.5 (encoding) / ~0.25 (inference), NOT 1: these are sub-components that
reach only a fraction of reliability -- unlike the whole-graph reinstatement (measure slide), whose
rho_sym ~ rel so it saturates the ceiling. Raw peaks: beta enc 0.29, beta inf 0.16.

Reads : data/sparsified_arc/enc_inf_arc_mst020/cohort_gate_vs_s.csv
        data/sparsified_arc/enc_inf_arc_mst020/rel_vs_s.csv   (split-half reliability rel(s))
Writes: data/preprint/figures/new_results_sec2/fig_encinf_scale_payoff.pdf
        data/outputs/figures/talk/slide15_encinf_scale_payoff.png   (transparent)
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.visuals.styles import use_lrg_style  # noqa: E402

# encoding / inference palette -- matches the branch-origin tree + chord network (one visual
# language across the slide): amber = encoding, teal = inference-specific.
AMBER = "#e0851a"
TEAL = "#1699ab"
NULLC = "0.70"

CSV = ROOT / "data" / "sparsified_arc" / "enc_inf_arc_mst020" / "cohort_gate_vs_s.csv"
REL = ROOT / "data" / "sparsified_arc" / "enc_inf_arc_mst020" / "rel_vs_s.csv"
PDF = ROOT / "data" / "preprint" / "figures" / "new_results_sec2" / "fig_encinf_scale_payoff.pdf"
PNG = ROOT / "data" / "outputs" / "figures" / "talk" / "slide15_encinf_scale_payoff.png"

BANDS = [("beta", r"$\beta$  (13–30 Hz)"), ("alpha", r"$\alpha$  (8–13 Hz)")]
SERIES = [("T_learn", "encoding", AMBER), ("T_infspec_pe", "inference-specific", TEAL)]
S_REPORT = 5.646


def _sig_columns(s):
    """Per-scale log-space cell edges, so a shaded column at scale s[i] tiles the axis and
    adjacent significant columns MERGE into one continuous window (isolated ones stay isolated)."""
    ls = np.log(np.asarray(s, float))
    mid = 0.5 * (ls[:-1] + ls[1:])
    lo = np.concatenate([[ls[0] - (mid[0] - ls[0])], mid])
    hi = np.concatenate([mid, [ls[-1] + (ls[-1] - mid[-1])]])
    return np.exp(lo), np.exp(hi)


def _runs(idx):
    """Group a sorted index array into (start, end) contiguous runs, inclusive."""
    runs = []
    for i in idx:
        if runs and i == runs[-1][1] + 1:
            runs[-1][1] = i
        else:
            runs.append([i, i])
    return runs


def _rel_on(rel_df, band, s_vals):
    """Per-scale split-half reliability rel(s), nearest-matched onto s_vals."""
    r = rel_df[rel_df.band == band].sort_values("s")
    rs, rv = r.s.values, r.rel_median.values
    return np.array([rv[int(np.argmin(np.abs(rs - s)))] for s in s_vals])


def panel(ax, g, rel_df, band, title):
    sub = g[g.band == band]
    s_all = np.sort(sub.s.unique())
    rel = _rel_on(rel_df, band, s_all)                  # reachable ceiling per scale (=1 after norm)
    # matched-strength null floor: surrogate median enveloped over the two plotted functionals,
    # expressed (like every curve here) as a FRACTION of the split-half reliability ceiling.
    surr = np.array([sub[(sub.s == ss) & (sub.functional.isin([f for f, _, _ in SERIES]))]
                     .surr_median.max() for ss in s_all]) / rel
    ax.axhline(0.0, color="0.82", lw=0.8, zorder=0)
    ax.fill_between(s_all, 0, np.clip(surr, 0, None), color=NULLC, alpha=0.30, lw=0, zorder=1)
    ax.plot(s_all, surr, color=NULLC, lw=1.0, zorder=1)
    ax.axvline(S_REPORT, color="0.80", lw=1.0, ls=(0, (2, 3)), zorder=0)

    counts = {}
    for func, lab, col in SERIES:
        x = sub[sub.functional == func].sort_values("s")
        ss, gp = x.s.values, x.gate_p.values
        obs = x.obs_median.values / _rel_on(rel_df, band, ss)   # concordance / reliability ceiling
        sig = gp < 0.05
        counts[lab] = int(sig.sum())
        # significance WINDOW: one translucent band per CONTIGUOUS run of significant scales, in the
        # curve's own colour, alpha baked into RGBA so a region is painted exactly once (abutting
        # per-scale spans otherwise composite and saturate). beta encoding -> one band spanning almost
        # the whole axis; alpha encoding -> a couple of slivers. Behind everything.
        lo, hi = _sig_columns(ss)
        band_rgba = (*to_rgb(col), 0.13)
        for a, b in _runs(np.where(sig)[0]):
            ax.axvspan(lo[a], hi[b], facecolor=band_rgba, edgecolor="none", lw=0, zorder=-1)
        ax.plot(ss, obs, "-", color=col, lw=2.6, zorder=3, solid_capstyle="round")
        ax.scatter(ss[sig], obs[sig], s=108, color=col, edgecolor="white", lw=1.8, zorder=5)
        ax.scatter(ss[~sig], obs[~sig], s=88, facecolor="white", edgecolor=col, lw=2.1, zorder=4)
        if func == "T_infspec_pe":                      # flag the ns finest scale (the "emergence")
            ax.annotate("ns", (ss[0], obs[0]), textcoords="offset points", xytext=(-2, -13),
                        ha="center", va="top", fontsize=9.5, color=col, fontweight="bold")

    ax.set_xscale("log")
    ax.set_title(title, fontsize=15, pad=6)
    ax.text(0.035, 0.965, rf"encoding  {counts['encoding']}/16",
            transform=ax.transAxes, ha="left", va="top", fontsize=11.5, color=AMBER,
            fontweight="bold")
    ax.text(0.035, 0.878, rf"inference-specific  {counts['inference-specific']}/16",
            transform=ax.transAxes, ha="left", va="top", fontsize=11.5, color=TEAL,
            fontweight="bold")
    ax.set_xlim(0.9, 220)
    ax.set_xticks([1, 10, 100])
    ax.set_xticklabels(["1", "10", "100"])
    ax.tick_params(labelsize=10)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def main():
    use_lrg_style()
    g = pd.read_csv(CSV)
    rel_df = pd.read_csv(REL)
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.5), sharey=True)
    for ax, (band, title) in zip(axes, BANDS):
        panel(ax, g, rel_df, band, title)
    axes[0].set_ylabel("concordance with persistence\n(fraction of split-half reliability)", fontsize=12)
    fig.text(0.5, 0.085, r"diffusion scale  $s=\tau\lambda_{\max}$   (fine $\rightarrow$ coarse)",
             ha="center", va="bottom", fontsize=12)

    handles = [
        Line2D([0], [0], color=AMBER, lw=2.6, marker="o", mfc=AMBER, mec="white", ms=9,
               label="encoding"),
        Line2D([0], [0], color=TEAL, lw=2.6, marker="o", mfc=TEAL, mec="white", ms=9,
               label="inference-specific"),
        Line2D([0], [0], color="0.4", lw=0, marker="o", mfc="0.4", mec="white", ms=9,
               label=r"significant ($p<0.05$, matched-strength)"),
        Line2D([0], [0], color="0.4", lw=0, marker="o", mfc="white", mec="0.4", mew=1.9, ms=9,
               label="ns"),
        Line2D([0], [0], color=NULLC, lw=6, alpha=0.4, label="matched-strength null"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=5, frameon=False, fontsize=9.6, handletextpad=0.5, columnspacing=1.3)
    fig.subplots_adjust(left=0.09, right=0.985, top=0.90, bottom=0.17, wspace=0.06)

    PDF.parent.mkdir(parents=True, exist_ok=True)
    PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF, transparent=True, bbox_inches="tight")
    fig.savefig(PNG, transparent=True, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # console: the numbers, for the record (raw rho + fraction-of-reliability peak)
    for band, _ in BANDS:
        for func, lab, _ in SERIES:
            x = g[(g.band == band) & (g.functional == func)].sort_values("s")
            nfire = int((x.gate_p < 0.05).sum())
            s1p = float(x.iloc[0].gate_p)
            reln = _rel_on(rel_df, band, x.s.values)
            raw_pk = float(x.obs_median.max())
            frac_pk = float((x.obs_median.values / reln).max())
            print(f"  {band:5s} {lab:18s}: {nfire}/16 scales  (finest s1 gate_p={s1p:.3f})  "
                  f"peak raw rho={raw_pk:.3f} -> {frac_pk:.2f} of reliability")
    print(f"\nwrote {PDF}\n      {PNG}")


if __name__ == "__main__":
    main()
