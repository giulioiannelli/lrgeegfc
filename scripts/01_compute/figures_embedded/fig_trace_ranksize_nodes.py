#!/usr/bin/env python3
"""fig_trace_ranksize_nodes — rank-size (Zipf) view of the per-node trace contribution.

Question (user, 2026-07-10): the cohort trace rho_sym is a cross-phase rank
correlation; decompose it per node (size_i) and rank-size plot size vs rank. Is it
Zipfian -- a few nodes driving the trace, most inert (a power law)?

Size = the per-node ADDITIVE trace contribution T_i = node_incidence_mean(c_sym)
(audit_154), where c_p = (rank(dD_task)-mid)(rank(dD_rest)-mid) and, exactly,
mean_i T_i is proportional to rho_sym. So T_i is each node's literal share of the
cohort trace: total trace = sum of node sizes. Cached: data/audit/
per_node_trace_decomposition_rhosym/per_node.csv (all 6 bands, n=10).

VERDICT (data, not vibes): the Zipf hypothesis is REFUTED for the real signal -- the
truth is the inversion. A strong trace is DELOCALISED: Pat_08/05/03 beta have ~88%
of nodes co-moving and the top-10 nodes carry <20%. Measured SIGN-INDEPENDENTLY
(Gini of |T_i|, so the trend is not a sign tautology), strong-trace beta sits at
Gini 0.32 vs resetters 0.51; cohort Spearman(rho_sym, Gini|T|) = -0.50 (p=4e-5).
[The naive Gini of the positive part gives -0.96 but is inflated -- it is slaved to
sign(mean T)=sign(rho) -- so it is NOT used here.] The Zipf-like CONCENTRATED profile
(a few idiosyncratic nodes dominate) is the signature of trace ABSENCE (resetters
Pat_15/10; null bands theta/delta). Delocalisation IS the trace; a heavy tail is what
noise looks like. Confirms + quantifies the delocalised-trace finding
(dendrogram_persistence_gate) and "beta delocalised single-contact" (heat-brain).

Panels
  A  beta rank-size per patient, LINEAR signed share (T_i / sum_j|T_j|), coloured by
     the patient's rho_sym: strong-trace = gentle positive plateau (delocalised);
     reset = a spike then a negative tail (concentrated).
  B  beta size distribution, LOG-BINNED pdf (log-spaced bins, counts / bin-width / N —
     the correct power-law view; a straight line = a power law). Pooled per group
     (strong-trace vs reset) with slope -1 / -2 guides: neither group is a straight
     line -> not a power law; the reset pool is only heavier-tailed, not Zipfian.
  C  Gini(contribution) vs rho_sym, all 60 patient x band cells, coloured by band:
     the monotone inversion (delocalised <=> trace).
  D  per-band Gini (slow->fast): beta/alpha least concentrated (real trace), theta/
     delta most concentrated (null).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

SRC = ROOT / "data/audit/per_node_trace_decomposition_rhosym/per_node.csv"
OUT = FIGURES_ROOT / "trace_ranksize_nodes.pdf"
BANDS = list(BRAIN_BANDS)


def gini_pos(x):
    """Gini of the non-negative part (concentration; 0 = uniform, ->1 = one node)."""
    x = np.sort(np.clip(np.asarray(x, float), 0, None))
    n = x.size
    if n == 0 or x.sum() == 0:
        return np.nan
    return (2 * np.sum(np.arange(1, n + 1) * x) / (n * x.sum())) - (n + 1) / n


def pr_over_n(x):
    """Effective participating fraction: (sum|x|)^2 / (N * sum x^2). 1 = uniform."""
    x = np.abs(np.asarray(x, float))
    s2 = np.sum(x ** 2)
    return (x.sum() ** 2 / s2) / x.size if s2 > 0 else np.nan


def logbin_pdf(x, nbins=13):
    """Log-binned probability density: log-spaced bins, counts / bin-width / N.
    Dividing by the (linear) bin width is what makes a power law a straight line;
    without it the heavy tail is artificially flattened. Geometric-mean centers."""
    x = np.asarray(x, float)
    x = x[x > 0]
    edges = np.logspace(np.log10(x.min()), np.log10(x.max()), nbins + 1)
    cnt, _ = np.histogram(x, bins=edges)
    width = np.diff(edges)
    ctr = np.sqrt(edges[:-1] * edges[1:])
    dens = cnt / width / x.size
    m = cnt > 0
    return ctr[m], dens[m], cnt[m]


def cell_stats(g):
    T = g.T_raw.to_numpy(float)
    rho = float(g.rho_sym.iloc[0])
    absum = np.sum(np.abs(T))
    share = T / absum if absum > 0 else T
    ssort = np.sort(share)[::-1]                     # signed share, descending (panel A)
    snorm = np.abs(T) / np.mean(np.abs(T))           # size / <|T|>, per-patient comparable (panel B)
    # concentration is measured SIGN-INDEPENDENTLY (Gini of |T|): the naive Gini of the
    # positive part is slaved to sign(mean T)=sign(rho) and gives a tautological -0.96.
    return dict(rho=rho, N=len(T), ssort=ssort, snorm=snorm,
                gini=gini_pos(np.abs(T)), prn=pr_over_n(T),
                frac_pos=float((T > 0).mean()))


def main():
    use_lrg_style()
    df = pd.read_csv(SRC)
    cells = {(p, b): cell_stats(g) for (p, b), g in df.groupby(["patient", "band"])}

    fig = plt.figure(figsize=(11.0, 8.6))
    gs = GridSpec(2, 2, figure=fig, hspace=0.30, wspace=0.26,
                  left=0.075, right=0.90, top=0.95, bottom=0.09)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axC, axD = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    # rho_sym -> colour (diverging; trace positive = blue, reset = red)
    rhos = [cells[(p, "beta")]["rho"] for p in df.patient.unique()
            if (p, "beta") in cells]
    norm = TwoSlopeNorm(vmin=min(rhos) - 1e-3, vcenter=0.0, vmax=max(rhos) + 1e-3)
    cmap = plt.get_cmap("RdBu")

    # ---- A: beta rank-size, linear signed share, coloured by rho_sym ----
    beta_cells = sorted(((p, cells[(p, "beta")]) for p in df.patient.unique()
                         if (p, "beta") in cells), key=lambda t: t[1]["rho"])
    for p, c in beta_cells:
        col = cmap(norm(c["rho"]))
        axA.plot(np.arange(1, c["N"] + 1), c["ssort"], color=col, lw=1.4,
                 alpha=0.9, solid_capstyle="round")
    axA.axhline(0.0, color="0.5", lw=0.7, ls=(0, (4, 4)))
    axA.set_xlabel("node rank")
    axA.set_ylabel(r"signed share  $T_i / \sum_j |T_j|$")
    axA.set_title(r"$\beta$  per-node contribution (linear)", fontsize=10)
    axA.text(0.97, 0.95, "strong trace: flat +plateau\n(delocalised)",
             transform=axA.transAxes, ha="right", va="top", fontsize=7.5, color="#2166ac")
    axA.text(0.97, 0.06, "reset: spike + negative tail\n(concentrated)",
             transform=axA.transAxes, ha="right", va="bottom", fontsize=7.5, color="#b2182b")
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=axA, pad=0.02, fraction=0.05)
    cb.set_label(r"$\rho_{\mathrm{sym}}$ (trace)", fontsize=8)

    # ---- B: beta size DISTRIBUTION, LOG-BINNED pdf (proper power-law view) ----
    pats_b = [p for p in df.patient.unique() if (p, "beta") in cells]
    strong = [p for p in pats_b if cells[(p, "beta")]["rho"] > 0.3]
    reset = [p for p in pats_b if cells[(p, "beta")]["rho"] < 0.0]
    slopes = {}
    for grp, col, lab in [(strong, "#2166ac", f"strong trace, n={len(strong)} (delocalised)"),
                          (reset, "#b2182b", f"reset, n={len(reset)} (concentrated)")]:
        if not grp:
            continue
        pooled = np.concatenate([cells[(p, "beta")]["snorm"] for p in grp])
        x, dens, cnt = logbin_pdf(pooled, nbins=13)
        axB.loglog(x, dens, "o-", color=col, ms=4.5, lw=1.4, label=lab, zorder=3)
        # descriptive OLS slope over populated log-bins (honest caveat: OLS-on-log-bins
        # is biased vs MLE; used only to show it is not a single straight power law)
        s = np.polyfit(np.log10(x), np.log10(dens), 1)[0]
        slopes[lab.split(",")[0]] = s
    # power-law slope guides (a straight line here = a power law)
    xg = np.array([0.12, 3.5])
    for a, y0 in [(-1.0, 0.6), (-2.0, 0.6)]:
        yg = y0 * (xg / xg[0]) ** a
        axB.plot(xg, yg, color="0.6", lw=0.8, ls="--", zorder=1)
        axB.text(xg[1] * 1.03, yg[1], f"$-{abs(a):.0f}$", fontsize=7, color="0.5", va="center")
    axB.set_xlabel(r"normalized size  $|T_i| / \langle |T| \rangle$")
    axB.set_ylabel("log-binned density")
    axB.set_title(r"$\beta$  size distribution (log-binned)", fontsize=10)
    axB.legend(fontsize=6.8, frameon=False, loc="lower left")

    # ---- C: Gini vs rho_sym, all cells, coloured by band ----
    for b in BANDS:
        xs = [cells[(p, b)]["rho"] for p in df.patient.unique() if (p, b) in cells]
        ys = [cells[(p, b)]["gini"] for p in df.patient.unique() if (p, b) in cells]
        axC.scatter(xs, ys, s=34, color=band_color(b), edgecolor="white",
                    linewidth=0.4, label=BRAIN_BAND_TEX_DICT.get(b, b), zorder=3)
    allr = [cells[k]["rho"] for k in cells]
    allg = [cells[k]["gini"] for k in cells]
    rho_s, p_s = spearmanr(allr, allg)
    axC.set_xlabel(r"$\rho_{\mathrm{sym}}$ (trace strength)")
    axC.set_ylabel(r"Gini of $|T_i|$ (sign-independent)")
    axC.set_title(f"concentration falls as trace grows  "
                  f"($\\rho_S={rho_s:.2f}$)", fontsize=10)
    axC.text(0.03, 0.05, "delocalised", transform=axC.transAxes, fontsize=7.5,
             color="0.3", va="bottom")
    axC.text(0.03, 0.95, "concentrated (Zipf-like)", transform=axC.transAxes,
             fontsize=7.5, color="0.3", va="top")

    # ---- D: per-band Gini, slow->fast ----
    for i, b in enumerate(BANDS):
        gv = np.array([cells[(p, b)]["gini"] for p in df.patient.unique()
                       if (p, b) in cells], float)
        gv = gv[np.isfinite(gv)]
        x = np.full(gv.size, i) + np.linspace(-0.16, 0.16, gv.size)
        axD.scatter(x, gv, s=26, color=band_color(b), edgecolor="white", linewidth=0.3, zorder=3)
        axD.plot([i - 0.28, i + 0.28], [np.median(gv)] * 2, color=band_color(b, 0.7), lw=2.4)
    axD.set_xticks(range(len(BANDS)))
    axD.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS])
    axD.set_ylabel(r"Gini of $|T_i|$")
    axD.set_title("per-band concentration (bar = median)", fontsize=10)

    handles, labels = axC.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.01),
               ncol=len(labels), frameon=False, fontsize=8)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight")
    plt.close(fig)

    # ---- cohort numbers for the writeup ----
    print(f"[fig] {OUT}")
    st = pd.DataFrame([dict(patient=p, band=b, rho=c["rho"], N=c["N"],
                            gini=c["gini"], prn=c["prn"], frac_pos=c["frac_pos"])
                       for (p, b), c in cells.items()])
    print("\nper-band medians (Gini|T|, PR/N, frac co-moving, rho_sym):")
    print(st.groupby("band").agg(rho=("rho", "median"), gini=("gini", "median"),
                                 prn=("prn", "median"),
                                 frac_pos=("frac_pos", "median")).round(3).to_string())
    strong = st[(st.band == "beta") & (st.rho > 0.3)]
    reset = st[(st.band == "beta") & (st.rho < 0)]
    print(f"\nbeta strong-trace (rho>0.3, n={len(strong)}): Gini|T| {strong.gini.median():.2f}, "
          f"PR/N {strong.prn.median():.2f}, frac co-moving {strong.frac_pos.median():.2f}")
    print(f"beta reset (rho<0, n={len(reset)}): Gini|T| {reset.gini.median():.2f}, "
          f"PR/N {reset.prn.median():.2f}, frac co-moving {reset.frac_pos.median():.2f}")
    print(f"Spearman(rho_sym, Gini|T|) all cells = {rho_s:.2f} (p={p_s:.1e}) "
          f"-- negative = delocalisation grows with trace (sign-independent metric)")
    print("\nlog-binned beta size-distribution OLS slopes (descriptive; a clean power "
          "law would be a single straight line ~ -1..-3):")
    for k, v in slopes.items():
        print(f"  {k:20s} slope {v:+.2f}")


if __name__ == "__main__":
    main()
