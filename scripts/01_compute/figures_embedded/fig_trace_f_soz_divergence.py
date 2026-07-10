#!/usr/bin/env python3
r"""fig:trace_f — over the seizure-onset zone, alpha and beta diverge (Results §1, para f).

CORE MESSAGE: the two carrier bands treat the diseased core oppositely, and the object
that diverges is the LRG cophenetic HIERARCHY itself — not an FC edge.

TOP (exemplar geometry, Pat_14). For every contact PAIR, the cophenetic distance is
tracked as it moves rest_pre -> task (x) and rest_pre -> rest_post (y), rank-scored
within its class. A pair on the green diagonal HELD its task-induced hierarchy move into
rest (concordant = trace). SOZ<->SOZ pairs (crimson) over the faint non-SOZ sea (grey):
in ALPHA the SOZ pairs collapse onto the diagonal (rho_sym = +0.43) — the seizure-onset
sub-network shares one held hierarchy move; in BETA the same pairs scatter (rho_sym =
-0.39) — beta does not hold the SOZ hierarchy.

BOTTOM (cohort forest, n=10, the load-bearing claim). Observed rho_sym per pair-class vs
its matched-strength surrogate null (open marker). ALPHA recruits the SOZ: SOZ<->SOZ
+0.41 sits far above a NEGATIVE null (weak contacts, yet strongest held trace; p=0.005,
pair-null p=0.000 concentrated). BETA is the mirror: SOZ<->SOZ +0.08 on its null
(p=0.31), the trace living instead in cortex<->cortex (+0.17). alpha pulls the diseased
core in; beta steers around it.

HONEST FRAMING. The per-pair co-movement MATRIX is node-polarity-streaked and shows no
SOZ block (that is why a clustermap fails); the effect is a WITHIN-CLASS rank statistic,
so the load-bearing panel is the forest with the matched-strength null. Pat_14 is one
exemplar (alpha SOZ +0.43 ~ cohort +0.41; its beta SOZ is a reset, cohort beta is
generic) chosen for SOZ coverage + the clearest alpha-concordant / beta-scattered split.

Reads : imcoh_abs FC + LRG cophenetic (task/post + rest_pre halves), Pat_14 alpha/beta;
        SOZ membership from per_node_trace_decomposition_rhosym/per_node.csv (is_epi);
        cohort forest from epi_wm_stratified_rhosym/pairclass_rhosym_cohort.csv.
Writes: data/preprint/figures/results_section1/fig_trace_f_soz_divergence.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from scipy.stats import rankdata, spearmanr

from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.metrics.node_localization import canonical_cophenet
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style, band_color
from lrg_eegfc.workflow.fc import load_fc_matrix

ROOT = setup_script_env()
use_lrg_style()

OUT = ROOT / "data/preprint/figures/results_section1/fig_trace_f_soz_divergence.pdf"
PN = ROOT / "data/audit/per_node_trace_decomposition_rhosym/per_node.csv"
FOREST = ROOT / "data/audit/epi_wm_stratified_rhosym/pairclass_rhosym_cohort.csv"
HALVES = CACHE_ROOT / "imcoh_halves_fc"

PATIENT, FC = "Pat_14", "imcoh_abs"
BANDS = ["alpha", "beta"]
C_SOZ = (0.80, 0.12, 0.14)                       # crimson seizure-onset contacts
C_SEA = (0.72, 0.72, 0.72)                       # faint non-SOZ sea
C_DIAG = (0.10, 0.49, 0.24)                      # green "held" diagonal
CLASSES = [("epi_epi", "SOZ – SOZ"), ("cross", "SOZ – cortex"),
           ("nonepi_nonepi", "cortex – cortex")]


def _rankz(x):
    return 2.0 * (rankdata(x) - 1.0) / (x.size - 1) - 1.0


def _cov_ellipse(ax, x, y, color, nstd=1.0):
    """Filled 1-sigma covariance ellipse = the SOZ pairs' CORE extent: a thin sliver
    along the alpha held-diagonal (concordant) vs a broad, near-round blob in beta."""
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    w, h = 2.0 * nstd * np.sqrt(np.maximum(vals, 1e-9))
    ax.add_patch(Ellipse((float(np.mean(x)), float(np.mean(y))), w, h, angle=ang,
                         facecolor=color, alpha=0.18, edgecolor=color, lw=1.4,
                         zorder=1.5, clip_on=True))


def _cophenet_moves(band):
    """Averaged-pre cophenetic shifts (rest->task, rest->post) + symmetric SOZ rho_sym."""
    Wt = load_fc_matrix(PATIENT, "task_test", band, FC)
    Wp = load_fc_matrix(PATIENT, "rest_post", band, FC)
    WA = np.load(HALVES / PATIENT / f"{band}_rest_pre_A_{FC}.npy")
    WB = np.load(HALVES / PATIENT / f"{band}_rest_pre_B_{FC}.npy")
    Dt, Dp, DA, DB = (canonical_cophenet(W) for W in (Wt, Wp, WA, WB))
    pre = 0.5 * (DA + DB)
    return Dt - pre, Dp - pre, Dt, Dp, DA, DB


def _soz_pair_mask(band):
    pn = pd.read_csv(PN)
    soz = (pn[(pn.patient == PATIENT) & (pn.band == "beta")]
           .sort_values("node_idx").is_epi.to_numpy(bool))
    N = len(soz)
    iu = np.triu_indices(N, 1)
    return soz[iu[0]] & soz[iu[1]]


def scatter_panel(ax, band, tag=None,
                  xlabel=r"hierarchy shift  rest$\rightarrow$task  (rank)",
                  ylabel=r"hierarchy shift  rest$\rightarrow$post  (rank)"):
    dtask, dpost, Dt, Dp, DA, DB = _cophenet_moves(band)
    soz_pair = _soz_pair_mask(band)
    other = ~soz_pair
    # canonical SYMMETRIC rho_sym restricted to SOZ pairs (matches the reported number)
    r1 = spearmanr((Dt - DA)[soz_pair], (Dp - DB)[soz_pair]).statistic
    r2 = spearmanr((Dt - DB)[soz_pair], (Dp - DA)[soz_pair]).statistic
    rho_soz = 0.5 * (r1 + r2)

    xs_soz, ys_soz = _rankz(dtask[soz_pair]), _rankz(dpost[soz_pair])
    ax.axhline(0, color="0.75", lw=0.6, zorder=0)
    ax.axvline(0, color="0.75", lw=0.6, zorder=0)
    ax.plot([-1, 1], [-1, 1], ls="--", lw=1.4, color=C_DIAG, zorder=1)
    ax.scatter(_rankz(dtask[other]), _rankz(dpost[other]), s=6, color=C_SEA,
               alpha=0.28, linewidths=0, zorder=2, rasterized=False)
    _cov_ellipse(ax, xs_soz, ys_soz, C_SOZ)                 # SOZ extent: thin in alpha, broad in beta
    ax.scatter(xs_soz, ys_soz, s=44, color=C_SOZ,
               edgecolor="white", linewidths=0.4, zorder=3)
    ax.set_xlim(-1.04, 1.04); ax.set_ylim(-1.04, 1.04)
    ax.set_aspect("equal")
    ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    col = band_color(band, shade=0.85)
    sym = r"$\alpha$" if band == "alpha" else r"$\beta$"
    # band symbol boxed at the top-left corner so it reads clearly above the point cloud
    # without colliding with the legend strip above; the SOZ rho_sym value is caption-only
    ax.text(0.04, 0.96, sym, transform=ax.transAxes, fontsize=23, va="top",
            ha="left", color=col, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none", alpha=0.72))
    if tag:
        ax.text(-0.20, 1.02, rf"$\mathbf{{{tag}}}$", transform=ax.transAxes,
                fontsize=18, va="bottom", ha="left", fontweight="bold")
    return rho_soz


def forest_panel(ax, sym_x=0.03):
    d = pd.read_csv(FOREST)
    d = d[d.stratify == "epi"]
    rows, y = [], 0
    yticks, ylabels = [], []
    band_bands = []
    for band in ["alpha", "beta"]:
        col = band_color(band)
        for cfg, lab in CLASSES:                         # top row = SOZ-SOZ
            r = d[(d.band == band) & (d.config == cfg)].iloc[0]
            rows.append((y, band, col, cfg, float(r.cohort_obs_rho_sym),
                         float(r.ms_surr_median), float(r.ms_wilcoxon_p)))
            yticks.append(y); ylabels.append(lab)
            y -= 1
        band_bands.append((band, col, y + 1))
        y -= 0.6                                          # gap between bands

    ax.axvline(0, color="0.6", lw=0.9, zorder=1)
    for yy, band, col, cfg, obs, null, p in rows:
        soz = cfg == "epi_epi"
        ax.plot([null, obs], [yy, yy], color=col, lw=2.6 if soz else 1.6,
                alpha=0.9 if soz else 0.55, zorder=2, solid_capstyle="round")
        ax.scatter([null], [yy], s=46, facecolor="white", edgecolor="0.45",
                   linewidths=1.3, zorder=3)
        ax.scatter([obs], [yy], s=150 if soz else 78, color=col,
                   edgecolor="white", linewidths=0.8, zorder=4)
        if p < 0.05:
            ax.text(obs + 0.05, yy, "$\\ast$", ha="left", va="center",
                    fontsize=15, color=col)
    ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=11)
    for band, col, ytop in band_bands:
        sym = r"$\alpha$" if band == "alpha" else r"$\beta$"
        ax.text(sym_x, ytop + 1.0, sym, transform=ax.get_yaxis_transform(),
                fontsize=22, va="center", ha="left", color=band_color(band, 0.85),
                fontweight="bold", clip_on=False)
    ax.set_ylim(y + 0.4, 0.9)
    ax.set_xlim(-0.09, 0.53)
    ax.set_xlabel(r"cophenetic co-movement  $\rho_{\mathrm{sym}}$  (cohort, $n=10$)")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    handles = [
        Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.45", ms=9,
               label="matched-strength null"),
        Line2D([0], [0], marker="o", ls="", mfc="0.35", mec="white", ms=10,
               label="observed"),
        Line2D([0], [0], marker="$\\ast$", ls="", color="0.35", ms=11,
               label=r"$p<0.05$ vs null"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=9.5,
              handletextpad=0.3, borderaxespad=0.4)


def draw(target):
    gs = target.add_gridspec(2, 2, height_ratios=[1.0, 0.58], hspace=0.13, wspace=0.24,
                             left=0.10, right=0.975, top=0.90, bottom=0.07)
    ax_a = target.add_subplot(gs[0, 0]); ax_b = target.add_subplot(gs[0, 1])
    ax_f = target.add_subplot(gs[1, :])

    ra = scatter_panel(ax_a, "alpha", tag="c")
    rb = scatter_panel(ax_b, "beta")

    handles = [
        Line2D([0], [0], marker="o", ls="", mfc=C_SOZ, mec="white", ms=10,
               label="SOZ–SOZ pair"),
        Line2D([0], [0], marker="s", ls="", mfc=C_SOZ, mec=C_SOZ, ms=11, alpha=0.32,
               label="SOZ extent (1σ)"),
        Line2D([0], [0], marker="o", ls="", mfc=C_SEA, mec="none", ms=9,
               label="non-SOZ pair"),
        Line2D([0], [0], ls="--", color=C_DIAG, lw=2, label="held (concordant) diagonal"),
    ]
    target.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.975),
                  ncol=4, frameon=False, fontsize=10, handletextpad=0.5, columnspacing=1.5)

    forest_panel(ax_f)


def main():
    fig = plt.figure(figsize=(9.6, 8.1))
    draw(fig)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
