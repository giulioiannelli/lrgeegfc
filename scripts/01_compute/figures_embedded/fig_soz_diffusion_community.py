#!/usr/bin/env python3
r"""fig — SOZ contacts co-cluster under the diffusion lens (epilepsy thread).

IDEA: lay every contact out by how heat diffuses between them on the |ImCoh| graph — a
force-directed layout of the heat-kernel affinity K = e^{-tau L} (combinatorial Villegas
Laplacian, tau = 1/lambda_max). The layout never sees the seizure-onset labels; we only
COLOUR the SOZ contacts afterwards. In every patient the SOZ contacts settle into a
coherent sub-region of the diffusion plane — i.e. they are more strongly diffusion-connected
to each other than to the rest of the brain: a diffusion community.

WHY THIS IS LEGITIMATE (not a proximity artifact): |ImCoh| is zero-lag-immune (Nolte 2004),
so instantaneous volume conduction / same-probe mixing cannot create the affinity — the
project therefore forbids a spatial-proximity null and requires only the strength-matched
null. Each panel carries its OWN strength-matched cohesion test: is the mean within-SOZ heat
affinity larger than for random contact sets MATCHED on node strength (so "SOZ are just
hubs" is controlled)? p from R strength-matched draws. The cohort-level, leave-one-out
validated version is audit_132 (heat-tau5, matched-strength AUC ~0.80 in delta).

Unsupervised layout; the eigenmap collapses to a blob and cannot show this — the force layout
can. delta band (strongest SOZ readout); rest_pre.

Reads : |ImCoh| FC (imcoh_abs), implant SOZ labels (build_epi_masks).
Writes: data/preprint/figures/results_section1/fig_soz_diffusion_community.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.lines import Line2D
from scipy.linalg import eigh

from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.fc import load_fc_matrix

ROOT = setup_script_env()
use_lrg_style()

OUT = ROOT / "data" / "preprint" / "figures" / "results_section1" / "fig_soz_diffusion_community.pdf"
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BAND, PHASE = "delta", "rest_pre"
KNN = 8                    # layout: keep each node's strongest diffusion links
TAU_FRAC = 1.0             # tau = TAU_FRAC / lambda_max (LRG default)
R_NULL = 2000              # strength-matched cohesion draws
C_SOZ, C_HEALTHY = "#d1352b", "#9aa0a6"


def _clean(W):
    W = np.asarray(W, float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, None)
    return 0.5 * (W + W.T)


def _within_mean(K, idx):
    if idx.size < 2:
        return np.nan
    sub = K[np.ix_(idx, idx)]
    return sub[np.triu_indices(idx.size, 1)].mean()


def _strength_matched_p(K, strength, y, rng, R=R_NULL):
    """p = P(random strength-matched set is as diffusion-cohesive as SOZ). Bins on strength
    quartiles; each draw keeps the SOZ per-bin counts, so 'SOZ are just hubs' is controlled."""
    obs = _within_mean(K, np.where(y)[0])
    edges = np.quantile(strength, [0.25, 0.5, 0.75])
    bins = np.digitize(strength, edges)
    per_bin = {b: (np.where(bins == b)[0], int(((bins == b) & y).sum()))
               for b in np.unique(bins)}
    ge = 1
    for _ in range(R):
        sel = np.zeros(len(y), bool)
        for idx_b, k in per_bin.values():
            if k:
                sel[rng.choice(idx_b, k, replace=False)] = True
        ge += _within_mean(K, np.where(sel)[0]) >= obs
    return obs / _within_mean(K, np.arange(len(y))), (ge / (R + 1))


def heat_affinity_and_layout(pat):
    W = _clean(load_fc_matrix(pat, PHASE, BAND, fc_method="imcoh_abs"))
    y = build_epi_masks(pat).epi_mask
    N = W.shape[0]
    w, V = eigh(np.diag(W.sum(1)) - W)          # combinatorial L
    tau = TAU_FRAC / w[-1]
    K = (V * np.exp(-tau * w)) @ V.T
    K = np.clip(K, 0.0, None)
    np.fill_diagonal(K, 0.0)
    A = np.zeros_like(K)                         # kNN sparsify for the layout only
    for i in range(N):
        j = np.argsort(K[i])[-KNN:]
        A[i, j] = K[i, j]
    A = np.maximum(A, A.T)
    G = nx.from_numpy_array(A)
    pos = nx.spring_layout(G, weight="weight", seed=1,
                           k=1.6 / np.sqrt(N), iterations=250)
    P = np.array([pos[i] for i in range(N)])
    P = (P - P.min(0)) / np.ptp(P, axis=0).clip(1e-9)     # unit square → panels align
    ratio, p = _strength_matched_p(K, W.sum(1), y, np.random.default_rng(7))
    return P, A, y, ratio, p


def main():
    print(f"fig — SOZ diffusion community ({BAND}, {PHASE}); strength-matched cohesion:")
    data = {}
    for pat in COHORT:
        P, A, y, ratio, p = heat_affinity_and_layout(pat)
        data[pat] = (P, A, y, ratio, p)
        print(f"  {pat:8s} nSOZ={int(y.sum()):2d}  within-SOZ affinity ×{ratio:4.1f}  "
              f"strength-matched p={p:.3f}", flush=True)

    fig, axs = plt.subplots(2, 5, figsize=(19.0, 8.0))
    for ax, pat in zip(axs.ravel(), COHORT):
        P, A, y, ratio, p = data[pat]
        G = nx.from_numpy_array(A)
        pos = {i: P[i] for i in range(P.shape[0])}
        nx.draw_networkx_edges(G, pos, ax=ax, edge_color="0.86", width=0.4, alpha=0.5)
        ax.scatter(P[~y, 0], P[~y, 1], s=22, c=C_HEALTHY, edgecolor="white",
                   linewidth=0.4, zorder=3)
        ax.scatter(P[y, 0], P[y, 1], s=62, c=C_SOZ, edgecolor="white",
                   linewidth=0.7, zorder=5)
        pstr = "p<0.001" if p < 0.001 else f"p={p:.3f}"
        ax.set_title(f"{pat}   SOZ ×{ratio:.1f}, {pstr}", fontsize=11.5,
                     color="0.2", pad=4)
        ax.set_xlim(-0.06, 1.06)
        ax.set_ylim(-0.06, 1.06)
        ax.set_box_aspect(1.0)
        ax.axis("off")

    handles = [
        Line2D([0], [0], marker="o", ls="", mfc=C_SOZ, mec="white", ms=11,
               label="seizure-onset (SOZ) contact"),
        Line2D([0], [0], marker="o", ls="", mfc=C_HEALTHY, mec="white", ms=9,
               label="other contact"),
        Line2D([0], [0], color="0.86", lw=1.4, label="strongest heat-diffusion links"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=3, frameon=False, fontsize=12, handletextpad=0.5, columnspacing=2.2)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.06,
                        wspace=0.06, hspace=0.14)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
