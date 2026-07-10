#!/usr/bin/env python3
r"""fig:branch-origin-tree — one patient's retained hierarchy, coloured by WHEN each branch
was written: encoding vs inference (§2, single-patient, dendrogram-grounded).

CORE MESSAGE. The consolidation scalars (encoding e, inference-specific f, persistence p)
are abstract. This grounds them in the actual tree: draw ONE patient's rest_post LRG
dendrogram — the hierarchy the brain KEEPS after reasoning — and colour every merge by the
stage that set its height:
    grey  = merge already at baseline (little task change)              -> nothing written
    AMBER = merge whose height was set mostly during ENCODING (shown premises)
    TEAL  = merge whose height was set mostly during INFERENCE (reasoning)
    branch thickness / vividness  ∝  how strongly that change PERSISTS into rest_post
The retained tree is then legible as a mosaic: amber branches = the premises it encoded,
teal branches = the relations it inferred and still holds offline (the abstraction).

Per merge ν of the rest_post tree (leaf-pair block P(ν) = leaves(left) × leaves(right),
all sharing the merge height in an ultrametric):
    h_X(ν) = mean_{(i,j)∈P(ν)} D_coph^X(i,j)         X ∈ {rest_pre, learn, test, rest_post}
    e = h_learn − h_pre   f = h_test − h_learn   p = h_rest_post − h_pre
    hue     = amber(|f|≪|e|) … teal(|f|≫|e|)    via inf_frac = |f| / (|e|+|f|)
    persist = clip(|p| / pctl75|p|, 0, 1)  ->  grey↔hue blend + linewidth
D_coph^X = cophenet(UPGMA(1/ρ̂(τ=1/λ_max))) for phase X (the LRG cophenetic distance;
combinatorial Laplacian, τ = 1/λ_max — the pipeline of record). Cophenetic-only: no
fcluster/ARI/Jaccard.

HONEST REGISTER. Illustrative single-patient decomposition (like the attractor), NOT a
gated test: full-phase cophenetic distances, descriptive attribution. The cohort claim is
the matched-strength gate T_infspec_pe (β p=0.010, 6/10) reported by the forest (Fig 3b).
Runs for ANY band / patient (default: a β tracer); off-trace bands show a mostly-grey tree,
which is itself the contrast.

Reads : phase FC via audit_63.load_phase_fc (imcoh_abs); LRG tree via audit_65.lrg_linkage.
Writes: data/preprint/figures/_drafts/fig_branch_origin_tree_{patient}_{band}_DRAFT.pdf
"""
from __future__ import annotations

import argparse
import sys

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import cophenet, dendrogram
from scipy.spatial.distance import squareform

from lrg_eegfc.visuals.styles import use_lrg_style

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc          # type: ignore
from audit_65_kc_matched_strength_surrogate import lrg_linkage       # type: ignore

OUTDIR = ROOT / "data/preprint/figures/_drafts"
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
BAND_SYM = {"delta": "δ", "theta": "θ", "alpha": "α", "beta": "β",
            "low_gamma": "γ_low", "high_gamma": "γ_high"}

AMBER = np.array([0.91, 0.55, 0.13])     # written during encoding
TEAL = np.array([0.10, 0.60, 0.67])      # written during inference
GREY = np.array([0.74, 0.76, 0.80])      # baseline / reverted


def coph_square(W):
    """N×N cophenetic-distance matrix of the LRG tree for one phase's FC."""
    Z = lrg_linkage(W)
    return squareform(cophenet(Z)), Z


def build_tree(pat, band):
    """Return rest_post linkage + per-phase cophenetic matrices (aligned leaves)."""
    Ws = {ph: load_phase_fc(pat, ph, band) for ph in PHASES}
    ns = {ph: W.shape[0] for ph, W in Ws.items()}
    if len(set(ns.values())) != 1:
        print(f"  !! {pat} {band}: phase N mismatch {ns} — skipping (cannot align pairs)")
        return None
    coph = {}
    Z_rp = None
    for ph, W in Ws.items():
        C, Z = coph_square(W)
        coph[ph] = C
        if ph == "rest_post":
            Z_rp = Z
    return Z_rp, coph


def attribute(Z, coph):
    """Per internal node: (e, f, p) from leaf-block cophenetic means across phases."""
    N = Z.shape[0] + 1
    members = {i: [i] for i in range(N)}
    efp = {}
    for k in range(N - 1):
        c1, c2 = int(Z[k, 0]), int(Z[k, 1])
        L, R = members[c1], members[c2]
        idx = np.ix_(L, R)
        h = {ph: coph[ph][idx].mean() for ph in PHASES}
        efp[N + k] = (h["task_learn"] - h["rest_pre"],
                      h["task_test"] - h["task_learn"],
                      h["rest_post"] - h["rest_pre"])
        members[N + k] = L + R
    return efp, members


def leaf_origin_weights(Z, N, efp):
    """Per-leaf persistence-weighted encoding/inference magnitudes, aggregated over every
    merge on the leaf's root-path: enc[l]=Σ|e|·|p|, inf[l]=Σ|f|·|p|, tot[l]=Σ|p|. Shared by
    the branch-origin node tint (network) and the e/f ribbon — one source of truth."""
    parent = {}
    for k in range(N - 1):
        parent[int(Z[k, 0])] = N + k
        parent[int(Z[k, 1])] = N + k
    enc = np.zeros(N); inf = np.zeros(N); tot = np.zeros(N)
    for l in range(N):
        nd = l
        while nd in parent:
            nd = parent[nd]
            e, f, p = efp[nd]
            w = abs(p)
            tot[l] += w; enc[l] += abs(e) * w; inf[l] += abs(f) * w
    return enc, inf, tot


def link_color(e, f, p, pscale):
    """Categorical origin (no muddy amber/teal midpoint): the dominant stage picks the hue,
    persistence sets vividness + width. Ambiguous *and* non-persistent -> stays grey."""
    side = AMBER if abs(e) >= abs(f) else TEAL
    pw = np.clip(abs(p) / pscale, 0.0, 1.0)
    col = (1 - pw**0.7) * GREY + pw**0.7 * side
    return col, 0.8 + 2.7 * pw, pw


def draw(pat, band):
    got = build_tree(pat, band)
    if got is None:
        return
    Z, coph = got
    N = Z.shape[0] + 1
    efp, members = attribute(Z, coph)
    pscale = np.percentile([abs(p) for _, _, p in efp.values()], 75) or 1e-9

    dend = dendrogram(Z, no_plot=True)
    xpos = {leaf: i for i, leaf in enumerate(dend["leaves"])}
    heights = np.sort(Z[:, 2])
    base = heights[0] * 0.8

    def cx(c):  # x of a node (leaf or internal), via subtree leaf mean
        return np.mean([xpos[l] for l in members[c]])

    def ch(c):  # bottom height of a child link
        return base if c < N else Z[c - N, 2]

    use_lrg_style()
    fig, ax = plt.subplots(figsize=(8.2, 5.4))

    links = []
    for k in range(N - 1):
        c1, c2 = int(Z[k, 0]), int(Z[k, 1])
        h = Z[k, 2]
        x1, x2 = cx(c1), cx(c2)
        e, f, p = efp[N + k]
        col, lw, pw = link_color(e, f, p, pscale)
        links.append((pw, x1, x2, ch(c1), ch(c2), h, col, lw))

    for pw, x1, x2, h1, h2, h, col, lw in sorted(links, key=lambda t: t[0]):  # vivid on top
        ax.plot([x1, x1], [h1, h], color=col, lw=lw, solid_capstyle="round", zorder=2 + pw)
        ax.plot([x2, x2], [h2, h], color=col, lw=lw, solid_capstyle="round", zorder=2 + pw)
        ax.plot([x1, x2], [h, h], color=col, lw=lw, solid_capstyle="round", zorder=2 + pw)

    ax.set_yscale("log")
    ax.set_ylim(base, heights[-1] * 1.05)
    ax.set_xlim(-1.5, N - 0.5)
    ax.set_xticks([])
    ax.set_ylabel(r"cophenetic height  (community distance $1/\hat\rho$)", fontsize=11)
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", length=0)

    handles = [Line2D([], [], color=AMBER, lw=3, label="written at encoding (premises)"),
               Line2D([], [], color=TEAL, lw=3, label="written at inference (reasoning)"),
               Line2D([], [], color=GREY, lw=3, label="baseline / reverted"),
               Line2D([], [], color=(0.4, 0.42, 0.46), lw=3.4,
                      label="thickness ∝ persistence into rest")]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.13),
              frameon=False, fontsize=9, ncol=2, handlelength=1.6, columnspacing=1.6)
    ax.text(0.006, 0.985, f"{pat} · {BAND_SYM.get(band, band)}", transform=ax.transAxes,
            ha="left", va="top", fontsize=10.5, color="#2b2f36", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.8))

    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / f"fig_branch_origin_tree_{pat}_{band}_DRAFT.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True, pad_inches=0.03)
    qa = "/tmp/claude-1000/-home-giulio-Documents-research-neural-networks-lrgeegfc/c2197e94-c496-4e4e-97cc-89e43564e7da/scratchpad"
    fig.savefig(f"{qa}/branchtree_{pat}_{band}.png", dpi=150, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    frac_inf = np.mean([abs(f) > abs(e) for e, f, _ in efp.values()])
    print(f"wrote {out.name}   inf-dominant merges: {frac_inf:.2f}   pscale={pscale:.3g}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--patients", nargs="+", default=["Pat_06"],
                    help="patient id(s); default a beta tracer (Pat_06)")
    ap.add_argument("--bands", nargs="+", default=["beta"],
                    help="band(s); default beta. any of "
                         "delta theta alpha beta low_gamma high_gamma")
    a = ap.parse_args()
    for pat in a.patients:
        for band in a.bands:
            draw(pat, band)


if __name__ == "__main__":
    main()
