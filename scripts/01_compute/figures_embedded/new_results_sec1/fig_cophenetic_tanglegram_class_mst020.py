#!/usr/bin/env python3
r"""fig:cophenetic-tanglegram-class (mst@0.20 rebuild) — the four cross-phase fates
(TRACE / ANCHOR / RESET / REORGANISED) each at a LARGE clade (~20-30 contacts), one fate per
figure, exemplar drawn from across bands x patients — rebuilt on the SPARSIFIED multiscale
cophenetic geometry (mst@0.20 backbone, single diffusion scale s = tau*lambda_max ~= 5.6).

WHAT CHANGED vs the dense original. The dense generator
(``figures_embedded/fig_cophenetic_tanglegram_class.py``) builds every per-phase cophenetic
matrix on the FULLY-CONNECTED graph at ``tau = 1/lambda_max`` (via ``coph_square`` ->
``lrg_linkage``). This fork swaps ONLY that single-scale tree/cophenetic step for the mst@0.20
backbone at the report scale ``s = S_REPORT`` (the honest joint alpha/beta scale from the
sparsified-arc recovery), using the validated shared layer ``new_results_sec1/_common.py``:

    dense:  coph_square(W)                      == squareform(cophenet(lrg_linkage(W))), Z
    here :  _coph_and_linkage(W, s=S_REPORT)    == squareform(cophenet(Z)), Z
            with  Z = linkage_at_scale(*eig_backbone(W), s)   # eig_backbone = mst@0.20 Laplacian

Everything downstream (clade enumeration, the a/b/c subset-rho fate scores, Bezier-ribbon
tanglegram, optimal_leaf_ordering, greedy untangling, dendrogram rows, contact-label beads,
honesty footnotes) is REUSED UNCHANGED from ``fig_cophenetic_tanglegram`` — it operates purely
on the cophenetic matrices, so the whole taxonomy machinery just sees the new graph. Because the
scan re-scores exemplars ON the mst@0.20 cophenetics, the chosen (patient, band, clade) will
generally DIFFER from the dense picks — the exemplar is selected on the SAME graph the figure
draws.

SELECTION (disclosed, seeded). For each fate, scan every (patient, band); take that cell's best
clade of size in [size-min, size-max] by the fate score (fig_cophenetic_tanglegram.scores); pick
the exemplar RANDOMLY among the top-k highest-scoring candidates (seeded, reproducible). The
chosen (patient . band . N . a,b,c) is printed and annotated on the figure.

HONEST REGISTER (inherits the dense figure's preamble). Illustrative, single-clade, subset
 rho^coph — NOT a gated test; trees are untangled by branch rotation so crossings are genuine.
The cohort claim stays the matched-strength gates (mst@0.20 recovery: beta scale-broad, alpha
mesoscale). A big clade tends toward higher rho (more pairs, more inertia), so a clean
"reorganised" (all-low) exemplar is the hardest to source — its printed rho values are reported
honestly whatever they are.

Output ONLY under data/preprint/figures/new_results_sec1/_drafts/ :
    fig_cophenetic_tanglegram_{fate}_{pat}_{band}_N{m}_mst020_DRAFT.pdf   (PDF only, fully vector)
"""
from __future__ import annotations

import argparse
import sys
import time

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import squareform

# shared mst@0.20 + s=5.6 backbone layer (this script lives in new_results_sec1/, so
# `import _common` resolves via the script-dir entry Python puts on sys.path[0]).
import _common as CM
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.utils.fc.heat_multiscale import linkage_at_scale
from lrg_eegfc.visuals.styles import use_lrg_style

# reuse the ENTIRE dense tanglegram machinery (DRY): only load_cophs + output are forked.
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
from fig_cophenetic_tanglegram import (                                # type: ignore  (DRY)
    best_clades, build_row, contact_names, legend_footnotes, n_groups,
    ROW_META, FATES, COHORT, UNIT_IN, SPAN_TOP, BOTTOM_IN, FIG_W,
)
from fig_branch_origin_tree import BAND_SYM                            # type: ignore

BANDS = CM.BANDS                                   # the six canonical bands
OUTDIR = CM.DRAFTDIR                               # data/preprint/figures/new_results_sec1/_drafts
S = CM.S_REPORT                                    # single report scale s = tau*lambda_max ~= 5.646


# ───────────────── mst@0.20 cophenetic geometry (drop-in for the dense load_cophs) ─────────────────
def _load_full_fc(pat, band, phase):
    """Full-phase imcoh_abs FC, normalised exactly as _common.load_phase.

    task_test / rest_post come straight from the canonical sparsified loader; rest_pre needs the
    FULL FC (load_phase only exposes the A/B rest_pre half-splits) via the unified loader, then the
    identical clean/clip/symmetrise so all three phases enter the backbone on the same footing.
    """
    if phase in ("task_test", "rest_post"):
        return CM.load_phase(pat, phase, band)
    W = np.asarray(load_fc_matrix(pat, phase, band, fc_method="imcoh_abs"), float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def _coph_and_linkage(W, s=S):
    """(N×N cophenetic matrix, linkage Z) of the mst@0.20 diffusion tree at scale ``s``.

    Drop-in for the dense ``coph_square(W)`` (fully-connected, tau=1/lambda_max): here the tree is
    the mst@0.20 backbone addressed by the dimensionless scale ``s``. One eigendecomposition, one
    linkage, cophenet derived from THAT same tree (provably consistent).
    """
    ev, V = CM.eig_backbone(W)                     # mst_union_top_fraction(W, 0.20) -> Laplacian eig
    Z = linkage_at_scale(ev, V, s)
    return squareform(cophenet(Z)), Z


def load_cophs(pat, band, task_phase):
    """Per-phase N×N cophenetic matrices (aligned leaves) + rest_post linkage on mst@0.20 at ``S``.

    Mirrors fig_cophenetic_tanglegram.load_cophs exactly, only the FC->cophenetic step is the
    sparsified multiscale one. None on phase-N mismatch (e.g. Pat_10 task-row drop), same as dense.
    """
    phases = ["rest_pre", task_phase, "rest_post"]
    C, Zpost, ns = {}, None, []
    for ph in phases:
        W = _load_full_fc(pat, band, ph)
        Cq, Z = _coph_and_linkage(W, S)
        C[ph] = Cq
        ns.append(W.shape[0])
        if ph == "rest_post":
            Zpost = Z
    if len(set(ns)) != 1:
        return None
    return C, Zpost


# ───────────────────────── fleet scan across patients × bands ─────────────────────────
def scan(patients, bands, task_phase, lo, hi):
    """Return {fate: [(score, pat, band, S, (a,b,c)), ...]} over every (patient, band)."""
    cand = {f: [] for f in FATES}
    combos = [(p, b) for p in patients for b in bands]
    t0 = time.perf_counter()
    for i, (pat, band) in enumerate(combos, 1):
        try:
            got = load_cophs(pat, band, task_phase)
        except Exception as e:                                        # noqa: BLE001
            print(f"  [{i}/{len(combos)}] {pat} {band}: skip ({e})", flush=True)
            continue
        if got is None:
            print(f"  [{i}/{len(combos)}] {pat} {band}: skip (phase N mismatch)", flush=True)
            continue
        best = best_clades(got[0], got[1], task_phase, lo, hi)
        for f in FATES:
            sc, Sc, tri = best[f]
            if Sc is not None and np.isfinite(sc):
                cand[f].append((float(sc), pat, band, Sc, tri))
        if i == 1:
            dt = time.perf_counter() - t0
            print(f"  [1/{len(combos)}] {pat} {band}  ({dt:.1f}s/combo → "
                  f"~{dt * len(combos):.0f}s total)", flush=True)
        else:
            print(f"  [{i}/{len(combos)}] {pat} {band}", flush=True)
    return cand


def pick(cand_list, seed, topk):
    """Random exemplar among the top-k highest-scoring candidates (seeded, reproducible)."""
    ranked = sorted(cand_list, key=lambda t: t[0], reverse=True)[:topk]
    rng = np.random.default_rng(seed)
    return ranked[int(rng.integers(len(ranked)))]


# ───────────────────────── single-fate render ─────────────────────────
def draw_class(fate, pat, band, Sclade, tri, task_phase):
    C = load_cophs(pat, band, task_phase)[0]
    names = contact_names(pat, C["rest_post"].shape[0])
    m = len(Sclade)
    span = m + SPAN_TOP + 1.4                                # extra top room for the class title
    fig_h = span * UNIT_IN + BOTTOM_IN
    use_lrg_style()
    fig = plt.figure(figsize=(FIG_W, fig_h))
    b_frac = BOTTOM_IN / fig_h
    gs = GridSpec(1, 1, top=0.92, bottom=b_frac, figure=fig)
    ax = fig.add_subplot(gs[0])
    # RESET reads clearest with rest_post in the CENTRE: the reverted pair (rest_pre ≈ rest_post)
    # then sits adjacent (gap 1, the similar/parallel one) and the task is the lone scrambled
    # departure on the right — "first similar, second scrambled".
    phase_order = ("rest_pre", "rest_post", task_phase) if fate == "reset" else None
    a, b, c, x1, x2, _ = build_row(ax, C, Sclade, task_phase, is_top=True,
                                   labels=[names[l] for l in Sclade], phase_order=phase_order)
    title, desc = ROW_META[fate]
    if fate == "reset":
        desc = "rest-post shown centre — rest-pre ≈ rest-post (reverted); task is the lone departure"
    fig.text(0.5, 0.982, title, ha="center", va="top", fontsize=15.5, fontweight="bold",
             color="#2b2f36")
    fig.text(0.5, 0.952,
             f"{pat} · {BAND_SYM.get(band, band)}   ·   N = {m} contacts   ·   "
             f"mst@0.20, s={S:.1f}   ·   {desc}",
             ha="center", va="top", fontsize=9, color="#6a7078", style="italic")
    legend_footnotes(fig, b_frac, n_groups(m))

    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / f"fig_cophenetic_tanglegram_{fate}_{pat}_{band}_N{m}_mst020_DRAFT.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True, pad_inches=0.04)   # PDF only, vector
    plt.close(fig)
    print(f"  {fate:11s} → {pat} {band} N={m}  a={a:.2f} b={b:.2f} c={c:.2f}  "
          f"crossings {x1}/{x2}   wrote {out.name}", flush=True)


def main():
    ap = argparse.ArgumentParser(description="single-fate large-N cophenetic tanglegram (mst@0.20)")
    ap.add_argument("--classes", nargs="+", default=list(FATES),
                    help=f"any of {list(FATES)}; default all four")
    ap.add_argument("--bands", nargs="+", default=BANDS)
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--task-phase", default="task_test")
    ap.add_argument("--size-min", type=int, default=20)
    ap.add_argument("--size-max", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0, help="random exemplar seed (per class)")
    ap.add_argument("--topk", type=int, default=6, help="sample among the top-k by score")
    a = ap.parse_args()

    t0 = time.perf_counter()
    print(f"[scan] mst@0.20  s={S:.3f}   {len(a.patients)}×{len(a.bands)} = "
          f"{len(a.patients) * len(a.bands)} (patient,band) combos, clade "
          f"{a.size_min}-{a.size_max}, task={a.task_phase}", flush=True)
    cand = scan(a.patients, a.bands, a.task_phase, a.size_min, a.size_max)

    for fate in a.classes:
        if not cand.get(fate):
            print(f"!! {fate}: no clade of size {a.size_min}-{a.size_max} found", flush=True)
            continue
        sc, pat, band, Sclade, tri = pick(cand[fate], a.seed + FATES.index(fate), a.topk)
        draw_class(fate, pat, band, Sclade, tri, a.task_phase)
    print(f"[done] {time.perf_counter() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
