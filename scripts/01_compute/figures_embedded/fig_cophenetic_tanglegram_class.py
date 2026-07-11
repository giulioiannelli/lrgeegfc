#!/usr/bin/env python3
r"""fig:cophenetic-tanglegram-class — ONE cross-phase fate at a LARGE clade (≥20 contacts),
one fate per figure, exemplar drawn from across bands × patients (§5.3 explainer, companion to
fig_cophenetic_tanglegram).

WHY A SECOND FIGURE. The combined figure packs all four fates into one patient at small N so the
taxonomy reads at a glance. This one does the opposite: it gives each fate its own big canvas
(≥20 contacts) so the crossing pattern is unmistakable, and it samples the exemplar from the
WHOLE cohort × six bands rather than one brain — evidence that each ρ^coph regime is a general
phenomenon, not a Pat_02 quirk. Four fates → four files:
    trace        rest_pre scrambled → task grouped → rest_post grouped   (written at task, held)
    anchor       grouped in all three phases                            (nothing moved)
    reset        grouped in rest_pre & rest_post, scrambled at task      (moved, then reverted)
    reorganized  scrambled everywhere                                    (changes every transition)

SELECTION. For each fate we scan every (patient, band), take that (patient, band)'s best clade of
size in [size-min, size-max] by the fate score (fig_cophenetic_tanglegram.scores), then pick the
exemplar RANDOMLY among the top-k highest-scoring candidates (seeded) — a strong but not
hand-tuned instance, and re-seeding gives a different valid one. The chosen (patient · band · N ·
a,b,c) is printed and annotated on the figure.

HONEST REGISTER (inherits the combined figure's preamble). Illustrative, single-clade, subset
ρ^coph — NOT a gated test; trees are untangled by branch rotation so crossings are genuine. The
cohort claim stays the matched-strength gates. A big clade tends toward higher ρ (more pairs, more
inertia), so a clean "reorganized" (all-low) exemplar is the hardest to source — its printed ρ
values are reported honestly whatever they are.

Reuses the entire machinery of fig_cophenetic_tanglegram (DRY); only the fleet-scan + single-row
layout are new.
Writes: data/preprint/figures/_drafts/fig_cophenetic_tanglegram_{fate}_{pat}_{band}_N{m}_DRAFT.pdf
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

from lrg_eegfc.visuals.styles import use_lrg_style

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
from fig_cophenetic_tanglegram import (                                # type: ignore  (DRY)
    load_cophs, best_clades, build_row, contact_names, legend_footnotes, n_groups,
    ROW_META, FATES, COHORT, UNIT_IN, SPAN_TOP, BOTTOM_IN, FIG_W, OUTDIR, QA,
)
from fig_branch_origin_tree import BAND_SYM                            # type: ignore

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]


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
            sc, S, tri = best[f]
            if S is not None and np.isfinite(sc):
                cand[f].append((float(sc), pat, band, S, tri))
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
def draw_class(fate, pat, band, S, tri, task_phase):
    C = load_cophs(pat, band, task_phase)[0]
    names = contact_names(pat, C["rest_post"].shape[0])
    m = len(S)
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
    a, b, c, x1, x2, _ = build_row(ax, C, S, task_phase, is_top=True,
                                   labels=[names[l] for l in S], phase_order=phase_order)
    title, desc = ROW_META[fate]
    if fate == "reset":
        desc = "rest-post shown centre — rest-pre ≈ rest-post (reverted); task is the lone departure"
    fig.text(0.5, 0.982, title, ha="center", va="top", fontsize=15.5, fontweight="bold",
             color="#2b2f36")
    fig.text(0.5, 0.952, f"{pat} · {BAND_SYM.get(band, band)}   ·   N = {m} contacts   ·   {desc}",
             ha="center", va="top", fontsize=9, color="#6a7078", style="italic")
    legend_footnotes(fig, b_frac, n_groups(m))

    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / f"fig_cophenetic_tanglegram_{fate}_{pat}_{band}_N{m}_DRAFT.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True, pad_inches=0.04)
    fig.savefig(f"{QA}/tanglegram_class_{fate}.png", dpi=150, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print(f"  {fate:11s} → {pat} {band} N={m}  a={a:.2f} b={b:.2f} c={c:.2f}  "
          f"crossings {x1}/{x2}   wrote {out.name}", flush=True)


def main():
    ap = argparse.ArgumentParser(description="single-fate large-N cophenetic tanglegram")
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
    print(f"[scan] {len(a.patients)}×{len(a.bands)} = {len(a.patients) * len(a.bands)} "
          f"(patient,band) combos, clade {a.size_min}-{a.size_max}, task={a.task_phase}", flush=True)
    cand = scan(a.patients, a.bands, a.task_phase, a.size_min, a.size_max)

    for fate in a.classes:
        if not cand.get(fate):
            print(f"!! {fate}: no clade of size ≥{a.size_min} found", flush=True)
            continue
        sc, pat, band, S, tri = pick(cand[fate], a.seed + FATES.index(fate), a.topk)
        draw_class(fate, pat, band, S, tri, a.task_phase)
    print(f"[done] {time.perf_counter() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
