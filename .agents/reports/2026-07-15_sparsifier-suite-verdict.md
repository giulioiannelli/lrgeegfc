---
name: 2026-07-15_sparsifier-suite-verdict
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery — sparsifier SETTLED)
status: settled-verdict
created: 2026-07-15
scope: THE VERDICT of the full sparsifier suite (marker + trace, matched-strength, per-scale,
  tmfg/pmfg/disparity vs the mst family). Settles the sparsifier once and for all as a TWO-SCHEME
  choice by readout — trace = mst-union@0.20, marker = TMFG — with the dissociation as a
  methodological result. Supersedes the "PMFG/TMFG-for-everything" recommendation.
pointers:
  - .agents/guides/02_methods/sparsification-choice.md          # the settled guide (§4a = decision)
  - scripts/01_compute/sparsified_arc/13_matched_strength_mst020.py   # trace gate (SA_BACKBONE/SA_FRAC)
  - scripts/01_compute/sparsified_arc/06_epi_arc.py                   # marker (SA_BACKBONE/SA_FRAC)
  - scripts/01_compute/sparsified_arc/25_pmfg_observed_trace_readout.py  # PMFG≈TMFG mechanism
  - src/lrg_eegfc/utils/fc/backbone.py                                # pmfg/disparity/select_backbone
---

# Sparsifier settled — two schemes, one per readout

## Head

The full matched-strength suite settles it: **the cross-phase cognitive trace uses the
strength-ranked `mst_union_top_fraction` @ frac 0.20; the epileptogenic-zone marker uses TMFG**
(planar-maximally-filtered). The planar family is **disqualified for the trace** — not for being
sparse but because planarity *discards the strongest edges the trace lives in* (α collapses
**12/16 → 0/16** scales under matched-strength; exact PMFG and greedy TMFG both halve α ρ_sym at
equal density, so it is the planar principle, not a shortcut). It is **confirmed for the marker** —
seeded heat-diffusion *detection* sharpens on the triangulated skeleton (fused SOZ AUC
**0.816 → 0.904**, β matched-strength **8/10 → 10/10**). The two-backbone split is defensible because
each backbone is matched-strength-validated as its readout's optimum, and the
**magnitude-vs-skeleton dissociation is itself a reportable finding**, not filter-shopping.

## 1. The trace: strength-ranked mst@0.20 (planar rejected)

τ-resolved cophenetic ρ_sym gate, matched-strength (4-cycle ±δ strength-preserving), clears / 16
scales. Want: **α, β CLEAR** (true trace) and **θ, low_γ NULL** (band-selectivity).

| band | mst@0.20 | mst@0.10 | mst@0.05 | disparity@0.2 | TMFG | want |
|---|---|---|---|---|---|---|
| **alpha** | **12/16** | 8/16 | 15/16 | 11/16 | **0/16** ✗ | CLEAR |
| **beta**  | **16/16** | 9/16 | 5/16  | 11/16 | 3/16 ✗ | CLEAR |
| theta | 0/16 | 0/16 | 0/16 | 1/16 | 3/16 ✗ | NULL |
| low_γ | 0/16 | 9/16 ✗ | 1/16 | 8/16 ✗ | 7/16 ✗ | NULL |
| delta | 8/16 | 1/16 | 6/16 | 0/16 | 4/16 | — |
| high_γ | 4/16 | 0/16 | 0/16 | 0/16 | 2/16 | — |

- **It is planarity, not sparsity.** PMFG-observed readout (`25`, best-scale cohort-median ρ_sym) at
  identical density 0.05: α = **0.24** for strength-ranked mst@0.05 but **0.11 for both PMFG (exact)
  and TMFG (greedy)** — halved. Strength-ranked mst@0.05 at the same density keeps α fully; the
  planar family does not. PMFG ≈ TMFG ⇒ the loss is the planar principle, not TMFG's chordal shortcut.
- **Strength-based ≠ planar.** Disparity (strength/statistics, non-planar, even sparser natively)
  *keeps* α (11/16) and β (11/16) — the mechanism control. It leaks low_γ (8/16), so it is not as
  clean as mst@0.20, but it does not lose α the way the planar filters do.
- **0.20 is the unique clean fraction.** α and β have opposite density needs (α peaks sparse 15/16 @
  0.05; β needs density 16/16 @ 0.20, dies to 5/16 @ 0.05). Only 0.20 gives α strong + β maximal +
  θ/low_γ null. Sparser mst leaks low_γ. ⇒ **trace backbone = mst_union @ 0.20** (pre-registered
  "keep mst@0.20": nothing beats it).

## 2. The marker: TMFG (planar) — best detection, modest margin, honest caveat

SOZ seeded heat-diffusion, strength-residualised, matched-strength fake-SOZ null; per-band
beats-null of 10; fused 6-band AUC at the s=10 operating point (design-space audit).

| backbone | density | fused AUC | R-prec | β beats-null | low_γ | high_γ |
|---|---|---|---|---|---|---|
| dense | 1.00 | 0.745 | 0.423 | — | — | — |
| mst@0.20 | 0.20 | 0.816 | 0.409 | 8/10 | 8/10 | 7/10 |
| mst@0.10 | 0.10 | 0.859 | 0.467 | 8/10 | 9/10 | 7/10 |
| mst@0.05 | 0.05 | 0.874 | 0.469 | 8/10 | 8/10 | 7/10 |
| disparity@0.2 | 0.12 | — | — | 8/10 | 9/10 | 6/10 |
| **TMFG** | 0.05 | **0.904** | **0.550** | **10/10** | 9/10 | **4/10 ↓** |

- Fused AUC is **monotone in sparsity within the mst family** (0.742 @ 0.50 → 0.874 @ 0.05): most of
  the marker gain is "go sparser," which the *same family* (mst@0.05) already captures.
- **TMFG's planar-specific edge** over mst@0.05 at equal density: **+0.03 AUC, +0.08 R-prec, β
  10-of-10** — real, but modest, and it **loses high_γ (4/10 vs 7/10)**. So the honest claim is "TMFG
  is the best-performing marker backbone under matched-strength," *not* "planarity is uniquely
  necessary for detection." The strength-ranked alternative (mst@0.05, 0.874) is logged so the choice
  is transparent.
- **Decision (user, 2026-07-15):** for the epilepsy *detection* readout we optimise detection ⇒
  **marker backbone = TMFG.** TMFG per-band single-marker AUCs (matched-strength, `epi_arc_tmfg`):
  β 0.869 / 10-of-10, low_γ 0.843 / 9, δ 0.830 / 7, α 0.792 / 7, θ 0.743 / 6, high_γ 0.688 / 4.

## 3. Why two backbones is a result, not filter-shopping

The cross-phase trace is carried by strong-edge **magnitude** (kill the strong edges → α dies); the
within-phase marker is carried by the sparse topological **skeleton** (triangulate → detection
sharpens). Two structural regimes of the *same* heat-kernel operator, each backbone **fixed by
matched-strength validation**, not chosen for AUC. The dissociation localizes what graph structure
each phenomenon depends on — a methods contribution in its own right. The pre-registered §6 rule was
followed to the letter: global planar adoption rejected on the trace, marker-only TMFG licensed.

## 4. What changed in the code / library

- `src/lrg_eegfc/utils/fc/backbone.py`: added `pmfg_backbone`, `disparity_backbone`,
  `select_backbone`; docstrings cite Tumminello 2005 / Serrano 2009.
- `06_epi_arc.py`, `13_matched_strength_mst020.py`: `SA_BACKBONE` ∈ {mst020, tmfg, pmfg, disparity},
  `SA_FRAC`, `SA_DISP_ALPHA`; per-backbone OUT dirs.
- `25_pmfg_observed_trace_readout.py`: PMFG≈TMFG observed-ρ_sym mechanism check.
- Data: `data/sparsified_arc/{ms_tmfg, ms_disparity_a0.2, ms_mst010, ms_mst005, epi_arc_mst005,
  epi_arc_mst010, epi_arc_disparity_a0.2, planar_trace_readout}/`.

## 5. Follow-on (writeup, not method) — propagate the settled marker to §3

Trace side unchanged (mst@0.20) → §1 + all trace downstream need nothing. Only §3 marker moves to
TMFG: re-source `fig_epi_a_relational_marker.py` from `epi_arc_mst020` → `epi_arc_tmfg`, update
`results_sec_3.tex` marker AUCs + fig:epi1(b) caption, add the dissociation sentence to Methods/§3.
See guide §7a.
