---
name: 2026-04-28_psi-tau-scan-verdict
type: handoff
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-28
updated: 2026-04-28
pointers:
  - data/audit/psi_tau_scan/psi_diagnostic.pdf
  - data/audit/psi_tau_scan/psi_grid.csv
  - .agents/reports/2026-04-28_residual-subspace-diagnostic.md
---

# Ψ(n;τ) and Ψ_L diagnostic — verdict

**Renormalization head.** The 3-patient × 6-band × 3-phase × 8-τ scan of the
Villegas 2025 partition stability index Ψ(n;τ) on imcoh_abs LRG dendrograms
shows that **84% of cells (362/432) put `argmax_n Ψ` at the deepest split
(n ≥ N−3, near singletons)**, **16% at trivial coarse cuts (n ≤ 5)**, and
**only 0.7% (3/432) at any intermediate cut**. The local stability index Ψ_L
(per-cluster log-Δ branch length) is uniformly weak: `Ψ_L p95 ∈ [0.1, 0.3]`
across all (band, phase), `n_psi_L > 0.5` has median 0 (max 2 per tree),
`n_psi_L > 1.0` is **zero everywhere**. The dendrograms are too tight in
log-Δ for stability-based selection to anchor anything. Ψ-as-selector is
dead on this data; Ψ_L cannot anchor a per-cluster trace claim.

## Inputs

- 3 patients (Pat_02, Pat_06, Pat_03), all 6 bands, 3 phases.
- 8 log-spaced τ-values per case on `τ ∈ [1/λ_max, 1/λ_gap]`.
- `D(τ) = 1/ρ(τ)` (canonical Villegas 2025 ultrametric); average-linkage
  dendrogram per τ; Ψ via existing
  `lrg_eegfc.visuals.lrg.compute_partition_stability_index` (paper formula).

## Headline numbers (`argmax_n` distribution, 432 cells)

| regime | range | count | % |
|---|---|---|---|
| **deep** | n ≥ 100 (≈ N−3) | 362 | 84% |
| **trivial coarse** | n ≤ 5 | 67 | 16% |
| **mid** | 6 ≤ n ≤ 99 | 3 | 0.7% |

`argmax_n Ψ` is bimodal (deep / trivial) — never at a meaningful intermediate
partition. Operationally useless as a partition selector.

## Ψ_L (per-cluster persistence)

| metric | typical range | reading |
|---|---|---|
| Ψ_L p95 | 0.10 – 0.30 (across all band, phase) | most stable cluster has < 1 decade of "stem" in log-Δ |
| n_psi_L > 0.5 (med, max) | 0, 2 | very few clusters stand out |
| n_psi_L > 1.0 (med, max) | 0, 0 | no cluster has order-of-magnitude stability |

Filtering candidate trace clusters by `Ψ_L > threshold` would yield zero
candidates in most cells.

## Structured patterns (such as they are)

- **δ task_test → n=3** in Pat_02 and Pat_06; **δ/θ task_test → n=2** in
  Pat_03. A coarse 2–3-cluster split that emerges only in task, only on the
  artefact-suspect bands. Chasing it would chase the slow-drift / EMG ghost
  already flagged in prior verdicts.
- **All non-extremal-band cases** land at the deep regime — i.e. Ψ argmax
  says "the most stable cut is one merge above the leaves", which is
  formally a "no signal" answer in this framework.

## What this implies for downstream design

The Villegas 2025 cut machinery is canonical and correct, but the
*application* requires dendrograms with separated mesoscales. Our imcoh_abs
FC trees are uniformly tight: log-merge-distances are nearly equispaced, no
gap stands out enough for Ψ to lock onto it. This is consistent with C(τ)
being a single broad peak on these networks (which is what dense weighted FC
graphs without strong spectral gaps produce). The "no clear mesoscale"
property is a property of the data + FC method + cohort, not a flaw of
LRG.

Downstream consequence: any framework that *requires* Ψ-stable clusters as a
filtering step is unusable here. This kill chained directly into the
residual-subspace alignment diagnostic
([`2026-04-28_residual-subspace-diagnostic.md`](2026-04-28_residual-subspace-diagnostic.md))
which avoided Ψ but failed for a different (rank-1 dominance) reason.

## Outputs preserved

- `data/audit/psi_tau_scan/psi_grid.csv` — 432 rows.
- `data/audit/psi_tau_scan/psi_diagnostic.pdf` — 4 pages (argmax_n, top1/top2
  ratio, Ψ_L count, Ψ heatmaps for Pat_02).
- `scripts/01_compute/diagnostics/diag_psi_tau_scan.py` — reproducible.

## Decision

- Drop Ψ-as-selector. Confirmed by experiment.
- Drop per-cluster Ψ_L filtering for trace candidates.
- Composite findings (Ψ scan + residual-subspace) push the deliverable
  toward Option B: α/β under H2c continuous controls as the working
  headline. Drafting plan in a separate handoff.
