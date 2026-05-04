---
name: 2026-04-29_critical-post-mortem
type: post-mortem
era: COHORT_N10
status: current
created: 2026-04-29
updated: 2026-04-29
pointers:
  - .agents/reports/2026-04-24_post-mortem-scalar-session.md
  - .agents/reports/2026-04-25_task-trace-audit-and-recovery.md
  - .agents/reports/2026-04-27_task-persistence-reconciliation.md
  - .agents/reports/2026-04-28_residual-subspace-diagnostic.md
  - .agents/reports/2026-04-28_psi-tau-scan-verdict.md
  - .agents/reports/2026-04-29_measure-correctness-audit.md
  - .agents/guides/02_methods/lrg-framework-guide.md
  - data/audit/cohort_coverage_matrix/triangulation_n10_imcoh_abs.csv
---

# Critical post-mortem — what worked, what didn't, where to go next

## Renormalization head

**Two months of work since the imcoh_abs reset (2026-04-15) tested at
least 13 distinct task-trace measures across five families
(continuous-distance, partition-divergence, subtree-identity,
per-leaf, spectral). Of those, exactly one cohort-wide claim survives
controls: a δ-band partition-divergence ridge at integer-k cuts.
Everything else is either uncontrolled, methodologically dead, or
silent under controls. The four critical lessons are: (a) every
spectrum-scalar quantity (entropy, specific heat, Ψ-as-selector,
rank-1 residual subspace) has been killed by our continuous spectrum;
(b) every cohort-aggregation strategy that pre-registers binary
thresholds masks the underlying gradient; (c) the dendrogram is
load-bearing but it compresses eigenvector information that the user
believes carries the actual diffusion-pattern reorganisation; (d) the
"task trace" claim has never been sharpened to a single testable
predicate — we have been target-shopping. The decision is to pivot
from the predicate-based geometric ladder to eigenvector-direct
spectral probes (E1–E3), keeping the LRG primitives but replacing
scalar-spectrum collapses with full-basis comparisons.**

---

## TL;DR (table of contents)

1. Six families of measures were tried since the imcoh_abs reset; thirteen distinct measures total.
2. Five things we know with high confidence.
3. Five things we thought we knew but don't.
4. Three open scientific questions.
5. Five methodological lessons learned.
6. Pivot direction: eigenvector-direct spectral probes (the E-series).

---

## 1. What was tried (chronological since the imcoh_abs reset, 2026-04-15)

**Abbreviations.** *imcoh_abs* = `<|ImCoh|>_f` (frequency-resolved imaginary coherence, magnitude-averaged). *FC* = functional connectivity. *VI* = variation of information. *KC* = Kendall–Colijn tree distance. *MRL* = module-retention landscape. *CBR* = community-boundary residual. *CNP* = cophenetic neighbourhood profile. *MSPC* = multiscale partition coherence. *Ψ* = partition-stability index per Villegas 2025. *PRL* = partition-resolution-locked.

| Family | Measure | Era | Status | Verdict |
|---|---|---|---|---|
| Continuous distance | H2c-shared (`Spearman ρ` on `(Δ_task, Δ_rest)` upper-triangle, shared `D_pre`) | 2026-04-22 | superseded | 6/6 cohort-positive but artefactual: shared baseline correlates Δ_task and Δ_rest via additive noise |
| Continuous distance | H2c-split-baseline / Run A (disjoint halves of `rest_pre`) | 2026-04-26 | current | 0–3/10 cohort-positive across bands; controlled silence |
| Continuous distance | H2c-cross-probe / Run B | 2026-04-26 | current | α / β / γ_l survive at 3-control consolidation; not load-bearing |
| Continuous distance | H2c-drift-floor / Run C (within-session halves) | 2026-04-25 | current | drift baseline; flat near zero across bands |
| Continuous distance | H2-RAW / H2-FROB (raw-D global similarity) | pre-reset | superseded | anatomy-dominated; failed |
| Partition-divergence | Δ_VI(k) / Δ_H(k) / Δ_NMI(k) / Δρ_H2d(k) | 2026-04-22 | current | δ ridge k=20–32 at ≥8/10; α/β/γ_h partial signal; θ/γ_l ergodic; CSV is what the v1 cohort-coverage matrix consumes |
| Partition-divergence | Δ_VI(h_rel) at fixed fractional dendrogram depth | 2026-04-26 | current | δ ridge dissolves at fixed h_rel — locked to per-patient resolution |
| Subtree-identity | MRL (hard-Jaccard task-anchored module count over h_rel grid) | 2026-04-25 | superseded | cohort-null (TAM cloud sits on diagonal at high Jaccard; trace zone empty) |
| Subtree-identity | CBR variants: legacy / task-anchored / containment / consensus-subtree | 2026-04-25 | superseded | classifier-failure analysis showed all four conflate scores; replaced by Cohesion-CBR |
| Subtree-identity | Strict-J (Jaccard ≥ 0.9) trace modules at integer-k cuts | 2026-04-26 | current | cohort-null: 1 patient × 1 band × 1 partition has any strict-J module (Pat_02 alpha) |
| Per-leaf | Cohesion-CBR (Jaccard + corner-distance affinity, soft) | 2026-04-25 | current | per-patient T-leaves visible; no cohort-wide test until cohort-coverage matrix |
| Per-leaf | CNP — cophenetic neighbourhood profile (per-leaf Spearman ρ on cophenetic vector) | 2026-04-25 | current | per-leaf affinity; uncontrolled cohort-wide |
| Per-leaf × multiscale | MSPC — per-(leaf, k) Jaccard cluster-mate sets | 2026-04-26 | current | per-(leaf, k) affinity; uncontrolled cohort-wide; 407k-row CSV |
| Spectral / scalar | Ψ(n;τ) and Ψ_L diagnostic (partition stability per Villegas 2025) | 2026-04-28 | dead | 84% of cells put `argmax_n Ψ` at the deepest split, 16% at trivial coarse cuts, 0.7% at any intermediate cut. Ψ-as-selector dead. |
| Spectral / scalar | Residual-subspace alignment `α_1` (rank-1, leading eigenvectors of `D^task − D^pre` vs `D^post − D^pre`) | 2026-04-28 | dead | mean `α_1 − null_p95 = −0.066` across 216 cells; below permutation null cohort-wide; β-band partial signal explained by near-rank-1 residual dominance, not task-specific |
| Spectral / scalar | Entropy curve / specific heat (S(τ), C(τ)) | (never built) | removed 2026-04-29 | non-informative for our continuous spectrum; permanently retired from the ladder |
| Triangulation | Cohort-coverage matrix v1 (8-rung × 6-band) | 2026-04-29 | current | δ partition-resolution-locked; α/β/γ_l/γ_h ergodic under strict gate; α/γ_h signal masked by ridge≥3 + dissenter inclusion |

Plus methodological / infrastructure:
- Phase 0 measure-correctness audit (2026-04-29): 15/15 measures `current`, 14/14 hand-recomputations pass, 18 superseded scripts archived.
- LRG framework guide (2026-04-29): canonical formulas + 5 verification snippets passing against live code; flagged our case as outlier wrt Villegas papers.
- Decision-rules pre-registration (2026-04-29): cohort threshold ≥8/10, ridge ≥3, BH-FDR q ≤ 0.05, dissenter handling, 7 triangulation predicates.

---

## 2. What we know with high confidence

**(K1) δ-band partition-divergence ridge at integer-k cuts is real.**
Cohort `frac_pos ≥ 0.8` for 13 contiguous k-values (k=20–32), peaks at 9/10–10/10 patients positive at k=26–28, Wilcoxon q=0.014 BH-FDR-corrected within rung. Robust under dropping dissenters (frac_pos_robust ≥ 0.875 at the same k). Survives the within-baseline drift floor (`Δ_VI > drift_dVI`). Methodological caveat: locked to per-patient partition resolution — dissolves at fixed fractional dendrogram depth.

**(K2) Strict-Jaccard subtree identity is null cohort-wide across every band.**
With J_min = 0.9 the trace-module cohort returns 1 row total (Pat_02 alpha) across n=10 × 6 bands × 2 phase-pairs. No cohort-wide named module. This *is* a finding: there is no shared anatomical / topological module that recurs across the cohort.

**(K3) The continuous-distance pair-coherent claim does not survive the split-baseline control.**
H2c-shared was 6/6 cohort-positive; H2c-Run-A (disjoint halves of rest_pre) is 0–3/10 across all bands. The original "task trace at the pair-distance level" claim was an artefact of shared-baseline noise correlation. Run B (cross-probe) recovers α/β/γ_l partially; not load-bearing for headlines.

**(K4) Spectrum-scalar quantities are dead on our data.**
Ψ-as-selector: `argmax_n Ψ` collapses to deepest split or trivial coarse cuts in 99.3% of cells. Ψ_L per-cluster local stability: `n_psi_L > 1.0` is zero across all (band, phase). Residual-subspace `α_1`: below permutation null in mean across 216 cells. Specific heat C(τ) and entropy S(τ): no discrete peaks (continuous spectrum). Every spectrum-scalar approach has failed in the same characteristic way: either the quantity is dominated by trivial extremes or it's mathematically trivial under our continuous-spectrum regime.

**(K5) The infrastructure works.**
LRG primitives verified against the canonical Villegas formulas (5 verification snippets pass). Halves cache supports within-baseline nulls. Phase 0 audit confirms every measure consumed by the matrix v1 is correctly implemented. Cohort metadata + ledger + audit verdicts gate downstream consumers.

---

## 3. What we thought we knew but don't

**(D1) "Task leaves a multiscale trace cohort-wide".**
This was the headline coming out of H2a (multiscale VI ≥8/10) at n=9. Under the 2026-04-29 controlled re-evaluation with `Δ_VI > drift_dVI` and ridge ≥ 3 and dissenters in the denominator, only δ survives at the integer-k partition rung; α and γ_h have partial signal that the strict gate trims. The original headline was **methodologically permissive** (no within-baseline null, no ridge condition, dissenter handling implicit). The current cohort-wide claim is much weaker.

**(D2) "Triangulation across rungs strengthens the verdict".**
The cohort-coverage matrix triangulation predicates conjugate per-rung verdicts under strict thresholds. In practice this *masks* gradients (α/γ_h's near-positive signal), embeds a hidden assumption that real trace must be both continuous AND modular AND topology-coherent (which by construction excludes purely-modular traces from "headline"), and creates label artefacts that look like ranked-strength but are actually rung-coverage statements. Triangulation is informative as a diagnostic map but it is **not a reliable cohort-wide trace detector**.

**(D3) "ImCoh-derived dendrograms preserve the spatial information".**
The dendrogram compresses eigenvector information into an ultrametric tree. We have been deriving every measure (L1, L4, L5_k, L5_hrel, L6, L7) from this single tree and conjuncting their verdicts. The user's intuition (now operationalised as the eigenvector-direct pivot) is that **diffusion-pattern reorganisation lives in the eigenvectors, not in the tree-discretised hierarchy**. Some of our cohort silence may be the dendrogram averaging out signal that the eigenvectors retain.

**(D4) "Per-leaf evidence (CBR, CNP, MSPC) is a cohort census".**
Per-leaf measures produce per-patient PDFs and per-(patient, band) leaf assignments, but until the cohort-coverage matrix v1 there was no cohort-level census joining them. Even now they are uncontrolled in v1 — the row-count sanity-check passes but the per-leaf 2D-corner-distance affinity formulas are not predicate-level verified.

**(D5) "Pre-registered thresholds are a substitute for sensitivity analysis".**
≥8/10 cohort, ridge ≥ 3, q ≤ 0.05 BH-FDR, eligibility 0.5 — all defensible, none derived. We pre-registered to avoid post-hoc target shopping, but we never ran the sensitivity sweep that would tell us which verdicts are stable across reasonable threshold choices. Until that sweep runs, every verdict carries an implicit "stable under these specific gates" qualifier.

---

## 4. What's still genuinely open

**(O1) Does the δ partition ridge reflect cognitive task processing or a non-cognitive confound?**
δ is the same band that captures slow drifts, drowsiness modulation, and slow ERP components. The δ PRL verdict is *real* in the sense of "modular reorganisation cohort-wide-survives-drift-floor-at-integer-k", but biologically interpreting it as task-induced requires a sham-task null (e.g. permute task labels or compare task-test to a passive resting block) that does not currently exist in the rebuild plan.

**(O2) Is the diffusion-pattern reorganisation visible only in eigenvectors?**
The dendrogram-based ladder has been silent or ergodic for α/β/γ_h under strict controls. The user's hypothesis is that the *eigenvector basis* of the Laplacian carries the trace even when the dendrogram does not. The E-series (E1 spectral subspace alignment, E2 diffusion-map embedding drift, E3 eigenvector-localisation shifts) tests this directly and is the immediate next priority.

**(O3) Is "cohort-wide" the right frame given implant heterogeneity?**
n=10 sEEG patients have *different* electrodes at *different* brain regions. Asking for shared dendrogram structure ignores that the underlying physical substrate differs. A "cohort-wide modular reorganisation" verdict is in some sense a topology-only claim (independent of where the modules sit anatomically). Open whether per-patient rich case studies + meta-analytic summary is methodologically more honest. *Not in scope for the E-series; flagged for paper-framing decision later.*

---

## 5. Methodological lessons learned

**(L1) Continuous spectrum kills every scalar-spectrum probe.**
Ψ, Ψ_L, S(τ), C(τ), τ\*, rank-1 residual `α_1` — all collapsed to triviality. Lesson: if the quantity is a scalar function of the eigenvalue distribution, it will be uninformative on our data. **Use the eigenvectors directly.**

**(L2) Pre-registered binary thresholds need a graded companion.**
The matrix v1's binary positive/silent verdicts hide the gradient (α near-positive, γ_h narrow-significant, β truly weak — read as three "ergodic" cells). Lesson: report cohort `frac_pos` distribution and ridge profile alongside binary verdict; never let the binary cut be the sole headline.

**(L3) Triangulation predicates encode hidden models of "what trace looks like".**
The headline-triangulated predicate requires both continuous AND modular AND topology-coherent agreement. A purely modular trace (like δ) is structurally excluded from "headline". Lesson: triangulation predicates are **methodological choices**, not neutral aggregators. Document them as such; do not treat them as derived consequences of the data.

**(L4) The "task trace" claim has been target-shopping.**
Five families of measures, thirteen distinct operationalisations, no single sharpened predicate. Lesson: **before adding any new measure, sharpen the claim**. What exactly do we want to demonstrate? Per-pair distance shift? Module identity persistence? Diffusion-mode rotation? Each is a different paper. Pick one (or two) and stick with it.

**(L5) Don't treat the dendrogram as the only LRG output.**
The dendrogram is a discretisation of a continuous object (the cophenetic ultrametric matrix from the density operator at τ=1/λ_max). It throws away the eigenvector basis and the τ-resolved kernel structure. Lesson: when a tree-based measure goes silent, ask whether the underlying matrix object (ρ̂(τ), K̂(τ)) might still carry the signal.

---

## 6. Pivot direction — eigenvector-direct spectral probes

The decision: **keep the LRG primitives (Laplacian, heat kernel propagator, density operator, communication distance, dendrogram) but pivot from spectrum-scalar collapses (S(τ), C(τ), Ψ) and tree-discretised conjunctions (the current 8-rung ladder) to direct eigenvector-basis probes (E1, E2, E3) that compare diffusion-pattern reorganisation across phases without compressing to a tree or to a scalar.**

Concrete plan owed (separate document next): `2026-04-29_eigenvector-direct-pivot-plan.md` defining

- **E1 — Spectral subspace alignment.** Top-k eigenvectors of `L̂` per phase as N×k subspace; principal-angle / Grassmann-distance comparison across phases.
- **E2 — Diffusion-map embedding drift.** Coifman-Lafon embedding via top-k eigenvectors weighted by `exp(−τλᵢ/2)`; per-node embedding-space distance across phases.
- **E3 — Eigenvector-localisation shifts.** Per-eigenvector inverse participation ratio; cohort distribution shift across phases.

Each E-rung has a within-baseline null (halves cache), pre-registered cohort threshold (taken from existing decision-rules), and a row in the cohort-coverage matrix v2.

**Critical distinction from the failed residual-subspace probe:** that
probe used rank-1 leading eigenvectors of *residuals* (`D^task − D^pre`,
`D^post − D^pre`) and was killed by near-rank-1 dominance making
alignment mathematically trivial. The E-series uses the **full top-k
eigenvector basis of the Laplacian directly** (not residuals), compared
phase-to-phase. Different math, different failure modes; not the same
trap.

What the E-series does NOT solve:
- The "task trace" claim is still unsharpened. (E-series sharpens it to "diffusion-mode reorganisation"; that's a specific, testable claim.)
- The implant heterogeneity / cohort-framing question (O3) remains open.
- The biological interpretation (O1) still needs a sham-task null.

What the E-series can deliver:
- A direct test of the user's hypothesis that the trace lives in the eigenvectors.
- A row in the cohort-coverage matrix that is genuinely independent from the dendrogram-based rungs (different math, different invariances).
- A clean "yes" or "no" on whether spectral subspace reorganisation is cohort-wide visible — something the dendrogram-based ladder cannot answer.

---

## 7. Outcome of this post-mortem

- The δ partition-resolution-locked verdict stays as a real result.
- The cohort-wide multi-band trace claim is downgraded — controls + thresholds remove most of it.
- The dendrogram-based ladder will not be extended further (no §7.3, §7.4, §7.5, §7.6 control scopes for the existing tree-based rungs).
- The §7.2 KC λ-sweep scope is paused pending E-series results — KC is a tree distance and may inherit the dendrogram-compression issue.
- The new direction is the E-series eigenvector-direct probes.
- The cohort-coverage matrix stays as the matrix v1.1 (with E-rungs added in v2 once they land), but the binary triangulation should grow a graded companion.

Detailed E-series plan in the paired document:
`.agents/plans/active/2026-04-29_eigenvector-direct-pivot-plan.md` (next).
