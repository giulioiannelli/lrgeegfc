---
name: task-persistence-investigation
type: guide
era: COHORT_N10
status: current
created: 2026-04-25
updated: 2026-04-25
pointers:
  - .agents/reports/2026-04-25_task-trace-audit-and-recovery.md
  - .agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md
  - .agents/reports/2026-04-24_h1-h4-vi-results.md
  - .agents/plans/active/2026-04-25_surface-multiscale-trace.md
---

# Task-persistence investigation

**Canonical home for every multiscale measure built to investigate task-induced
reorganization in the post-task resting state. New measures land here as
mathematically rigorous scope reports first; code follows. Mathematical
rigor is required so a careful reader can spot definition errors, edge cases,
and identifiability problems before any compute is launched.**

## Purpose

The standing scientific question (mirrored at the top of `CLAUDE.md`):

> Does `task_test` leave a band-specific, multiscale structural trace in
> `rest_post` LRG dendrograms that is cohort-wide (≥ 8/10 at n=10 under
> `imcoh_abs`)?
>
> Cohort returned to n=10 on 2026-04-25 with Pat_14 vendor-replacement.
> Canonical reformalization at the leaf-set level (regimes P/T/R/RA):
> [`2026-04-25_task-trace-canonical.md`](2026-04-25_task-trace-canonical.md).

Existing evidence is fragmented across multiple reports and operationalizations
(H1 VI(k), H2c continuous drift, H2d block-pair persistence, partition-level
cluster-permutation, H1-topo bipartition overlap, H2e split-half drift floor).
Each measure tests one slice of "task trace" — pairs, partitions, bipartitions,
distance shifts. None of them directly measures *which subtrees of the
hierarchy at which scales* are task-induced and persistent.

This folder collects the descriptive multiscale measures we develop to surface
that exact question. Explicitly **not** scalar-gate hypothesis tests at
n=9 + FDR q < 0.05 m=6 — the post-mortem
(`.agents/reports/2026-04-24_post-mortem-scalar-session.md`) closes that loop.

## File layout

```
task-persistence-investigation/
├── README.md                                        # this file (index)
└── YYYY-MM-DD_<measure-name>.md                     # one scope report per measure
```

Naming: `YYYY-MM-DD_<short-kebab-name>.md`. Date = day the scope was first
drafted; do not change on revisions, update the `updated` frontmatter field
instead.

## Required scope-report structure

Every report under this folder MUST contain, in order:

1. **Frontmatter** — `name, type=scope, era, status ∈ {draft, current, superseded}, created, updated, pointers`.
2. **Renormalization head** — 1–2 sentences. The juice. A reader skimming the
   first paragraph must come away knowing what the measure is and what
   question it answers.
3. **Notation** — symbols, indices, sets, and the mathematical objects
   the measure operates on. Define once; never overload.
4. **Definitions** — each predicate / function / aggregator stated as a
   formula with its domain and range. Cite literature when borrowing
   (Kendall–Colijn 2016, Robinson–Foulds 1981, Maris–Oostenveld 2007, etc.).
5. **Properties** — range, monotonicity, invariance, identifiability,
   complexity. State *what the measure cannot detect* (negative properties
   are as important as positive ones).
6. **Caveats & failure modes** — probe bias at coarse `h_rel`, edge effects
   at the root, ties in heights, parameter sensitivity, finite-N artifacts,
   anatomy-dominated baselines. Each caveat with a mitigation.
7. **Pseudocode** — language-agnostic, line-oriented, every loop bound
   explicit. Not Python — algorithmic.
8. **Visualization spec** — figure axes, colormaps, annotation rules,
   what a reader is supposed to *see*. Include reading rules
   ("vertical bands of color = nested persistence"; "horizontal contrast = band heterogeneity").
9. **Connection to prior tools** — a paragraph or table mapping this
   measure to existing operationalizations (H2c / H2d / H1-topo / etc.).
   What does it complement; what does it subsume; what does it not replace.
10. **Implementation plan** — script path, output paths, library entry points
    (must reuse `lrg_eegfc.utils.metrics.tree*`; no new helpers without a
    promotion-to-library justification per `coding-rules.md`).
11. **Open questions** — parameter defaults, scope decisions deferred.
    Marked explicitly so a future reader knows what's still load-bearing.

## Status semantics

- `draft` — written, not yet implemented or sanity-checked.
- `current` — implemented and feeding into a live analysis.
- `superseded` — replaced by a later scope; keeps the file (don't delete);
  add a top-of-body pointer to the successor.

## Active scope reports

| Date | Measure | Status | Tagline |
|------|---------|--------|---------|
| 2026-04-25 | [Task-trace canonical (P/T/R/RA reformalization)](2026-04-25_task-trace-canonical.md) | draft | **Definitional umbrella.** Names the four regimes (Persistence / Trace / Reset / Rearrange) at the leaf-set level and pins shared primitives (`J_min = 0.9`, `k_min = 3`, native heights) so MRL, CBR, and the audit-and-recovery surfacing tooling all speak the same vocabulary. T = `¬present_pre ∧ present_task ∧ present_post` is the central scientific claim; cohort-wide threshold `Π_T ≥ 8/10`. Residual claim, not global ("rest_post becomes task-like" is rejected by H2-RAW/H2-FROB). |
| 2026-04-25 | [Module-Retention Landscape (MRL)](2026-04-25_module-retention-landscape.md) | **superseded** | Discrete subtree-identity test using hard Jaccard thresholds. Returned cohort null (TAM `(J_pre, J_post)` cloud sits on the diagonal at high Jaccard; trace zone empty). Pathologies and reconciliation with CBR are in [the MRL ↔ CBR reconciliation](2026-04-25_mrl-vs-cbr-reconciliation.md). |
| 2026-04-25 | [MRL ↔ CBR reconciliation](2026-04-25_mrl-vs-cbr-reconciliation.md) | current | Post-mortem of the MRL arc (forward / forward−backward / TAM) and why it detected null while Cohesion-CBR shows per-patient traces on the same data: limbo-zone filtering, score conflation, scalar-field rendering. Frames the continuous-distance pivot. |
| 2026-04-25 | [CBR investigation (single-module-matching family)](2026-04-25_cbr-investigation.md) | active | Umbrella for the CBR variants drafted in the 2026-04-25 session: legacy / task-anchored / containment / consensus-subtree, plus the [classifier failure analysis](2026-04-25_cbr-classifier-failure-analysis.md) and [Cohesion-CBR](2026-04-25_cohesion-cbr.md) (Jaccard + corner-distance soft affinities, confidence-graded alpha — supersedes the four hard-threshold variants for the cohort census). |
| 2026-04-26 | [Continuous-trace matrix](2026-04-26_continuous-trace-matrix.md) | draft | Visual proof on the ultrametric distance matrix `D` itself — not the dendrogram. Per-pair distance shifts `Δ_task = D_test − D_pre` and `Δ_rest = D_post − D_pre`; sign-agreement map `σ(i,j)`; per-pair scatter; cohort `ρ` per band. Antidote to MRL's discretisation; visual companion to H2c (53/54 cells positive at n=9 prior). |
| 2026-04-29 | [Measure-correctness audit (Phase 0, BLOCKING)](2026-04-29_measure-correctness-audit.md) | draft | **Phase 0 of the task-trace rebuild.** Before any new measure or matrix-cell evaluation, audit every existing measure end-to-end — script↔scope-report fidelity (φ_fid), cohort/era correctness (φ_coh), FC-method order (φ_fc), library-helper reuse (φ_help), per-cell smoke reproducibility (φ_smoke), naming-convention compliance (φ_name) — and inventory every script under `scripts/01_compute/` for dead-duplicate archival. Two CSVs (`measure_audit_n10_imcoh_abs.csv`, `script_inventory_n10.csv`) plus a narrative report gate every downstream rebuild step. |
| 2026-04-29 | [Decision rules (pre-registration)](2026-04-29_decision-rules.md) | draft | **Pre-registered triangulation logic for the cohort-coverage matrix.** Per-rung positivity tests (L1–L7 + L5(h_rel)), per-band cohort verdicts `V(b, r) ∈ {positive, negative, silent, ineligible}`, dissenter handling (Pat_02/03 counted in denominator with "×"/"+" markers), and cross-rung triangulation predicates `T(b)` ∈ {headline-triangulated, partition-resolution-locked, continuous-only, anatomy-suspect, ergodic, uncontrolled-exploratory, inconsistent}. Locked before matrix evaluation runs. |
| 2026-04-29 | [Cohort-coverage matrix (data-flow / artefact)](2026-04-29_cohort-coverage-matrix.md) | draft | **Paired-PR sibling to the decision-rules scope.** Defines the per-cell CSV (one row per `patient × band × rung × scale_axis × scale_value`), the per-(band × rung) `band_verdict` CSV that materialises `V(b, r)`, the per-band `triangulation` CSV that materialises `T(b)`, and the band × rung heatmap PDF that supersedes `task_trace_band_k_n10_imcoh_abs.pdf` as the headline. Carries the consumer-list (which CSV/NPZ each rung pulls from), the eligibility gates, the dissenter map, and the v1-vs-v2 split (`uncontrolled` cells in v1; replaced after §7 controls land). Together with decision-rules.md fully specifies `audit_18_cohort_coverage_matrix.py`. |
| 2026-04-28 | [Functional tree distance (`δ(τ)` curves on `ρ` and `1/ρ`)](2026-04-28_functional-tree-distance.md) | **superseded** | Per-(patient, band, phase-pair) τ-curve on `D = 1/ρ` (rank `δ_S`) and `ρ` (height `δ_P`), integrated over `[1/λ_max, τ\*]`, gated by within-phase split-half noise floor. Implemented at n=10 with within-baseline null. Verdict: **3 trace cells at strict rank gate, 10 at permissive height gate, no band ≥ 3/10 — eighth measure to fail the cohort task-trace gate at n=10**. See [`2026-04-28_functional-tree-distance-verdict.md`](../../reports/2026-04-28_functional-tree-distance-verdict.md). |
| 2026-04-29 | [E1 — Spectral subspace alignment](2026-04-29_e1-spectral-subspace-alignment.md) | draft | **First eigenvector-direct rung in the pivot plan.** Grassmann chordal distance `d_chord(V_a, V_b) = sqrt(k − Σ σᵢ²)` between top-k non-trivial eigenspaces of L̂; trace `Δ^E1_k = d(V_test, V_post) − d(V_pre, V_post) < 0`; halves-cache null. Sees subspace rotation invisible to L1 (pair-distance) and L5_k (modular swap). k-grid `{2,3,5,8,13,21}`, ridge ≥ 2. Subspace-invariant: no eigenvector sign/order alignment needed (E2/E3 do). |
| 2026-05-04 | [Section 5 LRG-trace investigation plan](2026-05-04_lrg-trace-investigation-plan.md) | draft | **Methods review + empirical-pilot plan for manuscript Section 5.** Lifts the substrate triangle `T(p, b) = d(TT, RPost) − d(RPre, TT)` to the LRG hierarchy at τ=1/λ_max on three primitives — ρ(τ), D(τ), T(τ). Five-part deliverable: methods review (every repo-resident probe + literature candidates, with status), empirical pilot of the global per-band trace measure (5 candidates: D-rank triangle, CTM σ-aggregate, Grassmann triangle, KC λ-blend triangle, VI(k) multiscale band×k heatmap as a result-per-se), localization pilot (Cohesion-CBR re-anchored on rest_pre + triplet-consensus + MSPC f_trace), bidirectional cross-validation (LRG-side per-pair only, per user directive), pre-registered selection rule for the headline. Documents shelved/killed probes (FTD, RSA, E1 pairwise, Ψ_L) by elimination so the recommendation is grounded in the full repertoire. |
| 2026-05-06 | [RF(k) — clade persistence across phases](2026-05-06_rf-k-multiscale-measure.md) | **parked** | **Parked 2026-05-07, not in manuscript.** Hard + soft + per-clade strict 4-mode taxonomy implemented (audit_48 / 48b / 48c / 48d); cohort numbers 16/13/32/48 trace/reset/persist/rearrange across 60 cells; soft RF passes BH at α + β. Set aside because Jaccard's size-sensitivity prevents clean detection of graded clade reorganization — KC λ=0 (audit_36 + audit_46) and the audit_47-50 KC module-view 4-mode family handle the same scientific question more cleanly. Output kept on disk as a strict-identity sensitivity check. Verdict at `.agents/reports/2026-05-07_rf-clade-persistence-cohort-verdict.md`. |
| 2026-05-07 | [KC reset modules](2026-05-07_kc-reset-modules.md) | draft | **Symmetric inverse of the KC trace-module visualization (audit_47).** A reset module is a leaf set that exists as a coherent subtree in BOTH rest phases (rPre and rPost) but fragments during task. Anchored at rPre subtrees: gate 1 `J(L_pre, rPost) ≥ 0.65`, gate 2 `∀ tt: J(L_pre, tt) < 0.5`, gate 3 `|C_tt(L_pre)| ≥ 2·|L_pre|` (containment fragmentation guard, mirrors trace), gate 5 symmetric `(L_post, tt)` fragmentation. 3×3 dendrogram+network grid mirroring audit_47 with phase semantics swapped: coherent in pre and post, scattered in tt. `audit_48_kc_reset_module_view.py` planned. Decomposes the cohort scalar `T_d^KC < 0` into trace + reset components per patient. Does not detect height-only resets (same blind spot as audit_47 for height-only traces). |
| 2026-05-07 | [KC rearrangement modules](2026-05-07_kc-rearrangement-modules.md) | draft | **Counter-example to trace.** A rearrangement module is a leaf set coherent only in rest_post — fragmented in BOTH rest_pre and task_test. Anchored at rPost subtrees with both rPre and tt fragmentation gates: `∀ pre: J(L_post, pre) < 0.5` AND `|C_pre(L_post)| ≥ 2·|L_post|` AND symmetric for tt. By construction disjoint from trace (gate fails for tt-coherent), reset (gate fails for rPre-coherent), and anchor (gate fails for both). 3×3 visualization with rows = size strata (≥3 / ≥5 / ≥8) since λ-blend doesn't apply (no matched pair). `audit_49_kc_rearrangement_module_view.py` implemented 2026-05-07. First sweep: rearrangement is **abundant** (~10 modules per cell across all 6 bands, every patient) — most rPost coherence is emergent, not task-seeded. This complicates "rPost coherence is task-locked" — that claim is supported only by the SUBSET of rPost coherent subtrees that match a tt subtree (the trace population). Strict-emergent subcase of `terminology.md`'s `emergent` (also requires tt fragmentation). |
| 2026-05-07 | [KC anchor modules](2026-05-07_kc-anchor-modules.md) | draft | **Fourth class — the null.** An anchor module (user term: "persist") is a leaf set with a coherent matching subtree in ALL THREE phases. Triple-Jaccard gate: `J(L_pre, L_tt) ≥ 0.6 ∧ J(L_pre, L_post) ≥ 0.6 ∧ J(L_tt, L_post) ≥ 0.6`. Disjoint with trace, reset, and rearrangement by construction. 3×3 grid identical in shape to audit_47/48 (rows = λ regimes; height bars at all three phases when λ ∈ {0.5, 1}). `audit_50_kc_anchor_module_view.py` implemented 2026-05-07. Anchors are the BASELINE against which trace, reset, and rearrangement are deviations. High anchor count → brain mostly preserves module structure across rest→task→rest. The trace, reset, rearrangement, anchor partition is the complete cross-phase taxonomy from `terminology.md`. |

## Rule — every new task-trace tool lands here first

When developing a new descriptive multiscale measure for the task-trace question:

1. Write the scope report **first**. No code before scope.
2. Mathematical rigor is non-negotiable — every definition, every predicate,
   every aggregator stated as a formula. A reader must be able to spot
   errors from the math alone.
3. Cite primitives in `lrg_eegfc.utils.metrics.tree*` and
   `lrg_eegfc.utils.metrics.hypothesis`. No private re-implementations.
4. Cross-reference what the measure complements or supersedes among
   `2026-04-24_*` reports.
5. Update the **Active scope reports** table above with date, link, status, tagline.
6. State explicitly what the measure does **not** measure — the negative
   property is part of the contract.

This rule is mirrored in the project-level enforcement list:
`.agents/guides/04_rules/never-always-list.md` and the seed list at the top
of `CLAUDE.md` / `AGENTS.md`.

## Cross-references

- Existing evidence handoff: `.agents/reports/2026-04-24_multiscale-task-trace.md`
- VI(k) tables: `.agents/reports/2026-04-24_h1-h4-vi-results.md`
- Why scalar gates failed: `.agents/reports/2026-04-24_post-mortem-scalar-session.md`
- Active plan: `.agents/plans/active/2026-04-25_surface-multiscale-trace.md`
- Tree primitives: `src/lrg_eegfc/utils/metrics/tree.py`,
  `src/lrg_eegfc/utils/metrics/tree_distance.py`
- Hypothesis-test helpers (use, do not reinvent):
  `src/lrg_eegfc/utils/metrics/hypothesis.py`
- H2-method overview: `.agents/guides/02_methods/h2-metrics.md`
- ImCoh method note: `.agents/guides/02_methods/imcoh-guide.md`
- Probe-bias controls: `.agents/guides/02_methods/probe-bias-guide.md`
