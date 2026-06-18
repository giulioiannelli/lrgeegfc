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
| 2026-06-08/11 | **ARCHIVED** — CRC + coordinated cross-phase null (SB-CRC) | **archived 2026-06-12** | **Archived as opaque matrix-level nulls — see [`scripts/archive/2026-06_opaque-matrix-nulls/POSTMORTEM.md`](../../../scripts/archive/2026-06_opaque-matrix-nulls/POSTMORTEM.md).** Two coherency-matrix Haar-rotation surrogates built after matched-strength: CRC (on-manifold bracket; finding: matched-strength was adequate on the realizability axis, 30/30 cells) and SB-CRC (coordinated null to separate task-specific trace from a stable fingerprint). Retired at user direction: opaque (not one-sentence explainable), construction-dependent (verdict swung with backbone choice), and run onto real-data verdicts before synthetic validation. **Deeper reason:** separating task-specific persistence from a trait-level backbone has *no clean null — physical or matrix — without a task-free control session*; phase-randomized surrogates don't help (multivariate PR is a no-op for `\|ImCoh\|`, univariate destroys it). Keep matched-strength; state the control-session limitation in one sentence in the manuscript. Scope reports + scripts + data under the respective `archive/2026-06*` dirs. |
| 2026-06-05 | [Sampling-conditioned localization atlas](2026-06-05_localization-atlas.md) | **current** | **DONE — β trace localizes to MTL+OFC, epi-independent (un-buries the 2026-05-30 over-retraction).** Descriptive re-framing of trace localization. Drops the cohort permutation GATE (a structural false-negative at n=10 / ≤5-per-region) for sampling-conditioned positive rates `J⁺_g/K_g` = patients-where-trace-concentrates / patients-IMPLANTED — so Hip-β reads **4/5**, not 4/10. Four granularities (DK region → system MTL/limbic/… → lobe → hemisphere; pooling tests "different region per patient, same system"), both matched-strength TRACE probes (cophenet + cross-phase Grassmann, NOT the phase-avg anchor), epi-incl/excl split, depth-shaft adjacency control shown ALONGSIDE with a `deep_target` flag (electrode = structure for Hip/Amy, not a confound). No pre-registered gate; compute→plot→adjudicate. `audit_81_localization_atlas.py`. |
| 2026-06-05 | [Occult (non-labeled) epi-node marker](2026-06-05_occult-epi-node-marker.md) | **built-negative** | **Discovery reframe of the per-node marker.** Not "classify KNOWN epi" but "do NON-labeled contacts behave like epi beyond strength+shaft — occult candidates?". New per-node leave-node-out (LNO) influence `Δρ_i = ρ_full − ρ_{∖i}` on the cohort cophenetic trace, residualized WITHIN patient against strength + same-shaft-epi + along-shaft pos; bands α/β/γ_l (`audit_82`). **Result:** at β `Δρ` is a real strength-INDEPENDENT within-patient epi signal (AUC 0.63, 8/9, p=0.010, survives rank-partial) that the row-restricted F1 missed — but cross-patient LOPO is at chance (β PR 0.12 vs prev 0.11, lift 1.2×), the same ceiling as `audit_80`. No transferable occult marker; 154 candidates hypothesis-only (no surgical-outcome ground truth). |
| 2026-06-05 | [Per-node epi marker from trace behaviour](2026-06-05_epi-node-trace-marker.md) | draft | **Direction-2 marker scope (gated).** Lifts per-pair ρ_split (audit_76) + Grassmann mode-loading (audit_66) to a per-NODE feature vector, surrogate-z'd so it is strength-independent by construction, to flag epileptogenic contacts. Real only if it beats BOTH a node-strength and a probe-adjacency baseline in leave-one-patient-out (per-leaf coph ρ is global-position/strength dominated; epi clusters on shafts). Split into descriptive `audit_79` (F1 ρ_split_node, F2 coherence, F7 mode-load; this pass) + deferred classifier `audit_80` (LOPO + whole-shaft masking + candidate discovery). Gated on Direction-1 (audit_77/78) showing a per-band epi signature. |
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
| 2026-05-08 | [Epi cross-phase rigidity (TARR on E_p)](2026-05-08_epi-cross-phase-rigidity.md) | draft | **LRG-epilepsy Direction A.** Apply the trace / anchor / reset / rearrange (TARR) classification to the *fixed* leaf set E_p (epileptic node set per patient). Three KC distances on the induced subtree at λ ∈ {0, 0.5, 1}, classify per (patient, band, λ), z-score against a size-matched random-leaf-subset null (R=100). Cohort-positive at ≥ 8/9 patients with z>1.96 in any class × band × λ cell, BH-FDR q<0.05 cohort-Wilcoxon. Cross-probe variant `E_p^cp` reported alongside. Direct test of the textbook "epileptic networks are pathologically rigid across cognitive states" intuition. Either result publishable. `audit_54_epi_rigidity_compute.py` planned (slots 51-53 already taken). |
| 2026-05-08 | [Eigenmode-embedding motion](2026-05-08_eigenmode_embedding.md) | current | **Per-contact decomposition of the §5.4 Grassmann.** Procrustes-aligned eigenmode embedding `Y^φ(k) ∈ R^{N × k}` (columns `v_2 … v_{k+1}` of L̂); per-contact phase motion `Δ_i = ‖Y_i^a − Ỹ_i^b‖_2` after `R = argmin ‖Y^a R - Y^b‖_F` Procrustes; per-contact trace `T_i^Y = Δ^{tt→post} − Δ^{pre→tt}` mirrors §5 Result 2. Implemented n=10 × β at k=3. Pat_02 β shows TRACE-shape on probe Q (Q6/Q2/Q3): Δ_pre→tt ≈ 0.65–0.86 with Δ_pre→post ≈ 0.02–0.05 — the contacts move in task and snap back. Decomposes the global Grassmann scalar across rows; reveals localisation invisible to the cumulative cell. Multiscale extension k ∈ {2, 5, 10, 20} + universality / implant-bias discount deferred. `audit_59_eigenmode_embedding.py`. |
| 2026-05-08 | [Epi eigenmode localization (m^E_k on L̂)](2026-05-08_epi-eigenmode-localization.md) | draft | **LRG-epilepsy Direction C.** Mode projection mass `m^E_k = Σ_{i ∈ E_p} v_k(i)²` of L̂ eigenvectors restricted to E_p, binned on normalised eigenvalue λ̂ = λ/λ_max. Strength-stratified null (Q=5 quintiles, R=100 draws) corrects for heavy-tail-degree triviality. C1 (static cohort signature): does a cohort-consistent (≥ 8/9) localisation regime of width ≥ 3 contiguous bins survive BH-FDR? C2 (phase-trace extension): does the localisation regime shift across phases? Anderson-localisation analogy framed loosely. **Blocked** on the eigenvector cache extension from `2026-04-29_eigenvector-direct-pivot-plan.md` step 0. `audit_57_epi_eigenmode_compute.py` planned (slot 55 taken). |
| 2026-05-08 | [Band-agnostic LRG (band-pool ρ_split proxy)](2026-05-08_band_agnostic_lrg.md) | current | **Does the §5 cohort task-trace direction need band splitting?** Proxy-tier answer: NO. Stack the per-band Run-A Δ_task / Δ_rest pair vectors across all 6 bands per patient → single broadband ρ_split^pool. Cohort: 9/10 ρ^pool > 0, 10/10 ρ^pool > drift floor, paired Wilcoxon p = 0.0025 — direction count *exceeds* every per-band α/β/γ_l count individually (each 7–8/10). Median pool ρ = +0.148 lands between α (+0.115) and β (+0.222); band splitting therefore localizes magnitude, not direction. `audit_58_band_agnostic_lrg.py` implemented 2026-05-08. §11 of the scope specifies the canonical broadband |ImCoh| → LRG pipeline (≤ 30 min, parked) for the magnitude-tier confirmation. Pat_15 only proxy-negative; Pat_03 1024 Hz survives at +0.415. |
| 2026-05-08 | [Epi Grassmann embedding (E.align + E.resect)](2026-05-08_epi-grassmann-embedding.md) | draft | **LRG-epilepsy Direction E (user-flagged).** Two well-defined Grassmann scalars on the LRG eigenstructure of the FC graph. `ε^E_align = sqrt(min(k,|E_p|) − ‖V_k[E_p,:]‖_F²)` — chordal distance from the top-k eigenspace to the epi coordinate subspace `S_E ⊆ ℝ^N`; equivalent to the average fraction of top-k mode mass sitting on epi rows (`f^E_k = ‖V_k[E_p,:]‖_F²/k`). E.align is the cumulative-over-modes form of Direction C's per-mode `m^E_k`. `δ^E_resect = d_chord(Q_full[N_p], V_k^resect)` — chordal distance between the full-graph eigenspace restricted to non-epi rows (orthonormalised via QR) and the eigenspace freshly computed on the resected FC graph; this question is not asked by any other measure. Strength-stratified null (Q=5 quintiles, R=100 draws) shared with C. **No pre-registered acceptance gate** — we report the per-patient observation, the null reference, and the cohort distribution (median/IQR/per-patient) per (band, k, phase), then look at the plots and decide what's signal vs noise. `audit_61_epi_grassmann_compute.py` planned (slots 55–60 all taken). Promotes `chordal_distance`, `principal_angles`, `grassmann_to_coord_subspace`, `chordal_full_vs_resect` to new module `lrg_eegfc.utils.metrics.spectral` (≥ 3 callers across audit_37, audit_46, audit_61). Eigenvectors already cached in `IMCOH_LRG_CACHE`; no cache extension needed. |
| 2026-05-08 | [Anchor anatomy baseline + mutual-exclusivity tie-breaker](2026-05-08_taxonomy_mutual_exclusivity.md) | current | **§5.6 anchor anatomy caveat + canonical class-ordering.** Per-(patient, band) within-module same-probe-pair fraction `f_sp` for each of the four KC classes, compared to the per-patient cohort baseline `f_sp^baseline = Σ_q binom(n_q, 2) / binom(N, 2)` (random-leaf-draw geometry null). At β: anchor enrichment **6.22x** baseline (median 0.524 over n=43), trace 3.96x (median 0.333 over n=25), reset 4.28x, rearrange 2.25x. At γ_l/γ_h anchor median sits at 1.0 (entirely on a single probe, 11.87x baseline). The §5.6 prose framing "anchors are the cross-phase null" needs a probe-geometry caveat — anchors are **not** anatomy-neutral under |ImCoh|. Also formalizes the canonical mutual-exclusivity tie-breaker class-ordering **anchor > trace > reset > rearrange > diffuse** with `MIN_RESIDUAL = 3` floor. `audit_60_anchor_anatomy_baseline.py` + `audit_60b_anchor_anatomy_baseline_figure.py` implemented 2026-05-08. |
| 2026-05-30 | [Asymmetric pair-trace regression slope `s_TR`](2026-05-30_asymmetric-pair-trace-regression.md) | current | **Comparative methodology audit: does the asymmetric regression slope `s_TR = ⟨Δ_task, Δ_rest⟩/‖Δ_task‖²` beat the symmetric `ρ_split^coph`?** Four measures {`ρ_Spearman` (current), `ρ_Pearson`, `s_TR` slope, `R²`} × three primitives {`ρ̂(τ_max)`, `D(τ_max)`, `D_coph`} × 6 bands × n=10 under the shared R=200 matched-strength null (`audit_73`). **Verdict: NO — `s_TR` loses the β trace** (β 7/10 p=0.007 → 2/10 p=0.161 on `D_coph`) because magnitude-weighting reintroduces the strength-driven merge-height structure matched-strength nullifies (Pat_03 obs slope +0.612 vs own surrogate +0.626). Keep Spearman; `ρ_Pearson` is the defensible magnitude-aware companion (β 8/10 p=0.005); `R²` fires on all 6 bands incl. anti-trace θ (unusable). NAMING: slope is `s_TR`/`slope`, never `β` (reserved for the band). Reproduces audit_63 (Δ=4.6e-07) + preprint_05 (bit-exact). Report: `.agents/reports/2026-05-30_pair-trace-measures-methodology-audit.md`. |

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
