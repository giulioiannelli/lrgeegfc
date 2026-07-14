---
name: localization-cohort-statistic-audit
type: scope
era: IMCOH_ABS_COHORT_N10 (mst@0.20 recovery)
status: current
created: 2026-07-13
updated: 2026-07-13
status_note: C1–C5 DONE + verified (frac=0.10 sweep the only pending, non-decisive). Verdict below + in the report.
pointers:
  - .agents/reports/2026-07-13_CHECKPOINT-localization-scale.md          # the "delocalized" verdict this audits
  - scripts/01_compute/sparsified_arc/17_localization_arc_mst020.py      # the broken-statistic run
  - scripts/01_compute/audit/audit_83_localization_matched_strength.py   # the liberal-statistic (whole-graph) run
  - scripts/01_compute/audit/audit_162_localization_robustness_all_systems_rhosym.py  # whole-graph LOO closure
  - .agents/guides/task-persistence-investigation/2026-06-05_localization-atlas.md    # the original localization measure
---

# Localization cohort-statistic audit + controls plan (mst@0.20)

## Head

The mst@0.20 "localization is dead / β is delocalized" verdict rests on a **broken
cohort statistic** — a per-patient Wilcoxon signed-rank across the 5–8 patients that
sample each region, whose one-sided floor at n=5 is p=0.031, so a region sampled in
5 patients (OFC, insula, MTL) **cannot clear BH regardless of effect size**. But the
old whole-graph "β→OFC" it overturned rests on the *opposite* broken statistic — a
pooled cohort-median-vs-R-surrogate-medians test whose p-floor is decoupled from
patient count (1/1001), has **no cross-patient-consistency requirement**, and is
driven by 1–2 strong patients. **Neither verdict is trustworthy.** This scope defines
a localization test that is simultaneously *powered, consistency-requiring,
whole-grid-multiplicity-honest, LOO-robust, coverage-fair, and τ-swept*, and runs it
as the single arbiter. Honest prior (from a zero-recompute recombination of the
existing surrogate p-values): the arbiter will likely confirm a **weak, distributed
limbic–prefrontal concentration** (low_γ/PFC/MTL more than β/OFC) that does **not**
survive honest whole-grid multiplicity — i.e. genuine delocalization, but established
by a valid statistic rather than a DOA one.

## RESULTS (C1 bake-off + C3 duration + C4 raw-vs-coph — verified 2026-07-13)

**The wash-out was PART test-artifact, PART real — and correcting the artifact does
NOT rescue β.** The Wilcoxon in scripts 16/17 was genuinely broken (DOA at K=5), so its
"nothing clears BH anywhere" was under-powered. But under the powered, consistency-
favouring statistic (sign-consistency, whole-grid BH per band, LOO-robust), **every
localization that survives is in `low_gamma` — zero in β/α/δ.** Verified: independent
recompute matches the library; low_γ has 66/432 (15%) system cells beating the *entire*
matched-strength null vs 8–10/432 (~2%, pure noise floor) for β/α/δ.

- **β localization is GENUINELY DELOCALIZED (not a test artifact).** Trace/encoding/
  inference → OFC / cingulate / any system: dead under all four statistics + LOO on the
  non-degenerate backbone (β at the 2% noise floor). The pre-registered **β→OFC,
  encoding→OFC, inference→cingulate all fail.** The whole-graph β→OFC was itself
  LOO-fragile (audit_162), so nothing is lost — the delocalization verdict now rests on a
  *valid powered* statistic, not a DOA one.
- **`low_gamma` beats the null but is DIFFUSE, not focal** — 8/9 systems positive, top-ranked
  is occipital (K=3, coverage artifact); no coherent anatomical story. So low_γ is not
  delocalized like β, but gives **no focal cognitive home** either. The ONE specific low_γ
  concentration is **low_γ → SOZ** (disease marker, q≈.04; β→SOZ does not) — the epilepsy-line
  AUC .82 signature, a binary tissue split, not an anatomical system.
- **Statistic disagreement is diagnostic** (validates the 4-way design): Stouffer
  (omnibus) fires for β-PFC (single-patient-driven) while the median-based pooled +
  sign-consistency + LOO do not → those β "hits" are artifacts. The trustworthy
  (median-based) statistics agree, and they only clear in low_γ.
- **C3 (duration, NON-DESTRUCTIVE ratio regression):** every inference/encoding
  localization cell is duration-clean (|ρ(dur_ratio, loc)| ns, p≥0.25). No truncation used.
- **C4 (raw-vs-coph):** the cophenetic localization value-add is REAL but exists **only in
  low_γ** — 16 cells where coph concentrates OFC/PFC/cingulate/sensorimotor/limbic while
  raw edges are placeless; plus 26 raw-only cells (raw localizes elsewhere). Zero β
  value-add (β delocalized in both representations).

**Bottom line for the user's question:** yes, the localization tests had a real artifact
(Wilcoxon-at-n=5 hid the low_γ signal), but fixing it does not bring β→OFC back — β is
genuinely placeless on the non-degenerate backbone; the only anatomically-concentrated,
statistically-robust localization is `low_gamma` (disease band), distributed over
PFC/limbic/MTL/sensorimotor and the SOZ.

- **C5 (density):** pruning RULED OUT — at frac=0.50 (denser, toward the dense graph)
  β still has **0** LOO-robust survivors; low_γ strengthens (56 cells). β delocalization is
  density-robust. Full verdict: `.agents/reports/2026-07-13_localization-cohort-statistic-verdict.md`.

## 5-point critical preamble (mandatory, before any code)

1. **Claim under test.** On the non-degenerate mst@0.20 backbone, the per-pair
   cophenetic trace/encoding/inference concordance *concentrates* on specific
   anatomical systems (predicted: β→OFC, encoding→OFC, inference→cingulate) above a
   strength-preserving null, at some scale.
2. **Null.** Matched-strength (4-cycle ±δ, R≥200, then sparsify@0.20 → LRG) — the
   locked mandatory FC null. Preserves every node's total coupling; any concentration
   that is a pure function of node strength is reproduced.
3. **Strongest alternative the null must beat, by mechanism.** (a) *Node strength* —
   controlled by matched-strength (verified near-identity by audit_84: within-patient
   ρ(strength, concordance)=+0.16, residualization ≈ identity). (b) *A statistic
   artifact masquerading as biology* — the real adversary here, and the null does
   **not** control for it: the cohort combiner choice alone flips the verdict (30 vs 0
   BH cells on identical data). (c) *Multiple comparisons over 3 targets × 4 bands ×
   16 scales × ~9 systems = 1728 cells* — controlled only by whole-grid BH, never by
   the per-cell BH the shipped scripts use.
4. **What the null CANNOT reject.** Matched-strength says nothing about the cohort
   aggregation. A region can be strength-null yet still fail to localize because (i)
   the cohort statistic is DOA at its coverage (Wilcoxon@n=5), or over-fire because
   (ii) the statistic is liberal (pooled-median, 1-patient-driven). The **arbiter
   statistic**, not the null, is the object under audit here.
5. **Falsification.** If, under the arbiter statistic (powered + consistency-gated) at
   supersystem granularity (K=10, so no coverage floor) with whole-grid BH and LOO,
   *no* system/scale concentrates for any target, the trace is **genuinely
   delocalized** and the paper's localization leg is dropped (not merely "underpowered").
   If limbic/OFC/cingulate concentrate and survive LOO + whole-grid BH + (for
   inference) duration-ratio robustness, localization is *reinstated as a scale-resolved result*.

## Notation

- Patients `p ∈ P` (|P| = 10, COHORT n=10). Systems/units `u ∈ U` at granularity
  `g ∈ {system, supersystem, soz}`. `K_u` = # patients sampling `u`.
- Backbone `B_frac(W) = MST(W) ∪ top-frac edges`, default `frac = 0.20`.
- Dimensionless scale `s = τ·λ_max ∈ SGRID = logspace(0, log₁₀180, 16)`.
- Cophenetic condensed vector `c^{ph}_s = cophenetic_at_scale(eig(B_0.20(W^{ph})), s)`.
- Per-pair concordance targets (rank-based; audit_151/158/160 primitives):
  - trace `σ^T = ½[conc(c^{tt}−c^{A}, c^{post}−c^{B}) + conc(c^{tt}−c^{B}, c^{post}−c^{A})]`
  - encoding `σ^E` = same with `c^{learn}` in place of `c^{tt}`
  - inference `σ^I` = arm-symmetrized first-order **partial** concordance
    `conc_partial(c^{tt}−c^{learn}, c^{post}−c^{B} | c^{learn}−c^{A})`
- Per-unit observed score `m_u(p) = mean_{pairs incident to u} demean(σ)` (endpoint
  incidence, demeaned over the patient's kept pairs) — `node_incidence_mean(...,
  demean=True)` aggregated to `u`.
- Surrogate array `M_u(p) = {m_u(p; r)}_{r=1..R}` from the R matched-strength draws.
- Per-patient MS evidence: upper-tail `π_u(p) = (1 + #{r: M_u(p;r) ≥ m_u(p)})/(R+1)`
  and z-score `z_u(p) = (m_u(p) − mean_r M_u(p;r)) / sd_r M_u(p;r)`.

## Definitions — the four cohort statistics (all on IDENTICAL `m_u(p)`, `M_u(p)`)

For each `(target, band, scale, granularity, u)`:

- **(a) POOLED-MEDIAN (liberal reference; = audit_83/151/158/160/171).**
  `M_obs = median_p m_u(p)`; per realization `M_surr[r] = median_p M_u(p;r)`;
  `p_pool = (1 + #{r: M_surr[r] ≥ M_obs})/(R+1)`. Floor `1/(R+1)`, independent of `K_u`.
  *Failure mode: 1–2 strong patients drive the cohort median; no consistency required.*

- **(b) WILCOXON (DOA reference; = scripts 16/17).**
  `p_wilcox = wilcoxon(m_u(p) − median_r M_u(p;r), alternative='greater')`. Floor
  `1/2^{K_u}`; at K=5, `p_min=0.031` → BH-DOA. *Failure mode: structural under-detection
  at the anatomy of interest.*

- **(c) RANDOM-EFFECTS COMBINE (powered; omnibus).** Stouffer on `z_u(p)` with a
  between-patient variance inflation (or Lancaster/`combine_pvalues`), giving `p_re`.
  *Failure mode: omnibus — tests "signal somewhere in the cohort", 1-patient-drivable;
  use only cross-checked against (d).*

- **(d) SIGN-CONSISTENCY PERMUTATION (RECOMMENDED PRIMARY).** Cohort statistic
  `Θ_u = median_p z_u(p)` **gated by** a consistency requirement
  `frac_pos = #{p: π_u(p) < 0.5}/K_u`. Null distribution of `Θ_u` built from the
  matched-strength realizations themselves — `Θ_u^{null}[r] = median_p ẑ_u(p; r)`
  where `ẑ_u(p;r)` is realization `r` scored against the *rest* of `p`'s surrogate
  ensemble (leave-one-realization-out z), giving a proper R-resolution null of the
  *median-of-z* statistic. `p_sc = (1 + #{r: Θ^{null}[r] ≥ Θ_u})/(R+1)`, **reported
  only if `frac_pos ≥ ⅔`** (else marked "inconsistent"). This is powered (R-resolution,
  works at K=5) AND consistency-requiring (kills the 1-patient-driven liberal failure).

**Adjudication.** A cell is a *localization hit* iff it clears **(d)** under
**whole-grid BH** (family = all `target × band × scale × system` cells) AND is
**LOO-robust** (worst single-patient drop still `q < 0.05`) AND `frac_pos ≥ ⅔`.
(a)/(b)/(c) are reported alongside purely as sensitivity references, never as the gate.

## Properties

- **Range/identifiability.** `p_sc ∈ [1/(R+1), 1]`; achievable floor is R-limited (not
  K-limited), so OFC (K=5) is testable — the fix to (b). Consistency gate makes it
  non-liberal — the fix to (a).
- **Coverage-fairness.** At `g = supersystem`, `limbic` = OFC∪cingulate∪insula∪MTL is
  sampled in **K=10** (all patients) → Wilcoxon itself becomes viable (floor 0.001) and
  every statistic agrees; this is the load-bearing coverage-robust anchor.
- **What it CANNOT detect.** It cannot distinguish a genuine single-region hotspot from
  a distributed limbic territory better than the atlas granularity allows; and it cannot
  rescue a target whose per-pair concordance is truly flat across systems (a real
  delocalization reads as `frac_pos ≈ 0.5`, `p_sc ≈ 1`).
- **Complexity.** Same as script 17 (one eig per phase per realization, cophenetic at
  16 scales) plus O(R) per-cell for the LOO-realization z — dominated by the existing
  eigendecompositions. Must retain the per-patient surrogate unit-mean arrays to disk.

## Caveats & failure modes (each with mitigation)

- **Whole-grid multiplicity is the real gate.** Per-cell BH over ~9 systems is not
  honest when a "best scale" is then read off 16 scales. → Report both per-cell and
  whole-grid BH; the *verdict* uses whole-grid.
- **Coarse-scale reading.** Best cells often land at s=45–180. Cophenetic degeneracy is
  *ruled out* as a dilution mechanism (tie structure flat across scales; costs <1.2%
  rank variance; jitter-stable). But coarse scales are still where the global trace
  "collapses" (τ-sensitivity guide). → Read the full 16-scale curve; do not let a coarse
  best-scale be the headline; prefer fine→meso where the trace itself is τ-robust.
- **Inference duration confound.** Handled by the project-locked **duration-ratio
  regression** (non-destructive), NOT by truncation. The truncation / "length-matched"
  null (audit_113b/161, head-truncate task_test → task_learn) is **RETIRED-INVALID**
  (δ neg-control false-positives; a contiguous window injects a common-mode shift into
  `D_TT`). **No series is ever shortened.** Encoding is length-symmetric by construction
  (`e = c_learn − c_pre`, no test/learn asymmetry). → C3 = per-patient duration-ratio
  *sensitivity* of the inference localization; a length artifact would correlate with the
  task_test/task_learn duration ratio.
- **Raw-vs-coph.** "Localization is the cophenetic value-add" (audit_171: coph→cingulate,
  raw placeless) is the load-bearing claim and was pooled-median-based on the dense
  graph. → C4 re-establishes it on the backbone with statistic (d).
- **Backbone pruning.** frac=0.20 keeps only the strongest 20% edges; a weak long-range
  fronto-limbic substrate could be pruned. → C5 sweeps frac ∈ {0.10, 0.20, 0.50}.
- **Combiner independence.** (c)/(d) assume per-patient independence within a cell
  (plausible) but the SAME 200 surrogates recur across scales/targets (grid internally
  correlated). → whole-grid BH is conservative w.r.t. this; state it.

## Pseudocode (arbiter run, C1)

```
for target in {trace, encoding, inference}:
  for band in {delta, alpha, beta, low_gamma}:
    for p in P:
      for ph in phases(target):            # A,B,(learn),tt,post
        eig_obs[ph]      = eig(B_0.20(W[p,ph,band]))
        for r in 1..R:  eig_sur[ph][r] = eig(B_0.20(matched_strength(W[p,ph,band], rng)))
      for s in SGRID:
        c_obs   = {ph: cophenetic_at_scale(eig_obs[ph], s)}
        sigma   = concordance_target(target, c_obs)          # per-pair
        for g in {system, supersystem, soz}:
          m[p,s,g,:]   = unit_means(demean(sigma), incidence, labels[g])   # per unit
        for r in 1..R:
          c_r   = {ph: cophenetic_at_scale(eig_sur[ph][r], s)}
          M[p,s,g,:,r] = unit_means(demean(concordance_target(target,c_r)), ...)
      DUMP m[p], M[p] to disk                                # retain arrays
    for s, g, u with K_u>=5 (system) or K_u==10 (supersystem):
      z[p]      = (m[p]-mean_r M[p]) / sd_r M[p]
      pi[p]     = (1+#{M[p,r]>=m[p]})/(R+1)
      frac_pos  = #{pi[p]<0.5}/K_u
      Theta     = median_p z[p]
      Theta_null[r] = median_p LOO_realization_z(m,M,p,r)    # R-resolution null of median-z
      p_sc      = (1+#{Theta_null>=Theta})/(R+1)
      p_pool, p_wilcox, p_re = (a),(b),(c)                   # sensitivity refs
      loo_worst = max over p_drop of p_sc(P\{p_drop})
    q_grid_sc   = BH(all p_sc over target x band x scale x system)      # HONEST family
    q_cell_sc   = BH(p_sc within each (target,band,scale) over systems) # reference
  HIT iff q_grid_sc<0.05 AND loo_worst<0.05 AND frac_pos>=2/3
```

## Visualization spec

- **Per-target localization heatmap:** rows = systems (ordered limbic→neocortex),
  cols = 16 scales; cell color = signed `Θ_u` (diverging, band palette-neutral), a
  ring/marker on cells that are HITs (whole-grid q<0.05 ∧ LOO ∧ consistent). One panel
  per band; reading rule: *a horizontal band of same-color cells across scales = a
  scale-invariant anatomical concentration; a lone marked cell = a scale-specific,
  multiplicity-fragile hint.*
- **Statistic bake-off panel:** for the pre-registered targets (β/OFC, enc/OFC,
  inf/cingulate) and the limbic supersystem, a small-multiple of `p_pool` (liberal),
  `p_wilcox` (DOA), `p_sc` (arbiter) vs scale, with the whole-grid-q<0.05 line — makes
  the "which statistic you pick decides the verdict" explicit and honest.
- **Coverage bar:** `K_u` per system with the achievable-floor annotation, so the
  reader sees why n=5 systems were untestable under (b).

## Connection to prior tools

| Prior tool | Substrate | Cohort statistic | This audit's relation |
|---|---|---|---|
| audit_83/151 (trace) | dense | pooled-median R200 | (a) reference; liberal |
| audit_158 (enc→OFC), 160 (inf→cing) | dense | pooled-median R1000 | (a) reference; targets to re-test |
| audit_161 (inf truncation null) | dense | pooled-median R200 | **RETIRED-INVALID** (δ false-pos) — NOT ported; use ratio regression |
| audit_162 (all-systems LOO) | dense | pooled-median R1000 + LOO | fresh recompute = *nothing robust*; motivates LOO gate |
| audit_171 (raw-vs-coph) | dense | pooled-median R1000 | → **C4** port to backbone |
| audit_157/173 (scale sweep) | dense | Wilcoxon (global only) | no per-scale localization exists on the dense graph either |
| scripts 16/17 | mst@0.20 | Wilcoxon | (b) reference; DOA — the thing being fixed |

This audit **subsumes** the cohort-statistic choice of all of the above under one
principled test; it does **not** replace the per-pair concordance primitives (reused
verbatim) nor the matched-strength null (unchanged).

## Implementation plan (controls, prioritized)

- **C1 — arbiter bake-off (the linchpin).** New script
  `scripts/01_compute/sparsified_arc/23_localization_statistic_bakeoff_mst020.py`.
  Extends script 17's compute to retain per-patient surrogate arrays and emit all four
  statistics + whole-grid BH + LOO + `frac_pos`, at 3 granularities, all 16 scales, all
  3 targets, 4 bands. Output `data/sparsified_arc/localization_bakeoff_mst020/`.
  Library: reuse `utils.metrics.node_localization` (unit_means, rank_concordance,
  concordance_partial), `utils.fc.backbone`, `utils.fc.heat_multiscale`,
  `utils.metrics.surrogate`, `utils.metrics.hypothesis.bh_fdr`. Promote the four
  cohort-combiners to `utils/metrics/cohort_localization.py` (≥2 callers: C1 + C3/C4).
  **Optimize+time first** (numba surrogate already warm; time 2 cells, extrapolate,
  print `[i/N] ETA`). Est. ≈ script 17 runtime (~25–30 min) + disk for arrays.
- **C2 — coverage-fair anchor.** Folded into C1 as `g=supersystem` (limbic K=10). The
  cleanest verdict: does `limbic` concentrate under ALL four statistics at K=10?
- **C3 — duration-ratio sensitivity of inference localization (NON-DESTRUCTIVE, no
  truncation).** Locked: the truncation null is retired-invalid; **no series is
  shortened.** Regress per-patient inference→unit localization on the task_test/task_learn
  duration ratio across patients — a length artifact correlates, a real inference trace
  does not. Encoding (length-symmetric by construction) is the clean co-anchor.
- **C4 — raw-vs-coph value-add on backbone.** Port audit_171; does coph localize where
  raw is placeless, on the non-degenerate substrate, under statistic (d)?
- **C5 — backbone-density robustness.** Rerun C1's pre-registered cells at
  frac ∈ {0.10, 0.20, 0.50}; tests the "0.20 pruned the substrate" alternative.

## Open questions (deferred, load-bearing)

1. **Consistency threshold** `frac_pos ≥ ⅔` vs `≥ ½+1` — default ⅔; revisit if it
   kills a limbic K=10 signal that (a)/(c) both show.
2. **Random-effects vs fixed** for statistic (c) — default Stouffer + between-patient
   variance; a full mixed-effects model is overkill at K≤10.
3. **Whole-grid family definition** — per-band (16×9=144) vs full-grid (3×4×16×9=1728).
   Default: report both; verdict uses per-band whole-grid (the paper reads one band at
   a time), with the full-grid count stated as the conservative bound.
4. **If the arbiter confirms delocalization** — the paper's localization leg is
   replaced by the honest "trace is strong but placeless / distributed limbic–prefrontal"
   framing already supported by `node_ranksize_delocalisation` (β ~88% co-move).
