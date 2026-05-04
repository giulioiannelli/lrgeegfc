---
name: decision-rules
type: scope
era: COHORT_N10
status: draft
created: 2026-04-29
updated: 2026-04-29
pointers:
  - .agents/guides/task-persistence-investigation/README.md
  - .agents/guides/task-persistence-investigation/2026-04-29_measure-correctness-audit.md
  - .agents/reports/2026-04-27_task-persistence-reconciliation.md
---

# Decision rules — pre-registered triangulation logic for the cohort-coverage matrix

**Pre-registered, per-rung positivity / negativity / silent /
ineligible thresholds, dissenter-handling protocol, and triangulation
logic for the cohort-coverage matrix that is the rebuild's central
deliverable. These rules MUST be locked before any matrix-cell
evaluation runs (rebuild step 3 onwards) — once committed, no
threshold or logic in this document is tunable post-hoc on the basis
of what the matrix turns out to say. The rules apply across all 7
geometric rungs (L1–L7); for each rung they specify (a) the per-cell
positivity test, (b) the cohort threshold, (c) the band-level
verdict that aggregates per-cell results, and (d) how the
band-verdict feeds the cross-rung triangulation predicate that
delivers the headline claim.**

---

## Notation

- `P` — cohort, `|P| = 10` patients per `data/audit/cohort_metadata.csv`.
- `B = {δ, θ, α, β, γ_l, γ_h}` — six bands.
- `r ∈ R = {L1, L3, L4, L5(k), L5(h_rel), L6, L7}` — seven rungs of the
  geometric ladder per the rebuild plan §2. L0 closed by
  `2026-04-28_residual-subspace-diagnostic.md`. **L2 (entropy curve /
  specific heat) is permanently removed** — `S(τ)` and `C(τ)` are
  non-informative for our fully-connected weight-heterogeneous outlier
  case (continuous spectrum, no discrete `C(τ)` peaks; see
  `.agents/guides/02_methods/lrg-framework-guide.md` §6 + memory entry
  `lrg_outlier_case_fully_connected.md`).
- For each `(p, b, r)` cell: `Q_r(p, b)` is the per-cell statistic
  (rho, Δ_VI, T-module count, etc.) and `Q_r^{null}(p, b)` is the
  rung-matched within-baseline-fluctuation null (drift floor /
  half-baseline ρ / J_{rpre_A,rpre_B} / etc.).
- For rungs with a scale axis (L5(k), L5(h_rel), L7): `s` indexes the
  scale value (integer-k or h_rel bin or partition-k); cells become
  `(p, b, r, s)`.
- For per-leaf rungs (L4, L6, L7): the per-cell value is a *count of
  T-leaves* `n_T(p, b, k)` rather than a single scalar.

## Definitions

### Per-cell positivity `π(p, b, r [, s])`

Boolean. `1` iff the per-cell statistic exceeds its rung-matched null
in the trace direction, AND (for rungs with continuous statistics) is
strictly greater than zero.

| Rung | Statistic | Null | Positivity test |
|---|---|---|---|
| **L1** | per-patient Spearman ρ on `(Δ_task, Δ_rest)` upper-triangle, split-baseline | drift_floor `ρ_null_drift = Spearman(D^post_B − D^post_A, D^pre_B − D^pre_A)` | `ρ_split > 0 AND ρ_split > ρ_null_drift` |
| **L3** | per-patient `Δ_KC(λ) = d_KC(λ; T_test, T_post) − d_KC(λ; T_test, T_pre)` for `λ ∈ {0, 0.1, …, 1.0}` (test-tree closer to post than to pre = trace) | `d_KC(λ; T_rpre_A, T_rpre_B)` cohort distribution | per λ: `Δ_KC(λ) < 0` AND `|Δ_KC| > median(half-baseline KC)` |
| **L4** | per-patient T-module count at strict J_min=0.9 over `T_test ∪ T_learn` | half-baseline T-count using `Z_rpre_A` vs `Z_rpre_B` as the (pre, target) pair | `n_T(p, b) > n_T_null(p, b)` AND `n_T ≥ 1` |
| **L5(k)** | per-(p, b, k) `Δ_VI(k) = VI(c_pre, c_post) − VI(c_test, c_post)` | `drift_dVI(k) = VI(c_rpre_A, c_rpost_A) − VI(c_rpre_B, c_rpost_A)` | `Δ_VI(k) > 0 AND Δ_VI(k) > drift_dVI(k)` |
| **L5(h_rel)** | per-(p, b, h_rel) `Δ_VI(h_rel)` via `fcluster_at_h_rel` | `drift_dVI(h_rel)` analogue from halves | `Δ_VI(h_rel) > 0 AND Δ_VI(h_rel) > drift_dVI(h_rel)` AND cell is in the eligible h_rel zone (see §3) |
| **L6** | per-leaf trace-affinity excess `a_trace(ℓ) − a_trace_null(ℓ)` | half-baseline `ρ_{rpre_A, rpre_B}(ℓ)` corner-affinity | per-patient: `n_T-leaves(p, b) ≥ τ_leafcount(b)` where `τ_leafcount(b) = 95th percentile of n_T-leaves under half-baseline cohort distribution` |
| **L7** | per-(leaf, scale) `J_{tt, rpost}_k(ℓ) − J_{rpre_A, rpre_B}_k(ℓ)` corner affinity classified as TRACE | half-baseline `J_{rpre_A, rpre_B}_k(ℓ)` | per-(p, k): `n_T-leaves(p, b, k) ≥ τ_leafcount(b, k)` calibrated from half-baseline 95th percentile |

### Per-band cohort verdict `V(b, r)`

Computed by aggregating `π(p, b, r [, s])` across `p ∈ P`. For rungs
with a scale axis, also requires a contiguous-ridge condition.

| Aggregation | Definition |
|---|---|
| `n_pos(b, r [, s])` | `|{p ∈ P : π(p, b, r [, s]) = 1}|` (numerator) |
| `n_eligible(b, r [, s])` | `|{p ∈ P : cell is eligible per §3}|` (denominator; ineligible cells removed) |
| `frac_pos(b, r [, s])` | `n_pos / n_eligible` |
| `cohort_threshold` | `frac_pos ≥ 0.8` (= ≥ 8/10 with no ineligibles) |
| `wilcoxon_test` | one-sided paired Wilcoxon `Q_r(full) > Q_r^{null}` per band; FDR-BH `m = 6` bands per rung; require `q ≤ 0.05` |
| `ridge_length` | for L5(k), L5(h_rel), L7: longest contiguous run of `frac_pos ≥ 0.8` along the scale axis. Require `L ≥ 3` |
| **V(b, r)** ∈ `{positive, negative, silent, ineligible}` | as in §2.3 below |

```
V(b, r) =
    positive    if cohort_threshold ∧ wilcoxon_test ∧ (no scale axis OR ridge_length ≥ 3)
    negative    if frac_pos ≤ 0.2 across all eligible cells (anti-trace)
    silent      if 0.2 < frac_pos < 0.8 (under-powered or null)
    ineligible  if (n_eligible / |P|) < 0.5 (cohort comparison structurally degenerate)
```

### Triangulation `T(b)`

A band's headline label is a function of its `V(b, r)` across rungs.

| Label | Predicate |
|---|---|
| **headline-triangulated** | `V(b, L1) = positive` AND `V(b, L5(k)) = positive` AND (`V(b, L3) = positive at λ ≥ 0.5` OR `V(b, L7) = positive at any (b, k)`) |
| **partition-resolution-locked** | `V(b, L5(k)) = positive` AND `V(b, L5(h_rel)) ≠ positive` AND `V(b, L1) ≠ positive` |
| **continuous-only / non-modular** | `V(b, L1) = positive` AND `V(b, L5(k)) ≠ positive` AND `V(b, L4) ≠ positive` AND `V(b, L7) ≠ positive` |
| **anatomy-suspect** | `V(b, r) = positive for some r ∈ {L5(k), L5(h_rel), L7}` AND `V(b, r') = ineligible for r' in fine-scale rungs` (γ_h with `h_rel_min` up to 0.71) |
| **ergodic** | `V(b, L1) ∈ {silent, negative}` AND `V(b, L5(k)) ∈ {silent, negative}` AND `V(b, L5(h_rel)) ∈ {silent, negative}` AND `V(b, L4) ∈ {silent, negative}` AND `V(b, L7) ∈ {silent, negative}` |
| **uncontrolled / exploratory** | `V(b, r) = positive for some r` AND `V(b, r')` is contested or null is missing for the corresponding control rung |
| **inconsistent** | None of the above match (e.g. positive at L1 + positive at L5(k) but not at any other rung) — flag as a research finding requiring narrative explanation, not a verdict |

## Properties

- **Pre-registration is the load-bearing property.** Once this scope
  report is committed and step 1 of the rebuild proceeds, the
  thresholds (`cohort_threshold = 0.8`, `q ≤ 0.05` BH-FDR, `r_rb ≥ 0.6`,
  `ridge_length ≥ 3`, `eligibility_threshold = 0.5`,
  `τ_leafcount(b) = 95th percentile of half-baseline cohort distribution`)
  are immutable. Adjustments require a new dated scope report
  superseding this one.
- **Conservative aggregation.** A single verdict per (band, rung)
  cell — no scalar averaging across rungs. The triangulation logic
  preserves rung-specific information.
- **No silent dissenter exclusion.** Pat_02 and Pat_03 are *flagged*
  via "×"/"+" superscripts on cells but always counted in the
  denominator. A "robust-w/o-dissenter" sub-claim is reported alongside
  the strict cohort claim, never instead of it.
- **Multiple-comparison handling per rung, not across rungs.** BH-FDR
  applies to the 6 bands within each rung. Cross-rung triangulation
  uses no further correction (the predicates are designed to be
  *more conservative* by requiring conjunction across independent
  geometric reductions).
- **Scale-axis ridge condition is asymmetric.** A single positive
  cell is not enough for L5(k), L5(h_rel), or L7 — ridge length `L ≥ 3`
  is required. This guards against single-cell artefacts (k-collapse
  near singleton/giant zones).
- **Ineligibility is not silence.** A cell where the geometric
  comparison is structurally degenerate (e.g. fixed h_rel below the
  cohort `h_rel_min` for a band) is `ineligible`, not `silent`. The
  band verdict computes only over eligible cells; the ineligible-rate
  feeds the *anatomy-suspect* triangulation predicate.
- **What this does NOT specify.** Per-rung *implementation details*
  (which exact CSV column to use, which RNG seed for cluster-perm) are
  in the per-rung scope reports (`2026-04-29_kc-lambda-sweep.md`,
  etc.). This document specifies only the *decision logic over the
  implementation outputs*.

## Caveats and failure modes

| Caveat | Mitigation |
|---|---|
| **Missing nulls for L3/L6/L7 at time of step 3 v1 matrix.** Per the rebuild plan §10, the matrix v1 is built before §7.2–§7.6 nulls land. | Cells without a rung-matched null are marked `uncontrolled-exploratory`, NOT `positive`. The triangulation predicates explicitly require null-pass; uncontrolled cells contribute to *exploratory* labels, never to *headline-triangulated*. |
| **n=10 is small for BH-FDR with `m=6` bands.** Per-rung Wilcoxon with q ≤ 0.05 corresponds to raw p ≤ 0.0083 in the worst case. | Acceptable. The rebuild's headline rests on multi-rung triangulation, not on single-test FDR survival. The Wilcoxon condition is necessary but not sufficient. |
| **Pat_02 strengthens the cohort signal when dropped.** "Robust-w/o-dissenter" sub-claim might overstate. | Always report the strict cohort claim AS the headline; the "robust-w/o-dissenter" sub-claim appears in supplementary. Never let the dissenter-removed verdict be the headline. |
| **Ridge contiguity at integer-k vs h_rel is not the same metric.** L5(k) ridges in δ are contiguous over k=20–32 (L=13); L5(h_rel) ridges in δ dissolve. | The triangulation predicate explicitly distinguishes these via `partition-resolution-locked` vs `headline-triangulated`. Pre-registered, no mid-stream reframing. |
| **Eligibility threshold `0.5` (cell ineligible if < 50% patients can be evaluated) is arbitrary.** | Sensitivity analysis is owed: re-evaluate the matrix with eligibility ∈ {0.3, 0.5, 0.7} and report stability of band-verdicts in the final reconciliation. Not in this scope. |
| **τ_leafcount calibrated from half-baseline 95th percentile assumes the half-baseline is a valid null.** If the halves cache has structure (e.g. session-order drift exceeds within-session noise), the 95th percentile understates the null. | The L1 H2e split-half audit (`data/reports/imcoh_vi/h2e_split_half_rho_raw.csv`) showed `ρ_null_drift ∈ [-0.04, +0.03]` cohort-wide — drift is statistically indistinguishable from zero. The half-baseline calibration is therefore not contaminated by session drift. |
| **What if `T(b) = inconsistent` for a band?** | Report explicitly as "rung disagreement requiring mechanistic interpretation"; do NOT silently classify as `silent` or `exploratory`. The narrative report names the disagreeing rungs and references the principled-disagreement table (rebuild plan §3) for interpretation. |

## Pseudocode

```
INPUT  : per-rung CSVs from rebuild §6 reuse map (only "current"-verdict measures)
       : cohort_metadata.csv (dissenter flags)
       : per-rung null artefacts (drift floors, half-baselines)
OUTPUT : cohort_coverage_matrix_n10_imcoh_abs.csv
       : per-band triangulation verdict T(b) per §2.4

# Step 1: Per-cell positivity
for each measure m in current_measures:
    rung = m.rung
    for each (p, b [, s]) cell in m.artefact:
        Q     = cell value
        Q_null = lookup_null(rung, p, b [, s])
        eligible = check_eligibility(rung, p, b [, s])
        if not eligible: pi[p, b, r [, s]] = "ineligible"; continue
        passes_null = (Q > Q_null) for trace-direction
        passes_zero = (Q > 0) for continuous statistics
        pi[p, b, r [, s]] = 1 if (passes_null AND passes_zero) else 0

# Step 2: Per-band cohort verdict
for each (b, r):
    cells = pi[:, b, r [, all s]]
    n_eligible = count(cells != "ineligible")
    if n_eligible / |P| < 0.5: V[b, r] = "ineligible"; continue
    frac_pos = count(cells == 1) / n_eligible
    wilcoxon_q = bh_fdr(wilcoxon_one_sided(Q_full > Q_null per b)) for this rung
    if rung in {L5(k), L5(h_rel), L7}: ridge = longest_run(frac_pos >= 0.8 along s); else ridge = inf
    if frac_pos >= 0.8 AND wilcoxon_q[b] <= 0.05 AND ridge >= 3: V[b, r] = "positive"
    elif frac_pos <= 0.2: V[b, r] = "negative"
    else: V[b, r] = "silent"

# Step 3: Cross-rung triangulation
for each band b:
    if V[b, L1]=="positive" AND V[b, L5(k)]=="positive" AND (V[b, L3]=="positive" OR V[b, L7]=="positive"):
        T[b] = "headline-triangulated"
    elif V[b, L5(k)]=="positive" AND V[b, L5(h_rel)] != "positive" AND V[b, L1] != "positive":
        T[b] = "partition-resolution-locked"
    elif V[b, L1]=="positive" AND V[b, L5(k)] != "positive" AND V[b, L4] != "positive" AND V[b, L7] != "positive":
        T[b] = "continuous-only"
    elif V[b, L5(k)]=="positive" AND any(V[b, r] == "ineligible" for r in fine_scale_rungs):
        T[b] = "anatomy-suspect"
    elif all(V[b, r] in {"silent", "negative"} for r in {L1, L5(k), L5(h_rel), L4, L7}):
        T[b] = "ergodic"
    elif any(V[b, r] == "positive" for r) AND not (controls available for r):
        T[b] = "uncontrolled-exploratory"
    else:
        T[b] = "inconsistent"   # narrative explanation required

emit cohort_coverage_matrix CSV
emit triangulation summary {b: T[b]}
```

## Visualisation spec

The cohort-coverage matrix renders as a band × rung heatmap.

- **Rows:** 6 bands (δ, θ, α, β, γ_l, γ_h).
- **Columns:** 7 rungs (L1, L3, L4, L5(k), L5(h_rel), L6, L7). L2 is permanently removed.
- **Cell colour:** `frac_pos(b, r)` on a viridis colormap; `ineligible`
  cells rendered as light grey diagonal hatch (no colour).
- **Cell text:** `n_pos/n_eligible` as a small fraction (e.g. `9/10`).
  Pat_02 / Pat_03 dissenter glyphs ("×"/"+") in the corner if the cell
  is positive only when the dissenter is dropped.
- **Border:**
  - black solid: `V(b, r) = positive`
  - red dashed: `V(b, r) = negative`
  - thin grey: `V(b, r) = silent`
  - none: `V(b, r) = ineligible`
- **Right-side annotation column:** per-band triangulation label `T(b)`
  with colour coding (green=headline-triangulated, yellow=partition-resolution-locked,
  blue=continuous-only, orange=anatomy-suspect, grey=ergodic,
  red=inconsistent).
- **Caption:** explicitly states "cohort consensus is at *nominal* per-rung
  scale resolution, not at fractional dendrogram depth — see L5(h_rel)
  column for the fixed-fractional-depth view." Names the dissenters
  (Pat_02, Pat_03) and the anatomy-flag rule.

PDF only, full vector (no `set_rasterized`). No `fig.suptitle`.
≥3 bands × 3 rungs visible per facet.

## Connection to prior tools

This scope is **operative**, not exploratory: it specifies the rules
that consume the rebuild plan's §6 reuse map and the §7 new measures,
and emits the matrix that supersedes the existing `task_trace_band_k_n10_imcoh_abs.pdf`
headline figure.

| Prior tool | This scope's relationship |
|---|---|
| `2026-04-27_task-persistence-reconciliation.md` §8.8 per-band table | This scope's `T(b)` predicates re-derive the §8.8 verdicts from explicit per-rung evidence rather than from narrative synthesis. Cross-check: matrix v1 (rebuild step 5) must reproduce §8.8 verdicts when restricted to currently-controlled rungs. |
| `2026-04-25_task-trace-canonical.md` P/T/R/RA regimes | This scope's L4 positivity (`n_T(p, b) > n_T_null(p, b) AND n_T ≥ 1`) is the cohort aggregator over the canonical T-regime predicate. The canonical defines the leaf-set predicate; this scope defines the cohort decision. |
| `2026-04-26_continuous-trace-matrix.md` three-control consolidation | This scope's L1 positivity uses `ρ_split` (Run A) compared against `ρ_null_drift` (Run C) — exactly the consolidation that report defined. Run B (cross-probe) is not load-bearing for the L1 verdict but is reported alongside as anatomy-control evidence. |
| `2026-04-25_measure-ledger.csv` status column | This scope's evaluation eligibility = (ledger.status == "current") AND (audit_21 verdict == "current"). Both conditions must hold. |
| `2026-04-29_measure-correctness-audit.md` (Phase 0) | This scope is a downstream consumer. Phase 0 must be gate-green for every rung this scope evaluates. The audit's φ-checks are necessary preconditions for rules in this scope to apply. |

## Implementation plan

**Producer:** `scripts/01_compute/audit/audit_18_cohort_coverage_matrix.py`
(per rebuild plan §10 step 3). Reads:

- `data/reports/imcoh_vi/h2_partition_multiscale_raw.csv` → L5(k)
- `data/audit/dvi_split_baseline/dvi_split_baseline_n10_imcoh_abs.csv` → L5(k) null
- `data/audit/dvi_split_baseline/dvi_hrel_n10_imcoh_abs.csv` → L5(h_rel)
- `data/reports/imcoh_continuous_trace/per_cell_summary_split.csv` → L1
- `data/reports/imcoh_continuous_trace/controls_summary.csv` → L1 nulls
- `data/reports/imcoh_vi/h2e_split_half_rho_raw.csv` → L1 drift floor
- `data/audit/trace_modules/trace_subtrees_n10_imcoh_abs.csv` → L4
- `data/audit/per_patient_hierarchy_cnp/leaf_assignment.csv` → L6
- `data/audit/per_patient_hierarchy_mspc/multiscale_assignment.csv` → L7
- `data/audit/per_patient_hierarchy_cohesion/leaf_assignment.csv` → L4_aux
- `data/audit/cohort_metadata.csv` → dissenter flags
- `data/audit/dvi_split_baseline/dmin_dmax_inventory.csv` → eligibility gate

Outputs:
- `data/audit/cohort_coverage_matrix/cohort_coverage_matrix_n10_imcoh_abs.csv`
- `data/audit/cohort_coverage_matrix/triangulation_n10_imcoh_abs.csv`
  (one row per band: `b, T(b), V(b, L1), V(b, L5_k), V(b, L5_hrel), V(b, L4), V(b, L6), V(b, L7), V(b, L3)`)
- `data/outputs/figures/section6/cohort_coverage_matrix_n10_imcoh_abs.pdf`

**Library entrypoints (must reuse):**
`lrg_eegfc.utils.metrics.hypothesis.{wilcoxon_z, bh_fdr, rank_biserial}`,
`lrg_eegfc.utils.metrics.tree.{simpson_neff}` (for `n_eff` ineligibility check).

**Acceptance gate.** No matrix-cell evaluation runs until this scope
report is committed AND `2026-04-29_cohort-coverage-matrix.md`
(rebuild plan §7.7) is also committed. The two documents are written
as a paired-PR; commit both before audit_18.

## Open questions

1. **Eligibility threshold `0.5` is provisional.** A sensitivity sweep
   over `{0.3, 0.5, 0.7}` is owed in the final reconciliation. If
   verdict labels change under this sweep, document the dependency.
2. **L3 KC λ-sweep direction convention.** `Δ_KC < 0` means the test-tree
   is closer to post than to pre. Pre-register: trace ⇔ negative Δ_KC.
   Confirmed in `2026-04-29_kc-lambda-sweep.md` (when written).
3. **Per-leaf 2D corner-distance affinity formula consistency.**
   Cohesion-CBR (J_rpre, J_task), CNP (ρ_pre_post, ρ_task_post),
   MSPC (J_pp, J_tp) — three different 2D planes, same regime
   labels. Pre-register: each rung uses the formula in *its own*
   scope report; cross-rung label consistency is asserted in the
   triangulation predicate, not enforced at the per-cell statistic.
4. **What if the matrix v1 disagrees with `2026-04-27` §8.8 verdicts?**
   By construction, v1 cells have no controls for L3 / L4-control
   / L5(h_rel)-control / L6-control / L7-control. Disagreements
   are expected and constitute the *control-gap surface*. The
   reconciliation report (rebuild plan §10 step 19) must explicitly
   trace each disagreement either to "control gap, resolved by §7
   addition" or to "principled disagreement, captured by §3 of the
   rebuild plan".
5. **Triangulation under partial control coverage.** A band with
   `V(b, L1) = positive` (controlled) and `V(b, L5(k)) = positive`
   (controlled) and `V(b, L7) = uncontrolled-exploratory` should be
   labelled *headline-triangulated-pending-L7-control* rather than
   strictly *headline-triangulated*. This nuance is not in §2.4 above
   for compactness; the implementation should emit the
   `pending-<rung>-control` suffix when applicable.
