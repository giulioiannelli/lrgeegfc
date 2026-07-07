---
name: 2026-04-28_edge-vs-hierarchy-discriminability
type: scope
era: IMCOH_ABS × COHORT_N10
status: draft
created: 2026-04-28
updated: 2026-04-28
pointers:
  - 2026-04-25_task-trace-canonical.md
  - 2026-04-26_continuous-trace-matrix.md
  - 2026-04-28_functional-tree-distance.md
  - 2026-04-28_residual-subspace-trace.md
  - .agents/reports/2026-04-24_post-mortem-scalar-session.md
  - .agents/reports/2026-04-25_task-trace-audit-and-recovery.md
---

# Edge-vs-hierarchy substrate discriminability audit

**Renormalization head.** Sanity-check the load-bearing premise of the
LRG framework before any further measure: does the diffusion-induced
hierarchy `(M, T)` resolve coherent multi-edge phase shifts that the
raw FC matrix cannot resolve at the single-edge level? Operationally,
compare phase discriminability of three distances — edge-Frobenius
`d_E`, ultrametric-Spearman `d_H^S`, tree-Kendall-Colijn `d_H^{KC}` —
against a within-`rest_pre` split-half null, and read off the gap.
A gap of the form `n_+(d_H) ≫ n_+(d_E)` is the empirical signature of
"memory hidden in information-diffusion paths"; no gap means LRG is a
projection of an already-detectable signal; no signal anywhere means
the FC substrate isn't where memory lives.

## Why this audit exists

The post-mortem (`2026-04-24_post-mortem-scalar-session.md`) and the
audit-and-recovery handoff (`2026-04-25_task-trace-audit-and-recovery.md`)
both close the loop on scalar-gate hypothesis tests but leave the
substrate-level premise unaudited. Every measure currently active —
MRL, CBR variants, continuous-trace matrix, functional-tree distance
`δ(τ)`, residual-subspace alignment — *assumes* the LRG hierarchy
carries information not present in the FC matrix. This scope tests
that assumption directly. It is a Phase −1 substrate audit; it does
not propose a new task-trace measure.

## Notation

- Patients `p ∈ COHORT_N10`; bands `b ∈ {δ, θ, α, β, γ_l, γ_h}`;
  phases `Φ = {pre, task, post}` aliasing `{rest_pre, task_test, rest_post}`.
- Per `(p, b, φ)`: weighted FC matrix `A^φ ∈ ℝ^{N×N}` from `imcoh_abs`,
  restricted to giant-component intersection
  `V*(p, b) = ⋂_{φ} GC(|A^φ|)`. All matrices restricted to `V* × V*`.
- Laplacian `L^φ = D^φ − A^φ`; spectrum `(λ^φ_i, u^φ_i)` with
  `λ^φ_max = max_i λ^φ_i`.
- Heat-kernel density `ρ̂^φ(τ) = e^{−τ L^φ} / Tr(e^{−τ L^φ})`
  evaluated at `τ_φ = 1 / λ^φ_max` (codebase default; the only
  well-defined operating point per `lrg_outlier_case_fully_connected.md`).
- Communication / ultrametric distance
  `M^φ_ij = (1 − δ_ij) / ρ̂^φ_ij(τ_φ)` (Villegas 2025 §I).
- Single-linkage tree `T^φ = single_linkage(M^φ)`, the dendrogram
  object compared at the hierarchy level.
- Time-series `x^φ(t) ∈ ℝ^{N × T_φ}` for split-half null construction.

## Definitions

### Three distance functions on the substrate ladder

**Edge-level distance** (substrate-blind, treats edges iid):
```
d_E(φ_A, φ_B) = 1 − corr_P( triu(A^{φ_A}), triu(A^{φ_B}) )
```
Pearson on upper-triangular FC entries. No spectrum, no diffusion,
no hierarchy. Coherent multi-edge perturbations are not amplified.

**Hierarchy-level distance, scale-aware** (Spearman on ultrametric):
```
d_H^S(φ_A, φ_B) = 1 − corr_S( triu(M^{φ_A}), triu(M^{φ_B}) )
```
Sensitive to the rank-ordering of pairwise diffusion distances; robust
to overall scale of `M`. Picks up coherent reorderings of the
ultrametric induced by the heat kernel.

**Hierarchy-level distance, scale-blind** (Kendall–Colijn on tree):
```
d_H^{KC}(φ_A, φ_B) = KC( T^{φ_A}, T^{φ_B} )
```
Topology-only tree distance (Kendall & Colijn 2016) via
`lrg_eegfc.utils.metrics.tree_distance.kendall_colijn_distance` (verify
on import; promote if missing). Discards heights; pure dendrogram
topology.

The triple `(d_E, d_H^S, d_H^{KC})` forms a discriminability ladder:
edge → diffusion-distance → tree-only.

### Within-baseline split-half null

For each `(p, b)`: partition `x^{pre}(t)` into two non-overlapping
temporal halves `x^{pre_a}, x^{pre_b}` of equal length (random temporal
permutation, then split). Compute `A^{pre_a}, A^{pre_b}` via the same
`imcoh_abs` pipeline; restrict to `V*`; build `L, ρ̂, M, T`. Repeat
with `n_split = 50` independent random permutations.

Null distribution per distance:
```
N_d(p, b) = { d(pre_a, pre_b)_i : i = 1, ..., n_split }   d ∈ {E, S, KC}
```

### Discriminability statistic

For each `(p, b)`, distance `d`, phase pair `π ∈ {(task, pre), (post, pre)}`:
```
Z_d(p, b; π) = ( d_obs(π) − median(N_d(p, b)) ) / MAD(N_d(p, b))
```
robust z-score against the within-baseline floor.

### Cohort aggregation

```
n_+(b, d, π) = #{ p ∈ COHORT_N10 \ {Pat_03} : Z_d(p, b; π) > 2 }
```
Pat_03 (1024 Hz outlier, see `pat03_nperseg.md`) reported separately
as an out-of-pool diagnostic, not counted in `n_+`. Cohort-positive
verdict per the canonical threshold from `task-trace-canonical.md`:
`V(b, d, π) = positive` if `n_+(b, d, π) ≥ 8` (i.e. ≥ 8 of 9 in-pool).

## Properties

- **Range:** `Z_d ∈ ℝ`; `n_+ ∈ {0, ..., 9}`; verdict in `{positive, negative}`.
- **Cohort-aggregable:** scalar per `(p, b)` aggregates without leafset
  matching, anatomy maps, or partition selection (matches the canonical
  aggregation grade in `task-trace-canonical.md`).
- **Substrate-discriminating (the load-bearing observable):** the gap
  `Δn_+(b, π) = n_+(b, d_H^S, π) − n_+(b, d_E, π)`. Positive `Δn_+`
  with `n_+(d_H^S) ≥ 8` and `n_+(d_E) < 8` is the empirical signature
  of the user's hypothesis: "coherent multi-edge changes hidden in
  diffusion paths." Reported per band.
- **Decoupled from k-cuts and Ψ:** no flat-clustering at any `h_rel`,
  no Ψ-selector, no partition machinery. The 2026-04-28 Ψ τ-scan
  diagnostic showed Ψ has no purchase on `imcoh_abs` trees; this audit
  deliberately avoids that machinery.
- **What it does NOT measure:** *which* nodes / subtrees / scales carry
  the trace. This is a substrate-existence test, not a localization
  test. If positive, hand off to MRL / CBR / partition-multiscale /
  functional-tree-distance for localization. If negative for a given
  `(b, d, π)`, those localization tools are testing for signal that
  isn't there at the substrate's discriminability grade.
- **What it cannot distinguish:** edge-level signal that is genuinely
  weak vs edge-level signal that `d_E` underweights because of an
  ill-chosen edge norm (Pearson on `triu` may be unfairly weak; see
  Open Questions).

## Caveats and failure modes

| Caveat | Mitigation |
|---|---|
| `τ = 1/λ_max` is a single operating point per `lrg_outlier_case_fully_connected.md` — no τ sweep | Audit fixed at `τ_φ = 1/λ^φ_max` per phase; `τ`-scan is the scope of [`2026-04-28_functional-tree-distance.md`](2026-04-28_functional-tree-distance.md). |
| `rest_pre` non-stationarity inflates the null and biases toward null verdicts | Stationarity diagnostic per patient: variance ratio between halves + KL of channel amplitude distribution. Flag `(p, b)` if drift exceeds threshold. |
| Probe bias dominates coarse hierarchy structure (see `feedback_probe_bias_critical`) | Probe bias is in both observed AND null distances → it factors as a constant offset; `Z` insensitive. Robustness check: rerun on probe-bias-zeroed FC for one band. |
| Edge norm choice could be unfairly weak | Run `d_E ∈ {1−corr_P, Frobenius_norm, 1−corr_S}` and report all three; verdict only declared if at least two agree. |
| Split-half halves not independent (autocorrelation) | Use random temporal permutation before splitting (destroys time-order autocorrelation but preserves stationary FC); compare to chronological split as sensitivity check. |
| MAD = 0 (degenerate null, all splits identical) | Skip `(p, b)` and log; emit "ineligible" rather than divide-by-zero `Z`. |
| `n_split = 50` MAD unstable | Convergence check: `n_split ∈ {20, 50, 100, 200}` on Pat_06 × β. Promote default to 100 if MAD unstable below. |
| Pat_03 1024 Hz spectrum scale differs | Reported separately; not counted in `n_+`. Per cohort lock at n=10, in-pool denominator = 9 for `n_+`. |
| Single FC method | `imcoh_abs` only per current era; cross-method robustness is out of scope. |

## Pseudocode

```
INPUT:  P = COHORT_N10, B = {δ,θ,α,β,γ_l,γ_h}, n_split = 50
OUTPUT: discriminability_audit.csv with one row per (p, b, d, π)

FOR each (p, b) IN P × B:
    # 1. Build per-phase substrate
    FOR each φ IN {pre, task, post}:
        x_φ ← load_timeseries(p, φ)
        A_φ ← imcoh_abs(x_φ, b, nperseg_for_fs(fs(p)))
        keep_φ ← giant_component_nodes(A_φ)
    V* ← ⋂_φ keep_φ
    IF |V*| < 10: EMIT (p, b, *, *, "ineligible"); CONTINUE

    FOR each φ:
        A_φ ← restrict(A_φ, V*)
        L_φ ← laplacian(A_φ)
        (eigvals_φ, eigvecs_φ) ← eigh(L_φ)
        τ_φ ← 1 / max(eigvals_φ)
        ρ_φ ← rho_tau(eigvals_φ, eigvecs_φ, τ_φ)
        M_φ ← (1 − I) / ρ_φ
        T_φ ← single_linkage(condensed(M_φ))

    # 2. Observed phase-pair distances
    FOR π IN {(task, pre), (post, pre)}:
        (φ_A, φ_B) ← π
        d_E_obs   ← 1 − corr_P(triu(A_{φ_A}), triu(A_{φ_B}))
        d_HS_obs  ← 1 − corr_S(triu(M_{φ_A}), triu(M_{φ_B}))
        d_HKC_obs ← kendall_colijn(T_{φ_A}, T_{φ_B})

        # 3. Within-baseline split-half null (shared across distances)
        N_E, N_HS, N_HKC ← [], [], []
        FOR i IN 1..n_split:
            (x_a, x_b) ← random_split_half(x_pre, seed=i)
            A_a, A_b ← imcoh_abs(x_a, b, ...), imcoh_abs(x_b, b, ...)
            (restrict to V*; build L_a, L_b, ρ_a, ρ_b, M_a, M_b, T_a, T_b)
            APPEND N_E   ← 1 − corr_P(triu(A_a), triu(A_b))
            APPEND N_HS  ← 1 − corr_S(triu(M_a), triu(M_b))
            APPEND N_HKC ← kendall_colijn(T_a, T_b)

        # 4. Discriminability z-score
        FOR (d_obs, N_d, label) IN {(d_E_obs, N_E, "E"),
                                    (d_HS_obs, N_HS, "HS"),
                                    (d_HKC_obs, N_HKC, "HKC")}:
            IF MAD(N_d) == 0: Z ← NaN; flag "degenerate_null"
            ELSE: Z ← (d_obs − median(N_d)) / MAD(N_d)
            EMIT (p, b, label, π, d_obs, median(N_d), MAD(N_d), Z)

# 5. Cohort aggregation
FOR each (b, d_label, π):
    Z_pool ← {Z(p, b, d_label, π) : p ∈ COHORT_N10 \ {Pat_03}, Z finite}
    n_plus ← #{ z ∈ Z_pool : z > 2 }
    verdict ← "positive" if n_plus ≥ 8 else "negative"
    EMIT_BAND (b, d_label, π, n_plus, |Z_pool|, verdict)

# 6. Substrate gap
FOR each (b, π):
    Δn_plus ← n_plus(b, "HS", π) − n_plus(b, "E", π)
    EMIT_GAP (b, π, Δn_plus, hypothesis_signature ∈ {"hierarchy_only",
              "edge_only", "both", "neither"})
```

## Visualization spec

**Page 1 — discriminability ladder (the headline).**
Layout: 6 rows = bands, 2 columns = phase pairs `(task,pre)` and
`(post,pre)`. Each panel: grouped bar chart with three bars
`{d_E, d_H^S, d_H^{KC}}`, height = `n_+`. Horizontal dashed line at
`n_+ = 8` (cohort threshold). Reading rules:
- bars below 8 across the row → no substrate-level signal in that band;
- only `d_H^S` / `d_H^{KC}` bar above 8 → **hierarchy resolves what
  edges miss** (the user's hypothesis confirmed for that band);
- all three bars above 8 → signal exists at edge level too;
- only `d_E` above 8 → pathological (diffusion destroys signal),
  revisit τ.

**Page 2 — per-patient Z-scatter.** One panel per band, two markers
per patient (one per phase pair, distinguished by colour).
X-axis: `Z_E`. Y-axis: `Z_{HS}`. Identity line `y = x`. Reading rule:
points above identity = hierarchy more discriminative than edges; below
= edges more discriminative. Quadrant lines at `Z = 2` mark the
significance threshold per axis.

**Page 3 — within-baseline null distributions.** Grid: 6 bands × 3
distances. Each panel: violin / KDE of `N_d` pooled across patients
(coloured by patient), with vertical lines at `d_obs` for each phase
pair. Reading rule: where the obs lines fall in the null tail per
patient.

**Page 4 — substrate-gap summary table.** Single page, 6 × 2 cells
(band × phase pair). Each cell shows `Δn_+` and the verdict label
(`hierarchy_only`, `both`, `edge_only`, `neither`). The headline
artefact for the audit verdict report.

## Connection to prior tools

| This audit | Complements | Subsumes / replaces |
|---|---|---|
| `d_H^S` at `τ = 1/λ_max` | [`2026-04-28_functional-tree-distance.md`](2026-04-28_functional-tree-distance.md) (`δ(τ)` curves at snapshot + integral over `[1/λ_max, τ*]`) | None — pins the substrate question to a single τ; the τ-scan answers a different, downstream question. |
| Within-baseline split-half null | H2e split-half drift floor (existing scalar test) | Generalises H2e's null structure to three substrate levels; supersedes within-baseline normalisation hand-rolled inside individual measures. |
| `n_+ ≥ 8/9` cohort verdict | Π_T `≥ 8/10` from `task-trace-canonical.md` | Same threshold convention; substrate-existence rather than regime-classification. |
| `d_E` baseline | None prior — no current measure runs an edge-only head-to-head against the hierarchy | Provides the missing comparator for every existing measure that operates on `M` or `T` without first establishing that `M` / `T` add resolution. |

Does **not** replace MRL / CBR / partition-multiscale / continuous-trace
matrix / functional-tree distance / residual-subspace. This audit is
*upstream* of all of them: it tests whether they are looking for signal
on the right substrate. If positive, those tools localise the trace.
If negative, those tools are testing in the dark.

## Implementation plan

- **Audit script:** `scripts/01_compute/audit/audit_19_substrate_discriminability.py`.
- **Library entry points reused** (no new private helpers):
  - `lrg_eegfc.workflow.fc.load_fc_matrix` (per-phase FC).
  - `lrg_eegfc.utils.io.patient.load_timeseries` (raw time-series for split-half).
  - `lrg_eegfc.utils.fc.imcoh.imcoh_abs` (FC on split halves).
  - `lrg_eegfc.utils.lrg.<rho_tau>` (heat-kernel density from spectrum;
    confirm exact name on import).
  - `lrg_eegfc.utils.metrics.tree_distance.kendall_colijn_distance`
    (verify; if missing, promote per `coding-rules.md`).
  - `lrgsglib.nx_patches.funcs.get_giant_component`.
  - `scipy.cluster.hierarchy.linkage` (single linkage on condensed `M`).
- **New helper, library-promotion candidate:**
  `lrg_eegfc.utils.metrics.null.split_half_null(ts, distance_fn, n_split, seed)`
  — only promote if reused outside this script (e.g. by the
  functional-tree-distance scope, which needs the same null).
- **Outputs:**
  - `data/audit/substrate_discriminability/discriminability_audit.csv`
    (long format, one row per `(p, b, d, π)`).
  - `data/audit/substrate_discriminability/null_distributions.npz`
    (per `(p, b)` arrays of `N_E, N_{HS}, N_{HKC}`).
  - `data/audit/substrate_discriminability/diagnostic.pdf` (Pages 1–4).
  - Verdict report:
    `.agents/reports/2026-04-28_substrate-discriminability-verdict.md`
    with the verdict per band + cohort-level summary table + a
    1-paragraph implication for the rebuild plan.

## Pre-registered verdict semantics

Before running, fix:

- **Hypothesis confirmed (hierarchy is load-bearing):** for ≥ 1 band,
  `n_+(d_H^S) ≥ 8 ∧ n_+(d_E) < 8` for `π = (task, pre)` AND
  `π = (post, pre)`. Strongest if `d_H^{KC}` agrees with `d_H^S`.
- **Hierarchy is a projection of edge signal:** for ≥ 3 bands,
  `n_+(d_H^S) ≥ 8 ∧ n_+(d_E) ≥ 8`. LRG remains useful for
  interpretation but isn't required for detection.
- **No FC-substrate memory:** for all 6 bands and both phase pairs,
  `n_+(d_H^S) < 8 ∧ n_+(d_E) < 8`. Closes the static-FC LRG program
  for memory; pivots to temporal-dynamics or non-FC features.
- **Pathological:** `n_+(d_E) ≥ 8` while `n_+(d_H^S) < 8` for any band
  → flag, do not interpret as confirmation in either direction;
  diagnose τ choice or norm choice.

## Open questions

- **Edge norm (`d_E`).** Default `1 − corr_P(triu)` may underweight
  high-magnitude coordinated changes. Robustness panel: also report
  Frobenius and Spearman variants. Verdict declared only with ≥ 2-of-3
  agreement.
- **Tree distance choice.** KC default; weighted Robinson–Foulds
  (already in `tree_distance.py`) as second-tier robustness. Decide
  before runtime.
- **`τ` per phase vs shared.** Default `τ_φ = 1/λ^φ_max`
  (intrinsic-resolution per phase). Alternative: shared
  `τ = 1/min_φ(λ^φ_max)` (same physical diffusion time). Open whether
  the verdict survives the convention switch; included as sensitivity.
- **`n_split`.** Default 50; convergence check `{20, 50, 100, 200}` on
  Pat_06 × β. Promote default to 100 if MAD unstable.
- **Stationarity gate.** What variance-ratio / KL threshold flags a
  `(p, b)` as drift-contaminated? Pre-register threshold from a
  preliminary run on Pat_02 / Pat_06 baselines.
- **Phase 0 ordering.** [`2026-04-29_measure-correctness-audit.md`](2026-04-29_measure-correctness-audit.md)
  is the BLOCKING Phase 0 for downstream measure-cell evaluation. This
  substrate audit is conceptually *upstream* (Phase −1): if substrate
  discriminability is null, measure-correctness for downstream measures
  loses operational meaning. Open: should this audit gate the
  measure-correctness audit, run in parallel, or stay as an independent
  diagnostic? Default for now: independent, run before any new measure
  is greenlit; do not block Phase 0 of the existing rebuild.
- **What if results are mixed across bands?** E.g. `δ` shows
  hierarchy-only signature, `α` shows both, `θ/γ` show neither. The
  scope as written reports per-band verdicts without aggregating to a
  single global verdict — deliberately, to avoid the scalar-gate
  failure mode flagged in the post-mortem. Open question for the
  verdict report: how to communicate band-heterogeneous outcomes
  without re-introducing scalar aggregation.
