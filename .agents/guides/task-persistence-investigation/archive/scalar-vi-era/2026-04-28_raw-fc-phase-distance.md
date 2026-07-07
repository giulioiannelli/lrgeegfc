---
name: 2026-04-28_raw-fc-phase-distance
type: scope
era: IMCOH_ABS × COHORT_N10
status: draft
created: 2026-04-28
updated: 2026-04-28
pointers:
  - 2026-04-28_edge-vs-hierarchy-discriminability.md
  - 2026-04-26_continuous-trace-matrix.md
  - 2026-04-25_task-trace-canonical.md
  - .agents/reports/2026-04-24_post-mortem-scalar-session.md
---

# Raw-FC phase-distance audit (step 0, pre-LRG)

**Renormalization head.** The most basic question that must answer
"yes" before any LRG / hierarchy / multiscale tool is justified: do
the FC matrices `A^φ` themselves move detectably between phases? One
patient first (eyeball + three scalar distances against a within-
`rest_pre` split-half null), then cohort-wide counts. No diffusion,
no spectrum, no dendrogram, no τ. If the matrix doesn't budge, no
downstream cleverness saves it; if it does budge, we have quantified
the edge-level signal that every LRG-based measure must beat to
justify itself.

## Why this exists

Every active task-trace measure (MRL, CBR family, continuous-trace
matrix, functional-tree distance, residual-subspace, the substrate
ladder in [`2026-04-28_edge-vs-hierarchy-discriminability.md`](2026-04-28_edge-vs-hierarchy-discriminability.md))
implicitly assumes the FC matrices carry a phase-distinguishable
signal somewhere — they then ask which transformation surfaces it
best. This step-0 audit makes the implicit assumption explicit and
testable on the simplest possible substrate: the raw FC matrix `A`,
unfiltered, untransformed. It is meant to run first; everything
downstream conditions on its outcome.

## Notation

- Patients `p ∈ COHORT_N10`; bands `b ∈ {δ, θ, α, β, γ_l, γ_h}`;
  phases `Φ = {pre, task, post}` aliasing `{rest_pre, task_test, rest_post}`.
- Per `(p, b, φ)`: weighted FC matrix `A^φ ∈ ℝ^{N×N}` from `imcoh_abs`,
  restricted to the giant-component intersection
  `V*(p, b) = ⋂_{φ} GC(|A^φ|)`.
- Time-series `x^φ(t) ∈ ℝ^{N × T_φ}` for split-half null construction.
- Phase pairs `Π = {(pre, task), (pre, post), (task, post)}`.

## Definitions

### Three distance functions on the raw FC matrix

```
d_P(A, B) = 1 − corr_P( triu(A, k=1), triu(B, k=1) )       # Pearson
d_S(A, B) = 1 − corr_S( triu(A, k=1), triu(B, k=1) )       # Spearman (rank)
d_F(A, B) = ‖A − B‖_F / √( ‖A‖_F · ‖B‖_F )                # scale-normalised Frobenius
```

`d_P, d_S ∈ [0, 2]`; `d_F ∈ [0, ∞)` but with the chosen normalisation
typically lives in `[0, 2]`. Three distances are reported because
each is sensitive to a different aspect of the change: linear
co-variation (`d_P`), monotone re-ranking (`d_S`), absolute magnitude
shift (`d_F`). The verdict is reported per-distance and as a
robustness consensus (≥ 2-of-3 agreement).

### Within-`rest_pre` split-half null

For each `(p, b)`: random temporal permutation of `x^{pre}(t)`, then
split into two non-overlapping halves of equal length. Compute
`A^{pre_a}, A^{pre_b}` with the same `imcoh_abs` pipeline (same
`nperseg`, same band-averaging), restrict to `V*`, evaluate the three
distances. Repeat `n_split = 50` times.

```
N_d(p, b) = { d(A^{pre_a}_i, A^{pre_b}_i) : i = 1..n_split }   d ∈ {P, S, F}
```

### Robust z-score per phase pair

For each `(p, b)`, distance `d`, phase pair `π ∈ Π`:
```
Z_d(p, b; π) = ( d(A^{φ_A}, A^{φ_B}) − median(N_d(p, b)) ) / MAD(N_d(p, b))
```

### Triangle-style persistence test

A simple, model-free operationalisation of "post closer to task than
to pre" — no LRG required:
```
T_d(p, b) = d(A^{task}, A^{post}) − d(A^{pre}, A^{task})
```
- `T_d < 0`: `post` is closer to `task` than `pre` is to `task` →
  candidate persistence signature in raw FC.
- `T_d > 0`: `post` is closer to `pre` than to `task` → reset
  signature.
- `T_d ≈ 0`: ambiguous.

`T_d` is a per-patient scalar; cohort summary is the count
`#{p : T_d(p, b) < 0}` per band. **Important:** `T_d < 0` is necessary
for "raw FC carries a persistence trace" but not sufficient — it can
arise from generic non-stationarity. Pair with `Z_P((post, pre)) > 2`
to require the post-shift to also exceed within-baseline jitter.

### Cohort aggregation

```
n_+(b, d, π) = #{ p ∈ COHORT_N10 \ {Pat_03} : Z_d(p, b; π) > 2 }
```
Pat_03 (1024 Hz outlier, `pat03_nperseg.md`) reported separately as
out-of-pool diagnostic; in-pool denominator = 9. Cohort-positive per
the canonical threshold from `task-trace-canonical.md`:
`V(b, d, π) = positive` if `n_+(b, d, π) ≥ 8`.

## Properties

- **Range:** `Z_d ∈ ℝ`; `n_+ ∈ {0, ..., 9}`; `T_d ∈ ℝ`; verdicts
  `{positive, negative}` per `(b, d, π)`.
- **Cohort-aggregable** without partition selection, leafset matching,
  or anatomy mapping — a single scalar per `(p, b)` per distance.
- **Decoupled from LRG.** No spectrum, no `ρ̂`, no τ, no dendrogram,
  no flat-clustering, no Ψ. Failure mode of the LRG outlier-case
  (`lrg_outlier_case_fully_connected.md`) does not apply here.
- **What it cannot detect (the reason LRG might still matter):**
  reorganisations that preserve the marginal edge distribution while
  rerouting which-edge-talks-to-which along diffusion paths. If the
  per-edge values shift in a coordinated, low-rank way that integrates
  on the heat kernel but cancels in `d_P / d_S / d_F`, this audit
  returns null while the LRG hierarchy still detects signal. That is
  precisely the case the substrate-ladder scope tests, *conditional on
  this audit returning null*.
- **What it cannot distinguish:** generic non-stationarity (electrode
  drift, vigilance, recording setup change between phases) from
  task-induced reorganisation. The within-`rest_pre` split-half null
  controls only for jitter *within* `rest_pre`, not for systematic
  pre→post drift unrelated to the task.

## Caveats and failure modes

| Caveat | Mitigation |
|---|---|
| `rest_pre` non-stationarity inflates the null and biases toward null verdicts | Stationarity diagnostic: variance ratio between halves + KL of channel-amplitude distribution. Flag `(p, b)` if drift exceeds threshold; report flagged cells separately. |
| ~~Probe-bias edges~~ — under `imcoh_abs` the zero-phase-lag component is killed by construction (Nolte 2004); same-probe-zeroing is **not** a default control under this FC method (per updated `probe_bias_critical` memory) | No action; do not propose same-probe-zeroing here. |
| Different recording lengths between phases (e.g. `task_test` shorter than `rest_pre`) | Split-half null operates on `rest_pre` only — equal-length halves by construction. Phase-pair `d_obs` uses each phase at its native length; document length per `(p, φ)` in the output CSV for sanity. |
| `MAD = 0` (degenerate null, all splits coincident) | Skip `(p, b)`; emit `Z = NaN` with `flag = "degenerate_null"`. |
| `n_split = 50` MAD instability | Convergence check at `n_split ∈ {20, 50, 100, 200}` on Pat_06 × β. Promote default to 100 if MAD CV > 0.2. |
| Pat_03 1024 Hz spectrum-scale offset | Reported separately, not in `n_+`. |
| `nperseg` mismatch between phases | All phases of a patient use the same `nperseg = nperseg_for_fs(fs(p))` (4096 at 2048 Hz, 2048 at 1024 Hz Pat_03). |
| `imcoh_abs` only | Out of scope to vary FC method. Era convention. |

## Pseudocode

```
INPUT:  P = COHORT_N10, B = {δ,θ,α,β,γ_l,γ_h}, n_split = 50
OUTPUT: raw_fc_phase_distance.csv (one row per (p, b, d, π))
        raw_fc_triangle.csv         (one row per (p, b, d))

# === STEP 1: single-patient (default Pat_06; user-selectable) ===
# Same loop body as below, run for one p first; render Page 1+2 of the
# diagnostic PDF and pause for visual inspection before launching the
# cohort sweep.

FOR each (p, b) IN P × B:
    # Build per-phase FC
    FOR each φ IN {pre, task, post}:
        x_φ      ← load_timeseries(p, φ)
        A_φ_full ← imcoh_abs(x_φ, b, nperseg_for_fs(fs(p)))
        keep_φ   ← giant_component_nodes(A_φ_full)
    V*  ← ⋂_φ keep_φ
    IF |V*| < 10: EMIT (p, b, *, *, "ineligible"); CONTINUE
    FOR each φ:
        A_φ ← restrict(A_φ_full, V*)

    # Phase-pair observed distances
    FOR π IN {(pre, task), (pre, post), (task, post)}:
        (φ_A, φ_B) ← π
        d_P_obs ← 1 − corr_P(triu(A_{φ_A}), triu(A_{φ_B}))
        d_S_obs ← 1 − corr_S(triu(A_{φ_A}), triu(A_{φ_B}))
        d_F_obs ← ‖A_{φ_A} − A_{φ_B}‖_F / √(‖A_{φ_A}‖_F · ‖A_{φ_B}‖_F)

        # (z-score and emit per distance below, after null is built)

    # Within-pre split-half null (built once per (p, b), shared across π)
    N_P, N_S, N_F ← [], [], []
    FOR i IN 1..n_split:
        (x_a, x_b) ← random_split_half(x_pre, seed=i)
        A_a       ← restrict(imcoh_abs(x_a, b, nperseg), V*)
        A_b       ← restrict(imcoh_abs(x_b, b, nperseg), V*)
        APPEND N_P ← 1 − corr_P(triu(A_a), triu(A_b))
        APPEND N_S ← 1 − corr_S(triu(A_a), triu(A_b))
        APPEND N_F ← ‖A_a − A_b‖_F / √(‖A_a‖_F · ‖A_b‖_F)

    # Z-scores
    FOR π:
      FOR (d_obs, N_d, label) IN {(d_P_obs, N_P, "P"),
                                  (d_S_obs, N_S, "S"),
                                  (d_F_obs, N_F, "F")}:
          IF MAD(N_d) == 0: Z ← NaN; flag ← "degenerate_null"
          ELSE: Z ← (d_obs − median(N_d)) / MAD(N_d)
          EMIT_ROW (p, b, label, π, d_obs, median(N_d), MAD(N_d), Z, flag)

    # Triangle persistence scalar
    FOR label IN {P, S, F}:
        T_d ← d_obs(task,post; label) − d_obs(pre,task; label)
        EMIT_TRIANGLE (p, b, label, T_d)

# === STEP 2: cohort aggregation ===
FOR each (b, d_label, π):
    Z_pool ← {Z(p, b, d_label, π) : p ∈ COHORT_N10 \ {Pat_03}, Z finite}
    n_plus ← #{ z ∈ Z_pool : z > 2 }
    verdict ← "positive" if n_plus ≥ 8 else "negative"
    EMIT_BAND (b, d_label, π, n_plus, |Z_pool|, verdict)

FOR each (b, d_label):
    n_persist  ← #{ p \ {Pat_03} : T_d(p, b, d_label) < 0
                                   AND Z_d(p, b, (pre,post)) > 2 }
    EMIT_PERSIST (b, d_label, n_persist, verdict)
```

## Visualization spec

**Page 1 — single-patient FC matrices.** Layout: 6 rows = bands, 3
columns = phases. Each panel: `imshow(A^φ)` with shared colormap and
shared color limits per band (so cross-phase contrast is honest).
Channels ordered by anatomical region (Desikan–Killiany), not raw
index. Reading rule: visible block-structure shifts between columns
within a row = candidate signal in that band.

**Page 2 — single-patient phase-difference matrices.** Layout: 6 rows
= bands, 2 columns = `(A^task − A^pre)` and `(A^post − A^pre)`.
Diverging colormap centred at 0; per-band symmetric limits at the
99th percentile of `|Δ|`. Reading rule: spatially coherent
positive/negative regions = coordinated edge change; speckle = noise.

**Page 3 — single-patient null vs observed.** 6 rows × 3 columns
(bands × distance type). Histogram of `N_d(p, b)`; vertical lines for
the three observed `d(A^{φ_A}, A^{φ_B})`; annotation of `Z`. Reading
rule: where the obs lines fall in the null tail.

**Page 4 — cohort `n_+` bars.** Layout: 3 rows = phase pairs, 6
columns = bands. Each panel: three grouped bars `{d_P, d_S, d_F}`,
height = `n_+`, horizontal dashed line at 8. Reading rule: bars at /
above 8 across distances = robust phase difference at edge level for
that band.

**Page 5 — cohort triangle scatter.** 6 panels (one per band). X-axis:
`d(pre, task)`. Y-axis: `d(task, post)`. Each marker = one patient
(Pat_03 highlighted out-of-pool). Identity line `y = x`; markers
**below** identity = `T < 0` = candidate persistence. Reading rule:
density of below-identity markers per band = cohort-level FC-edge
persistence count.

## Connection to prior tools

| This audit | Complements | Subsumes / replaces |
|---|---|---|
| `d_P / d_S / d_F` on raw `A` against split-half null | `2026-04-26_continuous-trace-matrix.md` (operates on the ultrametric `D = 1/ρ`, not on `A`) | Provides the missing pre-LRG baseline; `continuous-trace-matrix` already assumes the diffusion transform is justified. |
| Step-0 substrate verdict | [`2026-04-28_edge-vs-hierarchy-discriminability.md`](2026-04-28_edge-vs-hierarchy-discriminability.md) | Strict prerequisite. The `d_E` rung of that scope is operationalised here in isolation; the hierarchy-vs-edge ladder runs only if step 0 returns positive somewhere. If step 0 returns null cohort-wide, the ladder pivots to the question "does LRG resolve a signal that `A` cannot resolve at all" — same scope, different load-bearing claim. |
| Triangle test `T_d < 0` | Cohesion-CBR per-patient persistence census | Edge-level analogue, no tree machinery; cheap sanity check before any tree-distance test. |
| Within-`rest_pre` split-half null | H2e split-half drift floor; null structure in [`2026-04-28_edge-vs-hierarchy-discriminability.md`](2026-04-28_edge-vs-hierarchy-discriminability.md) | Same null construction; this scope adopts and reports it for the raw-FC substrate only. |

Does **not** replace any LRG-based measure. This audit's purpose is
to set the floor; the LRG-based measures remain necessary if and only
if (a) step 0 returns positive *and* the LRG transform reveals
additional structure, or (b) step 0 returns null *and* the substrate
ladder shows the hierarchy resolves what edges miss.

## Implementation plan

- **Audit script:** `scripts/01_compute/audit/audit_20_raw_fc_phase_distance.py`.
- **Library entry points reused:**
  - `lrg_eegfc.workflow.fc.load_fc_matrix` (per-phase FC; cache-aware).
  - `lrg_eegfc.utils.io.patient.load_timeseries` (raw time-series).
  - `lrg_eegfc.utils.fc.imcoh.imcoh_abs` (FC on split halves).
  - `lrg_eegfc.config.const.nperseg_for_fs`.
  - `lrgsglib.nx_patches.funcs.get_giant_component`.
  - `numpy.corrcoef`, `scipy.stats.spearmanr`, `numpy.linalg.norm`.
- **New helpers (library-promotion candidates):**
  - `lrg_eegfc.utils.metrics.null.split_half_fc_null(ts, fc_fn, n_split, seed)`
    — promote only on second-use (consumed also by the substrate ladder
    scope; that's the second use, justifying promotion). Lives next to
    `hypothesis.py`.
  - No new distance helpers — `corr_P / corr_S / Frobenius` are
    one-liners.
- **Outputs:**
  - `data/audit/raw_fc_phase_distance/raw_fc_phase_distance.csv`
  - `data/audit/raw_fc_phase_distance/raw_fc_triangle.csv`
  - `data/audit/raw_fc_phase_distance/null_distributions.npz`
  - `data/audit/raw_fc_phase_distance/diagnostic.pdf` (Pages 1–5).
- **Verdict report:**
  `.agents/reports/2026-04-28_raw-fc-phase-distance-verdict.md`
  with the per-band verdict table, the triangle-persistence count,
  flagged cells (stationarity / probe-bias / degenerate-null), and a
  one-paragraph implication for the LRG program.

## Pre-registered verdict semantics

Locked before runtime:

- **Raw FC moves with phase, cohort-wide:** for ≥ 2 bands, ≥ 2 of 3
  distances `{d_P, d_S, d_F}` agree with `n_+ ≥ 8` for `(pre, task)`
  AND `(pre, post)`. → Edge-level signal exists; LRG must beat it.
- **Raw FC carries a persistence signature:** for ≥ 1 band,
  `n_persist ≥ 5` (≥ 5 of 9 in-pool). → Triangle test passes; pre-LRG
  evidence for task-trace.
- **Raw FC null:** for all bands and pairs, `n_+ < 8` across all
  three distances. → No edge-level signal at within-baseline grade;
  hands the question to the LRG substrate ladder
  ([`2026-04-28_edge-vs-hierarchy-discriminability.md`](2026-04-28_edge-vs-hierarchy-discriminability.md))
  with the sharpened claim "diffusion-induced hierarchy resolves what
  edges cannot."
- **Mixed:** verdict reported per-band; no scalar global aggregation,
  per the `feedback_dont_rerun_scalar_tests` lesson.

## Open questions

- **Single-patient default.** Suggest Pat_06 (mid-cohort, 2048 Hz, no
  known data quirks beyond the cohort-wide constraints). User
  override before run; alternative Pat_02 (canonical reference, but
  per `feedback_figure_variety` overused — prefer fresh cases for
  step-0 eyeballing).
- **Distance set.** `{d_P, d_S, d_F}` is the proposal. A fourth
  candidate is mean absolute edge difference (`L1` Frobenius). Decide
  before runtime; over-running distances is cheap, the question is
  reporting clarity.
- **Probe-bias robustness.** Run once with raw `A` and once with
  same-probe edges zeroed for one band (e.g. `β`) to confirm the
  `Z`-score is bias-insensitive. If `Z` shifts substantially, run
  the full audit on probe-zeroed FC.
- **Stationarity gate.** Variance-ratio threshold and channel-KL
  threshold for flagging drift-contaminated `(p, b)` cells.
  Pre-register from a preliminary run on Pat_02 / Pat_06 baselines.
- **`n_split` default.** 50 is the proposal; promote to 100 on MAD-CV
  evidence.
- **Triangle test threshold.** `T_d < 0` is the binary; finer
  pre-registration option is `T_d < − k · MAD(N_d)` for some
  `k ∈ {0.5, 1}` to demand the asymmetry exceed the within-baseline
  jitter, not just be negative. Decide before runtime.
- **Out-of-pool reporting.** Pat_03 reported separately by convention;
  open whether its 1024 Hz scaling makes the `Z`-score directly
  comparable. If not, report only sign of `T_d` and `n_+` membership
  qualitatively.
