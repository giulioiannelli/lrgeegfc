---
name: cross-phase-cophenetic-taxonomy
type: scope
era: IMCOH_ABS / COHORT_N10
status: draft
created: 2026-06-12
updated: 2026-06-12
status_note: >
  Implemented 2026-06-12 (audit_105/106/107). Trace guard 60/60 bit-exact. Two
  changes from this draft, both recorded below: (a) sigma is the RMS of Δbase,
  not median |Δbase| (median is 0 under quantization); (b) the coupled null
  FAILS its validation gate (§2.4) — coordinated strength-preserving null not
  viable, confirming audit_63 — so similarity channels use the geometry baseline
  only. Results: .agents/reports/2026-06-12_cross-phase-taxonomy.md.
scope: >
  Decompose the validated per-pair cophenetic cross-phase signal (the one that
  carries the rho_split trace) into four orthogonal-contrast channels —
  anchor / trace / reset / reorganize — and localize each. 100% rho^coph-native;
  no tree-cutting. Trace channel is identically the validated rho_split.
pointers:
  - .agents/guides/task-persistence-investigation/README.md            # scope-report template
  - .agents/guides/01_project/terminology.md                           # trace/anchor/reset/emergent defs
  - scripts/01_compute/audit/audit_63_split_baseline_surrogate.py       # FC->cophenet pipeline, obs_rho
  - scripts/01_compute/audit/audit_76_pair_displacement_taxonomy.py     # per-pair cloud, cell readings
  - scripts/01_compute/audit/audit_83_localization_matched_strength.py  # obs_trace/concordance/unit_means
  - scripts/01_compute/audit/audit_92_localization_R1000_surrogates.py  # R=1000 surrogate generation
  - src/lrg_eegfc/utils/surrogate/matched_strength.py                   # strength_preserving_shuffle
  - src/lrg_eegfc/utils/io/regions.py                                   # DK anatomy, SUPERSYSTEMS
  - perpair_cophenetic_trace_visualization.md                          # quantization memory (98% ties)
  - localization_audit_plan_2026_05_29.md                              # beta->OFC locked verdict
  - feedback_matched_strength_mandatory.md                             # null mandate
---

# Cross-phase cophenetic taxonomy — anchor / trace / reset / reorganize

**Head.** The per-pair cophenetic cross-phase signal that carries our validated
β/α **trace** (`ρ_split`) decomposes, with *no residual*, into a rigid backbone
(**anchor**, low cross-phase variance) plus a 2-D fluctuation spanned by two
orthonormal contrasts: the **persistence** contrast `φ₁ = (D_post − D_pre)/√2` —
which *is* what `ρ_split` measures — and the **excursion** contrast
`φ₂ = (2·D_task − D_pre − D_post)/√6`, the **reset** channel. Because a perfect
reset (post returns to pre) has `Δrest = 0`, it lives in the **null-space of
`ρ_split`**: the trace statistic is literally blind to it, so the reset/anchor
channels **cannot perturb the locked trace result** — that is algebra, not a
promise. Only the **trace channel is verified**; anchor / reset / reorganize are
**new, unverified** functionals and stay labelled so until they clear their own
nulls. This is the first **cut-free** version of the cross-phase taxonomy: it
*decomposes the cophenetic signal* and never re-introduces the retired
tree-cutting / Jaccard-subtree measures (whose reliability is unknown).

---

## 0. Critical preamble (5-point, mandatory)

1. **Claim.** The cross-phase cophenetic signal partitions into anchor / trace /
   reset / reorganize; and beyond the already-established trace, the **reset
   (excursion) channel carries a cohort-coherent, anatomically-localized
   task-transient component that is distinct from the trace** (different system
   than β→OFC, or an honest null).
2. **Null(s).** *Temporal:* a **coupled-cross-phase strength-preserving** surrogate
   (shared 4-cycle swap sequence across phases) — preserves each phase's node
   strength *and* the shared backbone, randomizes the specific cross-phase fine
   structure. *Spatial:* the **implant-geometry baseline** (count-matched random
   leaf draw; same-probe / within-DK-system concentration).
3. **Strongest alternative the null must kill.** The channels are **strength /
   anatomy artifacts**: (a) anchor = high-strength hub backbone that trivially
   merges early at every phase (the KC-β confound); (b) reset = generic
   task-vs-rest deflection reproducible by independent per-phase strength
   evolution; (c) reorganize = sampling noise; (d) the *independent-per-phase*
   matched-strength null **trivially confirms** anchor/reset because real
   same-patient phases are intrinsically correlated (RPre≈RPost, both rest).
4. **Does the null control it — by mechanism?** The **coupled** null preserves
   per-phase strength *exactly* and keeps the shared perturbation locations, so a
   channel that beats it is beyond *both* strength and the shared backbone — it
   defuses (a),(b),(d). It **cannot** by itself rule out that the task-transient
   is a generic arousal/time effect rather than task-content (we have no
   alternate-task null — see §5, §11). The geometry baseline defuses the "it's
   just where the electrodes are" alternative for the spatial claim only.
   **Validation gate:** the coupled surrogate's cross-phase ρ^coph must bracket
   the observed value (well above the independent-per-phase null). If it does not,
   the coupled null is mis-specified and reset/anchor revert to descriptive-only.
5. **Falsification / limitations.** Reset fails if it clears no band/system under
   the coupled null + geometry baseline (→ descriptive-only). The whole
   decomposition is *wrong/buggy* if the **trace share does not peak at β/α** or
   does not reproduce the locked `obs_rho` bit-exactly. Hard limits that remain:
   no per-pair labels (quantization, §5/§6); **emergent is not testable**
   (cut-free, membership needs cuts); task-*content* specificity is unfalsifiable
   with one task.

---

## 1. Notation

- `p ∈ COHORT_N10 = {Pat_02,03,05,06,07,08,10,13,14,15}`; band
  `b ∈ {δ,θ,α,β,γ_l,γ_h}`; phases `φ ∈ {pre, task, post}` (rest_pre, task_test,
  rest_post). Split baseline: `A,B` are the two independent halves of rest_pre.
- `N_p` contacts; `M_p = N_p(N_p−1)/2` unordered pairs; pair index `u`, endpoints
  `(i,j)`.
- `D^φ ∈ ℝ^{M_p}` — **cophenetic (ultrametric) condensed vector** of phase `φ`:
  FC matrix → LRG propagator at `τ = 1/λ_max` → communication distance
  `T_ρ = 1/ρ̂` → average-linkage → `scipy.cluster.hierarchy.cophenet`. This is the
  *exact* `a63.lrg_ultrametric_condensed` used by the locked trace
  (`audit_63`), reused verbatim — no fork. Quantized to ≤ `N_p−1` merge heights
  ⇒ ~98% tied (see `perpair_cophenetic_trace_visualization`).
- Per-pair phase triple: `𝐃_u = (D¹_u, D²_u, D³_u) := (D^A_u, D^task_u, D^post_u)`.
  (For the *trace* statistic we keep the cross-split `Δrest` against `D^B`; see §2.)
- Mean `m_u = (D¹_u+D²_u+D³_u)/3`; centered `(d̃¹,d̃²,d̃³) = 𝐃_u − m_u`.
- Validated axes (unchanged from `audit_63/76`):
  `Δtask_u = D^task_u − D^A_u`, `Δrest_u = D^post_u − D^B_u`,
  `Δbase_u = D^A_u − D^B_u`; noise scale
  `σ_{p,b} = sqrt(mean_u Δbase_u²)` (RMS). **[impl note]** the draft said
  `median_u |Δbase_u|`; under quantization ~98% of `Δbase` are exact ties (0) so
  the median is 0 and would clamp σ — RMS is robustly non-zero and is the natural
  amplitude scale for `var` (an energy).
- `ρ_split(p,b) = Spearman(Δtask, Δrest)` — the **locked trace statistic**;
  per-pair rank-concordance carrier `s_u` per `audit_83.concordance`.
- Anatomy: `sys(i) ∈` `regions.ANATOMICAL_SYSTEMS` (9) /
  `regions.SUPERSYSTEMS` (11); `probe(i)` the sEEG shaft; `(x,y,z)`, DK region
  from `implant_pat_NN.csv` via `regions.load_channel_regions` (never letter
  prefixes — `feedback_implant_anatomy_not_letters`).

## 2. Definitions

### 2.1 Orthonormal contrast basis (the decomposition)

On the centered triple, define two orthonormal contrasts:

```
φ₁(u) = (D³_u − D¹_u)/√2            persistence   (RPost vs RPre)
φ₂(u) = (2 D²_u − D¹_u − D³_u)/√6   excursion     (Task vs rest-average)
```

**Exact energy split** (proof: `φ₁²+φ₂² = d̃¹²+d̃²²+d̃³² = 3·var(𝐃_u)`):

```
var(𝐃_u) = (φ₁(u)² + φ₂(u)²) / 3
```

The cross-phase fluctuation of every pair is, with **zero residual**, persistence
energy `φ₁²` plus excursion energy `φ₂²`. `φ₁ ⟂ φ₂`.

### 2.2 The four channels (per pair, dimensionless in σ units)

```
anchor   w^A(u) = exp( − var(𝐃_u) / (2 σ²) )          ∈ (0,1], ≈1 ⟺ rigid
trace    s(u)   = audit_83 rank-concordance of (Δtask, Δrest)   [VALIDATED]
                  cell statistic = ρ_split(p,b)  (UNCHANGED)
reset    w^R(u) = (φ₂(u)/σ) · 1[ |Δtask_u| > σ ] · 1[ |Δrest_u| ≤ σ ]   sign-aware
reorganize       residual mover energy not reproduced by the coupled null (§2.4)
```

- **Consistency theorem (trace protection).** `Δrest_u = √2·φ₁(u) + Δbase_u`, i.e.
  the `ρ_split` "rest" axis is the persistence contrast plus baseline noise; and
  for a perfect reset (`D³=D¹`) `φ₁=0 ⇒ Δrest=Δbase` (noise only) ⇒ the pair
  contributes nothing to `ρ_split`. **Reset ⟂ trace and lies in the null-space of
  the trace statistic.** Computing `w^R` cannot change `ρ_split`. ∎
- **trace is φ₁-dominant, not pure φ₁** (a→b=c gives φ₁/φ₂ = √3), so we **keep
  `ρ_split` itself as the trace statistic** and only *add* the orthogonal `φ₂`
  channel; we never reconstruct the trace from φ₁.

### 2.3 Cohort decomposition (per patient × band → cohort distribution)

```
anchor mass     A_{p,b} = mean_u w^A(u)
mover set       𝓜 = { u : var(𝐃_u) > κ σ² }            (κ free, swept; default 1)
trace energy    E_T = Σ_{u∈𝓜} φ₁(u)²
reset energy    E_R = Σ_{u∈𝓜} φ₂(u)²
shares          (Ŝ_T, Ŝ_R) = (E_T, E_R)/(E_T+E_R)
reorganize      O_{p,b} = fraction of (E_T+E_R) within the coupled-null band (§2.4)
```

Report the **cohort distribution** of `(A, Ŝ_T·(1−O), Ŝ_R·(1−O), O)` per band —
never a pooled consensus scalar (`feedback_no_partition_metrics…`,
"never pool into a consensus scalar"). Guards: `Ŝ_T` peaks at β/α; `A` dominates
(Gratton backbone check).

### 2.4 Coupled-cross-phase null (temporal; the Full-scope construction)

```
draw one ordered quadruple list Q = [(a,b,c,d)_1 … (a,b,c,d)_{n_swaps}]   # SHARED
for φ in {A, B, task, post}:
    W̃^φ = strength_preserving_shuffle(W^φ, swaps=Q)   # same locations, feasible δ_φ
recompute D̃^φ → φ̃₁, φ̃₂, channel scores ; repeat R times (new Q per r)
```

Preserves per-phase node strength **exactly** and the shared perturbation
locations (⇒ shared backbone). One-sided rank p as in `audit_83.run_band`
(`p = (1+#{surr ≥ obs})/(R+1)`). **Validation gate (run first):**
`ρ^coph(W̃^A,W̃^B)` etc. must bracket observed and exceed the independent-per-phase
null; else fall back to geometry baseline only.

### 2.5 Geometry baseline (spatial; the kept KC-era borrow — combinatorial, not a tree measure)

For a per-channel per-system aggregate, the spatial null is a **count-matched
random leaf draw** from the implant: is the channel's endpoint set more
concentrated within `sys`/`probe` than `K`-of-`N_p` random contacts? Pure implant
combinatorics (`f_sp^baseline(p)=Σ_q C(n_q,2)/C(N_p,2)` per `audit_60`), **never
touches a dendrogram**. Loose sanity cross-check only against the retired KC
numbers (β anchor 6.22× same-probe).

## 3. Properties

- **Ranges.** `w^A∈(0,1]`; `φ₁,φ₂∈ℝ`; shares `∈[0,1]`, sum to 1 with `A`+`O`.
- **Exactness.** `var = (φ₁²+φ₂²)/3` is an identity (no cross terms, no residual).
- **Orthogonality / invariance.** `φ₁⟂φ₂`; both invariant to the per-pair mean
  `m_u` (anchor level) ⇒ the fluctuation channels do not see the DC backbone, and
  anchor does not see the fluctuation. Scale-equivariant in σ.
- **Trace protection (identifiability).** `ρ_split` is a function of φ₁ + noise
  only (§2.2); the reset channel is φ₂, orthogonal ⇒ the verified result is
  invariant under adding the new channels.
- **What it CANNOT detect (negative properties).**
  - **No per-pair labels.** Quantization (98% ties) ⇒ channels are only
    trustworthy as **cell/system aggregates**, never as "edge u is a reset edge"
    (`audit_76` finding: the trace is a diffuse bulk drift, not a separable
    subset).
  - **No emergent channel.** Emergent is a community-*membership* statement; it
    needs tree cuts and is excluded by the cut-free mandate
    (`terminology.md`: "edge-rank distances cannot test for emergent modules").
  - **No task-content specificity.** With a single task we cannot separate
    task-driven reset from arousal/time-on-task.
- **Complexity.** Per cell `O(M_p)` for the channels; null is `R` eigendecompositions
  per phase (reuse cached `load_or_compute_surrogate_eigs`); localization `O(M_p)`
  per system.
- **Literature placement.** anchor↔reorganize is the cophenetic-hierarchy,
  cross-state analog of network **flexibility / allegiance / recruitment /
  integration** (Mattar 2015; Braun 2015) — but with a **strength-matched null**
  that literature lacks; **anchor dominance** is the Gratton 2018 "stable backbone
  ≫ task variation" prediction tested in sEEG/LRG; **trace** sits in the
  post-task-persistence lineage (Tambini & Davachi 2013); **trace vs reset** is a
  **hysteresis** loop (Kim/Lee/Mashour 2018) — hysteretic vs reversible branch of
  the RPre→Task→RPost cycle.

## 4. Caveats & failure modes

| | Caveat | Mitigation |
|--|--|--|
| 1 | **Quantization** (≤N−1 heights, ~98% ties) ⇒ per-pair scores are coarse. | Aggregate to cell/system only; never per-pair labels. Mirror `audit_76`. |
| 2 | **Independent-per-phase null trivially confirms anchor/reset** (real RPre≈RPost). | Use the **coupled** null for anchor/reset temporal tests; validate the gate (§2.4). |
| 3 | **anchor = strength backbone** tautology (high-strength pairs merge early everywhere). | Coupled null preserves strength ⇒ anchor must beat strength; report strength corr; expect mostly confirmatory. |
| 4 | **reorganize = noise floor** by construction (null-residual). | Label it a *residual/null bin*, never a "finding". |
| 5 | **κ and σ-disk radius are free knobs.** | Sweep; report sensitivity. **No pre-registered acceptance gate** (`feedback_no_pre_registered_acceptance`). |
| 6 | **Coupled null may over- or under-preserve** cross-phase correlation. | Validation gate (§2.4) run *before* any p-value; fall back to geometry baseline if it fails. |
| 7 | **A/B baseline noise** enters `Δrest` (`+Δbase`). | Identical to the locked pipeline; σ scales it out; bit-exact `obs_rho` guard. |
| 8 | **Multiplicity** across systems×channels×bands. | BH-FDR over **systems within one band**, per the localization lesson (`localization_audit_plan_2026_05_29`). Channels reported separately. |
| 9 | **Probe pseudo-replication** in localization. | `audit_83 --shaft-collapse` + LOO-patient robustness required before any claim. |

## 5. Pseudocode

```
# ---------- per (patient p, band b) ----------
load W^A, W^B, W^task, W^post           # audit_63.load_phase_fc (half cache)
D1,D2,D3, DB = lrg_ultrametric_condensed(each)          # a63 verbatim
assert spearman(D2-D1, D3-DB) == locked obs_rho[p,b]    # BIT-EXACT trace guard
m   = (D1+D2+D3)/3
phi1= (D3-D1)/sqrt2 ;  phi2=(2*D2-D1-D3)/sqrt6
var = (phi1**2+phi2**2)/3
sig = median(abs(D1-DB))
wA  = exp(-var/(2*sig**2))
movers = var > kappa*sig**2
E_T = sum(phi1[movers]**2) ; E_R = sum(phi2[movers]**2)
record A=mean(wA), share_T=E_T/(E_T+E_R), share_R=E_R/(E_T+E_R)
# channel per-pair scores for localization
s_trace = audit_83.concordance(D2-D1, D3-DB)            # VALIDATED carrier
s_reset = (phi2/sig) * (abs(D2-D1)>sig) * (abs(D3-DB)<=sig)
s_anchor= -var/sig**2

# ---------- coupled temporal null ----------
for r in 1..R:
    Q = draw_quadruples(seed=r)                         # SHARED across phases
    for phase: Wt = strength_preserving_shuffle(W, swaps=Q)  # same locations
    recompute phi1,phi2 -> E_T^r, E_R^r, s_*^r, rho_coph^r
validate: median_r rho_coph(A,B) brackets observed ; > independent-null
p_channel = (1 + #{surr >= obs}) / (R+1)                # one-sided, audit_83 style

# ---------- spatial localization (per channel) ----------
for channel s in {s_reset, s_anchor}:                   # s_trace already done (OFC)
    unit_mean = audit_83.unit_means_from_s(s, sys, shaft_collapse=True)
    p_geom  = geometry_baseline_p(unit, channel)        # count-matched draw
    p_temp  = coupled_null_p(unit, channel)             # from R surrogates
    keep unit iff BH(p_geom)<q AND BH(p_temp)<q within band b
require: LOO-patient + shaft-collapse robust
```

## 6. Visualization spec

Activate `use_lrg_style()`; **PDF only, full vector, no suptitle, no watermark,
transparent=True, figure-level legends/colorbars** (plotting rules). Class colors
(avoid near-white — `feedback_no_near_white_cmaps`): anchor `#807dba`, trace
`#2c7fb8`, reset `#41ab5d`, reorganize `#fc9272`.

1. **Decomposition figure.** Per band, stacked composition of cohort-median
   `(A, share_T·(1−O), share_R·(1−O), O)` + per-patient spread (strip/box).
   Reading: anchor dominates everywhere; trace wedge peaks β/α; reset wedge is the
   question; reorganize = residual.
2. **Contrast-plane figure.** Per patient, `(Δtask/σ, Δrest/σ)` cloud (the
   `audit_76` portrait) with the φ₁ (diagonal=trace) and φ₂ (Δtask-axis=reset)
   directions overlaid; symlog, sign-colored. ≥3 patients shown.
3. **Reset anatomy figure.** Only if §2.4 gate passes: DK-system bars of reset
   `unit_mean` with geometry+coupled-null significance, β next to the OFC trace
   panel to show dissociation (or not).

## 7. Connection to prior tools

| Prior tool | Relation |
|--|--|
| `audit_63` (`load_phase_fc`, `lrg_ultrametric_condensed`, `obs_rho`) | **Source pipeline, reused verbatim.** Trace bit-exact guard. |
| `audit_76` (pair displacement) | This **supersedes/extends** its cell-level readings with the orthogonal-contrast decomposition + localization; inherits its "no per-pair hard partition" lesson. |
| `audit_83` (`obs_trace`,`concordance`,`unit_means_from_s`,`--shaft-collapse`) | **Localization engine, signal-agnostic** — swap input signal to `s_reset`/`s_anchor`. |
| `audit_92` + `matched_strength.py` | Surrogate eig cache + `strength_preserving_shuffle` → base for the **coupled** null (shared `Q`). |
| `regions.py` (`ANATOMICAL_SYSTEMS`,`SUPERSYSTEMS`,`load_channel_regions`) | DK anatomy for localization. |
| **Retired KC taxonomy** (`audit_47–52/60`, `2026-05-08_taxonomy_mutual_exclusivity`) | **NOT revived.** Unknown reliability; tree-cutting. Used only as a loose numeric sanity cross-check; the **geometry baseline** (`audit_60`'s `f_sp`) is borrowed as a *combinatorial anatomy null*, not a tree measure. |
| Locked trace verdict (`localization_audit_plan_2026_05_29`) | The trace channel = that result; this work sits *around* it, does not alter it. |

## 8. Implementation plan

- **New audit:** `scripts/01_compute/audit/audit_105_cross_phase_taxonomy.py`
  (decomposition + cohort shares + bit-exact trace guard). Imports `audit_63`
  pipeline verbatim; **no forks** (`audit` skill check before commit).
- **Coupled null:** add `coupled=True` path to the surrogate driver (shared `Q`),
  living next to `strength_preserving_shuffle` in
  `src/lrg_eegfc/utils/surrogate/matched_strength.py` (≥2 callers ⇒ library, not
  script). Cache under `data/cache/coupled_surrogate_lrg/`.
- **Localization:** reuse `audit_83` with `--signal {reset,anchor}` flag
  (parametrize its input `s`); outputs `data/audit/cross_phase_taxonomy/`.
- **Stats:** `lrg_eegfc.utils.metrics.hypothesis` (`wilcoxon_z`, `bh_fdr`,
  `boot_ci_mean`) — no private copies.
- **Figures:** `scripts/02_preprint/preprint_*` only if the explore phase warrants;
  default lands figures under `data/audit/cross_phase_taxonomy/figures/`.
- **Order:** (i) decomposition + trace guard + shares; (ii) coupled-null
  **validation gate**; (iii) reset/anchor localization *iff* gate passes;
  (iv) evaluate with user → decide preprint §5.x. No preprint commitment up front.

## 9. Open questions

1. **Coupled-null fidelity.** Does the shared-`Q` construction bracket observed
   cross-phase ρ^coph, or over-/under-preserve it? Alt: baseline-anchored drift
   (perturb a common `W0=rest_pre` by drifts matched to `‖W^φ−W0‖`). Decide by the
   validation gate, not a priori.
2. **Reset operationalization.** `φ₂` (excursion) vs an explicit "return"
   coefficient `1 − |Δrest|/|Δtask|`. φ₂ chosen for orthogonality+exact energy
   split; cross-check the two correlate before committing.
3. **κ / σ-disk.** Sweep `κ∈{0.5,1,2}`; report share sensitivity. No gate.
4. **Anchor null.** Is the geometry baseline + coupled null enough, or does anchor
   need an explicit strength-regression control (expected: anchor ≈ hubness)?
5. **Hysteresis quantifier.** Whether to add a scalar loop-area
   (`∮` over RPre→Task→RPost) as the physics-facing summary of trace-vs-reset, or
   keep the energy shares. Defer.
