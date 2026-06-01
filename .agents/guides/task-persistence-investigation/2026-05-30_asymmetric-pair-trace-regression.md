---
name: asymmetric-pair-trace-regression
type: scope
era: COHORT_N10
status: current
created: 2026-05-30
updated: 2026-05-30
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-26_continuous-trace-matrix.md
  - .agents/preprint/bands/01_beta.md
  - .agents/preprint/methods/methods_revision_2026-05-18_cophenet.md
  - scripts/01_compute/audit/audit_63_split_baseline_surrogate.py
  - scripts/02_preprint/preprint_05_allbands_matched_strength_raw_D.py
  - src/lrg_eegfc/utils/surrogate/matched_strength.py
  - src/lrg_eegfc/utils/metrics/hypothesis.py
  - src/lrg_eegfc/utils/metrics/tree.py
---

# Asymmetric pair-trace regression slope `s_TR` — comparative methodology audit on three LRG primitives

## Renormalization head

The current §5.3 / §3.2 manuscript measure `ρ_split^coph` is a **rank-symmetric** Spearman between `Δ_task` and `Δ_rest` on the cophenetic distance `D_coph`; it strips magnitude and direction-orders the two phase-difference vectors as if they were exchangeable. The asymmetric through-origin regression slope `s_TR = ⟨Δ_task, Δ_rest⟩ / ‖Δ_task‖²` — the OLS coefficient when post-shift is regressed on task-shift through the origin — preserves both magnitude and the task→rest direction-of-causation and gives the verbatim physical reading "the fraction of the task-induced per-pair shift that is recovered in post". This scope defines `s_TR` (and its rotation-invariant companion `R² = ρ_Pearson²`) alongside the current `ρ_Spearman` and `ρ_Pearson` across **three per-pair primitives** — raw density `ρ̂(τ_max)`, raw distance `D(τ_max) = 1/ρ̂`, and cophenet `D_coph = cophenet(UPGMA(D(τ_max)))` — and prescribes a **comparative methodology audit**, NOT a manuscript-replacement. The output is a 4-measure × 3-primitive verdict heatmap under matched-strength surrogacy; the verdict decides whether `s_TR` becomes a candidate methodological supplement or replacement for `ρ_split^coph` in a future revision.

> **Naming locked.** The slope is `s_TR` (math, `s_{\mathrm{TR}}`) / `slope` (code, CSV columns) / "asymmetric pair-trace slope" or "recovery slope" (prose). The bare letter `β` is reserved for the 13–30 Hz band and is NEVER used for the regression slope. See `.agents/guides/04_rules/never-always-list.md` and `feedback_no_beta_for_regression_slope`.

## 1. Critical preamble (5-point, mandatory per `feedback_critical_null_preamble.md`)

1. **Claim.** Under matched-strength surrogacy and split-baseline within-rsPre halves, the asymmetric per-pair slope `s_TR` on at least one of the three LRG primitives `{ρ̂(τ_max), D(τ_max), D_coph}` yields a cohort-paired one-sided Wilcoxon verdict structurally similar (in band selectivity, effect-size ratio, and patient-level agreement) to the current `ρ_split^coph` at the β band, while providing strictly more interpretable units ("fraction of task imprint recovered") and an asymmetric framing that pins task as the cause and post as the response.
2. **Null tested.** (i) Independent split-baseline halves rsPre_A and rsPre_B used to construct `Δ_task = X_task − X_{rsPre,A}` and `Δ_rest = X_rsPost − X_{rsPre,B}` — eliminates shared-baseline correlation between the two delta-vectors by design; (ii) matched-strength surrogate ensemble R=200 per (patient, band, phase) with seed 20260511, swap_factor 20, 4-cycle ±δ Maslov-Sneppen rewiring preserving each node's strength to 10⁻⁴ — eliminates the per-node strength-evolution alternative; (iii) per-patient observed-vs-surrogate one-sided p (`mean(surrogate >= observed)`) plus cohort-paired Wilcoxon (`observed > surrogate-median`); (iv) LOO sensitivity at the cohort level for descriptive single-patient-driver flag (NOT a gate, per `feedback_no_single_patient_p_driven.md`).
3. **Strongest alternatives.** (a) Per-node strength evolution alone (not pair-level reorganisation) reproduces the cohort `s_TR` direction; (b) a within-session monotone drift in signal quality inflates both Δ_task and Δ_rest with the temporal gap and produces `s_TR > 0` without any task-induced imprinting (the third-common-cause confound: `s_TR = 1`, `R² = 1` is achievable with NO actual task-driven trace); (c) for `s_TR` specifically, magnitude-asymmetric outlier pairs (a single Δ_task = high-magnitude tail) can dominate the inner product `⟨Δ_task, Δ_rest⟩` and the denominator `‖Δ_task‖²` differently than they dominate the Spearman ρ on ranks.
4. **Does the null cover it.** (i) Split-baseline eliminates shared-baseline noise correlation by construction. (ii) Matched-strength eliminates (a) — the surrogate ensemble preserves strength but destroys per-pair edge identity, so any pure-strength-driven `s_TR` reproduces in surrogates and gates out. (iii) Drift-floor null (cf. `audit_63_split_baseline_surrogate.py`) addresses (b) for `ρ_split^coph` but the **drift null for `s_TR` is NOT the drift null for ρ_Spearman** — the drift floor must be recomputed on `s_TR` using rsPre_A/B and rsPost_A/B halves (deferred to v2 if v1 verdict warrants; in v1 we report matched-strength only). (iv) The matched-strength null does **NOT** cover the magnitude-asymmetric outlier case (c) — a single high-leverage pair will dominate `s_TR` and `‖Δ_task‖²` and the matched-strength surrogate may also exhibit this, so the per-patient `s_TR` value should be inspected against the per-pair scatter for influential-point leverage (visualisation hook in §8). v1 will compute per-patient `s_TR` Cook's-distance-style robustness as a descriptive flag.
5. **Falsification + limitations.** Falsification: (i) `s_TR`-on-D_coph cohort-paired Wilcoxon p ≥ 0.05 at the β band (the band where `ρ_split^coph` is strongest); (ii) `s_TR`-on-D_coph β-band ratio-over-surrogate-median ≤ 5× (vs `ρ_split^coph`'s 23.7×); (iii) `s_TR` LOO max p > 0.05 at the β band with a single argmax patient (single-patient-driven verdict). Limitations: (a) `s_TR` is **not unit-free across primitives** — ρ̂, D, and D_coph have different dynamic ranges, so `s_TR` values are comparable within a primitive but not across; reading across the primitive axis must use ratio-over-surrogate or `R²` for normalised comparison. (b) The third-common-cause confound (point 3b) cannot be ruled out by any null defined on the same four phases; a no-task control session in the same patients would be required and is unavailable. (c) `R² = ρ_Pearson²` uses mean-centred Pearson while `s_TR` is through-origin — they answer slightly different geometric questions (cosine of mean-centred vectors vs slope of through-origin regression); §3 makes this explicit.

## 2. Notation

| Symbol | Domain | Meaning |
|:---|:---|:---|
| `p ∈ P` | `|P| = 10`, `lrg_eegfc.config.const.PATIENTS_4PHASE` | patient |
| `b ∈ B` | `|B| = 6` (δ, θ, α, β, γ_l, γ_h) | band |
| `φ ∈ Φ_split = {rest_pre_A, rest_pre_B, task_test, rest_post}` | 4 phases | split-baseline phase set |
| `N_p` | int | number of contacts for patient `p` (∈ [113, 122]) |
| `m_p = N_p (N_p − 1) / 2` | int | number of upper-triangle pairs |
| `iu(p)` | `m_p` index pairs | `np.triu_indices(N_p, k=1)` |

Per `(p, b, φ)` the LRG cache (`data/cache/imcoh_lrg/Pat_NN/{band}_{phase}_lrg_imcoh-abs.npz` for full phases; `data/cache/imcoh_lrg_halves/...` for the `rest_pre_{A,B}` splits) provides the eigendecomposition `(λ^(p,b,φ), V^(p,b,φ))` of the combinatorial Laplacian `L̂ = D̂ − W`. The three per-pair primitives are reconstructed deterministically from `(λ, V)`:

```
τ_max = 1 / λ_max
ρ̂(τ_max) = (V · diag(exp(−τ_max · λ)) · V^T) / Σ_k exp(−τ_max · λ_k)         ∈ ℝ^{N × N}
D(τ_max) = 1 / ρ̂(τ_max)                                                     ∈ ℝ^{N × N}_{>0}
D_coph    = cophenet( UPGMA( squareform(D(τ_max)) ) )                       ∈ ℝ^{N × N}_{≥0} (ultrametric)
```

We collectively denote any of the three by `X^{(p,b,φ)} ∈ ℝ^{N × N}`. Upper-triangle vectorisation: `vec(X^{(p,b,φ)}) := X^{(p,b,φ)}[iu(p)] ∈ ℝ^{m_p}`.

Per-pair phase-difference vectors (split-baseline):

```
Δ_task^{(p,b,X)} := vec(X^{(p,b,task_test)}) − vec(X^{(p,b,rest_pre_A)})    ∈ ℝ^{m_p}
Δ_rest^{(p,b,X)} := vec(X^{(p,b,rest_post)}) − vec(X^{(p,b,rest_pre_B)})    ∈ ℝ^{m_p}
```

Note: `rest_pre_A` and `rest_pre_B` are independent halves of rsPre by construction (`compute_imcoh_abs_halves` with `nperseg_for_fs(fs)//2`), so the two delta-vectors do not share a baseline noise term. This is the same convention as `audit_63_split_baseline_surrogate.py` and `preprint_05_allbands_matched_strength_raw_D.py`.

## 3. Definitions

For paired vectors `x, y ∈ ℝ^m`, define the four measures of cross-phase agreement:

### 3.1 Symmetric rank correlation (current measure)

```
ρ_Spearman(x, y) := Spearman( rank(x), rank(y) ) ∈ [−1, +1]
```

Rank-symmetric in `x ↔ y`. Throws away both magnitude and the absolute level of the deltas. This is `ρ_split^coph` when `(x, y) = (Δ_task, Δ_rest)` with `X = D_coph`, the current §5.3 / §3.2 quantity.

### 3.2 Symmetric magnitude-aware correlation

```
ρ_Pearson(x, y) := ⟨x − μ_x, y − μ_y⟩ / (σ_x · σ_y)   ∈ [−1, +1]
```

with `μ_x = mean(x)`, `σ_x² = ‖x − μ_x‖² / m`. Mean-centred. Symmetric in `x ↔ y`. Magnitude-aware but invariant under positive rescaling of either argument.

### 3.3 Asymmetric through-origin regression slope (new headline)

```
s_TR(x, y) := ⟨x, y⟩ / ⟨x, x⟩   ∈ ℝ
```

This is the OLS solution of `y ≈ s_TR · x` over the through-origin model (no intercept). **Asymmetric in `x ↔ y`** — the inputs play the roles of cause (`x = Δ_task`) and response (`y = Δ_rest`). Units of `s_TR` = (units of `y`)/(units of `x`). When `x` and `y` are both delta-distances on `D_coph`, `s_TR` is dimensionless and reads as **the fraction of the task-induced per-pair shift that is preserved in rsPost**:
- `s_TR = 1` ⇔ post shifts match task shifts pair-by-pair (full trace at the linear-model level)
- `s_TR ∈ (0, 1)` ⇔ partial recovery; e.g. `s_TR = 0.3` means rsPost recovered ~30 % of the per-pair shift task induced
- `s_TR = 0` ⇔ reset (post-shifts independent of task-shifts in the linear-model sense)
- `s_TR < 0` ⇔ anti-trace (post-shifts move opposite to task-shifts)
- `s_TR > 1` ⇔ amplification (post-shifts exceed task-shifts on the same direction; uncommon, biologically suspicious)

### 3.4 Rotation-invariant fit fraction (companion)

```
R²(x, y) := ρ_Pearson(x, y)²   ∈ [0, 1]
```

The fraction of `y`'s variance (mean-centred) explained by a linear fit `y ≈ a + s_centred · x`. Companion to `s_TR` but with a subtle distinction: `R²` here is mean-centred ("how aligned are the two delta-vectors as directions after centring"), whereas `s_TR` is through-origin. The two together fully characterise the linear-model trace at the cohort level; either alone is incomplete.

(Open question §11: whether to report a through-origin `R²_TR := ⟨x, y⟩² / (⟨x, x⟩ · ⟨y, y⟩) = ρ_uncentred²` as a fifth measure that is the geometrically natural companion to `s_TR`. v1 follows the other agent's recipe and reports the centred `R²` only; the through-origin variant is deferred unless v1 verdict warrants.)

### 3.5 Per-leaf decomposition of `s_TR`

For each contact `c ∈ {1, …, N_p}`, restrict the inner product and denominator to pairs that involve `c`:

```
s_TR^{(p, b, X)}(c) := ( Σ_{j ≠ c} Δ_task^{(p,b,X)}_{c,j} · Δ_rest^{(p,b,X)}_{c,j} )
                       / ( Σ_{j ≠ c} ( Δ_task^{(p,b,X)}_{c,j} )² )
                  ∈ ℝ
```

This is the **row-`c` restriction** of the global `s_TR`. Averaging `s_TR(c)` over `c` lands near (NOT exactly at) the global `s_TR` because of the per-row weighting by `Σ_j Δ_task²`. The per-leaf scalar admits the same "fraction recovered" reading as the global statistic, restricted to the contact's row of pair distances. It powers the per-leaf colored dendrogram visualisation (§8.3).

## 4. Properties

| # | Statement |
|:--|:---|
| 4.1 | Range / domain: `ρ_Spearman, ρ_Pearson ∈ [−1, 1]`; `R² ∈ [0, 1]`; `s_TR ∈ ℝ` (unbounded). |
| 4.2 | Asymmetry: `s_TR(x, y) ≠ s_TR(y, x)` in general (`s_TR(y, x) = ⟨x, y⟩ / ⟨y, y⟩`). The two are equal only when `‖x‖ = ‖y‖`. All three other measures are symmetric. |
| 4.3 | Mean-centring: `ρ_Pearson` and `R²` subtract means; `ρ_Spearman` is centring-irrelevant (ranks are centred by construction); `s_TR` does NOT subtract means — the through-origin interpretation requires the natural zero `Δ = 0` ⇔ "no shift", which would be destroyed by centring. |
| 4.4 | Anatomy subtraction: `Δ_task` and `Δ_rest` are subtraction-canonical (`X_task − X_{rsPre,A}` cancels the static anatomy component of `X` that appears in both phases). This is the same justification as in `2026-04-26_continuous-trace-matrix.md` §3.2. The cohort signal is on the *residual* per-pair shift. |
| 4.5 | Sensitivity to magnitude: `ρ_Spearman` is rank-only (insensitive to scale); `ρ_Pearson` is sensitive to scale of the centred vectors; `s_TR` and `R²` are both sensitive to magnitude but in geometrically distinct ways (`s_TR` to the inner product scaled by `‖Δ_task‖²`, `R²` to the cosine of the centred vectors). |
| 4.6 | Invariance under positive rescaling: `s_TR(α·x, y) = (1/α) · s_TR(x, y)` for `α > 0`; `s_TR(x, α·y) = α · s_TR(x, y)`. So `s_TR` is NOT rescale-invariant in either argument — units matter, comparability across primitives requires within-primitive ratio-over-null. |
| 4.7 | Identifiability: `s_TR` is the unique through-origin OLS coefficient given `‖Δ_task‖ > 0`. Undefined when `Δ_task ≡ 0` (no task-induced shift) — returned as NaN. |
| 4.8 | Complexity: all four measures `O(m_p)`. Cophenet recomputation is `O(N_p² log N_p)` per surrogate per phase, dominating runtime. |
| 4.9 | Reuses primitives: `load_or_compute_surrogate_eigs` from `lrg_eegfc.utils.surrogate.matched_strength`, `cophenet_matrix` from `lrg_eegfc.utils.metrics.tree`, `wilcoxon_z` and `loo_sensitivity` from `lrg_eegfc.utils.metrics.hypothesis`. One new helper: `regression_slope_through_origin(x, y) → (slope, r_squared)` to be added at the bottom of `hypothesis.py`. |
| 4.10 | NOT detected by `s_TR`: (a) emergent modules forming purely in rsPost (no task imprint at the pair level → `Δ_task ≈ 0` for the relevant pairs); (b) full-anatomy resets (`Δ_task` and `Δ_rest` both saturated but at opposite signs); (c) pair-level reorganisation that swaps within-`Δ_task` magnitudes but keeps the average — visible to per-leaf `s_TR(c)` and the per-pair scatter, invisible to the cohort-level scalar. |

## 5. Caveats & failure modes

### 5.1 Magnitude-asymmetric outlier dominance

`s_TR = ⟨Δ_task, Δ_rest⟩ / ‖Δ_task‖²` is sensitive to a single high-leverage pair in two ways simultaneously: that pair dominates the numerator AND inflates the denominator. The net effect on `s_TR` is NOT monotone in the leverage of a single pair. Mitigation: render the per-pair scatter (§8.2) and flag patients whose `s_TR` changes by > 30% under leave-one-pair-out (descriptive flag in `per_patient_per_band.csv`; v1 reports this only for the β band as a sanity check, full cohort if warranted).

### 5.2 Units mismatch across primitives

`s_TR` on ρ̂ is bounded ≈ ℝ (ρ̂ ∈ [0, 1] / N range, so `s_TR` is unbounded in principle but typically `O(1)`); `s_TR` on D is bounded by the dynamic range of the propagator distances (often `O(10²–10⁴)` on heterogeneous FC); `s_TR` on D_coph is bounded by the dendrogram merge-height range. Direct comparison of `s_TR` values across primitives is **meaningless**; only within-primitive ratio-over-surrogate-median, or the centred `R²` (rotation-invariant), are cross-primitive comparable. The comparison heatmap (§8.1) must enforce this in its colormap normalisation.

### 5.3 Pat_15 anti-anchor

Pat_15 is the cohort's anti-aligned anchor at the β LRG layer (right-hemisphere-only implant, 0 epi contacts; `.agents/preprint/bands/01_beta.md` §6). The cohort-paired Wilcoxon for `s_TR` will inherit this. LOO sensitivity will flag Pat_15 at the β band as expected; this is biology, not a methodological flaw.

### 5.4 Same-probe pair bias (deferred to v2)

Under `|ImCoh|` the zero-phase-lag component is killed by Nolte 2004 construction (`probe_bias_critical` memory), so the same-probe bias is muted relative to MSC. v1 does NOT compute a cross-probe-only variant; v2 adds it if v1 verdict warrants.

### 5.5 Third-common-cause confound (unfalsifiable in this design)

A confounder that drives both `Δ_task` and `Δ_rest` away from rsPre — e.g. arousal coupled to task block, hardware drift correlated with phase — yields `s_TR = 1, R² = 1` with NO real imprinting. The matched-strength null does NOT eliminate this. A no-task control session in the same patients would, but is unavailable. v1 reports `s_TR` with this limitation stated up front.

### 5.6 Bare letter `β` as naming hazard

The slope is `s_TR`, NEVER `β`. CSV columns: `*_slope` not `*_beta`. Figure axis labels: `s_TR` typeset, not `β`. Prose: "asymmetric pair-trace slope" or "recovery slope", never bare `β`. Locked 2026-05-30 per never-always-list.

## 6. Pseudocode

```
Input:
  P  = PATIENTS_4PHASE                  # n=10
  B  = ["alpha", "beta", "low_gamma", "delta", "theta", "high_gamma"]
  PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
  R = 200                               # matched-strength surrogates
  swap_factor = 20
  seed = 20260511                       # canonical
  fc_method = "imcoh_abs"

For each (p, b) ∈ P × B:
  # ---- observed primitives from cached eigendecompositions ----
  eigs_obs[φ] ← load_lrg_eig(p, b, φ) for φ ∈ PHASES         # 4 phases
  for φ in PHASES:
      ρ̂_obs[φ], D_obs[φ]  ← primitives_from_eig(eigs_obs[φ])  # N×N
      D_coph_obs[φ]        ← cophenet_matrix(linkage(squareform(D_obs[φ]), 'average'), condensed=False)
  N ← D_obs["rest_pre_A"].shape[0]; iu ← triu_indices(N, k=1)

  # ---- observed per-pair shift vectors per primitive ----
  for X in {ρ̂, D, D_coph}:
      dT[X] ← X_obs["task_test"][iu] − X_obs["rest_pre_A"][iu]
      dR[X] ← X_obs["rest_post"][iu] − X_obs["rest_pre_B"][iu]

  # ---- observed scalars (3 primitives × 4 measures = 12) ----
  for X in {ρ̂, D, D_coph}:
      obs[X]["spearman"]  ← spearmanr(dT[X], dR[X]).statistic
      obs[X]["pearson"]   ← pearsonr(dT[X], dR[X]).statistic
      obs[X]["slope"], obs[X]["rsq"] ← regression_slope_through_origin(dT[X], dR[X])

  # ---- surrogate ensemble per phase ----
  for φ in PHASES:
      W_φ ← load_FC_for_phase(p, b, φ)                          # halves or full
      evals_R[φ], evecs_R[φ] ← load_or_compute_surrogate_eigs(
          p, b, φ, W_φ, n_surr=R, swap_factor=swap_factor, seed=seed, fc_method=fc_method
      )                                                          # shapes (R, N), (R, N, N)

  # ---- per-surrogate scalars ----
  surr[X][m] ← empty(R) for X in {ρ̂, D, D_coph}, m in {spearman, pearson, slope, rsq}
  for r in 0..R-1:
      try:
          for φ in PHASES:
              ρ̂_r[φ], D_r[φ] ← primitives_from_eig(evals_R[φ][r], evecs_R[φ][r])
              D_coph_r[φ]    ← cophenet_matrix(linkage(squareform(D_r[φ]), 'average'), condensed=False)
          for X in {ρ̂, D, D_coph}:
              dT_r ← X_r["task_test"][iu] − X_r["rest_pre_A"][iu]
              dR_r ← X_r["rest_post"][iu] − X_r["rest_pre_B"][iu]
              surr[X]["spearman"][r] ← spearmanr(dT_r, dR_r).statistic
              surr[X]["pearson"][r]  ← pearsonr(dT_r, dR_r).statistic
              surr[X]["slope"][r], surr[X]["rsq"][r] ← regression_slope_through_origin(dT_r, dR_r)
      except (LinAlgError, NaN-in-eig):
          surr[X][m][r] ← NaN for all (X, m)

  # ---- emit per-patient-per-band row ----
  for X in {ρ̂, D, D_coph}, m in {spearman, pearson, slope, rsq}:
      finite ← surr[X][m][np.isfinite(...)]
      row[f"obs_{X}_{m}"]               ← obs[X][m]
      row[f"surr_p50_{X}_{m}"]          ← median(finite)
      row[f"surr_std_{X}_{m}"]          ← std(finite, ddof=1)
      row[f"obs_z_{X}_{m}"]             ← (obs[X][m] − mean(finite)) / std(finite)
      row[f"obs_p_one_sided_{X}_{m}"]   ← mean(finite >= obs[X][m])
  append row to per_patient_per_band

After cohort loop:
  for b in B, X in {ρ̂, D, D_coph}, m in {spearman, pearson, slope, rsq}:
      vec_obs  ← per_patient_per_band[b][f"obs_{X}_{m}"]
      vec_surr ← per_patient_per_band[b][f"surr_p50_{X}_{m}"]
      diff      ← vec_obs − vec_surr
      W, p_w    ← wilcoxon(diff, alternative='greater')
      n_above   ← (per_patient_per_band[b][f"obs_p_one_sided_{X}_{m}"] < 0.05).sum()
      loo       ← loo_sensitivity(diff, lambda v: wilcoxon(v, alternative='greater').pvalue)
      emit cohort_summary row(b, X, m, median(vec_obs), median(vec_surr),
                              median(vec_obs)/median(vec_surr), n_above, p_w,
                              loo.max_p, loo.argmax_patient)

Output:
  data/audit/pair_trace_measures_comparison/per_patient_per_band.csv   # 60 rows × 48 cols
  data/audit/pair_trace_measures_comparison/cohort_summary.csv         # 72 rows
  data/audit/pair_trace_measures_comparison/README.md                  # lab notes
```

### 6.1 Complexity

Per `(p, b)` cell: 4 phases × `O(N² log N)` cophenet × 201 realisations ≈ `O(800 · N² log N)`. At `N ≈ 115` that is `~ 10⁷` ops per cell. Total cohort `60 · 10⁷ = 6 · 10⁸` ops ≈ 60–90 min wall-clock on a single core.

### 6.2 NaN handling

A surrogate row whose `eigvals` or `eigvecs` contains NaN (rare; mostly numerical edge cases) is dropped — `finite = surr[np.isfinite(surr)]`. Per `feedback_brutal_honesty_no_sycophancy`, the script PRINTS the per-cell `n_surrogates_finite` and the report flags any cell with `n_surrogates_finite < 180` for inspection (< 90% surrogate yield).

## 7. Visualization spec

### 7.1 Cohort comparison heatmap (primary verdict figure)

`data/audit/pair_trace_measures_comparison/figures/comparison_heatmap.pdf`. PDF only, full vector.

Layout: 4 rows (measures: `ρ_Spearman, ρ_Pearson, s_TR, R²`) × 3 columns (primitives: `ρ̂, D, D_coph`). Each panel is a 6-band × n=10 grid of dots — y-axis = band (δ, θ, α, β, γ_l, γ_h), x-axis = patient ID, dot value = per-patient `(obs − surr_p50)` for that (measure, primitive). Dot colour = sign of difference on a saturated diverging palette (no `RdBu_r` if near-white interior wash-out; use `seismic` or custom palette per `feedback_no_near_white_cmaps`). Dot edge = darker if `obs_p_one_sided < 0.05`. **NO `n/10` reference lines** per `feedback_no_visual_threshold_lines`.

Below each panel: a single annotation line (no boxed text, no in-axes summary) giving per-band cohort-paired Wilcoxon p over the 6 bands, formatted `δ:p / θ:p / α:p / β:p / γ_l:p / γ_h:p`. Per `feedback_no_text_in_figures`, this annotation is the only text in the axes; full numbers go to the CSV and the report.

Reading rule: convergence at the β-band column on `D_coph` row (all four measures fire) supports the methodology migration. Divergence (e.g. `s_TR` fires on ρ̂ and D too where `ρ_Spearman` does not, indicating that `s_TR`-on-D is permissive in the same way raw D was permissive for `ρ_split^raw`) confirms `cophenet_methodology_rationale`.

### 7.2 Per-pair scatter at the β band (canonical-`s_TR` figure)

`data/audit/pair_trace_measures_comparison/figures/scatter_slope_betaband/<patient>.pdf`. One PDF per patient at β. PDF only, full vector.

Four panels per patient at β band:
- A — `Δ_task` heatmap (N × N), diverging palette, symmetric `vmin = −|max|, vmax = +|max|`, diagonal masked.
- B — `Δ_rest` heatmap, same palette and limits as A.
- C — sign-agreement `σ(i,j) = sgn(Δ_task) · sgn(Δ_rest)` heatmap (3-colour categorical: yellow `+1`, grey `0`, purple `−1`).
- D — `(Δ_task, Δ_rest)` scatter (cross-probe pairs only as the bulk of the cloud, same-probe pairs overlaid as a different marker so the reader can see the bias contribution). Through-origin regression line of slope `s_TR` drawn through (0, 0). Reference diagonal `y = x` drawn lightly. Annotation in the corner: `s_TR`, `R²`, `ρ_Spearman` numbers from the CSV. **NO `n/10` lines.** Per `feedback_no_text_in_figures`, this annotation is the only in-axes text.

Cohort overlay panel (separate PDF): all 10 patient regression lines overlaid on one axis with a thick patient-median line on top. Layout per `2026-04-26_continuous-trace-matrix.md` §6.2.

### 7.3 Per-leaf `s_TR(c)` dendrogram (β band, 4 exemplar patients)

`data/audit/pair_trace_measures_comparison/figures/per_leaf_slope_betaband/<patient>.pdf`. One PDF per patient ∈ {Pat_02, Pat_05, Pat_06, Pat_08} at β.

Per patient: rsPost dendrogram (linkage from `D_coph(rest_post)`) with leaves coloured by per-leaf `s_TR(c)` on a saturated diverging palette centred at 0 (avoid `RdBu_r` washing; prefer `coolwarm` or custom). Reuse the `link_color_func` pattern at `scripts/01_compute/audit/audit_47_kc_trace_network_view.py:340-372`. Subtree colour = parent merge-height-coloured by the subtree's mean `s_TR(c)` over leaves.

Y-axis: log-scale; `tmin = merge_heights[0] · 0.8`, `tmax = merge_heights[-1] · 1.05` per `feedback_dendrogram_limits`.

No suptitle per `feedback_no_suptitles`. Sidecar `.md` carries the caption.

## 8. Connection to prior tools

| Prior measure | Relation to this scope |
|:---|:---|
| `2026-04-26_continuous-trace-matrix.md` | Direct predecessor. Symmetric `ρ_Spearman` on `D` (single-scale, pre-cophenet). This scope generalises to 4-measure × 3-primitive grid and runs matched-strength on each cell, replacing the shared-baseline `ρ_Spearman` of the predecessor by the split-baseline version of all four measures. The predecessor's per-pair scatter visualisation is reused here for `Δ_task` vs `Δ_rest` (§7.2). |
| `audit_63_split_baseline_surrogate.py` | Matched-strength template for `ρ_split^coph` at α, β, γ_l. This scope's compute script `audit_73` extends the template to 4 measures × 3 primitives × 6 bands. The matched-strength surrogate cache (`data/cache/matched_strength_surrogate_lrg/`) is shared. |
| `preprint_05_allbands_matched_strength_raw_D.py` | All-bands raw-`D` `ρ_split` matched-strength run. This scope's compute script template (`audit_73`) is a direct extension of `preprint_05`'s structure, adding the cophenet primitive, the ρ̂ primitive, and the three new measures. |
| `audit_66/67` Grassmann subspace probe | Different probe family (subspace orientation, not per-pair geometry). Out of scope; the Grassmann analogue of `s_TR` (asymmetric subspace distance) is a v2 question. |
| KC tree distance (`audit_47/65`, retired 2026-05-18) | Retired probe family; not extended here. |

## 9. Implementation plan

### 9.1 Library (one helper)

`src/lrg_eegfc/utils/metrics/hypothesis.py` — append:

```python
def regression_slope_through_origin(x: np.ndarray,
                                   y: np.ndarray
                                   ) -> tuple[float, float]:
    """Through-origin slope and r-squared of y ≈ slope · x.

    slope     = <x, y> / <x, x>   (asymmetric in x ↔ y)
    r_squared = pearson(x, y) ** 2

    Returns (nan, nan) for empty / shape-mismatched / zero-norm inputs.
    """
```

Unit test at `tests/test_regression_slope_through_origin.py`:
- `slope(x, α·x) = α`, `r_squared = 1` for synthetic `α ∈ {-2, 0.5, 1, 3}`.
- `(slope, r_squared)` reproduces `numpy.linalg.lstsq` through-origin and `scipy.stats.pearsonr` to machine precision on random `(x, y)`.
- NaN / empty / zero-norm edge cases return `(nan, nan)`, not raise.

### 9.2 Compute script

`scripts/01_compute/audit/audit_73_pair_trace_measures_comparison.py`:

- Imports `PATIENTS_4PHASE` from `lrg_eegfc.config.const`, `load_or_compute_surrogate_eigs` from `lrg_eegfc.utils.surrogate.matched_strength`, `cophenet_matrix` from `lrg_eegfc.utils.metrics.tree`, `regression_slope_through_origin`, `wilcoxon_z`, `loo_sensitivity` from `lrg_eegfc.utils.metrics.hypothesis`.
- Phase tuple: `("rest_pre_A", "rest_pre_B", "task_test", "rest_post")`.
- Outputs:
  - `data/audit/pair_trace_measures_comparison/per_patient_per_band.csv` (60 rows × 48 numeric cols)
  - `data/audit/pair_trace_measures_comparison/cohort_summary.csv` (72 rows: band × primitive × measure)
  - `data/audit/pair_trace_measures_comparison/README.md` (lab notes; auto-generated)

### 9.3 Figure scripts

Three under `scripts/01_compute/figures_embedded/`, each calling `use_lrg_style()`:

1. `fig_pair_trace_comparison_heatmap.py` (§7.1)
2. `fig_pair_trace_scatter_slope_betaband.py` (§7.2)
3. `fig_per_leaf_slope_dendrogram_betaband.py` (§7.3)

### 9.4 Report

`.agents/reports/2026-05-30_pair-trace-measures-methodology-audit.md` — renormalisation head + 12-cell verdict heatmap + convergence/divergence cross-comparison + honest verdict + migration recommendation IF justified.

## 10. Open questions

1. **Through-origin R²_TR vs centred R²**. v1 follows the other agent's recipe (`R² = ρ_Pearson²`, centred). Geometrically the more natural companion to `s_TR` is `R²_TR = ⟨x, y⟩² / (⟨x, x⟩ · ⟨y, y⟩) = (uncentred Pearson)²`. Decision deferred to v2 if v1 verdict warrants — adds 12 more CSV columns.
2. **Cross-primitive `s_TR` comparability**. Within-primitive ratio-over-surrogate-median is the recommended comparison axis. The comparison heatmap (§7.1) normalises by this. Cross-primitive `R²` is rotation-invariant and IS comparable; cross-primitive `s_TR` is NOT.
3. **LOO at the patient level vs at the pair level**. v1 does patient-level LOO only (descriptive, not gating). Pair-level LOO (Cook's-distance-style for the regression) is a richer diagnostic for §5.1 (magnitude-asymmetric outlier dominance) but is deferred unless the per-patient scatter shows obvious leverage points.
4. **Drift-floor null for `s_TR`**. v1 reports matched-strength only. The drift-floor null on `s_TR` requires `s_TR( D_coph(rsPre_B) − D_coph(rsPre_A), D_coph(rsPost_B) − D_coph(rsPost_A) )` over 9 patients (no Pat_14 in H2e cache) and would land in v2.
5. **Manuscript framing**. If v1 verdict supports migration, the methods directive `.agents/preprint/methods/methods_revision_2026-05-30_asymmetric_slope.md` (new) crystallises the rename: `ρ_split^coph` becomes "rank-robust complement", `s_TR^coph` becomes "asymmetric pair-trace slope" as the headline. The current preprint cycle continues with `ρ_split^coph` until you explicitly approve the swap.

## 11. Status semantics

`draft` — written, code not yet implemented. Promote to `current` only after `audit_73` runs end-to-end, the (Pat_05, β, D_coph, ρ_Spearman) reproduction check passes within rounding, and the cohort_summary.csv is on disk. Promote to `superseded` only if v2 replaces this scope; in that case add a top-of-body pointer to the successor.

**Promoted to `current` 2026-05-30.** `audit_73` ran end-to-end (720 + 72 rows, 200/200 surrogates finite); the (Pat_05, β, `D_coph`, ρ_Spearman) reproduction passed (Δ = 4.6e-07 vs audit_63; raw-D bit-exact vs preprint_05). **Verdict: `s_TR` loses the β trace under matched-strength** (β 7/10 p=0.007 → 2/10 p=0.161 on `D_coph`) because magnitude-weighting reintroduces the strength-driven merge-height structure the null nullifies. Keep `ρ_split^coph`; `ρ_Pearson` is the defensible magnitude-aware companion; `s_TR` and `R²` are unsuitable as headlines. Full writeup: `.agents/reports/2026-05-30_pair-trace-measures-methodology-audit.md`.
