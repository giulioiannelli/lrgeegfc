---
name: coordinated-cross-phase-null
type: scope
era: COHORT_N10
status: draft
created: 2026-06-11
updated: 2026-06-11
pointers:
  - .agents/guides/task-persistence-investigation/2026-06-08_on-manifold-coherency-surrogate.md
  - src/lrg_eegfc/utils/surrogate/coherency_surrogate.py
  - src/lrg_eegfc/utils/surrogate/matched_strength.py
  - scripts/01_compute/audit/audit_97_crc_cohort_bracket.py
  - scripts/01_compute/audit/audit_63_split_baseline_surrogate.py
---

# Coordinated cross-phase null (SB-CRC) — a shared-backbone, deviation-randomized null that separates *task-specific* persistence from a *task-independent stable coupling backbone*

## Renormalization head

Both nulls we run today — matched-strength and CRC — randomize **each phase
independently**, so under the §5.3 split-baseline protocol both collapse to a
`ρ_split` floor of **≈ 0**. Clearing that floor proves only that *some*
cross-phase coupling-identity structure exists; it does **not** separate the two
explanations that structure could have: (i) **task-specific reorganization that
persists** into `rest_post` (the trace we claim), versus (ii) a **task-independent
stable coupling backbone** — the same higher-order geometry present in `rest_pre`
*and* `task` *and* `rest_post` (anatomy / individual fingerprint) — which, because
cophenetic distance is a **nonlinear** functional of the adjacency, induces a
**positive** `ρ_split` floor on its own, with zero task-specific persistence. The
**Shared-Backbone Coupling-Randomized Coherency (SB-CRC)** null fixes the
cross-phase common coherency `B(f)` verbatim across all four surrogate phases and
rotates only each phase's **deviation** `Δ_φ(f) = C_φ(f) − B(f)` by an independent
Haar real-orthogonal `O_φ`. This **lifts the floor off zero to exactly the
shared-backbone contribution** while making the task and rest deviations
conditionally independent given `B`. The test then decomposes the observed
`ρ_split` into a **shared-backbone part** (floor − 0) and a **task-specific part**
(observed − floor). If the observed clears the SB-CRC floor, the trace is
task-specifically persistent beyond the stable backbone; **if it does not, the
§5.3 cross-phase headline is substantially a stable-fingerprint effect and must be
reframed.** This null can bite the headline — that is the point.

> **Status: built + piloted; cohort running (2026-06-11).** Library helpers in
> `coherency_surrogate.py`; runner `audit_99_coordinated_null.py` (slot 98 was
> taken by the epi workstream). **Pilot (Pat_05): β floor lifts to +0.087, obs
> +0.463 → ~80% task-specific, clears the coordinated floor (p=0.000); α marginal
> (p=0.060) + grid-sensitive.** Substrate is CRC: the congruence makes `B` and the
> per-phase deviations cleanly separable and recombinable.
>
> **Empirical correction — PSD projection is MANDATORY, not cosmetic (2026-06-11
> cohort).** The deviation recombination `B + O Δ Oᵀ` is NOT PSD — every frequency
> bin strays negative. The pilot (Pat_05) suggested this was cosmetic (PSD-projected
> readout ρ=0.99); **the cohort proved that wrong.** For several patients (Pat_06,
> Pat_03, Pat_14) the non-PSD recombination drives the unit-diagonal renorm divisor
> `√diag` to ≈0 and the surrogate magnitude **blows up 1e2–1e8×** (Pat_06 β magR
> 3.3e5). Worse, even for *well-conditioned* patients (Pat_05, magR≈1) the raw
> off-manifold floor is **systematically too low**: projecting to the nearest PSD
> matrix raised Pat_05 β's floor +0.086 → +0.188 (Pat_06 β +0.029 → +0.312, which
> then *exceeds* its observed +0.263 → backbone-dominated, not task-specific). So the
> off-manifold construction was **anti-conservative everywhere** — it overstated the
> task-specific fraction. The fix (`psd_project=True`, now the library default, batched
> `eigh` clip-to-PSD before renorm) makes every surrogate a genuine realizable
> coherency, restores magnitude faithfulness (ratio 0.74–1.00), and gives the honest
> (higher) shared-backbone floor. **The first cohort run (raw, off-manifold) is void;
> the corrected run is canonical.** Lesson: a single-patient pilot cannot validate a
> null's calibration — cohort-wide magnitude/PSD checks are mandatory.

---

## 1. Critical preamble (5-point, mandatory per `feedback_critical_null_preamble.md`)

1. **Claim.** The §5.3 cross-phase cophenetic trace (`ρ_split^coph`; β
   load-bearing, α cophenetic-only) reflects **task-specific** reorganization that
   persists into `rest_post` — i.e. the alignment between the task-induced
   per-pair cophenetic shift and the `rest_post` per-pair shift **exceeds what a
   task-independent coupling backbone, shared across all phases, induces on its
   own**. The positive `ρ_split` is not merely the shadow of a stable individual
   connectivity fingerprint passing through the nonlinear LRG map.

2. **Null tested.** For each `(patient, band)`, build the four split-baseline
   phases' complex coherency stacks `C_φ(f)` on a **common frequency grid**
   (§5.1). Form the shared backbone `B(f) = Σ_φ w_φ C_φ(f)` (equal-condition
   weighting, §3.1) and the deviations `Δ_φ(f) = C_φ(f) − B(f)`. For each of
   `R = 200` draws and each phase φ, draw an **independent** Haar real-orthogonal
   `O_φ` (shared across the band) and form the surrogate coherency
   `C̃_φ(f) = renorm_unitdiag(B(f) + O_φ Δ_φ(f) O_φᵀ)`, then
   `W̃_φ = mean_f |Im C̃_φ(f)|`. Push each `W̃_φ` through the **identical**
   `cophenetic_condensed_from_adjacency` path and compute
   `ρ_split^surr = Spearman(D̃_task − D̃_preA, D̃_post − D̃_preB)`. Upper-tail
   `p = mean(ρ_split^surr ≥ ρ_split^obs)`; cohort paired one-sided Wilcoxon,
   `alternative='greater'`. The construction interpolates exactly between the
   independent-CRC null (put all content in `Δ`, `B=0` ⇒ floor ≈ 0) and the
   observed (put all content in `B`, `Δ=0` ⇒ floor = observed).

3. **Strongest alternative the null must control for.** A **task-independent
   stable coupling backbone.** The same coarse-scale communication geometry is
   present in every phase (it is the patient's anatomy / electrode placement /
   trait connectivity). Two mechanisms turn it into a positive `ρ_split` with **no
   task persistence**: (a) the cophenetic map is nonlinear, so a shared adjacency
   backbone does **not** cancel in the differences `D^task − D^preA`,
   `D^post − D^preB` (it would cancel only under an exact additive decomposition,
   which the Laplacian→propagator→ultrametric→cophenetic chain does not obey);
   (b) the backbone fixes **which pairs are labile** (near community boundaries)
   and which are rigid — a shared *heteroscedasticity* pattern — so the per-pair
   shifts in the task comparison and the rest comparison are rank-correlated purely
   through shared lability, independent of any task effect. The independent
   matched-strength and CRC nulls **destroy this backbone independently per phase**,
   so their floor (≈ 0) does not test it; they answer the wrong question
   ("structure vs. four independent draws") rather than the load-bearing one
   ("task-specific vs. shared-backbone").

4. **Does the null cover it — by mechanism, and what it CANNOT reject.**
   - **Covers it.** `B(f)` is held **identical** across all four surrogate phases,
     so the full cross-phase correlation it carries — including the nonlinear
     lability geometry (mechanisms a and b above) — is reproduced verbatim in the
     surrogate. Meanwhile the task and rest deviations are made conditionally
     independent given `B` (independent `O_task`, `O_post`). The resulting floor is
     therefore **exactly the `ρ_split` a shared backbone produces with no
     task-specific alignment**. Observed exceeding it isolates the task-specific
     part. The construction is on the CRC substrate, so each `W̃_φ` is a realizable,
     magnitude-faithful `|ImCoh|` matrix (the real-orthogonal congruence preserves
     `‖Im Δ_φ‖_F`; `B` contributes `|Im B|` verbatim) — the off-manifold and
     coupling-magnitude objections CRC already closes stay closed. ✓
   - **CANNOT reject — shared structure living in deviation second moments.**
     `B = mean` captures only the **first-moment** shared component. If the
     patient's stable structure also lives in a **shared covariance of the
     deviations** (a common lability sub-space across phases not captured by the
     mean), `B=mean` under-attributes it to the backbone ⇒ the floor is **too low**
     ⇒ the null is **anti-conservative** (overstates task-specificity). Mitigation
     and honest bracket: an **augmented backbone** `B⁺` carrying the leading shared
     low-rank deviation component, which raises the floor (§3.1, §11.1). Report
     both; a trace robust to both `B` levels is solid, a trace that clears `B=mean`
     but not `B⁺` is flagged.
   - **CANNOT reject — strength.** Like CRC, `O_φ` moves row-sums; SB-CRC does not
     control the hub/strength alternative. It inherits CRC's **bracket** status
     beside matched-strength. (Strength is already cleared independently;
     SB-CRC asks the orthogonal cross-phase question.)
   - **CANNOT reject — the third common cause (no-task control session).** SB-CRC
     is still a same-four-phase null. It separates task-specific from
     stable-backbone *within these phases*; it cannot exclude that a no-task
     control session would show the same `task→post` alignment for a reason
     unrelated to this task. Only an actual control session closes that.
   - **CANNOT reject — spatial embedding within the deviation.** `O_φ` mixes all
     channels, scrambling same-probe block structure **in the deviation**; the
     backbone's spatial structure is preserved (it is in `B`), but deviation-level
     probe geometry is not. Same anti-conservative-against-geometry caveat CRC has.

5. **Falsification + limitations.** *Falsification:* the observed `ρ_split` falls
   **inside** the SB-CRC distribution (upper-tail `p ≥ 0.05`) at β/α → the
   cross-phase trace is substantially explained by a task-independent stable
   backbone, and the §5.3 "task leaves a persistent multiscale trace" headline must
   be downgraded to "task-correlated structure consistent with a stable coupling
   backbone." *Limitations:* (a) `B=mean` weighting choice (equal-condition vs
   plain 4-mean, §3.1) and first-moment-only shared component (§4, mitigated by
   `B⁺`); (b) `B + O_φ Δ_φ O_φᵀ` is only **approximately** PSD (the `−O_φ B O_φᵀ`
   cross term can dip the spectrum slightly negative; reported as a min-eigenvalue
   residual, optional Higham projection §11.2 — note `W̃_φ ≥ 0` holds regardless,
   so LRG validity is never at risk); (c) the common-grid requirement forces the
   split-baseline halves onto the **full** `nperseg` (§5.1), shifting the
   within-pipeline observed slightly off the canonical `audit_63` value (both
   reported); (d) strength floats (bracket, not replacement).

---

## 2. Notation

| Symbol | Domain | Meaning |
|:---|:---|:---|
| `p ∈ P` | `|P| = 10` | patient (`COHORT_N10`) |
| `b ∈ B` | `{α, β, γ_l}` (primary) | band |
| `Φ` | `{rest_pre_A, rest_pre_B, task_test, rest_post}` | the four split-baseline phases |
| `N ≡ N_p` | int ∈ [113, 122] | contacts for patient `p` |
| `f ∈ F_b` | in-band bins on a **common grid** (§5.1) | `mask = (freqs ≥ f_min) & (freqs ≤ f_max)` |
| `C_φ(f)` | `ℂ^{N×N}` Hermitian PSD, unit diag | complex coherency of phase φ at bin `f` |
| `w_φ` | `Σ_φ w_φ = 1` | backbone weights (equal-condition default, §3.1) |
| `B(f)` | `ℂ^{N×N}` Hermitian PSD, unit diag | shared backbone `Σ_φ w_φ C_φ(f)` |
| `Δ_φ(f)` | `ℂ^{N×N}` Hermitian, **zero diag** | deviation `C_φ(f) − B(f)` |
| `O_φ` | `O(N)` real | Haar orthogonal, **independent per phase**, shared over `F_b`, fresh per surrogate |
| `C̃_φ(f)` | `ℂ^{N×N}` Hermitian, unit diag | surrogate coherency `renorm(B + O_φ Δ_φ O_φᵀ)` |
| `W̃_φ` | `ℝ_{≥0}^{N×N}`, `[0,1]`, zero diag | surrogate adjacency `mean_f |Im C̃_φ(f)|` |
| `D̃_φ` | `ℝ^{N(N−1)/2}` | surrogate condensed cophenetic vector at `τ = 1/λ_max` |
| `ρ_split` | `[−1, 1]` | `Spearman(D_task − D_preA, D_post − D_preB)` |

Everything downstream of `W̃_φ` (`L = D_W − W̃_φ → eigh → cophenetic_condensed_from_eigs`,
the `ρ_split` Spearman) is **identical** to `audit_63`/`audit_97`. SB-CRC changes
only the surrogate-generation step; it is plug-compatible with the existing
observed pipeline and reuses `cophenetic_condensed_from_adjacency`.

---

## 3. Definitions

### 3.1 Shared backbone `B(f)` and deviations `Δ_φ(f)`

`B(f)` is a convex combination of the four phase coherencies, hence Hermitian PSD
(convexity preserves PSD) with unit diagonal — itself a **valid coherency** (the
"average coherency"). Default **equal-condition** weighting (rest counts once, not
twice, since `pre_A`/`pre_B` are two noisy estimates of the *same* rest structure):

```
C_rest(f) = ½ (C_preA(f) + C_preB(f))
B(f)      = ⅓ ( C_rest(f) + C_task(f) + C_post(f) )           # w_preA = w_preB = 1/6, w_task = w_post = 1/3
Δ_φ(f)    = C_φ(f) − B(f)                                       # Hermitian, zero diagonal
```

`Δ_φ` has zero diagonal because `diag(C_φ) = diag(B) = 1`. Robustness knob: plain
four-phase mean `B = ¼ Σ_φ C_φ` (weights rest 1/2); reported as a sensitivity check
since it shifts the backbone toward rest.

**Augmented backbone `B⁺` (conservative companion, §11.1).** To guard against
shared structure hiding in deviation second moments, optionally fold the leading
**cross-phase-shared** deviation direction into the backbone: stack
`{Δ_φ(f)}_φ`, extract the rank-`k` component common to all phases (e.g. the top
singular subspace of the per-pair deviation matrix across phases), and add it to
`B`. `B⁺` raises the floor; the SB-CRC verdict is reported for both `B=mean`
(primary) and `B⁺` (conservative). v1 ships `B=mean`; `B⁺` is the first follow-up.

### 3.2 Haar rotation (Mezzadri 2007, real case)

`O_φ ∈ O(N)` real orthogonal, reusing `haar_orthogonal` from
`coherency_surrogate.py` (validated for CRC). A real congruence on the **complex**
Hermitian `Δ_φ` already yields a non-zero imaginary part
`Im(O Δ_φ Oᵀ) = O (Im Δ_φ) Oᵀ` and **preserves `‖Im Δ_φ‖_F`** — this is what keeps
the deviation's `|ImCoh|` magnitude faithful, exactly as in CRC.

### 3.3 Surrogate coherency (the construction)

For each `f ∈ F_b`, with one real orthogonal `O_φ` per phase (shared across the
band, fresh per surrogate):

```
C̃_raw,φ(f) = B(f) + O_φ Δ_φ(f) O_φᵀ                      # Hermitian; ‖Im Δ_φ‖_F preserved; B contributes |Im B| verbatim
d_φ(f)     = Re diag(C̃_raw,φ(f)) = 1 + Re diag(O_φ Δ_φ O_φᵀ)   # > 0 generically (guard)
C̃_φ(f)     = diag(d_φ)^{-1/2} · C̃_raw,φ(f) · diag(d_φ)^{-1/2}   # unit-diagonal congruence → valid coherency
```

`C̃_raw,φ` is **approximately** PSD: `B ⪰ 0` dominates and `Δ_φ` is a small
deviation, but the `−O_φ B O_φᵀ` cross term inside `O_φ Δ_φ O_φᵀ = O_φ C_φ O_φᵀ −
O_φ B O_φᵀ` can push the smallest eigenvalue slightly below 0. Reported as a
**PSD residual** `min_f λ_min(C̃_raw,φ(f))` in calibration; optional Higham PSD
projection (§11.2). Crucially, the downstream readout `W̃_φ = mean_f |Im C̃_φ|` is
**non-negative by the entrywise `|·|` regardless of PSD**, so the graph Laplacian
`L = D_W − W̃_φ` is PSD and the LRG propagator is always well-defined — the PSD
residual is an *on-manifold-fidelity* diagnostic, not a validity gate.

### 3.4 Surrogate adjacency and statistic (reused verbatim)

```
W̃_φ   = mean_{f∈F_b} |Im C̃_φ(f)|                  ∈ [0,1]^{N×N}, symmetric, zero-diag, ≥ 0
D̃_φ   = cophenetic_condensed_from_adjacency(W̃_φ)   # L = diag(W̃·1) − W̃ ; eigh ; ρ̂(τ_max) ; T=1/ρ̂ ; avg-linkage ; cophenet
ρ_split^surr = Spearman(D̃_task − D̃_preA, D̃_post − D̃_preB)
```

---

## 4. Properties

- **Floor lifts off zero** (the defining property). Because `B` is shared across all
  four surrogate phases and the cophenetic map is nonlinear, `ρ_split^surr` has a
  **positive** expectation = the shared-backbone contribution. Contrast: the
  independent-CRC / matched-strength nulls rotate `B+Δ` together per phase, so
  their backbone decorrelates across phases and their floor is ≈ 0. SB-CRC strictly
  generalizes both endpoints (§1.2).
- **ρ_split decomposition** (the deliverable). `ρ_split^obs ≈ (independent-CRC
  floor ≈ 0) + (SB-CRC floor − 0)  [shared-backbone] + (obs − SB-CRC floor)
  [task-specific]`. The test reads the **third** term; the **second** term is the
  newly-measured shared-backbone size.
- **On-manifold + magnitude-faithful** (inherited from CRC). Each `W̃_φ` is a
  realizable `|ImCoh|`; real-orthogonal congruence preserves `‖Im Δ_φ‖_F`, `B`
  contributes its own imaginary energy verbatim ⇒ calibration ratio ≈ 1 expected.
- **Non-degeneracy of the shared-`B` design.** Applying the *same* `O` to all
  phases' **deviations** (with `B` fixed) is **not** degenerate (unlike rotating
  `B+Δ` together, which is a consistent global relabeling that returns `T_obs`);
  but same-`O` would keep task/rest deviations aligned ⇒ floor too high ⇒
  over-conservative. **Independent `O_φ` is required** to make the deviations
  conditionally independent given `B`. (Resolves the §12 degeneracy note of the CRC
  scope.)
- **Strength NOT preserved; spatial deviation structure NOT preserved** (negative,
  inherited — bracket with matched-strength; backbone spatial structure IS
  preserved).
- **Anti-conservative if shared structure hides in deviation 2nd moments** (the one
  real threat to validity; mitigated by `B⁺`, §3.1/§11.1).
- **Complexity.** Same order as CRC × (4 phases simultaneously in memory):
  `O(R · 4 · |F_b| · N³)`. `build_coherency` in `audit_97` already holds all four
  phases, so memory is unchanged; runtime ≈ CRC cohort (~minutes/cell).

---

## 5. Construction details that differ from CRC

### 5.1 Common frequency grid (the one genuinely new requirement)

`B(f) = Σ_φ w_φ C_φ(f)` is a **per-frequency** average, so all four phases must
share the in-band frequency grid. The canonical split-baseline halves
(`compute_imcoh_abs_halves`) use `nperseg//2` (to keep the segment count, hence
CSD variance, comparable to the full phases), which puts the halves on a **coarser
grid** than `task`/`post`. SB-CRC therefore computes **all four phases at the full
`nperseg`** (`nperseg_for_fs(fs)`), splitting `rest_pre` at `T//2`. Consequence:
the halves carry ~half the Welch segments (noisier CSD) than the canonical
`nperseg//2` halves, so the **within-pipeline observed** `ρ_split` differs slightly
from the canonical `audit_63` value. Both are reported per cell; the floor and the
within-pipeline observed share this grid, so the test is internally exact, and the
canonical value is the cross-reference (sanity, like `audit_97`). If the shift is
material it is itself a grid-sensitivity finding → escalate to per-frequency
interpolation onto the canonical half-grid (§11.3).

### 5.2 Independent vs coordinated — what changes in code

`audit_97`'s `crc_null` calls `coupling_randomized_coherency(C_φ)` **independently
per phase** (rotates `B+Δ` together ⇒ floor ≈ 0). SB-CRC instead: compute `B`,
`Δ_φ` **once per cell from the four phases**, then per draw rotate **only `Δ_φ`**
per phase and recombine with the shared `B`. One new library function; the per-cell
orchestration, the cophenetic path, and `ρ_split` are reused.

---

## 6. Calibration (run before any verdict)

Per `(patient, band)`:
1. **Floor lift** — SB-CRC floor median vs the independent-CRC floor (≈ 0) on the
   **same common grid**. This is the measured shared-backbone size; the headline
   number. (Must be ≥ 0 and < observed for the construction to be sitting between
   its two endpoints — a sanity gate.)
2. **Backbone energy fraction** — `‖B‖_F² / mean_φ ‖C_φ‖_F²` and per-phase
   `‖Δ_φ‖_F / ‖C_φ‖_F` (how much is backbone vs deviation; sets how much room the
   randomization has).
3. **Magnitude match** — surrogate `mean(W̃[triu])` vs observed (expect ≈ 1, since
   `‖Im Δ_φ‖` preserved + `|Im B|` verbatim).
4. **PSD residual** — `min_f λ_min(C̃_raw,φ)` (how far off-manifold the recombination
   strays; expect small-negative; report worst over phases).
5. **Within-pipeline vs canonical observed** — `ρ_split` on the full-`nperseg` grid
   vs the cached `audit_63` value (grid-shift sanity, §5.1).

A failed calibration (floor lift < 0, magnitude ratio far from 1, or large
within-pipeline/canonical divergence) voids that cell's SB-CRC `p`.

---

## 7. Pseudocode

```
INPUT: patient p, band b, Φ = {preA, preB, task, post}, R = 200, seed, rng
# common grid: all four phases at full nperseg
fs ← sample_rate(p); nseg ← nperseg_for_fs(fs)
X_rest ← load_timeseries(p, "rest_pre"); T ← len
C_preA ← coherency_band(X_rest[:, :T//2], fs, b, nseg)        # (Fb, N, N) complex
C_preB ← coherency_band(X_rest[:, T//2:], fs, b, nseg)
C_task ← coherency_band(load_timeseries(p, "task"), fs, b, nseg)
C_post ← coherency_band(load_timeseries(p, "post"), fs, b, nseg)

# shared backbone + deviations (once per cell)
C_rest ← ½(C_preA + C_preB)
B      ← ⅓(C_rest + C_task + C_post)                          # equal-condition weights
Δ      ← { φ: C_φ − B  for φ in Φ }                            # zero-diagonal

# within-pipeline observed (full-grid) + canonical cross-ref
obs    ← rho_split({ φ: cophenetic(mean_f |Im C_φ|) for φ in Φ })
obs_canonical ← audit_63 cached value                          # sanity

FOR r in 1..R:
    D̃ ← {}
    FOR φ in Φ:
        O_φ ← haar_orthogonal(N, rng)                          # INDEPENDENT per phase
        Craw ← B + O_φ @ Δ[φ] @ O_φ.T                          # batched einsum over f
        d    ← max(Re diag(Craw), eps)
        Ctil ← Craw / sqrt(outer(d, d))                        # unit-diagonal congruence
        W̃    ← mean_f |Im Ctil|
        D̃[φ] ← cophenetic_condensed_from_adjacency(W̃)
        record PSD residual min_f λ_min(Craw)
    ρ_surr[r] ← spearman(D̃[task] − D̃[preA], D̃[post] − D̃[preB])

p ← mean(ρ_surr ≥ obs)                                          # upper tail
# companion (same machinery, B → 0): independent-CRC floor on the SAME grid
RETURN obs, obs_canonical, p, floor_median = median(ρ_surr),
       independent_crc_floor, calibration
```

---

## 8. Visualization spec

- **8.1 Decomposition panel** (headline). Per band, n=10 small multiples: a
  horizontal axis in `ρ_split`. For each patient draw **two** stacked null bars —
  the **independent-CRC floor** (≈ 0, reference, grey) and the **SB-CRC floor**
  (lifted, P5–P95, colour) — plus the observed marker. The visual story is the
  *gap structure*: `0 → SB-CRC floor` = shared-backbone part; `SB-CRC floor →
  observed` = task-specific part. A trace is **task-specific** when the observed
  sits in the right tail of the SB-CRC bar; it is **backbone-explained** when the
  observed sits **inside** the SB-CRC bar even though it cleared the independent
  floor.
- **8.2 Calibration panel** (§6): floor-lift (SB-CRC − independent), magnitude
  ratio, PSD residual, within-pipeline-vs-canonical observed — per band, dots, no
  in-axes text (print to stdout / companion `.md`).
- Style: `use_lrg_style()`, PDF-only full-vector, transparent, no suptitle, math
  `$i$,$j$` axis labels, ≥ 3 patients shown.

---

## 9. Connection to prior tools

| Tool | Preserves cross-phase | Randomizes | `ρ_split` floor | Separates |
|:---|:---|:---|:---|:---|
| Matched-strength (`audit_62/63/66`) | nothing (independent per phase) | adjacency, off-manifold | ≈ 0 | structure vs 4 indep. draws |
| CRC (`audit_96/97`) | nothing (independent per phase) | coupling identity, on-manifold | ≈ 0 | structure vs 4 indep. draws |
| **SB-CRC (this scope)** | **shared backbone `B` verbatim** | **deviation `Δ_φ` only, independent `O_φ`** | **> 0 (= backbone part)** | **task-specific vs stable backbone** |
| Split-baseline halves (`audit_63`) | independent rsPre halves | — | n/a | shared-baseline noise |

SB-CRC is the **coordinated cross-phase member** of the CRC family: same substrate,
same cophenetic pipeline, same `ρ_split` statistic; it differs only in **what is
held fixed across phases** (the backbone) vs **randomized** (the deviation). It
answers the question the independent nulls structurally cannot. It does **not**
replace them — matched-strength still owns the strength axis, CRC the
realizability axis; SB-CRC owns the **task-specific-vs-backbone** axis. The three
together bracket the trace on its three independent objections.

---

## 10. Implementation plan

- **Library** (`src/lrg_eegfc/utils/surrogate/coherency_surrogate.py`, extend):
  - `shared_backbone_deviation(C_by_phase, weights="equal_condition")
    -> (B_band, {phase: Δ_band})` — §3.1.
  - `deviation_rotated_coherency(B_band, Delta_band, rng, eps=1e-12)
    -> (W_surr, renorm_resid, psd_resid)` — one phase's surrogate, §3.3, batched
    einsum mirroring `coupling_randomized_coherency`.
  - export both via `surrogate/__init__.py`.
- **Runner (single dual-mode script)**:
  `scripts/01_compute/audit/audit_99_coordinated_null.py` (slot 98 was taken by
  `audit_98_epi_seed_prediction_pipeline.py`; 89–95 epi, 96/97 CRC). Defaults to the
  full n=10 cohort × {β, α, γ_l}, R=200, common full-`nperseg` grid; `--patients
  Pat_05 --bands beta alpha` reproduces the pilot. Emits within-pipeline obs,
  canonical obs, SB-CRC floor, independent-CRC floor, the decomposition, calibration,
  the per-band strip figure, `cohort_summary.csv` (paired Wilcoxon obs vs SB-CRC
  floor). Output `data/audit/coordinated_cross_phase_null/`.
- **Checkpoint with user** done after the pilot (β survives, ~80% task-specific) →
  cohort approved and running.
- **Cohort** verdict is DESCRIPTIVE (no pre-registered gate); the headline-relevant
  case is the *weaker* patients — backbone-dominated (floor ≈ obs) or task-specific?
- **No reuse violations**: `load_timeseries`/`welch_csd`/`cophenetic_*` from the
  library; stats from `utils.metrics.hypothesis`; paths from `config.paths`.

---

## 11. Open questions

1. **Backbone richness `B` vs `B⁺`.** Is `B = mean` enough, or does shared
   deviation-covariance demand the augmented `B⁺`? Ship `B=mean`; if a cell clears
   `B=mean` narrowly, run `B⁺` as the conservative bracket and report both.
2. **PSD projection.** Default OFF (report residual). Turn ON (Higham nearest-PSD on
   `C̃_raw,φ`) only if residuals are large enough to question the on-manifold claim;
   `W̃ ≥ 0` is never at risk either way.
3. **Common grid vs interpolation.** Full-`nperseg` halves (default, §5.1) vs
   interpolating `task`/`post` onto the canonical `nperseg//2` half-grid. Escalate
   to interpolation only if the within-pipeline/canonical observed shift is material.
4. **Backbone weighting.** Equal-condition (default) vs plain 4-mean vs
   task-excluded backbone `B = ½(C_rest + C_post)` (a "rest-only backbone" that asks
   whether task itself is even needed to explain the alignment). The last is a
   sharper, more adversarial backbone — worth a sensitivity row.
5. **Statistic.** Primary `ρ_split^coph` (manuscript headline); secondary triangle
   `T_d` (locked sign, positive = trace) and cross-phase Grassmann `T_G`. Lock to
   the §5 headline for apples-to-apples.
6. **`R` / seed.** `R = 200`, seed `20260611`, paralleling CRC; revisit if floor
   variance demands more draws.
