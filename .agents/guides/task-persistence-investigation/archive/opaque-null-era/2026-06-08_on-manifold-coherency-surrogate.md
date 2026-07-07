---
name: on-manifold-coherency-surrogate
type: scope
era: COHORT_N10
status: draft
created: 2026-06-08
updated: 2026-06-08
pointers:
  - src/lrg_eegfc/utils/surrogate/matched_strength.py
  - src/lrg_eegfc/utils/fc/coherence/imcoh.py
  - src/lrg_eegfc/utils/fc/coherence/_common.py
  - src/lrg_eegfc/utils/metrics/spectral.py
  - scripts/01_compute/audit/audit_63_split_baseline_surrogate.py
  - .agents/guides/task-persistence-investigation/2026-05-30_asymmetric-pair-trace-regression.md
  - .agents/guides/02_methods/imcoh-guide.md
---

# Coupling-Randomized Coherency (CRC) surrogate — an on-manifold, spectrum-matched null for `|ImCoh|` LRG statistics

## Renormalization head

The mandatory matched-strength null (`audit_62/63/65/66`) preserves each node's
strength but builds its surrogates by 4-cycle ±δ rewiring of the *already-computed*
`|ImCoh|` adjacency — which samples adjacency matrices that **no Hermitian PSD
cross-spectral density can generate** (off the `|ImCoh|` manifold) and distorts the
weight distribution. The **Coupling-Randomized Coherency (CRC)** surrogate is its
**complement**: for each in-band frequency it conjugates the *complex* coherency
`C(f)` by a Haar-random **real orthogonal** matrix `O` (`C̃ = O C(f) Oᵀ`,
randomizing *which* node-sets carry the coupling) — which preserves both the
coherency eigenvalue spectrum and the imaginary-coupling energy `‖Im C(f)‖_F`
(the *amount* of lagged coupling, hence the `|ImCoh|` magnitude) exactly — then
re-derives `|ImCoh|`. Every surrogate is therefore a genuine, realizable
`|ImCoh|` matrix with non-negative edges (valid Laplacian, well-defined
propagator) by construction. CRC controls the realizability + coupling-magnitude
alternatives that matched-strength cannot; it does **not** preserve node
strength, so it is a **bracket** beside matched-strength, never a replacement. A
cross-phase LRG trace that clears *both* nulls is robust to the strongest version
of each objection.

> **Status: current — cohort verdict in (2026-06-10).** Design pinned per the
> folder rule, prototyped (`audit_96`), cohort-run (`audit_97`). **Headline: CRC
> reproduces the matched-strength verdict at 30/30 cells** — the off-manifold
> concern about matched-strength is not load-bearing for the §5.3 cophenetic trace.

## Empirical status & construction correction (2026-06-08)

The construction below was originally drafted with a **complex** unitary on the
eigenvalue matrix (`C̃ = Q diag(Λ) Qᴴ`, §3.2–3.3 as written). The pilot
(`audit_96`, Pat_05 β) **rejected** that variant: a complex-Haar eigenbasis
*generates* imaginary mass out of real coherence, inflating the surrogate
`|ImCoh|` magnitude **4.2×** over observed (the data has `Im ≪ Re`; a "fair"
complex split hands the surrogate ~half the *total* coherency content as
imaginary). The **validated** construction is the **real-orthogonal congruence**
`C̃ = O C(f) Oᵀ` (`O ∈ O(N)` Haar): it preserves the eigenvalue spectrum (to
4e-15) **and** `‖Im C(f)‖_F` (to ratio 1.0000) exactly, so the surrogate
`|ImCoh|` magnitude matches observed (calibration ratio **1.12**, the 0.12 from
the unit-diagonal renorm). Read §3.2–3.3 with `Q`-complex replaced by `O`-real;
the complex variant survives only as a `complex_eigbasis` diagnostic in
`coupling_randomized_coherency`.

Two further empirical notes:
- **Renorm residual is order-1** (Pat_05 β ≈ 1.7). The unit-diagonal congruence
  distorts the spectrum more than "small" — magnitude survives it (ratio 1.12)
  but the *spectrum* is held only approximately. Property §4 "Λ preserved"
  is softened to "Λ and `‖Im C‖` preserved by the rotation, approximately after
  renorm." Open: a Higham-iterated unit-diagonal projection (§11) to tighten it.
- **Shared- vs split-baseline.** The `audit_96` pilot used a *shared*-baseline
  `ρ_split` (common `coph_pre`), which inflated **both** the observed *and* the
  CRC null (Pat_05 β: obs +0.66, CRC null mean +0.58, p=0.155 — a misleading
  "may-not-clear" flag). The §5.3 **split**-baseline protocol (`audit_63`,
  independent `pre_A`/`pre_B`) removes the shared term: the CRC null collapses to
  ≈0 (Pat_05 β CRC null median −0.013) and the observed (+0.49) clears both CRC
  and matched-strength.

### Cohort verdict (`audit_97`, n=10 × {α, β, γ_l}, R=200)

| band | obs med ρ | MS null med | CRC null med | MS Wilcoxon p | CRC Wilcoxon p | n>MS | n>CRC | CRC mag-ratio | verdict |
|---|---|---|---|---|---|---|---|---|---|
| α | +0.105 | +0.014 | +0.002 | 0.0049 | 0.0049 | 5/10 | 5/10 | 1.07 | clears_both |
| β | +0.221 | +0.012 | +0.001 | 0.0068 | 0.0068 | 7/10 | 7/10 | 1.07 | clears_both |
| γ_l | +0.083 | +0.001 | +0.001 | 0.116 | 0.116 | 5/10 | 5/10 | 1.02 | clears_neither |

**CRC reproduces the matched-strength clear/no-clear decision at 30/30 cells**
(per-patient `crc_p` vs `ms_p` agree 10/10 in every band); magnitude ratio
0.97–1.19 throughout (calibration holds). The on-manifold concern does **not**
overturn the §5.3 cophenetic trace — matched-strength was adequate *on the
realizability axis*. α/β cophenetic trace clears both nulls; γ_l is cohort-null
under both (split: strong in Pat_02/05/06/08, negative in Pat_07/13/14/15).

**Honest limitation (unchanged by CRC).** Both nulls use *independent per-phase*
randomization, so both collapse to ≈0 under split-baseline and both test the same
"is there cross-phase ρ beyond independent per-phase randomization" question. CRC
closes the *off-manifold* gap specifically; it does **not** test a *coordinated
cross-phase* (strength-trajectory-preserving) null — that gap, flagged in
`audit_63`'s own note, remains open for both. Outputs: `data/audit/crc_cohort_bracket/`.

---

## 1. Critical preamble (5-point, mandatory per `feedback_critical_null_preamble.md`)

1. **Claim.** The cross-phase LRG cophenetic trace (the §5.3 `ρ_split^coph` /
   triangle finding; β load-bearing, α cophenetic-only) reflects genuine
   *higher-order coupling geometry* — which node-sets block together at coarse
   communication scale — and survives a null that, unlike matched-strength,
   (i) keeps every surrogate a **realizable** `|ImCoh|` matrix and (ii) holds the
   **per-frequency coupling-magnitude spectrum** (the coherency eigenvalues) fixed.
   I.e. the trace is not an artifact of comparing real on-manifold matrices against
   off-manifold rewired ones, nor of the sheer amount/concentration of coupling per
   frequency.
2. **Null tested.** For each `(patient, band, phase)`, generate `R = 200` surrogate
   `|ImCoh|` matrices: per in-band frequency `f`, eigendecompose the complex
   Hermitian coherency `C(f) = U Λ Uᴴ` and replace `U` by a single Haar-random
   complex unitary `Q` (shared across the band's frequencies, fresh per surrogate,
   **independent per phase** — mirroring the per-phase independence of
   `load_or_compute_surrogate_eigs`): `C̃(f) = Q Λ(f) Qᴴ`, renormalized to unit
   diagonal by positive-diagonal congruence; then `W_surr = mean_f |Im C̃(f)|`. Push
   `W_surr` through the **identical** `_laplacian_eig → cophenetic_condensed_from_eigs`
   path used for the observed statistic. Cross-phase statistic computed with the
   same split-baseline halves convention as `audit_63`; upper-tail
   `p = mean(s_finite ≥ obs)`.
3. **Strongest alternatives the null must control for.**
   (a) **Realizability artifact** — the matched-strength 4-cycle shuffle samples
   adjacencies outside the image of the `|ImCoh|` map, so the observed value could
   read extreme partly *because surrogates are off-manifold*, not because of genuine
   structure.
   (b) **Coupling-magnitude artifact** — the trace could be driven purely by *how
   much* coupling exists at each frequency and *how concentrated* it is (the
   coherency eigenvalue spectrum / distance of `C` from the identity), independent of
   *which* nodes carry it.
4. **Does the null cover them — by mechanism, and what it CANNOT reject.**
   - Covers (a): each `C̃(f)` is a unitary conjugation of a Hermitian PSD matrix,
     then a positive-diagonal congruence to unit diagonal — hence Hermitian PSD with
     unit diagonal, i.e. **a genuine coherency**. `W_surr = mean_f |Im C̃|` is in the
     image of the `|ImCoh|` map by construction. Observed and surrogate are compared
     **on the same manifold**. ✓
   - Covers (b): the Haar conjugation `Q Λ Qᴴ` preserves `Λ(f)` *exactly*; the
     mandatory unit-diagonal renormalization perturbs it by a controlled amount set
     by the spread of `diag(Q Λ Qᴴ)` (reported as a calibration residual — must be
     `≪` the observed effect). So the per-frequency coupling magnitude/concentration
     is held (near-)fixed while eigenvector identity is fully randomized. ✓
   - **CANNOT reject — node strength.** Row-sums of `|ImCoh|` are a nonlinear
     functional of `C`; the Haar rotation moves them, so the surrogate strength
     sequence drifts off the observed. CRC does **not** control "the trace is a
     hub/strength effect" — that is exactly what matched-strength controls. **CRC and
     matched-strength are complementary brackets**: matched-strength fixes strength
     but breaks realizability + the coupling spectrum; CRC fixes realizability + the
     coupling spectrum but lets strength float. Clearing only one is not enough.
   - **CANNOT reject — spatial embedding.** `Q` mixes all channels uniformly, so any
     same-probe / proximity block structure present in *every* phase is scrambled —
     the same anti-conservative-against-geometry limitation matched-strength has
     (mitigation: probe-block-constrained rotation, §11). And CRC does **not**
     separate task from time-on-task drift (the within-baseline/drift null's job;
     CRC is run on the split-baseline halves exactly as `audit_63`).
5. **Falsification + limitations.** *Falsification:* the observed cross-phase
   cophenetic statistic falls **inside** the CRC distribution (upper-tail `p ≥ 0.05`)
   at the bands where matched-strength clears (β, α) → the on-manifold +
   coupling-spectrum explanation is sufficient and the higher-order-identity reading
   is not supported there. *Limitations:* (a) strength floats → read as a bracket,
   not a replacement; (b) the unit-diagonal renormalization perturbs `Λ` slightly
   (residual reported); (c) shared-`Q`-per-band imposes frequency-constant coupling
   directions (smoother in `f` than data) — conservative (wider null), but a modeling
   choice; per-`f` independent `Q` shrinks band-averaged magnitudes via the central
   limit theorem and is rejected as default (§11); (d) requires recomputing the
   complex coherency from the timeseries (the `imcoh` cache stores only `Im C`);
   (e) CRC preserves *total* coherency coupling content and lets the in-phase/lagged
   split vary under complex-Haar mixing, so the surrogate `|ImCoh|` magnitude level is
   set by total coupling, not by the observed `|ImCoh|` level — acceptable because the
   gated statistic is a *rank* correlation of per-pair shifts (magnitude-robust), but
   calibrated in §6/§8; (f) the third-common-cause (a no-task control session) is
   still unavailable — CRC, like every same-four-phase null, cannot rule it out.

---

## 2. Notation

| Symbol | Domain | Meaning |
|:---|:---|:---|
| `p ∈ P` | `|P| = 10` (`config.const.PATIENTS_4PHASE`) | patient |
| `b ∈ B` | `|B| = 6` (δ, θ, α, β, γ_l, γ_h) | band |
| `φ` | split-baseline phase set `{rest_pre_A, rest_pre_B, task_test, rest_post}` | phase |
| `N ≡ N_p` | int, ∈ [113, 122] | contacts for patient `p` |
| `X_p^φ` | `ℝ^{N × L}` | timeseries (Pat_10 task rows `[53,54,55]` dropped at load) |
| `f ∈ F_b` | in-band frequency bins | `mask = (freqs ≥ f_min) & (freqs ≤ f_max)` |
| `S(f)` | `ℂ^{N×N}`, Hermitian PSD | Welch CSD, `welch_csd(X, fs, nperseg_for_fs(fs))` |
| `P_i(f) = S_ii(f)` | `ℝ_{>0}` | per-channel auto-power |
| `C(f)` | `ℂ^{N×N}`, Hermitian PSD, unit diag | coherency `D^{-1/2} S(f) D^{-1/2}`, `D = diag(S)` |
| `Λ(f), U(f)` | `Λ ∈ ℝ_{≥0}^N`, `Σ_k Λ_k = N`; `U ∈ U(N)` | eigenvalues / eigenvectors of `C(f)` |
| `Q` | `U(N)` (complex) | Haar-random unitary, one per `(surrogate, phase)`, shared over `F_b` |
| `A_p^{b,φ}` | `ℝ_{≥0}^{N×N}`, `[0,1]`, zero-diag | **observed** `|ImCoh|` adjacency `= mean_{f∈F_b} |Im C(f)|` |
| `W_surr` | `ℝ_{≥0}^{N×N}`, `[0,1]`, zero-diag | **surrogate** adjacency `= mean_{f∈F_b} |Im C̃(f)|` |
| `X^{(p,b,φ)}` | `ℝ^{N×N}` | downstream LRG primitive (`ρ̂`, `D=1/ρ̂`, or `D_coph`) at `τ = 1/λ_max` |

The downstream reconstruction `(λ, V) → ρ̂(τ_max) → D → D_coph` and the cross-phase
statistic are **identical** to `2026-05-30_asymmetric-pair-trace-regression.md` §2 and
`cophenetic_condensed_from_eigs` in `matched_strength.py`. CRC changes only **how the
surrogate adjacency `W_surr` is generated**; every line downstream of `W_surr` is reused
verbatim, so the null is plug-compatible with the existing observed pipeline.

---

## 3. Definitions

### 3.1 Complex coherency (the on-manifold object)
```
C(f) := D(f)^{-1/2} S(f) D(f)^{-1/2},   D(f) = diag(S(f)) ∈ ℝ_{>0}^{N}
```
`C(f)` is Hermitian PSD (congruence of the PSD CSD), unit diagonal, eigenvalues
`Λ(f) ≥ 0` with `Σ_k Λ_k(f) = tr C(f) = N`. The observed adjacency is
`A = mean_{f∈F_b} |Im C(f)|` — *bit-identical* to `compute_imcoh` followed by
`band_average` then `np.abs` (the load-time magnitude transform). `Re C` is the
real coherence, discarded by `|ImCoh|`; CRC needs the full complex `C(f)`, hence the
recompute from timeseries.

### 3.2 Haar-random rotation (Mezzadri 2007)
Validated construction uses a **real orthogonal** `O ∈ O(N)`:
```
Z ~ N(0,1)^{N×N};   O0, R = qr(Z);   O = O0 · diag(sign(R_kk))   ∈ O(N), Haar
```
Because the congruence acts on the *full complex* `C(f)` (§3.3), a real `O`
already produces a non-zero imaginary part: `Im(O C Oᵀ) = O (Im C) Oᵀ ≠ 0`, and
its Frobenius norm is preserved. (The complex unitary `Q ∈ U(N)` — same recipe
with `Z = (G1 + i G2)/√2` — is needed *only* by the rejected `complex_eigbasis`
route, where it rotates the real eigenvalue matrix `diag(Λ)`; there a real `O`
would give `Im ≡ 0`. That route is rejected for magnitude inflation, §3.3.)

### 3.3 Coupling-randomized coherency (the construction — VALIDATED form)
For each `f ∈ F_b`, with one real orthogonal `O ∈ O(N)` shared across the band
(fresh per surrogate, per phase):
```
C̃_raw(f) = O · C(f) · Oᵀ                     # Hermitian PSD; Λ(f) AND ‖Im C(f)‖_F preserved exactly
d(f)      = Re diag(C̃_raw(f))               # > 0 generically (guard against 0)
C̃(f)      = diag(d(f))^{-1/2} · C̃_raw(f) · diag(d(f))^{-1/2}   # unit-diagonal congruence → valid coherency
```
A real-orthogonal congruence is a unitary similarity, so it preserves the
eigenvalues `Λ(f)` exactly **and** preserves the Frobenius norm of the imaginary
part `‖Im C(f)‖_F` (the ImCoh energy) exactly — this is what keeps the surrogate
`|ImCoh|` magnitude faithful (calibration ratio ≈ 1.1). The unit-diagonal
congruence then forces a valid coherency in one step; empirically it perturbs the
spectrum order-1 (residual ≈ 1.7) while leaving the magnitude faithful. Range:
`C̃(f)` Hermitian PSD, unit diagonal ⇒ genuine coherency ⇒ `|Im C̃_ij| ≤ 1`.

> **Rejected diagnostic (`complex_eigbasis`).** The original `C̃ = Q diag(Λ) Qᴴ`
> with complex Haar `Q ∈ U(N)` preserves `Λ` but, because a random *complex*
> eigenbasis splits the preserved off-diagonal content evenly into real/imaginary,
> it inflates `|ImCoh|` 4.2× when the data has `Im ≪ Re`. Kept in code only to
> reproduce this failure; never used as a null.

### 3.4 Surrogate adjacency
```
W_surr = mean_{f∈F_b} |Im C̃(f)|     ∈ [0,1]^{N×N}, symmetric, zero-diagonal, ≥ 0
```
Non-negativity is by the entrywise `|·|`, *identical* to the observed transform —
this is what guarantees `L = D_W − W_surr` is PSD and the LRG propagator
`ρ = e^{-τL}` is a valid density operator (the constraint clarified by the user:
non-negative edges, not `W` PSD).

### 3.5 Downstream statistic + p-value (reused verbatim)
`W_surr → (λ, V) = eigh(diag(W_surr·1) − W_surr) → cophenetic_condensed_from_eigs(λ, V)`
gives the surrogate cophenetic vector at `τ = 1/λ_max`. The cross-phase trace
statistic `T_obs` (primary: `ρ_split^coph`; secondary: triangle `T_d`, locked
sign `T_d = d(rsPre,task) − d(task,rsPost)`, positive = trace) is computed from the
surrogate phase-triplet exactly as for observed data. Cohort: per-patient upper-tail
`p_p = mean(s_finite ≥ T_obs,p)` and paired one-sided Wilcoxon
`T_obs > median(surrogate)`, `alternative='greater'`.

---

## 4. Properties

- **Non-negativity / LRG-validity** (positive). `W_surr ≥ 0` entrywise ⇒ valid graph
  Laplacian, well-defined propagator. Holds for *every* surrogate, no rejection step.
- **Realizability / on-manifold** (positive). Each surrogate derives from a genuine
  Hermitian PSD coherency ⇒ in the image of the `|ImCoh|` map. This is precisely the
  property matched-strength lacks.
- **Coupling-spectrum preservation** (positive, up to renorm). `Λ(f)` exact under `Q`,
  perturbed only by the unit-diagonal congruence; `‖C − I‖_F² = Σ_k(Λ_k − 1)²` (total
  coherency coupling content) preserved up to the same residual.
- **Strength NOT preserved** (negative — the contract). Row-sums drift; CRC says
  nothing about the hub/strength alternative. Must be bracketed with matched-strength.
- **Spatial embedding NOT preserved** (negative). Uniform channel mixing scrambles
  probe-block structure; anti-conservative against spatially-structured nuisance.
- **In-phase/lagged split NOT preserved** (negative). Complex-Haar mixing redistributes
  the `Re/Im` split of the preserved total coupling, so surrogate `|ImCoh|` *magnitude*
  is set by total coupling, not by observed `|ImCoh|` magnitude (calibrated, §6).
- **Identifiability** (negative edge case). When `C(f) ≈ I` (no coupling) `Λ ≈ 1`,
  `Q Λ Qᴴ ≈ I`, `Im C̃ ≈ 0` — surrogate *and* observed `|ImCoh| ≈ 0`, test correctly
  uninformative (no coupling to randomize).
- **Complexity.** `O(R · |F_b| · N³)` — one `eigh` per in-band frequency per surrogate,
  **heavier** than matched-strength (one `eigh` per surrogate). Mitigated by caching
  `W_surr`/its Laplacian eigendecomposition (§10) and by `|F_b|` being modest at the
  4096-sample `nperseg`.

---

## 5. Caveats & failure modes

| Failure mode | Mechanism | Mitigation |
|:---|:---|:---|
| **Band-averaged magnitude collapse** | per-`f` independent `Q` ⇒ `Im C̃(f)` independent across `f` ⇒ `mean_f|·|` concentrates (CLT) ⇒ artificially homogeneous, too-narrow null | **shared `Q` per band** as default; calibrate surrogate `|ImCoh|` entry + strength distributions against observed (§8) |
| **`Λ` perturbation by renorm** | unit-diagonal congruence is non-unitary | report `max_f ‖Λ̃(f) − Λ(f)‖ / N` and the `diag(C̃_raw)` spread; require `≪` observed effect |
| **Magnitude mismatch** | `Re/Im` split varies (Property §4) | calibrate; if surrogate `|ImCoh|` strength is badly off observed, escalate to CRC-v2 with an imaginary-energy constraint (§11) |
| **Rank-deficient CSD** | few Welch segments ⇒ some `Λ_k = 0` | harmless (zero eigenvalues rotate fine); Pat_03 1024 Hz handled at config via `nperseg_for_fs` |
| **Degenerate `Im C̃ ≈ 0`** | `C ≈ I` (no coupling) | guard; flag the cell as uninformative rather than NaN |
| **Off-band leakage** | wrong frequency mask | reuse `band_average`'s exact `(freqs ≥ f_min) & (freqs ≤ f_max)` mask |
| **Cost blow-up** | per-`f` `eigh` over `R=200` × cohort | cache surrogate Laplacian eigs (`load_or_compute_eigs_at_path` pattern); compute once, reuse for `ρ_split`/triangle/Grassmann |

---

## 6. Calibration (run before any verdict)

Before reading any p-value, confirm the null is well-posed by checking, per
`(patient, band, phase)`:
1. **Strength drift** — distribution of surrogate row-sums vs observed (expected to
   drift; quantify KS so the bracket vs matched-strength is honestly characterized).
2. **Magnitude match** — surrogate `mean(W_surr[triu])` vs observed `mean(A[triu])`
   (Property §4e); large mismatch → CRC-v2.
3. **`Λ` residual** — `max_f ‖Λ̃ − Λ‖/N` from the renormalization (§5).
4. **Non-collapse** — surrogate `|ImCoh|` entry-variance not `≪` observed (guards
   against the CLT shrinkage of §5).
These four are reported as a calibration panel, not gated, but a failed calibration
voids the verdict.

---

## 7. Pseudocode (language-agnostic)

```
INPUT: patient p, band b, phases Φ = {rsPre_A, rsPre_B, task, rsPost},
       R = 200, seed, rng = Generator(seed)
FOR each phase φ in Φ:
    X      ← load_timeseries(p, φ)                 # (N, L), Pat_10 rows [53,54,55] dropped
    fs     ← sample_rate(p)
    freqs, S ← welch_csd(X, fs, nperseg_for_fs(fs))   # S: (N,N,Fall) complex Hermitian PSD
    Pwr    ← real(diag_over_channels(S))           # (N, Fall)
    Fb     ← indices where freqs in band(b)
    PRECOMPUTE per f in Fb:  C[f] = S[:,:,f] / sqrt(outer(Pwr[:,f], Pwr[:,f]))
                             (Λ[f], _) = eigh(C[f])        # eigvecs not needed (Q replaces them)
    # observed adjacency (sanity vs cache): A = mean_f |Im C[f]|

    FOR r in 1..R:
        Q ← haar_unitary(N, rng)                   # complex; one Q for all f in Fb
        accum ← zeros(N, N)
        FOR f in Fb:
            Craw ← Q @ diag(Λ[f]) @ Q^H            # eigenvalues preserved exactly
            d    ← max(Re(diag(Craw)), eps)
            Ctil ← Craw / sqrt(outer(d, d))        # unit-diagonal congruence (PSD-preserving)
            accum ← accum + abs(Im(Ctil))
        W_surr ← accum / |Fb|                       # (N,N) ≥ 0, [0,1], zero-diag
        L      ← diag(rowsum(W_surr)) − W_surr
        (λ, V) ← eigh(L)
        Dcoph_surr[φ, r] ← cophenetic_condensed_from_eigs(λ, V)

# cross-phase statistic per surrogate r, identical to audit_63:
FOR r in 1..R:
    Δtask_r ← Dcoph_surr[task, r] − Dcoph_surr[rsPre_A, r]
    Δrest_r ← Dcoph_surr[rsPost, r] − Dcoph_surr[rsPre_B, r]
    T_surr[r] ← spearman(rank(Δtask_r), rank(Δrest_r))     # ρ_split^coph (primary)
T_obs ← same statistic on observed Dcoph
p     ← mean(T_surr_finite ≥ T_obs)            # upper-tail; positive = trace
RETURN T_obs, p, calibration(W_surr vs A)
```

---

## 8. Visualization spec

- **8.1 Bracket null panel** (headline). Per band, a row of n=10 small multiples (or a
  cohort strip): observed `T_obs` (vertical line) over **two** overlaid null densities —
  CRC (this scope) and matched-strength (`audit_63`) — same x-axis. Reading rule: a
  trace is **robust** when `T_obs` sits in the right tail of *both* densities; it is
  **strength-explained** when it clears CRC but sits inside matched-strength; it is
  **realizability/spectrum-explained** when it clears matched-strength but sits inside
  CRC. Colour CRC and matched-strength consistently across all band panels.
- **8.2 Calibration panel** (§6). Per band: scatter of surrogate-mean vs observed
  per-pair `|ImCoh|` (magnitude match, identity line); KS of surrogate vs observed
  strength (drift); `Λ` residual bar. No stats text in-axes — print to stdout, push to
  companion `.md` (per `feedback_no_text_in_figures`).
- Style: `use_lrg_style()`, PDF-only full-vector, transparent, no suptitle, math `$i$,$j$`
  axis labels, ≥ 3 patients/bands shown.

---

## 9. Connection to prior tools

| Tool | Preserves | Randomizes | Manifold | Controls |
|:---|:---|:---|:---|:---|
| **Matched-strength** (`audit_62/63/65/66`) | node strength (exact) | everything else; distorts weight dist | **off** (4-cycle rewire of adjacency) | hub/strength alternative |
| **CRC (this scope)** | coherency spectrum `Λ(f)` + realizability | eigenvector identity (who-couples); `Re/Im` split | **on** (genuine coherency) | realizability + coupling-magnitude alternatives |
| Within-baseline / split-half drift (`audit_74`) | within-phase noise structure | — | n/a | task-vs-time-on-task |
| Split-baseline halves (`audit_63`) | independent rsPre halves | — | n/a | shared-baseline noise correlation |

CRC **complements** matched-strength (orthogonal invariant: spectrum+realizability vs
strength) and **does not replace** it. It **reuses** `cophenetic_condensed_from_eigs`,
the split-baseline halves convention, and the `ρ_split^coph` statistic unchanged — only
the surrogate-generation step differs. It does **not** address task-vs-time (drift null)
or spatial nuisance (probe-block control, §11).

---

## 10. Implementation plan

- **Library** (`src/lrg_eegfc/utils/surrogate/coherency_surrogate.py`, sibling of
  `matched_strength.py`; promote on 2nd caller per `coding-rules.md`):
  - `haar_unitary(N, rng) -> (N,N) complex` (Mezzadri 2007).
  - `coupling_randomized_coherency(C_band, rng) -> W_surr` where `C_band` is the
    `(|F_b|, N, N)` complex coherency stack — the core of §3.3–3.4.
  - `complex_coherency_band(X, fs, band, nperseg) -> C_band` reusing
    `welch_csd` + the `compute_imcoh` normalization (factor the shared
    `Im(CSD)/sqrt(P_i P_j)` denominator so observed and surrogate agree bit-for-bit).
  - Cache wrapper paralleling `load_or_compute_surrogate_eigs` →
    `CRC_SURROGATE_LRG_CACHE = CACHE_ROOT / "crc_surrogate_lrg"`, filename
    `{band}_{phase}_R{R}_seed{S}_imcoh_abs.npz` (no `swap` token; add `Q-shared-band` tag).
  - Reuse `cophenetic_condensed_from_eigs`, `_laplacian_eig` from `matched_strength.py`;
    `chordal_distance` etc. from `spectral.py` for the Grassmann secondary.
- **Pilot script**: `scripts/01_compute/audit/audit_96_crc_surrogate_pilot.py`
  (next free slot; `audit_88` = signed-Laplacian epi). One patient × {β, α}, R=200,
  emit `T_obs`, CRC p, calibration, and the bracket panel vs the cached matched-strength
  null. Output `data/audit/crc_surrogate/`.
- **Cohort runner** (after pilot green-light): n=10 × {β, α, γ_l} × split-baseline,
  seed `20260608`, the bracket figure as deliverable.
- **No reuse violations**: data loading via `load_timeseries` / `load_fc_matrix`; stats
  via `utils.metrics.hypothesis`; no hardcoded paths (`config.paths`).

---

## 11. Open questions

1. **`Q` granularity.** Shared-per-band (default, conservative) vs per-`f` (CLT
   shrinkage, rejected) vs per-degenerate-eigen-block (maximal entropy that respects
   `Λ` multiplicities). Calibrate §6.4 before locking.
2. **Spatial control.** A probe-block-constrained `Q` (block-diagonal-ish over shafts,
   or a rotation that preserves the same-probe `|ImCoh|` block sums) would make CRC
   *also* control spatial nuisance — closing the one negative property it shares with
   matched-strength. Deferred; needs the probe-membership map.
3. **Combined strength-matched + on-manifold null.** Importance-reweight / accept the
   CRC ensemble toward the observed strength sequence for a single null that pins both.
   Worth it **only if the two brackets disagree** at β/α; otherwise the bracket reading
   (§8.1) is the cleaner story.
4. **Imaginary-energy constraint (CRC-v2).** If calibration (§6.2) shows surrogate
   `|ImCoh|` magnitude badly mismatched to observed, constrain the `Re/Im` split (e.g.
   match `‖Im C̃‖_F` per `f` to `‖Im C‖_F`) instead of letting complex-Haar set it.
5. **Gated statistic.** Primary `ρ_split^coph` (matches the manuscript headline);
   secondary triangle `T_d` and cross-phase Grassmann `T_G`. Lock to whatever the §5
   headline uses so the bracket is apples-to-apples.
6. **`R` and seed.** `R = 200`, seed `20260608` to parallel the matched-strength
   ensemble; revisit if the calibration variance demands more draws.

---

## 12. NEXT TASK — coordinated cross-phase null → **NOW SCOPED**

> **Resolved 2026-06-11: scoped as its own report,
> [`2026-06-11_coordinated-cross-phase-null.md`](2026-06-11_coordinated-cross-phase-null.md)
> (SB-CRC — shared-backbone, deviation-randomized null).** The design below is
> superseded by that report; kept here as the origin note. Headline of the new
> scope: hold the cross-phase common coherency `B` fixed across all four surrogate
> phases, rotate only each phase's deviation `Δ_φ` independently — this lifts the
> `ρ_split` floor **off zero to exactly the shared-backbone contribution** and
> decomposes the observed trace into a shared-backbone part and a task-specific
> part. **It can bite the §5.3 headline.** Built + piloted as
> `audit_99_coordinated_null.py` (slot 98 = epi); Pat_05 β ~80% task-specific,
> clears the coordinated floor (p=0.000); cohort running.

**Status 2026-06-10: CRC + matched-strength both DONE; they agree 30/30 (§ cohort
verdict). The ONE axis neither closes is the coordinated cross-phase null.** Per
the folder rule, this lands as its **own scope report** with a fresh 5-point
preamble **before** code.

**The gap, precisely.** Both matched-strength and CRC randomize each phase
*independently*. Note this already **preserves the per-node strength trajectory**
(node `i` keeps its real strength `s_i^φ` in every phase φ, since each phase's
surrogate matches that phase's strengths on the shared node set) — yet `ρ_split`
collapses to ≈0 because per-pair *coupling identity* decorrelates across phases.
Observed clearing them therefore proves **cross-phase coupling-identity structure
exists**, but does **not** separate the alternative we still cannot exclude:

> The cross-phase `ρ_split` could be carried by a **task-independent stable coupling
> backbone** (the same higher-order structure present in pre *and* task *and* post —
> anatomy/geometry), not by **task-specific** reorganization that persists. The
> independent nulls destroy that backbone independently per phase, so they don't test
> whether it alone explains `ρ_split`.

**DESIGN SUBTLETY (already resolved — do NOT rebuild the degenerate version).**
The naive "apply the *same* rotation `O` to all four phases" is **degenerate**:
`ρ_split` is **invariant under a consistent global node relabeling** (Spearman of
per-pair deltas is unchanged when the same permutation/rotation reindexes all
phases' pairs identically), so that null returns exactly `T_obs`. Useless.

**What the construction must actually do.** Preserve the cross-phase-*stable*
coupling structure while destroying the *task-specific deviation*, then test if
`T_obs` still exceeds it. Candidate constructions to scope:
- **Shared-backbone + randomized-deviation generative model.** Decompose each phase
  `C_φ = B + Δ_φ` (B = cross-phase common component, e.g. average coherency or its
  low-rank shared part); surrogate `C̃_φ = B + R_φ(Δ_φ)` where `R_φ` randomizes the
  deviation's coupling identity (e.g. CRC-rotate only `Δ_φ`, independent per phase)
  while `B` is shared verbatim. Preserves the stable backbone's cross-phase
  correlation; kills task-specific persistence. If `T_obs` clears this, the trace is
  task-specific beyond the stable backbone.
- **Partially-correlated rotation family** interpolating independent↔fully-coordinated
  (e.g. `O_φ = exp(θ·G_common + (1−θ)·G_φ)`), sweeping θ — but anchor on the
  generative-backbone version since the pure-rotation endpoints are degenerate/≈0.

**Why CRC is the right substrate.** Unlike matched-strength (where "same structural
permutation across phases" is ill-defined — `audit_63` note), CRC's congruence makes
a shared component `B` and per-phase deviations cleanly separable and recombinable on
the realizable manifold.

**Deliverable:** scope report → pilot (Pat_05 β, reuse `audit_97` scaffolding +
`coherency_surrogate.py`) → cohort. Reuse everything in `audit_97`; only the
surrogate-generation step changes. Outputs under
`data/audit/coordinated_cross_phase_null/`. Slot: `audit_99` (`audit_98` was already
taken by `audit_98_epi_seed_prediction_pipeline.py`; 96/97 = CRC). **Done — see the
resolution note at the top of this section and `2026-06-11_coordinated-cross-phase-null.md`.**
