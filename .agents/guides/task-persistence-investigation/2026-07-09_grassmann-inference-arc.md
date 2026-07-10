---
name: grassmann-inference-arc
type: scope
era: IMCOH_ABS / COHORT_N10
status: result
created: 2026-07-09
updated: 2026-07-09
pointers:
  - scripts/01_compute/audit/audit_165_grassmann_inference_arc.py
  - data/audit/grassmann_inference_arc/cohort_summary.csv
  - .agents/guides/task-persistence-investigation/2026-06-12_cophenetic-consolidation-arc.md
  - scripts/01_compute/audit/audit_66_grassmann_matched_strength_surrogate.py
  - scripts/01_compute/audit/audit_70_grassmann_cluster_extent.py
  - scripts/01_compute/audit/audit_152_consolidation_arc_rhosym.py
  - scripts/01_compute/audit/audit_113b_inference_lenmatched_null.py
  - .agents/preprint/directives/methods_grassmann_cluster_extent.md
---

# Grassmann subspace arc — a second, strength-independent witness for the inference-specific β trace

**Head.** The inference-specific β trace — *rest_post preferentially consolidates
the reorganization the brain adds when it must **infer** novel relations
(`task_test` beyond `task_learn`), not the mere encoding of premise pairs* — is
currently carried by a **single** probe: the cophenetic per-pair correlation
`T_infspec·e = ρ_sym(f, p · e)` (β obs med +0.091, surr med +0.005, Wilcoxon
p = 0.0098, LO-P15 0.020; audit_152). The whole-task β trace, by contrast, is
**doubly-witnessed** — cophenetic *and* Grassmann subspace geometry — which is
what makes it defensible. This scope defines the **Grassmann analog of the
inference arc**, so the flagship inference claim earns the same second,
structurally-independent witness from a probe that is blind to the linkage /
ultrametric machinery cophenetic depends on. Pure reuse: no new connectivity
estimator, and the matched-strength surrogate eigen-ensembles it needs
(`task_learn`, `rest_pre_B`, `task_test_lm_head`) **already exist** on disk from
the cophenetic-arc scripts (seed 20260511).

This is a **result-first** investigation. The deliverable is a yes/no: does an
eigenmode-subspace probe independently reproduce the β-only inference-specific
consolidation, and does it hold under length-matching?

> **VERDICT (2026-07-09, audit_165). NEGATIVE — clean, informative, on-thesis.**
> The Grassmann leading-*k* subspace probe does **not** second-witness the
> cophenetic inference-specific β trace. Standard `T_G^infspec` is **null in all
> six bands** (β `cluster_p_mass = 0.408`, `T_G* = 0.032`, 5/10 positive median,
> 4/10 above own surrogate; δ 0.756 / θ 0.612 / α 0.821 / γ_l 0.378 / γ_h 1.000;
> per-patient β medians scatter symmetrically over [−0.20, +0.31] with no cohort
> direction). **Not a power miss:** the probe resolves cross-phase structure —
> `T_G^onl` (whole-task) reproduces audit_70 to the decimal (**δ 38.07 / β 69.76 /
> γ_l 66.14**, all strong; cross-check max dev 1.7e-13) and `T_G^enc` (encoding
> echo) clears δ and β (0.005). **Reading:** `rest_post`'s leading-subspace
> relocation is toward *the task in general* (encoding), not specifically toward
> inference; the inference-minus-encoding refinement lives only in the fine
> cophenetic pair-identity / merge structure, below leading-eigenmode resolution.
> So the inference-specific β trace **stays singly-witnessed** (cophenetic ρ_sym
> p=0.0098, audit_152) — the hoped-for second witness did not materialise. The
> negative is **on-thesis**: a coarse/global spectral probe cannot see the
> inference signal, reinforcing that the *fine multiscale* (cophenetic) structure
> is necessary. Grassmann `enc` also fails to reproduce the cophenetic
> α-encoding channel (α `enc` p=0.552), tracking δ/β like `onl` — the subspace
> probe does not resolve the encoding/inference dissociation at all.
>
> **Caveat (iii) is INVERTED by the data — length-matching must NOT be applied to
> this probe.** `T_G^infspec_lm` fires "strong" (p=0.005) in **all six bands
> including θ** (masses 206–273, runs 62–104/111), and in θ/α the *observed* cohort
> median is **negative** while the verdict reads "strong" with 10/10 above
> surrogate — impossible for a real positive trace. Head-truncating `task_test`
> (fewer Welch segments) shifts its matched-strength surrogate baseline
> systematically negative, manufacturing a one-sided obs>surr artifact everywhere.
> The **raw** `T_G^infspec` is the interpretable arm and is length-robust *by
> construction* (the `task_test`-length bias enters `Δ(RP)` and `Δ̄(pre)` through
> the same `−d_G(·,TT)` term and cancels in the difference). If a length-matched
> arm is ever wanted, the surrogate must carry the truncation on both sides at the
> Δ level — the cached `task_test_lm_head` ensemble does not. Result memory:
> [[grassmann_inference_arc_negative_2026_07_09]].

---

## Critical preamble (5 points, before code)

**(1) Claim.** Across `rest_pre → task_learn → task_test → rest_post`, the
leading-*k* Laplacian eigenmode subspace of `rest_post` sits **preferentially
closer to the inference subspace (`task_test`) than to the encoding subspace
(`task_learn`)**, beyond the rest-baseline, in β and not in the null bands. Strong
form: the inference-specific subspace relocation `T_G^infspec` clears
matched-strength in β only, corroborating the cophenetic `T_infspec·e` dissociation
(α consolidates encoding; β consolidates encoding **+** inference-specific).

**(2) Null.** Matched-strength 4-cycle ±δ surrogate (R = 200, `SWAP_FACTOR = 20`,
seed 20260511), generated per phase, reusing the **exact** shared eigen-ensemble
(`load_or_compute_surrogate_eigs` → `data/cache/matched_strength_surrogate_lrg/`)
already built for the Grassmann whole-task probe (audit_66) and the cophenetic arc
(audit_103/110/111). Chordal distances recomputed on surrogate eigvecs; per-*k*
paired one-sided Wilcoxon (obs > surrogate mean); cluster-mass empirical
p = (1 + #{null ≥ obs}) / (R + 1). Matched-strength is mandatory before any cohort
claim (matched-strength rule).

**(3) Strongest plausible alternatives.**
  (i) **Encoding leakage** — `task_test` and `task_learn` subspaces share gross
      structure, so an apparent inference-specific relocation could be encoding
      bleeding through.
  (ii) **Distance-difference noise** — `Δ(r) = d_G(r,TL) − d_G(r,TT)` differences
      two noisy scalar distances; a weak `T_G^infspec` could be estimator noise,
      not absence of inference-specific persistence. Grassmann is *scalar per
      patient per k*, hence lower-powered than the per-pair cophenetic ρ.
  (iii) **Phase-duration / SNR asymmetry** — `task_test` is 1.4–2.5× longer than
      `task_learn`; longer records give better-estimated eigenvectors, so
      `d_G(·,TT)` and `d_G(·,TL)` differ in stability *independent of cognition*.
      Because `Δ` **directly differences** distances to TL and TT, this bias hits
      the Grassmann arm harder than the cophenetic arc.
  (iv) **Strength drift, not subspace reorganization** — the relocation could be
      per-node strength change rather than eigenmode-geometry change.
  (v) **Eigenvalue near-degeneracy** — near-degenerate Laplacian modes swap order
      across phases, inflating `d_G` at the affected *k* independent of structure.

**(4) Null's mechanical reach.** Matched-strength fixes the per-node strength
sequence exactly, so a surviving `T_G^infspec` is **not** explained by (iv). It
does **not** reach (i)–(iii) or (v). Those require, respectively: (i) the contrast
is **inference-minus-encoding by construction** — `Δ` differences closeness-to-TT
against closeness-to-TL, so shared task structure cancels and no post-hoc partial
is needed (the scalar geometry cannot support a stable n=10 partial Spearman
anyway); (ii) a **within-rest noise floor** — compare `|T_G^infspec|` to the
split-half `d_G(preA, preB)` band; abstain if inside it; (iii) the **length-matched
control** — recompute `U_k^{TT}` from the head-truncated `task_test`
(`_tt_fc_lenmatched`, n = `n_TL` samples, unchanged `nperseg`) and require β to
survive; (v) the **k-grid sweep + cluster-mass** (resilient to isolated bad *k*)
plus reporting the spectral gap at gated *k*.

**(5) Falsification + limitations.** "Inference-specific offline subspace
consolidation" requires, cohort-wide: `T_G^infspec(k) > 0` (per-*k* Wilcoxon,
one-sided) over a contiguous *k*-band with `cluster_p_mass < 0.05`, **AND**
β-selectivity (null bands do not clear), **AND** survival of length-matching.
It is **falsified** if `rest_post` is equidistant from TL and TT beyond baseline
(`T_G^infspec ≈ 0`), or if any apparent β effect collapses under length-matching,
or if the effect is not β-selective (fires in θ, the absent band).
**Limitations:** Grassmann is a *distance-relocation* probe — it answers "did the
subspace move toward inference," not "did the same pairs reorganize" (cophenetic's
question); the two can disagree, so **non-corroboration is ambiguous** (could be a
probe blind-spot or a power deficit), whereas **corroboration is strong** (two
independent geometries agree). Scalar-per-patient ⇒ no per-pair localization, so
this measure **cannot** localize the inference trace to OFC/cingulate — that needs
a separate per-node subspace-leverage construction (open question). n = 10.

---

## Notation

- `W_x ∈ ℝ^{N×N}` — `imcoh_abs` FC adjacency for phase `x` (zero-diag, clipped
  `[0,1]`, symmetrized), via `load_fc_matrix` / `load_phase_fc`.
- `L_x = D̂_x − W_x` — combinatorial Laplacian (Villegas canonical; `L_rw`, `L_sym`
  rejected, methods_grassmann §1). `eigh(L_x) → (λ, Φ)` ascending; trivial mode
  `φ_1` (λ_1 = 0) discarded.
- `U_x^{(k)} = Φ_x[:, 1:k+1] ∈ ℝ^{N×k}` — the *k* slowest non-trivial eigenmodes
  (leading invariant subspace).
- **Chordal (Grassmann) distance**, `k`-truncated (audit_66:235–256, audit_70:67–91):
  `d_G(x, y; k) = √( k − ‖ U_x^{(k)ᵀ} U_y^{(k)} ‖_F² )`,   domain `[0, √k]`;
  basis/rotation-invariant. Equivalent SVD form `√(k − Σ_i σ_i²)`, `σ_i` = singular
  values of `U_xᵀU_y` (= cosines of principal angles).
- Phases: `rest_pre` split into `preA`, `preB` (unbiased baseline, audit_63
  `ensure_half_fcs`); `task_learn` (TL); `task_test` (TT); `rest_post` (RP).
- *k*-grid: `K = {2, 3, …, 112}` (111 values), matching the whole-task probe.

## Definitions

**Inference-minus-encoding closeness of a reference phase `r`** (per *k*):
> `Δ(r; k) := d_G(r, TL; k) − d_G(r, TT; k)`   ∈ `[−√k, √k]`.
Positive ⇒ `r`'s subspace sits closer to the inference phase than to the encoding
phase.

**Arc functionals** (per *k*; sign convention **positive = trace**, locked
`T_d > 0`):

| Functional | Formula (per *k*) | Cophenetic analog |
|---|---|---|
| encoding echo | `T_G^enc(k)  = d_G(preA, TL; k) − d_G(TL, RP; k)` | `T_learn = ρ(e, p)` |
| inference-online | `T_G^onl(k)  = d_G(preA, TT; k) − d_G(TT, RP; k)` | `T_test  = ρ(g, p)` (= whole-task `T_G`) |
| **inference-specific** | `T_G^infspec(k) = Δ(RP; k) − Δ̄(pre; k)` | `T_infspec·e = ρ_sym(f, p·e)` **(HEADLINE)** |

where the **split-symmetrized baseline** mirrors `ρ_sym`:
> `Δ̄(pre; k) = ½ [ Δ(preA; k) + Δ(preB; k) ]`.

**Algebraic identity (motivation, not a second definition).** With the single-half
baseline `preA`, `T_G^onl − T_G^enc = Δ(RP) − Δ(preA)` exactly (`d_G` symmetric).
So `T_G^infspec` is simultaneously (a) the **difference-of-triangles**
inference-online-minus-encoding-echo, and (b) the **relative-closeness** shift of
RP toward inference-over-encoding beyond baseline — the same object from two
readings. The split-symmetrization replaces the single `preA` anchor by the mean
over both halves, exactly as `ρ_sym` averages the two arm assignments.

**Cohort aggregation over the *k*-grid** (reuse audit_70 verbatim):
- per-*k* paired one-sided Wilcoxon `H1: T_G^·(k) > 0` across the 10 patients
  (obs vs surrogate mean) → `p_k`.
- `cluster_mass(b) = Σ_{k : p_k < α_k} (−log10 p_k)`, `α_k = 0.05` (audit_70:161).
- normalized `T_G^*(b) = cluster_mass / (n_k · log10(R+1))` ∈ `[0,1]`,
  `n_k = 111`, `R = 200` (methods_grassmann §5b).
- empirical `cluster_p_mass = (1 + #{null_mass ≥ obs_mass}) / (R+1)`; phantom-
  surrogate null (each surrogate as obs, mean of the rest as reference).
- **gate (locked Decision 8, mass-only):** `strong ⇔ p_mass < 0.01`;
  `weak ⇔ 0.01 ≤ p_mass < 0.05`; `no trace ⇔ p_mass ≥ 0.05`; `cluster_p_LR`
  descriptive only. Strong-tier precondition `p_mass_loo_max < 0.05`.

## Properties / sanity contracts

- **Cross-check gate (anti-hallucination):** `T_G^onl(k)` computed here MUST equal
  the whole-task Grassmann `T_G(k)` of audit_66/70 per (patient, band, *k*) to
  ~1e-10 — identical inputs (`preA, TT, RP`), identical `chordal`. Script asserts
  and prints max deviation before any new number is trusted (mirrors audit_103's
  `T_test`↔audit_83 assertion).
- **Sign convention locked:** positive = trace, all three functionals.
- **Invariance:** `d_G` is invariant to within-subspace rotation and eigenvector
  sign/order *within a degenerate block* (chordal is a function of principal
  angles only) — so no eigenvector-sign bookkeeping is needed. It is **not**
  invariant to which modes fall inside the leading-*k* cut when eigenvalues are
  near-degenerate at the cut (caveat (v)).
- **Range / degeneracy:** `T_G^infspec(k) ∈ [−2√k, 2√k]`; grows with `k`, so the
  raw functional is **not** comparable across *k* — this is why aggregation is on
  the *significance* (`p_k`), not the raw magnitude (as in audit_70).
- **Negative property (contract):** blind to merge-height / ultrametric structure
  and to within-`U_k` reorganization; measures *subspace relocation*, not
  *pair-identity* reorganization. Cannot localize (scalar).
- **Complexity:** per (patient, band): 5 observed `eigh` (O(N³)) + one
  `AᵀB` per phase-pair per surrogate; the cumulative-Frobenius trick (audit_70)
  gives all 111 *k* from one `(N−1)×(N−1)` product. Dominated by cache reads;
  effectively free vs the cophenetic arc.

## Caveats & failure modes

| # | Failure mode | Mitigation (in-script) |
|---|---|---|
| (i) | encoding leakage TL↔TT | `Δ` differences TL vs TT → shared structure cancels by construction; report `T_G^enc` alongside to show the dissociation |
| (ii) | distance-difference noise; low power (scalar/patient) | within-rest floor `d_G(preA, preB; k)`; abstain if `\|T_G^infspec\|` inside floor; cluster-mass over 111 *k* recovers power; **flag non-corroboration as ambiguous, never as refutation** |
| (iii) | `task_test` 1.4–2.5× longer than `task_learn` biases `d_G(·,TT)` vs `d_G(·,TL)` | **mandatory** length-matched recompute: `U_k^{TT}` from `_tt_fc_lenmatched` (head-truncate to `n_TL`, unchanged `nperseg`); require β to survive; cached `*_task_test_lm_head_R200_…` surrogates already exist |
| (iv) | strength drift | matched-strength null fixes strength sequence exactly |
| (v) | eigen near-degeneracy at the *k*-cut | *k*-grid sweep + cluster-mass (resilient to isolated bad *k*); report spectral gap `λ_{k+1}−λ_k` at gated *k* |
| (vi) | root/whole-graph anatomy dominates coarse *k* | same-probe / epi controls are downstream; here rely on β-selectivity vs null bands as the specificity check |

## Pseudocode

```
K = [2..112];  COHORT = 10 patients;  BANDS = 6;  R = 200
for pat in COHORT:
  for band in BANDS:
    # ---- observed ----
    for x in [preA, preB, TL, TT, RP]:
        Phi[x] = eigvecs(L(load_phase_fc(pat, x, band)))      # eigh ascending, drop col 0
    TT_lm  = eigvecs(L(tt_fc_lenmatched(pat, band)))          # length-matched TT (caveat iii)
    for k in K:
        d = { (a,b): chordal(Phi[a], Phi[b], k) for needed (a,b) pairs }
        Δ_RP  = d[RP,TL] - d[RP,TT]
        Δ_pre = 0.5*((d[preA,TL]-d[preA,TT]) + (d[preB,TL]-d[preB,TT]))
        obs_infspec[k]  = Δ_RP - Δ_pre
        obs_enc[k]      = d[preA,TL] - d[TL,RP]
        obs_onl[k]      = d[preA,TT] - d[TT,RP]           # == whole-task T_G (assert vs audit_66)
        obs_floor[k]    = chordal(Phi[preA], Phi[preB], k)
        # length-matched variant reuses TT_lm in place of Phi[TT]
    # ---- surrogate (cached eigvecs, matched-strength) ----
    for r in 1..R:
         Psi[x][r] = surrogate_eigvecs(pat, band, x, r)       # load_or_compute_surrogate_eigs
        ... same k-loop → surr_infspec[k, r], surr_enc, surr_onl
cohort:
  for each functional, band:
    p_k       = wilcoxon_greater(obs[:, k] , surr_mean[:, k])   for k in K   # audit_70:140
    obs_mass  = cluster_mass(p_k, α=0.05)
    p_mass    = (1 + #{null_mass ≥ obs_mass}) / (R+1)           # phantom-surrogate
    gate      = strong/weak/none  by Decision-8 thresholds
    p_mass_loo_max = max over leave-one-patient-out               # strong-tier precondition
validate: assert max_k,pat | obs_onl(k) − audit66_TG(k) | < 1e-6
report:  per band × {enc, onl, infspec, infspec_lenmatched} → p_mass, T_G*, gate, LR, floor-fraction
```

## Visualization spec

- **Arc bar panel** (headline), per band cohort: `T_G^*` (normalized cluster-mass)
  for `enc` vs `onl` vs `infspec`, with the β column annotated `p_mass`. Reading
  rule: β-only tall `infspec` bar next to a broad `enc` bar (α + β) = the same
  dissociation the cophenetic arc shows, now in subspace geometry.
- **Per-*k* significance ribbon**, β: `−log10 p_k` vs `k ∈ [2,112]` for `infspec`
  (observed) with the `α_k = 0.05` line and the shaded within-rest floor band;
  the length-matched curve overlaid dashed. Reading rule: a contiguous run above
  the line that *survives the dashed length-matched overlay* = a real inference-
  specific subspace band; loss under dashing = a length artifact.
- **Cross-phase `d_G` 5×5 heatmap** (`preA, preB, TL, TT, RP`), cohort-median,
  β, at a representative gated `k` — the subspace trajectory geometry; mirror of
  the cophenetic 5×5.
- No pre-registered gate line drawn per patient (never/always list); mark every
  significant *k*-cell uniformly. `use_lrg_style`, PDF-only, transparent.

## Connection to prior tools

Extends the **whole-task Grassmann probe** (audit_66/70; `T_G` collapses the task
to `task_test`) to the **four-phase arc**, exactly as the cophenetic
consolidation arc (audit_103/152) extended the whole-task cophenetic trace
(audit_83). Same Laplacian, same chordal distance, same `T_G*` cluster-mass gate,
same matched-strength ensemble, same length-matched control (audit_113b). It is
the **subspace-geometry sibling** of `T_infspec·e = ρ_sym(f, p·e)`: complementary
(distance-relocation vs pair-identity correlation), not a replacement. Orthogonal
to the WM/epi node-stratifications (those re-slice *nodes*; this re-slices
*phases*, like the cophenetic arc). Does **not** subsume the cophenetic arc —
corroboration is the goal, and a genuine dual-probe result requires both.

## Implementation plan

- **Script:** `scripts/01_compute/audit/audit_XXX_grassmann_inference_arc.py`
  (next free audit index). Clone the audit_70 skeleton; swap STEP-1 (its
  `chordal_distances_for_k_grid` triangle) for the three arc functionals above;
  keep STEP-2/3 (per-*k* Wilcoxon, `cluster_mass`, phantom-surrogate `p_mass`,
  Decision-8 gate, LOO) verbatim.
- **Library reuse (no forks):** `chordal` / `topk_basis` / `laplacian_eig` —
  promote from audit_66 into `src/lrg_eegfc/utils/metrics/spectral.py` (this audit
  is their second caller ⇒ promotion required by coding-rules). `cluster_mass`,
  `wilcoxon_per_k_greater` likewise → `hypothesis.py` if not already there.
  Surrogate: `matched_strength.load_or_compute_surrogate_eigs`. FC + halves:
  `load_fc_matrix`, `ensure_half_fcs`, `_tt_fc_lenmatched` (audit_113b). **No new
  connectivity or surrogate code** — all ensembles cached (seed 20260511):
  `task_learn`, `rest_pre_B`, `task_test_lm_head` confirmed present for all bands
  at R=200 (β also R=1000).
- **Outputs:** `data/audit/grassmann_inference_arc/{per_patient,cohort_summary}.csv`
  + `README.md`; figure PDFs to `data/reports/…` (deferred until numbers read
  with the user).
- **Compute estimate:** cache-read-bound; ~5 observed `eigh` per (patient, band)
  + Frobenius cumsum. Expect < 2 min cohort×bands; time 1–2 cells and extrapolate
  before full launch (optimize-and-surface rule), print `[i/N] elapsed ETA`.

## Open questions

- **Inference contrast (load-bearing, pending user sign-off).** `T_G^infspec =
  Δ(RP) − Δ̄(pre)` bakes the encoding control into a *distance difference*. The
  documented alternative is a **principal-angle projection**: decompose the
  `preA → RP` subspace displacement onto the encoding axis (`preA→TL`) vs the
  inference axis (`TL→TT`) via principal vectors, and test the inference-axis
  component. Richer, but introduces tangent-space bookkeeping and its own
  identifiability caveats. Default = the closeness contrast (parameter-free,
  reuses `chordal` unchanged); the projection variant is deferred unless the
  closeness contrast is judged to under-isolate inference.
- **Localization.** Grassmann is scalar ⇒ no OFC/cingulate localization of the
  inference-specific component. A per-node subspace-**leverage** score (drop-node
  influence on `d_G`) would be a separate scope report if the arc result warrants
  anatomical follow-up.
- **k default for the trajectory heatmap** — pick the gated-band centroid *k*
  post-hoc; not a free parameter of the test (the test is *k*-swept).
- **Behavioral TC1** — as with the cophenetic arc, no behavior, so the trace
  cannot be linked to inference *success*.
