---
name: multiscale-notions-propagator-vs-grassmann
type: guide
era: COHORT_N10
status: current
created: 2026-07-04
updated: 2026-07-04
pointers:
  - .agents/guides/02_methods/lrg-framework-guide.md
  - .agents/guides/02_methods/linearity-and-higher-order-structure.md
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
  - data/audit/grassmann_cluster_extent/cohort_summary.csv
  - .agents/reports/2026-06-26_per-band-phenomenology-vision.md
---

# Two "multiscale" notions: the LRG propagator (cophenetic hierarchy) vs the Grassmann rank-sweep

## Renormalization head

We use **two** trace probes and call both "multiscale", but they are multiscale in
**different senses** and must not be presented as independent confirmations. The
**cophenetic ρ_split** (LRG propagator → dendrogram) is multiscale through a
**hierarchy of merge heights** and uses the **eigenvalue-weighted** diffusion
`e^{−τL}`. The **Grassmann `d_G(k)`** is a **single-scale subspace metric at fixed
rank k**; it becomes multiscale only through the **k-sweep + cluster-extent test**,
and it is **eigenvalue-blind** (it sees only the *span* of the top-k eigenvectors).
Both are windowings of the **same Laplacian spectrum** — a soft exponential window
(propagator, indexed by τ) vs a hard rank cutoff (Grassmann, indexed by k) — so they
are **correlated by construction, not independent**. Their **agreement** (β) is
robustness-to-windowing; their **dissociation** (α = per-pair-only, γ_l =
subspace-only) is the real result. **Recommended manuscript framing: "two spectral
views of the LRG structure — the per-pair hierarchy and the whole-network subspace",
NOT "two independent multiscale measures".**

---

## 1. What each probe actually is

**Cophenetic `ρ_split` (the "propagator" probe).**
Build the graph Laplacian `L`, form the diffusion affinity from the LRG propagator
`e^{−τL}` at `τ = 1/λ_max` (fine-scale operating point, τ-robust — see
[[tau_sensitivity_trace_2026_06_22]]), average-linkage → **dendrogram** → cophenetic
ultrametric `D_coph` (the merge height at which each pair joins). `ρ_split =
Spearman(D_tt − D_preA, D_post − D_preB)` over the per-pair condensed vectors.
- **Scale axis:** the **dendrogram hierarchy** — one cophenetic distance per pair
  already encodes *all* merge heights (fine pairs → coarse communities). Natively
  multiscale.
- **Uses:** eigenvectors (via clustering) **and** eigenvalues (via the `e^{−τλ}`
  diffusion weighting). Full spectrum, softly weighted.

**Grassmann `d_G(k)` (the "subspace" probe).**
Take the leading `k` eigenvectors of each phase → a `k`-plane on the Grassmann
manifold `Gr(k, N)`. Chordal distance between task and rest subspaces
`d_c(U,V) = ‖UUᵀ − VVᵀ‖_F / √2 = (Σ_i sin²θ_i)^{1/2}`, `θ_i` = principal angles.
Trace statistic `T_G` (lower tail = trace). Verdict via **cluster-extent
permutation over the k-sweep** (longest contiguous-significant run / cluster mass).
- **Scale axis:** the **spectral rank k** — small `k` = coarse global modes
  (Fiedler / communities), large `k` = finer modes added.
- **Uses:** only the **span** of the top-k eigenvectors. **Eigenvalue-blind** within
  the top-k (eigenvalues only *order* the modes); no diffusion weighting.

## 2. Is the Grassmann chordal distance "really" multiscale? — yes, but weakly

- The **chordal metric at fixed k is single-scale** (one subspace comparison, one
  rank). By itself it is NOT multiscale.
- The **`d_G(k)` k-sweep is multi-resolution**: sweeping `k` traverses coarse→fine
  spectral content, and the cluster-extent test asks "is the trace sustained across a
  *band* of scales." That aggregation is a legitimate multiscale test.
- **But it is not a renormalization.** `k` is a brick-wall spectral truncation with
  no semigroup / coarse-graining kernel. The RG-proper object is the propagator
  `e^{−τL}` (`e^{−τ₁L}e^{−τ₂L} = e^{−(τ₁+τ₂)L}`). Even our cophenetic *trace* is not
  literally an RG flow either — its multiscale-ness is the **linkage hierarchy** —
  but it sits closer to the RG because its affinity **is** the diffusion propagator.

**Verdict:** multiscale in the **multi-resolution / spectral-rank** sense — yes;
multiscale in the strict **RG / semigroup** sense — no. The propagator is the
eigenvalue-weighted, full-spectrum probe; the Grassmann is the coarser,
subspace-only, rank-truncated probe.

## 3. The unifying view (the "peace")

Both are **filter functions on the identical eigenbasis** of `L`:

| probe | filter on mode λ | scale index | scale mechanism | spectrum used |
|---|---|---|---|---|
| cophenetic ρ_split | `e^{−τλ}` (soft low-pass) | τ / merge height | diffusion + linkage hierarchy | eigvecs **and** eigvals |
| Grassmann `d_G(k)` | `𝟙[λ ≤ λ_k]` (hard cutoff) | k (integer rank) | rank sweep + cluster-extent | eigvecs only (span) |

Because a diffusion time `τ` effectively retains modes up to `λ ≈ 1/τ`, **τ and k are
monotonically linked** (large τ ↔ small k = coarse; small τ ↔ large k = fine). Same
coarse↔fine spectral axis, reached by soft vs hard windowing.

**Claim / do-not-claim:**
- ✅ Both live in the **LRG / Laplacian framework** — neither is raw FC
  ([[feedback_results_only_in_laplacian_framework]]).
- ✅ **Convergence** of the two = the trace is robust to *how you window the spectrum*
  (soft vs hard, τ vs k).
- ✅ **Dissociation** = the reorganization is expressed in one spectral aspect but not
  the other (see §4).
- ❌ Do **not** call them "independent" multiscale confirmations. They share the
  eigenbasis → correlated by construction. Agreement is corroboration of robustness,
  not two experiments.
- ❌ Do **not** imply the Grassmann k-sweep is a renormalization flow.

## 4. Empirical per-band picture (why this matters)

Per-band cohort verdict of each probe vs its own matched-strength null
(`✓` clears, `~` weak/trend, `·` none):

| band | raw per-edge | ρ_coph (per-pair hierarchy) | Grassmann `d_G(k)` (subspace) | reading |
|------|:---:|:---:|:---:|---|
| **β**   | ~ (p=.053) | **✓ (p=.005)** | **✓ (p=.005, LOO-robust)** | **both multiscale views → flagship** |
| **α**   | ✓ (.014) | **✓ (.002)** | · (.348) | **per-pair / hierarchy only** |
| **γ_l** | ~ (.053) | · (.116) | **✓ (.005)** | **subspace / whole-network only** |
| δ       | ✓ (.042) | · (.278) | ~ (.005 but LOO .055) | raw floor + fragile subspace |
| γ_h     | · (.188) | · (.246) | ~ (.060) | subset only |
| θ       | ~ (.138) | · (.722) | · (.159) | **absent in every view** |

- **β is the only band both multiscale probes clear** → the reorganization is both a
  **subspace rotation** (γ_l-type signature) **and** a **stable per-pair hierarchy**
  (α-type signature). This is *why* β is the flagship, from a second angle.
- **α = ρ_coph-only:** α reshuffles pairwise ultrametric relations but does **not**
  rotate the dominant eigen-subspace (fine/pairwise change, no global mode
  reorientation).
- **γ_l = Grassmann-only:** γ_l rotates the dominant subspace (a global/mesoscale mode
  change) but leaves **no** stable per-pair cophenetic hierarchy — consistent with the
  VERDICT_LEDGER "γ_l: strong trace, only Grassmann".

Sources: `matched_strength_surrogate_split_baseline/cohort_summary.csv` (ρ_coph),
`raw_vs_multiscale/` (raw, Q2), `grassmann_cluster_extent/cohort_summary.csv`
(Grassmann, `cluster_p_cluster_mass` + `cluster_p_mass_loo_max`).

## 5. Manuscript wording (recommended)

- **Preferred (option B):** "two spectral views of the LRG structure — the **per-pair
  hierarchy** (cophenetic ultrametric) and the **whole-network subspace** (Grassmann
  `d_G(k)`)". Sidesteps the RG-vs-resolution debate; unimpeachable.
- **Acceptable (option A):** keep the word "multiscale" but **define it once** —
  "multiscale = probed across a sweep of spectral scales (τ-diffusion hierarchy for
  ρ_coph; subspace rank k for Grassmann)".
- **Forbidden:** "two independent multiscale measures" (they share the eigenbasis);
  "the Grassmann k-sweep is a renormalization" (it is a spectral-rank truncation).

## 6. Open question / caveat

- The chordal `d_G(k)` is eigenvalue-blind within the top-k. A natural
  eigenvalue-weighted subspace metric (e.g. weighting principal angles by `e^{−τλ}`)
  would interpolate between the two probes and could test whether the α/γ_l
  dissociation is an artifact of the hard cutoff. **Not run.** Until then, read the
  dissociation as "the reorganization lives in different spectral windows", not as a
  proven mechanistic split.
- Figure `fig_per_band_phenomenology.py` panel D currently shows this comparison; its
  label should follow §5 (option B).
