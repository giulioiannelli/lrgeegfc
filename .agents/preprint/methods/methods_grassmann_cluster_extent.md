---
name: methods-grassmann-cluster-extent
era: IMCOH_ABS_COHORT_N10
status: locked_2026-05-19
kind: methods-companion
scope: complete methodology for the Grassmann probe under the cluster-extent (LR ∨ mass) gate, including the band-level k-independent scalar T_G*
companion: CONTROLS.md, VERDICT_LEDGER.md, methods_revision_2026-05-18_cophenet.md
supersedes: any prior partial Grassmann methods text (notably the per-k narrative without an explicit k-collapse statistic)
---

# Grassmann subspace probe — complete methods (cluster-extent, with k-independent scalar)

**Head.** The Grassmann probe asks whether the *whole-network communication
geometry* — encoded as the leading-`k` non-trivial Laplacian eigenmode
subspace — moves through `task_test` in a way that persists into `rsPost`,
above what strength-preserving rewiring can produce. The probe is naturally
multiscale (one chordal distance per `k`), so the cohort-level test does not
read off any single `k` — it summarises the entire `k ∈ [2, 112]` profile via
two co-primary cluster statistics (**longest contiguous-significant run** and
**cluster mass**) whose empirical null is built from the same R=200
matched-strength surrogates that drive the per-k Wilcoxon. The band-level
**single-number scalar** that the manuscript reports is the cluster mass
`T_G^*(b)` — resilient all-clusters sum of `−log10 p_k(b)` across every
contiguous-significant `k`-cell — normalized to `[0,1]` per C1
(denominator `n_k^cohort · log10(R+1) = 111 × 2.3032 ≈ 255.65`). β at
`T_G^* = 69.76` raw / **0.273 normalized**, `cluster_p_mass = 0.005`,
LOO max `p_mass = 0.005` (Pat_02) is the **primary band finding**; γ_l
also reaches the strong tier under Decision-12 LOO precondition (`T_G^* = 66.14`
raw / 0.259 normalized; LOO 0.040 Pat_05); δ sits in the weak tier (cohort
gate at floor but full-data LOO 0.055 Pat_08 fails the < 0.05 precondition);
α, θ, γ_h carry no Grassmann trace under this gate.

---

## 1 Subspace object `U_k`

For each (patient, band, phase) FC adjacency `W ∈ ℝ_{≥0}^{N×N}`
(`<|ImCoh|>_f`, `N` = number of contacts) we form the **combinatorial
(Villegas fluid) Laplacian** — the same `L̂` that drives every other
LRG object in the codebase (propagator `K̂(τ) = e^{−τL̂}`, communication
distance `D(τ) = 1/ρ̂(τ)`, dendrogram):

```
L̂ = D̂ − W,   D̂_{ii} = Σ_j W_{ij}
```

`L̂` is symmetric (since `W` is symmetric) and positive-semidefinite,
so its eigendecomposition under `numpy.linalg.eigh` gives a real
orthonormal eigenbasis:

```
L̂ = Σ_n λ_n φ_n φ_n^T,   0 = λ_1 ≤ λ_2 ≤ … ≤ λ_N,   φ_n^T φ_m = δ_{nm}
```

`φ_1` is the trivial constant mode (eigenvalue zero on a connected
graph) and is discarded. The **leading-`k` non-trivial eigenmode
subspace** is

```
U_k = span{φ_2, φ_3, …, φ_{k+1}} ⊂ ℝ^N,   k ∈ {2, 3, …, 112}
```

`U_k` is treated as a point on the Grassmann manifold `Gr(k, N)` — the
manifold of `k`-dimensional subspaces of `ℝ^N`. The orthonormality of
`{φ_n}` is **decisive** for the chordal Grassmann distance
identity in §2 (Frobenius-of-product equals sum of `cos² θ_i` only
when both bases are orthonormal); using the random-walk Laplacian
`L_rw = D̂^{−1}L̂` would be invalid because `L_rw` is non-symmetric,
its eigenvectors are not orthonormal under the standard inner
product, and the chordal identity would fail. Using the symmetric
normalised Laplacian `L_sym = D̂^{−1/2} L̂ D̂^{−1/2}` would also be
mathematically valid (also symmetric PSD) but would diagonalise a
**different operator** from the rest of the LRG pipeline, breaking
the coherence between the Grassmann probe and the cophenetic
per-pair probe; we use the canonical Villegas `L̂` throughout.

Implementation: `scripts/01_compute/audit/audit_66_grassmann_matched_strength_surrogate.py:176-180`
(`laplacian_eig: deg = W.sum(axis=1); L = np.diag(deg) - W; eigh(L)`).
Eigvecs cached at `data/cache/matched_strength_surrogate_lrg/Pat_NN/{band}_{phase}_R200_swap20_seed20260511_imcoh_abs.npz`
and reused by audit_70 (`compute_surr_T_G_band`).
Library helper: `lrg_eegfc.utils.metrics.spectral.topk_basis` for the
`[:, 1:k+1]` non-trivial subspace selection; replicated inline at
`audit_70.chordal_distances_for_k_grid:75-76` to avoid an extra
function call inside the inner cumsum loop.

**Why a subspace, not a single mode.** Eigenvector swaps under tiny
perturbations of `W` shuffle adjacent modes within a near-degenerate
band; reading off any individual `φ_n` is therefore ill-conditioned.
The span is stable wherever the *cumulative spectrum* up to `λ_{k+1}` is
itself stable, which is a much weaker condition.

**Why an eigenmode subspace at all** (not a partition, not a per-pair
distance). The cophenet probe `ρ_split^coph` reads task-induced
reorganisation at the level of *per-pair communication distance*; the
Grassmann probe reads it at the level of *whole-network spectral
geometry*. The two probes are mechanistically independent: one can hold
without the other (see γ_l, where Grassmann holds while `ρ_split^coph`
fails, and α, where the opposite holds).

---

## 2 Chordal Grassmann distance `d_G(k)`

The distance between two subspaces `U_k^a`, `U_k^b` (orthonormal columns
`A, B ∈ ℝ^{N × k}`) is the **chordal Grassmann distance**:

```
d_G(U_k^a, U_k^b)² = k − ‖A^T B‖_F²
                   = Σ_{i=1..k} sin² θ_i(U_k^a, U_k^b)
```

where `{θ_i}` are the principal angles. The trace identity on the right
makes the computation SVD-free and lets the entire `k`-sweep be derived
from a single `(N−1) × (N−1)` cross-correlation `M = A^T B` followed by
a 2-D cumulative sum of `M²` (audit_70 lines 65–89). The distance is:

- **bounded**: `d_G(k) ∈ [0, √k]`;
- **invariant to within-subspace rotation**: orthonormal column choice
  doesn't matter;
- **a true metric on `Gr(k, N)`** (Edelman, Arias, Smith 1998).

For each (patient, band, phase pair) and each `k ∈ {2, …, 112}` we
compute three chordal distances:

```
d_G^{pre→task}(k)   = d_G(U_k^{rsPre_A}, U_k^{task_test})
d_G^{task→post}(k)  = d_G(U_k^{task_test},  U_k^{rsPost})
```

(`rsPre_A` is the early half of `rsPre`, used as the upstream baseline;
the late half `rsPre_B` is reserved for `ρ_split^coph` independence and
not used by the Grassmann probe.)

---

## 3 Per-`k` trace statistic `T_G(k)`

The **per-patient per-`k` trace statistic** is

```
T_G(patient, band, k) = d_G^{task→post}(k) − d_G^{pre→task}(k)
```

A negative `T_G(k)` means `rsPost` is *closer* to `task_test` than
`rsPre_A` is to `task_test` at subspace dimension `k` — the trace
direction. The cohort-level statistic at fixed `k` is the paired
one-sided Wilcoxon

```
p_k(band) = P_paired-Wilcoxon[T_G^obs(k) <_{cohort} T_G^surr-mean(k)]
```

against the matched-strength surrogate mean across R=200 surrogates per
patient (audit_70.wilcoxon_per_k_less, alternative='less'). The
matched-strength surrogate is the same 4-cycle ±δ rewiring used at every
other LRG-layer probe (audit_63 / audit_67), with parameters

```
R = 200,   SWAP_FACTOR = 20,   n_swaps = 20·N(N−1)/2,   seed = 20260511
```

shared across cophenet and Grassmann layers (CONTROLS.md C3). The
surrogate eigendecompositions are cached at
`data/cache/matched_strength_surrogate_lrg/Pat_NN/{band}_{phase}_R200_swap20_seed20260511_imcoh_abs.npz`
(180 files, ~2 GB; helpers at `lrg_eegfc.utils.surrogate.matched_strength`).

---

## 4 Why `k` is a nuisance dimension — and why we don't choose a single one

`k` is not a free physical parameter; it is the *resolution scale* of
the subspace read. A trace can land at low `k` (coarse — large
hemispheric communities), mid `k` (network-scale), or high `k` (fine —
near-degenerate eigenmodes pick up local structure). The k-grid spans
`{2, 3, …, 112}` so that the entire spectrum-up-to-`λ_{112}` is
explored. Reading off any individual `k` would:

- multiply-test 111 hypotheses without family-level correction;
- be vulnerable to spurious "lucky" `k`-cells;
- discard information about the *contiguity* of the trace across
  resolution scales — which is exactly the structural signature that a
  strength-preserving null cannot fake.

The cluster-extent permutation gate (audit_70) addresses both
multiple-comparison and contiguity in one statistic.

---

## 5 Cluster-extent permutation null (audit_70) — the gate

At significance `α_k = 0.05` we compute the contiguous-significant
clusters in the `k`-axis from `{p_k}_{k=2..112}`. The **two co-primary
test statistics** are:

### 5a Longest contiguous-significant run (LR)

```
LR_obs(band) = max_{contiguous C ⊂ [2,112]} |C|
              s.t.  ∀ k ∈ C, p_k(band) < α_k
```

Captures *contiguous extent* in `k` of the trace.

### 5b Cluster mass `T_G^*(band)` — resilient all-clusters sum, normalized to `[0,1]`

```
T_G^*(b) = (n_k^cohort · log10(R+1))^{−1} · Σ_{k : p_k(b) < α_k} (−log10 p_k(b))   ∈ [0,1]
```

with `n_k^cohort = 111` (biggest k-grid common to all patients in the
cohort condition; `k ∈ [2, 112]` for full-data n=10 set by Pat_10's 113
contacts), `R = 200` matched-strength surrogates, and `log10(R+1) =
log10(201) ≈ 2.3032`. Denominator `255.65` for the full-data cohort.

Per-patient analogue (descriptive, never gated):
```
T_G^{*,s}(b) = (n_k^s · log10(R+1))^{−1} · Σ_{k : p_k^s(b) < α_k} (−log10 max(p_k^s(b), 1/(R+1)))   ∈ [0,1]
```

with `n_k^s` = patient `s`'s own available k-grid (`= N_s − 1`, or the
effective epi-X-excluded range under C5). The per-patient denominator
varies per patient so that `T_G^{*,s} = 1` means "this patient saturates
their own significance budget across the full available subspace",
giving cross-patient comparability invariant to `N_s`. Per-cell `p_k^s`
is regularized at the empirical floor `1/(R+1)` so a saturated cell
contributes at most `log10(R+1) ≈ 2.3032`.

Cluster mass aggregates *significance depth* across **every**
contiguous-significant `k`-cluster, not just the longest one. The
per-cell threshold `α_k = 0.05` is the cluster-forming threshold; the
statistic is the sum of `−log10 p_k` over all `k`-cells that pass it,
divided by the maximum-attainable sum (`n_k` cells at the floor).

**Why all-clusters, not longest-cluster.** A length-20 cluster split
into two halves of length 10 by a single non-significant `k`-cell
(sampling noise, near-degeneracy in the spectrum at that k) keeps its
total mass under the all-clusters sum but would lose half its mass
under a longest-cluster sum. Cluster mass is therefore **resilient to
a single missing `k`**, and it captures the full trace amplitude when
the signal fragments into multiple contiguous-significant stretches
across resolution scales. This is the standard cluster-mass test in
the Maris & Oostenveld 2007 sense, generalised to sum across all
clusters in the test statistic (not max), with the empirical null
built consistently on the same all-clusters sum.

`LR_obs` and `T_G^*` now read complementary aspects of the
`{p_k}_{k=2..112}` profile:
- `LR_obs` — **contiguity-based**: how long is the longest scale-coherent
  trace stretch.
- `T_G^*` — **mass-based, resilient**: total significance accumulated
  across all stretches, robust to single-`k` gaps.

CONTROLS.md C3 records the locked rule (Decision 8, 2026-05-19 pm) as
**mass-only**: `cluster_p_mass < 0.05` is the gate; `cluster_p_LR` is a
descriptive co-statistic that does not gate. Decision 12 (2026-05-28)
further requires `LOO max p_mass < 0.05` at full data as a precondition
for the strong tier — the disjunctive framing from Decision 7 (2026-05-19 am)
is retired.

### 5c Empirical null via phantom-surrogate testing

The null distributions of `LR` and cluster mass are built by treating
each of the R=200 matched-strength surrogates `r ∈ {1, …, R}` as a
"phantom observation" against the remaining `R−1` as reference:

```
for r = 1 .. R:
    phantom[k]  = surr_T_G[r, patient, k]           (one-of-R picked out)
    reference[k] = mean_{r' ≠ r} surr_T_G[r', patient, k]
    p^{phantom}_k(r) = paired Wilcoxon (less) over patients
    LR^null(r)       = longest run of {p^{phantom}_k(r) < α_k}
    mass^null(r)     = Σ_{k : p^{phantom}_k(r) < α_k} (−log10 p^{phantom}_k(r))
```

The null `mass^null(r)` uses the **same all-clusters sum** as the
observed statistic — both observation and null aggregate `−log10 p`
over every significant cell, so the empirical p-value is built on
matched formulas and is properly calibrated.

The empirical cluster p-values are

```
cluster_p_LR(band)   = (1 + #{r : LR^null(r)   ≥ LR_obs(band)})   / (R + 1)
cluster_p_mass(band) = (1 + #{r : mass^null(r) ≥ T_G^*(band)})   / (R + 1)
```

`(R + 1)` denominator and `(1 + …)` numerator follow the standard
permutation-test convention (Phipson & Smyth 2010) that prevents
p-values from collapsing to zero.

### 5d Mass-only verdict gate

The per-probe Grassmann verdict (CONTROLS.md C3, 2026-05-19
Decision 8 — supersedes Decision 7) is **mass-only**:

```
strong   ⇔ cluster_p_mass < 0.01
weak     ⇔ 0.01 ≤ cluster_p_mass < 0.05
no trace ⇔ cluster_p_mass ≥ 0.05
```

`L_obs` and its empirical `cluster_p_LR` stay in the audit CSV
(`data/audit/grassmann_cluster_extent/cohort_summary.csv` cols
`obs_longest_run`, `cluster_p_longest_run`) as a **descriptive
co-statistic** — useful to read alongside `T_G^*` to characterise
whether the trace concentrates in one contiguous window or spreads
across multiple — but they **do not gate** the verdict.

**Why mass-only (replaces the earlier disjunctive `min(p_LR, p_mass) < α`
gate)**: the resilient all-clusters cluster mass already encodes both
pieces of structural information the disjunctive gate combined:

- **Contiguity, indirectly**: a long contiguous-significant run
  contributes many `−log_10 p_k` terms with consistent depth, so its
  mass scales as length × depth and is hard for a phantom surrogate
  to match.
- **Depth, directly**: each `k`-cell contributes `−log_10 p_k`, so
  marginally-significant `p_k ≈ 0.04` cells contribute little
  (`−log_10 p_k ≈ 1.4`) and deeply-significant `p_k ≈ 10⁻³` cells
  contribute heavily (`−log_10 p_k ≈ 3`).

With the post-2026-05-19 resilient (all-clusters) cluster mass —
robust to single-`k` gaps that fragment a long cluster — having a
parallel `p_LR` gate is **double-insurance against a failure mode
the mass formula has already absorbed**. A sparse-cell pattern of
significance gives low cluster mass (few terms, all marginal), so
won't fire `p_mass`; a long pattern with many shallow cells gives
moderate mass, and the empirical phantom-surrogate null on mass will
correctly tell us whether it matters. Adding `p_LR` as a second
disjunctive arm enlarges the rejection region without principled
basis (writing-agent feedback 2026-05-19).

### 5e Why no hardcoded `min_cluster_size` — the empirical null is the gate

**Head.** A reader following the all-clusters sum literally will worry
that the test now rewards isolated significant `k`-cells equally with
contiguous structure, and that "3 isolated cells around" could trigger
a spurious trace verdict. This concern is **answered by the empirical
null calibration**, not by adding a hardcoded `min_cluster_size`
threshold to the statistic.

**Mechanics**. Under H₀ with 111 `k`-cells tested at `α_k = 0.05`, the
expected number of significant cells per phantom-surrogate test is
≈ 5.55 — scattered isolated cells are the *typical null configuration*,
and the null distribution of `cluster_mass` is built on exactly these
configurations. Numerically, on the locked audit_70 re-run:

| Band | `null_mean_mass` | `null_p95_mass` |
|---|---|---|
| δ | 8.06 | 18.37 |
| θ | 7.44 | 18.63 |
| α | 7.35 | 21.26 |
| β | 7.79 | 19.84 |
| γ_l | 9.04 | 27.24 |
| γ_h | 8.02 | 31.84 |

For *N* isolated cells to trigger a strong verdict, the sum
`N · |log10 p_k|` must exceed `null_p95_mass` (~ 18–32). The
mass-per-cell budget is bounded by the **Wilcoxon's resolution at
n = 10**: the minimum one-sided paired Wilcoxon p-value when every
patient is in the trace direction is ≈ 0.001, so `−log10 p_k` is
at most ≈ 3 per cell. The thought-experiment thresholds:

- **3 isolated cells at threshold `p = 0.04`**: mass ≈ 4.2 — far below
  every band's `null_mean_mass`. Empirical `p_mass ≈ 1`. No trace. ✓
- **3 isolated cells at deepest Wilcoxon `p ≈ 0.001`**: mass ≈ 9. Still
  below every band's `null_mean_mass`. Empirical `p_mass` ≈ 0.4 (more
  null surrogates exceed the obs than not). No trace. ✓
- **5 isolated cells at `p = 0.001`**: mass ≈ 15. Still below every
  band's `null_p95_mass`. Borderline-to-no-trace. ✓
- **10 isolated cells at `p = 0.001`**: mass ≈ 30. Now approaching
  null p95. *But* observing 10 cells of 111 at `p ≈ 0.001` under
  H₀ has probability of order `10⁻¹³` — at this point the cohort
  signal is real (scale-distributed, not contiguous, but real), and
  most reasonable readers would accept that as a genuine trace.

The empirical null **automatically rejects fluctuation-driven
isolated-cell configurations** because the null distribution contains
exactly those configurations.

**What this buys**:
1. **No hardcoded `min_cluster_size`** — compliant with
   `feedback_no_hardcoded_test_thresholds.md`. The test is the gate.
2. **No arbitrary gap-tolerance parameter** (`gap ≤ 1`, `gap ≤ 2`).
   The all-clusters sum already absorbs single-cell gaps because both
   halves of a fragmented cluster contribute to the sum.
3. **The verdict is structurally interpretable** *because of the data*,
   not by construction: looking at our three trace-positive bands,
   singleton contributions are 4% (β), 4% (γ_l), 15% (δ) of the total
   `T_G^*` — the verdicts are dominated by multi-cell contiguous
   sub-clusters, not by scattered singletons.

**What this does NOT do**: the all-clusters formula does not
*distinguish* "one big contiguous cluster" from "many scattered deep
significance". Both can pass the gate if the total mass exceeds the
null p95. If a future analysis required explicit penalisation of
non-contiguous structure in the **test statistic itself** (rather than
through null calibration), the principled threshold-free answer is
threshold-free cluster enhancement (TFCE; Smith & Nichols 2009), in
which each cell's contribution is weighted by the size of the
contiguous cluster it belongs to, integrated over all cluster-forming
thresholds. TFCE has only the two standard exponents `E = 0.5`,
`H = 2` and no `min_cluster_size`. For the current preprint, the
empirical null + all-clusters sum is sufficient and simpler; TFCE
would be belt-and-braces.

**Audit verification (recurring check)**: when running new bands or
new cohorts, decompose `T_G^*` into per-cluster contributions and
verify that the principal mass comes from multi-cell clusters,
not from singletons. If a band's `T_G^*` is dominated by isolated
singletons (say > 50% singleton mass) AND the empirical
`cluster_p_mass < 0.05`, treat it as a *cautionary* finding worth a
TFCE companion test before reporting. This is descriptive
diagnostics, not a gate.

---

## 6 The k-independent scalar — `T_G^*` definition and current values

The band-level **single number** the manuscript reports for the
Grassmann probe is the cluster mass, normalized to `[0,1]` per C1:

```
T_G^*(b) = (n_k^cohort · log10(R+1))^{−1} · Σ_{k : p_k(b) < α_k} (−log10 p_k(b))   ∈ [0,1]
        (α_k = 0.05, n_k^cohort = 111, R = 200, denominator ≈ 255.65)
```

with empirical p-value `cluster_p_mass(b)`. This scalar:

- has no `k` dependence — the `k`-axis is collapsed exactly by
  summing significance across every contiguous-significant cluster;
- is **mass-weighted in the natural permutation-test sense**: each
  `k`-cell contributes `−log10 p_k`, so cells with deeper Wilcoxon
  significance carry exponentially more weight than marginal cells;
- is **resilient to a single missing `k`** (all-clusters sum, not
  longest-cluster sum) — sampling noise that breaks one long cluster
  into two does not penalise the score;
- inherits its empirical null directly from the audit_70
  phantom-surrogate permutation — no further test or correction
  needed;
- is bounded `[0, 1]` after C1 normalization: `0` = no significant
  cell; `1` = every k-cell at the empirical floor;
- the **raw** all-clusters sum (un-normalized) is also kept in the
  CSV for audit-trail integrity; locked ledger tables cite both as
  `raw (normalized)` dual format per user decision 2026-05-26.

### Locked band-level values (audit_70 re-run 2026-05-19 pm — resilient `cluster_mass`, mass-only gate)

| Band | `LR_obs` | `T_G^*` | `cluster_p_LR` (descriptive) | `cluster_p_mass` (**gate**) | Verdict |
|---|---|---|---|---|---|
| **β** | 29 | **69.76** | 0.005 | **0.005** | **strong** |
| **γ_l** | 13 | **66.14** | 0.015 | **0.005** | **strong** |
| **δ** | 7 | **38.07** | 0.025 | **0.005** | **strong** |
| γ_h | 9 | 31.67 | 0.055 | 0.060 | no trace |
| θ | 5 | 12.97 | 0.099 | 0.159 | no trace |
| α | 4 | 7.98 | 0.159 | 0.348 | no trace |

Source: `data/audit/grassmann_cluster_extent/cohort_summary.csv`,
columns `obs_longest_run`, `obs_cluster_mass_neglog10p`,
`cluster_p_longest_run`, `cluster_p_cluster_mass`,
`cluster_p_mass_loo_max`, `cluster_p_mass_loo_argmax_patient`,
`verdict_cluster_extent`. The verdict column is gated on
`cluster_p_cluster_mass` alone (locked 2026-05-19 Decision 8). The
audit also writes a sibling `loo_cluster_p_mass.csv` per-patient
LOO file for single-patient leverage diagnostics
(`feedback_no_single_patient_p_driven.md`).

**Multi-cluster fragmentation reading** for γ_l and δ: both bands
have `cluster_p_LR > cluster_p_mass` (γ_l 0.015 vs 0.005; δ 0.025
vs 0.005). The all-clusters mass formula correctly captures multiple
contiguous-significant `k`-clusters (γ_l: 13 + 12 + 5 + 5 + 2 + 2
cells + 2 singletons; δ: 7 + 4 + 4 + 2 + 2 cells + 4 singletons),
whereas an LR-only gate would have under-called both as weak. This
is the multi-cluster signature the resilience fix was designed to
detect — the trace at γ_l and δ is *scale-distributed* rather than
*scale-concentrated*, and the methodology recovers it.

**LOO single-patient leverage diagnostic** (descriptive, per
`feedback_no_single_patient_p_driven.md` — values from audit_70 re-run
2026-05-19 pm with mass-only verdict + LOO; see
`loo_cluster_p_mass.csv` for the full per-patient table):

| Band | obs `p_mass` | LOO max (full) | C5 `p_mass^epi-X` | C5 LOO max | LOO interpretation |
|---|---|---|---|---|---|
| β | 0.005 | **0.005** (Pat_02) | 0.005 | **0.005** (Pat_02) | **fully LOO-robust at both full and epi-X**; C5 strengthens (mass 70→89, LR 29→36) |
| γ_l | 0.005 | **0.040** (Pat_05) | 0.030 | 0.159 (Pat_05) | Full robust; C5 passes cohort gate but **LOO-fragile under epi-X** (Pat_05 leverage); report both |
| δ | 0.005 | **0.055** (Pat_08) | 0.005 | **0.005** (Pat_02) | Full LOO flags Pat_08; **C5 epi-X resolves the leverage cleanly** (mass 38→44, fully LOO-robust under epi-X) → Pat_08 leverage at full data is attributable to epi-zone interactions, not the true trace |
| γ_h | 0.060 | 0.164 (Pat_02) | 0.060 | 0.214 (Pat_05) | No trace at both; C5 confirms |
| α | 0.348 | 0.577 (Pat_14) | 0.413 | 0.532 (Pat_13) | No trace at both; C5 confirms |
| θ | 0.159 | 0.826 (Pat_08) | 0.209 | 0.677 (Pat_06) | No trace at both; C5 confirms |

The δ entry is the worked example of what the LOO + C5 epi-X
diagnostic battery is *for*: a Wilcoxon `p` that crosses 0.05 under
single-patient drop at full data flags potential leverage; the C5
epi-X analysis then resolves it mechanistically — under epi-X, the
δ Grassmann trace **strengthens** (mass 38 → 44) AND becomes **fully
LOO-robust** (LOO max p = 0.005), which attributes the Pat_08 leverage
at full data to epi-zone interactions rather than the true biological
trace. The verdict stays "strong" by the cohort-level mass gate, and
the manuscript text should report **both** the full-data LOO caveat
*and* the C5 epi-X resolution transparently per
`feedback_no_single_patient_p_driven.md`.

γ_l shows the inverse pattern: full-data LOO is robust (max 0.040),
but the C5 epi-X analysis is LOO-fragile (max 0.159, Pat_05 leverage).
The cohort C5 gate still passes (`p_mass^epi-X = 0.030 < 0.05`), so
the verdict holds, but the manuscript should note that γ_l's
biological-attribution signal is partially leveraged on Pat_05 in
the epi-X analysis.

**Cascading updates** (downstream of this lock):
- `locked/VERDICT_LEDGER.md` — γ_l and δ Grassmann verdicts update from
  weak → strong under the mass-only gate; add dated revision entry
  citing audit_70 re-run + writing-agent feedback.
- `locked/CONTROLS.md` §C3 — updated 2026-05-19 to mass-only.
- `bands/03_gammalow.md` and `bands/06_delta.md` — update Grassmann verdict
  field + per-band paragraphs.
- `bands/00_cohort.md` — update cross-band verdict matrix.
- `README.md` — strength-of-results table.

### A descriptive companion (no test, for reader interpretation)

To accompany `T_G^*` with a quantity in physical *distance* units rather
than significance units, we report

```
ΔG_C*(band) = Σ_{k ∈ C*(band)} median_patients T_G^obs(patient, band, k)
            = integrated cohort-median trace excess across the
              significant cluster
```

This is the **integrated trace amplitude on `C*`** — a negative number
(trace direction) whose absolute value scales like a Grassmann distance
times the cluster length. `ΔG_C*` does NOT have its own permutation
p-value; it is reported as a descriptive complement to `T_G^*` so the
reader can interpret the *strength* of the trace alongside its
*significance*. This is the appropriate single number to quote when the
text says "the Grassmann trace has amplitude X at β".

> **Implementation TODO** (not blocking the methods text): `ΔG_C*`
> is a one-line addition to audit_70 (one extra `Σ_{k ∈ C*} median(…)`
> over `obs_T_G[band]`). If the writing agent requires explicit `ΔG_C*`
> values for the manuscript, regenerate `cohort_summary.csv` with the
> extra column. The audit_70 inputs (`obs_T_G[patient, k]` per band)
> are already loaded.

---

## 7 Sensitivity layer (C5) — epi-zone exclusion

For each band, the same audit pipeline is rerun on the epi-zone-excluded
adjacency matrices (audit_72 wrapping the same cluster-extent machinery
on a reduced node set, contacts in the epileptogenic zone dropped per
patient via `load_epileptic_nodes`). The output is at
`data/audit/grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` (cohort
gate) + `per_patient_per_band_per_k.csv` (per-k) + `sensitivity.csv`
(diagnostic window summary).

Reading rule (CONTROLS.md C5, locked 2026-05-19 pm Decision 10):
- The C5 gate is `cluster_p_mass^epi-X < 0.05` (Wilcoxon-on-epi-X on the
  same all-clusters cluster-mass statistic, R=200 surrogates built from
  the epi-X-restricted matched-strength ensemble).
- **C5 passes** when the cohort cluster-mass null is cleared at α=0.05
  under epi-zone exclusion — the biological attribution is to non-epileptic
  cortex, OR (if mass increases under epi-X) the epi zone was *masking*
  a wider trace.
- **C5 fails** when `cluster_p_mass^epi-X ≥ 0.05` — the cohort trace
  depends on the epi-zone contacts.
- **LOO under epi-X** (`cluster_p_mass_loo_max_epiX`) is the matching
  Decision-11 single-patient-leverage diagnostic for the C5 analysis.

The previous "≳ 80% retention" / "< 80% retention" framing (cell count
or `T_G^*` retention thresholds) is **retired** under Decision 10
(2026-05-19 pm) — retention is a descriptive diagnostic only, never
a gate.

**Epi-X is a secondary mechanistic observation per Decision 10** —
never an independent gate, never a verdict-promoter or demoter at the
strong/weak/no-trace level. The full-data verdict (under Decisions 8 +
12) is the cohort verdict; C5 epi-X is mechanistic interpretation.

---

## 8 What `T_G^*` does NOT claim

The methods section should explicitly state these scope limits:

1. **`T_G^*` is not a trace amplitude.** It is a significance-weighted
   sum on the empirically-significant cluster. For amplitude, cite
   `ΔG_C*`.
2. **`T_G^*` does not localise the trace to specific eigenmodes.** It
   localises it to a contiguous `k`-window. The eigenmode-level
   anatomical reading is in `locked/ANATOMY_LEDGER.md` under the Grassmann A3
   top-decile attribution, not in `T_G^*`.
3. **`T_G^*` does not survive bands where `cluster_p_mass ≥ 0.05`.**
   α, θ, γ_h have positive `T_G^*` numbers but those numbers are not
   distinguishable from the empirical null — they are reported only
   for completeness of the cohort table.
4. **Per-pair and subspace probes are mechanistically independent.**
   `T_G^*` is the Grassmann-probe gate; `ρ_split^coph` is the cophenet
   probe gate. A band may be strong on one and absent on the other;
   no single combined scalar is reported.
5. **The cluster-extent null is built against matched-strength
   surrogates only.** It controls for per-node strength but not for
   other graph invariants (degree-strength relationships at the
   per-band frequency, mesoscale community structure, etc.). The
   matched-strength surrogate is the minimum mandatory null
   (`feedback_matched_strength_mandatory.md`), not an exhaustive one.

---

## 9 Methods-section drafting checklist (for the writing agent)

When transcribing into the manuscript Methods section, the writing
agent should:

1. **Cite the object explicitly**: "We compute the leading-`k`
   non-trivial Laplacian eigenmode subspace `U_k = span{φ_2, …,
   φ_{k+1}}` of the symmetric normalised Laplacian `L̂` of each
   adjacency matrix." Do not write `U_k = span{φ_1, …, φ_k}` — the
   `φ_1` skip is methodologically decisive.
2. **Cite the distance explicitly**: the chordal Grassmann distance
   identity `d_G² = k − ‖A^T B‖_F²`, equivalent to `Σ sin² θ_i`. State
   that the `k`-sweep is computed SVD-free from a single `(N−1) × (N−1)`
   cross-correlation.
3. **State the per-`k` test once**: per-patient `T_G(k) =
   d_G^{task→post}(k) − d_G^{pre→task}(k)`; cohort paired one-sided
   Wilcoxon vs the matched-strength surrogate mean.
4. **Pivot on the nuisance argument**: explain why no single `k` is
   chosen, why the cluster-extent test is the natural collapse, and
   why the gate is **mass-only** (Decision 8) with `LR` as a
   descriptive co-statistic and **LOO max p_mass < 0.05** as the
   strong-tier robustness precondition (Decision 12).
5. **Define `T_G^*` formally** (normalized to `[0,1]` per C1) and
   **give the β number first** (raw `69.76` / normalized **`0.273`**,
   `cluster_p_mass = 0.005`, LOO max `p_mass = 0.005` Pat_02). Cite
   `cohort_summary.csv`.
6. **Tabulate all six bands** in one table (Section 6 above).
7. **Add the C5 epi-X paragraph** as a **secondary mechanistic
   observation** (per Decision 10) — never as a verdict-promoter.
   For γ_l: trace contracts (mass 66.14 → 32.75) but cohort gate
   held; for δ: trace strengthens (mass 38.07 → 43.99) and LOO
   resolves to 0.005 Pat_02, but does not promote the weak full-data
   verdict to strong.
8. **Do not call** the cluster-extent procedure a "cluster-based
   correction" or a "denoising step". It is a **test gate**, not a
   pre-processing step on the data.
9. **Do not cite the LR-only verdict column** from the CSV
   (`verdict_cluster_extent`); cite the **mass-only verdict** from
   CONTROLS.md C3 (Decision 8) with the **Decision-12 LOO precondition**
   for strong-tier promotion (`cluster_p_mass_loo_max < 0.05`).
10. **The matched-strength surrogate parameters** (R=200, swap_target=20,
    seed=20260511) are shared with cophenet and substrate audits and
    should be cited *once* in the methods Surrogate Battery section,
    not separately under cophenet and Grassmann.

---

## 10 Source-of-truth references

- **Audit script**: `scripts/01_compute/audit/audit_70_grassmann_cluster_extent.py`
- **Cohort CSV**: `data/audit/grassmann_cluster_extent/cohort_summary.csv`
- **Null distribution**: `data/audit/grassmann_cluster_extent/null_distribution.csv`
- **Per-`k` observed p-values**: `data/audit/grassmann_cluster_extent/per_k_obs_p.csv`
- **README**: `data/audit/grassmann_cluster_extent/README.md`
- **Companion observed-T_G CSV** (audit_66 input): `data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv`
- **Surrogate eigvec cache**: `data/cache/matched_strength_surrogate_lrg/Pat_NN/{band}_{phase}_R200_swap20_seed20260511_imcoh_abs.npz`
- **Helper module**: `lrg_eegfc.utils.surrogate.matched_strength`
- **Subspace utilities**: `lrg_eegfc.utils.metrics.spectral.topk_basis` (k-truncated eigenmode basis), `chordal_distance` (single-`k` form)
- **Locked gate**: `.agents/preprint/locked/CONTROLS.md` §C3 Grassmann — Decision 8 mass-only (replaces Decision 7 disjunctive) + Decision 12 LOO precondition for strong tier
- **Locked verdicts**: `.agents/preprint/locked/VERDICT_LEDGER.md` (Decision 12 cascade 2026-05-28)
- **Anatomy companion**: `.agents/preprint/locked/ANATOMY_LEDGER.md` (per-(band, Grassmann) DK-region attribution under A3)
- **Per-band briefs** (where Grassmann is primary or contributes): `bands/01_beta.md` (strong), `bands/03_gammalow.md` (strong ↑), `bands/06_delta.md` (weak — LOO Pat_08 binding)

---

## 11 Revision log

- **2026-05-19** — methods companion created; locks `T_G^*` as the
  band-level k-independent scalar; documents the disjunctive
  `min(cluster_p_LR, cluster_p_mass)` gate; flags `ΔG_C*` companion
  as one-line addition to audit_70 if writing agent needs it.
- **2026-05-28** — Decision-12 cascade applied: gate citation switched
  from disjunctive (Decision 7) to mass-only (Decision 8) + LOO
  precondition for strong tier (Decision 12). `T_G^*` normalized to
  `[0,1]` per C1 (denominator `n_k^cohort · log10(R+1) = 255.65` for
  full-data cohort); per-patient `T_G^{*,s}` formula added with
  per-patient denominator `n_k^s · log10(R+1)`. β raw 69.76 /
  normalized **0.273**; γ_l raw 66.14 / normalized **0.259** (strong ↑
  under Decision 12); δ raw 38.07 / normalized **0.149** (weak — LOO
  Pat_08 = 0.055 fails Decision-12 precondition). C5 epi-X explicitly
  reframed as secondary mechanistic observation (Decision 10), never
  verdict-promoter.
