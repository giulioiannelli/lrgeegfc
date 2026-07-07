---
name: eigenmode-embedding-motion
type: scope
era: IMCOH_ABS × COHORT_N10
status: current
date: 2026-05-08
created: 2026-05-08
updated: 2026-05-08
scope: eigenmode-embedding-motion
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-29_e1-spectral-subspace-alignment.md
  - .agents/guides/task-persistence-investigation/2026-05-08_epi-eigenmode-localization.md
  - .agents/reports/archive/2026-05/2026-05-05_result-2-lrg-beta-trace.md
  - .agents/guides/02_methods/lrg-framework-guide.md
---

# Eigenmode-embedding motion across phases (Procrustes-aligned per-contact shift)

## Renormalization head

**The leading-k Laplacian eigenvectors of the |ImCoh| graph place each
contact at a point Y_i(k) ∈ R^k; cross-phase motion ‖Y_i^φ - Y_i^ψ‖_2
after orthogonal-Procrustes alignment of the eigenbases tells us *which
contacts move and how far* through the LRG-relevant low-dimensional
embedding. The Pat_02 β first-pass at k=3 shows the motion is sharply
localised: three contacts on probe Q (Q6, Q2, Q3) account for virtually
all the embedding shift, and Q6/Q2 return to their rest_pre positions in
rest_post (Δ_pre→tt ≈ 0.65–0.86, Δ_pre→post ≈ 0.02–0.05) — i.e. a
TRACE-shaped, contact-localised reorganisation visible only in the
embedding (it is invisible to global Grassmann d_chord, which compresses
all directions in one scalar). This is the per-contact decomposition of
the §5.4 Grassmann chordal distance.**

## Notation

(Inherits LRG primitives from `2026-04-29_e1-spectral-subspace-alignment.md`
and `lrg-framework-guide.md`.)

- `V_p` — node set for patient `p`, size `N_p`. `i ∈ V_p` indexes a
  contact in the canonical FC ordering (matches `load_channel_labels`).
- `Â_φ^p` — symmetric weighted FC matrix (|ImCoh|, default
  `fc_method='imcoh_abs'`) at phase φ ∈ {pre, tt, post}, band b. We
  drop the `^p, b` superscripts where unambiguous.
- `D̂_φ` — diagonal of row sums of `Â_φ`.
- `L̂_φ = D̂_φ − Â_φ` — symmetric weighted Laplacian.
- Eigendecomposition `L̂_φ v_k^φ = λ_k^φ v_k^φ`, ordered ascending,
  `0 = λ_1 < λ_2 ≤ … ≤ λ_N`. The trivial mode `v_1` (constant) is
  always discarded.
- `Y^φ(k) ∈ R^{N × k}` — the **eigenmode embedding at order k**: column
  `j` of `Y^φ(k)` is `v_{j+1}^φ` (j = 1..k).
- `Y_i^φ(k) ∈ R^k` — row `i` of `Y^φ(k)`, the embedding point of
  contact `i` at phase φ.

## Definitions

### Eigenmode embedding

```
Y^φ(k) := [v_2^φ | v_3^φ | … | v_{k+1}^φ]    ∈ R^{N × k}
Y_i^φ(k) := (v_2^φ(i), …, v_{k+1}^φ(i))      ∈ R^k
```

Each contact `i` has unit-magnitude variance across modes only on
average; `‖Y_i^φ‖_2` may vary widely with node strength. **Magnitude
is not normalised** — that's deliberate: a hub contact is *informationally*
larger and should sit further from the origin in the embedding, exactly
as in standard spectral-graph embeddings (Belkin–Niyogi 2003;
Coifman–Lafon 2006).

### Orthogonal-Procrustes alignment

Eigenvectors are defined up to sign; degenerate eigenvalues introduce
free rotations within the degenerate subspace. To compare embeddings
between phases, anchor on rest_pre and rotate the other phases:

```
R^{pre→tt}    := argmin_{R ∈ O(k)}  ‖Y^pre R − Y^tt‖_F
R^{pre→post}  := argmin_{R ∈ O(k)}  ‖Y^pre R − Y^post‖_F
```

solved by `scipy.linalg.orthogonal_procrustes`.

We then form the **aligned embeddings in the rest_pre reference frame**:

```
Ỹ^tt    := Y^tt   R^{tt→pre}      = Y^tt   (R^{pre→tt})^T
Ỹ^post  := Y^post R^{post→pre}    = Y^post (R^{pre→post})^T
```

Per-pair motion is computed *in the shared reference frame*:

```
Δ_i^{pre→tt}  (k)  := ‖Y_i^pre  − Ỹ_i^tt‖_2
Δ_i^{tt→post}(k)  := ‖Ỹ_i^tt   − Ỹ_i^post‖_2
Δ_i^{pre→post}(k) := ‖Y_i^pre  − Ỹ_i^post‖_2
```

The triangle inequality holds:
`Δ_i^{pre→post} ≤ Δ_i^{pre→tt} + Δ_i^{tt→post}` for every contact.
Equality holds iff `Y_i^pre, Ỹ_i^tt, Ỹ_i^post` are collinear in R^k.

> **Subtle alignment choice.** We do **not** Procrustes-align tt→post
> independently. Doing so would put `Ỹ^tt` and `Ỹ^post` in different
> reference frames and `Δ_i^{tt→post}` would no longer be coherent with
> `Δ_i^{pre→tt}` and `Δ_i^{pre→post}`. The chained alignment via the
> rest_pre anchor is the only choice that preserves the triangle
> inequality at the per-contact level.

### Per-contact trace score (analogue of Result 2 T_d)

By symmetry with the §5.3 / KC λ-blend trace formulation
`T_d = d(tt, post) − d(pre, tt)`:

```
T_i^Y(k) := Δ_i^{tt→post}(k) − Δ_i^{pre→tt}(k)
```

A trace-shaped contact has `T_i^Y < 0`: it moved during task and stayed
near the task position in rest_post. A reset-shaped contact has
`T_i^Y > 0`. An anchor-shaped contact has both small.

## Properties

- **Range.** `Δ_i ∈ [0, 2]` since each row of `Y^φ(k)` has
  `‖Y_i‖_2 ≤ 1` (unit-norm columns; entries bounded by 1 in absolute
  value). Empirically the dominant mass at k=3 is `< 1`.
- **Sign invariance.** `Δ_i` is invariant under sign flip of any single
  eigenvector because Procrustes searches over the full orthogonal
  group `O(k)` (sign flips ∈ O(k)).
- **Permutation invariance.** Procrustes also handles permutation
  within degenerate subspaces — at non-degenerate spectra (verify
  `gap_min`) this is a no-op.
- **Triangle.** As above:
  `Δ_i^{pre→post} ≤ Δ_i^{pre→tt} + Δ_i^{tt→post}` (per contact, by
  Euclidean triangle inequality).
- **Connection to Grassmann chordal distance.** The cumulative scalar
  `d_chord(Y^pre, Y^tt) = √(k − Σ σ_j^2)` (where σ_j are singular
  values of `(Y^pre)^T Y^tt`) **is** the principal-angles aggregate
  over the same subspaces; the per-contact `Δ_i` decomposes that
  scalar across rows. This is the key complementary relation: §5.4
  Grassmann gives the *global* embedding rotation, this measure gives
  the *per-contact* contribution to that rotation.
- **What this measure does NOT detect.**
  - **Pure scaling within a single mode.** If `Y^tt = αY^pre` with
    `α ≠ 1` along one column, Procrustes restricts to `O(k)` and
    cannot absorb the scale; instead it would partially rotate. A
    pure-scaling distortion shows up as cohort-wide motion, not
    localised — degree-correction nulls would catch this.
  - **High-order embedding effects.** k=3 is sensitive to coarse
    block structure; finer features (e.g. probe-level micro-clusters)
    appear only at higher k. Multiscale extension k ∈ {2, 5, 10, 20}
    is the natural follow-up.
  - **Phase-wise reordering of degenerate modes.** If λ_2 ≈ λ_3 in
    one phase but not in another, the column ordering can flip and
    the Procrustes alignment recovers it — but if the degeneracy is
    genuine (gap < 1e-6), the embedding *itself* is not unique and
    any motion measure is undefined for that contact.
- **Complexity.** Per (patient, band, k): three eigendecompositions
  (already cached) + two Procrustes solves (k×k SVD, O(k^3)) + N
  row-norms. Negligible: <100 ms per cell.

## Caveats & failure modes

| | Caveat | Mitigation |
|--|--|--|
| 1 | Procrustes anchored on rest_pre is asymmetric — choosing tt or post as anchor gives slightly different rotations. | Document the anchor; the cumulative quantity `d_chord(Y^a, Y^b)` is anchor-invariant; per-contact rows are not. The `T_i^Y` trace score is anchor-symmetric in the sense that it depends only on `Y_i^φ` displacements after a single self-consistent reference is chosen. |
| 2 | Pat_03 1024 Hz outlier — eigenvalue spread ~2× compared to 2048 Hz patients; this affects which modes occupy the first three slots. | Flag in cohort summary; report Pat_03 separately if a peak/null shows up only there. |
| 3 | Degeneracy at small λ_k (rare in our weight-heterogeneous full-rank graphs but possible). | Include a `gap_min` column in the per-contact CSV; flag (patient, phase) cells with `gap_min < 1e-6`. |
| 4 | Same-probe contiguous contacts on a sEEG probe have correlated FC by construction (probe-bias residual). They tend to embed close to each other and may co-move — that's a *feature* of the geometry, not a measurement bias, but it complicates the "which contacts move" interpretation. | Report probe identity alongside `Δ_i`; cohort-aggregate motion *per probe* is the natural next slice. |
| 5 | Magnitude of `Δ_i` depends on k; reporting raw values without specifying k is meaningless. | The CSV always carries `k`; the figure header always specifies k. |
| 6 | Embedding axes have no intrinsic meaning (the eigenvectors are determined only up to the orthogonal group acting on each eigenspace). The 2D scatter axes are the *Procrustes-aligned* basis, which fixes them up to a sign/permutation but not a unique geometric embedding. | Reading rule: "look at relative motion, not absolute coordinates"; multiple 2D pairs (panels (i)/(ii)/(iii)) help see consistent vs scatter-driven displacement. |
| 7 | Pat_10 channel mask (113 channels vs 116-row implants) and Pat_03 nperseg (1024 Hz) are honoured by the cached LRG result; per-contact rows align with `load_channel_labels(pat)`. | Defensive label fallback in script (uses `iN` placeholders if length mismatch). |

## Pseudocode

```
input:  fc_method = 'imcoh_abs', cohort, band, k
output: per-(p, b, contact): Δ_i^{pre→tt}, Δ_i^{tt→post}, Δ_i^{pre→post}

for p in cohort:
    for phase in {pre, tt, post}:
        load LRG result for (p, b, phase, fc_method='imcoh_abs')
        Y^phase  = eigenvectors[:, 1:k+1]                # N x k

    R_tt   = argmin_{R in O(k)} || Y^pre R - Y^tt   ||_F     # Procrustes
    R_post = argmin_{R in O(k)} || Y^pre R - Y^post ||_F
    Y_tilde^tt   = Y^tt   R_tt^{-1}     # alignment INTO the rest_pre frame
    Y_tilde^post = Y^post R_post^{-1}

    for i in 1..N:
        Δ_i^{pre→tt}    = || Y_i^pre   - Y_tilde_i^tt   ||_2
        Δ_i^{tt→post}   = || Y_tilde_i^tt - Y_tilde_i^post ||_2
        Δ_i^{pre→post}  = || Y_i^pre   - Y_tilde_i^post ||_2

    emit row(p, band=b, contact=i, label=ch_labels[i],
              dist_rspre_taskt=Δ_i^{pre→tt}, ...)
```

## Visualization spec

**First-pass figure** (Pat_02 β at k=3) — single PDF, 5 panels:

- **Top row (3 panels)**: 2D scatters of the embedding for the three
  phases overlaid:
  - panel (i): `v_2` vs `v_3`
  - panel (ii): `v_2` vs `v_4`
  - panel (iii): `v_3` vs `v_4`
  - colour by phase (blue = rest_pre, red = task_test, green = rest_post)
  - marker size ∝ max(Δ_i^{pre→tt}, Δ_i^{tt→post}) — directs the eye to
    movers
  - one shared figure-level legend (top centre, no per-axis legend)
- **Bottom row, left/centre (2/3 width)**: bar plot of top-15 movers,
  ranked by max(Δ_i^{pre→tt}, Δ_i^{tt→post}), with paired
  red (`pre→tt`) and green (`tt→post`) bars and channel-label x-ticks.
- **Bottom row, right (1/3 width)**: cohort teaser — one pair of bars
  per patient (mean over contacts of `pre→tt` and `tt→post`).

**Reading rules.**
- A scattered point in the 2D scatter panels with **no nearby same-colour
  partner of the other phases** = a contact that moved.
- If `red` and `green` cluster together but `blue` is far → trace-shaped
  motion (changed during task, stayed there).
- If `red` is far from both `blue` and `green` (which cluster) → reset.
- A bar with `Δ_pre→tt ≈ Δ_tt→post` and `Δ_pre→post ≈ 0` (use the CSV
  for the third number) = trace.
- A patient with both `pre→tt` and `tt→post` bars near zero in the
  cohort teaser = global anchor (no embedding reorganisation).

**Multiscale extension.** k ∈ {2, 5, 10, 20} as a 4×N small-multiples
grid in the future; for now k=3 only.

## Connection to prior tools

| Prior tool | Relation |
|--|--|
| §5.4 Grassmann chordal distance (`audit_27`, `audit_37`, `audit_46`) | **Direct cumulative aggregator.** `d_chord(Y^pre, Y^tt) = √(k − Σ σ_j^2)` is the global subspace rotation; `Σ_i Δ_i^2 = ‖Y^pre − Ỹ^tt‖_F^2 = 2(k − tr(Σ))` so `‖Y^pre - Ỹ^tt‖_F` is a 1-1 function of `d_chord` after Procrustes alignment. The per-contact `Δ_i` decomposes the Grassmann scalar across rows. |
| §5.3 per-pair correlation (`audit_55_ctm_per_pair_scatter`) | Independent decomposition of the *same global signal*. §5.3 collapses the embedding into a single Spearman ρ on pair distances `D_{ij}`; this measure decomposes it across rows of the embedding matrix. The two views are complementary: §5.3 sees pair-level reorganisation; this sees node-level reorganisation. |
| `2026-04-29_e1-spectral-subspace-alignment.md` (E1) | E1 is the cohort-level Δ-of-d_chord triangle. This is its per-contact decomposition. E1 detects subspace rotation; this localises it. If E1 is cohort-positive at (band, k), this measure should show the same direction at the cohort-mean of contacts. |
| `2026-05-08_epi-eigenmode-localization.md` (Direction C) | Both rest on the Laplacian eigenmodes. Direction C measures **mode-mass on a fixed clinical node set E_p**: a static, per-band, per-mode quantity. This measure is **per-contact embedding shift across phases**: a dynamic, per-band, per-contact quantity. They probe different facets (mode-mass vs node-coordinate-shift). |
| §5 LRG dendrogram-based probes (KC, RF, MSPC, Cohesion-CBR) | Hierarchical-tree based; this is embedding-based at fixed τ = 1/λ_max. Tree probes integrate over the full Laplacian; embedding-based at small k uses only the leading non-trivial spectrum. The two are not duplicative. |
| Result 2 (β LRG trace) | β was the loaded band in §5; the first-pass cell uses Pat_02 β. Predicted: contacts whose `T_i^Y < 0` should overlap with the audit_53/47 β trace-leaf set. **Cross-check pending.** |

What this measure **subsumes**: nothing in the current repo (per-contact
embedding-shift was not previously measured).

What this measure **complements**: §5.3 (pair-level reorganisation),
§5.4 Grassmann (global subspace rotation), Direction C
(mode-mass-on-E_p localisation).

What this measure **does not replace**: the cohort-level Grassmann
triangle test stays as the global rung; this is its row-decomposition.

## Implementation plan

### Status

- **Implemented** in `scripts/01_compute/audit/audit_59_eigenmode_embedding.py`
  (lapbrain conda env), n=10 cohort × 1 band (β) × k=3 first-pass run.
- **CSV:** `data/audit/eigenmode_embedding/per_contact_motion.csv`
  (1177 rows = 117+122+118+115+116+120+113+119+119+118 contacts).
- **Figure:** `data/reports/notes_verification_2026-05-08/figures/eigenmode_motion_pat02_beta.pdf`.

### Library reuse

- `lrg_eegfc.workflow.lrg.load_lrg_result` — cached eigenvectors (the
  E1-pivot step-0 cache extension landed already; eigenvectors of
  shape `(N, N)` are present for all `imcoh_abs` cells per
  `audit_26_eigenvector_cache_check.py`).
- `lrg_eegfc.utils.io.load_channel_labels` — per-patient channel labels
  (FC ordering); promoted to library 2026-05-08 (used here as the 10th
  caller).
- `scipy.linalg.orthogonal_procrustes` — direct dependency, no library
  helper needed.
- `lrg_eegfc.config.const.BRAIN_BAND_TEX_DICT` — LaTeX band labels.

No new helpers in `scripts/`. No new module promotion required for the
first-pass cell; if k-sweep + cohort-level statistics expand,
`embedding_motion(Y_a, Y_b, anchor)` would be the library entry point
under `lrg_eegfc.utils.metrics.spectral`.

### Multiscale extension (deferred)

- k ∈ {2, 5, 10, 20} as a sweep — same script, `--k` argument already
  accepts arbitrary k.
- Cohort-level statistics: per-contact cohort-Wilcoxon on `T_i^Y`,
  BH-FDR within band; report per-contact trace counts (analogue of
  Result 2 but at the contact granularity).
- Within-baseline null: split-half `rest_pre_A/B` analogue from the
  halves cache (already exists, see `audit_27` for the pattern).

### Universality / implant-bias-discount (the user's brief)

The user specifically asked: *"once discounted the implant bias reveal
universality across patients"*. The cohort teaser bar shows that mean
contact motion varies 5× across patients (Pat_07: 0.018; Pat_13: 0.116
for `pre→tt`). To discount implant bias we need to:

1. **Anatomy-stratified normalisation.** Group contacts by Desikan-Killany
   region per `implant_pat_NN.csv`, then compare *region-level mean
   motion* across patients. If the region-level pattern is preserved
   when patient identity is shuffled, that's universality.
2. **Probe-bias correction.** Same-probe contacts co-move. Report
   `Δ_i` separately for cross-probe vs same-probe neighbour pairs of
   contacts (the embedding's local geometry). Defer.
3. **Strength-stratified contact ranking.** High-strength contacts
   contribute more to the embedding by construction. Rank `Δ_i` within
   strength quintile (mirrors Direction C's degree-corrected null) so
   the comparison is across patient-matched contacts.

These are **deferred** from the first pass; the goal of the first pass
is to verify the math runs and produce a Pat_02 β image. Universality
analysis is a follow-up scope.

## Open questions

1. **Anchor choice for Procrustes.** rest_pre is the anatomical baseline
   and the natural anchor for our trace-direction tests
   (`d(tt, post) − d(pre, tt) < 0`). An alternative is "midpoint
   anchor" via generalised Procrustes (rotate all three to a common
   centroid). Defer; rest_pre anchor is the canonical default.
2. **k = 3 vs higher k.** Three is the smallest non-degenerate
   embedding for visualisation in 2D pairs; the spectrum has many more
   informative modes. The k ∈ {2, 5, 10, 20} sweep is the natural
   follow-up — the script already accepts `--k`.
3. **Cohort-positive criterion at the per-contact level.** A statistic
   like "fraction of contacts with `T_i^Y < 0` exceeds a within-baseline
   null" is the right cohort gate. Defer pre-registration.
4. **Probe-aware vs probe-blind aggregation.** Same-probe co-motion is
   real and biologically interesting (the probe samples a small tissue
   volume). Report both views; do not pre-select.
5. **Consistency with Result 2 β trace cohort.** The β trace cohort
   from §5 was 10/10 vs within-baseline null. Predicted overlap with
   contacts having `T_i^Y < 0` at β. Cross-check pending —
   `audit_47_kc_trace_network_view` β-trace-leaf set vs the
   per-contact `T_i^Y < 0` mask is the natural next probe.
6. **Effect of Procrustes orthogonal vs general linear.** Using
   `O(k)` rather than `GL(k)` is intentional (preserves Euclidean
   distance and the row-norm interpretation). A general-linear
   alignment would absorb the per-mode scaling and reduce `Δ_i` —
   but the resulting motion would no longer be a Euclidean
   distance in a fixed metric. Stay with `O(k)`.

## First-pass empirical numbers (Pat_02 β, k=3)

```
top-5 movers (Pat_02 β k=3):
channel  Δ_pre→tt  Δ_tt→post  Δ_pre→post
   Q6     0.866    0.819      0.054     ← TRACE-shape: large pre→tt and tt→post, ~0 pre→post
   Q2     0.643    0.655      0.017     ← TRACE-shape
   Q3     0.290    0.520      0.231
   H1     0.071    0.133      0.202
   H2     0.064    0.103      0.165
```

- **Q6 and Q2** — task moves them dramatically, then they return almost
  exactly to their rest_pre embedding position. This is the per-contact
  signature of a TRACE: not "the contact ends up somewhere new", but
  "the contact moves in task and snaps back, leaving a *transient*
  reorganisation". Whether this corresponds to the §5 β trace
  cohort-level finding is the natural cross-check.
- **Q3** — moved during task but did NOT snap back (Δ_pre→post = 0.231).
  Mixed pattern.
- **Probe localisation.** The top-3 movers (Q6, Q2, Q3) are all on
  probe Q. The next two (H1, H2) are on probe H. Probe-bias residuals
  in the embedding are real and the natural anatomy-stratified
  follow-up should treat probe identity as a covariate.
- **Cohort means** (mean contact motion):
  - `pre→tt`: Pat_07 (0.018) < Pat_02 (0.027) < … < Pat_08 (0.101) <
    Pat_13 (0.089).
  - `tt→post`: Pat_14 (0.026) < Pat_02 (0.038) < … < Pat_13 (0.116).
  - 5× spread across patients — implant-bias discount is required
    before universality claims.

These numbers are the first observation; sign and magnitude predictions
for the cohort universality test will be drafted after the implant-bias
discount step.
