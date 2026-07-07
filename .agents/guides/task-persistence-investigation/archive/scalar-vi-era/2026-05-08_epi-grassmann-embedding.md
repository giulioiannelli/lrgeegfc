---
name: epi-grassmann-embedding
type: scope
era: IMCOH_ABS × COHORT_N10
status: draft
created: 2026-05-08
updated: 2026-05-08
pointers:
  - .agents/plans/active/2026-05-08_lrg-epilepsy-research-directions.md
  - .agents/guides/task-persistence-investigation/2026-04-29_e1-spectral-subspace-alignment.md
  - .agents/guides/task-persistence-investigation/2026-05-08_epi-eigenmode-localization.md
  - .agents/reports/2026-05-08_direction-a-induced-subtree.md
  - .agents/reports/archive/2026-05/2026-05-08_trace-minus-epi-resection.md
  - .agents/guides/02_methods/lrg-framework-guide.md
---

# Epi Grassmann embedding — subspace alignment of epi nodes with the top-k eigenspace

## Renormalization head

**Treat the top-k Laplacian eigenspace `V_k = span(v_2…v_{k+1})` as a
point on the Grassmann manifold `Gr(k, N)` and ask one of two questions
about the epi node set `E_p`. (i) *Alignment.* How close is `V_k` to
the epi coordinate subspace `S_E = span{e_i : i ∈ E_p}` ⊆ ℝ^N? Measured
by the chordal Grassmann distance
`ε^E_align(p, b, φ, k) = sqrt(min(k, |E_p|) − ‖V_k[E_p, :]‖_F²)`.
Small ⇔ the global modes of the FC graph live preferentially on epi
nodes. (ii) *Resection.* How close is the *full-graph* eigenspace
restricted to non-epi rows to the eigenspace you'd compute from
scratch on the resected FC graph (epi rows/cols deleted)? Measured by
`δ^E_resect(p, b, φ, k) = d_chord(Q_full[N_p], V_k^resect)`. Small ⇔
removing epi nodes leaves the global mode structure intact; large ⇔
epi nodes carry information that doesn't exist in the resected graph.
Static, per-phase measure — no cross-phase aggregation needed for the
primary signature. Direction E is the cumulative-over-modes companion
to Direction C's per-mode `m^E_k`; the resection variant is genuinely
independent of C and of any prior tool. Eigenvectors already cached in
`IMCOH_LRG_CACHE`; no cache extension required.**

---

## Plain-language summary

**The eigenvectors of the FC Laplacian are the natural patterns of
coordinated activity that the connectivity supports.** Each pattern
("mode") assigns a value to every brain contact saying how much of
that pattern lives at that contact. The first few modes (the
"top-k") are the broadest, most coarse-grained ones. We have a
marked subset of contacts: the epileptic nodes `E_p`. Two questions:

**E.align — do the broad modes preferentially sit on epi contacts?**
For each top-k mode, sum the squared values at epi contacts; that's
a number in [0, 1] saying "what fraction of this mode lives on epi".
Average over the top-k modes → one summary number `f^E_k ∈ [0, 1]`
per (patient, band, k, phase). The Grassmann form of the same scalar
is `ε^E_align = sqrt(min(k, |E_p|) − k · f^E_k)`. **This is the
cumulative-over-modes form of Direction C's per-mode mass; same
content, different resolution.**

**E.resect — do the broad modes change structure if we remove the
epi contacts from the FC graph?** Build the FC graph on all contacts,
take its top-k modes; build the FC graph on only the non-epi contacts
(delete epi rows/cols of the FC matrix), take *its* top-k modes; for
the non-epi contacts, you now have two patterns. Compute the angle
between them. Small ⇔ removing epi barely changed anything; large ⇔
removing epi created a different mode set. **This is genuinely
independent of Direction C, Direction A, and audit_56.**

**Both variants are static (one number per phase). No cross-phase
trace structure is forced into the primary computation.** If the
static plots show something interesting per phase, a phase
comparison is a follow-up.

**Reporting.** Per (band, k, phase): plot the cohort distribution of
`ε^E_align` and `δ^E_resect` against the strength-stratified null —
median, IQR, per-patient values, null reference. **No pre-baked
acceptance gate.** We look at the picture together and decide what
counts as signal and what counts as noise.

---

## TOC

1. **Two scalars per (patient, band, phase, k).** `ε^E_align`
   (Grassmann distance from top-k eigenspace to epi coordinate
   subspace) and `δ^E_resect` (Grassmann distance between full-graph
   modes and resected-graph modes on the non-epi rows).
2. **`ε^E_align` ≡ cumulative form of Direction C's per-mode mass.**
   `f^E_k = (1/k) Σ_{j=2}^{k+1} m^E_j`. Same primitive, different
   resolution. Reported jointly because the Grassmann reading is
   user-flagged ("embedding in Grassmann space").
3. **`δ^E_resect` requires one fresh eigendecomposition per
   (p, b, φ).** ~5 s per patient. Genuinely new question.
4. **Strength-stratified null (Q=5 quintiles, R=100 draws).** Same
   null protocol as Direction C; protects against the trivial
   "high-strength nodes accumulate mode mass" confound.
5. **No acceptance gate.** Cohort distribution + null reference are
   the deliverable. Signal vs noise judged from the plot, post hoc.
6. **Cross-phase comparison deferred** until static plots are
   inspected. The epilepsy frame is not trace-locked
   (`feedback_epilepsy_not_trace_locked.md`).

---

## 1. Notation

(Inherits the LRG primitives from `2026-04-29_e1-spectral-subspace-alignment.md`
and `lrg-framework-guide.md`.)

- `V_p` — node set for patient `p`, size `N_p` (100 ≤ N_p ≤ 122 across
  the n=9 epi-annotated cohort).
- `E_p ⊆ V_p` — epi node set (channel-matched via
  `load_epileptic_nodes`). Same operational definition as Directions
  A and C. `N_p^E = |E_p|`. Pat_15 has no epi annotation → drops.
- `N_p^N = V_p \ E_p` — non-epi node set, size `N_p − |E_p|`.
- `Â_φ^p,b` — symmetric weighted FC matrix at phase
  `φ ∈ {pre, tt, post}`, band `b`. (Diagonal handling per the LRG
  pipeline.)
- `D̂_φ^p,b` = diag of row sums of `Â_φ^p,b`.
- `L̂_φ^p,b` = `D̂_φ^p,b − Â_φ^p,b` — symmetric weighted Laplacian.
  Spectrum `0 = λ_1 < λ_2 ≤ … ≤ λ_{N_p}`; eigenvectors `v_1, …, v_{N_p}`
  unit-norm; `v_1` is the trivial constant mode.
- `V_k^φ,p,b ∈ ℝ^{N_p × k}` — top-k non-trivial eigenvectors as
  column block: `[v_2, v_3, …, v_{k+1}]`. Orthonormal (V_k.T V_k = I_k).
- `S_E ⊆ ℝ^{N_p}` — *epi coordinate subspace*: `span{e_i : i ∈ E_p}`
  where `e_i` is the i-th standard basis vector. Dimension `|E_p|`.
- `Â_resect_φ^p,b ∈ ℝ^{N_p^N × N_p^N}` — FC submatrix with epi rows
  and columns deleted. Then `L̂_resect = D̂_resect − Â_resect` where
  `D̂_resect` is the diag of row sums of the resected adjacency. (Note:
  this differs from "delete rows/cols of the original L̂" because the
  diagonal degrees change.)
- `V_k^resect,φ,p,b ∈ ℝ^{N_p^N × k}` — top-k non-trivial eigenvectors
  of `L̂_resect_φ^p,b`.
- `Q[E]` — for any matrix `M`, `M[E, :]` is the row-restriction to
  rows indexed by set `E`.

## 2. Definitions

### 2.1 Grassmann chordal distance (recap)

Following `2026-04-29_e1-spectral-subspace-alignment.md` §2.1: given
two orthonormal bases `V_a ∈ ℝ^{D × p}` and `V_b ∈ ℝ^{D × q}` (in a
common ambient ℝ^D), let `σ_1 ≥ … ≥ σ_{min(p,q)} ≥ 0` be the singular
values of the cross-Gram `V_a.T V_b ∈ ℝ^{p × q}`. Each `σ_i = cos θ_i`
is the cosine of the i-th principal angle. The chordal distance is

```
d_chord(V_a, V_b) = sqrt( min(p, q) − Σ_i σ_i² )    ∈ [0, sqrt(min(p, q))]
```

`d_chord = 0` ⇔ one subspace is contained in the other (all `σ_i = 1`).
`d_chord = sqrt(min(p, q))` ⇔ the subspaces are mutually orthogonal.

For non-orthonormal column bases, orthonormalize via QR first.

### 2.2 Variant 1 — alignment with epi coordinates (E.align)

The epi coordinate subspace `S_E = span{e_i : i ∈ E_p}` has the
canonical orthonormal basis `E ∈ ℝ^{N_p × |E_p|}` with
`E[i, j] = 1` iff `i` is the `j`-th element of `E_p`, else 0.

Cross-Gram of `V_k^φ,p,b` with `E`:

```
V_k.T E = (V_k[E_p, :]).T   ∈ ℝ^{k × |E_p|}
```

Let `σ_E,1 ≥ … ≥ σ_E, m` (where `m = min(k, |E_p|)`) be its singular
values. Equivalently, `σ_E,i² = i-th eigenvalue of G_E := V_k[E_p, :].T V_k[E_p, :] ∈ ℝ^{k × k}`.

**Definition (E.align).**

```
ε^E_align(p, b, φ, k) = d_chord(V_k^φ,p,b, S_E)
                      = sqrt( min(k, |E_p|) − Σ_i σ_E,i² )
                      = sqrt( min(k, |E_p|) − ‖V_k^φ,p,b[E_p, :]‖_F² )
```

Range: `[0, sqrt(min(k, |E_p|))]`.

**Equivalent reading — epi fraction.** Define

```
f^E_k(p, b, φ) = ‖V_k[E_p, :]‖_F² / k   ∈ [0, 1]
```

By unit-norm of each `v_j` (`j = 2..k+1`):
`Σ_{i ∈ V_p} v_j(i)² = 1`, so
`Σ_{j=2}^{k+1} Σ_{i ∈ E_p} v_j(i)² = ‖V_k[E_p, :]‖_F²`. Therefore
`f^E_k = (1/k) Σ_{j=2}^{k+1} m^E_j` where `m^E_j = Σ_{i ∈ E_p} v_j(i)²`
is Direction C's per-mode projection mass.

**Identity:**

```
ε^E_align² = min(k, |E_p|) − k · f^E_k                  (when k ≤ |E_p|)
```

So `ε^E_align` and `f^E_k` are equivalent reparametrisations of the
same scalar — `ε^E_align` is the Grassmann form, `f^E_k` is the
fraction form. Both reported.

### 2.3 Variant 2 — resection stability (E.resect)

Compute the resected eigendecomposition `L̂_resect_φ^p,b v_k^resect = λ_k^resect v_k^resect`,
ordered ascending, drop the trivial mode, take the top k:
`V_k^resect,φ,p,b ∈ ℝ^{N_p^N × k}` (orthonormal).

Restrict the full eigenspace to non-epi rows:
`V_k^full[N_p^N, :] ∈ ℝ^{N_p^N × k}`. **This is in general not
orthonormal** (it is a row-block of an orthogonal matrix). Compute the
QR factorization `V_k^full[N_p^N, :] = Q_full · R` where
`Q_full ∈ ℝ^{N_p^N × r}` is orthonormal with
`r = rank(V_k^full[N_p^N, :]) ≤ k`. Generically `r = k` for our
weighted FC graphs (no exact eigenvector–coordinate alignments).

**Definition (E.resect).**

```
δ^E_resect(p, b, φ, k) = d_chord(Q_full, V_k^resect)
                       = sqrt( min(r, k) − Σ_i σ_R,i² )
```

where `σ_R,i` are singular values of `Q_full.T V_k^resect ∈ ℝ^{r × k}`.

Range: `[0, sqrt(min(r, k))]`. Generically `[0, sqrt(k)]`.

**Reading.** `δ^E_resect ≈ 0`: the modes you compute on the resected
graph are the same modes you'd see on the non-epi rows of the full
graph — epi nodes do not change the global mode structure on the rest.
`δ^E_resect ≈ sqrt(k)`: removing epi creates a fundamentally different
mode set; the epi nodes are structurally critical to the global
diffusion geometry.

### 2.4 Strength-stratified null

Both variants share the same null draw protocol (mirrors Direction C
§2.4):

For each patient `p`, partition `V_p` into `Q = 5` strength quintiles
on per-(phase, band) node strength `s_i = Σ_j Â_φ^p,b[i, j]`. Let
`q_p^i ∈ {0, …, Q−1}` denote node `i`'s quintile and `n_p^q` the
number of epi nodes in quintile `q`.

A null draw `S_r` is produced by sampling, in each quintile,
`n_p^q` nodes without replacement from
`{i ∈ V_p : q_p^i = q}` (or from non-epi nodes in that quintile if
the user opts for the stricter "exclude E_p from the pool"
variant — see §4 caveat 5).

For each null draw `S_r`, compute the variant-specific scalar:

- E.align null: `ε^{S_r}_align(p, b, φ, k) = sqrt(min(k, |S_r|) − ‖V_k[S_r, :]‖_F²)`
- E.resect null: re-eigendecompose `L̂_resect_S_r` (Laplacian of FC
  with rows/cols at `S_r` removed) and compute
  `δ^{S_r}_resect(p, b, φ, k)`.

`R = 100` null draws per (p, b, φ). Per-cell z-score:

```
z^E_align(p, b, φ, k) = (ε^E_align − μ̂(ε^{S_r}_align)) / σ̂(ε^{S_r}_align)
z^E_resect(p, b, φ, k) = (δ^E_resect − μ̂(δ^{S_r}_resect)) / σ̂(δ^{S_r}_resect)
```

**Direction reading (descriptive, not a test).** The "natural"
prediction is *small* `ε^E_align` (epi rows carry the global modes)
and *large* `δ^E_resect` (resection disrupts the modes). The data may
or may not show this. Both signs are reported per cell with no
pre-registered direction-test gate — what we actually find is what we
report.

### 2.5 What we report

No pre-registered acceptance gate. For each (band, k, phase) cell we
emit, per variant:

- **Per-patient observation:** `ε^E_align(p, b, k, φ)` and
  `δ^E_resect(p, b, k, φ)` (raw scalars, both reported).
- **Per-patient null reference:** mean and standard deviation of the
  R = 100 strength-stratified null draws for the same cell.
- **Per-patient z-score:** `(observed − null_mean) / null_std`.
  Reported alongside the raw scalar; not used as a gate.
- **Cohort summary:** median, Q1, Q3, min, max of the per-patient
  z-scores and of the per-patient raw observations across the n=9
  patients. Pat_03 (1024 Hz outlier) reported separately.
- **Footnote stats:** cohort-Wilcoxon two-sided p-value and BH-FDR q
  per band over m = 6 (k) × 3 (phase) = 18 cells, both variants.
  These are *reported*, not used to gate the lead description.

Decision about what counts as signal vs noise is made post hoc from
the figures and tables — not pre-baked into the script.

### 2.6 Phase-trace extension (E.trace, deferred)

Optional follow-up: if the static plots show something interesting,
compute the per-patient triangle

```
T_E_align(p, b, k) = z^E_align(post) − z^E_align(pre) − [z^E_align(tt) − z^E_align(pre)]
```

across phases (and the analogue for `δ^E_resect`), and inspect the
cohort distribution. Not part of the primary deliverable. Per
`feedback_epilepsy_not_trace_locked.md`, trace is one connectivity
pattern among many; not the default frame.

## 3. Properties

- **Range.**
  - `ε^E_align ∈ [0, sqrt(min(k, |E_p|))]`. Generically
    `sqrt(k)` upper bound when `k ≤ |E_p|`.
  - `δ^E_resect ∈ [0, sqrt(k)]` generically.
- **Invariance.**
  - Both are subspace-functionals: invariant under sign flips and
    intra-block rotations of the eigenbasis at degenerate eigenvalues
    (same as E1 §3 — this is the cleanest spectral primitive,
    sign-alignment-free).
  - Invariant under permutation of `E_p` (set-valued).
  - **Not** invariant under reweighting `Â → c · Â` for `c > 0` because
    the spectrum scales but the eigenvectors don't change — so
    eigenvectors *are* invariant, hence the measure is invariant under
    positive global rescaling of FC.
- **What E.align cannot detect.**
  - Localisation patterns split across small disjoint subsets of
    `E_p` (averaging over modes loses sub-cluster structure; addressed
    by Direction C per-mode resolution).
  - Higher-moment localisation (e.g. heavy-tailed amplitude
    distribution within a mode); addressed by `IPR^E_k` (Direction C
    secondary).
  - Phase information (the eigenvectors are real-valued; sign matters
    for individual entries but `m^E_k` uses the second moment).
- **What E.resect cannot detect.**
  - Mode rotations *within* the non-epi span (same `Q_full`, different
    rotation): δ is by construction subspace-invariant.
  - Eigenvalue-magnitude shifts that don't rotate the subspace (the
    scalar measures *direction*, not *speed*; complementary measure
    is the eigenvalue-shift `Δλ_k = λ_k^full − λ_k^resect`,
    deferred).
- **What both cannot detect** (shared with Direction C).
  - Effects driven by *which specific* modes carry the localisation.
    The measure aggregates the top-k subspace; it doesn't tell you
    "the 4th mode peaks on epi". Use Direction C `m^E_k` for that.
- **Sample size & multiple comparisons.** n=9 patients × 6 bands × 6
  k-values × 3 phases = 972 z-scores per variant. Per-band BH-FDR
  on m=18 (k × phase) cells gives the secondary table; primary lead
  is regularity-based (no FDR gate).
- **Complexity.**
  - E.align: per (p, b, φ): one SVD on a `k × |E_p|` matrix (k ≤ 21,
    |E_p| ~10-40). Negligible (<10 ms). With R = 100 null draws ×
    9 × 6 × 3 = 162 cells = 16 200 SVDs. Total < 30 s.
  - E.resect: per (p, b, φ): one eigendecomposition of an
    `(N_p − |E_p|) × (N_p − |E_p|)` matrix (~85 × 85), plus
    `R + 1` SVDs. Eigendecomposition ~50 ms; total ≈ 8 s × 162 cells
    ≈ 22 minutes. With nulls: ~36 minutes wall-clock if serial,
    or ~5 minutes if parallelised over patients.

## 4. Caveats & failure modes

| | Caveat | Mitigation |
|--|--|--|
| 1 | E.align is mathematically a *cumulative* form of Direction C's `m^E_k`. A regularity at low-k in E.align reflects the same per-mode structure that C resolves. | Be explicit in the report: don't double-claim. E.align answers "does the top-k subspace as a whole sit on epi"; C answers "at which modes does the mass concentrate". Direction E adds *resection* as the genuinely new question. |
| 2 | Eigenvalue degeneracy at the k vs k+1 cutoff makes the included eigenvector ambiguous. | k-grid `{2, 3, 5, 8, 13, 21}` (sparse) reduces collision probability; sensitivity check at k → k+1 should not flip the verdict. Filter cells where `min(λ_{k+2} − λ_{k+1}, λ_{k+1} − λ_k) < 1e-6 · λ_max`. |
| 3 | Heavy-tailed degree distribution gives high-strength nodes more eigenmode mass for trivial reasons. | Strength-stratified null (Q = 5 quintiles) is mandatory. Uniform-random null reported only as a sanity check showing the problem. |
| 4 | The resected Laplacian uses degrees recomputed on the resected adjacency, *not* the principal-submatrix of the full L̂. | Pseudocode is explicit: `L̂_resect = diag(rowsums(Â_resect)) − Â_resect`. Document in the script header. |
| 5 | Strength-stratified null with epi nodes included in the pool dilutes the effect (epi can be sampled into S_r). | Default: pool excludes E_p (canonical "non-epi controls" matched on strength). Sensitivity variant: pool includes E_p (looser; report alongside). |
| 6 | Pat_15 has no epi annotation. | Drops out (n = 9, same as Direction A and C). |
| 7 | Same-probe / cross-probe split (analogous to Direction A `E_p^cp`). | Optional cross-probe-restricted variant `E_p^cp` reported alongside `E_p`. |
| 8 | Pat_03 outlier (1024 Hz, larger eigenvalue spread). | Log Pat_03 separately; flag in cohort plot with distinct marker. |
| 9 | The QR step in E.resect for `V_k^full[N_p, :]` can rank-drop if `V_k^full` is anomalously aligned with epi coordinates (`r < k`). | Detect and report cells where `r < k`; treat as a *positive* alignment finding (epi rows are necessary to span the top-k subspace) rather than a numerical failure. |
| 10 | `ε^E_align` and `δ^E_resect` are not directly comparable (different scales, different upper bounds). | Don't aggregate them into a single score. Report side-by-side tables and per-band heatmaps; reader interprets each independently. |
| 11 | n=9 is small. The IQR floor is the right guardrail at this n. | Floor + cohort-direction count + bootstrap CI on the cohort median. Per-patient transparency in the spaghetti plot. |

## 5. Pseudocode

```
input:  IMCOH_LRG cache (per p, b, φ → eigenvalues, eigenvectors)
        Â_φ^p,b for each (p, b, φ) (FC matrix; needed for E.resect)
        E_p (epi node indices) per p
        K_GRID = [2, 3, 5, 8, 13, 21]
        R = 100 null draws
        Q = 5 strength quintiles
output: data/audit/epi_grassmann/E_per_patient.csv
        data/audit/epi_grassmann/E_cohort.csv
        figures: cohort heatmap, regularity scatter, resect-overlay,
                 (optional) phase-trace

# ---- helper (promoted to library) ----
function chordal_distance(A, B):
    # A: D × p (orthonormal), B: D × q (orthonormal)
    M = A.T @ B               # p × q
    sigma = singular_values(M)
    p_eff = min(rows(A.cols), cols(B))   # = min(p, q)
    return sqrt( p_eff − sum(sigma**2) )

function grassmann_chordal_to_subspace_E(V_k, E_idx):
    # V_k: N × k orthonormal; E_idx: subset of {0..N-1}
    block = V_k[E_idx, :]         # |E| × k
    fro2  = sum(block ** 2)
    m     = min(k, len(E_idx))
    return sqrt( m − fro2 ), fro2 / k       # eps_align, f_E_k

function chordal_full_vs_resect(V_k_full, V_k_resect, N_idx):
    # V_k_full: N × k; V_k_resect: |N_idx| × k
    Q_full, _ = qr(V_k_full[N_idx, :])      # |N_idx| × r, r ≤ k
    return chordal_distance(Q_full, V_k_resect), rank(Q_full)

# ---- per-(p, b, φ) ----
for each p in COHORT_N9_EPI:
    for each band b in BRAIN_BANDS_NAMES:
        for each phase φ in {pre, tt, post}:
            lrg = load_lrg_result(p, φ, b, 'imcoh_abs')
            V_full = lrg.eigenvectors          # N × N
            λ_full = lrg.eigenvalues
            A      = load_fc_matrix(p, φ, b, 'imcoh_abs')
            E_idx  = load_epileptic_nodes(p)
            N_idx  = [i for i in 0..N-1 if i not in E_idx]
            s      = A.sum(axis=1)
            quint  = quantile_bin(s, Q)

            # Resected eigendecomposition (E.resect only)
            A_R = A[N_idx, :][:, N_idx]
            L_R = diag(A_R.sum(axis=1)) − A_R
            λ_R, V_R = eig(L_R)  # ascending
            # drop trivial v1 (constant), take next k columns

            for each k in K_GRID:
                # ---- E.align observed ----
                eps_E,  f_E  = grassmann_chordal_to_subspace_E(V_full[:, 1:k+1], E_idx)
                # ---- E.resect observed ----
                delta_E, r   = chordal_full_vs_resect(
                                 V_full[:, 1:k+1], V_R[:, 1:k+1], N_idx)

                # ---- nulls (R draws, both variants share draws) ----
                eps_null   = []
                delta_null = []
                for r in 1..R:
                    S_r = strength_stratified_draw(quint, n_per_q[E_idx], pool=N_idx)
                    eps_S, _  = grassmann_chordal_to_subspace_E(V_full[:, 1:k+1], S_r)
                    eps_null.append(eps_S)
                    # null for resect: re-eigendecompose with S_r as "epi" set
                    A_R_S = A[V_p \ S_r, :][:, V_p \ S_r]
                    L_R_S = diag(A_R_S.sum(axis=1)) − A_R_S
                    λ_R_S, V_R_S = eig(L_R_S)
                    delta_S, _ = chordal_full_vs_resect(
                                   V_full[:, 1:k+1], V_R_S[:, 1:k+1], V_p \ S_r)
                    delta_null.append(delta_S)

                z_align  = (eps_E   − mean(eps_null))   / std(eps_null)
                z_resect = (delta_E − mean(delta_null)) / std(delta_null)

                emit row (p, b, φ, k, eps_E, f_E, delta_E, r,
                          z_align, z_resect)

# ---- cohort aggregation (descriptive only; no acceptance gate) ----
for each (b, k, φ):
    for variant in {align, resect}:
        obs_vec = [observed(p, b, k, φ, variant) for p in COHORT_N9_EPI]
        z_vec   = [z(p, b, k, φ, variant)        for p in COHORT_N9_EPI]
        emit row (b, k, φ, variant,
                  median(obs_vec), Q1(obs_vec), Q3(obs_vec),
                  median(z_vec),   Q1(z_vec),   Q3(z_vec),
                  n_pos_obs   = #{p : obs(p) below null mean},
                  n_pos_z     = #{p : z(p) < 0},
                  wilcoxon_p_two_sided(z_vec),
                  bh_q_within_band)
```

## 6. Visualisation spec

PDF only, full vector. No `fig.suptitle`. Adjacency-axis labels use
math `$i$`, `$j$` per `feedback_no_watermark_default.md`. Watermark
opt-in via `watermark=False`.

**Figure 1 — Cohort z-score heatmap (one variant per page; two-page PDF).**

For each variant `v ∈ {align, resect}`, a 6 × 3 grid (rows = bands,
cols = phases). Each cell is a small heatmap:

- Rows of heatmap = patients (9 rows + cohort-median row at the bottom).
- Columns = k ∈ K_GRID (6 columns, log spacing on x-axis labels).
- Cell value = `z_v(p, b, k, φ)`.
- Colormap: diverging `RdBu_r` centred at 0, range `[-3, 3]`.
- No regularity / regime annotations. We read the heatmap directly.

**Figure 2 — Cohort spaghetti per band, one panel per phase (3 panels horizontal).**

For each variant, x = `k`, y = `z_v`. One faint line per patient
(n=9 lines), thick line for cohort median, shaded IQR. Horizontal
line at 0 (null mean). One panel per phase. Pat_03 plotted in a
distinct style (1024 Hz outlier).

**Figure 3 — Raw observed scalars (no z-scoring).**

For each variant, mirror of Figure 2 with y = raw `ε^E_align` /
`δ^E_resect` and the per-patient null mean shown as a faint
dashed reference per patient. Lets us see the raw scale, not just
the z-distance from null.

**Figure 4 — `E_p` vs `E_p^cp` overlay (sensitivity panel, optional).**

Mirror of Figure 2 but with `E_p` (solid) and `E_p^cp` (dashed) on
the same axes per band per phase. Reveals whether any pattern is
cross-probe-driven (long-range epi network) or same-probe-driven
(focal epi clustering).

**Figure 5 (optional, deferred).**

Per-patient phase trajectory at any (band, k) cell that looks
interesting in the static plots: x = phase, y = `z_v`, line per
patient. Computed only after we look at Figures 1–3 together.

## 7. Connection to prior tools

| Prior tool | Relation |
|--|--|
| `2026-04-29_e1-spectral-subspace-alignment.md` (E1) | E1: cross-phase Grassmann distance between `V_k^pre`, `V_k^tt`, `V_k^post` (full-graph eigenspaces). Direction E: cross-region Grassmann distance between `V_k` and the epi coordinate subspace (E.align), or between full vs resected eigenspaces (E.resect). Same primitive `d_chord`; different leaf-set partitioning. E1 is global cross-phase; E is epi cross-region. Both can run independently; cross-rung interpretation is layered if both are cohort-positive. |
| `2026-05-08_epi-eigenmode-localization.md` (Direction C) | E.align = cumulative-over-modes form of Direction C. Identity `f^E_k = (1/k) Σ_{j=2}^{k+1} m^E_j`. So a regularity in Direction C at low normalised eigenvalue λ̂ is essentially the same finding as a regularity in E.align at low k. **Direction E.align does not subsume Direction C** — C resolves per mode (which mode peaks where), E.align integrates. **Direction E.resect is genuinely new** — does not reduce to any C-style measure. |
| `audit_37_e1_grassmann_triangle.py` | Has the `chordal()` primitive in script-private form. Promote to `lrg_eegfc.utils.metrics.spectral.chordal_distance` on this PR (≥ 3 callers: audit_37 + audit_46 + audit_61/E). |
| `audit_46_grassmann_principal_angles.py` | Decomposes the cross-phase chordal distance into per-angle contributions. Direction E does not need this decomposition for the lead claim, but if a regime is cohort-positive, the per-angle decomposition becomes a useful follow-up panel (deferred). |
| `lrg_eegfc.workflow.lrg.load_lrg_result` | Eigenvectors already cached (`r.eigenvectors`, shape `(N, N)`); no cache extension required. Confirmed 2026-05-08 on `Pat_02/alpha_rest_post_lrg_imcoh-abs.npz`. |
| `lrg_eegfc.utils.io.patient.load_epileptic_nodes` | Same `E_p` loader as Direction A and C. |
| `2026-05-08_epi-cross-phase-rigidity.md` (Direction A, audit_54) | Independent measures on independent objects (induced KC tree distances vs Grassmann subspace distance). Joint reading: both cohort-positive in the same band ⇒ epi has both hierarchical *and* spectral signatures; Direction E positive alone ⇒ spectral signature without hierarchical reorganisation; A positive alone ⇒ hierarchical signature without spectral. Pre-registered cross-rung interpretation table (deferred). |
| `2026-05-08_trace-minus-epi-resection.md` (Direction A', audit_56) | A' did *tree-level* resection (induced subtree on `V \ E_p` of the cached linkage). Direction E.resect does *FC-level* resection (delete rows/cols of A, recompute LRG eigendecomposition from scratch). Stronger version of the same conceptual question. A' answered "is the β trace non-epi-carried at the dendrogram level"; E.resect can answer "do epi nodes carry information that doesn't exist on the resected eigenspace". Possible diagnostic table: A' β-stable + E.resect β-positive ⇒ β trace is non-epi-tree-carried *and* epi nodes hold spectral info beyond the tree representation. |
| `2026-04-29_eigenvector-direct-pivot-plan.md` step 0 (cache extension) | Step 0 was scoped to extend the LRG NPZ to include eigenvectors. **Already complete in the imcoh_abs cache** (verified). The eigenmode-localization scope's "BLOCKING" note about cache extension is stale. |

What this measure **subsumes:** nothing in the current repo (cross-region
Grassmann embedding on a clinical node subset has not been done).

What this measure **complements:** Direction C (per-mode resolution) and
Direction A (hierarchical position). The three together test whether
epi nodes have a fingerprint at hierarchical, spectral-mode, and
subspace-rotation levels respectively.

What this measure **does not replace:** the global E1 trace test stays
as the cohort-wide cross-phase spectral-rotation rung; this measure is
epi cross-region.

## 8. Implementation plan

### Slot allocation

- **Compute script:** `scripts/01_compute/audit/audit_61_epi_grassmann_compute.py`
  (slot 61 is the next genuinely free slot as of 2026-05-08:
  55=ctm_per_pair_scatter, 56=anatomy_distribution + trace_minus_epi
  (double), 57=cohort_n_audit_figure, 58=band_agnostic_lrg,
  59=eigenmode_embedding, 60=anchor_anatomy_baseline + 60b — see §9
  open question 8).
- **Figures script:** `scripts/01_compute/audit/audit_61b_epi_grassmann_figures.py`.

### Library reuse (no new helpers in `scripts/`)

- `lrg_eegfc.workflow.lrg.load_lrg_result` for cached eigenpairs.
- `lrg_eegfc.workflow.fc.load_fc_matrix` for `Â_φ^p,b` (needed for
  E.resect re-eigendecomposition).
- `lrg_eegfc.utils.io.patient.load_epileptic_nodes` for `E_p`.
- `lrg_eegfc.utils.io.patient.load_channel_labels` for cross-probe
  variant `E_p^cp`.
- `lrg_eegfc.utils.metrics.hypothesis.{wilcoxon_z, bh_fdr,
  rank_biserial, boot_ci_mean}` for the secondary stats footnote.
- `lrg_eegfc.config.const.{PATIENTS_LIST, BRAIN_BANDS_NAMES, PHASE_LABELS}`.

### Library promotions (this PR)

**New module:** `lrg_eegfc.utils.metrics.spectral` with three callers
already (audit_37 already has its own private `chordal()`,
audit_46 already has `principal_angles()` + `chordal_from_angles()`,
new audit_61 needs both). Per the ≥ 2-caller library-promotion rule,
this module promotion is overdue.

API:

```python
def principal_angles(V_a: np.ndarray, V_b: np.ndarray) -> np.ndarray:
    """Principal angles in radians via SVD of V_a.T @ V_b."""

def chordal_distance(V_a: np.ndarray, V_b: np.ndarray) -> float:
    """Chordal Grassmann distance √(min(p,q) − Σ σᵢ²) on orthonormal bases.
    Both bases must be in a common ambient ℝ^D; V_a is D×p, V_b is D×q."""

def grassmann_to_coord_subspace(
    V_k: np.ndarray, idx: Sequence[int]
) -> tuple[float, float]:
    """Chordal distance from span(V_k) to coord-subspace span{e_i : i ∈ idx}.

    Returns (eps_align, f_idx) where
        eps_align = √(min(k, |idx|) − ‖V_k[idx, :]‖_F²)
        f_idx     = ‖V_k[idx, :]‖_F² / k    ∈ [0, 1]
    """

def chordal_full_vs_resect(
    V_k_full: np.ndarray, V_k_resect: np.ndarray, retained_idx: Sequence[int]
) -> tuple[float, int]:
    """δ^E_resect: chordal distance between V_k_full restricted to
    retained rows (orthonormalized via QR) and V_k_resect.

    Returns (delta_resect, rank_full_block).
    """
```

audit_37 and audit_46 should be updated in the same PR to import from
the library (separate commits within the PR for the migration to keep
the diff readable).

### Outputs

```
data/audit/epi_grassmann/
├── E_per_patient.csv           # rows: (p, b, φ, k, variant, value, z, …)
├── E_cohort.csv                # rows: (b, φ, k, variant, median, iqr, n_pos, regularity, p_wilcoxon, q_bh)
├── E_regimes.csv               # rows: (b, φ, variant, k_start, k_end, …) for regularity-passing regimes
├── figures/
│   ├── fig_01_cohort_heatmap_align.pdf
│   ├── fig_01_cohort_heatmap_resect.pdf
│   ├── fig_02_spaghetti_align.pdf
│   ├── fig_02_spaghetti_resect.pdf
│   ├── fig_03_variant_comparison.pdf
│   └── fig_04_cp_overlay.pdf      (optional)
```

### Scope partitioning

E.align (Variant 1) is the **first deliverable** — fast, cache-only,
~50 lines using existing helpers.

E.resect (Variant 2) is layered after E.align lands and the library
promotion is in place. Adds the eigendecomposition loop; ~30 minutes
wall-clock cohort run.

E.trace (cross-phase extension, §2.6) only runs if the static
signature (E.align or E.resect) is cohort-positive in any (band, k,
phase) cell. Per the trace-not-the-default rule.

### Sanity gate (BLOCKING before cohort run)

Mirror E1's SBM gate (E1 §8) at smaller scope:

- 4-block SBM, 25 nodes/block, intra-block weight 0.8, inter-block 0.1
  + Gaussian noise σ=0.05.
- Designate `E_p` = block 1 (25 nodes). Compute `ε^E_align`. **Pass:**
  at k=3 (= n_blocks − 1), `ε^E_align(SBM, E_p=block1)` is much smaller
  than `ε^E_align(SBM, E_p=random25)` (the Fiedler / mode-2 vector should
  be block-aligned).
- Reflexive: `chordal_distance(V, V) < 1e-5`.
- For E.resect: removing block 1 from a 4-block SBM should leave a
  3-block SBM whose top-3 eigenspace is well-defined; verify
  `δ^E_resect(SBM, E_p=block1) > δ^E_resect(SBM, E_p=random25)` at k=2.

CLI flag `--sanity-only` runs just the gate.

### CLI

```
python scripts/01_compute/audit/audit_61_epi_grassmann_compute.py [--sanity-only] [--align-only | --resect-only] [-v]
```

### Cost

- E.align cohort run: < 30 s.
- E.resect cohort run: ~ 5–10 minutes (parallelised).
- Figures: ~ 30 s.

## 9. Open questions

1. **Sign reading.** Both directions of departure from the null are
   reported per cell (no pre-registered "predicted direction").
   `ε^E_align < null` means modes lean toward epi; `> null` means
   modes lean away from epi. Either is informative; we look at the
   cohort distribution post hoc.

2. **Pool inclusion for the strength-stratified null.** Default
   excludes `E_p` from the pool (canonical "non-epi controls").
   Including `E_p` in the pool (broader sample, looser reference)
   is a sensitivity variant. Both are reported in `E_per_patient.csv`.

3. **k-grid choice.** `{2, 3, 5, 8, 13, 21}` per the eigenvector
   pivot plan and E1. Open question whether epilepsy-relevant scales
   are the same. A finer log-grid at small k (2, 3, 4, 5, 6, 8, 13,
   21) is a sensitivity check but expensive for E.resect (one
   eigendecomposition per k for each null draw).

4. **Cross-band reporting.** Per-band cohort summaries are the
   primary table. Cross-band patterns (same k showing structure in
   multiple bands) are themselves potentially interesting; inspected
   from the figures, not pre-registered.

5. **Pat_03 separate handling.** Per the standing rule, Pat_03
   (1024 Hz outlier) is included but flagged. Specifically: cohort
   summary reported both with and without Pat_03 in
   `E_cohort.csv`.

6. **Resection-graph degeneracy.** If `r < k` in the QR step (the
   non-epi rows of `V_k^full` rank-drop), E.resect is technically
   undefined for the missing dimensions. Default: report the
   reduced-rank chordal distance and flag the cell with `r < k` in
   `E_per_patient.csv`. Rank-drop is itself a finding (epi rows were
   necessary for the full span).

7. **Comparison with eigenvalue-shift `Δλ_k`.** E.resect measures
   *direction* changes, not *speed* changes. A complementary scalar
   is `Δλ_k = (λ_k^full − λ_k^resect) / λ_max^full` (the relative
   eigenvalue shift). Reported as a per-cell column in
   `E_per_patient.csv`.

8. **Slot assignments for Direction C and Direction E.** As of
   2026-05-08, audit slots 55–60 are all taken
   (55=ctm_per_pair_scatter, 56=anatomy_distribution +
   trace_minus_epi double-occupied, 57=cohort_n_audit_figure,
   58=band_agnostic_lrg, 59=eigenmode_embedding,
   60=anchor_anatomy_baseline + 60b figures). Direction E lands at
   audit_61 (next genuinely free slot). Direction C's scope file
   (`2026-05-08_epi-eigenmode-localization.md`) still claims
   audit_57, which is taken — that scope's slot claim is stale.
   When Direction C is implemented, it will need a fresh slot
   (audit_62 or later); cross-link the renumber back here.
