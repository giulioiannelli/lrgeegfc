---
name: epi-eigenmode-localization
type: scope
era: IMCOH_ABS × COHORT_N10
status: draft
created: 2026-05-08
updated: 2026-05-08
pointers:
  - .agents/plans/active/2026-05-08_lrg-epilepsy-research-directions.md
  - .agents/plans/active/2026-04-29_eigenvector-direct-pivot-plan.md
  - .agents/guides/task-persistence-investigation/2026-04-29_e1-spectral-subspace-alignment.md
  - .agents/reports/2026-05-07_epileptic-n10-revisit.md
  - .agents/guides/02_methods/lrg-framework-guide.md
---

# Epi eigenmode localization — IPR / projection mass on E_p

## Renormalization head

**For each Laplacian eigenmode v_k of the FC graph at fixed τ = 1/λ_max,
compute the fraction of mode mass concentrated on the epi node set
m^E_k = Σ_{i ∈ E_p} v_k(i)² ∈ [0, 1]. Plotted against the normalised
eigenvalue λ̂_k = λ_k / λ_max ∈ [0, 1], a cohort-consistent peak in
m^E_k at some λ̂ range identifies a "characteristic communication
scale" at which the epileptic network's modes do not propagate to
the rest of the brain — direct analogue of Anderson localisation in
disordered solid-state systems. Direction C is split into a static
signature (C1: does the peak exist in any single phase, cohort-wide?)
and an optional phase-trace extension (C2: does the localisation peak
shift across rest_pre → task → rest_post?). The degree-corrected null
is mandatory: in graphs with heavy-tailed degree distributions,
high-strength nodes accumulate eigenmode mass without any "trapping"
mechanism — the null must subtract out this trivial localisation.
Depends on the eigenvector cache extension scoped in the
2026-04-29 eigenvector pivot plan; cannot run until the cache is
extended.**

## Notation

(Inherits LRG primitives from `2026-04-29_e1-spectral-subspace-alignment.md`
and `lrg-framework-guide.md`.)

- `V_p` — node set for patient `p`, size `N_p`.
- `E_p ⊆ V_p` — epi node set (channel-matched). Same as Direction A.
- `Â_φ^p` — symmetric weighted FC matrix (|ImCoh|, with the bias
  diagonal handled per the LRG pipeline) at phase `φ ∈ {pre, tt, post}`,
  band `b`.
- `D̂_φ^p` — diagonal of row sums of `Â_φ^p`.
- `L̂_φ^p` = `D̂_φ^p − Â_φ^p` — symmetric weighted Laplacian.
- Eigendecomposition: `L̂_φ^p v_k = λ_k v_k`, ordered
  `0 = λ_1 < λ_2 ≤ … ≤ λ_{N_p}`. The trivial mode `v_1` (constant)
  is always discarded.
- `λ̂_k = λ_k / λ_max` ∈ (0, 1] — normalised eigenvalue. The
  normalisation makes spectra comparable across patients despite
  per-patient `λ_max` differing.
- `‖v_k‖_2 = 1` (unit-norm convention).
- `s_i = Â_φ^p[i, :].sum()` — node strength of node `i`.

## Definitions

### Mode projection mass on E_p

```
m^E_k(φ, b, p) = Σ_{i ∈ E_p} v_k^φ,b,p(i)²
```

Range: `m^E_k ∈ [0, 1]`. By unit-norm,
`Σ_{i ∈ V_p} v_k(i)² = 1`, so `m^E_k` is the fraction of the mode's
squared amplitude concentrated on the epi set.

Random-baseline expectation under "uniform" null (no localisation):
`E[m^E_k | random] ≈ |E_p| / N_p`.

### Inverse Participation Ratio (alternative, optional)

```
IPR_k = Σ_i v_k(i)^4    (range (1/N, 1])
IPR^E_k = (Σ_{i ∈ E_p} v_k(i)²)² / Σ_{i ∈ E_p} v_k(i)^4
```

`IPR_k` is the standard solid-state localisation measure (large
when `v_k` is concentrated on few nodes). `IPR^E_k` is its
restriction to E_p — a "localisation strength among epi nodes
given that the mode is on the epi set". Defer `IPR^E_k` to the
secondary analysis; primary measure is `m^E_k`.

### Eigenvalue binning

Define a fixed normalised-eigenvalue grid:

```
B = 50 bins on [0, 1]: bin b = [b/B, (b+1)/B], b = 0..B-1
```

For each (patient, phase, band), bin the modes:
`Bin(λ̂_k) = floor(λ̂_k · B)`.

Per-bin median projection mass:

```
m̄^E_b(φ, b, p) = median_{k : Bin(λ̂_k) = b} m^E_k(φ, b, p)
```

(Use median, not mean — robust to single-mode outliers.)

### Degree-corrected null

Because high-strength nodes carry disproportionate eigenmode mass
in heavy-tailed graphs, a uniform random subset is not the right
null. Instead: stratified shuffling on node strength.

For each patient, partition `V_p` into `Q = 5` strength quintiles
based on per-phase, per-band node strength `s_i`. Let `q_p^i` be
the quintile of node `i`. Let `n_p^q` = `|{i ∈ E_p : q_p^i = q}|`
— number of epi nodes per quintile.

Strength-stratified null draw `S_r`: in each quintile, sample
`n_p^q` nodes without replacement from `{i ∈ V_p : q_p^i = q}`.
Concatenate to get `S_r` of size `|E_p|` with the same
node-strength profile as `E_p`.

Random null:
```
m^{S_r}_k = Σ_{i ∈ S_r} v_k(i)²
```
Aggregate per bin: `m̄^{S_r}_b`.

Per-bin z-score:
```
z̄^E_b(φ, b, p) = (m̄^E_b − μ̂(m̄^{S_r}_b)) / σ̂(m̄^{S_r}_b)
```
where `μ̂` and `σ̂` are sample mean and std over `R = 100` null draws.

### Static cohort signature (C1)

For each (band, normalised-eigenvalue bin `b`), aggregate per-patient
z-scores over the cohort:

- **Cohort-Wilcoxon** one-sided (`z̄^E_b > 0`) on per-patient z-scores.
- **Cohort consistency**: count patients with `z̄^E_b > 1.96`
  (≥ 8 / 9 ⇒ cohort-positive bin).

A **localisation regime** is a contiguous range of bins
`[b_start, b_end]` all marked cohort-positive. Bin width 0.02; a
localisation regime of width ≥ 0.06 (≥ 3 contiguous bins) is a
"resolved peak"; narrower regimes are flagged as "punctate".

C1 default: pool over phases (use phase-averaged `m̄^E_b(p, b)` per
patient); phase-resolved peaks are C2.

### Phase-trace extension (C2)

For each (band, bin `b`), compute the per-patient localisation triangle:

```
ΔL^E_pre→tt   = z̄^E_b(tt) − z̄^E_b(pre)
ΔL^E_tt→post  = z̄^E_b(post) − z̄^E_b(tt)
T_loc^E_b     = ΔL^E_tt→post − ΔL^E_pre→tt   (mirrors KC trace formulation)
```

Cohort-Wilcoxon one-sided on `T_loc^E_b` per (band, bin) — significant
in trace direction means the localisation strength changes during
task and the change persists into rest_post.

C2 is layered on C1: only run C2 in (band, bin) cells where C1 is
cohort-positive (else underpowered).

## Properties

- **Range.** `m^E_k ∈ [0, |E_p|/N_p]` for delocalised modes, up to
  `1` if v_k is fully concentrated on E_p.
- **Invariance.** `m^E_k` is invariant under sign flip of `v_k`
  (we use squared amplitude). It is **not** invariant under
  permutation of eigenvector ordering at degenerate eigenvalues —
  for any patient with degenerate eigenvalues (rare in our
  full-rank weighted graphs but possible at tiny λ_k near zero),
  flag and exclude from the analysis.
- **Identifiability of localisation regimes.** A regime is
  identifiable iff it spans ≥ 3 contiguous bins (width ≥ 0.06 in
  normalised eigenvalue) and survives the degree-corrected null at
  z > 1.96 in ≥ 8 / 9 patients. Narrower / weaker regimes are
  ambiguous between localisation and degree-confound residuals.
- **Sample size.** n=9 patients × 6 bands × 50 bins × 3 phases =
  8 100 z-scores. After cohort-Wilcoxon per (band, bin) cell, BH-FDR
  within band (m = 50) gives the multiple-comparison correction.
- **Negative properties (what this measure does NOT detect).**
  - **Localisation patterns split across small disjoint subsets of
    E_p.** If E_p has two clusters and modes localise on one or
    the other but not both, `m^E_k` averages and may miss the
    structure. Mitigation: run a within-E_p substructure analysis
    if cohort-positive (cluster E_p anatomically and compute m^E_k
    per cluster).
  - **Localisation that is not visible at the second moment** (`v_k(i)²`).
    Higher-moment localisation (e.g. a few large-amplitude
    contributions among many small ones) may be better captured by
    `IPR^E_k` — deferred to secondary analysis.
  - **Cross-band coherence of localisation regimes.** A regime that
    appears at the same `λ̂` in multiple bands is not detected by
    this measure (each band is an independent dendrogram); a
    follow-up cross-band aggregation is layered if cohort-positive.
- **Complexity.** Per (patient, phase, band): one eigendecomposition
  (already cached after the eigenvector pivot plan extension) + R
  null draws × N_p mass evaluations. With R = 100, N ≈ 115, n = 9
  patients, 6 bands, 3 phases: ≈ 2 × 10^6 mass evaluations. <
  1 minute.

## Caveats & failure modes

| | Caveat | Mitigation |
|--|--|--|
| 1 | Eigenvector cache extension (one-line LRG NPZ change) is required and currently not done. | Block C until the cache extension lands per the eigenvector pivot plan. |
| 2 | Degenerate eigenvalues (near `λ_2`) produce non-unique eigenvectors. | Filter modes with `min(λ_{k+1} − λ_k, λ_k − λ_{k-1}) < 1e-6`; report fraction filtered per patient. |
| 3 | Sign convention of `v_k` is not unique; `m^E_k` is sign-invariant by construction. | Use squared amplitude. |
| 4 | `λ̂_k` normalisation can compress spectra non-linearly across patients with very different `λ_max`. | Report per-patient λ_max distribution; if cohort spread > 5×, also report on the unnormalised log-eigenvalue axis. |
| 5 | Heavy-tailed degree distribution can produce trivial "localisation" on high-strength nodes. | Strength-stratified null (Q = 5 quintiles) is mandatory; uniform-random null is reported only as a sanity check showing the problem. |
| 6 | Pat_15 has no epi annotation. | Drops out (n = 9, same as Direction A). |
| 7 | Same-probe / cross-probe split (analogous to Direction A `E_p^cp`). | Optional cross-probe-restricted variant `E_p^cp` reported alongside `E_p`. |
| 8 | The Anderson-localisation analogy is loose — in graphs, "localisation" is more degree-driven than disorder-driven. | Frame finding as "epi-projected eigenmode mass exceeds degree-corrected null", not as a literal Anderson transition. |
| 9 | Per-band variation in spectrum shape (lower-frequency bands often have flatter spectra). | Per-band binning + per-band cohort-Wilcoxon is independent across bands. |
| 10 | Pat_03 outlier (1024 Hz, larger eigenvalue spread). | Log Pat_03 separately; flag in cohort plot. |

## Pseudocode

```
input:  Â_φ^b^p (FC matrix), E_p, R = 100, B = 50, Q = 5
output: per-(p, b, φ): m̄^E_b for b = 0..B-1
        per-patient z̄^E_b
        cohort: mask of cohort-positive bins per band

# 1. Per-patient eigendecomposition (one-shot per phase, band)
L = D - A    # weighted Laplacian
λ, V = eig(L)
sort_idx = argsort(λ); λ = λ[sort_idx]; V = V[:, sort_idx]
λ_max = λ[-1]
λ̂_k = λ / λ_max          # normalised
modes_keep = [k for k in 1..N-1 if min_gap(λ, k) > 1e-6]

# 2. Strength quintiles (degree-corrected null)
s = A.sum(axis=1)
quintile = quantile_bin(s, Q)

# 3. Per-mode m^E_k
m_E = [sum(V[i, k]^2 for i in E_p) for k in modes_keep]

# 4. Bin into normalised-eigenvalue grid
bin_id = floor(λ̂_k * B) for each k in modes_keep
m̄^E_b = median(m_E[k] for k where bin_id[k] == b)  for b = 0..B-1

# 5. Strength-stratified null
m̄^null_b = empty array (R, B)
for r in 1..R:
    S_r = empty
    for q in 0..Q-1:
        n_q = count(i in E_p with quintile[i] == q)
        S_r += sample without replacement of size n_q from
              {i : quintile[i] == q and i not in S_r}
    m_S = [sum(V[i, k]^2 for i in S_r) for k in modes_keep]
    m̄^null_b[r, :] = bin_median(m_S, bin_id, B)

μ_null = mean(m̄^null_b, axis=0)
σ_null = std(m̄^null_b, axis=0)
z̄^E_b = (m̄^E_b - μ_null) / σ_null

# 6. Cohort aggregator (over patients for fixed (b, b_idx))
for band in BANDS:
    for b_idx in 0..B-1:
        z_vec = [z̄^E_{b_idx}(p) for p in cohort]
        p_one = wilcoxon_z(z_vec, alternative='greater')
        consistent = sum(z_vec > 1.96) >= 8
    BH-FDR over b_idx per band

# 7. Resolve regimes
for band:
    contiguous_pos_runs = [(start, end) where all bins in [start, end]
                            are BH-significant AND consistent]
    regimes = filter(width >= 3) of contiguous_pos_runs
```

## Visualization spec

**Figure 1 — Cohort heatmap of z̄^E_b across the bands.**

6 × 1 grid (rows = bands), each row is a heatmap:

- **Rows of heatmap** = patients (9 rows + cohort-median row).
- **Columns** = normalised eigenvalue bin `b` ∈ [0, 1].
- **Cell value** = `z̄^E_b(phase-pooled, p)`.
- **Colormap**: diverging RdBu_r centred at 0, range [-3, 3].
- **Annotation**: top of each band row, mark cohort-positive bins
  with a black bar above the heatmap. Width of the bar = regime width.

**Reading rule**: a vertical band of red across all 9 rows of a band
panel = a cohort-consistent localisation regime. Scattered red =
single-patient effects. White = null.

**Figure 2 — Cohort-median line per band.**

Single panel: x = `λ̂` ∈ [0, 1], y = cohort-median `z̄^E_b`,
one line per band (6 lines). Shaded IQR. Horizontal line at 1.96.
Reads at a glance: which bands have peaks, where, how high.

**Figure 3 — `E_p` vs `E_p^cp` (cross-probe restricted) overlay.**

Mirror of Figure 2 but with `E_p` (solid) and `E_p^cp` (dashed) on
the same axes. Difference reveals whether localisation is
cross-probe-driven (long-range epi network) or same-probe-driven
(focal epi clustering).

**Figure 4 (optional, C2 only) — Phase-trace per cohort-positive bin.**

For each cohort-positive (band, bin) cell from C1, a small panel:
3-phase boxplot of `z̄^E_b(φ, p)` across patients. Trace direction
is `z̄^E_b(post) > z̄^E_b(pre)` and similar to `z̄^E_b(tt)`. Stars:
BH-FDR `q < 0.05` cohort-Wilcoxon trace direction.

## Connection to prior tools

| Prior tool | Relation |
|--|--|
| `2026-04-29_e1-spectral-subspace-alignment.md` (E1) | E1 measures rotation of the *full top-k subspace* between phases. C measures *localisation of individual modes* on E_p. Different math, different question. E1 is a global cross-phase distance; C is a per-mode mass restricted to a specific node set. |
| `2026-04-29_eigenvector-direct-pivot-plan.md` (E1/E2/E3 rungs) | Direction C is a fourth eigenvector-direct rung scoped from epilepsy rather than task-trace. Shares the cache extension dependency. |
| `lrg_eegfc.workflow.lrg.compute_lrg_analysis` | Eigendecomposition lives here. The cache extension (eigenvector pivot plan step 0) is the prerequisite. |
| `lrg_eegfc.utils.io.patient.load_epileptic_nodes` | Same E_p loader as Direction A. |
| `2026-05-08_epi-cross-phase-rigidity.md` (Direction A) | Independent measures on independent objects (induced subtree distances vs eigenmode mass). Joint interpretation in the manuscript: a cohort-positive Direction A class + a cohort-positive Direction C regime in the same band is corroborative; flat A + sharp C is "no hierarchical signature, but communication-mode signature"; sharp A + flat C is "hierarchical signature without spectral fingerprint". |
| Anderson localisation literature | Loose analogy. Cited in framing only; we do not claim a localisation transition. |
| `lrg_eegfc.utils.metrics.hypothesis.wilcoxon_z`, `bh_fdr` | Direct reuse for cohort tests. |

What this measure **subsumes**: nothing in the current repo
(eigenmode-localisation analysis on a clinical node subset has not
been done).

What this measure **complements**: Direction A (hierarchical
position) and Direction E1 (subspace rotation). The three
together test whether epi nodes have a fingerprint at hierarchical,
spectral-mode, and subspace-rotation levels respectively.

What this measure **does not replace**: the global E1 trace test
(if completed) stays as the cohort-wide spectral-rotation rung;
this measure is epi-specific.

## Implementation plan

### Prerequisite (BLOCKING)

The eigenvector cache extension scoped in
`2026-04-29_eigenvector-direct-pivot-plan.md` step 0 must be
implemented first. The cached LRG NPZ currently stores `eigenvalues`
and `eigenvectors` for some pipelines but not consistently across
the |ImCoh| LRG cache. Required: confirm `IMCOH_LRG_CACHE` NPZs
contain `eigenvectors` and `eigenvalues`; if not, extend
`compute_lrg_analysis` and re-cache n=10 × 6 bands × 3 phases.

(Per the eigenvector pivot plan, this is a one-line change; the
re-cache is ~480 small files and runs in minutes.)

### Implementation

- **Script:** `scripts/01_compute/audit/audit_57_epi_eigenmode_compute.py`
- **Library reuse (no new helpers in scripts/):**
  - `lrg_eegfc.workflow.lrg.load_lrg_result` for cached eigenpairs.
  - `lrg_eegfc.utils.io.patient.load_epileptic_nodes` for `E_p`.
  - `lrg_eegfc.utils.metrics.hypothesis.wilcoxon_z`, `bh_fdr`.
- **Helper to add (library, not script):**
  `lrg_eegfc.utils.metrics.spectral.epi_projection_mass(V, E_p)` —
  computes `m^E_k` for all modes. Promote to library because Direction
  C C1 + C2 + downstream variants are 3 callers.
- **Outputs:**
  - `data/audit/epi_eigenmode/M_per_patient.csv` — one row per
    (patient, band, phase, bin) with `z̄^E_b`.
  - `data/audit/epi_eigenmode/M_cohort.csv` — one row per
    (band, bin) with cohort-Wilcoxon p, BH q, count, regime flag.
  - `data/audit/epi_eigenmode/figures/fig_01_cohort_heatmap.pdf`
  - `data/audit/epi_eigenmode/figures/fig_02_cohort_lines.pdf`
  - `data/audit/epi_eigenmode/figures/fig_03_cp_overlay.pdf`
  - `data/audit/epi_eigenmode/figures/fig_04_phase_trace.pdf` (C2 only)
- **Figures script:**
  `scripts/01_compute/audit/audit_57b_epi_eigenmode_figures.py`.

### Scope partitioning

C1 (static signature) is the **first deliverable** and gates whether
C2 runs. C1 alone is publishable as "epileptic networks have a
characteristic spectral localisation regime in band X at λ̂ ∈ [a, b]".

C2 (phase-trace extension) is conditional on C1 being cohort-positive
in at least one (band, bin) cell. Without C1, C2 is underpowered.

## Open questions

1. **Bin width.** B = 50 bins gives 0.02 width on `λ̂` ∈ [0, 1]. Too
   fine and per-bin medians are noisy; too coarse and regimes are
   indistinguishable. 50 is a defensible default; sensitivity check
   at B ∈ {25, 50, 100} in the first sweep.

2. **Strength quintiles vs continuous-strength matching.** Q = 5 is
   coarse; a continuous-matching null (KDE on `s_i` then Mahalanobis
   resampling) is stricter but more complex. Defer; quintile null
   is the canonical first pass.

3. **Median vs mean per-bin aggregation.** Median is robust but
   loses information when regimes are bimodal within a bin. Both
   are reported in the per-patient CSV; figures use median.

4. **Phase pooling for C1.** Per-phase z-scores can be averaged
   over phases (the default) or computed per-phase and combined via
   meta-analysis (Stouffer's method). Defer.

5. **Cross-band aggregation.** If multiple bands have peaks at the
   same `λ̂`, that is itself a cohort-level finding (cross-band
   coherence of epi localisation). Defer; cross-band tests are
   layered after per-band tests close.

6. **Pre-registration of acceptance gate.** Cohort-positive (C1)
   = ≥ 8 / 9 patients with `z̄^E_b > 1.96`, regime width ≥ 3
   contiguous bins, BH-FDR `q < 0.05` cohort-Wilcoxon, in any one
   band. C2 trace direction = ≥ 8 / 9 patients with
   `T_loc^E_b < 0`, BH-FDR `q < 0.05`.

7. **Anderson-localisation framing in the manuscript.** Use as
   loose analogy in the discussion only. The empirical claim is
   "epi-projected eigenmode mass exceeds degree-corrected null in
   regime X" — not a transition.

8. **Visualisation choice if no regime survives.** If C1 returns
   null at the regime-level but per-bin signal exists, fall back
   to a "narrowest cohort-positive bin" report. The negative
   result form is a publishable null: "epileptic networks do not
   show a spectral localisation fingerprint beyond what
   degree-matched controls give".
