---
name: continuous-trace-matrix
type: scope
era: COHORT_N10
status: current
created: 2026-04-26
updated: 2026-04-26
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-25_mrl-vs-cbr-reconciliation.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cohesion-cbr.md
  - .agents/reports/2026-04-24_multiscale-task-trace.md
  - scripts/01_compute/hypothesis_tests/h2c_ultrametric_drift.py
  - scripts/01_compute/hypothesis_tests/continuous_trace_matrix.py
  - scripts/01_compute/hypothesis_tests/continuous_trace_controls.py
  - scripts/01_compute/hypothesis_tests/h2e_split_half.py
  - src/lrg_eegfc/workflow/lrg.py
---

# Continuous-trace matrix — visual proof on the ultrametric distance matrix, before the dendrogram throws information away

**The ultrametric distance matrix `D^{p,b,φ} ∈ ℝ^{N×N}_{≥0}` is the
object the LRG dendrogram is built from. Every dendrogram is a
discretisation of `D` that loses the per-pair magnitudes of distance
shifts. The cohort-wide trace claim is a continuous statement about
those shifts: pairs `(i, j)` whose distance changed during task tend
to change in the same direction during the post-rest. This scope
specifies the descriptive, visual companion to the existing H2c test
(per-patient Spearman ρ on `Δ_task` vs `Δ_rest` in the upper triangle):
per-patient distance-shift heatmaps + sign-agreement map + per-pair
scatter; cohort overlay of per-patient ρ. No thresholds. No dendrogram
discretisation. The proof is the matrix.**

## 1. Notation

Indices and sets:

| Symbol | Domain | Meaning |
|:---|:---|:---|
| `p ∈ P` | `|P| = 10` (`COHORT_N9` + restored `Pat_14`) | patient |
| `b ∈ B` | `|B| = 6` | band |
| `φ ∈ Φ = {pre, learn, test, post}` | | phase (snake_case canonical) |
| `L_p` | `|L_p| = N_p` | leaf (contact) index set |

Per `(p, b, φ)`:

- `D^{p,b,φ} ∈ ℝ_{≥0}^{N_p × N_p}` — ultrametric distance matrix from
  the LRG. Symmetric, zero diagonal. Stored condensed in
  `LRGResult.ultrametric_matrix`; recovered as a square matrix via
  `scipy.spatial.distance.squareform`.

Upper-triangle indexing: `(i, j)` with `i < j`, total of
`m_p := N_p (N_p − 1) / 2` pairs per patient.

Distance-shift vectors (length `m_p` each):

```
Δ_task^{p,b}(i, j) := D^{p,b,test}(i, j) − D^{p,b,pre}(i, j)
Δ_rest^{p,b}(i, j) := D^{p,b,post}(i, j) − D^{p,b,pre}(i, j)
```

Both real-valued. Positive `Δ_task` = pair moved apart during task;
negative = moved closer.

## 2. Definitions

### 2.1 Per-pair co-shift

```
S^{p,b}(i, j) := Δ_task^{p,b}(i, j) · Δ_rest^{p,b}(i, j) ∈ ℝ
```

`S > 0` ⇔ both shifts have the same sign (pair moves consistently in
task and post). `S < 0` ⇔ opposite signs (pair moves in task and
reverses in post). The signed product encodes "trace direction" at
the pair level; `|S|` encodes magnitude.

### 2.2 Per-pair sign-agreement

```
σ^{p,b}(i, j) := sgn(Δ_task) · sgn(Δ_rest)   ∈ {−1, 0, +1}
```

(Convention: `sgn(0) = 0`, so pairs with no movement in either phase
contribute 0.) The sign-agreement strips magnitude information; useful
for the binary "same direction yes/no" reading at the pair level.

### 2.3 Per-patient continuous trace correlation

```
ρ^{p,b} := Spearman( vec_upper(Δ_task^{p,b}), vec_upper(Δ_rest^{p,b}) )   ∈ [−1, +1]
```

This is the existing H2c statistic per patient (`h2c_ultrametric_drift.py`).
`ρ > 0` = trace direction; cohort-wide claim: most patients have
`ρ > 0` per band, with `mean ρ` significantly above 0 against the
H2e split-half null.

### 2.4 Cohort-pooled summary

For each band `b`, two cohort-level objects:

- **Per-patient distribution**: `{ρ^{p,b} : p ∈ P}` — 10 values per band.
- **Pooled per-pair scatter**: union of `(Δ_task^{p,b}(i,j), Δ_rest^{p,b}(i,j))`
  across all `p ∈ P` and `(i, j)` upper-triangle pairs. Total
  `Σ_p m_p ≈ 10 · 115^2 / 2 ≈ 66 000` points per band. Visual density
  on the (Δ_task, Δ_rest) plane shows the trace as a positive-diagonal
  cloud.

## 3. Properties

| # | Statement |
|:--|:----------|
| 3.1 | `Δ_task, Δ_rest ∈ ℝ` per pair; `S ∈ ℝ`; `σ ∈ {−1, 0, +1}`; `ρ ∈ [−1, +1]`. |
| 3.2 | `Δ` is **anatomy-subtracted by construction**: the static structural component of `D` cancels under subtraction (it appears in both `D^pre` and `D^test`), leaving only the shift. The continuous trace claim is therefore on the *residual* signal — exactly where the prior project found it (see `2026-04-24_multiscale-task-trace.md` §7c). |
| 3.3 | Asymmetric in time: `Δ_task` measures pre→test, `Δ_rest` measures pre→post. The trace claim is `corr(Δ_task, Δ_rest) > 0`. The mirror question (do pairs that moved during *learn* also move in post?) uses `Δ_learn := D^learn − D^pre`. |
| 3.4 | Multiscale by virtue of `D` itself: `D` already encodes the multiscale structure of the LRG hierarchy (the cophenetic / ultrametric distance for pair `(i,j)` is the merge height of their LCA in the dendrogram, summing both depth and breadth of subtree shared structure). Looking at `Δ` instead of pre/post separately localises *which scale* changed — a small `Δ` for a pair with high `D^pre` is a deep-scale shift; a small `Δ` for a pair with low `D^pre` is a fine-scale shift. |
| 3.5 | Scale-comparable across patients via `D` normalisation: `D` ranges `[0, dmax]` per (p, b, φ). For cohort scatter pooling, normalising each patient's `Δ` by `dmax(D^{p,b,pre})` makes the magnitudes comparable (see §6.4). |
| 3.6 | NOT a discrete-module claim. This scope cannot say "subtree X persisted from task into post" — for that, see `2026-04-25_cohesion-cbr.md` (per-leaf classification with soft affinities). It says "the *pattern of distance shifts* is consistent between task and post" — a residual geometric claim. |
| 3.7 | NOT a hypothesis test. The cohort `ρ` distribution is descriptive. The corresponding inferential test exists already at `h2c_ultrametric_drift.py` with the H2e split-half null. This scope produces the *visual* companion to that test. |
| 3.8 | Reuses primitives — `load_lrg_result(...).ultrametric_matrix` (condensed), `scipy.spatial.distance.squareform`, `scipy.stats.spearmanr`. No new helpers. |

## 4. Caveats & failure modes

### 4.1 Anatomy retained at the pair level

Even after subtracting `D^pre`, pairs on the same anatomical neighbourhood may continue to share *correlated noise* between `D^test` and `D^post` — same probe, same drift, same line noise. This inflates `ρ` slightly even with no task-induced trace. Mitigation: report `ρ` against the H2e split-half drift null already computed for H2c (`data/reports/imcoh_vi/h2e_split_half.csv` if available). Visual mitigation: split scatter into *same-probe pairs* vs *cross-probe pairs* and confirm the trace direction holds in both.

### 4.2 Same-probe bias on `D` magnitudes

Same-probe pairs have systematically smaller `D` (probe-bias guide). They also have systematically smaller `|Δ|` because their distances vary in a narrow range. The pair-pooled scatter will have most same-probe points clustered near the origin. Cross-probe pairs carry the load-bearing cohort signal. Mitigation: colour same-probe pairs separately in the per-pair scatter; report pooled ρ both with and without same-probe pairs in the cohort plot.

### 4.3 Pat_03 1024 Hz outlier

`D` magnitudes for Pat_03 are systematically inflated (probe-bias guide and per-patient quirks). Normalising `Δ` by `dmax(D^pre)` per (p, b) (§6.4) makes the cohort scatter pooling defensible. Without normalisation, Pat_03 dominates the scatter visually.

### 4.4 Pair-pool independence violation

The `m_p ≈ 6 600` upper-triangle pairs per (p, b) are not statistically independent — pairs share leaves. The Spearman ρ is fine as a *descriptive* effect size, but inference on the pooled scatter would over-state precision. Per-patient ρ is the correct inference unit. The visual scatter shows the trace; the cohort dot of `{ρ^{p,b}}` shows whether it generalises.

### 4.5 N varies across patients

`Pat_10` is 113 channels; the rest 115–116. Cohort pooling on the scatter mixes patients with different `m_p`. This is fine for visual density but should be flagged — patients with more pairs contribute proportionally more points. Per-patient ρ aggregation does not have this issue.

### 4.6 Edge cases on `sgn`

`Δ_task = 0` is rare in continuous-valued distance shifts (pre-existing pairs that didn't move), but it does happen (especially same-probe pairs at a coherence floor). `σ(i, j) = 0` for any pair with either shift exactly zero. These pairs do not contribute to the sign-agreement reading.

### 4.7 What this measure cannot detect

- *Modules forming and dissolving*: a coherent group of leaves clustering together in task and dissipating in post. The continuous distance shift vanishes if the average pair-distance returns to baseline; the discrete CBR / Cohesion-CBR analysis catches this.
- *Asymmetric loosening vs tightening*: `Δ_task < 0` (pair came together) and `Δ_rest > 0` (pair moved apart) gives `S < 0` → no trace at this pair. A pair-level reset is invisible to the *correlation* (which averages out cancellations) but visible to the per-pair scatter (off-diagonal cells).

## 5. Pseudocode

```
Input:
  P, B, Φ
  fc_method = "imcoh_abs"
  same_probe_mask(p) — boolean N_p × N_p of same-probe pairs
                       (load_epileptic_nodes / channel_labels.csv)

For each (p, b) ∈ P × B:
  D_pre  ← squareform(load_lrg_result(p, "rest_pre",  b, fc_method).ultrametric_matrix)
  D_test ← squareform(load_lrg_result(p, "task_test", b, fc_method).ultrametric_matrix)
  D_post ← squareform(load_lrg_result(p, "rest_post", b, fc_method).ultrametric_matrix)
  N_p ← D_pre.shape[0]
  iu  ← np.triu_indices(N_p, k=1)             # upper triangle
  Δ_task ← (D_test - D_pre)[iu]
  Δ_rest ← (D_post - D_pre)[iu]
  ρ        ← scipy.stats.spearmanr(Δ_task, Δ_rest).statistic
  σ        ← sgn(Δ_task) * sgn(Δ_rest)        # in {-1, 0, +1}
  S        ← Δ_task * Δ_rest                  # signed product
  same_probe ← same_probe_mask(p)[iu]
  emit row(s):
    per-cell ρ summary, per-pair (Δ_task, Δ_rest, σ, S, same_probe)

Output:
  data/reports/imcoh_continuous_trace/per_cell_summary.csv
    columns: patient, band, N_p, m_p, rho, frac_pos_sigma, frac_neg_sigma,
             rho_no_same_probe
  data/reports/imcoh_continuous_trace/per_pair.parquet
    columns: patient, band, i, j, dD_task, dD_rest, sigma, S, same_probe
    (parquet for size; ~70k rows × 60 cells = 4M rows total)

For each band b:
  cohort_rho[b] ← {ρ^{p,b} : p ∈ P}
  pooled_pairs[b] ← concat per-patient (dD_task, dD_rest, same_probe)
                     with optional per-patient normalisation
                     dD_task / dmax(D_pre) and dD_rest / dmax(D_pre).

Render figures (§6).
```

### 5.1 Complexity

Per `(p, b)`: `O(N_p^2)` to build the upper-triangle vectors;
`O(m_p log m_p)` for Spearman ρ. Total per cohort: `60 cells · ~10^4
ops` ≈ a few seconds. Trivial.

The pooled scatter can have ~600 000 points per band (n_pairs × n_patients ≈ 6 600 × 10 × 6 bands) — needs hexbin or 2D-histogram rendering, not raw scatter.

## 6. Visualization spec

Two figure families.

### 6.1 Per-patient × band figure (one PDF per (patient, band))

4 panels in a row, sharing aspect:

- **Panel A — `Δ_task` heatmap (N×N)**. Diverging colormap (`RdBu_r`),
  centered at 0, symmetric `vmin = −max|Δ_task|, vmax = +max|Δ_task|`.
  Diagonal masked. `imshow(..., rasterized=True)`. Title:
  *Δ_task = D_test − D_pre*.
- **Panel B — `Δ_rest` heatmap (N×N)**. Same colormap and limits as A.
  Title: *Δ_rest = D_post − D_pre*. **Reading rule**: cells that are
  red in A should mostly be red in B (and blue ↔ blue). Visual trace.
- **Panel C — sign-agreement map `σ(i, j)` (N×N)**. Three-colour
  categorical: same-direction (yellow), opposite-direction (purple),
  zero (grey). Diagonal masked. Annotation in corner: fraction
  same-direction `Π_+ = (σ > 0).mean()`, fraction opposite
  `Π_− = (σ < 0).mean()`.
- **Panel D — per-pair scatter `Δ_task` vs `Δ_rest`**. Hex-bin or
  2D histogram (`m_p ≈ 6 000` points). Diagonal `y = x` line, `x = 0`
  and `y = 0` axes. Annotation: `ρ`. **Reading rule**: a positively
  oriented elliptical cloud = trace direction; a circular blob
  centred at origin = no trace; a negatively oriented ellipse =
  anti-trace.

No `fig.suptitle`. Sidecar `.md` carries the caption per band.

PDF only, all heatmaps/hexbins rasterised at `dpi=200` per the
project rule. One PDF per (patient, band): `data/reports/imcoh_continuous_trace/figures/per_cell/<patient>_<band>.pdf`.

### 6.2 Cohort figure (one PDF per band)

Three panels:

- **Panel A — per-patient `ρ` dotplot for this band**. 10 dots,
  cohort median + IQR, reference at 0, optional H2e null shading.
  Reading: dots above 0 = trace, below = anti-trace.
- **Panel B — pooled-cohort hex-bin** of normalised `(Δ_task /
  dmax_pre, Δ_rest / dmax_pre)` across all patients. ~60 000 points;
  hexbin with log density. Diagonal line. Annotation: cohort-pooled
  Spearman ρ on the pooled vectors and per-patient median.
- **Panel C — sign-agreement bar chart**: per patient, two stacked
  bars (yellow `Π_+`, purple `Π_−`, grey `Π_0`). Visual cohort
  unanimity readout.

`data/reports/imcoh_continuous_trace/figures/cohort_<band>.pdf`. PDF only, rasterised hexbin.

### 6.3 Cross-band cohort figure (one PDF, all bands)

Single panel: `ρ^{p,b}` arranged by band on the x-axis, dot per
patient, cohort median + IQR per band. Same layout as
`H3_imcoh_abs.pdf`-style band-by-band reading. Replicates the H2c
result visually. PDF: `data/reports/imcoh_continuous_trace/figures/cohort_all_bands.pdf`.

## 7. Connection to prior tools

| Prior measure | Relation to this scope |
|:---|:---|
| **H2c** (`h2c_ultrametric_drift.py`) | Computes the per-patient `ρ` already; reports cohort-level Wilcoxon. This scope is its **visual companion** — the figures show *what the ρ summarises*. The CSV produced by H2c can be reused directly for Panel A of the cohort figure. |
| **H2e** (`h2e_split_half.py`) | The split-half drift floor providing `ρ_null ≈ 0` per band — overlay on the cohort dotplot if available. |
| **MRL family** (this folder, superseded) | MRL operates on *discrete subtree identity* and lost the continuous-shift signal by thresholding. This scope is the antidote — the same data, pre-discretisation. |
| **Cohesion-CBR** (`audit_12`) | Per-leaf discrete classification on dendrograms. Complementary, not redundant: Cohesion-CBR shows *which leaves* exhibit each pattern; this scope shows *how distances between contact pairs shift*. Both can be in the same paper. |
| **H2-RAW / H2-FROB** (`2026-04-24_multiscale-task-trace.md` §7c) | Measure global similarity post→test vs post→pre on raw `D`. They fail because anatomy dominates baseline `D`. The Δ formulation removes the static anatomy component — same data, better-conditioned question. |

## 8. Implementation plan

### 8.1 Compute script

`scripts/01_compute/hypothesis_tests/continuous_trace_matrix.py` — new, ~200 lines.

CLI:
```
python scripts/01_compute/hypothesis_tests/continuous_trace_matrix.py \
    [--patients Pat_02,...] \
    [--bands delta,theta,...] \
    [--fc-method imcoh_abs] \
    [--no-same-probe-split] \
    [-v]
```

Imports: `load_lrg_result`, `squareform`, `scipy.stats.spearmanr`,
`load_epileptic_nodes` (for same-probe mask), `wilcoxon_z`/`bh_fdr` if
running the cohort inferential layer (re-use H2c's CSV instead by
default).

Outputs:
- `data/reports/imcoh_continuous_trace/per_cell_summary.csv` —
  one row per (patient, band): N_p, m_p, ρ, frac_pos_sigma,
  frac_neg_sigma, ρ_no_same_probe, dmax_pre.
- `data/reports/imcoh_continuous_trace/per_pair.parquet` —
  long-format per-pair frame for the figures.
- `data/reports/imcoh_continuous_trace/cohort_summary_<band>.csv`
  per band (10 patient ρ values, cohort median, IQR).

### 8.2 Figure scripts

Three:

- `scripts/01_compute/figures_embedded/fig_continuous_trace_per_cell.py`
  — per-(patient, band) 4-panel figure (§6.1). Loops the cohort.
  Output: `data/reports/imcoh_continuous_trace/figures/per_cell/<patient>_<band>.pdf`.
- `scripts/01_compute/figures_embedded/fig_continuous_trace_cohort.py`
  — per-band cohort 3-panel figure (§6.2). Output:
  `data/reports/imcoh_continuous_trace/figures/cohort_<band>.pdf`.
- `scripts/01_compute/figures_embedded/fig_continuous_trace_all_bands.py`
  — cross-band cohort dotplot (§6.3). Output:
  `data/reports/imcoh_continuous_trace/figures/cohort_all_bands.pdf`.

### 8.3 Library promotion (deferred)

Per `coding-rules.md`, helpers stay script-local in v1. Promotion
candidates if a second consumer arises:
- `pair_distance_shifts(D_pre, D_post) → (Δ, iu)` to
  `lrg_eegfc.utils.metrics.pair_shifts`.
- `same_probe_mask(patient)` to `lrg_eegfc.utils.io.patient`.

### 8.4 Tests

`tests/test_continuous_trace.py`:

- Identity: `D_pre = D_test = D_post` → `Δ_task = Δ_rest = 0` →
  ρ undefined (handle as `np.nan`).
- Sign-agreement: synthetic `Δ_task = (a, b, c)` and
  `Δ_rest = (a/2, −b, c)` → expected `σ = (+1, −1, +1)`.
- Symmetry: `Δ_task` and `Δ_rest` are vectors of equal length; ρ is
  symmetric in arguments.

## 9. Open questions

1. **Pooled cohort scatter normalisation** — divide each patient's `Δ`
   by `dmax(D^pre)` before pooling? Recommended (§3.5) but this
   discards a small amount of information (patient-level magnitude
   variation is lost). Default: yes-normalised in the cohort hex-bin,
   raw in the per-cell heatmaps.
2. **Same-probe split** — always split into same-probe vs cross-probe
   pairs in Panel D and the cohort hex-bin? Recommended yes; helps
   defeat the probe-bias caveat. Default: on.
3. **Δ_learn alongside Δ_rest** — should we also produce per-cell
   figures with `Δ_learn = D^learn − D^pre` to confirm the trace
   direction is consistent across both task phases? Cheap to add;
   makes the per-cell figure 6 panels instead of 4.
4. **Null overlay on cohort dotplot** — pull H2e split-half null from
   the existing CSV if present; otherwise leave null absent and
   document it. Recommended yes if the file exists.
5. **Per-leaf back-projection (hybrid)** — could be a future scope.
   Aggregate `|S(i,j)|` per leaf `i` (sum or mean over `j ≠ i`) to
   get a per-leaf "trace strength" score. This bridges to the
   per-leaf colouring style of Cohesion-CBR while staying continuous.
   Out of scope here; flagged for follow-up.

## 10b. Results — descriptive run (2026-04-26)

Cohort `n = 10` under `imcoh_abs`, full-phase `D` matrices, shared `D_pre`
baseline (i.e. the original `H2c` formulation, descriptive companion).

| band | median ρ | mean ρ | n_+ / 10 |
|------|---------:|-------:|---------:|
| δ    | +0.452 | +0.477 | 10/10 |
| θ    | +0.354 | +0.379 | 10/10 |
| α    | +0.485 | +0.452 | 10/10 |
| β    | +0.577 | +0.528 | 10/10 |
| γ_l  | +0.388 | +0.490 | 10/10 |
| γ_h  | +0.509 | +0.472 | 9/10  |

59/60 cells positive. Visually convincing per-cell figures (Pat_06 α
sanity check ρ=+0.725; per-band cohort PDFs at
`data/reports/imcoh_continuous_trace/figures/`). Reproduces the prior
`H2c` 53/54-cell positive finding at the now-restored cohort `n = 10`.

The descriptive run is not the final claim — it shares `D_pre` between
`Δ_task` and `Δ_rest`, which inflates ρ via shared-baseline noise.
Quantification of that artefact and the three-control consolidation
that follows it are in §11.

## 11. Three-control consolidation (2026-04-27)

The shared-baseline ρ in §10b has three independent failure modes that
must be neutralised before the cohort claim is publishable.

### 11.1 Run A — independent half-baselines

`D_pre_A`, `D_pre_B` from disjoint halves of the rest_pre time series
(LRG cache at `data/cache/imcoh_lrg_halves/`, populated by
`h2e_split_half.py`). Halved Welch `nperseg` keeps segment count
comparable. Define::

    Δ_task = D_test − D_pre_A,  Δ_rest = D_post − D_pre_B

`D_pre_A ⊥ D_pre_B` by construction → shared-baseline-noise correlation
removed. Population baseline `μ_pre` still cancels in expectation.
Implementation: `continuous_trace_matrix.py --mode split-baseline`.

### 11.2 Run B — probe-bias split

Same-probe pairs are anatomically near each other and trivially
correlated (probe-bias guide). Cross-probe pairs are the volume-
conduction-immune, load-bearing subset. Run A's per-pair vectors are
masked via `lrg_eegfc.utils.probe.build_probe_mask` and ρ is recomputed
on each subset.

### 11.3 Run C — drift floor (from H2e, same noise regime as Run A)

`ρ_null_drift = Spearman(D^pre_B − D^pre_A, D^post_B − D^post_A)` — pure
within-session drift on the resting recordings, no task involvement.
Pulled from `data/reports/imcoh_vi/h2e_split_half_rho_raw.csv`. Run A
and Run C live in the same halved-data noise regime, so paired
ρ_split > ρ_null_drift is a fair test.

### 11.4 Outcome — band-specific survival

Cohort medians (`n = 10`; null_drift `n = 9`, no Pat_14 in H2e):

| band | shared | split | drift | cross-probe | same-probe | split>drift | Wilcoxon p |
|------|-------:|------:|------:|------------:|-----------:|------------:|-----------:|
| δ    | +0.452 | +0.031 | +0.059 | +0.032 | +0.040 | 7/10 | 0.222 |
| θ    | +0.354 | −0.049 | −0.083 | −0.048 | −0.003 | 6/10 | 0.254 |
| **α**    | +0.485 | **+0.115** | −0.016 | **+0.105** | +0.198 | **8/10** | **0.008** |
| **β**    | +0.577 | **+0.222** | −0.043 | **+0.223** | +0.278 | **8/10** | **0.014** |
| **γ_l**  | +0.388 | **+0.140** | −0.020 | **+0.143** | +0.173 | **9/10** | **0.011** |
| γ_h  | +0.509 | −0.014 | −0.065 | −0.016 | +0.074 | 5/10 | 0.399 |

**Pass:** α, β, low_γ. All three bands clear the three pass criteria
(§10.1–10.3 below).
**Fail:** δ, θ, high_γ. Cohort split ≈ 50/50 after the artefact is
removed; the original signal was almost entirely shared-baseline
contamination.

The Run A drop ρ_shared − ρ_split = 0.25–0.52 quantifies how much of
the original `H2c`-style ρ was the formula talking to itself.

Outputs:
- `data/reports/imcoh_continuous_trace/per_cell_summary_split.csv`
- `data/reports/imcoh_continuous_trace/controls_summary.csv` — per-cell, all 5 ρ
- `data/reports/imcoh_continuous_trace/controls_band_stats.md` — per-band table
- `data/reports/imcoh_continuous_trace/run_A_baseline_split_summary.md`
- `data/reports/imcoh_continuous_trace/figures/baseline_split_overlay.pdf`
- `data/reports/imcoh_continuous_trace/figures/controls_overview.pdf`

Scripts:
- `scripts/01_compute/hypothesis_tests/continuous_trace_matrix.py`
  (`--mode {shared,split-baseline}`)
- `scripts/01_compute/hypothesis_tests/continuous_trace_controls.py`
- `scripts/01_compute/figures_embedded/fig_continuous_trace_baseline_split.py`
- `scripts/01_compute/figures_embedded/fig_continuous_trace_controls.py`

### 11.5 Final claim

> Task-induced reorganization of the LRG ultrametric distance matrix
> persists in the post-task rest in **α, β, and low_γ bands**, in
> **7–9 of 10 patients**. The residual signal exceeds the within-
> session drift floor (paired Wilcoxon `p ≤ 0.014`, ≥ 8/10 patients
> above the floor) and is not driven by same-probe anatomy
> (cross-probe ρ retains the cohort sign and magnitude).

Bands `δ, θ, high_γ` do **not** survive the controls. Their `H2c`
positive ρ was largely a shared-baseline + drift artefact.

## 10. Pass criterion (descriptive)

The continuous-trace-matrix passes the surface-evidence test if:

1. **Per-cell figures**: at least 3 patients per band exhibit a
   visually clear positive-diagonal cloud in Panel D and visibly
   correlated `Δ_task`/`Δ_rest` heatmaps in A/B.
2. **Cohort figure (per band)**: cohort median ρ > 0 in every band;
   ≥ 7/10 patients positive in every band (replicating the H2c
   53/54 finding now at n=10).
3. **Cross-band figure**: visible band-by-band consistency without
   the artificial flatness that the discrete-MRL forward score
   produced.
4. **Same-probe split** stays positive (or near-zero in cross-probe
   only) — confirms the trace is not a probe-bias artefact.

If 1–3 hold, the figures are publishable as the visual companion to
the H2c result. Together with the per-patient discrete proof from
Cohesion-CBR (`per_patient_hierarchy_cohesion/`), this gives the two
load-bearing visual statements about the trace:

- *Cohort, continuous*: pairs that moved during task move the same
  way in post (this scope).
- *Per-patient, discrete*: leaves bunch in task and post, scatter in
  pre (Cohesion-CBR).
