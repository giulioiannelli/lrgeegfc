> **⚠ QUANTITATIVE_STALE (2026-04-15)** — numbers in this report were
> computed under the pre-reset ImCoh mislabelling (stored `|ImCoh|²`
> under the name `imcoh`). Qualitative conclusions survive (rankings
> preserved under sqrt); regenerate absolute values against
> `fc_method="imcoh_abs"` before quoting. See
> `~/.claude/projects/.../memory/imcoh_taxonomy.md`.

# ImCoh Process Report: From Same-Probe Bias to Solution

This document records the full journey from discovering the same-probe MSC
bias to the ImCoh solution and its validated LRG results. It is intended as
a comprehensive reference for any agent working on this project.

---

## 1. Discovery of the Problem (March 2026)

During the epileptic node investigation, we noticed that LRG communities at
coarse scales (n=10-30) did not reflect functional brain organization but
instead followed the physical geometry of the sEEG electrode probes.

### Key observations

- **Same-probe MSC is 2-8x higher than cross-probe MSC.** Contacts on the
  same electrode probe (e.g., A1-A12) show artificially inflated coherence
  because they are physically adjacent (~3-5 mm spacing) and share both
  volume conduction and common-reference noise.

- **LRG communities at coarse scales follow probe geometry.** At n=10
  communities, same-community/same-probe enrichment is 3-8x. At n=30 this
  can reach 4-8x. The dendrogram merges same-probe contacts first because
  they have the shortest ultrametric distance.

- **Any community-level metric is confounded.** Participation coefficient,
  community cosine, metastability, allegiance -- all inherit the probe bias
  at scales n < ~15-20.

### How it was detected

The epileptic node investigation repeatedly produced results where nodes on
the same probe as the epileptic focus appeared "special" in the LRG
hierarchy. Visual inspection of 3D brain connectome plots (coloured by
community) immediately revealed that communities were coloured clusters along
electrode trajectories, not functionally meaningful groups.

Quantification via same-community/same-probe enrichment ratio confirmed the
bias systematically across all 5 patients. Variation of Information (VI)
between probe labels and community labels was misleading -- it can appear
"high" even when communities clearly follow probe geometry (due to many-to-one
mappings). The enrichment ratio is the correct diagnostic.

---

## 2. Failed Approaches

Four debiasing approaches were investigated before arriving at ImCoh.

### 2a. Zero same-probe edges

Set all MSC values between contacts on the same probe to zero.

- **Why it fails:** Ad hoc. Destroys genuine within-probe local circuit
  dynamics. Would not survive peer review -- reviewers would correctly ask
  why real connectivity was discarded rather than corrected.

### 2b. Percentile rescaling

Rescale same-probe edge weights to match the cross-probe distribution by
percentile mapping.

- **Why it fails:** Entirely ad hoc with no physical or statistical
  justification. The rescaled values do not correspond to any meaningful
  coherence quantity. Would not survive peer review.

### 2c. Bipolar re-referencing

Subtract consecutive contacts on the same probe (V_{k+1} - V_k) before
computing MSC. This is standard in clinical sEEG practice.

- **Why it fails: MAKES THE BIAS WORSE.** The same-probe ratio goes from
  ~5x (common reference) to 8-17x (bipolar). This happens because adjacent
  bipolar channels share a physical contact: bipolar channel k = V_{k+1} -
  V_k and channel k+1 = V_{k+2} - V_{k+1} both contain V_{k+1} with
  opposite sign. This shared contact creates a new, stronger source of
  spurious coherence that was not present in the common-reference montage.

### 2d. Distance regression

Fit a linear model of MSC vs. inter-contact distance and subtract the
predicted component.

- **Why it fails:** The relationship between MSC and distance is not linear
  (it is highly nonlinear with a sharp drop at same-probe distances). A
  linear model is insufficient and no principled nonlinear model exists.

---

## 3. The Physics of the Problem

The same-probe MSC inflation has two physical sources:

1. **Common reference (V_G2):** All channels are recorded as V_i - V_G2.
   The reference electrode G2's activity appears in every channel, creating
   a shared signal component that inflates coherence between all pairs
   (but especially same-probe pairs which are already correlated).

2. **Volume conduction:** Electric fields propagate through brain tissue at
   effectively the speed of light. Physically adjacent contacts (same probe,
   ~3-5 mm apart) record overlapping neural populations.

Both sources share one critical physical property: they are **instantaneous**.
Instantaneous mixing produces **zero phase lag** between the mixed components.
This means the cross-spectral density (CSD) contribution from these sources
is **purely real** (the imaginary part is exactly zero).

MSC uses |S_ij|^2 / (S_ii * S_jj), which includes both real and imaginary
parts of the CSD. It therefore captures the full volume conduction artifact.

---

## 4. The Solution: Imaginary Coherence (ImCoh)

### The key insight

The imaginary part of coherency, introduced by Nolte et al. (2004, >3000
citations), exploits the physics described above:

1. Volume conduction is instantaneous.
2. Instantaneous mixing produces zero phase lag.
3. Zero phase lag means Im(S_ij) = 0 for the volume-conducted component.
4. Therefore, Im(S_ij)^2 / (S_ii * S_jj) is exactly zero for any volume-
   conducted signal.
5. Any nonzero ImCoh must reflect genuine neural interaction with a nonzero
   time delay.

This is not an approximation or statistical correction. It is an exact
consequence of the physics.

### Implementation

ImCoh was implemented as a `metric` parameter in the existing MSC computation
functions. No separate computation path was needed:

```python
# Same function, different metric parameter
freqs, Coh = compute_msc_welch(X, fs=2048, nperseg=4096, metric="imcoh")

# Or via the pipeline
adj_matrices = coherence_fc_pipeline(X, fs=2048, nperseg=4096, metric="imcoh")
```

Computation time is identical to MSC (~35 seconds per patient/phase) because
the Welch CSD estimation is shared and only the final derivation step changes
(taking Im(S_ij)^2 instead of |S_ij|^2).

### Same-probe bias results

| Metric | Same-probe ratio range | Typical |
|--------|----------------------|---------|
| MSC    | 1.6-6.0x             | 3-5x    |
| ImCoh  | 0.6-2.0x             | ~1x     |

The residual ImCoh ratio (especially ~2x in beta for Pat_07) reflects genuine
lagged neural coupling between nearby contacts, not artifact.

---

## 5. The Caching Bug

### What happened

After computing ImCoh matrices and running them through `compute_lrg_analysis()`,
the resulting dendrograms appeared identical to the MSC-derived ones. This was
initially interpreted as a deep result: "LRG is invariant to MSC vs ImCoh
because `compute_normalized_linkage()` normalizes by max distance."

### What actually happened

`compute_lrg_analysis()` uses the `fc_method` parameter to construct cache
file paths. When called with `fc_method="msc"` (the default) while passing
an ImCoh matrix as input, it found and returned the existing MSC cache file
instead of computing fresh from the ImCoh matrix. The ImCoh matrix was never
actually processed.

### How it was caught

Exhaustive verification of all 120 cases (5 patients x 4 phases x 6 bands)
using raw lrgsglib functions (bypassing the caching layer) showed that
120/120 cases produce DIFFERENT dendrograms -- different merge pairs,
different merge heights, and community overlap at n=5 of only 4-15%.

### How to avoid

1. Use `fc_method="imcoh"` for proper cache path separation
2. Use `use_cache=False` to force recomputation
3. Or bypass `compute_lrg_analysis()` entirely:
   ```python
   from lrgsglib.core import compute_normalized_linkage, compute_ultrametric_distance
   dist = compute_ultrametric_distance(imcoh_matrix)
   linkage = compute_normalized_linkage(dist)
   ```

The `fc_method="imcoh"` option is now accepted by `compute_lrg_analysis()`
for proper cache separation going forward.

---

## 6. Results: ImCoh Strengthens the Main Findings

With the caching bug fixed and correct ImCoh LRG results computed, the
hypothesis testing was repeated:

| Hypothesis | Description | MSC cells | ImCoh cells | Change |
|------------|-------------|-----------|-------------|--------|
| H1 | Task stability | 112 | 193 | +72% |
| H2a | Task trace in rest | 33 | 107 | +224% |
| H2b | Task approach | 31 | 19 | -39% |
| H3 | Within < cross | 32 | 38 | +19% |

"Cells" refers to unanimous cells in the unanimity map (patient x band x
scale combinations where all patients agree on the direction of the effect).

Key findings:

- **H1 strengthens massively.** Task phases are more stable under ImCoh,
  presumably because removing the dominant volume-conduction component
  reveals the true task-related functional reorganization.
- **H2a strengthens dramatically (+224%).** The task trace in resting-state
  connectivity is much more visible when volume conduction is removed.
- **H2b narrows.** Task approach survives only in beta band.
- **High_gamma shows no signal** in any hypothesis under ImCoh.

---

## 7. Band-Specific H2a Breakdown

H2a (task trace: rsPre vs rsPost differ, indicating that the task leaves a
lasting imprint) shows strong band specificity under ImCoh:

| Band | Unanimous+ cells | k range | All 5 patients positive? |
|------|-------------------|---------|--------------------------|
| beta | 36 | k=16-58 | Yes -- strongest signal |
| low_gamma | 31 | k=21-58 | Yes |
| alpha | 29 | k=17-58 | No (Pat_02 barely negative, -0.02) |
| delta | 6 | — | No (Pat_02 negative) |
| theta | 5 | — | No (Pat_02 and Pat_08 negative) |
| high_gamma | 0 | — | No signal |

Beta is the dominant band for task-trace effects, with all 5 patients
showing positive effects across a wide range of scales. Low_gamma is the
second strongest. High_gamma is completely absent.

H2b (task approach: rsPre approaches task) survives only in beta band
(16 cells, k=13-58).

---

## 8. Files Produced

### ImCoh FC matrices (120 files)
```
data/imcoh_cache/{patient}/{band}_{phase}_imcoh_*.npy
```
One file per patient/phase/band combination. 5 patients x 4 phases x 6
bands = 120 files.

### ImCoh LRG results (120 files)
```
data/imcoh_lrg_cache/{patient}/{band}_{phase}_lrg_imcoh.npz
```
Computed using raw lrgsglib functions (bypassing the caching bug). Contains
linkage matrix, ultrametric distance matrix, and community assignments.

### VI profiles
```
data/imcoh_vi/vi_raw_profiles.csv
data/imcoh_vi/hypothesis_contrasts.csv
```
Variation of Information profiles across scales for all patient/band/phase
combinations, and the derived hypothesis contrast statistics.

### Unanimity maps
```
data/imcoh_unanimity/unanimity_map_H1.pdf
data/imcoh_unanimity/unanimity_map_H2a.pdf
data/imcoh_unanimity/unanimity_map_H2b.pdf
data/imcoh_unanimity/unanimity_map_H3.pdf
```
Visual maps showing which (band, scale) cells are unanimous across patients.

### Scalar contrasts
```
data/imcoh_unanimity/scalar_contrasts.csv
data/imcoh_unanimity/unanimity_summary.md
```
Numerical summary of hypothesis cell counts and per-band breakdowns.

---

## 9. Code Changes

### `src/lrg_eegfc/utils/fc/msc/msc.py`
- Added `metric` parameter to `compute_msc_welch()`.
- Accepts `"msc"` (default, legacy), `"imcoh"`, or `"wpli"`.
- The Welch CSD estimation is shared; only the final derivation step differs.

### `src/lrg_eegfc/utils/fc/msc/__init__.py`
- Propagated `metric` parameter to `coherence_fc_pipeline()`.
- All downstream functions pass through the metric choice.

### `src/lrg_eegfc/workflow/lrg.py`
- Accepts `fc_method="imcoh"` for proper cache path separation.
- Previously only accepted `"msc"` and `"corr"`.

### `src/lrg_eegfc/utils/fc/msc/sparsify.py`
- Added `zero_same_probe_edges()` -- zeros MSC between same-probe contacts.
- Added `rescale_same_probe_edges()` -- percentile rescaling (deprecated).
- Added `distance_regression_rescale()` -- distance regression (deprecated).
- These are the failed approaches documented in section 2, kept for reference.

### `src/lrg_eegfc/utils/io/patient.py`
- Added `bipolar_rereference()` -- bipolar re-referencing of raw signals.
- Added `parse_seeg_label()` -- parses probe name and contact number from
  channel label.
- Added `load_epileptic_nodes()` -- loads epileptic node annotations.

---

## 10. Surrogates and ImCoh

Circular-shift surrogates are **WRONG** for ImCoh.

### Why

Circular shifts preserve the power spectrum but scramble the phase
relationships between channels. However, they do not specifically preserve
the zero-lag structure. The shifted signals develop artificial phase lags
at all frequencies. ImCoh detects these artificial lags as "genuine"
coupling, producing a systematically inflated null distribution.

### Consequence

If surrogates are used for sparsification (keeping only edges where MSC >
95th percentile of null), the threshold is too high and real ImCoh signal
is discarded.

### Recommendation

- Use **dense ImCoh** (sparsify="none") for all LRG analyses.
- If sparsification is needed for other purposes, use the **disparity
  filter** or **ECM** (edge-centric modularity), neither of which requires
  surrogates.
- Do NOT use circular-shift surrogates with ImCoh.

---

## 11. Recommendations

1. **Use ImCoh as the primary FC metric** for all LRG-based analyses where
   same-probe bias matters (which is essentially all of them).

2. **Keep MSC as the default in the codebase** until the full analysis
   pipeline is replicated end-to-end with ImCoh. This avoids breaking
   existing workflows.

3. **When calling `compute_lrg_analysis()` with ImCoh,** ALWAYS either:
   - Use `fc_method="imcoh"` for proper cache separation, OR
   - Use `use_cache=False` to force recomputation, OR
   - Bypass the function entirely with raw lrgsglib functions.

4. **The wPLI implementation exists** (`metric="wpli"`) but has performance
   issues (Python loop over segments). Use ImCoh unless wPLI is specifically
   needed for its theoretical advantages at small phase lags.

5. **High_gamma produces no signal** under ImCoh in any hypothesis. This
   band can likely be excluded from future analyses to save computation time.

6. **Beta is the dominant band** for task-trace (H2a) and task-approach (H2b)
   effects under ImCoh. Focus resources on beta and low_gamma for the main
   findings.

---

## Related Files

- ImCoh guide: `.agents/guides/IMCOH_GUIDE.md`
- Probe bias guide: `.agents/guides/PROBE_BIAS_GUIDE.md`
- MSC method guide: `.agents/guides/MSC_METHOD_GUIDE.md`
- Memory: `imcoh_lrg_invariance.md` (corrected)
- Memory: `probe_bias_critical.md`
