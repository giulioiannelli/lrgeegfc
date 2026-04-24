---
name: imcoh-guide
type: guide
era: IMCOH_ABS
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# ImCoh (Imaginary Part of Coherency) Guide

## Why MSC Is Contaminated

Magnitude-squared coherence (MSC) measures the total linear coupling between
two signals at each frequency. It captures **all** sources of coherence
indiscriminately, including:

- **Volume conduction:** physically adjacent sEEG contacts record the same
  neural population's activity, creating spurious coherence.
- **Common reference:** all channels are referenced to G2, so G2's activity
  appears in every channel, inflating coherence between all pairs.
- **Shared far-field potentials:** distant strong sources (e.g., muscle,
  cardiac) appear identically at nearby contacts.

These contamination sources all share one physical property: they are
**instantaneous** (speed-of-light propagation through brain tissue).
Instantaneous mixing produces **zero phase lag** between the mixed
components, which means the cross-spectral density (CSD) contribution
from these sources is **purely real** (zero imaginary part).

MSC uses `|S_ij|^2 / (S_ii * S_jj)`, which includes both real and
imaginary parts of the CSD. It therefore captures the full volume
conduction artifact.

In our data, this manifests as same-probe MSC being 2-8x higher than
cross-probe MSC, dominating the network structure and confounding all
downstream analyses (LRG hierarchies, community detection, metastability).

## Why ImCoh Fixes This

The imaginary part of coherency (ImCoh), introduced by Nolte et al. (2004),
exploits the physics of volume conduction to eliminate it:

1. Volume conduction is instantaneous (speed of light through tissue).
2. Instantaneous mixing produces zero phase lag between source and sensor.
3. Zero phase lag means the cross-spectral density contribution is purely
   real: `Im(S_ij) = 0` for the volume-conducted component.
4. Therefore, `Im(S_ij) / √(S_ii · S_jj)` is **exactly zero** for any
   volume-conducted signal.
5. Any nonzero ImCoh **must** reflect genuine neural interaction with a
   nonzero time delay (i.e., real functional connectivity).

This is not an approximation or a statistical correction. It is an exact
consequence of the physics. ImCoh is immune to volume conduction by
construction.

**Key insight:** ImCoh does not remove volume conduction from the data --
it uses a metric that is mathematically blind to it. The underlying CSD
computation is identical; only the final step (taking the imaginary part
instead of the magnitude) changes.

## Implementation (post-2026-04-15 reset)

ImCoh lives in the new `src/lrg_eegfc/utils/fc/coherence/imcoh.py` module
and shares the Welch/CSD backend with MSC via `coherence/_common.py`. The
canonical signed quantity per Nolte 2004 is what is computed and cached.
Magnitude (`|ImCoh|`) and squared magnitude (`|ImCoh|²`) are derived at
load time by the workflow loader.

### fc_method taxonomy

| `fc_method` | Formula | Range | On disk | Notes |
|---|---|---|---|---|
| `imcoh` | `Im(S_ij) / √(S_ii · S_jj)` | `[-1, 1]` | **yes** (freq-resolved) | Nolte 2004 canonical. **Cannot feed LRG** (raises `ValueError`). |
| `imcoh_abs` | `<\|ImCoh\|>_f` | `[0, 1]` | no | Ewald 2012 / Bastos & Schoffelen 2016 convention; default for Section 2. |
| `imcoh_sq` | `<\|ImCoh\|²>_f` | `[0, 1]` | no | Ewald 2012 "squared imaginary coherence"; matches pre-reset archive. |

### Why freq-resolved caching

Band-averaging signed ImCoh values allows positive/negative phase-lag
bins within a band to cancel (Jensen's inequality): `abs(mean(signed)) ≠
mean(|signed|)` and `mean(signed)² ≠ mean(signed²)`. To compute each of
the three literature conventions correctly, the loader must apply the
transform per-frequency-bin **before** band-averaging. The cache stores
the freq-resolved `(N, N, F_band)` signed tensor, and the loader
dispatches:

```python
# pseudocode in workflow/fc.py:_load_imcoh
freq_resolved = np.load(cache_path)              # (N, N, F_band), signed
if coh_transform == "identity":
    return freq_resolved.mean(axis=-1)           # imcoh
if coh_transform == "abs":
    return np.abs(freq_resolved).mean(axis=-1)   # imcoh_abs
if coh_transform == "sq":
    return (freq_resolved ** 2).mean(axis=-1)    # imcoh_sq
```

### Recommended usage

```python
from lrg_eegfc.workflow.fc import load_fc_matrix

# Band-averaged |ImCoh| for Section 2 network analysis
A = load_fc_matrix("Pat_02", "rest_pre", "beta", "imcoh_abs")
# Signed ImCoh for directional / phase-lag analysis
A_signed = load_fc_matrix("Pat_02", "rest_pre", "beta", "imcoh")
```

### Low-level (compute the freq-resolved tensor)

```python
from lrg_eegfc.utils.fc.coherence import compute_imcoh
freqs, ImCoh = compute_imcoh(X, fs=2048, nperseg=4096)  # (N, N, F), signed
```

## Test Results: Same-Probe Bias (post-reset)

Recomputed 2026-04-15 from the canonical signed cache. Pat_02, rest_pre,
alpha band, nperseg=4096:

| `fc_method` | Same-probe mean | Cross-probe mean | Ratio |
|---|---:|---:|---:|
| `msc`       | 0.251 | 0.046 | **5.49×** |
| `imcoh_abs` | 0.069 | 0.061 | **1.13×** |
| `imcoh_sq`  | 0.009 | 0.007 | 1.35× |

`msc` retains its ~5.5× inflation (volume conduction). `imcoh_abs`
(`|ImCoh|`, the literature magnitude convention) achieves a 1.13× ratio
— effectively no bias. `imcoh_sq` is 1.35× because squaring amplifies
differences and keeps a small residual ratio.

Note: the pre-reset "ImCoh = 1.35×" number in older reports was
actually `imcoh_sq`. The post-reset `imcoh_abs` ratio of 1.13× is the
correct number to quote when the paper narrative says `|ImCoh|`.

The same-probe bias drops from 5.5x to 1.35x -- essentially eliminated.
The residual 1.35x ratio likely reflects genuine local connectivity between
nearby contacts, not artifact.

## wPLI: Also Implemented, Use With Caution

The weighted Phase Lag Index (wPLI, Vinck 2011) is also implemented in
`compute_msc_welch(X, fs, metric="wpli")`. wPLI is also immune to volume
conduction and has some theoretical advantages over ImCoh (better noise
properties, less sensitivity to volume conduction leakage at small phase lags).

However, the current wPLI implementation has **performance issues**: it
requires per-segment cross-spectrum computation (a Python loop over segments),
making it significantly slower than MSC or ImCoh (which accumulate CSD across
segments in a single vectorized einsum).

**Recommendation:** Use ImCoh for now. wPLI can be explored later if ImCoh
proves insufficient, but performance optimization would be needed first.

## Cache Convention

ImCoh matrices are cached in a separate directory tree that mirrors the MSC
cache structure:

```
data/imcoh_cache/Pat_XX/{band}_{phase}_msc_sparsify-{method}_nperseg-{value}.npy
```

Note: the filename still contains `_msc_` for historical compatibility with
the cache-path generation code. The directory name (`imcoh_cache` vs
`msc_cache`) is what distinguishes the metric.

The MSC cache remains at:
```
data/msc_cache/Pat_XX/{band}_{phase}_msc_sparsify-{method}_nperseg-{value}.npy
```

## LRG Results: ImCoh Produces Different Dendrograms

> **CRITICAL:** MSC and ImCoh produce fundamentally different LRG hierarchies.
> An earlier analysis incorrectly concluded they were identical. That result
> was caused by a caching bug (see below).

### Exhaustive verification (120/120 cases different)

All 120 patient/phase/band combinations were checked. In every case, the
ImCoh-derived LRG dendrogram differs from the MSC-derived one:

- **Merge pairs are different** (different nodes merge at each step)
- **Merge heights are different** (different ultrametric distances)
- **Community overlap at n=5 is 4-15%** (essentially unrelated structures)

ImCoh is NOT a rescaling of MSC. The hierarchical structure is fundamentally
different because ImCoh removes the dominant volume-conduction signal that
drives the MSC hierarchy.

### The `compute_lrg_analysis()` caching bug

The function `compute_lrg_analysis()` uses `fc_method` to construct cache
paths. When called with `fc_method="msc"` (the default) while passing an
ImCoh matrix, it matches and returns the existing MSC cache instead of
computing fresh from the ImCoh matrix. This was initially misinterpreted
as "LRG is invariant to MSC vs ImCoh."

**How to avoid:**

1. Use `fc_method="imcoh"` for proper cache separation
2. Use `use_cache=False` to force recomputation
3. Or bypass the function entirely and use raw lrgsglib functions:
   ```python
   from lrgsglib.core import compute_normalized_linkage, compute_ultrametric_distance

   dist = compute_ultrametric_distance(imcoh_matrix)
   linkage = compute_normalized_linkage(dist)
   ```

The correct ImCoh LRG results are stored in `data/imcoh_lrg_cache/`.

### Hypothesis results: ImCoh strengthens main findings

| Hypothesis | MSC cells | ImCoh cells | Change |
|------------|-----------|-------------|--------|
| H1 (task stability) | 112 | 193 | +72% |
| H2a (task trace) | 33 | 107 | +224% |
| H2b (task approach) | 31 | 19 | -39%, survives ONLY in beta |
| H3 (within < cross) | 32 | 38 | +19% |

### Band-specific H2a findings

H2a (task trace in resting state) is strongest in specific bands:

| Band | Unanimous+ cells | k range | All patients positive? |
|------|-------------------|---------|----------------------|
| beta | 36 | k=16-58 | Yes (strongest) |
| low_gamma | 31 | k=21-58 | Yes |
| alpha | 29 | k=17-58 | No (Pat_02 barely negative, -0.02) |
| delta | 6 | — | No (Pat_02 negative) |
| theta | 5 | — | No (Pat_02, Pat_08 negative) |
| high_gamma | 0 | — | No signal in ANY hypothesis |

H2b (task approach) survives only in beta (16 cells, k=13-58).

High_gamma shows no signal in any hypothesis with ImCoh.

### ImCoh results files

- ImCoh FC matrices: `data/imcoh_cache/{patient}/{band}_{phase}_imcoh_*.npy`
- ImCoh LRG: `data/imcoh_lrg_cache/{patient}/{band}_{phase}_lrg_imcoh.npz`
- VI profiles: `data/imcoh_vi/vi_raw_profiles.csv`, `hypothesis_contrasts.csv`
- Unanimity maps: `data/imcoh_unanimity/unanimity_map_H*.pdf`
- Scalar contrasts: `data/imcoh_unanimity/scalar_contrasts.csv`
- Summary: `data/imcoh_unanimity/unanimity_summary.md`

## Surrogates and ImCoh

Circular-shift surrogates are **WRONG** for ImCoh. Circular shifts introduce
artificial phase lags into the surrogate data, which ImCoh detects as
"genuine" coupling. This means the surrogate null is systematically too high,
making sparsification overly conservative and potentially destroying real
signal.

**Recommendation:** Use dense ImCoh (sparsify="none") for LRG analysis.
If sparsification is ever needed, use the disparity filter or ECM (no
surrogates).

## References

1. **Nolte G, Bai O, Wheaton L, Mari Z, Vorbach S, Hallett M.** (2004)
   "Identifying true brain interaction from EEG data using the imaginary
   part of coherency." *Clinical Neurophysiology*, 115(10):2292-2307.
   - The original ImCoh paper. Over 3000 citations.
   - Key result: Im(coherency) is exactly zero for volume conduction.

2. **Vinck M, Oostenveld R, van Wingerden M, Battaglia F, Pennartz CMA.** (2011)
   "An improved index of phase-synchronization for electrophysiological data
   in the presence of volume-conduction, noise and sample-size bias."
   *NeuroImage*, 55(4):1548-1565.
   - Introduces wPLI as an improvement over PLI with better noise properties.

3. **Bastos AM, Schoffelen JM.** (2016)
   "A tutorial review of functional connectivity analysis methods and their
   interpretational pitfalls." *Frontiers in Systems Neuroscience*, 9:175.
   - Comprehensive review comparing FC metrics including MSC, ImCoh, PLI, wPLI.
   - Clear recommendation to use ImCoh or wPLI for electrophysiology data.

## Related Files

- Implementation: `src/lrg_eegfc/utils/fc/msc/msc.py` (the `metric` parameter)
- Pipeline: `src/lrg_eegfc/utils/fc/msc/__init__.py` (`coherence_fc_pipeline`)
- Probe bias guide: `.agents/guides/probe-bias-guide.md`
- MSC method guide: `.agents/guides/msc-method-guide.md`
- Process report: `.agents/guides/2026-04-15_imcoh-process.md`
