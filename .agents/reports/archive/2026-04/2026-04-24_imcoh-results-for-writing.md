---
name: imcoh-results-for-writing
type: report
era: IMCOH_ABS
status: superseded
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

> **⚠ SUPERSEDED (2026-04-24)** — This document was the writing-agent
> handoff for the n=5 / n=7 / pre-cohort-expansion phase of the project.
> It has been **superseded** by
> [`2026-04-24_multiscale-task-trace.md`](2026-04-24_multiscale-task-trace.md),
> which covers the current n=9/n=10 cohort with the reorganized H2
> hypothesis framework (H2c universal drift direction, H2d band-selective
> block persistence, H2a demoted). Numbers in this file also predate the
> 2026-04-15 ImCoh reset and reflect the pre-reset `|ImCoh|²` mislabel.
> **Do not quote numbers from this file.** Kept for historical continuity
> of the MSC → ImCoh transition narrative.

# ImCoh vs MSC Results — Comprehensive Comparison for Writing Agent
## [SUPERSEDED — see 2026-04-24_multiscale-task-trace.md]

This document provides the complete, verified results comparing the original
MSC-based analysis with the ImCoh-based reanalysis. Use this to decide which
figures to reproduce, which claims to update, and which narrative to adopt.

---

## 1. What Changed and Why

**Problem:** MSC (magnitude-squared coherence) captures all coupling including
instantaneous (zero-phase-lag) components. Volume conduction and common
reference noise produce zero-lag correlations, inflating same-probe MSC by
2-8× and contaminating the hierarchical community structure.

**Fix:** ImCoh (imaginary part of coherency, Nolte 2004) uses only the
imaginary part of the cross-spectral density. Volume conduction is
instantaneous → zero phase lag → purely real CSD → Im(CSD) = 0.
ImCoh is mathematically guaranteed to be immune.

**Implementation:** `metric="imcoh"` parameter in `compute_msc_welch()`.
Same computation time (~35s per patient/phase). All channels preserved.

**Verification:** Same-probe/cross-probe MSC ratio drops from 2-6× to 0.6-2×
across all 5 patients. The LRG dendrograms are fundamentally different
(0/116 merge pairs identical, 4-15% community overlap at n=5).

---

## 2. Unanimity Map Comparison (Full Hierarchy k=2 to N-1)

### Raw cell counts

```
                     --- MSC (k=2..61) ---     |     --- ImCoh (k=2..N-1) ---
band               H1   H2a  H2b+  H2b-    H3 |    H1   H2a  H2b+  H2b-    H3
---------------------------------------------------------------------------
delta              17     1     1     0     6 |    72    16     5     1     6
theta               4     3     2     0     1 |   110     5     0     9    46
alpha              17    29    17     0     1 |    75    32     3     3    11
beta               41     1     0     0    14 |    82    82    21     0     0
low_gamma          12     0    13     0     5 |    73    59     8     1    20
high_gamma         21     1     0     0     4 |    32     0     0     0    28
---------------------------------------------------------------------------
TOTAL             112    35    33     0    31 |   444   194    37    14   111
```

**NOTE:** MSC was computed to k=N/2 only; ImCoh goes to k=N-1. Even
accounting for the extended range, the changes in dominant bands are real
(verified by checking the k=2..N/2 subset of ImCoh separately).

### ImCoh at k=2..N/2 only (for fair comparison with MSC)

```
              H1    H2a   H2b+  H2b-   H3
delta         35      6     0     1     0
theta         54      5     0     9     16
alpha         24     29     3     2     3
beta          32     36    16     0     0
low_gamma     48     31     0     1     18
high_gamma     0      0     0     0     1
TOTAL        193    107    19    13     38
```

---

## 3. What Survived, What Changed, What Died

### H1 — Task Stability: VI(TL,TT) < mean(other pairs)

| Status | Detail |
|--------|--------|
| **SURVIVED and STRENGTHENED** | 112 → 444 cells (at full k) |
| MSC dominant band | **beta** (41 cells) |
| ImCoh dominant band | **theta** (110 cells) |
| Narrative change | **YES — H1 shifts from beta to theta** |
| All bands positive? | YES (delta through high_gamma all have cells) |

**Key finding:** Task stability is now primarily a **theta-band** phenomenon
(92% of all k-levels are unanimous). Beta is second (68%). This is more
consistent with the neuroscience literature (theta oscillations are linked
to working memory maintenance and cognitive stability).

### H2a — Task Trace: VI(Pre,Post) > VI(TT,Post)

| Status | Detail |
|--------|--------|
| **SURVIVED and DRAMATICALLY STRENGTHENED** | 35 → 194 cells |
| MSC dominant band | **alpha** (29 cells, only significant band) |
| ImCoh dominant band | **beta** (82 cells, 68% coverage) |
| Second band | **low_gamma** (59 cells, 49% coverage) |
| Alpha still present? | Yes (32 cells), but no longer dominant |
| Narrative change | **YES — H2a shifts from alpha to beta/low_gamma** |
| All 5 patients positive? | beta: YES, low_gamma: YES, alpha: NO (Pat_02 = -0.02) |

**Band-specific detail (ImCoh, full k=2..N-1):**

| Band | Unanimous+ cells | k-range | Coverage | All 5 patients + mean? |
|------|----------------:|---------|----------|:----------------------:|
| **beta** | **82** | 16-116 | **68%** | **YES** |
| **low_gamma** | **59** | 21-116 | **49%** | **YES** |
| alpha | 32 | 17-116 | 27% | NO (Pat_02 = -0.02) |
| delta | 16 | 8-116 | 13% | NO (Pat_02 = -0.10) |
| theta | 5 | 45-49 | 4% | NO (Pat_02 = -0.16) |
| high_gamma | 0 | — | 0% | NO |

### H2b — Task Approach: VI(Pre,TT) > VI(TT,Post)

| Status | Detail |
|--------|--------|
| **WEAKENED but survives in beta** | 33 → 37 cells (net) |
| MSC dominant band | alpha (17) + low_gamma (13) |
| ImCoh dominant band | **beta only** (21 cells) |
| Narrative change | **YES — alpha and low_gamma H2b disappear** |
| Theta REVERSAL | 9 cells unanimous NEGATIVE |

**Band-specific detail (ImCoh, full k=2..N-1):**

| Band | Unanimous+ | Unanimous− | Interpretation |
|------|----------:|----------:|----------------|
| **beta** | **21** | 0 | Survives, clean |
| low_gamma | 8 | 1 | Weak positive |
| delta | 5 | 1 | Weak positive |
| alpha | 3 | 3 | Mixed, lost |
| **theta** | **0** | **9** | **REVERSED** |
| high_gamma | 0 | 0 | Nothing |

### H3 — Within < Cross: mean(VI_cross) > mean(VI_within)

| Status | Detail |
|--------|--------|
| **SURVIVED and STRENGTHENED** | 31 → 111 cells |
| MSC dominant band | beta (14) |
| ImCoh dominant band | **theta** (46) + **high_gamma** (28) + **low_gamma** (20) |
| Narrative change | **YES — broadens from beta-only to theta-dominated** |

---

## 4. The Old Narrative vs The New Narrative

### MSC narrative (now superseded)

> "Alpha and beta show a dissociation: alpha carries the task trace (H2a)
> while beta provides task stability (H1). The task approach effect (H2b)
> is present in alpha and low gamma."

### ImCoh narrative (correct, volume-conduction-free)

> "**Beta is the universal task-reorganization band**, carrying both the
> task trace (H2a, 82 cells, all patients positive) and task approach
> (H2b, 21 cells, only band that survives). **Theta is the dominant
> task-stability band** (H1, 110 cells, 92% of all k-levels). The
> alpha-band task trace effect (29-32 cells) is real but secondary to
> beta. The apparent alpha dominance in the MSC analysis was inflated
> by volume conduction, which preferentially affects alpha-band coherence
> due to the strong posterior alpha rhythm generating zero-lag
> correlations across nearby contacts."

### Why the shift makes neuroscience sense

1. **Beta for task reorganization:** Beta oscillations (13-30 Hz) are
   associated with motor planning, decision-making, and top-down control.
   Task-induced reorganization of beta connectivity is well-documented
   in the motor and cognitive literature.

2. **Theta for cognitive stability:** Theta oscillations (4-8 Hz) are
   linked to working memory, hippocampal-cortical communication, and
   sustained attention. Theta coherence stability during task execution
   is a hallmark of successful cognitive engagement.

3. **Alpha inflation by volume conduction:** The posterior alpha rhythm
   (8-13 Hz) is the strongest oscillatory signal in the brain and
   generates large electric fields. Volume conduction of alpha is
   therefore disproportionately strong compared to other bands, inflating
   alpha MSC between nearby contacts. ImCoh removes this inflation.

4. **High gamma artifact:** The MSC analysis showed 21 H1 cells in
   high_gamma. In ImCoh, high_gamma H1 persists (32 cells, expanded with
   full k range) but H2a/H2b are zero — the high_gamma task trace was
   pure artifact.

---

## 5. Consistency Checks Performed

| Check | Result |
|-------|--------|
| Same-probe MSC ratio (MSC) | 2-6× across patients |
| Same-probe MSC ratio (ImCoh) | 0.6-2× (bias removed) |
| Bipolar re-referencing tested? | Yes — MAKES BIAS WORSE (5× → 8-17×) |
| Percentile rescaling tested? | Yes — ad hoc, non-physical, rejected |
| Zero same-probe edges tested? | Yes — ad hoc, loses information, rejected |
| LRG dendrograms different? | YES — 120/120 cases different, 4-15% community overlap |
| Full k-range tested? | YES — k=2 to N-1 (~116 levels) |
| Surrogates needed? | NO — circular-shift surrogates are WRONG for ImCoh |
| All 5 patients included? | YES (Pat_02, Pat_03, Pat_05, Pat_07, Pat_08) |
| All 4 phases included? | YES (rsPre, taskLearn, taskTest, rsPost) |
| All 6 bands tested? | YES (delta through high_gamma) |
| Pat_03 (1024 Hz) handled? | YES — nperseg=2048 |

---

## 6. Figures That Need Reproduction

### Must reproduce (results changed)

1. **Unanimity maps** — H1, H2a, H2b, H3 — all look different now
2. **Radar charts** — scalar contrasts per band are different
3. **Mean contrast heatmap** (4-panel: H1, H2a, H2b, H3)
4. **VI(k) profiles by band** — the per-band curves are different
5. **Adjacency matrix heatmaps** — show MSC vs ImCoh side-by-side
6. **Eigenvalue spectra** — Laplacian eigenvalues differ
7. **Entropy/susceptibility curves** — S(τ) and C(τ) from ImCoh LRG
8. **Brain connectome plots** — community structure is different

### Can keep (methodology unchanged)

- Description of the LRG method itself
- Description of the VI metric and unanimity criterion
- Patient demographics and recording details
- Electrode implantation descriptions

### Must add (new content)

- Section on volume conduction and ImCoh justification
- Same-probe ratio comparison table (MSC vs ImCoh)
- Discussion of the alpha→beta narrative shift
- References: Nolte 2004, Vinck 2011, Bastos & Schoffelen 2016

---

## 7. Data File Locations

| Data | Path |
|------|------|
| ImCoh FC matrices | `data/imcoh_cache/{patient}/{band}_{phase}_imcoh_*.npy` |
| ImCoh LRG results | `data/imcoh_lrg_cache/{patient}/{band}_{phase}_lrg_imcoh.npz` |
| ImCoh VI profiles (full) | `data/imcoh_vi/vi_raw_profiles_full.csv` |
| ImCoh contrasts (full) | `data/imcoh_vi/hypothesis_contrasts_full.csv` |
| ImCoh unanimity maps | `data/imcoh_unanimity/unanimity_map_H*.pdf` |
| ImCoh scalar contrasts | `data/imcoh_unanimity/scalar_contrasts.csv` |
| MSC VI profiles (original) | `data/wp0_metric_exploration/task3_multiscale/vi_raw_profiles.csv` |
| MSC LRG results | `data/lrg_cache/{patient}/{band}_{phase}_lrg_msc.npz` |
| Visualization gallery | `data/imcoh_figures/gallery/` |
| Process report | `.agents/guides/2026-04-15_imcoh-process.md` |

### Each .npz LRG file contains:
```
ultrametric_matrix, linkage_matrix, entropy_tau, entropy_1_minus_S,
entropy_C, optimal_threshold, patient, phase, band, fc_method, n_nodes
```

### VI profiles CSV columns:
```
patient, band, k, pair, vi
```

### Hypothesis contrasts CSV columns:
```
hypothesis, band, k, patient, contrast, sign
```

---

## 8. Code Entry Points

```python
# Compute ImCoh FC
from lrg_eegfc.utils.fc.msc import coherence_fc_pipeline
bands = coherence_fc_pipeline(ts, fs, metric="imcoh", sparsify="none", n_surrogates=0)

# Load cached ImCoh matrix
A = np.load(f"data/imcoh_cache/{patient}/{band}_{phase}_imcoh_sparsify-none_nperseg-{nperseg}.npy")

# Compute LRG on ImCoh (BYPASS cache bug)
from lrgsglib.core import entropy, compute_laplacian_properties, compute_normalized_linkage
# ... use raw functions, NOT compute_lrg_analysis() which has a caching bug
# OR use compute_lrg_analysis(..., fc_method="imcoh", use_cache=False)

# Load cached ImCoh LRG
data = np.load(f"data/imcoh_lrg_cache/{patient}/{band}_{phase}_lrg_imcoh.npz")
Z = data['linkage_matrix']
```

---

## 9. Open Questions for the Writing Agent

1. **Should the paper present ONLY ImCoh results?** Or MSC + ImCoh comparison?
   Recommendation: ImCoh as primary, MSC as supplementary showing "before correction."

2. **How to frame the alpha→beta shift?** As a correction of a volume-conduction
   artifact, or as a new finding? Recommendation: frame as "volume conduction
   preferentially inflates alpha-band coherence, masking the true beta-dominant
   task reorganization signal."

3. **The k-range extension (k=2..N-1 vs k=2..N/2):** The original MSC analysis
   stopped at N/2. Should the ImCoh analysis also stop at N/2 for comparability,
   or use the full range? Recommendation: show both — N/2 for fair comparison,
   full range as the complete picture.

4. **Pat_02 in H2a alpha:** Pat_02 has a mean contrast of -0.02 in alpha H2a,
   breaking unanimity. Is this worth discussing? The other 4 patients are
   positive. For beta and low_gamma, all 5 patients are positive.

5. **Theta H2b reversal:** 9 cells of unanimous NEGATIVE in theta H2b. This means
   in theta, the task-test state is FARTHER from rsPost than rsPre is. This is
   the opposite of the approach hypothesis. Worth discussing as a band-specific
   dissociation.
