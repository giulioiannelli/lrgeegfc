---
name: writing-agent-briefing
type: report
era: IMCOH_SQ
status: superseded
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

> **⚠ QUANTITATIVE_STALE (2026-04-15)** — every numerical value in this
> document was computed under the pre-reset mislabelling, in which
> `fc_method="imcoh"` actually stored `|ImCoh|²`. After the 2026-04-15
> reset the three ImCoh quantities are correctly defined: `imcoh`
> (signed Nolte 2004), `imcoh_abs` (`|ImCoh|`, Ewald 2012), `imcoh_sq`
> (`|ImCoh|²`, Ewald 2012). Qualitative conclusions (MSC-vs-ImCoh
> contrast, probe-bias reduction, community structure, unanimity counts,
> enrichment ranks) are preserved because they depend on rankings, which
> are invariant under the monotonic `sqrt` transform linking
> `|ImCoh|` and `|ImCoh|²`. **Numerical values in this file must be
> regenerated from `imcoh_abs` before quoting.** See
> `~/.claude/projects/.../memory/imcoh_taxonomy.md`.

# Writing Agent Briefing: MSC → ImCoh Transition

This document summarises the shift from Magnitude-Squared Coherence (MSC) to
Imaginary Coherence (ImCoh) as the primary functional connectivity metric in
the LRG-sEEG pipeline. All numbers have been independently verified (see
`2026-04-15_imcoh-verification.md`). Use this as the single reference when
drafting or revising manuscript sections related to FC methodology and results.

---

## 1. The Problem with MSC

MSC (`|S_ij|² / (S_ii · S_jj)`) captures **all** coupling, including
instantaneous (zero-phase-lag) components arising from:

1. **Common reference noise.** All sEEG channels are recorded as V_i − V_G2.
   The shared reference G2 injects a common signal into every channel.
2. **Volume conduction.** Electric fields propagate effectively instantaneously
   through brain tissue. Physically adjacent contacts (~3–5 mm on the same
   electrode probe) record overlapping neural populations.

Both sources produce zero-phase-lag correlations → purely real
cross-spectral density → inflated MSC.

### Quantified impact

- Same-probe MSC is **2–8× higher** than cross-probe MSC across all patients.
- LRG communities at coarse scales (n ≤ 10–30) follow electrode probe
  geometry rather than functional brain organisation.
- Same-community/same-probe enrichment reaches 3–8× at n = 10.
- All community-level metrics (participation coefficient, community cosine,
  metastability, allegiance) are confounded.

### How it was detected

During the epileptic-node investigation, 3D brain connectome plots coloured
by LRG community revealed that communities traced electrode trajectories.
Quantitative enrichment ratios confirmed the bias systematically across all
5 patients. Variation of Information (VI) between probe labels and community
labels was misleading (many-to-one mappings mask the bias); the enrichment
ratio is the correct diagnostic.

---

## 2. Failed Debiasing Approaches

| Approach | Outcome |
|----------|---------|
| Zero same-probe edges | Ad hoc. Destroys genuine within-probe local-circuit dynamics. Would not survive peer review. |
| Percentile rescaling | No physical or statistical justification. Rescaled values do not correspond to any meaningful coherence quantity. |
| Bipolar re-referencing | **Makes bias worse** (5× → 8–17×). Adjacent bipolar channels share a physical contact (V_{k+1}), creating a new, stronger source of spurious coherence. |
| Distance regression | MSC–distance relationship is highly nonlinear with a sharp drop at same-probe distances. No principled nonlinear model exists. |

---

## 3. The Solution: Imaginary Coherence

### Physics

1. Volume conduction is instantaneous.
2. Instantaneous mixing → zero phase lag.
3. Zero phase lag → Im(S_ij) = 0 for the volume-conducted component.
4. Therefore Im(S_ij)² / (S_ii · S_jj) is **exactly zero** for any
   volume-conducted signal.
5. Any nonzero ImCoh reflects genuine neural interaction with a nonzero
   time delay.

This is an exact physical consequence, not an approximation or statistical
correction.

### Key reference

Nolte G, Bai O, Wheaton L, Mari Z, Vorbach S, Hallett M (2004). Identifying
true brain interaction from EEG data using the imaginary part of coherency.
*Clinical Neurophysiology*, 115(10), 2292–2307. (>3 000 citations)

Additional references: Vinck et al. 2011 (wPLI), Bastos & Schoffelen 2016
(review of FC metrics and volume conduction).

### Same-probe ratio after ImCoh

| Metric | Same-probe ratio range | Typical |
|--------|:----------------------:|:-------:|
| MSC    | 1.6–6.0×               | 3–5×    |
| ImCoh  | 0.6–2.0×               | ~1×     |

Residual ImCoh ratio (e.g. ~2× in beta for Pat_07) reflects genuine lagged
short-range coupling, not artifact.

### Technical notes

- Same Welch CSD estimation; only the final derivation step changes
  (`Im(S_ij)²` instead of `|S_ij|²`). Computation time identical (~35 s per
  patient/phase).
- Circular-shift surrogates are **wrong** for ImCoh: they scramble phase
  structure, creating artificial lags that ImCoh detects as "real" coupling.
  The null distribution is systematically inflated. Use **dense ImCoh**
  (no surrogate-based sparsification).
- The LRG dendrograms are fundamentally different: 120/120 patient–phase–band
  cases produce different merge sequences, with only 4–15% community overlap
  at n = 5.

---

## 4. Hypothesis Testing: MSC vs ImCoh

The pipeline tests four hypotheses using Variation of Information (VI) between
LRG community partitions at each hierarchical scale k. "Unanimous cells" are
(band, k) positions where all 5 patients agree on the sign of the contrast.

### 4.1 Cell counts

**Full k-range (k = 2 to N − 1, ~116 levels):**

| Hypothesis | Description | MSC | ImCoh | Change |
|------------|-------------|:---:|:-----:|:------:|
| H1  | Task stability: VI(TL,TT) < mean(other pairs) | 112 | 444 | +296% |
| H2a | Task trace: VI(Pre,Post) > VI(TT,Post) | 35 | 194 | +454% |
| H2b+| Task approach: VI(Pre,TT) > VI(TT,Post) | 33 | 37 | +12% |
| H2b−| Reversed approach | 0 | 14 | new |
| H3  | Within < cross-condition VI | 31 | 111 | +258% |

**Fair comparison at k = 2 to N/2 only (matching MSC range):**

| Hypothesis | MSC | ImCoh (N/2) | Change |
|------------|:---:|:-----------:|:------:|
| H1  | 112 | 193 | +72% |
| H2a | 35  | 107 | +206% |
| H2b+| 33  | 19  | −42% |
| H3  | 31  | 38  | +23% |

Even restricted to the same k-range, ImCoh strengthens H1, H2a, and H3.
H2b narrows (survives only in beta).

### 4.2 Band-level breakdown (ImCoh, full k-range)

**H1 — Task stability:**

| Band | Cells | Dominant? |
|------|:-----:|:---------:|
| theta | 110 | **YES** (92% of k-levels) |
| beta | 82 | second |
| alpha | 75 | |
| low_gamma | 73 | |
| delta | 72 | |
| high_gamma | 32 | fine-scale only (k > 84) |

**H2a — Task trace:**

| Band | Cells | k-range | All 5 patients positive? |
|------|:-----:|---------|:------------------------:|
| **beta** | **82** | 16–108 | **YES** |
| **low_gamma** | **59** | 21–108 | **YES** |
| alpha | 32 | 17–61 | NO (Pat_02 = −0.05) |
| delta | 16 | 8–71 | NO |
| theta | 5 | 45–49 | NO |
| high_gamma | 0 | — | NO |

**H2b — Task approach:**

| Band | Unanimous + | Unanimous − | Note |
|------|:-----------:|:-----------:|------|
| **beta** | **21** | 0 | Only clean survivor |
| low_gamma | 8 | 1 | Weak |
| delta | 5 | 1 | Weak |
| alpha | 3 | 3 | Mixed → lost |
| **theta** | **0** | **9** | **REVERSED** (k = 28–36) |
| high_gamma | 0 | 0 | — |

**H3 — Within < cross:**

| Band | Cells | Dominant? |
|------|:-----:|:---------:|
| **theta** | **46** | **YES** |
| high_gamma | 28 | |
| low_gamma | 20 | |
| alpha | 11 | |
| delta | 6 | |
| beta | 0 | (Pat_07 and Pat_08 negative → unanimity impossible) |

### 4.3 Patient-level H2a beta contrasts

| Patient | Mean contrast | Frac positive k-levels |
|---------|:------------:|:----------------------:|
| Pat_02 | +0.262 | 96% |
| Pat_03 | +0.259 | 83% |
| Pat_05 | +0.187 | 93% |
| Pat_07 | +0.192 | 88% |
| Pat_08 | +0.182 | 90% |

Group mean: **+0.217 ± 0.036** (all positive, robust).

Reference VI distances (beta, averaged across k):
- VI(Pre, Post) ≈ 1.1 (largest: resting states most different)
- VI(TT, Post) ≈ 0.9 (task-test closer to rsPost)
- VI(TL, TT) ≈ 0.7 (smallest: task phases most similar = H1)
- Δ_H2a ≈ 0.2 → ~20% of VI(Pre, TT), a meaningful effect size

---

## 5. The Narrative Shift

### MSC narrative (superseded)

> Alpha and beta show a dissociation: alpha carries the task trace (H2a)
> while beta provides task stability (H1). The task-approach effect (H2b) is
> present in alpha and low gamma.

### ImCoh narrative (correct, volume-conduction-free)

> **Beta is the universal task-reorganisation band**, carrying both the task
> trace (H2a: 82 unanimous cells, all 5 patients positive) and task approach
> (H2b: 21 cells, only surviving band). **Theta is the dominant
> task-stability band** (H1: 110 cells, 92% of k-levels) and dominates
> rest/task discrimination (H3: 46 cells). The alpha-band task trace (32
> cells) is real but secondary; its apparent dominance under MSC was an
> artifact of volume conduction, which preferentially inflates alpha-band
> coherence due to the strong posterior alpha rhythm generating large
> zero-lag electric fields across nearby contacts.

### Theta–beta dissociation

|  | Theta | Beta |
|--|:-----:|:----:|
| H1 (stability) | **110** | 82 |
| H2a (trace) | 5 | **82** |
| H2b (approach) | 0 / **9 reversed** | **21** / 0 |
| H3 (within < cross) | **46** | 0 |

- **Theta** = stability + discrimination (H1 + H3), negligible trace
- **Beta** = trace + approach (H2a + H2b), strong stability but no discrimination

The dissociation is clean for H2a, H2b, and H3. Both bands contribute to H1.

Interpretation: beta connectivity reorganises during task performance (H2a),
which disrupts rest/task discrimination in beta (H3 = 0). Theta connectivity
remains stable (H2a ≈ 0), preserving rest/task discrimination (H3 = 46).

### Why this makes neuroscience sense

1. **Beta (13–30 Hz)** is associated with motor planning, decision-making,
   and top-down control. Task-induced reorganisation of beta connectivity is
   well-documented in the motor and cognitive literature.
2. **Theta (4–8 Hz)** is linked to working memory, hippocampal–cortical
   communication, and sustained attention. Theta coherence stability during
   task execution is a hallmark of successful cognitive engagement.
3. **Alpha inflation by volume conduction:** the posterior alpha rhythm
   (8–13 Hz) is the strongest oscillatory signal in the brain and generates
   the largest volume-conducted fields. MSC alpha was therefore
   disproportionately inflated relative to other bands.
4. **High-gamma task trace was artifact:** MSC showed 21 H1 cells in
   high_gamma; ImCoh shows 0 H2a/H2b cells. The high-gamma "signal" was
   volume conduction.

---

## 6. Probe Bias Under ImCoh

Same-community/same-probe enrichment ratio by scale:

| Patient | n = 3 | n = 5 | n = 10 | n = 15 | n = 20 | n = 30 |
|---------|:-----:|:-----:|:------:|:------:|:------:|:------:|
| Pat_02 | 1.0 | 1.2 | 1.3 | 1.7 | 1.7 | 2.5 |
| Pat_03 | 1.0 | 1.0 | 1.1 | 1.1 | 1.2 | 1.3 |
| Pat_05 | 1.1 | 1.1 | 1.2 | 1.5 | 1.5 | 3.0 |
| Pat_07 | 1.0 | 1.0 | 1.1 | 1.5 | 1.8 | 2.4 |
| Pat_08 | 1.0 | 1.2 | 2.4 | 2.8 | 3.2 | 3.8 |

ImCoh dramatically reduces probe bias at coarse scales (n ≤ 10) compared to
MSC (3–8× at n = 10). Residual bias at fine scales (n ≥ 15) reflects genuine
lagged short-range coupling — ImCoh guarantees it is not volume conduction.

---

## 7. Flags and Caveats

1. **H2a beta cells are not contiguous.** The 82 cells span k = 16–108 with
   gaps where 1–2 patients dip below zero. Report total count but note
   non-contiguity.
2. **Probe bias persists at fine scales** (n ≥ 15) for some patients
   (Pat_08: 2.4× at n = 10). This is genuine lagged coupling, not artifact.
3. **Pat_02 is the weakest patient for H2a** — negative mean contrast in
   alpha, delta, and theta; positive only in beta and low_gamma.
4. **Beta H1 = 82** means the theta–beta dissociation is imperfect for H1
   (both bands show strong task stability). The dissociation is clean only
   for H2a/H2b/H3.
5. **High-gamma H1 = 32** appears only at fine scales (k > 84). Likely noise
   amplification near the ImCoh noise floor (~0.001). Report cautiously or
   note as fine-scale-only.
6. **Theta H2b reversal (9 cells, k = 28–36):** in theta, task-test is
   *farther* from rsPost than rsPre is — the opposite of approach. All 5
   patients negative, contrasts −0.04 to −0.39 (genuine, not marginal).
   Worth discussing as a band-specific dissociation.
7. **Pat_03 (negative control):** recorded at 1024 Hz (others at 2048 Hz),
   MSC 3× higher than any other patient. ImCoh normalises Pat_03's behaviour.
   Excluding Pat_03 from H2a beta → 84 cells with 4/4 unanimous (stronger).

---

## 8. Manuscript Recommendations

### Framing

- Present **ImCoh as the primary analysis**. Show MSC results as supplementary
  material demonstrating the "before correction" state.
- Frame the alpha → beta shift as: "volume conduction preferentially inflates
  alpha-band coherence, masking the true beta-dominant task-reorganisation
  signal."
- The theta–beta dissociation is the central finding. It is more
  interpretable than the previous alpha–beta narrative and better aligned
  with the established literature.

### k-range

- Show both N/2 (for fair comparison with MSC) and full k-range (complete
  picture). The qualitative conclusions are the same; the full range
  amplifies them.

### Sections to add

1. **Volume conduction and ImCoh justification.** Explain the physics
   (zero-lag → purely real CSD → Im = 0). Cite Nolte 2004.
2. **Same-probe ratio comparison table** (MSC vs ImCoh).
3. **Discussion of the alpha → beta narrative shift.**
4. **Surrogate note:** explain why circular-shift surrogates are wrong for
   ImCoh and why dense (unsparsified) ImCoh is used.

### Sections to keep unchanged

- LRG method description
- VI metric and unanimity criterion definition
- Patient demographics and recording details
- Electrode implantation descriptions

### Figures to reproduce with ImCoh

1. Unanimity maps (H1, H2a, H2b, H3)
2. Scalar contrast radar charts / heatmaps
3. VI(k) profiles by band
4. Adjacency matrix heatmaps (MSC vs ImCoh side-by-side)
5. Entropy/susceptibility curves S(τ), C(τ)
6. Brain connectome plots (community structure)

---

## 9. Data Locations

| Data | Path |
|------|------|
| ImCoh FC matrices | `data/cache/imcoh/{patient}/{band}_{phase}_imcoh_*.npy` |
| ImCoh LRG results | `data/cache/imcoh_lrg/{patient}/{band}_{phase}_lrg_imcoh.npz` |
| ImCoh VI profiles | `data/reports/imcoh_vi/vi_raw_profiles_full.csv` |
| ImCoh contrasts | `data/reports/imcoh_vi/hypothesis_contrasts_full.csv` |
| Unanimity maps (PDF) | `data/reports/imcoh_unanimity/unanimity_map_H*.pdf` |
| Scalar contrasts | `data/reports/imcoh_unanimity/scalar_contrasts.csv` |
| MSC VI profiles | `data/reports/metric_exploration/task3_multiscale/vi_raw_profiles.csv` |
| MSC LRG results | `data/cache/lrg/{patient}/{band}_{phase}_lrg_msc.npz` |

### LRG .npz contents

```
ultrametric_matrix, linkage_matrix, entropy_tau, entropy_1_minus_S,
entropy_C, optimal_threshold, patient, phase, band, fc_method, n_nodes
```

### VI profiles CSV columns

```
patient, band, k, pair, vi
```

### Hypothesis contrasts CSV columns

```
hypothesis, band, k, patient, contrast, sign
```

---

## 10. Related Documentation

| Document | Purpose |
|----------|---------|
| `.agents/reports/2026-04-15_imcoh-process.md` | Full discovery journey with technical details |
| `.agents/reports/2026-04-24_imcoh-results-for-writing.md` | Detailed MSC vs ImCoh comparison tables |
| `.agents/reports/2026-04-15_imcoh-verification.md` | Independent verification of all numbers |
| `.agents/guides/msc-method-guide.md` | MSC methodology reference |
| `.agents/guides/probe-bias-guide.md` | Same-probe bias analysis |
