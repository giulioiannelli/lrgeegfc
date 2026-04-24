#!/usr/bin/env python3
"""Curate report figures and generate technical report.

Selects 20 key figures from data/figures/metric_exploration/ subfolders,
copies them with meaningful names into report_figures/, and generates
a comprehensive TECHNICAL_REPORT.md for a writing agent.

Usage:
    python scripts/py/curate_report_figures.py [--dry-run]
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.paths import FIGURES_ROOT
SRC = FIGURES_ROOT / "metric_exploration"
DST = FIGURES_ROOT / "metric_exploration" / "report_figures"

# ═══════════════════════════════════════════════════════════════════════
# FIGURE MANIFEST — (target_filename, source_relative_to_SRC)
# ═══════════════════════════════════════════════════════════════════════
FIGURES = [
    # ── Section A: Overview & Methods ─────────────────────────────────
    ("A01_all_pairs_nvi_per_band.pdf",
     "continuous_all_pairs/all_pairs_per_band.pdf"),
    ("A02_k_at_h_profile.pdf",
     "continuous_all_pairs/k_at_h_profile.pdf"),
    ("A03_nvi_heatmaps_per_pair.pdf",
     "continuous_all_pairs/heatmaps_per_pair.pdf"),

    # ── Section B: H1 — Task Stability ───────────────────────────────
    ("B01_H1_task_stability_normalized.pdf",
     "final_hypotheses/H1_task_stability.pdf"),
    ("B02_H1_task_stability_multiscale.pdf",
     "multiscale_hypotheses/H1_task_stability.pdf"),

    # ── Section C: H2 — Task Trace ───────────────────────────────────
    ("C01_H2_definitive_heatmap.pdf",
     "continuous_h2/definitive_h2_heatmap.pdf"),
    ("C02_H2_per_band_detail.pdf",
     "continuous_h2/h2a_per_band_detail.pdf"),
    ("C03_H2_vi_raw_all_pairs.pdf",
     "continuous_h2/vi_raw_all_pairs.pdf"),
    ("C04_H2_task_trace_normalized.pdf",
     "final_hypotheses/H2_task_trace.pdf"),
    ("C05_H2_scale_specific_unanimity.pdf",
     "scale_specific_h2/h2a_unanimity_heatmap.pdf"),

    # ── Section D: H3 — Within vs Cross ──────────────────────────────
    ("D01_H3_within_vs_cross.pdf",
     "final_hypotheses/H3_within_vs_cross.pdf"),
    ("D02_H3_multiscale_contrasts.pdf",
     "continuous_all_pairs/contrast_heatmaps.pdf"),

    # ── Section E: H4 — Frequency Gradient ───────────────────────────
    ("E01_H4_frequency_gradient.pdf",
     "final_hypotheses/H4_frequency_gradient.pdf"),
    ("E02_H4_ranking_agreement.pdf",
     "h4_gradient/h4_ranking_agreement.pdf"),

    # ── Section F: Outlier Diagnostics ───────────────────────────────
    ("F01_outlier_Pat02_theta.pdf",
     "outlier_diagnostics/outlier_Pat_02_theta.pdf"),
    ("F02_outlier_Pat03_high_gamma.pdf",
     "outlier_diagnostics/outlier_Pat_03_high_gamma.pdf"),
    ("F03_outlier_Pat08_delta.pdf",
     "outlier_diagnostics/outlier_Pat_08_delta.pdf"),

    # ── Section G: Supplementary ─────────────────────────────────────
    ("G01_all_hypotheses_summary.pdf",
     "final_hypotheses/all_hypotheses_summary.pdf"),
    ("G02_pre_post_all_6_patients.pdf",
     "continuous_all_pairs/pre_post_all_patients.pdf"),
    ("G03_alternative_metrics_h2.pdf",
     "alternative_metrics/h2_unanimity_heatmap.pdf"),
]


# ═══════════════════════════════════════════════════════════════════════
# TECHNICAL REPORT
# ═══════════════════════════════════════════════════════════════════════
REPORT = r"""# Multiscale Reorganization of sEEG Brain Networks: Technical Report

> This document contains all definitions, methods, numerical results, and figure
> references needed to write the Results and Methods sections of the paper. Every
> symbol is defined before first use. Figure codes (e.g. **[A01]**) refer to the
> PDFs in this folder.

---

## 1. Definitions and Methodology

### 1.1 Experimental Design

Six patients (Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08) underwent
stereo-EEG (sEEG) recordings across experimental phases:

| Patient | Phases available |
|---------|-----------------|
| Pat_02  | rest_pre, task_learn, task_test, rest_post |
| Pat_03  | rest_pre, task_learn, task_test, rest_post |
| Pat_05  | rest_pre, task_learn, task_test, rest_post |
| Pat_06  | rest_pre, rest_post |
| Pat_07  | rest_pre, task_learn, rest_post |
| Pat_08  | rest_pre, task_learn, task_test, rest_post |

- **rest_pre**: resting state before task
- **task_learn**: learning phase of the task
- **task_test**: testing phase of the task
- **rest_post**: resting state after task

Hypothesis testing uses the 4 patients with all 4 phases (Pat_02, Pat_03,
Pat_05, Pat_08). Pat_06 and Pat_07 are included in the Pre-Post comparison
only.

Six frequency bands were analyzed:

| Band | Symbol | Frequency range |
|------|--------|----------------|
| delta | $\delta$ | 1–4 Hz |
| theta | $\theta$ | 4–8 Hz |
| alpha | $\alpha$ | 8–13 Hz |
| beta | $\beta$ | 13–30 Hz |
| low_gamma | $\gamma_l$ | 30–70 Hz |
| high_gamma | $\gamma_h$ | 70–150 Hz |

For each (patient, phase, band) combination, functional connectivity was
computed using Magnitude Squared Coherence (MSC, nperseg=4096), and a
hierarchical community structure was obtained via the Laplacian Renormalization
Group (LRG) method.

### 1.2 LRG Dendrograms and Cophenetic Height h

The LRG produces a dendrogram (scipy linkage matrix Z) for each (patient,
phase, band) triplet. Each dendrogram encodes a hierarchy of ~117 nodes
(sEEG contacts).

**Cophenetic height h**: the threshold on the ultrametric distance matrix at
which the dendrogram is cut. Cutting at height h yields a partition into k(h)
communities.

- h near 0: fine partition (many small communities, k large)
- h near 1: coarse partition (few large communities, k small)

The LRG diffusion distances produce a log-scaled hierarchy, meaning most merge
events occur at small h. We therefore sample h on a logarithmic grid:

$$h \in \text{geomspace}(0.003,\; 0.995,\; 400)$$

This gives uniform resolution across all hierarchical scales on a log axis.

**Reference k(h) values** (mean across all 144 trees): **[A02]**

| h    | Mean k |
|------|--------|
| 0.01 | ~98    |
| 0.05 | ~55    |
| 0.10 | ~35    |
| 0.20 | ~20    |
| 0.30 | ~12    |
| 0.50 | ~5     |
| 0.80 | ~2     |

### 1.3 Variation of Information (VI)

Given two partitions P and Q of n nodes:

$$\text{VI}(P, Q) = H(P) + H(Q) - 2\,I(P; Q)$$

where:
- $H(P) = -\sum_i p_i \ln p_i$ is the Shannon entropy (in nats) of partition P
- $p_i = |C_i| / n$ is the fraction of nodes in cluster $C_i$
- $I(P; Q) = \sum_{i,j} p_{ij} \ln \frac{p_{ij}}{p_i \cdot q_j}$ is the
  mutual information between P and Q
- $p_{ij} = |C_i^P \cap C_j^Q| / n$

**Properties:**
- True metric: symmetric, satisfies triangle inequality, VI(P,P) = 0
- Range: $\text{VI} \in [0, \ln n]$ where n is the number of nodes
- Decomposition: $\text{VI}(P,Q) = H(P|Q) + H(Q|P)$ — the information lost
  plus information gained when moving from P to Q
- VI = 0 means identical partitions; larger values mean more different

**Why VI over ARI:** The Adjusted Rand Index (ARI) is a similarity metric
bounded in $[-1, 1]$ that does not satisfy the triangle inequality and does not
scale naturally with partition size. VI is a proper distance metric with
information-theoretic interpretation. We cross-validated VI against 7
alternative metrics (cophenetic Pearson/Spearman correlation, Baker's Gamma,
normalized L1/L2 on ultrametric, top-k merge agreement, weighted ARI); all
produced consistent results for our hypotheses. **[G03]**

### 1.4 Normalized VI (NVI)

Raw VI scales with the number of clusters: fine partitions (large k) naturally
have higher entropy, leading to larger VI values. To compare across scales:

$$\text{NVI}(h) = \frac{\text{VI}(h)}{\ln(\bar{k}(h))}$$

where $\bar{k}(h)$ is the mean number of clusters at height h across all trees.
This normalization by the entropy capacity at each scale makes NVI values
comparable: NVI $\approx$ 1.0 means "as different as random partitions"
regardless of whether k = 5 or k = 50.

NVI $\in [0, \sim 2]$ in practice.

### 1.5 Signed Unanimity

To assess cross-patient agreement for each hypothesis contrast at each scale:

$$U(h) = \frac{1}{N}\sum_{i=1}^{N} \text{sign}\big(\Delta_i(h)\big)$$

where $\Delta_i(h)$ is the contrast value for patient $i$ at height $h$, and
$N$ is the number of patients with available data for that contrast.

- $U = +1$: all patients show positive contrast (unanimous support)
- $U = -1$: all patients show negative contrast (unanimous rejection)
- $U = 0$: patients split evenly

In heatmaps, $U$ is displayed as color (blue = all positive, red = all
negative). Contiguous regions where $|U| = 1$ (4/4 patients agree) are marked
with black borders.

### 1.6 Hypothesis-Specific Contrasts

All six phase pairs are denoted:

| Pair | Label | Type |
|------|-------|------|
| task_learn $\leftrightarrow$ task_test | TL-TT | within-task |
| rest_pre $\leftrightarrow$ rest_post | Pre-Post | within-rest |
| rest_pre $\leftrightarrow$ task_learn | Pre-TL | cross |
| rest_pre $\leftrightarrow$ task_test | Pre-TT | cross |
| task_learn $\leftrightarrow$ rest_post | TL-Post | cross |
| task_test $\leftrightarrow$ rest_post | TT-Post | cross |

**H1 — Relative Gain (Task Stability):**

$$G(h) = \frac{\overline{\text{VI}}_{\text{others}}(h) - \text{VI}_{\text{TL-TT}}(h)}{\overline{\text{VI}}_{\text{others}}(h)}$$

where $\overline{\text{VI}}_{\text{others}}$ is the mean VI of the 5 non-TL-TT
pairs. $G = 0.3$ means "TL-TT is 30% closer than the average pair." $G > 0$
supports H1.

**H2a — Relative Contrast (Task Persistence):**

$$C_a(h) = \frac{\text{VI}(\text{Pre,Post}; h) - \text{VI}(\text{TT,Post}; h)}{\frac{1}{2}\big(\text{VI}(\text{Pre,Post}; h) + \text{VI}(\text{TT,Post}; h)\big)}$$

Interpretation: Is rest_post more similar to task_test (small VI_TT-Post) than to
rest_pre (large VI_Pre-Post)? $C_a > 0$ means task structure **persists** in
rest_post. $C_a < 0$ means brain **recovers** to pre-task state.

**H2b — Relative Contrast (Task Proximity):**

$$C_b(h) = \frac{\text{VI}(\text{Pre,TT}; h) - \text{VI}(\text{TT,Post}; h)}{\frac{1}{2}\big(\text{VI}(\text{Pre,TT}; h) + \text{VI}(\text{TT,Post}; h)\big)}$$

Interpretation: Is rest_post closer to task_test than rest_pre was? $C_b > 0$ means
rest_post has moved toward the task configuration compared to baseline.

**H3 — Relative Gap (Within vs Cross):**

$$\text{Gap}(h) = \frac{\overline{\text{VI}}_{\text{cross}}(h) - \overline{\text{VI}}_{\text{within}}(h)}{\frac{1}{2}\big(\overline{\text{VI}}_{\text{cross}}(h) + \overline{\text{VI}}_{\text{within}}(h)\big)}$$

where within = {TL-TT, Pre-Post} and cross = {Pre-TL, Pre-TT, TL-Post,
TT-Post}. Gap > 0 means same-type phase pairs are more similar than cross-type.

**H4 — Kendall's W (Frequency Gradient):**

$$W = \frac{12S}{n^2(m^3 - m)}$$

where $S$ is the sum of squared deviations of rank sums from their mean,
$n$ = number of judges (patients), $m$ = number of items (bands). W = 0 means
no agreement on which bands reorganize most; W = 1 means perfect agreement.

### 1.7 Computational Details

- **FC method**: Magnitude Squared Coherence (MSC), nperseg = 4096
- **Sparsification**: soft thresholding based on surrogate testing (200
  surrogates)
- **LRG**: Laplacian Renormalization Group, producing scipy-format linkage
  matrices
- **VI computation**: direct from partition labels using natural logarithm
- **Height grid**: 400 log-spaced points in [0.003, 0.995]
- **Number of nodes**: ~117 sEEG contacts per patient (varies slightly)

---

## 2. Hypothesis Results

### 2.1 H1: Task Stability — SUPPORTED

**Claim:** The two task phases (task_learn, task_test) produce the most similar
hierarchical network structure, i.e., TL-TT should have the lowest VI among all
6 pairs.

**Result:** Supported on average with a mean relative gain of ~34% across all
bands and scales (TL-TT is 34% closer than the average pair). **[B01]**

**Multiscale detail** (continuous h, using signed unanimity across 4 patients):
**[B02]**

| Band | Unanimous h range | k range | Notes |
|------|------------------|---------|-------|
| $\beta$ | [0.02, 0.34] | 12–87 | Broadest continuous range |
| $\alpha$ | [0.05, 0.14] + [0.16, 0.18] + [0.20, 0.21] + [0.35, 0.40] | 10–64 | Multiple fragmented windows |
| $\theta$ | [0.11, 0.20] + [0.32, 0.34] | 12–40 | Two windows |
| $\gamma_l$ | [0.02, 0.09] + [0.15, 0.16] + [0.22, 0.28] | 16–90 | Fine + meso |
| $\gamma_h$ | [0.04, 0.05] | 63–71 | Narrow window |
| $\delta$ | — | — | No unanimous range |

**Interpretation:** Beta shows remarkably robust task stability from very fine
(k~87) to meso (k~12) scales — the task phases maintain nearly identical
hierarchical structure across a wide range of community resolutions. The effect
is consistent but fragmented for other bands.

### 2.2 H2: Task Trace — SUPPORTED (Band-Dependent, Scale-Dependent)

**Claim:** rest_post carries an imprint of the task, visible as rest_post being more
similar to task_test than to rest_pre at certain scales. The nature of this imprint
(persistence vs recovery) depends on frequency band.

**Result:** The central finding of this analysis. The H2 effect is both
**band-dependent** and **scale-dependent** — it lives primarily at fine
hierarchical scales (h < 0.3, k > 15). **[C01]**

#### H2a contrast: VI(Pre,Post) > VI(TT,Post)?

| Band | Verdict | h range | k range | # patients |
|------|---------|---------|---------|------------|
| $\alpha$ | **PERSIST** | [0.09, 0.09], [0.12, 0.14], **[0.17, 0.26]** | 16–48 | 4/4 unanimous |
| $\beta$ | **PERSIST** | **[0.02, 0.04]** | 69–87 | 4/4 unanimous |
| $\gamma_l$ | **RECOVER** | **[0.01, 0.04]**, [0.06, 0.08] | 50–98 | 4/4 unanimous |
| $\delta$ | No signal | — | — | Not unanimous |
| $\theta$ | No signal | — | — | Not unanimous |
| $\gamma_h$ | No signal | — | — | Not unanimous |

**Key interpretation:**
- **Alpha persistence** (h = 0.17–0.26, k $\approx$ 16–26): The alpha-band
  community structure established during the task remains present in rest_post.
  This is the strongest and most robust effect, occurring at the mesoscale
  (medium-sized communities). All 4 patients show the same pattern.
- **Beta persistence** (h = 0.02–0.04, k $\approx$ 69–87): Task structure
  persists at very fine scale in beta. This represents local circuit-level
  reorganization.
- **Low gamma recovery** (h = 0.01–0.08, k $\approx$ 50–98): The low gamma
  community structure in rest_post reverts to the pre-task configuration, not
  retaining the task pattern. The brain "forgets" the task in this band at fine
  scale.
- **Delta, theta, high gamma**: No unanimous signal. Patient-to-patient
  variability prevents a clear conclusion. These bands may show weak or mixed
  effects.

**Supporting figures:** Per-band patient traces **[C02]**, raw VI curves showing
the curve separations that generate the contrast **[C03]**, normalized discrete
view **[C04]**, scale-specific breakdown into macro/meso/micro **[C05]**.

### 2.3 H3: Within-Type vs Cross-Type — SUPPORTED

**Claim:** Phase pairs of the same type (rest-rest or task-task) are more
similar than cross-type pairs (rest-task).

**Result:** Supported for 5 out of 6 frequency bands. **[D01]**

Mean relative gap (mean_cross − mean_within) / mean is positive for delta,
theta, alpha, beta, and low_gamma. High gamma shows mixed results.

The multiscale contrast heatmap confirms this holds across scales, not just as
a mean effect. **[D02]**

### 2.4 H4: Frequency Gradient — REJECTED

**Claim:** There exists a consistent ordering of frequency bands by
reorganization strength (e.g., low frequencies reorganize more than high, or
vice versa).

**Result:** Rejected. **[E01]**, **[E02]**

- **Kendall's W = 0.043** (weak agreement; threshold for moderate is W > 0.3)
- Each patient has a different "top reorganizing band":
  - Pat_02: theta/delta reorganize most
  - Pat_03: theta/alpha reorganize most
  - Pat_05: high_gamma/low_gamma reorganize most
  - Pat_08: delta/beta reorganize most
- The rank heatmap **[E01]** shows shuffled colors across patients — no column
  is consistently dark or light.
- Pairwise Spearman correlations between patients are near zero or negative.

**Conclusion:** Reorganization strength is **patient-specific**, not
frequency-driven. There is no universal frequency gradient.

---

## 3. Per-Patient Findings

### 3.1 Pat_02

- **4 phases**: rest_pre, task_learn, task_test, rest_post
- **Outlier behavior**: Theta band — TL-TT pair shows mean VI = 0.877, making
  it the worst (most dissimilar) pair rather than the best. At k = 2–4 (coarse
  scale), the binary split differs completely between task_learn and task_test.
  **[F01]**
- **Explanation**: The theta dendrogram structure undergoes a genuine coarse
  reorganization between the two task phases. At fine scale (k > 10), TL-TT
  becomes more similar, consistent with other patients.
- **Contribution to hypotheses**: Supports H2 (alpha persist) and H3
  (within < cross). Partially violates H1 at coarse scale for theta.

### 3.2 Pat_03

- **4 phases**: rest_pre, task_learn, task_test, rest_post
- **Outlier behavior**: High gamma — Pre-Post mean VI = 0.964, the most extreme
  rest-to-rest reorganization in the dataset. **[F02]**
- **Explanation**: The high gamma dendrogram transitions from a star-like
  topology (one large hub community) in rest_pre to a modular structure in rest_post.
  This is a genuine biological reorganization, not a metric artifact.
- **Contribution**: Supports H2 and H3 normally for other bands.

### 3.3 Pat_05

- **4 phases**: rest_pre, task_learn, task_test, rest_post
- **No outlier behavior**: Consistent contributor across all hypotheses.
- **Notable**: Ranks high-frequency bands as most reorganized (opposite to
  Pat_02/03), contributing to the rejection of H4.

### 3.4 Pat_06

- **2 phases only**: rest_pre, rest_post (no task phases recorded)
- **Cannot contribute** to H1, H2, or H3 (requires task phases)
- **Pre-Post comparison** included in the 6-patient overview **[G02]**: NVI
  curves fall within the normal range of 4-phase patients, suggesting similar
  baseline variability.

### 3.5 Pat_07

- **3 phases**: rest_pre, task_learn, rest_post (missing task_test)
- **Can contribute** to Pre-Post and Pre-TL comparisons but not to contrasts
  involving task_test (H2a, H2b).
- **Individual profile** shows patterns consistent with the 4-phase group.

### 3.6 Pat_08

- **4 phases**: rest_pre, task_learn, task_test, rest_post
- **Outlier behavior**: Delta band — Pre-Post VI = 0.144 (most similar of all
  pairs), while TL-TT = 0.527. **[F03]**
- **Explanation**: Delta resting-state network structure is ultra-stable across
  the experiment for this patient. The two rest recordings produce nearly
  identical dendrograms at all scales. This makes Pat_08 an extreme supporter
  of H3 for delta (within ≪ cross).

---

## 4. Scale Interpretation

### 4.1 Why Fine Scale Shows More Effect

The predominant finding is that reorganization effects (especially H2) live at
**fine hierarchical scales** (h < 0.3, k > 15). This is because:

1. **Coarse scale convergence**: At h > 0.5 (k < 5), dendrograms collapse into
   2–3 large communities. With so few clusters, most partitions are
   identical — VI approaches 0 trivially. There is no meaningful signal to
   detect.

2. **Information capacity**: Fine partitions carry more information about
   network structure. At k = 50, there are ~50! possible partitions of 117
   nodes into 50 groups, while at k = 2, there are only ~2^116 binary
   partitions. Fine scale is inherently more discriminative.

3. **Neurobiological interpretation**: The task manipulation (learning and
   testing) affects **local circuit-level** community structure rather than
   gross hemisphere-level organization. The brain's macroscopic parcellation
   is robust to cognitive state changes, but the fine-grained sub-community
   structure reflects the specific computations performed during the task.

### 4.2 Scale Ranges

| Scale | h range | k range | Neurobiological level |
|-------|---------|---------|----------------------|
| Macro | > 0.5 | 2–5 | Hemisphere / lobe divisions |
| Meso | 0.15–0.5 | 5–25 | Functional networks (DMN, salience, etc.) |
| Micro | < 0.15 | 25–100+ | Local circuits, sub-communities |

The H2 task trace effects span the micro-to-meso boundary:
- Alpha persist: meso (k ~ 16–26)
- Beta persist: micro (k ~ 69–87)
- Low gamma recover: micro (k ~ 50–98)

The H1 task stability effect (beta) extends from micro well into meso scale.

---

## 5. Methodological Notes

### 5.1 Choice of Distance Metric

VI was selected as the primary metric after systematic comparison with 7
alternatives **[G03]**:

| Metric | Type | H2a unanimous bands | H2b unanimous bands |
|--------|------|--------------------|--------------------|
| mean VI | distance | Best or tied-best | Best or tied-best |
| Cophenetic Pearson | similarity | Fewer | Fewer |
| Baker's Gamma | similarity | Fewer | Fewer |
| Normalized L1 | distance | Comparable | Fewer |
| Normalized L2 | distance | Comparable | Fewer |
| Top-k merge agreement | similarity | Fewer | Fewer |
| Weighted ARI | similarity | Fewer | Fewer |

All metrics produce directionally consistent results, but VI achieves the
highest number of unanimous bands for H2, making it the most sensitive choice.

### 5.2 Log-Spaced Height Grid

The LRG diffusion process generates merge distances that are approximately
log-uniformly distributed. A linear grid (e.g., linspace) would:
- Oversample the coarse scale (h > 0.5) where signal is weak
- Undersample the fine scale (h < 0.1) where most effects live

The geometric spacing `geomspace(0.003, 0.995, 400)` allocates equal resolution
per octave across the hierarchy. Results are robust to grid choice (verified
with 200 and 300 points producing identical conclusions).

### 5.3 NVI Normalization

Without normalization, raw VI at fine scale (k = 50) can reach
$\ln(50) \approx 3.9$ nats, while at coarse scale (k = 3) it is bounded by
$\ln(3) \approx 1.1$ nats. Any aggregate metric (e.g., integral over h) would
be dominated by fine-scale contributions.

NVI = VI / $\ln(\bar{k})$ divides out this scale dependence, yielding values
in [0, ~2] that are directly comparable. NVI $\approx$ 0 means identical
partitions; NVI $\approx$ 1 means approximately random-level disagreement.

### 5.4 Patient Inclusion Criteria

- **H1, H2, H3 testing**: 4 patients with all 4 phases (Pat_02, 03, 05, 08)
- **H4 ranking**: same 4 patients
- **Pre-Post comparison**: all 6 patients (Pat_06 and Pat_07 added)
- **Unanimity threshold**: $|U| = 1$ requires all patients to agree. With
  N = 4, this is a conservative criterion — a single dissenting patient
  prevents the region from being flagged as unanimous.

---

## 6. Figure Manifest

| Code | Filename | Description |
|------|----------|-------------|
| **[A01]** | `A01_all_pairs_nvi_per_band.pdf` | All 6 NVI(h) curves for each band — the foundational raw data |
| **[A02]** | `A02_k_at_h_profile.pdf` | Mapping from cophenetic height h to number of clusters k |
| **[A03]** | `A03_nvi_heatmaps_per_pair.pdf` | Band x h heatmaps for each of the 6 phase pairs |
| **[B01]** | `B01_H1_task_stability_normalized.pdf` | H1: relative gain showing TL-TT is 34% closer on average |
| **[B02]** | `B02_H1_task_stability_multiscale.pdf` | H1: multiscale view with beta as broadest unanimous band |
| **[C01]** | `C01_H2_definitive_heatmap.pdf` | H2: signed unanimity heatmap — alpha persist, low_gamma recover |
| **[C02]** | `C02_H2_per_band_detail.pdf` | H2: per-band panels with individual patient traces |
| **[C03]** | `C03_H2_vi_raw_all_pairs.pdf` | H2: raw VI(h) for 4 key pairs showing curve separations |
| **[C04]** | `C04_H2_task_trace_normalized.pdf` | H2: normalized H2a/H2b at discrete k values |
| **[C05]** | `C05_H2_scale_specific_unanimity.pdf` | H2: macro/meso/micro scale breakdown |
| **[D01]** | `D01_H3_within_vs_cross.pdf` | H3: bar chart with patient dots, 5/6 bands positive |
| **[D02]** | `D02_H3_multiscale_contrasts.pdf` | H3 + all contrasts: 6-panel heatmap at all scales |
| **[E01]** | `E01_H4_frequency_gradient.pdf` | H4: rank heatmap showing patient-specific ordering |
| **[E02]** | `E02_H4_ranking_agreement.pdf` | H4: Kendall's W = 0.043, Spearman correlation matrix |
| **[F01]** | `F01_outlier_Pat02_theta.pdf` | Outlier: Pat_02 theta coarse-scale divergence |
| **[F02]** | `F02_outlier_Pat03_high_gamma.pdf` | Outlier: Pat_03 high_gamma star-to-modular transition |
| **[F03]** | `F03_outlier_Pat08_delta.pdf` | Outlier: Pat_08 delta ultra-stable rest |
| **[G01]** | `G01_all_hypotheses_summary.pdf` | Combined 4-panel hypothesis overview |
| **[G02]** | `G02_pre_post_all_6_patients.pdf` | Pre-Post NVI for all 6 patients |
| **[G03]** | `G03_alternative_metrics_h2.pdf` | Cross-validation of VI against 7 alternative metrics |

---

## 7. Summary of Key Numbers

| Quantity | Value |
|----------|-------|
| Number of patients (hypothesis testing) | 4 |
| Number of patients (Pre-Post) | 6 |
| Number of frequency bands | 6 |
| Number of sEEG nodes | ~117 per patient |
| Height grid | 400 log-spaced points in [0.003, 0.995] |
| FC method | MSC, nperseg = 4096 |
| Number of surrogates | 200 |
| H1 mean relative gain | ~34% |
| H1 broadest band | beta, h = [0.02, 0.34] |
| H2 alpha persistence | h = [0.17, 0.26], k ~ 16–26 |
| H2 beta persistence | h = [0.02, 0.04], k ~ 69–87 |
| H2 low_gamma recovery | h = [0.01, 0.08], k ~ 50–98 |
| H3 bands supporting | 5/6 |
| H4 Kendall's W | 0.043 (weak agreement) |
| Outlier Pat_02 theta TL-TT VI | 0.877 |
| Outlier Pat_03 high_gamma Pre-Post VI | 0.964 |
| Outlier Pat_08 delta Pre-Post VI | 0.144 |
"""


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Curate report figures and generate technical report")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be done without doing it")
    args = parser.parse_args()

    # ── Validate sources ──────────────────────────────────────────────
    missing = []
    for target, source_rel in FIGURES:
        src = SRC / source_rel
        if not src.exists():
            missing.append((target, str(src.relative_to(ROOT))))

    if missing:
        print(f"WARNING: {len(missing)} source file(s) not found:")
        for t, s in missing:
            print(f"  {t} <- {s}")
        print()

    # ── Create output directory ───────────────────────────────────────
    if not args.dry_run:
        DST.mkdir(parents=True, exist_ok=True)

    # ── Copy figures ──────────────────────────────────────────────────
    copied = 0
    for target, source_rel in FIGURES:
        src = SRC / source_rel
        dst = DST / target
        if not src.exists():
            continue
        if args.dry_run:
            print(f"  COPY  {src.relative_to(ROOT)}")
            print(f"     -> {dst.relative_to(ROOT)}")
        else:
            shutil.copy2(src, dst)
            copied += 1
            print(f"  Copied {target}")

    # ── Write report ──────────────────────────────────────────────────
    report_path = DST / "TECHNICAL_REPORT.md"
    if args.dry_run:
        print(f"\n  WRITE {report_path.relative_to(ROOT)}"
              f" ({len(REPORT):,} chars)")
    else:
        report_path.write_text(REPORT, encoding="utf-8")
        print(f"\n  Report written: {report_path.relative_to(ROOT)}"
              f" ({len(REPORT):,} chars)")

    # ── Summary ───────────────────────────────────────────────────────
    prefix = "DRY RUN: " if args.dry_run else ""
    print(f"\n{prefix}Copied {copied}/{len(FIGURES)} figures"
          f" to {DST.relative_to(ROOT)}/")
    if not args.dry_run:
        print(f"Total files in output: {len(list(DST.iterdir()))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
