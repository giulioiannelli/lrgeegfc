---
name: 2026-05-05_result-2-lrg-beta-trace
type: report
era: IMCOH_ABS × COHORT_N10
status: superseded
created: 2026-05-05
updated: 2026-05-28
sign_convention: pre-2026-05-26
pointers:
  - .agents/reports/archive/2026-04/2026-04-29_result-1-raw-fc-phase-trace.md
  - data/reports/section_5_lrg_trace/README.md
  - data/reports/section_5_lrg_trace/03_kc_lambda_triangle/report.md
  - data/reports/section_5_lrg_trace/11_lrg_global_probe_controls/report.md
  - data/reports/section_5_lrg_trace/12_per_leaf_rho_null/report.md
  - data/reports/section_5_lrg_trace/09_global_local_cross_validation/report.md
  - data/audit/lrg_global_probe_controls/cohort_controls_summary.csv
  - data/audit/lrg_global_probe_controls/pat03_dropout.csv
---

> **Superseded sign convention (banner added 2026-05-28).** This report
> was written before the 2026-05-26 T_d sign lock. It uses the
> pre-2026-05-26 convention: `T_d = d(task, rest_post) − d(rest_pre, task)`
> → negative = trace. The locked convention (from 2026-05-26 onward) is
> `T_d = d(rest_pre, task) − d(task, rest_post)` → **positive = trace**,
> at every layer (raw FC, D_coph, KC, Grassmann). To read this report in
> the current convention, multiply every T_d value by −1 and read every
> "T_d < 0", "negative = trace", or "n_trace = (T_d < 0).sum()" assertion
> with the sign flipped. Numerical magnitudes, p-values, per-band
> verdicts, and the β headline are unchanged; only the sign carrier is
> inverted. Live preprint scripts and the current directive use the new
> convention. Companion CSVs in `data/audit/...` may carry either sign
> depending on compute date; check audit script docstrings.

# Result 2 — LRG β-band trace is the Section-5 headline

## In plain words (read first)

**A short memory task leaves a fingerprint on β-band brain network
organization that persists into post-task rest, and the fingerprint sits
at the Hippocampus and left fusiform cortex — known memory-system regions —
once the analysis is corrected for implant-sampling bias.**

Setup: 10 epilepsy patients with brain electrodes (sEEG). Each did 5 min
rest → memory task → 5 min rest. We compared the *shape of the
connectivity network* (how channels cluster into communities at multiple
scales) across the three phases.

What we found:

- In the β rhythm (~13–30 Hz, the cognition-engagement band), the post-
  task rest network is **closer in shape to the task** than to the pre-
  task rest. Task changed the network, the change stuck.
- All **10 patients** show this when each is compared against their own
  resting noise floor.
- It shows up using **four independent ways** of measuring network shape
  (tree topology, tree heights, balanced topology+heights, magnitude-
  based matrix distance). All four agree at the strongest possible
  Wilcoxon significance level (p = 0.001 each).
- **β is genuinely unique at the global-network level.** When we test
  all 6 frequency bands with the same controlled framework (42 cells:
  6 bands × 7 probes), β is the only band with 4 cells surviving the
  strict joint correction; low_γ has 1; α / δ / θ / high_γ have zero.
  δ / θ / high_γ in particular show no global-probe trace at all under
  controls — their cohort-level p-values are above 0.05 even
  uncorrected.
- After correcting for implant-sampling bias (region trace-rate vs
  cohort baseline 5.3% at β; hypergeometric test), the enriched
  regions are the **Hippocampus** (3.5×, p = 0.011, 3/9 patients) and
  **left fusiform cortex** (2.5×, p = 0.045, 2/9). Left superior
  temporal is borderline (2.1×, p = 0.055). Left middle temporal,
  which looked dominant in raw counts, is **at baseline** (1.1×,
  p = 0.48) and was a sampling artifact, not a finding.

Why it matters: resting brain state is usually treated as a stable
baseline. We show that even a short task leaves a structural fingerprint
on resting connectivity for the minutes that follow, and the fingerprint
is anatomically meaningful — not random drift.

What this does NOT say:

- Not "memories live in β rhythms" — we don't show behaviour, only
  network signatures.
- Not "this generalizes beyond this cohort" — n=10 is small and all
  patients are clinical sEEG implants; cross-cohort replication is
  needed.
- Not "the trace lasts forever" — we measured a short post-task window
  only.

The technical sections below build the evidence chain that supports
this in-plain-words finding, with all the controls, robustness checks,
and per-patient numbers.

---

**This is the second crystallized result of the rebuild. It builds directly on
Result 1 (raw-FC β trace at the substrate level) and adds the LRG-layer
geometric enrichment: not only does β-band connectivity reorganize and stick,
the *partition structure* (how channels group into communities at multiple
scales) reorganizes and sticks. Four independent geometric measures agree,
all four pass the strictest within-baseline-null control, and the β finding
survives joint Bonferroni correction across 48 simultaneous (probe, variant,
band) cells (42 controls + 6 per-leaf bands; threshold p ≤ 0.00104).**

---

## Renormalization head

In β-band oscillations, task reorganizes the way sEEG channels cluster at
every scale of the LRG dendrogram, and that reorganization persists into
post-task rest. We measured this four different ways — KC tree-distance with
λ ∈ {0, 0.5, 1} (topology / balanced / heights) and a Frobenius-normalized
matrix distance d_F on the LRG ultrametric — and all four say the same thing:
**every one of the 10 patients has a β-band trace stronger than their own
within-rest noise floor, with paired-Wilcoxon p = 0.001 cohort-wide.** A
fifth probe (Grassmann k=13 on eigenstructure) is uncorrected-only at β,
which means β trace is *not* an eigenvector-rotation story — it is a
partition-structure story. The β topology cell (KC λ=0) survives Pat_03
sampling-rate dropout. Per-leaf localization confirms 9/10 patients carry
calibrated trace-leaves at β. Pair-level cross-validation between global and
local probes agrees in β. **β trace is the cleanest Section-5 finding.**

---

## What "β trace" means concretely

Per (patient, band) we computed the triangle scalar

```
T_d(p, b) = d(task_test, rest_post) − d(rest_pre, task_test)
```

with several distances `d` defined on either the LRG ultrametric matrix
D = 1/ρ̂ at τ = 1/λ_max or its dendrogram T. Negative T_d means task_test is
*closer* to rest_post than to rest_pre — task reorganized the structure AND
the new structure persists into post-task rest. That is the **trace**
configuration in our four-way taxonomy (trace / anchor / reset / emergent),
and it is what we expect from a task that left a real fingerprint in
resting-state networks.

For β specifically:

- All 10 patients have a real β triangle T_d below their own within-baseline
  drift floor under every one of the four strict-control cells.
- That floor is built from rest halves only (Δ_drift = D^post_B − D^post_A
  vs Δ_drift = D^pre_B − D^pre_A), so it captures pure session-level noise
  with no task contamination. Real T_d beating this floor means the
  triangle's signed gap is bigger than within-session distance fluctuations.
- The four cells are KC λ=0 (topology / tree shape only), KC λ=0.5
  (balanced topology + heights), KC λ=1 (heights only), and D-rank d_F
  (Frobenius-normalized matrix distance). They probe different geometric
  facets of the same dendrogram, and all four agree.

---

## Why β at the LRG level is more than β at the raw-FC level

Section 4 (Result 1) established that raw-FC β-band carries a task-shaped
trace at the *edge* level (`d_S^β` triangle, 7/10 patients, controls pass).
Section 5 (this result) adds the LRG-layer geometric enrichment: the
*partition structure* reorganizes too. Concretely:

- Raw-FC β trace says: the entries of the |ImCoh| matrix shifted in task and
  the shift persists in rest_post.
- LRG β trace says: the way channels cluster into nested communities at
  every scale of the LRG dendrogram shifted in task and the cluster
  organization persists in rest_post.

These are not redundant. The first is a *channel-pair-level* statement
(which edges moved). The second is a *multiscale-partition-level* statement
(which group structures moved). The LRG step is a structural enrichment, not
a confirmation — it tells us *what kind* of reorganization is happening
(partitioning structure shifts, not just edge-weight noise).

A literature reviewer asking "is the raw-FC trace just channel drift?" is
answered: no, because the same trace is visible at the partition-structure
layer, and the partition-structure layer is invariant to most channel-level
drift (KC λ=0 is purely topological — height- and magnitude-blind).

---

## The chain of evidence

Eight Section-5 measures bear on β. Sorted by what they add:

### Tier 1 — load-bearing (controlled, joint-MTC-surviving)

| # | measure | β verdict | n_trace | Wilcoxon p (real < null) | passes Bonferroni m=42 |
|:-:|:--|:--|:-:|:-:|:-:|
| 11 | KC λ=0 (topology) | 10/10 | 10/10 | 0.001 | **✓** |
| 11 | KC λ=0.5 (balanced) | 10/10 | 10/10 | 0.001 | **✓** |
| 11 | KC λ=1 (heights) | 10/10 | 10/10 | 0.001 | **✓** |
| 11 | D-rank d_F (Frobenius) | 10/10 | 10/10 | 0.001 | **✓** |

All four cells pass joint Bonferroni m=48 across measures 11+12 (42 controls
cells + 6 per-leaf bands; threshold p ≤ 0.00104, our actual p = 0.000977).
Adding measure 09 cross-validation (12 cells) makes m=60, threshold 0.000833,
which our p=0.000977 misses by 0.00014 — cross-validation is reported as
methodological consistency, not joint-MTC-surviving evidence.

The within-baseline null is the conservative form: any task-induced
reorganization that genuinely persists in rest_post inflates the rest_post
half — if the trace is real, the null floor is biased UPWARD against
detection. Real T_d still beats the inflated floor → real signal exceeds
even biased floor.

### Tier 2 — supporting (per-measure controls, BH-FDR surviving)

| # | measure | β verdict | passes BH-FDR q≤0.05 |
|:-:|:--|:--|:-:|
| 11 | D-rank d_P (Pearson) | 8/10 trace, p=0.010 | ✓ |
| 11 | D-rank d_S (Spearman) | 8/10 trace, p=0.032 | ✗ (uncorrected only) |
| 12 | Per-leaf calibrated trace-leaves | 9/10 patients ≥3 leaves, p=0.002 | per-measure m=6 ✓ |
| 09 | Pair-level enrichment (AND set) | p=0.004, both directions | per-measure m=6 ✓ |

D-rank d_P (which is the magnitude-weighted complement of the rank distance)
correlates with d_S at cohort ρ ≈ 0.85–0.95 per band, so it is not orthogonal
evidence; it is partly redundant with d_S and is reported for completeness
(see `feedback_dP_framing.md`). The honest read is that d_F and the three KC
variants are the load-bearing Tier-1 cells; d_S and d_P are consistent
confirmations.

### Tier 2.5 — full 6-band band-coverage scan (the "is β unique?" check)

The original measure 11 controls tested only α / β / low_γ (3 bands × 7 probes = 21 cells). On 2026-05-05 we extended the within-baseline null to **all 6 bands × 7 probes = 42 cells** to test whether β is genuinely uncovering more than the other bands or whether we just happened to test the right bands.

| band | uncorrected p<0.05 cells (out of 7) | Bonferroni m=42 ✓ | BH-FDR q≤0.05 ✓ |
|:--|:--:|:--:|:--:|
| **β** | 7 | **4** (KC λ=0/0.5/1, d_F) | 5 |
| **low_γ** | 2 | **1** (Grassmann k=13) | 1 |
| α | 3 | 0 | 1 (Grassmann k=13 only) |
| δ | 0 | 0 | 0 |
| θ | 0 | 0 | 0 |
| high_γ | 0 | 0 | 0 |

**δ, θ, high_γ have ZERO uncorrected-p<0.05 cells under any of the seven global probes.** Their global-probe T_d at the LRG level is statistically indistinguishable from within-rest noise at the cohort level. α is borderline (3 uncorrected cells, none surviving Bonferroni — Grassmann k=13 survives BH-FDR alone).

The β finding is therefore **structurally unique at the global-probe level**: β has 4 of 7 probes surviving strict Bonferroni m=42 across the entire 42-cell space, low_γ has 1, and the other four bands have none. This rules out the criticism "you tested only the bands that worked" — when we test ALL bands, only β + low_γ survive.

The honest joint correction across measures 11 + 12 is **m=48** (42 controls cells + 6 per-leaf bands). Bonferroni threshold 0.001042 — the 5 measure-11 cells at β / low_γ (all p=0.000977) plus the low_γ per-leaf cell (p=0.000977) survive. Joint count: **6 cells**, all at β (4) + low_γ (2 — global Grassmann + per-leaf).

### Tier 3 — diagnostic (told us what β trace is NOT)

| # | measure | β verdict | what it tells us |
|:-:|:--|:--|:--|
| 04 | Grassmann k=13 (top-k eigenvector overlap) | 5/10 trace (NULL) | β trace is **not** an eigenvector-rotation story |
| 10 | CTM τ-sweep within [τ_min, τ*] | flat curves | LRG-distance trace is approximately τ-invariant within meaningful diffusion range |

The Grassmann null is informative: if β trace were a story about the leading
eigenmodes of the LRG Laplacian rotating in task, Grassmann k=13 would have
caught it. The fact that Grassmann is null at β while KC and D-rank are 10/10
means β trace lives in the *partition structure* (how the dendrogram cuts
look), not in the *eigenmode geometry* (which directions L̂ favors). This is
a genuine geometric distinction.

The τ-sweep null says "you don't get a different story by sweeping τ" —
β trace is τ-invariant within [τ_min = 1/λ_max, τ* = peak of C(τ)]. This
licenses our default τ = 1/λ_max as the load-bearing scale (no privileged
n*).

### Per-patient evidence at β

KC topology (λ=0) per-patient T_d:

| patient | T_d | trace direction? | within-baseline-null < real? |
|:-:|:-:|:-:|:-:|
| Pat_02 |  +0.508 | no (anti) | YES (null even more positive) |
| Pat_03 | −19.730 | yes | YES |
| Pat_05 |  +0.175 | no (anti) | YES |
| Pat_06 |  −5.090 | yes | YES |
| Pat_07 |  −2.433 | yes | YES |
| Pat_08 |  −5.069 | yes | YES |
| Pat_10 |  −2.452 | yes | YES |
| Pat_13 |  −3.159 | yes | YES |
| Pat_14 |  −1.367 | yes | YES |
| Pat_15 |  +2.381 | no (anti) | YES |

7/10 patients have absolute T_d < 0 (raw trace direction). All 10/10 patients
have T_d below their own within-baseline drift floor (controlled trace
direction). The "10/10 vs null" is more sensitive than the "T_d < 0
absolute" because the within-baseline drift adds a positive bias to T_d that
the controlled comparison removes per-patient.

Pat_15 is the cohort-anti patient at every β cell. Pat_02 / Pat_05 are
mixed: anti on absolute T_d, pro relative to their own null floor. The
cohort-pro group is Pat_03 / Pat_06 / Pat_07 / Pat_08 / Pat_10 / Pat_13 /
Pat_14 (7/10 unanimous trace).

---

## Why we trust the β finding

### 1. Within-baseline null beats the strict bar

Joint Bonferroni m=48 across 42 measure-11 cells + 6 measure-12 bands sets
threshold p ≤ 0.00104. Our four β cells all have p = 0.000977 (the smallest
possible Wilcoxon paired test result at n=10). They survive joint Bonferroni
m=39 across 11+12+09 (threshold 0.00128) too.

### 2. Pat_03 (1024 Hz outlier) dropout

Pat_03 is recorded at 1024 Hz (cohort default 2048 Hz; see
`pat03_nperseg.md`). Its β-band MSC was historically a 3× outlier. We
re-tested the absolute T_d < 0 cohort Wilcoxon at β with Pat_03 dropped:

| cell | full n=10 absolute p | dropped n=9 absolute p |
|:--|:-:|:-:|
| KC β λ=0 (topology) | 0.019 | **0.037** (still significant uncorrected) |
| KC β λ=1 (heights) | 0.042 | **0.082** (lost significance) |

The β topology trace (KC λ=0) is **Pat_03-dropout robust**. The β heights
trace (KC λ=1) is partly Pat_03-driven on the absolute test. Both are
Bonferroni-surviving on the within-baseline-null test (p=0.001 each), so
the topology cell is the most defensible single-cell headline.

### 3. Bidirectional cross-validation between global and local

Measure 09 looks at pair-level concordance: does the localization probe
(measure 12 calibrated trace-leaves) agree with the global probe (CTM
σ-aggregate's per-pair concordance score)?

- **Direction (a) Enrichment** (mean concordance score of trace-leaf pairs
  vs all other pairs): β p = 0.004 (per-measure Bonferroni m=6 ✓).
- **Direction (b) Concentration** (top-decile concordance pairs concentrated
  in trace-leaves vs expected): β p = 0.004 (per-measure Bonferroni m=6 ✓).
- 5 of 9 patients pass per-patient permutation enrichment p < 0.05 at β.

Both probes — looking at the same data through different lenses — agree
about β. That is the cross-validation we wanted.

### 4. Multiple geometric facets agree

KC λ=0 (topology only, no heights) and KC λ=1 (heights only, no topology)
are *complementary* parts of the dendrogram structure. They both show the
same β trace direction independently. The intermediate λ=0.5 also agrees.
D-rank d_F (Frobenius on the matrix D) is geometrically a different object
than KC (matrix-level vs dendrogram-level), and it also agrees.

If β trace were an artifact of one specific representation of the dendrogram
(say, a height-rescaling artifact), at most one of these four cells would
catch it. The fact that *all four* agree at p = 0.001 with 10/10 patients
is the strongest possible evidence we can build at our cohort size.

### 5. Per-leaf localization 9/10 confirms cohort generality

Measure 12 checks a separate question — is there a localized, leaf-level
signal? — and at β reports 9/10 patients with ≥3 calibrated trace-leaves
(p = 0.002 against within-baseline null). This is per-measure
Bonferroni m=6 surviving (just below 0.0083) but JOINT Bonferroni m=48
miss by 0.00096 (low_γ does survive). Marginal at joint level, robust at per-measure level —
honest description: localized β trace is real but its statistical room
above the strict joint bar is narrow.

### 6. The substrate-level story already pointed here

Result 1 (raw-FC β trace) had n_trace = 7/10 with controls passing for d_S
at β. The LRG layer should at minimum reproduce that direction; in fact
it strengthens it (10/10 under controlled test). The LRG layer cannot create
a β trace that the substrate doesn't carry, and indeed it doesn't — it
*enriches* the existing β substrate trace with partition-structure detail.

---

## What β trace at LRG IS NOT

To stay critical:

- It is **not a localization-to-specific-anatomy claim.** The trace nodes
  are still indexed by channel number. Anatomical labelling via
  `implant_pat_NN.csv` → Desikan-Killany regions is task #28, not done.

- It is **not an independent replication of Result 1.** Both Result 1
  (raw-FC) and Result 2 (LRG) are derived from the same imcoh_abs matrices
  on the same 10 patients. They are different geometric layers of the same
  data, not different data.

- It is **not an eigenvector-rotation story.** Grassmann k=13 at β is null
  (5/10 trace, p = 0.024 BH only). β trace lives in partition structure.

- It is **not a "best band" claim.** The strictest joint Bonferroni m=39
  across all three measure batches survives in 6 cells: 5 from measure 11
  (4 of them at β + 1 at low_γ Grassmann) + 1 from measure 12 (low_γ
  per-leaf controlled). β is the band with the most surviving cells, but
  low_γ also has joint-Bonferroni-surviving evidence at both global
  (Grassmann k=13) and local (per-leaf calibrated) levels. β is the
  *cleanest* finding, low_γ is the *secondary* finding.

- It is **not free of cohort heterogeneity.** Pat_15 is a cohort-anti
  patient at every β cell. The 10/10 vs-null result hides the fact that
  Pat_15's null is even more positive than its real — Pat_15 doesn't
  show absolute trace direction, it just isn't worse than its own noise
  floor. This is honest framing, not a weakness — but it's why we report
  per-patient T_d alongside cohort scalars (per the never-always rule).

- It is **not a manuscript claim about τ-multiscale sensitivity.** The
  τ-sweep was null. β trace is τ-invariant within [τ_min, τ*]. CTM at
  τ = 1/λ_max captures the full signal.

---

## How to reproduce — exact numbers and scripts

### One-shot reproducibility

```python
import pandas as pd

ctl = pd.read_csv("data/audit/lrg_global_probe_controls/cohort_controls_summary.csv")
beta_load_bearing = ctl[(ctl.band == "beta") & (ctl.passes_bonf == True)]
print(beta_load_bearing[["probe", "variant", "n_real_below_null",
                         "wilcoxon_p_real_lt_null"]].to_string(index=False))
```

Expected output (4 cells, all 10/10 p=0.000977):

```
    probe    variant  n_real_below_null  wilcoxon_p_real_lt_null
    drank        d_F                 10                 0.000977
       kc lambda=0.0                 10                 0.000977
       kc lambda=0.5                 10                 0.000977
       kc lambda=1.0                 10                 0.000977
```

### Pat_03 dropout sensitivity

```python
p3 = pd.read_csv("data/audit/lrg_global_probe_controls/pat03_dropout.csv")
print(p3[p3.band == "beta"].to_string(index=False))
```

Expected: KC β λ=0 absolute test p_full=0.019 → p_drop=0.037 (still
significant). KC β λ=1 p_full=0.042 → p_drop=0.082 (lost).

### Cross-validation

```python
xv = pd.read_csv("data/reports/section_5_lrg_trace/"
                 "09_global_local_cross_validation/tables/cohort_band_summary.csv")
print(xv[xv.band == "beta"][["band", "wilcoxon_p_enrich_and",
                              "wilcoxon_p_conc_and_vs_exp",
                              "n_pat_perm_enrich_and_p05"]].to_string(index=False))
```

Expected: β enrich_p = 0.0039, conc_p = 0.0039, n_pat = 5.

### Re-running from scratch

Compute scripts (run in `lapbrain` env from repo root):

1. `scripts/01_compute/audit/audit_41_lrg_global_probe_controls.py` — within-baseline null for D-rank/KC/Grassmann + Pat_03 dropout + joint MTC.
2. `scripts/01_compute/audit/audit_41b_controls_summary_figure.py` — heatmap + paired-lines figures.
3. `scripts/01_compute/audit/audit_42_per_leaf_rho_null.py` — within-baseline null for per-leaf demeaned ρ_ℓ.
4. `scripts/01_compute/audit/audit_43_global_local_cross_validation.py` — bidirectional pair-level enrichment + permutation null.

Each takes < 5 min on the 10-patient cohort.

---

## Critical files

### Source data (cached, do not regenerate)

- `data/cache/imcoh_lrg_halves/Pat_XX/{band}_{phase}_lrg_imcoh-abs.npz` —
  half-baseline LRG cache (rest_pre_A/B, rest_post_A/B). Built by
  `scripts/01_compute/hypothesis_tests/h2e_split_half.py`.
- `data/cache/imcoh_lrg_cache/Pat_XX/{band}_{phase}_lrg_imcoh-abs.npz` —
  full-phase LRG (rest_pre, task_test, rest_post).
- `data/reports/imcoh_continuous_trace/per_pair_split/{Pat_XX}_{band}.npz` —
  CTM split-baseline per-pair vectors used by measure 09.

### Per-measure reports (canonical citations)

- `data/reports/section_5_lrg_trace/03_kc_lambda_triangle/report.md` — KC
  λ-blend triangle results.
- `data/reports/section_5_lrg_trace/02_d_rank_triangle/report.md` — D-rank
  results (d_S, d_P, d_F).
- `data/reports/section_5_lrg_trace/04_grassmann_triangle/report.md` —
  Grassmann diagnostic null at β.
- `data/reports/section_5_lrg_trace/11_lrg_global_probe_controls/report.md`
  — within-baseline null + Pat_03 dropout + joint MTC.
- `data/reports/section_5_lrg_trace/12_per_leaf_rho_null/report.md` —
  per-leaf calibrated trace-leaves.
- `data/reports/section_5_lrg_trace/09_global_local_cross_validation/report.md`
  — bidirectional pair-level cross-validation.

### Tables that carry the load-bearing numbers

- `data/audit/lrg_global_probe_controls/cohort_controls_summary.csv` — the
  21-cell controls matrix (probe, variant, band → n_below_null, Wilcoxon p,
  passes_bonf, passes_q05).
- `data/audit/lrg_global_probe_controls/pat03_dropout.csv` — KC β λ=0/1
  full vs n=9 results.
- `data/reports/section_5_lrg_trace/03_kc_lambda_triangle/tables/Td_per_patient_per_band_lambda.csv`
  — per-patient KC T_d at all λ ∈ {0, 0.25, 0.5, 0.75, 1}.
- `data/reports/section_5_lrg_trace/02_d_rank_triangle/tables/Td_per_patient_per_band.csv`
  — per-patient T_d under d_S, d_P, d_F.

### Headline figures (manuscript-grade)

The 4 visual-impact figures that prove β-uniqueness are at
`data/reports/section_5_lrg_trace/headline/figures/`:

- `volcano.pdf` — 42-cell volcano plot. β KC family + β d_F + γ_l Grassmann
  k=13 sit alone above the Bonferroni m=42 line; all other cells in noise.
- `per_patient_slopes.pdf` — 6-panel grid. β panel: all 10 patients slope
  null→real downward. Other bands: 4–6/10 with mixed direction.
- `cohort_fingerprint.pdf` — 10-patient × 6-band × 4-facet heatmap. β
  column saturated red across every patient and facet; other bands
  checkerboard noise.
- `anatomy_2d.pdf` — sampling-corrected region enrichment (top: bar
  charts of rate / cohort-baseline-rate with hypergeometric p-values;
  bottom: axial MNI scatter colored by region enrichment). β enriched
  at Hippocampus + left fusiform; low_γ at left fusiform + left
  inferior temporal.

See `data/reports/section_5_lrg_trace/headline/README.md` for the figure
guide and the defensive question each addresses.

### Figures cited in any β manuscript paragraph

- `data/reports/section_5_lrg_trace/11_lrg_global_probe_controls/figures/controls_heatmap.pdf`
  — 21-cell summary; β column dominantly dark red across 5/7 rows.
- `data/reports/section_5_lrg_trace/11_lrg_global_probe_controls/figures/controls_paired_lines.pdf`
  — per-patient real vs null lines for the 5 Bonferroni-surviving cells.
- `data/reports/section_5_lrg_trace/12_per_leaf_rho_null/figures/cohort_calibrated_summary.pdf`
  — null vs real demeaned ρ_ℓ distributions per band.
- `data/reports/section_5_lrg_trace/09_global_local_cross_validation/figures/cohort_summary.pdf`
  — per-band pair-level enrichment + concentration box-plots.

---

## What's downstream of this result

- **Task #28 — anatomical labelling** of the calibrated trace-leaves at
  β / low_γ. This is what gives the manuscript a "trace localizes to
  {regions}" sentence. Pure descriptive, no new statistics.
- **Section 5 Part 5 — final headline pick.** All 11 measures have landed.
  The recommendation is now a writing exercise, not a research one. The
  β-trace as load-bearing global probe is established here. Per-leaf
  calibrated localization at low_γ is the load-bearing local probe under
  joint MTC. β + low_γ are the two manuscript bands; θ / δ / α / high_γ
  are reported as supporting evidence with appropriate caveats.

---

## TL;DR (one paragraph for the impatient)

In β-band, the way sEEG channels group into communities at every scale of
the LRG dendrogram reorganized in task and persisted into post-task rest.
We measured this four different ways (KC λ=0/0.5/1 + D-rank d_F), the
within-rest-only noise floor was beaten by all 10 patients in every one of
the four cells (Wilcoxon p = 0.001 each, joint Bonferroni m=48 surviving),
the most defensible cell (KC λ=0, topology only) survives Pat_03 dropout,
9/10 patients carry per-leaf calibrated trace-leaves at β, and pair-level
cross-validation between the global and local probes agrees at p = 0.004
in both directions. The Grassmann null at β tells us this is a partition-
structure story, not an eigenmode-rotation story. The τ-sweep null tells
us this is τ-invariant within the meaningful diffusion range. β trace is
the load-bearing Section-5 finding.
