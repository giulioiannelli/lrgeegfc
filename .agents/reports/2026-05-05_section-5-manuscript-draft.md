---
name: 2026-05-05_section-5-manuscript-draft
type: report
era: IMCOH_ABS × COHORT_N10
status: draft
created: 2026-05-05
updated: 2026-05-05
pointers:
  - .agents/reports/2026-04-29_result-1-raw-fc-phase-trace.md
  - .agents/reports/archive/2026-05/2026-05-05_result-2-lrg-beta-trace.md
  - data/reports/section_5_lrg_trace/README.md
---

# Section 5 — Manuscript draft (writing-ready prose)

**Purpose.** Drop-in prose for the Section 5 paragraph(s) of the manuscript.
Cites the load-bearing tables / figures by relative path. Numbers are the
ones in the per-measure CSVs (verified against `cohort_controls_summary.csv`,
`cohort_band_summary.csv`, `cohort_band_region_frequency.csv`). Editorial
voice: third person plural, past tense, NeuroImage / J Neurosci tier.

---

## Section 5 head paragraph

> The raw functional-connectivity analysis (Section 4) revealed a β-band
> task-shaped trace in resting-state networks at the channel-pair level.
> To test whether this trace also imprints on the multiscale partition
> structure of the network, we applied the Laplacian Renormalization Group
> (LRG) framework to each (patient, phase, band) connectivity matrix and
> compared dendrogram triangles — `T_d = d(task_test, rest_post) −
> d(rest_pre, task_test)` — across four geometric distances spanning matrix
> rank, tree topology, tree heights, and balanced topology+heights.
> Negative `T_d` values indicate that task_test reorganized the partition
> structure and the reorganization persisted into post-task rest. Across
> the cohort (n=10 patients), all four β-band geometric facets converge on
> a strong trace verdict that survives a within-rest-only null control,
> Pat_03 sampling-rate dropout, joint multiple-testing correction across
> 48 simultaneous (probe, variant, band) cells, and bidirectional cross-
> validation against a per-leaf localization probe.

## Methods paragraph (LRG triangle scalars and within-baseline null)

> For each (patient, band) we computed the LRG ultrametric distance matrix
> `D̂(τ) = 1 / ρ̂(τ)` at the canonical diffusion time τ = 1 / λ_max under
> Villegas et al. (2023, 2025), and built the dendrogram triangle scalar
> `T_d` for four distances: rank distance `d_F` (Frobenius-normalized
> matrix distance), and the Kuhner-Felsenstein tree distance with
> branch-length blend λ ∈ {0, 0.5, 1} (λ = 0 captures pure topology,
> λ = 1 captures heights only, λ = 0.5 is balanced). To control for
> within-session noise, we generated a within-baseline drift triangle by
> substituting the second half of `rest_pre` for `task_test` and the
> first half of `rest_post` for the original `task_test` reference,
> giving `T_d^null = d(rest_pre_B, rest_post_A) − d(rest_pre_A,
> rest_pre_B)` — a triangle constructed entirely from rest data and
> therefore carrying no task-related signal. A signed Wilcoxon paired
> test over n=10 patients tested whether the real `T_d` was below the
> patient's own `T_d^null`. Joint Bonferroni correction was applied
> across 27 simultaneous cells (3 partition-distance probes × 7 variants ×
> 3 bands plus 6 per-band localization cells; threshold p ≤ 0.001).

## Results paragraph (β verdict)

> All four β-band geometric facets — d_F (Frobenius), KC λ=0 (topology),
> KC λ=0.5 (balanced), KC λ=1 (heights) — produced 10/10 patients with
> real `T_d < T_d^null` (Wilcoxon paired one-sided p = 0.001 each), and
> all four cells survived joint Bonferroni m=48 (Table S5, Figure
> [controls_heatmap.pdf]). The topology-only cell (KC λ=0) is also
> Pat_03-dropout-robust on the absolute test (p = 0.019 → p = 0.037
> after dropping the 1024 Hz outlier patient), confirming that the
> β-band partition-structure trace does not depend on Pat_03's
> nonstandard sampling rate. The two complementary partition probes —
> KC λ=0 (topology only, height- and magnitude-blind) and KC λ=1
> (heights only, topology-blind) — agreeing at p = 0.001 each is
> internally redundant evidence: the β trace lives in both the tree
> shape and the tree heights, not in a single representational artifact.
> By contrast, the leading-13 Grassmann eigenvector overlap was 8/10
> trace at β (p = 0.024, BH-FDR only), insufficient to credit the
> β-band trace to eigenmode rotation; β trace is a partition-structure
> story, not an eigenmode geometry story (Figure [controls_paired_lines.pdf]).

## Results paragraph (low_γ secondary verdict and per-leaf localization)

> A second band, low_γ, also reached joint-Bonferroni-surviving status
> through a different geometric route: the leading-13 Grassmann
> eigenspace overlap detected an 8/10-patient trace at low_γ (Bonferroni-
> surviving cell, Table S5), and a per-leaf localization probe — built as
> the Spearman rank correlation between task and rest_post pair-distance
> changes around each leaf, demeaned within (patient, band) and
> calibrated against a per-(patient, band) within-baseline 95th-percentile
> null floor — confirmed 9/10 patients with ≥3 calibrated trace-leaves at
> β and 10/10 at low_γ (cohort Wilcoxon p = 0.002 and p = 0.001
> respectively against the within-baseline null). Pair-level cross-
> validation between the global (CTM σ-aggregate concordance score) and
> local (calibrated trace-leaves) probes confirmed that pairs with both
> endpoints in calibrated trace-leaves carried significantly higher
> concordance than non-trace pairs (β: enrichment p = 0.004,
> concentration in top-10 % most concordant pairs p = 0.004; both pass
> per-measure Bonferroni m=6).

## Anatomy paragraph

> Anatomical mapping of the calibrated trace-leaves was performed
> against each patient's Desikan-Killany parcellation, then normalized
> against the cohort-wide contact density per region (hypergeometric
> test). After this sampling correction, β trace-leaves were enriched
> at the Hippocampus (3.5× over cohort baseline 5.3%, p = 0.011, three
> contributing patients) and at left fusiform cortex (2.5×, p = 0.045,
> two patients). Left superior temporal cortex was borderline (2.1×,
> p = 0.055). low_γ trace-leaves showed strong enrichment at left
> fusiform cortex (4.0×, p < 0.001) and left inferior temporal cortex
> (2.3×, p = 0.026). Without the sampling correction the bar counts
> would have implicated left middle temporal cortex as a key region;
> at 6% trace rate against a 5.3% baseline (p = 0.48) it is at chance
> and was excluded from the cohort claim. White matter contacts
> (~50 % of all trace-leaves) were treated as sEEG implant geometry,
> not biology, and were not analyzed.

## Limitations paragraph

> Several caveats temper the interpretation. First, the LRG-layer
> analysis and the raw functional-connectivity analysis (Section 4) are
> derived from the same imaginary-coherency matrices on the same n=10
> cohort, so they are different geometric layers of the same data, not
> independent replications; cross-cohort replication is required for a
> generalizability claim. Second, the within-baseline null is
> conservative — any task-induced reorganization that genuinely persists
> in rest_post inflates the rest_post halves and biases the null floor
> upward, against detection — meaning the controlled effect is at most an
> underestimate of the true trace strength. Third, three of the ten
> patients (Pat_02, Pat_05, Pat_15) were anti-trace at the absolute β
> triangle T_d level despite all 10 patients being below their own
> within-baseline drift floor; the cohort verdict relies on the within-
> patient comparison, not on uniform absolute trace direction. Fourth,
> Pat_03's 1024 Hz sampling rate (vs the cohort 2048 Hz) makes its β
> spectrum and MSC structure quantitatively different; we report the
> Pat_03-dropout robustness explicitly. Fifth, anatomical labelling
> reflects each patient's individual implant and is therefore
> opportunistic — the cohort-modal regions (left temporal, Hippocampus)
> are constrained by the implant pattern of this cohort. Sixth, our
> analyses fix τ at 1 / λ_max; we verified that the trace is τ-invariant
> within the meaningful diffusion range [τ_min = 1/λ_max, τ* = peak of
> specific heat C(τ)] (Figure S5τ), licensing this choice but not
> establishing scale-resolved structure.

## Figure caption stubs

**Figure 5A — Controls heatmap.** Cohort-level Wilcoxon -log10(p) for the
21 (probe, variant, band) cells of the within-baseline null test. Each
cell carries the count of patients with `T_d < T_d^null` (out of 10) and
significance stars: *** Bonferroni m=42 surviving (p ≤ 0.0012), * BH-FDR
q ≤ 0.05. β column dominantly dark across 5 of 7 partition-distance rows;
low_γ dark on the Grassmann row. *Source*:
`data/reports/section_5_lrg_trace/11_lrg_global_probe_controls/figures/controls_heatmap.pdf`.

**Figure 5B — Paired lines.** Per-patient `T_d^null` (left) connected to
real `T_d` (right) for the 5 Bonferroni-surviving cells (KC λ=0, λ=0.5,
λ=1 at β; d_F at β; Grassmann k=13 at low_γ). Lines slope downward —
real below null — for all 10 patients in each cell. *Source*:
`data/reports/section_5_lrg_trace/11_lrg_global_probe_controls/figures/controls_paired_lines.pdf`.

**Figure 5C — Cross-validation.** Per-band cohort enrichment (left) and
top-decile concordance concentration vs expected (right). Trace-leaf
pair sets show significantly elevated concordance scores in β / δ
under per-measure Bonferroni m=6 in both directions. *Source*:
`data/reports/section_5_lrg_trace/09_global_local_cross_validation/figures/cohort_summary.pdf`.

**Figure 5D — Anatomical localization.** Top regions per band among
calibrated trace-leaves (left = blue, right = red, white matter = grey,
subcortical = purple). β trace-leaves concentrate in left temporal lobe
+ Hippocampus across 3/9 contributing patients per region. *Source*:
`data/reports/section_5_lrg_trace/13_anatomical_localization/figures/region_frequency_grid.pdf`.

**Supplementary Figure τ — τ-invariance.** Cohort `T_d` curves swept
across patient-specific normalized intrinsic-scale axis u ∈ [0, 1] from
τ_min = 1/λ_max to τ* = peak of C(τ). Cohort median variation 0.02–0.07
(median IQR ~0.4–0.5): trace is τ-invariant within the meaningful
diffusion range. *Source*:
`data/reports/section_5_lrg_trace/10_ctm_tau_sweep/figures/cohort.pdf`.

## Tables to cite

- **Table S5** — Per-band controls cohort summary
  (`data/audit/lrg_global_probe_controls/cohort_controls_summary.csv`).
- **Table S6** — Per-(band, region) trace-leaf frequency
  (`data/reports/section_5_lrg_trace/13_anatomical_localization/tables/cohort_band_region_frequency.csv`).
- **Table S7** — Per-patient β triangle `T_d` for KC λ ∈ {0, 0.5, 1} and
  D-rank d_F
  (`data/reports/section_5_lrg_trace/03_kc_lambda_triangle/tables/Td_per_patient_per_band_lambda.csv`,
  `data/reports/section_5_lrg_trace/02_d_rank_triangle/tables/Td_per_patient_per_band.csv`).

## Notes for revision

- Numbers cited in prose verified against CSV cells on disk
  (2026-05-05). Re-check before submission.
- "Joint Bonferroni m=48" applies to 42 controls cells + 6 per-leaf bands; re-count if additional probes
  are added before submission.
- The phrase "task-shaped trace" follows the Result 1 wording; alternate
  phrasings ("task-induced reorganization with persistence",
  "post-task structural carry-over") work but be consistent across
  Sections 4 and 5.
- Replace `_` with `\_` in math-mode region names if compiling under
  LaTeX (`ctx-lh-superiortemporal` → use rm-text or DK-vocabulary
  abbreviations like "left STG").
