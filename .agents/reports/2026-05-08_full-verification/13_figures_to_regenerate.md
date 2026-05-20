---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
scope: per-figure regeneration brief for the coding agent — cohort fixes (n=6 → n=10), not visual quality
---

# Figures to regenerate — cohort fix queue

**Head.** §2 figures in the current PDF are still rendered at **n=6** (Pat_02, 03, 05, 06, 07, 08), pre-2026-04-25 cohort. The configuration constant `PATIENTS_4PHASE` in `src/lrg_eegfc/config/const.py:134-137` is now n=10 (Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15), and the `_shared.py` module in `scripts/01_compute/figures_for_notes/_shared.py:65` resolves `ALL_PATIENTS` to the same n=10. Re-running each script produces n=10 outputs. §4 figures (`data/audit/raw_fc_phase_distance/`) are already at n=10, no regeneration needed. §5 figures (`data/audit/section5_v2_round3_redo/figures/`) are already at n=10, no regeneration needed. **Action: re-run the §2 generators and update captions** (the manuscript text says "all six patients" / "for all six patients" in several places — those need rewording too once the figures are at n=10).

## §2 figures — regeneration queue (CRITICAL)

All paths are absolute under `/home/giulio/Documents/research/neural_networks/lrgeegfc/`.

### Fig. 1 (Same-shaft enrichment ratio MSC vs |ImCoh|, β rsPre)

- **Manuscript ref**: §2.2.2 introductory figure, p.4 of PDF.
- **Current cohort in PDF**: n=6 (Pat_02, 03, 05, 06, 07, 08).
- **Target cohort**: n=10.
- **Current path**: `data/outputs/figures/section2/fig_01/` (or `fig_01_probe_bias.pdf` at top-level depending on layout).
- **Generator**: `scripts/01_compute/figures_for_notes/fig_01_probe_bias.py`.
- **Action**: re-run the script. Output should now have 10 rows in the heatmap (one per patient, sorted by patient ID).
- **Caption update**: "Same-shaft to cross-shaft mean weight ratio for MSC (left) and `|ImCoh|` (right), computed at the edge-weight level for **all ten** patients..." (was "all six patients").

### Fig. 2 (Community-scale dependence of same-shaft enrichment, β band)

- **Manuscript ref**: §2.2.3 figure on p.5.
- **Current cohort**: n=6.
- **Target**: n=10.
- **Current path**: `data/outputs/figures/section2/fig_D/fig_D1_enrichment_heatmap_redesigned.pdf` and the bias-reduction panels (`fig_D2_bias_reduction_beta_rsPre.pdf` etc.).
- **Generator**: `scripts/01_compute/figures_for_notes/fig_D_enrichment.py` (uses `ALL_PATIENTS` from `_shared.py`).
- **Action**: re-run via `python scripts/01_compute/figures_for_notes/fig_D_enrichment.py --phase rest_pre`.
- **Caption update**: per-scale bias-reduction now over n=10. Cohort mean / SD will shift slightly.

### Fig. 3 (Cross-spectral profiles, all bands, all patients, rsPre)

- **Manuscript ref**: §2.3.1, p.7.
- **Current cohort**: n=6 (panel grid 2x3).
- **Target**: n=10 (will need a 2x5 or 3x4 grid layout).
- **Current path**: `data/outputs/figures/section2/fig_C/fig_C2_spectral_distribution_msc_rsPre.pdf` and `fig_C2_spectral_distribution_imcoh_rsPre.pdf`.
- **Generator**: `scripts/01_compute/figures_for_notes/fig_C2_spectral_distribution.py`.
- **Action**: re-run; the layout function will need to be updated to handle 10 panels (current code was written assuming 6 panels = 2x3 grid). Look for the `subplots(2, 3, ...)` call in the generator and replace with a flexible layout (e.g. `2x5` or `3x4` with masked extra cells).
- **Caption update**: "Cross-spectral profiles for **all ten** patients (rsPre phase). Each panel..." (was "all six patients").

### Fig. 4 (KDE of β-band edge weights, same-shaft vs cross-shaft, all patients, rsPre)

- **Manuscript ref**: §2.3.2, p.8.
- **Current cohort**: n=6.
- **Target**: n=10.
- **Current path**: `data/outputs/figures/section2/fig_C/fig_C3_weight_overlay_all_patients_beta.pdf`.
- **Generator**: `scripts/01_compute/figures_for_notes/fig_C_distributions.py` (or the section3 variant `gen_fig_*` in `figures_for_notes/section3/`).
- **Action**: re-run for all bands of interest (β confirmed, also α / θ / γ_l if used). Color cycle will need 10 distinct colors — verify `PATIENT_COLORS` in `_shared.py:103` uses tab10 which gives 10 distinct colors.
- **Caption update**: "Kernel density estimates of β-band edge weights for same-shaft (dashed) and cross-shaft (solid) pairs, rsPre phase, **all ten** patients overlaid by color."

### Fig. 5 / Fig. 6 (Band-resolved FC matrices Pat_02 / Pat_05)

- **Manuscript ref**: §2.3.3, p.9–10.
- **Current cohort**: per-patient (Pat_02 and Pat_05 individually).
- **Target**: per-patient (no cohort change). **No regeneration needed** — these are single-patient figures and Pat_02, Pat_05 are in both n=6 and n=10 cohorts.
- **Action**: none.

### Fig. 7 (Band-resolved FC matrices for all six patients, β rsPre, MSC vs ImCoh)

- **Manuscript ref**: §2.3.3, p.11.
- **Current cohort**: n=6.
- **Target**: n=10 (will need a 2x5 grid instead of 2x6 = 1x6 + 1x6).
- **Current path**: `data/outputs/figures/section2/fig_B/fig_B_beta_rsPre_MSC_vs_ImCoh.pdf` (and `fig_B_beta_rest_pre_MSC_vs_ImCoh.pdf` from the rest_pre alias).
- **Generator**: `scripts/01_compute/figures_for_notes/fig_AB_adjacency.py`.
- **Action**: re-run; the panel grid will go from 2×6 to 2×10. This may need a multi-row layout (2 rows × 5 columns × 2 methods stacked) for legibility. Alternative: split into two figures (n=5 each) or use a smaller per-patient panel size.
- **Caption update**: "Band-resolved FC matrices for **all ten** patients..." (was "all six patients"). Patient channel counts row will have 10 entries.

### Fig. 8 (Force-directed network embeddings, β rsPre)

- **Manuscript ref**: §2.3.4, p.11.
- **Current cohort**: 3 selected patients (Pat_02, Pat_05, Pat_08).
- **Target**: 3 selected patients (no cohort change). **No regeneration needed for cohort.** But see `14_figures_to_redesign.md` for selection-rationale issues.
- **Current path**: `data/outputs/figures/section2/fig_E/fig_E2_spring_Pat_02_Pat_05_Pat_08_beta_rest_pre.pdf` (and the `_rsPre.pdf` alias).
- **Generator**: `scripts/01_compute/figures_for_notes/fig_E_network.py`.
- **Optional action**: discuss whether the 3-patient choice (Pat_02, 05, 08) should be expanded to include one new-cohort patient (e.g. Pat_10 or Pat_13) for cohort-disclosure transparency. **Default**: no change. The 3-patient layout is for visual contrast, not cohort representation.

### Fig. 9 (Sorted Laplacian eigenvalues for Pat_05 across bands and phases)

- **Manuscript ref**: §2.3.4, p.12.
- **Current cohort**: 1 patient (Pat_05). **No regeneration for cohort.**
- **Current path**: `data/outputs/figures/section2/fig_G/fig_G1_eigenvalue_spectrum_Pat_05_MSC_vs_ImCoh.pdf` (and the `Pat_02`, `Pat_08` siblings).
- **Generator**: `scripts/01_compute/figures_for_notes/fig_G_eigenvalue.py`.
- **Action**: none for cohort. See `14_figures_to_redesign.md` for label-sizing.

### Fig. 10 (Normalized entropy + susceptibility for Pat_08 α band)

- **Manuscript ref**: §2.3.4, p.12.
- **Current cohort**: 1 patient. No cohort regeneration.
- **Action**: none.

## §3 figures — verification only

### Fig. 11 (Dendrograms for 3 patients, MSC vs ImCoh)

- **Manuscript ref**: §3.2.1, p.14.
- **Current cohort**: 3 patients (Pat_02 β, Pat_05 θ, Pat_08 α).
- **Target**: same. **No cohort regeneration needed.**
- **Current path**: `data/outputs/figures/section3/fig_H/fig_H1_dendrograms_MSC_vs_ImCoh_3patients.pdf` and `for_writing_agent/` mirror.
- **Generator**: `scripts/01_compute/figures_for_notes/section3/gen_fig_H_dendrograms.py`.
- **Action**: none for cohort. See `14_figures_to_redesign.md` for leaf-label legibility.

### Fig. 12 (LRG community partitions on glass-brain)

- **Manuscript ref**: §3.2.2, p.15.
- **Current cohort**: 3 patients (Pat_02, 03, 05).
- **Target**: same. **No regeneration**.
- **Generator**: `scripts/01_compute/figures_for_notes/section3/gen_fig_I_brain_communities.py`.

## §4 figures — already n=10, verification only

All figures in `data/audit/raw_fc_phase_distance/` are at n=10. They are produced by `scripts/01_compute/audit/audit_25_raw_fc_phase_distance.py` (or the `audit_28_drift_triangle_null.py` companion).

### Fig. 12 §4 (Four-phase geometry on rank distance dS, per band)

- **Path**: `data/audit/raw_fc_phase_distance/diagnostics.pdf` (multi-panel).
- **Cohort**: n=10. ✓
- **Action**: confirm Pat_03 1024 Hz outlier marking is visible (orange triangle).

### Fig. 13 §4 (Chord-diagram four-phase distance geometry)

- **Path**: separate panel in same script output. Verify it's part of `diagnostics.pdf` or a separate chord-figure output.
- **Cohort**: n=10. ✓

### Fig. 14 §4 (Trace scatter on rank distance dS, per band)

- **Path**: same audit_25 output.
- **Cohort**: n=10. ✓

### Fig. 15 §4 (Per-band T_d distributions across distances)

- **Path**: same.
- **Cohort**: n=10. ✓

### Fig. 16 / Fig. 17 §4 (Convergence dS vs dP / dS vs dF)

- **Path**: same.
- **Cohort**: n=10. ✓

### Fig. 18 §4 (Symmetric cross-baseline scatter on rank dS)

- **Path**: `data/audit/raw_fc_phase_distance/cross_baseline_scatter.pdf`.
- **Cohort**: n=10. ✓

### Fig. 19 §4 (Two phase-pair examples, Pat_14 low-γ + Pat_02 high-γ)

- **Path**: per-patient subdirs under `data/audit/raw_fc_phase_distance/Pat_*/`.
- **Cohort**: example panels — no cohort change.

## §5 figures — already n=10, verification only

All figures in `data/audit/section5_v2_round3_redo/figures/` and the kc_leaf / ctm_per_pair sibling directories are at n=10.

### Fig. 20 (Ψ profile on |ImCoh| substrate at Pat_02 β rsPre)

- **Path**: `data/audit/section5_v2_round3_redo/figures/psi_irrelevance.pdf`.
- **Cohort**: 1 example. ✓
- **Generator**: `scripts/01_compute/audit/audit_round3_section5_redo.py` (psi_irrelevance helper).

### Fig. 21 (Ψ argmax position histogram across cohort)

- **Path**: `data/audit/section5_v2_round3_redo/figures/psi_argmax_histogram.pdf`.
- **Cohort**: n=10 × 6 bands × 3 phases = 180 cells. ✓

### Fig. 22 (Geometric decomposition of KC tree-distance, per band)

- **Path**: `data/audit/section5_v2_round3_redo/figures/kc_decomposition.pdf`.
- **Cohort**: n=10. ✓

### Fig. 23 (Per-leaf decomposition of KC λ=0 topology trace)

- **Path**: `data/audit/section5_v2_round3_redo/figures/kc_topology_per_pair_signed.pdf` (or `kc_topology_eye_proof.pdf` for the alternative composite).
- **Cohort**: 4 patients × 3 phases. ✓

### Fig. 24 (Controlled cohort ρ_split at the LRG ultrametric layer)

- **Path**: `data/audit/section5_v2_round3_redo/figures/ctm_headline.pdf`.
- **Cohort**: n=10. ✓

### Fig. 25 (Per-pair scatter of Δ_task vs Δ_rest at three exemplar cells)

- **Path**: `data/audit/ctm_per_pair_scatter/figures/manuscript_ctm_per_pair_scatter.pdf`.
- **Cohort**: 3 exemplar cells. ✓

### Fig. 26 (Multiscale partition trace via VI(k))

- **Path**: `data/audit/section5_v2_round3_redo/figures/` (filename to be confirmed; cohort heatmap + per-patient strip).
- **Cohort**: n=10. ✓

### Fig. 27 (Multiscale spectral trace via Grassmann principal angles)

- **Path**: `data/audit/section5_v2_round3_redo/figures/grassmann_principal_angles.pdf`.
- **Cohort**: n=10. ✓

### Fig. 28 (Anatomical localization of CTM trace-leaves)

- **Path**: `data/audit/anatomy_2d_disclosed/figures/` (from audit_50_anatomy_55_figure.py + audit_51_anatomy_2d_disclosed.py).
- **Cohort**: n=10. ✓

### Fig. 29–34 (§5.6 taxonomy class composites)

- **Path**: `data/audit/section5_v2_round3_redo/figures/taxonomy_class_*.pdf` (per-class) + `taxonomy_pat*_beta_composite.pdf` (per-patient β).
- **Cohort**: n=10 (per-patient panel). ✓
- **Generators**: `scripts/01_compute/audit/audit_51a_taxonomy_cohort_stack.py` + `audit_52_section5_class_showcases.py`.

## Caption rewrites needed when regenerated figures land

For each §2 figure that goes from n=6 to n=10, search the .tex source for the phrase "all six patients" / "for all six patients" / "across all six patients" / "across the six patients" / "the six patients" and rewrite to "all ten patients" / "for all ten patients" / "across all ten patients" / "across the ten patients" / "the ten patients". The phrase appears at:
- Fig. 1 caption
- Fig. 2 caption (same-shaft enrichment, β band)
- Fig. 3 caption
- Fig. 4 caption
- Fig. 7 caption (β rsPre, MSC vs ImCoh)
- §2.3.1 spectral profiles prose paragraph

Cohort-mean numerical values (e.g. "cohort mean 3.2× under MSC", "cohort mean 0.99× under |ImCoh|") will shift slightly when going to n=10. The qualitative pattern (MSC ~ 1.4×–6× same-shaft enrichment; |ImCoh| ~ 1× indistinguishable) holds at n=10 per the n=10 data already in `data/audit/raw_fc_phase_distance/cohort_summary.csv` and the §2.2.2 prose paragraph; only the headline numbers should be re-cited from the regenerated figures.

## Action — coding agent

Order of operations:

1. Re-run §2 generators (`fig_01_probe_bias.py`, `fig_D_enrichment.py`, `fig_C2_spectral_distribution.py`, `fig_C_distributions.py`, `fig_AB_adjacency.py`).
2. Verify the panel-grid layout for n=10 (current code may break on 2x3 → needs 2x5 or 3x4).
3. Update `_shared.py:65` to drop the redundant `| {"Pat_06"}` (since PATIENTS_4PHASE already contains Pat_06):
   ```python
   ALL_PATIENTS = list(PATIENTS_4PHASE)  # n=10
   ```
4. Verify Pat_03 1024 Hz outlier marking is preserved (orange triangle / distinct legend color).
5. Re-render the captions where "all six" → "all ten".
6. Confirm output cohort matches by reading row count from regenerated PDFs (or re-run the underlying CSV to verify n=10 patient rows).
