---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
scope: visual-quality redesign queue (NOT cohort fixes); referee-audit issues collected from the 260508 PDF
---

# Figures to redesign — visual-quality queue

**Head.** Issues here are **visual quality** — labels, fonts, colorbars, overplotting — not cohort or numerical content. These complement `13_figures_to_regenerate.md`. None of the items below changes the underlying data; all are edits to the generator scripts to improve readability for a referee.

## Taxonomy class composites — missing in-figure patient/band labels

**Affected figures.** §5.6, Figs. 29 / 30 / 31 / 32 / 33 / 34. Files `taxonomy_class_trace.pdf`, `taxonomy_class_reset.pdf`, `taxonomy_class_rearrange.pdf`, `taxonomy_class_anchor.pdf`, `taxonomy_class_diffuse.pdf`, plus the full Pat_13 β taxonomy composite `taxonomy_pat13_beta_composite.pdf` shown as Fig. 34 in the manuscript (and the per-patient siblings `taxonomy_pat02_beta_composite.pdf` etc.).

**Path.** `data/audit/section5_v2_round3_redo/figures/`.

**Issue.** The exemplar (patient, band) for each class figure is currently noted only in the figure caption text ("Fig. 29 shows trace at Pat_07 β"). Within the figure itself, the panel headers say "rest pre / task / rest post" but do not label the patient or band. A reader looking at the figure standalone (e.g. a slide deck) cannot tell which patient/band is shown.

**Fix.** Add an in-figure label in the upper-left corner of the leftmost panel (or as a top-strip) of each class composite:
- `taxonomy_class_trace.pdf` → "Pat_07 · β · 6 modules · 32/116 class leaves" (already shown in the top-right of the existing PDF, verify it survives the redesign).
- Pat_07 β composite already has this — confirm consistency for other class figures.

**Generator.** `scripts/01_compute/audit/audit_52_section5_class_showcases.py` and `audit_51a_taxonomy_cohort_stack.py`.

**Action.** Add an in-figure annotation block (likely at `ax.set_title` or `ax.text` for the first panel of each row) that contains:
- patient ID (e.g. "Pat_07")
- band (e.g. "β")
- module count
- class-leaf count and total leaf count

Currently the existing PDF headers DO have these (Fig. 29 PDF shows "Pat_07 beta 6 modules · 32/116 class leaves" in upper right). Confirm this is preserved in the regeneration. If missing in any of the 5 class composites, add it.

## Dendrogram leaf labels — unreadable at fig_H1 scale

**Affected figure.** Fig. 11. File `fig_H1_dendrograms_MSC_vs_ImCoh_3patients.pdf`.

**Path.** `data/outputs/figures/section3/fig_H/fig_H1_dendrograms_MSC_vs_ImCoh_3patients.pdf` and `for_writing_agent/` mirror.

**Issue.** Dendrograms have ~117 leaves (Pat_02), ~118 (Pat_05), ~120 (Pat_08). Leaf labels at the bottom are typeset at ~3pt and overlap. They are colored by shaft identity (mentioned in caption) but not legible at default print scale. The caption says "Leaves are coloured by shaft identity" — true, but the legend mapping shaft → color is also missing.

**Fix.** Two options:
1. Drop leaf labels entirely; rely on the colored leaf bars at the dendrogram bottom (already present). Add a small colored-bar legend on the right side mapping shaft prefix → color.
2. Keep leaf labels but rotate 90°, decimate (every 5th label), and increase font size to ~6pt.

Option 1 is cleaner for a publication figure. Option 2 is appropriate for a supplementary figure.

**Generator.** `scripts/01_compute/figures_for_notes/section3/gen_fig_H_dendrograms.py`.

**Action.** Apply Option 1. Verify the legend (shaft prefix list + color swatch) does not occlude the dendrogram axis.

## KC heatmap — colorbar label missing

**Affected figure.** Fig. 22. File `kc_decomposition.pdf`.

**Path.** `data/audit/section5_v2_round3_redo/figures/kc_decomposition.pdf`.

**Issue.** The right panel ("per-band imprint counts n_imprint/10 on each axis") is bar chart — no colorbar issue there. But the left panel scatter on the (`-T_KC(λ=0)`, `-T_KC(λ=1)`) plane uses red/grey markers without a colorbar key in the panel itself. The markers code `p ≤ 0.05` (red filled) vs `p > 0.05` (white filled with marker outline, grey edge), per the legend at top-right. Legend is present but small. **Colorbar issue is for the per-band axis-strip below**: it has a green/red imprint-direction encoding but no axis-strip colorbar label is visible.

**Wait** — re-checking Fig. 22 in the PDF (`p.28`): the figure has two panels (left = scatter; right = bar chart). The colorbar/legend labeling looks OK in the existing PDF. **Possibly false positive.** Confirm with the writing agent before action.

**If confirmed**: add a small axis-strip color key (red = trace direction p ≤ 0.05; green = anti-trace; grey = ns) under the scatter panel, or include in the existing top-right legend.

**Generator.** `scripts/01_compute/audit/audit_46_kc_section5_controls.py` (or audit_round3_section5_redo).

## CTM headline — colorbar label missing

**Affected figure.** Fig. 24. File `ctm_headline.pdf`.

**Path.** `data/audit/section5_v2_round3_redo/figures/ctm_headline.pdf`.

**Issue.** The legend at the bottom shows `ρ_split` (controlled, blue), `ρ_drift` (grey), `ρ_xprobe` (orange) — three boxes per band. Per-patient lines are colored green (`ρ_split > 0`) or red (anti). The Wilcoxon p-values are annotated above the blue box per band. **Colorbar label issue**: the green/red patient-line color does not have a key inside the figure (only mentioned in caption: "colored green when `ρ_split > 0` (trace direction) and red otherwise"). For a referee, the red lines are the visible-but-unkeyed signal.

**Fix.** Add a small inset legend in the upper-right or upper-left of the figure stating "patient: pro (`ρ_split > 0`) — green; patient: anti (`ρ_split ≤ 0`) — red" with green and red line swatches.

**Generator.** `scripts/01_compute/audit/audit_round3_section5_redo.py` (CTM headline panel).

**Action.** Add inset legend.

## VI smallmult — inconsistent header fonts

**Affected figure.** Fig. 26. File `vi_smallmult_polished.pdf`.

**Path.** `data/audit/section5_v2_round2/figures/vi_smallmult_polished.pdf`.

**Issue.** The cohort-heatmap top panel and the per-patient bottom panels use different header font sizes for the band labels and patient IDs. The patient IDs in the bottom panels are typeset at ~7pt (legible) but the band labels in the top panel are at ~10pt and bolded. The visual asymmetry suggests two different rendering passes were stitched together.

**Fix.** Normalize header font size and weight across the cohort heatmap and the per-patient panels:
- Cohort heatmap top panel band labels: 8pt regular.
- Per-patient panel patient IDs: 8pt regular.
- Per-patient panel band labels (vertical axis): 7pt regular.

**Generator.** `scripts/01_compute/audit/audit_round2_section5_v2_fixes.py` or the `vi_smallmult_polished` helper in audit_round2_section5_v2.py.

**Action.** Normalize header / axis-label font sizes and weights.

## Manuscript KC leaf topology — heatmap fonts unreadable

**Affected figure.** §5.5 (anatomy figure or §5.2 KC composite — manuscript reference unclear). File `manuscript_kc_leaf_topology.pdf`.

**Path.** `data/audit/kc_leaf_topological_memory/figures/manuscript_kc_leaf_topology.pdf`.

**Issue.** The composite layout has multiple sub-panels (heatmap + dendrogram + leaf score curve). Heatmap row labels (anatomical region names like "ctx-lh-fusiform") are typeset at ~5pt and unreadable at print scale.

**Fix.**
- Increase heatmap row-label font size to ≥7pt.
- Truncate Desikan-Killiany region names to short forms ("L fus", "L par opc", etc.) and add a footnote with the full names if space requires.

**Generator.** `scripts/01_compute/audit/audit_53b_manuscript_composite.py` or `audit_53_kc_leaf_topological_memory.py`.

**Action.** Apply font-size / region-label-truncation fixes.

## Manuscript CTM per-pair scatter — dense overplotting

**Affected figure.** Fig. 25. File `manuscript_ctm_per_pair_scatter.pdf`.

**Path.** `data/audit/ctm_per_pair_scatter/figures/manuscript_ctm_per_pair_scatter.pdf`.

**Issue.** Each of the three panels (Pat_06 β, Pat_05 α, Pat_08 θ) contains ~6500–7140 contact pairs as scatter points. With same-probe pairs in light gray and cross-probe pairs in dark blue, the plotting area is mostly black with a faint diagonal pull at α and β. The diagonal pull is the visible signal, but the black saturation makes it hard to see the cohort-trend overlay.

**Fix.**
- Switch from scatter to **hexbin** (`matplotlib.pyplot.hexbin`) with a dual-color encoding (same-probe: log-density gray; cross-probe: log-density blue, overlaid).
- Or reduce alpha to 0.15 for the dense scatter, retain the marker for ≤ 1000 randomly-sampled cross-probe pairs at full alpha for visibility.
- Add the `y = x` identity line in red (currently dashed gray, low contrast against the black saturation).

**Generator.** `scripts/01_compute/audit/audit_55_ctm_per_pair_scatter.py`.

**Action.** Implement hexbin variant; verify the trace direction (diagonal pull along `y = x`) remains visible after the density transform.

## Eigenvalue spectrum (Pat_05 / Pat_02 / Pat_08) — inconsistent label sizing

**Affected figure.** Fig. 9 + supplementary siblings. Files `fig_G1_eigenvalue_spectrum_Pat_*_MSC_vs_ImCoh.pdf`.

**Path.** `data/outputs/figures/section2/fig_G/`.

**Issue.** Across the three sibling files (Pat_05, Pat_02, Pat_08), x-axis tick labels and y-axis labels render at slightly different point sizes. The triangular markers (`λ_2` red, `λ_max` grey) are inconsistent in marker size between panels. This looks like the figures were rendered at different times with different rcParams.

**Fix.**
- Set rcParams once in the generator before plotting (use `apply_pub_style()` from `_shared.py:117`).
- Verify all three sibling files render with identical font sizes and marker sizes.
- Confirm the triangular markers (`λ_2`, `λ_max`) are at the same point size across the three patient panels.

**Generator.** `scripts/01_compute/figures_for_notes/fig_G_eigenvalue.py`.

**Action.** Run `apply_pub_style()` at the top of the script; re-render the three patient siblings; visually confirm consistency.

## Force-directed network embeddings — same-probe color saturation

**Affected figure.** Fig. 8. File `fig_E2_spring_Pat_02_Pat_05_Pat_08_beta_rsPre.pdf` (and `_rest_pre.pdf` alias).

**Path.** `data/outputs/figures/section2/fig_E/`.

**Issue.** Same-shaft edges are colored by shaft identity using the matplotlib tab10 cycle, which has limited contrast for some pairs (e.g. orange tab10[1] vs pink tab10[6] are both in the warm half). With ~10–20 shafts per patient, visually-distinct shaft coloring runs out around shaft 8.

**Fix.**
- Use a perceptually-uniform sequential or qualitative palette with ≥ 20 distinct hues (e.g. matplotlib's `glasbey` or `seaborn.color_palette("tab20")`).
- Or use a single color for same-shaft edges + intensity encoding for shaft index.
- Cross-shaft edges in black is fine (current).

**Generator.** `scripts/01_compute/figures_for_notes/fig_E_network.py`.

**Action.** Switch palette; verify all patient panels render with distinct shaft colors.

## §4 figures — diagnostic plots in the public PDF

**Affected figure.** Fig. 19 §4 (the two phase-pair examples).

**Path.** `data/audit/raw_fc_phase_distance/Pat_*/report.pdf` (per-patient) — but the publication version in the manuscript is composed differently.

**Issue.** The §4 Fig. 19 in the PDF shows two rows: top = Pat_14 low-γ TL→RPost (topological reorganization); bottom = Pat_02 high-γ TT→RPost (magnitude redistribution). Each row has 5 panels: matrix `A^φ`, matrix `A^φ'`, log-log scatter, network `A^TL`, network `A^RPost`. Font sizes within the matrix axes (probe labels) are too small to read. Probe color legend below ("node = electrode (color = sEEG probe), cross-probe edge (gray, weight-scaled), same-probe edge (probe color, full alpha), identity, power-law fit") is one line of dense text.

**Fix.**
- Drop the probe labels entirely from the matrix axes (axes are anonymous `i`, `j` per the never/always rule on adjacency labels).
- Reformat the legend as a multi-line block with one item per line.

**Generator.** `scripts/01_compute/audit/audit_25_raw_fc_phase_distance.py` or its companion.

**Action.** Apply the relabel and legend reformat.

## Cross-probe signs (Fig. supplementary in §5)

**Affected figure.** Supplementary file `cross_probe_signs.pdf`.

**Path.** `data/audit/section5_v2_round3_redo/figures/cross_probe_signs.pdf`.

**Issue.** Used in §5.3 prose to show the cross-probe restriction outcome. Verify the figure exists and is referenced correctly in the manuscript text. If not currently a numbered manuscript figure, it can stay as a supplementary file.

**Action.** Confirm whether this figure is supposed to be numbered in the manuscript. If not, no action.

## Action priority

1. **High** — Dendrogram leaf labels in Fig. H1 (publication legibility).
2. **High** — CTM per-pair scatter overplotting (referee will look at this hard).
3. **High** — Manuscript KC leaf topology heatmap fonts.
4. **Medium** — Taxonomy class composites in-figure labels.
5. **Medium** — VI smallmult header font normalization.
6. **Medium** — Eigenvalue spectrum label sizing across siblings.
7. **Low** — Force-directed embedding palette.
8. **Low** — §4 Fig. 19 diagnostic relabel.
9. **Low** — KC heatmap / CTM headline colorbar / legend additions (probable false-positives — confirm first).

## Action — coding agent

For each redesign item:
1. Read the corresponding generator script.
2. Apply the fix (font-size, palette, legend, etc.).
3. Re-render to the same path.
4. Verify visually.
5. Note any cohort dependencies (most redesign items are cohort-agnostic, but a few — e.g. the spring layout — depend on n=3 vs n=10 patient choice).

Do not couple cohort-regeneration (`13_figures_to_regenerate.md`) and visual-redesign passes in the same run unless the affected figure is in both queues. Run cohort-regeneration first; then redesign on the n=10 outputs.
