---
name: talk-slide-11-pipeline
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 11
status: draft
updated: 2026-07-11
canva: page 11 · ~55%
---

# Slide 11 — The network-analysis pipeline

1. TITLE
The network-analysis pipeline

2. MAIN CONCEPT
- Recap in one panel (slides 7 & 9, shown by the figure sequence, NOT re-derived): timeseries → per-band |ImCoh| FC → diffusion → cophenetic dendrogram (the deterministic linkage ℒ of the communication distance D(τ)).
- The one genuinely new step is the cross-phase comparison: rank-correlate two phases' cophenetic distances to ask how similar their hierarchies are. Illustrated by two trees with a node that kept its place.
- The trace schema uses the three phases (rest_pre → task → rest_post): a trace = the task-induced reorganization that persists into rest_post. The estimator is ρ — a Spearman of the cophenetic displacements, symmetrized over a split rest_pre baseline.

3. ON-SLIDE TEXT
Bullet list (on-slide):
- recap (slides 7 & 9): timeseries → |ImCoh| → diffusion → dendrogram (linkage ℒ)
- compare two phases: rank-correlate their cophenetic distances
- a trace across rest_pre → task → rest_post: does the task change persist? → ρ

Formulas — compile in Canva's LaTeX tool. ℒ = the deterministic (UPGMA / average-linkage) hierarchical clustering; D(τ) = the communication distance (slide 7); A,B = the two split-halves of rest_pre.

(1) cophenetic distance = merge-height of the lowest common ancestor in the average-linkage tree
D^{\mathrm{coph}}_{ij}=h_{\mathcal{L}[D]}\!\bigl(\mathrm{LCA}(i,j)\bigr)

(ℒ[D] = the UPGMA / average-linkage tree of the distance D — ℒ already shown; h(·) = merge height of a node; LCA(i,j) = lowest common ancestor of leaves i, j. So D^coph_ij = the height at which i and j first join in the dendrogram.)

(2) trace estimator across rest_pre → task → rest_post (A, B = the two split-halves of rest_pre)
\rho=\tfrac{1}{2}\Bigl[\rho_{S}\bigl(\Delta^{\mathrm{task}}_{A},\,\Delta^{\mathrm{rest}}_{B}\bigr)+(A\leftrightarrow B)\Bigr]

\Delta^{\mathrm{task}}_{X}=D^{\mathrm{coph}}_{\mathrm{task}}-D^{\mathrm{coph}}_{X},\qquad
\Delta^{\mathrm{rest}}_{X}=D^{\mathrm{coph}}_{\mathrm{post}}-D^{\mathrm{coph}}_{X},\qquad X\in\{A,B\}

((A↔B) = the same term with A and B swapped; the ½ averages the two arm assignments of the split rest_pre baseline — removes the arbitrary A/B labeling. Δtask = how the task moved the hierarchy; Δrest = how rest_post moved it; ρ asks whether the two agree.)

4. SPEECH
Here's the whole machine on one slide — and most of it you've already seen. From the timeseries we build the per-band imaginary-coherence network, run diffusion on it, and cluster the resulting communication distances into a dendrogram — the multiscale hierarchy, through a deterministic linkage. What's new here is the comparison. To compare two phases we simply rank-correlate their cophenetic distances — how deep in the tree each pair of contacts meets — which tells us how similar two hierarchies are; picture two trees with a node that kept its place. And to detect a trace we use all three phases: does the reorganization the task induced still show up at rest afterwards? Our estimator, ρ, is exactly that — a rank correlation between the task change and the rest change.

Careful: no τ-sweep — the multiscale content is the dendrogram's own hierarchy read at τ = 1/λmax (don't imply we scan τ). ρ is the trace estimator (three phases: task change vs rest change). Present the estimator here, NOT a verdict — the control that turns ρ into a verdict comes on slides 14–15 (kept off this slide on purpose).

5. FIGURES
- THE PIPELINE SEQUENCE — 5 clean single panels (Pat_05 β), compose left→right with arrows in Canva. BUILT 2026-07-11, gen `scripts/07_figures/gen_pipeline_sequence.py` (`--patient Pat_05 --band beta`):
  1. example timeseries — data/outputs/figures/talk/pipeline_seq_1_timeseries.png (a few sEEG contacts, ~2 s).
  2. |ImCoh| FC matrix — data/outputs/figures/talk/pipeline_seq_2_imcoh_matrix.png (generic FCᵢⱼ colorbar).
  3. communication distance D(τ)=(1−δ)/K̂(τ) — data/outputs/figures/talk/pipeline_seq_3_distance.png (the matrix that FEEDS the clustering; photographic negative of the FC — strong coupling = small distance; heatmap at τ=1/λmax, log colour). Shows D not K on purpose: at τ_min K̂ is ~monotone in the FC (Spearman 0.98) so K would read as a near-copy of panel 2; D is the distinct object UPGMA runs on.
  4. cophenetic dendrogram — data/outputs/figures/talk/pipeline_seq_4_dendrogram.png (UPGMA tree, log-height).
  5. cross-phase comparison — two dendrograms (rest_pre vs rest_post) with a RETAINED clade highlighted (5 parallel amber ribbons = kept its place; grey crossing ribbons = reorganised); the visual intuition for ρ (a pair that keeps its place = agreement between phases). — data/outputs/figures/talk/pipeline_seq_5_rhocoph_tanglegram.png. NOTE: the baked annotation is a pairwise rest_pre↔rest_post rank agreement (illustrative only); the slide's single estimator symbol is ρ (the symmetric three-phase form) — do NOT print a separate ρ^coph symbol on-slide.
- (Alternative single-image pipeline) the 3D HERO — data/outputs/figures/talk/pipeline_brain_tree_3d.png. NOTE: shared with slide 08 (the bet). If slide 08 uses it as the vision/teaser, slide 11 should LEAD with the step-by-step SEQUENCE above (the build), not the same still.
- The two formulas render as on-slide LaTeX (compiled in Canva), not Lane-F PDFs.

6. REFERENCES
- LRG operator + communication distance + UPGMA dendrogram: Villegas, P., Gili, T., Caldarelli, G., Gabrielli, A. (2023), "Laplacian renormalization group for heterogeneous networks", Nature Physics 19, 445–450; Villegas, P., Gabrielli, A., Poggialini, A., Gili, T. (2025), Physical Review Research 7, 013065.
- |ImCoh| FC: Nolte et al. 2004 / Ewald et al. 2012 (introduced slide 9).

7. CANVA STATUS
Page 11 · ~55% (content locked 2026-07-11; figures BUILT). Have: title; the recap concept; TWO compile-ready LaTeX formulas (D^coph via ℒ · the trace estimator ρ) — minimal, matched to code, NO re-derived ladder (it's on slides 7/9), NO null on this slide (deferred to 15–16 per user); speech; the 5-panel pipeline sequence (Pat_05 β) gathered in talk/. ⚠ 2026-07-11: CORRECTED + TRIMMED the formulas — RE-RENDER in Canva. (1) D^coph was `ℒ[D(τ)]` (type-wrong: ℒ builds a TREE, not a distance matrix); now `D^coph_ij = h_{ℒ[D]}(LCA(i,j))` — the LCA merge-height in the average-linkage tree ℒ[D], matching scipy `cophenet` (`tree.py` cophenet_matrix). Single line, ℒ reused (already shown), τ dropped. (2) trace estimator: DROPPED the standalone non-symmetric ρ^coph formula (per user); the estimator is now just `ρ` (no `_sym` subscript) in the ½-averaged split-half form matching audit_150 `rho_sym_split`, with the second arm written `(A↔B)` instead of spelled out; Spearman denoted `ρ_S`. ρ^coph symbol removed from on-slide text/concept/speech too. ⚠ Fixed the name/header 10→11 numbering lag. ⚠ 3D hero shared with slide 08 → here lead with the step-by-step SEQUENCE (build), not the same still. Remaining to 100%: presenter composes the 5-panel sequence + 3 formulas in Canva.
