---
name: talk-slide-11-pipeline
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 11
status: draft
updated: 2026-07-13
canva: page 11 · ~55% (⚠ content updated 2026-07-13 — sparsification + τ-sweep added; RE-DO on deck)
---

# Slide 11 — The network-analysis pipeline

1. TITLE
The network-analysis pipeline

2. MAIN CONCEPT
- Recap in one panel (slides 5, 7 & 9, shown by the figure sequence, NOT re-derived): timeseries → per-band |ImCoh| FC → **sparsify to a backbone** → **diffusion swept over τ** → cophenetic dendrogram (the deterministic linkage ℒ of the communication distance D(τ)).
- Two steps make the multiscale real (both are what changed since the early draft — say them plainly):
  - **Sparsification** — keep each network's maximum spanning tree plus its strongest ~20% of edges (`mst@0.20`). On the fully-connected |ImCoh| graph the propagator is *degenerate* at the fine scale (the heat geometry is ≈ the raw weights), so diffusion buys nothing; the backbone makes it non-degenerate and gives τ something to resolve. β is robust to *which* backbone; α needs one that keeps the strong edges — that's why MST-plus-strong-edges, not the parameter-free percolation skeleton.
  - **τ is swept, not fixed** — diffusion time τ is the scale knob (slide 7). Sweeping it turns the one dendrogram into a multiscale *scanner*: the β reorganization is scale-invariant (clears at every τ), α is mesoscale-tuned. The band-specificity is *in* that scale dependence — read at a single τ it isn't there.
- The one genuinely new comparison step: rank-correlate two phases' cophenetic distances to ask how similar their hierarchies are. Illustrated by two trees with a node that kept its place.
- The trace schema uses the three phases (rest_pre → task → rest_post): a trace = the task-induced reorganization that persists into rest_post. The estimator is ρ — a Spearman of the cophenetic displacements, symmetrized over a split rest_pre baseline, read across τ.

3. ON-SLIDE TEXT
Bullet list (on-slide):
- recap: timeseries → |ImCoh| → **sparsify (MST + top-20%)** → **diffusion e^{−τL̂}** → dendrogram (linkage ℒ)
- **τ swept = multiscale scanner** (β at every scale · α at the mesoscale)
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
Here's the whole machine on one slide — and most of it you've already seen. From the timeseries we build the per-band imaginary-coherence network, then thin it to its backbone — the spanning tree plus its strongest fifth of edges — which is exactly what lets diffusion see structure the dense graph hides. We run diffusion on that backbone and *sweep the diffusion time* to scan every scale, and cluster the resulting communication distances into a dendrogram — the multiscale hierarchy, through a deterministic linkage. What's new here is the comparison. To compare two phases we simply rank-correlate their cophenetic distances — how deep in the tree each pair of contacts meets — which tells us how similar two hierarchies are; picture two trees with a node that kept its place. And to detect a trace we use all three phases: does the reorganization the task induced still show up at rest afterwards? Our estimator, ρ, is exactly that — a rank correlation between the task change and the rest change, read across scales.

Careful: τ IS swept now — it is the scale scanner, and the multiscale content lives in the scale-dependence (β clears at every scale, α at the mesoscale). The earlier "read at a single τ = 1/λmax, no sweep" framing is RETIRED: on the dense fully-connected graph that operating point was degenerate (heat geometry ≈ raw weights), which is why the sparsification step is now real, not cosmetic — the MST-plus-strongest-20% backbone is what makes the propagator non-degenerate. ρ is the trace estimator (three phases: task change vs rest change). Present the estimator here, NOT a verdict — the matched-strength control that turns ρ into a verdict comes on slides 14–15 (kept off this slide on purpose). Matched-strength is now the ONLY null (drift retired — a directional task makes a drift null degenerate with the trace).

5. FIGURES
⚠ PIPELINE UPDATE 2026-07-13 (figures need a rebuild — deferred to the figure lane / task #13): the sequence must now show (a) a **sparsification panel** — the |ImCoh| matrix → its `mst@0.20` backbone (dense heatmap → sparse backbone / MST+top-20% overlay), inserted between panels 2 and 3; and (b) the diffusion step as a **τ-sweep**, not a single τ (a short D(τ) morph or 2–3 τ snapshots), so the visual matches "τ = the scanner". Until rebuilt, the 5-panel sequence below is the OLD dense/single-τ build.
- THE PIPELINE SEQUENCE — 5 clean single panels (Pat_05 β), compose left→right with arrows in Canva. BUILT 2026-07-11, gen `scripts/07_figures/gen_pipeline_sequence.py` (`--patient Pat_05 --band beta`):
  1. example timeseries — data/outputs/figures/talk/pipeline_seq_1_timeseries.png (a few sEEG contacts, ~2 s).
  2. |ImCoh| FC matrix — data/outputs/figures/talk/pipeline_seq_2_imcoh_matrix.png (generic FCᵢⱼ colorbar).
  [NEW] 2b. sparsified backbone — the |ImCoh| graph thinned to `mst@0.20` (MST + strongest 20%); the step that makes the propagator non-degenerate. (to build)
  3. communication distance D(τ)=(1−δ)/K̂(τ) — data/outputs/figures/talk/pipeline_seq_3_distance.png (the matrix that FEEDS the clustering; photographic negative of the FC — strong coupling = small distance; heatmap at τ=1/λmax, log colour). Shows D not K on purpose: at τ_min K̂ is ~monotone in the FC (Spearman 0.98) so K would read as a near-copy of panel 2; D is the distinct object UPGMA runs on. NOW: show as a τ-sweep, not a single τ.
  4. cophenetic dendrogram — data/outputs/figures/talk/pipeline_seq_4_dendrogram.png (UPGMA tree, log-height).
  5. cross-phase comparison — two dendrograms (rest_pre vs rest_post) with a RETAINED clade highlighted (5 parallel amber ribbons = kept its place; grey crossing ribbons = reorganised); the visual intuition for ρ (a pair that keeps its place = agreement between phases). — data/outputs/figures/talk/pipeline_seq_5_rhocoph_tanglegram.png. NOTE: the baked annotation is a pairwise rest_pre↔rest_post rank agreement (illustrative only); the slide's single estimator symbol is ρ (the symmetric three-phase form) — do NOT print a separate ρ^coph symbol on-slide.
- (Alternative single-image pipeline) the 3D HERO — data/outputs/figures/talk/pipeline_brain_tree_3d.png. NOTE: shared with slide 08 (the bet). If slide 08 uses it as the vision/teaser, slide 11 should LEAD with the step-by-step SEQUENCE above (the build), not the same still.
- The two formulas render as on-slide LaTeX (compiled in Canva), not Lane-F PDFs.

6. REFERENCES
- LRG operator + communication distance + UPGMA dendrogram: Villegas, P., Gili, T., Caldarelli, G., Gabrielli, A. (2023), "Laplacian renormalization group for heterogeneous networks", Nature Physics 19, 445–450; Villegas, P., Gabrielli, A., Poggialini, A., Gili, T. (2025), Physical Review Research 7, 013065.
- |ImCoh| FC: Nolte et al. 2004 / Ewald et al. 2012 (introduced slide 9).

7. CANVA STATUS
⚠ 2026-07-13 PIPELINE UPDATE (mst@0.20 recovery): the slide now carries TWO new pipeline steps — **sparsification** (MST + strongest 20%, because the dense propagator is degenerate at the fine scale) and **τ swept as the multiscale scanner** (β scale-invariant, α mesoscale). The old "no τ-sweep, read at τ=1/λmax" framing is RETIRED throughout (concept/on-slide/speech). Matched-strength is now the ONLY null (drift retired). On-slide bullets + speech updated in-repo; figure sequence needs a backbone panel + τ-sweep morph (task #13). RE-RENDER the on-slide bullets in Canva.
Page 11 · ~55% (content locked 2026-07-11; figures BUILT). Have: title; the recap concept; TWO compile-ready LaTeX formulas (D^coph via ℒ · the trace estimator ρ) — minimal, matched to code, NO re-derived ladder (it's on slides 7/9), NO null on this slide (deferred to 15–16 per user); speech; the 5-panel pipeline sequence (Pat_05 β) gathered in talk/. ⚠ 2026-07-11: CORRECTED + TRIMMED the formulas — RE-RENDER in Canva. (1) D^coph was `ℒ[D(τ)]` (type-wrong: ℒ builds a TREE, not a distance matrix); now `D^coph_ij = h_{ℒ[D]}(LCA(i,j))` — the LCA merge-height in the average-linkage tree ℒ[D], matching scipy `cophenet` (`tree.py` cophenet_matrix). Single line, ℒ reused (already shown), τ dropped. (2) trace estimator: DROPPED the standalone non-symmetric ρ^coph formula (per user); the estimator is now just `ρ` (no `_sym` subscript) in the ½-averaged split-half form matching audit_150 `rho_sym_split`, with the second arm written `(A↔B)` instead of spelled out; Spearman denoted `ρ_S`. ρ^coph symbol removed from on-slide text/concept/speech too. ⚠ Fixed the name/header 10→11 numbering lag. ⚠ 3D hero shared with slide 08 → here lead with the step-by-step SEQUENCE (build), not the same still. Remaining to 100%: presenter composes the 5-panel sequence + 3 formulas in Canva.
