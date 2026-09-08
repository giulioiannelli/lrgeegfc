---
name: talk-slide-11-pipeline
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym · mst@0.20 backbone)
slide: 11
status: draft
updated: 2026-07-16
canva: page 11 · ~55% (⚠ figures rebuilt 2026-07-13 — 7-panel sparsification sequence; RE-COMPOSE in Canva)
---

# Slide 11 — The network-analysis pipeline

1. TITLE
The network-analysis pipeline

2. MAIN CONCEPT
- One example patient, one band (Pat_05 β), left→right through the whole machine: brain electrodes → sEEG timecourses → |ImCoh| connectivity → **sparse backbone** → diffusion / specific heat → cophenetic tree → the trace estimator ρ.
- The genuinely NEW step is the **sparsification** (`mst@0.20`). The |ImCoh| graph is fully connected — every contact couples to every other — and a fully-connected graph has essentially **one scale**: the propagator is degenerate at the fine scale (its heat geometry ≈ the raw weights), so diffusion buys nothing and C(τ) collapses to a single peak. We keep each network's maximum spanning tree plus its strongest ~20 % of edges: a connected, cycle-rich backbone. On THAT graph the propagator e^{−τL̂} telescopes and C(τ) grows a mesoscale **ladder** — the multiple scales the LRG is built to read. (Speaker context: β is robust to *which* backbone; α needs one that keeps the strong edges — hence MST-plus-strong-edges, not the parameter-free percolation skeleton.)
- **τ is the scale knob and it is swept.** Sweeping diffusion time turns the one tree into a multiscale scanner; the band-specificity lives in that scale-dependence, not at any single τ. The C(τ) panel (with the propagator inset) is what shows this on the slide.
- The comparison step: rank-correlate two phases' cophenetic distances to ask how similar their hierarchies are. The trace schema uses all three phases (rest_pre → task → rest_post): a trace = the task-induced reorganization that persists into rest_post. The estimator is ρ — a Spearman of the cophenetic displacements, symmetrized over a split rest_pre baseline.

3. ON-SLIDE TEXT
Bullet list (on-slide):
- timeseries → |ImCoh| connectivity
- sparsify → connected cycle-rich backbone (dense graph = one scale; backbone = many)
- diffuse & sweep τ → specific heat C(τ) resolves the scales → cophenetic tree
- compare rest_pre → task → rest_post → ρ

Formulas — compile in Canva's LaTeX tool. ℒ = the deterministic (UPGMA / average-linkage) hierarchical clustering; the backbone is the maximum spanning tree ∪ the strongest edges (`mst@0.20`); τ = 1/λmax at the fine scale.

(1) cophenetic distance = merge-height of the lowest common ancestor in the average-linkage tree
D^{\mathrm{coph}}_{ij}=h_{\mathcal{L}[D]}\bigl(\mathrm{LCA}(i,j)\bigr)

(ℒ[D] = the UPGMA tree of the communication distance D — ℒ already shown; h(·) = merge height of a node; LCA(i,j) = lowest common ancestor of leaves i, j. So D^coph_ij = the height at which i and j first join in the dendrogram.)

(2) trace estimator across rest_pre → task → rest_post
\rho=\tfrac12\bigl[\rho_{S}\bigl(D_{\mathrm{task}}-D^{A}_{\mathrm{pre}},\ D_{\mathrm{post}}-D^{B}_{\mathrm{pre}}\bigr)+(A\!\leftrightarrow\!B)\bigr]

(every D is the cophenetic D^coph; **D^A_pre, D^B_pre = the two split-halves of rest_pre** — that is what A and B are; ρ_S = Spearman; (A↔B) = the same term with the two rest_pre halves swapped, and the ½ averages them so no half is privileged. First slot = how the task moved the hierarchy; second slot = how rest_post moved it; ρ asks whether the two agree. In the real analysis ρ is read across the τ-scan, not at one scale.)

4. SPEECH
Here's the whole machine on one slide, on one example patient. From the raw sEEG we build the per-band imaginary-coherence network — and now the step you haven't seen yet. That network is fully connected: every contact couples to every other, and a fully-connected graph has essentially one scale — the diffusion has nothing multiscale to resolve. So we sparsify: keep the maximum spanning tree plus the strongest fifth of the edges — a connected, cycle-rich backbone. On THAT backbone the diffusion propagator, e to the minus tau L, genuinely telescopes: as we sweep the diffusion time the specific heat C of tau grows a mesoscale ladder — several peaks, several scales — instead of one collapse; the inset shows the propagator spreading as tau increases. The last panel reads that diffusion geometry into a cophenetic tree — the hierarchy of the network at scale. How we turn two of those trees, one per phase, into a single trace number — ρ — is the next slide.

Careful: the sparsification is WHY the LRG is non-degenerate here — on the fully-connected graph the trace is identical with and without diffusion (audit_174), so the backbone is not cosmetic; it is what makes the propagator multiscale (audit_175). τ IS swept — it is the scale scanner, and the band-specificity lives in the scale-dependence. Present ρ as the estimator, NOT a verdict — the matched-strength control that turns ρ into a verdict comes on slides 14–15 (kept off this slide on purpose). Matched-strength is now the ONLY null (drift retired — a directional task makes a drift null degenerate with the trace).

5. FIGURES
- THE PIPELINE SEQUENCE — 6 clean single panels (Pat_05 β; panel 4 in 2 ring-order variants, panel 6 with an optional tensor v2), compose left→right with arrows in Canva. The terminal ρ^coph symbol is NOT here (moved to the ρ-measure slide). Rebuilt 2026-07-13 for the sparsification arc, gen `scripts/07_figures/gen_pipeline_sequence.py` (`--patient Pat_05 --band beta`). All PNG (transparent, Canva-friendly) in data/outputs/figures/talk/:
  1. brain electrodes — pipeline_seq_1_brain.png (3D translucent pial brain, all 10 implants as mm-space spheres coloured per patient, ONE left-hemisphere view — the slide-03 cohort-overlay style, single perspective; the establishing "here is our sEEG coverage" shot. `--brain patient` renders the example patient alone if a single-implant lead-in is preferred).
  2. sEEG timecourses — pipeline_seq_2_timeseries.png (a handful of raw contacts, ~2 s, common gain).
  3. |ImCoh| FC matrix — pipeline_seq_3_imcoh_matrix.png (generic FCᵢⱼ colorbar, log colour; the DENSE graph).
  4. mst@0.20 backbone network — circular **hierarchical-edge-bundled chord** (graph-tool/cairo): the 1384 backbone edges Holten-bundled THROUGH the LRG communication hierarchy into cool magenta arcs; beads shaft-coloured; ONLY the ring + bundled arcs (no dendrogram overlay, colorbars, or labels). TWO ring-order versions (both emitted automatically by one run; `panel4_backbone(..., ring_order=)`):
     - `pipeline_seq_4_backbone_network_shaft.png` — beads **sorted by shaft**, one contiguous colour segment per sEEG probe (anatomical reading; busier centre).
     - `pipeline_seq_4_backbone_network_auto.png` — beads in the **drawing algorithm's own order** = the LRG hierarchy leaf order (`leaves_list`), so communities sit contiguously and the bundles tighten into a clean branching HEB (the "cool graph-tool chord" look). Beads still shaft-coloured, so you see how each probe scatters across communities.
     (Same family as the results-section `fig_hierarchy_chord_network` drafts. Choose per slide taste — `_auto` is the cleaner structure shot, `_shaft` the anatomically-labelled one.)
  5. specific heat + entropy — pipeline_seq_5_specific_heat.png (WIDE landscape: C(τ) and Ŝ(τ) vs s = τ·λmax on the LEFT with 4 timescales marked; the propagator K̂(τ)=e^{−τL̂} at those 4 τ as a **2×2 colour-framed block on the RIGHT** — each matrix's frame colour = its vertical mark on C(τ), reading TL,TR,BL,BR = ascending s; NO s= titles (frame colour is the key, the imshow colour carries the value). Shows K spreading/coarsening with τ — the multiscale ladder the sparsification unlocks and the visual for "τ = the scanner").
  6. cophenetic tree — pipeline_seq_6_dendrogram.png (UPGMA tree of the backbone communication distance, **normalised** ultrametric height, with ONE illustrative horizontal cut at ψ_n(τ) colouring the partitions — a teaching device for "the LRG carves communities at scale τ"; `--cut-k` sets the illustrative partition count, default 12. The real cross-phase comparison uses the FULL cophenetic distances, not this cut. NB Pat_05 β is delocalised — genuinely ~two dominant communities — so the cut is hand-placed to read as a partition, not derived from `optimal_threshold`, which is degenerate on this floor-tailed tree).
  6b. (ALT / v2) dendrogram TENSOR — pipeline_seq_6_dendrogram_tensor.png (`--tensor-cards`, default 4). The tree at 4 **increasing** diffusion times τ (≈0.3, 5, 51, 268 for Pat_05 β), spanning the FULL diffusion range from finest (`Ŝ≈0.98`, many communities) to the collapse (`Ŝ≈0.02`, ALL nodes in one cluster = the ground state), each drawn 2D on oblique-projected slices receding up-and-right — a "tensor" of the hierarchy. Front (lower-left, opaque) = finest τ (rich tree, ~9 illustrative partitions); each deeper slice is coarser (fewer partitions), the LAST (largest-τ) slice being the single-cluster collapse the LRG renormalises to (labelled `Ŝ→0`). This is the textbook multiscale **coarsening**: more diffusion = more merging = fewer clusters (9→6→3→1). Deeper slices colour-matched to the panel-5 propagator insets (`MARK_CMAP`), faded with depth. Height axis = the symbol **𝒟** (`\mathcal{D}`, the dendrogram merge height — NOT `D^coph`, to avoid confounding the tree height with the cophenetic-distance matrix in the ρ formula); cards labelled by **order-of-magnitude τ** (τ∼10⁻¹, 10⁰, 10¹, →∞ — the last is the collapse; the numbers are scale markers, not precise values). Deck laid ~45° (landscape/horizontal, wide). NB β is delocalised (natural community count ~2 at every resolved scale) — the stepped cut is an illustrative teaching cut, not a balanced-community claim; and the count *coarsens* with τ (it cannot grow while ending in the single-cluster S→0 collapse). Purely illustrative "τ is the scanner" beat — use INSTEAD of panel 6 for a single striking image, or beside it. Stack order is flippable (collapse-in-front) if the slide wants partitions reading as growth-into-page. Same generator (`gen_pipeline_sequence.py`).
- The trace-estimator symbol (ρ^coph) is NOT on this slide — it lives on the dedicated **ρ-measure slide**. (`panel7_rho()` retained in the generator but not emitted by the pipeline run.)
- The two formulas render as on-slide LaTeX (compiled in Canva), not Lane-F PDFs.
- (Alternative single-image pipeline) the 3D HERO — data/outputs/figures/talk/pipeline_brain_tree_3d.png. NOTE: shared with slide 08 (the bet). If slide 08 uses it as the vision/teaser, slide 11 should LEAD with the step-by-step SEQUENCE above (the build), not the same still.
- The retained-clade tanglegram (ρ^coph visual intuition) is NO LONGER a pipeline panel — it lives as its own explainer figure (`fig_cophenetic_tanglegram.py`, commit a557223) if a "what does ρ mean" beat is wanted separately.

6. REFERENCES
- LRG operator + communication distance + UPGMA dendrogram: Villegas, P., Gili, T., Caldarelli, G., Gabrielli, A. (2023), "Laplacian renormalization group for heterogeneous networks", Nature Physics 19, 445–450; Villegas, P., Gabrielli, A., Poggialini, A., Gili, T. (2025), Physical Review Research 7, 013065 (the fully-connected → single-C(τ)-peak degeneracy and the multiscale-ladder regime).
- Sparse backbone (MST ∪ strongest edges, connected + spanning): Tumminello et al. 2005 (PNAS); Massara, Di Matteo & Aste 2016 (TMFG). audit_175 established `mst@0.20` for this cohort.
- |ImCoh| FC: Nolte et al. 2004 / Ewald et al. 2012 (introduced slide 8).

7. CANVA STATUS
⚠ 2026-07-13 (late): ρ^coph glyph REMOVED from this slide (it belongs on the dedicated ρ-measure slide); pipeline is now 6 panels ending at the tree. Panel 4 chord emitted in TWO ring orders (`_shaft` anatomical + `_auto` algorithm/leaf-order, the clean HEB). Panel 6 tensor v2 reworked: 4 GRADUAL trees at increasing τ (labelled by τ, not s), all inside the resolution window (no star-collapse), partitions GROW front→back (≈4→9) via a stepping cut — "communities emerge with τ, never many→1 in one step" (β delocalisation means one giant clade persists; the count is the illustrative cut).
Page 11 · ~55% (content locked 2026-07-13; 6-panel sequence rebuilt + polished). Have: title; the 7-panel example pipeline (Pat_05 β) with the sparsification + specific-heat panels; TWO compile-ready LaTeX formulas (D^coph via ℒ · the trace estimator ρ) — one line each, matched to code (`heat_multiscale.rho_sym`, `tree.cophenet_matrix`), NO null on this slide (deferred to 14–15). ⚠ 2026-07-13: NARRATION CHANGED — sparsification (`mst@0.20`) is now a first-class pipeline step and τ is the swept scanner. Panels went 5→7. ⚠ 2026-07-13 PANEL POLISH (4 fixes): (1) panel 1 brain = the 3D translucent-pial cohort overlay (all 10 implants, per-patient colour), ONE left view — replaced the flat nilearn glass brain; (2) panel 4 backbone = Fruchterman–Reingold spread + faint grey web (killed the same-probe RED edges and the hairball ball); (3) panel 6 tree = **normalised** ultrametric height + ONE illustrative ψ_n(τ) cut colouring the partitions (hand-placed, `--cut-k 12`; `optimal_threshold` is degenerate here because Pat_05 β is delocalised — genuinely ~2 dominant communities); (4) panel 7 glyph = **ρ^coph** (was ρ). ⚠ FORMULA (2) is ONE line with A/B written as D^A_pre, D^B_pre (the two rest_pre split-halves). NB the on-slide estimator formula is ρ (symmetric); the terminal pipeline glyph is ρ^coph (the cross-phase cophenetic correlation) — intentional per user. RE-RENDER the two LaTeX boxes + RE-COMPOSE the 7-panel row in Canva. Remaining to 100%: presenter composes the 7-panel sequence + 2 formulas in Canva.
