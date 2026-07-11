---
name: talk-slide-04-functional-connectivity
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 4
status: locked
updated: 2026-07-11
canva: page 5 · ~55%
---

# Slide 4 — Functional connectivity

1. TITLE
Functional connectivity

2. MAIN CONCEPT
- The network-scientist's route to the relational features of neural activity is functional connectivity: for each pair of contacts, a coupling estimate on their time series gives a weighted graph — FCᵢⱼ = coupling(xᵢ, xⱼ) → W. (The estimator is kept generic here; the specific choice, imaginary coherence, arrives on slide 9.)
- It is a genuinely functional object — related to, but not reducible to, anatomy: strong coupling appears even between regions with no direct structural link (Honey et al. 2009). We do NOT teach a "structure vs function" dichotomy — only that we work with a functional object, potentially independent of the wiring.
- Network neuroscience is a mature field with fundamental results on this graph: brain networks are small-world (Bassett & Bullmore 2006), built around a rich club of hubs (van den Heuvel & Sporns 2011), and organized into modular communities (Sporns & Betzel 2016) — balancing integration and segregation.
- "Communities" already hints at the multiscale organization we push on next (slide 5). Our contribution is what we do next to that graph.

3. ON-SLIDE TEXT
FCᵢⱼ = coupling(xᵢ, xⱼ)  →  weighted graph W
coupling ≠ wiring
small-world · rich-club · communities

4. SPEECH
How does a network scientist get at relations between brain areas? Functional connectivity. For every pair of contacts we take their time series and compute a coupling estimate — how much they fluctuate together — and that pairwise coupling is a weighted graph. Here it is for four of our patients, across the four phases of the experiment — the same functional object every time. It is genuinely functional: it does not just retrace the anatomy — strong coupling shows up even between regions with no direct structural link. On that graph, decades of network neuroscience have found real organization — small-world, a rich club of hubs, communities that trade off integration against segregation. A mature, productive toolkit. Our contribution is what we do next to that graph.

Careful: the coupling symbol on the formula is GENERIC — do NOT write imCoh or name β here (the specific estimator comes on slide 9). On this slide the matrix is just "a functional network".

5. FIGURES
- FC matrix mosaic — 4 patients × 4 phases (rest_pre → task_learn → task_test → rest_post) — data/outputs/figures/talk/n4_delta_imcoh_abs_generic_log_per_row.pdf. Method-neutral BY CONSTRUCTION: colorbar reads a generic "FCᵢⱼ" (no imCoh / no band named), no probe-letter ticks; dense structure reads straight off the dataset slide and quietly previews the four-phase structure the talk exploits. (Underlying data = δ imcoh_abs — honest in the filename, NOT shown on-figure.) Generator (uses the new colorbar_label override + use_lrg_style fix): `.agents/guides/05_plotting/fc_templates/mosaic_patient_phase.py --patients Pat_02,Pat_05,Pat_08,Pat_13 --band delta --fc-method imcoh_abs --tick-labels generic --colorbar-label '$\mathrm{FC}_{ij}$' --out-dir data/outputs/figures/talk`.
- (Replaced 2026-07-11) the sfdp matrix↔network PDF (Pat_05_beta_rpre_…_g2.pdf) — too sparse, read as "minimal" (user) → dropped for the mosaic above.
- (Backup / not on this slide) 3D cortex connectome — data/outputs/figures/talk/fig_brain_connectome_Pat_05_n10_rsPre_beta.png. Available for later use.
- (Dropped) structure-vs-function schematic — the point ("functional, need not mirror anatomy") is carried inline (on-slide text + Honey 2009); no Canva schematic.
- (Optional Canva tag) a small "small-world · rich-club · communities" caption strip beside the mosaic, so the fundamental results are visible, not just spoken.

6. REFERENCES
Umbrella review (canonical FC / graph-analysis intro):
- Bullmore, E. & Sporns, O. (2009), "Complex brain networks: graph theoretical analysis of structural and functional systems", Nature Reviews Neuroscience 10(3), 186–198.
Fundamental results pointed to on-slide (name the finding, tag 1–2 in Canva):
- Small-world: Bassett, D. S. & Bullmore, E. (2006), "Small-world brain networks", The Neuroscientist 12(6), 512–523. [Lane R: pin locators/DOI; optionally the 2017 "…Revisited", Neuroscientist 23(5), 499–516.]
- Rich-club hubs: van den Heuvel, M. P. & Sporns, O. (2011), "Rich-club organization of the human connectome", Journal of Neuroscience 31(44), 15775–15786. [Lane R: pin DOI.]
- Modular communities: Sporns, O. & Betzel, R. F. (2016), "Modular brain networks", Annual Review of Psychology 67, 613–640. [Lane R: pin DOI.]
- Function ≠ structure (the "functional object, potentially unrelated to wiring" point): Honey, C. J. et al. (2009), "Predicting human resting-state functional connectivity from structural connectivity", PNAS 106(6), 2035–2040. [Lane R: pin DOI.]

7. CANVA STATUS
Page 5 · ~55% (LOCKED 2026-07-11, mosaic + formal-object framing). Figure LOCKED: FC matrix mosaic — 4 patients × 4 phases, generic FCᵢⱼ colorbar (talk/n4_delta_imcoh_abs_generic_log_per_row.pdf). Text final (formal object FCᵢⱼ = coupling(xᵢ,xⱼ) → W · coupling ≠ wiring · the three fundamental results). Speech final. ⚠ Deck page still has raw markdown pasted in — replace with the compressed 3-line block; drop the 3D-brain + structure/function placeholders; optionally add the "small-world · rich-club · communities" caption strip + 1–2 citation tags. Keep the coupling symbol generic (no imCoh). Re-sync exact page in the one-pass Canva re-map (deck shifted ~+5 after the slide-2 split). Remaining to 100%: presenter pastes text + places the mosaic (+ optional tags) in Canva; Lane R pins the four new DOIs.
