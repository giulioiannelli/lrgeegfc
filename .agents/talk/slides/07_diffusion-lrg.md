---
name: talk-slide-06-diffusion-lrg
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 7
status: draft
updated: 2026-07-11
canva: page 7 · ~40%
---

# Slide 6 — Revealing multiscale structure with diffusion

1. TITLE
Revealing multiscale structure with diffusion

2. MAIN CONCEPT
- The lens is as fundamental as any process in nature: diffusion (the LRG framework, Villegas–Gabrielli).
- Picture a quantity spreading on the network like heat in space: dx/dt = −L̂ x, with L̂ = D̂ − Â the graph Laplacian.
- Let it flow, and diffusion reveals structure: the diffusion time τ is a scale knob — run it briefly to see local structure, longer to see global structure — so one parameter sweeps the network's intrinsic topological scales. (τ is the diffusion-time parameter, NOT the neural time axis — that's the frequency bands, slide 10. Don't let "time" here read as the time axis.)
- The propagator K̂(τ) = e^{−τL̂} gives a communication distance D_ij(τ) ~ 1/K_ij(τ); sweeping τ zooms out and resolves nested, multiscale organization.
- Higher-order *graph* features come for free: e^{−τL̂} = Σ_n (−τL̂)ⁿ/n! sums walks of every length, so multi-step paths and communities emerge from the *same pairwise graph* — slide 5's "no extra model" road paying off. (Higher-order in the GRAPH, from *linear* diffusion — NOT simplicial interactions, NOT nonlinear signal statistics. See imcoh-2nd-order memory / linearity-and-higher-order-structure guide.)
- Diffusion resolves the TOPOLOGY axis of the multiscale brain — the relational scales, from pairwise edges to nested communities. It completes the trio: the frequency bands give the time axis (slide 10), the sEEG contacts give the space axis (slide 3), diffusion gives the topology axis.

3. ON-SLIDE TEXT
heat on a network: dx/dt = −L̂ x
propagator e^{−τL̂} = Σ (−τL̂)ⁿ/n!  (sums walks of all lengths)
communication distance D_ij(τ) ~ 1/K_ij(τ)
τ = the zoom knob → nested structure
diffusion = the topology axis   (bands = time · sEEG = space)

4. SPEECH
Our tool is diffusion — heat spreading on the network, governed by the Laplacian. Let it flow and something remarkable happens: how long you let it run sets the scale it resolves — briefly, you see local structure; longer, whole communities appear. So from the propagator we read a communication distance between contacts, and by turning that one knob — the diffusion time — we zoom out and watch nested structure emerge. And the higher-order structure comes for free: the propagator sums walks of every length, so multi-step paths and communities are built in — from the same pairwise graph, no simplices to model by hand. That's the "no extra model" road from two slides back, cashed out. Finally, how it fits the bigger picture: diffusion resolves the topological axis — the relational scales. Put it with the frequency bands for the time axis and the sEEG contacts for the space axis, and our approach covers all three axes of the multiscale brain.

5. FIGURES
- Video: diffusion τ-sweep on an FC graph (fc_diffusion_tausweep / lrg_diffusion_zoom).
- Dendrogram cut snapshots — data/outputs/figures/lrg_diffusion_zoom/lrg_dendrogram_cut_snapshots.pdf

6. REFERENCES
- LRG framework — diffusion, Laplacian, propagator, entropy (the tool this slide introduces): Villegas, P., Gili, T., Caldarelli, G., Gabrielli, A. (2023), "Laplacian renormalization group for heterogeneous networks", Nature Physics 19, 445–450.
- Multi-scale communities + the communication distance D_ij(τ) our dendrogram actually uses: Villegas, P., Gabrielli, A., Poggialini, A., Gili, T. (2025), "Multi-scale Laplacian community detection in heterogeneous networks", Phys. Rev. Research 7, 013065.

7. CANVA STATUS
Page 7 · ~40%. Have: good concept text, the reference, presenter notes. Missing: the diffusion video/figure, compressed on-slide text.
