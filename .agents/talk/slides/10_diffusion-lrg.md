---
name: talk-slide-10-diffusion-lrg
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 10
status: draft
updated: 2026-07-16
canva: page 7 · ~45%. FOLDED IN 14b (2026-07-14 regroup): the "what the scale parameter means"
  concept + its micro→meso→macro coarse-graining figure now live here — a standalone scale-concept
  slide was under-used once localization was dropped, but the concept belongs exactly where diffusion
  is introduced. The "functional, not spatial" reach curve is a speaker aside + Q&A backup (it leans
  on |ImCoh|, introduced slide 08), NOT a main panel.
---

# Slide 10 — Revealing multiscale structure with diffusion

1. TITLE
Revealing multiscale structure with diffusion — one parameter that coarse-grains the network
(Canva title is "Laplacian as a useful *cagata pazzesca*" — KEEP AS-IS: a deliberate in-joke callback to another conference speaker's talk, per user 2026-07-16. Do NOT "fix" / sanitize it.)

2. MAIN CONCEPT
- The lens is as fundamental as any process in nature: diffusion (the LRG framework, Villegas–Gabrielli).
- Picture a quantity spreading on the network like heat in space: dx/dt = −L̂ x, with L̂ = D̂ − Â the graph Laplacian.
- Let it flow, and diffusion reveals structure: the diffusion time τ is a scale knob — run it briefly to see local structure, longer to see global structure — so one parameter sweeps the network's intrinsic topological scales. (τ is the diffusion-time parameter, NOT the neural time axis — that's the frequency bands, slide 9. Don't let "time" here read as the time axis.)
- The propagator K̂(τ) = e^{−τL̂} gives a communication distance D_ij(τ) ~ 1/K_ij(τ); sweeping τ zooms out and resolves nested, multiscale organization.
- **What the knob does, made concrete (the main figure).** One implant, read at three diffusion scales: **MICRO** (many single-population clusters) → **MESO** (~10-node groups) → **MACRO** (~2 blocks). As τ grows the local populations merge into progressively coarser communities — that *is* the coarse-graining, and it is the whole meaning of "multiscale."
- Higher-order *graph* features come for free: e^{−τL̂} = Σ_n (−τL̂)ⁿ/n! sums walks of every length, so multi-step paths and communities emerge from the *same pairwise graph* — slide 5's "no extra model" road paying off. (Higher-order in the GRAPH, from *linear* diffusion — NOT simplicial interactions, NOT nonlinear signal statistics. See imcoh-2nd-order memory / linearity-and-higher-order-structure guide.)
- Diffusion resolves the TOPOLOGY axis of the multiscale brain — the relational scales, from pairwise edges to nested communities. It completes the trio: the frequency bands give the time axis (slide 9), the sEEG contacts give the space axis (slide 3), diffusion gives the topology axis. (That relational scale *traces* a physical extent too — local clusters out to the whole implant as you zoom — but it is a **functional / communication** scale, delocalised even at its finest, NOT a spatial zoom. So the topology axis maps only loosely onto space; that honest nuance is quantified later, once |ImCoh| is on the table.)

3. ON-SLIDE TEXT
heat on a network: dx/dt = −L̂ x
propagator e^{−τL̂} = Σ (−τL̂)ⁿ/n!  (sums walks of all lengths)
communication distance D_ij(τ) ~ 1/K_ij(τ)
τ = the scale parameter → coarse-grains the network:  micro → meso → macro
  single populations → ~10-node groups → ~2 blocks   (as τ grows)
diffusion = the topology axis   (bands = time · sEEG = space)  —  a relational, functional scale, not a spatial zoom

4. SPEECH  (~100 s — the coarse-graining mechanism, then the effective spatial scale read off community size)
Our lens is diffusion — heat spreading along the network's connections, governed by the graph Laplacian. [gesture: equations, left] The propagator, the heat kernel, tells you how much reaches one contact from another after a diffusion time tau.

Let it flow [gesture: FC matrix / play], and how long you let it run sets how coarse the communities you read. On one implant [gesture: micro-meso-macro brain panel]: briefly, dozens of little single-population clusters; longer, they merge into a handful of mid-sized groups, then two big blocks — micro to meso to macro, which is what "multiscale" actually means (illustrative on the dense graph, not a data claim). And tau is the diffusion time, a scale parameter, not the neural time axis — that's the frequency bands, from a moment ago.

And the higher-order structure comes for free. [gesture: propagator series] The propagator sums walks of every length, so multi-step paths and communities emerge from the same pairwise graph, by linear diffusion — a higher-order graph feature, never nonlinear, no simplices built by hand. That's the "no extra model" road from earlier, cashed out.

From the propagator we read a communication distance between contacts, and linkage turns it into a hierarchy. [gesture: dendrogram] The tree is the coarse-graining written out — cut it low for many small communities, high for a few large ones.

Now, does that scale mean anything physical? [gesture: right-hand curve] Here's the honest answer, and it's read straight off community size: as the scale grows, the co-diffusing communities grow, and we can measure how much of the brain each one physically spans — that's this curve, the average extent of the community a contact belongs to, against the diffusion scale. To the left of this line you're still reading raw, single edges; cross it, where real diffusion communication begins, and the groups are already spread over the implant — a few centimetres across for most bands — not compact spots. Push it further and they merge out to the whole implant, about nine centimetres. So the scale does earn a physical footprint — local groups out to global — but through community size, as a communication ruler, never a millimetre resolution zoom. These are functional groups, not proximity clusters — precise from how we just measured connectivity, a couple of slides back.

So diffusion resolves the topology axis — from pairwise edges up to nested communities. [gesture: brain+tree 3D scene] With the bands giving time and the contacts giving space, our lens now covers all three axes of the multiscale brain.

Careful:
- The micro/meso/macro coarse-graining panel is ILLUSTRATIVE mechanism on the DENSE graph (only the dense graph gives a balanced 3-regime split; the sparse backbone is core-dominated) — tag it "illustrative", not a data claim.
- The effective-spatial-scale curve is `scripts/07_figures/talk_fig_spatial_reach.py` → `data/outputs/figures/talk/fig_spatial_reach.png` (NOT the preprint `new_results_sec1` asset — that is a DIFFERENT, max-diameter figure). Its ℓ(s) = per-contact MEAN physical DIAMETER of the co-diffusing community (connected components of K ≥ 1/N). It FLOORS at the electrode pitch (3.5 mm, singletons := pitch) at fine s and saturates at the implant span (~91.5 mm) by s ≈ 2–3. So on THIS curve the finest point IS a single contact (raw regime) — do NOT say "the finest community is ~15 mm, never a single electrode" over it (that 15 mm belongs to the OTHER max-diameter figure). Anchor delocalisation at the s = 1 dashed line (raw↔communication boundary): there the reach is already ~25–40 mm for most bands (cohort-median δ39.6/α40.3/θ34.5/γl28.6/β25.2/γh9.9 mm) — delocalised the moment it is communication. Bands differ in the RATE/onset of merging, not the endpoints (all floor at 3.5 mm, all saturate at ~91 mm).
- The vertical dashed line is s = 1 = τ_min = 1/λmax (raw↔communication boundary), NOT s_report (5.6).
- "Functional, not spatial" leans on |ImCoh| geometry (median backbone edge ~35 mm, functional not proximity), introduced slide 08 — forward-reference it here ("precise once connectivity is on the table"), do not fully justify it on this slide.
- τ is the diffusion-time / scale parameter, NOT the neural time axis (that's the bands, slide 9).
- Higher-order in the GRAPH from LINEAR diffusion — never "nonlinear", never simplicial.
- No localization / anatomy here (that beat is retired).

5. FIGURES
- **MAIN — the coarse-graining panel (folded in from 14b, BUILT).** `data/outputs/figures/talk/fig_scale_coarse_grains_space.png` (gen `scripts/07_figures/talk_fig_scale_coarse_grains_space.py`, transparent, grey glass shell, Canva-ready). One fixed implant (Pat_05, β, rest_pre) read at three diffusion scales: MICRO (~46 single-population clusters, avg ~2–3 nodes) · MESO (12 balanced ~10-node groups) · MACRO (2 blocks). Coloured contacts + strongest intra-community edges bowed inward (contacts on top); size-rank colour so the biggest community keeps its hue across scales. Illustrative mechanism (dense graph — the only graph that yields a balanced 3-regime coarse-graining), NOT a data claim. No dendrogram, no text on the asset.
- **MAIN (2nd panel, PROMOTED to on-slide per user 2026-07-16) — the effective spatial scale from community size.** `data/outputs/figures/talk/fig_spatial_reach.png` (gen `scripts/07_figures/talk_fig_spatial_reach.py`, cache `_reach_cohort_cache.npz`). Per-band cohort-median ℓ(s) = mean physical DIAMETER of the co-diffusing community (connected components of the heat kernel K = e^{−τL} ≥ 1/N on the mst@0.20 backbone) vs diffusion scale s = τλmax, with three glass-brain coarse-graining insets and a vertical dashed line at s = 1. Floors at the pitch (3.5 mm) at fine s, saturates at the implant span (~91.5 mm) by s ≈ 2–3; bands differ in the RATE of merging (s=1 reach δ39.6/α40.3/θ34.5/γl28.6/β25.2/γh9.9 mm). This is the visual home of the "s → an effective physical scale, via community size" beat. ⚠ NOT `data/preprint/figures/new_results_sec1/fig_spatial_reach.png` — that is a DIFFERENT (max-diameter) figure that floors at ~15 mm; never swap their numbers.
  - ⚠ HONESTY CALL FOR THE USER: this on-slide talk figure floors at the pitch (single contact) by construction (singleton := pitch; its docstring literally says "the finest resolved scale is the single contact"), so its fine end visually reads as "single-contact resolution" — the exact thing the functional-not-mm lock guards against. Two clean options: (a) KEEP it, and narrate delocalisation from the s = 1 boundary (current speech does this — below s=1 = raw single edges, above = communication, already ~25–40 mm); or (b) SWAP to the max-diameter version (floors ~15 mm, never a single electrode) for a figure that itself never touches the single contact — at the cost of the |ImCoh| dependency being more explicit before slide 08. Descriptive geometry only — no null, no p-value.
  - ⚠ Deck-order tension (flagged, user's call): the "functional, not spatial" justification leans on |ImCoh| (slide 08) which the audience has not met yet; the speech forward-references it rather than fully justifying it here.
- SECONDARY / optional — the diffusion τ-sweep **video** (`fc_diffusion_tausweep` / `lrg_diffusion_zoom`) + dendrogram cut snapshots (`data/outputs/figures/lrg_diffusion_zoom/lrg_dendrogram_cut_snapshots.pdf`): the motion version of the same coarse-graining. Use if the slide wants animation; otherwise the static panel carries it.

6. REFERENCES
- LRG framework — diffusion, Laplacian, propagator, entropy (the tool this slide introduces): Villegas, P., Gili, T., Caldarelli, G., Gabrielli, A. (2023), "Laplacian renormalization group for heterogeneous networks", Nature Physics 19, 445–450.
- Multi-scale communities + the communication distance D_ij(τ) our dendrogram actually uses: Villegas, P., Gabrielli, A., Poggialini, A., Gili, T. (2025), "Multi-scale Laplacian community detection in heterogeneous networks", Phys. Rev. Research 7, 013065.
- Spatial-reach geometry (backup figure): `scripts/01_compute/sparsified_arc/26_diffusion_spatial_reach_mst020.py` (descriptive, no null).

7. CANVA STATUS
Page 7 · ~45%. FOLDED IN 14b (2026-07-14 regroup): the standalone "what the scale means" slide is retired here (→ `archive/2026-07/`) — a full slide for the scale concept was under-used once localization was dropped, but the concept + its micro→meso→macro figure belong exactly where diffusion is introduced. Have: concept text, references, presenter notes, the BUILT coarse-graining panel. Missing on deck: place `fig_scale_coarse_grains_space.png` as the main visual, compress the on-slide text, (optional) the τ-sweep motion version. The "functional, not spatial" reach curve is a Q&A backup (imcoh-dependent), not a main panel. NB the same coarsening reappears as the dendrogram *tensor* on slide 11 (panel 6b) — 07 shows it on the brain (concept), 11 on the tree (mechanics); keep both (reinforcing, two modalities) or thin if it reads repetitive on the deck.
