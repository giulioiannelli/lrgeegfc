---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
section: 3
---

# Section 3 — fix list

**Head.** §3.1 LRG recap is mathematically clean and matches `lrg_canonical_formulas.md` (memory). Two issues at §3.2: (i) the
"MSC n*≈23-24 vs |ImCoh| n*≈20-22 conventional cuts" comparison is apples-to-oranges
because §3.2 itself states the |ImCoh| cuts are "conventional cuts chosen for visual
comparability with MSC, not Ψ-optimal cuts in the strict sense" — the contrast as
written suggests both are Ψ-best when only the MSC ones are; (ii) the brain-space
Fig. 12a/12b uses three-patient exemplars (Pat_02, 03, 05) without disclosing the
selection. The §3.2 cohort-count framing also inherits the n=6 vs n=10 contradiction
flagged in §2 (§3.2.2 prose says "across the cohort", which the §4/§5 cohort reads as
n=10, but §3.2 figures still draw from the n=6 subset).

## Confirmed (no action)

- §3.1 Eq. (4): `ρ̂(τ) = e^{−τL} / Tr(e^{−τL})` matches NOTES_MSC Eq. (4) and `lrg_canonical_formulas.md` (codebase implementation in `lrg_eegfc.utils.lrg.propagator`). ✓
- §3.1 Eq. (5): communication distance `D_ij(τ) = (1 − δ_ij)/K_ij(τ)` with `K(τ) = e^{−τL}` is the Villegas 2023 ultrametric construction. ✓
- §3.1 Eq. (6): `Ψ(n; τ) = N · [log10 Δ_n − log10 Δ_{n+1}]` matches NOTES_MSC Eq. (8) and the `compute_optimal_threshold` implementation per `lrg_psi_optimal_threshold.md`. ✓
- §3.1 closing: "Inter-phase comparison uses the Variation of Information evaluated across all dendrogram cut levels k = 2, ..., ⌊N/2⌋, combined with a strict cross-patient unanimity criterion" — matches NOTES_MSC §5.4 and the cohort-coverage-matrix work. ✓
- §3.2.1 dendrogram structure observation: MSC vertically stratified with bursts at distinct heights, |ImCoh| quasi-uniform merge distribution, no preferred scale — matches the §5.1 ψ-irrelevance finding (`_psi_argmax_pinning.md`: 89/180 = 49% strict argmax-at-extremum) and `lrg_outlier_case_fully_connected` memory.
- §3.2.1 closing: "the MSC dendrogram's apparent preferred scale n* ≈ 23-24 is the algorithmic signature of the same-shaft bias and dissolves once |ImCoh| is adopted" — load-bearing claim, supported by the Ψ-irrelevance discussion (89/180 cells pin to a Ψ extremum on |ImCoh|).
- §3.2.2 brain-space same-shaft enrichment ratios: MSC n=5 mean 1.40×, n=20 mean 4.21×; |ImCoh| n=5 mean 1.05×, n=20 mean 1.75× → cohort-mean reduction ≈ 2.4× at n=20. Consistent with §2.2.3.

## Numerical corrections (action: writing agent)

- **§3.2.1**: `n*≈23-24` for MSC and `n*≈20-22` for |ImCoh|. The latter is described as "conventional cuts chosen for visual comparability with MSC, not Ψ-optimal cuts in the strict sense" — this disclaimer is correct but the contrast as written ("a well-defined preferred scale at n*≈23-24" vs "a larger number of smaller communities n*≈20-22") reads as a like-for-like comparison. Recommend rewording: "The cuts shown for |ImCoh| (n* ≈ 20-22) are conventional reading scales chosen for visual comparability with the MSC Ψ-best cuts (n* ≈ 23-24); they are not themselves Ψ-best, because as developed in §5.1 the Ψ stability index has no interior maximum on |ImCoh| substrates." (Already in §3.2.1 prose later — can be tightened by lifting it into the headline contrast sentence.)
- **§3.2 cohort-count framing**: "Same-shaft enrichment at n=5 [..] mean 1.40×; at n=20 [..] mean 4.21×" — confirm these are computed over the same n=6 pre-replacement subset as §2 (Fig. 1, 2). If so, add a sentence: "computed on the n=6 pre-replacement subset (Pat_02, 03, 05, 06, 07, 08); cohort regeneration at n=10 is forthcoming." OR regenerate at n=10 (preferred per §2 fix list).
- **§3.2.2 last sentence**: "Multi-community structure under |ImCoh| becomes informative only at n ≥ 20, consistent with the flat hierarchical landscape of the |ImCoh| dendrogram" — supported by Fig. 12b at n=5 collapse to single dominant cluster + few singletons. ✓

## Framing rewrites (action: writing agent)

- **§3.2.1 lead sentence** (currently "Under MSC the branch structure is vertically stratified ..."): consider opening with the headline: "The MSC and |ImCoh| dendrograms differ structurally rather than by rescaling. The MSC linkage is vertically stratified into distinct merge bursts (Fig. 11, top), producing a Ψ-best partition at n* ≈ 23-24; the |ImCoh| linkage is quasi-uniform with no Ψ-best interior cut (Fig. 11, bottom; full Ψ-irrelevance discussion in §5.1)." This makes the apple-to-apple status clear up front.
- **§3.2.1 second-paragraph open**: "Under |ImCoh| the picture changes qualitatively" → keep, but the very next clause ("Merges are distributed quasi-uniformly along the height axis: branch lengths are short and no preferred scale dominates") would benefit from a concrete number from `_psi_argmax_pinning.md`: "89 of 180 cohort cells (49%) pin Ψ to one of the two dendrogram extrema; 169/180 (94%) within five indices of a boundary".
- **§3.2.2 brain-space exemplar disclosure**: the figure caption mentions "for Pat_02, Pat_03, and Pat_05" but the §3.2.2 text doesn't disclose why these three. Add: "We show Pat_02, Pat_03, and Pat_05 as exemplars of three contrasting implant geometries (left-hemisphere temporal-dominant, right-hemisphere temporal, and bilateral); the same-shaft enrichment numbers are computed over the full §2.3 cohort, not the three-exemplar set." If the enrichment numbers are computed over the n=6 subset, see numerical-correction note above.

## Caveats to add (action: writing agent)

- **§3.2.1 — Ψ-failure reference forward**: at the end of §3.2.1 ("Once the bias is removed by switching to |ImCoh|, the gap closes, and as we develop in §5.1, Ψ no longer selects an interior scale anywhere across the cohort"), add a one-line forward to the actual count: "(89 of 180 (patient, band, phase) cells pin to a Ψ extremum; only 2/180 sit more than ten merge indices interior; full diagnostic in §5.1)." This makes the "no interior cut" claim concrete and saves the reader from §5.1 hunting.
- **§3.2.2 cohort-count mismatch**: if the enrichment numbers are from n=6, state this explicitly (otherwise see Figure actions below).

## Figure actions (coding agent)

- **Fig. 11 (dendrograms, three exemplars Pat_02 β, Pat_05 θ, Pat_08 α)**: keep the three-exemplar choice (it's pedagogically targeted at three different bands × three different patients). Optional: regenerate the |ImCoh| panels with n* set to actual Ψ-argmax (likely a leaf-pair-side extremum) to make the Ψ-irrelevance visually explicit. The current n*=20-22 conventional cuts work fine for the §3.2 visual purpose if the prose disclaimer is lifted into the figure caption.
- **Fig. 12a / 12b (brain-space n=5 vs n=20)**: three-patient exemplar (Pat_02, 03, 05). Disclose exemplar choice in caption: "Three exemplars chosen to span three implant geometries (Pat_02 left-temporal, Pat_03 right-hemisphere, Pat_05 bilateral); the same-shaft enrichment values quoted are cohort-mean averages over the full n=10 cohort, see Tab X" — if and only if the values are cohort-wide. If they are restricted to the three exemplars, state that.
- If §2 figures are regenerated at n=10, **regenerate Fig. 11 dendrogram exemplars at the latest cache state** to ensure the |ImCoh| substrate is the post-replacement Pat_14 cache.

## Deferred / questions

- Should the §3 dendrogram exemplar set (Pat_02 β, Pat_05 θ, Pat_08 α) be expanded to include at least one of the post-replacement patients (Pat_10, 13, 14, 15)? The current trio happens to all be from the n=6 subset; including a Pat_14 panel would help bridge the n=6 / n=10 visual gap. Optional and is a polish question for §3.
- §3.2.1 cites a "20–22 conventional cuts" range without disclosing which patient gets which n*. Pat_02 β = 20, Pat_05 θ = 22, Pat_08 α = 22. Since §3.2.1 is methodological-comparison and not statistical, consistency is fine, but a one-line disclosure ("n* per panel: 20, 22, 22 respectively") would help.
