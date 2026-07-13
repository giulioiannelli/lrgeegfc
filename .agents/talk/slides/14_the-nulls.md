---
name: talk-slide-14-the-nulls
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 14
status: draft
updated: 2026-07-13
canva: page 12 · ~60% (merge of old 14 who's-right + 15 the-nulls; ⚠ 2026-07-13 drift null RETIRED — matched-strength is the sole referee; RE-DO the null cartoon on deck)
---

# Slide 14 — Who's right? The nulls decide

1. TITLE
Who's right? — the nulls decide

2. MAIN CONCEPT
- Read the SAME connectivity three ways and you get three different answers about the trace:
  - raw edges (pairwise): the low-contrast view — most bands crowd the threshold; only β (and low-γ) barely cross, with no clean selective peak.
  - Grassmann (spectral subspace): δ, β, γ_low light up.
  - cophenetic (multiscale hierarchy): α and β — and it disagrees with Grassmann.
  - β is the one band all three pick out — but they still disagree about every other band, so they can't all be right about the rest.
- They contradict, so the referee can't be another measure — it has to be a null. The meaningful one is **matched-strength**: rewire while preserving each contact's total connection strength → kills any "trace" that is just where connectivity happens to be strong (the classic mistake in this literature).
- Why only one null (retired the second on 2026-07-12): a "drift" null — does time alone within a baseline reproduce it? — sounds fair but is *degenerate* here. The task is **directional**, so a genuine trace is monotonic by construction, and the drift test ends up asking the same question as the alternative it was meant to exclude. Matched-strength is the single bar every claim in the talk must clear.
- The next slide runs matched-strength on all three reads — and then asks the harder question: is the surviving trace really *multiscale*, or could a simpler measure have found it?

3. ON-SLIDE TEXT
same data, three reads:
raw edges → low-contrast: only β, γ_low barely cross
Grassmann (spectral) → δ, β, γ_low
cophenetic (multiscale) → α, β
β = the one band all three pick out — so who's right?
the referee is a null, not another measure:
matched-strength → removes strength-driven structure (the only meaningful null)
(drift retired — a directional task makes it degenerate with the trace)
a real trace must clear matched-strength → is it also multiscale? (next)

4. SPEECH
Here's the problem. Look at the same data three ways and you get three stories. The raw edges are the low-contrast view — most bands just crowd the threshold, only beta and low-gamma barely tip across, so on their own they make the weakest claim about band structure. The spectral, Grassmann view lights up a few bands — delta, beta, low-gamma. Our cophenetic hierarchy says alpha and beta — and it doesn't even agree with Grassmann. The one thing all three share is beta. They contradict each other, so they can't all be right, and the tie-breaker can't be yet another measure — it has to be a null model. The one that bites is matched-strength: it reshuffles the network but keeps every contact's total strength, so anything that was just "the strongly connected regions" disappears. You might expect a second, drift null — does time alone reproduce it? — but here that's a trap: the task is directional, so a genuine trace is monotonic by construction, and the drift test collapses onto the very alternative it was meant to rule out. So we hold every claim to the one bar that actually bites, matched-strength — and on the next slide we run it, and then ask the sharper question: is the trace that survives really multiscale, or could a plain pairwise or spectral measure have found it?

Careful: the cophenetic read-out lights α AND β at the matched-strength gate — do NOT say "cophenetic → β-only" here. β-selectivity is a downstream result (β is the band that also localizes and stays lit across every scale), not the raw coph verdict. β = the intersection of all three measures, which is the honest "who's right?" hook. And do NOT say raw "sees no structure / all bands move together" — in the figure raw lights β AND low-γ (two filled discs, MS p = 0.024 / 0.032). Frame raw as the LOW-CONTRAST / least-selective view (most bands hug threshold, no clean peak), never as blind. DRIFT IS RETIRED (2026-07-12) — do not mention a "second null" or a "double bar"; matched-strength is the sole referee.

5. FIGURES
- Three-method per-band comparison — raw / Grassmann / cophenetic side by side ("who's right") — data/outputs/figures/talk/fig_whos_right_three_methods.pdf (3×6 method×band verdict matrix; disc size = −log10 p_MS, filled = clears the null; β the only band lit by all three). gen scripts/01_compute/figures_embedded/fig_whos_right_three_methods.py
- Null schematic — matched-strength, what it destroys (build in Canva). (drift retired 2026-07-12 — do NOT draw a second/drift null panel.)
- Cophenetic null triangle — data/preprint/figures/all_bands/fig_bands_null_triangle_coph.pdf

6. REFERENCES
- None (matched-strength + drift lineage is ours).

7. CANVA STATUS
Page 12 · ~60% (repo content FINAL). MERGED 2026-07-11 from old slide 14 (who's right) + old slide 15 (the nulls) — the "nulls" half of the 14-15-16 collapse; the survivors table is now slide 15. Repo-side DONE: three-method panel generated + verified against the locked MS summaries (raw β/γ_low, Grassmann δ/β/γ_low, coph α/β; β lit by all three), on-slide text + speech finalized, raw framing reconciled with the figure (low-contrast, not "no structure"). Remaining = presenter's Canva job: place the three-method panel + null-triangle, draw the matched-strength "what it destroys" cartoon in Canva (single null — NO drift panel, retired 2026-07-12), add presenter notes. ⚠ 2026-07-13: on-slide text + speech updated to the single-null (matched-strength) framing; RE-RENDER bullets.
