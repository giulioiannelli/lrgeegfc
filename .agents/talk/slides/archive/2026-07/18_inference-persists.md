---
name: talk-slide-18-inference-persists
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 18
status: draft
updated: 2026-07-13
canva: pages 17–18 · ~30% (encoding + inference persistence, and their different cortex)
---

# Slide 18 — Encoding and inference both persist — in different cortex

1. TITLE
Both persist — different cortex: encoding → OFC, inference → cingulate

2. MAIN CONCEPT
- HEAD (the juice): what post-task rest keeps is not just the pairs the patient was *shown* (encoding) but the order the patient *inferred* (inference-specific) — and the two components sit in different cortex, premises toward orbitofrontal, inferred order toward cingulate.
- ENCODING persists. In β it clears at every scale — 16/16, p_meso ≈ .003 (scale-invariant). In α it persists too (3/16 scales; τ_min .042, meso .032), and there ordinary tools already register it — raw FC .032, clustering .010. So encoding shows the same split the whole trace does: β genuinely multiscale, α resolvable by simple measures.
- INFERENCE-specific persists — and the correction to the old story: it is **NOT β-only**. β 7/16 scales (p = .005) AND α 7/16 (p = .024). The order the brain reasoned out — encoding controlled away — survives into rest, an offline abstraction of a learned relational map (interpretation, not a behavioural claim). [The earlier δ inference-specific effect is a partial-correlation artifact — δ is null for both encoding and the whole trace — and is CUT.]
- The multiscale value-add on the inference cell is **localization, not detection**. Detection is tied by construction: raw FC already detects β inference (p = .010), because the cophenetic distance is built from the same edges. What the hierarchy adds is *where* — it pools many weak edges into a coherent anatomical read that the raw edges cannot. Encoding → OFC (q = .010 / .040 with/without SOZ); inference → cingulate, away from OFC (q = .030 / .040). ⚠️ These addresses are whole-graph (pre-`mst@0.20`) trends, re-derivation on the backbone pending — present as trend-level, not closed.
- Not a recording-length effect: the β inference signal does not scale with test/learn duration (ρ = +0.10, p = .78).
- Companion band result: a focal **low-γ ENCODING** trace in the cingulate that whole-brain averaging hides (ρ_sym +0.25, q = .035, 8/8 sampled ⚠️ whole-graph) — a memorizing copy of the learning phase, distinct from the β/α inference persistence.

3. ON-SLIDE TEXT
encoding and inference both persist — in different cortex
encoding — β 16/16 scales (p_meso .003) · α 3/16 (raw-visible: raw .032, clustering .010)
inference — NOT β-only: β 7/16 (.005) · α 7/16 (.024)
detection tied (raw sees β inference .010) → the value-add is localization
encoding → OFC (q .010/.040) · inference → cingulate (q .030/.040)   ⚠ whole-graph trend
not a recording-length effect (ρ = +0.10, p = .78) · null = matched-strength, read per-scale

4. SPEECH
The trace splits into two components. Encoding is the reorganization while the brain sees the adjacent pairs; inference-specific is the extra reorganization for the order it had to work out, with encoding controlled away. Both survive into rest. Encoding persists strongly — in β it clears at every single scale, sixteen of sixteen, around p .003; in α it persists too, and there ordinary tools already see it. The inferred order also persists, and here is the correction to the older version of this talk: it is not beta alone. Beta at seven of sixteen scales, p .005, and alpha at seven of sixteen, p .024. So rest holds not just the pairs the brain saw, but the order it inferred.
Now, what does the multiscale read buy us on inference? Not detection — raw connectivity already sees the beta inference signal at p .010, because the hierarchy is built from the same edges. What it buys is location. It places the two components in different cortex: the seen pairs anchor toward orbitofrontal, the inferred order sits away from it, toward the cingulate. I'll present those addresses honestly as trends — they come from the whole-graph pipeline and re-derivation on the sparsified backbone is still pending. And it isn't just a recording-length artifact: the inference signal doesn't scale with how long the phases ran, correlation about plus point-one, p .78.

Careful: inference is α AND β, never "β-only" — and the δ inference-specific effect is a partial-correlation artifact, CUT (never mention it). Don't say "only multiscale sees inference" / "invisible to simple methods" — raw FC detects β inference at .010; the value-add is LOCALIZATION, not detection (say "localization, not detection"; never "coph detects better"). The OFC/cingulate addresses are WHOLE-GRAPH trends, re-derivation pending — lead with them as trend-level, not a closed result. OFC here = the ENCODING anchor only; never say "the β *trace* → OFC" (the β trace is delocalized). Null = matched-strength only — no drift, no "second null", no "survives both nulls"; the inference signal stands on matched-strength + duration-immunity, not on drift. Read per-scale — never quote a single +0.091 / p=0.0098 scalar. Keep the low-γ-encoding→cingulate result distinct from the β/α inference. No behavioural data.

5. FIGURES
- Encoding-vs-inference dissociation (main figure: encoding → OFC vs inference → cingulate) — data/preprint/figures/results_section1/fig_encoding_inference_dissociation.pdf  ⚠️ whole-graph localization; present as trend. VERIFY it does not encode "inference = β-only"; REBUILD if it does.
- Inference → cingulate — data/outputs/figures/talk/fig_arc_f_inference_cingulate.pdf  ⚠️ whole-graph, trend-level.
- Low-γ encoding → cingulate (companion band) — data/outputs/figures/talk/fig_arc_e_lowgamma_cingulate.pdf  ⚠️ whole-graph, trend-level.
- REBUILD or RETIRE: data/outputs/figures/talk/fig_reasoning_bloom_inference_raw.pdf — it shows the β wedge glowing ALONE, which now contradicts the settled result (inference persists in α AND β, 7/16 each). Either rebuild so both α and β fire, or drop it.
- (optional) OFC + cingulate on one brain, per-scale inference counts as a small β/α bar (build in Canva).

6. REFERENCES
- Behrens 2018; Wilson 2014; Schuck 2016 (cognitive map); audit_171 (localization, ours).

7. CANVA STATUS
⚠ 2026-07-13 REFOCUSED (settled three results): retitled off "(β only)"; inference is now α AND β (β 7/16 .005, α 7/16 .024), δ inference-specific CUT; per-scale counts replace the +0.091/p=0.0098 scalar; encoding persistence added (β 16/16 .003, α 3/16 raw-visible); OFC/cingulate marked whole-graph trends; drift analogy removed. RE-RENDER on deck.
Pages 17–18 · ~30%. Two beats on one slide (both components persist, then their addresses) — prune on-slide text hard in Canva. Missing: dissociation + cingulate figures placed with the whole-graph trend caveat, a rebuilt (or retired) inference figure that shows α+β not β-only, compressed on-slide text, presenter notes.
