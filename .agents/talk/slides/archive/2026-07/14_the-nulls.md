---
name: talk-slide-14-selectivity
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 14
status: draft
updated: 2026-07-13
canva: page 12 · REFOCUSED 2026-07-13 (SELECTIVITY / controls-ladder; matched-strength sole null) — retitled off "who's right?", Grassmann demoted, "it must be multiscale" removed; RE-RENDER on deck
---

# Slide 14 — The null names the bands: selectivity, not exclusivity

1. TITLE
The null names the bands — selectivity, not exclusivity

2. MAIN CONCEPT
- **Head (the juice):** the one null that bites — **matched-strength** — doesn't crown a single "right" read; it exposes what each read is *for*. Send a whole ladder of read-outs through that one gate on the *same* sparsified graph and the payoff is **selectivity**: only the multiscale cophenetic hierarchy *names the two carrier bands* (α, β). Raw edges detect the reorganization in four bands — β among them, p = 0.024 — but can't say which band carries the trace; a genuinely spectral read is dead in every band. **The edge is which-bands, not who-sees-it.**
- The referee has to be a null, not another measure (slide 13 set that up). There is exactly **one** meaningful null: **matched-strength** — rewire while preserving each contact's total connection strength, so anything that was merely "where connectivity happens to be strong" disappears (the classic mistake in this literature). Drift is *retired* (2026-07-12): a **directional** task makes a trace monotonic by construction, so a drift null collapses onto the very alternative it was meant to exclude. Matched-strength is the single bar every claim in the talk clears.
- **The selectivity ladder** — same graph, same null, seven read-outs:
  - **raw FC (pairwise edges)** → fires in **four** bands (δ .024 · α .019 · β .024 · low-γ .032). It *detects* the reorganization — β included — but registers movement everywhere it can, so it cannot separate a carrier band from a non-carrier. **Non-selective, not blind.**
  - **node strength** → degenerate: the strength-matched surrogate reproduces it exactly, so it can never clear its own null (p = 1).
  - **clustering coefficient** → catches **α only**, misses β.
  - **graph geodesic** → fires the **wrong** bands (δ, low-γ), misses both carriers.
  - **effective resistance (spectral)** → the sharp foil: a global, multi-step Laplacian read-out that *should* be the strongest rival — and it is **dead in every band** (best cell β p = 0.116, nothing clears). Global and spectral is not enough.
  - **cophenetic (multiscale, τ-swept)** → the **only** band-selective read: **α and β and nothing else**, holding from the fine scale to the mesoscale (β 0.014 → 0.001; α 0.024 → 0.007).
- **The take-away is information-DEPTH, framed as selectivity.** Ordinary pairwise FC gives a scalar verdict per band — "did this band lean back toward the task?" — and several bands say yes. That scalar registers *that* something moved; it cannot say *which* structure carries it. The multiscale read is a structured object, so it can: it names the two carriers cleanly where the pairwise, node, topological, and spectral read-outs all fail to. **Not "only multiscale sees the trace" — the multiscale read is the one that selects the bands.**
- Grassmann is a **complementary whole-graph companion** (it also detects β, whole-graph — not on the sparsified backbone), *not* a co-equal read on this ladder; the spectral foil here is effective resistance. (*Which scale* each carrier lives at — β flat, α peaked — is the next slide.)

3. ON-SLIDE TEXT
same sparsified graph · one null (matched-strength) · which read *names* the bands?

raw edges (pairwise)              → 4 bands: δ α β γ_low   — detects, doesn't select (β p = .024)
node strength                     → degenerate: null preserves it (p = 1)
clustering                        → α only
graph geodesic                    → wrong bands: δ, γ_low (misses α, β)
effective resistance (spectral)   → dead in every band
cophenetic (multiscale, τ-swept)  → α, β only  ✓ selective  (β .014→.001 · α .024→.007)

⇒ the edge is SELECTIVITY, not who-sees-it — only the multiscale read names the carriers
(Grassmann = complementary whole-graph companion, β only — not on this ladder)

4. SPEECH
On the last slide the three reads disagreed, and I said the tie-breaker can't be yet another measure — it has to be a null. Here it is, and there's only one that actually bites: matched-strength. It reshuffles the network but keeps every contact's total connection strength, so anything that was merely "where connectivity is strong" vanishes. Now watch what happens when we push a whole ladder of read-outs through that one gate, on the very same graph. The raw edges fire in four bands — delta, alpha, beta, low-gamma — so they clearly detect the reorganization; beta itself clears at p = 0.024. But four bands is not an answer to "which band carries the trace." Node strength is degenerate — the null preserves it exactly, so it can't clear anything. Clustering catches only alpha. The graph geodesic fires the wrong bands, delta and low-gamma, and misses both carriers. And the sharp one: effective resistance — a genuinely spectral, global, multi-step read that ought to be the strongest rival — is dead in every band. Only the multiscale cophenetic hierarchy names two bands and stops there: alpha and beta. That's the whole point. The edge isn't that the multiscale read is the only thing that sees a trace — the raw edges see beta. It's that it's the only read that *selects* which bands carry it. Selectivity, not exclusivity. And which scale each of those two carriers lives at — that's next.

Careful: matched-strength is the SOLE null (drift RETIRED 2026-07-12 — a directional task makes a drift null degenerate with the trace). Do NOT mention a "second null", a "drift" column, or "survives both nulls". Do NOT say "only multiscale sees the trace" / "invisible to simple methods" / "multiscale-exclusive" — raw DETECTS β at p = .024; frame raw as NON-SELECTIVE (fires 4 bands), never as blind. Do NOT retitle this "it must be multiscale" (that is the retired exclusive-detection framing). The cophenetic read lights α AND β — do NOT say "coph → β-only" here. Grassmann is a whole-graph companion (β only), NOT a co-equal read and NOT next to resistance ("two spectral things") — the spectral foil on this ladder is effective resistance. No localization here (no β→OFC — that is the encoding anchor, slide 18 only), no inference, and NEVER rank bands by the ×-null ratio. Do NOT re-open "who's right?" as a fresh hook — slide 13 owns that tension; this slide RESOLVES it.

5. FIGURES
- **MAIN — the selectivity / controls-ladder matrix.** Descriptor × band, one matched-strength gate; a cell is ringed when it clears α = 0.05. Use the **whole-task-trace (`T_test`) panel**: raw FC = 4 rings (δ/α/β/low-γ), node strength empty (p = 1), clustering = α only, geodesic = δ/low-γ, resistance = empty, cophenetic (τ_min + meso) = α/β only. This is the picture of "only the multiscale read is band-selective." Source: `data/sparsified_arc/figures/controls_ladder/controls_ladder_gate.pdf` (2×2 over functionals — crop the top-left `T_test` panel for the talk); numbers `data/sparsified_arc/controls_ladder/cohort_gate.csv`; gen `scripts/01_compute/sparsified_arc/fig_controls_ladder.py`. ⚠ REBUILD a single-panel talk cut (T_test only) into `data/outputs/figures/talk/` (PNG fine) — keep it **DISTINCT from slide 15** (which becomes the ρ_sym(τ) scale-shape curves).
- Matched-strength null cartoon — "what it destroys" (build in Canva). **Single null — NO drift/second panel** (retired 2026-07-12).
- (optional supporting) cophenetic null triangle — `data/outputs/figures/talk/fig_bands_null_triangle_coph.pdf` (whole-graph; use only if space allows).
- **RETIRE for this slide:** `fig_whos_right_three_methods.pdf` and `fig_before_nulls_three_measure.pdf` — both are the 3-co-equal-reads (raw/Grassmann/cophenetic) synthesis, which puts Grassmann on par with cophenetic and reads as exclusive-detection. Grassmann is now a demoted whole-graph companion, so these panels contradict the selectivity framing. Do not place them here.

6. REFERENCES
- Matched-strength null and the controls-ladder audit (`14_controls_ladder_mst020`) are ours. Effective resistance = Laplacian pseudo-inverse distance, a standard graph quantity (resistance distance, Klein & Randić 1993) — used here as the spectral foil, no external result imported.

7. CANVA STATUS
⚠ 2026-07-13 REFOCUSED (settled three results): slide is now the SELECTIVITY pillar — the controls-ladder run through the sole matched-strength null (drift retired), landing on "only the multiscale read *names* the carrier bands (α/β)"; retitled off "who's right?" (slide 13 owns that hook), Grassmann demoted to a whole-graph companion, and the retired "it must be multiscale" / exclusive-detection framing removed. RE-RENDER on deck.

Page 12 · repo content FINAL for the refocus. Numbers locked against `controls_ladder/cohort_gate.csv` (raw δ.024/α.019/β.024/low-γ.032 = 4 bands; strength p=1; clustering α .019; geodesic δ.042/low-γ.014; resistance dead, best β .116; coph shown at the SETTLED trace verdict α .024→.007 / β .014→.001 — ladder cohort_gate.csv has α .019→.0068, a near-duplicate harmonized to the locked source of truth). Remaining = presenter's Canva job: place the single-panel selectivity matrix (T_test cut), draw the matched-strength "what it destroys" cartoon (single null — NO drift panel), retitle the page, add presenter notes. Coordination flag: this slide now carries the controls-ladder that currently also sits on slide 15 — **slide 15 must be rebuilt to the ρ_sym(τ) scale-shape job** (its current "Who survives — it must be multiscale" title violates the guardrails and duplicates this slide); flag left for the slide-15 rewrite agent.
