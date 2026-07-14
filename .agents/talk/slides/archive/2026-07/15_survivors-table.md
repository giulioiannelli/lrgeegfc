---
name: talk-slide-15-scale-shape
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 15
status: draft
updated: 2026-07-13
canva: NEW — no Canva page yet
---

# Slide 15 — Grade the trace by scale: only β is genuinely multiscale

1. TITLE
Grade the trace by scale — only β is present at every scale

2. MAIN CONCEPT
- **Head (the juice):** selectivity told us *which* bands fire (α and β); the τ-sweep tells us *how* they fire — and only β fires at **every** scale. Because the multiscale read is a **structured object**, ρ_sym is a *curve* over diffusion scale, not a single number. That curve **grades each firing band by its scale-shape**: β is **scale-invariant** (present at every scale) = genuinely multiscale, the deepest trace; α is a **single mesoscale peak** = single-scale, real but resolvable; δ fires only at scattered coarse scales = shallow, not claimed. A scalar verdict has no scale axis and cannot grade any of this — this is the payoff of the sweep.
- The measurement: for each band, cross-phase ρ_sym(τ) is swept over the dimensionless scale s = τ·λ_max, and **each scale is held to its own matched-strength null** (the sole null). We read the gate **per-scale**, never best-scale.
- **β — scale-invariant (genuinely multiscale).** Clears the matched-strength gate at **16/16 scales**, gate p **.014 at τ_min → .001 at the mesoscale**. The curve sits above its null everywhere; there is no characteristic scale — the trace is present at all of them. This is the *deep* trace.
- **α — single-scale (mesoscale).** Clears at **12/16 scales**, gate p **.024 → .007**, with the effect **peaking at s≈5** and fading at the fine and coarse ends. A real trace with a single characteristic scale — the sort of thing an ordinary measure can also catch (that is the point of the next section on encoding vs inference).
- **δ — fires-but-shallow (not claimed).** Significant only at scattered coarse scales; **fails both τ_min (.216) and the mesoscale (.116)**. It rides the coarse-scale collapse rather than holding a genuine hierarchy, so we do **not** claim it as a trace.
- **θ, low-γ — silent (0/16).** Nothing fires at any scale. Their silence is what makes the grading trustworthy: the coarse scales do *not* light everything up, so β's presence at every scale is a real property, not an artifact.
- **The one-liner:** only the τ-sweep can grade this. A single-scale measure would have called α and β the same band; the scale axis is what separates the genuinely-multiscale carrier (β) from the single-scale one (α). Information-depth = the *scale-structure* of the trace, not who detects it.

3. ON-SLIDE TEXT
grade each firing band by its ρ_sym(τ) curve — one matched-strength null, read per-scale
β  → clears at EVERY scale  (16/16,  p .014 → .001)  = scale-invariant · genuinely multiscale
α  → one mesoscale peak      (12/16,  p .024 → .007,  peaks s≈5)  = single-scale
δ  → scattered coarse scales  (fails τ_min .216, meso .116)  = fires-but-shallow, not claimed
θ, low-γ → silent (0/16)
⇒ only β is present at every scale — the deepest trace

4. SPEECH  (~80 s)
We know α and β carry the trace — that was selectivity. But "does a band fire" is a yes/no. The multiscale read gives us more than a yes/no: ρ_sym is a *curve* over diffusion scale, so we can ask a sharper question — *how* does each band fire? For every scale we hold the curve against its own matched-strength null, and we read it scale by scale, never cherry-picking the best one. Beta clears at every single scale — all sixteen — sharpening from p around fourteen-thousandths at the finest scale to a thousandth at the mesoscale. It has no characteristic scale: the trace is simply present everywhere you look. That's what "genuinely multiscale" means — the deepest kind of trace. Alpha is different: it clears at twelve of sixteen scales, but the effect is a single bump peaking around scale five and fading at both ends — one characteristic scale, real but resolvable. Delta only flickers at a few coarse scales and fails both the finest scale and the mesoscale, so it's shallow — we don't claim it. And theta and low-gamma are silent everywhere, which is exactly what tells us the coarse scales aren't just lighting everything up. So the grading is real: only beta is present at every scale. A single-scale measure would have called alpha and beta the same — the scale axis is what separates them.

Careful:
- The null is **matched-strength ONLY** — drift was RETIRED 2026-07-12 (a directional task makes a drift null degenerate with the trace). Do NOT reintroduce a drift column, a "second null", or a "survives both nulls" framing.
- Read **per-scale**. Never say "best scale" or collapse the curve to one number — the whole point is the *shape* of ρ_sym(τ).
- Do NOT say "it must be multiscale" or "only the multiscale read sees the trace" (retired overclaim — raw FC detects β at p=.024). The claim here is **scale-shape grading**: β is present at every scale, α at one. Frame it as *characterization*, not *exclusive detection*.
- Do NOT rank bands by any ×-null ratio; the gate p per scale is the read.
- δ is "fires-but-shallow, not claimed" — do NOT call it a trace, and keep it OFF the epilepsy story here (δ is the epilepsy band on slide 20; the band-dissociation belongs on slide 19, not this one).

5. FIGURES
- **Primary (the payoff):** `data/sparsified_arc/figures/trace_arc/trace_vs_tau_summary.pdf` — one panel per band: cohort-median observed ρ_sym(s) with IQR (band colour) over the matched-strength surrogate p50–p95 envelope (grey), firing scales (gate p<0.05) ticked. This IS the scale-shape figure: β's curve stays above its null across all τ (flat/broad), α's curve bumps at the mesoscale and falls off, δ only crosses at coarse scales, θ/low-γ never clear. gen `scripts/01_compute/sparsified_arc/fig_trace_arc.py` from `data/sparsified_arc/ms_mst020/cohort_gate.csv` (+ `trace_arc/cohort_gate_vs_s.csv`).
  - ⚑ **May need a talk-polished rebuild:** for a single glance, a 2-panel highlight showing **β (flat across τ)** next to **α (peaked at the mesoscale)** — same data, β-blue / α-purple from `visuals.styles.band_color`, PNG is fine for the talk. The 6-panel summary is busy for a slide.
- **Companion (whole landscape in one image):** `data/sparsified_arc/figures/trace_arc/trace_gate_heatmap.pdf` — band × scale grid coloured by −log10(cohort gate p): β = a solid row across all scales, α = a mesoscale block, δ = a coarse-only smear, θ/low-γ = empty. Good backup if the curves read as too busy.
- **RETIRED from this slide:** the controls-ladder / selectivity table (`data/sparsified_arc/figures/controls_ladder/controls_ladder_gate.pdf`) — it MOVED to **slide 14** (selectivity). Do not place it here; this slide is scale-shape only.

6. REFERENCES
- None new. The matched-strength null and the τ-swept cophenetic ρ_sym are ours (`13_matched_strength_mst020`). The dimensionless scale s = τ·λ_max is the LRG diffusion time (Villegas et al. 2023, *Nat. Phys.*; Villegas et al. 2025, *PRR*) — cite on the pipeline slide, not here.

7. CANVA STATUS
⚠ 2026-07-13 REFOCUSED (settled three results): slide re-tasked from the controls-ladder ("it must be multiscale") to the **SCALE-SHAPE pillar** — the τ-sweep grades each firing band, landing on "only β is present at every scale". The controls-ladder / selectivity content MOVED to slide 14. RE-RENDER on deck.

NEW — no Canva page yet. Build: the ρ_sym(τ) panel (β flat across τ vs α peaked at the mesoscale) + the grading text block (β 16/16 .014→.001 scale-invariant · α 12/16 .024→.007 single-scale · δ shallow/not-claimed · θ,low-γ silent) → punchline "only β is genuinely multiscale — present at every scale". Numbers are the settled-doc canonical set (`.agents/preprint/established_results/2026-07-13_settled-three-results.md`, R1 rows 1.2/1.3/1.7) and `data/sparsified_arc/ms_mst020/cohort_gate.csv`.
