---
name: talk-slide-15-survivors-table
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 15
status: draft
updated: 2026-07-13
canva: NEW — no Canva page yet
---

# Slide 15 — Who survives: it must be multiscale

1. TITLE
Who survives — it must be multiscale

2. MAIN CONCEPT
- The trace clears matched-strength. But is it really *multiscale* — or could a simpler, single-scale measure have found the same thing? Run a **ladder of read-outs on the same sparsified graph**, each through its own matched-strength gate, and see which ones actually *carve the bands*:
  - **raw FC (pairwise edges)** — fires in *four* bands (δ/α/β/low-γ): it registers the reorganization but cannot tell a trace band from a non-trace band. Non-selective.
  - **node strength** — degenerate: the strength-matched surrogate preserves it exactly, so it can never clear its own null (p = 1).
  - **clustering coefficient** — catches α only, misses β.
  - **graph geodesic** — fires the *wrong* bands (δ, low-γ); misses both carriers.
  - **effective resistance (spectral)** — a global, multi-step Laplacian read-out, and yet *multiscale-blind* — **dead in every band**. This is the sharp one: a genuinely spectral method that *should* be the strongest rival finds nothing.
  - **cophenetic (multiscale, τ-swept)** — the ONLY read that is band-selective: α and β and nothing else, sharpening from the fine scale to the mesoscale (β p 0.014 → 0.001; α 0.019 → 0.007).
- Only the multiscale hierarchy separates the bands. Pairwise, network, geodesic, AND spectral tools all fail — so the trace is not just "where connectivity is strong", not just a graph statistic, and not just a spectral pattern: **it must be multiscale**, and the thing that carves the bands is the scale-sweep in τ.

3. ON-SLIDE TEXT
same graph, one null (matched-strength) — which read carves the bands?
raw edges              → fires 4 bands (non-selective)
node strength          → degenerate (p = 1)
clustering             → α only
geodesic               → wrong bands (δ, low-γ)
spectral (resistance)  → dead everywhere
cophenetic (multiscale, τ-swept) → α, β  ✓ selective
⇒ only the multiscale read carves the bands → it must be multiscale

4. SPEECH  (~75 s)
The trace clears matched-strength. But here's the question that actually matters for this talk: is it really *multiscale*, or could some ordinary measure have found the same thing? So we line up a ladder of read-outs on the very same network and hold each to the same strength-matched null. The raw edges light up in four bands — they see a reorganization but can't separate a band that carries the trace from one that doesn't. Node strength is degenerate — the null preserves it, so it can't clear anything. Clustering catches only alpha. The graph geodesic fires the wrong bands entirely. And then the sharpest test: the effective-resistance distance — a genuinely spectral, global read-out — is *dead in every band*, because it's still blind to scale. Only the cophenetic hierarchy, read across the sweep of diffusion times, carves the bands cleanly: alpha and beta, and it sharpens as you go to the mesoscale. Pairwise, network, geodesic, spectral — all fail exactly where the multiscale read holds. That's the sentence: it must be multiscale.

Careful: the null is matched-strength ONLY (drift RETIRED 2026-07-12 — a directional task makes a drift null degenerate with the trace). Do NOT reintroduce a "drift" column or a "survives both nulls" framing. The killer comparison is the SPECTRAL one — effective resistance = the Laplacian pseudo-inverse distance: it is multi-step and global yet multiscale-blind, so its failure is what licenses "not just spectral — multiscale". Grassmann (slide 13) is a *different*, complementary spectral read (global modes) that DOES survive matched-strength for β; keep it OFF this slide to avoid "two spectral things" — this slide's spectral foil is resistance.

5. FIGURES
- The controls-ladder table (build in Canva — descriptor × band, matched-strength gate; only the cophenetic rows are band-selective; the resistance row is empty). Source numbers: `data/sparsified_arc/controls_ladder/cohort_gate.csv` (functional = T_test), gen `scripts/01_compute/sparsified_arc/14_controls_ladder_mst020.py`; preview figure `scripts/.../fig_controls_ladder.py`.
- (optional) the cophenetic↔Grassmann confirmation matrix — data/outputs/figures/talk/fig_coph_grassmann_confirmation_matrix.pdf — only if slide 13 didn't already use it.

6. REFERENCES
- None (controls-ladder audit `14_controls_ladder_mst020` + matched-strength null are ours; effective resistance = Laplacian pseudo-inverse, standard).

7. CANVA STATUS
NEW — no Canva page yet. ⚠ 2026-07-13 REBUILT for the mst@0.20 recovery: the old "3 measures × 2 nulls (matched-strength + drift) → cophenetic β sole survivor" table is RETIRED (drift discarded — task is directional). The slide is now the CONTROLS LADDER — raw / strength / clustering / geodesic / spectral-resistance / cophenetic, one matched-strength null — landing on "only the multiscale read carves the bands". Build the ladder table + the punchline here; numbers in the OVERVIEW directive (`.agents/preprint/directives/writing_directive_2026-07-13_sparsified-recovery-OVERVIEW.md`) and `controls_ladder/cohort_gate.csv`.
