---
name: talk-slide-15-survivors-table
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 15
status: draft
updated: 2026-07-11
canva: NEW — no Canva page yet
---

# Slide 15 — Who survives: it must be multiscale

1. TITLE
Who survives — it must be multiscale

2. MAIN CONCEPT
- Run the two nulls on all three measures. Who is robust to whom:
  - raw FC: ~half killed by matched-strength, and fails drift (raw whole-task is entirely drift) → out.
  - Grassmann (spectral): survives matched-strength, but is drift-vulnerable (does not clear the fair drift null, p ≈ 0.12, trends) → out.
  - cophenetic β (multiscale): survives matched-strength (23.7×) AND survives drift (C2, p = 0.0137) → sole survivor.
- Only the multiscale hierarchy carries a trace robust to both nulls. The pairwise and spectral tools fail — so the trace is not just structure, and not just time: it must be multiscale.

3. ON-SLIDE TEXT
                matched-strength   drift
raw FC          ~half killed       fails      → out
Grassmann       survives           vulnerable → out
cophenetic β    survives (23.7×)   survives   → SOLE survivor
⇒ it must be multiscale

4. SPEECH
Now run both nulls on all three measures. Raw edges: matched-strength already halves them, and drift finishes them — raw whole-task is basically just drift. Grassmann is tougher — it beats matched-strength — but it does not clear the fair drift null; it's drift-vulnerable. Only the cophenetic, multiscale trace in β survives both: strength-matched by a factor of twenty-plus, and drift-clean. So the two standard tools — pairwise and spectral — fail exactly where the multiscale one holds. That's the sentence: it must be multiscale.

Careful: say Grassmann is drift-*vulnerable* (fails the fair full-duration drift null, trends), not "washed out" — it does pass matched-strength. The dual-null survivor is cophenetic β specifically.

5. FIGURES
- The survivors table (build in Canva — method × null).
- Cophenetic ↔ Grassmann confirmation matrix — data/outputs/figures/talk/fig_coph_grassmann_confirmation_matrix.pdf

6. REFERENCES
- None (audit_167/169 + C2 drift null are ours).

7. CANVA STATUS
NEW — no Canva page yet. The old "gate" result has no standalone Canva page; build the table + the punchline here.
