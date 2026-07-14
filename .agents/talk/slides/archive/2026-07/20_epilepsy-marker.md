---
name: talk-slide-20-epilepsy-marker
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 20
status: draft
updated: 2026-07-13
canva: page 19 · ~30% (cuttable; ⚠ 2026-07-13 REFOCUSED to settled R3 — updated AUCs + two-populations + LOPO detector)
---

# Slide 20 — The same operator marks epileptogenic tissue

1. TITLE
The same operator marks the seizure zone

2. MAIN CONCEPT
- Head (the juice): the *same* sparsified diffusion propagator that carries the cognitive trace — no retuning, read within a single recording as seed affinity — localizes epileptogenic tissue as a **strength-independent co-diffusing community**. Cognitive lens, clinical dividend, for free.
- SOZ is **relational, not per-contact**: the seizure-onset contacts co-diffuse as a group above their node-strength baseline — it is not "where connectivity is strong."
- Seed-affinity AUC, strongest in δ: **δ 0.83 · low-γ 0.82 · β 0.745 · α 0.61**. Each band beats a **matched-strength fake-SOZ null** (the sole null): δ 7/10, low-γ 8/10, β 8/10 patients above their own null.
- **Band dissociation** (the R2→R3 bridge): cognition lives in **α/β**, epilepsy in **δ/low-γ/β**, and **β is the bridge**. δ marks the seizure zone but carries no cognitive trace. Over the SOZ the cognitive bands split too: **β spares the SOZ** (the trace lives off epileptogenic tissue), **α recruits it** (α SOZ–SOZ +0.41, p=.005).
- τ (the sweep) lifts **ranking, not precision**: multiscale beats single-scale AUC in 7–8/10 patients, but top-of-list precision comes from **fusing the bands**, not from τ — single-band prec@5 = .40, **fused prec@5 = 60% (~7× base rate)**. So τ buys ordering; band fusion buys the hit-rate.
- **Two populations**: 8 patients are the co-diffusing community (community-level AUC up to **0.99**); 2 are right-hemisphere hub cases where the community read fails (0.49 / 0.57) and node strength rescues them. A **calibrated leave-one-patient-out detector** lands at **median AUC 0.87, 9/10 above chance**.
- **Robust to sparsification**: the marker reproduces dense → `mst@0.20` (δ .80→.83, low-γ .74→.82, β .69→.745) — a property of the operator, not of how the graph was thinned.
- Scope: this reads clinically-labelled **SOZ, not surgical outcome** — a **triage marker, not a detector**.

3. ON-SLIDE TEXT
Same operator, read within one recording → the seizure zone
SOZ = a strength-INDEPENDENT co-diffusing community (relational)
seed-affinity AUC: δ 0.83 · low-γ 0.82 · β 0.745 · α 0.61
beats matched-strength fake-SOZ null: δ 7/10 · low-γ 8/10 · β 8/10
β spares the SOZ · α recruits it (α SOZ–SOZ +0.41, p=.005)
cognition = α/β · epilepsy = δ/low-γ/β · β = the bridge
calibrated LOPO detector: median AUC 0.87 (9/10 > chance)
(SOZ labels, not surgical outcome — a triage marker)

4. SPEECH
And the same tool, with nothing retuned, pays a second dividend — a clinical one. Read inside a single recording as diffusion affinity from a seed, the seizure-onset contacts light up as their own co-diffusing community. Crucially it's strength-independent: it beats a matched-strength null that scrambles the SOZ but keeps every node's connectivity, so this is genuine relational geometry, not just "where the graph is strong." The strongest band is delta, area under the curve 0.83, then low-gamma 0.82, then beta 0.745 — and each beats its own fake-SOZ null in most patients. Notice the band dissociation that ties the whole talk together: cognition lives in alpha and beta, epilepsy in delta, low-gamma and beta — beta is the bridge. And over the seizure zone the cognitive bands even split: beta spares it, alpha recruits it. Sweeping tau reorders the list but doesn't sharpen the top of it — precision comes from fusing the bands, sixty percent at the top five, about seven times chance. There are really two populations: eight patients are this diffusion community, up to 0.99 at the community level; two are right-hemisphere hubs where plain strength does the work. Calibrated leave-one-patient-out, the detector sits at a median AUC of 0.87, nine of ten above chance. One operator, one null, a cognitive lens that doubles as a clinical one.

Careful (what NOT to say): scope is clinically-labelled SOZ, NOT surgical outcome — this is a triage marker, not a diagnostic. Say "the SOZ is a strength-INDEPENDENT co-diffusing community" (relational), never a per-contact classifier. The clinical bands are δ/low-γ/β; do NOT call δ a cognitive band (δ carries no trace — that's the dissociation, not a contradiction). τ = ranking, NOT precision; precision is from band FUSION — don't credit τ with the hit-rate. The null is matched-strength ONLY — no "drift", no "second null", no ranking bands by ×-null ratio. Read AUCs as reported (per-band), not best-scale. Cuttable if you run long.

5. FIGURES
- SOZ = strength-independent diffusion community — data/outputs/figures/talk/fig_soz_diffusion_community.pdf (still valid — relational-community layout)
- β spares / α recruits the SOZ — data/outputs/figures/talk/fig_trace_f_soz_divergence.pdf (still valid)
- Relational marker + per-band AUC vs fake-SOZ null — data/reports/results_section3/fig_epi_a_relational_marker.pdf ⚠ REBUILD: script docstring + panel still carry the DENSE AUCs (δ .80 / low-γ .74 / β .69); re-render on `mst@0.20` to the settled δ .83 / low-γ .82 / β .745.
- Two populations (community vs right-hemi hub) — data/reports/results_section3/fig_epi_c_two_populations.pdf (verify it shows community AUC↑.99 + the 2 hub cases .49/.57)
- Calibrated LOPO detector — data/reports/results_section3/fig_epi_b_calibrated_detector.pdf (verify median AUC .87, 9/10 > chance)
- Scope / ceiling (τ = ranking not precision; fused prec@5 = 60%) — data/reports/results_section3/fig_epi_d_scope_ceiling.pdf
- PNG is fine for the talk; export whichever compound reads cleanest at slide size.

6. REFERENCES
- Ours (audit_132 all-contacts propagator SOZ marker; audit_89 diffusion community; mst@0.20 recovery arc 2026-07-12).

7. CANVA STATUS
⚠ 2026-07-13 REFOCUSED (settled three results): updated AUCs to δ .83 / low-γ .82 / β .745, added strength-independent relational community + fake-SOZ null (δ7/low-γ8/β8), τ=ranking-not-precision + band-fusion prec@5 60%, two populations (community↑.99 + 2 right-hemi hubs), LOPO detector median .87, dense→mst@0.20 robustness, and the α/β·δ-low-γ-β·β-bridge dissociation. RE-RENDER on deck.
Page 19 · ~30% ("Coda N3"). Have: the concept + scope text. Missing: the four section-3 figures (rebuild fig_epi_a to mst@0.20 AUCs first), compressed on-slide text, presenter notes. Cuttable if running long.
