---
name: talk-slide-R3-beta-to-ofc
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-10
updated: 2026-07-10
slide: R-3
part: III — Results
duration: ~75 s
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
  - localization_audit_plan_2026_05_29.md
---

# R-3 · β → orbitofrontal cortex

## Slide placeholder (copy into Canva)

**Main point.** The β trace does not sit everywhere — it **concentrates in
orbitofrontal cortex**, the *only* brain system that survives the matched-strength
enrichment gate (BH q ≈ 0.01, all four conditions), while **PFC / sensorimotor are
depleted**. OFC is exactly where the I-2 bet said a cognitive map should consolidate
(Wilson 2014, Schuck 2016) — the payoff. **A hotspot, not a container.**

**Concepts to land.**
- **OFC is the only robust enrichment survivor.** Per-system matched-strength null
  (**R = 1000**) across the **9 a-priori systems** within β: OFC clears at
  **BH q ≈ 0.010–0.015 in all four conditions** (contact / shaft-collapsed ×
  epi-include / exclude), **leave-one-out robust** (survives dropping the β-anti
  Pat_10). **PFC and sensorimotor are the depleted pole** (robustly *anti*-enriched).
  A cohort enrichment gate, not a hand-picked region list.
- **Hotspot, not container.** OFC is **bilateral**; within-OFC only **~59 % of
  contacts are pro-trace**; OFC is implanted in **5/10 patients (4 positive)**.
  Read it as *concentration above baseline*, never "the trace lives only in OFC."
  The one marginal leave-one-out (q ≈ 0.06) is the known n = 5 coverage fragility.
- **Strength-independent.** The OFC survivors are **low-strength, non-hub** contacts
  (negative strength deviation) — the address is not a byproduct of where
  connectivity is strongest. Same matched-strength logic as M-4, now made spatial.
- **The bet pays off.** OFC is the cognitive-map / task-space hub (Wilson 2014,
  Schuck 2016). A relational map, consolidated offline, landing in OFC is the
  anatomically expected address — **structure + anatomy + literature converging**,
  not a brain–behaviour correlation.

**Figures / visuals.**
- **Load-bearing — systems-enrichment glass brain** (three views). Contacts coloured
  by verdict: **OFC = green (enriched, robust)**, **PFC / sensorimotor = orange
  (depleted, robust)**, pale = single-patient-driven, grey = n.s. / white matter.
  Green concentrates ventral-frontal and **bilateral**; orange over dorsal PFC / SM.
  Path: `data/preprint/figures/results_section1/fig_trace_b_ofc_localization.pdf`
  (polished alt: `data/preprint/figures/beta/anatomy/fig_beta_anatomy_brain.pdf`).
  *Currency note:* §1 Fig 2 was folded into a compound on 2026-07-10, so this standalone
  on disk is pre-refactor (message-correct, but rebuild `fig_trace_b_ofc_localization.py`
  for a fresh cut before dropping into Canva).
- **Flashy talk asset (optional) — rotatable 3D β-vs-α brain.** Contacts as spheres
  sized by |ρ_sym| (orange = carrier, blue = anti, grey = the "sea"), **OFC green
  halo**, green links = co-movements reaching OFC; header reads **β movers reach OFC
  ×1.5 vs α ×1.1**. Node-driven (arcs are illustration — see *Keep honest*). Paths:
  `data/preprint/figures/_drafts/fig_beta_ofc_pairglow_3d_DRAFT.html` (interactive) ·
  `..._DRAFT.png` (still). Builder:
  `scripts/01_compute/figures_embedded/fig_beta_ofc_pairglow_3d.py`.

**References.**
- Wilson, Takahashi, Schoenbaum & Niv (2014). Orbitofrontal cortex as a cognitive
  map of task space. *Neuron* 81:267–279. https://doi.org/10.1016/j.neuron.2013.11.005
- Schuck, Cai, Wilson & Niv (2016). Human orbitofrontal cortex represents a cognitive
  map of state space. *Neuron* 91:1402–1412. https://doi.org/10.1016/j.neuron.2016.08.019

---

## Keep honest (content constraints, not styling)

- **The β → OFC address is system/node-level matched-strength enrichment (q ≈ 0.01),
  NOT per-pair arcs.** OFC-incidence among top-mover pairs enriches only
  **×1.46 (β) vs ×1.12 (α)** and is patient-variable (OFC-frac 0.16 ± 0.21) — so the
  green arcs in the 3D figure are **illustration**, never the evidence. If pressed,
  cite the systems enrichment (the glass brain), not the arcs.
- **Hotspot, not container** (bilateral, ~59 % within-OFC pro, 5/10 implanted) —
  concentration above baseline, not "only in OFC." (Structure-file accuracy flag 5.)
- **Hippocampus / MTL is a sub-threshold hint** (matched-strength p ≈ 0.4,
  LOO-fragile, Pat_02-driven) — real but **not the headline**. Don't upgrade it to a
  co-localization.
- **Shaft-collapse robust** → the OFC concentration is *not* electrode-adjacency
  inflation. A distributed "paralimbic ring" reading failed shaft-collapse and was
  retracted; the verdict flipped three times before shaft-collapse settled it on OFC.
  Don't reopen the ring.
- **α does not localize** (no FDR-surviving system) — the OFC address is **β-specific**
  here. Don't imply α shares it. (α's lack of an anatomical home is the band-taxonomy
  story, not R-3.)
- Numbers are **ρ_sym-era**; **no behavioural correlation** — the "cognitive map"
  reading is licensed by anatomy + literature, never by performance.
