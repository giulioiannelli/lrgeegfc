---
name: talk-slide-M2-frequency-bands
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-09
updated: 2026-07-10
slide: M-2
part: II — Methods
duration: ~40 s
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
---

# M-2 · Why frequency bands

## Slide placeholder (copy into Canva)

**Main point.** Different rhythms index different processes, so
**band-specificity is both the signal and a control**: a real effect is
frequency-selective; an artifact is band-blind. And in multiscale terms,
splitting into bands **disentangles the temporal axis** of brain organization —
the time-scale complement to the spatial / topological multiscale the LRG
hierarchy reads. We analyse six bands separately rather than broadband — which
would wash the signal band into the null.

**Concepts to land.**
- Different **rhythms** index different processes → analyse **per band**, not
  broadband.
- **Band-specificity is signal *and* control:** a real trace is
  **frequency-selective**; an artifact shows up **everywhere** (band-blind).
- **Band-splitting = the temporal axis of the multiscale.** Betzel & Bassett's
  multi-scale organization has several axes; **one is time**. Frequency bands
  *are* temporal scales, so separating them disentangles the multiscale along its
  **temporal** dimension — the complement to the spatial / topological multiscale
  the LRG diffusion hierarchy reads (I-4 / I-5).
- The six bands:
  **δ 0.5–4 · θ 4–8 · α 8–13 · β 13–30 · γ_low 30–80 · γ_high 80–300 Hz**
  (γ_high is **intracranial**, *not* the textbook 30–80).
- **Broadband would wash it out** — averaging the signal band (β) together with a
  null band (θ).
- Foreshadow (honest): at the **raw-edge** level the drift is positive in *every*
  band, θ included (band-blind) — specificity only **emerges** through the
  multiscale read (→ R-1).

**Figures / visuals.**
- A **6-band row of FC matrices** — same brain, six bands, the structure differs
  by band:
  `data/outputs/figures/fc_templates/row_per_band/Pat_05_rest_pre_imcoh_abs_chnames_log_shared.pdf`
- The **shaded band regions in the M-1 cross-spectral figure do double duty** here
  (they *are* the six bands).

**References.**
- Multiscale / temporal axis: Betzel & Bassett 2017, *NeuroImage* — "Multi-scale
  brain networks" (one axis is temporal → bands).
  https://doi.org/10.1016/j.neuroimage.2016.11.006 (also cited at I-4).
- Band definitions are otherwise standard; **γ_high 80–300 Hz is an intracranial
  convention** (high-gamma / ripple range), stated as ours.

---

## Keep honest (content constraints, not styling)

- **γ_high = 80–300 Hz is intracranial, explicitly NOT the textbook 30–80.** Say
  so on the slide.
- **Band-specificity is *earned*, not raw.** Be honest that raw edge-level drift
  is **band-blind** (positive everywhere incl. θ); the frequency selectivity
  appears only through the multiscale hierarchy (R-1). Don't claim
  band-specificity at the raw-FC level.
- **Don't preview which bands win** — the gate and the ρ_sym verdicts are R-2.
