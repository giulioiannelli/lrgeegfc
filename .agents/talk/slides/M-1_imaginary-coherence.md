---
name: talk-slide-M1-imaginary-coherence
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-09
updated: 2026-07-09
slide: M-1
part: II — Methods
duration: ~55 s
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
---

# M-1 · Why imaginary coherence

## Slide placeholder (copy into Canva)

**Main point.** The killer methods slide: intracranial FC is dominated by
**volume conduction** — instantaneous, zero-lag signal leakage. **Imaginary
coherence is blind to it by construction** (a zero-lag mix is purely real → its
imaginary part is exactly zero), so the hierarchy reflects brain coupling, not
electrode geometry.

**Concepts to land.**
- **The problem:** intracranial FC is swamped by **volume conduction** — one
  source seen at the same instant by nearby contacts (zero-lag leakage), not a
  true interaction.
- **The fix:** **imaginary coherence** ignores the zero-lag part **by
  construction** — an instantaneous mix is purely real, so its imaginary part is
  exactly zero. Only genuine **time-lagged** coupling survives.
- **Standard coherence (MSC)** uses the full cross-spectrum → **captures the
  artifact** → **same-probe** coupling inflated, driving the hierarchy by
  electrode geometry, not brain function.
- **What we feed the pipeline:** `imcoh_abs = ⟨|ImCoh(f)|⟩` over in-band bins —
  non-negative (needed for the Laplacian); absolute value **per bin first**, then
  band-average.
- **The number:** same-vs-cross-probe coupling ratio drops **5.5× (MSC) → 1.13×
  (imcoh_abs)** — the geometry bias is essentially gone.

**Figures / visuals.**
- **How ImCoh is (across frequency) + the bands we average over** — the
  cross-spectral distribution (the "imcoh notes" figure: per-patient, all pairs,
  same-shaft vs off-shaft, δ θ α β γ γ_h shaded). MSC vs |ImCoh| makes the
  volume-conduction point directly — **MSC carries a big same-shaft spike;
  |ImCoh| suppresses it**:
  - (a) MSC: `data/outputs/figures/section2/fig_C/fig_C2_spectral_distribution_msc_rsPre.pdf`
  - (b) |ImCoh|: `data/outputs/figures/section2/fig_C/fig_C2_spectral_distribution_imcoh_rsPre.pdf`
  - cleaner single-patient variant: `…/fig_C2_spectral_example_Pat_05_rsPre.pdf`
- **The band-average per pair (mechanism companion)** — signed ImCoh(f) → |·| per
  bin → in-band mean = the one number we feed the pipeline:
  `data/outputs/figures/talk/imcoh_band_average.pdf` (built fresh; script
  `scripts/01_compute/figures_embedded/fig_imcoh_band_average_demo.py`).
- **Downstream contrast (dendrogram level)** — the bias propagates to the tree:
  `data/outputs/figures/section3/fig_H/fig_H1_dendrograms_MSC_vs_ImCoh_3patients.pdf`.

**References.**
- Nolte et al. 2004, *Clin. Neurophysiol.* — imaginary part of coherency.
  https://doi.org/10.1016/j.clinph.2004.04.029
- Bastos & Schoffelen 2016, *Front. Syst. Neurosci.* — FC methods & pitfalls.
  https://doi.org/10.3389/fnsys.2015.00175

---

## Keep honest (content constraints, not styling)

- **ImCoh removes zero-lag leakage — it is not a panacea.** It also **attenuates
  genuine zero-lag coupling** (a known, accepted cost). Frame it as *the
  volume-conduction control*, not "true connectivity."
- **imcoh_abs is second-order** (a linear / pairwise measure). The higher-order
  structure comes from the **LRG graph read**, not from ImCoh being "nonlinear."
  Never call it nonlinear.
- The **|·|-then-average** order matters (Jensen: `|mean| ≠ mean|·|`); imcoh_abs
  is a magnitude, kept non-negative for the Laplacian.
- **Same-probe bias is THE artifact to defeat** here; the 5.5× → 1.13× ratio is
  the headline evidence the control works.
