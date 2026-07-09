---
name: talk-slide-I4-higher-order-and-multiscale
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-09
updated: 2026-07-09
slide: I-4
part: I — Introduction
duration: ~50 s
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
---

# I-4 · Higher-order — hypergraphs vs the forgotten multiscale

## Slide placeholder (copy into Canva)

**Main point.** Beyond the single pairwise edge lies higher-order structure. The
loud frontier builds **hypergraphs**; the under-explored route is the
**quantitative multiscale** one — organization at many scales at once. That's our
lane — and here's a glimpse of what it finds.

**Concepts to land.**
- Real interactions aren't only pairs — there is **higher-order** structure
  (groups, not just links). The visible push: **hypergraphs / simplicial
  complexes**.
- The forgotten route: **quantify the multiscale** — structure living at many
  scales **simultaneously**, with **no privileged resolution**. Loud on
  hypergraphs, quiet on multiscale.
- Light bridge: even ordinary **pairwise** data hides higher-order structure when
  you read it **across scales** — diffusion carries higher powers of the operator
  (one line; the mechanism lands in I-5). We do **not** build hypergraphs; we
  extract higher-order **graph** features from a pairwise operator.
- **This is our lane** — cognition may hide in these multiscale features,
  invisible to the single edge and to the fixed-*k* partition.
- **Flash-forward (light, don't explain):** a glimpse of what the multiscale read
  actually yields — a band-specific trace that lands in one cortical system. Pays
  off in Results.

**Figures / visuals.**
- Multiscale in one picture — the **same** network at two scales (5 vs 20
  communities):
  `data/outputs/figures/section3/fig_I/fig_I2_brain_communities_ImCoh_n5_n20.pdf`
- **Main-result teaser** (flash-forward — show, don't unpack):
  - multiscale earns its keep (raw vs multiscale):
    `data/preprint/figures/results_section1/fig_trace_c_raw_vs_multiscale.pdf`
  - the payoff glimpse (β trace → OFC on the brain):
    `data/reports/260707_meeting/figures/fig_trace_b_ofc_localization.pdf`
  - optional band map:
    `data/reports/rho_sym_band_map/fig_rho_sym_band_map.pdf`

**References.**
- Higher-order (Crossref-verified): Battiston 2020, *Phys. Rep.*
  (https://doi.org/10.1016/j.physrep.2020.05.004); Battiston 2021, *Nat. Phys.*
  (https://doi.org/10.1038/s41567-021-01371-4); Santoro et al. 2024,
  *Nat. Commun.* (https://doi.org/10.1038/s41467-024-54472-y).
- Multiscale — the forgotten quantitative route: Betzel & Bassett 2017,
  *NeuroImage* — "Multi-scale brain networks."
  https://doi.org/10.1016/j.neuroimage.2016.11.006

---

## Keep honest (content constraints, not styling)

- **Our "higher-order" ≠ hypergraphs.** It is higher *powers of a pairwise
  operator* → higher-order **graph** features ("emergent from pairwise
  connectivity"). Never call it "nonlinear," never imply we build a hypergraph.
- **"Most forgot multiscale" is rhetorical positioning** — say *under-explored
  relative to the hypergraph push*, not literally "everyone forgot."
- **Teasers stay light — don't explain them here.** If any number appears use the
  ρ_sym values (β p=0.032, 16.6× null); **OFC is a hotspot, not a container**;
  **no behavioral data**. Full unpack is Results (R-1 … R-7).
