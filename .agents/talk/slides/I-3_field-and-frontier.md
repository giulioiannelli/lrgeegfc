---
name: talk-slide-I3-network-neuroscience-and-fc
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-09
updated: 2026-07-09
slide: I-3
part: I — Introduction
duration: ~40 s
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
---

# I-3 · Network neuroscience & functional connectivity

## Slide placeholder (copy into Canva)

**Main point.** Why we build networks at all: complex systems — the brain among
them — are understood through their web of **interactions**, and **functional
connectivity** is how neuroscience builds that web from brain signals (which
regions *fluctuate together*), giving network science a handle on brain
organization.

**Concepts to land.**
- **Network science, in general:** a complex system is units + interactions; the
  network view exposes organization — **hubs, communities, pathways** — you can't
  see one unit at a time. Same toolkit across physics, biology, society.
- **Why functional connectivity:** in the brain you rarely watch the wiring
  *work*; instead you infer **coupling from signals** — which regions **fluctuate
  together** — and get a **functional network**.
- **What it buys you:** read brain organization straight from recordings, compare
  states / conditions, find **modules and hubs**, relate structure to cognition
  and disease.
- One image: a **coupling matrix ↔ a graph** — the whole idea at a glance.
- Network neuroscience is a **large, mature** field built on exactly this.

**Figures / visuals.**
- FC in one image — coupling **matrix ↔ network**:
  `data/outputs/figures/network_templates/matrix_plus_network/Pat_05_beta_rpre_imcoh_abs_lrg_sfdp_cmap_g2.pdf`
  (present **method-neutrally** — just "a functional network").
- Brain-on-cortex connectome exists as **interactive HTML**
  (`data/outputs/figures/brain_connectome_multiscale/Pat_05/`) — static export
  available on request.

**References.**
- Bullmore & Sporns 2009, *Nat. Rev. Neurosci.* — complex brain networks (the
  field's canonical reference). *DOI to verify — not yet in our Crossref list.*
- Bassett & Sporns 2017, *Nat. Neurosci.* — network neuroscience.
  https://doi.org/10.1038/nn.4502
- Optional "network science in general" nod: Barabási, *Network Science* (2016),
  or Newman, *Networks*.

---

## Keep honest (content constraints, not styling)

- **FC = statistical coupling, not wiring and not causation.** Say "fluctuate
  together," never "connected by fibres" or "A drives B." Structural and
  effective/causal connectivity are different things.
- Present the exemplar **method-neutrally** — don't name ImCoh (M-1) or LRG (M-5)
  yet, and don't imply one FC measure *is* the field.
- Keep this slide purely foundational — the higher-order / multiscale framing is
  the next slide (I-4).
