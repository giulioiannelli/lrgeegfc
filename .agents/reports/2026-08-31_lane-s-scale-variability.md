---
name: 2026-08-31_lane-s-scale-variability
kind: report
era: PAPER_FINALIZATION (Wave 0, lane W0-S)
status: draft
created: 2026-08-31
scope: Does the cross-phase trace vary across diffusion scales, or is the incumbent readout blind to scale? Pre-registered before any number, run on the locked substrate knob-integrated over the plateau, everything referenced to a matched-strength null and gated by the locked cohort gate. RESULTS PENDING — the master grid is still computing.
pointers:
  - .agents/guides/task-persistence-investigation/2026-08-31_scale-local-trace-readouts.md
  - scripts/01_compute/paper_final/w0s_01_scale_locality_grid.py
  - data/paper_final/lane_s_scale/
---

# W0-S — is the scale axis flat, or is the readout blind to it?

*Draft. The master grid is still running; every results section below is empty until it lands, and nothing here may be cited yet.*

## The question, and why it was worth asking

Two W0-C numbers sit badly together. Across the swept scale axis the hierarchy changes enormously — the tree resolves 118 components at the finest scale and 4 at the coarsest, and the physical communication reach grows to ~92% of the implant span by `s ~ 2`. Yet the cross-phase statistic barely moves: the 28-point sweep is worth 1.2–1.9 independent tests, with mean cross-scale correlation +0.70 to +0.93. A statistic that stays nearly constant while the object it measures changes thirty-fold is a suspect statistic, so the project's conclusion "there is no multiscale structure" might really be "our readout cannot see scale".

The suspect mechanism, which this lane tests rather than assumes: `rho_sym` is a Spearman over all ~7 000 contact pairs at every scale, the cophenetic distance between two contacts is dominated by where they sit in the global tree, and that global ordering is largely preserved as the diffusion coarsens — so every scale is scored on nearly the same information and whatever is scale-specific is diluted by the bulk.

## Pre-registration

Written before any number and not editable in response to one: `.agents/guides/task-persistence-investigation/2026-08-31_scale-local-trace-readouts.md`. It fixes the notation, the two stratifications, the three candidate readouts, the three criteria a readout must satisfy to replace the incumbent, and the condition under which the dilution hypothesis is declared wrong.

<!-- RESULTS SECTIONS PENDING -->
