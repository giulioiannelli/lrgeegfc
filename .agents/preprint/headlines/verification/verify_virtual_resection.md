---
name: verify-virtual-resection
era: IMCOH_ABS_COHORT_N10
status: needs_new_audit
kind: verification-brief
headline: N3 (.agents/preprint/headlines/03_epileptogenic_markers.md)
owner_agent: epi-marker-trace-analysis
updated: 2026-06-22
---

# Verify — virtual resection of the epileptogenic diffusion community

**Head.** Mechanistic validation for N3: if the δ-band diffusion community
(`audit_89`/`audit_94`) is the epileptogenic backbone, **computationally removing
it should disrupt the network's diffusion organization more than removing a
size/strength-matched set of healthy nodes**. This is the "virtual resection"
(epilepsy Direction B, scoped 2026-05-08, never executed). It needs a new audit.

## 5-point critical preamble

1. **Claim.** Removing the predicted epileptogenic community (and/or clinical SOZ
   contacts) from the graph degrades a diffusion-organization metric (e.g.
   global diffusion efficiency / modularity of the propagator / Fiedler structure)
   **more** than matched control removals.
2. **Null.** Δ(metric) from resecting the SOZ community ≤ Δ from resecting
   **size-matched random** sets AND **strength-matched healthy** sets AND a generic
   **node-count decimation** of the same K.
3. **Strongest alternative.** (A) **Node-count/strength confound** — removing any
   K high-strength hubs perturbs the network; SOZ ≈ hubs in some patients
   (Pat_10/15). (B) **Trivial size effect** — bigger removals do more. (C) The
   metric is dominated by a few edges unrelated to epileptogenicity.
4. **Does the null control it.** vs (A): **strength-matched** healthy removal +
   the decimation control (`feedback_decimation_control_for_subset_exclusion`) — the
   SOZ effect must exceed strength-matched, not just random. vs (B): match K
   exactly. vs (C): test multiple diffusion metrics + per-patient, report which
   survive.
5. **Falsification + limits.** Falsified if SOZ resection ≈ strength-matched
   healthy resection. **Hard limit:** "virtual" only — **no surgical outcome
   (Engel/ILAE) to validate against** (unavailable); this validates the *community
   = backbone* claim, not clinical resection benefit. Two-population split
   (community vs hub patients) will likely show: validation strong in community
   patients, weak in hub patients — report per-patient.

## Steps

1. Define resection target per patient: (i) clinical SOZ contacts; (ii) the
   data-driven diffusion community (`audit_89`/`audit_94`).
2. Choose diffusion-organization metric(s) on ρ̂(τ): e.g. global diffusion
   efficiency, propagator modularity, spectral gap / Fiedler value, communicability.
3. Recompute metric after removing the target; compute Δ.
4. Build the controls: K-matched random removal, **strength-matched healthy**
   removal, decimation null (R≥200).
5. Cohort Wilcoxon on (Δ_SOZ − Δ_matched), LOO-max; per-patient table.

## Data / scripts

- Diffusion community: `audit_89`/`audit_94` outputs (see
  `../epi-marker-analysis.md` provenance; memory
  `epi_propagator_diffusion_community_2026_06_08`).
- Archived prior note: `.agents/reports/archive/2026-05/2026-05-08_trace-minus-epi-resection.md`.
- **New audit required** (none exists) — build on `build_operators` (audit_101) +
  the decimation control (`audit_85`). |ImCoh| substrate (proximity-immune).

## Pass / fail

- **Pass:** Δ_SOZ > Δ_strength-matched (LOO-robust) on ≥1 principled metric. →
  community = epileptogenic backbone; strong N3 mechanistic figure.
- **Fail:** Δ_SOZ ≈ Δ_strength-matched → the community is not a special backbone
  beyond hubness; report honestly, do not claim resection validation.

## Caveats

- No outcome labels → "virtual" claim only.
- Hub patients (Pat_10/15) expected to fail; per-patient mandatory.

## What to return

Verdict + the metric(s) that survive, the strength-matched comparison, per-patient
table, and a schematic (community removal → fragmentation). Fill below.

> **Result (agent fills):** …
