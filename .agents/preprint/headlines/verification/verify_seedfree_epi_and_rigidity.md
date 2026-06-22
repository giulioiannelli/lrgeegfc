---
name: verify-seedfree-epi-and-rigidity
era: IMCOH_ABS_COHORT_N10
status: live_audit_118_119_120_121
kind: verification-brief
headline: N3.3 (.agents/preprint/headlines/03_epileptogenic_markers.md) + candidate sub-result
owner_agent: epi-marker-trace-analysis
updated: 2026-06-22
---

# Verify — seed-free SOZ markers + cross-phase routing rigidity

**Head.** Two live epilepsy directions that could strengthen N3. **(A) Seed-free
markers** (`audit_118/119/120`): can the epileptogenic community be found *without*
known SOZ seeds? — the clinically valuable, harder case (prior label-free work,
`audit_94`, failed cross-patient selection at chance). **(B) Cross-phase routing
rigidity** (`audit_121`): is epileptic-zone diffusion routing more *rigid* across
phases (rest/task) than healthy routing? — a candidate new sub-result. Both are
running; adjudicate honestly.

## Part A — seed-free SOZ marker

### 5-point preamble
1. **Claim.** A seed-free read-out (path routing / routing-module structure,
   `audit_118/119/120`) ranks/narrows SOZ contacts above chance — at least
   *within* patient, ideally *across* patients.
2. **Null.** Seed-free ranking ≤ chance AND ≤ a node-strength baseline; cross-
   patient *selection* of the pathological module ≤ chance.
3. **Strongest alternative.** (A) **Hubness** — seed-free structure just recovers
   high-strength hubs. (B) **The audit_94 failure mode** — multiple physiological
   modules look equally "pathological" with zero labels, so cross-patient
   selection is at chance even if within-patient narrowing works.
4. **Does the null control it.** Strength-residualize; report **within-patient**
   vs **cross-patient** separately (the honest split); off-shaft / leave-one-shaft
   evaluation (proximity is removed by |ImCoh| but shaft-clustering of labels is
   not).
5. **Falsification + limits.** Falsified (for the cross-patient claim) if selection
   is at chance. Expected honest outcome: seed-free is **less precise** than the
   seed-based 6-band detector (N3.2); frame as "narrows, does not localize."

### Steps / data
- Run/inspect `audit_118_epi_seedfree_marker.py`, `audit_119_epi_path_routing.py`,
  `audit_120_epi_seedfree_routing_module.py` (2026-06-22, live).
- Compare to seed-based detector `data/audit/epi_propagator_detector/` (`audit_117`).
- Metrics: within-patient AUC (strength-residualized), cross-patient selection
  accuracy, off-shaft precision; strength baseline + label-shuffle null.

### Pass / fail
- **Pass (modest):** within-patient seed-free narrowing beats strength + chance →
  N3.3 becomes a real "seed-free narrows SOZ" result (with the precision-loss
  caveat).
- **Fail:** cross-patient at chance and within-patient ≈ strength → report the
  negative; N3 stays seed-based.

## Part B — cross-phase routing rigidity

### 5-point preamble
1. **Claim.** Diffusion routing through/around the epileptic zone changes **less**
   across phases (rest_pre→task→rest_post) than healthy-tissue routing — epileptic
   routing is **rigid**.
2. **Null.** Epi rigidity ≤ matched-strength surrogate AND ≤ strength-matched
   healthy nodes.
3. **Strongest alternative.** Hubs are rigid regardless of pathology
   (strength/hubness confound); rigidity is just low temporal SNR in some contacts.
4. **Does the null control it.** Strength-matched healthy comparison + matched-
   strength surrogate; per-band; LOO-max.
5. **Falsification + limits.** Falsified if epi rigidity ≈ strength-matched
   healthy. Limits: n=10; "rigidity" definition must be pre-registered (e.g. cross-
   phase variance of a per-node routing vector) before fishing.

### Steps / data
- Run/inspect `audit_121_epi_crossphase_routing_rigidity.py` (2026-06-22, live).
- Define routing rigidity (cross-phase stability of a propagator routing measure);
  compare epi vs strength-matched healthy; matched-strength null; per band; LOO.

### Pass / fail
- **Pass:** epi routing rigid beyond strength (LOO-robust) → candidate new N3
  sub-result (mechanistic: the seizure network is a fixed conduit).
- **Fail:** ≈ strength → drop.

## Caveats (both parts)

- Clinical-label not outcome (no Engel/ILAE).
- Two-population split (community vs hub patients) — per-patient reporting.
- Matched-strength / strength-residualization mandatory; hubness is the recurring
  confound in all epi-marker work.

## What to return

Two verdicts (A seed-free: narrows/within-only/fails; B rigidity: pass/fail) with
the strength controls, within-vs-cross-patient split, LOO-max, and whether either
becomes an N3 sub-result. Fill below.

> **Result A (agent fills):** …
> **Result B (agent fills):** …
