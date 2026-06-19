---
name: cross-phase-taxonomy-results
type: report
era: IMCOH_ABS / COHORT_N10
status: current
date: 2026-06-12
created: 2026-06-12
updated: 2026-06-12
scope: >
  Results of the anchor/trace/reset/reorganize taxonomy built as a cut-free
  decomposition of the rho^coph signal. Steps 1-3 (decomposition, coupled-null
  gate, geometry-baseline localization). Trace verified; reset/anchor unverified.
pointers:
  - .agents/guides/task-persistence-investigation/2026-06-12_cross-phase-cophenetic-taxonomy.md  # scope
  - src/lrg_eegfc/utils/metrics/cross_phase.py                          # decomposition lib
  - src/lrg_eegfc/utils/surrogate/matched_strength.py                  # coupled_surrogate_cophenet
  - scripts/01_compute/audit/audit_105_cross_phase_taxonomy.py           # step 1
  - scripts/01_compute/audit/audit_106_coupled_null_validation.py       # step 2 (gate)
  - scripts/01_compute/audit/audit_107_taxonomy_localization.py         # step 3 (geometry baseline)
  - scripts/01_compute/audit/audit_110_taxonomy_matched_strength.py     # step 4 VERIFICATION (kills reset/anchor)
  - scripts/01_compute/audit/audit_109_taxonomy_brain_figure.py         # brain figure
  - data/audit/cross_phase_taxonomy/                                    # all outputs
  - localization_audit_plan_2026_05_29.md                              # locked beta->OFC trace
---

# Cross-phase taxonomy — anchor / trace / reset / reorganize (results)

**Head (updated 2026-06-18).** The per-pair cophenetic signal that carries the
validated β trace decomposes — exactly, with no residual — into a rigid **anchor**
backbone plus two orthogonal fluctuation contrasts, **trace** (persistence φ₁, =
ρ_split) and **reset** (excursion φ₂). The decomposition is sound (reproduces the
locked trace bit-exact in 60/60 cells; trace energy peaks at β; anchor dominates —
the Gratton stable-backbone picture). The localization looked like a **3-way β
dissociation** (trace→OFC, reset→lateral-temporal, anchor→insula) under the
geometry baseline — **but that was retracted by the mandatory matched-strength
test (audit_110, §5): only trace → OFC survives the strength null.** reset and
anchor were **node-strength artifacts** of the geometry baseline (which controls
contact count, not strength) — the KC-style failure the matched-strength rule
exists to catch. **Publishable result: none that is new** — trace→OFC is the
already-locked one. What's left is a clean *framing* (the flow decomposition) and
the honest node-flow conclusion that **only the persistent component is
anatomically specific beyond connectivity strength; reversible/rigid flow is
strength-organized**. The coupled-cross-phase null also failed its own gate
(confirming audit_63), so no temporal null is available either.

---

## 1. What was built (all cut-free, ρ^coph-native; no tree-cutting)

- **Scope report** `…/2026-06-12_cross-phase-cophenetic-taxonomy.md` (11-section,
  5-point preamble) — written before code.
- **Library** `utils/metrics/cross_phase.py` — `phase_contrasts` (φ₁/φ₂/var, exact
  split `var=(φ₁²+φ₂²)/3`), `cross_phase_channels`, `cohort_share_row`.
- **`coupled_surrogate_cophenet`** in `utils/surrogate/matched_strength.py` — the
  shared-swap-sequence coordinated null (+ independent variant).
- **audit_105** decomposition + trace guard; **audit_106** coupled-null validation
  gate; **audit_107(/b)** geometry-baseline localization + figures.

## 2. Decomposition (audit_105) — VERIFIED trace channel

**Trace guard: 60/60 cells reproduce the locked `obs_rho` bit-exact** (tol 1e-6).
The decomposition cannot and does not perturb the trace (reset ⟂ trace, φ₂ in the
null-space of ρ_split — a perfect reset has Δrest=0).

Cohort composition (median over n=10; `comp_*` is mass×energy, sums to 1):

| band | anchor | trace | reset | among-mover share_T | mover% | r(anchor,strength) |
|---|---|---|---|---|---|---|
| δ | 0.873 | 0.066 | 0.057 | 0.62 | 13% | +0.19 |
| θ | 0.857 | 0.039 | 0.038 | 0.46 | 14% | +0.09 |
| α | 0.918 | 0.024 | 0.059 | 0.42 | 8% | +0.05 |
| **β** | 0.603 | **0.306** | 0.108 | **0.57** | 40% | +0.19 |
| low-γ | 0.828 | 0.032 | 0.067 | 0.38 | 17% | +0.38 |
| high-γ | 0.192 | 0.036 | 0.185 | 0.50 | 81% | +0.26 |

Readings (both scope-report guards pass):
- **Trace energy peaks at β** (0.306 vs ≤0.066 elsewhere) — consistent with the
  locked verdict that β is the trace band.
- **Anchor dominates every band** (0.60–0.92) except **high-γ (0.19, 81% movers)**
  — high-γ is the noise band: everything moves, no persistence direction
  (share_T 0.50, balanced) → correctly flagged as reorganize/noise, not trace.
- **β movement is trace-dominated** (share_T 0.57); α/low-γ movement is
  reset-dominated (0.42/0.38) — when the β tree moves it tends to *persist*; when
  the α/low-γ tree moves it tends to *revert*.
- **Strength diagnostic:** anchor is only weakly–moderately strength-correlated
  (cohort r 0.05–0.38), so anchor is not a pure hubness tautology, but strength
  contributes (strongest at low-γ).

**Caveat (σ fragility).** The anchor/mover split uses a per-patient RMS noise
scale σ (median |Δbase| is 0 under quantization — switched to RMS). σ is
heterogeneous (β mover% ranges 0.7%→99.7% across patients), so the *absolute*
anchor mass is patient-fragile. The *band ordering* (β trace peak) and the
*among-mover direction* (share_T) are robust; the localization (§4, no σ
threshold) is robust by construction.

## 3. Coupled-cross-phase null (audit_106) — FAILS its gate (honest negative)

The coupled construction shares one 4-cycle swap sequence across phases
(per-phase feasible δ). **Validation gate (run first): FAIL.** At the canonical
20·E swaps the coupled surrogate preserves only ~31% of the observed cross-phase
backbone similarity (`pre_post` coupled/observed ratio 0.31; coupled–independent
margin +0.03 — negligible). Full matched-strength randomization decorrelates the
phases regardless of shared swap locations. This **confirms audit_63's statement**
that a coordinated cross-phase surrogate "isn't well-defined in the
strength-preserving family." Consequence: anchor/reset are tested with the
**geometry baseline** only (the chosen spatial null); they get **no temporal
significance**, and are reported **unverified**. (Sanity preserved: independent-null
ρ_split median ≈ 0, matching audit_63.)

## 4. Localization (audit_107, geometry baseline) — the dissociation

Per-pair channel signal → demeaned per-system endpoint mean (audit_83 engine) →
count-matched node-label-permutation null (R=500) → one-sided enrichment p,
BH-FDR over systems within (band, signal). **Sanity: trace → OFC reproduced**
(β, q=0.040), validating the pipeline.

β systems (epi-include, non-shaft), top by p (q = BH within band×signal):

| channel | β top system(s) | reading |
|---|---|---|
| trace (concordance) | **OFC** q=0.040, occipital q=0.040 | reproduces locked β→OFC |
| persist (φ₁²) | cingulate q=0.15, **MTL** q=0.15 | the MTL "hint" (limbic) |
| reset (φ₂²) | **lateral_temporal** q=0.020 | reversible component |
| reset−trace (φ₂²−φ₁²) | **lateral_temporal** q=0.010, **sensorimotor** q=0.010 | reset-specific |
| anchor (−var) | **insula** q=0.020, sensorimotor q=0.030 | rigid backbone |

**The dissociation:** at β the persistent trace (OFC), the reversible reset
(lateral temporal / sensorimotor), and the rigid anchor (insula) live in
**anatomically distinct systems**.

**Cross-band anchor regularity:** anchor localizes to **limbic** systems in every
trace band — insula (β, q=0.020), **MTL (α q=0.040, low-γ q=0.040)**. The rigid
backbone is anatomically limbic — the "anchor = anatomical component" hypothesis,
consistent with the retired-KC anchor-anatomy 6.2× (different method, same theme).
At α/low-γ the trace and reset channels are mostly non-significant (consistent with
β being THE trace band).

### Robustness — shaft-collapse (decisive for any *distributed* claim)

The β dissociation **survives shaft-collapse** (each electrode shaft collapsed to
one observation before aggregating — the audit_83 control that killed the locked
"distributed paralimbic ring"). Non-shaft → shaft-collapsed (R=200) best-system q:

| β channel | non-shaft | shaft-collapsed | verdict |
|---|---|---|---|
| trace → OFC | q=0.040 | q=0.100 (OFC stays top) | reproduces locked OFC (locked uses R=1000 MS; here R=200 geometry) |
| reset → lateral_temporal + sensorimotor | q=0.020 | **q=0.050** | **survives** |
| reset−trace → lateral_temporal + sensorimotor | q=0.010 | **q=0.050** | **survives** |
| anchor → insula | q=0.020 | **q=0.050** | **survives** |

So reset (lateral temporal / sensorimotor) and anchor (insula) are **not
single-shaft artifacts** — multi-shaft anatomical concentration. They land
**right at q≈0.05** (shaft-collapse cuts the effective sample to distinct shafts
and R=200 coarsens the floor), so the dissociation is robust but **not
overwhelming**; R=1000 + LOO-patient would firm it up before any headline.

## 5. Honest verdict — UPDATED 2026-06-18 by the matched-strength test (audit_110)

**The dissociation did NOT survive the mandatory matched-strength null. Only
trace → OFC stands.** audit_107 localized reset/anchor with the GEOMETRY baseline
only (controls contact count, not node strength). audit_110 reran the localization
of every channel against the STRENGTH-PRESERVING surrogate (the locked null behind
β→OFC, R=200 cached; R=1000 available). Result:

| β channel | geometry baseline (audit_107) | **matched-strength (audit_110)** | verdict |
|---|---|---|---|
| **trace → OFC** | p=0.008, q=0.040 | **p=0.005, q=0.025, LOO≤0.005, str−0.51** | **SURVIVES** (low-strength, LOO-robust) |
| reset → lateral_temporal | p=0.002, q=0.020 | **p=0.104, q=0.55** | **FAILS** — strength artifact |
| reset_dom → lat_temp / sensorimotor | p=0.002–0.010 | p=0.10 / 0.18, q≈0.50 | FAILS |
| anchor → insula | p=0.002, q=0.020 | **p=0.31, q=0.63** | **FAILS** — strength artifact (as predicted) |
| persist → MTL/cingulate | p=0.02–0.03 | p≈0.13, q=0.66 | FAILS |

- **Verified (and only this):** **trace → OFC** = the locked β→OFC result.
  Reproduces under matched-strength (p=0.005, **LOO-robust**, strength_dev −0.51 =
  OFC contacts are *below* average strength, so not a hub artifact). The
  decomposition is exact (reproduces the locked trace 60/60). This is the
  **existing** result, not new.
- **FALSIFIED — reset / anchor dissociation.** The geometry-baseline significance
  was **node-strength**: lateral-temporal's reset excursion and insula's anchor
  rigidity are reproduced by strength-matched surrogate graphs. This is exactly the
  KC-style failure mode the matched-strength rule exists to catch (the geometry
  baseline didn't reach the strength alternative). The trace/reset *anatomical
  dissociation* is **retracted** — not publishable. (anchor failing was expected:
  −var ≈ rigid backbone ≈ strength.)
- **Null by construction:** reorganize = the noise floor (high-γ exemplar).
- **What this leaves:** the per-pair flow decomposition (anchor/trace/reset) is a
  clean *framing* of the existing OFC result; the **only strength-independent,
  anatomically-specific node-flow signature is the persistent (trace) component →
  OFC**. The reversible (reset) and rigid (anchor) components are **strength-
  organized, not region-specific**. That is the honest node-flow characterization
  (see §6.5) — a nuanced/partly-negative result, not the 3-way dissociation.
- **Cognitive reading** (task = transitive inference; OFC = cognitive map): the
  *persistent* component consolidates in OFC. The earlier "reversible reset sits in
  lateral-temporal / rigid anchor in limbic" reading is **withdrawn** (strength
  artifact). Only "the cognitive map (persistent flow) is retained in OFC" stands.

## 6. Brain figure (audit_109) — the dissociation, on the brain

The catchy anatomical figure the taxonomy was built for (prior 3D / trajectory
attempts failed — phase-states are near-equidistant, so the trace is a subtle
rank-level pull, not a geometric collapse; see `cross_phase_taxonomy_2026_06_12`).
Pools all 10 patients' contacts into approximate MNI space (the per-patient
DK-anchored affine in `visuals.spatial_coords` — adequate for a cohort glass
brain, not mm-precise) and lights up each contact by its anatomical SYSTEM's
**incident-edge enrichment** — how strongly the channel concentrates on the
EDGES incident to that system (the audit_107 edge→endpoint incidence,
geometry-baseline `−log10 p`); white-matter / non-enriched contacts are faint
grey context; systems passing **BH q<0.05** get a black ring. **The taxonomy is
per-PAIR, never per node:** a contact's colour is its system's incidence
enrichment, *not* a node classification — a system/node has three independent
incidence scores (trace / reset / anchor), not one label, which is why the
figure is three separate maps. This is the same edge→endpoint step behind the
locked β→OFC verdict, so "trace → OFC" means *edges incident to OFC carry the
trace* (relational, "concentration not container"), not "OFC nodes are trace
nodes". It is a *result-at-system-granularity* map (the granularity audit_107
verified), **not** a per-contact field — single-contact cophenetic values are
too noisy to read (audit_76/108).

- `figures/brain_dissociation_beta.pdf` — 3 rows (trace / reset / anchor), β.
  The rows glow in three distinct territories. **⚠ SUPERSEDED for the reset/anchor
  rows** — that dissociation FAILED matched-strength (§5, audit_110); the figure
  shows the *geometry-baseline* glow, which is strength-confounded for reset/anchor.
  Only the **trace → OFC** row is a verified result. Do not present the 3-row
  dissociation as a finding; if used, keep the trace row only (or relabel
  reset/anchor "strength-baseline, not strength-controlled").
- `figures/brain_trace_bands.pdf` — trace channel, α / β / low-γ. Only **β**
  rings/glows (OFC + occipital); α and low-γ are empty grey context →
  **β-specificity** of the localized trace.

**Read with the locked caveats.** trace → **OFC** is the verified/locked result
(K=5 implanted, survives shaft-collapse, matched-strength R=1000, BH q≈0.01);
**occipital glows equally here but is K=3** and geometry-baseline-only — a
secondary hint, not co-equal with OFC (per `localization_audit_plan_2026_05_29`).
reset (lateral-temporal) and anchor (insula/sensorimotor) are **new and
UNVERIFIED** — geometry baseline only (the legend says so), no temporal null
(audit_106 coupled null failed its gate). The figure is the cleanest statement
of the *exploratory* dissociation, not a verification.

## 7. Next (if the user wants to proceed)

- shaft-collapse + LOO-patient robustness on the reset dissociation before any
  claim beyond "hint".
- decide preprint placement (candidate §5.x decomposition + the brain figure)
  — exploration first.
- optional: a genuine coordinated null via the shared-backbone-deviation /
  SB-CRC route (`coherency_surrogate.py`) rather than coupled rewiring.
