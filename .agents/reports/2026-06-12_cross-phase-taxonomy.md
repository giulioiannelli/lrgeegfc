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
  - scripts/01_compute/audit/audit_107_taxonomy_localization.py         # step 3
  - data/audit/cross_phase_taxonomy/                                    # all outputs
  - localization_audit_plan_2026_05_29.md                              # locked beta->OFC trace
---

# Cross-phase taxonomy — anchor / trace / reset / reorganize (results)

**Head.** The per-pair cophenetic signal that carries the validated β trace
decomposes — exactly, with no residual — into a rigid **anchor** backbone plus
two orthogonal fluctuation contrasts, **trace** (persistence φ₁, = ρ_split) and
**reset** (excursion φ₂). Three findings: (1) the decomposition **reproduces the
locked trace bit-exact in all 60 cells** and the trace energy **peaks at β**, with
anchor dominating every band (the Gratton "stable backbone" picture) except the
high-γ noise band; (2) at β the channels **localize to distinct systems** — trace
→ **OFC** (reproduced), reset → **lateral temporal**, anchor → **insula** — and
**anchor localizes to limbic systems across bands** (insula β, MTL α/low-γ), the
"anatomical backbone" hypothesis; (3) the attempted **coupled-cross-phase
matched-strength null fails its validation gate**, confirming audit_63's claim
that a coordinated strength-preserving cross-phase null is not viable, so the
similarity channels rest on the geometry baseline. **Only the trace channel is
verified**; reset and anchor are new, **unverified** (geometry baseline controls
contact count, not strength or a temporal null).

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

## 5. Honest verdict

- **Verified:** the trace channel = the locked β→OFC result; the decomposition is
  exact and reproduces it 60/60.
- **New + unverified (geometry baseline only, no temporal null):**
  - **reset → lateral temporal + sensorimotor at β** — a genuine trace/reset
    *anatomical dissociation*: distinct from the OFC trace, and it **survives
    shaft-collapse** (q≈0.05). Reset is otherwise modest in magnitude (cohort
    reset energy ≤0.11 except the high-γ noise band) and α/low-γ reset is
    non-significant. The single most promising *new* result; firm up with
    R=1000 + LOO before any headline.
  - **anchor → limbic across bands** — insula (β, survives shaft-collapse
    q≈0.05), MTL (α/low-γ, q=0.04). The "anatomical backbone" hypothesis holds,
    but partly a strength/backbone story (geometry baseline does not control
    strength; cohort r(anchor,strength) 0.05–0.38).
- **Null by construction:** reorganize = the noise floor (high-γ is the exemplar:
  anchor 0.19, everything moves, no direction).
- **Framing:** RPre→Task→RPost as a hysteresis loop — trace = hysteretic
  (non-returning) branch, reset = reversible branch, anchor = rigid, reorganize =
  noise. Literature: flexibility/allegiance (Mattar/Braun) with the strength null
  those lack; Gratton stable-backbone (anchor dominance, confirmed); Tambini
  post-task persistence (trace); Kim/Lee/Mashour hysteresis (trace vs reset).

- **Cognitive reading (task = transitive inference; OFC = cognitive map, see
  `task_paradigm_transitive_inference`).** The dissociation lines up cleanly: the
  **persistent trace sits in OFC** — the learned relational hierarchy / cognitive
  map is what *consolidates* offline into rest_post; the **reversible reset sits in
  lateral-temporal / sensorimotor** — the perceptual/motor machinery that runs the
  task and then reverts; the **rigid anchor is limbic** (insula/MTL) — the
  state-invariant backbone. So the taxonomy reads as: the cognitive map is retained
  (trace) while the task's perceptual scaffolding is released (reset). This is a
  *post hoc* interpretation, not a tested claim.

## 6. Next (if the user wants to proceed)

- shaft-collapse + LOO-patient robustness on the reset dissociation before any
  claim beyond "hint".
- decide preprint placement (candidate §5.x decomposition) — exploration first.
- optional: a genuine coordinated null via the shared-backbone-deviation /
  SB-CRC route (`coherency_surrogate.py`) rather than coupled rewiring.
