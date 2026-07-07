---
name: headline-n3-epileptogenic-markers
era: IMCOH_ABS_COHORT_N10
status: current
kind: headline
scope: N3 — the same Laplacian propagator that carries the cognitive trace also localizes epileptogenic tissue. Relational δ diffusion community → seed-based 6-band calibrated P(SOZ) detector. Seed-free markers are a documented NEGATIVE (5 confirmations), folded OUT of the Results (2026-07-07). Detailed audit history in epi-marker-analysis.md.
owner_agent: epi-marker-trace-analysis chat
updated: 2026-07-07
---

# N3 — The same diffusion propagator localizes epileptogenic tissue: interpretable seed-based SOZ markers

> One operator, a second read-out (CORE). Read ρ̂(τ)=e^{−τL}/Z as inter-contact
> diffusion *affinity* instead of a cross-phase distance, and the seizure network
> shows up. Full audit-by-audit history → [`epi-marker-analysis.md`](epi-marker-analysis.md).
> Numbers in cached CSVs (§F).

> **Update 2026-07-07 (paper reconciliation).** The **seed-free** thread (former **N3.3**) is a
> **documented negative** — five independent tests (`audit_80/81`, `audit_94`, `audit_118`,
> `audit_120`, `audit_121`) all sit at the strength or label-shuffle floor; only electrode depth
> (an implantation-placement prior) transfers. It is **folded OUT of the Results**. Paper
> paragraphs renumbered **R3.3 = two-populations, R3.4 = scope**; headline **N3.4→N3.3, N3.5→N3.4**.

## §A — Result in plain language

The same lens that finds the cognitive trace also finds the **seizure-onset zone
(SOZ)**. For a neuroscience (not engineering) audience, **lead with the finding,
not the detector**: epileptogenic tissue is **not identifiable node-by-node** — it
is organized **relationally**, as a strength-independent, spatially-irreducible
**diffusion community**. That is the brain-organization claim; the calibrated
detector below is its *application*. The honest arc, in one breath:

- A **node-by-node** "is this contact epileptic?" marker **does not exist** beyond
  trivial confounds (within a patient it reduces to hubness + electrode depth;
  across patients nothing transfers).
- But read **relationally**, the epileptic contacts form a **strength-independent,
  spatially-irreducible diffusion community** — a group that *co-diffuses* — and
  that community is a genuine structural finding that survives strength- and
  space-matched nulls.
- Turned into a **marker**: given a few known SOZ contacts as seeds, the diffusion
  picture points at the others, **including contacts far from the seeds, on
  different electrodes** (spatial proximity is removed by construction, because
  the |ImCoh| substrate is volume-conduction immune).
- Pushed to a **deployable detector** (current state): fuse the propagator
  features across **all six bands** and train **leave-one-patient-out** → a
  **calibrated probability P(SOZ)** for every contact. This is a **triage /
  shortlist** tool (good ranking, modest precision), not a from-scratch clinical
  localizer.
- Seed-free (no known SOZ) is a **documented negative** on this data — five independent
  tests; only implantation depth (a placement prior) transfers, not connectivity. The
  deliverable is the **seed-based** detector above.

**Why it matters:** a *mechanistic, multiscale, strength-orthogonal* SOZ read-out
from the same framework as the cognitive result — cognition and clinic from one
operator. New to the SOZ-marker field (which is dominated by black-box or
single-feature methods).

## §B — Technical statement (per subheadline)

**N3.1 — Relational δ diffusion community marks distant SOZ.** A per-node read-out
of how strongly a contact co-diffuses with known SOZ recovers **distant** SOZ
(off-shaft, proximity removed), **δ-band**, stable across rest and task, with a
label-shuffle null at p≈0. The all-contacts version was a proximity tautology and
was **withdrawn**; the honest metric is leave-one-shaft-out off-shaft discovery.

**N3.2 — Seed-based 6-band compound detector (canonical, current).** Propagator-
only features (heat-affinity at multiple scales, PPR, Katz, communicability,
commute/resistance, diffusion-distance, seed-relative segregation), **no
coordinates**, per band, concatenated across 6 bands, L2-logistic, k=3 seeds,
LOPO → calibrated P(SOZ). Headline: **mean LOPO AUC ≈ 0.81, precision@5 ≈ 60%
(~7× chance), 9/10 patients above chance**, clean label-shuffle null. The lift is
**cross-band fusion**, not richer single-band features. A gradient-boosted variant
reaches higher AUC but **no precision gain** and overfit risk on n=10 → **rejected
in favour of the interpretable logistic**. (⚠ A first draft scored 0.97 via a
**label leak** — segregation referenced the full SOZ set, not the seeds; caught,
fixed, shuffle-null certifies the clean pipeline.)

**N3.3 — The two-population split.** Detector performance is bimodal: a **strong-
community** group (several patients at AUC ≈ 0.9–0.99, precision@5 up to 100%) and
a **hub** group (Pat_10, Pat_15 — right-hemisphere implants where SOZ = hubs, not
a community) where it is marginal. Report per-patient, not just the cohort mean —
the cohort median can lie.

**N3.4 — Honest scope.** Validated against **clinical SOZ labels, not resection +
surgical outcome (Engel/ILAE)**. **Seed-based** (N3.2; cannot localize from zero —
seed-free is a closed negative). Precision is a **triage** number, **ceiling-bound** (~60%) by target
rarity × ranker quality — no filter rescues it (white-matter filter even hurts:
~41% of SOZ are WM-labelled). Occult candidates (unmarked contacts ranked high —
e.g. hippocampus/amygdala, medial-OFC) are **hypotheses**, none outcome-validated.

## §C — Critical issues & powerful strengths

**Strengths:** mechanistic + interpretable (logistic on named diffusion features);
strength-orthogonal and proximity-immune by construction; cross-band fusion is a
real, defensible lift; calibrated probabilities; clean nulls after the leak fix;
finds *distant* SOZ, which is the clinically hard case.

**Lead-with weaknesses:**
- **Clinical-label not outcome** — the honest ceiling on every claim.
- **Seed-based** is the current canonical; the clinically valuable **seed-free**
  case is unproven (N3.3 live).
- **Two-population split** — hub patients (Pat_10/15) are genuine failures, not
  noise.
- **Precision ceiling ~60%** is structural (rare target); don't imply a
  standalone localizer.
- **n=10** — small; the GBM-vs-logistic choice is partly an overfitting-defense.

## §D — To-dos & verifiables (owner: epi-marker-trace-analysis chat)

- [ ] **Seed-free markers** (`audit_118/119/120`): characterize honestly — how
  much precision is lost vs seed-based; is cross-patient selection still at chance?
- [ ] **Cross-phase routing rigidity** (`audit_121`): is epileptic routing more
  *rigid* across phases than healthy routing? (potential new sub-result.)
- [ ] Lock the **canonical detector config** description against
  `data/audit/epi_propagator_detector/README.md` (keep the leak-fix + shuffle-null
  story; they are the credibility).
- [ ] Decide framing: **N3 as a "method generality" headline** (one operator →
  clinic) vs a standalone clinical paper. Keep it folded here for now.
- [ ] Outcome-label feasibility (Engel/ILAE) — the upgrade from triage to claim.

## §E — Figure / representation ideas

- Diffusion-community schematic: seeds → co-diffusion affinity → distant SOZ lit
  up on a different shaft (the "proximity removed" point).
- Per-patient detector strip (the two-population split) — **not** a single cohort
  number.
- Calibration curve for P(SOZ) + precision@k.
- Occult-candidate brain (clearly labelled HYPOTHESES).
- (When ready) seed-based vs seed-free side-by-side precision.

## §F — Provenance (CSV · script · timestamp)

- Detector (canonical) → `data/audit/epi_propagator_detector/{detector_lopo_per_patient.csv,
  detector_node_predictions.csv,occult_candidates_graymatter.csv,detector_ablation.csv,
  feature_importance.csv}` · `audit_117_epi_propagator_detector.py` · 2026-06-22;
  README `data/audit/epi_propagator_detector/README.md`.
- Relational community + distant marker (history) → see
  [`epi-marker-analysis.md`](epi-marker-analysis.md) §§ + reports
  `.agents/reports/2026-06-08_epi-propagator-soz-marker.md`,
  `2026-06-12_propagator-distant-soz-marker.md`,
  `2026-06-16_epi-distant-marker-walkthrough.md`; memory
  `epi_propagator_diffusion_community_2026_06_08`.
- **Seed-free (LIVE)** → `audit_118_epi_seedfree_marker.py`,
  `audit_119_epi_path_routing.py`, `audit_120_epi_seedfree_routing_module.py`,
  `audit_121_epi_crossphase_routing_rigidity.py` · 2026-06-22 (running).
- Literature positioning → `.agents/reports/2026-06-11_epi-soz-marker-literature-positioning.md`.

## §G — Missing parts / open

- **Seed-free verdict** (N3.3) — gates whether N3 can claim zero-label
  localization.
- **Surgical-outcome validation** (Engel/ILAE) — the upgrade path.
- A clean statement of **why δ** (and which bands the 6-band fusion actually
  leans on — `feature_importance.csv`).
