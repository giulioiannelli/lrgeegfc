---
name: propagator-distant-soz-marker
era: IMCOH_ABS_COHORT_N10
date: 2026-06-12
status: current
kind: report
scope: VERIFIED — the LRG heat-kernel propagator (slow τ, δ band) marks DISTANT epileptic nodes (off-shaft, proximity discarded, strength-controlled); label-shuffle null p=0 + nested-LOPO held-out 0.72. Supersedes the proximity-confounded audit_92/98 framing.
supersedes: the all-contacts marker numbers in .agents/reports/2026-06-08_epi-propagator-soz-marker.md §3-4 (those were proximity-confounded)
---

# The propagator marks distant epileptic nodes (δ) — verified

**Head.** Using the LRG framework directly: the heat-kernel propagator **ρ(τ)=e^{−τL}/Z
at slow diffusion time, in the δ band**, given SOZ on some electrodes, **ranks SOZ on
*other* electrodes above healthy contacts** — the non-trivial *distant* discovery, where
spatial proximity is useless by construction. Median leave-one-shaft-out AUC **0.72,
8/10 patients**, strength-residualised. It is verified two ways: a **label-shuffle null
gives p=0.000** (random same-size contact sets collapse to chance), and **nested
leave-one-patient-out holds at 0.72** (the marker is selected on 9 patients, scored on the
10th). The heat kernel — *the* propagator — beat communicability, PageRank, Katz, and
effective-resistance in a 16-operator sweep. This is a **small, true** claim that replaces
the **big, proximity-confounded** numbers from 2026-06-08.

---

## 1. Why the metric is leave-one-SHAFT-out (proximity is a rigged baseline)

SOZ contacts are labelled because they sit in a clinically-marked epileptogenic **area**,
so "a contact near a known SOZ is itself SOZ" is **baked into the labelling**, not a
discovery. A trivial "flag the nearest contacts" baseline therefore beats any network
marker ~2× on all-contacts recall (audit_99) — but that is the labelling tautology, not
performance. The only honest, non-trivial task is **distant** discovery:

> **Leave-one-shaft-out (LOSO):** hide *all* SOZ on one electrode; using SOZ on the *other*
> electrodes as seeds, score every contact by propagator affinity to the seeds, and ask —
> among contacts on electrodes that contain **no seed** — do the hidden shaft's SOZ rank
> above healthy? Proximity cannot help: targets and controls are all far from the seeds.

Headline metric = **strength-residualised** LOSO AUC (so hubness is also controlled).
Proximity is **discarded entirely**, not used as a baseline.

## 2. The result (δ band, heat kernel e^{−τL}, slow τ)

Per-patient strength-residual LOSO distant-discovery AUC (`marker_library_per_patient.csv`):

| patient | heat_t5 | note |
|---|---|---|
| Pat_08 | 0.95 | responder |
| Pat_14 | 0.91 | responder |
| Pat_05 | 0.84 | responder |
| Pat_06 | 0.82 | responder |
| Pat_13 | 0.76 | responder |
| Pat_03 | 0.68 | |
| Pat_02 | 0.66 | |
| Pat_07 | 0.57 | marginal |
| Pat_10 | 0.41 | ✗ hub-patient |
| Pat_15 | 0.21 | ✗ hub-patient |

**Median 0.72, mean 0.68, 8/10 > 0.5.** The 2 failures are the hub-patients (SOZ are hubs,
not a co-diffusing community) — the same two-population split seen across the whole
investigation. heat_t3/t4/t5 all give median 0.72–0.75 (robust across the slow-τ range);
normalized-Laplacian heat (`heatN_t4`) matches (0.74, 8/10). Slow diffusion **reaching
far** is the mechanistically pre-specified operator for distant discovery — it is not
cherry-picked, and it beat every non-heat operator.

## 3. Verification (audit_102)

- **(A) Label-shuffle null — p = 0.000** (500 draws, each hero marker). Relabel n_epi
  random contacts (spanning ≥2 shafts) as fake-SOZ and rerun LOSO: null median ≈ **0.487**
  (p95 ≈ 0.56) vs real **0.72–0.75**. Random same-size sets are NOT distant-recoverable →
  the recovery is specific to the **true SOZ being a co-diffusing community** (the
  node-level expression of the audit_89 C3/C5 finding). Kills "graph artifact" / "any set
  works".
- **(B) Nested LOPO — held-out median 0.716** (mean 0.671, 7/10 > 0.5). Best marker
  selected on 9 patients, scored on the held-out 10th; `heat_t5` is selected in **9/10**
  folds. Kills the 64-cell multiple-comparison concern: the same marker wins out of sample.

## 4. The 16-operator sweep (audit_101) — heat kernel wins

All operators from L = D − W of the |ImCoh| graph, seed-affinity = mean over seeds of the
operator column. Families: heat kernel e^{−τL} (6 τ), normalized-Laplacian heat (6 τ),
personalized PageRank (α∈{0.5,0.85,0.95}), Katz, communicability expm(W/ρ), negative
effective-resistance (commute distance via L⁺), diffusion-distance; + strength baseline
and seed-free intrinsics (heat-diagonal, slow-mode participation, Fiedler). **In δ the top
6 cells are all heat/heatN variants** (median 0.72–0.75); communicability/PPR/Katz/
resistance and all intrinsics sat below. β/low-γ/α carry no distant-discovery signal
(median ~0.50–0.60, sign-test n.s.) — **the marker is δ-specific.**

## 5. Honest scope
- **δ band only.** Other bands do not carry distant discovery.
- **Seed-based.** Needs known SOZ on ≥2 electrodes to find SOZ on others; not from-scratch,
  cannot run a zero-label patient.
- **8/10 patients;** fails on the 2 hub-patients (Pat_10, Pat_15).
- **Clinical-label validated, not surgical outcome** (no resection/Engel ground truth).

## 6. Supersession — what this corrects
The 2026-06-08 marker numbers (audit_92 masking AUC ~0.81, audit_98 deployment, "beats
strength") were computed on an **all-contacts** metric and are **proximity-confounded** —
audit_99 showed a trivial nearest-contact baseline beats them ~2×. Those numbers are
**withdrawn as a performance claim**. The δ slow-heat-kernel LOSO result here is the clean
replacement: proximity removed (off-shaft), strength removed (residual), null- and
selection-verified. **Lesson (persist):** for any seed-based epileptic-node discovery,
evaluate **leave-one-shaft-out / off-shaft only** — the all-contacts metric is a labelling
tautology.

## 7. Deployment — honest precision, the compound, P(SOZ) (added 2026-06-16; audit_103/104)

Turning the verified ranker into a tool, evaluated OFF-SHAFT throughout:

- **Precision is real but concentrated.** Off-shaft SOZ prevalence ≈ 4%; cohort-median
  **lift@5 ≈ 2.2× (responders), recall@10 ≈ 0.34**. Strongly **bimodal**: Pat_08/14/05 get
  9–15× lift, several patients ≈0 at the top-5 even with decent AUC (precision@k with tiny
  targets-per-shaft is brutal). A good ranker, a concentrated shortlist — not a clean catch.
- **Compound: switch, don't blend.** Propagator and strength are complementary (propagator =
  responders, strength = hub-patients). A **hard seed-regime SWITCH** — use strength when the
  seeds are hub-like (mean seed strength-percentile r ≥ 0.70), else the pure propagator —
  lifts AUC **0.716 → 0.750, 9/10**, *preserving* the strong responders (Pat_08 stays 0.95)
  and rescuing the one recoverable hub-patient (Pat_15 0.21 → 0.93). A continuous **blend is
  misleading** — it raises cohort *lift* but **damages the best responders** (Pat_08 0.95 →
  0.56); the cohort median hid it. Regime detector imperfect: Pat_10 unrecoverable (both
  signals fail), borderline r ≈ 0.72 patients ambiguous.
- **P(SOZ) is emittable, calibrated-at-top, modest.** A leave-one-patient-out logistic on
  [affinity, strength] gives a probability that tracks the diagonal at the high end (P≈0.13 →
  13% observed) and transfers for responders (held-out AUC median 0.79) but not hubs; absolute
  values are small (top candidates ~0.1–0.35) and Brier ≈ the 4%-prevalence floor. It is a
  triage score, not a per-node verdict.
- **Effective-marker path:** adopt the switch; filter white-matter contacts; more seeds across
  more shafts; the real ceiling is the two-population structure + the absence of
  surgical-outcome ground truth (precision here is an upper bound on true occult discovery).
- Figures: `preprint_39` (core marker) + `preprint_40` (deployment). Audits: `audit_103`
  (precision + candidate shortlists), `audit_104` (compound + probability). READMEs in
  `data/audit/epi_marker_{precision,compound}/`. Walkthrough:
  `.agents/reports/2026-06-16_epi-distant-marker-walkthrough.md`.

## Provenance (numbers in CSVs; scripts + 2026-06-12)
- 16-operator sweep: `scripts/01_compute/audit/audit_101_epi_marker_library.py` →
  `data/audit/epi_marker_library/marker_library_{per_patient,cohort}.csv`.
- Verification (null + LOPO): `scripts/01_compute/audit/audit_102_epi_marker_verify.py` →
  `data/audit/epi_marker_library/verify_{nulls,lopo}.csv`.
- All-contacts vs off-shaft discovery (why proximity is discarded):
  `scripts/01_compute/audit/audit_99_epi_discovery_yield.py` →
  `data/audit/epi_discovery_yield/`.

## Next (optional)
- Combine the δ propagator across τ-scales (area-under-diffusion-curve) — does it lift 8/10
  toward 10/10, or just refit the 2 hub-patients (it should NOT rescue them)?
- Figure: per-patient distant-discovery AUC + the label-shuffle null band, δ heat kernel.
- The mechanism link: distant-discovery responders = the audit_89 C3/C5 community
  responders (same ~5–8 patients); confirm the overlap explicitly.
