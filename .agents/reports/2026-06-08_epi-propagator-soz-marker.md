---
name: epi-propagator-soz-marker
era: IMCOH_ABS_COHORT_N10
date: 2026-06-08
status: current
kind: report
scope: Part A (epi diffusion community, C1–C5 all closed incl. new C5 spatial null) + Part B (interpretable per-node SOZ marker — relational features, masking recovery, cross-patient classifier, occult discovery)
supersedes_forward_plan: .agents/preprint/directives/investigation_directive_2026-06-08_epi-propagator-soz-marker.md
---

# Epileptic diffusion community + interpretable SOZ marker — results

**Head.** Read relationally, the LRG propagator ρ(τ)=e^{−τL}/Z shows epileptic
contacts forming a **strength-independent, spatially-irreducible, band-specific
diffusion community** — and that community can be turned into a **per-node SOZ
marker that recovers held-out SOZ within patient and transfers across patients**.
Part A is now airtight: the epi block survives the full C1–C5 battery, including the
new **C5 spatial-matched null** (δ/β/low-γ + α survive; high-γ a clean null). Part B
delivers the marker: strength-orthogonal relational features separate epi from non
at AUC 0.74–0.86; a masking experiment recovers hidden SOZ at AUC 0.72–0.82
(beating the hubness baseline); an **interpretable leave-one-patient-out logistic**
predicts SOZ cross-patient at AUC 0.68–0.81 (permutation-p = 0) where node strength
is at chance; and the calibrated P(SOZ) ranks **occult (unlabelled) candidates** with
plausible anatomy (e.g. Pat_06 medial-OFC, P=0.87, 3.5 mm from a labelled SOZ).
**Honest scope:** the marker is *seed-based* — it needs a few labelled SOZ per
patient to define the community; it is not a from-scratch detector, and occult
candidates are hypotheses (no resection/outcome ground truth in this cohort).

Numbers live in CSVs (no tables here), each with its generating script + a
2026-06-08 timestamp. Pointers in §Provenance.

---

## TOC (plain-English)
1. **Part A is closed** — the epi diffusion community survives strength *and* space.
2. **The relational features are the marker** — affinity to the SOZ community, not
   node strength; strength-orthogonal by matched-strength-z.
3. **It recovers hidden SOZ within patient** — mask half the SOZ, find them back.
4. **It transfers across patients** — an interpretable logistic, calibrated P(SOZ),
   where hubness fails cross-patient.
5. **It proposes occult SOZ candidates** — ranked, with proximity + region + reasons.
6. **What it is not** — seed-based, not from-scratch; candidates are hypotheses.
7. **Outliers, bands, and the honest caveats.**

---

## 1. Part A — the diffusion community survives C1–C5

The measure (audit_89): partition node *pairs* by epileptic membership and read the
propagator as block means ρ_EE / ρ_EN / ρ_NN. The question is whether ρ_EE is
elevated **beyond strength and beyond anatomy**. Control battery (cite by number):

- **C1** matched-strength absolute (ρ_EE vs strength-preserving rewiring) — passes
  but **inflated** (rewiring destroys all community structure); context only.
- **C2** matched-strength contrast (ρ_EE − ρ_NN) — epi block enriched more than the
  healthy block. Strong δ/α/β, moderate θ/low-γ, null high-γ.
- **C3** strength-matched random-subset null on the real graph — epi vs arbitrary
  equal-size, equal-strength groups. **Decisive — cannot be hubness.** Strong
  δ/low-γ/β, moderate α/θ, null high-γ.
- **C4** cross-shaft restriction — positive where tested; statistical control only,
  not a figure (the full epi–epi block stays in the visual).
- **C5 (NEW, audit_90) — strength+spatial-matched random-subset null.** Random
  groups matched to the epi set on *both* the strength-quintile profile (exact, by
  construction) *and* spatial spread (radius of gyration in implant x,y,z, matched
  to ~1e-6 by within-quintile swaps). **This closes the anatomy confound.**

**C5 verdict** (cohort n=9, max over τ; numbers in `spatial_null_cohort.csv`): every
community band survives — δ +3.05 (p.020), low-γ +3.38 (p.002), α +1.81 (p.020),
β +1.59 (p.006), θ +2.37 (p.020); **high-γ null** (+1.03, p.455, 1/9). The
directive's gate ("all of δ/β/low-γ survive → full claim") is met. The epi diffusion
community is **neither hubness (C3) nor spatial clustering (C5)**.

Transversal robustness: positive across all 6 τ; high-γ a clean internal negative
at every control. Literature framing unchanged — "epileptogenic zone internally
hyperconnected" is known; the new part is **strength-independent + spatially-
irreducible + multiscale + band-resolved**.

## 2. The marker is relational, not node strength (audit_91)

The whole investigation's lesson (audit_85/86/87 node strength & recovery; audit_88
signed Laplacian): **node-intrinsic / strength features fail beyond hubness.** The
win is relational. Per node i and the known epi set E, three interpretable,
strength-orthogonal (matched-strength-z) features from ρ(τ):

- `f_aff` mean ρ(i, E∖i) — affinity to the SOZ community.
- `f_seg` ρ̄(i,E) − ρ̄(i,N) — interface segregation (the V3 signal: audit_89 ρ_EN is
  *suppressed* beyond strength, so SOZ couple to SOZ more than to healthy tissue).
- `f_part` share of i's diffusion routed to SOZ.

Within-patient leave-self-out epi-vs-non separation (τ-mean MS-z; numbers in
`relational_separation.csv`): δ 0.863, α 0.827, β 0.818, low-γ 0.744, θ 0.707 — all
p<0.05, far above the strength baseline (~0.60–0.74). **high-γ 0.318 (p.94) — null.**
The interface feature `f_seg` separates as strongly as `f_aff` (V3 settled).

## 3. Within-patient masking recovery (audit_92, VII1)

The non-tautological predictive test: hide 50% of each patient's labelled SOZ,
learn only "affinity to the *visible* SOZ", recover the hidden SOZ among healthy
contacts (200 splits) + leave-one-epi-out. Held-out SOZ never inform their own
feature. Scores: `resid` (affinity residualised on strength — strength-orthogonal)
vs `strength` (hubness). Numbers in `masked_recovery_cohort.csv`:

- `resid` recovery AUC: δ 0.823, low-γ 0.795, θ 0.748, β 0.721, α 0.616 — all
  p<0.05. **high-γ 0.450 (p.85) — null.**
- `resid` beats `strength` in every non-null band (Δ +0.06 → +0.21), decisively at
  δ (paired p.020), directionally elsewhere (n=9 paired Wilcoxon underpowered).

The dissociation is clean per-patient: Pat_02's hidden SOZ are recovered by `resid`
(0.80) not strength (0.40, below chance); Pat_10's by strength (0.80) not `resid`
(0.63) — Pat_10's epi are hubs, not a community (the C3/C5 non-responder).

## 4. Cross-patient classifier + calibrated P(SOZ) (audit_93, VI2/VI3)

An **interpretable logistic** on `z_f_aff_mean` + `z_f_seg_mean`, leave-one-patient-
out (a patient's own nodes never train its scores). Numbers in `classifier_lopo.csv`:

- LOPO ROC-AUC: δ 0.806, β 0.766, low-γ 0.752, α 0.689, θ 0.677 — **permutation-p =
  0** (within-patient label shuffle). **high-γ 0.482 (perm-p 1.0) — null.**
- **Strength fails cross-patient** (AUC 0.40–0.52, ~chance) — reproducing audit_80's
  cross-patient failure for node-intrinsic markers; the **relational** features are
  what transfer. `strength+relational ≈ relational` (strength adds nothing, even
  hurts low-γ) — odds ratios: affinity and segregation both > 1.
- P(SOZ) is reasonably **calibrated** (`classifier_calibration.csv`; β top bin
  predicts 0.39 / observes 0.45; known epi median P 0.21 vs non-epi 0.067).

This is the surprising, strong result: a per-node SOZ probability that **generalises
to patients it was never trained on**, where hubness does not.

## 5. Occult-candidate discovery (audit_93, VII2/VII4)

Apply the LOPO model to a patient's **unlabelled** contacts, rank by cross-band
P(SOZ) over the airtight bands (δ/β/low-γ). Each candidate carries proximity to the
nearest labelled SOZ, anatomical region, and the feature values that fired it.
Numbers in `occult_candidates.csv`; top hypotheses:

- **Pat_06** medial-orbitofrontal nodes (P=0.87, 0.72) — **3.5–3.8 mm** from a
  labelled SOZ, flagged in all three airtight bands; high affinity (z≈+3.8) and
  segregation (z≈+2.4). Most plausible occult extension in the cohort.
- **Pat_08** superior-temporal cluster (P=0.66–0.68, ~16 mm).
- **Pat_13** hippocampal/amygdalar contacts (P≈0.2–0.4) among its already-large SOZ.

The OFC + lateral-temporal pattern resonates with the independent trace-localization
result (β → orbitofrontal + lateral-temporal) — a cross-consistency, not proof.

## 6. What the marker is **not** (honest scope)

- **Seed-based / semi-supervised.** The features are affinity to the *known* SOZ;
  the marker needs ≥ a few labelled SOZ per patient to define the community. It is
  the directive's intended "given seeds, recover the rest + rank occult" tool, **not
  a from-scratch detector**. It cannot run on a zero-label patient (Pat_15).
- **The per-node separation operationalises the audit_89 community finding** — it is
  the node-level expression of "epi co-diffuse", not independent new evidence. The
  genuinely *predictive* claims are the **masking recovery** (held-out SOZ) and the
  **cross-patient transfer** (LOPO).
- **Occult candidates are hypotheses.** No resection/outcome ground truth exists in
  this cohort. White-matter candidates (e.g. one Pat_02 hit) are lower-confidence
  (cf. WM-exclusion result). Proximity + cross-band agreement are plausibility
  filters, not validation.

## 7. Outliers, bands, caveats

- **Outliers (don't pool).** Pat_10 = consistent non-responder (epi are hubs, not a
  community); Pat_07 + Pat_03 weak/underpowered (n_epi 7/6); Pat_06/Pat_10 respond
  only at low-γ. Effect carried by Pat_08 (strongest), Pat_05, Pat_13, Pat_14,
  Pat_02. Per-patient geometry in `epi_geometry_per_patient.csv` (Pat_10/Pat_15 are
  right-hemisphere; the cohort is mostly left — not itself an explanation).
- **Bands.** δ/β/low-γ airtight (C3+C5); α survives C2+C5; θ moderate; **high-γ null
  everywhere** (the method returns null when it should). θ being a per-node SOZ
  separator (AUC 0.70) while it is an "anti-trace" in the cross-phase story is a
  reminder that the **θ "anti-trace" label needs nuance** (flagged for
  `bands/04_theta.md` + `locked/VERDICT_LEDGER.md`).
- **Every audit opened with the 5-point critical preamble**; matched-strength is the
  mandatory null at every layer.

## 8. Label-free cross-patient localization — narrows, cannot select (audit_94, 2026-06-10)

The sharp question: with zero labels in a target patient, can single-patient
knowledge localize its SOZ? Node-level cannot (audit_80; strength at chance here).
The only place it could work is the **community** level — the transferable invariant
is "the SOZ is a strength-independent diffusion community". audit_94 tests it: build
a **label-free** beyond-strength co-diffusion matrix A^z (per-pair ρ z vs the
matched-strength surrogate — the surrogate is the strength reference, no labels),
detect communities **unsupervised** (spectral clustering on A^z⁺), and either train a
cross-patient community classifier or apply a simple "most-coherent community = SOZ"
rule. Two clearly separated outcomes (`community_features.csv`,
`localization_*.csv`):

- **Within-patient narrowing WORKS (label-free).** Unsupervised clustering surfaces
  one community (~15 nodes of ~115) that captures **~50–62 % of the SOZ at ~4–5×
  enrichment** over prevalence (δ precision 0.48 / recall 0.53 / lift 5.0×; β
  0.42 / 0.62 / 4.5×; low-γ 0.38 / 0.57 / 4.2×). **5–6 of 8–9 patients** have a
  label-free community with SOZ-precision ≥ 0.40. The SOZ *is* recoverable as a
  coherent community with zero labels. Pat_10 fails (0.20) — its SOZ is not a
  community (epi-as-hubs), the expected false-negative.
- **Cross-patient SELECTION FAILS (label-free).** Picking *which* of the ~8
  communities is the pathological one, with zero target labels, is **at chance** —
  the trained classifier and the simple "most-coherent" rule both put the SOZ
  community at #1 in ~12 % of patients (chance 12.5 %), top precision ≈ prevalence
  (lift < 1), p ≈ 0.9–1.0, every band. Multiple *physiological* communities are
  equally beyond-strength-coherent; nothing label-free marks the pathological one.

**Verdict (honest).** Fully autonomous, zero-label cross-patient SOZ detection does
**not** work — confirming the marker is fundamentally **seed-based**. What the
community lens adds is **search-space reduction**: label-free, it narrows the SOZ
from ~115 contacts to one ~15-contact community at ~4–5× enrichment, capturing the
majority of the SOZ — but a label or clinical prior is still needed to *select* which
coherent community is pathological. The "unit of transfer = community" hypothesis was
right that it helps within-patient recoverability; it does not unlock cross-patient
selection. Caveat: K=8 communities fixed a priori (a free parameter); the at-chance
selection is robust to that, the narrowing magnitude is not finely tuned.

## 9. "Mark few, discover many?" — yes, as a shortlist; not clean discovery (audit_95, 2026-06-10)

VII1 (§3) hid 50 % of the SOZ — heavy seeding. The clinically exciting claim is
stronger: mark only a FEW seeds (k=2,3,5,8) and recover the MANY remaining SOZ.
audit_95 sweeps k, with a ranking metric (AUC of remaining-SOZ vs healthy) and the
HONEST hard metric — precision@top-m (flag the top-m contacts, m = remaining SOZ;
fraction truly SOZ). Numbers in `seed_curve_cohort.csv`:

- **δ (best): even k=2 works** — AUC 0.73, precision@top-m 0.29 (prevalence 0.07,
  **2.6× lift**); rising to AUC 0.85 at k=8. From 2–3 seeds the remaining SOZ are
  recovered usefully.
- **low-γ:** AUC 0.64→0.81 (k=2→8); precision 0.34–0.39 (**4–7× lift**) — strong
  enrichment even at k=2.
- **β:** AUC 0.63→0.74; precision 0.21–0.38 (**3.8–4.7× lift**).
- **α: fails at few seeds** — AUC 0.57–0.65, lift ~1.4× (≈ prevalence).
- Beats the node-strength baseline at every k (strength precision 0–0.16).

**Honest framing.** From a few seeds, the marker produces a **prioritized shortlist**
in which remaining SOZ are **3–5× enriched** and recovered at AUC 0.7–0.85 (δ best,
then low-γ/β; α not). But **absolute precision is ~30 %** — flag the top-m and most
are still false positives. So this is a genuine **search-narrowing / triage** tool
(a 30 %-precision shortlist beats random review 3–5×), **not** a clean "find all the
other SOZ" detector, and it is band-dependent (δ/low-γ/β, not α). That distinction
is the difference between an honest contribution and an overclaim.

---

## Provenance (numbers in CSVs; scripts + 2026-06-08 timestamp)
- Few-seed curve (audit_95, 2026-06-10):
  `scripts/01_compute/audit/audit_95_epi_seed_curve.py` →
  `data/audit/epi_marker_relational/seed_curve_{cohort,per_patient}.csv`.
- Label-free community localization (audit_94, 2026-06-10):
  `scripts/01_compute/audit/audit_94_epi_community_localization.py` →
  `data/audit/epi_community_localization/{community_features,localization_per_patient,
  localization_cohort}.csv`.
- C5 spatial null: `scripts/01_compute/audit/audit_90_epi_propagator_spatial_null.py`
  → `data/audit/epi_propagator_blocks/{spatial_null_cohort,spatial_null_per_patient,
  epi_geometry_per_patient}.csv`, `README_spatial_null.md`.
- Relational features: `audit_91_epi_relational_features.py` →
  `data/audit/epi_marker_relational/{node_relational_features,relational_separation}.csv`.
- Masking recovery: `audit_92_epi_masked_recovery_relational.py` →
  `…/masked_recovery_{cohort,per_patient}.csv`.
- Classifier + discovery: `audit_93_epi_soz_classifier_discovery.py` →
  `…/{classifier_lopo,classifier_calibration,occult_candidates}.csv`.
- Figures: `scripts/02_preprint/preprint_36_epi_propagator_soz_marker.py` →
  `data/preprint/figures/all_bands/fig_epi_{propagator_controls,soz_marker,
  occult_discovery}.pdf`.
- Part A finding (audit_89): `data/audit/epi_propagator_blocks/propagator_blocks_*.csv`.

## Next (optional)
- Prospective check if any cohort patient has resection/outcome metadata (would turn
  occult candidates from hypotheses into a testable prediction).
- Decide the manuscript home (dedicated epilepsy subsection under `.agents/preprint/`).
