---
name: epi-propagator-soz-marker
era: IMCOH_ABS_COHORT_N10
status: completed_2026-06-08
results_report: .agents/reports/2026-06-08_epi-propagator-soz-marker.md
kind: directive
scope: established propagator diffusion-community epi finding (audit_89, C1–C5) + forward directive to build an interpretable per-node SOZ classifier/predictor from diffusion-lens features (epi↔epi community AND epi↔non interface)
---

# Epilepsy section — propagator diffusion-community finding + SOZ-marker directive

> **STATUS 2026-06-08 — COMPLETE.** Part A closed (C5 spatial null passed:
> δ/β/low-γ + α survive, high-γ null; audit_90). Part B delivered: relational
> features (audit_91, AUC 0.74–0.86 strength-orthogonal), within-patient masking
> recovery (audit_92, AUC 0.72–0.82 > hubness), interpretable cross-patient LOPO
> classifier with calibrated P(SOZ) (audit_93, AUC 0.68–0.81, perm-p 0, strength at
> chance), and occult-candidate discovery (audit_93). Marker is **seed-based** (needs
> a few labelled SOZ/patient), occult candidates are **hypotheses**. Full results +
> honest scope: **`.agents/reports/2026-06-08_epi-propagator-soz-marker.md`**. The
> phases below are retained as the original plan of record.

**Head.** Read relationally, the LRG propagator ρ(τ)=e^{−τL}/Z shows that
**epileptic contacts form a strength-independent, band-specific diffusion
community** — they co-diffuse with each other beyond what node strength forces
and beyond strength-matched random groups (audit_89; δ/β robust, low-γ strong on
the decisive control, α on the contrast, **high-γ a clean null**). This is the
first epi per-contact signal in the whole investigation that survives its nulls —
unlike hubness, depth, the signed Laplacian, and node-level propagator features,
all of which failed beyond hubness. **Next section:** turn this diffusion-lens
signature — the epi↔epi community *and* the epi↔non interface — into an
**interpretable** per-node classifier that assigns each contact a calibrated
**P(SOZ | node)**, validate it by masking on high-epi patients, and use it to
propose **occult (unlabeled) SOZ candidates** in low-epi patients. **Post-compact
start point: Part B, Phase V.**

---

## Part A — Established finding (validated 2026-06-08)

### Provenance (numbers live in the CSVs, not in this file)
- Build script: `scripts/01_compute/audit/audit_89_epi_propagator_blocks.py` (2026-06-08).
- Numbers: `data/audit/epi_propagator_blocks/propagator_blocks_cohort.csv`
  (per band × τ) and `…/propagator_blocks_per_patient.csv` (per patient × band × τ).
  Human summary: `data/audit/epi_propagator_blocks/README.md`.
- Cohort n=9 (Pat_15 has 0 epi → structurally undefined, **not** a dropout).
- Substrate `imcoh_abs`, phase `rest_post`, LRG propagator on `L = D − W`,
  τ-grid of 6 points `1/λmax → 10/λmax`. Surrogate ensemble = the cached
  matched-strength eigs (R=200, 4-cycle ±δ, swap 20, seed 20260511) — the same
  mandatory null used by the trace battery.

### What the measure is
Partition node **pairs** by epileptic membership and read the propagator as block
means: ρ_EE (epi–epi), ρ_EN (epi–nonepi interface), ρ_NN (healthy–healthy). The
question is whether epi tissue is an internally hyper-communicating sub-network
**beyond strength** — a relational property node-strength cannot produce.

### Verdict (prose; cite the CSV rows for values)
- **Robust on BOTH epi-specific controls, broad across all 6 τ — the safest
  claims: δ and β.**
- **low-γ:** strong on the decisive specificity control (C3), weaker on the
  contrast (C2).
- **α:** strong on the contrast (C2), only moderate on specificity (C3, ~half
  the cohort).
- **θ:** moderate on both, scale-dependent (grows at coarse τ).
- **high-γ:** **null on both controls at every τ** — the internal negative that
  proves the method returns null when it should (so the positives are not a
  universal artifact).

### Control battery — C1–C5 (epi-propagator)
Numbered analogously to `locked/CONTROLS.md` (trace). Cite by number in any
writeup. Column names in brackets index the CSV.
- **C1 — matched-strength surrogate, absolute** (ρ_EE vs strength-preserving
  rewiring). PASSES but **INFLATED**: rewiring destroys *all* community
  structure, so any coherent subset beats it, and the raw value rises with
  frequency (sparser graph). **Context only — never the headline.** [`z_EE`]
- **C2 — matched-strength contrast** (ρ_EE − ρ_NN vs same surrogate): epi block
  enriched *more* than the healthy block. Epi-specific. Strong δ/α/β, weaker
  θ/low-γ, null high-γ. [`z_EE_minus_NN`]
- **C3 — strength-matched random-subset null on the REAL graph** (epi vs
  arbitrary equal-size, equal-strength groups; strength-quintile-matched,
  R_SUB=400). **Decisive specificity control — cannot be hubness** (strength is
  matched). Strong δ/low-γ/β, moderate α/θ, null high-γ. [`z_rand_EE`]
- **C4 — cross-shaft restriction** (epi pairs on *different* electrodes): rules
  out on-shaft proximity. Positive where tested. **Statistical control only — NOT
  a figure panel** (see issues). [`z_EEx`]
- **C5 — strength+spatial-matched random-subset null: CLOSED 2026-06-08
  (audit_90).** Random groups matched to the epi set on BOTH strength-quintile
  profile (exact) AND spatial spread (radius of gyration in implant x,y,z, matched
  to ~1e-6 by within-quintile swaps). **PASSED:** δ +3.05 (p.020), low-γ +3.38
  (p.002), α +1.81 (p.020), β +1.59 (p.006), θ +2.37 (p.020); high-γ null (p.455).
  The finding is now airtight — neither hubness nor spatial clustering. Numbers in
  `data/audit/epi_propagator_blocks/spatial_null_cohort.csv`. [`z_spat_EE`]

Transversal robustness (not a numbered control): the signal is positive across
**all 6 τ** (not a cherry-picked scale), and high-γ is a clean internal negative.

### Flagged finalization issues (HQ — "quartier generale")
General preprint-finalization directives recorded here so they survive compaction:
- **No number tables in `.md`.** Reference the cohort/per-patient **CSV** + the
  generating **script** + the **timestamp** — never hardcode values. (General
  rule; propagate to `WRITING_GUIDE.md` on user confirmation. Applied in this
  file.)
- **Statistical suite framed as C1–C5** (above) — complete battery for this
  finding.
- **Outliers (do not pool; report per-patient).** `Pat_10` is the consistent
  cohort non-responder (δ/β/α near or below 0); `Pat_03` is underpowered
  (n_epi=6, marginal β/α); `Pat_07` drops out in low-γ. Effect carried by
  `Pat_08` (strongest), `Pat_13`, `Pat_05`, `Pat_02`, `Pat_14`, `Pat_06`.
  Investigate Pat_10 / Pat_07 coverage & SOZ extent (V2 below).
- **Cross-probe (C4 / "X-epi") is a control, not a visualization.** Figures keep
  the **full** epi–epi block (all pairs incl. cross-shaft); do **not** carve out
  a separate cross-probe panel — "let the X-epi be there." C4 stays a statistical
  check only.
- **θ "anti-trace" is more complicated than a single label** (TRACE section, not
  epilepsy). Flag for `bands/04_theta.md` revision; do not carry the bare
  "θ = anti-trace" claim into the manuscript without nuance. Recorded here so it
  is not lost across compaction; belongs in `locked/VERDICT_LEDGER.md` +
  `bands/04_theta.md`.
- **Literature framing:** "the epileptogenic zone is internally hyperconnected"
  is already known. Our genuinely new part is the **strength-independent,
  multiscale, band-resolved** diffusion characterization. Frame as *sharpening* a
  known phenomenon, not discovering it.

### Manuscript home
No existing band brief covers this (`bands/` is the cross-phase **trace** story;
this is the orthogonal **epilepsy** axis). Decide with the user whether the
epilepsy results get a dedicated subsection file under `.agents/preprint/`; until
then **this directive is the source of truth** for the finding.

---

## Part B — Next-section directive: interpretable per-node SOZ classifier

**Goal (user, 2026-06-08).** Turn the diffusion-lens signature into a
**predictor**: assign every contact a calibrated **P(SOZ | node)** from
**interpretable** propagator features (NOT a black box), validate by masking known
SOZ on high-epi patients, then **propose occult SOZ candidates** in low-epi
patients with quantified probability. Use **both** the epi↔epi community **and**
the epi↔non **interface** — the SOZ→non linking pattern is as important because it
*separates* SOZ from the rest. **Do not stop until the marker predicts held-out
SOZ above nulls with calibrated probability AND emits ranked occult candidates
with interpretable reasons.**

Run all phases in order; iterate until results are satisfying.

### Phase V — finish validating the measure (gate before any classifier claim)
- **V1 (C5): spatial-matched random-subset null.** Draw random groups matched to
  the epi set's spatial spread (use implant `(x,y,z)` from `implant_pat_NN.csv`,
  per `feedback_implant_anatomy_not_letters.md` — not channel letters). Close the
  anatomy confound. Outcome decides scope: all of δ/β/low-γ survive → full claim;
  only β survives → β-led claim.
- **V2: per-patient outlier analysis.** Why Pat_10 (δ/β/α) and Pat_07 (low-γ)
  don't respond — implant coverage, SOZ size/laterality, n_epi power. Report;
  never drop (`feedback_no_patient_dropout.md`).
- **V3: interface (EN) characterization as a SEPARATE signal.** Is ρ_EN
  distinctive — do SOZ nodes link to non-SOZ with a characteristic diffusion
  fingerprint (SOZ as diffusion sources/sinks; EN suppressed or directionally
  structured)? This boundary signal is what *separates* SOZ from non-SOZ and
  feeds feature f2 below. Test vs C1/C3 nulls.

### Phase VI — interpretable per-node classifier
- **VI1: small set of INTERPRETABLE per-node diffusion features**, all
  strength-orthogonal (matched-strength surrogate-z, so beyond hubness):
  - `f1` affinity-to-SOZ-community: mean ρ(node, known-epi) — diffusion coupling
    into the SOZ block.
  - `f2` interface/segregation: contrast of ρ to epi vs ρ to non
    (e.g. ρ̄_epi − ρ̄_non, or a diffusion participation/boundary score) — encodes
    the EN linking pattern (V3).
  - `f3` within-community diffusion centrality restricted to the SOZ block.
  - `f4` (optional) multiscale shape of f1–f3 across τ.
  - Per band, in the bands where the community is real (δ/β/low-γ; test α; skip
    high-γ).
- **VI2: interpretable model** — logistic regression or monotone scoring on
  f1–f4. **Report weights / odds-ratios** (interpretability is a hard
  requirement). Output **calibrated P(SOZ | node)** (Platt/isotonic; show
  calibration curve). Must beat the **strength-only baseline** (node-level
  strength failed beyond hubness in audit_86/87 — the relational features are the
  bet).
- **VI3: nulls on the CLASSIFIER itself** (not just features): the per-node score
  must beat matched-strength + strength-matched-random expectation.

### Phase VII — validation by masking, then discovery
- **VII1: masking recovery on high-epi patients** (Pat_13 n=30; Pat_02/05 n=14;
  Pat_14 n=12): hide a fraction of labeled SOZ, train on visible + healthy,
  recover hidden SOZ. Report ROC-AUC / PR-AUC + calibration; leave-one-epi-out
  for fine-grained recovery. This is the rigorous within-patient validation (the
  diffusion-affinity feature is the one that beat its nulls).
- **VII2: discovery in low-epi patients** — apply the calibrated classifier to
  **unlabeled** contacts → rank → propose top candidates as possible occult SOZ
  with P(SOZ) + uncertainty. **No ground truth exists** for unlabeled contacts
  (no resection/outcome in cohort) → frame as **hypotheses**: check anatomical /
  spatial plausibility (proximity to known SOZ, plausible regions) and cross-band
  agreement (a candidate flagged in δ AND β AND low-γ is stronger).
- **VII3: probability calibration + cohort summary** — per-candidate P(SOZ) with
  CI; expected occult-SOZ count per patient; rank-stability under bootstrap.
- **VII4: interpretability report** — for each candidate, WHY it fired (affinity
  to SOZ community vs interface pattern), so it is explainable to a clinician, not
  a black-box flag.

### Cross-cutting constraints
- **Interpretable, not black box** (logistic / monotone; report weights).
- **Matched-strength control mandatory at every layer** (features and classifier
  strength-orthogonal) — `feedback_matched_strength_mandatory.md`.
- **Per-band**; lead with δ/β/low-γ (real community), test α, skip high-γ (null).
- **Use BOTH epi↔epi (community) AND epi↔non (interface)** — interface required
  (it separates SOZ from non-SOZ).
- **Terminology:** dataset labels are clinically-marked epileptic contacts; "SOZ"
  per user usage. Never state a proposed candidate *is* SOZ — report **P(SOZ) as
  a hypothesis** (no ground truth for occult nodes).
- **Brutal honesty / 5-point critical preamble before each new audit**
  (`feedback_critical_null_preamble.md`); state limitations first.

### Infra to reuse (don't reinvent)
- `audit_89` block machinery: `block_means`, `strength_matched_subsets`,
  `subset_offdiag_means`.
- Matched-strength eigs: `load_or_compute_eigs_at_path`,
  `es.surr_eig_path` / `es.cell_rng` (`_epi_stratify.py`); cached rest_post,
  seed 20260511.
- `build_epi_masks` (FC-order epi mask + probes); `load_phase_fc` (imcoh_abs).
- `diffusion_from_eigs` (audit_86) for per-node propagator features (ρ_ii, Hdiff,
  T_comm).
- **Lesson from the negatives** (audit_85/86/87 node strength/recovery; audit_88
  signed Laplacian): node-level / strength features fail beyond hubness. The win
  is **RELATIONAL** (block / affinity / interface). Build the classifier on
  relational features, not node-strength.
- History report: `.agents/reports/2026-06-05_propagator-subspace-epi-recovery.md`
  (the negatives + depth-confound retraction). This directive is the forward plan.

### Done-when
A per-band, interpretable classifier that (1) recovers held-out SOZ above
matched-strength + random nulls with calibrated probability on high-epi patients,
(2) emits ranked occult-SOZ candidates with P(SOZ) + CI + interpretable reasons on
low-epi patients, and (3) has the spatial confound (C5) closed.
