---
name: propagator-subspace-epi-recovery
type: report
era: IMCOH_ABS_COHORT_N10
status: current
date: 2026-06-05
created: 2026-06-05
updated: 2026-06-05
companion:
  - data/audit/epi_marker/
sources:
  - audit_85  epi_propagator_recovery            # within-patient template recovery (raw)
  - audit_86  epi_propagator_matched_strength     # FAIR baseline + matched-strength certification
  - audit_87  epi_masked_recovery                 # hide-and-seek validation + within-shaft control (DECISIVE)
pointers:
  - .agents/guides/task-persistence-investigation/2026-06-05_epi-node-trace-marker.md   # scope §11
  - .agents/reports/2026-06-05_epilepsy-headline-and-occult-node-marker.md              # parallel occult report (consistent)
  - .agents/preprint/bands/02_alpha.md
  - .agents/preprint/bands/01_beta.md
---

# Within-patient propagator-subspace recovery of epileptic contacts — certified

**Head (FINAL, after three controls).** There is **no dynamics-based epileptic-
contact marker** — not in the LRG propagator ρ(τ) subspace, not in the signed-
ImCoh magnetic Laplacian. Within a patient, epileptic contacts *are* recoverable,
but the recoverable signal decomposes entirely into two **geometric / known-
biology** axes: **node strength (hubness)** at the gross level, and **along-shaft
depth** at the fine (within-shaft) level — epileptic contacts sit at
characteristic electrode depths. When we hide half the known epileptic contacts
and search within the same shafts, plain **contact depth recovers them at AUC
0.72–0.74 in every band, and *beats* both the propagator (0.60–0.67) and the
signed-Laplacian (0.60–0.71) features; adding the dynamics to depth improves
nothing (all p ≥ 0.33).** Cross-patient transfer is null throughout. This report
is the full audit trail of THREE successive "positive" results that each
dissolved under the right control — (1) a β "breakthrough" that was an unfair
strength baseline, (2) a "shaft-independent" within-shaft recovery that was
along-shaft depth — kept in full so the dead ends are not re-walked. The honest
deliverable: epileptic contacts are hubs at characteristic depths; our LRG /
diffusion / directional machinery adds nothing a clinician's electrode map does
not already give.

> **Brutal-honesty banner.** The first internal pass (audit_85) reported
> "propagator β recovery 0.720 > strength 0.622" and was framed as a major
> breakthrough. That comparison scored the strength baseline by *template
> distance*, which cripples a single monotone feature. Under the **fair**
> strength baseline (within-patient Mann-Whitney AUC) strength is **0.738 at β**
> and **beats** the propagator *on the full contact set*. The β headline was
> retracted. A second "positive" — the within-shaft hide-and-seek of §9 — then
> **also fell** to the along-shaft DEPTH confound (§10): contact depth recovers
> hidden epi within-shaft better than the dynamics, which add nothing on top of
> it. Sections 3–9 are the chronological audit trail; **§10 is the final
> verdict: no dynamics marker; epileptic contacts are hubs at characteristic
> depths (geometry).** Two retractions, both caught by the next control — that is
> the process working, not a result.

---

## 1. The question and the reframe

The preprint's epilepsy section had reached a clean negative on a *cross-patient*
epileptic-node marker: per-node identity features (audit_79/80) and per-node
leverage (audit_81) are all at chance in leave-one-patient-out (LOPO), because
absolute FC strength lives on a patient-specific scale and no global threshold
transfers. The parallel occult-node report
(`2026-06-05_epilepsy-headline-and-occult-node-marker.md`) reached the same
verdict from a third angle (leave-node-out influence on the cohort trace).

The reframe that motivated this work: "rediscover an unlabelled epileptic node"
is **not** a cross-patient classification problem. It is a **within-patient**
anomaly/template problem — the realistic clinical scenario where a clinician has
marked *some* contacts and wants the rest. Self-normalising within each patient
removes the FC-scale confound that killed the cross-patient approach.

The feature space is the **LRG propagator** the user asked for: the information-
diffusion density matrix `ρ(τ) = e^{−τL}/Z` (the canonical LRG object,
`lrgsglib.compute_laplacian_properties`), summarised per node as

- **return probability** `ρ_ii(τ)` (diffusion localisation at the node),
- **diffusion-row entropy** `H_i(τ) = −Σ_j p_ij log p_ij`, `p_ij = K_ij/Σ_j K_ij`
  (how far the node's diffusion spreads),
- **communication distance** `T_ρ = 1/ρ_ij` summarised as mean / min over the
  node's row (how central the node is in the ultrametric),

over a per-patient τ-grid `geomspace(1/λ_max, 10/λ_max, 6)`, plus the per-pair
trace features (audit_79) and leverage (audit_81), standardised **within
patient**.

---

## 2. Within-patient recovery (audit_85) — works, but the baseline was unfair

Method: per (patient, band), standardise features within patient; **leave one
epileptic node out**, build the epi template = mean feature vector of the *other*
epileptic contacts, score every node by similarity (−distance) to the template,
and measure the recovery AUC = P(held-out epi ranks above a non-epi).

audit_85 reported within-patient recovery clearly above chance (full subspace
β 0.715, θ 0.679, low-γ 0.637, α 0.629; up to 0.85 in individual patients) and a
"propagator beats strength" gap. **The flaw:** the strength baseline was scored
with the *same template-distance* metric. Template distance flags nodes *near the
epi mean* and penalises nodes *more extreme* than the epi mean — which is exactly
wrong for a monotone feature like strength, where "more is more epileptic." So
the strength baseline (0.622 at β) was artificially depressed, manufacturing a
gap that is not real.

---

## 3. Certification (audit_86) — fair strength + matched-strength surrogate

Two controls, both mandatory before any "beats strength / strength-orthogonal"
claim:

1. **Fair strength baseline** — within-patient **Mann-Whitney AUC** of node
   strength (epi rank high), the natural monotone scoring.
2. **Matched-strength surrogate-z** — recompute every propagator feature on the
   node's R=200 matched-strength surrogate ensemble (the cached rest_post 4-cycle
   ±δ eigendecompositions, seed 20260511) and z-score the observation against it.
   The surrogate fixes each node's strength **exactly**, so any epi separation
   left in the surrogate-z feature is **strength-orthogonal by construction**.

**Recovery AUC (within-patient, cohort median over patients):**

| band | fair strength (MW) | propagator raw | propagator **strength-removed (MS-z)** | MS-z vs chance (Wilcoxon) | reading |
|---|---|---|---|---|---|
| δ | 0.639 | 0.606 | 0.548 | p 0.125 (ns) | hubness only |
| θ | 0.672 | 0.589 | 0.646 | p 0.020 | MS-z ≈ strength |
| **α** | 0.600 | 0.622 | **0.667** | **p 0.049** | **MS-z > strength** ✓ |
| β | **0.738** | 0.720 | 0.606 | p 0.034 | **hubness dominates** |
| low-γ | 0.610 | 0.613 | 0.608 | p 0.027 | MS-z ≈ strength |
| high-γ | 0.545 | 0.585 | 0.541 | p 0.102 (ns) | weak / hubness |

Two facts from this table:

- **β is hubness.** Fair strength (0.738) beats the propagator in every form
  (raw 0.720, strength-removed 0.606). The audit_85 β headline does not stand.
- **The strength-orthogonal propagator signal is real but modest, and it is at
  α.** The strength-removed propagator recovers α epi at **0.667**, exceeding the
  fair strength baseline (0.600); it clears its random-set permutation null at
  θ/α/β/low-γ but only *exceeds fair strength* at α.

Per-feature MS-z separation (within-patient, oriented AUC) departs from 0.5 — at
β `Tcomm_min` 0.75 and `Hdiff` ≈ 0.71, at α `Hdiff` ≈ 0.62 — but **no single
propagator feature is individually significant at n = 9**; the signal lives in
the combined subspace, not in any one descriptor.

---

## 4. The genuine finding — a strength-orthogonal α propagator signal

At α, the part of the propagator that *cannot* be explained by node strength
still recovers epileptic contacts above the fair strength baseline (0.667 vs
0.600, p = 0.049). This is modest and borderline, but it is **matched-strength
certified** and it is **mechanistically coherent**: the independent epi-
stratification result (`2026-06-05_epilepsy-headline-and-occult-node-marker.md`
§1.3) shows that **α uniquely recruits epileptic tissue into the task trace** —
the α epi↔epi pair class carries the trace strongly (ρ_split +0.404, p = 0.006,
6/9 patients), whereas β *spares* the epileptic core. So α is precisely the band
where epileptic contacts have a distinctive *communication-geometry* signature
beyond being hubs, and that is exactly where the strength-removed propagator
recovers them. The two analyses, built independently, point the same way.

This is the defensible, publishable claim: **a modest, band-specific (α),
strength-orthogonal propagator-geometry signature of epileptic contacts,
detectable within patient.** It is a supporting/mechanistic result, not a
clinical marker.

---

## 5. What does NOT hold (retractions)

- **No β breakthrough.** Within-patient β epileptic recovery is node strength
  (hubness). Retract any "propagator beats strength at β."
- **No cross-patient marker.** Unchanged from audit_80/81 and the occult report:
  nothing transfers across patients; an unlabelled patient cannot be scored from
  a model trained on others.
- **Discovery candidates remain hypothesis-only.** The within-patient candidate
  lists (`propagator_candidates.csv`; closest unlabelled contacts to the epi
  template) are dominated by hubness and have no surgical-outcome ground truth in
  this cohort. They cannot be called predicted foci.

---

## 6. Bounds and caveats

- **Within-patient only.** Needs ≥ 3 labelled epileptic contacts to build the
  template (the partial-label scenario). Pat_15 (0 epi) is excluded.
- **n = 9 patients with epi.** Per-feature tests are underpowered; the α result
  (p = 0.049) is borderline and should be reported as such, with the α-recruits-
  epi coherence as the corroborating argument rather than the p-value alone.
- **AUC ~ 0.6–0.7 is useful, not clinical-grade.** Heterogeneous across patients.
- **The matched-strength surrogate** preserves node strength exactly but drifts
  the weight distribution (KS ≈ 0.36); it is the more conservative null on paired
  cross-phase tests (see `matched_strength_null_defense_2026_05_28`).

---

## 7. Files

| path | what |
|---|---|
| `scripts/01_compute/audit/audit_85_epi_propagator_recovery.py` | within-patient template recovery (raw; renamed from a colliding audit_82) |
| `scripts/01_compute/audit/audit_86_epi_propagator_matched_strength.py` | fair strength baseline + matched-strength surrogate-z certification |
| `scripts/02_preprint/preprint_34_epi_recovery.py` | honest certification figure |
| `data/audit/epi_marker/node_propagator_features.csv` | per-node propagator ρ(τ) features |
| `data/audit/epi_marker/recovery_per_band.csv` | audit_85 recovery (raw, with the unfair baseline — kept for provenance) |
| `data/audit/epi_marker/propagator_ms_features.csv` | per-node raw + matched-strength surrogate-z propagator features |
| `data/audit/epi_marker/propagator_ms_certify.csv` | the certified table (fair strength, raw, MS-z) |
| `data/audit/epi_marker/README_ms_certify.md` | audit_86 self-documenting README |
| `data/preprint/figures/all_bands/fig_epi_recovery.pdf` | (a) three-way recovery (b) strength-orthogonal margin |

Scope: `.agents/guides/task-persistence-investigation/2026-06-05_epi-node-trace-marker.md`
§11 (RESULT banner + CERTIFICATION correction).

---

## 8. Honest verdict

Within-patient epileptic-contact recoverability is **hubness-dominated**
(strongest at β). On top of that there is a **modest, matched-strength-certified,
strength-orthogonal propagator-geometry signature that is α-specific**
(strength-removed recovery 0.667 > fair strength 0.600, p = 0.049), coherent with
α's recruitment of epileptic tissue. This is a genuine mechanistic result worth a
sentence or two in the α subsection — **not** a cross-patient biomarker, **not** a
β headline. The earlier breakthrough framing was a baseline-scoring artifact and
has been retracted here in full.

**Possible next step (only if pursued):** a within-class matched-strength null
for the α epi↔epi trace and a per-feature α propagator analysis at larger R to
firm up the borderline p = 0.049 — but the cohort is n = 9 epi-bearing patients,
so the ceiling on significance is low. No clinical claim is available without
surgical-outcome ground truth.

---

## 9. Hide-and-seek validation with shaft control (audit_87) — the decisive test

The user's validation design, done properly: hide a random 50% of each patient's
known epileptic contacts (≥6 epi required), learn the "epileptic direction" in
feature space from the **visible** half via a Fisher linear discriminant (fair to
a single monotone feature and to a multi-feature subspace alike — removing the
template-distance handicap of §2), and recover the hidden half. 200 random splits
per patient; recovery AUC = P(hidden epi ranks above a healthy contact).

Two negative pools:
- **all** healthy contacts;
- **within-shaft**: healthy contacts on the **epi-bearing shafts only**, so shaft
  position is constant and the scorer must find the epi *contact*, not the epi
  *shaft*. This is the decisive control (the shaft-proximity confound that the
  earlier passes had not addressed).

Recovery AUC (median over the 9 epi-bearing patients), strength-removed
(matched-strength surrogate-z) propagator vs fair node strength:

| band | strength (all) | propagator (all) | **strength (within-shaft)** | **propagator (within-shaft)** | propagator within-shaft > chance | prop − strength (within-shaft) |
|---|---|---|---|---|---|---|
| δ | 0.633 | 0.689 | **0.495** | **0.610** | p 0.010 | **+0.081** |
| θ | 0.621 | 0.683 | 0.590 | 0.632 | p 0.002 | +0.030 |
| α | 0.545 | 0.601 | 0.687 | 0.633 | p 0.002 | −0.037 |
| β | 0.741 | 0.707 | 0.595 | 0.642 | p 0.020 | −0.041 |
| low-γ | 0.605 | 0.708 | 0.540 | 0.666 | p 0.020 | −0.011 |
| high-γ | 0.546 | 0.623 | **0.499** | **0.607** | p 0.014 | **+0.054** |

**What this establishes:**

1. **A real, strength- AND shaft-independent marker.** The strength-removed
   propagator recovers hidden epileptic contacts above chance **in every band,
   within shaft** (p ≤ 0.02). Holding shaft constant and removing strength, it
   still finds the epi contact (AUC ≈ 0.61–0.67). Six independent bands pointing
   the same way is not a fluke. This is the cleanest evidence in the whole
   investigation that the epileptic-contact signal is a genuine per-contact
   communication-geometry property, not hubness and not electrode placement.
2. **It matters most exactly where hubness fails.** Within shaft, node strength
   drops to chance at δ (0.495) and high-γ (0.499) — contacts on one shaft have
   similar strength — and there the propagator beats it outright (Δ +0.08, +0.05).
   Where strength is still informative within-shaft (α 0.687, the
   epi-recruitment band), the propagator ties/trails it.
3. **The earlier "modest α only" (§3–4, from the template-distance leave-one-out
   in audit_86) understated it.** With a fair scorer and the shaft control, the
   strength-orthogonal propagator signal is present in all bands. α remains
   special only in that strength *also* works there.

**Bounds (unchanged).** AUC ≈ 0.61–0.67 = useful, not clinical-grade;
within-patient (needs labelled epi to learn the direction); no cross-patient
transfer; n = 9 epi-bearing patients. Discovery candidate lists remain
hypothesis-only (no surgical-outcome ground truth).

Artifacts: `data/audit/epi_marker/{masked_recovery_per_band.csv,
masked_recovery_per_patient.csv, README_masked_recovery.md}`;
script `scripts/01_compute/audit/audit_87_epi_masked_recovery.py`.

**Revised verdict (supersedes §5/§8's "hubness, modest α"):** there IS a real,
matched-strength-certified, shaft-independent within-patient propagator marker of
epileptic contacts, robust across all six bands and most valuable for the
fine-grained within-shaft call where hubness is uninformative. It is a
within-patient discovery aid (find more epi given some), not a cross-patient
biomarker, and not yet clinical-grade.

---

## 10. The along-shaft DEPTH confound — overturns §9; final verdict

§9's within-shaft recovery looked like a real, shaft-independent marker. It is
not. The control §9 missed: **along-shaft depth.** Epileptic foci sit at
characteristic electrode depths (deep/mesial contacts), so *within* a shaft,
epileptic and healthy contacts differ in **position**, and any feature that
varies smoothly along the shaft (propagator localisation, directional flow)
partially encodes that position — which is why they showed within-shaft recovery.

Adding `contact_index_norm` (along-shaft depth) as a baseline settles it
(within-shaft recovery AUC, median over patients):

| band | depth alone | propagator (MS-z) | signed-Laplacian | propagator>depth? | dynamics add to depth? |
|---|---|---|---|---|---|
| δ | **0.74** | 0.61 | 0.71 | p 0.97 (no) | p ≥ 0.82 (no) |
| θ | **0.73** | 0.62 | 0.60 | p 0.92 (no) | p ≥ 0.75 (no) |
| α | **0.73** | 0.66 | 0.67 | p 0.88 (no) | p ≥ 0.46 (no) |
| β | **0.73** | 0.65 | 0.60 | p 0.85 (no) | p ≥ 0.63 (no) |
| low-γ | **0.74** | 0.67 | 0.64 | p 0.79 (no) | p ≥ 0.33 (no) |
| high-γ | **0.72** | 0.61 | 0.62 | p 0.92 (no) | p ≥ 0.63 (no) |

- **Depth alone (0.72–0.74) beats both the propagator and the signed-Laplacian
  features in every band.** The dynamics do **not** beat depth (all p ≥ 0.79).
- **Adding the dynamics to depth improves nothing** in any band (all p ≥ 0.33;
  several combinations are slightly worse from added noise).
- Therefore the §9 within-shaft signal was **along-shaft depth**, a geometric /
  known-anatomy feature (same class as the audit_80 electrode-geometry result),
  not a propagator or directional dynamics marker.

**Final verdict (supersedes §8 and §9).** Within-patient epileptic-contact
recoverability decomposes into **node strength (hubness, gross level) + along-
shaft depth (fine level)** — both geometric / already known to the implanting
team. **Neither the LRG propagator subspace nor the signed-ImCoh magnetic
Laplacian adds a marker beyond these.** Cross-patient transfer is null. There is
no novel dynamics-based epileptic-node marker in this cohort with these measures.
The signed-ImCoh / magnetic-Laplacian direction (audit_88,
`data/audit/epi_signed_laplacian/`) is **Stage-1 negative**: its within-shaft
recovery is depth, so the **Stage-2 sign-randomising surrogate is NOT warranted.**

Check script: the depth-baseline comparison is a standalone verification (run
inline; numbers above). audit_87/88 within-shaft tables remain on disk as the
*pre-depth-control* record, now superseded by this section.
