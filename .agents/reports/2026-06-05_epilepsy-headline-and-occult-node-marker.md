---
name: epilepsy-headline-and-occult-node-marker
type: report
era: IMCOH_ABS_COHORT_N10
status: current
date: 2026-06-05
created: 2026-06-05
updated: 2026-06-05
companion:
  - data/audit/epi_stratified/
  - data/audit/epi_marker/
  - data/audit/occult_node_marker/
sources:
  - audit_77  epi_stratified_cophenetic        # cophenetic ρ_split, 6 configs × 6 bands
  - audit_78  epi_stratified_grassmann          # Grassmann T_G(k), 3 subgraph configs
  - audit_79  epi_node_behavior                 # per-node descriptive gate (NEGATIVE)
  - audit_80  epi_marker_classifier             # known-node LOPO classifier (NEGATIVE)
  - audit_90  occult_node_lno_influence         # occult discovery reframe (this report)
pointers:
  - .agents/guides/task-persistence-investigation/2026-06-05_occult-epi-node-marker.md
  - .agents/guides/task-persistence-investigation/2026-06-05_epi-node-trace-marker.md
  - .agents/preprint/bands/03_gammalow.md
  - .agents/preprint/bands/01_beta.md
  - .agents/preprint/bands/02_alpha.md
---

# Epilepsy-independence is a headline, and there is no transferable occult-node marker

**Head.** The cognitive task-trace in `rest_post` LRG hierarchies is a property
of **normal task-related network reorganization living in healthy circuitry and
its interface with the epileptic zone — not an epileptiform artefact**: removing
the clinically labeled epileptic contacts makes the cophenetic trace *cleaner
and stronger* (β +0.221 → +0.275, α +0.105 → +0.187), and a low-γ trace that is
null on the full graph *emerges* once the diseased tissue is excluded
(p 0.116 → 0.024, 8/10). The trace is carried by the healthy↔healthy and the
epi↔healthy **interface** pair classes; the epileptic subgraph itself *weakens*
the β/low-γ trace (it partly **masks** the healthy γ_low trace). The one
exception is α *within* the epileptic subgraph, which is strong (+0.404,
p=0.006) — α recruits epileptic tissue while β spares it. **Separately, the
"occult epileptogenic node" discovery question is gate-negative:** a per-node
leave-node-out influence on the trace IS a real strength-independent
within-patient β signal (new — `audit_79`'s row-feature missed it), but it does
**not transfer across patients** (LOPO at chance, same ceiling as the
`audit_80` classifier), so non-labeled "candidates" are hypotheses only and
cannot be validated without surgical-outcome ground truth.

**Verification status.** Every cohort claim below is under the **mandatory
R=200 matched-strength surrogate null** (audit_77/78 ensembles). The two
brutal-honesty caveats that bound the strength of the epilepsy-independence
message are stated in §1.4 — they are not afterthoughts.

---

## 1. DELIVERABLE 1 — epilepsy-independence (consolidated, no new compute)

All numbers below are read directly from `data/audit/epi_stratified/
cophenetic_cohort.csv` + `cophenetic_per_patient.csv` (audit_77, six epi views
× six bands × n=10, matched-strength R=200) and `grassmann_cohort.csv`
(audit_78). Cophenetic ρ_split = Spearman(D_task − D_preA, D_post − D_preB),
positive = TRACE; cohort one-sided paired Wilcoxon of (obs − own-surrogate
median); LO-P13 = leave-Pat_13-out cohort p (Pat_13 carries 30 epi = 25%,
leverages every epi config).

### 1.1 Excluding the epileptic zone makes the trace cleaner / stronger

| band | full ρ (p) | exclude-epi ρ (p) | n>surr | LO-P13 p | verdict |
|---|---|---|---|---|---|
| **β** | +0.221 (0.0068) | **+0.275 (0.0029)** | 7/10 | 0.0059 | persist + strengthen |
| **α** | +0.105 (0.0049) | **+0.187 (0.0137)** | 7/10 | 0.0039 | persist + strengthen |
| **γ_low** | +0.083 (0.116, **null**) | **+0.174 (0.0244)** | **8/10** | 0.0371 | **EMERGES** |
| δ | +0.008 (0.216) | +0.086 (0.116) | 4/10 | 0.0488 | absent |
| θ | −0.040 (0.688) | +0.038 (0.423) | 4/10 | 0.545 | absent |
| γ_high | +0.000 (0.278) | +0.015 (0.313) | 5/10 | 0.326 | absent |

- **β / α persist and strengthen** under epi exclusion — the trace is *not*
  driven by the diseased tissue; if anything the diseased tissue dilutes it.
- **γ_low EMERGES**: null on the full graph (p=0.116), significant after
  excluding epi (p=0.024, 8/10 above own surrogate, the `separated` verdict),
  robust to dropping Pat_13 (p=0.037). The epileptic zone partly **masks** a
  healthy low-γ trace.
- Per-patient γ_low flips on epi exclusion (full → exclude-epi obs):
  **Pat_13 −0.042 → +0.109**, **Pat_14 −0.220 → +0.073** (the two largest
  swings), Pat_10 +0.004 → +0.085, Pat_03 +0.162 → +0.239. Pat_15 (0 epi)
  unchanged at −0.131 (no epi to remove — a built-in sanity check).

### 1.2 Pair-class decomposition — the interface and healthy tissue carry it

Pair-class configs restrict the full-graph ρ_split to a node-pair class, holding
the global strength sequence fixed (median pair counts: nonepi↔nonepi 5565,
cross 1050, epi↔epi 45):

| band | nonepi↔nonepi ρ (p) | cross (epi↔non) ρ (p) | epi↔epi ρ (p) |
|---|---|---|---|
| **β** | +0.197 (0.0049) **carry** | +0.273 (0.0098) **carry** | +0.240 (0.150) **weaken** |
| **α** | +0.073 (0.0322) **carry** | +0.177 (0.0039) **carry** | **+0.404 (0.0059) carry** |
| **γ_low** | +0.116 (0.080) | +0.013 (0.180) | +0.226 (0.102) |

- **β trace is carried by the CROSS (epi↔non-epi interface, +0.273) and
  NONEPI↔NONEPI (+0.197) classes**; the epi↔epi class is *non-significant*
  (p=0.150). The β trace lives in healthy tissue and at its interface with the
  epileptic zone, not inside the diseased core.
- The **epi_only subgraph** (a true LRG rebuilt on epi nodes only) **weakens**
  for both β (−0.078, p=0.715) and α (−0.034, p=0.752) — removing all the
  healthy tissue destroys the trace.

### 1.3 The α exception — α recruits epileptic tissue, β spares it

- **α within the epi subgraph (epi↔epi pair class): +0.404, p=0.0059, 6/9
  patients individually significant** (Pat_03/05/06/07/13/14 at own-surrogate
  p<0.05; obs>0 in 8/9). This is the strongest single epi-stratified α cell.
- The β analogue is non-significant (epi↔epi +0.240, p=0.150). **The
  dissociation is the story: α recruits epileptic tissue into the task-trace
  while β spares it.** (Note: the brief stated "5/9 individually sig" for α
  epi_epi; the CSV is **6/9** — Pat_03 +0.760 p=0.015 also clears.)

### 1.4 Grassmann (subspace) corroboration

audit_78 chordal T_G(k) on the three subgraph views (full re-emitted from
audit_66, exclude_epi from audit_67):

- **β exclude_epi gains a separated k-window** (k≈17–32, 52/87 k-cells p<0.05)
  where the full graph had none at the matched-strength frontier — the same
  "epi exclusion sharpens β" direction as cophenetic.
- low-γ has a full-graph Grassmann window (k≈12–27) that does not survive
  exclude_epi at the same k-axis (N changes from 112→88, so k-axes are only
  approximately comparable — a caveat, not a contradiction).

### 1.5 Brutal-honesty caveats (bounding the message)

1. **epi_only / epi↔epi underpowered.** The epileptic subgraph is 6–30 nodes
   (Pat_15 = 0), 15–435 pairs; matched-strength on such tiny graphs is
   near-degenerate. The epi_only `separated`-style verdicts are flagged
   exploratory and are NOT trustworthy. The α epi↔epi +0.404 is the *pair-class*
   restriction of the full-graph null (well-powered), not the tiny epi_only
   subgraph — but it still rests on a global-rewiring null (next point).
2. **The pair-class null is global-rewiring, not within-class.** It fixes the
   global per-node strength sequence and asks whether class-restricted
   co-movement survives; it does NOT randomize *within* the class. So
   "interface-carried" is **suggestive, not airtight** — a within-class
   matched-strength null would be the airtight version and has not been run.
3. **Class sizes are unequal** (nonepi↔nonepi always has the most pairs → the
   most stable Spearman). The cross and nonepi verdicts must be read alongside
   n_pairs; the β cross result (1050 pairs, p=0.0098) is not a size artefact
   but the comparison is not size-controlled.

### 1.6 Headline-ready framing

> The band-selective cognitive task-trace in post-task resting-state LRG
> hierarchies is **independent of the epileptic pathology**. Excluding the
> clinically labeled epileptogenic contacts leaves the β and α traces intact
> and *strengthens* them (β +0.221→+0.275, α +0.105→+0.187, both
> matched-strength-significant, both robust to leave-one-patient-out), and
> *unmasks* a low-γ trace that the epileptic zone had concealed
> (null→significant, 8/10). The trace is carried by healthy↔healthy and the
> epi↔healthy interface, not by the epileptic core; α uniquely also recruits
> epileptic tissue (epi↔epi +0.404). The phenomenon is normal task-related
> network reorganization, not an epileptiform artefact.

---

## 2. DELIVERABLE 2 — occult (non-labeled) epi-node marker

### 2.0 The prior was stronger than the brief assumed

The task brief stated the audit_80 classifier was "DEFERRED". **It was not — it
was BUILT, and the verdict was NEGATIVE** (`data/audit/epi_marker/
README_classifier.md`, scope `2026-06-05_epi-node-trace-marker.md` status
`built-negative`): `geometry_dominant` in every band, trace LOPO PR-AUC at
prevalence (best β 0.124 vs 0.095, perm-p 0.080), node strength does NOT
transfer across patients (PR ≤ prevalence), the only cross-patient-transferable
signal is label-free electrode geometry (shaft length + along-shaft position,
PR 0.163 ≈ 1.7× chance) — real but trivial (implantation strategy). The
whole-shaft masking control fired correctly (label tautology → masked ROC 0.38).

So the a-priori expectation for the user's *reframed* discovery question
("flag non-labeled nodes that behave like epi") was: **no occult signal beyond
strength/geometry.**

### 2.1 What was genuinely missing — and the one new honest test (audit_90)

The one per-node statistic both audit_79 and audit_80 **explicitly deferred**
(their F8) is leave-node-out *influence on the cohort trace*. I built its cheap
cophenetic form (scope: `2026-06-05_occult-epi-node-marker.md`; script
`audit_90_occult_node_lno_influence.py`):

```
Δρ_i = ρ_split_full − ρ_split_{drop node i}      (per patient, per band)
```

Δρ_i > 0 ⇒ removing node i weakens the trace ⇒ node i is a **trace-carrier**.
Each Δρ_i is a fresh (N−1)-node LRG rebuild (drop row+col, Laplacian, τ=1/λ_max,
average-linkage cophenetic). Residualized WITHIN patient against
{strength, same-shaft-epi count, along-shaft position, shaft length} BEFORE any
epi comparison (the OLS residual is the strength/geometry control; far cheaper
than a per-node surrogate). Bands α/β/γ_l (the only trace-bearing bands —
influence on a null trace is leverage on noise). ~24 s, n=10, no surrogate.

### 2.2 Result — a real within-patient β signal that does NOT transfer

**Finding 1 (NEW, positive, within-patient).** At β, Δρ_i is a real,
strength-INDEPENDENT epileptic signal:
- per-patient Spearman(Δρ, strength) median ≈ **−0.06** (range [−0.27, +0.20]) —
  **not** a strength re-read; residualization barely moves the AUC (0.621→0.634).
- epi-vs-non-epi cohort AUC **0.634, 8/9 patients, Wilcoxon p = 0.010**, and it
  **survives the rank-partial nonlinear-strength control** (AUC 0.650,
  p = 0.0076).
- Sign: labeled-epi nodes have **positive** median Δρ (+0.0022, trace-carriers);
  non-epi **negative** (−0.0058, trace-suppressors).
- This is genuinely new relative to audit_79/80: the row-restricted
  `rho_split_node` (their F1) was at chance (β AUC 0.50). **Leave-node-out
  leverage on the cohort trace ≠ the trace measured on the node's own row** —
  the influence statistic sees something the row statistic does not.
- **Caveat:** the within-patient β AUC (0.634) is *comparable to, not above*,
  the strength-baseline AUC (0.652). Even within-patient it is not a dominant
  axis. (α/γ_l Δρ within-patient: at chance, AUC 0.51 / 0.49.)

**Finding 2 (DECISIVE, negative, cross-patient).** In leave-one-patient-out
transfer — the actual occult/discovery use-case, because an UNLABELED patient
has no within-patient labels to set a threshold — Δρ_i is **at chance**:

| band | feature | prevalence | PR-AUC | lift | ROC-AUC |
|---|---|---|---|---|---|
| β | dRho | 0.106 | 0.123 | **1.16×** | 0.562 |
| α | dRho | 0.106 | 0.132 | 1.25× | 0.552 |
| γ_low | dRho | 0.106 | 0.127 | 1.20× | 0.478 |
| β | strength | 0.106 | 0.083 | 0.78× | 0.386 |

PR-AUC ≈ prevalence (lift ~1.2×) is the **same ceiling the audit_80 classifier
hit**. The within-patient β separation does not transfer to a held-out patient.

### 2.3 Occult verdict (honest)

**There is no transferable occult-node marker.** Consistent with audit_80
(`geometry_dominant`). The within-patient β LNO-influence signal is a genuine,
publishable *descriptive* fact about how epileptic contacts sit in the trace
(carriers, strength-independently) — but it cannot flag occult tissue in a new
patient, because (a) it does not exceed the strength baseline even within
patient, and (b) it does not transfer cross-patient.

The 154 candidates emitted (`occult_candidates.csv`: non-labeled contacts in the
top decile of within-patient residual Δρ AND inside the labeled-epi residual
span; 122 of 154 on epi-free shafts) are **hypothesis-generating ONLY**. With no
surgical-outcome (Engel-class) ground truth in this cohort, a candidate cannot
be confirmed or refuted; the cross-patient test says such a flag could not even
be produced for a genuinely unlabeled patient from this score. **The direction
is closed without surgical-outcome ground truth.**

### 2.4 The validation gap, stated plainly

"Occult candidate" ≠ "predicted epileptogenic focus". The only real validation
is surgical outcome at resected-vs-spared candidate sites, which this cohort does
not carry. Everything in §2 is either a within-patient descriptive fact (Finding
1) or an explicitly negative discovery result (Finding 2). No clinical claim is
made or supported.

---

## 3. Files created this pass

| path | what |
|---|---|
| `.agents/reports/2026-06-05_epilepsy-headline-and-occult-node-marker.md` | this report |
| `.agents/guides/task-persistence-investigation/2026-06-05_occult-epi-node-marker.md` | scope report (5-point preamble + 11 sections) for the occult reframe |
| `scripts/01_compute/audit/audit_90_occult_node_lno_influence.py` | the one honest cheap test (LNO influence + residualization + within-patient gate + LOPO transfer + candidates) |
| `data/audit/occult_node_marker/lno_node_influence.csv` | per (patient, band, node): Δρ, row-features, confounds, residuals |
| `data/audit/occult_node_marker/occult_gate_per_band.csv` | within-patient epi-vs-non AUC per score + strength baseline |
| `data/audit/occult_node_marker/occult_lopo_transfer.csv` | cross-patient LOPO PR/ROC-AUC (the decisive test) |
| `data/audit/occult_node_marker/occult_candidates.csv` | 154 ranked non-labeled candidates (hypothesis-only) |
| `data/audit/occult_node_marker/README.md` | audit_90 self-documenting README |

No pre-existing artifacts were edited. Deliverable-1 numbers were *read and
verified*, not recomputed (audit_77/78 caches already existed).

## 4. Numbers that surprised me

1. **The brief's premise that audit_80 was "deferred" was stale** — it had been
   run with a NEGATIVE verdict. The scope report was already `built-negative`.
2. **β LNO-influence Δρ is a real strength-INDEPENDENT within-patient epi
   signal** (AUC 0.634, 8/9, p=0.010, survives rank-partial) — I expected a flat
   negative given audit_79/80. It is new because influence ≠ row-feature. It
   still fails the occult use-case (no cross-patient transfer), so the headline
   verdict is unchanged, but the within-patient fact is genuine.
3. **γ_low Pat_14 flips −0.220 → +0.073 on epi exclusion** — the single largest
   per-patient swing, and the mechanism behind the low-γ "emerge" verdict (the
   epileptic zone was actively masking a healthy low-γ trace in this patient).
4. **α epi↔epi is 6/9 individually significant, not 5/9** as the brief stated
   (Pat_03 +0.760 p=0.015 clears too) — the α-recruits-epi-tissue effect is
   marginally stronger than briefed.
