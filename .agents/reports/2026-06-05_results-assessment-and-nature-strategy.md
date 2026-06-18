---
name: results-assessment-and-nature-strategy
type: report
era: IMCOH_ABS_COHORT_N10
status: current
created: 2026-06-05
updated: 2026-06-05
pointers:
  - .agents/preprint/locked/VERDICT_LEDGER.md
  - .agents/preprint/locked/ANATOMY_LEDGER.md
  - .agents/preprint/responses/2026-06-01_internal-referee-review.md
  - data/audit/epi_stratified/README.md
  - data/audit/anatomy_localization_wilcoxon/README.md
  - data/audit/per_patient_localization/README.md
  - .agents/reports/2026-06-05_epilepsy-headline-and-occult-node-marker.md
  - .agents/reports/2026-06-05_rhocoph-grassmann-dissociation.md
  - .agents/reports/2026-06-05_eigenmode-localization.md
---

# Results assessment + Nature-Neuroscience presentation strategy

## Head

The science is in much better shape than the preprint. There is **one genuinely
strong, fully-defended finding** — a β-band hierarchical task-trace that survives
the mandatory matched-strength null on **two** independent probes — **one solid
partner** (α, cophenetic only), and **two reframings that are currently buried but
are the best Nature-style material we have**: the trace is **spatially distributed**
(not focal) and it **lives in clinically-healthy tissue** (it *strengthens* when the
epileptic zone is removed). The preprint reads as weak not because the science is
weak but because it is **mid-retraction-cascade and method-first**. The single real
cap on the journal ceiling is that **there is no behavioral or cognitive anchor**
anywhere in the project (the task paradigm itself is undocumented) — that is what
separates "strong specialty paper" from "Nature Neuroscience."

This report is the honest inventory + the presentation strategy. Decisions locked
today (2026-06-05) are in §7.

---

## 1. Strength schematic — what we actually hold

Verdict vocabulary: **SURVIVES** = clears the matched-strength surrogate (the only
bar that counts). **EMERGES** = null on the full graph, clears once the epileptic
zone is excluded. **REFRAME** = a retraction turned into a positive message.
**WEAK / NULL** = honest.

| Finding | Probe(s) | Key number | Tier | The one caveat |
|---|---|---|---|---|
| **β hierarchical trace** | ρ^coph **+** Grassmann | coph p=0.005 / 23.7× / 7-of-10; Grassmann p_mass=0.005, k=27–55 (29 cells), LOO-robust | **★ HEADLINE** | drift-control margin is the thinnest (C2 p=0.014) |
| **α hierarchical trace** | ρ^coph only | p=0.002 / 8.3× / 5-of-10 | **★ Secondary** | lower patient consensus; no global-mode (Grassmann) trace |
| **Band-selectivity** (θ clean null) | both | θ coph p=0.72 (anti); Grassmann p=0.14 | **★ Supporting** | it is the *control* — present it as one |
| **Trace is DISTRIBUTED** | both | 0 DK regions >5/10 pts; per-pt 0–1/10 every band | **★ REFRAME → message** | confounded with sparse sEEG sampling (§6) |
| **Trace is EPILEPSY-INDEPENDENT** | ρ^coph | β +0.221→+0.275 (p=0.003); α +0.105→+0.187; γ_l EMERGES p=0.024 | **★ REFRAME → message** | epi-only subgraph underpowered; global-rewiring null |
| **LRG-emergent** | — | raw FC fails MS gate (β p=0.053); only the hierarchy clears | **Supporting argument** | — |
| α-recruits-epi / β-spares-epi | ρ^coph | α epi-core +0.404 p=0.006, 5/9 pts sig | **Secondary obs.** | 45-pair subgraph, small |
| γ_low | Grassmann (full) + ρ^coph (epi-X) | Grassmann p_mass=0.005, k=12–23; coph EMERGES under epi-X | **Secondary** | Grassmann LOO-fragile under epi-X (0.159) |
| δ | Grassmann only | p_mass=0.005 but LOO 0.055 (Pat_08) | **Honest weak** | single-patient-leveraged at full data |
| γ_high | both | p≈0.055–0.060 | **Null (near-miss)** | one permutation from flipping |
| θ | both | clean negative | **Null (the benchmark)** | — |
| **All anatomy / localization** | — | RETRACTED | **✗ Dead** | cascade still unfinished in brief frontmatter |

---

## 2. The strongest result, stated cleanly

**β (13–30 Hz) is the only band with the trace on both LRG probes, and both clear
matched-strength:**

- **ρ^coph (per-pair cophenetic — the backbone measure):** cohort paired Wilcoxon
  **p = 0.00488**, effect **23.7× surrogate**, median +0.221, 7/10 patients above
  their own surrogate. Passes all five controls (C1 split 0.005, C2 drift 0.014,
  C3 matched-strength 0.005, C4 cross-probe non-degradation, C5 epi-X 0.003).
- **Grassmann (leading-subspace chordal — the spectral cross-check):**
  cluster-mass **p = 0.005**, all-clusters mass T_G* = 69.76, a 29-cell contiguous
  significant run at k = 27–55, fully LOO-robust (drop Pat_02 → 0.005), strengthens
  under epi-exclusion (mass → 89.04).

**α (8–13 Hz)** is the clean secondary: ρ^coph **p = 0.00195**, 8.3×, median +0.105
— cophenetic only (no Grassmann trace). **θ is the clean negative** (coph p=0.72 in
the *anti* direction) and is the band-specificity benchmark that kills "any cognitive
engagement perturbs all connectivity."

Everything else is honestly weak and belongs in a dissociation table, not as
findings: γ_low and δ survive on Grassmann only (δ single-patient-leveraged);
γ_high is a near-miss null; raw FC does not clear the drift floor in any band.

---

## 3. The LRG-emergent point (sell it)

β is the **weakest band in raw connectivity** — raw FC fails the matched-strength
gate at β (p = 0.053) and fails the drift floor in *every* band (best is γ_low at
p ≈ 0.057). The trace becomes significant **only after** each matrix is summarized
as a multiscale hierarchy. That is the methodological justification, handed to us:
the LRG step is *what makes the signal detectable*, not decoration. Frame the
persistence as living in *how the network reorganizes hierarchically*, not in which
edges changed.

---

## 4. The two reframings (currently buried — promote them)

### 4A. Distributed, not focal
Every anatomical localization claim was correctly retracted: no Desikan-Killiany
region is sampled by more than 5/10 patients; per-patient localization is 0–1/10 in
every band; the apparent "hotspots" are electrode-shaft spatial autocorrelation
(region↔shaft NMI 0.65), not anatomy. The **positive reading**: the β/α trace is a
**whole-network topological reorganization**, not edge-strengthening at one site —
and at β, 8/10 patients are 61.7% sign-consistent across the *entire* pairwise graph.
(Honest caveat in §6: delocalization is partly forced by sparse sEEG coverage.)

### 4B. Epilepsy-independent (the strongest clinical sentence we have)
When epileptic contacts are excluded, the trace **strengthens**: β +0.221 → +0.275
(p = 0.003), α +0.105 → +0.187, and the epileptic-core subgraph *weakens*. **This
pre-empts the first reviewer question** ("isn't this just pathological-network
reorganization?") and converts the cohort's central liability (epilepsy patients)
into a control: if the effect were an artefact of pathological hypersynchrony or
interictal spiking, removing the epileptic zone would weaken it — instead it cleans
up. The trace lives in healthy circuitry and its interface with the epileptic zone,
and the epileptic zone partly *masks* it. **This was sitting in YAML frontmatter as
a robustness footnote. As of today it is a primary message** (see §7).

---

## 5. Numbers we under-noticed

- **α recruits the epileptic tissue; β spares it.** Within the epileptic subgraph,
  α shows a *strong* trace (+0.404, p = 0.006, 5/9 patients individually significant)
  while β weakens there. A clean, two-method-corroborated band dissociation about how
  different rhythms engage pathological vs healthy circuitry.
- **Two patients carry the β trace purely in subspace rotation.** Pat_10 and Pat_14
  are flat/anti on ρ^coph but strongly positive on Grassmann (z ≈ +5). And `onlyC = 0`
  in every band — every cophenetic-positive patient is also Grassmann-positive. The
  two probes are **nested in sensitivity yet dissociate per-band** (Grassmann 9/10 vs
  coph 7/10 at β).
- **γ_low emerges on epi-exclusion** — full-graph p = 0.116 (null) → exclude-epi
  p = 0.024 (8/10), with Pat_13 (−0.042→+0.109) and Pat_14 (−0.220→+0.073) *flipping
  sign* when epi contacts are dropped. The epileptic zone was actively *inverting* a
  healthy-tissue γ_low trace. This contradicts the old γ_low brief sentence ("no
  per-pair cophenetic trace") — updated today.

---

## 6. The honest ceiling (what a skeptical Nat Neuro reviewer hits first)

1. **No cognition, no behavior, no mechanism.** The paper claims a *cognitive* task
   leaves a trace but cannot say what the task was, what was learned/tested, or
   whether trace strength relates to anything the brain accomplished. No accuracy/RT/
   recall anywhere → no brain–behavior correlation. This demotes it from a
   cognitive-neuroscience finding to a **connectivity-dynamics** observation.
   **Highest-value thing to chase: recover the task paradigm + any behavioral measure.**
2. **n = 10, heterogeneous implants — delocalization confounded with sampling.** No
   DK region is sampled by >5/10 patients, so cohort anatomy is *untestable by
   construction*. "Distributed" is partly a sampling statement. Blunt it with the
   within-patient whole-graph sign-consistency (which needs no shared anatomy).
3. **A modest, abstract structural-distance effect.** The headline is a Spearman
   correlation of merge-order changes (median ρ ≈ +0.22 β / +0.11 α) on a derived
   ultrametric object — not a power change or coupling strength. β (both probes, 23.7×)
   is genuinely solid; the weak bands invite "12 cells tested, 2 survive cleanly."
4. **LRG is non-standard and runs in an admitted outlier regime** (fully-connected,
   continuous-spectrum; standard scale-identification fails). Matched-strength guards
   it, but it is an attack surface for a methods reviewer.

**One-line verdict:** the β-band hierarchical post-task trace is real, well-controlled,
epilepsy-independent, and methodologically novel — but it is a modest, abstract
structural-distance effect in a small clinically-constrained cohort with no behavioral
anchor. As-is, it is a strong **specialty-journal** paper (eLife / NeuroImage / Network
Neuroscience / Brain). For Nat Neuro it is two ingredients short: a behavioral correlate,
and evidence that "distributed" is biology not sampling.

---

## 7. Decisions locked today (2026-06-05, PI directive)

1. **Backbone = ρ^coph.** The per-pair cophenetic trace is the **primary, novel,
   multiscale-aware** measure the paper introduces. **Grassmann is the secondary
   spectral cross-check** ("usual spectral analysis"), valuable mainly to test for a
   **dissociation** from ρ^coph.
2. **The α dissociation is a finding.** α has a ρ^coph trace but **no** Grassmann
   trace = the task reorganizes *fine relational structure* that persists, **without**
   reorganizing the *dominant collective modes*. (Direction: characterize across all
   bands + figure — in progress.)
3. **Epi-exclusion is elevated from a secondary footnote to a primary interpretive
   lens.** The headline is not just "the trace survives epi-exclusion" but "the trace
   is a property of healthy network reorganization that the epileptic zone *masks*."
   We are *also* investigating whether these measures can flag **occult (non-labeled)
   epileptic nodes** — honest prior is negative (audit_79 per-node gate failed vs the
   strength baseline); in progress.
4. **γ_low verdict updated** to surface the epi-X cophenetic emergence (audit_77):
   full graph null → exclude-epi p = 0.024 (8/10). Ledger + brief edited today.

## 8. Three directions — outcomes (completed 2026-06-05)

**1. Epilepsy headline — STRONG POSITIVE.** Consolidated, all numbers verified
against CSVs. β +0.221→+0.275 (p=0.003, 7/10, LO-Pat_13 0.006); α +0.105→+0.187
(p=0.014); low-γ EMERGES (p=0.024, 8/10). Pair-class decomposition: carried by the
**cross** (interface, β +0.273 p=0.010) + **healthy↔healthy** (β +0.197 p=0.005);
the **epi↔epi block is non-significant** (β p=0.150). α exception: α epi↔epi +0.404
(p=0.006, **6/9** patients individually sig). → `2026-06-05_epilepsy-headline-and-occult-node-marker.md`.
*Occult-node marker — NEGATIVE for transfer.* Cross-patient LOPO at chance
(ROC 0.562). No transferable marker. New within-patient nuance: β leave-node-out
influence is strength-independent and separates epi from non-epi (AUC 0.634, p=0.010,
8/9) — the trace carries epi-relevant per-node structure that does not generalize.
154 candidate nodes = hypotheses only (no surgical ground truth).

**2. ρ^coph × Grassmann dissociation — CLEAN; validates the backbone decision.**
4-cell map: **BOTH** (β), **COPH-ONLY** (α — flagship local-but-no-global-modes),
**GRASSMANN-ONLY** (δ, γ_l full-graph), **NEITHER** (θ, γ_h). γ_l crosses into BOTH
under epi-exclusion (cophenet emerges, Grassmann contracts — opposite directions).
Honest caveat: per-patient nesting (`onlyC=0`) is partly a gating artefact (per-pair
Wilcoxon vs curve-level cluster-mass FWE); the **cohort dissociation is the defensible
statement**. → `2026-06-05_rhocoph-grassmann-dissociation.md` + dissociation-map figure
(`preprint_34`).

**3. Eigenvector localization — CLEAN DOUBLE NEGATIVE.** Q1 (do localized non-leading
modes explain the α dissociation?) FALSIFIED 10/10: leading modes are the *most*
localized, not extended; the trace rides extended modes, **but the within-baseline
control shows the same** → a generic FC-reproducibility property, not a trace channel.
**The α-dissociation mechanism remains OPEN** (the obvious explanation is ruled out).
Q2 (epi eigenmode trapping, Anderson analogy) CLEAN NEGATIVE — consistent with the
spatial-delocalization verdict. One generic structural fact: FC-Laplacian modes are
2–10× more localized than the strength null (all 240 cells, strongest γ) — real but
not a trace/epilepsy signature. Promoted library helpers (`participation_number`,
`node_set_mode_mass`, `subspace_displacement`). → `2026-06-05_eigenmode-localization.md`.

**Concurrency note (2026-06-05).** Scripts created outside these three agents appeared
in the same window — `audit_80_epi_marker_classifier` (18:21, predates my agents;
built-NEGATIVE, geometry-dominant), `audit_81_localization_atlas` (re-examines the
2026-05-30 anatomy retraction at coarser system/lobe granularity conditioned on implant
sampling — *confirms* the region-level retraction, offers a descriptive coarser reading),
`audit_81_epi_node_leverage`. A parallel session is actively climbing audit numbers
(it took 82–86 during this work). Per PI directive, my two scripts were moved CLEAR of
that range — `audit_90_occult_node_lno_influence` and `audit_91_eigenmode_localization`
— leaving the 80s block to the parallel session (their files untouched). The parallel
`localization_atlas` / `epi_propagator` work is not integrated or verified here.
