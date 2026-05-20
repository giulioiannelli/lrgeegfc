---
name: writing-directive-methods-beta-latex-verification
era: IMCOH_ABS_COHORT_N10
status: directive
kind: writing-agent-directive
date: 2026-05-19
target: manuscript Methods (§ssec:methods_lrg + §ssec:methods_compare) + β results (§ssec:results_beta)
source_of_truth:
  - .agents/preprint/locked/CONTROLS.md (Decisions 8 + 9 + 10 + 11, locked 2026-05-19 pm)
  - .agents/preprint/locked/VERDICT_LEDGER.md (post-Decision 8 verdicts)
  - .agents/preprint/locked/ANATOMY_LEDGER.md
  - .agents/preprint/methods/methods_revision_2026-05-18_cophenet.md
  - .agents/preprint/methods/methods_grassmann_cluster_extent.md
  - .agents/preprint/methods/methods_section_review_2026-05-19.md
  - .agents/preprint/directives/writing_directive_2026-05-19_beta_post_methods_revision.md
  - .agents/preprint/directives/methods_directive_2026-05-19_TG_normalization.md
  - .agents/preprint/responses/feedback_response_2026-05-19_gate_refactor.md
  - .agents/preprint/bands/01_beta.md
supersedes: none (this directive is a cumulative verification pass on top of prior writing/methods directives)
---

# Methods + β LaTeX verification pass — what to patch (2026-05-19)

## Head

After full cross-check of the supplied Methods (§ssec:methods_lrg +
§ssec:methods_compare with four subsubsections) and β results
(§ssec:results_beta — six paragraphs + per-patient table) against
the locked Decisions 8/9/10/11, the locked verdict and anatomy
ledgers, and the prior writing/methods directives, the manuscript is
substantially aligned. **Three must-fix items in Methods** and
**three in β results** remain. They are all single-paragraph or
single-sentence patches and do not require text reflows. Two open
data items (the `T_G^*` numerical-value discrepancy and the BH
family-size choice for matrix distances on `D_coph`) need a user
decision before the writing agent re-touches the affected
sentences. The β subsection's structure is sound; once the patches
land, β is methods-aligned, the per-patient table is locked, and
the same structure transposes cleanly to α / γ_l / θ / γ_h / δ via
the per-band briefs already in `bands/02_alpha.md` …
`bands/06_delta.md`.

## Status snapshot (updated 2026-05-19 evening after user decisions)

| Section | Status | Must-fix items | Editorial review items | Open data items |
|---|---|---|---|---|
| Methods §ssec:methods_lrg | aligned | 0 | 0 | 0 |
| Methods §sssec:methods_compare_rawfc | aligned | 0 | 0 | 0 |
| Methods §sssec:methods_compare_ctm | partial (C4 gate stale) | 1 (M1: C4) | 0 | 0 (D2 RESOLVED) |
| Methods §sssec:methods_compare_grassmann | aligned | 0 | 1 (S4: audit citations) | 0 (D1 RESOLVED) |
| Methods §sssec:methods_compare_stats | partial (C5 + LOO stale) | 2 (M2: C5, M3: LOO) | 1 (S4: audit citations) | 0 |
| β results head | partial (γ_l + δ vocabulary stale) | 1 (R1) | 0 | 0 |
| β results per-pair paragraph | aligned | 0 | 0 | 0 |
| β results per-pair table | partial (bold criterion hardcoded) | 1 (S1: null-based bold) | 0 (S2 RESOLVED — keep substrate) | 0 |
| β results per-pair interpretation | aligned | 1 (R3: "order of magnitude" leftover) | 0 | 0 |
| β results Grassmann paragraph | aligned | 0 | 1 (S3: epi-X specificity) | 0 |
| β results anatomy paragraph | aligned | 1 (R2: medial-temporal phrasing) | 0 | 0 |
| β results synthesis | aligned | 0 | 0 | 0 |
| Cross-band synthesis (forward-looking) | framing correction | 0 (writing agent guidance in Part VI) | 1 (cascading edits to `bands/00_cohort.md`) | 0 |

**Total (post user-decision update)**: 7 must-fix patches (3 methods
+ 4 β results — including S1 promoted from SHOULD to MUST per user
directive), 3 editorial review items (S3, S4, Part VI cascading), 0
open data items.

---

## User decisions (2026-05-19 evening — locked)

The user resolved the four pending items and corrected one framing
issue:

- **D1 → (A) all-clusters formula stands.** `T_G^*(b) = Σ_{k : p_k <
  α_k} (−log10 p_k)` over **all** significant `k`-cells (not the
  longest contiguous cluster). Locked normalized values: **β =
  0.273** (raw 69.76), **γ_l = 0.259** (raw 66.14), **δ = 0.149**
  (raw 38.07). The Methods Eq.~\eqref{eq:methods_TGstar} in the
  LaTeX already cites the all-clusters sum; no methods patch
  needed. **Pending action (user / non-writing-agent)**: regenerate
  `data/audit/grassmann_cluster_extent/cohort_summary.csv` by re-
  running the current (corrected) audit_70 script so the stored
  `obs_cluster_mass_neglog10p` value for β matches the methods
  equation (69.76 instead of the stored 52.97). After the re-run,
  every cohort-row reference to β `T_G^*` becomes the all-clusters
  value, and the `bands/01_beta.md` brief + `methods_grassmann_cluster_extent.md`
  §6 are already aligned with this.

- **D2 → (A) m = 6 for `d_S` alone stands.** `d_P` is permanently
  dropped from matrix distances on `D_coph`. The current LaTeX
  §sssec:methods_compare_stats already has m = 6 for `d_S` on
  `D_coph` and is correct as written. **Cascading action (non-
  writing-agent)**: update `methods/methods_revision_2026-05-18_cophenet.md`
  Statistical Inference subsection to align m = 6 (currently says
  m = 12 with `d_S` + `d_P`).

- **S1 → bold criterion changes from hardcoded to null-based.** The
  per-patient `n_sig,k ≥ 20` heuristic is **dropped** (no hardcoded
  per-patient threshold). Per-patient significance on `T_G^{*,p}`
  becomes an **empirical p-value against the patient's own null
  distribution of `T_G^{*,p}`**, in parallel with the cohort-level
  cluster-mass permutation gate of Eq.~\eqref{eq:methods_cluster_p}.
  The patches are now in **Part III S1 (promoted to MUST-FIX)**.
  **Pending action (user / non-writing-agent)**: extend audit_70 to
  also emit per-patient `T_G^{*,p,null}` distributions (R = 200
  phantom-surrogate values per patient via per-patient phantom
  permutation) and per-patient empirical `p_mass^p`, written to
  `per_patient_per_band_per_k.csv` or a sibling CSV. The bold
  cells flow from the new column.

- **S2 → keep substrate column in the per-patient table.** R6's
  "drop substrate column" directive is **superseded by user
  decision**; the substrate column (col 8 of `tab:beta_per_patient`)
  is retained as a comparison reference. The current LaTeX caveat
  "does not gate LRG-trace classification" + the trace-direction
  sign-flip per R6 stand. No LaTeX patch needed for S2.

**Cross-band synthesis framing correction (user 2026-05-19
evening).** The "three-layer cohort table" terminology was being used
in two distinct senses across the briefs, and the user has clarified
which is the cross-band synthesis and which is a sensitivity
argument:

- **The three analytical layers** of the LRG pipeline are
  **substrate (raw FC) → per-pair multiscale (D_coph) → subspace
  (Grassmann)**. These are the three layers of the **cross-band
  cohort synthesis**: each band's verdict reads at each of the three
  layers (β positive at all; α at substrate + D_coph; γ_l/δ at
  substrate + Grassmann; γ_h/θ at substrate only).
- **The raw FC / raw `D(τ_max)` / cophenet `D_coph`** three-layer
  table from `bands/00_cohort.md` §2 and `bands/01_beta.md`
  "Headline three-layer cohort table" is a **sensitivity argument
  internal to the D_coph probe's adoption** — it shows that raw FC
  and raw `D(τ_max)` are operationally indistinguishable at
  matched-strength gating, and that the cophenet wrap is what
  delivers band-resolution at the per-pair multiscale layer. It is
  **not** the cross-band synthesis; it is **methodology for adopting
  D_coph** at one specific layer.

Part VI of this directive is updated to reflect this distinction.
**Cascading action (non-writing-agent)**: `bands/00_cohort.md` §2
should be reframed to label the raw FC / raw D / cophenet table as
"Sensitivity check: cophenet vs raw D(τ_max) at the D_coph layer"
rather than as "Three-layer cohort matched-strength gating", and a
new cross-band synthesis section organized around the three
analytical layers (substrate / D_coph / Grassmann) should be added.
The writing agent's cross-band synthesis paragraph should follow the
analytical-layer organization, not the sensitivity-argument table.

**Updated totals**: 7 must-fix patches (3 methods + 4 β results),
3 editorial review items (S3 + S4 + Part VI cascading edits), 0
open data items remaining (D1 + D2 resolved; both have downstream
user/non-writing-agent action items but no further LaTeX patches
hinge on them).

---

# Part I — Methods section: must-fix patches

## M1. C4 cross-probe gate (§sssec:methods_compare_ctm, last paragraph)

**Issue.** The closing paragraph of §sssec:methods_compare_ctm currently states:

> "The cohort claim at the per-pair multiscale layer requires \(\rhosplit\) to sit above \(\rhodrift\) under a one-sided paired Wilcoxon, with \(\rhoxprobe\) preserving the cohort-median sign and **patient count**."

The "patient count" clause is the **retired** C4 gate. Decision 9
(locked 2026-05-19 pm, `locked/CONTROLS.md` §C4) retired patient-count
thresholds (`n_trace_xprobe ≥ 6/10`) in favor of a paired Wilcoxon
non-degradation test plus sign-agreement, on `feedback_no_hardcoded_test_thresholds.md`
+ `feedback_patient_counts_never_the_gate.md` grounds. The methods
section must cite the locked gate.

**Source of truth.** `locked/CONTROLS.md` §C4 Decision 9.

**Replacement text** (full sentence, in place):

> "The cohort claim at the per-pair multiscale layer requires \(\rhosplit\) to sit above \(\rhodrift\) under a one-sided paired Wilcoxon, and the cross-probe restriction \(\rhoxprobe\) to satisfy a paired one-sided non-degradation test: the paired Wilcoxon on \((\rhosplit - \rhoxprobe)\) under \(H_{1}: \rhosplit > \rhoxprobe\) must fail to reject at \(\alpha = 0.05\) (no significant degradation under cross-probe restriction), with \(\mathrm{sign}(\rhoxprobe^{\rm med}) = \mathrm{sign}(\rhosplit^{\rm med})\)."

**Severity**: MUST. Decision 9 is locked; the LaTeX as written cites a retired gate.

---

## M2. C5 epi-X gate (§sssec:methods_compare_stats, last paragraph)

**Issue.** The closing paragraph of §sssec:methods_compare_stats currently states:

> "A single systematic robustness regime is applied throughout: an epileptogenic-zone exclusion that drops contacts in each patient's clinically identified epileptogenic zone and recomputes each probe on the residual contact set. The cohort claim at any band is reported as robust only when the directional sign of the per-patient triangle scalar is preserved under the regime **and the one-sided Wilcoxon \(p\)-value remains below 0.05**."

This is the **superseded** C5 form (≳80% retention + ambiguous
"Wilcoxon"). Decision 10 (locked 2026-05-19 pm, `locked/CONTROLS.md`
§C5) replaced this with two probe-specific Wilcoxon-on-epi-X gates:

- **Grassmann (all bands)**: re-run the cluster-mass permutation test on the epi-X surrogate ensemble, gate at `cluster_p_mass^epi-X(b) < 0.05`.
- **Cophenet (α only)**: one-sample one-sided Wilcoxon on per-patient `ρ_split^epi-X` under `H_1: ρ_split^epi-X > 0`, gate at `wilcoxon_one_sided_p < 0.05`.

The current LaTeX collapses these two distinct gates into a single
"one-sided Wilcoxon" sentence, which is technically wrong for the
Grassmann layer (the Grassmann C5 statistic is the cluster-mass
permutation `p_mass^epi-X`, not a Wilcoxon).

**Source of truth.** `locked/CONTROLS.md` §C5 Decision 10;
`methods/methods_section_review_2026-05-19.md` M4 RESOLVED entry.

**Replacement text** (full paragraph, in place):

> "A single systematic robustness regime is applied throughout: an epileptogenic-zone exclusion that drops contacts in each patient's clinically identified epileptogenic zone per~\cite{xxx?} and recomputes each probe on the residual contact set. The robustness gate is probe-specific. At the Grassmann layer, the cluster-mass permutation test of \eqref{eq:methods_TGstar}--\eqref{eq:methods_cluster_p} is re-run on the epi-zone-excluded surrogate ensemble, and the band is reported as robust when \(p_{\rm mass}^{\rm epi-X}(b) < 0.05\). At the per-pair cophenetic layer the corresponding test is a one-sample one-sided Wilcoxon on the per-patient \(\rhosplit^{\rm epi-X}\) under \(H_{1}: \rhosplit^{\rm epi-X} > 0\), gated at \(p < 0.05\); this is reported for \(\alpha\), where it is the band that the cophenet probe identifies as a non-epileptic-cortex trace partly masked by epi-zone contacts in the full-cohort analysis. For every other band the epi-X analysis at the cophenetic layer is not run, and the verdict from controls 1--4 stands."

**Severity**: MUST. Decision 10 is locked; the current LaTeX form
does not correctly describe the Grassmann C5 gate (which is a
cluster-mass permutation, not a Wilcoxon).

---

## M3. LOO max-p sensitivity diagnostic (§sssec:methods_compare_stats, new paragraph)

**Issue.** Methods section currently does not state the LOO max-p
diagnostic that is mandatory under Decision 11
(`locked/CONTROLS.md`). The diagnostic is descriptive — never gates
— but flags single-patient-leveraged verdicts (the canonical
worked example is δ Grassmann full-data, where LOO max
`p_mass = 0.055` crosses the gate when Pat_08 is dropped; the C5
epi-X analysis subsequently resolves this and the full-data
leverage is attributed to epi-zone interactions rather than the
true biological trace).

This must be stated in Methods so the Results text's δ-Grassmann
LOO caveat (in the δ band paragraph, when that paragraph is
written) has a methods anchor.

**Source of truth.** `locked/CONTROLS.md` Decision 11;
`feedback_no_single_patient_p_driven.md`.

**Replacement** (new paragraph, appended at the end of §sssec:methods_compare_stats after the epi-X paragraph):

> "Every cohort-level Wilcoxon-based and cluster-mass-permutation-based gate is paired with a leave-one-out (LOO) sensitivity diagnostic: the cohort test is recomputed \(n = 10\) times, each excluding one patient, and the maximum resulting \(p\)-value (LOO max \(p\)) is reported alongside the cohort \(p\) and the patient whose exclusion drives it. LOO max \(p\) is descriptive only and never gates a verdict; verdicts where LOO max \(p\) crosses 0.05 are flagged in the Results text as single-patient-leveraged. The pairing is reported for every gate in the four trace controls (C1 split, C2 drift, C3 matched-strength on both probes, C4 cross-probe) and the C5 epi-X gate."

**Severity**: MUST. Decision 11 mandates the LOO diagnostic; without
it, the δ Grassmann LOO caveat in the δ band Results paragraph (when
written) has nowhere to anchor.

---

# Part II — β results: must-fix patches

## R1. Verdict vocabulary for γ_l + δ in head paragraph

**Issue.** The β head paragraph (paragraph 1 of §ssec:results_beta)
currently states:

> "the remaining bands carrying at most a partial signature: \(\alpha\) at the per-pair level only, \(\gamma_{\mathrm{low}}\) and \(\delta\) at the subspace level only and **at a weak cluster-extent gate**, and \(\gamma_{\mathrm{high}}\) and \(\theta\) none"

Post-Decision 8 (mass-only Grassmann gate, locked 2026-05-19 pm),
γ_l and δ both have `cluster_p_mass = 0.005`, which clears the
strong-trace threshold of `cluster_p_mass < 0.01` in
`locked/CONTROLS.md` §C3. The locked
verdicts in `locked/VERDICT_LEDGER.md` are **strong trace, only
Grassmann** for both bands — the upgrade from weak → strong on the
mass-only gate is explicit in the ledger's 2026-05-19 pm Decision 8
entry.

The phrase "at a weak cluster-extent gate" contradicts the locked
ledger. Drop it; the per-pair vs subspace distinction is the
substantive cross-band point, not the strong/weak gradient.

**Source of truth.** `locked/VERDICT_LEDGER.md` (γ_l + δ rows,
revised 2026-05-19 pm Decision 8); `locked/CONTROLS.md` §C3
Decision 8.

**Replacement text** (drop the "at a weak cluster-extent gate"
qualifier):

> "the remaining bands carrying at most a partial signature: \(\alpha\) at the per-pair level only, \(\gamma_{\mathrm{low}}\) and \(\delta\) at the subspace level only, and \(\gamma_{\mathrm{high}}\) and \(\theta\) none"

**Severity**: MUST. Locked verdicts are strong for γ_l and δ at the
Grassmann probe; the head paragraph must not contradict the ledger.

---

## R2. Anatomy phrasing — medial-temporal at both probes

**Issue.** The β anatomy paragraph (paragraph 5 of §ssec:results_beta)
currently states:

> "The two probes therefore read complementary anatomical sub-networks of the same \(\beta\) reorganization: **Hippocampus and parahippocampal cortex appear at both**, and the insula appears at both with opposite hemispheres (a Desikan--Killiany label match, not a contact-set match)."

The "Hippocampus and parahippocampal cortex appear at both" clause is
**misleading** as the per-DK-label cross-probe overlap is:

- Cophenet probe: **left parahippocampal cortex** (DK
  `ctx-lh-parahippocampal`); the Hippocampus does **not** appear.
- Grassmann probe: **Hippocampus** (DK `Hip`); parahippocampal cortex
  does **not** appear.

These are two distinct DK labels (one cortical, one subcortical) in
the same medial-temporal-lobe family. The sentence as written reads
to a careful neuroanatomy reviewer as "both regions appear at both
probes", which is false at the DK-label level.

**Source of truth.** `locked/ANATOMY_LEDGER.md` β cophenet vs β
Grassmann region lists.

**Replacement text** (replace the clause about medial-temporal
overlap):

> "The two probes therefore read complementary anatomical sub-networks of the same \(\beta\) reorganization: the medial-temporal-lobe family appears at both (left parahippocampal cortex on the cophenet probe; the Hippocampus on the subspace probe), and the insula appears at both with opposite hemispheres (a Desikan--Killiany label match, not a contact-set match)."

**Severity**: MUST. The current wording will be flagged by a careful
reviewer at neuroanatomy level. Fixing makes the same cross-probe
complementarity claim correctly.

---

## R3. "More than an order of magnitude" leftover in interpretation paragraph

**Issue.** The interpretation paragraph (paragraph 3 of
§ssec:results_beta) currently states:

> "The magnitude is modest in absolute terms but **more than an order of magnitude above the matched-strength noise floor** of the same statistic, which is the operative comparison: ..."

R5 of `writing_directive_2026-05-19_beta_post_methods_revision.md`
removed the "more than an order of magnitude" phrasing from the
per-pair results paragraph in favor of "≈ 23× ratio". The
interpretation paragraph (one paragraph later) was not patched in the
same pass and retained the leftover "more than an order of magnitude"
qualitative phrasing.

It now reads inconsistently with paragraph 2 (which uses "≈ 23×
ratio"). Pick one phrasing; the locked directive R5 prefers the
explicit ratio.

**Source of truth.** `writing_directive_2026-05-19_beta_post_methods_revision.md`
R5.

**Replacement text** (replace the phrase):

> "The magnitude is modest in absolute terms but ≈ 23× the matched-strength noise floor of the same statistic, which is the operative comparison: ..."

(Alternative qualitative form: "but well above the matched-strength
noise floor of the same statistic". Either works; the directive
prefers the explicit ratio for consistency with paragraph 2.)

**Severity**: SHOULD (consistency with R5; not a verdict-level
correction).

---

# Part III — Should-fix / editorial review

## S1. Per-patient bold criterion: drop hardcoded `n_sig,k ≥ 20`, use null-based gate **(promoted to MUST-FIX per user directive)**

**Issue.** The per-patient table caption currently defines:

> "Bold values in cols.\ 3 and 6 indicate per-patient significance against the matched-strength surrogate (\(z \ge 1.96\) on \(\rhosplit\); **\(n_{\mathrm{sig},k} \ge 20\) on \(T_G^{*,p}\), well above the null expectation of \(\approx 5.5\) cells**)."

The `z ≥ 1.96` rule for `ρ_split^coph` is null-based and is correct.
The `n_sig,k ≥ 20` rule for `T_G^{*,p}` is **hardcoded**, which the
user has retired per `feedback_no_hardcoded_test_thresholds.md`. The
per-patient `T_G^{*,p}` should instead use a **per-patient
empirical p-value against the patient's own null distribution of
`T_G^{*,p}`**, in direct parallel with the cohort-level cluster-mass
permutation gate of Methods Eq.~\eqref{eq:methods_cluster_p} and
with the `z ≥ 1.96` `ρ_split^coph` rule that the table already uses
for col. 3.

**Source of truth.** User directive 2026-05-19 evening (no
hardcoded thresholds on patient results; the null distribution we
already control against is what determines significance).

**Replacement (table caption, in place — replaces the bold-criterion
sentence)**:

> "Bold values in cols.\ 3 and 6 indicate per-patient significance against the patient's own \(R = 200\) matched-strength surrogate ensemble: \(z \ge 1.96\) on \(\rhosplit\), and empirical \(p\text{-value} < 0.05\) on \(T_G^{*,p}\) computed as \((1 + \#\{r : T_G^{*,p,\rm null}_{r} \ge T_G^{*,p,\rm obs}\}) / (R + 1)\), where \(T_G^{*,p,\rm null}_{r}\) is the per-patient phantom-surrogate cluster mass obtained by treating surrogate \(r\) as the observation against the remaining \(R - 1\) surrogates for that patient — the per-patient analogue of the cohort gate of Methods Eq.~\eqref{eq:methods_cluster_p}."

**Optional Methods addition** (one sentence at the end of
§sssec:methods_compare_grassmann after the cohort `p_mass`
discussion):

> "A per-patient analogue \(p_{\rm mass}^{p}(b) = (1 + \#\{r : T_G^{*,p,\rm null}_{r}(b) \ge T_G^{*,p}(b)\}) / (R + 1)\) is computed by phantom-surrogate permutation at the single-patient level and is reported in the per-patient table as a descriptive single-patient sensitivity diagnostic; the cohort gate at the Grassmann layer remains the band-level \(p_{\rm mass}(b) < 0.05\) on \(T_G^*\)."

**Implementation note (non-writing-agent action item)**: The audit_70
script needs a small extension to emit per-patient phantom-surrogate
cluster-mass null distributions and per-patient empirical `p_mass^p`
per band. Output should land in
`data/audit/grassmann_cluster_extent/per_patient_p_mass.csv` or as a
new column in `per_patient_per_band_per_k.csv`. Until this CSV is
produced, the bold cells in col. 6 of the table cannot be filled
deterministically. Expected practical outcome: 9/10 patients bold
(Pat_15 with `T_G^{*,p} = 0.006` will not clear; all others have
substantial per-patient cluster mass).

**Per-patient "LRG-trace classification" rule** (separate from the
bold-cell criterion): the rule "a patient is LRG-trace if at least
one of cols.\ 3 or 6 is bold" is retained — it reads the cohort
structure cleanly and the closing synthesis paragraph relies on it.
The classification is descriptive only and does not gate any cohort
verdict.

**Severity**: MUST (per user directive removing hardcoded
patient-level thresholds).

---

## S2. Substrate column kept in per-patient table — **RESOLVED**

**Status.** R6's "drop substrate column" directive is **superseded
by user decision 2026-05-19 evening**: keep col 8 (substrate trace
direction) as a comparison reference. No LaTeX patch needed; the
current `tab:beta_per_patient` 9-column layout stands.

The current caveat in the caption ("does not gate LRG-trace
classification") is correct framing and should remain. The Notes
column's "mixed" annotations on Pat_10 + Pat_14 (substrate-vs-LRG
disagreement) are descriptive context for the reader and should
remain.

**Severity**: RESOLVED. No action required.

---

## S3. Epi-X specificity in Grassmann paragraph

**Observation.** The Grassmann paragraph (paragraph 4) currently
states the epi-X result qualitatively only:

> "Removing the epileptic-zone contacts and rebuilding the analysis on the non-pathological tissue preserves the cluster-mass gate and strengthens the cohort signal, so the rotation of the leading communication modes is carried by physiological cortex."

Under R1 of
`writing_directive_2026-05-19_beta_post_methods_revision.md`
(drop `L_obs` references from results), the absence of "29 → 36
cells" or "k ∈ [21, 56]" is correct.

But the locked C5 cluster-mass numbers are still citable: the C5
epi-X cluster-mass clears the gate at `p_mass^epi-X = 0.005`, and
the **normalized cluster mass strengthens from `T_G^* = 0.273`
(full data) to `T_G^*^epi-X = 0.348`** (using the all-clusters
formula; from `bands/01_beta.md` β C5 entry: mass 69.76 → 89.04 raw,
or 0.273 → 0.348 normalized). These are the load-bearing C5
numbers under the locked battery and could anchor the qualitative
"strengthens" claim.

**Caveat**: D1 below (`T_G^*` discrepancy CSV vs Methods equation)
applies; until resolved, the writing agent should cite `p_mass^epi-X
= 0.005` only and skip the specific `T_G^*` strengthening number.

**Recommendation.** Add the `p_mass^epi-X = 0.005` citation:

> "Removing the epileptic-zone contacts and rebuilding the analysis on the non-pathological tissue preserves the cluster-mass gate at \(p_{\rm mass}^{\rm epi-X} = 0.005\) and strengthens the cohort signal, so the rotation of the leading communication modes is carried by physiological cortex."

**Source of truth.** `locked/VERDICT_LEDGER.md` β Grassmann C5 row;
`data/audit/grassmann_epi_exclusion/c5_wilcoxon_cohort.csv`.

**Severity**: SHOULD (R1-compliant; adds specificity without
introducing `L_obs`).

---

## S4. Cite audit_71 and audit_72 in Methods

**Observation.** Two new audit scripts produced the locked C4 +
C5 CSVs but the Methods section does not name them. Adding a
citation pattern (in line with the existing `audit_70` reference)
gives the reader a single source for the cluster-extent + cross-probe
+ epi-X audits.

**Recommendation.** Add one footnote or parenthetical citation per
audit when each gate is introduced:

- In §sssec:methods_compare_ctm at the C4 gate sentence (after M1
  patch): cite `audit_71` (`data/audit/ctm_triangle/c4_wilcoxon_cohort.csv`).
- In §sssec:methods_compare_stats at the C5 paragraph (after M2
  patch): cite `audit_72` (`data/audit/{grassmann_epi_exclusion,alpha_epi_exclusion}/c5_wilcoxon_cohort.csv`).

These citations match the cluster-extent citation pattern (the LaTeX
already implicitly cites `audit_70` via the equation labels). Keep
the citations to once per gate; do not litter the methods text.

**Severity**: OPTIONAL.

---

# Part IV — Open data items — **RESOLVED**

## D1. `T_G^*` formula — **RESOLVED: (A) all-clusters stands**

**User decision (2026-05-19 evening)**: (A) — methods equation
\eqref{eq:methods_TGstar} stands as the all-clusters sum over every
contiguous-significant `k`-cell. Locked normalized values: β = 0.273
(raw 69.76), γ_l = 0.259 (raw 66.14), δ = 0.149 (raw 38.07). The
`methods_grassmann_cluster_extent.md` §6 locked band-level table
already reflects this; the LaTeX Methods equation already reflects
this.

**Action item (non-writing-agent)**: regenerate
`data/audit/grassmann_cluster_extent/cohort_summary.csv` by re-running
the corrected audit_70 script (already in the corrected all-clusters
form per the script's `cluster_mass()` definition at lines 154–167).
After the re-run, `obs_cluster_mass_neglog10p` for β stores 69.76,
γ_l 66.14, δ 38.07. The β verdict `p_mass = 0.005` is unchanged
(both 52.97 and 69.76 sit at the empirical-null floor 1/(R+1) of the
R=200 ensemble). γ_l + δ verdicts are upgraded weak → strong under
the mass-only gate Decision 8 (already reflected in the locked
VERDICT_LEDGER post-2026-05-19 pm).

**LaTeX consequence**: writing agent may continue citing `p_mass =
0.005` only for β (no inline raw `T_G^*` value) until the
cohort_summary CSV is regenerated. After regeneration, the
normalized β `T_G^* = 0.273` becomes citable if the writing agent
wants to anchor the cohort-level cluster-mass scalar in the prose.

---

## D2. BH family m for matrix distances on `D_coph` — **SUPERSEDED 2026-05-20**

**Status:** RETIRED. The cross-band BH-FDR sentence is dropped in
full from methods (both `ρ_split^coph` `m = 6` and `d_S` `m = 6`
clauses). `d_S on D_coph` is no longer reported. See
`METHODS_AUDIT_ISSUES.md` §B2 (revised 2026-05-20) and
`writing_directive_2026-05-20_methods_audit_application.md` §B2 for
the broader policy + replacement sentence. Three-point check governing
all future corrections: `feedback_no_unmotivated_bh_fdr.md`.

~~**User decision (2026-05-19 evening)**: (A) — m = 6 for `d_S` alone
stands. `d_P` is permanently dropped from matrix distances on
`D_coph`. The LaTeX §sssec:methods_compare_stats with m = 6 for
`d_S` is correct as written.~~

~~**Action item (non-writing-agent)**: update
`methods/methods_revision_2026-05-18_cophenet.md` Statistical
Inference subsection (currently says "m = 12 (six bands × two
distances `d_S`, `d_P`)") to reflect m = 6 for `d_S` alone, with the
Spearman-as-natural-match-for-cophenet justification carried over
from the current LaTeX §sssec:methods_compare_ctm rationale
paragraph.~~

~~**LaTeX consequence**: no change.~~

---

# Part V — What is already correct (sanity audit)

Major checks that pass — no patch needed:

## Methods §ssec:methods_lrg

- Combinatorial Laplacian `L̂ = D̂ − W` explicitly named; random-walk Laplacian explicitly ruled out [matches `methods_grassmann_cluster_extent.md` §1].
- `τ' = 1/λ_max` canonical scale stated [matches `bands/01_beta.md` and `methods_revision_2026-05-18_cophenet.md`].
- Continuous-spectrum rationale for retiring τ-sweep + spectral-gap selector stated [matches `lrg_outlier_case_fully_connected.md` memory].
- Two coupled objects (`D_coph` + `U_k`) explicitly framed; bare `D(τ)` is "not read directly" [matches methods_revision binding directive].
- `D(τ') = (1 − δ_ij) / ρ̂_ij(τ')` correct (uses density matrix `ρ̂`, not propagator `K̂`) — M2 RESOLVED [matches `methods_section_review_2026-05-19.md`].
- `D_coph` ultrametric "by construction (cophenetic image of any agglomerative linkage)" — correct framing, not the Villegas (i)-(iii) conditions [matches methods_revision anti-pattern #4].
- `U_k = span{φ_2, …, φ_{k+1}}` — discards trivial mode `φ_1` [matches `methods_grassmann_cluster_extent.md` §1].
- Two probes "not derivable from each other" — correct mechanistic independence framing.

## Methods §ssec:methods_compare + §sssec:methods_compare_rawfc

- Phase-triangle convention `T_d(p, b) < 0 ⇔ trace` stated.
- Substrate per-pair shifts `Δ^raw_task`, `Δ^raw_rest` with independent half-baselines explicit.
- `ρ^raw_split` defined as cohort-level statistic.
- `d_F` retained "as a drift diagnostic for global signal-quality variation across the recording session rather than as a trace probe" [matches `feedback_dP_framing.md` + `result_1_raw_fc_phase_trace.md` memories].

## Methods §sssec:methods_compare_ctm (per-pair on D_coph)

- `d_S = 1 − ρ_S` matrix-level Spearman + justification of Spearman-as-natural-match for D_coph.
- `ρ_split^coph = ρ_S(Δ_task, Δ_rest)` definition with independent half-baselines.
- `T_ctm = −ρ_split^coph` sign relation.
- Methodological rationale paragraph on cophenet vs raw D (single-scale vs multiscale; continuous spectrum makes τ-sweep degenerate).
- `ρ_drift^coph` rest-only counterpart correctly defined with both rsPre + rsPost halves.

## Methods §sssec:methods_compare_grassmann (Grassmann)

- `U^Φ_k` and chordal Grassmann distance `d_G² = k − ‖A^T B‖²_F = Σ sin²θ_i` correct.
- SVD-free `k`-sweep via cumulative sum of `|A^T B|²`.
- `k ∈ {2, …, 112}` grid with `N − 1` smallest-N rationale.
- `T_G(k) < 0` trace direction explicit — m5 RESOLVED.
- Normalized cluster mass `T_G^*` ∈ [0, 1] with denominator `K · log10(R+1)` — methods_directive_2026-05-19_TG_normalization.md applied; M1 RESOLVED.
- Sparse-cell + long-pattern + isolated-lucky-cell rationale paragraph (§5e of `methods_grassmann_cluster_extent.md`).
- Per-patient `T_G^{*,p}` ∈ [0, 1] defined in Methods (not just table caption).
- `L_obs` defined as descriptive companion, "does not gate the verdict".
- Phantom-surrogate permutation null construction explicit.
- `p_mass(b) < 0.05` cohort gate on cluster mass alone — Decision 8 applied.
- "T_G^* is a significance-weighted scalar, not a trace amplitude" — methods scope-limits acknowledged.
- Cross-probe-mechanistic-independence paragraph.

## Methods §sssec:methods_compare_stats

- `n = 10` cohort stated — m8 RESOLVED.
- `R = 200` surrogates stated — m7 RESOLVED.
- 4-cycle ±δ surrogate description with "preserve per-node strength to 10⁻⁶" + "keep every edge weight within the original [0,1] interval; the marginal edge-weight distribution drifts within these constraints" — M3 RESOLVED.
- BH family m = 6 for `ρ_split^coph` (six bands) stated.
- BH family m = 6 for `d_S` on `D_coph` (six bands) stated; **note D2 above**.
- BH families "corrected independently" — m9 RESOLVED.
- Grassmann probe gated by cluster-extent permutation, no separate within-family BH.

## β results

- Head paragraph correctly identifies β as the only band with trace at both probes (modulo R1 above).
- Per-pair paragraph cites C1 p=0.005, C2 p=0.014, C4 cohort median +0.223, C3 ratio ≈23× + n_above 7/10 + p=0.005; numbers cross-verified against `bands/01_beta.md` and CSVs.
- Per-pair table: per-patient numbers cross-verified against
  `writing_directive_2026-05-19_beta_post_methods_revision.md` R6
  number-fidelity table; cohort row matches CSV.
- R3 Pat_03 dropout reference removed.
- R4 Pat_15 framed as "single-patient sensitivity".
- R5 "≈23× ratio" instead of "more than an order of magnitude" in paragraph 2 (modulo R3 above for paragraph 3).
- R6 positive=trace convention applied to per-patient table.
- Grassmann paragraph cites `p_mass = 0.005` only (R1 compliant; does not cite `L_obs` or `p_LR`).
- Grassmann paragraph: k≈40 descriptive context "≈3× surrogate; 9/10 below own surrogate" (R2 applied).
- Anatomy paragraph: 7 cophenet regions + 7 Grassmann regions match `locked/ANATOMY_LEDGER.md` β rows.
- Anti-pattern scan: no `D(τ)` for the LRG probe, no "denoising", no τ-sweep, no KC, no VI, no "diffuse" near β anatomy, no "Hippocampus + left fusiform" KC-era claim, no Pat_03 dropout reference.

---

# Part VI — Template for the remaining bands (α / γ_l / θ / γ_h / δ)

Once Methods + β are finalized via Parts I–II above, the same
structural template transposes to the five remaining band sections.
The per-band briefs `bands/02_alpha.md` … `bands/06_delta.md` carry
the locked numbers; the writing agent's task is prose, not
re-derivation.

## Common per-band paragraph structure

For trace-positive bands (α, γ_l, δ):

1. **Lead paragraph** — verdict tag + role within the cross-band
   hierarchy (e.g., "second strongest trace, per-pair only" for α).
2. **Per-pair (cophenetic) paragraph** — if trace-positive at this
   probe, cite C1, C2, C3, C4 with paired-Wilcoxon p-values per the
   locked CONTROLS battery, and the C5 epi-X gate where it has been
   run (α only).
3. **Grassmann paragraph** — if trace-positive at this probe, cite
   `p_mass` (the cohort gate per Decision 8) and the C5 epi-X
   `p_mass^epi-X` if available. Use the descriptive k-anchor
   convention (representative `k`; cohort-median observed `T_G(k)`
   in trace-negative magnitude vs surrogate median; n_below cohort
   count) per writing_directive R2.
4. **LOO caveat where applicable** — δ Grassmann full-data Pat_08
   leverage (LOO max `p_mass = 0.055`), resolved under C5 epi-X to
   LOO max `p_mass^epi-X = 0.005` fully robust; γ_l Grassmann C5
   Pat_05 leverage under epi-X (LOO max `p_mass^epi-X = 0.159`).
5. **Anatomy paragraph** — region list + gate citation (A1+A3 for
   cophenet probes, A3 alone for Grassmann probes); cross-probe
   complementarity if both probes are trace-positive.
6. **Synthesis** — per-patient probe consensus (e.g., "X of ten
   patients carry the trace direction at this probe").

For trace-negative bands (θ, γ_h):

1. **Lead paragraph** — verdict tag (no trace, both probes) +
   role within cross-band picture (cleanest negative reference for
   θ; closest miss for γ_h).
2. **Single short paragraph** — cite per-probe `p` values that fail
   the gate (e.g., γ_h cophenet C3 p=0.246, Grassmann `p_mass`=0.060).
3. **No anatomy paragraph** — per the locked ANATOMY_LEDGER, anatomy
   is not run for trace-negative bands.

## Band-by-band anchors (one-line cheatsheet — pulled from
   `bands/02_alpha.md` … `bands/06_delta.md`)

### α (8–13 Hz) — strong trace, only D_coph

- **Lead**: "α is the second strongest trace in the cohort; per-pair only, no Grassmann signal."
- **Per-pair C1**: p=0.010. C2: p=0.007. C3: p=0.002, ratio 8.3×, n_above 5/10 (borderline).
- **C4**: paired p=0.461 ≥ 0.05 (no degradation), sign-agree.
- **C5 epi-X cophenet (load-bearing for α)**: one-sample Wilcoxon p=0.0098; obs_rho_median^epi-X = +0.187; LOO max p=0.0195 (Pat_03).
- **Narrative point**: C5 epi-X **strengthens** decisively (ratio 8.3× → 27.7×; n_above 5/10 → 7/10) → α is a non-epi-cortex trace partly masked in the full cohort by epi-zone contacts.
- **Anatomy (cophenet probe)**: 11 named DK regions (bilateral cingulate + parahippocampal + medial OFC + caudal middle frontal + postcentral + precuneus + superior parietal); identical under C5 epi-X.
- **No Grassmann anatomy** (Grassmann is no trace at α; 4-cell run within null).

### γ_l (30–80 Hz) — strong trace, only Grassmann

- **Lead**: "γ_l carries a subspace-only trace; no per-pair signal."
- **Per-pair**: C3 fails (p=0.116) → cophenet no trace.
- **Grassmann**: `p_mass = 0.005` (mass-only gate, locked Decision 8); LOO max `p_mass = 0.040` (Pat_05) robust.
- **C5 epi-X**: `p_mass^epi-X = 0.030` passes cohort gate; LOO max `p_mass^epi-X = 0.159` (Pat_05) **LOO-fragile under epi-X**; T_G^*^epi-X contracts (raw mass 66 → 33).
- **LOO caveat to report**: "the γ_l Grassmann trace clears the cohort cluster-mass gate under both full data (`p_mass = 0.005`) and epi-X (`p_mass^epi-X = 0.030`); under epi-X the cohort signal contracts and is partly leveraged on Pat_05 (LOO max p = 0.159), to be read alongside the cohort verdict."
- **Anatomy (Grassmann probe)**: 6 DK regions (left fusiform + left inferior/middle/superior temporal + right paracentral + right pars triangularis); **left fusiform retraction from β** (left fusiform is anatomically at γ_l Grassmann, NOT at β).

### θ (4–8 Hz) — no trace

- **Lead**: "θ is the cleanest negative reference; no trace at either probe."
- **Cophenet**: C1 p=0.688, C2 p=0.278, C3 p=0.722 (anti-direction ratio −10.2×, n_above 2/10).
- **Grassmann**: `p_mass = 0.144`, `cluster_p_LR = 0.099`.
- **Use**: band-specificity benchmark for the β + α claims.

### γ_h (80–300 Hz) — no trace (closest miss)

- **Lead**: "γ_h is the closest miss in the cohort; just outside the gate at both probes."
- **Cophenet**: C3 p=0.246 (fails).
- **Grassmann**: `p_mass = 0.060` (just outside `< 0.05`); 9-cell observed contiguous run at k=19..27 sits at the null 95th percentile of 8.0.
- **C5 epi-X**: under epi-X, 8 new emergent cells appear at neighboring k, hinting at physiological-attribution scenario (per `audit_67_epi_exclusion_verdict.md`); but the cluster-mass cohort gate still fails (`p_mass^epi-X = 0.060`).
- **Use**: discussion-worthy borderline case; demoted from "weak" (under retired 8-cell hardcoded threshold) to "no trace" (under cluster-extent permutation) on 2026-05-19 Decision 6.

### δ (0.53–4 Hz) — strong trace, only Grassmann (LOO caveat resolved under C5)

- **Lead**: "δ carries a subspace-only trace; cohort-level cluster-mass gate clear, with a worked-example LOO + C5 epi-X resolution."
- **Per-pair**: C3 fails (p=0.278; ratio 0.98×).
- **Grassmann**: `p_mass = 0.005` (mass-only gate); **full-data LOO max `p_mass = 0.055` (Pat_08) crosses 0.05**.
- **C5 epi-X**: cluster-mass **strengthens** under epi-X (raw mass 38 → 44; `p_mass^epi-X = 0.005`); **LOO max `p_mass^epi-X = 0.005` (Pat_02) fully robust**.
- **Narrative point** (the worked example): the full-data Pat_08 leverage is attributable to epi-zone interactions, NOT the true biological trace. Under epi-X both the cohort signal strengthens and the single-patient leverage disappears.
- **Anatomy (Grassmann probe, FULL k=57..63)**: 6 DK regions (Amy + caudal anterior cingulate + medial OFC + left fusiform + left inferior parietal + left bankssts) — overlaps published δ anchor anatomy near epi zones (`memory/epileptic_imcoh_universal.md` 1.55× cross-probe).
- **Anatomy (Grassmann probe, C5 epi-X k=33..39)**: 6 DK regions (left superior parietal + right postcentral + right rostral middle frontal + left inferior temporal + left fusiform + left inferior parietal). **Only 2/6 regions shared** (fusiform + inferior parietal).
- **Use**: the δ Grassmann probe is a **mixture of two phenomena** — an anchor-anatomy signal from epi-zone contacts (full) and a separate physiological parietal-cortex signal (C5 epi-X). Report both and the dissociation.
- **Separately**: the C4 cross-probe sign agreement (`ρ_xprobe = +0.032`, `sign(ρ_xprobe) == sign(ρ_split)`) is a descriptive note about the known δ anchor anatomy (1.55× cross-probe ratio in `memory/epileptic_imcoh_universal.md`), NOT a trace claim — C3 fails so the per-pair cophenet probe shows no trace.

## Cross-band synthesis (for §results_cohort or Discussion section)

**Locked framing (user 2026-05-19 evening)**: the cross-band
synthesis is organized around the **three analytical layers** of the
LRG pipeline — **substrate (raw FC) → per-pair multiscale (D_coph) →
subspace (Grassmann)**. Each band registers (or does not) at each of
the three layers, and the combination of all three gives the
cross-band cohort map.

The cross-band layer-by-band matrix (locked):

| Band | Substrate (raw FC) | D_coph (per-pair multiscale) | Grassmann (subspace) | Coverage tag |
|---|---|---|---|---|
| **β** | borderline (p=0.053, 6/10) | **strong** (p=0.005, ratio 23.7×, 7/10) | **strong** (p_mass=0.005; C5 strengthens) | **strong trace, both probes** |
| **α** | passes (p=0.014, 8/10) | **strong** (p=0.002; C5 epi-X strengthens 8.3× → 27.7×) | no trace (p_mass=0.10) | **strong trace, only D_coph** |
| **γ_l** | borderline (p=0.053, 7/10) | no trace (p=0.116) | **strong** (p_mass=0.005, LOO max 0.040) | **strong trace, only Grassmann** |
| **δ** | passes (p=0.042, 7/10) | no trace (p=0.278) | **strong** (p_mass=0.005, Pat_08 leverage resolved under C5) | **strong trace, only Grassmann** |
| **γ_h** | borderline (p=0.19, 7/10) | no trace (p=0.246) | no trace (p_mass=0.060, closest miss) | **no trace** |
| **θ** | borderline (p=0.14, 7/10) | no trace (p=0.722) | no trace (p_mass=0.144) | **no trace** |

**Headline cross-band methodological argument**:

- Raw FC alone gives 6–8/10 cohort agreement across every band; it
  is **band-agnostic** at the substrate level. Substrate is the
  necessary but not sufficient layer — the per-pair signal exists in
  every band's raw `|ImCoh|` matrix, but without further analysis it
  cannot identify which bands carry **task-specific** reorganization.
- The **D_coph (per-pair multiscale)** probe filters the substrate
  signal through the dendrogram-induced merge-height continuum and
  reads cross-phase persistence at the per-pair hierarchical-scale
  level. Of the six bands, only β and α survive: D_coph is
  **band-selective** at the per-pair layer.
- The **Grassmann (subspace)** probe reads cross-phase rotation of
  the leading-k Laplacian eigenmode subspaces — a structurally
  distinct question from the per-pair coalescence-scale shifts.
  Grassmann is also **band-selective**: β + γ_l + δ pass; α + γ_h +
  θ do not. The Grassmann band-selectivity is **complementary to**
  D_coph's, not redundant with it.
- The combination of all three analytical layers gives the
  cross-band cohort map. **β is the only band that registers at both
  LRG probes**; **α uniquely registers at D_coph only** (its
  per-pair multiscale signal is partly masked by epi-zone contacts
  in the full cohort and reveals under C5 epi-X); **γ_l + δ uniquely
  register at Grassmann only** (their subspace rotation is not
  reducible to per-pair coalescence-scale shifts); **γ_h + θ
  register at neither LRG probe**.

The two LRG probes (D_coph + Grassmann) **read distinct facets of
the same diffusion communication geometry**, and their cross-band
band-selectivity profiles are complementary, not redundant. The
load-bearing claim is that **the structural-memory phenomenon is
multiscale and multi-facet at the LRG layer**; the cross-band
synthesis demonstrates which bands carry which facet.

---

**A separate methodological argument (not the cross-band synthesis)**:
the **cophenet vs raw `D(τ_max)`** sensitivity table from
`bands/00_cohort.md` §2 — currently mis-labeled as "Three-layer
cohort matched-strength gating" — is an **internal sensitivity
check at the D_coph layer**, NOT the cross-band synthesis. It shows
that:

- Raw FC and raw `D(τ_max)` are operationally indistinguishable at
  matched-strength gating (6–8/10 cohort agreement in every band on
  both); the LRG diffusion at the finest `τ` does not add
  band-resolution.
- The cophenet wrap on `D(τ_max)` is what produces band-resolution
  at the D_coph layer (bands without multiscale tree-shift demote;
  β + α survive).

This is the **methodological argument for adopting cophenet over
raw `D(τ_max)`** at the per-pair multiscale layer — useful as a
methods appendix or as a sensitivity paragraph in the D_coph
methods section. **It is not the cross-band cohort synthesis.**

**Cascading edit (non-writing-agent)**: `bands/00_cohort.md` §2
currently labels the raw FC / raw D / cophenet table as "Three-layer
cohort matched-strength gating", which conflates this internal
D_coph-adoption sensitivity argument with the cross-band cohort
synthesis. Reframe the table label to "Sensitivity check:
band-resolution emerges at the cophenet wrap, not at raw
`D(τ_max)`" and add a separate cross-band synthesis section
organized around the substrate / D_coph / Grassmann three-layer
analytical structure above. `bands/01_beta.md` "Headline three-layer
cohort table" requires the same relabeling.

---

# Part VII — Order of operations for the writing agent

With D1 + D2 + S2 resolved by user decision (Part "User decisions"
above) and S1 promoted from SHOULD to MUST, the writing agent
applies the patches in the following order:

1. **Methods M1** (C4 paired-Wilcoxon gate replacement) — single sentence in §sssec:methods_compare_ctm.
2. **Methods M2** (C5 epi-X Wilcoxon gate replacement) — single paragraph in §sssec:methods_compare_stats.
3. **Methods M3** (LOO max-p sensitivity diagnostic insertion) — new short paragraph appended to §sssec:methods_compare_stats.
4. **β results R1** (drop "weak cluster-extent gate" qualifier) — single phrase in head paragraph.
5. **β results R2** (medial-temporal phrasing) — single sentence in anatomy paragraph.
6. **β results R3** (replace "more than an order of magnitude") — single phrase in interpretation paragraph.
7. **β results S1** (drop hardcoded `n_sig,k ≥ 20` per-patient bold criterion, replace with null-based empirical p-value < 0.05) — single sentence in table caption + optional one-sentence Methods addition.

After these **seven** patches the manuscript is methods-aligned and
internally consistent with the locked Decisions 8/9/10/11 +
user-decision additions. The remaining editorial items (S3 epi-X
specificity, S4 audit citations, Part VI cascading edits) can be
addressed in a follow-up pass.

The five remaining band sections (α / γ_l / θ / γ_h / δ) should then
be drafted per the Part VI template — anchored on the **substrate
+ D_coph + Grassmann three-layer analytical structure** (NOT the
raw FC / raw D / cophenet sensitivity-table framing) — with the
per-band briefs (`bands/02_alpha.md` … `bands/06_delta.md`) as the
locked source of all numbers.

## Pending non-writing-agent action items (user / methods agent / data side)

These do not block the writing agent's seven LaTeX patches but should
be tracked separately:

- **(D1)** Regenerate `data/audit/grassmann_cluster_extent/cohort_summary.csv` by re-running audit_70 with the current (all-clusters) cluster_mass() definition.
- **(S1)** Extend audit_70 to emit per-patient `T_G^{*,p,null}` distributions and per-patient empirical `p_mass^p` per band, into a new column or sibling CSV. Required before the bold cells of `tab:beta_per_patient` col. 6 can be filled deterministically under the null-based criterion.
- **(D2)** Update `methods/methods_revision_2026-05-18_cophenet.md` Statistical Inference subsection: m = 12 → m = 6 for matrix distances on `D_coph`.
- **(Part VI cascading)** Reframe `bands/00_cohort.md` §2 + `bands/01_beta.md` "Headline three-layer cohort table" to clarify: (a) the raw FC / raw `D(τ_max)` / cophenet table is an internal D_coph-adoption sensitivity argument, not the cross-band synthesis; (b) add a separate cross-band synthesis section organized around the substrate / D_coph / Grassmann three analytical layers.

## Source-of-truth references

- `locked/CONTROLS.md` (Decisions 8 + 9 + 10 + 11, 2026-05-19 pm).
- `locked/VERDICT_LEDGER.md` (post-Decision 8 verdicts).
- `locked/ANATOMY_LEDGER.md`.
- `methods/methods_grassmann_cluster_extent.md` (`T_G^*` + cluster-extent).
- `methods/methods_revision_2026-05-18_cophenet.md` (binding methods directive).
- `methods/methods_section_review_2026-05-19.md` (M1-M4 + m5-m9 + n10-n14).
- `directives/methods_directive_2026-05-19_TG_normalization.md` (normalized `T_G^*`).
- `directives/writing_directive_2026-05-19_beta_post_methods_revision.md` (R1-R6).
- `directives/writing_directive_2026-05-19_beta_red_paragraphs.md` (red-paragraph un-redding).
- `responses/feedback_response_2026-05-19_gate_refactor.md` (cascade summary).
- `bands/01_beta.md` (β brief — source of all β numbers).
- `bands/02_alpha.md` … `bands/06_delta.md` (per-band briefs for upcoming sections).

## Revision history

- **2026-05-19** (initial) — Directive issued. Verification pass on
  Methods + β LaTeX after Decisions 8-11 +
  methods_directive_TG_normalization +
  writing_directive_beta_post_methods_revision cascade. Six must-fix
  patches (3 methods, 3 β results), four editorial review items,
  two open data items needing user decision. Part VI template
  prepares the writing agent for α / γ_l / θ / γ_h / δ.
- **2026-05-19** (evening — user decisions locked) — All four
  pending items resolved: D1 → (A) all-clusters formula stands;
  D2 → (A) m = 6 for `d_S` alone stands; S1 → drop hardcoded
  `n_sig,k ≥ 20` bold criterion in favor of per-patient
  null-based empirical p_mass^p < 0.05 (promoted from SHOULD to
  MUST per user directive); S2 → keep substrate column in
  per-patient table (R6 superseded). One framing correction:
  the cross-band synthesis is organized around the three
  analytical layers (substrate / D_coph / Grassmann), not the
  cophenet-vs-raw-D sensitivity table from `bands/00_cohort.md`
  §2 (the latter is an internal D_coph-adoption sensitivity
  argument). Part VI updated. Updated totals: 7 must-fix
  patches, 3 editorial review items, 0 open data items.
  Cascading actions for non-writing-agent: (i) regenerate
  audit_70 cohort_summary.csv with the current corrected
  all-clusters script; (ii) extend audit_70 to emit per-patient
  T_G^{*,p,null} distributions + p_mass^p column; (iii) update
  methods_revision_2026-05-18_cophenet.md m=12 → m=6; (iv)
  reframe `bands/00_cohort.md` §2 + `bands/01_beta.md` "three-layer
  cohort table" sections to clarify sensitivity-argument vs
  cross-band-synthesis distinction.
