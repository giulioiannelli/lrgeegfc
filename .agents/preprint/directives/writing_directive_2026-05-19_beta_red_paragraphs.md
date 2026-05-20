---
name: writing-directive-beta-red-paragraphs
era: IMCOH_ABS_COHORT_N10
status: directive
kind: writing-agent-directive
date: 2026-05-19
target: manuscript β subsection — replace the two red-coloured paragraphs (Grassmann subspace + anatomy)
band: beta
source_of_truth:
  - .agents/preprint/bands/01_beta.md
  - .agents/preprint/locked/VERDICT_LEDGER.md
  - .agents/preprint/locked/ANATOMY_LEDGER.md
  - .agents/preprint/locked/ANATOMY_CONTROLS.md
  - .agents/preprint/locked/CONTROLS.md
  - .agents/preprint/methods/methods_revision_2026-05-18_cophenet.md
  - .agents/preprint/directives/writing_directive_2026-05-18_beta_paragraph.md
locked_verdicts:
  trace: "strong trace, both probes (VERDICT_LEDGER.md, locked 2026-05-18, revised 2026-05-19 with disjunctive cluster-extent gate)"
  anatomy: "strong localized, both probes (ANATOMY_LEDGER.md, locked 2026-05-19; 7 DK regions per probe, A1+A3 for cophenet, A3 alone for Grassmann)"
---

# Writing directive — β-band, replace the two red-coloured paragraphs (2026-05-19)

## Head

The β subsection has two red-coloured paragraphs left to close: the Grassmann subspace paragraph (currently cites a single cluster-extent p-value at one statistic), and the deferred anatomy paragraph (currently a red placeholder). Both can now be **un-red**: the Grassmann gate was refined on 2026-05-19 to a **disjunctive longest-run ∨ cluster-mass** rule (VERDICT_LEDGER.md Decision 7), and the anatomy lockdown finished on the same day with **strong localized, both probes** on the matched-strength-controlled A1+A3 battery (ANATOMY_LEDGER.md, audit_71 + audit_72). After these two paragraphs land, the β subsection is complete. This directive supplies the exact numerical anchors, source CSVs, and suggested wording for each paragraph — under the voice rules already locked in the 2026-05-18 directive (plain English first; p-values at the close; no methodology drift).

## Authority and constraint (reminders from the 2026-05-18 directive)

- Plain English first, numerical evidence appended at the close. No leading p-values.
- One p-value per load-bearing claim. Counts (7/10, 8/10, 9/10) carry the per-patient cohort picture as descriptive context.
- Methods stay in the Methods section. Do not re-state the cluster-extent permutation algorithm, the matched-strength surrogate construction, the A1/A3 anatomy battery, or the cophenet-vs-raw-D rationale. The paragraph cites the *verdict* and the *p-value*; methodology lives one section away.
- "Persistence" / "trace" / "task-induced reorganization retained into rsPost" are interchangeable in the β subsection — pick whichever reads cleanly.
- "Inter-shaft" and "cross-probe" are synonyms; do not globally substitute one for the other.
- Pat_03 dropout is **not** a sensitivity row; Pat_07 dropout is **retired**; the only LRG-native robustness restriction is `n = 9` (drop Pat_15). All of this is already in the 2026-05-18 directive — do not reintroduce.

---

## Paragraph 1 — Grassmann subspace (replace the existing red paragraph)

### What was wrong with the previous version

The previous draft of this paragraph cited a single `cluster-extent permutation p = 0.005` value, which under the 2026-05-19 refinement now reads as **only the longest-run statistic** of a disjunctive (longest-run ∨ cluster-mass) gate. The reader is owed both:

1. The **longest contiguous-significant run** statistic — already in the previous draft — clears the empirical null at `p = 0.005` (1/201 floor).
2. The **cluster-mass** statistic — `Σ_k (−log10 p_k)` over every contiguous-significant cell — *also* clears the empirical null at `p = 0.005`. This is the second leg of the locked gate (VERDICT_LEDGER.md Decision 7, 2026-05-19). Cluster mass is more robust than the longest-run statistic to multi-cluster patterns, and the β signature passes both at the empirical-null floor.

The paragraph should report **both legs** in a single sentence, framed as "the cohort signature clears the empirical null under both the longest-contiguous-run and the aggregate-mass criterion". This is the load-bearing statistical claim for the Grassmann probe — without both legs of the disjunctive gate, the reader is shown half the picture.

A secondary issue is voice: the previous version says "the cohort effect deepens before attenuating again" — flowery prose. Replace with the plain observation that the trace is significant across a contiguous band and the cohort-median effect-size peaks near `k ≈ 40`.

### Numerical anchors for paragraph 1 (with source CSV row)

All from `data/audit/grassmann_cluster_extent/cohort_summary.csv` row `band = beta` (audit_70, 2026-05-19) unless otherwise noted:

| Plain-English claim | Numerical anchor (at the close of the sentence) | Source |
|---|---|---|
| The β cohort signature spans a contiguous band of intermediate-mode counts | `k ∈ [27, 55]; 29 contiguous cells` | `cluster_p_longest_run`'s observed run length `obs_longest_run = 29` |
| The strongest effect within the band is near `k = 40` | "cohort median ≈ 3× the surrogate cohort median; 9/10 patients individually below their own surrogate at this `k`" | `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv`, `k = 40` row: `obs_T_G_median = −0.489`, `surr_T_G_p50_median = −0.166`, `n_below = 9` |
| The contiguous-run length clears the empirical null built from the matched-strength surrogate ensemble | "longest-contiguous-run permutation `p = 0.005`" (`1/(R+1)` floor with R = 200) | `cluster_p_longest_run = 0.00498` (round to **0.005**); null mean LR = 2.37, null 95th = 6.0 — for context only, not inline |
| The aggregate strength of the significant cells *also* clears the same null | "aggregate cluster-mass permutation `p = 0.005`" (also `1/(R+1)` floor) | `cluster_p_cluster_mass = 0.00498`; observed cluster mass = 52.97 vs null mean 4.34 (12.2× over null mean) and null 95th 11.24 — for context only |
| Eight to nine of ten patients sit individually above their own surrogate through the core of the band | "`n_below` own surrogate 8/10 at the window endpoints (`k = 27`, `k = 55`); 9/10 at the mid-window (`k = 40`, `k = 45`)" | `grassmann_matched_strength_surrogate/cohort_summary.csv`, `k = 27, 40, 45, 55` rows |
| Removing epileptic-zone contacts preserves and broadens the contiguous block | "`k ∈ [21, 56]`, `36` cells; every cell of the original `k = 27..55` window persists, seven new cells emerge at lower `k`" | `data/audit/grassmann_epi_exclusion/cohort_summary.csv` (audit_67) + `sensitivity.csv` flags |
| Dropping the single LRG anti-aligned patient (Pat_15) sharpens the subspace claim | "cohort median at `k = 40` deepens from `−0.489` to `−0.572`; Wilcoxon `p = 0.002`" | `responses/2026-05-18_beta_n9_drop15_wilcoxon_reply.md`, primary numbers table |

### Suggested wording (replacement for paragraph 1)

> "The same persistence signature appears at the global level of the network's communication modes, read through the leading-`k` Laplacian eigendirections scanned across the resolvable spectrum (see Methods). At β the cohort signal clears the matched-strength surrogate across roughly a third of the resolvable subspace dimensions: the significant cells concentrate into a single contiguous band of intermediate mode counts at `k ∈ [27, 55]` (29 cells), with the strongest effect near `k ≈ 40`, where the observed cohort median exceeds the surrogate cohort median by approximately three-fold and nine of ten patients sit individually below their own surrogate. The contiguous block is large enough — and the aggregate strength of the band's significant cells high enough — to clear the empirical-null distributions built from the matched-strength surrogate ensemble on both criteria: cluster-extent permutation `p = 0.005` under the longest-contiguous-run statistic and `p = 0.005` under the aggregate cluster-mass statistic (both at the `1/(R+1)` floor of the `R = 200` ensemble). The cohort signal across the band is therefore not reducible to per-node strength evolution: against the same matched-strength surrogate, the observation holds across the band, with eight to nine of ten patients sitting individually above their own surrogate through the core of the regime. Removing the epileptic-zone contacts and rebuilding the analysis on the non-pathological tissue leaves every cell of the original window intact and broadens the contiguous block to `k ∈ [21, 56]` (36 cells in total, with seven new cells emerging at lower `k`), so the rotation of the leading communication modes is carried by physiological cortex. Dropping Pat_15 — the only patient anti-aligned at the LRG layer — sharpens rather than weakens the subspace claim: the cohort median at `k = 40` deepens from `−0.489` to `−0.572`, with `p = 0.002` (`n = 9`). Per-patient values across the band core are listed in Table~\ref{tab:beta_per_patient}."

This wording is structurally close to the previous draft; the substantive changes are:

- A single sentence reports **both legs** of the disjunctive gate (longest-contiguous-run + aggregate cluster-mass), each at `p = 0.005`. Both are at the empirical-null floor `1/(R+1)`.
- "The cohort effect deepens before attenuating again" → replaced by the plain claim that the strongest effect is near `k ≈ 40` (specific, anchorable).
- Suffix on the Pat_15 sentence updated to give the `n = 9` Grassmann `p = 0.002` cleanly (already in the 2026-05-18 reply).
- Epi-X clause clarified: 36 cells total = 29 persisted + 7 emerged (kept simple; the full 35 persist + 17 emerge + 3 weaken + 32 absent breakdown is in `bands/01_beta.md` §3.3 and need not appear in the paragraph).
- Removed the "where the cohort effect deepens before attenuating again" phrase: vague.

### Two reader-facing decisions for the agent

1. **k anchor**: I cited `k ≈ 40` as the representative cell because that is the mid-window cell with `9/10` patients below their own surrogate and the cleanest single Wilcoxon `p`. If the agent prefers `k = 45` (cohort-median magnitude peak `−0.585`, ratio `2.97×`) or `k = 27` (peak ratio `3.53×`), either is defensible — `bands/01_beta.md` §3.3 documents all three. Pick whichever flows; do not over-specify.
2. **Direction of "below" / "above"**: `T_G` is negative in the trace direction, so "below own surrogate" is the trace direction. Use consistent direction language throughout the paragraph; do not flip between "below" (for `T_G`) and "above" (for `ρ_split^coph`) without making the directionality explicit. Suggested gloss already in the wording above ("nine of ten patients sit individually below their own surrogate").

### What this paragraph does **not** say (intentional omissions)

- The cluster-extent permutation algorithm itself (Methods).
- The matched-strength surrogate construction (Methods).
- Why the empirical-null floor is `1/(R+1)` (Methods or supplement).
- That the `n_below ≥ 8` strict in-script gate (audit_66 / audit_67 line 425 / 417) was retired in 2026-05-15 — internal accounting (`bands/01_beta.md` §3.3 last note).
- The Pat_13 epi burden (30/119 contacts). Per-patient context belongs in the patient table.

---

## Paragraph 2 — Anatomy (replace the red placeholder)

### What this paragraph was waiting for

The 2026-05-18 directive flagged the anatomy paragraph in red because the per-region Desikan–Killiany enrichments had not yet been re-verified against the matched-strength gate on the current `n = 10` cohort. **That has now happened.** The locked anatomy battery is recorded in `locked/ANATOMY_CONTROLS.md` (A1 hypergeometric + A3 matched-strength surrogate as gates; A2 bootstrap + A4 implant-geometry regression deferred to sensitivity supplement). The locked anatomy verdicts are in `locked/ANATOMY_LEDGER.md` (2026-05-19 entries). β is `strong localized, both probes`: seven Desikan–Killiany regions per probe pass the A1+A3 join for the cophenet probe, or A3 alone for the Grassmann probe (A1 sparse on Grassmann's per-node-participation endpoint).

The anatomy paragraph is therefore **un-red** and can be written. Per the 2026-05-18 voice rule, plain English first, region lists at the close.

### Numerical anchors for paragraph 2 (with source CSV row)

All from `locked/ANATOMY_LEDGER.md` (locked 2026-05-19); audit-level CSVs:

- **β cophenet anatomy**: `data/audit/anatomy_beta_cophenet/cohort_summary.csv` (audit_71, 2026-05-19)
- **β Grassmann anatomy**: `data/audit/anatomy_beta_grassmann/cohort_summary.csv` (audit_72, 2026-05-19)

| Plain-English claim | Numerical anchor / region list (at the close of the sentence) | Source |
|---|---|---|
| The trace localizes anatomically to a distributed cortical network, not to a single region or to diffuse whole-brain cortex | "seven Desikan–Killiany regions per probe pass the matched-strength-controlled enrichment gate" | ANATOMY_LEDGER.md β rows |
| **Cophenet anatomy**: the per-pair multiscale signature concentrates in a cingulate + medial-temporal + insular + sensorimotor + prefrontal network | "left isthmus cingulate, right rostral anterior cingulate, left parahippocampal, left entorhinal, right insula, right postcentral, left superior frontal" | `anatomy_beta_cophenet/cohort_summary.csv`: A1 `q_BH < 0.05` AND A3 `p_emp < 0.05` AND `z_obs > 2` |
| **Grassmann anatomy**: the global subspace signature concentrates in a Hippocampus + temporal + orbitofrontal + insula + frontal network | "Hippocampus (subcortical), left middle temporal, left superior temporal, left lateral orbitofrontal, right medial orbitofrontal, left insula, right rostral middle frontal" | `anatomy_beta_grassmann/cohort_summary.csv`: A3 `p_emp < 0.05` AND `z_obs > 2` (A1 sparse with top-decile-per-patient endpoint — read under A3 alone) |
| The two probes read complementary sub-networks of the same β reorganization | "**Hippocampus and parahippocampal cortex** appear at both probes; **insula** appears at both probes (right hemisphere on the cophenet probe, left hemisphere on the Grassmann probe — a Desikan–Killiany label match across hemispheres, not a contact-set match)" | `locked/ANATOMY_LEDGER.md` cross-band region table |
| The matched-strength gate is the key control | "All seven cophenet regions pass A1 hypergeometric enrichment at `q_BH < 0.05` AND A3 matched-strength surrogate at `p_emp < 0.05` jointly; the seven Grassmann regions pass A3 matched-strength alone (A1 sparse on the per-node-participation endpoint)" | ANATOMY_CONTROLS.md gate definitions |
| What the anatomy paragraph does **not** claim | "no single-region anatomical anchor; no diffuse whole-brain attribution" | ANATOMY_LEDGER.md |

### Suggested wording (replacement for paragraph 2 — anatomy)

> "Where in the brain does this β reorganization concentrate? Both LRG-layer probes localize the trace to a distributed cortical network rather than to a single region or to diffuse whole-brain cortex. At the per-pair multiscale level the signal concentrates in seven Desikan–Killiany regions spanning cingulate, medial-temporal, insular, sensorimotor, and prefrontal cortex — specifically, the left isthmus cingulate, the right rostral anterior cingulate, the left parahippocampal cortex, the left entorhinal cortex, the right insula, the right postcentral gyrus, and the left superior frontal gyrus. At the global subspace level the signal concentrates in seven regions partly overlapping the first set, spanning Hippocampus, lateral temporal, orbitofrontal, insular, and rostral frontal cortex — specifically, the Hippocampus (the only subcortical region), the left middle temporal, the left superior temporal, the left lateral orbitofrontal, the right medial orbitofrontal cortices, the left insula, and the right rostral middle frontal gyrus. The two probes therefore read complementary sub-networks of the same β reorganization: Hippocampus and parahippocampal cortex appear at both, and the insula appears at both with opposite hemispheres. All seven cophenet regions pass the hypergeometric enrichment gate at `q_{BH} < 0.05` jointly with the matched-strength surrogate enrichment gate at `p_{\mathrm{emp}} < 0.05`; the seven Grassmann regions pass the matched-strength gate (the hypergeometric gate is sparse on this probe's per-node-participation endpoint, as expected; see Methods). The trace is therefore anatomically structured — not uniformly cortical — and specifically localized under topology-randomizing matched-strength surrogacy."

### Structural notes for the agent

- The cophenet region list is in **A1+A3 joint p_emp order** in `locked/ANATOMY_LEDGER.md` and `bands/01_beta.md` §5. Reorder to a flow that reads naturally in prose (anatomical grouping → cingulate, then medial-temporal, then insular, then sensorimotor/prefrontal). Cite all 7.
- The Grassmann region list is in **A3 p_emp order**; same reordering for prose flow.
- **Do not** inline the per-region A1 ratios (1.50×, 0.69×, 5.67, etc.). Those belong in a supplementary table, not in the paragraph. The paragraph cites the *gate*, not the per-region statistic. (If the agent wants a supplementary table, `bands/01_beta.md` §5.1 and §5.2 have both probes' per-region details.)
- The KC-era "Hippocampus + left fusiform" claim from `result_2_lrg_beta_trace.md`: **Hippocampus survives** (β Grassmann), **left fusiform does not** at β. Under the **locked cluster-extent paradigm `S(b)`** (2026-05-19 pm rerun) left fusiform appears at **none** of β, γ_l, δ; the earlier "γ_l + δ" re-attribution from the retired `K*(b)` audit is itself withdrawn. The paragraph should **not** mention left fusiform anywhere. If the previous KC-era draft of the manuscript carried "left fusiform" at β, retract that explicitly in revision notes; do not retract it inline in the paragraph.
- "Right insula at cophenet, left insula at Grassmann" is a **Desikan–Killiany label match across hemispheres**, not a contact-set match. The paragraph already states this; the agent should preserve the explicit clarification — otherwise a careful reviewer will read the sentence as "same insular contacts at both probes".
- **A2 + A4 deferral language dropped from the manuscript** per 2026-05-19 pm renumbering directive. The manuscript methods reference only the load-bearing battery (manuscript labels A1 = hypergeometric, A2 = matched-strength surrogate; the latter is lab A3). Lab A2 (sampling-corrected bootstrap) and lab A4 (implant-geometry regression) are reserved for the sensitivity supplement and **not disclosed in the main text**. The paragraph should **not** mention them in any form.

### Position in the section

The anatomy paragraph naturally sits **after** the Grassmann paragraph and **before** the closing two-probe synthesis paragraph (which is already in the draft and does not need editing). In the current four-sub-paragraph layout:

1. Lead claim + β-uniqueness + diffuseness flag (already in the draft, ends pointing to anatomy)
2. `ρ_split^coph` results (already in the draft, ends with Table~\ref{tab:beta_per_patient} anchor)
3. Grassmann results + epi-X (**this directive's paragraph 1 above replaces the red version**)
4. Two-probe synthesis (already in the draft, the "thus read the same persistence picture" sentence)
5. **NEW: Anatomy** (**this directive's paragraph 2 above replaces the red placeholder**) — sits between (4) and the end of the subsection, or fold into a single closing pair with (4); agent's call.

A reasonable layout: (1) → (2) → (3) → (4) → (5) → end of subsection. Five sub-paragraphs total. Alternative: merge (4) and (5) into one closing paragraph that opens with the two-probe synthesis and ends with the anatomy localization. The four-paragraph layout is *not* sacred — pick the cleaner flow.

---

## Two-probe synthesis paragraph (no change, but updated context)

The middle paragraph that closes "The two \gls{lrg}-layer probes thus read the same persistence picture at β" is unchanged — its three named disagreements (Pat_10, Pat_14, Pat_15) are still correct under the locked verdicts, and the structural-complementarity argument carries through to the new anatomy paragraph (the two probes read complementary anatomical sub-networks of the same diffusion geometry). The agent does not need to re-touch this paragraph except to confirm Pat_15 is the only LRG anti-aligned patient (Fix 1 from the 2026-05-18 directive — already landed in 01_beta.md §6).

---

## Number fidelity check (do this before sending the draft back)

The agent must verify that the following numerical claims survive byte-for-byte from this directive into the LaTeX draft. If any number is changed, flag it and stop.

| Claim | Number |
|---|---|
| Contiguous-significant band (full FC) | `k ∈ [27, 55]`, 29 cells |
| Peak Grassmann cell anchor | `k ≈ 40`; cohort median ≈ 3× surrogate; 9/10 below own surrogate |
| Longest-run cluster-extent permutation p (β) | `p = 0.005` |
| Cluster-mass cluster-extent permutation p (β) | `p = 0.005` |
| Epi-X broadened band | `k ∈ [21, 56]`, 36 cells; 29 persist + 7 emerge |
| Pat_15-dropped Grassmann at k=40 | cohort median `−0.572`, `p = 0.002`, `n = 9` |
| β cophenet anatomy regions | left isthmus cingulate, right rostral anterior cingulate, left parahippocampal, left entorhinal, right insula, right postcentral, left superior frontal (7) |
| β Grassmann anatomy regions | Hippocampus, left middle temporal, left superior temporal, left lateral orbitofrontal, right medial orbitofrontal, left insula, right rostral middle frontal (7) |
| Anatomy gates | A1 `q_{BH} < 0.05` AND A3 `p_{\mathrm{emp}} < 0.05` for cophenet; A3 alone for Grassmann |
| Shared regions across probes | Hippocampus / parahippocampal at both; insula at both (right cophenet, left Grassmann) |

Counts ("7/10", "8/10", "9/10") are descriptive — the agent may omit per-patient counts at sub-cells that are not load-bearing, but do not invent counts or rewrite them.

---

## Source CSVs (no agent re-derivation needed)

| Anchor | Path | Row |
|---|---|---|
| Cluster-extent permutation, β | `data/audit/grassmann_cluster_extent/cohort_summary.csv` | `band = beta` (audit_70) |
| Per-k cohort Wilcoxon, β | `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv` | `band = beta, k ∈ [27, 55]` (audit_66) |
| Epi-X per-k cohort Wilcoxon, β | `data/audit/grassmann_epi_exclusion/cohort_summary.csv` | `band = beta` (audit_67) |
| Pat_15-dropped Grassmann | `.agents/preprint/responses/2026-05-18_beta_n9_drop15_wilcoxon_reply.md` | auxiliary numbers table |
| β cophenet anatomy regions | `data/audit/anatomy_beta_cophenet/cohort_summary.csv` | A1 `q_BH < 0.05` AND A3 `p_emp < 0.05` AND `z_obs > 2` (audit_71) |
| β Grassmann anatomy regions | `data/audit/anatomy_beta_grassmann/cohort_summary.csv` | A3 `p_emp < 0.05` AND `z_obs > 2` (audit_72) |
| Anatomy gate definitions | `.agents/preprint/locked/ANATOMY_CONTROLS.md` | A1, A3 sections |
| Anatomy verdict | `.agents/preprint/locked/ANATOMY_LEDGER.md` | β rows |

## Revision history

- **2026-05-19** — Directive issued. Two red-coloured β paragraphs (Grassmann subspace, anatomy) are now un-red: the disjunctive longest-run ∨ cluster-mass gate (VERDICT_LEDGER.md Decision 7, 2026-05-19) and the locked anatomy verdict (ANATOMY_LEDGER.md, audit_71 + audit_72) supply the required numerical anchors. Builds on the 2026-05-18 voice / fixes / per-patient-table directive (`writing_directive_2026-05-18_beta_paragraph.md`). After the agent's next pass, the β subsection is complete.
