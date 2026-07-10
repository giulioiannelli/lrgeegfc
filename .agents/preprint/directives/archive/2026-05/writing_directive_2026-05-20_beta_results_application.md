---
name: writing-directive-beta-results-application
era: IMCOH_ABS_COHORT_N10
status: open
kind: writing-directive
scope: results-section LaTeX edits for the β subsection, derived from the 2026-05-20 cross-check against the post-audit Methods + locked audit CSVs
created: 2026-05-20
companion: writing_directive_2026-05-20_methods_audit_application.md (Methods edits, already in flight), METHODS_AUDIT_ISSUES.md (master spec), bands/01_beta.md (lab brief)
target: writing agent producing manuscript Results §`ssec:results_beta`
---

# β Results-section edits — writing-agent directive

**Head.** Every cohort-level number in the β paragraph is verified against the locked audit CSVs and matches to the cited decimal. The required edits are five corrections (R1–R5), three clarifications (C1–C3), four style/voice items (S1–S4), three figure label/caption placeholders (F1–F3), and one structural flow item (Fl1). The numerical claims are sound; the work is copy-editing the β paragraph against the now-finalized Methods conventions and Nature Neuroscience descriptive style. Return the revised LaTeX block + a one-paragraph change log.

## What to read first

1. **This directive** — edit list + verified-numbers table.
2. **`writing_directive_2026-05-20_methods_audit_application.md`** — the in-flight Methods edits the β paragraph is harmonized against.
3. **`bands/01_beta.md`** — lab β brief; the load-bearing per-patient numbers and anatomy region lists come from here.
4. **`locked/CONTROLS.md`** §§C1–C5 + **`locked/VERDICT_LEDGER.md`** β row — verdict source.

## Numerical verification — all claims hold

Every numerical claim in the β paragraph was cross-checked against the source CSV row and matches to the cited decimal.

| # | β paragraph claim | CSV row | Result |
|---|---|---|---|
| 1 | ρ_split^raw = +0.258 | `raw_fc_matched_strength/cohort_summary.csv` β | 0.2577 ✓ |
| 2 | 21× noise floor (+0.012) | same, surr_median = 0.0121 | 21.3× ✓ |
| 3 | 6/10 above surrogate | same, n_above_surrogate | 6/10 ✓ |
| 4 | substrate p = 0.053 | same, paired_wilcoxon_p | 0.0527 → 0.053 ✓ |
| 5 | ρ_split = +0.222 (within-baseline) | `ctm_triangle/cohort_summary.csv` β rho_split_median | 0.2222 ✓ |
| 6 | cohort Wilcoxon p = 0.005 | same, wilcoxon_split_gt_0_p | 0.0049 ✓ |
| 7 | 8/10 above drift, p = 0.014 | same, n_above_drift + wilcoxon_split_gt_drift_p | 8/10 + 0.0137 ✓ |
| 8 | ρ_xprobe = +0.223 | same, rho_xprobe_median | 0.2228 ✓ |
| 9 | non-degradation Wilcoxon p = 0.385 | `ctm_triangle/c4_wilcoxon_cohort.csv` β paired_p | 0.3848 ✓ |
| 10 | ρ_split = +0.221 (matched-strength), 23× ratio (+0.221 vs +0.009), 7/10, p = 0.005 | `matched_strength_surrogate_split_baseline/cohort_summary.csv` β | obs 0.2206 / surr 0.0093 = 23.7× ; 7/10 ; p = 0.0049 ✓ |
| 11 | n=9 (drop Pat_15) ρ_split = +0.221→+0.230, p = 0.010 | recomputed Wilcoxon on `per_patient_per_band.csv` β | median 0.2302, p = 0.0098 ✓ |
| 12 | p_mass = 0.005 | `grassmann_cluster_extent/cohort_summary.csv` β cluster_p_cluster_mass | 0.0050 ✓ |
| 13 | "roughly a third" of resolvable subspace dimensions individually significant | `grassmann_cluster_extent/per_k_obs_p.csv` β | 40/111 = 36% ≈ 1/3 ✓ |
| 14 | k≈40 cohort obs T_G ~3× surr median in magnitude | `grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv` β k=40 | 0.489 / 0.166 = 2.94 ✓ |
| 15 | 9/10 above own surrogate at k=40 | same, obs_z column | 9/10 ✓ |
| 16 | n=9 (drop Pat_15) k=40 +0.489→+0.572, p = 0.002 | recomputed | 0.572 ; p = 0.00195 ✓ |
| 17 | β epi-X Grassmann cluster mass survives | `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` β | p_mass = 0.0050, c5_pass=True ✓ |
| 18 | per-pair (cophenet) anatomy: 7 DK regions | `anatomy_beta_cophenet/cohort_summary.csv` passes_both excl. Unk | 7 cortical regions ✓ |
| 19 | Grassmann anatomy: 7 DK regions | `anatomy_beta_grassmann_clusterext/cohort_summary.csv` passes_A3 | 7 regions ✓ |
| 20 | medial-temporal overlap (lh-parahippocampal ↔ Hip) + bilateral insula | both anatomy CSVs | ✓ |
| 21 | Per-patient Table 1 values (all 10 patients, ρ_split + z) | `matched_strength_surrogate_split_baseline/per_patient_per_band.csv` β | matched to the cited 3-decimal precision ✓ |

**Action**: do not change any number. The numerical content is correct; the edits below are framing, terminology, and style.

## Corrections (R1–R5) — apply

### R1 — Drop Pat_03 "f_s = 1024 Hz" from the per-patient Table 1 Notes column

- **Location**: `\patient{03}` row of `tab:beta_per_patient`, last column ("Notes").
- **Current**: `\(f_s=\) \qty{1024}{\hertz}`.
- **Fix**: Replace with `---`. Pat_03 is a full cohort member treated identically at the analysis layer; the sampling-rate difference is absorbed at the config layer (`nperseg_for_fs(fs)`, `FS_OVERRIDES`) and does not propagate to per-figure or per-row annotation. The previous "outlier / 1024 Hz" framing was retired 2026-05-18 (see project memory `feedback_pat03_no_dropout.md`).
- **Rationale**: An entry in the Notes column reads to a reviewer as "this patient is treated specially in some way". The locked policy is that Pat_03 is treated identically. Drop the annotation.

### R2 — Reframe "preserves the cluster-mass gate and strengthens the cohort signal" (epi-X claim)

- **Location**: paragraph immediately before the Grassmann heatmap figure: "Removing the epileptic-zone contacts and rebuilding the analysis on the non-pathological tissue preserves the cluster-mass gate and strengthens the cohort signal, so the rotation of the leading communication modes is carried by physiological cortex."
- **Issue**: The gate `p_mass` value is at the empirical-null floor `1/(R+1) = 0.0050` in both the full-cohort and epi-X analyses; "strengthens the cohort signal" overstates the change at the gate level (it cannot drop below 0.0050). What changes under epi-X is the **k-cell support** of `T_G^*`: in the full-cohort analysis the contiguous run is ~29 cells; under epi-X the contiguous extent broadens (per `bands/01_beta.md` brief and `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` β `obs_longest_run_epiX = 36`).
- **Fix direction**: Reword to "Restricting the analysis to non-epileptic cortex preserves the cluster-mass gate at the empirical-null floor (p_mass^epi-X = 0.005) and broadens the support of T_G* across the resolvable spectrum, so the rotation of the leading communication modes is carried by physiological cortex." Cite `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` β.
- **Rationale**: Accurate description of what changes (k-extent) and what stays fixed (floor-level p) under epi-X.

### R3 — Demote "the strongest alternative is per-node strength evolution" — methods-section voice in Results

- **Location**: paragraph "Eight of ten patients carry ρ_split…", second-to-last sentence: "The strongest alternative is per-node strength evolution: on a dense weighted substrate, simply preserving each node's strength can reproduce part of the cohort signal, and indeed at the substrate level it does. At the LRG layer it does not …"
- **Issue**: "The strongest alternative is" is Methods-section framing (the 5-point critical preamble vocabulary; appears explicitly in the β brief §3.1 preamble). In a Results paragraph it reads as out-of-register — Results state findings, Methods state nulls.
- **Fix direction**: Convert to a descriptive contrast without the "strongest alternative" framing — e.g., "Per-node strength evolution alone — preservation of each node's strength sequence under topological rewiring — reproduces part of the substrate-level signal but does not reproduce the cophenetic-layer trace: against the matched-strength surrogate, the observed cohort median exceeds the surrogate cohort median by an ≈ 23× ratio (+0.221 vs +0.009), with seven of ten patients individually above their own surrogate (p = 0.005)." Numbers are unchanged.
- **Rationale**: The factual content survives; the voice shifts from "we enumerate the alternatives we ruled out" (Methods) to "the observed effect exceeds the matched-strength surrogate by Xx" (Results).

### R4 — Reframe the "+0.222 vs +0.221" parenthetical at the end of the Table-1 paragraph

- **Location**: paragraph immediately preceding the Pat_15-sensitivity sentence: "The cohort row of Table~\ref{tab:beta_per_patient}, col. 3 reports the matched-strength split-baseline realization of ρ_split (+0.221); the within-baseline realization cited in the text is the closely related +0.222, the two differing only by the choice of null built from the same halved-data noise budget."
- **Issue**: This is methods-side reconciliation language inside a Results paragraph; it interrupts the narrative flow and the reader does not need the explanation at this point — the two numbers are identical to the cited 3-decimal precision modulo the null construction (within-baseline vs matched-strength) and the difference is reconcilable by reading the Methods.
- **Fix direction**: Move the reconciliation note into a Table-1 footnote (or absorb into footnote (a) on `T_G^{*,s}` and add a separate footnote (d): "The cohort-row ρ_split value (+0.221) is the matched-strength realization; the within-baseline realization +0.222 cited in the text differs only by the choice of null and is reconciled in Methods §sssec:methods_compare_ctm.")
- **Rationale**: Caption-level reconciliation is the conventional place for methodological cross-references.

### R5 — Reframe "LRG anti-aligned patient" — define on first use or substitute

- **Location**: appears three times in the β paragraph and once in Table-1 Pat_15 Notes column: "the only LRG anti-aligned patient — right-hemisphere-only implant, no epileptic contacts" (body) and "right-hemisphere only; LRG anti-aligned" (table Notes).
- **Issue**: "LRG anti-aligned" is internal lab jargon (informally: "the only patient where the cohort-paired Wilcoxon would change sign if this patient's contribution were doubled" or similar). The phrase does not appear in the locked Methods and a reviewer encountering it cold will not know what it means.
- **Fix direction**: At first use replace with a definition-by-context — "Pat_15, the single patient whose cophenetic-layer cohort contribution is opposed to the cohort consensus at β (z = +1.27 against own matched-strength surrogate vs the cohort median z = +6 to +9 in the trace direction, Table~\ref{tab:beta_per_patient}; consistent with the right-hemisphere-only implant geometry, no epileptic-zone contacts)" — then refer to "Pat_15" alone in subsequent mentions. In the Table-1 Notes column use "right-hemisphere only" alone (drop "LRG anti-aligned").
- **Rationale**: Self-defining language with a Methods-traceable basis is the Nature-Neuroscience-style alternative to in-house shorthand.

## Clarifications (C1–C3) — apply

### C1 — Disambiguate `T_G*` (band-level, eq. methods_TGstar) vs `T_G^{*,s}` (per-patient, eq. methods_TGstar_perpatient)

- **Location**: paragraph describing the Grassmann subspace probe, the sentence "the cluster mass aggregating these significant cells clears the empirical-null distribution built from the same surrogate ensemble (cluster-mass permutation p_mass = 0.005, …)".
- **Issue**: This sentence (correctly) reports `p_mass` from the band-level `T_G*`. Table 1 then reports `T_G^{*,s}` cohort median = 0.653 — a numerically distinct quantity built from per-patient empirical p-values rather than from the cohort-level Wilcoxon. The two quantities are both ∈ [0,1] but they are not the same number for the same band, and a reviewer scanning the band-level normalized value `T_G*(β) = 0.273` (from `grassmann_cluster_extent/cohort_summary.csv` `obs_cluster_mass_neglog10p_norm`) against the Table-1 cohort row `T_G^{*,s}^cohort_median = 0.653` may flag an inconsistency.
- **Fix direction**: Where the band-level value is cited in body text, cite the band-level `T_G*(β) = 0.273` explicitly (one sentence, e.g., "the band-level cluster mass `T_G*(β) = 0.273` of eq. (methods_TGstar) clears the empirical-null floor at p_mass = 0.005"). In Table 1 the `T_G^{*,s}` column is correctly labeled — keep its footnote (a) clarifying the per-patient definition, and add one sentence to that footnote: "`T_G^{*,s}` is the per-patient analogue of the band-level `T_G*(β) = 0.273` (eq. methods_TGstar); the cohort median of `T_G^{*,s}` is not equal to `T_G*(β)`."
- **Rationale**: The two quantities answer different questions (band-level cluster-mass gate vs per-patient descriptive cluster mass). Surfacing both with an explicit pointer prevents a reviewer comparison error.

### C2 — Mark the n=9 sensitivity as a single biology-driven exception, not as a general dropout pattern

- **Location**: paragraph "The cohort effect is stable across the LRG-native robustness restriction (n=9 cohort excluding Pat_15, the only LRG anti-aligned patient — right-hemisphere-only implant, no epileptic contacts)".
- **Issue**: The Methods (post-edit) gates everything at n=10 and uses LOO max-p as the only paired sensitivity. An n=9 cohort dropping Pat_15 is permitted by the locked patient-dropout policy as a single biology-driven exception (right-hemisphere-only implant), but it should be explicitly framed as such — not as a generic "robustness restriction".
- **Fix direction**: Reword to "A single-patient sensitivity restricting the cohort to n=9 by dropping Pat_15 — the one patient whose β cophenet contribution opposes the cohort consensus and whose implant covers only the right hemisphere — leaves the cohort median in the same narrow band (+0.221 to +0.230) and the gate intact (p = 0.010, n=9). Pat_15 is also the patient whose exclusion drives the LOO max-p on the Grassmann band-level cluster mass at the manuscript value p_mass^LOO = 0.005 (Table~\ref{tab:beta_per_patient}; see Methods §sssec:methods_compare_stats for the LOO diagnostic)."
- **Rationale**: Pat_15 is the single biology-driven n=9 exception per the locked policy. Frame it as such, not as a general sensitivity check.

### C3 — Sign-convention note for the audit-trail traceability (footnote, not body text)

- **Issue (internal, not visible to the manuscript reader)**: The CSV `grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv` stores `obs_T_G` with the opposite sign from the Methods convention (CSV uses "lower-is-trace", Methods uses "T_d > 0 = trace"). The β paragraph correctly reports values in the Methods convention (e.g., cohort obs T_G at k=40 = +0.489; CSV stores −0.489). This is a documentation issue for downstream verification, not a manuscript issue.
- **Action**: No edit to the manuscript text. Add a one-line internal comment in the LaTeX source — `% Note: CSV stores T_G with opposite sign convention; paragraph reports Methods-convention sign (multiply CSV obs_T_G by -1).` This protects any downstream evaluator running EVALUATION_PROTOCOL.md from raising a sign mismatch.
- **Rationale**: Audit-trail traceability. No reader-facing change.

## Style and voice (S1–S4) — apply

### S1 — Remove or soften emphatic / judgemental adjectives

The user's directive: Nature Neuroscience descriptive style; no judgemental terms ("load-bearing" forbidden — does not appear in this paragraph, but related framings do).

Items to soften:

| Current | Suggested replacement |
|---|---|
| "the strongest and most geometrically coherent task-induced reorganization" (opening sentence) | "the only band whose post-task signature clears matched-strength surrogacy at both LRG-layer probes, with the largest cohort-level effect size at the per-pair multiscale probe" |
| "thus a coherent multiscale reshaping of information-flow geometry" | "a multiscale reorganization of the information-communication geometry" |
| "The LRG step amplifies this borderline substrate signal into a strong cohort verdict at both probes" | "Wrapping the substrate through the LRG construction shifts the cohort claim above the matched-strength surrogate at both probes" |
| "the cognitive episode does not rewrite the resting communication geometry, it perturbs it in a structured, detectable, pair-specific way" | "the cognitive episode does not overwrite the resting communication geometry but leaves a structured, pair-specific shift detectable across the cohort" |
| "is the operative comparison" | "is the relevant comparison" |
| "The post-task resting state is not a return to the rsPre baseline; it inherits, in its hierarchical communication geometry, the specific pair-level signature of the cognitive episode that preceded it." | "Post-task resting connectivity does not return to the pre-task baseline; in its hierarchical communication geometry it retains a pair-specific signature of the preceding cognitive episode." |

### S2 — Tighten the "persistence picture" / "emerges" language

- "the persistence picture that emerges in this band" → "the cross-probe picture at β"
- "the same persistence signature appears at the global level" → "The same retention pattern appears at the global subspace level"
- "the two LRG-layer probes thus read the same persistence picture at β" → "the two LRG-layer probes therefore agree on the β retention pattern"

Bare "persistence" remains acceptable per the softened 2026-05-18 terminology rule (see `feedback_trace_terminology.md`); the issue here is the literary "picture / emerges" framing, not the noun.

### S3 — Replace "Taken together, these checks identify" with a declarative finding sentence

- Current: "Taken together, these checks identify a persistent reorganization of the per-pair hierarchical communication geometry: …"
- Replace: "The per-pair multiscale layer identifies a structured reorganization of the cophenetic communication geometry: …"
- Rationale: "Taken together, these checks identify" is a hedging meta-statement; the declarative form is the Nature-Neuroscience-style equivalent.

### S4 — Active voice for the cohort-level finding

- Current: "the post-task signature clears the matched-strength surrogate" (passive-construction with a non-human agent). Acceptable.
- Current: "At the LRG layer it does not — against the matched-strength surrogate (see Methods), the observed cohort median exceeds the surrogate cohort median by an ≈ 23× ratio". Acceptable (the agent is the cohort median).
- Action: no global voice change; the paragraph alternates active/passive at acceptable cadence. Apply S1–S3 edits as listed and leave voice unchanged otherwise.

## Figure label + caption placeholders (F1–F3) — apply

All three inline figures currently carry `\caption{Caption}` and `\label{fig:placeholder}`. Replace with descriptive labels + tentative placeholder captions per the user directive ("placeholders for needed updates to the caption for when they will be definitive"). The figures themselves will be reviewed and embellished after the Results section stabilizes.

### F1 — `fig_beta_rho_split.pdf`

- **Label**: `\label{fig:beta_rho_split}`
- **Placeholder caption** (tentative — to be replaced once the figure is finalized): "Cohort-level per-pair multiscale trace correlation `ρ_split` at β across the n=10 patient set. Top: per-patient observed `ρ_split` against the patient-matched matched-strength surrogate distribution (R = 200, 4-cycle ±δ rewiring; Methods §sssec:methods_compare_stats). Bottom: cohort-paired Wilcoxon against the within-baseline drift floor `ρ_drift` and the cross-probe restriction `ρ_xprobe`. **\textit{Figure placeholder — final caption pending updated figure.}**"
- **Body reference**: add `(Fig.~\ref{fig:beta_rho_split})` at the end of the sentence "Eight of ten patients also lie above their own drift-floor null …".

### F2 — `fig_beta_grassmann_heatmap_visual.pdf`

- **Label**: `\label{fig:beta_grassmann_heatmap}`
- **Placeholder caption**: "Per-`k` cohort effect of the β Grassmann subspace probe `T_G(k)` across `k ∈ {2, …, 112}`. Color encodes the per-`k` cohort-paired Wilcoxon `p_k(β)` against the matched-strength surrogate ensemble; the significance support set `S(β)` of the cluster-mass statistic `T_G*(β)` (Methods eq. methods_TGstar) is highlighted at the per-cell threshold `α_k = 0.05`. **\textit{Figure placeholder — final caption pending updated figure.}**"
- **Body reference**: add `(Fig.~\ref{fig:beta_grassmann_heatmap})` at the end of the sentence "The strongest per-`k` effect inside this regime sits near `k ≈ 40` …".

### F3 — `fig_beta_anatomy_brain.pdf`

- **Label**: `\label{fig:beta_anatomy_brain}`
- **Placeholder caption**: "Anatomical localization of the β trace at both LRG-layer probes. Left: Desikan–Killiany regions clearing the joint A1+A2 enrichment gate at the per-pair multiscale probe `ρ_split` (Methods §sssec:methods_anatomy). Right: Desikan–Killiany regions clearing the A2 enrichment gate at the Grassmann probe `T_G*`. Region overlap on the medial-temporal-lobe family and the bilateral insula is indicated. **\textit{Figure placeholder — final caption pending updated figure.}**"
- **Body reference**: add `(Fig.~\ref{fig:beta_anatomy_brain})` at the end of the anatomy paragraph.

## Structural flow (Fl1) — apply

### Fl1 — Reorder the Pat_15 n=9 sensitivity sentence

- **Current order**: "[Grassmann heatmap figure]" → "As a single-patient sensitivity, dropping Pat_15 — the only LRG anti-aligned patient — sharpens rather than weakens the subspace claim: the cohort median at `k = 40` rises from +0.489 to +0.572 (p = 0.002, n = 9). Per-patient cluster-mass values are listed in Table~\ref{tab:beta_per_patient}." → "[shared-null construction paragraph]" → "[per-patient Table]".
- **Issue**: The Pat_15-sensitivity sentence sits awkwardly between the Grassmann figure and the table; it touches both the cophenet n=9 sensitivity (cited earlier in the paragraph) and the Grassmann k=40 sensitivity. The shared-null construction paragraph immediately after it is methods-recap and breaks narrative flow.
- **Fix direction**: After applying C2 above (which merges the n=9 cophenet sensitivity and the k=40 Grassmann sensitivity into a single Pat_15 paragraph framed by the locked patient-dropout policy), the structural flow becomes: cophenet section → cophenet figure (F1) → Grassmann section → Grassmann figure (F2) → unified Pat_15 sensitivity paragraph (per C2) → Table 1 + footnotes → anatomy paragraph → anatomy figure (F3) → closing two-probes-agree paragraph. The "shared-null construction" paragraph (lines describing patient-level z and empirical p) becomes one sentence absorbed into the Table-1 introduction.
- **Rationale**: Single coherent placement for the Pat_15 sensitivity; tighter narrative; methods-recap content cited not described.

## What you must NOT do

- Do **not** change any cited number — every claim is CSV-verified.
- Do **not** introduce KC, VI(k), τ-sweep, K*(b), or any retired probe terminology.
- Do **not** invert the trace sign convention: positive `T_d` (or `ρ_split`) = trace direction throughout, per the post-audit Methods (B1 fix in the Methods directive).
- Do **not** report raw cluster-mass values (e.g., 69.76); use the normalized `T_G*(β) = 0.273` ∈ [0,1] when the band-level scalar is cited.
- Do **not** introduce cross-band BH-FDR or any other cross-band correction — per the locked decision dropped 2026-05-20 (`feedback_no_unmotivated_bh_fdr.md`).
- Do **not** add `\setcaptionlinewidth` or other per-figure tweaks until the figures themselves are finalized — the placeholder captions are tentative.
- Do **not** propose new analyses or new controls.

## Deliverable

Revised β subsection LaTeX (the entire `\subsection{β: a structurally coherent task signature in resting connectivity}` block) with:

- R1–R5 corrections applied.
- C1–C3 clarifications applied (C3 is a LaTeX-comment-only edit).
- S1–S4 style edits applied.
- F1–F3 figure labels + tentative placeholder captions + body references applied.
- Fl1 paragraph reordering applied.

Plus a one-paragraph change log naming each applied edit by ID.

After return, the user pastes the revised LaTeX back to Claude Code for verification under `EVALUATION_PROTOCOL.md` (number-level CSV cross-check, anti-pattern scan, verdict consistency). Iterate if findings; otherwise the β paragraph is locked and we move to the figure embellishment + caption finalization pass, then to α / γ_l / δ / γ_h / θ Results subsections in turn.
