---
name: writing-directive-beta-paragraph-revision
era: IMCOH_ABS_COHORT_N10
status: directive
kind: writing-agent-directive
date: 2026-05-18
target: manuscript β subsection (\ssec:results_beta)
band: beta
source_of_truth:
  - .agents/preprint/bands/01_beta.md
  - .agents/preprint/locked/VERDICT_LEDGER.md
  - .agents/preprint/locked/CONTROLS.md
  - .agents/preprint/methods/methods_revision_2026-05-18_cophenet.md
locked_verdict: "strong trace, both probes (VERDICT_LEDGER.md locked 2026-05-18, revised 2026-05-19 for cluster-extent gate)"
---

# Writing directive — β-band results paragraph revision (2026-05-18)

## Head

The drafted β paragraph is structurally close to ready for the Nature Neuroscience preprint but needs a focused set of edits before the next round. Voice rule from the author: each result is stated in plain English first; numerical evidence is appended at the end of the claim, not led with. This directive lists the mandatory fixes, the additions, and the deferrals (per-patient table to build, anatomy paragraph to flag in red as pending re-verification). Nothing in this directive describes methodology — methods stay in the Methods section.

## Authority and constraint

- Every number must come from a row in `data/audit/<probe>/cohort_summary.csv` (or the per-patient companion). The brief at `.agents/preprint/bands/01_beta.md` transcribes these and is the lookup table for the agent.
- The locked verdict is `strong trace, both probes`. Do not re-derive; document.
- Plain English first, then the statistics at the close of each claim. Do not lead with p-values and do not over-quote them — one p-value per load-bearing claim is enough, with descriptive counts (8/10, 7/10) carrying the per-patient cohort picture.
- The methodology-comparison story (cophenet vs raw D; why the cophenet step is the multiscale-resolution operator) lives in the Methods section; do not re-state it in this paragraph.

---

## Mandatory fixes

### Fix 1 — Pat_07 attribution at the LRG layer

**Current draft sentence:**
> "and when the two anti-aligned implant geometries Pat_07 and Pat_15 are dropped (remaining cohort of eight, signal preserved)"

**Problem.** Pat_07 is *pro* at `ρ_split^coph` with z > 3 (brief §6). Only anti at the substrate level. Calling Pat_07 "anti-aligned" in the LRG section is incorrect.

**At the LRG layer:**
- Pat_15 is the only LRG anti-aligned patient (right-hemisphere-only implant, 0 epi contacts).
- Pat_07 is anti at the substrate (raw FC `d_S`) only — pro at `ρ_split^coph`.

**Required edit.** Rewrite the named-patient sentence so that the LRG-layer anti-aligned anchor (Pat_15) is identified correctly, and the robustness restriction drops **only Pat_15** — not Pat_07. The legacy `n=8` restriction that also dropped Pat_07 was a substrate-layer construct (Pat_07 is anti at raw-FC `d_S` with `S = +0.005`, near-zero) and is **retired** at the LRG layer because Pat_07 is solidly pro at both LRG probes (`ρ_split^coph` z = +4.47; Grassmann `T_G(k=40)` z = −3.90). Suggested wording:

> "The result is preserved under the LRG-native pro-cohort restriction (n = 9 cohort excluding Pat_15, the only LRG anti-aligned patient — right-hemisphere-only implant, no epileptic contacts): trace direction and surrogate gap intact."

Exact patient counts and p-values to attach at end of the sentence per Fix 4.

**Patient-dropout policy (updated 2026-05-18, second update of the day).** The project rule is: **don't drop patients in analysis**, except for two narrow exceptions — (a) a single patient anti-aligned at the cohort level (biology-driven, e.g. Pat_15 right-hemisphere-only at the LRG layer); (b) genuinely problematic data (vendor corruption, sampling-rate handled at config layer — none currently active in the cohort). Pat_03 dropout was retired earlier today (config-layer-only sampling-rate handling). The Pat_07 dropout is **retired now** because Pat_07 is pro at LRG; only the substrate sign is anti (and only marginally, `S = +0.005`). Do **not** reintroduce a Pat_07-dropout row in any robustness panel, figure caption, or per-patient table. The single acceptable robustness dropout at β is `n=9` (drop Pat_15).

### Fix 2 — Cite the audit_70 cluster-extent permutation for Grassmann

**Problem.** The current Grassmann passage describes the contiguous-run length (29 cells at k ∈ [27, 55]) and the per-k effect-size ratios, but does not cite the cluster-extent permutation null. The locked Grassmann gate (revised 2026-05-19, `locked/VERDICT_LEDGER.md` Decision 6) is the empirical-null distribution of the longest contiguous-significant run built from R = 200 matched-strength phantom surrogates: `data/audit/grassmann_cluster_extent/cohort_summary.csv`. β cluster-extent p = 0.005.

**Required edit.** Append the cluster-extent p-value at the end of the contiguous-band claim, in plain language. Suggested rewording (numerical anchor at the close):

> "the significant dimensions are not scattered across the spectrum but concentrate into a single contiguous band of intermediate mode counts (k ∈ [27, 55]), with the strongest effect at k ≈ 40 where the observed cohort median exceeds the surrogate cohort median by approximately three-fold — a contiguous block large enough to clear the empirical-null distribution of contiguous-run length built from the matched-strength surrogate ensemble (cluster-extent permutation p = 0.005)."

This is the load-bearing statistical claim for the Grassmann probe; without it, the Grassmann passage has no cohort-level p-value.

### Fix 3 — Add the β-uniqueness sentence (result-level, not methodology)

**Problem.** The paragraph reads as if β is the only band under discussion. A Nature Neuroscience reader scanning this subsection cannot tell that no other band passes both LRG probes under the locked matched-strength gate. This is the headline result and must appear in the first paragraph.

**Required addition.** Insert one sentence near the top of the section, immediately after the structural-coherence framing. Suggested:

> "Across the six bands tested, β is the only one whose post-task signature clears the matched-strength surrogate at both LRG-layer probes — the per-pair multiscale measure and the global subspace measure described below; the remaining bands carry at most a partial signature (α at the per-pair level only; γ_low and δ at the subspace level only and at a weak cluster-extent gate; γ_high and θ none)."

Do NOT, per author direction, describe the methodological reason (cophenet vs raw D as the discriminating operator) — that is Methods.

### Fix 4 — p-values attached *at the end* of plain-English claims

**Voice rule (author directive 2026-05-18).** Each statistical claim leads with the plain-English statement of what is true; the numerical evidence is appended at the close of the sentence. Do not lead with `p < 0.05`. Do not over-quote — one p-value per load-bearing claim is sufficient; counts (8/10, 7/10) carry the per-patient cohort picture as descriptive context.

**Numerical anchors the next draft must hit, with source CSV:**

| Plain-English claim | Numerical anchor (end of sentence) | Source CSV |
|---|---|---|
| Cohort median ρ_split^coph in trace direction | "+0.222; 8/10 patients individually positive; cohort-paired Wilcoxon p = 0.005" | `data/audit/ctm_triangle/cohort_summary.csv` |
| Same eight patients above drift-floor null | "8/10 above drift; cohort Wilcoxon p = 0.014" | same |
| Cross-probe restriction does not shift the cohort effect | "cohort median +0.223, 8/10 patients" | same (`rho_xprobe_median`) |
| ρ_split^coph above matched-strength surrogate | "+0.222 vs surrogate +0.009 (≈23× ratio); 7/10 above own surrogate; cohort Wilcoxon p = 0.005" | `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` |
| <!-- BH-FDR q = 0.027 over m = 6 bands row RETIRED 2026-05-20 per writing-agent feedback (cross-band BH-FDR dropped; per-band controls already gate). See feedback_no_unmotivated_bh_fdr.md + METHODS_AUDIT_ISSUES.md §B2. --> | | |
| LRG-native robustness restriction (drop Pat_15, the only LRG anti-aligned patient) preserves the trace | "n = 9, p = 0.010" (`ρ_split^coph` matched-strength); Grassmann `T_G(k=40)` strengthens to `p = 0.002`; drift-floor `p = 0.027` | `.agents/preprint/responses/2026-05-18_beta_n9_drop15_wilcoxon_reply.md` |
| Grassmann contiguous band clears matched-strength | "k ∈ [27, 55], cluster-extent p = 0.005" | `data/audit/grassmann_cluster_extent/cohort_summary.csv` |
| Peak Grassmann effect at k ≈ 40 | "observed cohort median ≈ 3× surrogate; 9/10 above own surrogate" | `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv` |
| Epi-zone exclusion preserves and broadens Grassmann | "every cell of the original window persists; the contiguous block extends to k ∈ [21, 56], 36 cells in total" | `data/audit/grassmann_epi_exclusion/cohort_summary.csv` |

**Do not invent numbers.** If a control was not run at β (e.g. C5 epi-X cophenet — audit_68 covered α only), do not cite a value. Silence is acceptable.

### Fix 5 — Number fidelity to source CSVs

- Cohort median `ρ_split^coph`: the brief carries +0.2222 (within-baseline, `ctm_triangle/cohort_summary.csv`) and +0.2206 (matched-strength, `matched_strength_surrogate_split_baseline/cohort_summary.csv`). Write **+0.222** uniformly for the within-baseline cohort median; **+0.221 vs surrogate +0.009** for the matched-strength contrast.
- Cross-probe cohort median: **+0.223** (`rho_xprobe_median` = +0.2228, rounded).
- Grassmann effect at k ≈ 40: actual obs/surr ratio = 2.94×, **not** "factor of three" strictly. Write "approximately three-fold" or "≈ 3×" — do not write "more than three-fold".
- Peak-ratio cell across the Grassmann band is k = 27 (ratio 3.53×); peak-magnitude cell is k = 45 (cohort median −0.585). The k ≈ 40 anchor in the current draft is fine as a representative mid-window cell; alternatively quote k = 45 for the magnitude peak. The agent may pick whichever reads cleaner — do not over-specify both.

---

## Voice and stability framing

### Stable-under-restriction framing (author directive 2026-05-18, point 2)
The cohort effect should be stated **without correction first**, then explicitly framed as stable across the sensitivity battery. The robustness checks (cross-probe restriction, LRG-native `n=9` dropping the single anti-aligned Pat_15) all return the cohort effect within a narrow band around +0.22 to +0.23. The agent should state this stability explicitly rather than recite each check as a separate caveat. Suggested closing clause for the per-pair sub-paragraph:

> "Across every robustness restriction — inter-shaft pairs only and the LRG-native `n=9` dropping the single anti-aligned patient — the cohort median sits in the same narrow band (≈ +0.22 to +0.23), so the per-pair trace is stable rather than carried by any one substructure or by any one patient."

**Pat_03 dropout is NOT a robustness check.** Pat_03 at 1024 Hz is fully integrated at the config layer (`nperseg_for_fs(fs)`); it is treated identically to every other patient at the analysis layer. Do not reintroduce Pat_03 dropout as a robustness row, sensitivity check, or per-patient marker.

This is the load-bearing summary; it replaces a per-check enumeration with a single stability claim.

### Statistical reporting density (author directive 2026-05-18, point 1)
Cite p-values, but do not overstate. Each load-bearing claim gets one p-value at its close (table in Fix 4). BH-FDR q may go in a parenthetical or footnote. Per-patient z-scores go in the supplementary per-patient table (see Per-patient table below); do not list them inline.

### Persistence terminology — softened (author directive 2026-05-18, point 4)
The earlier strict rule against bare "persistence" is **relaxed**. When the subsection context makes the meaning unambiguous (i.e. this is the task-trace section and the reader knows from the Methods that we mean task-induced retention), "persistence" reads cleanly. The explicit trace / anchor / reset / emergent taxonomy is still required where the reader could otherwise read "anchor" — cross-phase taxonomy tables, mixed-band paragraphs, captions of cross-phase summary figures. For this β subsection, "persistence", "trace", and "task-induced reorganization that is retained into rsPost" can be used interchangeably; pick whichever flows better.

### "Inter-shaft" / "cross-probe" (author directive 2026-05-18, point 9)
Synonyms. Do not globally substitute one for the other. Use whichever reads naturally in each sentence.

### Methods-vs-results discipline (author directive 2026-05-18, point 10)
Three sentences in the current draft drift into methodology; trim to a single nominal phrase. The methodology lives in the Methods section.

| Current (methods drift) | Trimmed (results-only, pointer to Methods) |
|---|---|
| "we quantify persistence as the Spearman correlation `ρ_split` between the task-induced and rest-induced shifts of the cophenetic communication distance `D_coph`, each shift referenced against an independent half of `rsPre` so that no shared baseline can inflate the coupling by construction" | "the per-pair correlation `ρ_split^coph` between task- and rest-induced shifts of the cophenetic communication distance `D_coph`, constructed against independent `rsPre` halves (see Methods)" |
| "Reading the leading-`k` eigendirections of the Laplacian — the slow, large-scale directions along which diffusion is most coherent — as a `k`-dimensional subspace, we ask whether the `rsPost` subspace sits closer to the `taskT` subspace than the `rsPre` subspace does, and we scan `k` across the full resolvable spectrum" | "for the global subspace probe — the leading-`k` Laplacian eigendirections, scanned across the resolvable spectrum (see Methods) —" |
| "against a matched-strength surrogate ensemble that preserves each node's strength exactly while randomizing which pair carries the weight" | "against the matched-strength surrogate (see Methods)" |

The "by construction" / "to eliminate" / "we ask whether" / "we scan" justification clauses are Methods-side framings; trim them.

---

## Additions

### Per-patient summary table (new artefact, agent to construct)

**Author directive 2026-05-18, point 5.** The paragraph names individual patients (Pat_03, Pat_07, Pat_15) but the reader has no per-patient context for the cohort claim. A manuscript-grade per-patient table should accompany the paragraph (decision on main text vs supplementary deferred — agent to recommend).

**Required columns:**

| Patient | n_chan | ρ_split^coph (obs) | surrogate median | z vs surrogate | sign at Grassmann k = 40 | substrate `T_d^(d_S)` sign | Notes |
|---|---|---|---|---|---|---|---|
| Pat_02 | … | … | … | … | pro / anti | pro / anti | β cohort backbone |
| Pat_03 | … | … | … | … | … | … | acquired at 1024 Hz, handled at config layer only (no special treatment at analysis layer); pro at LRG |
| Pat_05 | … | … | … | … | … | … | β cohort backbone |
| Pat_06 | … | … | … | … | … | … | β cohort backbone |
| Pat_07 | … | … | … | … | … | … | **anti at substrate only, pro at LRG (z > 3)** |
| Pat_08 | … | … | … | … | … | … | β cohort backbone |
| Pat_10 | … | … | … | … | … | … | task rows [53, 54, 55] dropped at load |
| Pat_13 | … | … | … | … | … | … | heaviest epi burden (30/119); trace-positive at epi-X manuscript window |
| Pat_14 | … | … | … | … | … | … | task_test vendor-replaced 2026-04-25 |
| Pat_15 | … | … | … | … | … | … | **right-hemi-only implant, 0 epi; LRG anti-aligned anchor** |

**Data sources for the table:**
- `ρ_split^coph` per patient + surrogate median + z: `data/audit/matched_strength_surrogate_split_baseline/per_patient.csv` (or equivalent per-patient artefact; verify with `audit_63_split_baseline_surrogate.py` output structure).
- Grassmann sign at k = 40: `data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv`.
- Substrate `T_d^(d_S)`: `data/audit/lrg_phase_distance/Td_per_patient_per_band.csv`.

**Anchor in the paragraph.** Add a single reference sentence at the close of the per-pair sub-paragraph: "Per-patient values are listed in Table S?." Same for Grassmann.

### Anatomy paragraph — **RED FLAG, DO NOT INLINE NUMBERS** (author directive 2026-05-18, point 8)

**Status.** The anatomical-distribution paragraph is currently a single phrase in the introduction sentence ("we do not observe ... a specific localization"). The brief at `.agents/preprint/bands/01_beta.md` §5 documents uncorrected enrichment leads (Hippocampus 3.49×, left fusiform 2.48×, left superior temporal 2.13×) and a cohort B_hemi correlation (ρ_S = +0.697, p = 0.025 uncorrected). None of these survive multi-region or BH-FDR correction.

**Author concern.** The user is not confident that the anatomical tests have been re-verified against the matched-strength gate, or re-run on the current n = 10 cohort. The anatomy section is therefore **deferred** for this draft.

**Agent action.** Insert a LaTeX placeholder block where the anatomy paragraph would otherwise sit, formatted as a clearly visible red flag in the source so it is impossible to miss in review:

```latex
% =====================================================================
% RED FLAG — ANATOMY PARAGRAPH DEFERRED (2026-05-18)
% ---------------------------------------------------------------------
% The anatomical-distribution paragraph for the β band is intentionally
% NOT included in this draft. The underlying per-region Desikan–Killiany
% enrichment + cohort B_hemi correlation analyses have not been
% re-verified against the locked matched-strength gate on n = 10, and
% the leads documented in .agents/preprint/bands/01_beta.md §5 (Hippocampus
% 3.49x, left fusiform 2.48x, B_hemi rho_S = +0.697) do not survive
% Bonferroni / BH-FDR. DO NOT inline anatomy numbers here until the
% user confirms the analysis is current and stable.
%
% Source-of-truth references:
%   .agents/preprint/bands/01_beta.md §5
%   data/audit/lrg_localization_anatomy/cohort_summary.csv
%   data/audit/implant_geometry/cohort_correlations.csv
% =====================================================================
```

Once the user confirms that the anatomy analyses are current and stable, the agent will be able to draft a single sentence along the lines of:

> "The trace is anatomically diffuse: no single Desikan–Killiany region carries the cohort signal under multi-region correction; uncorrected enrichment is seen at the Hippocampus and the left fusiform, with a suggestive cohort correlation with left-hemisphere implant fraction that does not survive correction over the implant-geometry contrasts."

Until then: red-flag block, no inlined numbers.

---

## Section structure (preserve the 4-paragraph layout, add Table S? anchor)

The current four-sub-paragraph layout is sound:

1. **Lead claim + β-uniqueness + diffuseness flag** — add the β-uniqueness sentence (Fix 3); the diffuseness sentence stays as a single phrase pointing forward to the anatomy paragraph that will be added once verified.
2. **`ρ_split^coph` results** — per-pair multiscale measure, drift-floor + cross-probe + matched-strength + robustness, closing with the stability framing.
3. **Interpretation of the per-pair finding** — keep largely as written; this is the cleanest sub-paragraph in the draft.
4. **Grassmann results + epi-X** — add the cluster-extent permutation p (Fix 2); keep the multiscale-block framing and the physiological-cortex close.

Add a single anchor sentence at the close of paragraphs 2 and 4 pointing to the per-patient Table S?.

Anatomy paragraph: deferred, red-flag block only.

---

## Source-of-truth references

- `.agents/preprint/bands/01_beta.md` — full β brief; every numerical claim in the paragraph must trace to a row here.
- `.agents/preprint/locked/VERDICT_LEDGER.md` — locked verdict; do not re-derive.
- `.agents/preprint/locked/CONTROLS.md` — locked 5-control battery.
- `.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md` — binding methods directive (KC / VI(k) / τ-sweep retired; D_coph adopted).
- `data/audit/ctm_triangle/cohort_summary.csv` — C1 split, C2 drift, C4 cross-probe.
- `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` — C3 matched-strength for `ρ_split^coph`.
- `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv` — Grassmann full-FC per-k.
- `data/audit/grassmann_cluster_extent/cohort_summary.csv` — audit_70 cluster-extent permutation null (β cluster p = 0.005).
- `data/audit/grassmann_epi_exclusion/cohort_summary.csv` — Grassmann epi-X (β: 29/29 persist + 7 emerge → 36-cell window).

## Revision history

- 2026-05-18 — Initial directive from author review of the drafted β paragraph. Five mandatory fixes (Pat_07 attribution, audit_70 cluster-extent p, β-uniqueness sentence, p-values at end of claims, number fidelity to source CSVs). Per-patient table requested. Anatomy paragraph deferred with red-flag block. Persistence terminology rule softened (now permitted when context is unambiguous).
- 2026-05-18 (later that day) — Pat_03 outlier/dropout framing **fully retired** on user direction. Pat_03 is now a full cohort member treated identically to every other patient at the analysis layer; sampling-rate difference (1024 Hz vs 2048 Hz) is absorbed at the config layer only (`nperseg_for_fs(fs)`, `FS_OVERRIDES`). All Pat_03-dropout robustness rows, "Pat_03-fragile" framing in anatomy, "Pat_03 (1024 Hz outlier)" tags, and figure marker distinctions are dropped. Rule landed in `CLAUDE.md`, `AGENTS.md`, `.agents/guides/04_rules/never-always-list.md`, `.agents/guides/01_project/agent-playbook.md`, `.agents/guides/03_implementation/data-layout.md`.
