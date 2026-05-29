---
name: preprint-evaluation-protocol
era: IMCOH_ABS_COHORT_N10
status: locked_2026-05-19
kind: writing-agent-output-evaluation-protocol
companion: HANDOFF_INDEX.md
---

# Evaluation protocol — checking writing-agent LaTeX output (locked 2026-05-19)

**Head.** When the user pastes LaTeX from the writing agent into a Claude Code session, run this protocol in order. Each step has a specific check, a specific failure mode, and a specific corrective action. Report findings as a side-by-side delta vs the briefs: cite the LaTeX line / paragraph + the ledger or brief entry that contradicts it. Do **not** silently fix issues — flag them so the user can decide whether to send the LaTeX back to the writing agent or to override.

## Read order for evaluation

1. The pasted LaTeX (the input).
2. `HANDOFF_INDEX.md` (writing-agent contract).
3. `locked/VERDICT_LEDGER.md` + `locked/ANATOMY_LEDGER.md` (locked verdicts).
4. `locked/CONTROLS.md` + `locked/ANATOMY_CONTROLS.md` (locked batteries).
5. Per-band briefs referenced in the LaTeX section being checked.

Each step below references which artifact is authoritative for the check.

## Step 1 — Number-level CSV cross-check (most important)

**Goal**: every numerical claim in the LaTeX must trace back to a CSV row cited in a per-band brief.

**Procedure**:
- Extract every number from the LaTeX (Wilcoxon p, ratio, n_above, cluster_p, region enrichment, k-window bounds, etc.).
- For each, find the originating per-band brief (`bands/01_beta.md` … `bands/06_delta.md`) and the cited CSV row.
- Compare the LaTeX number to the brief's number to the CSV value.

**Failure modes**:
- LaTeX number disagrees with brief (typo / drift).
- LaTeX cites a number not present in any brief (writing agent invented it OR pulled it from memory / older artifact).
- LaTeX number rounded to fewer digits than the brief.

**Authoritative source**: brief tables → CSVs in `data/audit/*/cohort_summary.csv`.

**Reporting**: "LaTeX claims X = a.bc; brief says X = a.bd (`01_beta.md:LL`); CSV says X = a.bd1 (`<csv>:row`). Action: correct LaTeX to a.bd or a.bd1."

## Step 2 — Object name compliance

**Goal**: every per-pair distance citation is `D_coph` (or its long form); every subspace citation is on `U_k`; every τ is `τ_max = 1/λ_max`.

**Procedure**: grep the LaTeX for:
- `D(\tau)`, `D(τ)`, "propagator distance", "raw D" — if present in Results or Methods sections, must be in the **substrate / raw D layer** of the three-layer cohort table only. If in the main text describing the primary LRG probe, replace with `D_coph` or `D_coph = cophenet(UPGMA(D(τ_max)))`.
- "ultrametric matrix", "linkage", "Z" — if present, must reference `D_coph` correctly.
- Grassmann distance citations — must specify `U_k` and the k-range.
- "diffusion time", "LRG time" without `τ_max = 1/λ_max` — replace.

**Failure mode**: writing agent confuses `D(τ)` (raw propagator distance) with `D_coph` (cophenetic ultrametric) — see `feedback_never_confuse_D_with_cophenet.md`. The cached field `lrg.ultrametric_matrix` IS `D_coph`, not raw `D(τ)`.

**Reporting**: "LaTeX para X writes `D(τ)` describing the per-pair LRG probe — should be `D_coph`. The raw propagator distance only appears in the three-layer cohort table as the substrate layer."

## Step 3 — Verdict consistency

**Goal**: every per-band verdict in the LaTeX matches the locked verdict in `locked/VERDICT_LEDGER.md` + `locked/ANATOMY_LEDGER.md`.

**Procedure**: for each band, locate the LaTeX verdict statement (e.g., "β shows a strong trace at both probes"). Compare to:
- Trace verdict from `locked/VERDICT_LEDGER.md` § (one of: strong trace, both probes / strong trace, only D_coph / weak trace, only Grassmann / no trace).
- Anatomy verdict from `locked/ANATOMY_LEDGER.md` § (one of: strong localized / weak localized / not localized).

**Failure modes**:
- LaTeX upgrades a "weak" verdict to "strong" or vice versa.
- LaTeX writes "strong trace, both probes" for α (which is "strong trace, only D_coph").
- LaTeX writes "diffuse, not localized" for β (which is "strong localized, both probes" — retired KC-era framing).
- LaTeX writes "Hippocampus + left fusiform" for β anatomy (retracted; under the locked cluster-extent paradigm `S(b)`, left fusiform appears at **none** of β / γ_l / δ).

**Reporting**: "LaTeX claims X for band Y; ledger says Z. Reference: `VERDICT_LEDGER.md:LL` or `ANATOMY_LEDGER.md:LL`."

## Step 4 — Control battery citation

**Goal**: every cohort claim in the Results section cites which controls support it.

**Procedure**: for each cohort claim, check that at least the locked C3 matched-strength is cited (mandatory per `feedback_matched_strength_mandatory.md`). For D_coph traces, C1+C2+C4 should also be cited where relevant. For Grassmann traces, the cluster-extent permutation p-value (LR + mass) should be cited. For anatomy claims, **manuscript labels A1 (hypergeometric) + A2 (matched-strength surrogate; lab A3) must be cited**; for Grassmann-only anatomy, manuscript A2 alone is sufficient (with the caveat that A1 was sparse). Lab labels (A1 / A3) remain in the briefs and ledgers for traceability; the manuscript LaTeX uses A1 / A2 per the 2026-05-19 pm renumbering directive.

**Failure modes**:
- LaTeX writes "p < 0.05" without specifying which null.
- LaTeX cites only the cohort sign agreement (e.g., "7/10 patients in the same direction") without the matched-strength control.
- LaTeX cites old controls (Pat_03 dropout, coverage-matched permutation as a gate) — those are descriptive-only or retired.
- LaTeX cites Bonferroni at m=48 for anatomy without explaining that A2 (lab A3 matched-strength surrogate) is the gate, not Bonferroni.

**Reporting**: "LaTeX para X cites p=0.005 without naming the null. Per `locked/CONTROLS.md` C3, this is paired one-sided Wilcoxon vs R=200 matched-strength surrogates (`feedback_matched_strength_mandatory.md`). Action: add 'matched-strength' or cite C3 explicitly."

## Step 5 — Anti-pattern scan (full `HANDOFF_INDEX.md` checklist)

Grep the LaTeX for every anti-pattern in `HANDOFF_INDEX.md` "Anti-pattern checklist". Report any match.

**Specific patterns to grep**:
- "diffuse" near "β" or "beta" → retracted (β anatomy is strong localized, 7+7 regions).
- "Hippocampus + left fusiform" → retracted (Hippocampus survives at β Grassmann; left fusiform appears nowhere under the cluster-extent paradigm `S(b)`).
- "KC" / "T_KC" / "Kuhner-Felsenstein" → KC is retired 2026-05-18.
- "VI(k)" / "variation of information" → VI is retired 2026-05-18.
- "τ-sweep" / "tau sweep" / "across multiple diffusion times" — retired 2026-05-18; LRG spectrum is gap-less; `τ_max = 1/λ_max` is the canonical scale.
- "Pat_03 dropout" / "Pat_03 outlier" / "Pat_03 acquired at 1024 Hz (excluded)" — retired 2026-05-18.
- "persistence" without "trace / anchor / reset / emergent" disambiguation in cross-phase contexts.
- "denoising" near "cophenet" — wrong; cophenet does band resolution.
- "brain-wide" / "uniformly distributed across cortex" for any trace-positive band — retracted; the trace is localized to a band-specific cortical network.
- "volume + topology" / "orthogonal" framing of `d_P` vs `d_S` — retracted.
- "n_below ≥ 8" / "frac_consistent ≥ 7/10" as a gate — retracted.

**Reporting**: list each match with line number + the retracting reference (memory or ledger).

## Step 6 — Figure-path coverage

**Goal**: every figure reference in the LaTeX corresponds to a file at the cited path.

**Procedure**: `ls` the cited paths under `data/preprint/figures/<band>/`. Check that PDF exists (note: per `feedback_no_png_duplicates.md`, only PDF should exist; if the LaTeX cites a `.png`, flag it).

**Failure mode**: writing agent assumes a figure exists that hasn't been generated yet; cites a path that's empty.

**Reporting**: "LaTeX figure ref \\includegraphics{data/preprint/figures/...} at `LL` — file does not exist. Generate figure first."

## Step 7 — Multiple-comparison correction consistency (revised 2026-05-20)

**Goal**: every BH-FDR (or Bonferroni / Holm / any correction) cited in the LaTeX must pass the three-point check from `feedback_no_unmotivated_bh_fdr.md` — corrected `p`/`q` gates a verdict, family is a coordinated unit of inference, per-test gate doesn't already address multiple-testing.

**Authoritative**:
- **Per-pair `ρ_split^coph`**: NO cross-band BH-FDR. Each band's verdict is gated independently by C1 / C2 / C3 / C4 / C5 from `locked/CONTROLS.md`; the six per-band claims are not a coordinated cross-band family. (Cross-band BH at `m = 6` was retired 2026-05-20 per writing-agent feedback.)
- **Grassmann cluster-extent**: NO cross-band BH (cluster-extent permutation is the family-level gate across the `k`-grid per band).
- **A1 anatomy hypergeometric**: BH-FDR across DK regions per (band, probe) IS the gate — multi-region family is the coordinated unit of inference per (band, probe). Family size `m ≈ N_DK_regions_with_coverage_per_band`. KEEP.
- **A4 implant-geometry regression (deferred)**: BH-FDR across covariates per (band, probe). KEEP if/when A4 is run.

**Failure modes**:
- LaTeX cites BH at `m=6` across bands on `ρ_split^coph` or any LRG probe (retired 2026-05-20).
- LaTeX cites BH at `m=48` for cohort claims (KC-era hardcoded; retired).
- LaTeX claims Bonferroni at any LRG-probe layer (controls are matched-strength surrogate + cluster-extent + A1+A3, none Bonferroni-based).
- LaTeX cites any correction without specifying (a) the corrected `p`/`q` is the gate, (b) the family is coordinated, (c) the per-test gate doesn't already absorb the concern.
- LaTeX drops anatomy A1 BH-FDR — that one IS decisive and must remain.

## Step 8 — Methodological-language scan (binding directive)

**Authoritative**: `methods/methods_revision_2026-05-18_cophenet.md` is the binding methods directive.

**Procedure**: grep the LaTeX Methods section for:
- The definition of `D_coph` (must be explicit: `D_coph = cophenet(UPGMA(D(τ_max)))`).
- The choice `τ_max = 1/λ_max` (stated as the canonical scale; reason: continuous spectrum, no gap).
- The cluster-extent permutation (LR + mass, disjunctive gate).
- The matched-strength algorithm (4-cycle ±δ, R=200, swap_target=20, seed=20260511).
- The substrate (`imcoh_abs = <|ImCoh|>_f`).

**Failure modes**:
- LaTeX writes `K̂(τ_max) = exp(−τ_max · L̂)` but does not state that `τ_max = 1/λ_max`.
- LaTeX writes about τ-sweep / multiscale via τ-sweep — retired.
- LaTeX writes "denoising" instead of "band resolution" for cophenet.
- LaTeX writes a different surrogate seed or different R.

## Step 9 — Consistency across paragraphs

**Goal**: the LaTeX doesn't contradict itself.

**Procedure**: cross-check:
- The Results verdict for each band matches the Discussion framing.
- The β anatomy in §Results matches the β anatomy in §Discussion.
- The cross-band synthesis paragraph (in Discussion) doesn't contradict the per-band Results paragraphs.

**Failure mode**: writing agent produces a slightly different framing in two places (e.g., "β trace is diffuse" in Discussion vs. "β trace localizes to a cortical network" in Results).

## Step 10 — Final delta report

Produce a single document at the end with:
1. **Number disagreements**: LaTeX value vs. brief value vs. CSV value, with line numbers.
2. **Object-naming corrections**: each `D(τ)` → `D_coph` correction; each missing `U_k` specification.
3. **Verdict contradictions**: each LaTeX verdict that contradicts `locked/VERDICT_LEDGER.md` or `locked/ANATOMY_LEDGER.md`.
4. **Control battery gaps**: each cohort claim without proper C3 + (where applicable) cluster-extent / A3 citation.
5. **Anti-pattern matches**: each anti-pattern hit with the retracting reference.
6. **Figure path gaps**: each cited figure path that doesn't exist.
7. **BH-FDR / methodological-language gaps**: each m-mismatch or wrong-vocabulary instance.
8. **Internal inconsistencies**: each paragraph-vs-paragraph disagreement.

End with: "Send back to writing agent with the above corrections" OR "LaTeX is consistent with locked ledgers and briefs; ready for user review."

## Reporting style

For every finding cite:
- LaTeX line / paragraph reference.
- Brief / ledger entry that contradicts it (file:line).
- Specific corrective action.

**Do not silently fix issues.** Flag them. The user decides whether to send back to writing agent or override.

**Do not propose new analyses.** If a verdict is contested by a finding, the answer is "the locked ledger stands; the writing agent must rephrase." Verdicts change only via a dated revision in the ledger, which requires a new audit.

## What this protocol explicitly does NOT do

- Does not rewrite the LaTeX. Output is a delta report, not a corrected manuscript.
- Does not re-run audits. The CSVs cited in briefs are the substrate.
- Does not propose new figures. If a figure path is empty, flag it; do not generate the figure as part of the eval pass.
- Does not invent verdicts. If the LaTeX makes a claim not in any ledger, the corrective action is "remove the claim", not "audit it".

## Revision history

- **2026-05-19** — Locked. Protocol covers the 5-control trace battery + 4-control anatomy battery + the disjunctive Grassmann cluster-extent gate. All anti-patterns from `HANDOFF_INDEX.md` are enumerated for the grep pass.
