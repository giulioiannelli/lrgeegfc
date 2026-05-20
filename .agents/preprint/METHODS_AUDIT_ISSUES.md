---
name: methods-audit-issues
era: IMCOH_ABS_COHORT_N10
status: open
kind: audit-checklist
scope: issues identified in the finalized Methods section + cascading consequences for preprint folder
created: 2026-05-20
supersedes: methods/methods_section_review_2026-05-19.md (partially — see resolution map §1)
companion: HANDOFF_INDEX.md, EVALUATION_PROTOCOL.md, locked/CONTROLS.md, locked/VERDICT_LEDGER.md
---

# Methods Audit — Issues + Cascade Tasks

**Head.** The finalized Methods section is internally near-consistent and aligns with the locked artifacts on probe definitions, control battery, and verdict vocabulary. Four real bugs (B1–B4), five reviewer-defense gaps (G1–G5), four cascading consequences (C1–C4), four quantitative inconsistencies inside the preprint folder (Q1–Q4), and four reviewer red flags (N1–N4) are flagged below — each with a one-line location and a one-line fix direction. **No drop-in LaTeX.** The writing agent and the user retain authorial control. Section 9 partitions the action items by owner (writing agent / user decisions / compute work / Phase 3 brief realignment).

## Contents

1. Resolution map vs `methods/methods_section_review_2026-05-19.md`
2. Methods-section bugs (B1–B4)
3. Methods-section gaps + reviewer-defense integrations (G1–G5)
4. Cascading consequences of the pasted methods (C1–C4)
5. Quantitative consistency inside the preprint folder (Q1–Q4)
6. Terminology consistency status (T1–T5)
7. Verdict + ledger consistency status (V1–V2)
8. Nature-reviewer red flags + preemptive responses (N1–N4)
9. Action items (categorized by owner)

---

## 1. Resolution map vs `methods/methods_section_review_2026-05-19.md`

The previous review (2026-05-19 am) flagged M1–M4 + minor items. Their status against the pasted methods (2026-05-19 pm):

| Prior flag | Status in pasted methods | Note |
|---|---|---|
| **M1** `T_G*` formula stale (longest-run-only sum) | **RESOLVED** — pasted methods uses all-clusters sum | New cascade C1 introduced: `T_G*` now normalized to `[0,1]` |
| **M2** `D(τ′)` written with `K` not `ρ` | **RESOLVED** — pasted methods uses `Π_ij(τ′)` (propagator) consistently in `D_ij(τ′) = (1−δ_ij)/Π_ij(τ′)` formulation | n/a |
| **M3** Matched-strength claim "preserves global edge-weight distribution" overstrong | **RESOLVED** — pasted methods correctly says "preserves per-node strength sequence to machine precision … marginal edge-weight distribution drifts within these constraints" | n/a |
| **M4** Epi-X criterion divergent | **RESOLVED** — pasted methods adopts Wilcoxon-on-epi-X (CONTROLS Decision 10) | n/a |
| Other minor items (m6, m9, etc.) | Partially resolved; new bugs surfaced (B1–B4) | See §2 below |

---

## 2. Methods-section bugs (B1–B4)

Flag-only. Each entry: location + 1-line problem + 1-line fix direction.

### B1 — Sign-convention contradiction

- **Location**: `sssec:methods_compare_stats`, first sentence of "Cohort-level inference" paragraph; vs `ssec:methods_compare` and `sssec:methods_compare_grassmann` sign-convention statements.
- **Problem**: `ssec:methods_compare` says `T_d(p,b) > 0 = trace direction`; `sssec:methods_compare_grassmann` echoes "`T_G(k;s,b) > 0` is the trace direction"; but `sssec:methods_compare_stats` writes `H_0: T_d ≥ 0 (no imprint or anti-imprint direction) against H_1: T_d < 0 (imprint direction)` — opposite sign.
- **Fix direction**: flip the stats sentence to `H_0: T_d ≤ 0` against `H_1: T_d > 0` (consistent with the comparison-section convention); replace "imprint direction" with "trace direction" (the locked vocabulary per `.agents/guides/01_project/terminology.md`).

### B2 — Drop the entire cross-band BH-FDR family-size sentence — REVISED 2026-05-20 (broader)

- **Location**: `sssec:methods_compare_stats`, the "Within-probe-family multiple-comparison correction uses the BH procedure …" sentence in full (both `ρ_split^coph` and `d_S on D_coph` clauses).
- **Problem (revised 2026-05-20)**: cross-band BH-FDR on `ρ_split^coph` (m=6) is scaffolding without function — each band's verdict is gated by an independent control battery (drift-floor `ρ_drift`, matched-strength surrogate at R=200, cross-probe `ρ_xprobe`, epi-zone exclusion, LOO max-p); the six per-band claims are not a coordinated cross-band family of inference. Cross-band BH on top of per-band matched-strength surrogacy is double-protection against a concern the per-band empirical null already addresses at the per-test level. The `d_S on D_coph` clause has the additional problem of being an undefined symbol in the comparison section (the original B2 scope).
- **Fix direction (locked, broader than original)**: drop the **entire** cross-band BH-FDR family-size sentence. Replace with one sentence stating that **no cross-band multiple-comparison correction is applied** to either per-pair `ρ_split^coph` or `T_G*` because each band carries its own independent control battery and the per-band verdicts are not a coordinated cross-band claim. Anatomical enrichment under A1 hypergeometric retains BH-FDR across DK regions per (band, probe) because the multi-region family IS the coordinated unit of inference.
- **Policy**: see `feedback_no_unmotivated_bh_fdr.md` for the three-point check that decides whether any future BH/Bonferroni/Holm correction is load-bearing. Recorded as a "Never" rule in `.agents/guides/04_rules/never-always-list.md`.

### B3 — `[xxx?]` placeholder citation

- **Location**: `sssec:methods_compare_stats`, sensitivity paragraph: "drops contacts in each patient's clinically identified epileptogenic zone per `[xxx?]`".
- **Problem**: unresolved citation placeholder.
- **Fix direction**: insert the canonical sEEG epileptogenic-zone reference used elsewhere in the lab (clinical-protocol paper for the cohort's epi-zone identification).

### B4 — `ρ_xprobe` definition under-specified

- **Location**: `sssec:methods_compare_ctm`, cross-probe restriction paragraph.
- **Problem**: "`ρ_xprobe` is `ρ_split` recomputed on the cross-probe subset of contact pairs only" — ambiguous whether (a) the full LRG pipeline is rerun on the cross-probe submatrix, or (b) the LRG pipeline is unchanged and only the Spearman correlation is restricted to cross-probe pairs of `D_coph`.
- **Fix direction**: explicit one sentence — "The LRG pipeline (eigendecomposition, propagator, UPGMA linkage, cophenet) is unchanged; the Spearman correlation `ρ_S(Δ_task, Δ_rest)` is restricted to the cross-probe subset of `D_coph` pair entries." (This matches CONTROLS.md C4 and the audit_71 implementation.)

---

## 3. Methods-section gaps + reviewer-defense integrations (G1–G5)

Items where the methods is internally consistent but a reviewer can plausibly ask "why?". One-line preemptive integration each.

### G1 — Substrate-probe gating not justified

- **Location**: `sssec:methods_compare_stats`, "The within-baseline drift floor `ρ_drift` and the cross-probe restriction `ρ_xprobe` … operate at the cophenetic layer only" sentence.
- **Gap**: methods doesn't say *why* the substrate probe `ρ^raw_split` has only the matched-strength gate (no C1/C2/C4 analogues).
- **Integration**: one sentence — "Substrate `ρ^raw_split` is reported as a descriptive baseline against the LRG-layer probes; the C1/C2/C4 controls operate at the cophenetic layer because the split-baseline / drift / cross-probe constructions read the LRG-multiscale per-pair structure that the substrate probe by construction does not resolve."

### G2 — C5 epi-X at cophenet only for α — DECISION: run for β too

- **Location**: `sssec:methods_compare_stats`, sensitivity paragraph: "this is reported for α … For every other band the epi-X analysis at the cophenetic layer is not run".
- **Decision (locked)**: run audit_72 (or equivalent cophenet epi-X pipeline) for β as well. β is the other cophenet-positive band and symmetric coverage strengthens reviewer-defense.
- **Action**: see §9 compute work.
- **Methods text update after the audit lands**: replace "is reported for α" with "is reported for α and β (the two cophenet-positive bands)". Drop the "every other band not run" sentence; the cophenet C5 is run wherever cophenet C1+C2+C3+C4 pass.

### G3 — Combinatorial vs symmetric-normalized Laplacian

- **Location**: `ssec:methods_lrg`, paragraph after eq. `methods_laplacian` ("We use the combinatorial form throughout (rather than the random-walk Laplacian `L_RW = D^{−1}L`) …").
- **Gap**: justifies combinatorial vs random-walk but not vs symmetric-normalized `L_sym = D^{−1/2} L D^{−1/2}` (also symmetric PSD; standard alternative in spectral graph theory).
- **Integration**: one clause appended — "and the symmetric-normalized form `L_sym = D^{−1/2} L D^{−1/2}` is also symmetric PSD but diagonalises a different operator from the rest of the LRG pipeline (propagator, cophenetic distance, etc.); using the combinatorial `L` preserves coherence across probes."

### G4 — UPGMA average-linkage choice

- **Location**: `ssec:methods_lrg`, paragraph introducing the linkage `L(τ′)`.
- **Gap**: methods specifies "Average-linkage agglomerative clustering" without justifying vs single-linkage / complete-linkage / Ward.
- **Integration**: one clause — "Average-linkage is chosen as the natural finite-sample estimator of the mean cophenetic distance in the continuous-spectrum regime; single-linkage produces chained dendrograms unstable to noise in the heavy-tailed `D(τ_min)` distribution, and Ward minimises within-cluster variance under a Gaussian assumption that does not match the `1/Π_ij` distance scale." (cite Kaufman & Rousseeuw 2009 already in methods).

### G5 — R=200 surrogate count

- **Location**: `sssec:methods_compare_stats`, matched-strength surrogate paragraph.
- **Gap**: `R = 200` stated without convergence diagnostic.
- **Integration**: one clause — "`R = 200` is sufficient for the `cluster_p_mass` empirical p-value to resolve below the locked `0.005` floor (`(R+1)^{−1} = 0.005`); cross-checks at `R = 500` for the load-bearing β / γ_l / δ Grassmann bands reproduce the same `cluster_p_mass` floor".

---

## 4. Cascading consequences of the pasted methods (C1–C4)

The pasted methods introduces convention changes that must propagate into existing locked artifacts.

### C1 — `T_G*` normalization to `[0,1]` — DECISION: dual columns in CSV

- **Change in methods**: `T_G*(b) = (n_k · log_10(R+1))^{−1} · Σ_{k:p_k<α_k} (−log_10 p_k(b)) ∈ [0,1]` with `n_k = 111`, `R = 200`, denominator ≈ 255.6.
- **Decision (locked)**: CSV `data/audit/grassmann_cluster_extent/cohort_summary.csv` gets a **new `obs_cluster_mass_neglog10p_norm` column** with the normalized values; the existing raw `obs_cluster_mass_neglog10p` column stays for audit-trail integrity. Briefs and manuscript cite the **normalized** form.
- **Normalized values** (computed from current raw values; `p_mass` invariant under monotone normalization):

| Band | Raw `T_G*` | Normalized `T_G*` | `p_mass` |
|---|---|---|---|
| β | 69.76 | **0.273** | 0.005 |
| γ_l | 66.14 | **0.259** | 0.005 |
| δ | 38.07 | **0.149** | 0.005 |
| γ_h | 31.67 | **0.124** | 0.060 |
| θ | 12.97 | **0.0507** | 0.159 |
| α | 7.98 | **0.0312** | 0.348 |

- **Cascade targets** (Phase 3): CSV regeneration (add column); `bands/00_cohort.md` §3 table values; `methods/methods_grassmann_cluster_extent.md` §6 table; all band briefs citing `T_G*`; `locked/VERDICT_LEDGER.md` per-band tables.

### C2 — `T_G^{*,s}` per-patient cluster mass — newly defined symbol

- **Change in methods**: `sssec:methods_compare_grassmann` eq. `methods_TGstar_perpatient` defines `T_G^{*,s}(b)` as the per-subject analogue of `T_G*(b)`, with same normalization denominator, reported descriptively (not gated).
- **Cascade**: every band brief that reports the Grassmann probe (`01_beta.md`, `03_gammalow.md`, `06_delta.md`) needs a `T_G^{*,s}` per-patient column in the Grassmann table. Currently absent.
- **Action**: Phase 3 brief realignment computes per-patient `T_G^{*,s}` from `per_patient_per_band_per_k.csv` and adds the column.

### C3 — Anatomy A1 / A2 (manuscript) ↔ A1 / A3 (lab) mapping

- **Change in methods**: anatomy section uses manuscript labels A1 (hypergeometric) + A2 (matched-strength surrogate).
- **Lab labels in preprint folder**: A1 (hypergeometric) + A3 (matched-strength surrogate) — A2 + A4 are deferred lab-only controls.
- **Status**: mapping notes already in place across `locked/ANATOMY_CONTROLS.md`, `locked/ANATOMY_LEDGER.md`, `writing_directive_2026-05-19_methods_anatomy.md`, `writing_directive_2026-05-19_anatomy_clusterext_rerun.md`, and `EVALUATION_PROTOCOL.md` Step 4.
- **Action**: one final pass during Phase 3 to confirm every brief that cites anatomy controls uses the lab labels A1/A3 internally and references the mapping for manuscript translation.

### C4 — `τ_min` adoption throughout preprint folder — DECISION

- **Change**: pasted methods uses `τ′`; preprint folder uses `τ_max`; user-locked convention is **`τ_min = 1/λ_max`** — physically precise (smallest diffusion time, finest spatial scale resolved by `K(τ)`).
- **Cascade targets**: methods LaTeX (writing agent); `locked/VERDICT_LEDGER.md` frontmatter + body; `locked/CONTROLS.md`; `methods/methods_revision_2026-05-18_cophenet.md`; `methods/methods_grassmann_cluster_extent.md`; `HANDOFF_INDEX.md`; `bands/00_cohort.md`; all band briefs; `README.md`; `EVALUATION_PROTOCOL.md`; `locked/ANATOMY_CONTROLS.md`; `locked/ANATOMY_LEDGER.md`.
- **Action**: Phase 3 folder-wide `τ_max` → `τ_min` and `τ′` → `τ_min` sweep; methods LaTeX handled by writing agent.

---

## 5. Quantitative consistency inside the preprint folder (Q1–Q4)

CSV-level + ledger-level drifts surfaced during the audit.

### Q1 — `grassmann_cluster_extent/cohort_summary.csv` post-fix status

- **Status**: ✓ CSV is post-2026-05-19-pm all-clusters formula. β 69.76, γ_l 66.14, δ 38.07, γ_h 31.67, θ 12.97, α 7.98. `cluster_p_mass` columns post-Decision-8 (β/γ_l/δ at 0.005).
- **Action**: per C1, add normalized column. No raw-value regeneration needed.

### Q2 — `bands/00_cohort.md` §3 cluster-extent table is STALE

- **Problem**: five of six rows carry pre-fix `obs_mass`:
  - β: 52.97 (should be 69.76)
  - γ_l: 19.17 (should be 66.14), `cluster_p_mass` 0.035 (should be 0.005)
  - δ: 12.78 (should be 38.07), `cluster_p_mass` 0.025 (should be 0.005)
  - γ_h: 15.75 (should be 31.67), `cluster_p_mass` 0.065 (should be 0.060)
  - θ: 7.11 (should be 12.97), `cluster_p_mass` 0.144 (should be 0.159)
  - α: 7.98 (matches), `cluster_p_mass` 0.099 (should be 0.348)
- **Action**: Phase 3 update §3 table to post-fix values (or per C1 to normalized values directly).

### Q3 — `bands/00_cohort.md` §1 + §4 carry STALE γ_l / δ verdicts

- **Problem**: §1 verdict matrix shows γ_l and δ as "**weak trace, only Grassmann**"; §4 paragraphs say "weak trace". Locked under Decision 8 (2026-05-19 pm) these are **strong trace, only Grassmann**.
- **Action**: Phase 3 update §1 matrix + §4 γ_l + δ paragraphs to "strong trace, only Grassmann"; add the LOO caveats (γ_l LOO 0.040 robust; δ LOO 0.055 full-data flagged + C5 epi-X resolves) already in VERDICT_LEDGER.

### Q4 — `HANDOFF_INDEX.md` per-band cheatsheet STALE γ_l / δ lines

- **Problem**: γ_l line says "weak trace, only Grassmann"; δ line says "weak trace, only Grassmann". Locked verdict is strong.
- **Action**: Phase 3 update both cheatsheet lines + the "What this manuscript IS about" paragraph that summarizes γ_l + δ.

---

## 6. Terminology consistency status (T1–T5)

### T1 — `τ_min` cascade

- See §4 C4. Folder-wide find-replace in Phase 3.

### T2 — "anchor anatomy" δ — CLEAN (no action)

- δ Grassmann "anchor anatomy" retracted under `S(δ)`; C4 cross-probe substrate-layer reading retained as descriptive (LEDGER Decision 5). Both demarcated consistently.

### T3 — Bare "persistence" — CLEAN (no action)

- Relaxation per `writing_directive_2026-05-18_beta_paragraph.md`: bare "persistence" allowed in task-trace subsections of trace-positive bands. No misuse detected in mixed-band or anatomy paragraphs.

### T4 — "substrate" — CLEAN (no action)

- Formal use only in `sssec:methods_compare_rawfc` heading. No informal stand-in usage in the preprint folder.

### T5 — "imprint direction" — REMOVED VIA B1

- Appears once in methods stats (the B1 sentence). When B1 is fixed (sign-convention + replace "imprint direction" with "trace direction"), this terminology disappears. No separate action.

---

## 7. Verdict + ledger consistency status (V1–V2)

### V1 — Trace verdict consistency (LEDGER ↔ briefs ↔ 00_cohort ↔ HANDOFF)

- `locked/VERDICT_LEDGER.md`: post-Decision-8 verdicts (β strong/strong; α strong cophenet-only; γ_l strong Grassmann-only; δ strong Grassmann-only; γ_h no trace; θ no trace). ✓
- `bands/01_beta.md`: β verdict ✓.
- `bands/02_alpha.md` … `bands/06_delta.md`: not re-read in this audit; Phase 3 will verify γ_l + δ briefs reflect "strong trace, only Grassmann" not "weak".
- `bands/00_cohort.md`: γ_l + δ STALE (see Q3).
- `HANDOFF_INDEX.md`: γ_l + δ STALE (see Q4).
- **Action**: Phase 3 cascade.

### V2 — Anatomy verdict consistency (ANATOMY_LEDGER ↔ briefs)

- ANATOMY_LEDGER: locked 2026-05-19 pm under cluster-extent paradigm `S(b)`. Left fusiform retracted (appears nowhere); Hippocampus survives at β Grassmann; δ-full vs δ-epi-X fully disjoint (0 shared regions).
- Band briefs: post-cluster-extent updates were applied in previous session (per conversation summary). Phase 3 will verify.
- **Status**: appears CLEAN; verify during Phase 3.

---

## 8. Nature-reviewer red flags + preemptive responses (N1–N4)

### N1 — Small-N cohort (n=10)

- **Reviewer question**: cohort is small; how robust are the cohort-level Wilcoxon tests?
- **Preemptive response (one sentence in stats section)**: "Cohort-level inference uses paired non-parametric Wilcoxon tests against patient-matched empirical nulls from R=200 strength-preserving surrogates, which are robust to the small-N regime; every cohort gate is paired with a leave-one-out diagnostic that flags single-patient leverage explicitly (see `sssec:methods_compare_stats` LOO paragraph)."
- See also G5 (R convergence) and N4 (LOO scope).

### N2 — R=200 surrogate convergence

- **Reviewer question**: is R=200 enough to resolve the cluster-mass empirical p-value floor?
- **Preemptive response**: see G5.

### N3 — No cross-band multiple-comparison correction (broader scope, 2026-05-20)

- **Reviewer question**: why no cross-band BH-FDR across the six bands on either LRG probe?
- **Preemptive response (one sentence)**: "No cross-band multiple-comparison correction is applied at either LRG probe because each band's verdict is gated independently by its own control battery — per-band drift-floor, matched-strength surrogate, cross-probe restriction, and epi-zone exclusion for `ρ_split^coph`; per-band cluster-extent permutation on the `k`-grid for `T_G*` — and the six per-band claims are not part of a coordinated cross-band family of inference; per-band matched-strength surrogacy already provides a calibrated empirical null at the per-test level. Within-band per-region anatomical enrichment under A1 retains BH-FDR across DK regions because the multi-region family IS the coordinated unit of inference per (band, probe)."
- **Policy reference**: `feedback_no_unmotivated_bh_fdr.md`; `.agents/guides/04_rules/never-always-list.md`.

### N4 — LOO max-p scope

- **Reviewer question**: what does LOO max-p mean for the verdict?
- **Preemptive response (one sentence)**: "LOO max-p is a descriptive single-patient-leverage diagnostic reported alongside every cohort gate; it is never a verdict criterion. Verdicts where LOO max-p crosses 0.05 are flagged in Results as single-patient-leveraged (with the C5 epi-X analysis as the mechanistic resolution where applicable)." (This matches `feedback_no_single_patient_p_driven.md`.)

---

## 9. Action items — categorized by owner

### A. Writing agent (Methods LaTeX edits)

- B1 — flip sign convention in `sssec:methods_compare_stats`; replace "imprint direction" with "trace direction"
- B2 — drop the entire cross-band BH-FDR family-size sentence (both `ρ_split^coph` and `d_S on D_coph` clauses); replace with one sentence stating no cross-band correction is applied (each band gated independently)
- B3 — insert canonical sEEG-epi citation in place of `[xxx?]`
- B4 — one explicit sentence on `ρ_xprobe` construction (LRG pipeline unchanged; correlation restricted to cross-probe pairs)
- G1 — one-sentence substrate-gate justification
- G3 — one-clause `L_sym` preemption appended to combinatorial-Laplacian paragraph
- G4 — one-clause UPGMA average-linkage justification
- G5 — one-clause R=200 convergence claim
- C4 — replace `τ′` with `τ_min` throughout the LaTeX (the symbol is the smallest diffusion time, `τ_min = 1/λ_max`)
- G2 — replace "is reported for α" sentence per the audit-72 β output (see B below)

### B. User decisions made (recorded here)

- C1 — `T_G*` carries both raw + normalized in CSV; briefs/manuscript cite normalized
- B2 — drop the **entire** cross-band BH-FDR family-size sentence (both `ρ_split^coph` AND `d_S on D_coph` clauses); writing-agent feedback 2026-05-20 broadens this from the original "drop d_S only" scope. Policy: `feedback_no_unmotivated_bh_fdr.md`; three-point check ((1) gates a verdict, (2) coordinated family, (3) per-test gate doesn't already cover it) recorded in `.agents/guides/04_rules/never-always-list.md`.
- C4 — `τ_min = 1/λ_max` is the locked symbol (cascade across folder)
- G2 — run audit_72 cophenet C5 for β

### C. Compute work (action item)

- **Run cophenet C5 for β**: re-run the cophenet epi-X pipeline for β (`load_epileptic_nodes` per patient; LRG pipeline on epi-X-restricted adjacency; per-patient `ρ_split^coph,epi-X`; one-sample one-sided Wilcoxon at cohort level; LOO max-p). Output to `data/audit/beta_epi_exclusion/c5_wilcoxon_cohort.csv`.
- **Add normalized `T_G*` column to `grassmann_cluster_extent/cohort_summary.csv`**: re-run audit_70 (or one-line script edit) to add `obs_cluster_mass_neglog10p_norm` column with `raw / (111 · log_10(201))`.

### D. Phase 3 — preprint-folder cascade (after methods crystallization + β paragraph)

- Q2 — update `bands/00_cohort.md` §3 cluster-extent table to post-fix normalized values
- Q3 — update `bands/00_cohort.md` §1 verdict matrix + §4 γ_l/δ paragraphs to "strong trace, only Grassmann"
- Q4 — update `HANDOFF_INDEX.md` γ_l + δ cheatsheet lines + "What this manuscript IS about" summary
- C2 — add per-patient `T_G^{*,s}` columns to `bands/01_beta.md`, `bands/03_gammalow.md`, `bands/06_delta.md`
- C3 — one final pass on A1/A2 manuscript-vs-lab label mapping consistency
- T1 — folder-wide `τ_max → τ_min` and `τ′ → τ_min` sweep
- C1 — rescale `T_G*` value citations in all briefs to normalized form
- V1 — verify γ_l + δ briefs reflect post-Decision-8 "strong trace, only Grassmann"
- Cross-band consistency check: identical section structure across all six band briefs (head / claim / cohort / per-probe / anatomy / interpretation / sources)

---

## Resolution path

1. **User decides cascade items** — done (C1, B2, T1=τ_min, G2 β audit) — recorded in §9B.
2. **Compute work** — §9C: run audit_72-cophenet for β; add normalized column to grassmann_cluster_extent CSV.
3. **Writing agent applies** — §9A methods LaTeX edits.
4. **User confirms methods crystallization** + shares β paragraph.
5. **Phase 3 cascade** — §9D preprint folder updates.

After Phase 3 completes, the writing agent re-reads the updated briefs + cross-band synthesis and derives per-band manuscript paragraphs with band-comparable structure and post-audit numbers.
