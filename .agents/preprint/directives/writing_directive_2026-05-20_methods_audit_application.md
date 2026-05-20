---
name: writing-directive-methods-audit-application
era: IMCOH_ABS_COHORT_N10
status: open
kind: writing-directive
scope: methods-section LaTeX edits derived from the 2026-05-19/20 methods audit
created: 2026-05-20
companion: METHODS_AUDIT_ISSUES.md (master spec), methods/methods_revision_2026-05-18_cophenet.md, methods/methods_grassmann_cluster_extent.md
target: writing agent producing manuscript Methods section LaTeX
---

# Methods-section edits — writing-agent directive

**Head.** Apply nine edits to the Methods LaTeX. Four are bugs (B1–B4), five are reviewer-defense gaps (G1–G5), and four are convention changes (C1–C4). The audit doc `METHODS_AUDIT_ISSUES.md` is the master spec; this directive is its writing-agent-facing distillation. Do not re-derive verdicts. Do not propose new analyses. Apply the edits, return the revised LaTeX.

## What to read first (in order)

1. **This directive** — the edit list and the locked conventions.
2. **`METHODS_AUDIT_ISSUES.md`** — full audit context for each item; §9A is the writing-agent action block; §9B records the user decisions; §1–§8 are the rationale.
3. **`methods/methods_revision_2026-05-18_cophenet.md`** — binding methods directive for the cophenet wrap.
4. **`methods/methods_grassmann_cluster_extent.md`** — locked Grassmann probe methodology including the all-clusters `T_G*` formula and the cluster-extent permutation null.

## Locked conventions (do not deviate)

- `τ_min = 1/λ_max` is the canonical symbol throughout the LaTeX. Replace every `τ′`, `τ_max`, "diffusion time", "LRG time" with `τ_min`. Physical reasoning: `λ_max` is the largest Laplacian eigenvalue and `τ_min = 1/λ_max` is the smallest diffusion time, corresponding to the finest spatial scale resolved by `Π(τ)`.
- `T_G*(b)` reported in the manuscript is the **normalized** form `∈ [0,1]`: `T_G*(b) = (n_k · log_10(R+1))^{-1} · Σ_{k:p_k(b)<α_k} (−log_10 p_k(b))` with `n_k = 111`, `R = 200`, denominator `≈ 255.6`. Lab CSV carries both raw and normalized columns; cite normalized values.
- Anatomy controls in the manuscript use compressed labels **A1** (hypergeometric per-region) + **A2** (matched-strength surrogate). The lab-side label for matched-strength is A3; the deferred lab controls A2 (sampling-corrected) and A4 (implant-geometry) are NOT in the manuscript. The mapping is documented in `locked/ANATOMY_CONTROLS.md`.

## The nine edits

### Bugs (4) — fix these

**B1 — Sign-convention contradiction**
- Section: `sssec:methods_compare_stats`, first sentence of "Cohort-level inference".
- Current text: "under `H_0: T_d ≥ 0` (no imprint or anti-imprint direction) against `H_1: T_d < 0` (imprint direction)".
- Fix direction: flip to `H_0: T_d ≤ 0` against `H_1: T_d > 0`, and replace "imprint direction" with "trace direction". This matches the sign convention `T_d > 0 = trace` stated in `ssec:methods_compare` and `sssec:methods_compare_grassmann`.

**B2 — Drop the entire cross-band BH-FDR family-size sentence** (REVISED 2026-05-20 — broader than original)
- Section: same paragraph as B1.
- Current text: "Within-probe-family multiple-comparison correction uses the BH procedure at `q = 0.05`, with family sizes `m = 6` for the controlled per-pair multiscale correlation `ρ_split` (six bands) and `m = 6` for the matrix-level Spearman distance `d_S` on `D_coph` (six bands); the two families are corrected independently as they test distinct statistics on `D_coph` (controlled per-pair split-baseline correlation versus single-scale matrix-level distance)."
- Fix direction: **delete the entire sentence**. Replace with one sentence: "No cross-band multiple-comparison correction is applied at either LRG probe; each band's verdict is gated independently by its own control battery (drift-floor `ρ_drift`, matched-strength surrogate at R=200, cross-probe `ρ_xprobe`, epi-zone exclusion for `ρ_split^coph`; per-band cluster-extent permutation on the `k`-grid for `T_G*`), and the six per-band claims are not a coordinated cross-band family of inference. Within-band per-region anatomical enrichment under A1 (Section `sssec:methods_anatomy`) retains BH-FDR across DK regions per (band, probe) where the multi-region family is the coordinated unit of inference."
- **Rationale**: writing-agent feedback 2026-05-20 — cross-band BH on `ρ_split^coph` (m=6) added no verdict-level information because each band is gated by independent controls and the per-band matched-strength surrogate is already a calibrated empirical p; layering BH on top is scaffolding without function. The `d_S on D_coph` clause has the additional problem of being undefined in the comparison section (original B2 scope; subsumed here). Three-point check before any future correction: (1) does the corrected `p`/`q` gate a verdict, (2) is the family a coordinated unit of inference, (3) does the per-test gate already address the multiple-testing concern. If any check fails, the correction is scaffolding and should not appear.

**B3 — `[xxx?]` epileptogenic-zone citation**
- Section: `sssec:methods_compare_stats`, sensitivity paragraph: "drops contacts in each patient's clinically identified epileptogenic zone per `\cite{xxx?}`".
- Fix direction: insert the canonical sEEG-epi reference used elsewhere in the lab's prior work for the cohort's epileptogenic-zone identification protocol.

**B4 — `ρ_xprobe` definition explicit**
- Section: `sssec:methods_compare_ctm`, cross-probe restriction paragraph.
- Current text: "`ρ_xprobe(p, b)` is `ρ_split` recomputed on the cross-probe subset of contact pairs only".
- Fix direction: one explicit sentence — "The LRG pipeline (eigendecomposition, propagator at `τ_min`, average-linkage UPGMA, cophenet) is unchanged; the Spearman correlation `ρ_S(Δ_task, Δ_rest)` is restricted to cross-probe entries of `D_coph`. This matches the audit_71 cohort implementation."

### Gaps (5) — preempt reviewers

**G1 — Substrate-probe gating justification**
- Section: `sssec:methods_compare_stats`, after the sentence "The within-baseline drift floor `ρ_drift` and the cross-probe restriction `ρ_xprobe` … operate at the cophenetic layer only."
- Fix direction: append one sentence — "Substrate `ρ^raw_split` is reported as a descriptive baseline against the LRG-layer probes; the C1/C2/C4 controls operate at the cophenetic layer because the split-baseline / drift / cross-probe constructions read the LRG-multiscale per-pair structure that the substrate probe by construction does not resolve."

**G2 — C5 epi-X cophenet for β (audit-driven)**
- Status: a new audit run (cophenet C5 for β) is queued; CSV will land at `data/audit/beta_epi_exclusion/c5_wilcoxon_cohort.csv`. **Wait for the CSV before applying this edit.**
- Section: `sssec:methods_compare_stats`, sensitivity paragraph.
- Current text: "is reported for α, where it is the band that the cophenet probe identifies as a non-epileptic-cortex trace partly masked by epi-zone contacts in the full-cohort analysis. For every other band the epi-X analysis at the cophenetic layer is not run, and the verdict from controls 1–4 stands."
- Fix direction (apply once the β CSV lands): "is reported for the two cophenet-positive bands α and β. For all other bands the cophenet probe has no positive trace, so an epi-zone-exclusion sensitivity is not informative; the verdict from controls 1–4 stands." Cite both `alpha_epi_exclusion/c5_wilcoxon_cohort.csv` and `beta_epi_exclusion/c5_wilcoxon_cohort.csv`.

**G3 — Combinatorial vs symmetric-normalized Laplacian**
- Section: `ssec:methods_lrg`, paragraph after eq. `methods_laplacian`.
- Current text justifies the combinatorial form vs the random-walk Laplacian.
- Fix direction: append one clause — "and the symmetric-normalized form `L_sym = D^{−1/2} L D^{−1/2}` is also symmetric PSD but diagonalises a different operator from the rest of the LRG pipeline (propagator, cophenetic distance), so the combinatorial `L` is used throughout to preserve cross-probe coherence."

**G4 — UPGMA average-linkage choice**
- Section: `ssec:methods_lrg`, paragraph introducing the linkage on `D(τ_min)`.
- Fix direction: append one clause to the "Average-linkage agglomerative clustering" sentence — "Average-linkage is chosen as the natural finite-sample estimator of mean cophenetic distance in the continuous-spectrum regime; single-linkage produces chained dendrograms unstable to noise in the heavy-tailed `D(τ_min)` distribution, and Ward minimises a within-cluster variance under a Gaussian assumption that does not match the `1/Π_ij` distance scale (Kaufman & Rousseeuw 2009)."

**G5 — R=200 surrogate convergence claim**
- Section: `sssec:methods_compare_stats`, matched-strength surrogate paragraph.
- Fix direction: append one clause — "`R = 200` is sufficient for the cluster-mass empirical `p_mass` to resolve below the locked `0.005` floor (`(R+1)^{-1} = 0.005`); the load-bearing β, γ_l, and δ Grassmann bands sit at this floor."

### Convention changes (4)

**C1 — `T_G*` normalization to [0,1]**
- Already in the pasted methods text (eq. `methods_TGstar`). The locked normalized values for the manuscript Results / Table are:

| Band | Normalized `T_G*` | `p_mass` |
|---|---|---|
| β | **0.273** | 0.005 |
| γ_l | **0.259** | 0.005 |
| δ | **0.149** | 0.005 |
| γ_h | 0.124 | 0.060 |
| θ | 0.0507 | 0.159 |
| α | 0.0312 | 0.348 |

- Action: use these normalized values everywhere `T_G*` is quoted. The raw values (β 69.76 etc.) belong to the audit-trail only and should NOT appear in the manuscript.

**C2 — `T_G^{*,s}` per-patient cluster mass**
- Already defined in the pasted methods (eq. `methods_TGstar_perpatient`). No edit needed in Methods. The per-patient table in Results / Supplement is populated in Phase 3 (Claude side, not writing-agent side).

**C3 — Anatomy labels A1 / A2 (manuscript) — already locked**
- No edit needed. The pasted `sssec:methods_anatomy` already uses A1 + A2; the lab-side A1 / A3 mapping is documented in `locked/ANATOMY_CONTROLS.md`.

**C4 — `τ_min` symbol throughout**
- Already covered in "Locked conventions" above. Find/replace `τ′` → `τ_min` throughout the LaTeX. The user has corrected the previous `τ_max` / `τ′` convention to `τ_min` (smallest diffusion time, finest scale). Methods narrative: "We adopt `τ_min = 1/λ_max`, the smallest diffusion time on the natural resolution window `[τ_min, τ*]`; this is the finest spatial scale resolved by the heat propagator `Π(τ)`."

## What you must NOT do

Anti-patterns from `HANDOFF_INDEX.md` apply.

- Do not re-derive verdicts — they're locked in `locked/VERDICT_LEDGER.md`.
- Do not propose new analyses or new controls.
- Do not introduce KC, VI(k), τ-sweep, K*(b), Pat_03-dropout — all retired.
- Do not change the matched-strength surrogate parameters (`R=200`, swap_target=20, seed=20260511) or the cohort definition (n=10, Pat_02..Pat_15 per VERDICT_LEDGER frontmatter).
- Do not write the Results section based on this directive — this directive is Methods-only. Results revision comes separately after the user shares the β paragraph for Phase 3 cascade.

## Deliverable

Revised Methods LaTeX with all nine edits applied (G2 deferred until the β CSV lands; apply the other eight now). Return as a single LaTeX block + a one-paragraph change log naming the edits applied.

After you return the revised LaTeX, the user pastes it back to Claude Code for `EVALUATION_PROTOCOL.md` verification (number-level CSV cross-check, anti-pattern scan, verdict consistency). Iterate if findings.
