---
name: methods-section-review-2026-05-19
era: IMCOH_ABS_COHORT_N10
status: mostly_superseded
status_updated: 2026-05-28
superseded_by: METHODS_AUDIT_ISSUES.md
kind: review-checklist
scope: factual + consistency review of the draft Methods LaTeX (Welch → ImCoh → LRG → comparison probes → cluster-extent → stats + robustness) against the locked artifacts and the actual code
companion: methods_revision_2026-05-18_cophenet.md, methods_grassmann_cluster_extent.md, CONTROLS.md, VERDICT_LEDGER.md, METHODS_AUDIT_ISSUES.md
---

> **Note (2026-05-28).** This file is **mostly superseded** by
> [`METHODS_AUDIT_ISSUES.md`](../METHODS_AUDIT_ISSUES.md), which carries
> the live methods-section audit items + C1 normalization lock + Decision-12
> cascade. M1, M4, m6, m9 here are marked RESOLVED/SUPERSEDED inline; the
> §54-55 "pre-fix vs post-fix" passage has a forward-pointer footnote to
> the C1 normalized values. Kept verbatim as historical record of the
> 2026-05-19 review pass — do **not** cite live numbers from this file;
> route to `METHODS_AUDIT_ISSUES.md` instead.

# Methods section review — what shall be corrected (2026-05-19)

**Head.** The draft Methods is structurally sound and the LRG-pipeline framing matches the locked artifacts; the changes below are factual corrections, not framing edits. Four are MAJOR (the `T_G^*` formula is stale post the 2026-05-19 pm fix; the communication distance is written with `K` instead of `ρ`; the matched-strength surrogate description overstates what is preserved; the epi-X robustness criterion diverges from the locked C5 reading). Five are MODERATE (under-specified numbers / thresholds: `T_G < 0` trace direction, C4 cross-probe gate, `R = 200`, cohort `n = 10`, BH family sizes). The rest are minor polish. After these fixes, every numerical or symbolic statement in Methods will be traceable to either a locked artifact or a script line.

**UPDATE 2026-05-19 pm**: writing-agent feedback prompted three further
gate refactors that supersede the earlier M4 / m6 wording in this
review: (i) Grassmann verdict gate becomes **mass-only**
(`cluster_p_mass < 0.05`); the disjunctive `min(p_LR, p_mass)` rule is
retired (CONTROLS.md §C3 Decision 8). (ii) C4 cross-probe gate becomes
a **paired Wilcoxon non-degradation test**; the `n_trace_xprobe ≥ 6/10`
threshold is retired (CONTROLS.md §C4 Decision 9). (iii) C5 epi-X gate
becomes a **Wilcoxon on epi-X data** (Grassmann: cluster-mass on epi-X
eigvec cache; cophenet α: one-sample Wilcoxon on per-patient
`obs_rho^epi-X`); the ≳80% retention rule is retired (CONTROLS.md §C5
Decision 10). All three gates now consistently avoid hardcoded
patient-count or magnitude thresholds, per
`feedback_no_hardcoded_test_thresholds.md` and writing-agent feedback
2026-05-19. Every cohort Wilcoxon also gets a LOO max-p diagnostic
column per `feedback_no_single_patient_p_driven.md` — descriptive,
never a gate. Audit scripts: `audit_71_c4_wilcoxon_cohort.py`,
`audit_72_c5_wilcoxon_cohort.py`; outputs in
`data/audit/ctm_triangle/c4_wilcoxon_cohort.csv`,
`data/audit/{grassmann_epi_exclusion,alpha_epi_exclusion}/c5_wilcoxon_cohort.csv`.

---

## MAJOR — must fix before submission

### M1 — `T_G^*` formula in eq. `(methods_TGstar)` is OBSOLETE

**Current**: `T_G^*(b) = Σ_{k ∈ C^*(b)} (−log_10 p_k(b))`, with `C^*` defined
as "the longest contiguous run of `k`-cells satisfying `p_k(b) < α_k`".

**Correct**: `T_G^*(b) = Σ_{k : p_k(b) < α_k} (−log_10 p_k(b))` — the
sum runs over **all** significant `k`-cells across **every**
contiguous-significant cluster.

**Why**: `audit_70.cluster_mass` was corrected this session
(2026-05-19 pm) from the longest-run-only sum to the resilient
all-clusters sum. The all-clusters formula is robust to a single
non-significant `k` that fragments a long cluster (a length-20 cluster
split into 10+10 keeps total mass; the longest-run sum halves it).
After the fix, γ_l and δ both upgrade weak → strong on the disjunctive
gate (`cluster_p_mass = 0.005` for both). The pre-fix locked values
in the manuscript draft (β 52.97, γ_l 19.17, δ 12.78) are obsolete;
re-run values are β 69.76, γ_l 66.14, δ 38.07.

> **Forward pointer (added 2026-05-28).** Per C1 lock (`METHODS_AUDIT_ISSUES.md:116-129`)
> these raw values are subsequently normalized to `[0,1]`: β → **0.273**,
> γ_l → **0.259**, δ → **0.149**. Per Decision-12 cascade
> (`VERDICT_LEDGER.md:519+`, locked 2026-05-28), δ is **demoted back to
> "weak"** because full-data LOO max `p_mass = 0.055 (Pat_08)` fails
> the < 0.05 LOO precondition; the "weak → strong" upgrade above held
> under Decision 8 mechanical rule but is retracted under Decision 12.
> γ_l remains strong (LOO 0.040 Pat_05 passes). The raw values quoted
> above are kept verbatim as historical record of the pre-fix → post-fix
> formula change; live citations should use normalized + Decision-12
> verdicts.

**No `min_cluster_size` threshold is added** despite the all-clusters
sum including isolated significant `k`-cells. The empirical null
calibration absorbs the singleton-driven concern: under H₀ with 111
k-cells at `α_k = 0.05` the expected null configuration is ≈ 5.55
scattered significant cells, so the null distribution of `cluster_mass`
already sees the configurations the concern points at. Quantitatively
(`data/audit/grassmann_cluster_extent/cohort_summary.csv`): null p95
mass is 18–32 across bands; 3 isolated cells at the threshold p give
mass ~4 (far below null mean); 3 cells at the deepest possible
Wilcoxon p (≈ 0.001 at n=10) give mass ~9 (still below null mean).
The empirical `p_mass` for any such configuration returns ≈ 1. The
test correctly refuses to call 3 isolated cells a trace — through
the null, not through a hardcoded statistic-level threshold. This is
compliant with `feedback_no_hardcoded_test_thresholds.md`: the test
is the gate.

**Verdict-anatomy diagnostic for our data** (decomposition of `T_G^*`
into per-cluster contributions): β singleton mass 4%, γ_l 4%, δ 15%
of total. Load-bearing mass comes from multi-cell contiguous
sub-clusters, not isolated singletons — verdicts are structurally
defensible *because of the data*, not by hardcoded construction.

**Where**: eq. `(methods_TGstar)` and the preceding sentence
introducing `C^*`. Keep `L_obs = |C^*|` as the contiguity-based
co-primary statistic — that one is still the longest run.

**Reference**: `.agents/preprint/directives/methods_grassmann_cluster_extent.md`
§5b (formula), §5c (null formula), §5e (no-hardcoded-min-size
rationale + null-calibration math + verdict-anatomy diagnostic), §6
(locked table with re-run values); memory
`feedback_cluster_mass_all_clusters.md`,
`grassmann_null_calibration_protects.md`.

---

### M2 — Eq. `(methods_distm)`: `D(τ')` is `1/ρ(τ')`, not `1/K(τ')`

**Current**: `D_ij(τ') = (1 − δ_ij) / K_ij(τ')`.

**Correct**: `D_ij(τ') = (1 − δ_ij) / ρ_ij(τ')`.

**Why**: the code at `lrgsglib/src/lrgsglib/utils/lrg/spectral.py:121-124`
computes:
```python
num = expm(-tau * L)              # K(τ)
den = np.trace(num)               # Z(τ)
rho = num / den                   # ρ(τ) = K(τ)/Z(τ)
Trho = 1.0 / rho                  # D(τ) = 1/ρ(τ) = Z(τ)/K(τ)
```
So `D = 1/ρ`, not `1/K`. The two differ by the scalar `Z(τ) = Tr[e^{−τL}]`.
The cophenetic dendrogram is invariant under this multiplicative
rescaling (UPGMA is monotone-rescale invariant), so all reported
numbers are unchanged. But the statement as written is inconsistent
with both the code and Villegas 2025 Eq. 1 (which uses `ρ`).

**Where**: eq. `(methods_distm)` and the sentence introducing it
(*"the propagator induces a per-pair communication distance..."*).

**Reference**: `.agents/guides/02_methods/lrg-framework-guide.md`
§2.5 (the framework guide has the same error in its formula listing
and should also be cleaned up in the same pass).

---

### M3 — Matched-strength surrogate: "preserves global edge-weight distribution" is too strong

**Current**: *"...repeated 4-cycle ±δ edge-weight rewirings that
preserve, by construction, the per-node strength sequence **and the
global edge-weight distribution** of the original |ImCoh| graph."*

**Correct (suggested wording)**: *"...repeated 4-cycle ±δ edge-weight
rewirings that preserve, by construction, the per-node strength
sequence (verified to 10⁻⁶ per surrogate) and keep each edge weight
within the original [0, 1] interval; the marginal edge-weight
distribution drifts within these constraints."*

**Why**: the 4-cycle ±δ rewiring at
`scripts/01_compute/audit/audit_67_raw_fc_matched_strength.py:155-186`
preserves only:
1. Per-node strengths to 10⁻⁶ (verified by `verify_strengths` lines 189-193);
2. Global edge-weight SUM (= 2 · Σ_i s_i, consequence of (1));
3. Edge values in `[0, w_max]` (via the `lo`/`hi` bounds at lines 167-176).

It does NOT preserve the histogram of edge values — each ±δ is drawn
fresh from `Uniform(lo, hi)` per swap, so the distribution drifts
under repeated application.

**Where**: §`sssec:methods_compare_stats`, the matched-strength
surrogate description.

---

### M4 — Epi-X robustness criterion (RESOLVED 2026-05-19 pm)

**Resolved via CONTROLS.md §C5 Decision 10**: the ≳80% retention
rule was retired (arbitrary hardcoded threshold). C5 now uses a
Wilcoxon-based gate:

- **Grassmann (all bands)**: re-run audit_70 cluster-mass test on
  the epi-X surrogate eigvec cache. C5 passes iff
  `cluster_p_mass^epi-X < 0.05`. Audit at
  `scripts/01_compute/audit/audit_72_c5_wilcoxon_cohort.py`,
  CSV at `data/audit/grassmann_epi_exclusion/c5_wilcoxon_cohort.csv`.
- **Cophenet (α only)**: one-sample Wilcoxon on per-patient
  `obs_rho^epi-X` under `H_1: rho_split^epi-X > 0`. C5 passes iff
  `wilcoxon_one_sided_p < 0.05`. CSV at
  `data/audit/alpha_epi_exclusion/c5_wilcoxon_cohort.csv`.

LaTeX action: replace the M4 paragraph in §`sssec:methods_compare_stats`
with this Wilcoxon-based C5 description. Cite both audit_72 CSVs and
the LOO max-p diagnostic per
`feedback_no_single_patient_p_driven.md`.

**Reference**: `.agents/preprint/locked/CONTROLS.md` §C5 (updated 2026-05-19
Decision 10).

---

## MODERATE — fill in missing numbers / thresholds

### m5 — Eq. `(methods_T_G)` trace direction is not explicit

`T_G(k; p, b)` is defined but never explicitly tied to `T_G < 0` as
the trace direction. §`ssec:methods_compare` already establishes
`T_d < 0` ⇔ trace via eq. `(methods_triangle)`, so the connection is
implicit — but a single sentence after eq. `(methods_T_G)` would
prevent reader confusion: *"\(T_G(k; p, b) < 0\) is the trace
direction at scale `k`, consistent with the phase-triangle convention
of eq. (methods_triangle)."*

### m6 — C4 cross-probe gate (RESOLVED 2026-05-19 pm)

**Resolved via CONTROLS.md §C4 Decision 9** (option (b) selected
per writing-agent feedback): the `n_trace_xprobe ≥ 6/10` patient-count
threshold was retired. C4 now uses a paired one-sided Wilcoxon on
`(rho_split − rho_xprobe)` per patient under `H_1: rho_split >
rho_xprobe`, **failing to reject** (`paired_wilcoxon_p ≥ 0.05`)
indicates no significant degradation under cross-probe restriction;
this is combined with `sign(rho_xprobe_median) == sign(rho_split_median)`.

Audit at `scripts/01_compute/audit/audit_71_c4_wilcoxon_cohort.py`;
CSV at `data/audit/ctm_triangle/c4_wilcoxon_cohort.csv`. All 6 bands
pass the new C4 gate. LOO max-p diagnostic included per
`feedback_no_single_patient_p_driven.md`.

LaTeX action: replace the m6 paragraph in §`sssec:methods_compare_ctm`
with the Wilcoxon-based C4 description. Cite the audit_71 CSV.

### m7 — `R` (number of surrogates) is referenced but not stated

`R` appears in the LaTeX without a value. Locked: **R = 200**, swap
target = **20 per node** (so `n_swaps_per_surrogate = 20 · N(N−1)/2`),
seed = **20260511**, shared across audits 63 / 66 / 67 / 70. Cite
once in §`sssec:methods_compare_stats`.

### m8 — Cohort size `n = 10` is never stated

The 10-patient cohort (Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15)
is never named in the Methods section. State it once in
§`ssec:methods_fc` (*"for each of the \(n = 10\) patients in the
cohort"*) or §`sssec:methods_compare_stats`. Without it the reader
cannot interpret the Wilcoxon power claims.

### m9 — BH family sizes worth a sentence — **SUPERSEDED 2026-05-20**

**Status:** RETIRED by the 2026-05-20 BH-drop directive
(`METHODS_AUDIT_ISSUES.md` §B2 revised; writing-agent feedback
2026-05-20). The cross-band BH-FDR on `ρ_split^coph` (`m = 6`) is
dropped entirely from methods; `d_S on D_coph` is dropped from
methods entirely (undefined symbol). The only BH-FDR that remains
is the anatomy A1 hypergeometric across DK regions per (band, probe).
See `feedback_no_unmotivated_bh_fdr.md` for the three-point check
that governs this and any future correction. Original m9 text below
preserved for audit trail.

~~The LaTeX says BH family `m = 6` for `ρ_split^coph` (six bands) and
`m = 6` for `d_S` on `D_coph` (six bands). These are two
independent families, both at q = 0.05. Worth one sentence stating
they are independent (`d_S` is a single-scale matrix-level test;
`ρ_split^coph` is the controlled per-pair split-baseline test). The
Grassmann probe explicitly does NOT get a separate BH correction
because the cluster-extent permutation is the family-level gate —
already correctly stated in the LaTeX. Just confirm the two
cophenet-family sizes against `locked/VERDICT_LEDGER.md`.~~

---

## MINOR — internal consistency / polish

### n10 — k-grid upper bound

LaTeX *"the upper bound matching \(N - 1\) of the smallest-N patient"*
is verified by
`scripts/01_compute/audit/audit_66_grassmann_matched_strength_surrogate.py:107`:
*"min(N) across cohort = 113 (Pat_10) → max k = 112"*. ✓ Statement
accurate. Optional: name Pat_10 explicitly, but the abstract phrasing
is acceptable.

### n11 — `D̄` notation and its relation to `D_coph`

The LaTeX symbol for the per-pair communication distance (renders
something like `D̄_ij`) is introduced, immediately fed to UPGMA, and
then superseded by `D_coph`. The §`ssec:methods_lrg` paragraph
already says *"The raw distance \(\distm(\tau')\) is not read
directly"* ✓. No edit needed — just confirm no later equation
references `D̄` outside the linkage construction.

### n12 — `rspost` halves required by `ρ_drift`

Eq. `(methods_rho_drift_coph)` uses
`D_coph^rspost_B − D_coph^rspost_A`. This means **rspost is split
into halves** (not just rspre as in the substrate audit_67).
**Verify** that a script generates `rspost_A`, `rspost_B` half-FCs
and runs the full LRG pipeline on them; the substrate-level
`audit_67_raw_fc_matched_strength.py:122-132` only splits rspre. If
no rspost-half cache exists, the `ρ_drift` equation in the LaTeX
cannot currently be evaluated and either (a) the C2 control needs
to point at a different audit that does compute rspost halves, or
(b) `ρ_drift` is computed in a different way than what eq.
`(methods_rho_drift_coph)` states. **1-minute grep needed**.

### n13 — Halved-data noise budget for `ρ_drift`

The LaTeX says `ρ_drift` *"captures the within-session drift
component of the per-pair multiscale correlation"*. Consider one
sentence about the **noise-budget equivalence** with `ρ_split`
(`ρ_drift` is built from half-data shifts on both sides, matching
the half-data noise floor of `Δ_task = D_coph^taskt − D_coph^rspre_A`
and `Δ_rest = D_coph^rspost − D_coph^rspre_B`). Pre-empts the
reviewer's concern that `ρ_drift` might be noisier because it sees
only rest data.

### n14 — Phase-set glossary entries

The LaTeX uses `\txtacr{rspre}`, `\txtacr{taskl}`, `\txtacr{taskt}`,
`\txtacr{rspost}`. These correspond to codebase
`{rest_pre, task_learn, task_test, rest_post}`. Cross-check that the
`\txtacr{}` macro is applied **uniformly**: I noticed `rspost`
appearing both with and without macro wrapping in different
equations. Single-pass cleanup.

---

## What is already correct (no edits needed)

1. **Bands table** matches `BRAIN_BANDS` in `src/lrg_eegfc/config/const.py:191`
   exactly, including γ_h = 80-300 Hz (memory `brain_bands_definition.md`).
2. **Welch nperseg** `= fs × 2` gives `Δf = 0.5 Hz` at both 2048 Hz
   and 1024 Hz — matches `nperseg_for_fs(fs)` at `const.py:235-241`.
3. **`<|ImCoh|>_f`** eq. `(methods_imcoh_bandavg)` is the canonical
   `imcoh_abs` operator (memory `imcoh_taxonomy.md`).
4. **Combinatorial Laplacian `L = D − A`** matches
   `audit_66_grassmann_matched_strength_surrogate.py:176-180` and
   the canonical Villegas fluid Laplacian
   (`.agents/guides/02_methods/lrg-framework-guide.md:53`). The
   explicit ruling-out of `L_RW` is correct and addresses the
   writing-agent's earlier flag.
5. **`τ' = 1/λ_max`** as the canonical scale matches
   `spectral.py:120` (`tau = 1/max(spectrum)`).
6. **τ-sweep retired** rationale (continuous spectrum, no entropy
   plateaus, dendrogram replaces gap-based scale ID) matches
   `methods_revision_2026-05-18_cophenet.md` and memory
   `lrg_outlier_case_fully_connected.md`.
7. **Cophenetic image** as multiscale carrier — correct framing.
8. **Phase-triangle scalar** `T_d = d(taskt, rspost) − d(rspre, taskt)`
   with `T_d < 0` as trace direction matches
   `result_1_raw_fc_phase_trace.md` memory.
9. **Substrate `ρ^raw_split`** construction matches
   `audit_67_raw_fc_matched_strength.py:198-204` exactly.
10. **Matrix-level `d_S`/`d_P`/`d_F`** at substrate as descriptive
    sensitivity (not load-bearing) matches `feedback_dP_framing.md`
    on never framing `d_P` as orthogonal to `d_S`, and
    `result_1_raw_fc_phase_trace.md` on `d_F` as drift diagnostic.
11. **Chordal Grassmann identity** `d_G² = k − ‖A^T B‖_F²
    = Σ sin²θ_i` matches `audit_70.chordal_distances_for_k_grid:65-89`
    (SVD-free via cumsum).
12. **k ∈ {2, ..., 112}** matches `audit_70.K_GRID` and the
    audit_66 `min(N) = 113` rationale.
13. **`T_G(k)` eq. `(methods_T_G)`** uses `rest_pre_A` (early half)
    for the upstream baseline — matches
    `audit_70.compute_surr_T_G_band:109-116`.
14. **Phantom-permutation null** construction matches
    `audit_70:217-228`.
15. **Disjunctive gate `min(p_LR, p_mass) < 0.05`** matches
    `locked/CONTROLS.md` §C3 Decision 7.
16. **One-sided paired Wilcoxon under `H_1: T_d < 0`** matches every
    cohort-test script in the repo.

---

## Suggested fix order (for the writing agent)

1. **M1** — replace the `T_G^*` formula and re-quote the cohort
   numbers (β 69.76, γ_l 66.14, δ 38.07 instead of 52.97 / 19.17 / 12.78).
2. **M2** — swap `K_ij` → `ρ_ij` in eq. `(methods_distm)`.
3. **M3** — soften the matched-strength surrogate description.
4. **M4** — align epi-X criterion with CONTROLS.md C5 (or revise
   CONTROLS.md).
5. **m5-m9** — fill in `T_G < 0` direction, C4 threshold, `R = 200`,
   `n = 10`, BH family clarification.
6. **n10-n14** — minor polish at the writing agent's discretion.

After 1-5, the Methods section should pass `EVALUATION_PROTOCOL.md`
without further changes at the methodology layer.

---

## Source-of-truth references

- `.agents/preprint/directives/methods_grassmann_cluster_extent.md` — locked
  Grassmann methodology + post-fix `T_G^*` definition + locked
  band-level values.
- `.agents/preprint/locked/CONTROLS.md` — locked 5-control battery (C1-C5).
- `.agents/preprint/locked/VERDICT_LEDGER.md` — locked per-band verdicts.
- `.agents/preprint/directives/archive/2026-05/methods_revision_2026-05-18_cophenet.md` — binding
  methods directive (KC retired, VI(k) retired, τ-sweep retired,
  `D_coph` adopted).
- `scripts/01_compute/audit/audit_67_raw_fc_matched_strength.py` —
  substrate per-pair `ρ_split^raw` (M3 source).
- `scripts/01_compute/audit/audit_70_grassmann_cluster_extent.py` —
  cluster-extent permutation, post-fix resilient `cluster_mass` (M1 source).
- `lrgsglib/src/lrgsglib/utils/lrg/spectral.py:121-128` — code
  definition of `D(τ) = 1/ρ(τ)` (M2 source).
- Memory `feedback_cluster_mass_all_clusters.md` — methodological
  justification for M1.
- Memory `methods_grassmann_TG_star.md` — locked `T_G^*` definition
  and re-run values.
- Memory `laplacian_choice_grassmann.md` — combinatorial vs RW vs
  symmetric-normalised Laplacian disambiguation.
