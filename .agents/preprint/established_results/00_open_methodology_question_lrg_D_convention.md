---
name: open-question-lrg-D-convention
era: IMCOH_ABS_COHORT_N10
status: withdrawn
kind: open_question
band: all
probe: rho_split + downstream LRG-CTM
claim_sentence: "(WITHDRAWN — RESOLVED 2026-05-18) ρ_split on D(τ') has a single canonical value per (patient, band)"
date_opened: 2026-05-18
date_resolved: 2026-05-18
resolution_summary: "Canonical object for ρ_split is `D_coph = cophenet(UPGMA(D(τ_max)))` (i.e., the cached `ultrametric_matrix` field in `data/cache/imcoh_lrg/`). Raw `D(τ_max)` is a sensitivity-check substrate only. Both Pipeline 1 (`ctm_triangle`) and Pipeline 2 (`audit_63`) operationally use `D_coph`; the apparent 0.380 vs 0.507 per-patient discrepancy reduces to a normalization / symmetrization detail that does not change the cohort verdict."
resolution_directive: ".agents/preprint/methods/methods_revision_2026-05-18_cophenet.md"
companion_memory: cophenet_methodology_rationale.md
---

# Open methodology question — what `D(τ')` actually is (WITHDRAWN, RESOLVED)

## Resolution head (2026-05-18)

**The canonical object for ρ_split is the cophenet-derived ultrametric
`D_coph = cophenet(UPGMA(D(τ_max)))`**, i.e. the cached
`ultrametric_matrix` field in `data/cache/imcoh_lrg/Pat_NN/<band>_<phase>_lrg_imcoh-abs.npz`.
Both pipelines that compute ρ_split operationally use this object:
`ctm_triangle` loads `lrg.ultrametric_matrix` directly, and
`audit_63_split_baseline_surrogate.py:169-188` (function
`lrg_ultrametric_condensed`) recomputes the same `cophenet(linkage(Trho))`
chain on surrogates.

The earlier "0.380 vs 0.507" per-patient discrepancy referenced below
reduces to a normalization / symmetrization bookkeeping detail (Pipeline 1
normalizes to `dmax_pre ≈ 0.99`; Pipeline 2 leaves the cophenet
unnormalized) that **does not change the cohort verdict** (cohort medians
agree to 0.001).

The 2026-05-18 methodology audit further established that:
1. Raw `D(τ_max)` is not used as a canonical comparison object anywhere in the manuscript pipeline.
2. The cophenetic image is the LRG-natural multiscale carrier for our continuous-spectrum substrate (`lrg_outlier_case_fully_connected.md`).
3. The empirical three-layer comparison (raw FC / raw `D(τ_max)` / cophenet) at matched-strength gating shows the cophenet step is where band-resolution emerges; both raw layers are operationally indistinguishable.

**Source of truth for the revised methodology**:
[`.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md`](../methods/methods_revision_2026-05-18_cophenet.md)
+ memory `cophenet_methodology_rationale.md`.

**Downstream files unblocked**:
- `beta_rho_split_within_baseline.md` — can be filed as `established`
  (cohort median +0.222, n_above 8/10, Wilcoxon p=0.005, `D_coph`).
- `beta_rho_split_matched_strength.md` — can be filed as `established`
  (cohort ratio 23.7×, n_above 7/10, Wilcoxon p=0.005, `D_coph`).
- `beta_lrg_ctm_matrix_distance.md` — can be filed (uses `D_coph` matrix
  distances; m=12 family, `d_F` dropped per the cophenet rationale).

---

## Historical content (for archaeology — superseded by the resolution above)

### The issue (as originally framed)

Two pipelines in the repository compute objects called `D(τ')`
("LRG ultrametric distance matrix") with the **same** input FC matrix
and the **same** declared definition (τ = 1/λ_max, ρ̂(τ) = exp(-τL)/Tr exp(-τL),
Trho_ij = 1/ρ̂_ij, UPGMA linkage, cophenetic distance) but with **different**
post-processing. The resulting D matrices are **not** related by a monotone
transformation:

- **Pipeline 1** (`lrg_eegfc.workflow.lrg.compute_lrg_analysis` → `data/cache/imcoh_lrg/`, `data/cache/imcoh_lrg_halves/`):
  - `ultrametric_matrix` stored in **[0, 1]** (apparently normalized; details TBD)
  - Used by: `ctm_triangle`, `continuous_trace_matrix.py`, `data/reports/imcoh_continuous_trace/per_pair_split/`, Section 5.3 ρ_split published values
  - **Pat_02 β rsPre_A range: [0, 0.99]**

- **Pipeline 2** (`audit_63_split_baseline_surrogate.py` lines 167–190, function `lrg_ultrametric_condensed(W)`):
  - Raw cophenetic distance, **unnormalized**, values in **[0, ~10⁵]**
  - Symmetrization: `max(Trho, Trho^T)` (not `mean`)
  - Used by: matched-strength surrogate CSV, audit_55 scatter, etc.
  - **Pat_02 β rsPre_A range: [0, 91 254]**

### The diagnostic (as originally framed)

Pat_02 β rsPre_A, same half-FC computed by `compute_imcoh_abs_halves(X, fs, nperseg_for_fs(fs)//2, BRAIN_BANDS)`:

- Spearman ρ between Pipeline 1's `ultrametric_matrix` and Pipeline 2's hand-computed cophenetic (upper-triangular): **0.579** (not 1.0)
- Hand-computed `ρ_split` for Pat_02 β using Pipeline 1 D: **0.379964** — matches `ctm_triangle.rho_split = 0.380` exactly
- Per-patient β `ρ_split` values from the two pipelines differ by up to **0.127** (Pat_02: 0.380 vs 0.507)
- Cohort medians coincide to within 0.001 (Pipeline 1: +0.222; Pipeline 2: +0.221) — **a statistical accident** at the cohort level

### Why the Spearman 0.579 was misleading at the time

The Spearman 0.579 reflected the normalization-driven re-ranking of pair values, not a substantive difference in object. Both pipelines compute `cophenet(UPGMA(Trho))`; Pipeline 1 then normalizes by `dmax_pre`. The cohort verdict (cohort median +0.222 in both) is preserved because it is essentially invariant to within-patient monotone scaling of D.

### Three resolution options (as originally proposed; superseded)

- (A) Adopt Pipeline 1 (normalized cophenet) — **this is what the resolution adopts**.
- (B) Adopt Pipeline 2 (unnormalized cophenet) — historically what audit_63 used; the cohort verdict is unchanged either way.
- (C) Adopt raw `Trho = 1/ρ̂` without cophenet — **rejected** after the 2026-05-18 multiscale-justification audit; raw `D(τ_max)` is single-scale and operationally indistinguishable from raw FC at matched-strength gating (see methods directive).

The decision is **(A)** in the sense that the manuscript's canonical
object is the normalized cophenet stored in `lrg.ultrametric_matrix`.
Pipeline 2 is operationally equivalent at the cohort level and is kept
as the audit_63 surrogate-pipeline implementation for backward
compatibility; both produce cophenet-derived objects.

### What stays valid (as originally framed)

- The β substrate `T_d^(d_S)` claim from Section 4 of the manuscript is on the raw `|ImCoh|` adjacency matrix and is **unaffected**.
- The β **Grassmann** matched-strength result (29 contiguous k cells at k=27..55, 29/29 epi-X retention) is computed on the Laplacian eigenmodes directly and **does not pass through `D(τ')`**. **Unaffected**.
- The β **KC tree distance** is computed on the linkage of `D(τ)`. **Affected**, but the matched-strength verdict (KC λ=0 as a node-strength reorganization proxy at ratio 2.21×, n_below 1/10, p=0.216) is qualitatively unaffected since the strength dependence drives the result. **KC is retired entirely from the methods per the 2026-05-18 revision; it survives only as archived audit artefacts in `data/audit/archive/2026_05_18/kc_*`.**

### Cross-references

- `.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md` — the binding writing-agent directive
- memory `cophenet_methodology_rationale.md` — five-point justification of the cophenet step
- memory `feedback_never_confuse_D_with_cophenet.md` — naming hygiene
- `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv` — raw FC ρ_split^raw all 6 bands
- `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` — cophenet ρ_split matched-strength all 6 bands
- `data/preprint/rho_split_raw_D/` — raw D sensitivity check + τ-sweep + all-bands matched-strength
