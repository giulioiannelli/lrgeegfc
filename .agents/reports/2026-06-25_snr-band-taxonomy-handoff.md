---
name: snr-band-taxonomy-handoff
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-06-25
updated: 2026-06-25
pointers:
  - .agents/reports/2026-06-25_per-node-trace-anatomy-and-heterogeneity.md
  - .agents/guides/task-persistence-investigation/2026-06-25_per-node-trace-decomposition.md
  - .agents/preprint/headlines/02_encoding_vs_inference.md
  - .agents/reports/2026-06-18_inference-mark-handoff.md
  - scripts/01_compute/audit/audit_83_localization_matched_strength.py
  - scripts/01_compute/audit/audit_144_per_node_trace_decomposition.py
  - scripts/01_compute/audit/audit_63_split_baseline_surrogate.py
  - data/audit/per_node_trace_decomposition/per_node.csv
  - data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
---

> **Head — resume token: "continue the SNR / band-taxonomy work (run 1→2→3)".**
> The per-node trace investigation produced two pivots the PI has accepted. **(1)
> The per-patient heterogeneity is a measurement-DETECTABILITY axis** — trace
> magnitude is gated by task-vs-baseline SNR (β ρ=+0.78, pooled +0.70, LOO-robust),
> confirmed as reliability not biology (the matched-strength surrogate is itself
> SNR-ordered at ~10× smaller magnitude). **(2) Under the SNR lens the bare trace is
> BROADBAND, not β-specific** — among high-SNR patients it clears matched-strength in
> β/γ_l (5/5), γ_h (4/5), δ/α (3/5); **θ is the only genuinely trace-free band**.
> **DECISIONS:** promote **N2 (encoding/inference dissociation)** to the central core
> — it makes β special for *cognitive content* (inference-specific consolidation),
> not bare-trace exclusivity, which *resolves* the broadband finding instead of being
> threatened by it. Recover the band narration as a **three-tier consistency
> taxonomy**: consistent (β, γ_l) / patient-specific (δ, α, γ_h) / absent (θ).
> **NEXT: run 1→2→3 below, augmented with the PI's two methodological requirements
> (Q1 full multi-phase SNR; Q2 multiscale-must-beat-raw).** Do NOT touch the locked
> ledgers until the taxonomy is locked.

---

## 0. Resume contract

After compaction, say "continue the SNR / band-taxonomy work". Then:
1. Read this file + `2026-06-25_per-node-trace-anatomy-and-heterogeneity.md` (the result)
   + the N2 home `.agents/preprint/headlines/02_encoding_vs_inference.md`.
2. Execute steps **1 → 2 → 3** (§3), folding in **Q1** (§4) and **Q2** (§5).
3. Honor constraints (§7). Python: `/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python`
   (the `conda run` wrapper errors on a g++ hook — call the interpreter directly).

## 1. What is DECIDED / LOCKED (do not relitigate)

- **N2 promoted to central core.** The encoding/inference dissociation (task_learn vs
  task_test; β inference-specific consolidation arc, duration-robust, →OFC,
  mesoscale-verified) is the flagship. It is SAVED and solid (memory
  [[arc_inference_consolidation_2026_06_18]], [[inference_mark_localization_2026_06_18]];
  handoff `2026-06-18_inference-mark-handoff.md`). β is special for **cognitive
  content**, not bare-trace exclusivity.
- **Heterogeneity = SNR-detectability** (not anatomy/implant/epi — all those were
  tested and null; see the result report §4). Apparent non-tracers are SNR-limited,
  not trace-free. This DEFENDS the cohort claim.
- **Bare trace is broadband, θ-exempt.** β-exclusivity was a detectability artifact.
- **Anatomy result (ρ^coph, β):** carried by OFC/MTL/cingulate (OFC locked via
  audit_83 q=0.009); **sensorimotor is the only net-anti system** (within-patient
  p=0.006) — SPECIFICITY not antagonism (closing test 0/6 beat random).
- **Library:** `rank_concordance`, `canonical_cophenet`, `node_incidence_mean`
  promoted to `node_localization.py`; audit_83 delegates (bit-identical, OFC lock
  unchanged).

## 2. The narrative spine (the thing we are building toward)

- **Central core:** N2 encoding/inference dissociation (cognitive result; β carries
  the inference-specific component, every other band null on it).
- **Supporting result A — the offline trace is broadband, θ-exempt, SNR-gated.**
  Methodological frame: detectability explains the per-patient heterogeneity AND why
  the trace looked β-exclusive before.
- **Supporting result B — the recovered band narration: a three-tier CONSISTENCY
  taxonomy** (NOT presence). Bands differ in how *cohort-consistent* the trace is:

  | tier | bands | reading |
  |---|---|---|
  | consistent | **β, γ_l** | present AND same structure across patients (β→OFC cognitive; γ_l→limbic) |
  | patient-specific | δ, α, γ_h | real where SNR allows, localizes differently per patient |
  | absent | θ | no trace even at high SNR — the one true negative |

  This unifies the gate report's presence-vs-consistency axis with the SNR finding:
  the "present-but-split" bands split *because* of SNR; their trace is patient-specific.

## 3. THE PLAN — run 1 → 2 → 3

### Step 1 — per-band localization LOCK (converts the taxonomy from hypothesis to result)
The discriminator for "consistent vs patient-specific" is **cohort-consistent
localization**, tested with the VALIDATED audit_83 machinery (cohort-median per-system
trace + matched-strength + shaft-collapse + LOO), run on **all 6 bands** — NOT the
crude n=5 pairwise-Spearman check used in discussion (that even mis-read β as
inconsistent vs the locked OFC result, so it is unreliable; discard it).
- **Prerequisite:** audit_83 reads the matched-strength surrogate eig cache, which
  exists at R=200/seed=20260511 only for {α, β, γ_l}. **Generate the cache for
  {δ, θ, γ_h}** first via `lrg_eegfc.utils.surrogate.load_or_compute_surrogate_eigs`
  (canonical R=200, SWAP=20, seed=20260511, 4 phases rest_pre_A/B, task_test,
  rest_post). Then run `audit_83 --band <b>` per band.
- **Also (sensorimotor-specificity lock):** add a **lower-tail** matched-strength test
  to audit_83 (M_obs < surrogate, per system) so "trace lives in the cognitive-map
  core, absent in sensorimotor" is matched-strength-grade, comparable to the OFC
  upper-tail. Additive column only (do not change existing columns → OFC lock safe).
- **Output verdict per band:** is there a cohort-consistent system localization
  (→ consistent tier) or not (→ patient-specific tier)? β→OFC already locked; decide
  γ_l, δ, α, γ_h. This locks the taxonomy.

### Step 2 — the detectability figure
Per-band `ρ_split` vs SNR scatter (patients labelled) + the obs-vs-surrogate
separation (high-SNR clears surrogate, low-SNR doesn't), showing broadband + θ-exempt.
This is the figure that defends the cohort claim against "half your patients don't
show it". PDF only, `use_lrg_style()`, no suptitle. **Build AFTER Q1 (§4) so the SNR
axis on the figure is the corrected multi-phase SNR, not the rsPre-only proxy.**

### Step 3 — headline-structure discussion
Once the taxonomy is locked: N2 central; "broadband trace, θ the exception,
consistency-graded by band" as the framing headline. Discuss with PI; then (only then)
cascade into the locked ledgers / preprint. **Cascade is ON HOLD until here.**

## 4. Q1 (PI) — full multi-phase SNR, not rsPreA/B alone

**The concern is valid.** Current `d_noise = 1−Spearman(D_preA, D_preB)` only estimates
**rest_pre** reliability. ρ_split's attenuation is governed by ALL four phases:
`observed ρ ≈ true ρ × √(rel(Δ_task)·rel(Δ_rest))`, where `rel(Δ_task)` needs
`rel(D_tt)` and `rel(D_preA)`, `rel(Δ_rest)` needs `rel(D_post)` and `rel(D_preB)`.
**Fix:** split-half **each** phase — task_test → A/B, rest_post → A/B (rest_pre halves
exist) — get per-phase cophenetic reliability `rel_phase = Spearman(D_phaseA,
D_phaseB)`, build the full attenuation model, and recompute SNR.
- **Prerequisite:** task_test + rest_post half-FCs are **NOT cached** (only rest_pre).
  Compute via `compute_imcoh_abs_halves` (the audit_63 helper; Welch nperseg//2),
  cache alongside `data/cache/imcoh_halves_fc/`.
- **The scientific question it answers (the real one):** after correcting for full
  per-phase reliability, is there RESIDUAL between-patient variability (biology) or is
  the heterogeneity ENTIRELY measurement reliability? If entirely reliability →
  strong "universal trace, detectability-limited" claim. If residual → a real
  biological axis remains to characterize. Either is publishable; we need to know which.
- Re-verify the SNR→ρ headline with the corrected SNR before it goes in the figure.

## 5. Q2 (PI) — the multiscale tool MUST beat the raw per-edge comparison

**Locked rule** ([[feedback_results_only_in_laplacian_framework]]): raw FC is the
baseline, never the result. If a raw per-edge trace tells the same story, the LRG
machinery is redundant. **On the SAME high-SNR patients**, compute the raw-edge analog
`raw_ρ = Spearman(A_tt − A_preA, A_post − A_preB)` over `triu(FC)` and compare to
cophenetic ρ_split on three axes:
1. **Magnitude** — does the hierarchy persist MORE than its edges (cophenetic > raw)?
2. **Story** — is the band-consistency taxonomy (β/γ_l consistent vs δ/α/γ_h
   patient-specific vs θ absent) **multiscale-specific** or already in raw edges? If
   raw reproduces it, it is NOT our result.
3. **SNR gating** — is SNR→trace a raw-FC property or multiscale-specific?
**Win condition** ([[imcoh_lrg_second_order_not_nonlinear]]): the multiscale tool
surfaces higher-order GRAPH structure (community / multi-step / hierarchy persistence,
OFC localization, the inference component) that pairwise edges do not. The raw
comparison is the FOIL that proves the multiscale adds value — put it in the
figure/narrative, not just a private check. Run this on the locked taxonomy bands.

## 6. Key numbers (so the next session does not re-derive)

- **High-SNR patients (median SNR across bands):** Pat_06 (6.1), Pat_05 (2.0),
  Pat_02 (1.6), Pat_03 (1.4), Pat_08 (1.0). **Low-SNR:** Pat_15, Pat_14, Pat_07,
  Pat_13, Pat_10 (≤0.94). SNR is largely a patient-level (recording-quality) property.
- **Per-band, clears own matched-strength surrogate (obs_rho > surr_p95):**

  | band | full cohort | high-SNR-5 | low-SNR-5 | snr↔ρ (n=10) |
  |---|---|---|---|---|
  | δ | 4/10 | 3/5 | 1/5 | +0.87 |
  | θ | 2/10 | 2/5 | 0/5 | +0.37 |
  | α | 5/10 | 3/5 | 2/5 | +0.33 |
  | β | 7/10 | 5/5 | 2/5 | +0.78 |
  | γ_l | 5/10 | 5/5 | 0/5 | +0.72 |
  | γ_h | 4/10 | 4/5 | 0/5 | +0.96 |

  (CAVEAT: high-SNR is n=5, Pat_06 an SNR outlier. Direction clear; exact ranking not.
  "clears surrogate" is the attenuation-robust metric — trust it over raw ρ.)
- **β per-system trace propensity (gray, patient-demeaned):** OFC +0.37, MTL +0.28,
  cingulate +0.16 (carriers); sensorimotor −0.28 (only net-anti). γ_l: MTL-led,
  occipital most-anti. α: flat/incoherent.
- **SNR control:** snr↔observed-ρ +0.78, snr↔surrogate-ρ +0.84 but surrogate ρ ∈
  [−0.004, +0.091] (~10× smaller than observed [−0.09, +0.51]) ⇒ detectability, not
  artifact. Encoding partial (d_task | d_noise, full clean baseline) +0.61.

## 7. Constraints (unchanged + two new methodological musts)

1. Patient is not an axis unless explained — now explained (SNR-detectability), so OK.
2. **BAND result stays** — recovered as the consistency taxonomy (§2 tier table).
3. **ρ^coph ⊥ Grassmann distinct** — this whole thread is ρ^coph; Grassmann per-node
   (`grassmann_trace_contributions`) untouched, stays the distinct subspace story.
4. Machinery-dependence is INTERNAL, never in the paper.
5. **Prioritize OLD results** — N2 (now central) + the β trace + OFC.
6. Decimation control mandatory for any "exclude subset → trace changes" claim.
7. Brutal honesty, no confidence laundering.
8. **NEW (Q1):** SNR must use full per-phase reliability, not rsPreA/B alone.
9. **NEW (Q2):** every multiscale trace claim must be shown to BEAT the raw per-edge
   comparison on the same patients (raw is baseline, never the result).
10. **Run SNR/reliability gating EARLY** in any future per-patient analysis (the lesson
    of this session — it should have been the first check).

## 8. Artifacts produced this session

- Scope: `task-persistence-investigation/2026-06-25_per-node-trace-decomposition.md`.
- Code: `audit_144_per_node_trace_decomposition.py`; library helpers in
  `node_localization.py` (audit_83 delegates, bit-identical).
- Data: `data/audit/per_node_trace_decomposition/{per_node,per_patient_summary,
  characterization_*}.csv`.
- Report: `2026-06-25_per-node-trace-anatomy-and-heterogeneity.md` (anatomy + SNR
  explanator §4b). Memory: [[per_node_trace_anatomy_2026_06_25]].
- Gate report `2026-06-25_cophenetic-gate-presence-vs-consistency.md` §6 superseded.
- **No locked-ledger edits. Cascade ON HOLD until the taxonomy is locked (Step 1).**
