---
name: r2-mesoscale-encoding-rhosym-migration
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-07
supersedes_claims:
  - "R2.2 mesoscale/τ-sweep numbers are ρ_split-era (audit_103c/d)"
  - "R2.3 encoding→OFC localization numbers are ρ_split-era (audit_110)"
pointers:
  - scripts/01_compute/audit/audit_157_arc_scale_null_rhosym.py
  - scripts/01_compute/audit/audit_158_encoding_localization_rhosym.py
  - data/audit/consolidation_arc_rhosym/arc_scale_cohort_verdict.csv
  - data/audit/inference_localization_rhosym/encoding_localization_rhosym_{include,exclude}.csv
---

# R2.2 (mesoscale) + R2.3 (encoding→OFC) — ρ_sym migration

## Head

The two remaining ρ_split-era legs of the flagship — the **mesoscale scale-signature**
(R2.2) and the **encoding→OFC anchor** (R2.3) — are re-estimated under **ρ_sym**
(`audit_157`, `audit_158`), so no ρ_split number lands in R2. **Both verdicts hold
like-for-like.** All correctness anchors PASS bit-exact.

## audit_157 — arc τ-sweep + mesoscale null under ρ_sym (R2.2)

Ports `audit_103c/d` (τ-sweep + mesoscale matched-strength null) to the symmetric
two-arm concordance of `audit_152`, at each τ_mult ∈ {1.0, 2.61, 6.81}. Reuses the
cached seed-20260511 surrogate eigendecompositions (cophenetic reformed at τ; **no
regeneration**). Anchors: arm1 == `audit_103d` ρ_split at every τ (n=180, max|Δ|=1.1e-16);
sym at τ=1 == `audit_152` fine-scale ρ_sym (n=60, max|Δ|=1.1e-16).

**β `T_infspec_pe` (inference-specific | encoding), matched-strength Wilcoxon p:**

| τ/λmax | obs med | surr med | wilcox p | LO-P15 | note |
|---|---|---|---|---|---|
| 1.00 | +0.091 | +0.005 | 0.0098 | 0.020 | fine-scale lock (== audit_152) |
| **2.61** | **+0.126** | +0.009 | **0.0049** | **0.010** | **strongest (mesoscale)** |
| 6.81 | +0.121 | +0.006 | 0.0186 | 0.037 | still clears |

Every **control band** (α, δ, θ, low-γ, high-γ) is **null at every τ** (all p ≥ 0.25;
smallest control p at the mesoscale = high-γ 0.246). Neg-controls staying null at the
mesoscale rules out the coarse-graining artifact (the retired `audit_103b` δ-false-positive
lesson; δ here stays null, p=0.42). Verdict: inference-specific consolidation is
**multiscale, mesoscale-favouring** (strongest at τ≈2.6, not exclusive) — ρ_sym.

## audit_158 — encoding→OFC localization under ρ_sym (R2.3)

Ports `audit_110`'s **encoding** target to `audit_151`'s symmetric per-node concordance
(task-side phase = `task_learn`). R=1000 (beta R=1000 task_learn surrogates cached from
`audit_111`), both epi modes, BH over the same 10-system family as `audit_155`. Anchor:
arm1 encoding M_obs == `audit_110` R=1000 ρ_split (max|Δ|=0.0, all 10 systems).

| condition | OFC M_obs | MS_p | BH q_up | strength_dev | verdict |
|---|---|---|---|---|---|
| epi-include | +3.77e5 | 0.004 | **0.040** | −0.54 | **CARRIER (only survivor)** |
| epi-exclude | +5.57e5 | 0.001 | **0.010** | −0.60 | **CARRIER (only survivor)** |

**PFC DEPLETED** both modes (q_lo=0.010). OFC is the **same hotspot** the test-phase trace
concentrates in (R1.2 / `audit_155` R=1000, q_up=0.010 both modes) → **OFC anchored by both
task phases**. The symmetrization slightly softens OFC (arm1/audit_110 q=0.010 → ρ_sym
q=0.040 include) but it still clears BH. Encoding uses only the **shorter learn phase** →
carries no test/learn length asymmetry → the encoding anchor is duration-immune.

## What R2.2 / R2.3 say (ρ_sym)

- **R2.2:** β inference-specific consolidation clears matched-strength at every scale, is
  strongest at the mesoscale (τ≈2.6, p=0.005), controls null throughout → mesoscale-favouring
  (not exclusive); the mesoscale = scale of multi-step relational integration.
- **R2.3:** the encoding (learning-phase) β trace concentrates in OFC (only carrier;
  q=0.010 excl / 0.040 incl; below-average strength; PFC depleted), the same hotspot as the
  test-phase trace → OFC anchored by both phases; duration-immune (shorter learn phase).

Framework purity: cophenetic only. Both migrations reuse cached surrogates (no regeneration),
verify arm1 == ρ_split bit-exact, and do NOT edit the ρ_split originals (new scripts/dirs).

---

# R2.4–R2.7 — ρ_sym migration (added 2026-07-07)

The rest of the flagship's supporting cast, migrated to ρ_sym. All arm1 cross-checks PASS
bit-exact; no ρ_split number in the paper.

## R2.4 — duration control under ρ_sym (band dissociation + length regression)

Band dissociation is `audit_152` (already ρ_sym). Duration regression recomputed 2026-07-07:
per-patient ρ_sym `T_infspec_pe` (`consolidation_arc_rhosym/arc_per_patient.csv`) vs test/learn
length ratio (`inference_localization/length_control.csv`), Spearman, n=10.

| band | ρ(component, length ratio) | p | note |
|---|---|---|---|
| **β** | **+0.10** | **0.78** | the inference band — **does NOT track length** (was +0.25 ρ_split) |
| α | +0.55 | 0.098 | length-tracking, but inference-null band |
| high-γ | +0.70 | 0.025 | length-tracking, but inference-null band |

⚠️ **Outline correction:** "only α tracks length" is FALSE under ρ_sym — high-γ tracks it harder.
Honest framing: length-tracking sits in **non-tracing** bands (α, high-γ); β (the inference band)
is duration-flat → length and the trace are disjoint by rhythm.

## R2.5 — low-γ cingulate memory (`audit_159`, reruns `audit_112`)

Symmetric within-system ABSOLUTE test. Anchor: arm1 == `audit_112` `within_system_trace`
(max|Δ|=9.7e-17). Low-γ **encoding → cingulate**, the only BH survivor:

| epi | ρ_sym (arm1) | Wilcoxon p | BH q | n_pos |
|---|---|---|---|---|
| include | +0.25 (+0.39) | 0.0039 | **0.035** | 8/8 |
| exclude | +0.39 (+0.51) | 0.0039 | **0.035** | 8/8 |

Encoding-specific (same cingulate edges: inference null p=0.93, standard p=0.07). β-cingulate
encoding is null (p=0.24/0.073, `audit_158`) → cingulate carries **memory at low-γ**. **Solid.**

## R2.6 — inference → cingulate (`audit_160` full-length + `audit_161` length-matched)

`audit_160` (reruns `audit_110`, R=1000): arm1 == `audit_110` inference_pe (max|Δ|=5.8e-11).
**β inference → cingulate CLEARS the full-length null** — sole concentration, q_up=0.030 include /
0.040 exclude, near-average coupling. **Stronger than the ρ_split "BH-borderline hint".**

`audit_161` (reruns `audit_113b` length-matched control, R=200): arm1 == `audit_113b`
(max|Δ|=5.8e-11). Truncating `task_test` to `task_learn` length, **cingulate does NOT survive**:
BH q=0.15 include (raw p=0.030 but cingulate/insula/occipital tie ~q=0.15, none clears) / q=0.30
exclude (occipital leads q=0.10; cingulate p=0.060 q=0.30). → the full-length concentration is
**length-assisted**, replicating the ρ_split downgrade. **Verdict: cingulate = DIRECTIONAL LEAD,
not an established inference location.** R2.6 reports BOTH nulls.

## R2.7 — cingulate multiplexing (synthesis, no new compute)

2×2 (ρ_sym): cingulate = **memory at low-γ** (encoding q=0.035, inference null) and **inference at
β** (inference q=0.030/0.040 full-length, encoding null). One region, content by frequency —
**firm for the memory half, provisional for the inference half** (length-assisted).

## New scripts / outputs (R2.4–R2.7)

- `audit_159_within_system_trace_rhosym.py` → `inference_localization_rhosym/within_system_trace_rhosym_{include,exclude}.csv`
- `audit_160_inference_localization_rhosym.py` → `inference_localization_rhosym/inference_pe_localization_rhosym_{include,exclude}.csv`
- `audit_161_inference_lenmatched_null_rhosym.py` → `inference_localization_rhosym/lenmatched_null_rhosym_R200_{include,exclude}.csv`
- Paragraphs: `results_paragraphs/R2.{4,5,6,7}_*.tex`.
