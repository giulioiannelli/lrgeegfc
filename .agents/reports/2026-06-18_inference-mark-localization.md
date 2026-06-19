---
name: inference-mark-localization
type: report
era: IMCOH_ABS / COHORT_N10
status: result
created: 2026-06-18
updated: 2026-06-18
pointers:
  - scripts/01_compute/audit/audit_110_inference_mark_localization.py
  - scripts/01_compute/audit/audit_110b_inference_localization_figure.py
  - scripts/01_compute/audit/audit_111_gen_task_learn_r1000_surrogates.py
  - scripts/01_compute/audit/audit_112_within_system_trace.py
  - .agents/guides/task-persistence-investigation/2026-06-18_inference-mark-localization.md
  - .agents/reports/2026-06-18_consolidation-arc-handoff.md
---

# Where the β inference-mark sits — it dissociates from OFC into the cingulate

**Head.** Within the β band the offline trace **splits by what it carries**: the
standard trace and the **encoding** component (the memorised pairs) over-accumulate
in **orbitofrontal cortex** (OFC), but the **inference-specific** component — the
figured-out relations, with encoding statistically removed — over-accumulates in
the **CINGULATE**. This is an *accumulation hotspot on a distributed trace*
(per-patient demeaned, [[feedback-localization-is-overexpression-not-container]]),
not a container. **At R=1000 the dissociation is BH-significant on both sides:**
encoding/standard → OFC (q=0.010), inference → cingulate (**q=0.050 epi-include**,
borderline q=0.060 epi-exclude). The cingulate inference hotspot is robust across
all four conditions (contact/shaft × epi-in/out), strength-independent (dev +0.04),
K=8/10, 5/8 patients, leave-one-out-robust median, and survives a length-equalized
duration control (audit_113, below). (Earlier R=200 framing of
"suggestive, q≈0.10–0.20" was a floor artifact — the true tail sat below the
1/201 R=200 resolution; R=1000 resolved it to p=0.005.) The corroborating evidence
is a **double dissociation** — the cingulate carries *inference* in β but
*encoding* in α/low-γ — which independently re-derives the consolidation arc's band
dissociation ([[arc-inference-consolidation-2026-06-18]]) from anatomy alone.

This localizes the arc result: β does not just consolidate "the inference" as an
abstract component — that component has an anatomical home, the cingulate, distinct
from where the memorised material consolidates (OFC).

**Band picture (a focal trace can hide under a null whole-brain trace).** Cohort
whole-brain trace clears only in α (encoding+standard) and β (all three); low-γ is
null brain-wide. Yet the **absolute within-system trace** (audit_112, no demeaning,
Wilcoxon obs>surrogate on system-incident pairs) finds a **genuine focal low-γ
encoding trace in the cingulate** (ρ=+0.40, 8/8 patients, q=0.035; epi-exclude
ρ=+0.51, 7/8, q=0.070 — bigger effect, power-limited) — a real trace that
whole-brain averaging dilutes to nothing. α concentrates (encoding) in
temporal-limbic cortex (MTL/lateral-temporal/cingulate, large but per-region
underpowered). **Caveat:** the absolute test is structurally underpowered for
K≤5 systems (OFC ρ=+0.56 cannot reach BH at K=5), so "only cingulate clears it" is
partly the cingulate being best-sampled (K=8); OFC's evidence is the more powerful
demeaned localization. The recurring hub across bands is the **cingulate**
(encoding in α/low-γ, inference in β).

---

## What was done (plain)

The consolidation arc (audit_103) decomposed the β offline trace into an
**encoding** part (`e = D_taskLearn − D_preA`) and an **inference-specific** part
(`f = D_taskTest − D_taskLearn`), and found β uniquely keeps an offline mark of
the inference part beyond encoding (`T_infspec·e` p=0.0068, β-only). This asked
the matching anatomical question, **reusing the locked OFC-localization machinery
verbatim** (audit_83 rank concordance → per-system per-patient-demeaned endpoint
aggregation → matched-strength R surrogate), with the per-pair vector swapped from
the standard trace to encoding / inference. Four targets, all through the same
null: `standard = concordance(g,p)`, `encoding = concordance(e,p)`,
`inference_raw = concordance(f,p)`, `inference_pe = concordance_partial(f,p|e)`.

## Evidence — the β dissociation (demeaned localization, R=1000)

Per-system matched-strength p at **R=1000** (positive = accumulation above the
distributed baseline). `inference_pe` = the headline (inference controlling
encoding); BH q across the 9 systems within β:

| target | top β system | MS p (epi-incl) | BH q (epi-incl) | epi-excl q | shaft-collapse (R200) |
|---|---|---|---|---|---|
| standard | **OFC** | 0.002 | **0.010 ✓** | 0.010 ✓ | p=0.025 |
| encoding | **OFC** | 0.001 | **0.010 ✓** | 0.010 ✓ | p=0.030 |
| inference_pe | **cingulate** | 0.005 | **0.050 ✓** | 0.060 (borderline) | p=0.020 / 0.020 |

- **Cingulate β inference-mark clears BH at R=1000** (q=0.050 epi-include,
  borderline q=0.060 epi-exclude) — the R=200 "q≈0.15" was a 1/201-floor artifact.
  Strength_dev +0.04 (not a hub), K=8/10 (better sampled than OFC's K=5), 5/8
  patients positive, cohort median positive under every leave-one-out drop,
  survives the shaft-collapse control that killed the earlier "distributed ring".
  OFC drops to p=0.13–0.15 for the inference component (the dissociation).

### The double dissociation (the strong evidence) — cingulate role-flip by band

| band | cingulate ENCODING | cingulate INFERENCE_pe |
|---|---|---|
| α     | p=0.005–0.010 ✓ (q≤0.075) | p=0.89–1.0 (strongly **anti**) |
| β     | p=0.12 (n.s.)             | **p=0.010–0.020 ✓** |
| low-γ | p=0.005 ✓ (q=0.025)       | p=0.74–0.82 (anti) |

The cingulate is an **encoding** hub in α/low-γ but flips to an **inference** hub in
β — exactly the arc's dissociation (α = encoding-only; β = encoding + inference),
recovered independently by anatomy. Encoding/standard sit in OFC throughout; α
encoding additionally engages lateral-temporal + MTL (q=0.033).

## Absolute within-system trace (audit_112) — a focal trace under a null whole-brain trace

The demeaned localization tests *relative* concentration (system above the patient's
own mean). To test whether a system carries an **absolute** trace — the steelman of
"α/low-γ are not clear overall but locally strong" — audit_112 restricts the cohort
split-baseline trace to system-incident pairs (Wilcoxon obs > matched-strength
surrogate median, no demeaning), per band×target×system.

Cohort whole-brain trace (arc null) clears only **α** (standard p=0.005, encoding
p=0.014) and **β** (all three); **low-γ is null brain-wide** (p=0.12 / 0.12 / anti).
Yet within-system:

| band | target | system | obs ρ | n_pos | Wilcoxon p | BH q |
|---|---|---|---|---|---|---|
| low-γ | encoding | **cingulate** | +0.40 (excl +0.51) | **8/8** (excl 7/8) | 0.004 | **0.035 ✓** (excl 0.070) |
| α | standard | PFC | +0.13 | 7/7 | 0.008 | 0.070 · |
| α | encoding | cingulate / MTL / lat-temp | +0.35–0.38 | 7/8, 5/5 | 0.03–0.06 | 0.11–0.12 · |
| β | encoding | MTL / lat-temp / OFC | +0.27–0.56 | 4–5/5 | 0.03–0.09 | 0.14–0.28 |
| β | inference_pe | cingulate | +0.11 | 5/8 | 0.16 | 0.28 |

- **The standout: a genuine focal low-γ ENCODING trace in the cingulate** (ρ=+0.40,
  8/8, q=0.035) where the whole-brain low-γ trace is null — vindicates the
  "focal-trace-hidden-by-averaging" reading. **Survives epi-exclusion** (ρ rises to
  +0.51; q drifts to 0.070 only from lost power) → not epileptic-contact driven.
  Most other low-γ localization hotspots (insula/OFC/sensorimotor/PFC) do **not**
  survive this absolute test → they were demeaning structure.
- **Power caveat (decisive for reading this table):** the absolute Wilcoxon is
  structurally underpowered for K≤5 systems — OFC has a *large* β encoding trace
  (ρ=+0.56) but cannot reach BH at K=5 (min achievable q≈0.28). So "only cingulate
  clears" is partly the cingulate being best-sampled (K=8); OFC's evidence is the
  demeaned localization, not this test. The β *inference* cingulate result is a
  **concentration** (small absolute ρ=+0.11), not a large absolute trace — the
  inference component is small in magnitude.

## Cross-checks (anti-hallucination, all passed)

- standard β/system/OFC `M_obs` reproduces the locked audit_83 value **exactly**
  (+6.133e5) → this pipeline == the canonical OFC localizer.
- `inference_pe` per-pair partial reproduces the arc's formula partial Spearman to
  **5.6e-17** (machine precision); `standard`/`inference_raw` sums reproduce
  scipy Spearman to ~2e-3 (tie-correction). The localization is tied to the
  validated arc numbers.

## Duration / length-equalized control (audit_113) — content, not data-amount

`task_test` is **1.35–2.53× longer than `task_learn` in all 10 patients** (median
1.57×), so the inference vector `f = D_taskTest − D_taskLearn` pairs a `task_test`
cophenetic estimated from ~1.6× more Welch segments against a noisier `task_learn`
one — a systematic data-amount asymmetry that could inflate the cingulate inference
hotspot independent of content. audit_113 truncates `task_test` to `task_learn`'s
sample count (3 window placements head/center/tail), recomputes ONLY that phase's
canonical cophenetic, and reruns the exact audit_110 localization. Two built-in
checks pass: the full-length recompute reproduces the cached cophenetic **exactly**
(Spearman ρ=1.000000, max|Δ|≈1.5e-4 = float32 rounding, on all 3 validation
patients spanning the ratio range); and the **encoding** target (no `task_test`
dependence) is bit-identical full-vs-matched, as it must be.

**Verdict: the β cingulate inference hotspot SURVIVES length-matching — content-
driven, not a duration artifact.** Of the 8 patients that sample the cingulate, all
**8/8 keep the sign** of their per-patient demeaned inference value; the cohort
median stays positive at **0.61× full** (+4.36e5→+2.67e5 epi-include;
+4.90e5→+2.97e5 epi-exclude); per-patient full-vs-matched values correlate at
**ρ=0.74–0.76**; and a paired Wilcoxon finds matched **indistinguishable from full
(p≈0.95)**. The 5/8-positive cohort split is reproduced exactly. The OFC standard
hotspot likewise survives (0.83–0.87× full, sign 5/5, ρ=0.90–1.0); OFC encoding is
identical by construction. **Direction of bias:** truncation discards data →
noisier `D_TT'` → biases the concordance toward zero, so the ~40% point-estimate
attenuation is the expected (statistically non-significant) noise penalty, not a
collapse; survival under the noisier matched length is conservative evidence.

**Stage limit:** observed-statistic robustness only — it does NOT regenerate the
matched-strength surrogate on truncated `task_test`, so it yields no new BH q at
matched length. Given the observed value barely moves (ρ_fm=0.74, paired p≈0.95,
8/8 sign), a surrogate-on-truncated stage would be expected to reproduce q≈0.05;
it is available if a referee insists. Data: `length_control.csv` (540 rows).

## Honest limits

- **β cingulate inference clears BH only at the boundary** (q=0.050 epi-include,
  q=0.060 epi-exclude) — established but marginal, not the decisive q=0.010 the
  OFC encoding/standard hotspots reach. The corroboration is the 4-condition
  consistency + the double dissociation, not just the corrected p.
- **Duration controlled, not surrogate-re-nulled (audit_113).** `task_test` runs
  1.5–2.5× longer than `task_learn`; the cingulate inference (and OFC) hotspots
  survive length-matching `task_test` (8/8 sign-preserved, matched indistinguishable
  from full p≈0.95), excluding the data-amount confound — but this is the observed-
  statistic stage; the matched-strength surrogate was not regenerated on truncated
  data.
- **System scale only.** Region-level cingulate sub-regions are all K≤2 (implant
  ceiling) and individually n.s. — no anterior-vs-posterior cingulate claim; the
  signal is the pooled cingulate **system** (same as the OFC verdict's system
  scale).
- **No behavior** → this is the offline residue of the inference-specific
  *reorganization*, not inference *success* (deferred behavioral TC1).
- A low-γ → PFC inference cell appears at contact level (p=0.005) but weakens
  under shaft-collapse (p=0.010) and localizes a cohort inference component the
  arc found non-significant → flagged, not interpreted.

## Figures (data/audit/inference_localization/figures/)

- `brain_inference_dissociation_beta.pdf` — 3 glass-brain rows (standard /
  encoding / inference), β: the OFC→OFC→cingulate migration; grey = distributed
  implant footprint, glow = system accumulation, ring = BH q<0.05.
- `systems_inference_dissociation.pdf` — per-system normalised accumulation
  (distributed field → 0, hotspots pop) + the cingulate's α/β/low-γ role-flip
  crossover.

## Next steps

1. ✅ DONE — β `inference_pe` finalized at R=1000 (q=0.050 incl / 0.060 excl);
   ✅ DONE — within-system absolute trace (audit_112): low-γ cingulate encoding
   focal trace (q=0.035, epi-robust). Refresh the figure glow from the R=1000 CSV.
2. (Deferred, needs behavior) does the cingulate β inference-mark predict
   inference success / the symbolic-distance effect? — behavioral TC1.
3. ✅ DONE — length-equalized duration control (audit_113): `task_test` is
   1.35–2.53× longer than `task_learn` (the confound was real, NOT moot), but the
   cingulate inference and OFC hotspots survive truncating `task_test` to
   `task_learn` length (8/8 cingulate-sampling patients sign-preserved; matched
   indistinguishable from full, p≈0.95) → duration confound excluded (observed
   stage). Low-γ cingulate-encoding focal trace (no whole-brain trace) remains a
   candidate result in its own right — worth a dedicated within-region writeup if
   pursued.
