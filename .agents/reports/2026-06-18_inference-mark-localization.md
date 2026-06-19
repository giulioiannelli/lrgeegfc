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

# Where the β inference-mark sits — directionally cingulate, but the localization does NOT survive the duration control (memorizing→OFC does)

**Head.** Within the β band the offline trace **splits by what it carries**: the
standard trace and the **encoding** component (the memorised pairs) over-accumulate
in **orbitofrontal cortex** (OFC), but the **inference-specific** component — the
figured-out relations, with encoding statistically removed — over-accumulates in
the **CINGULATE**. This is an *accumulation hotspot on a distributed trace*
(per-patient demeaned, [[feedback-localization-is-overexpression-not-container]]),
not a container. **At full length both sides were BH-significant (R=1000):** encoding/standard → OFC
(q=0.010), inference → cingulate (q=0.050 epi-include, q=0.060 epi-exclude — already
the marginal piece). **A duration control then SPLIT these two verdicts** (`task_test`
runs 1.35–2.53× longer than `task_learn` in every patient — a real, one-directional
data-amount asymmetry feeding `f = D_taskTest − D_taskLearn`; audit_113/113b below):

- **Memorizing → OFC SURVIVES.** Truncating `task_test` to `task_learn`'s length and
  re-running the matched-strength null, the **standard** trace's OFC hotspot still
  clears (q=0.025 incl / 0.050 excl), and **encoding** is duration-invariant by
  construction (`e = D_taskLearn − D_preA` has no `task_test` term → bit-identical;
  its full-length q=0.010 stands). The memorizing localization is duration-robust.
- **Inference → cingulate does NOT survive.** At matched length the cingulate
  inference component sits **inside** its strength null (p=0.124 incl / 0.179 excl,
  q=0.42/0.76) and is no longer even the top-ranked inference system (insula/MTL/
  occipital rank above it; none clear BH). The cingulate stays the **directional**
  best guess (audit_113 stage-1: 8/8 cingulate-sampling patients keep the sign,
  per-patient values track full at ρ=0.74), but its STATISTICAL significance is not
  separable from `task_test`'s longer recording.

**Honest verdict: inference → cingulate is a directional hint, NOT an established
q<0.05 localization once the duration confound is controlled; only memorizing → OFC
is duration-robust.** The double dissociation's *encoding* arms (the cingulate
carrying encoding in α/low-γ — both ⊥ `task_test`) are unaffected; its β-*inference*
arm is downgraded to directional. **Upstream flag:** the consolidation arc itself
([[arc-inference-consolidation-2026-06-18]]) is built on the same `f = D_TT − D_TL`
and has NOT yet been duration-tested — the whole-brain "β consolidates an
inference-specific component" headline now needs the same length-matched control
before it can be called duration-robust (open: a stage-2-arc).

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
  **⚠ FULL-LENGTH result — superseded for the inference row.** The duration control
  (audit_113b, below) puts inference → cingulate back inside its null at matched
  length (p=0.124 incl, q=0.42); the standard/encoding OFC rows survive. Read this
  table as the full-length analysis, NOT the duration-controlled verdict.

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

## Duration / length-equalized control (audit_113 + audit_113b) — it SPLITS the verdict

`task_test` is **1.35–2.53× longer than `task_learn` in all 10 patients** (median
1.57×), so `f = D_taskTest − D_taskLearn` pairs a `task_test` cophenetic from ~1.6×
more Welch segments against a noisier `task_learn` one — a one-directional data-amount
asymmetry that could inflate the inference localization independent of content. (The
prior handoff called this control "moot by band/target specificity"; the durations
show it is **not** — the confound is real.) Two stages.

**Stage 1 — observed robustness (audit_113).** Truncate `task_test` to `task_learn`'s
sample count (3 placements head/center/tail), recompute ONLY that phase's canonical
cophenetic, rerun the audit_110 localization. Built-in checks pass: full-length
recompute reproduces the cached cophenetic exactly (Spearman ρ=1.000000, max|Δ|≈1.5e-4
= float32 rounding; 3 validation patients), and **encoding** (no `task_test` term) is
bit-identical full-vs-matched. Result: the cingulate inference **DIRECTION** is
preserved — 8/8 cingulate-sampling patients keep the sign, per-patient values track
full at ρ=0.74, cohort median 0.61× full. A paired Wilcoxon finds matched
indistinguishable from full (p≈0.95) — but at K=8 that test is underpowered: it shows
the point estimate did not significantly *drop*, NOT that it still beats a null. Data
`length_control.csv`.

**Stage 2 — matched-length strength null (audit_113b), the decisive test.**
Regenerate the matched-strength surrogate ON the truncated `task_test` (canonical
4-cycle ±δ, swap 20, seed 20260511; non-canonical cache path so the canonical
ensemble is untouched), keep the cached full-length surrogates for the four unchanged
phases, recompute per-system matched-strength p + BH q at matched length (R=200, head
placement). **The two sides split:**

| target (β, matched length) | KEY system | matched-strength p | BH q | survives? |
|---|---|---|---|---|
| standard | OFC | 0.005 / 0.005 | **0.025 / 0.050 ✓** | YES |
| encoding | OFC | 0.010 (R=200 floor) | full-length q=0.010 (⊥`task_test`) | YES (by construction) |
| inference_pe | cingulate | **0.124 / 0.179** | **0.42 / 0.76 ✗** | **NO** |

- **Memorizing → OFC SURVIVES the matched-length null.** Standard-trace OFC still
  clears (q=0.025 incl / 0.050 excl); encoding is duration-invariant (its full-length
  R=1000 q=0.010 is unchanged — the q=0.100 shown at R=200 is only the floor, not a
  real move). The pipeline produces significance when warranted (OFC standard +
  occipital both clear) → the inference failure is NOT a dead test.
- **Inference → cingulate does NOT survive.** At matched length the cingulate
  inference component is **inside** its strength null (p=0.124 incl / 0.179 excl),
  ranks only 3rd/2nd among inference systems (insula p=0.025, MTL p=0.040 sit above
  it), and **no system clears BH for inference**. The full-length q=0.050 was
  **duration-assisted** — its significance depended on `task_test`'s extra data
  sharpening `D_TT`. R won't rescue it (p=0.12 is far from the 1/(R+1) floor); nor
  will placement (the strongest placement's observed ≈0.6× full cannot reach the BH
  threshold, which needs ≈the full-length value). Data
  `lenmatched_null_R200_{include,exclude}.csv`.

**Reading stages 1+2 together:** the cingulate is the inference component's
*directional* home (sign-robust 8/8, ρ=0.74 with full), but its significance is not
separable from the longer test recording — so **inference → cingulate is a directional
hint, not an established localization.** Whether the full-length q=0.050 was a pure
duration artifact or a real-but-data-hungry signal cannot be distinguished here;
either way the clean q<0.05 claim does not hold at matched length. (R=1000 would
sharpen the numbers but cannot change this negative; optional.)

## Honest limits

- **β inference → cingulate does NOT survive the duration control (the headline
  caveat).** Full length it was already marginal (q=0.050/0.060); the matched-length
  strength null (audit_113b) puts it inside its null (p=0.124/0.179, q≥0.42) — its
  significance was duration-assisted. It remains directionally robust (8/8 sign,
  ρ=0.74) but is a **directional hint, not an established localization**. Only
  memorizing → OFC is duration-robust (standard OFC q=0.025–0.050 at matched length;
  encoding ⊥`task_test`).
- **The consolidation arc shares this confound (open).** The arc's
  `T_infspec = concordance(f, p | e)` is built on the same `f = D_TT − D_TL`; it has
  NOT been duration-tested. The whole-brain "β consolidates an inference-specific
  component" headline needs the same length-matched null (a stage-2-arc) before it is
  duration-robust — flagged, not yet run.
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
3. ✅ DONE — duration control, BOTH stages. Stage 1 (audit_113, observed): direction
   robust (8/8 sign). Stage 2 (audit_113b, matched-length strength null, R=200): the
   two sides split — **memorizing → OFC survives** (standard q=0.025/0.050; encoding
   ⊥`task_test`), **inference → cingulate does NOT** (p=0.124/0.179, q≥0.42, not even
   top inference system). The inference→cingulate localization is downgraded to a
   directional hint. The confound was real, not moot.
4. ⚠ OPEN (raised by #3) — duration-test the consolidation ARC itself (audit_103,
   same `f = D_TT − D_TL`): a stage-2-arc whole-brain matched-length null. Until run,
   the "β consolidates an inference-specific component" headline is duration-uncontrolled.
5. Optional: R=1000 matched-length null (sharpens audit_113b numbers, cannot change
   the negative — cingulate p=0.12 is far from the floor).
6. (Deferred, needs behavior) behavioral TC1 — does the β cingulate inference value
   predict inference success? Low-γ cingulate-encoding focal trace (audit_112, ⊥
   `task_test`, unaffected by all this) remains a candidate result on its own.
7. ⚠ Figures `brain_inference_dissociation_beta.pdf` / `systems_inference_dissociation.pdf`
   still show the inference→cingulate BH ring from the full-length result — regenerate
   without the inference-row significance ring (or annotate "full length") to match the
   duration-controlled verdict.
