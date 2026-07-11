---
name: ms-relative-scale-sweep-verdict
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-11
updates: .agents/reports/2026-07-07_r1-tau-multiscale-reexamination.md
  (executes the there-deferred audit_159; upgrades the α lead from
   observed-only/marginal to matched-strength-clean + LOO-robust)
pointers:
  - scripts/01_compute/audit/audit_172_ms_relative_scale_sweep.py
  - scripts/01_compute/audit/audit_172b_scale_sweep_figure.py
  - data/audit/tau_sweep_ms_gate/cohort_gate_s.csv
  - data/audit/tau_sweep_ms_gate/audit_172_scale_sweep_verdict.pdf
  - .agents/guides/task-persistence-investigation/2026-07-11_ms-relative-scale-sweep.md
  - scripts/01_compute/audit/audit_121_tau_sweep_cophenetic_trace.py
memory: tau_sensitivity_trace_2026_06_22
---

# MS-relative scale sweep — verdict: α has a characteristic mesoscale; β is scale-broad; the rest is collapse

## Head

We ran the referee audit_121 lacked: the canonical **ρ_sym gate (audit_150)
swept across the diffusion scale `s ≡ τ·λ_max` with a matched-strength null
rebuilt at every s** (audit_172; `s = 1` reproduces the certified gate
bit-for-bit, max|Δ| = 1×10⁻¹⁶). **The α-band trace is genuinely scale-tuned:
it strengthens at a mesoscale `s ≈ 1.5` (cohort gate `p 0.024 → 0.0098`,
z_median 2.2 → 3.6, leave-one-out worst-case p = 0.0195), on a broad significant
plateau `s ∈ [0.5, 2.3]`, and it *weakens* past the Fiedler scale — so it is not
collapse.** β clears at *every* collapse-free scale (strongest `s ≈ 1.86`,
p = 0.003) — **scale-broad**, a robustness reframe, not a new gate. **δ and
low_γ "clear" only as s climbs toward collapse** (monotone z, strengthening past
Fiedler, δ being the known SOZ-anatomy band): **rejected as task-traces** by the
neg-control discipline. θ and high_γ stay null throughout. Net: there is a real
**band ↔ scale structure — β broad, α mesoscale-peaked, the rest scale-less** —
and reading α at τ_min *under-reads it*. The Methods claim that "the τ-sweep is
degenerate … no hierarchy remains to read" is false for α (true enough for β).

## What was run

audit_172 generalises the audit_150 gate from `s = 1` to a 16-point s-grid
(`s ∈ [0.5, 12]`, dense near 1–3), all 10 patients × 6 bands. Each phase matrix
is eigendecomposed once; `D_coph(s)` is reformed from the eigenpairs at every s;
the trace statistic is the canonical ρ_sym; the null is the canonical
matched-strength 4-cycle ±δ swap (R = 200), rebuilt at every s. Seeds and swap
order are audit_150's exactly, so the **`s = 1` column reproduces
`per_patient_per_band.csv` obs and surrogate mean to machine precision** — the
anchor that certifies every `s ≠ 1` reading as real scale physics, not code
drift. 5-point preamble + design:
`.agents/guides/task-persistence-investigation/2026-07-11_ms-relative-scale-sweep.md`.

## Results — `gate_p_sym(s)`, collapse-free window (s < Fiedler s_F) marked `*`

```
band        s_F |  0.85   1.00   1.51   1.86   2.29   2.82  | coarse (≥ s_F)
alpha       3.8| .024*  .024*  .010*  .014*  .014*  .053* | .097 .161  weakens
beta        3.2| .019*  .032*  .005*  .003*  .014*  .032* | .032 .010  broad
delta       3.7| .053*  .080*  .065*  .024*  .042*  .032* | .007 .014  STRENGTHENS→
low_gamma   2.3| .042*  .080*  .138*  .080*  .097   .042  | bouncy, no plateau
theta       3.4| .461*  .784*  .754*  .688*  .577*  .615* | null ✓
high_gamma  1.8| .097*  .080*  .116*  .097   .116   .080  | null ✓
```

Leave-one-out (worst-case p over the 10 drops) at each candidate's collapse-free
best s: **α @ s=1.51 → 0.0195 (robust)**, β @ s=1.86 → 0.0059 (robust),
α @ s=1.00 → 0.0488 (barely). δ @ s=5.24 → 0.0137 (robust — but coarse, see below).

## The discriminator — genuine (α) vs collapse (δ, low_γ)

`z_median(s)` separates them by **shape**, not by significance:

- **Scale-tuned trace = PEAK then decline** inside the collapse-free window. α:
  z 2.2 → **3.6 @ s≈1.5** → decline; gate weakens past Fiedler. A real
  characteristic scale.
- **Collapse / anatomy channel = MONOTONE climb** toward collapse. δ:
  z 0.7 → 1.9 → 2.9 → 4.6, gate *strengthening* as s → collapse; low_γ, high_γ
  the same. δ is the band with a strong anatomical slow mode (the SOZ marker
  lives at δ, coarse τ). Matched-strength does **not** rescue this into a trace:
  MS destroys the anatomy geometry, so a phase-invariant **anchor** surviving the
  ρ_sym differencing at coarse s (where the differencing degrades toward
  collapse) still beats the strength-matched surrogate without being task-driven.
  Rejected per the scope's neg-control discipline; would need a dedicated
  within-phase / anatomy control before any δ claim.

## Verdict

| band | call | evidence |
|---|---|---|
| **α** | **PASS — genuine mesoscale trace** | broad plateau s∈[0.5,2.3]; peak gate 0.0098 @ s≈1.5; z peak 3.6; LOO 0.0195; weakens past Fiedler |
| β | scale-broad (reframe) | clears every collapse-free s; strongest 0.003 @ s≈1.86; LOO 0.006 |
| δ | REJECT (collapse/anatomy) | monotone z-climb; strengthens into collapse; known SOZ slow-mode band |
| low_γ | no clean trace | bouncy, monotone-z; consistent with audit_121 "not rescued at any τ" |
| θ, high_γ | null | stay non-significant at every s (clean negative controls) |

## Manuscript implications (PI decision — hardened prose NOT edited here)

- The **"τ-sweep is degenerate / no hierarchy remains to read"** language
  (methods.tex §Decomposing…, §LRG framework) is **contradicted for α** and
  should soften to something like: *the flagship β trace is scale-broad and read
  at τ_min; the α trace is mildly mesoscale-favouring (Fig. Sx), the only band
  with a scale it does not share with its anatomy.* β stays at τ_min regardless.
- This upgrades the r1-reexam (2026-07-07) status: the α mesoscale is **no longer
  observed-only/marginal** — it is matched-strength-clean and LOO-robust. Whether
  it enters the paper (as an α "scale though no place" point, strengthening the
  N1→N2 arc) is a writing call, not a null call.
- **Do not** headline δ/low_γ coarse-τ significance — it is collapse/anatomy.

## Part B — encoding vs inference over scale (audit_173)

Same construction applied to the decomposition functionals (`T_learn` = encoding
echo; `T_infspec_pe` = inference-specific | encoding), 5-phase arc, over a dense
s-grid `s ∈ [0.5, 8]`. Generalises audit_157 (3 τ) exactly as audit_172 ⊃
audit_150, reusing the **cached seed-20260511 surrogate eigenpairs**. Anchors:
arm1 at s ∈ {1.0, 2.610, 6.813} == audit_103d (ρ_split); sym at s=1 == audit_152
(ρ_sym fine) — **both to 1×10⁻¹⁶**.

`cohort gate p(s)` (collapse-free `*` = s < Fiedler):

```
T_infspec_pe   s=0.5   1.0    2.1    2.6    3.8    8.0
beta          .005*  .010*  .003*  .005*  .014   .003     ← β-ONLY, every scale
(all others)   null    null   null   null   null   null    (0.25–0.98)

T_learn        s=0.5   0.85   1.0    2.6    3.8    8.0
alpha         .042*  .010*  .014*  .024*  .097   .246     ← FINE, dies at coarse
beta          .010*  .032*  .032*  .024*  .161   .032     ← fine/meso
delta         .024*  .014*  .042*  .014*  .010   .003     ← STRENGTHENS → collapse (REJECT)
```

**Verdict:**
- **Inference-specific (T_infspec_pe) is β-only at every scale** — clears at τ_min
  (0.010) and everywhere else (best 0.003 @ s≈2.1); **no other band ever clears**.
  Scale-broad, mildly mesoscale-favouring → confirms audit_157 with the full curve
  and vindicates the τ_min reading.
- **Encoding echo (T_learn) is α+β and fine-scale-favouring** — α strongest at
  s≈0.85 (0.010), **weakening monotonically to null past Fiedler** (0.246 @ s=8);
  β fine/meso.
- **Scale dissociation:** encoding lives fine and dissolves at coarse; inference
  persists broadly. A component-level echo of Part A's α-mesoscale / β-broad split
  (*what was encoded is finest-scale; what was computed sits slightly coarser and
  persists*).
- **δ encoding strengthens with s** (0.014→0.007→0.003) — same collapse/anatomy
  signature as Part A; **rejected**.
- **No paper change:** β inference and α/β encoding all clear at τ_min; the same
  components are carried by the same bands across the collapse-free window.

Figure: `data/audit/arc_scale_sweep_rhosym/audit_173_arc_scale_verdict.pdf`.

## Provenance

- audit_172 (engine) + audit_172b (figure), 2026-07-11, R=200, |s|=16, 553 s.
- audit_173 (arc, engine) + audit_173b (figure), 2026-07-11, R=200 cached
  seed-20260511 eigs reform-only, |s|=17, 943 s. Anchors arm1==103d, sym@1==152.
- Outputs: `data/audit/tau_sweep_ms_gate/{per_patient_per_band_s,cohort_gate_s}.csv`,
  `audit_172_scale_sweep_verdict.pdf` (+ `.md`), `README.md`.
- s=1 anchored bit-exact to audit_150 (BASE_SEED 20260706 + cell_idx, same swap).
- Estimator ρ_sym; null matched-strength 4-cycle ±δ (SWAP_FACTOR 20), reformed at
  each s from cached eigenpairs. Cohort gate Wilcoxon(obs − surr_p50) one-sided.
- Does NOT modify audit_121/150; new pipeline, new output dir.
