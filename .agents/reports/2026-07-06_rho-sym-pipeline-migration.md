---
name: rho-sym-pipeline-migration
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-06
supersedes_claims:
  - "bare ρ_split is the trace estimator of record (audit_63) — replaced by ρ_sym"
pointers:
  - scripts/01_compute/audit/audit_149_estimator_robustness_regate.py
  - scripts/01_compute/audit/audit_150_rho_sym_gate.py
  - scripts/01_compute/audit/audit_151_localization_rhosym.py
  - data/audit/rho_sym_gate/
  - data/audit/localization_atlas_rhosym/
---

# ρ_sym pipeline migration — COMPACTION HANDOFF (resume the cascade here)

## Head

The cophenetic trace estimator is now **ρ_sym** (canonical), replacing bare
**ρ_split**. ρ_sym = ½[ρ(A→task,B→rest) + ρ(B→task,A→rest)] — it averages the two
equally-valid split-half arm assignments, removing the **arbitrary-half artifact**
that made ρ_split sign-unstable for near-zero patients (20/60 cells flip under A/B
swap; audit_149). **Both flagship pillars are re-confirmed under ρ_sym and did not
move:** the cohort gate (α/β CLEAR, rest fail) and β→OFC localization (carrier
q=0.025, sensorimotor/PFC depleted). New audits only — **audit_63/83 were NOT
edited** (user rule). All surrogate compute is now **numba-JIT'd (98× faster,
bit-identical)**.

---

## Two standing rules (locked 2026-07-06, memories saved)

1. **ρ_sym is the canonical estimator; forget bare ρ_split** (keep it only as a
   supplement column). [[feedback_rho_sym_canonical_estimator]]
2. **Always numba for heavy compute** (JIT the swap loop; rng arrays pre-generated
   outside → bit-identical; no `cache=True`; warm-compile before Pool fork).
   [[feedback_numba_for_surrogates]]
3. **Never edit old audits — new audit number per step, new output dir.**

## What's DONE (verdicts estimator-invariant — ENTIRE CASCADE COMPLETE 2026-07-06)

| step | script | output | result |
|---|---|---|---|
| estimator-robustness | `audit_149` | `data/audit/estimator_regate/` | split=sym gate, ZERO flips; full-baseline REJECTED (Pat_02 catastrophic cancellation 0.775−0.767) |
| **cohort gate** (ρ_sym) | `audit_150` | `data/audit/rho_sym_gate/` | **α p=.024, β p=.032 CLEAR**; δ/θ/low-γ/high-γ fail. split ref reproduces published (α.003/β.005) |
| **β→OFC localization** (ρ_sym) | `audit_151` | `data/audit/localization_atlas_rhosym/` | **OFC carrier q=.025** (both epi); **sensorimotor+PFC depleted q≤.033**; cingulate/occipital/lat-temp secondary carriers |
| **consolidation arc** (ρ_sym) | `audit_152` | `data/audit/consolidation_arc_rhosym/` | **β T_infspec_pe p=.0098 (LO-P15 .020), β-ONLY** (all other bands min p=.246); **α/β T_learn p=.014/.032** (learning-phase own trace); T_test α/β p=.032. arm1=ρ_split reproduces audit_83 full-graph to 1e-16 (CROSS-CHECK PASS). N2 climax estimator-proof. |
| **cross-phase taxonomy** (ρ_sym) | `audit_153` | `data/audit/cross_phase_taxonomy_rhosym/` | composition holds: **β trace-dominant** (comp_trace .263, mover_frac .374 vs anchor .63; every other band anchor-dominated, α anchor .905); anchor≠hubness (β r(anchor,strength)=.14). trace-guard **60/60 bit-exact** (arm1 ρ_split vs locked audit_63); ρ_sym vs audit_150 max\|Δ\|=1e-16 |
| **per-node decomposition** (ρ_sym) | `audit_154` | `data/audit/per_node_trace_decomposition_rhosym/` | carrier/anti split holds; **node-level ρ_sym↔ρ_split ρ=.96/.95/.98** (β/α/γₗ), sign-concordant 9/10 β; anti-node property preserved (β is_epi anti-enriched p=.021, MTL 2nd); non-tracers = anti-driven (Pat_10 44 anti=reset, Pat_15 16 anti=null) |

**Nothing moved.** Every verdict is bit-for-bit the same taxonomy under ρ_sym; only the
arbitrary-half artifact is gone. The compute cascade (Section A below) is COMPLETE;
old audits 63/83/103/105/144 untouched (audit_83_localization shows a pre-existing
2026-06-25 library-promotion refactor, NOT this session).

Per-patient reporting rule: report **ρ_sym + split-uncertainty ½|ρ_AB−ρ_BA|**; mark
`|ρ_sym| < 1 SE` **"undetermined"** (near-zero labels = estimator noise, not biology
— borderline third of cohort; e.g. Pat_13/Pat_15 β).

## CASCADE STATUS

**A. Downstream ρ_split-based audits → new ρ_sym audits — ✅ DONE 2026-07-06.**
All three ported with the sym concordance `s_sym = ½[concordance(A→task, B→rest) +
concordance(B→task, A→rest)]`, cached surrogates reused (no regeneration), old audits
untouched. Results in the table above — every verdict held.
  - ✅ **consolidation_arc** → `audit_152` (f=inference-specific contrast is arm-invariant;
    only its persistent axis symmetrizes — a weaker de-bias than T_test/T_learn, noted).
  - ✅ **cross_phase_taxonomy** → `audit_153` (symmetric composition = avg of both arms;
    kept the bit-exact arm1 ρ_split guard vs locked audit_63 as the pipeline anchor).
  - ✅ **per_node_trace_decomposition** → `audit_154` (audit_151 template; reuses R=200 cache).
  - ⬜ (optional robustness, NOT blocking) β→OFC LOO + shaft-collapse under ρ_sym at R=1000
    — audit_151 already confirmed the verdict at R=200 (q=.025); R=1000/LOO/shaft polish
    would migrate the ANATOMY_LEDGER headline number like-for-like (currently R=1000 q=.009–.013).

**B. Ledgers / manuscript (swap ρ_split → ρ_sym numbers; PI / headline-chat domain):**
  - `VERDICT_LEDGER.md`, `ANATOMY_LEDGER.md` — update trace/localization numbers to
    the ρ_sym values above (verdict tags UNCHANGED).
  - `2026-06-26_per-band-phenomenology-vision.md` — estimator = ρ_sym; add the
    "undetermined" per-patient tier; keep taxonomy.
  - Add the **estimator-robustness supplement** (audit_149 split=sym invariance;
    full-baseline rejected) — the referee-facing defense that the trace is not an
    arbitrary-half artifact.
  - [[interpatient_variability_resolution_2026_07_04]] heterogeneity framing stays
    (emergent/global, not coverage) + now "near-zero labels are estimator-undetermined".

## Key numbers to carry forward
- Gate (ρ_sym, R=200): α p=.024, β p=.032 CLEAR; δ .080, θ .784, low-γ .080, high-γ .080 (all fail).
- β→OFC: carrier q=.025 (include & exclude); sensorimotor depleted q_lo=.017/.033; PFC depleted q_lo=.017/.025.
- numba shuffle: 882 ms → 9 ms (98×), max|nb−py|=0.0.
- Infra: `_swap_loop` numba pattern + MP Pool in audit_150/151; cached surrogates reused in audit_151.
