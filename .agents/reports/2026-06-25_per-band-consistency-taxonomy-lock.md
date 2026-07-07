---
name: per-band-consistency-taxonomy-lock
type: report
era: "IMCOH_ABS × COHORT_N10"
status: complete
created: 2026-06-25
pointers:
  - scripts/01_compute/audit/audit_83_localization_matched_strength.py
  - scripts/01_compute/audit/audit_83b_per_band_taxonomy.py
  - data/audit/localization_atlas/per_band_taxonomy_verdict.csv
  - data/audit/localization_atlas/matched_strength_allbands_*.csv
  - .agents/reports/2026-06-25_snr-band-taxonomy-handoff.md
  - .agents/reports/2026-06-25_per-node-trace-anatomy-and-heterogeneity.md
---

> **HEAD.** The per-band consistency taxonomy is **LOCKED** with the validated
> audit_83 machinery (cohort-median per-system trace × matched-strength × both
> epi modes × shaft-collapse, all 6 bands): **consistent {β→OFC, γ_l→PFC} |
> patient-specific {α, δ, γ_h} | absent {θ}**. It reproduces the hypothesised
> 3-tier structure — but now EARNED, not from the crude n=5 check that had
> mis-read it. Two honest refinements emerged: (1) **every band's carrier
> survives shaft-collapse**, so the tier is set by whether a *net cohort trace*
> exists to localise (the gate), NOT by localization robustness; (2) **θ and γ_h
> DO localize** (θ→MTL, γ_h→parietal, both shaft-robust) but on a band with ~zero
> net cohort trace — a spatial redistribution, not a trace — so they are not
> "consistent". β→OFC is fully locked (prior LOO); **γ_l→PFC still owes a
> leave-one-patient-out** before it is as hard as β.

---

## 1. The locked verdict

`data/audit/localization_atlas/per_band_taxonomy_verdict.csv`:

| band | tier | trace_level | gate (n>p95 / med ρ / Wilcoxon p) | robust carrier | robust depleted |
|---|---|---|---|---|---|
| **β** | **consistent** | cohort | 7/10 / +0.221 / 0.005 | **OFC** | PFC, sensorimotor |
| **γ_l** | **consistent** | cohort | 5/10 / +0.083 / 0.116 | **PFC** | occipital |
| **α** | patient-specific | cohort | 5/10 / +0.105 / 0.002 | — (none survives) | — |
| **δ** | patient-specific | subset | 4/10 / +0.008 / 0.278 | — | — |
| **γ_h** | patient-specific | subset | 4/10 / +0.000 / 0.246 | (parietal — net-null*) | OFC, insula, sensorimotor |
| **θ** | absent | none | 2/10 / −0.040 / 0.722 | (MTL — net-null*) | — |

`*` = `loc_without_cohort_trace`: a system that clears the matched-strength
localization test on a band whose *net* cohort trace is null. Real spatial
structure, but not a trace (see §3).

## 2. What each tier means (and what is locked)

- **consistent** — a cohort-level trace exists AND concentrates in one system
  that survives matched-strength under BOTH epi modes AND shaft-collapse.
  - **β → OFC** (q=0.025 contact both-epi; q=0.025 shaft-collapse; **+ sensorimotor
    DEPLETED**, new lower-tail q=0.025). Fully locked — matches the standing
    OFC result (q=0.009 at R=1000, bilateral, LOO + shaft-collapse robust). The
    sensorimotor *depletion* is the new matched-strength-grade confirmation of
    the per-node "sensorimotor is the only net-anti system" finding.
  - **γ_l → PFC** (q=0.050 all four conditions; occipital depleted). New lock at
    matched-strength + shaft-collapse + both-epi. **LOO still owed** — until then
    γ_l-consistency is one rung below β.
- **patient-specific** — a trace exists for several patients but does NOT
  localize to any system surviving the filters.
  - **α** — clearest case: a real *cohort* trace (gate p=0.002) with **no**
    surviving carrier. The α trace is real but delocalised — each patient places
    it differently. (Matches the per-node "α flat/incoherent" finding.)
  - **δ** — subset-level trace (4/10), no robust carrier.
  - **γ_h** — subset-level (4/10) with NO cohort-median shift (median ρ=0.000);
    its parietal lean is real but rides a net-null band, so it is subset-specific,
    not consistent.
- **absent** — **θ**: no net cohort trace (2/10, median ρ<0, gate p=0.72). The
  one genuinely trace-free band at the cohort level.

## 3. The two honest refinements (do not bury)

1. **All carriers are shaft-robust.** The pre-registered worry was that θ→MTL /
   γ_h→parietal were single-shaft artifacts. They are NOT — every contact-level
   carrier survives shaft-collapse (β-OFC 0.025→0.025, γ_l-PFC 0.050→0.050,
   γ_h-parietal 0.050→0.050, θ-MTL 0.025→0.050, δ-cingulate 0.100→0.100). So the
   taxonomy split is **not** "which localizations are real" — they all are — but
   "which bands have a *net cohort trace* for the localization to be a trace OF."
2. **θ and γ_h localize without a net trace.** θ-MTL and γ_h-parietal clear the
   matched-strength *localization* test while the band's cohort ρ_split is ≈0.
   Mechanically: the per-pair concordance is spatially concentrated (some systems
   above-baseline, others below) but sums to ~zero net persistence. This is a
   genuine finding — **the multiscale localizer detects reproducible spatial
   structure even where the net trace is null** — and it is a point FOR the
   method, but it is NOT a cohort trace and must not be reported as one. θ stays
   "the absent band"; γ_h stays patient-specific.

The clean organising principle is two orthogonal axes: **(net cohort-trace
strength)** β > α≈γ_l > δ≈γ_h > θ — the detectability/SNR axis — and
**(consistent anatomical home)** β,γ_l,δ,γ_h,θ yes / α no. "Consistent" tier =
strong on BOTH; "absent" = null on the first.

## 4. Caveats

- **LOO owed for γ_l→PFC** (audit_83 has no LOO column; β→OFC already has it from
  the prior R=1000 lock). Run before publishing γ_l as hard-consistent.
- **Reconcile with the broadband/SNR table** (`per_patient_per_band.csv`, handoff
  §6): θ is 2/10 cohort and 2/5 even among high-SNR patients — the weakest band
  on BOTH the detectability axis and here. "Absent (net cohort trace)" and
  "sporadic at best individually" are the same band; no contradiction.
- **Thresholds are a-priori, not tuned to the result:** trace_level=cohort iff
  one-sided gate p<0.05 OR ≥5/10 patients clear own p95; subset iff 3–4/10; none
  iff ≤2/10. Carrier bar = the locked β→OFC bar (BH q<0.05 across a-priori
  systems within band, both epi, contact-level; q<0.10 under shaft-collapse).
  The match to {β,γ_l | α,δ,γ_h | θ} is genuine support, not circular fitting —
  the underlying graded numbers are in the table so borderlines (γ_l gate p=0.116;
  γ_h parietal-lean) are visible.

## 5. Method

- `audit_83 --all-bands` and `--all-bands --shaft-collapse` (both write epi
  include + exclude) produced the 4 condition CSVs in
  `data/audit/localization_atlas/matched_strength_allbands*`. The lower-tail
  matched-strength column (`matched_strength_p_lower`) was added additively — the
  existing upper-tail column is bit-identical, so the OFC q=0.009 lock is
  unchanged. Surrogate cache: the R=200/seed20260511 split ensemble was found
  ALREADY cached for all 6 bands (the handoff's "generate δ/θ/γ_h first"
  prerequisite was unnecessary — no cache generation was needed).
- `audit_83b_per_band_taxonomy.py` synthesises the 4 conditions × the audit_63
  cohort gate into the per-band verdict (5-point preamble in its docstring).

## 6. Cascade implications (for the headline rewrite)

- The band narration is RECOVERED as this consistency taxonomy — supporting
  result B in the spine (`2026-06-25_snr-band-taxonomy-handoff.md` §2). It is no
  longer a placeholder; the directive
  `writing_directive_2026-06-25_headline-restructure-broadband-snr.md` can be
  updated to the locked table.
- R1's band story: broadband+θ-exempt presence (detectability) is ORTHOGONAL to
  this consistency grading. β is special on BOTH (strong trace + OFC home) — the
  bridge to R2 (β = cognitive content).
- **Cascade into locked ledgers still ON HOLD** until the headline structure is
  agreed with the PI (Step 3). γ_l-LOO is the one open compute.
