---
name: per-band-phenomenology-vision
type: report
era: "IMCOH_ABS × COHORT_N10"
status: current
created: 2026-06-26
supersedes_claims:
  - "2026-06-25_per-band-consistency-taxonomy-lock.md: γ_l→PFC 'consistent' (now LOO-FAILED → strong-subset)"
  - "per_node_trace_anatomy_2026_06_25 memory: 'heterogeneity = detectability/SNR' (now: residual biology, reliability explains ~0% in β)"
pointers:
  - data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
  - data/audit/localization_atlas/per_band_taxonomy_verdict.csv
  - data/audit/localization_atlas/carrier_loo.csv
  - data/audit/multiphase_snr/per_band_summary.csv
  - data/audit/ofc_coverage/ofc_coverage_vs_trace.csv
  - .agents/reports/2026-06-25_per-band-consistency-taxonomy-lock.md
  - .agents/reports/2026-06-25_multiphase-snr-reliability.md
  - .agents/reports/2026-06-25_raw-vs-multiscale-trace.md
---

> **HEAD.** Only **β** is a fully-locked cohort result: 7/10 patients trace, the
> group gate clears at p=0.005, it concentrates in **OFC** (LOO-robust, every
> drop p=0.005) and is actively *depleted* in sensorimotor/PFC. **α** is also
> cohort-significant (p=0.002) but has **no anatomical home** — each patient
> places it differently. **Low-γ** is a strong 5-patient subset that **fails the
> cohort gate** (p=0.12) and whose PFC home **fails leave-one-out** — not locked.
> **δ** and **high-γ** are patient-specific subsets (4/10, no net cohort trace).
> **θ is absent** (2/10, median ρ<0). The apparent "broadband" trace is a
> per-PATIENT property of ~5 strong tracers (Pat_06 traces all 6 bands, Pat_15
> traces none), **not** a cohort property — and that per-patient spread is **real
> biology**: it is NOT explained by measurement reliability (β: ~0%) nor by OFC
> electrode coverage (refuted, ρ=−0.10).

---

> **⚙ ESTIMATOR UPDATE (2026-07-06) — this report predates the ρ_sym migration; read
> its numbers with this note.** The cophenetic trace estimator is now **ρ_sym** =
> ½[ρ(D_task−D_preA, D_post−D_preB) + ρ(D_task−D_preB, D_post−D_preA)], the mean of both
> split-half arm assignments — it removes the arbitrary-half artifact of bare ρ_split
> (which sign-flipped 20/60 near-zero cells; `audit_149`,
> `feedback_rho_sym_canonical_estimator`). **Every verdict in this report holds**; the gate
> p-values become **β p=0.032, α p=0.024 (both CLEAR)**; δ/θ/low-γ/high-γ still fail
> (`audit_150`, `data/audit/rho_sym_gate/`). These are larger than the split values in the
> HEAD (β 0.005, α 0.002) at the **same** R=200 — not a grid effect but the estimator
> correctly pulling ill-conditioned near-zero patients toward zero; strong tracers unchanged,
> zero verdict flips (`audit_149`). β→OFC (q=0.025, sensorimotor/PFC depleted; `audit_151`),
> the β-only inference-specific arc (T_infspec_pe p=0.0098; `audit_152`), the anchor/trace
> composition (`audit_153`), and the carrier/anti per-node split (node-level ρ=0.96;
> `audit_154`) all reproduce. **New per-patient reporting tier:** report ρ_sym alongside
> its split-uncertainty ½|ρ_AB−ρ_BA|, and mark **|ρ_sym| < 1 SE "undetermined"** — the
> near-zero third of the cohort (e.g. Pat_13/Pat_15 β) whose sign is estimator noise, not
> biology. The "broadband is per-patient, real biology, not detectability" framing (HEAD)
> is unchanged and now also carries the estimator-invariance guarantee. Synthesis:
> `.agents/reports/2026-07-06_rho-sym-pipeline-migration.md`.
>
> **📖 READING KEY — never report a bare ρ; pair it with three confound-free companions**
> (`× null` effect size, `% held` = of movers, persistence-vs-reset, `% struct` = trace share
> of total cross-phase energy): **β** ρ=0.20 · **16.6× null** · **54% held** · **26% struct**
> (structurally substantial → localizes to OFC); **α** ρ=0.10 · 14× null · 41% held · **3%
> struct** (statistically real but a thin whisper → no anatomical home). Guards: high ×null
> alone ≠ trace (low-γ 48× but fails the gate), and % held is only meaningful once the gate
> clears (δ 61% = noise). Full key + how-to-say-it: `2026-07-06_rho-sym-panoramic-and-methodology.md`
> (Reading key) + `2026-07-06_critical-appraisal-of-the-measure.md`.

---

## 0. The one master table

Per-patient `ρ_split` (cophenetic trace; >0 = trace), individual matched-strength
clearance (✓ = own p<0.05, **a** = significantly *anti*, · = null), and the
cohort verdict. Source: `per_patient_per_band.csv` + `cohort_summary.csv`.

| Patient | δ | θ | α | β | γ_l | γ_h | #bands traced |
|---|---|---|---|---|---|---|---|
| Pat_02 | ✓ +.42 | · −.00 | · −.04 | ✓ +.51 | ✓ +.76 | a −.17 | **3** |
| Pat_03 | ✓ +.31 | a −.17 | ✓ +.36 | ✓ +.37 | ✓ +.16 | ✓ +.51 | **5** |
| Pat_05 | · +.04 | ✓ +.09 | · +.08 | ✓ +.49 | ✓ +.85 | ✓ +.92 | **4** |
| Pat_06 | ✓ +.79 | ✓ +.78 | ✓ +.83 | ✓ +.21 | ✓ +.59 | ✓ +.52 | **6** |
| Pat_07 | ✓ +.16 | · +.09 | · +.08 | ✓ +.23 | · −.06 | a −.15 | **2** |
| Pat_08 | a −.12 | a −.13 | ✓ +.38 | ✓ +.50 | ✓ +.37 | ✓ +.31 | **4** |
| Pat_10 | · −.02 | · −.08 | ✓ +.13 | a −.09 | · +.00 | · +.01 | **1** |
| Pat_13 | a −.15 | · −.08 | · +.03 | ✓ +.21 | · −.04 | · −.01 | **1** |
| Pat_14 | · −.07 | · +.04 | ✓ +.21 | · −.05 | a −.22 | a −.14 | **1** |
| Pat_15 | · −.05 | a −.10 | · +.04 | · +.08 | · −.13 | a −.47 | **0** |
| **#trace /10** | **4** | **2** | **5** | **7** | **5** | **4** | |
| **median ρ** | +.008 | −.040 | +.105 | **+.221** | +.083 | +.000 | |
| **gate Wilcoxon p** | .278 | .722 | **.002** | **.005** | .116 | .246 | |
| **cohort gate** | ✗ | ✗ | **✓** | **✓** | ✗ | ✗ | |
| **robust home** | — | (MTL*) | — | **OFC** | PFC† | (parietal*) | |
| **tier** | patient-spec | **absent** | cohort/no-home | **consistent** | strong-subset | patient-spec | |

`*` localizes on a band with ~zero net cohort trace ("localization without a
trace", §4). `†` clears full-cohort but **fails LOO** (§3).

---

## 1. The cohort axis (group evidence)

Two and only two bands pass the cohort gate (paired Wilcoxon, obs vs
surrogate-median, one-sided):

- **β: p=0.005**, median ρ=+0.221, 7/10 trace. Largest effect size of any band.
- **α: p=0.002**, median ρ=+0.105, 5/10 trace. *Smallest* gate p — but a smaller,
  more uniform shift (no patient goes strongly anti; the strongest negative is
  Pat_02 at −0.04).

Everything else fails the gate: low-γ p=0.12, high-γ p=0.25, δ p=0.28, θ p=0.72.
The β/α split from the rest is genuine and quantitative — not a threshold
artifact (the underlying medians are +0.22, +0.11 vs +0.08, +0.01, +0.00, −0.04).

**Why low-γ/δ/high-γ fail despite 4–5 clearers:** strong anti-tracers sink the
paired test. γ_l has Pat_14 at −0.22; γ_h has Pat_15 at −0.47 and Pat_02 at −0.17.
β, by contrast, has only mild negatives (Pat_10 −0.09, Pat_14 −0.05) against six
large positives — so it survives the group test where the others don't.

## 2. The per-patient axis (variability IS the signal)

There is a clean, monotone **"number of bands traced" gradient** from 6 down to 0:

> Pat_06 (6) → Pat_03 (5) → Pat_05, Pat_08 (4) → Pat_02 (3) → Pat_07 (2) →
> Pat_10, Pat_13, Pat_14 (1) → Pat_15 (0)

- **Pat_06 is the universal tracer** — clears all six bands, with huge ρ in the
  low bands (δ +0.79, θ +0.78, α +0.83).
- **Pat_15 is the universal non-tracer** — zero bands, strongly anti in high-γ
  (−0.47). (Consistent with its right-hemisphere-only implant, the standing
  single-patient β exception.)
- **β is the great equalizer:** it is the *only* band traced by patients who
  trace nothing else (Pat_13: β +0.21 and nothing; Pat_07: only δ and β). That is
  why β is the cohort flagship — its support is broad *and* reaches the weak
  patients.

**The "broadband" appearance is per-patient, not cohort.** The five patients who
light up many bands (Pat_06/03/05/08/02) are exactly the five "high-SNR" patients.
So the trace looks broadband because a handful of strong tracers carry it across
bands — at the cohort level only α and β survive. *Recommended framing:* "β and α
are cohort-consistent; γ_l/δ/γ_h appear in a strong-tracer subset; θ is absent" —
**not** "the trace is broadband." (This supersedes the earlier "broadband,
θ-exempt LOCKED" framing.)

## 3. The localization axis (anatomical home)

Per-system matched-strength (audit_83, BH across a-priori systems, both epi modes,
contact + shaft-collapse) + LOO (audit_83c):

- **β → OFC. LOCKED.** q=0.025 (contact, both epi), shaft-robust, and
  **LOO-robust: every leave-one-out drop stays p=0.005, q≤0.05** (`carrier_loo.csv`).
  **Sensorimotor and PFC are significantly *depleted*** (new lower-tail q=0.025) —
  the trace doesn't just concentrate in OFC, it actively avoids sensorimotor.
  This matches the standing R=1000 OFC lock (q=0.009, bilateral).
- **low-γ → PFC. NOT robust.** Clears the full cohort (q=0.025–0.05) but **fails
  LOO**: dropping any of Pat_05/06/07/10/14/15 pushes it to q>0.05 (worst:
  −Pat_14, p=0.035, q=0.17). Six of ten single-patient drops break it ⇒ PFC is
  **not** a cohort-stable home. Downgrade γ_l from "consistent" to **strong
  subset with a single-patient-fragile localization**. (This closes the one open
  item from the 2026-06-25 taxonomy lock.)
- **α → no home.** Real cohort trace (gate p=0.002), but **no system survives**
  the localizer under any condition. The α trace is genuine and delocalised — each
  patient places it somewhere different.
- **δ → no robust home** (cingulate is borderline, exclude-only).
- **θ → MTL/lateral-temporal**, **γ_h → parietal**: both localize (shaft-robust)
  **but on bands with ~zero net cohort trace** — see §4.

So **β is the only band that is strong on BOTH axes** (cohort trace + stable
anatomical home). That single fact is the spine of the whole result.

## 4. "Localization without a trace" (a method point, not a result)

θ-MTL and γ_h-parietal clear the per-system matched-strength *localization* test
while their band's net cohort ρ_split ≈ 0 (θ median −0.04, γ_h median +0.00). The
per-pair concordance is spatially concentrated (some systems above baseline,
others below) but **sums to zero net persistence**. This shows the multiscale
localizer detects reproducible spatial structure even where the net trace is null
— a point *for* the method's sensitivity — but it is **not a cohort trace** and
must never be reported as one. θ stays "the absent band"; γ_h stays patient-specific.

## 5. Is the per-patient spread just measurement noise? (Q1 — NO)

The decisive question for the cohort claim: is "who traces" an artifact of how
cleanly each patient was measured? **No.** (`multiphase_snr/per_band_summary.csv`,
`ofc_coverage_vs_trace.csv`.)

- **Split-half reliability explains almost none of it.** Spearman(ρ_split,
  within-phase reliability) is near zero in every band and essentially **zero for
  β (ρ=0.055, R²=0.0006)**. Across bands reliability accounts for <20% of the
  spread (δ 12%, θ 19%, α 5%, β 0.06%, γ_l 4%, γ_h 3%).
- **The killer counterexample:** Pat_15 has the **highest** β split-half
  reliability (0.842) yet the **lowest** β trace (+0.083). Most cleanly measured,
  least trace. The spread is not noise.
- **Coverage is refuted too.** Spearman(β ρ_split, OFC-electrode fraction) =
  **−0.10** (flat/negative). **Pat_08** is the 2nd-strongest β tracer (+0.50) with
  **zero OFC electrodes**; Pat_03 (+0.37) and Pat_13 (+0.21) also have none. A
  patient does **not** need to implant the OFC hub to show the β trace —
  consistent with "OFC = hotspot, not container."

**Caveat on the SNR gradient.** The band-count gradient (§2) *does* track the
task-vs-baseline SNR ranking — but that SNR is partly circular (a strong
reorganization inflates both the task contrast *and* ρ_split). The clean
measurement-noise axis (split-half) is flat. So the honest reading: the spread is
**residual biology, real and currently unexplained** — not detectability, not
coverage. It does **not** undermine the cohort β claim, which stands on the group
gate + matched-strength + LOO independently of why individuals differ.

## 6. Multiscale vs raw (Q2 — the trace only becomes interpretable under the filter)

The raw per-edge trace is *larger* in magnitude and present *everywhere* — including
θ, where there is no real structure. That ubiquity is precisely its **triviality**:
raw can't separate signal from noise. Passed through the multiscale (cophenetic)
filter, the trivial blob is **refined** — θ collapses to ~0, β sharpens onto OFC,
sensorimotor empties out. The interpretable anatomy (§3) is a property of the
multiscale measure, not of the raw edges. (Report:
`2026-06-25_raw-vs-multiscale-trace.md`; head-to-head with spectral clustering:
`audit_141_spectral_clustering_headtohead.py`.)

## 7. The current vision in one paragraph

A task (transitive inference) leaves an **offline cophenetic trace** in rest_post.
At the cohort level it is carried by **β** (the flagship: 7/10, p=0.005, concentrates
in **OFC**, avoids sensorimotor, LOO-robust) and, more weakly and without an
anatomical home, by **α** (5/10, p=0.002, delocalised). **Low-γ** rides a strong
5-patient subset but is not cohort-locked and its PFC home is single-patient-fragile;
**δ** and **high-γ** are patient-specific subsets; **θ** is absent. Underneath the
cohort verdict is a real, biological **per-patient gradient** — some patients trace
across many bands, some not at all — that is **not** measurement noise and **not**
implant coverage, and remains unexplained. β is the one band strong on every axis,
which is why it anchors the manuscript and why the encoding/inference dissociation
(N2) is built on it.
