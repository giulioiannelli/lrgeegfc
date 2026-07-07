---
name: per-node-trace-anatomy-and-heterogeneity
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-06-25
updated: 2026-06-25
pointers:
  - .agents/guides/task-persistence-investigation/2026-06-25_per-node-trace-decomposition.md
  - .agents/reports/2026-06-25_trace-heterogeneity-handoff.md
  - scripts/01_compute/audit/audit_144_per_node_trace_decomposition.py
  - data/audit/per_node_trace_decomposition/per_node.csv
  - .agents/reports/archive/2026-06/  # localization_audit_plan (audit_83 OFC lock)
---

> **Head.** The bottom-up per-node decomposition of `ρ_split` did exactly what the
> PI asked — found, per patient, which nodes carry vs oppose the trace — and it
> delivers **one real result and one honest non-result**. (1) **The β trace is
> anatomically specific**: it is carried by the limbic / cognitive-map core
> (**OFC** [locked], MTL, cingulate — all net-carrier) and is the **only** opposed
> by **sensorimotor cortex**, the single cortical system whose nodes lean
> *anti-trace* (carrier 0.10 / anti 0.20; multivariate within-patient p=0.006 after
> depth+strength+patient-FE). This deepens the flagship OFC result by naming its
> anatomical complement: the trace is a property of the cognitive-map system, not
> the whole brain. (2) **The per-patient heterogeneity is NOT explained by any node
> property** — not system coverage, tissue, epilepsy, hemisphere, depth, strength,
> baseline split-half stability, or data size. Every cross-patient correlation is
> weak (|r| ≤ 0.5, none significant at n=10) and several contradict. The non-tracers
> fail to consolidate the trace *even in their carrier-system tissue* (their anti
> nodes scatter across systems, including OFC/MTL), so the heterogeneity is a
> patient-**state** phenomenon, not a patient-**sampling** one.
>
> **UPDATE (the explanator, §4b): the heterogeneity is a measurement-DETECTABILITY
> axis, not an anatomical one.** The per-patient trace magnitude is gated by the
> **task-vs-baseline signal-to-noise ratio** `snr = d_task/d_noise` — how far the
> task reorganizes the hierarchy (`d_task = 1−Spearman(D_tt, D_pre)`) relative to the
> resting baseline's split-half noise floor (`d_noise = 1−Spearman(D_preA, D_preB)`).
> Cross-patient **ρ(ρ_split, snr) = +0.78 (β), +0.72 (γ_l), +0.70 pooled (n=30),
> LOO-robust [+0.70,+0.88]**. This is a *detectability/reliability* effect, not
> biology: the matched-strength surrogate trace is itself snr-ordered (+0.84) at ~10×
> smaller magnitude. Where snr is high the real trace is large and clearly above its
> surrogate (Pat_02/05/08 obs ≈0.5 vs surr ≈0.05); where snr ≈ 1 neither observed nor
> surrogate shows structure — **the trace is undetectable, not demonstrably absent**.
> So the heterogeneity is a detectability continuum consistent with a trace present
> across the cohort whose visibility is gated by recording/encoding SNR — which
> *defends* the cohort claim (apparent non-tracers are SNR-limited, not
> counter-evidence) and explains why no spatial factor separated the groups (the
> separating variable is patient/recording-level SNR, not anatomy). The publishable
> layer is the cohort trace + its anatomy; the heterogeneity is framed as
> SNR-gated detectability.

## 1. What was built (scope-first, library-first)

- **Scope report** (rule-mandated, written before code):
  `.agents/guides/task-persistence-investigation/2026-06-25_per-node-trace-decomposition.md`
  (11-section + 5-point critical preamble).
- **Library promotions** into `lrg_eegfc.utils.metrics.node_localization` (general
  names, no manuscript scope; audit_83's local copies now delegate — verified
  **bit-identical**, conservation residual rel 9e-17, so the locked β→OFC result is
  unchanged): `rank_concordance`, `canonical_cophenet`, `node_incidence_mean`.
- **Script** `scripts/01_compute/audit/audit_144_per_node_trace_decomposition.py`
  → `data/audit/per_node_trace_decomposition/{per_node,per_patient_summary,characterization_*}.csv`.

## 2. The decomposition is faithful (the mechanism)

Per-node `T_i = mean_j c_ij` (endpoint-incidence mean of the rank-concordance);
`mean_i T_i ∝ ρ_split` (conservation). Each node is labelled **carrier / neutral /
anti** against its **own matched-strength surrogate** (R=200, the audit_63/83
ensemble). Cross-patient, the carrier/anti balance *is* the trace:

| band | corr(ρ_split, carrier_frac) | corr(ρ_split, anti_frac) |
|---|---|---|
| β | **+0.88** | −0.79 |
| α | +0.90 | −0.85 |
| γ_l | +0.95 | −0.81 |

The non-tracers are non-tracers because they are anti-heavy: Pat_10 β (10 carrier /
22 anti, ρ=−0.09), Pat_14 β (12 / 18, ρ=−0.05), Pat_14 γ_l (14 / 64, ρ=−0.22).
Strong tracers are nearly all-carrier (Pat_08 β 66 / 0). This is the conservation
identity made visible — and the platform for asking *what the anti nodes share*.

## 3. RESULT — the β trace is anatomically specific (carrier core vs sensorimotor antagonist)

Per-system trace propensity (gray nodes, patient-demeaned continuous polarity
`π̃`, pooled cohort; carrier_frac − anti_frac in the last column):

| system | n | π̃ (demeaned) | carrier | anti | net |
|---|---|---|---|---|---|
| **OFC** | 58 | **+0.37** | 0.41 | 0.05 | +0.36 |
| **MTL** | 50 | **+0.28** | 0.58 | 0.06 | +0.52 |
| **cingulate** | 36 | +0.16 | 0.42 | 0.06 | +0.36 |
| lateral_temporal | 169 | −0.02 | 0.39 | 0.03 | +0.36 |
| occipital | 26 | −0.03 | 0.46 | 0.04 | +0.42 |
| insula | 22 | −0.03 | 0.32 | 0.14 | +0.18 |
| PFC | 175 | −0.16 | 0.31 | 0.10 | +0.21 |
| parietal | 42 | −0.24 | 0.29 | 0.07 | +0.21 |
| **sensorimotor** | 51 | **−0.28** | 0.10 | 0.20 | **−0.10** |

- The carrier top is **OFC / MTL / cingulate** — the limbic / cognitive-map core,
  recapitulating the locked audit_83 β→OFC localization (OFC matched-strength
  q=0.009–0.013) and adding MTL/cingulate as node-level carrier-leaning (consistent
  with audit_83's paralimbic hint; MTL was LOO-fragile at the cohort-median rung, so
  treat OFC as the locked carrier and MTL/cingulate as corroborating).
- **Sensorimotor is the only net-anti system** (carrier 0.10 < anti 0.20; most
  negative π̃). Control: multivariate logistic `anti ~ sensorimotor + depth +
  strength + C(patient)` gives sensorimotor coef +2.83, within-patient permutation
  **p=0.006** — it is anti BEYOND depth, strength, and per-patient base rate.
- **Honest bound (no overclaim):** this is a node-level *leaning*, NOT a causal
  suppressor. The closing test (exclude all sensorimotor, recompute ρ_split vs
  size-matched random-node decimation) gives median Δρ=+0.014, **0/6 patients beat
  random** — removing sensorimotor does not raise the cohort trace beyond node-count.
  So the defensible claim is *specificity* ("the trace lives in the cognitive-map
  core and not in primary sensorimotor cortex"), not *antagonism* ("sensorimotor
  suppresses the trace").
- γ_l gradient differs (MTL-led, occipital most-anti −0.66); α gradient is flat
  (MTL +0.07 top) — consistent with α being the pairwise-only secondary band. The
  band structure is preserved (constraint 2).

## 4. NON-RESULT — the heterogeneity is not node-property-explained (exhaustive)

Every candidate explanator of *who traces* was tested and failed. Cross-patient
(n=10) correlations of β ρ_split against coverage / quality factors:

| candidate | metric | corr with ρ_split | verdict |
|---|---|---|---|
| sensorimotor coverage | frac sensorimotor (gray) | −0.44 | n.s. (n=10 needs |r|>0.63) |
| carrier-core coverage | frac OFC+MTL+cingulate | **+0.07** | flat — coverage of carrier tissue does NOT predict tracing |
| epilepsy (node level) | epi→anti logistic, patient-FE | p=0.389 (β); α: epi never anti | null — epi nodes are NOT the anti nodes (epi is carrier-leaning if anything) |
| epilepsy (coverage) | frac epi vs ρ (α) | −0.76 | confounded coverage artifact — dies at the node level (row above) |
| hemisphere / implant side | frac R | −0.42 | null — Pat_14 is L-only and a non-tracer |
| along-shaft depth | anti−carrier depth gap | (controlled) | anti slightly deeper cohort-wide, but sensorimotor survives it; not the driver |
| baseline split-half stability | Spearman(D_preA,D_preB) | +0.42 (coph) / **−0.44 (raw FC)** | self-contradictory → no clean baseline-quality story |
| data size / hubness | N_nodes / strength CV | +0.28 / −0.49 | n.s. |

The decisive observation: the non-tracers' anti nodes are **scattered across systems
including the carrier core** (Pat_14 β anti = PFC 5, WM 5, MTL 3, lat-temp 3, SM 2;
Pat_10 β anti = PFC 7, WM 5, OFC 3, SM 3, …). In a non-tracer, normally-carrier
tissue (OFC, MTL) goes anti. So the failure is **not** "this patient sampled the
wrong tissue" — it is "in this patient the trace did not consolidate anywhere,
including the cognitive-map core." That is a patient-state / data-provenance fact
(note Pat_14's task_test was vendor-replaced; Pat_10 is R-only with dropped rows),
not an implant-anatomy fact, and at n=10 it is not separable from noise.

## 4b. THE EXPLANATOR — task-vs-baseline SNR (detectability), not anatomy

The variability is not spatial; it is a **measurement-detectability** axis. Define,
per patient and band, on the canonical cophenetic distances:

- `d_task  = 1 − Spearman(D_tt, D_pre)`   — how far the task moved the hierarchy (encoding magnitude)
- `d_noise = 1 − Spearman(D_preA, D_preB)` — the resting baseline split-half noise floor
- `snr     = d_task / d_noise`             — task signal above the baseline noise

| relation (cross-patient, n=10; pooled n=30) | value |
|---|---|
| ρ(ρ_split, **snr**) — β / γ_l / α | **+0.78** / +0.72 / +0.33 |
| ρ(ρ_split, snr) — pooled band-demeaned (n=30) | **+0.70** |
| β LOO range (drop each patient) | [+0.70, +0.88] |
| β ρ(ρ_split, d_task \| d_noise) — encoding, noise held | +0.61 |
| β ρ(ρ_split, d_noise \| d_task) — noise floor, signal held | −0.59 |

Non-tracers sit at snr ≈ 1.0–1.3 (task barely clears the baseline noise): Pat_15
encodes weakly (`d_task` 0.16), Pat_10 reorganizes strongly but on a noisy baseline
(`d_task` 0.71, `d_noise` 0.53) — different failures, both captured by snr ≈ 1.
Strong tracers sit at snr 1.7–2.8.

**Honesty control (matched-strength surrogate):** snr also rank-orders the
*surrogate* ρ (+0.84), but at ~10× smaller magnitude — surrogate ρ ∈ [−0.004,
+0.091] vs observed [−0.09, +0.51]. So this is fundamentally a **detectability /
reliability** axis, NOT a biological one: where snr is high the real trace is large
and clearly separated from its surrogate (Pat_02/05/08 obs ≈0.5 vs surr ≈0.05); where
snr ≈ 1 neither observed nor surrogate shows structure (**undetectable, not
demonstrably absent**). The substantive (non-attenuation) half is the encoding term
`d_task` (+0.61 partial, computed against the *full* clean rest_pre, not the noisy
halves): patients whose task reorganized the hierarchy more leave a stronger trace —
the cross-patient echo of the N2 encoding result.

**Why this is the right answer to "what explains the variability":** the separating
variable is **patient/recording-level SNR**, which is exactly why no *spatial*
(node/system/tissue/epi/hemisphere) factor separated the groups in §4 — they were the
wrong axis. It reframes the heterogeneity as a **detectability continuum** consistent
with a trace present across the cohort whose visibility is gated by how cleanly each
patient's task reorganization rises above baseline noise. This **defends** the
cohort-level claim: apparent non-tracers are SNR-limited, not counter-evidence.
(Source: `data/audit/per_node_trace_decomposition/` + inline analysis 2026-06-25.)

## 5. Implication for the manuscript (consistent with the PI's constraints)

- **Constraint 1 satisfied honestly:** "it's the patient" is not a result, and we
  now know *why* — the heterogeneity is **SNR-gated detectability** (§4b), not
  biology or implant. This is publishable framing: present the trace as cohort-level
  with **detectability gated by task-vs-baseline SNR**, so apparent non-tracers are
  signal-limited, not counter-evidence. The consistency among the patients where the
  trace is *measurable* IS the result; the spread is detectability noise.
- **Constraint 5 served:** the new anatomy (carrier core + sensorimotor specificity)
  strengthens the OLD flagship (OFC cognitive-map trace) by showing the trace is
  cognitive-system-specific, not brain-wide.
- **Constraint 3 respected:** this is the ρ^coph probe only; the Grassmann per-node
  analog (`grassmann_trace_contributions`) is untouched and remains the distinct
  subspace story.
- **The gate report's §6** (`2026-06-25_cophenetic-gate-presence-vs-consistency.md`,
  Pat_06 coverage) is now fully superseded: coverage explains neither Pat_06 nor the
  cohort heterogeneity.

## 6. Status + next steps

The heterogeneity question is **answered**: it is SNR-gated detectability (§4b), not
a node/anatomical/implant property (§4). Two deliverables remain to make this
publication-grade:

1. **Detectability figure** (the explanator): per-band scatter `ρ_split` vs `snr`
   (patients labelled), beside the obs-vs-surrogate separation showing high-snr
   patients clear their surrogate and low-snr patients don't. This is the figure that
   defends the cohort claim against "half your patients don't show it".
2. **Sensorimotor-specificity lock**: run the anatomical-specificity through the
   audit_83 machinery as a **lower-tail** matched-strength test (M_obs < surrogate,
   per system) with shaft-collapse + LOO — directly comparable to the OFC upper-tail
   lock — so "trace lives in the cognitive-map core, absent in sensorimotor" is a
   matched-strength result, not just a within-cohort gradient.

Optional / lower priority: characterize `d_noise` (which patients have noisy
baselines and why — recording length / artifact load) to make the SNR axis concrete,
but this is a data-audit thread, not load-bearing for the framing.
