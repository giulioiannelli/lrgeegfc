---
name: preprint-anatomy-ledger
era: IMCOH_ABS_COHORT_N10
status: locked_2026-05-19
kind: anatomy-verdict-lockdown
supersedes: kc_era_anatomy_artifacts
companion: ANATOMY_CONTROLS.md
parent: VERDICT_LEDGER.md
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
fc_method: imcoh_abs
tau: 1/lambda_max
---

# Anatomy verdict ledger — per-(band, probe) trace localization (locked 2026-05-19)

> ## ✅ FINAL AUDITED VERDICT 2026-06-10 — β trace is brain-wide but CONCENTRATES in ORBITOFRONTAL CORTEX (matched-strength R=1000, BH q≈0.01, shaft- & LOO-robust). Concentration above a distributed baseline, not exclusive presence; bilateral, peaks at the *system* scale. Supersedes the 2026-06-05 R=200 banner below; corrects the α claim.
>
> The 2026-06-05 exploration was carried through a **six-rung control gate + a
> self-audit that found and fixed two real flaws**. Final state
> (`data/audit/localization_atlas/README.md`, 2026-06-10):
>
> - **β trace concentrates in orbitofrontal cortex (OFC system)** above its
>   brain-wide baseline (the trace exists cohort-wide; this measures above-baseline
>   *concentration*, per-patient demeaned — even within OFC the lean is ≈0.59,
>   ~41% anti). Matched-strength surrogate
>   (**R=1000**, fresh independent ensemble = seed-robustness check), rank
>   concordance, **canonical observed trace** (recomputed to match the surrogate
>   construction — fix below): **BH q = 0.009–0.013 across the 9 a-priori
>   anatomical systems within β, in ALL FOUR conditions** (contact/shaft-collapsed
>   × epi-include/exclude). 4/5 implanted patients positive, **leave-one-out
>   robust** (LOO-min median ≫ 0, survives dropping the −1.5M anti patient Pat_10),
>   **shaft-collapse robust** (multi-shaft in the two strongest patients), in
>   **low-strength / non-hub** tissue (strength dev −0.51). The concentration
>   **peaks at the *system* scale** — it is **bilateral** (no L/R lateralization,
>   p≈0.48/0.76) and **washes out** when pooled to the whole frontal lobe (p≈0.29)
>   or split by hemisphere; at single-region granularity it is undersampled (fails
>   BH over 53 regions). This is the **only** band/region surviving the full gate.
> - **Hippocampus / MTL — sub-threshold** (matched-strength p ≈ 0.4 at R=1000,
>   leave-one-out fragile, Pat_02-driven). Real, strength-independent, but
>   underpowered at n=5 — a hint, NOT the localization. (Supersedes the
>   2026-06-05 Hip-β reinstatement banner: the matched-strength-significant β
>   localization is **orbitofrontal**, not hippocampal.)
> - **"Distributed paralimbic ring" — RETRACTED.** Ring-minus-OFC clears at
>   contact level but **fails shaft-collapse** (mode-unstable, q=0.07–0.32); the
>   pooled "limbic" supersystem clears only because OFC is in it (and itself
>   fails shaft-collapse-include). The localization is **OFC-specific**, not a ring.
> - **α — NO FDR-surviving localization** (CORRECTS the 2026-06-05 "α → frontal
>   operculum" line). pars opercularis raw p=0.005 but the region family is 53 DK
>   regions → BH q=0.066–0.26; no a-priori system clears. β localizes at the coarse
>   *system* scale (family 9, survives); α only at a fine single region (family 53,
>   fails).
> - **Occipital** clears every control incl. shaft + R=1000 (q=0.006–0.045) but on
>   only **K=3** patients → genuine low-coverage **secondary**, not headline.
> - **low-γ** epi-mode-unstable; **Grassmann subspace** trace does not localize
>   (single-probe, edge-cophenetic localization).
>
> **Self-audit fixes (2026-06-08/10):** (1) the observed trace was being read from
> the LRG-pipeline ultrametric while the surrogate used the canonical
> `cophenetic_condensed_from_eigs`; they differ by a per-phase ρ-trace scale that
> perturbs cross-phase ranks (Spearman 0.92–0.99). `audit_83.obs_trace` now
> recomputes the observed trace canonically (matches the surrogate AND the locked
> §5.3 audit_63); OFC verdict unchanged, null now exact. (2) α demoted (above).
> Multiplicity family is **a-priori systems within β only** — NOT bands (β is *the*
> band; α/γ_l are separate questions), NOT nested granularities, NOT pooled epi modes.
>
> **Status:** the β-OFC verdict is **audited and unattackable on method**; the only
> residual is the irreducible n=5 OFC coverage. Locked-verdicts table row for β is
> updated below. Source: `data/audit/localization_atlas/README.md`. Cascade into
> `00_cohort.md` / `01_beta.md` landed 2026-06-10; **commit pending user green-light**.
>
> ---

> ## ⚠️ VERDICT CORRECTION 2026-06-05 — the 2026-05-30 retraction over-reached on "no anatomy anywhere"; Hip-β (cophenet TRACE) is a real, modest, sampling-conditioned localization
>
> Re-examination (user-prompted, 2026-06-05) keeps the **multi-region-list** retractions
> below but **withdraws the blanket "spatially delocalized / no anatomy at cohort OR
> single-patient level"** claim. Two errors in the 2026-05-30/06-01 reading:
> 1. **Denominator.** "3–5/10 patients = marginal" uses the full cohort. The localization
>    denominator is patients **implanted in the region**. The 5/10 ceiling is implant
>    coverage, not a per-region weakness — so the ">5/10 or thin" bar retracts every region
>    by construction. A cohort permutation gate on n=10 / ≤5-per-region is a structural
>    false-negative machine (the audit's own §"methodological trap" shows both natural
>    statistics failing oppositely). Non-significance there ≠ absence.
> 2. **Anchor/trace conflation.** The "Hip = marginal hint" line leaned on the β **Grassmann**
>    Hip (phase-averaged ANCHOR). The genuine signal is the β **cophenet TRACE**.
>
> **Corrected Hip-β cophenet (trace):** of the 5 patients implanted in hippocampus, **4/5**
> show the β trace concentrating in Hip above their own demeaned baseline (Pat_02 +0.046,
> Pat_08 +0.016, Pat_03 +0.007, Pat_13 +0.006; Pat_14 −0.003 only). Cohort demeaned median
> **+0.0073 > shuffle q95 +0.0042**, **perm_p = 0.0046** (best-sampled passer of 53 regions;
> exact one-sided Wilcoxon p = 0.0625 borderline; robust to dropping Pat_02). A **real,
> modest hippocampal concentration of the β trace**, reported descriptively + sampling-
> conditioned. Caveats stand: half the cohort isn't implanted there; no single patient has
> Hip as #1 hotspot; prespecified-hippocampal framing (not 53-way fishing) is what carries it.
> **Status:** multi-region lists stay retracted; "delocalized" framing withdrawn; a broader
> sampling-conditioned + system-level (MTL/limbic) + epi-X-resolved localization exploration
> is **underway (2026-06-05)** before any re-lock. Source: `data/audit/anatomy_localization_wilcoxon/README.md` (2026-06-05 banner). See Revision history.
>
> **⚠️ FINAL UPDATE 2026-06-05 (verified — 3 nested controls + adversarial review):** the
> exploration ran (`audit_81` atlas + `audit_83` matched-strength + `audit_84` strength-
> residualization). Matched-strength verdict for the **β** trace (rank-concordance,
> strength-preserving R=200): **localizes to the ORBITOFRONTAL SYSTEM (medial+lateral OFC
> pooled), MS_p=0.005 include & exclude-epi, leave-one-out robust, 4/5 patients** + lateral
> temporal (superiortemporal, bankssts) at fine granularity. **α → frontal operculum**
> (parsopercularis MS_p=0.005). **HIPPOCAMPUS/MTL is a genuine but SUB-THRESHOLD concentration**
> (MS_p≈0.10): strength-INDEPENDENT (survives residualization; within-patient strength↔trace
> ρ≈0.14; MTL is low-strength) — a real moderate-effect hint underpowered at n=5, **NOT a
> strength artifact and NOT the headline**. Grassmann subspace + low-γ MTL do NOT localize
> under matched strength. **This supersedes the Hip-β reinstatement above:** the matched-
> strength-significant β localization is **orbitofrontal**, not hippocampal; Hip is a
> suggestive hint. Multi-region lists stay retracted; "delocalized" withdrawn (OFC clears).
> Adversarial review (agent a3dc4b5d) confirmed code clean + OFC LOO-robust + Hip-not-artifact.
> Source: `data/audit/localization_atlas/README.md`. **Cascade into verdict table / 00_cohort /
> 01_beta + commit: pending user green-light (finding flipped twice; PI should see final verdict).**
>
> ---

> ## ⚠️ RETRACTION 2026-05-30 — the "strong localized" anatomy verdicts do NOT survive a signed, sampling-aware localization test
>
> A signed, threshold-free localization audit (`diag_anatomy_localization_wilcoxon.py`;
> report `data/audit/anatomy_localization_wilcoxon/README.md`) replaced audit_71/72's
> direction-blind `|Δ|` + top-decile + endpoint-count enrichment with the signed score
> `s_ij = dD_task·dD_rest` (cophenet) / per-node participation deviation (Grassmann),
> per-patient demeaned (removes base rate), tested per region by a **floor-free
> permutation rank** of the across-patient median vs a within-patient region-label
> shuffle null. **Decisive finding: no DK region is sampled by more than 5/10 patients,
> and the locked region-lists rest on 1–4 patients each** (β `ctx-rh-insula` = Pat_10
> alone — the β-*anti* patient; several well-sampled locked regions are *anti*-localized).
> Requiring a region to be adequately sampled (n≥5, the only regions where a cohort claim
> is testable) AND localized AND in the locked list gives:
>
> | cell | locked region passing on n≥5 patients | verdict |
> |---|---|---|
> | β cophenet | 0/7 | **RETRACT** the DK list |
> | α cophenet | 0/11 | **RETRACT** (diffuse) |
> | β Grassmann | only Hip (n=5) crosses the test — marginal, not a localization | **RETRACT** the DK list |
> | γ_l Grassmann | 0/7 | **RETRACT** the DK list |
> | δ Grassmann (full) | 0/4 | **RETRACT** (well-sampled regions *anti*-localized) |
> | δ Grassmann (epi-X) | 0/3 | **RETRACT** (no locked region localized; mask bug also fixed) |
>
> **No region in any band or probe reaches a defensible cohort anatomical localization.**
> The Hippocampus (β) is the single region that crosses the permutation threshold, but on
> only **3–5/10 patients** (3 individually significant at cophenet; 5 sampling it at
> Grassmann) — the *same* thin-sampling zone the other regions were retracted for. By the
> consistent standard it is a **marginal hint, NOT an established localization**, and
> carries no claim. (It is also an *anchor*-flavored phase-averaged Grassmann quantity, not
> a cross-phase trace.) The per-(band, probe) "strong localized" subsections below are
> **superseded** and retained only as historical record. See the 2026-05-30 revision entry.
> (β/α cophenet were user-approved; β Grassmann + γ_l/δ retractions extend the identical
> finding from the same run — the earlier "soften to Hippocampus" framing was itself an
> inconsistency and is withdrawn.)
>
> ### ⚠️ Completing test 2026-06-01 — PER-PATIENT localization is ALSO null (trace is delocalized)
>
> The natural fallback — "maybe the trace localizes *within* each patient, just to a
> different region per patient (heterogeneous implants)" — was tested directly and **also
> fails**. For the cophenet trace, **0–1 of 10 patients** have a DK region beating their own
> within-patient implant-label-shuffle null in **every band** (Binomial p ≥ 0.40), and this
> holds at the higher-powered **lobe** granularity too (β 0/10, α 1/10). Every patient's
> strongest region *is* different (distinct argmax 9–10/10), but that is the **signature of
> no localization** (noise scattering the argmax over heterogeneous implants), not
> idiosyncratic localization. A between-region η² "concentration" that an adversarial check
> raised (5–8/10, α/β/γ) is **electrode-shaft spatial autocorrelation, not anatomy** — the
> anatomy-free shaft partition reproduces it equally (region↔shaft NMI = 0.65). The MAX test
> is calibrated *and* powered (it detects a planted hotspot at ≈0.9 SD; observed hotspots are
> below each patient's own detection floor). **Verdict: the verified β/α trace is spatially
> DELOCALIZED — no anatomical anchor at the cohort OR the single-patient level.** Report:
> `data/audit/per_patient_localization/README.md`; verified workflow `wf_ddbfbe43-da7`.

**Head.** This ledger is the **single source of truth** for where each trace lives anatomically, under the locked anatomy control battery (`ANATOMY_CONTROLS.md`, A1 hypergeometric + A3 matched-strength surrogate; A2 sampling-corrected bootstrap and A4 implant-geometry regression deferred to sensitivity supplement). All audits are now run; verdicts below are **locked**. KC-era anatomy artifacts (`lrg_localization_anatomy/`, the "Hippocampus + left fusiform" memory claim) are **retired** under the 2026-05-18 trace-side lockdown and **not citable** for preprint claims.

## Locked verdicts

| Band | Trace verdict (CONTROLS) | Probe(s) audited | Anatomy verdict (ANATOMY_CONTROLS) | Localization audit 2026-05-30 |
|---|---|---|---|---|
| **β** | strong trace, both probes | cophenet + Grassmann | ~~strong localized, both probes (7+7 named DK regions)~~ | **REINSTATED 2026-06-10 → CONCENTRATES in ORBITOFRONTAL CORTEX (cophenet).** Multi-region lists stay retracted; the brain-wide β cophenet trace **concentrates** (per-patient demeaned, above baseline — not exclusive) in the **OFC system**: matched-strength R=1000, BH q=0.009–0.013 (9 a-priori systems within β, all 4 contact/shaft × incl/excl), 4/5 implanted, LOO- & shaft-robust, low-strength, **bilateral** (peaks at *system* scale; washes out at lobe/hemisphere). Hip/MTL = sub-threshold hint (p≈0.4); occipital = K=3 secondary; paralimbic-ring & Grassmann = no concentration. See `localization_atlas/README.md`. |
| α | strong trace, only D_coph | cophenet (full + C5 epi-X) | ~~strong localized, only D_coph (11 named regions)~~ | **RETRACTED — no FDR-surviving localization** (2026-06-10: pars opercularis fails BH over 53 regions q=0.07–0.26; no system clears. Earlier "diffuse; 0/11" and the 2026-06-05 "frontal operculum" line both withdrawn.) |
| γ_l | strong trace, only Grassmann ↑ | Grassmann | ~~strong localized, only Grassmann (7 named DK regions)~~ | **RETRACTED (0/7 cohort-supported)** |
| δ | weak trace, only Grassmann (LOO Pat_08 fails Decision-12 precondition) | Grassmann (full + C5 epi-X) | ~~strong localized, only Grassmann (4 + 3 named regions)~~ | **RETRACTED (0/4 full, 0/3 epi-X; well-sampled regions anti-localized; epi-X mask bug fixed)** |
| θ | no trace | — | n/a | — |
| γ_h | no trace | — | n/a | — |

**Aggregate reading (REVISED 2026-05-30)**: the earlier claim that the trace
"is localized to band-specific distributed cortical networks" is **not supported**
once direction-blindness and the top-decile threshold are removed and per-region
patient sampling is accounted for. The locked DK region-lists are dominated by
1–4 patient findings (no region exceeds 5/10 patients). **No region — including
the Hippocampus — reaches a defensible cohort anatomical localization**: Hip is
the single region that crosses the statistical threshold, but on only 3–5/10
patients, the same thin-sampling zone every retracted region sits in, so it is a
marginal hint and not a claim. The anatomy is, at this implant coverage, **not
cohort-resolvable** — neither the clean distributed-network story previously
locked nor diffuse-brain-wide; it is *undersampled*. Per-patient, the trace lands
on largely patient-specific regions (α: zero cross-patient overlap; β: near-zero,
only Hip recurring in 3 patients), and that per-patient localization is itself
weak (~1.2–1.3 significant regions/patient vs ~0.8 by chance).

## Per-(band, probe) verdict detail

### β cophenet anatomy — ~~strong localized~~ **RETRACTED 2026-05-30**

> **RETRACTED** (signed localization audit). 0/7 of these regions are both
> adequately sampled (n≥5) and cohort-localized. Per-region patient support:
> insula **n=1 (Pat_10, β-anti)**, parahippocampal n=1, postcentral n=2,
> rostralanteriorcingulate n=2, isthmuscingulate n=2, superiorfrontal n=3,
> entorhinal n=5 (the only well-sampled one — and it **fails** the signed test).
> The only well-sampled cophenet-localized region is `Hip` (n=5), which is NOT
> in this list. The enrichments below were driven by direction-blind `|Δρ|` +
> top-decile selection over thinly-sampled regions. Retained as historical record.

Top-decile per-pair `|Δρ_split^coph|` pairs. Seven named DK regions pass A1+A3 jointly (`q_BH<0.05` for A1; `p_emp<0.05` AND `obs_z>2` for A3):
- ctx-lh-isthmuscingulate (A1: 2.33×, p_hyper=2.07e-08; A3: z=6.09, p_emp=0.005)
- ctx-lh-superiorfrontal (1.90×, 4.95e-14; z=7.36, 0.005)
- ctx-rh-insula (2.33×, 7.09e-19; z=5.47, 0.005)
- ctx-lh-parahippocampal (2.13×, 1.03e-06; z=3.59, 0.015)
- ctx-rh-postcentral (1.82×, 6.75e-11; z=2.96, 0.020)
- ctx-lh-entorhinal (1.44×, 1.18e-04; z=2.20, 0.035)
- ctx-rh-rostralanteriorcingulate (2.52×, 1.49e-41; z=2.05, 0.040)

Network: cingulate + medial-temporal + lateral-temporal + sensorimotor + prefrontal.

Source: `data/audit/anatomy_beta_cophenet/cohort_summary.csv` (audit_71, 2026-05-19).

### β Grassmann anatomy — ~~strong localized~~ **RETRACTED 2026-05-30**

> **RETRACTED** (signed localization audit). Of the 7 regions, only **`Hip`**
> crosses the test, and only on **3–5/10 patients** (5 sample it at Grassmann;
> 3 individually significant at β cophenet) — the same thin-sampling zone that
> disqualified every other region, so it does **not** support a cohort
> localization claim either. It is recorded as the single *closest-to-signal*
> region — a **marginal hint, not a localization**. `ctx-lh-insula` (n=4) is
> **anti-localized** (median deviation −0.09); the rest are thin (n=2–3). Note the
> Grassmann quantity is **phase-averaged participation** (an *anchor* — where
> leading modes live — not a cross-phase trace). Retained as historical record.

Top-decile per-node participation in `U_k`, k=27..55. Seven named DK regions pass A3 alone:
- **Hip** (A1: 1.50×; A3: z>>2, p_emp=0.005)
- ctx-lh-insula (0.69×; >>2, 0.005)
- ctx-lh-middletemporal (0.88×; 5.67, 0.005)
- ctx-lh-lateralorbitofrontal (1.72×; >>2, 0.005)
- ctx-rh-medialorbitofrontal (1.02×; >>2, 0.005)
- ctx-lh-superiortemporal (0.91×; 28.28, 0.005)
- ctx-rh-rostralmiddlefrontal (2.36×; 3.58, 0.030)

Network: Hippocampus + temporal cortex + orbitofrontal + insula (left, opposite hemisphere to cophenet's right insula) + rostral middle frontal.

A1 below 1.0× for several A3-passing regions — interpretation: these regions are A3-significant because matched-strength surrogate enrichments cluster strictly below the observation, not because the observed enrichment is itself absolutely high. Standard reading under A3-only.

Source: `data/audit/anatomy_beta_grassmann/cohort_summary.csv` (audit_72, 2026-05-19).

### α cophenet anatomy — ~~strong localized~~ **RETRACTED 2026-05-30 (diffuse)**

> **RETRACTED** (signed localization audit). 0/11 cohort-supported; the cohort
> statistic is **diffuse** (perm-null p ≈ 0.16–0.26, not significant) even before
> the sampling lens. All 11 locked regions are n≤3; no adequately-sampled (n≥5)
> region localizes at all. The C5 epi-X "reproduces identically" claim is moot
> (the underlying localization does not hold). Retained as historical record.

Top-decile per-pair `|Δρ_split^coph|` pairs. Eleven named DK regions pass A1+A3 jointly:
- ctx-lh-caudalanteriorcingulate (A1: 2.06×, p_hyper=4.07e-06; A3: z=4.67, p_emp=0.005)
- ctx-lh-parahippocampal (3.80×, 7.53e-27; z=9.76, 0.005)
- ctx-lh-rostralanteriorcingulate (3.04×, 3.66e-21; z=14.26, 0.005)
- ctx-rh-medialorbitofrontal (1.77×, 1.22e-17; z=3.22, 0.005)
- ctx-rh-caudalmiddlefrontal (2.90×, 2.49e-22; z=6.23, 0.005)
- ctx-rh-postcentral (3.47×, 6.70e-67; z=4.47, 0.005)
- ctx-lh-caudalmiddlefrontal (1.49×, 4.10e-04; z=3.02, 0.015)
- ctx-rh-caudalanteriorcingulate (2.44×, 2.85e-09; z=3.45, 0.020)
- ctx-rh-precuneus (1.47×, 2.92e-04; z=2.60, 0.025)
- ctx-rh-superiorparietal (1.64×, 7.94e-06; z=3.64, 0.025)
- ctx-rh-posteriorcingulate (1.99×, 7.39e-11; z=3.13, 0.045)

Network: bilateral cingulate (anterior caudal + anterior rostral + posterior) + left parahippocampal + right medial OFC + bilateral caudal middle frontal + right postcentral + right precuneus + right superior parietal.

**C5 epi-X reproduces identically** — all 11 named DK regions pass A1+A3 at the same enrichments and p_emp values when epi-zone contacts are removed per patient. No new region emerges; no region drops out. Confirms α anatomy is **wholly non-epi-cortex driven**.

Sources: `data/audit/anatomy_alpha_cophenet/cohort_summary.csv` + `data/audit/anatomy_alpha_cophenet_epiX/cohort_summary.csv` (audit_71 + audit_71 --epi-x, 2026-05-19).

### γ_l Grassmann anatomy — ~~strong localized~~ **RETRACTED 2026-05-30**

> **RETRACTED** (signed localization audit, cluster-extent S(γ_l) k-set). 0/7
> cohort-supported: locked regions are n≤4 except `ctx-lh-middletemporal` (n=5),
> which **fails** the signed test; no adequately-sampled region localizes. The
> A3-only enrichments below reflect phase-averaged participation concentration in
> thinly-sampled regions, not cohort localization. Retained as historical record.

### γ_l Grassmann anatomy — **strong localized** (A3 alone; occipito-temporal + frontal + medial-OFC; locked under `S(γ_l)` cluster-extent paradigm 2026-05-19 pm)

Top-decile per-node participation in `U_k` aggregated (unweighted average) over
`S(γ_l) = {k : p_k(γ_l) < α_k}` — the 41-cell support of the cluster-mass
statistic `T_G^*`. The retired `K*(γ_l) = [12, 23]` (longest-contiguous-significant
window, 13 cells) is **not used**. Seven named DK regions pass A3 alone:

- ctx-lh-lateraloccipital (A1: 2.43×; A3: z>>2 very large, p_emp=0.005)
- ctx-lh-middletemporal (1.59×; 3.78, 0.005)
- ctx-lh-rostralmiddlefrontal (2.37×; 16.30, 0.005)
- ctx-lh-superiortemporal (1.82×; 7.16, 0.005)
- ctx-rh-medialorbitofrontal (1.54×; 3.36, 0.005)
- ctx-rh-parstriangularis (3.89×; 9.01, 0.005)
- ctx-lh-cuneus (2.43×; 5.25, 0.040)

Network: **occipito-temporal + frontal + medial-OFC** (left lateral occipital
+ left cuneus + left middle/superior temporal + left rostral middle frontal +
right pars triangularis + right medial OFC).

**Cross-band note (KC-era retraction).** The earlier `K*(γ_l)` audit had picked
up left fusiform under the retired contiguous-window aggregation. **Left
fusiform does NOT appear under `S(γ_l)`** — fully retracted from the γ_l
anatomy. The KC-era "Hippocampus + left fusiform" β claim is partially
retracted (Hip retained at β Grassmann; left fusiform appears at no band
under `S(b)` everywhere).

Source: `data/audit/anatomy_low_gamma_grassmann_clusterext/cohort_summary.csv`
(audit_72 --cluster-extent, 2026-05-19 pm). The retired
`anatomy_low_gamma_grassmann/cohort_summary.csv` is superseded.

### δ Grassmann anatomy — ~~strong localized~~ **RETRACTED 2026-05-30**

> **RETRACTED** (signed localization audit, cluster-extent S(δ)/S^epiX(δ) k-sets,
> with the epi-X mask BUGFIXED). Full: 0/4 cohort-supported — the only passing
> locked region is `ctx-lh-inferiorparietal` **n=1**, while the well-sampled
> `ctx-lh-superiortemporal` (n=4) and `ctx-lh-inferiortemporal` (n=5) are
> **anti-localized** (−0.10, −0.03). Epi-X: 0/3 — **no** locked region localizes;
> `ctx-rh-rostralmiddlefrontal` is anti-localized. NB the locked δ epi-X was
> computed with a no-op epi mask (audit_72 `np.isin(int, str)` bug, fixed
> 2026-05-30); it never excluded epileptic nodes. The "fully disjoint networks"
> framing below is moot (neither network is cohort-localized). Historical record.

### δ Grassmann anatomy — **strong localized; full and epi-X read FULLY DISJOINT NETWORKS** (locked under `S(δ)` / `S^epiX(δ)` cluster-extent paradigm 2026-05-19 pm)

Top-decile per-node participation in `U_k` aggregated (unweighted average) over
`S(δ) = {k : p_k(δ) < α_k}` (23 cells) for the full cohort and
`S^epiX(δ)` (22 cells, derived from per-k cohort Wilcoxon on
`grassmann_epi_exclusion/per_patient_per_band_per_k.csv`, span
`k = [2, 34..52, 87, 88]`) for the C5 epi-X analysis. The retired
`K*(δ) = [57, 63]` (7 cells) and `K*^epiX(δ) = [33, 39]` (7 cells)
contiguous windows are **not used**. The dissociation strengthens from
"2/6 shared" (under retired `K*`) to **0/3 shared** (under locked `S(b)`)
— δ Grassmann is now read as a fully disjoint mixture of two distinct
phenomena.

**Full-cohort (over `S(δ)`), 4 named DK regions A3-pass**:
- ctx-lh-inferiortemporal (A1: 1.34×; A3: z>>2 very large, p_emp=0.005)
- ctx-lh-inferiorparietal (3.89×; >>2 very large, 0.005)
- ctx-rh-parstriangularis (0.65×; >>2 very large, 0.005)
- ctx-lh-superiortemporal (0.61×; 14.11, 0.010)

Network: **temporal + parietal + frontal** (left inferior temporal + left
superior temporal + left inferior parietal + right pars triangularis). The
"anchor anatomy" interpretation (Amygdala + cingulate + medial OFC +
fusiform + bankssts) of the retired `K*(δ)` analysis is **not supported**
under `S(δ)` — replaced by a parietal-temporal-frontal network.

**C5 epi-X (over `S^epiX(δ)`), 3 named DK regions A3-pass**:
- ctx-lh-superiorparietal (A1: 9.73×; A3: z>>2 very large, p_emp=0.005)
- ctx-rh-rostralmiddlefrontal (1.47×; 2.71, 0.005)
- ctx-lh-superiorfrontal (0.88×; 14.11, 0.010)

Network: **left superior parietal + right rostral middle frontal + left
superior frontal**. Fusiform, inferior parietal, postcentral, and inferior
temporal that appeared in the retired `K*^epiX(δ)` analysis are **not
present** under `S^epiX(δ)`; left superior frontal emerges as new.

**Shared regions between `S(δ)` and `S^epiX(δ)`: 0/3**. The two networks are
fully disjoint under the cluster-extent paradigm (was 2/6 shared under the
retired contiguous windows; dissociation strengthens to "completely
non-overlapping").

**Reading**: full-cohort δ Grassmann (over `S(δ)`) reads a temporal-parietal-frontal
cortical pattern; C5 epi-X (over `S^epiX(δ)`) reveals a distinct superior-parietal
+ frontal pattern. The two windows share no DK regions under the locked
methodology. The preprint should report both networks and note the strengthened
dissociation. The δ cross-probe anchor anatomy (1.55× cross-probe ratio, known
biology per `memory/epileptic_imcoh_universal.md`) is a **separate, descriptive
substrate-layer observation** (LEDGER Decision 5), independent of the Grassmann
anatomy result above.

Sources: `data/audit/anatomy_delta_grassmann_clusterext/cohort_summary.csv` +
`data/audit/anatomy_delta_grassmann_epiX_clusterext/cohort_summary.csv` (audit_72
--cluster-extent, full + --epi-x --k-list ..., 2026-05-19 pm). The retired
`anatomy_delta_grassmann{_epiX}/cohort_summary.csv` are superseded.

## Cross-band anatomy comparison

> **SUPERSEDED 2026-05-30.** The cross-band region motifs below rest on the
> per-(band, probe) region-lists that the signed localization audit retracted
> (see Head banner + 2026-05-30 revision). **No motif survives the sampling-aware
> test** — including the Hippocampus, which crosses the test on only 3–5/10
> patients (a marginal hint, not a localization). Treat the overlap table as
> historical; do not cite any multi-band motif as a cohort finding.

### Region overlaps across trace-positive bands (locked under `S(b)` cluster-extent paradigm, 2026-05-19 pm)

| Region | β cophenet | β Grassmann | α cophenet | γ_l Grassmann | δ Grassmann (full) | δ Grassmann (epi-X) |
|---|---|---|---|---|---|---|
| Hippocampus / parahippocampal | parahippocampal | **Hip** | parahippocampal | — | — | — |
| Insula (any hemisphere) | right | left | — | — | — | — |
| Cingulate (any subregion) | isthmus + ACC (rh) | — | ACC bilateral + PCC (rh) | — | — | — |
| Medial / lateral OFC | — | bilateral | medial (rh) | medial (rh) | — | — |
| Fusiform | — | — | — | — | — | — |
| Inferior parietal | — | — | — | — | left | — |
| Superior parietal / postcentral | postcentral (rh) | — | superior + postcentral (rh) | — | — | superior (lh) |
| Temporal cortex (mid / sup / inf) | — | middle + superior (lh) | — | middle + superior (lh) | inferior + superior (lh) | — |
| Amygdala | — | — | — | — | — | — |
| Bankssts | — | — | — | — | — | — |
| Pars (orbitalis / triangularis) | — | — | — | triangularis (rh) | triangularis (rh) | — |
| Rostral / caudal middle frontal | — | rostral (rh) | caudal bilateral | rostral (lh) | — | rostral (rh) |
| Superior frontal | left | — | — | — | — | left |
| Precuneus | — | — | right | — | — | — |
| Paracentral | — | — | — | — | — | — |
| Lateral occipital | — | — | — | left | — | — |
| Cuneus | — | — | — | left | — | — |

**Pattern (under locked `S(b)` cluster-extent paradigm)**:
- **Parahippocampal / Hippocampus** appears at β (both probes) and α cophenet — beta-and-alpha medial-temporal involvement (unchanged from KC era).
- **Medial OFC** appears at β Grassmann + α cophenet + γ_l Grassmann — three-band medial-OFC motif (newly visible under `S(γ_l)` after cluster-extent rerun).
- **Temporal cortex (middle / superior)** appears at β Grassmann + γ_l Grassmann + δ Grassmann (full) — three-band temporal motif. Inferior temporal at δ Grassmann (full) only.
- **Pars triangularis (right)** appears at γ_l Grassmann + δ Grassmann (full) — gamma-low and delta share a right inferior frontal node.
- **Parietal cortex** (inferior + superior + postcentral) appears at β cophenet (postcentral) + α cophenet (superior + postcentral + precuneus) + δ Grassmann full (inferior) + δ Grassmann epi-X (superior) — a recurring parietal motif preserved across paradigms.
- **Cingulate** appears at β cophenet + α cophenet only (no longer at δ Grassmann under `S(δ)` — the retired `K*(δ)` caudal ACC entry retracts).
- **Insula** appears only at β (both probes, opposite hemispheres) — unchanged.

**KC-era retractions under `S(b)`**:
- **Left fusiform** appears at **NO band** — fully retracted from γ_l Grassmann, δ Grassmann (full), and δ Grassmann (epi-X). The KC-era "Hippocampus + left fusiform" β claim is partially retracted (Hip retained at β Grassmann; left fusiform appears nowhere under the locked cluster-extent paradigm).
- **Amygdala** retracts from δ Grassmann (full) — the "anchor anatomy" interpretation (Amy + cingulate + medial OFC + fusiform + bankssts) is not supported under `S(δ)`. The δ cross-probe anchor anatomy (1.55× cross-probe ratio, known biology per `memory/epileptic_imcoh_universal.md`) is a separate descriptive substrate-layer observation, not a Grassmann anatomy claim.
- **Bankssts** retracts from δ Grassmann (full).
- **Right caudal anterior cingulate** retracts from δ Grassmann (full).
- **Right paracentral** retracts from γ_l Grassmann.

## Manuscript ↔ lab label mapping (locked 2026-05-19 pm)

The manuscript methods §sssec:methods_anatomy uses **A1 (hypergeometric) +
A2 (matched-strength surrogate)** in compressed form. In this ledger and in
all lab-side artifacts the labels **A1 / A3** are kept for traceability with
the locked battery. Mapping:

| Manuscript label | Lab label | Test |
|---|---|---|
| **A1** | A1 | Hypergeometric per-region enrichment |
| **A2** | **A3** | Matched-strength surrogate enrichment |
| — | A2 | Sampling-corrected bootstrap (sensitivity supplement, not in manuscript) |
| — | A4 | Implant-geometry cohort regression (sensitivity supplement, not in manuscript) |

All region-level entries below (z values, p_emp values, "passes A3" tags) are
expressed in lab labels. Manuscript LaTeX uses the manuscript convention
without renaming the underlying tests. Full mapping rationale in
`ANATOMY_CONTROLS.md`.

## Anti-revisitation clause

A verdict in this ledger is **locked**. Re-evaluation requires:
1. A new dated audit run producing fresh `data/audit/anatomy_<band>_<probe>[_epiX]/cohort_summary.csv`.
2. A dated revision entry in this ledger citing the new audit.
3. Cascading updates in the per-band briefs.

KC-era anatomy memory entries (`result_2_lrg_beta_trace.md` "Hippocampus + left fusiform" claim) are not citable. Hippocampus survives at β Grassmann; **left fusiform appears nowhere under the locked cluster-extent paradigm `S(b)` — neither at β, γ_l, nor δ** (audit_72 --cluster-extent rerun 2026-05-19 pm). The earlier "γ_l + δ" re-attribution from the retired `K*(b)` audit is itself withdrawn.

## What this ledger does not lock

- **Mechanistic interpretation** of why a region enriches (left to the writing agent / neuroscientific discussion).
- **Descriptive δ anchor anatomy** at C4 cross-probe (1.55× ratio, `memory/epileptic_imcoh_universal.md`) — known biology, logged in `VERDICT_LEDGER.md` Decision 5 as separate descriptive observation, not part of the δ trace localization (which is at the Grassmann probe, audited above).
- **A2 sampling-corrected bootstrap and A4 implant-geometry regression** — deferred to sensitivity supplement. The locked anatomy verdict rests on A1+A3 (cophenet probes) and A3 alone (Grassmann probes).

## Revision history

- **2026-06-12** — **Framing clarification (no new audit; same verdict).** Reframed
  the β verdict from "localizes to OFC" to "**brain-wide trace that concentrates in
  OFC** above a per-patient-demeaned baseline" — OFC is a reproducible *hotspot*, not
  a container (the trace exists cohort-wide; even within OFC the directional lean is
  only ≈0.59). Also surfaced the **scale-dependence** from the existing atlas data:
  the concentration **peaks at the a-priori *system* scale** and is **bilateral**
  (no L/R lateralization, hemisphere p≈0.48/0.76), **washes out** at the whole-frontal-
  lobe scale (p≈0.29), and is **undersampled at single-region granularity** (fails BH
  over 53 regions — the same wall α hits). Hemisphere/lobe were descriptive
  (label-shuffle, rung ①) only; the matched-strength gate covered
  region/system/supersystem/paracore, and the gated coarse units (limbic supersystem,
  paralimbic core) point back to OFC. Top banner + verdict row + `00_cohort.md` /
  `01_beta.md` heads + `localization_atlas/README.md` updated; the 2026-06-10 entry
  below retains its original "localizes" wording as historical record.

- **2026-06-10** — **β anatomy REINSTATED → orbitofrontal cortex (new audited verdict; satisfies the anti-revisitation clause).** The 2026-06-05 sampling-conditioned exploration was carried through a six-rung control gate + a self-audit. New dated audits: `audit_83_localization_matched_strength.py` (now canonical observed trace + `--shaft-collapse` + `--R`), `audit_92_localization_R1000_surrogates.py` (R=1000 ensemble), `diag_localization_per_patient_decomposition.py`; CSVs + report under `data/audit/localization_atlas/`. **β cophenet trace localizes to the OFC system**: matched-strength R=1000, rank concordance, BH q=0.009–0.013 across the 9 a-priori systems within β in all four conditions; 4/5 implanted positive; leave-one-out robust; **shaft-collapse robust** (multi-shaft in the 2 strongest patients); low-strength (non-hub). **Self-audit fixes:** (1) observed trace was read from the LRG-pipeline ultrametric while the surrogate used the canonical `cophenetic_condensed_from_eigs` — they differ by a per-phase ρ-trace scale that perturbs cross-phase ranks (Spearman 0.92–0.99); `obs_trace` now recomputes canonically (matches surrogate + locked §5.3 audit_63); OFC verdict unchanged. (2) **α demoted** to no FDR-surviving localization (pars opercularis fails BH over 53 regions; corrects the 2026-06-05 "frontal operculum" line). **Retractions held:** multi-region DK lists stay retracted; "distributed paralimbic ring" retracted (fails shaft-collapse); Hip/MTL = sub-threshold hint (p≈0.4, LOO-fragile); occipital = K=3 low-coverage secondary; Grassmann subspace + low-γ do not localize. Multiplicity family = a-priori systems within β only (NOT bands, NOT nested granularities, NOT pooled epi modes). Top banner updated; β/α rows in the locked-verdicts table updated. Adversarial review (a3dc4b5d) + self-audit (2026-06-08/10). Cascade into `00_cohort.md` / `01_beta.md` landed same day; commit pending user green-light.

- **2026-06-05** — **VERDICT CORRECTION (interpretation, not a new audit).** Re-examination
  of the *same* 2026-05-30 audit CSVs (`data/audit/anatomy_localization_wilcoxon/`) plus a
  fresh per-patient Hip-β cophenet recomputation found the retraction correct about the
  multi-region lists but over-aggressive on "no anatomy anywhere." (1) The localization
  denominator is patients **implanted in** a region, not /10; the 5/10 ceiling is implant
  coverage, so the ">5/10 or thin" bar is a structural false-negative on n=10. (2) The
  buried genuine signal is **Hip in β COPHENET (trace `s_ij=dD_task·dD_rest`)**, distinct
  from the β Grassmann Hip (phase-avg ANCHOR) the retraction leaned on. Corrected reading:
  **4/5 hippocampus-implanted patients** show the β trace concentrating in Hip above their
  own baseline (median +0.0073 > q95 +0.0042, perm_p 0.0046; exact Wilcoxon 0.0625
  borderline). **Multi-region lists stay retracted; the "spatially delocalized" blanket is
  withdrawn**, replaced by a descriptive, sampling-conditioned frame. Broader localization
  exploration (sampling-conditioned atlas + MTL/limbic system pooling + epi-X-resolved)
  underway before any re-lock. Banner at top of this ledger + audit README 2026-06-05 banner.

- **2026-05-19** — Ledger locked. Five audit runs landed: β cophenet (audit_71), β Grassmann (audit_72), α cophenet (audit_71, full + C5 epi-X), γ_l Grassmann (audit_72), δ Grassmann (audit_72, full k=57..63 + C5 epi-X k=33..39). All four trace-positive bands locked as **strong localized**. KC-era "Hippocampus + left fusiform" claim partially retracted (Hip retained at β Grassmann; fusiform shifts to γ_l + δ).

- **2026-05-19 (pm)** — **Methodology fix landed**. Grassmann anatomy regenerated under the locked all-clusters paradigm with per-node participation aggregated over `S(b) = {k : p_k(b) < α_k}` (the support of `T_G^*`), unweighted across `k`. The retired aggregation over `K*(b)` (longest-contiguous-significant window) is replaced everywhere. **Per-band cluster-extent diff** (locked CSVs at `data/audit/anatomy_<band>_grassmann{_epiX}_clusterext/`):
    - **β Grassmann**: `|S(β)| = 40` cells (k ∈ [21, 90] non-contiguous, was K*=[27,55] 29 contiguous). A3-passing region set is **identical 7/7** to the retired ledger: Hip, ctx-lh-insula, ctx-lh-middletemporal, ctx-lh-lateralorbitofrontal, ctx-rh-medialorbitofrontal, ctx-lh-superiortemporal, ctx-rh-rostralmiddlefrontal. ctx-rh-rostralmiddlefrontal strengthens (z=3.58 → 5.21, p_emp=0.030 → 0.010). **β anatomy verdict unchanged**.
    - **γ_l Grassmann**: `|S(γ_l)| = 41` cells, was K*=[12,23] 13 contiguous. A3-passing region set **shifted to 7 regions**, only 3/6 overlap with the locked ledger. Retained: ctx-lh-middletemporal, ctx-lh-superiortemporal, ctx-rh-parstriangularis. **Dropped**: ctx-lh-inferiortemporal, ctx-lh-fusiform, ctx-rh-paracentral. **Added**: ctx-lh-lateraloccipital, ctx-lh-rostralmiddlefrontal, ctx-rh-medialorbitofrontal, ctx-lh-cuneus. The "temporal-cortex dominant" narrative weakens — the new network is **occipito-temporal + frontal + medial-OFC**. **The KC-era "left fusiform" claim at γ_l is fully retracted** (was carried under K*; not in S(γ_l)).
    - **δ-full Grassmann**: `|S(δ)| = 23` cells, was K*=[57,63] 7 contiguous. A3-passing region set **shrinks to 4 regions**, only 1/6 overlaps. Retained: ctx-lh-inferiorparietal. **Dropped**: ctx-lh-fusiform, ctx-rh-medialorbitofrontal, ctx-rh-caudalanteriorcingulate, **Amy**, ctx-lh-bankssts. **Added**: ctx-lh-inferiortemporal, ctx-rh-parstriangularis, ctx-lh-superiortemporal. **The "anchor anatomy" interpretation (Amy + cingulate + fusiform) is not supported under S(δ)** — replaced by a parietal-temporal-frontal network. The δ cross-probe 1.55× anchor-anatomy known-biology layer (`memory/epileptic_imcoh_universal.md`) is a separate descriptive observation and is not affected.
    - **δ-epi-X Grassmann**: `|S^epiX(δ)| = 22` cells (k = [2, 34..52, 87, 88], derived from per-k Wilcoxon on `data/audit/grassmann_epi_exclusion/per_patient_per_band_per_k.csv`), was K*=[33,39] 7 contiguous. A3-passing region set **shrinks to 3 regions**, only 2/6 overlap. Retained: ctx-lh-superiorparietal, ctx-rh-rostralmiddlefrontal. **Dropped**: ctx-lh-fusiform, ctx-lh-inferiorparietal, ctx-rh-postcentral, ctx-lh-inferiortemporal. **Added**: ctx-lh-superiorfrontal. The parietal-dominant narrative is preserved (superior parietal kept).
    - **δ full-vs-epiX dissociation strengthens under S(b)**: the two networks now share **0/3 named regions** (locked ledger had 2/6 shared via fusiform + inferiorparietal). The "mixture of two distinct phenomena" claim becomes stronger.

    α cophenet and β cophenet verdicts are unaffected (cophenet anatomy uses top-decile per-pair `|Δρ_split^coph|`, no `k`-aggregation involved).

    The Per-(band, probe) verdict-detail subsections above remain in this ledger as **historical reference under K*(b)** until cascaded edits land in `bands/01_beta.md` and the β results LaTeX. The cluster-extent CSVs at `data/audit/anatomy_<band>_grassmann{_epiX}_clusterext/cohort_summary.csv` are the **new source of truth** for any per-band anatomy citation.

- **2026-05-30** — **Anatomy localization RETRACTED via a signed, sampling-aware test.**
  New dated audit (satisfies the anti-revisitation clause): script
  `scripts/01_compute/diagnostics/diag_anatomy_localization_wilcoxon.py`, report +
  per-region CSVs under `data/audit/anatomy_localization_wilcoxon/`.

    **Why the locked verdicts were wrong.** audit_71/72 enrichment is (i)
    **direction-blind** (`|Δρ|` / participation magnitude mixes trace and
    anti-trace) and (ii) gated by a **top-decile** threshold, then tallied as
    region **endpoint counts** — a hypergeometric/surrogate enrichment that can be
    driven by a handful of high-magnitude, thinly-sampled pairs. The new test uses
    the **signed** score `s_ij = dD_task·dD_rest` (cophenet) / per-node
    participation deviation (Grassmann), **per-patient demeaned** (removes the
    cohort base rate so it measures localization, not trace existence), aggregated
    per region by the **across-patient median**, and scored by a **floor-free
    permutation rank** against a within-patient region-label-shuffle null. No
    hardcoded thresholds gate the verdict (the permutation null is the gate).

    **Decisive fact = sampling.** No DK region is sampled by more than **5/10**
    patients in any cell; only 4 regions reach n=5 anywhere (Hip, entorhinal,
    inferiortemporal, middletemporal). The locked region-lists are dominated by
    1–4 patient findings. Requiring adequate sampling (n≥5) AND localization AND
    locked-membership: **β cophenet 0/7, α cophenet 0/11, β Grassmann only Hip
    crosses (and only on 3–5/10 patients), γ_l Grassmann 0/7, δ Grassmann full
    0/4, δ epi-X 0/3.** β `ctx-rh-insula` "localizes" on **Pat_10 alone — the
    β-anti patient**; δ's well-sampled locked regions (`superiortemporal` n=4,
    `inferiortemporal` n=5) are **anti-localized**.

    **Verdicts.** Retract **all six** cells' anatomy localization claims (β/α
    cophenet, β Grassmann, γ_l Grassmann, δ Grassmann full + epi-X). **No region —
    including the Hippocampus — reaches a defensible cohort localization**: Hip is
    the single region that crosses the test, but on only 3–5/10 patients, the same
    thin-sampling zone every retracted region sits in, so it is recorded as a
    *marginal hint, not a localization* (and is an anchor-flavored phase-averaged
    Grassmann quantity, not a trace). β/α cophenet were user-approved 2026-05-30;
    β Grassmann + γ_l/δ retractions extend the identical finding from the same run.
    **Correction (same day):** an earlier draft "softened β Grassmann to
    Hippocampus only" — that was an inconsistency (3–5/10 patients is exactly the
    coverage that disqualified the retracted regions) and is withdrawn; β Grassmann
    is retracted like the rest.

    **Two collateral findings.** (1) `audit_72` epi-X mask was a **no-op bug**
    (`np.isin(np.arange(N), string_labels)` → excludes 0 nodes; verified Pat_02
    0/117 vs correct 14) — the locked δ Grassmann epi-X never excluded epileptic
    nodes; **fixed** 2026-05-30 (label-based mask, audit_71 pattern). (2) The
    Grassmann anatomy quantity (`participation_phase_avg`) averages phases rather
    than differencing them — it is an **anchor** ("where leading modes live"),
    not a cross-phase trace.

    **Adversarial review.** The signed test was independently reviewed
    (statistician + code-reviewer + adjudicator); the review caught and the script
    fixed a `−log10 p` mass floor-dilution (false-negative) and a base-rate
    confound. Implementation (region groupby, index alignment incl. Pat_10 at 113,
    null calibration) verified sound.

    **Cascade TODO** (not yet done): per-band briefs (`bands/01_beta.md`,
    `bands/02_alpha.md`, `bands/03_gammalow.md`, `bands/06_delta.md`),
    `bands/00_cohort.md`, `HANDOFF_INDEX.md`, and `VERDICT_LEDGER.md` anatomy
    references. LaTeX is out of scope unless the user requests it.
