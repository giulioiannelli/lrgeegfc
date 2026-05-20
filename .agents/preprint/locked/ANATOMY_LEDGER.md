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

**Head.** This ledger is the **single source of truth** for where each trace lives anatomically, under the locked anatomy control battery (`ANATOMY_CONTROLS.md`, A1 hypergeometric + A3 matched-strength surrogate; A2 sampling-corrected bootstrap and A4 implant-geometry regression deferred to sensitivity supplement). All audits are now run; verdicts below are **locked**. KC-era anatomy artifacts (`lrg_localization_anatomy/`, the "Hippocampus + left fusiform" memory claim) are **retired** under the 2026-05-18 trace-side lockdown and **not citable** for preprint claims.

## Locked verdicts

| Band | Trace verdict (CONTROLS) | Probe(s) audited | Anatomy verdict (ANATOMY_CONTROLS) |
|---|---|---|---|
| **β** | strong trace, both probes | cophenet + Grassmann | **strong localized, both probes** (7+7 named DK regions) |
| α | strong trace, only D_coph | cophenet (full + C5 epi-X) | **strong localized, only D_coph; epi-X reproduces identically** (11 named regions both) |
| γ_l | weak trace, only Grassmann | Grassmann | **strong localized, only Grassmann** (6 named DK regions, temporal-cortex-dominant) |
| δ | weak trace, only Grassmann | Grassmann (full + C5 epi-X) | **strong localized, only Grassmann; full and epi-X read DIFFERENT NETWORKS** (6 + 6 named regions, only 2 shared) |
| θ | no trace | — | n/a |
| γ_h | no trace | — | n/a |

**Aggregate reading**: the trace, where it exists, is **localized to band-specific distributed cortical networks** — never single-region, never diffuse-brain-wide. The 5-control trace battery and the 4-control anatomy battery agree: when a trace exists, it concentrates anatomically under matched-strength surrogacy.

## Per-(band, probe) verdict detail

### β cophenet anatomy — **strong localized**

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

### β Grassmann anatomy — **strong localized** (A3 alone; A1 sparse)

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

### α cophenet anatomy — **strong localized; identical under C5 epi-X**

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

### γ_l Grassmann anatomy — **strong localized** (A3 alone; temporal-cortex-dominated)

Top-decile per-node participation in `U_k`, k=12..23. Six named DK regions pass A3 alone:
- ctx-lh-inferiortemporal (A1: 1.01×; A3: z=3.05, p_emp=0.005)
- ctx-lh-fusiform (0.26×; >>2, 0.005)
- ctx-lh-middletemporal (1.95×; 5.00, 0.005)
- ctx-lh-superiortemporal (2.13×; 11.06, 0.005)
- ctx-rh-paracentral (1.62×; >>2, 0.005)
- ctx-rh-parstriangularis (3.24×; 10.13, 0.005)

Network: **temporal cortex dominant** (left fusiform + left inferior/middle/superior temporal) + right paracentral + right pars triangularis.

**Cross-band note**: KC-era "left fusiform" claim at β is retracted — left fusiform passes A1+A3 here at γ_l Grassmann, not at β.

Source: `data/audit/anatomy_low_gamma_grassmann/cohort_summary.csv` (audit_72, 2026-05-19).

### δ Grassmann anatomy — **strong localized, two non-overlapping networks across full and epi-X windows**

The full-cohort window (k=57..63) and C5 epi-X window (k=33..39) are non-overlapping in k (`VERDICT_LEDGER.md` Decision 6) AND non-overlapping in anatomy (only 2/6 named regions shared). This is the strongest evidence that δ Grassmann is a **mixture of two distinct phenomena**.

**Full-cohort (k=57..63), 6 named DK regions A3-pass**:
- ctx-lh-fusiform (A1: 1.05×; A3: z=12.62, p_emp=0.005)
- ctx-lh-inferiorparietal (7.78×; 5.07, 0.005)
- ctx-rh-medialorbitofrontal (1.02×; 4.08, 0.005)
- ctx-rh-caudalanteriorcingulate (9.73×; 2.84, 0.005)
- **Amy** (4.42×; 2.33, 0.030)
- ctx-lh-bankssts (4.86×; 2.42, 0.045)

Network: Amygdala + cingulate + OFC + fusiform + inferior parietal + bankssts — overlaps published δ "anchor anatomy" (low-frequency synchronization near epi zones, `memory/epileptic_imcoh_universal.md`).

**C5 epi-X (k=33..39), 6 named DK regions A3-pass**:
- ctx-lh-fusiform (A1: 0.79×; A3: z>>2, p_emp=0.005)
- ctx-lh-inferiorparietal (3.89×; >>2, 0.005)
- ctx-lh-superiorparietal (9.73×; >>2, 0.005)
- ctx-rh-postcentral (1.95×; 3.03, 0.005)
- ctx-rh-rostralmiddlefrontal (1.47×; 3.52, 0.005)
- ctx-lh-inferiortemporal (1.34×; 5.69, 0.035)

Network: parietal-dominant (left inferior + left superior + right postcentral) + left fusiform + rostral middle frontal + left inferior temporal — distinct from full-cohort's anchor-anatomy network.

**Shared regions**: only fusiform + inferior parietal (2/6). The remaining 4 regions in each window are non-overlapping.

**Reading**: full-cohort δ Grassmann is dominated by epi-zone anchor anatomy; C5 epi-X reveals a separate, physiologically interpretable parietal-cortex network. The preprint should report both networks and the dissociation.

Sources: `data/audit/anatomy_delta_grassmann/cohort_summary.csv` + `data/audit/anatomy_delta_grassmann_epiX/cohort_summary.csv` (audit_72 + audit_72 --epi-x, 2026-05-19).

## Cross-band anatomy comparison

### Region overlaps across trace-positive bands

| Region | β cophenet | β Grassmann | α cophenet | γ_l Grassmann | δ Grassmann (full) | δ Grassmann (epi-X) |
|---|---|---|---|---|---|---|
| Hippocampus / parahippocampal | parahippocampal | **Hip** | parahippocampal | — | — | — |
| Insula (any hemisphere) | right | left | — | — | — | — |
| Cingulate (any subregion) | isthmus + ACC (rh) | — | ACC bilateral + PCC (rh) | — | caudal ACC (rh) | — |
| Medial / lateral OFC | — | bilateral | medial (rh) | — | medial (rh) | — |
| Fusiform | — | — | — | **left** | **left** | **left** |
| Inferior parietal | — | — | — | — | left | left |
| Superior parietal / postcentral | postcentral (rh) | — | superior + postcentral (rh) | — | — | superior (lh) + postcentral (rh) |
| Temporal cortex (mid / sup / inf) | — | middle + superior (lh) | — | inferior + middle + superior (lh) | — | inferior (lh) |
| Amygdala | — | — | — | — | **Amy** | — |
| Bankssts | — | — | — | — | left | — |
| Pars (orbitalis / triangularis) | — | — | — | triangularis (rh) | — | — |
| Rostral / caudal middle frontal | — | rostral (rh) | caudal bilateral | — | — | rostral (rh) |
| Superior frontal | left | — | — | — | — | — |
| Precuneus | — | — | right | — | — | — |
| Paracentral | — | — | — | right | — | — |

**Pattern**:
- **Cingulate** is the most multi-band region — appears in β cophenet, α cophenet, δ Grassmann (full). Across-band cingulate involvement is consistent with cingulate's role as a multi-band integrator.
- **Parahippocampal / Hippocampus** appears at β (both probes) and α cophenet — beta-and-alpha medial-temporal involvement.
- **Fusiform** appears at γ_l Grassmann + δ Grassmann (both windows) — NOT at β (retracting the KC-era β fusiform claim).
- **Insula** appears only at β (both probes, opposite hemispheres).
- **Amygdala** appears only at δ Grassmann full-cohort (consistent with the δ anchor anatomy).
- **Parietal cortex** (inferior + superior + postcentral) appears at β cophenet (postcentral) + α cophenet (superior + postcentral + precuneus) + δ Grassmann epi-X (superior + inferior + postcentral) — a recurring parietal motif.

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

- **2026-05-19** — Ledger locked. Five audit runs landed: β cophenet (audit_71), β Grassmann (audit_72), α cophenet (audit_71, full + C5 epi-X), γ_l Grassmann (audit_72), δ Grassmann (audit_72, full k=57..63 + C5 epi-X k=33..39). All four trace-positive bands locked as **strong localized**. KC-era "Hippocampus + left fusiform" claim partially retracted (Hip retained at β Grassmann; fusiform shifts to γ_l + δ).

- **2026-05-19 (pm)** — **Methodology fix landed**. Grassmann anatomy regenerated under the locked all-clusters paradigm with per-node participation aggregated over `S(b) = {k : p_k(b) < α_k}` (the support of `T_G^*`), unweighted across `k`. The retired aggregation over `K*(b)` (longest-contiguous-significant window) is replaced everywhere. **Per-band cluster-extent diff** (locked CSVs at `data/audit/anatomy_<band>_grassmann{_epiX}_clusterext/`):
    - **β Grassmann**: `|S(β)| = 40` cells (k ∈ [21, 90] non-contiguous, was K*=[27,55] 29 contiguous). A3-passing region set is **identical 7/7** to the retired ledger: Hip, ctx-lh-insula, ctx-lh-middletemporal, ctx-lh-lateralorbitofrontal, ctx-rh-medialorbitofrontal, ctx-lh-superiortemporal, ctx-rh-rostralmiddlefrontal. ctx-rh-rostralmiddlefrontal strengthens (z=3.58 → 5.21, p_emp=0.030 → 0.010). **β anatomy verdict unchanged**.
    - **γ_l Grassmann**: `|S(γ_l)| = 41` cells, was K*=[12,23] 13 contiguous. A3-passing region set **shifted to 7 regions**, only 3/6 overlap with the locked ledger. Retained: ctx-lh-middletemporal, ctx-lh-superiortemporal, ctx-rh-parstriangularis. **Dropped**: ctx-lh-inferiortemporal, ctx-lh-fusiform, ctx-rh-paracentral. **Added**: ctx-lh-lateraloccipital, ctx-lh-rostralmiddlefrontal, ctx-rh-medialorbitofrontal, ctx-lh-cuneus. The "temporal-cortex dominant" narrative weakens — the new network is **occipito-temporal + frontal + medial-OFC**. **The KC-era "left fusiform" claim at γ_l is fully retracted** (was carried under K*; not in S(γ_l)).
    - **δ-full Grassmann**: `|S(δ)| = 23` cells, was K*=[57,63] 7 contiguous. A3-passing region set **shrinks to 4 regions**, only 1/6 overlaps. Retained: ctx-lh-inferiorparietal. **Dropped**: ctx-lh-fusiform, ctx-rh-medialorbitofrontal, ctx-rh-caudalanteriorcingulate, **Amy**, ctx-lh-bankssts. **Added**: ctx-lh-inferiortemporal, ctx-rh-parstriangularis, ctx-lh-superiortemporal. **The "anchor anatomy" interpretation (Amy + cingulate + fusiform) is not supported under S(δ)** — replaced by a parietal-temporal-frontal network. The δ cross-probe 1.55× anchor-anatomy known-biology layer (`memory/epileptic_imcoh_universal.md`) is a separate descriptive observation and is not affected.
    - **δ-epi-X Grassmann**: `|S^epiX(δ)| = 22` cells (k = [2, 34..52, 87, 88], derived from per-k Wilcoxon on `data/audit/grassmann_epi_exclusion/per_patient_per_band_per_k.csv`), was K*=[33,39] 7 contiguous. A3-passing region set **shrinks to 3 regions**, only 2/6 overlap. Retained: ctx-lh-superiorparietal, ctx-rh-rostralmiddlefrontal. **Dropped**: ctx-lh-fusiform, ctx-lh-inferiorparietal, ctx-rh-postcentral, ctx-lh-inferiortemporal. **Added**: ctx-lh-superiorfrontal. The parietal-dominant narrative is preserved (superior parietal kept).
    - **δ full-vs-epiX dissociation strengthens under S(b)**: the two networks now share **0/3 named regions** (locked ledger had 2/6 shared via fusiform + inferiorparietal). The "mixture of two distinct phenomena" claim becomes stronger.

    α cophenet and β cophenet verdicts are unaffected (cophenet anatomy uses top-decile per-pair `|Δρ_split^coph|`, no `k`-aggregation involved).

    The Per-(band, probe) verdict-detail subsections above remain in this ledger as **historical reference under K*(b)** until cascaded edits land in `bands/01_beta.md` and the β results LaTeX. The cluster-extent CSVs at `data/audit/anatomy_<band>_grassmann{_epiX}_clusterext/cohort_summary.csv` are the **new source of truth** for any per-band anatomy citation.
