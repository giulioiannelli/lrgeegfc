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
| γ_l | strong trace, only Grassmann ↑ | Grassmann | **strong localized, only Grassmann** (7 named DK regions, occipito-temporal + frontal + medial-OFC) |
| δ | weak trace, only Grassmann (LOO Pat_08 fails Decision-12 precondition) | Grassmann (full + C5 epi-X) | **strong localized, only Grassmann; full and epi-X read FULLY DISJOINT NETWORKS** (4 + 3 named regions, **0 shared** under `S(b)`) |
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

- **2026-05-19** — Ledger locked. Five audit runs landed: β cophenet (audit_71), β Grassmann (audit_72), α cophenet (audit_71, full + C5 epi-X), γ_l Grassmann (audit_72), δ Grassmann (audit_72, full k=57..63 + C5 epi-X k=33..39). All four trace-positive bands locked as **strong localized**. KC-era "Hippocampus + left fusiform" claim partially retracted (Hip retained at β Grassmann; fusiform shifts to γ_l + δ).

- **2026-05-19 (pm)** — **Methodology fix landed**. Grassmann anatomy regenerated under the locked all-clusters paradigm with per-node participation aggregated over `S(b) = {k : p_k(b) < α_k}` (the support of `T_G^*`), unweighted across `k`. The retired aggregation over `K*(b)` (longest-contiguous-significant window) is replaced everywhere. **Per-band cluster-extent diff** (locked CSVs at `data/audit/anatomy_<band>_grassmann{_epiX}_clusterext/`):
    - **β Grassmann**: `|S(β)| = 40` cells (k ∈ [21, 90] non-contiguous, was K*=[27,55] 29 contiguous). A3-passing region set is **identical 7/7** to the retired ledger: Hip, ctx-lh-insula, ctx-lh-middletemporal, ctx-lh-lateralorbitofrontal, ctx-rh-medialorbitofrontal, ctx-lh-superiortemporal, ctx-rh-rostralmiddlefrontal. ctx-rh-rostralmiddlefrontal strengthens (z=3.58 → 5.21, p_emp=0.030 → 0.010). **β anatomy verdict unchanged**.
    - **γ_l Grassmann**: `|S(γ_l)| = 41` cells, was K*=[12,23] 13 contiguous. A3-passing region set **shifted to 7 regions**, only 3/6 overlap with the locked ledger. Retained: ctx-lh-middletemporal, ctx-lh-superiortemporal, ctx-rh-parstriangularis. **Dropped**: ctx-lh-inferiortemporal, ctx-lh-fusiform, ctx-rh-paracentral. **Added**: ctx-lh-lateraloccipital, ctx-lh-rostralmiddlefrontal, ctx-rh-medialorbitofrontal, ctx-lh-cuneus. The "temporal-cortex dominant" narrative weakens — the new network is **occipito-temporal + frontal + medial-OFC**. **The KC-era "left fusiform" claim at γ_l is fully retracted** (was carried under K*; not in S(γ_l)).
    - **δ-full Grassmann**: `|S(δ)| = 23` cells, was K*=[57,63] 7 contiguous. A3-passing region set **shrinks to 4 regions**, only 1/6 overlaps. Retained: ctx-lh-inferiorparietal. **Dropped**: ctx-lh-fusiform, ctx-rh-medialorbitofrontal, ctx-rh-caudalanteriorcingulate, **Amy**, ctx-lh-bankssts. **Added**: ctx-lh-inferiortemporal, ctx-rh-parstriangularis, ctx-lh-superiortemporal. **The "anchor anatomy" interpretation (Amy + cingulate + fusiform) is not supported under S(δ)** — replaced by a parietal-temporal-frontal network. The δ cross-probe 1.55× anchor-anatomy known-biology layer (`memory/epileptic_imcoh_universal.md`) is a separate descriptive observation and is not affected.
    - **δ-epi-X Grassmann**: `|S^epiX(δ)| = 22` cells (k = [2, 34..52, 87, 88], derived from per-k Wilcoxon on `data/audit/grassmann_epi_exclusion/per_patient_per_band_per_k.csv`), was K*=[33,39] 7 contiguous. A3-passing region set **shrinks to 3 regions**, only 2/6 overlap. Retained: ctx-lh-superiorparietal, ctx-rh-rostralmiddlefrontal. **Dropped**: ctx-lh-fusiform, ctx-lh-inferiorparietal, ctx-rh-postcentral, ctx-lh-inferiortemporal. **Added**: ctx-lh-superiorfrontal. The parietal-dominant narrative is preserved (superior parietal kept).
    - **δ full-vs-epiX dissociation strengthens under S(b)**: the two networks now share **0/3 named regions** (locked ledger had 2/6 shared via fusiform + inferiorparietal). The "mixture of two distinct phenomena" claim becomes stronger.

    α cophenet and β cophenet verdicts are unaffected (cophenet anatomy uses top-decile per-pair `|Δρ_split^coph|`, no `k`-aggregation involved).

    The Per-(band, probe) verdict-detail subsections above remain in this ledger as **historical reference under K*(b)** until cascaded edits land in `bands/01_beta.md` and the β results LaTeX. The cluster-extent CSVs at `data/audit/anatomy_<band>_grassmann{_epiX}_clusterext/cohort_summary.csv` are the **new source of truth** for any per-band anatomy citation.
