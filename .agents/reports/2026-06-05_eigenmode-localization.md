---
name: eigenmode-localization
era: IMCOH_ABS_COHORT_N10
status: current
date: 2026-06-05
companion: data/audit/eigenmode_localization/
scope_report: .agents/guides/task-persistence-investigation/2026-05-08_epi-eigenmode-localization.md
sources:
  - audit_91  eigenmode_localization                       # this analysis
  - audit_63  matched_strength_surrogate_split_baseline    # cached null eigvecs
  - audit_66  grassmann_cluster_extent                     # leading-subspace probe
  - audit_77/78/79  epi_stratified_*                       # delocalization prior
pointers:
  - .agents/reports/2026-06-01_trace-concordance-vs-blind-fc.md
  - memory/epi_stratified_trace_2026_06_05.md
  - memory/dendrogram_persistence_gate_2026_06_04.md
---

# Does the task-trace live in localized eigenmodes, and do they trap on epileptic nodes?

**Head.** Two scientifically motivated hypotheses, both **negative**, plus one
strong structural fact. The a-priori headline — that α's per-pair cophenetic
trace lives in *localized, non-leading* eigenmodes while the *leading extended*
subspace stays put (which would mechanistically explain why α has a cophenetic
trace but no leading-subspace Grassmann trace) — is **doubly falsified**:
the leading (slowest-λ) modes are the *most* localized, not the most extended,
and the modes that actually move task→rest_post are the *extended* ones, not
the localized ones (cohort Spearman(PR, displacement) = +0.66 to +0.70,
positive in 10/10 every band). Eigenmode localization does **not** explain the
α dissociation. The epileptic-trap hypothesis (Q2) is a **clean negative**: the
epileptic node set carries no more mass in localized modes than a
strength-matched null, and no more than a *random node set of its size* — fully
consistent with the established spatial-delocalization verdict. The one
positive: observed FC-Laplacian eigenmodes are 2–10× **more localized** than a
degree/strength-matched graph (participation number PR ≈ 2.4–6.5 of N ≈ 113–122
vs surrogate PR ≈ 5–40; ratio < 1 in **all 240** patient×band×phase cells), but
this is a generic property of every band and phase, not band- or trace-specific,
so it does no work for the task-trace question.

All nulls are the mandatory degree/strength-preserving matched-strength
surrogate (R = 200, eigendecomposition already cached); the observed
eigenvectors are one `eigh(L = D − W)` per cell, identical in construction, so
observed and null are directly comparable. Build:
`scripts/01_compute/audit/audit_91_eigenmode_localization.py`. Companion CSVs +
3 vector PDFs: `data/audit/eigenmode_localization/`.

---

## 1. What was asked (two connected questions)

- **Q1 (priority).** Does the trace live in **localized** (high inverse
  participation) vs **extended** modes? The a-priori most-plausible story:
  α's localized non-leading modes carry the per-pair cophenetic trace while
  its dominant extended subspace is unchanged ⇒ the leading-subspace Grassmann
  probe is null at α even though the per-pair cophenetic probe fires. This
  would *mechanistically explain the α dissociation*.
- **Q2.** Do localized modes concentrate on the **epileptic** node set
  (graph analogue of Anderson localization — the epi zone as an eigenmode
  trap)? Connects to the parallel occult-epi-marker investigation.

The 5-point critical preamble (claim / null / strongest alternative / does the
null mechanically control it + what it cannot reject / falsification) was
written into the scope doc **before** any code
(`2026-05-08_epi-eigenmode-localization.md`).

## 2. Measures (library-first, all matched-strength-nulled)

- **Participation number** `PR_k = (Σ_i v_k(i)²)² / Σ_i v_k(i)⁴` per nontrivial
  mode (small = localized; range [1, N]). New library primitive
  `lrg_eegfc.utils.metrics.spectral.participation_number`.
- **Cross-phase mode displacement** `disp_k = 1 − max_j |⟨v_k^post, v_j^pre⟩|`
  per rest_post mode (0 = unchanged, large = reorganized). New library
  primitive `spectral.subspace_displacement`.
- **Epi mode-mass** `m^E_k = Σ_{i∈E} v_k(i)²`. New library primitive
  `spectral.node_set_mode_mass`. Q2 uses the epi mass of the **most-localized
  decile** of modes (`m^E_loc`, the Anderson-trap statistic) and the single
  largest epi-mass mode (`m^E_max`). The mean over ALL modes is identically
  `|E|/N` for any orthonormal basis (rows are unit-norm) and is **vacuous** —
  flagged and not used.

## 3. Results

### 3.1 Structural localization — strong, but generic (not the answer)

Observed eigenmodes are far more localized than the strength-matched null in
**every** cell (rest_post shown; all phases identical in pattern):

| band | median PR ratio (obs/null) | n localized p<.05 |
|---|---|---|
| δ | 0.294 | 10/10 |
| θ | 0.297 | 10/10 |
| α | 0.290 | 10/10 |
| β | 0.260 | 10/10 |
| low-γ | 0.155 | 10/10 |
| high-γ | 0.109 | 10/10 |

Observed PR 2.4–6.5; surrogate PR 5–40 (N = 113–122). Ratio < 1 in **all 240**
patient×band×phase cells; lower-tail surrogate p ≤ 0.005 (floor) cohort-wide.
The matched-strength null preserves every node's strength exactly, so this is
**not** a degree artifact — it is genuine concentration of mode mass. But it is
present in all bands and all phases at similar strength (and strongest in the
γ bands, which carry no trace), so it is a generic property of weighted FC
Laplacians, **not** a task-trace signature. Figure 1.

### 3.2 Q1 — FALSIFIED (the a-priori headline)

Two independent refutations, both 10/10 in every band:

**(i) The leading modes are the MOST localized, not extended.** The 10 slowest
nontrivial modes (exactly the block the leading-subspace Grassmann probe lives
on) have lower PR than the bulk in **10/10** patients, robust over K∈{5,10,20}.
The hypothesis required leading = *extended*; the data give the opposite.

**(ii) The trace is carried by the EXTENDED modes.** Per patient,
Spearman(PR_rest_post, displacement_pre→post) is **positive in 10/10 every
band** (median +0.66 α, +0.68 β, up to +0.70 γ; range [+0.42, +0.78], never
negative). Extended modes (high PR) reorganize more across phases; localized
modes are stable. Concretely the extended-mode median displacement is ~1.3–1.5×
the localized-mode displacement (α: 0.55 vs 0.43; β: 0.56 vs 0.41).

| band | lead<bulk (K=10) | median ρ(PR,disp) | ρ>0 (extended carry) | ρ<0 (localized carry) |
|---|---|---|---|---|
| α | 10/10 | +0.660 | 10/10 | 0/10 |
| β | 10/10 | +0.681 | 10/10 | 0/10 |
| low-γ | 10/10 | +0.704 | 10/10 | 0/10 |
| high-γ | 10/10 | +0.703 | 10/10 | 0/10 |

**Within-baseline control.** The PR-displacement link is just as strong
across the two split-half rest_pre estimates (rest_pre_A→rest_pre_B, no task:
median ρ +0.67 α / +0.67 β, 10/10 each) as across the task contrast. So
"extended modes move more" is a **generic structural property** of the FC
eigenbasis (extended modes are less reproducible across *any* two FC
estimates), not a task-trace channel. Localized modes carry no differential
task-induced change — there is no localized-mode trace channel to find.

**The α Grassmann-null / cophenet-positive dissociation is therefore NOT a
localized-vs-extended mode split.** Eigenmode localization does not explain it.
The dissociation must come from elsewhere — most likely that the per-pair
cophenetic probe is sensitive to small reorganizations *across the whole
spectrum* (including the many extended modes that each move a little), while the
leading-subspace chordal Grassmann probe only sees rotation of the fixed top-k
block; at α the cross-phase movement is spread thinly over many modes rather
than concentrated in the leading block, so the per-pair probe integrates it and
the leading-subspace probe misses it. Figure 2.

### 3.3 Q2 — CLEAN NEGATIVE (no epileptic eigenmode trap)

The Anderson-trap statistic — epi mass of the most-localized modes `m^E_loc`
vs the strength-matched null — does not fire in any band (cohort Wilcoxon
p = 0.18–0.99; n trapping ≤ 3/9 patients):

| band | n (epi) | n trap_loc p<.05 | cohort Wilcoxon p_loc | n max-special vs random p<.05 | median p (max vs random) |
|---|---|---|---|---|---|
| δ | 9 | 1 | 0.715 | 1 | 0.438 |
| θ | 9 | 0 | 0.180 | 0 | 0.522 |
| α | 9 | 0 | 0.986 | 2 | 0.149 |
| β | 9 | 1 | 0.980 | 3 | 0.224 |
| low-γ | 9 | 2 | 0.951 | 3 | 0.189 |
| high-γ | 9 | 3 | 0.715 | 3 | 0.269 |

**Disambiguation of one apparent positive.** The single-largest-epi-mass mode
`m^E_max` *does* exceed the strength null (7–9/9 patients, p≈0.002 every band).
This is **not** a trap — it is the localization-extreme artifact: because the
observed modes are localized (§3.1), the maximum of *any* per-mode mass
statistic beats the extended-null maximum. The decisive control is a
**random same-size node set on the observed eigenvectors**: the epi set's
`m^E_max` is *not* exceptional among random sets of size |E| (median upper-p
0.15–0.52, only 0–3/9 patients < 0.05, last two table columns). Epi behaves
like a random node set of its size. No epileptic trapping beyond strength —
fully consistent with the established spatial-delocalization verdict
(per-patient localization null 0–1/10; cohort anatomy retracted; apparent
hotspots = electrode-shaft autocorrelation). Figure 3.

The surrogate does not control sEEG-shaft spatial autocorrelation, but because
Q2 is negative *even before* any geometric control, the negative is conclusive
(a positive would have needed whole-shaft masking before any claim).

## 4. Honest bottom line

- **Q1 (priority): NO.** Localization does not explain the α dissociation. The
  a-priori hypothesis is refuted in the strongest possible direction — the
  leading modes are the most localized and the trace rides the extended modes.
- **Q2: NO.** No epileptic eigenmode trap beyond strength; epi = random set of
  its size on the observed modes.
- **Structural fact (verified, strength-independent): YES but unhelpful.** FC
  Laplacians are strongly localized vs degree-matched graphs, generically
  across bands/phases — interesting physics, but not a task-trace or epilepsy
  signature.

This is a clean negative on both scientific questions. Nothing here changes the
β-trace headline or the §5.3 cohort verdict; it closes the eigenmode-
localization direction (Direction C of the epilepsy plan) as not productive for
either the trace mechanism or the epi marker, and removes "localized modes
carry the α trace" from the menu of candidate mechanisms for the α dissociation.

## 5. Files

- `scripts/01_compute/audit/audit_91_eigenmode_localization.py` — build script
  (5-point preamble in docstring + scope doc).
- `src/lrg_eegfc/utils/metrics/spectral.py` — new library primitives
  `participation_number`, `node_set_mode_mass`, `subspace_displacement`.
- `data/audit/eigenmode_localization/`
  - `localization_per_cell.csv` — per (patient, band, phase)
  - `q1_displacement_per_band.csv` — per-patient PR-displacement Spearman
  - `q1_cohort_per_band.csv`, `q2_epi_mass_per_band.csv` — cohort tables
  - `README.md`
  - `figures/fig_01_pr_obs_vs_null.pdf` — observed vs null PR per band
  - `figures/fig_02_q1_pr_vs_displacement.pdf` — Q1 PR-displacement
  - `figures/fig_03_q2_epi_mass.pdf` — Q2 epi mass of localized modes
- `.agents/guides/task-persistence-investigation/2026-05-08_epi-eigenmode-localization.md`
  — scope doc, now status `implemented` with outcome summary + 5-point preamble.
