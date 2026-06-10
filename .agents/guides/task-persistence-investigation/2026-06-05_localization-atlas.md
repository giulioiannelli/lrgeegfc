---
name: localization-atlas
type: scope
era: IMCOH_ABS_COHORT_N10
status: current
created: 2026-06-05
updated: 2026-06-05
result: data/audit/localization_atlas/README.md (β trace localizes to MTL+OFC, epi-independent; α diffuse; Grassmann probe null; low-γ epi-carried)
pointers:
  - data/audit/anatomy_localization_wilcoxon/README.md
  - data/audit/per_patient_localization/README.md
  - .agents/preprint/locked/ANATOMY_LEDGER.md
  - src/lrg_eegfc/utils/metrics/node_localization.py
  - src/lrg_eegfc/utils/io/regions.py
---

# Sampling-conditioned localization atlas (descriptive; multi-granularity; epi-resolved)

**Head.** A *descriptive* re-framing of trace localization that drops the cohort
permutation **gate** (a structural false-negative machine at n=10 / ≤5 patients
per region) and instead reports, for every anatomical unit at four granularities
(DK region → system → lobe → hemisphere), the **sampling-conditioned positive
rate** `J⁺_g / K_g` = (patients in whom the trace concentrates in unit `g` above
their own baseline) / (patients **implanted** in `g`), plus the effect-size
distribution. The denominator is *patients implanted there*, not the full cohort,
so Hip-β reads **4/5**, not 4/10. Run for both matched-strength-verified TRACE
probes (cophenet + cross-phase Grassmann), with and without epileptic-node
exclusion, and the depth-shaft adjacency control reported **alongside** (not used
to dissolve deep-target findings). This is exploration, not a verdict: no
pre-registered acceptance gate; compute → plot → adjudicate with the user.

## 0. Why this measure (and what was wrong before)

`diag_anatomy_localization_wilcoxon.py` (cohort) and `diag_per_patient_localization.py`
(within-patient) both made the **permutation gate the verdict**. At n=10 with a
≤5/10 per-region sampling ceiling, that gate is biased to "diffuse/retract": its
own README documents the two natural cohort statistics failing in *opposite*
directions (−log10 p mass → false negative; floor-free count → thin-region false
positive). Treating non-significance there as "no localization" is unsound. Two
concrete errors propagated into the 2026-05-30/06-01 retraction:

1. **Wrong denominator.** "3–5/10 patients = thin = marginal" uses the full cohort
   as denominator. No DK region is sampled by **>5/10** anywhere — an implant
   *ceiling*, not a per-region weakness — so a "n≥5 or it's thin" bar disqualifies
   *every* region by construction.
2. **Anchor/trace conflation.** The "Hip = marginal" line leaned on the β
   **Grassmann** quantity, which is **phase-averaged** participation (an *anchor*).
   The genuine signal is the β **cophenet TRACE** `s_ij = dD_task·dD_rest`.

This atlas fixes both: sampling-conditioned denominator, trace probes only,
descriptive-first.

### 5-point critical preamble (per `feedback_critical_null_preamble.md`)

1. **Claim.** Given a verified cohort trace (β/α; matched-strength-locked at §5.2/5.4),
   the trace concentrates in specific anatomy in the patients sampled there — and
   may do so at a coarser *system* granularity even when no single DK region recurs.
2. **Null (demoted to a column, not the gate).** Within-patient node→unit
   label-shuffle (preserves implant unit-size multiset + value distribution;
   randomizes which node carries which label). Reused verbatim from the existing
   localization probes.
3. **Strongest alternative.** (a) The per-unit positive rate is a coin-flip
   artifact (with K=5, 4/5 positives has binomial p≈0.19 under a fair coin — weak
   on its own). (b) Apparent concentration is **electrode-shaft spatial
   autocorrelation**, not anatomy (adjacent contacts measure near-identical signal).
4. **Does the design cover it — by mechanism.** (a) The rate is **never read
   alone**: it is paired with the effect-size distribution, the demoted
   permutation p, and a stated prior (a *prespecified* hippocampal hypothesis is
   not 53-way fishing). (b) The shaft control (`η²` region-vs-shaft, NMI) is
   computed and shown **side by side** — but with a `deep_target` flag: for
   deliberately-targeted deep structures (Hip, Amy, Thal) on a depth shaft the
   electrode **is** the structure, so shaft-level concentration there is the
   localization, not a confound; the shaft control discounts only *cortical-surface*
   adjacency. (User decision 2026-06-05: report both readings.)
5. **Falsification / limits.** A unit's localization claim is **descriptive and
   conditional**: it holds only among implanted patients and says nothing about
   un-implanted ones; it does **not** re-test trace existence (demeaning removes
   the base rate). A unit with high `J⁺/K` but near-zero effect sizes, or one whose
   rate collapses on epi-exclusion while its support was epileptic contacts, is
   **not** a clean anatomical localization. The atlas surfaces candidates; it does
   not lock them.

## 1. Notation

- Cohort `P`, `|P| = 10`. Bands `b ∈ {δ, θ, α, β, γ_l, γ_h}`. Trace probes
  `π ∈ {coph, grass}` (cophenet edge-trace, cross-phase Grassmann node-trace —
  **both** matched-strength-verified directions; the phase-avg Grassmann *anchor*
  is excluded by design).
- Epi mode `ε ∈ {incl, excl}` (excl = drop epileptic-labelled contacts via the
  fixed label-based `epi_keep_mask`).
- Granularity `G ∈ {region, system, lobe, hemi}`; a *unit* `g` is a value at that
  granularity (e.g. `g = Hip` at `region`, `g = MTL` at `system`).
- Per patient `p`, probe `π`, band `b`: signed contributions
  `S_p = (x_c, ν_c)_{c=1..C_p}` where `x_c ∈ ℝ` is a signed trace score and
  `ν_c ∈ {1..N_p}` is its endpoint node index.
  - **coph**: per upper-triangle pair `(i,j)`, `s_ij = dD_task·dD_rest`; each pair
    contributes `x = s_ij` to **both** endpoints `ν = i` and `ν = j`
    (`C_p = 2·\binom{N_p}{2}`).
  - **grass**: per node `i`, `s_i = (p_i^{task}−p_i^{pre})·(p_i^{post}−p_i^{pre})`
    over the band k-set; `ν = i`, `C_p = N_p`.
- `u_G : {1..N_p} → units` maps each node to its unit at granularity `G` (from
  `load_channel_regions`: `region`, derived `system`, `lobe`, `hemisphere`).

## 2. Definitions

**Per-patient demeaning** (isolates *localization* from trace *existence*):
`x̃_c = x_c − x̄_p`, `x̄_p = mean over anatomical contributions of patient p`.
A globally positive trace lifts every unit's raw mean; after demeaning only units
exceeding the patient's **own** average register as concentrated. `x̄_p` is
label-independent ⇒ a per-patient constant under the unit-label shuffle, so the
null stays calibrated. (Identical to the existing localization probes.)

**Per-(patient, unit) statistics** (anatomical contributions only;
`Wm/Unk/unknown` dropped):
- mean `m_{p,g} = mean{ x̃_c : u_G(ν_c) = g }` ∈ ℝ. `>0` ⇒ trace concentrates in
  `g` above patient baseline.
- standardized strength `z_{p,g} = (m_{p,g} − μ⁰_{p,g}) / σ⁰_{p,g}`, where
  `(μ⁰, σ⁰)` are the mean/SD of `m_{p,g}` under the within-patient label shuffle.
  Descriptive effect size; **not** thresholded.
- demoted per-patient permutation p (one-sided greater), kept as a column.

**Sampling-conditioned atlas cell** — for each `(π, b, ε, G, g)`:
- `K_g = #{ p : patient p has ≥1 anatomical contribution with u_G(ν)=g }`
  (patients *implanted* in `g`).
- `J⁺_g = #{ p sampling g : m_{p,g} > 0 }`. **Positive rate `r_g = J⁺_g / K_g`.**
- cohort effect: `M_g = median_{p sampling g} m_{p,g}`, `Z_g = median z_{p,g}`.
- `p_g` = floor-free permutation rank of `M_g` vs the shuffle null **(demoted;
  reported, never the verdict)**.
- flags: `deep_target` (g is/contains a depth-target subcortical structure:
  Hip, Amy, Thalamus), `n_pos`, `n_neg`, sign-consistency.

**Per-patient profile** — for each `(p, π, b, ε, G)`: the argmax unit
`g*_p = argmax_g m_{p,g}`, its `m`, `z`, and the **system** it belongs to. Lets us
read "different DK region per patient, same system" directly (tally how often
`g*_p` lands in the same *system* across patients even when the DK `region`
differs).

**Shaft both-readings** (reuse `η²`, `perm_p_eta2`, `nmi`, `shaft_of`): per
`(π, b, ε)`, the between-unit `η²` for `region` vs `shaft`, with NMI collinearity
— plus, for `deep_target` units, an explicit "electrode = structure" note so the
shaft control is not read as dissolving them.

## 3. Anatomical systems (the new grouping)

`system_of(region)` (added to `regions.py`, general-purpose):
- **MTL** = {Hip, Amy, ctx-{lh,rh}-entorhinal, ctx-{lh,rh}-parahippocampal} —
  medial temporal lobe. (Cohort sampling: Hip 5, entorhinal 5, Amy 3,
  parahippocampal 1 — *union* raises `K` well above any single region; the central
  test of the "same system, different region per patient" hypothesis.)
- **limbic** = MTL ∪ cingulate (all subregions) ∪ insula.
- **lat_temporal** = superior/middle/inferior temporal + bankssts + fusiform +
  temporalpole + transversetemporal.
- **OFC** = {ctx-{lh,rh}-{medial,lateral}orbitofrontal}.
- **dlPFC/PFC**, **sensorimotor**, **parietal**, **occipital** — coarse cortical
  systems (reuse `_REGION_GROUPS` where they already exist; add only the missing).
- `lobe` and `hemi` come straight from `load_channel_regions` (existing).

Systems are anatomy, defined a priori from the DK atlas — **not** data-driven, so
they introduce no selection circularity.

## 4. Properties

- **Range.** `r_g ∈ [0,1]`; `m_{p,g} ∈ ℝ`; `z_{p,g} ∈ ℝ`.
- **Sampling-conditioned by construction.** `K_g` is the denominator everywhere;
  un-implanted patients never enter a unit's rate.
- **Invariances.** `m`, `z`, `p` are demean-invariant (shift cancels). `r_g`
  invariant to monotone rescaling of `x` (sign-based).
- **Pooling monotonicity.** Coarser `G` ⇒ larger `K_g` (more support) but lower
  spatial resolution; the atlas reports all four so the resolution/support
  trade-off is explicit, never hidden.
- **What it does NOT measure (contract).**
  1. Trace *existence* — removed by demeaning; existence is the matched-strength
     rung (locked β/α). The atlas is *conditional on* a trace.
  2. Anything about patients not implanted in `g`.
  3. Anatomy *beyond electrode geometry* for deep targets — for a single depth
     shaft in Hip, "anatomy" and "shaft" are the same object (reported, not
     resolved).
  4. A significance verdict — `r_g` and `M_g` are descriptive; `p_g` is demoted.

## 5. Caveats & failure modes

| Caveat | Mitigation |
|---|---|
| `r_g = J⁺/K` is weak alone (K=5, 4/5 ⇒ binom p≈0.19) | Always paired with `M_g`, `Z_g`, demoted `p_g`, and a stated prior; never reported as a standalone verdict. |
| Cophenet endpoint double-counting (pairs → 2 contributions) ⇒ non-independent node contributions | Descriptive `m` unaffected; the permutation null shuffles node labels preserving this structure (reused verbatim). |
| Thin units (K=1–2) | Reported with `K` annotated; not gated out, but flagged low-support; systems pooling is the primary remedy. |
| Shaft autocorrelation inflates apparent regional concentration | `η²` region-vs-shaft + NMI shown side by side; `deep_target` flag distinguishes "electrode = structure" (keep) from cortical-surface adjacency (discount). |
| Hip/Amy carried by epileptic contacts (MTL is a common SOZ) | `ε ∈ {incl, excl}` split; a rate that survives `excl` is clean anatomy, one that collapses was epi-driven (ties to audit_77: α/β strengthen under exclude_epi). |
| Pat_15 (RH-only, β-anti) and Pat_10 (β-anti) drag rates | Kept in denominator (no dropout per policy); flagged in per-patient profiles; sign-consistency column surfaces their effect. |
| `system` definition is a researcher choice | Defined a priori from DK, documented here; sensitivity to MTL membership (±Amy) reported descriptively if it matters. |

## 6. Pseudocode

```
for π in {coph, grass}:
  for ε in {incl, excl}:
    atlas_rows = []; profile_rows = []
    for b in BANDS:
      per_patient = {}
      for p in COHORT:
        (x, ν, region_vec) = contributions(π, p, b)        # lib probe
        if ε == excl: (x, ν, region_vec) = drop_epi(x, ν, region_vec, p)
        x̃ = x − mean(x over anatomical contributions)
        per_patient[p] = (x̃, ν, unit_maps(p))              # unit_maps: G→ node→unit
      for G in {region, system, lobe, hemi}:
        units = ∪_p set(unit_maps(p)[G])  minus NON_ANATOMICAL
        for g in units:
          samplers = [p : g in unit_maps(p)[G] image]
          K = len(samplers)
          m = [ mean(x̃_p[ unit_maps(p)[G][ν]==g ]) for p in samplers ]
          z = [ standardize(m_p vs within-patient shuffle null) for p in samplers ]
          Jpos = count(m_p > 0)
          atlas_rows += { π,b,ε,G,g, K, Jpos, rate=Jpos/K,
                          median_m=median(m), median_z=median(z),
                          perm_p=floorfree_rank(median(m), shuffle_null),   # demoted
                          deep_target=is_deep(g), n_pos, n_neg }
        for p in samplers-of-all:                            # per-patient profile
          g* = argmax_g m_{p,g}
          profile_rows += { p,b,π,ε,G, top_unit=g*, m=m_{p,g*}, z=z_{p,g*},
                            system=system_of(g*) }
    write atlas_rows, profile_rows
    # shaft both-readings (reuse existing η² machinery), per b:
    for b in BANDS: write eta2_region, eta2_shaft, nmi, deep_target_units
```

`drop_epi` uses the fixed label-based `epi_keep_mask` (audit_71 pattern).
`floorfree_rank` and the shuffle are the existing `perm_p` logic.

## 7. Visualization spec

1. **Atlas heatmap** (headline). Rows = units (grouped: DK regions, then systems,
   then lobes, then hemispheres, with separators); columns = 6 bands; one figure
   per `(π, ε)`. Cell **fill = `r_g`** (positive rate) on a saturated sequential
   cmap (turbo / cividis — **never** a near-white-interior cmap per
   `feedback_no_near_white_cmaps`); cell **annotation = `J⁺/K`**. Cell **border
   thickness ∝ |Z_g|** (effect size) so a high rate with tiny effect looks visibly
   weak. `deep_target` units marked with a left-margin tick. Reading rule: *a
   bright, thick-bordered cell with a high J/K and large K is a localization
   candidate; a bright thin-bordered cell with K=1 is noise.*
2. **Per-patient hotspot dot-plot.** For a chosen band (β first), x = patients,
   y = top-unit `m` with the unit label; color = the unit's **system**. Reading
   rule: *same color across patients = same-system localization even when labels
   differ* (the user's hypothesis).
3. **Epi split panel.** Hip/MTL `r_g` and `M_g` at `incl` vs `excl`, per band —
   bars side by side. Reading rule: *bar holds or rises under excl ⇒ clean / epi-
   independent; collapses ⇒ epi-carried.*
4. **Shaft both-readings.** Existing `η²` region-vs-shaft scatter + NMI, annotated
   with which units are `deep_target`. Reading rule: *region η² ≈ shaft η² ⇒
   cortical-surface concentration is adjacency; but deep-target units are read as
   electrode = structure regardless.*

All figures: `use_lrg_style()`, PDF-only, full vector, transparent, no suptitle,
no in-axes stat text (push to README / stdout).

## 8. Connection to prior tools

| Prior tool | Relation |
|---|---|
| `diag_anatomy_localization_wilcoxon` (cohort gate) | **Reinterprets, does not replace.** Same signed `s_ij`/`s_i`, same demeaning, same shuffle null; demotes the gate to a column and adds sampling-conditioned rates + systems + epi split. |
| `diag_per_patient_localization` (within-patient MAX) | Complements: the MAX-vs-own-shuffle is conservative and gave 0/10; the atlas reads the *descriptive* per-patient hotspot + system membership instead of pass/fail. |
| `diag_per_patient_concentration_control` (η² shaft) | **Reuses** `η²/perm_p_eta2/nmi/shaft_of` (to be promoted to lib); keeps the shaft reading but adds the `deep_target` distinction and shows it side-by-side rather than as a dissolver. |
| audit_77/78/79 (epi-stratified trace) | The `ε` split connects here: audit_77 found α/β trace *strengthen* under exclude_epi carried by cross+nonepi pairs — the atlas asks whether the *hippocampal* concentration specifically is epi-carried or epi-independent. |
| ANATOMY_LEDGER (locked) | The atlas feeds the 2026-06-05 verdict correction; multi-region lists stay retracted, Hip-β reinstated descriptively. |

## 9. Implementation plan

- **Script:** `scripts/01_compute/audit/audit_81_localization_atlas.py`
  (5-point preamble in docstring; `--probe`, `--epi-mode`, `--band`, `--n-perm`,
  `--granularity`). Reuses `cophenet_trace_contributions` /
  `grassmann_trace_contributions` and `load_channel_regions`.
- **Library promotion** (per `coding-rules.md`, ≥2 callers): move `shaft_of`,
  `eta2`, `perm_p_eta2`, `nmi` from `diag_per_patient_concentration_control.py`
  into `lrg_eegfc.utils.metrics.node_localization` (general names; the control
  script and audit_81 both import them). Add `system_of(region)` +
  `ANATOMICAL_SYSTEMS` and a `system` column to `load_channel_regions` in
  `lrg_eegfc.utils.io.regions` (general anatomy, no manuscript scope).
- **Figures:** `scripts/02_preprint/preprint_NN_localization_atlas.py` (after the
  CSVs land and we've looked at them — do not pre-build).
- **Outputs:** `data/audit/localization_atlas/`
  - `atlas_{π}_{ε}.csv`, `per_patient_{π}_{ε}.csv`,
    `shaft_both_readings_{π}_{ε}.csv`, `summary.csv`, `README.md`.
- **Env:** `lapbrain`. Read-only on caches (no FC/LRG recompute).

## 10. Open questions (deferred to post-compute adjudication)

- MTL membership: include Amy? include Thalamus? (Report ±Amy sensitivity if it
  changes the MTL rate.)
- Effect-size standardization: `z` vs the within-patient shuffle SD is the default;
  consider also `m / IQR_p` as a robust alternative if shuffle SD is unstable at
  small K.
- Whether to add a `consensus` system-recurrence statistic (how often the
  per-patient argmax lands in the *same* system) as a headline number — decide
  after seeing the per-patient profiles.
- No acceptance gate is pre-registered (per `feedback_no_pre_registered_acceptance`):
  we compute, plot, and decide signal vs noise with the user.
