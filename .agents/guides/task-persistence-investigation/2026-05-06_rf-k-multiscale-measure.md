---
name: 2026-05-06_rf-k-multiscale-measure
type: guide
era: IMCOH_ABS x COHORT_N10
status: parked
created: 2026-05-06
updated: 2026-05-07
pointers:
  - .agents/reports/archive/2026-05/2026-05-06_section-5-critical-review.md
  - .agents/reports/2026-05-07_rf-clade-persistence-cohort-verdict.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-anchor-modules.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-rearrangement-modules.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-reset-modules.md
  - data/reports/section_5_lrg_trace/14_rf_clade_persistence/
  - data/audit/kc_trace_network_view/
  - data/audit/kc_reset_module_view/
  - data/audit/kc_rearrangement_module_view/
  - data/audit/kc_anchor_module_view/
  - scripts/01_compute/audit/audit_48_rf_k_clade_persistence.py
  - scripts/01_compute/audit/audit_48b_rf_k_soft_jaccard.py
  - scripts/01_compute/audit/audit_48c_rf_real_vs_null_scatter.py
  - scripts/01_compute/audit/audit_48d_trace_clades_on_dendrograms.py
---

> **PARKED 2026-05-07.** This measure (hard RF, soft RF, per-clade strict
> 4-mode taxonomy) is set aside as a Section-5 candidate. We document
> the definitions and findings on disk in case it becomes useful later
> as a strict-identity sensitivity check, but the probe is **not in the
> manuscript** and not scheduled for further development.
>
> Reason: RF reads clade leafsets via Jaccard, which is symmetric and
> size-sensitive. Neither hard threshold counting nor soft mean-best-J
> cleanly captures graded clade reorganization (a 10-leaf task clade
> nested inside a 15-leaf rs_pre subtree gives J=0.67 — looks like
> presence, but the leaves ARE actually clustered together in a larger
> module; vs J=0.50 against a 5-leaf subset — looks like partial
> presence, but the leaves are partially scattered). The metric shape
> doesn't match the question.
>
> KC λ=0 (audit_36 cohort scalar + audit_46 within-baseline-null
> controls + audit_47-50 module-view 4-mode catalogs) is the
> load-bearing measure for Section 5.4. KC operates on per-leaf-pair
> MRCA depths (integer counts of internal nodes), is naturally graded,
> and finds β at 10/10 vs null with p=0.000977 cleanly. The audit_47-50
> module catalogs use Jaccard with size-matched + containment-factor
> gates calibrated to the issue and surface clade-with-drift modules
> visually. RF doesn't add information beyond what KC + audit_47-50
> already establish more cleanly.
>
> Full reasoning and parked-state numbers (16/13/32/48 trace/reset/
> persist/rearrange across 60 cells) at
> `.agents/reports/2026-05-07_rf-clade-persistence-cohort-verdict.md`.

# RF(k) — clade persistence across phases as a multiscale topology probe

## What is missing from section 5

Section 5 currently carries five LRG probe families. KC (Kendall-Colijn at lambda in {0, 0.5, 1}) reads the leaf-pair common-ancestor depth basis as a single scalar per tree pair. D-rank reads the upper-triangular entries of D-hat(tau). VI(k) reads cluster-assignment label-space mismatch at one cut depth at a time. Grassmann reads the cumulative principal-angle rotation between leading eigenmode subspaces at a chosen spectral cutoff. CTM reads per-pair correlation between phase-shift maps in distance space.

None of these probes asks the cleanest mechanistic question we have for the trace: do entire branches of the dendrogram (clades, defined as leaf sets descending from one internal node) persist across phases? A clade that exists in the rsPre tree at some height and exists in the taskt tree at some height with the same leaf composition represents a coherent subnetwork that the LRG identifies as a communication unit at both phases. If the same clade also persists into the rsPost tree, it is a preserved communication branch across the rest-task-rest sequence: a "trace" in the most direct mechanistic sense.

The new measure introduced here is RF(k), where k indexes the cut height in the dendrogram. It computes clade persistence across phases at every level of the hierarchy and aggregates into a single per-(patient, band) trace scalar that supports cohort Wilcoxon testing and joins the multiscale companion probes of section 5.4.

## Definition

For a single (patient, band) cell, let T_phi denote the dendrogram at phase phi in {rsPre, taskt, rsPost}. Each T_phi has N-1 internal nodes; each internal node defines a clade by its leaf descendant set. Let C_phi(k) be the set of clades obtained by cutting T_phi at the k-th merge level (i.e. above merge index k from the root, the tree partitions into k+1 clades; the cut level k indexes the partition coarseness, similar to VI(k)).

For two phases phi_A and phi_B, define the clade persistence rate at cut level k:

P(phi_A, phi_B; k) = (1 / |C_{phi_A}(k)|) sum over C in C_{phi_A}(k) of indicator[max over C' in C_{phi_B}(k) of Jaccard(C, C') >= theta]

with theta the matching threshold. Two thresholds are reported: theta_strict = 0.85 (almost-identical leaf sets) and theta_partial = 0.70 (mostly-identical leaf sets, tolerant of a few leaf flips). Both produce one persistence rate per (phi_A, phi_B, k, threshold) cell.

The trace scalar at cut level k is

T_RF(k; p, b) = P(taskt, rsPost; k) - P(rsPre, taskt; k)

with T_RF(k) > 0 the trace direction (more clades persist from taskt into rsPost than from rsPre into taskt at this cut level). Note the sign convention is opposite to KC, VI, and Grassmann triangle scalars: for those, trace direction is negative. We invert the sign in plots for visual consistency with the other section 5 probes (positive on the y-axis = trace).

The cohort triangle scalar at the band level, integrating across cut levels, is

T_RF^cohort(p, b) = mean over k of T_RF(k; p, b)

restricted to a meaningful k range that excludes (i) the trivial root cut (k=1, where both trees have one clade containing all leaves and persistence is always 100%) and (ii) the singleton-dominated fine cuts (k > N - 5, where most clades are single leaves and Jaccard becomes degenerate). The acceptable range is approximately k in {2, ..., N - 5}.

For cohort statistics, use the per-patient T_RF^cohort(p, b) as the patient observation; n=10 patients per band; one-sided paired Wilcoxon T_RF^cohort > 0 against zero (or against a within-baseline drift-floor null if available; see Controls below).

## Why this probe earns its place in section 5

The probe reads strictly topology: leaf set composition of clades, not merge heights. It is therefore a topology-only probe like KC lambda=0 but with two key differences. First, it operates at every cut level of the dendrogram simultaneously rather than aggregating over all leaf-pair common-ancestor depths into one L2 scalar. This makes it multiscale in the same way VI(k) is multiscale, but at the clade level rather than the partition-assignment level. Second, the comparison rule is set-based with a tolerance threshold, which is more interpretable than KC's depth-vector L2 norm: a clade either persists across phases (above threshold) or it does not.

The mechanistic interpretation is direct: persistent clades are subnetworks identified as coherent communication units at both phases. Trace direction means more such subnetworks survive the task-to-rest transition than survived the rest-to-task transition. This is the cleanest "trace as preserved communication branch" reading the section can produce.

## How RF(k) relates to existing probes

RF(k) is correlated with but not redundant with KC lambda=0, VI(k), and CTM:

- KC lambda=0 measures graded common-ancestor depth shifts at the leaf-pair level. RF(k) measures binary clade persistence at the clade-set level. KC penalizes a one-leaf shift in a clade by changing one pair's depth slightly; RF(k) under threshold theta penalizes the same shift not at all if Jaccard remains above threshold and fully if it drops below.
- VI(k) measures cluster-assignment label-space mismatch given a single cut. RF(k) measures clade leaf-set mismatch across cuts. VI(k) is sensitive to leaf reassignments between clusters; RF(k) is sensitive to clade dissolution or formation.
- CTM measures per-pair correlation in distance shifts under three controls. RF(k) operates one geometric layer up, on the dendrogram derived from those distances.

If RF(k) converges with VI(k), Grassmann, and CTM on alpha/beta/low-gamma as the three trace bands, the convergence strengthens the section 5.4 multiscale companion narrative with a fourth independent probe that asks the cleanest mechanistic question.

## Controls

Two control regimes apply, in parallel to the CTM control structure:

1. Pat_03 dropout. Recompute T_RF^cohort restricting to nine patients. Report whether the band-level Wilcoxon survives.
2. Within-baseline-null. Build a null T_RF^null using rsPre split halves only: clade persistence between rsPre_A and rsPre_B serves as the drift floor; compare T_RF^real - T_RF^null cohort-wide via paired Wilcoxon. This is the analog of CTM rho_drift.

If both controls survive at uncorrected p < 0.05 with within-probe BH-FDR at m=6, RF(k) becomes a controlled cohort claim alongside CTM. If only the absolute Wilcoxon survives, RF(k) sits in section 5.4 as a multiscale companion at the same level as VI(k) and Grassmann.

## Output structure

The measure produces:

- `data/reports/section_5_lrg_trace/14_rf_clade_persistence/tables/Td_per_patient_per_band.csv`: rows = (patient, band, threshold), columns = (T_RF_cohort, T_RF_per_k as compact JSON, k_range_used, n_clades_total).
- `data/reports/section_5_lrg_trace/14_rf_clade_persistence/tables/cohort_summary.csv`: rows = (band, threshold), columns = (T_RF_median, n_trace, wilcoxon_p_oneside, BH_q_within_m6, pat03_dropout_p).
- `data/reports/section_5_lrg_trace/14_rf_clade_persistence/figures/rf_k_cohort_heatmap.pdf`: cohort-level heatmap n_trace(k)/10 across (band, k) cells, analogous to the VI(k) cohort heatmap.

The per-clade matching results (which clades persist between which phase pairs) are saved as a separate output that drives the visualization figure (deliverable 4 below).

## Pipeline placement

RF(k) is added to section 5.4 as the third multiscale companion probe alongside VI(k) and Grassmann. The section 5.4 prose currently reads two probes; it expands to three. The other section 5 subsections (5.1, 5.2, 5.3, 5.5, 5.6) are not affected.

## Verification expectations

Before the measure can be cited in the manuscript, the following must hold:

1. The Jaccard threshold sensitivity is reported at theta = 0.70 and theta = 0.85; both columns appear in the cohort table; the band ordering should be similar at the two thresholds.
2. The cohort persistence rates at both thresholds should be intermediate (neither saturating at 100% nor collapsing to 0%) at the meaningful k range. If most clades persist across phases at all thresholds, the measure has no resolution; if none do, the threshold is too strict.
3. The cohort verdict at alpha, beta, low-gamma should be in the trace direction (positive T_RF^cohort) under at least one of the two thresholds, with a one-sided Wilcoxon p < 0.10 if the convergence with CTM/VI(k)/Grassmann is genuine. If the verdict is not in the trace direction, this is itself a finding: it means clade persistence is not the right level of abstraction for what KC/VI/CTM are detecting.
4. The disagreement-bearing band (where RF(k) splits from the multiscale consensus) should be identified explicitly and reported.

## What this is not

- Not a replacement for KC. KC produces a single graded scalar per tree pair; RF(k) produces a stack of binary persistence rates indexed by k.
- Not a replacement for the cohort-level CTM controlled claim. RF(k) joins the section 5.4 companion view; the controlled cohort claim of section 5 remains CTM in section 5.3.
- Not a localization measure. RF(k) is band-level. The localization (per-leaf catalog, anatomical enrichment) lives in section 5.5.

---

# Soft variant — primary headline (added 2026-05-07)

## Definition

Replace the threshold-count rule with a continuous mean-best-Jaccard:
```
P_soft(phi_A, phi_B; k) = (1 / |C^phi_A(k)|) *
                          sum_{C in C^phi_A(k)} max_{C' in C^phi_B(k)} J(C, C')
```
Range: P_soft ∈ [0, 1]. No threshold theta. Triangle and cohort scalar are
identical in form to the hard version:
```
T_RFsoft(k; p, b) = P_soft(taskt, rsPost; k) - P_soft(rsPre, taskt; k)
T_RFsoft^cohort(p, b) = mean over k in [2, N - 5] of T_RFsoft(k; p, b)
```
Sign convention unchanged: T_RFsoft > 0 = trace direction.

## Why the soft variant is preferred

The hard threshold (J ≥ theta) systematically rejects clades with graded
reorganization in [0.5, 0.85]. β's task-induced reorganization lives
exactly there — clades reshuffle ~20–50% of their leaves between phases
without dissolving below J=0.5 and without staying above J=0.7. Hard RF
counts those zero; soft RF counts them at their actual best-Jaccard value.

Empirically (see findings below): hard RF has β n_above_null=7/10 (p=0.053,
q_BH(m=6)=0.10, fails control). Soft RF has β n_above_null=9/10 (p=0.0068,
q_BH(m=6)=0.021, passes control). α is controlled under BOTH variants.

# Empirical findings (2026-05-07)

## Cohort verdict (n=10, |ImCoh|_abs, τ=1/λ_max)

Two-tailed within-baseline-null Wilcoxon (paired one-sided, real > null);
within-probe BH-FDR at m=6:

| band | hard θ=0.70 n_above_null | hard p_null | hard q_BH | **soft n_above_null** | **soft p_null** | **soft q_BH** | soft p_drop_null |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| δ | 6/10 | 0.246 | 0.25 | 6/10 | 0.246 | 0.25 | 0.41 |
| θ | 7/10 (anti) | 0.032 | 0.10 | 8/10 (anti) | 0.019 | 0.037 | 0.037 |
| α | **9/10** | **0.0049** | **0.029** | **9/10** | **0.0020** | **0.012** | **0.0039** |
| **β** | 7/10 | 0.053 ✗ | 0.105 ✗ | **9/10** | **0.0068** | **0.021** | **0.014** |
| low_γ | 7/10 | 0.080 | 0.12 | 7/10 | 0.042 | 0.063 | 0.064 |
| high_γ | 7/10 | 0.138 | 0.17 | 7/10 | 0.116 | 0.14 | 0.21 |

**Controlled trace bands under soft RF**: α + β (BH-FDR survives at m=6,
Pat_03-dropout robust). low_γ uncorrected significant. δ, high_γ null.
θ "above-null" pass is in the ANTI direction (T_RFsoft median negative)
— drift-floor effect, not a trace claim.

## Per-clade dendrogram visualization (4-mode taxonomy)

`audit_48d` enumerates clades anchored at the appropriate phase per mode,
checking presence/absence via bestJ = max Jaccard against ALL internal-node
leafsets of the target tree (NOT same-cut-level only — this is critical
to distinguish PERSIST from RESET). Strict gates: presence ≥ 0.65,
absence ≤ 0.30. Top-3 disjoint clades per (p, b, mode) rendered on the 4
phase dendrograms with distinct colors.

Four mutually exclusive modes:

| mode | anchor | predicate | visual |
|:--|:--|:--|:--|
| **trace** | task_test | bestJ(C,post)≥0.65 ∧ bestJ(C,pre)≤0.30 | scatter→cluster→persist |
| **reset** | rs_pre | bestJ(C,post)≥0.65 ∧ bestJ(C,task)≤0.30 | cluster→scatter→cluster |
| **persist** (= anchor) | rs_pre | bestJ(C,task)≥0.65 ∧ bestJ(C,post)≥0.65 | cluster across all phases |
| **rearrange** | rs_post | bestJ(C,pre)≤0.30 ∧ bestJ(C,task)≤0.30 | scatter→scatter→cluster |

Cohort coverage (cells with ≥ 1 matching clade / 10):

| band | trace | reset | persist | rearrange |
|:--|:--:|:--:|:--:|:--:|
| δ | 3 | 1 | 5 | **10** |
| θ | 3 | 3 | 4 | **9** |
| α | 3 | 3 | 5 | **9** |
| **β** | 3 | 4 | **9** | 7 |
| low_γ | 2 | 2 | 7 | 8 |
| high_γ | 2 | 0 | 2 | 5 |
| total / 60 | **16** | **13** | **32** | **48** |

**Key reframe**: β is dominantly PERSIST (9/10), not strict TRACE (3/10).
The cohort soft-RF result (β at 9/10 above null) is asymmetric drift
WITHIN persistent modules, not strict task-induced new modules. REARRANGE
is near-universal (rs_post emergent structure exists in 8-10/10 patients
in 5/6 bands). Strict TRACE is rare (16/60 total, no band > 3/10) — the
canonical 2026-04-25 strict-J=0.9 cohort-null is now contextualized as
"strict emergent-then-persist clades are sparse across the cohort, but
cohort-level graded asymmetry exists at α + β + low_γ".

## Convergence with the substrate and other LRG probes

| probe | α | β | low_γ |
|:--|:--:|:--:|:--:|
| substrate d_S (raw FC) | 8/10 trace | 7/10 trace | 7/10 trace |
| KC λ=0 (per-pair MRCA) | passes BH | **passes BH (10/10 vs null)** | passes BH |
| soft RF (per-clade graded) | **9/10 above null, q=0.012** | **9/10 above null, q=0.021** | 7/10 above null, q=0.063 |
| **strict per-clade trace** (4-mode, mutually exclusive) | 3/10 cells | **3/10 cells** | 2/10 cells |
| **strict per-clade persist** (= anchor) | 5/10 cells | **9/10 cells** | 7/10 cells |

α is convergent on graded probes (substrate, KC, soft RF) but rare at
the strict per-clade trace level (3/10). β is the same shape: graded
trace cohort claim is real (soft RF, KC), but at the strict per-clade
level it is dominantly PERSIST (9/10). The β cohort claim should be
worded as **"asymmetric drift within persistent β modules"** rather
than "task-induced β modules persist into rest." The mechanistic
narrative under strict gates is PERSIST, not TRACE.

## Methodological lesson

Threshold-based clade-set probes (RF at θ ≥ 0.7, MRL at J=0.9, Trace-Modules
at strict-J) systematically miss graded clade reorganization. The 2026-04-25
MRL cohort-null and the 2026-04-25 Trace-Modules strict-J=0.9 cohort-null
are now reinterpreted: they did not measure "no clade-level reorganization
exists", they measured "no clade-level reorganization crosses our chosen
hard threshold." The continuous soft variant resolves both — and the
per-clade dendrogram visualization at strict gates (still produces 29/60
cells with ≥ 1 trace clade across the cohort, with β leading at 8/10)
shows the visual narrative the hard cohort scalar can't.

The remaining open question for the manuscript: is α + β the headline,
or just β (since β has the cleanest per-pair KC + per-clade soft RF +
per-clade dendrogram convergence)?
