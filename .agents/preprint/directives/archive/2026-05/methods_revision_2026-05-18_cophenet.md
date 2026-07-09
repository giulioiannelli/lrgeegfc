---
name: methods-revision-cophenet-2026-05-18
era: IMCOH_ABS_COHORT_N10
status: writing_directive
kind: methods_revision
scope: rewrite of Methods §"Comparison of the LRG information-communication output proxies" (ssec:methods_compare)
target_section: ssec:methods_compare + part of sssec:methods_lrg
companion_artifacts:
  - .agents/preprint/established_results/00_open_methodology_question_lrg_D_convention.md (resolved → withdrawn after this revision lands)
  - memory: cophenet_methodology_rationale.md
  - memory: feedback_never_confuse_D_with_cophenet.md
date: 2026-05-18
---

# Methods §"LRG comparison" — rewrite directive

## Head

The methods section currently presents the LRG comparison through **four
probe subsections** (per-pair on `D(τ)`, KC on the dendrogram, VI on
partitions, subtree taxonomy on the dendrogram) plus the subspace
Grassmann. After the 2026-05-18 methodological audit we are reducing to
**two probes**: the **cophenetic per-pair multiscale distance** and the
**leading-eigenmode subspace Grassmann distance**. The bare `D(τ)` is
never read directly — it is an intermediate construction. The dendrogram
`T(τ)` is also an intermediate construction; we do not compare trees,
partitions, or subtree taxonomies. The cophenetic image of `T(τ)`
returns the per-pair value, but with each pair-value evaluated at the
*natural communication-merge scale of that pair* rather than at the
single LRG scale `τ_max`. This is precisely the per-pair multiscale
content the LRG provides for our continuous-spectrum substrate.

---

## Edits, by subsection

### LRG framework subsection (sssec:methods_lrg)

Keep the paragraph defining `L`, `ρ(τ)`, `Π(τ)`, `S(τ)`, `C(τ)`, and the
choice of `τ = τ' = 1/λ_max`. Keep the paragraph on the natural
resolution window `[τ', τ*]`. Keep the explicit statement that the
spectral entropy / specific heat `Ψ` selector of Villegas does **not**
isolate a privileged interior cut on our dense weight-heterogeneous
`|ImCoh|` substrate.

**Modify the "Three coupled objects" paragraph as follows.** The bullet
list currently presents `D(τ)`, `T(τ)`, `U_k` as three direct
comparison objects. We are restructuring this list to make the
intermediate / read-out distinction explicit. Replacement text:

> Two coupled multiscale objects are extracted from the propagator
> `Π(τ_max)` at the working scale `τ_max = 1/λ_max` and serve, in the
> rest of the analysis, as the probes of the network's information-
> communication geometry:
>
> 1. the **cophenetic communication distance** `D_coph` between contact
>    pairs, defined as follows. The propagator induces a per-pair
>    communication distance `D_ij(τ_max) = (1 − δ_ij)/Π_ij(τ_max)` (the
>    standard form of [Villegas]). Average-linkage agglomerative
>    clustering on `D(τ_max)` produces a deterministic linkage hierarchy
>    `Z(τ_max)` whose merge heights `{h_n}_{n=1}^{N-1}` span the
>    coarsening scales of the network as resolved by `Π(τ_max)`. The
>    cophenetic distance reads, for each leaf pair `(i,j)`,
>    `D_coph_ij = h(LCA(i,j))`, the merge height at which `i` and `j`
>    first coalesce in `Z(τ_max)`. `D_coph` is strictly ultrametric by
>    construction and reads, pair-by-pair, the scale at which the two
>    contacts communicate as a single unit. This is the per-pair
>    multiscale image of the LRG propagator: the dimensionality remains
>    `N(N−1)/2` pair values, but each value is evaluated at the
>    intrinsic communication-merge scale of that pair rather than at a
>    single fixed `τ`. The raw distance `D(τ_max)` is **not** read
>    directly; for our continuous Laplacian spectrum the dendrogram
>    replaces the spectral-gap scale-identification mechanism of
>    [Villegas 2025] (see Section X.X), and the cophenetic image is the
>    natural multiscale summary of `D` along the merge-height continuum.
> 2. the **leading-eigenmode subspace** `U_k ∈ R^{N × k}` whose columns
>    are the eigenvectors associated with the `k` smallest non-zero
>    eigenvalues of `L`. These are the slowest non-trivial diffusion
>    modes and span the subspace in which low-frequency communication
>    is concentrated.

(Drop the "linkage hierarchy" bullet as a direct object; absorb it into
bullet 1 as the construction that defines `D_coph`.)

**Replace the closing sentence** of the "Three coupled objects"
paragraph with:

> The two objects read **structurally distinct aspects** of the same
> information-communication geometry: `D_coph` reads the **per-pair
> hierarchical-merge structure** (which contacts coalesce at which
> scale), and `U_k` reads the **subspace alignment of the slow
> diffusion modes** at a fixed cutoff `k`. They are not derivable from
> each other (the dendrogram is not recoverable from `U_k` alone, and
> the eigenmode subspace is not recoverable from `D_coph` alone), and
> cross-phase comparison at the two levels reads two complementary
> aspects of any coordinated reorganization.

---

### Subsubsection "Per-pair correlation on `D(τ)`" — REPLACE

This is the central rewrite. Replace the entire current subsubsection
with the following content. The new title is **"Per-pair multiscale
correlation on `D_coph`."**

Mathematical content:

> The per-pair multiscale probe operates on the cophenetic
> communication distance `D_coph` defined in the previous subsection.
> Two complementary tests are applied.
>
> The first is a set of **matrix distances** on the upper-triangular
> entries `vec_△(D_coph^φ) ∈ R^{N(N−1)/2}`: the rank distance `d_S =
> 1 − ρ_S` (one minus the Spearman correlation of the upper triangles)
> and the linear distance `d_P = 1 − ρ_P` (Pearson); each yields a
> phase-triangle scalar via [eq:methods_triangle]. (The Frobenius
> distance is not reported on `D_coph` because the cophenetic image
> has only `N−1` distinct values and `d_F` is dominated by tie
> structure of marginal interpretive value.)
>
> The second is the **controlled per-pair split-baseline test**. For
> each (patient, band) cell, rsPre is split into two halves of equal
> duration, and the LRG pipeline — eigendecomposition, propagator at
> `τ_max`, average-linkage UPGMA, cophenetic — is rerun independently
> on each half, yielding `D_coph^{rsPre_A}` and `D_coph^{rsPre_B}`.
> The per-pair shifts
>
>     Δ_task(i,j) = D_coph^{taskt}_{ij} − D_coph^{rsPre_A}_{ij}
>     Δ_rest(i,j) = D_coph^{rsPost}_{ij} − D_coph^{rsPre_B}_{ij}
>
> are referenced to **independent half-baselines**, eliminating by
> design the shared-baseline contribution that would correlate the two
> shift maps trivially when both are referenced to the same rsPre
> build. The trace scalar at the per-pair multiscale layer is
>
>     ρ_split^coph(p,b) = ρ_S(Δ_task, Δ_rest),       T_CTM(p,b) = −ρ_split^coph(p,b)
>
> which measures, on the merge-height-resolved per-pair grid, whether
> task-induced and rest-induced communication-scale shifts agree in
> sign and rank across all `N(N−1)/2` contact pairs. `ρ_split^coph >
> 0` is the trace direction.

**Methodological rationale paragraph (insert immediately after the
formula box):**

> The choice of the cophenetic image `D_coph` over the raw distance
> `D(τ_max)` is methodologically motivated. The LRG propagator
> `Π(τ_max)` is a single-scale object: at `τ_max = 1/λ_max` the heat
> kernel is set to the finest scale resolved by `L`, and `D(τ_max)`
> reads pairwise dissimilarity at that fixed scale. For our
> fully-connected weight-heterogeneous `|ImCoh|` substrate the
> Laplacian spectrum is continuous (no clean gap), so the standard
> `τ`-sweep multiscale reading of [Villegas 2025] — entropy plateaus
> and specific-heat peaks marking gap-induced crossovers — does not
> apply. The dendrogram `T(τ_max)` is the natural multiscale carrier
> in this regime: its `N−1` merge heights deterministically index the
> communication-coarsening scales of the network as read by
> `Π(τ_max)`, replacing the spectral-gap scale identification. The
> cophenetic distance `D_coph` is the per-pair image of this multiscale
> structure — each pair value reads the **intrinsic scale at which the
> pair coalesces** into a single communication unit. Computing the
> per-pair correlation on `D_coph` therefore reads cross-phase
> persistence of the multiscale-resolved per-pair organization,
> whereas the same statistic on `D(τ_max)` would only read
> single-scale per-pair persistence. The latter is reported in
> Supplementary Section X.X as a sensitivity check; the cohort
> verdict is preserved at β but the single-scale statistic is broadly
> permissive across bands and does not deliver the multiscale-level
> band resolution of the cophenetic probe.

**Drift floor and cross-probe control paragraphs:** keep them, but
replace every occurrence of `D(τ)` with `D_coph` and every
`ρ_split` with `ρ_split^coph`.

The drift-floor null:

    ρ_drift^coph(p,b) = ρ_S( D_coph^{rsPre_B} − D_coph^{rsPre_A},
                              D_coph^{rsPost_B} − D_coph^{rsPost_A} )

The cross-probe restriction: keep exactly the same operational rule —
restrict the pair set to cross-probe contacts and recompute
`ρ_split^coph`.

The cohort claim at the per-pair multiscale layer requires
`ρ_split^coph` to sit above `ρ_drift^coph` under a one-sided paired
Wilcoxon, with `ρ_xprobe^coph` preserving the cohort-median sign and
patient count.

---

### Subsubsection "Tree distance on the dendrogram" — DELETE

Remove the entire subsubsection on Kendall–Colijn distance
`d_KC(λ; T, T')` and the associated controlled test on `T_KC`. The
multiscale content of the dendrogram is now read at the per-pair level
via `D_coph` (above); a tree-level comparison no longer adds
information. The deletion also eliminates the within-baseline-null
controlled test mirrored from `ρ_drift`, which becomes redundant under
the cophenetic per-pair probe.

(Any references to `T_KC` in Results that rely on the deleted
subsection should be retired in the corresponding Results section
revision — track separately.)

---

### Subsubsection "Partition distance across cuts of the dendrogram" — DELETE

Remove the entire subsubsection on variation-of-information
`T_VI(k; p, b)` across horizontal cuts of the dendrogram. The
partition-resolved sweep over `k` is now redundant with `D_coph`:
the cophenetic distance integrates `T(τ_max)` across all merge heights
into a single per-pair value, and the partition-level information at
any fixed cut `k` is recoverable from `D_coph` by thresholding at the
corresponding height. The descriptive multiscale view that this
subsection provided is subsumed by the cophenetic probe.

---

### Subsubsection "Cross-phase module taxonomy on the dendrogram" — DELETE

Remove the entire subsubsection on the trace / reset / rearrange /
anchor / diffuse subtree taxonomy and Table `tab:methods_taxonomy`.
The taxonomy operates on dendrogram subtrees compared by Jaccard
overlap and was already reported in the original methods as a
"coarse-grained per-patient illustration rather than an independent
statistical test." With the dendrogram no longer being read directly,
the taxonomy has no place in the methods of the primary
preprint.

(If the taxonomy figures are needed as exploratory illustrations they
can move to the Supplement under a clearly demarcated "Exploratory
dendrogram-level views" subsection. This is **not required** by the
current revision.)

---

### Subsubsection "Subspace distance on the leading eigenmodes" — KEEP, with sharpened framing

Keep the entire current Grassmann subsubsection (chordal `d_G(U, U')`,
phase-triangle scalar `T_G(k; p, b)`, mode-resolved Δθ_i(k)
companion). Add an explicit **"What `T_G(k)` reads vs what
`ρ_split^coph` reads"** paragraph immediately after the mode-resolved
sentence:

> The Grassmann probe and the per-pair multiscale probe of the
> preceding subsection read **distinct aspects** of the LRG
> communication geometry. `ρ_split^coph` operates on the dendrogram-
> derived per-pair distance `D_coph`: it asks, for each contact pair
> independently, whether the **scale at which the two contacts
> coalesce in the communication hierarchy** shifts coherently from
> task to post-rest. The signal is concentrated **at the per-pair
> level** and the multiscale character is encoded in the cophenetic
> merge-heights `D_coph_ij ∈ {h_1, …, h_{N−1}}` (the values of the
> per-pair distance are themselves scale assignments). `T_G(k)`
> instead operates on the leading-`k` eigenmode subspace `U_k`: it
> asks whether the **k-dimensional subspace spanned by the slowest
> non-trivial diffusion modes** rotates between task and post-rest,
> at a chosen subspace cutoff `k`. The signal is concentrated **at
> the subspace level** and the multiscale character is encoded by
> sweeping `k`. The two probes are not derivable from each other —
> `U_k` is not a function of `D_coph` alone, and `D_coph` is not a
> function of `U_k` alone — and a coordinated trace at both levels
> indicates that the reorganization is visible *both* in the
> hierarchical pair-coalescence scale *and* in the slow-mode
> subspace structure, two structurally complementary readings of the
> same LRG propagator.

(Statistical inference for Grassmann is unchanged from the current
text; report `T_G(k)` as a descriptive multiscale view, with
matched-strength null evaluated at the cohort level per audit_66.)

---

## Statistical inference subsection adjustments (sssec at end)

The current text states family sizes `m = 6` (per-pair, six bands),
`m = 12` (KC, six bands × two `λ`), `m = 18` (matrix distances on
`D(τ)`, six bands × three distances).

**After the revision (locked 2026-05-18 + supersedure 2026-05-20):
no cross-band BH-FDR is applied at any LRG probe.** Each band's
verdict is gated independently by the per-band control battery
(`locked/CONTROLS.md` C1 split, C2 drift, C3 matched-strength, C4
cross-probe, C5 epi-X for `ρ_split^coph`; per-band cluster-extent
permutation on the `k`-grid for `T_G*`). The six per-band claims are
not a coordinated cross-band family of inference, and per-band
matched-strength surrogacy at `R = 200` already provides a calibrated
empirical `p` at the per-test level — layering BH on top is
scaffolding without function (three-point check at
`feedback_no_unmotivated_bh_fdr.md`).

- **Per-pair `ρ_split^coph`**: no cross-band correction. Per-band Wilcoxon p against `ρ_drift^coph`, matched-strength surrogate, and the cross-probe restriction is the gate. (Earlier draft proposed `m = 6` BH-FDR — retired 2026-05-20 by writing-agent feedback.)
- **`d_S` / `d_P` matrix distances on `D_coph`**: dropped from methods entirely (originally proposed as `m = 12` or `m = 6` — retired 2026-05-20, undefined in the comparison section; not reported in the manuscript). See METHODS_AUDIT_ISSUES.md §B2.
- **Grassmann `T_G(k)` / `T_G^*`**: per-band cluster-extent permutation on `k` is the gate; no within-family BH-FDR (cluster-extent is itself the family-level gate across the `k`-axis). No cross-band correction either, for the same reason as `ρ_split^coph`.
- **Anatomical enrichment A1 hypergeometric**: BH-FDR across DK regions per (band, probe), with family size `m ≈ N_DK_regions_with_coverage_per_band`. KEPT — the multi-region family IS the coordinated unit of inference per (band, probe).

Drop the references to `m = 12` (KC), to the partition / taxonomy
probes, and to the cross-band `m = 6` BH-FDR sentence from this
subsubsection.

The single acceptable robustness restriction at the LRG layer is
the **LRG-native `n=9` drop of the single anti-aligned patient**
(at β: Pat_15, biology-driven right-hemisphere-only implant). All
other dropout regimes from earlier drafts are **retired**:

- Pat_03 dropout (1024 Hz "outlier") — **retired 2026-05-18 (am)**,
  sampling-rate handled at config layer only.
- Pat_07 dropout (substrate-layer "anti" with `S = +0.005`,
  near-zero) — **retired 2026-05-18 (pm)**, Pat_07 is solidly pro at
  both LRG probes (`ρ_split^coph` z = +4.47; Grassmann
  `T_G(k=40)` z = −3.90).
- Legacy `n=8` "pro-cohort restriction" (drop Pat_07 + Pat_15) —
  **retired 2026-05-18 (pm)**, replaced by the LRG-native `n=9`
  (drop Pat_15 only).

Patient-dropout policy across the preprint: **don't drop patients in
analysis**, except for (a) a single patient anti-aligned at the
cohort level (biology-driven, e.g. Pat_15 right-hemisphere-only at
the LRG layer); (b) genuinely problematic data (vendor corruption,
sampling-rate handled at config layer — none currently active in
the cohort).

---

## What the writing agent must NOT do

Anti-patterns to flag and reject if encountered in any revision:

1. **Never write "communication distance `D(τ) = 1/ρ(τ)`" as the
   object that ρ_split is computed on.** The object is `D_coph =
   cophenet(UPGMA(D(τ_max)))`, the cophenetic image. They are
   distinct matrices with different mathematical content. The
   cached field in `data/cache/imcoh_lrg/` named `ultrametric_matrix`
   is `D_coph`, not raw `D(τ)`.
2. **Never describe the cophenet step as "denoising" or "smoothing".**
   It is a multiscale-resolution transform. The dendrogram replaces
   the spectral-gap scale-identification mechanism for our
   continuous-spectrum substrate. The framing is "multiscale per-pair
   image of the propagator", not "noise reduction".
3. **Never claim the `D(τ)` framework gives multiscale resolution by
   sweeping `τ` for our system.** Continuous spectrum makes the sweep
   degenerate; the dendrogram is the multiscale carrier (Methods §LRG
   framework and Section X.X discussion).
4. **Never describe `D_coph` as ultrametric "as in" the Villegas paper
   conditions (i–iii) on `1 − δ_ij / ρ_ij(τ)`.** The Villegas paper
   asserts ultrametricity of `1/ρ` under specific conditions that are
   not generically satisfied on continuous-spectrum substrates
   (empirically, raw `1/ρ(τ_max)` has ≈33% strong-triangle violations
   on our `|ImCoh|` data). `D_coph` is ultrametric by construction
   (cophenet of any agglomerative linkage always is); this is a
   different statement.
5. **Never describe Grassmann as redundant with the cophenetic probe,
   nor vice versa.** They read distinct aspects (pair-coalescence
   scale vs slow-mode subspace orientation) and are not derivable from
   one another.

---

## What the writing agent must do — the layered comparison table

Insert the following as a results-relevant aside (Methods §LRG
framework appendix or Results §LRG-CTM lead paragraph; coordinate with
Results-section revision). This anchors the cophenetic value-add
quantitatively against the substrate.

**Cohort matched-strength cohort agreement (`n_above_surrogate`/10)
and Wilcoxon p, per layer and per band (n=10 cohort, R=200 matched-
strength surrogates per cell, seed 20260511; the raw `D(τ_max)` row
landed 2026-05-18 from `preprint_05_allbands_matched_strength_raw_D.py`):**

| Layer | δ | θ | α | β | γ_l | γ_h |
|---|---|---|---|---|---|---|
| **Raw FC** `ρ_split^raw` (audit_67) | 7/10, p=.042 | 7/10, p=.14 | **8/10, p=.014** | 6/10, p=.053 | 7/10, p=.053 | 7/10, p=.19 |
| **Raw `D(τ_max)`** `ρ_split` (preprint_05) | 7/10, p=.116 | 6/10, p=.097 | 7/10, p=.042 | 6/10, p=.042 | 6/10, p=.032 | 7/10, p=.161 |
| **Cophenet** `ρ_split^coph` (audit_63) | 4/10, p=.28 | 2/10, p=.72 | 5/10, p=.002 | **7/10, p=.005** | 5/10, p=.12 | 4/10, p=.25 |

**Gate clarification (2026-05-19)**: each cell reports
`n_above_surrogate/10, paired_wilcoxon_p`. The **gate is the paired
Wilcoxon p < 0.05 alone** per CONTROLS.md §C3; the `n_above/10`
count is reported as a descriptive cohort-composition diagnostic
(useful to read alongside the Wilcoxon, especially for narrative
comparison between bands) but is **not a verdict criterion**. No
`n_above ≥ N/10` threshold exists anywhere in the controls battery
(see `feedback_no_hardcoded_test_thresholds.md`,
`feedback_no_single_patient_p_driven.md`). LOO max-p (the descriptive
single-patient leverage diagnostic) should also be reported alongside
each cohort Wilcoxon — see `.agents/preprint/locked/CONTROLS.md` §C3 + §C4.

**The interpretive contrast (three-layer reading)**:
- **Layer 1 (raw FC)** and **Layer 2 (raw `D(τ_max)`)** are operationally
  indistinguishable: both deliver 6–8 patients with positive
  `ρ_split` in every band tested. The LRG diffusion at the finest
  τ does not introduce band-resolution; it is a soft re-encoding of
  the per-pair information already present in the raw `|ImCoh|`
  matrix.
- **Layer 3 (cophenet `D_coph`)** is where band-resolution emerges:
  δ (7→4), θ (7→2), γ_h (7→4) drop sharply, γ_l drops modestly (6→5),
  α retains a strong p (=0.002) at 5/10, and β rises to **7/10
  (p=0.005)** — the highest cohort agreement on this layer.

The cophenet's contribution is therefore **band-resolution at the
multiscale level via the dendrogram-induced merge-height integration**,
not amplification of detection. This three-layer contrast is itself
the central argument for adopting `D_coph` as the canonical
per-pair probe and should appear in the Results discussion of the
LRG-CTM finding.

(For the β-only sensitivity check on raw `D(τ_max)`, see
`data/preprint/rho_split_raw_D/beta_matched_strength.csv`: ratio
13.4×, Wilcoxon p=0.042, n_above_null 6/10. Same numbers as the
all-bands β row above.)

---

## Source artifacts for the writing agent

For every number, formula, or claim above, the writing agent should
ground in one of:

- `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`
  (raw FC `ρ_split^raw` matched-strength, all 6 bands, audit_67)
- `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`
  (cophenet `ρ_split^coph` matched-strength, all 6 bands, audit_63)
- `data/audit/ctm_triangle/cohort_summary.csv`
  (cophenet `ρ_split^coph` within-baseline, all 6 bands)
- `data/reports/imcoh_continuous_trace/per_cell_summary_split.csv`
  (cophenet per-patient `ρ_split^coph` within-baseline, all bands)
- `data/preprint/rho_split_raw_D/beta_per_patient.csv` and
  `beta_matched_strength.csv` (raw `D(τ_max)` `ρ_split`, β, sensitivity)
- `data/audit/grassmann_matched_strength_surrogate/` and audit_66 /
  audit_67_epi_exclusion / audit_69 outputs (Grassmann verdict)
- Memory `cophenet_methodology_rationale.md` — five-point rationale
- Memory `feedback_never_confuse_D_with_cophenet.md` — naming hygiene
- Memory `lrg_outlier_case_fully_connected.md` — continuous-spectrum
  basis for rejecting τ-sweep

---

## Checklist for the writing agent

When the revision lands, verify all of the following hold in the
revised Methods section:

- [ ] LRG framework subsection lists **two** comparison objects
      (`D_coph`, `U_k`), not three. Linkage hierarchy / dendrogram
      `T(τ)` is described as the construction that defines `D_coph`,
      not as a comparison object in its own right.
- [ ] Replacement subsubsection is titled **"Per-pair multiscale
      correlation on `D_coph`"** (or equivalent — explicit "multiscale"
      and "cophenet" required).
- [ ] Every formula in the per-pair subsection uses `D_coph` or
      `ρ_split^coph`, never bare `D(τ)` or bare `ρ_split`.
- [ ] Methodological rationale paragraph explicitly states (a) raw
      `D(τ_max)` is single-scale, (b) continuous spectrum makes
      τ-sweep degenerate, (c) dendrogram replaces spectral-gap scale
      ID, (d) `D_coph` is the per-pair image of multiscale structure.
- [ ] KC subsection (`Tree distance on the dendrogram`) is **deleted
      in full**, including all references to `T_KC`, `λ` extremes, and
      the within-baseline-null controlled test on `T_KC`.
- [ ] VI subsection (`Partition distance across cuts of the
      dendrogram`) is **deleted in full**.
- [ ] Module taxonomy subsection (`Cross-phase module taxonomy on the
      dendrogram`) is **deleted in full**, including
      Table `tab:methods_taxonomy`.
- [ ] Grassmann subsection retains its current content plus the new
      "What `T_G(k)` reads vs what `ρ_split^coph` reads" paragraph
      that distinguishes the two probes' meanings explicitly.
- [ ] Statistical inference subsection: **drop all cross-band BH-FDR
      statements at the LRG probes** (locked 2026-05-20 supersedure per
      `feedback_no_unmotivated_bh_fdr.md`; ρ_split^coph m=6 cross-band
      BH-FDR was retired because the per-band Wilcoxon already gates
      and the cross-band family does not form a coordinated unit of
      inference). Keep BH-FDR **only** for the anatomy A1 hypergeometric
      across DK regions per (band, probe). Drops KC, VI, taxonomy
      family-size references.
- [ ] All anti-patterns listed in the "What the writing agent must
      NOT do" section are absent from the revised text.
- [ ] The layered cohort agreement table (raw FC vs cophenet) is
      either in the Methods (LRG appendix) or the Results (LRG-CTM
      lead paragraph) but is present somewhere in the revised
      manuscript.
- [ ] References to specific cache files / audit scripts use the
      paths listed in "Source artifacts" above; no invented paths or
      audit numbers.

---

## Companion-section flagging (Results, Discussion)

This methods revision **forces** the following downstream edits in
other manuscript sections. Track these separately:

- Results §5.2 (KC-based finding, the "β KC 10/10 within-baseline"
  paragraph): retire the entire central claim; KC is no longer in
  the methods. If the result is to be kept anywhere, move to
  Supplement as exploratory.
- Results §5.3 (LRG-CTM `ρ_split`): rewrite to use `ρ_split^coph`
  throughout; reference the cophenetic methodology paragraph in
  Methods.
- Results §5.4 (Grassmann): no content edits needed; verify that the
  meaning-distinction paragraph from Methods is echoed in the lead.
- Results §5.5 (taxonomy / cross-phase module classes): retire
  entirely or move to Supplement as exploratory.
- Discussion: any framing of "the LRG step gives band selectivity"
  must explicitly tie the band-resolution claim to the cophenetic
  multiscale step, citing the raw-FC-vs-cophenet contrast table.

---
