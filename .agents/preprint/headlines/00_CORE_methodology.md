---
name: headline-core-methodology
era: IMCOH_ABS_COHORT_N10
status: current
kind: headline
scope: CORE — the multiscale Laplacian framework from which N1–N3 emerge. Two probes (ρ^coph per-pair multiscale + Grassmann global-mode), the C1–C5 null battery, the multiscale-beats-global-spectral argument, band-resolution-not-amplification. This is the methods spine, not a neuro headline.
owner_agent: preprint-general-questioning (null model)
updated: 2026-06-22
---

# CORE — A multiscale Laplacian diffusion read-out of sEEG connectivity

> Methods spine. The three neurophysiological headlines (N1 trace, N2 content,
> N3 epi) all hang off this. Numbers live in cached CSVs (§F) — this file argues
> the method, it does not tabulate results.

## §A — In plain language

Brain functional connectivity is a weighted network, and the usual ways to
summarize it throw away most of its structure: pick one threshold, pick one
scale, or compress everything into a handful of "global" components (PCA,
spectral clustering, the leading eigenmodes). Those summaries are **scale-blind**
and **hard to interpret** — they tell you *that* something moved, not *which
relationships at which scale*.

We instead let the network **diffuse**. The Laplacian density operator
ρ̂(τ) = e^{−τL̂}/Z (the Laplacian Renormalization Group construction) describes
how signal spreads through the graph over a continuum of scales τ, and from it we
build the full **hierarchy** of how tightly every pair of contacts is bound
together across *all* scales at once. We then read this hierarchy two ways:

- a **per-pair, multiscale** distance (`ρ^coph`) — interpretable and localizable:
  it says, for each pair of brain sites, how their multiscale "communication
  distance" changed; and
- a **global-mode** distance (`d_G`, Grassmann) — the rotation of the leading
  diffusion eigenmodes, i.e. essentially the information a spectral-clustering /
  PCA pipeline would keep.

The punchline of the method is that **these two read-outs are not redundant** —
each band shows up on one, the other, or both — which proves that a global /
spectral summary alone *misses* real, multiscale, per-pair reorganization. And
every claim is judged against a **strength-matched surrogate** (a network with
the same per-node connection strengths but scrambled structure): the referee
that has killed several earlier headlines.

## §B — Technical statement

### Substrate (one choice, defended)
FC = band-averaged magnitude of imaginary coherence `<|ImCoh|>_f` (`imcoh_abs`;
Nolte 2004; Ewald 2012; Bastos & Schoffelen 2016). Immune by construction to
zero-phase-lag artifacts — same-shaft volume conduction and common-reference
inflation — so spatial proximity cannot leak into any result (no spatial null
needed). Raw |ImCoh| is the **comparison baseline only, never a result**
(`feedback_results_only_in_laplacian_framework`).

### The operator and the two probes
- ρ̂(τ) = e^{−τL̂}/Z, with the canonical scale τ = 1/λ_max. **τ-robustness:** a
  τ-sweep (`audit_121`) vindicates this choice — the genuine trace is fine-scale
  and τ-stable; the apparent "gains" at coarse τ are a collapse artifact (a
  placebo rises to match past the Fiedler scale). So the result is not a
  τ-cherry-pick.
- **`ρ^coph`** — cophenetic distance of UPGMA over `D(τ_max)`: a per-pair number
  integrating all N−1 dendrogram merge-height scales. Per-pair → **localizable**
  and **scale-resolved**.
- **Trace estimator: `ρ_sym`** (of record since 2026-07-06; supersedes bare
  `ρ_split`). The cross-phase trace on `ρ^coph` uses a **split** rest_pre baseline —
  `Spearman(D_task − D_preA, D_post − D_preB)` with independent halves A, B — so
  shared baseline noise cannot inflate it. Bare **`ρ_split`** hard-coded *which* half
  fed *which* arm, an arbitrary choice that flipped the sign of near-zero patients
  (20/60 patient×band cells under A↔B swap; `audit_149`). **`ρ_sym = ½[ρ(A→task,B→rest)
  + ρ(B→task,A→rest)]`** averages the two equally-valid assignments → invariant to the
  half-labelling, same data, same matched-strength surrogate. Every verdict is
  **estimator-invariant (0/6 bands flip)**; ρ_sym is deliberately *more conservative*
  (β gate p=0.032, α 0.024 at R=200) because it stops a lucky half pushing ill-conditioned
  near-zero patients above their surrogate. Per-patient reporting: `ρ_sym ± ½|ρ_AB−ρ_BA|`,
  with `|ρ_sym| < 1 SE` labelled **undetermined**. The full (shared) baseline was tested
  and **rejected** — it reintroduces shared-error inflation (Pat_02 obs 0.775 vs surrogate
  median 0.767). `ρ_split` retained as a supplement column only. (`audit_150`;
  `.agents/reports/2026-07-06_rho-sym-panoramic-and-methodology.md`.)
- **`d_G(k)`** — chordal (Grassmann) distance between the leading-k Laplacian
  eigenmode subspaces, swept over k. Whole-network, leading-mode → this is the
  **global/spectral** view. (Separate measure; unaffected by the ρ_split→ρ_sym change.)

### The C1–C5 null battery (frame completely; reference, don't re-derive)
- **C1** within-baseline split-half null (ρ^coph): trace > within-`rest_pre`
  noise?
- **C2** drift-floor null (ρ^coph): trace > slow within-session drift?
- **C3** **matched-strength surrogate — MANDATORY, both probes**: trace beyond
  degree-preservation? *The referee.* (ρ^coph: split-baseline surrogate;
  Grassmann: cluster-extent permutation on matched-strength surrogates,
  Decision-8 mass gate + Decision-12 LOO precondition.)
- **C4** cross-probe (cross-shaft) restriction (ρ^coph): paired Wilcoxon, does
  restricting to cross-shaft pairs degrade the trace? *(Control only — do not
  build figures around it; see README §5.)*
- **C5** epi-zone exclusion (both probes, **secondary/mechanistic**): does the
  trace survive/strengthen with SOZ contacts dropped? Always paired with the
  **node-count decimation control** (`audit_85`) — rebuilding LRG on any K-node
  submatrix inflates ρ^coph, so "exclude X → strengthens" is a node-count
  confound until decimation-controlled (`feedback_decimation_control_for_subset_exclusion`).

### Central methodological claims (these belong in the Methods/Results spine)

**Claim 1 — Multiscale per-pair geometry beats global-spectral summaries
(non-redundancy).** The locked band × probe verdicts dissociate:
- **β** carries a trace on **both** probes;
- **α** carries a trace on **`ρ^coph` only** — a per-pair multiscale persistence
  with *no* detectable leading-mode rotation. **A global/spectral analysis would
  miss α entirely.**
- **γ_low and δ** carry a trace on **Grassmann only** — the global subspace
  persists while the per-pair multiscale geometry is reorganized (no persistent
  cophenet trace). A global summary reads "persistence" and **misses that the
  multiscale per-pair structure has been reshuffled**.
- **β** = both; **θ, γ_high** = neither.

Because neither probe subsumes the other, the leading-eigenmode / PCA /
spectral-clustering view is provably **insufficient** for this question:
important multiscale information is only disentangled by tracing **per-pair paths
through the propagator**. This is the methodological reason the paper exists, and
it is the cited justification inside N1 (β both probes) and the framing of N2/N3.

**Head-to-head confirmation — the external baseline (`audit_143`, 2026-06-22).**
To meet the obvious referee demand (*"run an actual spectral-clustering / PCA
baseline"*), we put a literal **leading-k spectral embedding** (the NJW-normalized
representation k-means clusters) through the **same `ρ_split` statistic** as the
cophenetic distance — swapping *only* the per-pair distance (cophenetic = all N−1
scales; embedding = leading-k modes; with raw-FC and single-τ diffusion arms).
Cophenetic reproduced its locked α/β `ρ_split` per-patient **bit-exact** (the
correctness gate). Findings, honest limit first:
- **β is a tie** — the leading-k embedding recovers β as well as cophenetic; the
  method-superiority claim is **never** made on β magnitude.
- **α is cophenetic-only across *both* spectral read-outs** — *neither* the per-pair
  embedding *nor* the subspace Grassmann recovers α, while the full multiscale
  cophenetic does (BH-clears within band). A standard spectral-clustering / PCA
  analysis would **not see the α trace at all** — the sharpest single proof the
  multiscale step is necessary.
- **Selectivity is the cleaner edge** — cophenetic fires on *exactly* the two signal
  bands (α, β); the spectral read-outs either miss α (Grassmann) or **over-detect
  into the null band** (the embedding fires on θ, where cophenetic and Grassmann are
  both null) — so the team's Grassmann is itself a *better-behaved* spectral
  baseline than the textbook embedding. A diffusion-τ arm localizes the advantage in
  the **multiscale UPGMA aggregation**, not the heat kernel: single-τ diffusion gets
  α only weakly and leaks into γ_low; UPGMA sharpens α *and* adds the selectivity —
  so it is *multiscale*, not *diffusion*, doing the work.
  Numbers: `data/audit/spectral_distance_swap/cohort_summary.csv`.

**Claim 2 — Band resolution, not amplification.** A three-layer contrast (raw FC
→ raw `D(τ_max)` → cophenet `ρ^coph`) shows raw FC already "detects" change in
essentially every band at 6–8/10 cohort agreement; the cophenet/multiscale wrap
**selectively demotes** δ/θ/γ_high to obs ≈ surrogate while preserving β (and α).
The multiscale step contributes **which bands carry structured, hierarchical
persistence**, not raw detection power.

## §C — Critical issues & powerful strengths

**Strengths.** Interpretable + multiscale + localizable (per-pair); a single
mandatory referee (matched-strength) applied uniformly; volume-conduction-immune
substrate; τ-robust; the probe dissociation is itself a positive result for the
method, not just a caveat.

**Critical issues (lead with these).**
- The **matched-strength surrogate** is the entire credibility of the programme.
  It must be defended as adequate (on-manifold CRC cross-check reproduces it;
  `crc_onmanifold_surrogate_2026_06_10`). The **null-model design is the live
  topic** of the preprint-general-questioning chat — keep it there.
- **Whole-brain concordance is corruptible; the demeaned per-system read-out is
  not** (lesson from N2.5, 2026-06-22). A contiguous-window *truncation* null
  injected a common-mode structural shift that an **un-demeaned** whole-brain
  concordance read as fake signal (it false-positived a dead-null band, δ); the
  **per-patient demeaned** read-out was immune — which is why the demeaned
  localization null stayed clean while the truncation arc-null did not. Reusable
  rule: prefer demeaned per-unit read-outs + artifact-free controls (e.g. a
  length-ratio regression) over un-demeaned global concordance whenever a null
  perturbs global structure. Worth one Methods sentence.
- **Grassmann is an anchor-flavored, phase-quantity in some uses** — be careful
  the global probe is a genuine *cross-phase* trace, not a phase-averaged anchor
  (this conflation caused an anatomy retraction). State exactly which Grassmann
  quantity each claim uses.
- **θ "anti-trace" is more complicated than a sign** — see README §5; do not
  assert it; characterize layer-dependence + outlier leverage.
- **Per-patient leverage / outliers** — under ρ_sym the β spread is cleanly tiered
  (`audit_154`): 6 stable tracers, 2 **undetermined** (Pat_13, Pat_15 — `|ρ_sym|<1 SE`),
  2 stable resets (Pat_10, Pat_14). Pat_15 is now undetermined, **not** β-anti — revisit
  its sanctioned-dropout framing. Report LOO-max; no single-patient "strong" tags.

## §D — To-dos & verifiables (owner: preprint-general-questioning / null model)

- [x] **(DONE 2026-06-22, `audit_143`)** Non-redundancy argument written into
  Claim 1: the internal band × probe dissociation **plus** the external leading-k
  spectral-embedding head-to-head (α cophenetic-only across both spectral read-outs;
  selectivity; β tie stated honestly). Drive any panel from
  `data/audit/spectral_distance_swap/cohort_summary.csv`, not hardcoded.
- [ ] **Defend the matched-strength null** against the obvious referee attacks
  (degree-preservation sufficiency; independence per phase vs coordinated
  cross-phase — the audit_63 coordinated-null gap is open).
- [ ] **θ characterization** — is the cophenet anti-direction robust or
  outlier/layer-driven? Resolve before any θ statement.
- [ ] **τ-robustness** writeup from `audit_121` (fine-scale, collapse-artifact at
  coarse τ); decide if it's a Methods figure or supplement.
- [~] Explicit **interpretability contrast** vs PCA/spectral-clustering — the
  *significance* half is **DONE** (`audit_143` α head-to-head); the *localization*
  vignette still open: a worked example where the global modes look unchanged but
  ρ^coph localizes the α reorganization to a region (α is the natural candidate).

## §E — Figure / representation ideas

- **Schematic**: FC matrix → ρ̂(τ) diffusion → dendrogram (per-pair `ρ^coph`) and
  leading eigenmode subspace (`d_G`) side by side — "two read-outs of one
  operator."
- **Band × probe dissociation matrix** (β both / α coph-only / γ_low,δ
  Grassmann-only / θ,γ_high none) — the single most sellable methods panel; it
  *is* Claim 1. (Drive from CSV, not hardcoded.)
- **Three-layer band-resolution** strip (raw FC → raw D(τ) → cophenet) showing
  the selective demotion.
- A **global-modes-look-the-same-but-per-pair-moved** α vignette.
- **Three-way head-to-head** (`audit_143`): cophenetic vs leading-k spectral
  embedding vs Grassmann, per band — α caught by cophenetic *only*; the selectivity
  contrast (cophenetic = α,β only; embedding leaks into θ). Drive from
  `spectral_distance_swap/cohort_summary.csv`. A strong companion to the dissociation
  matrix — it answers "did you compare to actual spectral clustering?" in one panel.
- Do **not** make a C4 cross-probe figure; **do** keep C5 X-epi panels (README §5).

## §F — Provenance (CSV · script · timestamp)

- Probe dissociation / per-band verdicts → `../locked/VERDICT_LEDGER.md` (locked
  2026-05-18, rev. through 2026-05-28); Grassmann gate
  `data/audit/grassmann_cluster_extent/cohort_summary.csv` ·
  `audit_70_grassmann_cluster_extent.py` · 2026-05-26.
- C3 ρ^coph matched-strength (ρ_split ref, supplement) →
  `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` ·
  `audit_63_split_baseline_surrogate.py` · 2026-05-15.
- **ρ_sym estimator (of record)** — gate `data/audit/rho_sym_gate/`
  (`audit_150`, β p=0.032/α p=0.024); estimator-robustness `data/audit/estimator_regate/`
  (`audit_149`, 0/6 flips, full-baseline rejected); downstream migrations
  `data/audit/{localization_atlas_rhosym,consolidation_arc_rhosym,cross_phase_taxonomy_rhosym,per_node_trace_decomposition_rhosym}/`
  (`audit_151`–`155`) · 2026-07-06. Surrogate shuffle now numba-JIT (98× faster,
  bit-identical). Panoramic + supplement:
  `.agents/reports/2026-07-06_rho-sym-{panoramic-and-methodology,pipeline-migration,estimator-robustness-supplement}.md`.
- C1/C2/C4 →
  `data/audit/ctm_triangle/cohort_summary.csv` · `audit_33_ctm_triangle.py` ·
  2026-05-26 (C4: `c4_wilcoxon_cohort.csv` · `audit_71_c4_wilcoxon_cohort.py` ·
  2026-05-19).
- Three-layer band resolution →
  `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv` (2026-05-11),
  `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`,
  + the C3 CSV above; synthesized in `../bands/00_cohort.md` §2.
- τ-robustness → `audit_121_tau_sweep_cophenetic_trace.py` · 2026-06-22 (+ figure
  `audit_121b_tau_sweep_figure.py`); memory `tau_sensitivity_trace_2026_06_22`.
- Null adequacy cross-check → `crc_onmanifold_surrogate_2026_06_10` (audit_96/97).
- **Spectral head-to-head (external baseline, non-redundancy)** →
  `data/audit/spectral_distance_swap/cohort_summary.csv` ·
  `audit_143_spectral_distance_swap_headtohead.py` · 2026-06-22 (cophenetic vs
  leading-k embedding vs raw-FC vs single-τ diffusion, all on the same `ρ_split`;
  cophenetic arm reproduces audit_63 α/β per-patient bit-exact). Precursor
  `audit_141_spectral_clustering_headtohead.py` (v1, triangle+ARI/NMI — superseded:
  wrong statistic).

## §G — Missing parts / open

- The **coordinated cross-phase null** (audit_63 gap) — matched-strength is
  per-phase independent; a referee may ask for a jointly-coordinated null.
- ~~A crisp, citable PCA/spectral-clustering baseline run on the same data~~
  **DONE 2026-06-22 (`audit_143`)** — the leading-k spectral-embedding head-to-head
  closes this: α is cophenetic-only across *both* spectral read-outs (embedding +
  Grassmann), β is a tie, and selectivity is the edge. What remains open is only the
  **interpretability/localization vignette** (§D) — showing the global modes look
  unchanged while ρ^coph localizes the α reorganization to a region.
