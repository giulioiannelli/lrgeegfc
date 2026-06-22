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
- **`d_G(k)`** — chordal (Grassmann) distance between the leading-k Laplacian
  eigenmode subspaces, swept over k. Whole-network, leading-mode → this is the
  **global/spectral** view.

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
- **Grassmann is an anchor-flavored, phase-quantity in some uses** — be careful
  the global probe is a genuine *cross-phase* trace, not a phase-averaged anchor
  (this conflation caused an anatomy retraction). State exactly which Grassmann
  quantity each claim uses.
- **θ "anti-trace" is more complicated than a sign** — see README §5; do not
  assert it; characterize layer-dependence + outlier leverage.
- **Per-patient leverage / outliers** (Pat_15/10/08/02) — report LOO-max; no
  single-patient "strong" tags.

## §D — To-dos & verifiables (owner: preprint-general-questioning / null model)

- [ ] Write the **non-redundancy argument** as a standalone Methods-Results
  paragraph: tabulate (in the CSV, cite here) the band × probe dissociation and
  state the "global-spectral misses α / misreads γ_low,δ" conclusion explicitly.
- [ ] **Defend the matched-strength null** against the obvious referee attacks
  (degree-preservation sufficiency; independence per phase vs coordinated
  cross-phase — the audit_63 coordinated-null gap is open).
- [ ] **θ characterization** — is the cophenet anti-direction robust or
  outlier/layer-driven? Resolve before any θ statement.
- [ ] **τ-robustness** writeup from `audit_121` (fine-scale, collapse-artifact at
  coarse τ); decide if it's a Methods figure or supplement.
- [ ] Explicit **interpretability contrast** vs PCA/spectral-clustering: a worked
  example where the global modes look unchanged but ρ^coph localizes a real
  per-pair reorganization (α is the natural candidate).

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
- Do **not** make a C4 cross-probe figure; **do** keep C5 X-epi panels (README §5).

## §F — Provenance (CSV · script · timestamp)

- Probe dissociation / per-band verdicts → `../locked/VERDICT_LEDGER.md` (locked
  2026-05-18, rev. through 2026-05-28); Grassmann gate
  `data/audit/grassmann_cluster_extent/cohort_summary.csv` ·
  `audit_70_grassmann_cluster_extent.py` · 2026-05-26.
- C3 ρ^coph matched-strength →
  `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` ·
  `audit_63_split_baseline_surrogate.py` · 2026-05-15.
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

## §G — Missing parts / open

- The **coordinated cross-phase null** (audit_63 gap) — matched-strength is
  per-phase independent; a referee may ask for a jointly-coordinated null.
- A crisp, citable **PCA/spectral-clustering baseline** run on the same data, to
  make Claim 1 a head-to-head rather than an argument from the probe dissociation
  alone.
