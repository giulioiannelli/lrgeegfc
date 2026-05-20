---
name: 2026-05-08_lrg-epilepsy-research-directions
type: plan
era: IMCOH_ABS × COHORT_N10
status: active
created: 2026-05-08
updated: 2026-05-08
pointers:
  - .agents/reports/2026-05-07_epileptic-n10-revisit.md
  - .agents/reports/2026-05-08_direction-a-induced-subtree.md
  - .agents/reports/2026-05-08_trace-minus-epi-resection.md
  - .agents/reports/2026-05-05_result-2-lrg-beta-trace.md
  - .agents/guides/01_project/terminology.md
  - .agents/guides/02_methods/lrg-framework-guide.md
  - .agents/plans/active/2026-04-29_eigenvector-direct-pivot-plan.md
---

# LRG-epilepsy research directions — connectivity-pattern map

## Renormalization head

**The trace / anchor / reset / rearrange (TARR) taxonomy was the
natural starting point given Section 5's loaded β-trace finding.
For epilepsy specifically, that was the *wrong* default: epileptic
networks have many connectivity-pattern signatures and only a few of
them are about cross-phase reorganisation. The plan now indexes
directions by the graph-theoretic primitive they probe — spectral
(eigenvectors and eigenvalues of the LRG Laplacian L̂), path-integral
(heat kernel ρ̂(τ) and its derivatives), or distance-distribution
(cophenetic / resistance / commute time on the FC graph). The
2026-05-08 round closed Direction A (trace on epi-induced subtree:
band-resolved phenomenology with β-dissociation and δ-heights RESET)
and Direction A' / audit_56 (β trace is non-epi-carried, robust to
resection; α/λ=0 sharpens). The user's explicit follow-ups —
Grassmann embedding of epi nodes and Laplacian-driven information
flow — now sit at the top of the queue alongside Direction C
(eigenmode localisation). Trace is one connectivity pattern among
many, not the umbrella; the plan reflects that.**

## TL;DR

1. **Trace-frame closed for epilepsy.** A + audit_56 done; reported
   as band-resolved phenomenology, not a "discovery" but a clean
   disambiguation that includes β-dissociation, δ-heights RESET,
   α-topology TRACE.
2. **Connectivity-pattern directions go beyond trace.** Three
   families: spectral (E, C), path-integral (F, G, B),
   distance-distribution (H, D).
3. **Top of queue (user-flagged):** Grassmann embedding (E),
   Laplacian information flow (G), eigenmode localisation (C).
4. **Each direction tests a *different question* about how the
   epi nodes sit in the network**, not a different operationalisation
   of the same question.
5. **Acceptance criterion stays the same**: cohort regularity passing
   the IQR > 20% × |median| floor (per
   `feedback_iqr_vs_cohort_slope.md` /
   `feedback_regularity_over_bh_null.md`); BH-FDR is a footnote, not
   a gate.
6. **Order is question-driven, not feasibility-driven.** Tackle the
   directions whose answers would actually change how we think about
   the epileptic network.

---

## What's done (2026-05-08)

### A — Cross-phase rigidity (TARR on epi-induced subtree)

`audit_54_epi_rigidity_compute.py` + `audit_54b_*_figures.py`.
Verdict at
`.agents/reports/2026-05-08_direction-a-induced-subtree.md`:

- **β-band dissociation** (RESET at the same band where Result 2
  shows the cohort TRACE).
- **δ-heights RESET** (75% IQR — largest single effect).
- **α-topology TRACE** (55% IQR).
- audit_48 fig_05 α-topology pair-mask result (8/9 p=0.006)
  attenuates to 5/9 under induced subtree but direction preserved.

### A' — Virtual resection on the LRG tree (audit_56)

`audit_56_trace_minus_epi.py`. Verdict at
`.agents/reports/2026-05-08_trace-minus-epi-resection.md`:

- **β trace is non-epi-carried** — robust to resection per-patient,
  cohort Δmed ≤ 0.2 at all three λ values.
- **α/λ=0 trace strengthens** with resection (-0.81 → -1.29; IQR
  fraction 19% → 40%).
- **δ/λ=0 direction flips** — global RESET driven by epi-rest
  interface; non-epi part flat.

These two together close the trace-frame for epilepsy at the
cohort level.

---

## Connectivity-pattern directions — three families

Each direction below has:
- a one-line *hypothesis* (what's the question?),
- a *primitive* (what graph-theoretic object is computed?),
- a *cohort signature* (what would a positive finding look like?),
- *connection to existing tools* (so we know the implementation is
  cheap),
- *feasibility* (cache reuse, library reuse, blocking dependencies).

When ready to implement, each direction gets its own rigorous scope
report under
`.agents/guides/task-persistence-investigation/` (despite the name
— that folder is the canonical home for any new descriptive
multiscale measure on this project, per the never-always rule).

### Family 1 — Spectral (eigenvectors / eigenvalues of L̂)

#### **E — Grassmann embedding of epi vs full eigenspace** *(user-flagged)*

- **Hypothesis.** The k-dimensional subspace spanned by the top-k
  L̂ eigenvectors restricted to the epi rows (a |E_p|×k matrix
  treated as a k-dim subspace embedded in |E_p|-dimensional space)
  is *not* a random projection of the full top-k eigenspace. The
  principal angles between V_k[E_p, :] and V_k[full] (after
  appropriate normalisation) reveal whether the epi nodes "carry"
  the modes (small angles) or "interrupt" them (large angles).
- **Primitive.** Grassmann distance via principal angles. Two
  variants:
  - *Direct.* `d_chord(V_k[E_p, :], V_k[N_p, :]) = sqrt(k - Σ σᵢ²)` —
    Grassmann chordal distance between the epi-row block and the
    non-epi-row block. Tests whether epi and non-epi rows span
    similar k-dimensional subspaces (the same modes "live on"
    both regions) or different ones (the epi rows project onto
    a different sub-space than the non-epi rows).
  - *Restricted-to-full.* `d_chord(V_k_resect, V_k_full)` — same
    measure but between the eigendecomposition of the resected
    FC graph and the full FC graph. Tests how much the epi
    nodes alter the global eigenstructure.
- **Cohort signature.** Per (patient, band, k, phase): one
  Grassmann distance. Cohort regularity: cohort-median Grassmann
  distance vs random-leaf-subset null at matched size.
- **Connection.** E1 spectral subspace alignment scope at
  `.agents/guides/task-persistence-investigation/2026-04-29_e1-spectral-subspace-alignment.md`
  has the cross-phase version; this is the cross-region version.
  Same Grassmann primitive, different leaf-set partitioning.
- **Feasibility.** Eigenvectors already cached (audit_48 reads
  `r.eigenvalues` and `r.eigenvectors` successfully). Pure
  linear algebra; ~50 lines. **Quickest spectral direction.**

#### **C — Eigenmode localisation on E_p (m^E_k / IPR)** *(scope written)*

- **Hypothesis.** There exists a normalised-eigenvalue range
  `λ̂ = λ/λ_max ∈ [a, b]` where the epi-projection mass
  `m^E_k = Σ_{i ∈ E_p} v_k(i)²` is anomalously high vs a
  degree-corrected null — a "characteristic communication scale"
  at which epi modes don't propagate to the rest of the brain
  (Anderson-localisation analogy).
- **Primitive.** `m^E_k` per mode; binned on `λ̂` grid;
  strength-stratified null (Q=5 quintiles, R=100 draws).
- **Status.** Scope at
  `.agents/guides/task-persistence-investigation/2026-05-08_epi-eigenmode-localization.md`
  (audit_57 to be assigned; old slot 55 taken by
  `audit_55_ctm_per_pair_scatter.py`).
- **C1** (static) is independent; **C2** (cross-phase shift) is
  conditional on C1 cohort-positive.
- **Feasibility.** Cache already in use; ~3 days dev.

### Family 2 — Path-integral (heat kernel ρ̂(τ) and derivatives)

#### **G — Laplacian information flow / heat flux** *(user-flagged)*

- **Hypothesis.** The rate of communication between an epi node
  and the rest of the network is *qualitatively different* from
  what node degree alone predicts. Operationalised via the
  derivative of the heat kernel: `dρ̂_ij/dτ` at small τ measures
  the initial diffusion velocity from i to j. Per epi node, two
  fluxes to characterise:
  - *Internal flux* — `Σ_{j ∈ E_p \ {i}} dρ̂_ij/dτ` (diffusion
    velocity within the epi sub-network).
  - *External flux* — `Σ_{j ∈ N_p} dρ̂_ij/dτ` (diffusion velocity
    from epi i to non-epi).
  - *Balance ratio* — `External / (Internal + External)`. A balance
    < 0.5 means information stays inside the epi cluster; > 0.5
    means it leaks out.
- **Primitive.** `dρ̂/dτ = -L̂ K̂(τ)` evaluated at τ = τ_min,
  τ_min · 2, τ_min · 4 etc. This is well-defined and pure linear
  algebra.
- **Cohort signature.** Per (patient, band, τ): cohort-median
  Balance ratio for epi nodes vs non-epi nodes. If `Balance(epi) <
  Balance(non-epi)` consistently, epi nodes are "communication
  sinks" (information enters but doesn't leave).
- **Connection.** ρ̂ already cached. The τ-derivative is one
  matrix multiplication per τ. Maps directly to the Bartolomei
  PZ literature ("propagation zone" = where seizures pass through;
  the high-Balance-ratio non-epi nodes are PZ candidates).
- **Feasibility.** ~2 days dev. **Theoretically rich and clinically
  interpretable.**

#### **F — Multiscale communicability centrality**

- **Hypothesis.** The total communicability of a node,
  `C_i(τ) = Σ_j ρ̂_ij(τ)`, distinguishes hub nodes (high C) from
  peripheral nodes (low C). At each τ, the rank of an epi node
  in the cohort-wide centrality distribution is a "multiscale hub
  signature". Are epi nodes systematically more central at some
  τ regime?
- **Primitive.** `C_i(τ) = (Σ_j K̂_ij(τ)) / sqrt(K̂_ii K̂_jj)` —
  normalised. Per-patient, rank epi vs non-epi.
- **Cohort signature.** Per (patient, band, τ): mean centrality
  rank of E_p (relative to non-epi). Cohort-median rank vs
  random-leaf-subset null.
- **Connection.** Same eigenvector cache as C, E. One τ-sweep per
  patient.
- **Feasibility.** ~2 days dev. Smaller scope than G.

#### **B — True virtual resection on FC** *(carried over from earlier plan)*

- **Hypothesis.** The drop in non-epi-to-non-epi communicability
  Δ*T*(τ) when the epi rows/columns are deleted from `A` (then
  recomputing `L̂_resect` and `ρ̂_resect`) is anomalously large at
  some τ regime — epi as critical relay.
- **Primitive.** `T(G; τ) = Σ_{i,j ∈ N_p} ρ̂_ij(τ)`; resect:
  `T(G_resect; τ) - T(G; τ)`.
- **Cohort signature.** Cohort z-score of epi-resection ΔT vs
  size-matched random-leaf-subset resection. Strong z = epi is
  critical relay; near-zero z = epi is redundant.
- **Connection.** audit_56 did the *tree-level* resection (induced
  subtree on V \ E_p of the cached linkage). This direction does
  the *FC-level* resection (delete rows/cols of A, recompute LRG
  from scratch). Stronger version of the same question; tests
  what audit_56 can't (the eigenvectors are different on the
  resected graph).
- **Feasibility.** Recompute LRG for n=9 × 6 bands × 3 phases × R
  random-resection nulls. Heavier; ~4 days dev. Limited by
  missing surgical-outcome metadata for clinical interpretation.

### Family 3 — Distance-distribution (cophenetic / resistance / commute time)

#### **H — Effective resistance / commute time on epi pairs**

- **Hypothesis.** The distribution of effective resistances
  `R_ij = (e_i - e_j)^T L̂^+ (e_i - e_j)` for epi-epi pairs is
  characteristically different from non-non or random size-matched
  pairs. (Effective resistance ↔ commute time on the random walk
  defined by L̂.)
- **Primitive.** Pseudo-inverse of L̂; `R_ij` matrix is `r_ii + r_jj
  - 2 r_ij` where `r = L̂^+`. Cheap once the eigendecomposition
  is in hand.
- **Cohort signature.** Per (patient, band): KS or earth-mover
  distance between epi-epi `R` distribution and non-non
  distribution. Or simpler: cohort-median ratio.
- **Connection.** Resistance distance is the standard Markov-chain
  distance on a graph. Distinct from the LRG ρ̂(τ) at any single
  τ (resistance integrates over all τ). Different geometric
  primitive than what we've used so far.
- **Feasibility.** Pseudo-inverse via existing eigenvector cache;
  ~1 day dev. Smallest and quickest.

#### **D — ρ̂-leakage / propagation-zone (PZ) candidates**

- **Hypothesis.** For each non-epi node i, the fraction of its
  communicability flowing into the epi set
  `ℓ_i(τ) = Σ_{j ∈ E_p} ρ̂_ij(τ) / Σ_{j ≠ i} ρ̂_ij(τ)`
  identifies PZ candidates — non-epi nodes communicatively
  "close" to the epi set.
- **Primitive.** Per-node `ℓ_i`; ranked top-k.
- **Cohort signature.** Anatomical clustering of top-k PZ
  candidates with E_p neighbours (Desikan-Killany overlap)
  beyond chance.
- **Status.** Original Direction D in this plan; preserved here.
- **Feasibility.** ~2 days dev; most preliminary direction.

---

## Recommended order (question-driven)

1. **E — Grassmann embedding** (user-flagged, fastest spectral
   direction, ~50 lines). Tests *do epi nodes carry the global
   modes or not*. Either answer is theory-relevant.
2. **G — Laplacian heat flux** (user-flagged, theoretically
   richest, clinically interpretable via PZ link). Tests *does
   information leak out of the epi sub-network or stay trapped*.
3. **C — Eigenmode localisation** (scope written, ready to
   implement). Tests *do epi modes have a characteristic
   communication scale at which they don't propagate*.
4. **F — Multiscale centrality** (smaller scope; descriptive
   per-patient figure).
5. **H — Effective resistance** (smallest, complementary).
6. **B — True virtual resection on FC** (heaviest, defer until
   surgical metadata available).
7. **D — PZ leakage** (most preliminary; layered on top of G's
   findings if positive).

The order is question-driven, not feasibility-driven: tackle
directions whose answers would change how we think about the
epileptic network's connectivity. E + G + C give us a 3-pronged
spectral/path-integral attack within ~10 days dev total.

## Acceptance criterion (shared across all directions)

Per `feedback_regularity_over_bh_null.md`:

- **Lead with the cohort shape.** Per-patient distribution + median
  + IQR + per-band signature.
- **IQR floor:** report as a regularity if `|median| ≥ 0.20 × IQR`
  (Eric's rule). Cohort consistency ≥ 8/9 patients in the same
  direction is a separate gate.
- **BH-FDR / Wilcoxon as footnote.** Don't lead with "null" or "all
  cells significant".
- **Theoretical reading attached** to every regularity: what does
  the per-band / per-τ / per-k pattern *mean* in the context of
  epilepsy phenomenology?

## Slot allocation

Audit slots taken: 50, 51a/51b, 52, 53, 54 (Direction A),
54b (figures), 55 (CTM), 56 (Direction A' resection).

**Next free slots:** 57+ for Directions E, G, C, F, H, B, D as
they get implemented.

## What this plan deliberately does NOT do

- Re-litigate the trace frame. A + A' close it.
- Promise clinical biomarkers. Effect sizes are modest; n=9 is
  small.
- Promise surgical-outcome prediction. Metadata not available.
- Run all directions. The user's two flagged ones (E, G) are
  the priority; C is queued; F/H are quick follow-ups; B/D are
  optional.

## Open questions

1. **For E (Grassmann)**: which k? `k = 2, 3, 5, 8, 13, 21` per
   the eigenvector pivot plan, or a different scale-aware grid
   for epilepsy? Probably the same — but the question is band-
   specific, so each band may have a different "natural" k.
2. **For G (heat flux)**: which τ regime? At τ → 0 the flux is just
   the adjacency matrix (trivial); at τ → ∞ the flux is uniform
   (also trivial). The interesting regime is intermediate; expect
   τ ∈ [τ_min, τ_min · 10].
3. **For C (eigenmode IPR)**: is the strength-stratified null
   strict enough? An anatomy-stratified null (sample non-epi
   nodes from same Desikan-Killany regions as E_p) would be
   stricter and answer "is this purely anatomical?"
4. **Does it make sense to do Direction A again with different
   ε_KC?** Probably not — the IQR-passing regularities are already
   identified; tuning ε_KC won't add new information. Defer
   indefinitely.
