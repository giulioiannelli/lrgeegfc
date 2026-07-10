---
name: methods-neurophysiological-interpretation
era: IMCOH_ABS_COHORT_N10
status: current
kind: methods
scope: neurophysiological reading of the three-rung probe ladder (raw FC ρ_split, cophenet D_coph ρ_split, Grassmann T_G\*) and per-pair vs global persistence; per-band physical interpretation under the n=10 cohort
sources:
  - .agents/preprint/bands/00_cohort.md (cross-band synthesis, technical)
  - .agents/preprint/bands/01_beta.md … 06_delta.md (per-band locked briefs)
  - .agents/preprint/directives/archive/2026-05/methods_revision_2026-05-18_cophenet.md (cophenet methodology)
  - .agents/preprint/directives/methods_grassmann_cluster_extent.md (Grassmann cluster-mass)
  - data/audit/raw_fc_matched_strength/cohort_summary.csv
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
  - data/audit/grassmann_cluster_extent/cohort_summary.csv
  - data/audit/rho_split_raw_propagator/cohort_summary.csv (raw D(τ_max) ρ_split comparison)
---

# Neurophysiological interpretation — three-rung probe ladder and per-pair vs global persistence

## Head

Across three probes the cohort gives a four-class band map. **β passes everywhere — the only band with full multiscale trace.** α lives at the hierarchical partition level alone, γ_l at the eigenmode-subspace level alone, δ only at the eigenmode-subspace level and with a single-patient caveat. θ and γ_h are null. The two LRG probes (cophenet `D_coph` and Grassmann `T_G\*`) are not redundant — they measure genuinely different aspects of multiscale persistence, and the cohort splits cleanly across the two.

This companion translates the probe geometry into a neurophysiological reading: what each rung is sensitive to, what it can't tell you, what mechanism candidates fit the per-pair vs global dissociation, and the per-band physical interpretation under the locked verdicts.

---

## 1. What each probe measures, and what it can't tell you

### 1.1 Raw FC ρ_split — the edge level

Spearman correlation between Δ_task and Δ_rest, where both Δs are vectorized edge-by-edge over the FC adjacency. This is "did the same contact pairs that the task pulled around stay pulled around in rest_post" at the level of every connectivity edge.

The matched-strength surrogate (4-cycle ±δ rewiring) produces ρ_split clouds with high variance because edge-level Spearman is sensitive to bulk magnitude noise that degree-preserving rewiring doesn't fully eliminate. The signal is there but the test rarely calls it clean.

Cohort (matched-strength C3 gate, n=10, R=200):

  - β +0.258 / p=0.053 / LOOmax=0.102 — borderline
  - α +0.158 / p=0.014 / LOOmax=0.027 — passes nominally
  - γ_l +0.126 / p=0.053 / LOOmax=0.102 — borderline
  - δ +0.111 / p=0.042 / LOOmax=0.082 — fragile under LOO
  - θ / γ_h null

At this rung one band looks decisive (α), three look "almost", two look null. The verdicts feel uncertain because they are.

### 1.2 Cophenet D_coph ρ_split — the per-pair multiscale level

Apply LRG, take the propagator distance `D(τ_max) = 1/ρ̂(τ_max)`, feed to UPGMA hierarchical clustering, then take the cophenetic matrix — the ultrametric distance induced by the merge tree. Compute ρ_split on `D_coph`.

**Why this step is principled in the continuous-spectrum case.** The dendrogram-cophenet sequence does two things at once: (a) it replaces gap-based scale identification (unavailable in our continuous-spectrum LRG outlier regime — see `directives/archive/2026-05/methods_revision_2026-05-18_cophenet.md`) with a merge-tree multiscale identifier, and (b) it discards the absolute magnitude of `D(τ_max)` and keeps only the order in which nodes merge.

The magnitude-discarding aspect is what sharpens the matched-strength contrast. Raw `D(τ_max)` ρ_split for β is +0.220, almost identical to cophenet's +0.221, but the raw-D matched-strength surrogate cloud has variance comparable to the observed signal — paired p only reaches 0.042. Matched-strength rewiring scrambles the merge order dramatically (degree-preserving rewiring is bad at preserving who clusters with whom across the scale ladder); so the surrogate cloud in cophenet space collapses around zero. The same observed median becomes a much tighter rejection: p=0.005 instead of p=0.042.

Cohort (matched-strength C3 gate):

  - α +0.105 / p=0.002 / LOOmax=0.004 — decisive
  - β +0.221 / p=0.005 / LOOmax=0.010 — decisive (both C5 epi-X gates pass)
  - γ_l +0.083 / p=0.116 — null
  - δ +0.008 / p=0.278 — null
  - θ −0.040 / p=0.722 — null (slight anti-alignment)
  - γ_h ≈ 0 / p=0.246 — null

The cophenet step concentrates the claim on α + β at the per-pair multiscale rung.

### 1.3 Grassmann T_G\* — the global eigenmode-subspace level

For each `k` from 2 to 112, take the smallest k non-zero Laplacian eigenvectors. These span a k-dimensional invariant subspace `U_k = span{φ_2, …, φ_{k+1}}` of the Laplacian — the principal global communication modes at progressively finer scales. Measure the chordal distance between the rsPre subspace and the task subspace, vs the task subspace and rsPost subspace; integrate the per-k significant cells into a cluster-mass scalar `T_G\* = Σ_{k: p_k<α} −log₁₀ p_k`.

This is a fundamentally different reading from cophenet. **Cophenet asks: "do specific node pairs preserve their merge order across phases?". Grassmann asks: "does the principal k-dimensional communication subspace rotate consistently from rsPre→task→rsPost?".** The first lives in node-pair space (which contacts co-cluster); the second lives in subspace rotation (the basis spanning dominant communication modes can rotate without any specific pair preserving rank). The two probes are complementary — a band can show one without the other.

Cohort (cluster-extent permutation, audit_70, 6-band refresh 2026-05-26):

  - β `T_G^* = 0.273` (raw 69.76) / `p_mass = 0.005` / LOOmax=0.005 Pat_02 — decisive, both C5 gates pass
  - γ_l `T_G^* = 0.259` (raw 66.14) / `p_mass = 0.005` / LOOmax=0.040 Pat_05 — strong under Decision 12 (LOO precondition passes; C5 epi-X p_mass^epi-X=0.030 secondary)
  - δ `T_G^* = 0.149` (raw 38.07) / `p_mass = 0.005` / LOOmax=0.055 Pat_08 — **weak** under Decision 12 (LOO precondition fails; C5 epi-X strengthens to 0.220 but secondary, not verdict-promoter)
  - α `T_G^* = 0.0312` (raw 7.98) / `p_mass = 0.348` — null
  - γ_h `T_G^* = 0.124` (raw 31.67) / `p_mass = 0.060` — borderline (fails Decision-8 cohort gate)
  - θ null

---

## 2. Per-pair vs global persistence — physical interpretation

**Head.** Per-pair persistence asks "do the same contacts cluster together?". Global persistence asks "does the dominant communication pattern point in the same direction?". A band can keep one and break the other.

### 2.1 Per-pair (cophenet) — discrete, identity-bearing

The cophenet distance `D_coph(i,j)` is defined for every pair of contacts — it's the height at which `i` and `j` first merge in the UPGMA dendrogram built from `D(τ_max)`. Stack all pairs into a vector of length `N(N−1)/2` and take Spearman across phases.

  - **Sensitive to:** the identity of which node groups join up at which scale. If contacts A, B, C cluster together at a coarse scale in rsPre, task, and rsPost, the per-pair Spearman picks it up.
  - **Blind to:** the magnitude of communication along any particular direction — cophenet has thrown the scalar axis away by passing through a ranking step (the dendrogram only remembers merge order).

Persistence at this rung means: **the same contact groupings, at the same relative scales, recur across phases.** A discrete, identity-bearing reading.

### 2.2 Global (Grassmann) — continuous, gauge-invariant

The k-dimensional subspace `U_k` is the "principal hyperplane" on which the dominant communication modes at coarse scale k live, inside the N-dim space of all possible patterns over the contacts. Chordal distance compares two such subspaces as objects in a Grassmann manifold — it measures how much you'd have to *rotate* the rsPre k-plane to align it with the task k-plane. Crucially, this distance is invariant under rotations *inside* the k-plane: you can mix eigenmode 1 with eigenmode 2 freely without moving the subspace.

  - **Sensitive to:** where in N-space the bulk of the communication points.
  - **Blind to:** which specific node combinations span it — two physically different sets of contacts can span the same k-plane.

Persistence at this rung means: **the principal communication subspace stays oriented the same way across phases**, regardless of which specific contacts implement it.

### 2.3 Why they dissociate

A band can preserve "who clusters with whom" without preserving "what direction the bulk communication points in", and vice versa. The two metrics live on genuinely different objects:

  - Cophenet `ρ_split` lives on a vector space of upper-triangular pair-rank entries — every node is named.
  - Grassmann chordal distance lives on `Gr(k, ℝ^N)` — the manifold of k-dimensional subspaces of ℝ^N. Eigenvector identity is washed out by within-subspace rotation; only the subspace orientation matters.

The first is an *identity-preserving* metric; the second is a *gauge-invariant* metric.

### 2.4 One-line analogy

Per-pair persistence: do the same people stay in the same friend groups across time? Global persistence: does the social space the whole group occupies stay oriented the same way, regardless of who-is-friends-with-whom? β keeps both. α keeps the friend groups while the social space drifts. γ_l keeps the social space while the friend groups reshuffle.

---

## 3. Neurophysiological mechanism candidates for global-only persistence

A global-only trace (Grassmann passes, cophenet null) means the network kept the *shape* of how it was communicating, not the *who-with-whom*. The same coarse vibrational geometry persists, but the specific anatomical implementation can re-route. That's the signature of a state — neuromodulatory, top-down, distributed — rather than a memory engram pinning down particular contacts.

### 3.1 Neuromodulatory carryover

Acetylcholine, noradrenaline, dopamine all act by setting region-wide gain/excitability — they shape the principal modes of communication without directly wiring specific A↔B pairs. A task that recruits neuromodulatory release leaves a post-task state in which the modulator's influence on the global mode geometry persists for tens of seconds to minutes, even as specific spike-pair correlations decorrelate. **Global-only persistence is the natural readout of "the cortex came out of task in a different neuromodulatory state".**

### 3.2 Distributed / redundant cortical computation

If a task recruits a many-contact distributed representation, the *function* (which dominant modes carry information) can be invariant across re-implementations using different specific contacts. A motor-planning representation in γ_l could be implemented at one set of pairs during task and a partially overlapping set after, while the *direction in N-space* of dominant activity stays the same. This is the population-code framing — the readout doesn't care which specific neurons fire, only that the population activity points the same way.

### 3.3 Top-down state without bottom-up engram

In predictive-coding language, the global mode geometry is what the cortex is *expecting / modeling*; the per-pair merge structure is the *evidence routing*. A trace at the global level alone is: "the model state survived; the evidence path did not" — the cortex came out of task with the same internal hypothesis about what's going on, but the specific sensorimotor routing through which that hypothesis was being evaluated has decayed or rearranged.

### 3.4 Mechanism candidates for per-pair-only persistence

Contrast: a per-pair-only trace (cophenet passes, Grassmann null — α in our cohort) suggests a *structural scaffold* with phase-dependent *content*. The cleanest reading is a stable thalamocortical / cortico-cortical architecture whose anatomical implementation is fixed by network structure (the same loops, the same connectivity backbone), but what the activity does on that scaffold is task-conditioned. The structural backbone persists; the global communication direction doesn't.

### 3.5 Mechanism candidates for both (β)

When both probes agree, the most parsimonious reading is that the task-driven reorganization sets up a self-sustaining pattern at *both* levels — the specific contact pairings the task pulled together stay pulled together, *and* the principal communication direction of the whole network stays rotated toward the task state. That's the strongest possible signature of multiscale trace.

---

## 4. Per-band neurophysiological reading

### 4.1 β (13–30 Hz, motor / cognitive control)

The only band where all three probes agree, including both C5 epi-X gates (cophenet p=0.003, Grassmann p=0.005). The same multiscale reorganization the task induces is preserved across the dendrogram merge order **and** the principal eigenmode subspace **and** the bulk edge-level Spearman, and stays so when epileptogenic-zone contacts are excised. Median Spearman is identical between raw `D(τ_max)` and cophenet (+0.220 vs +0.221), but only cophenet hits decisive significance.

**Physical reading.** β is the canonical motor / cognitive control band. The cohort shows that whatever task-driven reorganization occurs in β at task_test is then carried into rest_post both in the hierarchical merge order **and** in the principal eigenmode subspace. Robust to Pat_15 anti-alignment and to removing epileptogenic contacts. β is the central trace claim of the manuscript.

### 4.2 α (8–13 Hz, posterior / attentional)

Trace at the cophenet rung (p=0.002), at raw FC (p=0.014), C5 epi-X cophenet passes (p=0.010), but Grassmann is null (p=0.348).

**Physical reading.** In α, the same node groups keep co-clustering across phases — there's a stable partition structure — but the principal α-mode subspace doesn't rotate coherently. Plausibly a thalamocortical α scaffold with a fixed who-clusters-with-whom structure across phases, while the dominant α mode shape stays orthogonal to the across-phase rotation Grassmann would detect. Reads as "α as gating / idling architecture": the architecture is the same, the gating pattern differs.

### 4.3 γ_l (30–80 Hz, local processing + cross-region binding)

Strong trace at Grassmann (`T_G^* = 0.259` raw 66.14, `p_mass = 0.005`, LOO 0.040 Pat_05 passes Decision-12 precondition), C5 epi-X Grassmann secondary (p_mass^epi-X = 0.030, contracts to `T_G^*^epi-X = 0.163` raw 32.75, LOO under epi-X 0.159 Pat_05 fragile — interpreted as signal recruiting cortex straddling epi-zone boundary). Cophenet null (p=0.116). Raw FC borderline (p=0.053).

**Physical reading.** γ_l's trace is at the dominant communication mode angle — the principal subspace rotates coherently across phases — but no specific node groups preserve their merge order in the dendrogram. The principal binding pattern (what the cortex is grouping into coherent percepts/actions) persists, but it's being remapped to slightly different contact pairings. Reads as a binding-state carryover whose anatomical implementation reshuffles patient-to-patient — a binding instruction "co-activate this k-D pattern" survives, while the anatomical instance of who-binds-with-whom drifts.

### 4.4 δ (0.5–4 Hz, slow oscillation)

Grassmann **weak** under Decision 12 (`T_G^* = 0.149` raw 38.07, `p_mass = 0.005` clears Decision-8 gate at floor but full-data LOO 0.055 Pat_08 fails Decision-12 < 0.05 precondition — single-patient leveraged at full data). C5 epi-X strengthens decisively (`T_G^*^epi-X = 0.220` raw 43.99, p_mass^epi-X = 0.005, LOO under epi-X = 0.005 Pat_02 fully robust — secondary mechanistic observation per Decision 10, interpretable as Pat_08 full-data leverage being epi-zone-driven, but never a verdict-promoter). Cophenet null. Raw FC fragile-positive (p=0.042 / LOOmax 0.082).

**Physical reading.** Cohort-level δ Grassmann signal is real but the Wilcoxon is being pulled by one patient; should be reported cautiously per the no-single-patient rule. The C5 pass is reassuring (the δ trace isn't a pure epileptogenic-zone artifact) but the cohort-level claim should note Pat_08-leveraging. If the signal is true cohort-wide, the same global-only mechanism candidates (§3) apply.

### 4.5 θ (4–8 Hz)

Null at every probe. Median ρ_split^coph slightly negative (−0.040). No multiscale trace in this cohort under the matched-strength gate. The cleanest negative reference in the manuscript and a band-specificity benchmark for β and α.

### 4.6 γ_h (80–300 Hz, broadband + HFO/ripple)

Borderline at Grassmann (p=0.060), null at cophenet (p=0.246) and raw (p=0.188). C5 epi-X Grassmann test fails (p=0.060). No conclusive trace; the broadband + HFO regime may be too noise-dominated for the pipeline at this cohort size.

---

## 5. The 2D probe map

|                       | Grassmann passes        | Grassmann null              |
|-----------------------|-------------------------|-----------------------------|
| **cophenet passes**   | **β** — full multiscale trace | **α** — partition-level trace only |
| **cophenet null**     | **γ_l, δ\*** — eigenmode-subspace trace only | θ, γ_h — no trace          |

\*δ with single-patient LOO caveat.

The four-class split is clean and the two LRG probes give genuinely complementary information. β is the only band where the task-driven reorganization persists into rsPost at every angle of the analysis ladder. α and γ_l each tell you something different about how the post-task brain state retains the task imprint: α keeps the partition of who-clusters-with-whom, γ_l keeps the rotation of dominant communication modes. They are not the same neural mechanism, and the probes can tell them apart.

---

## 6. Bigger picture

The per-pair reading tells you what *specific groupings* the brain still recognizes after task; the global reading tells you what *kind of state* the brain is in. They can dissociate because neural computation is implemented redundantly across many specific anatomical paths, and because the state space of cortex is bigger than its anatomical specificity.

A clinically / biologically interesting consequence: if you want to argue that *something specific* about the task is in the brain's post-task state, you need the per-pair rung. If you want to argue that the brain is *in a different state* after the task, the global rung is what you cite. β supports both; α supports the first; γ_l supports the second.

---

## 7. Limitations / caveats

- All claims gated at C3 (matched-strength surrogate, 4-cycle ±δ rewiring, R=200 per cell). C5 (epi-X) gates pass for β at both probes; for α at cophenet; for γ_l at Grassmann; for δ at Grassmann (but with Pat_08 LOO caveat). The α Grassmann null is a robust null at the cluster-mass test, not absence-of-signal at any subspace metric.
- The δ Grassmann verdict is single-patient leveraged (LOO p=0.055 when Pat_08 dropped). Flag in writeups per `feedback_no_single_patient_p_driven.md`.
- Raw FC ρ_split is reported as substrate-rank descriptive only, not as a stand-alone trace claim — the matched-strength surrogate cloud is too loose at this rung to discriminate the candidate alternative.
- The neurophysiological mechanism candidates (§3) are interpretation, not test outcomes. The cohort data is consistent with neuromodulatory / distributed-redundancy / top-down readings of global-only persistence, but the experiment doesn't disambiguate which is operative.
- γ_h spans 80–300 Hz (broadband + lower HFO/ripple regime) — see `memory/brain_bands_definition.md`. The borderline Grassmann verdict in γ_h could be either signal at the cohort size limit or pipeline noise; not resolved at n=10.

---

## Revision history

- **2026-05-26** — Initial methods companion. Consolidates the three-rung probe interpretation, per-pair vs global persistence neurophysiological reading, and per-band physical interpretation. Cross-references `bands/00_cohort.md` §2 (three-layer methodological argument) and the per-band briefs (§4 of each).
