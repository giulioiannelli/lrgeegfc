---
name: pairwise-descriptor-ladder
type: scope
era: IMCOH_ABS / COHORT_N10
status: result
created: 2026-07-09
updated: 2026-07-09
pointers:
  - scripts/01_compute/audit/audit_166_pairwise_descriptor_ladder.py
  - src/lrg_eegfc/utils/metrics/graph_descriptors.py
  - data/audit/pairwise_descriptor_ladder/cohort_summary.csv
  - scripts/01_compute/audit/audit_150_rho_sym_gate.py
  - scripts/01_compute/audit/audit_152_consolidation_arc_rhosym.py
  - src/lrg_eegfc/utils/surrogate/matched_strength.py
  - .agents/reports/2026-07-06_raw-vs-cophenet-resolving-power.md
  - .agents/guides/task-persistence-investigation/2026-07-09_grassmann-inference-arc.md
---

# Pairwise-descriptor ladder — do the usual low-level network scalars reproduce the multiscale trace?

**Head.** Our discriminating task-trace results — band-selective whole-task
consolidation (α/β clear, θ silent), a **β-only** inference-specific refinement
(`T_infspec·e` p=0.0098), and β→OFC localization — currently exist **only** after
the LRG cophenetic transform. The obvious referee worry, and the user's, is that
this ornate machinery merely re-encodes something a first-year network-science
scalar already sees: **node strength, weighted clustering, a centrality**. This
scope defines a **ladder of classical pairwise descriptors run through the
identical ρ_sym estimator and the identical matched-strength null**, ordered by
structural order (raw edges → node scalars → cophenetic), scored on the four
signatures that make the cophenetic result a result. The deliverable is a
**verdict table**: which rungs reproduce the discrimination, which are null, which
only look positive until matched-strength strips the strength out. The claim we
are trying to *break* is "you need the multiscale comparison"; the ladder is built
to break it if it is breakable.

Pure reuse of the ρ_sym + MS harness; the only new code is a `graph_descriptors`
library module (there is currently **no** classical-descriptor primitive in the
codebase — strength lives inline as `W.sum(1)`) and the swap-the-representation
audit driver.

This is a **result-first** investigation. The yes/no: does any low-level pairwise
descriptor reproduce the cophenetic discrimination (band-selectivity **and**
inference-β-specificity **and** →OFC), and does any apparent hit survive
matched-strength?

> **VERDICT (2026-07-09, audit_166). The named disaster is AVERTED; the unique
> multiscale contribution is BAND-SELECTIVITY, not inference-visibility.**
> Cophenetic cross-check PASS (α/β T_test .032/.032, θ silent .784, β-only
> `T_infspec·e` .0098 — reproduces audit_150/152 exactly, harness validated).
> Ladder (gate p, whole-task T_test | β `T_infspec·e`):
> - **Local classical scalars are BLIND.** `strength` p=1.000 all bands
>   (degenerate under MS **by construction** — MS preserves strength exactly, so
>   the strength-trace is reproduced and cannot clear the null; the cleanest
>   possible "strength is NOT causing our results"). `clustering_onnela` (0.99
>   strength-slaved) clears nothing, β inference .161. `pagerank` (1.00 slaved)
>   nothing, β inference .920. → the *usual quantities the user named* recover none
>   of the trace.
> - **Band-selectivity is UNIQUE to cophenetic.** Only the cophenetic rung recovers
>   α (.032) while staying silent on the null band θ (.784). Raw FC is BROAD (≈.05
>   in every band, θ .053 not silent); `eigcent` is β-only (misses α). No pairwise
>   scalar reproduces the α-recovery + θ-silence signature. This is the defensible
>   "you need the multiscale comparison" (band-selectivity IS the LRG contribution,
>   [[cophenet_methodology_rationale]]).
> - **HONESTY FLAG — the β inference-specific trace is a CONVERGENCE phenomenon,
>   NOT exclusively multiscale.** Reproduced by raw FC edges (.010, 8/10),
>   eigenvector centrality (.002; thin — 9-10/10 signed but 3/10 individually), and
>   closeness (.024), as well as cophenetic (.010, 6/10). β is the convergence cell
>   where raw and multiscale agree (ρ 0.85–0.95, [[feedback_dP_framing]]). So "only
>   the multiscale comparison sees the inference-specific trace" is an
>   OVERSTATEMENT. The inference trace is **fine-grained** (edge + hierarchy level),
>   not coarse-subspace (Grassmann null, audit_165) — that is the reconciled
>   statement. This **partly walks back the N2.7 "specifically multiscale" framing**
>   folded into `headlines/02_encoding_vs_inference.md` from the Grassmann result;
>   flagged for user reframing, not silently rewritten.
> - **Gradient (the real story):** local scalars → blind; fine/global probes (raw
>   edges, eigcent, closeness) → see the β convergence signal; multiscale → adds
>   band-selectivity + (downstream) β→OFC. Result memory:
>   [[pairwise_descriptor_ladder_2026_07_09]].
> - **DRIFT ADDENDUM (2026-07-10, WITHDRAWN 2026-07-10).** I briefly read the
>   raw/eigcent/closeness inference detection as session DRIFT via a windowed drift null
>   (audit_167/168). **Withdrawn (user-corrected):** a windowed drift arc is NOT a valid
>   null for the CONDITIONAL functional `T_infspec·e` (it conditions on encoding `e`);
>   the sham arc, built entirely from pre-task rest where no signal can exist, returns a
>   *significantly positive* conditional trace on its own (β drift>0 p=0.007) — the null
>   is broken, not the trace. The "convergence phenomenon" stands as a **matched-strength**
>   finding (part-3 above), not drift. Whole-task β IS drift-clean (locked **C2** 0.0137 +
>   audit_167 0.042); α construction-dependent — those whole-task statements remain valid.
>   Rule: conditional trace → validate with MS, not a drift window. See
>   `.agents/preprint/supplementary/S2_drift_controls.md`, [[feedback_drift_null_mandatory]].

---

## Critical preamble (5 points, before code)

**(1) Claim.** No classical low-level pairwise-network descriptor, turned into the
*same* ρ_sym rank correlator under the *same* matched-strength null, reproduces the
discriminating signature of the cophenetic trace. Concretely, per rung we test the
four cophenetic signatures: (a) whole-task band-selectivity — recovers the α/β
cohort gate while staying silent on the clean-null band θ; (b) the **β-only**
inference-specific refinement `T_infspec·e`; (c) survival of matched-strength;
(d) β→OFC localization (only evaluated for a rung that clears (a)–(c) in β). The
strong claim: **only the cophenetic rung passes all four**; every classical rung is
either null or an artifact that matched-strength kills.

**(2) Null.** Per rung, per (patient, band): matched-strength 4-cycle ±δ
strength-preserving surrogate (R=200, `SWAP_FACTOR=20`, `W_MAX=1.0`), regenerated
per phase with a numba shuffle (bit-identical to the pure-Python engine,
`strength_preserving_shuffle`), fresh seed 20260709 (independent per-cell rng =
`BASE_SEED + cell_idx`). **Every rung is computed on the SAME regenerated
surrogate adjacencies**, so the ladder is internally apples-to-apples. Surrogate
p = upper tail `mean(surr ≥ obs)`; cohort gate = one-sided Wilcoxon(obs − surr_p50)
`alternative='greater'` (audit_150 `cohort_gate`). Matched-strength is the minimum
mandatory null for any FC cohort claim (matched-strength rule).

**(3) Strongest plausible alternatives the null should control for.**
  (i) **Strength-slaving.** The single most dangerous alternative: our whole story
      is an elaborate re-encoding of **node strength**, the one quantity the MS
      null preserves exactly. Weighted clustering ≈ strength at r≈0.96 on dense FC;
      several centralities are strongly strength-slaved. If the trace were a
      strength effect, a strength-slaved descriptor would carry it AND MS would
      kill it.
  (ii) **A non-strength low-level descriptor genuinely reproduces the trace.** A
      centrality / closeness that encodes real topology beyond strength could both
      survive MS and localize to OFC — this would (partially) *falsify* the
      "multiscale is necessary" claim. This is why we run the FULL menu, not just
      strength.
  (iii) **Presence ≠ discrimination.** Raw FC is *positive in all six bands*
      (audit_67); a rung can show "a trace" while failing to be band-selective or
      inference-specific. The verdict must be scored on discrimination, not on the
      mere existence of a positive ρ.

**(4) Does the null actually control for it? (the hinge)**
  - Alternative (i) is controlled **by construction**: MS preserves each node's
    strength to float precision, so **node strength itself is degenerate under the
    null** — the surrogate reproduces the observed strength-trace exactly (the
    per-phase strength vector is identical to observed for every surrogate), giving
    p≈1. Node strength is *definitionally incapable* of producing an
    MS-surviving trace. Any descriptor that is a deterministic function of strength
    inherits this: its trace is reproduced by the surrogate and cannot clear MS.
    This is the cleanest possible answer to "is strength already causing our
    results — no."
  - Alternative (ii) is **not** auto-covered by MS — a non-strength descriptor can
    survive. It is covered by *running it*: the full menu gives every simple
    pairwise quantity a fair, identical chance. The design can therefore falsify us.
  - Alternative (iii) is covered by the **four-criterion score** (band-selectivity
    + inference-specificity are discrimination tests, not presence tests) and by
    including the raw-FC rung as the explicit "broad, non-selective positive"
    reference.

**(5) Falsification + limitations.** The claim is **falsified** if any classical
rung (a) survives MS, (b) is band-selective (recovers α, silent on θ), AND
(c) reproduces the β-only inference-specific refinement — or, weaker, if a rung
survives MS in β and localizes to OFC. Any such rung must be reported as a pairwise
scalar that (partially) recovers the trace, softening "LRG is necessary."
**Limitations:** node-level ρ_sym is a Spearman over N≈100–130 nodes; descriptors
with heavy ties (integer k-core, degree on a fully-connected graph) are
ρ_sym-degenerate and are **excluded and named, not silently dropped**. Descriptors
that are near-constant on a dense graph (local efficiency: every node's
neighborhood is the whole graph) are degenerate and excluded. Global summary
scalars (one value/phase) cannot enter the ρ_sym ladder and are relegated to a
weaker per-patient-displacement tier (Tier B). n=10; the cohort Wilcoxon p-grid is
coarse.

---

## Notation

- `W_x ∈ ℝ^{N×N}` — `imcoh_abs` FC adjacency for phase `x` (zero-diag, clipped
  `[0,1]`, symmetrized), via `load_phase` (audit_150:62-68).
- `s_i(W) = Σ_j W_ij` — node strength (the MS-preserved quantity).
- Phases: `rest_pre` split into `A`, `B` (unbiased baseline halves,
  `ensure_half_fcs`); `task_learn` (TL); `task_test` (TT); `rest_post` (RP).
- A **representation** is a map `φ : W ↦ v ∈ ℝ^m` producing a per-phase vector on a
  common index set: pair-level (`m = N(N−1)/2`, raw edges / cophenetic) or
  node-level (`m = N`, the classical scalars).
- **ρ_sym estimator** (audit_150:115-119, copied verbatim, representation-agnostic):
  > `ρ_sym(vA, vB, vT, vP) = ½[ Spearman(vT−vA, vP−vB) + Spearman(vT−vB, vP−vA) ]`

  Sign convention **positive = trace** (locked `T_d>0`). It consumes any four
  equal-length vectors; only `φ` changes across rungs.

## Definitions — the rungs

Each rung is a representation `φ`. The estimator, null, and cohort test are held
identical (audit_150 / audit_152 machinery).

| rung | `φ(W)` | level | note |
|---|---|---|---|
| **R0 raw FC** | `triu(W)` edge weights | pair | bottom rung; ρ_sym re-run of audit_67 (never done under ρ_sym) |
| R1 strength | `s_i = Σ_j W_ij` | node | **MS-degenerate anchor** (p≈1 by construction) |
| R2 clustering (Onnela) | `Ĉ_i = diag(Â³)_i /(k(k−1))`, `Â=(W/W_max)^{1/3}`, `k=N−1` | node | strength-slaving showcase (r≈0.96) |
| R3 eigenvector centrality | leading eigenvector of `W` (abs) | node | |
| R4 weighted PageRank | stationary of `α D⁻¹W + (1−α)/N` | node | α=0.85 |
| R5 closeness centrality | `(N−1)/Σ_j d_ij`, `d = ` shortest path on `1/W` | node | integration scalar |
| R6 betweenness centrality | Brandes on `1/W` distances | node | **conditional**: include only if per-graph cost < ~30 ms; else exclude-and-log |
| **R_top cophenetic** | `cophenetic_condensed_from_adjacency(W)` (τ=1/λ_max) | pair | flagship; library helper, no fork of `ultra` |

**Excluded, named (not silently dropped):** weighted degree / k-core (ties →
ρ_sym-degenerate on a fully-connected graph); local efficiency (near-constant on a
dense graph → degenerate); weighted modularity Q (partition-stochastic, LRG-adjacent).

**Two targets** (both computed per rung, mirroring audit_150 + audit_152):
- **Target 1 — whole-task trace** `T_test = ρ_sym(φ(A), φ(B), φ(TT), φ(RP))`, per
  band. The band-selectivity score (recovers α, silent on θ).
- **Target 2 — the consolidation arc** (audit_152 functionals, per band):
  `e = φ(TL)−φ(A)`, `g = φ(TT)−φ(A)`, `f = φ(TT)−φ(TL)`, `p = φ(RP)−φ(B)`,
  symmetrized over arms; `T_learn = ρ(e,p)` (encoding), `T_infspec = ρ(f,p)`,
  **`T_infspec·e = partial_ρ(f, p | e)`** (the β-only headline). The
  inference-specificity score.

## Scoring rubric (per rung)

Four discrimination criteria + one diagnostic:

| criterion | cophenetic (target to beat) | test |
|---|---|---|
| C1 recovers α | α gate p=0.024 | whole-task Wilcoxon gate at α < 0.05 |
| C2 silent on θ | θ gate p=0.78 (clean null) | θ gate p **not** < 0.05 |
| C3 survives MS | α/β clear | any band clears the MS gate |
| C4 β→OFC | q 0.010–0.015 (R=1000) | localization, only if C1–C3 pass in β |
| C5 inference-β-only | `T_infspec·e` β p=0.0098, rest > 0.24 | β < 0.05 AND all other bands null |
| **DIAG strength-slaving** | — | cohort-median Spearman(φ(W), s(W)) per phase; PLUS the MS-survival column |

**Verdict logic.** Expected: R1 MS-degenerate (p≈1); R2–R6 null or MS-killed
(strength-slaved); R0 broad/non-selective positive (recovers α like audit_67 but
fails C2/C5); only R_top passes C1–C5 + C4. If any classical rung passes C3 **and**
(C4 or C5), report it honestly — the falsification path is open on purpose.

## Properties / sanity contracts

- **Cross-check gate (anti-hallucination).** R_top (cophenetic) MUST reproduce the
  audit_150 **verdict** (α/β CLEAR, θ fail) and the audit_152 β-only
  `T_infspec·e`; seed differs from the flagship (20260709 vs 20260706/20260511) so
  this reproduces the *verdict*, not the bit value (audit_150 preamble pt 4). Script
  asserts α and β clear for R_top before any classical number is trusted.
- **R1 degeneracy contract.** Under MS the per-phase strength vector equals the
  observed to float precision, so `surr_rho[strength] ≈ obs_rho[strength]` for every
  surrogate ⇒ p≈1 (up to float noise). If R1 shows a "significant" trace, the
  reconstruction or the null is mis-wired — R1 is a wiring test.
- **Representation-length invariance.** ρ_sym is a Spearman; it is invariant to the
  entry count. Comparability across rungs is on the four criteria (p-value pattern),
  never on raw ρ magnitude (stated to preempt "apples vs oranges").
- **Sign convention locked:** positive = trace, all functionals.
- **Complexity.** Dense-graph descriptors reduce to numpy/scipy linear algebra:
  strength O(N²); Onnela `Â³` O(N³); eigenvector-centrality `eigh(W)` O(N³);
  pagerank sparse-solve O(N³); closeness Floyd-Warshall O(N³). All ~ms at N≈130.
  Betweenness (Brandes) is the only super-cubic term ⇒ conditional. Hot loop:
  `for r in R: {numba shuffle 4 phases; descriptor bank; ρ_sym}`.

## Caveats & failure modes

| # | Failure mode | Mitigation (in-script) |
|---|---|---|
| (i) | descriptor is a pure function of strength → MS-degenerate | that IS the finding; DIAG reports the strength correlation + MS-survival to make it explicit, not a bug |
| (ii) | a non-strength descriptor survives MS + localizes | run the full menu; report honestly (falsification path) |
| (iii) | ρ_sym-degenerate descriptor (ties / near-constant on dense graph) | exclude + name (k-core, local efficiency, degree); never silently drop |
| (iv) | presence mistaken for discrimination | four-criterion score; raw-FC rung as the broad-positive reference |
| (v) | reconstruction float noise perturbs strength | Mode A regenerates W directly (numba shuffle) — no eig round-trip; strengths exact |
| (vi) | betweenness cost blows up runtime | time one graph first; include only if < ~30 ms, else exclude-and-log (no silent cap) |

## Pseudocode

```
BANDS=6; COHORT=10; R=200; SEED=20260709; PHASES=[A,B,TL,TT,RP]
DESC = {raw_fc, strength, clustering_onnela, eigcent, pagerank, closeness, [betweenness], cophenetic}
for (idx, pat, band) in cells:
    W = {x: load_phase(pat, x, band) for x in PHASES}
    # observed
    V = {d: {x: phi_d(W[x]) for x in PHASES} for d in DESC}
    for d in DESC:
        obs.T_test[d]      = rho_sym(V[d][A], V[d][B], V[d][TT], V[d][RP])
        obs.arc[d]         = sym_functionals(V[d])      # T_learn, T_infspec, T_infspec_pe
        diag.slave[d]      = median_x Spearman(V[d][x], strength(W[x]))
    # surrogate (Mode A: regenerate all phases, same draws for every descriptor)
    rng = default_rng(SEED + idx)
    for r in 1..R:
        Ws = {x: numba_shuffle(W[x], n_swaps, rng) for x in PHASES}   # 4-cycle ±δ
        Vs = {d: {x: phi_d(Ws[x]) for x in PHASES} for d in DESC}
        for d: surr.T_test[d][r] = rho_sym(...);  surr.arc[d][r] = sym_functionals(Vs[d])
    for d: p.T_test[d]   = mean(surr.T_test[d] >= obs.T_test[d])          # upper tail
           p.arc[d]      = mean(surr.arc[d] >= obs.arc[d])
cohort:
    for d, band: gate_p[d,band] = wilcoxon(obs - surr_p50, greater)      # C1..C3
    assert R_top clears alpha and beta, fails theta                       # cross-check gate
    score each rung on C1..C5 + DIAG → verdict table
report: data/audit/pairwise_descriptor_ladder/{per_patient,cohort_summary,strength_slaving}.csv + README
```

## Visualization spec

- **Ladder verdict heatmap** (headline): rows = rungs (R0…R_top, structural order
  bottom→top), columns = the four criteria + the six-band gate; cell = pass/fail or
  −log10 p. Reading rule: a solid bottom band of failures with only the top row
  (cophenetic) lit across C1–C5 = "only the multiscale comparison discriminates."
- **Strength-slaving scatter** (diagnostic): per descriptor, cohort-median
  Spearman(φ, strength) on the x-axis vs MS-survival (does the trace clear) on the
  y-axis; annotate clustering at r≈0.96. Reading rule: descriptors cluster at high
  strength-correlation + no-MS-survival (bottom-right) except cophenetic.
- **Band profile lines**: `−log10 gate_p` vs band for R0 (raw, broad), R_top
  (cophenetic, α/β-peaked), and one representative classical scalar (flat/null).
- `use_lrg_style`, PDF-only, transparent, no pre-registered gate lines; band palette
  from `visuals.styles`.

## Connection to prior tools

Sibling of the **Grassmann inference arc** (2026-07-09): both ask "does a
*non-cophenetic* probe reproduce the trace." Grassmann probes coarse subspace
geometry and answered NO for the inference-specific component (audit_165); this
probes the *classical low-level* scalars and asks the same of the whole-task and
inference targets. Together they bracket the cophenetic result from above
(spectral subspace) and below (pairwise scalars). Reuses the ρ_sym gate
(audit_150), the consolidation arc (audit_152), the matched-strength engine
(`matched_strength.py`), and the raw-FC baseline (audit_67, now re-run under
ρ_sym as R0). Does **not** replace the cophenetic result — it is the negative
control fence around it.

## Implementation plan

- **Library (new, general names):** `src/lrg_eegfc/utils/metrics/graph_descriptors.py`
  — `node_strength(W)` (unifies scattered `W.sum(1)`, 3+ callers), plus
  `weighted_clustering_onnela`, `eigenvector_centrality`, `pagerank`,
  `closeness_centrality`, `raw_edges` (triu), `[betweenness_centrality]`. Pure
  numpy/scipy, no networkx in the hot loop. Names are general graph concepts (never
  manuscript-local). R_top reuses `cophenetic_condensed_from_adjacency` (library) —
  no `ultra` fork.
- **Audit driver:** `scripts/01_compute/audit/audit_166_pairwise_descriptor_ladder.py`
  (next free index — confirm). Clone audit_150 skeleton (`load_phase`, numba
  `_swap_loop`/`shuffle`, `cohort_gate`, MP Pool, live `[i/N]` progress); add the
  descriptor bank and the audit_152 arc functionals (`_arm_functionals`,
  `_partial_rho`). Mode A regeneration, seed 20260709.
- **Compute estimate:** time 1–2 cells before full launch (optimize-and-surface
  rule); extrapolate; numba warm-compiled in parent before fork. Expect a few
  minutes at R=200 × 60 cells with the pure-numpy descriptor bank; betweenness is
  the only risk and is gated on a timing check.
- **Outputs:** `data/audit/pairwise_descriptor_ladder/{per_patient_per_band,
  cohort_summary,strength_slaving}.csv` + `README.md`; verdict figure PDF deferred
  until numbers are read with the user.

## Open questions

- **Raw-FC arc.** R0 gets the whole-task ρ_sym re-run for free; whether to also
  compute the raw-FC *inference-specific* arc (T_infspec·e on edges) — expected
  null/broad — is a cheap add and included by default (it is the honest bottom-rung
  comparison for C5).
- **Betweenness inclusion** decided at runtime by the timing gate; if excluded,
  logged in the README (no silent truncation).
- **Global Tier B** (efficiency, char-path, transitivity, assortativity as
  per-patient displacement Wilcoxon) is a secondary add; build only if the node-level
  ladder leaves the "even the whole-graph summary" question open.
- No behavior (TC1): the verdict is about *representational* recoverability, not
  cognitive linkage.
