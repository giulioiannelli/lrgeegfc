---
name: sparsification-choice
type: methods-guide
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery → sparsifier SETTLED)
status: SETTLED 2026-07-15 — two-scheme by readout (trace = mst-union@0.20; marker = TMFG). Suite run complete.
created: 2026-07-13
updated: 2026-07-15
scope: PRINCIPLED choice of the graph sparsifier/filter for the LRG diffusion pipeline. Written
  BEFORE committing to stop the empirical bake-off; the full suite (marker + trace, matched-strength,
  per-scale, tmfg/pmfg/disparity vs mst-family) then SETTLED it. The pre-registered rule (§6) and the
  suite outcome (§4a) are the record. Read before touching any backbone in the arc or marker pipelines.
pointers:
  - .agents/reports/2026-07-15_sparsifier-suite-verdict.md               # THE VERDICT (all tables, this decision)
  - .agents/reports/2026-07-13_multiscale-marker-exploration-scope.md    # where TMFG beat mst@0.20 (marker)
  - src/lrg_eegfc/utils/fc/backbone.py                                   # select_backbone / tmfg / pmfg / disparity / mst_union
  - scripts/01_compute/sparsified_arc/13_matched_strength_mst020.py      # trace gate (SA_BACKBONE/SA_FRAC)
  - scripts/01_compute/sparsified_arc/06_epi_arc.py                      # marker (SA_BACKBONE/SA_FRAC)
  - scripts/01_compute/sparsified_arc/25_pmfg_observed_trace_readout.py  # PMFG≈TMFG mechanism confirmation
  - memory/villegas2025_sparsification_grounding_2026_07_11.md           # dense→degenerate, sparsify→multiscale
---

# Choosing the sparsifier — SETTLED: two schemes, one per readout

## Head

**SETTLED (2026-07-15).** The full matched-strength suite disqualified the planar family for the
*cognitive trace* and confirmed it for the *epilepsy marker* — so the pipeline uses **two backbones,
each the validated optimum for its readout, and the dissociation is itself a methodological result**:

- **Cross-phase cognitive trace (§1) → strength-ranked `mst_union_top_fraction` @ frac 0.20.** The
  trace lives in the *strongest* connections. A planar filter that drops them destroys it: under
  matched-strength, α collapses **12/16 → 0/16** scales on TMFG. Strength-preservation is the
  principle; 0.20 is the unique fraction where α+β are *jointly* clean (α needs sparse-strong edges,
  β needs the denser mesoscale skeleton, θ/low_γ stay null).
- **Within-phase epileptogenic-zone marker (§3) → TMFG** (planar-maximally-filtered, Massara 2016).
  Seeded heat-diffusion *detection* sharpens on the triangulated planar skeleton: fused SOZ AUC
  **0.816 → 0.904**, β matched-strength **8/10 → 10/10**, hub patient Pat_15 0.49 → 0.88.

The **dissociation is the finding**, not a compromise: the cross-phase trace is carried by strong-edge
*magnitude*, within-phase detection by the sparse topological *skeleton* — two structural regimes of
the same diffusion operator, each backbone fixed by matched-strength (not by shopping). This flips the
2026-07-13 recommendation (which was PMFG/TMFG-for-everything, decision-pending): planarity is **wrong
for the trace** and only **modestly better for the marker** — see the honest caveats in §4a/§5. The
original recommendation and its five-criteria derivation (§1–§3 below) are kept as the reasoning trail;
§4a records what the suite actually found and why the conclusion moved.

---

## 1. Why we sparsify at all (the principle, not a preference)

Villegas 2025 (PRR; our LRG grounding) makes it non-optional: on a dense, weight-heterogeneous
graph the heat-kernel propagator \(K=e^{-\tau L}\) sits in the Wigner regime — the specific-heat
\(C(\tau)\) has a **single peak**, the geometry is degenerate, and at \(\tau_{\min}\) diffusion is
indistinguishable from the raw weights (audit_174). **Multiscale structure only exists once the
graph is sparsified** (audit_175: the \(C(\tau)\) ladder turns on for cycle-rich backbones at
density ≈ 0.10–0.20; near-tree MST/TMFG-poor fails). So the sparsifier is a **theory requirement**,
and its job is defined: *move the propagator out of the degenerate regime while preserving the
strongest relational structure, with the least arbitrariness.*

## 2. The five criteria (derived from the requirement above)

1. **Sparse** — leave the dense/Wigner regime. Prefer an edge budget that scales O(N), not O(N²).
2. **Cycle-rich / non-trivial** — NOT a tree (a tree is trivially ultrametric — no mesoscale; the
   diffusion needs triangles/cliques). audit_175.
3. **Parameter-free** — no threshold, no fraction, no target edge count that a reviewer can call
   tuned. (This is the criterion `mst@0.20` fails — the 0.20 was never justified.)
4. **Connectivity-preserving** — never fragment; keep the strongest relational skeleton (ideally
   contain the MST).
5. **Established & domain-appropriate** — a named, cited method, ideally already used on
   correlation-like / brain networks.

## 3. Method map

| method | sparse | cycle-rich | param-free | connected | established (correlations) | verdict |
|---|---|---|---|---|---|---|
| global threshold | ✓ | ✓ | ✗ (θ) | ✗ fragments | weak | reject (arbitrary θ) |
| **MST** | ✓✓ | ✗ **tree** | ✓ | ✓ | ✓ | reject (trivial: no cycles) |
| **mst@f (current)** | ✓ | ✓ | ✗ **(f)** | ✓ | ✗ ad-hoc | reject (free parameter f) |
| effective-resistance (Spielman–Srivastava 2011) | ✓ | ✓ | ✗ (target m) | ✓ | ✓ | reject (*approximates dense spectrum → preserves the degeneracy we flee*) |
| **disparity filter** (Serrano 2009 PNAS) | ✓ | ✓ | ✗ **(α)** | ✗ can fragment | ✓✓ *"multiscale backbone"* | **B — principled rival, statistical, has α** |
| **PMFG** (Tumminello 2005 PNAS) | ✓ 3(N−2) | ✓✓ cliques | ✓ | ✓ **⊇ MST** | ✓✓ built for corr. matrices | **A — principled parent** |
| **TMFG** (Massara 2016 J.Complex.Net) | ✓ 3(N−2) | ✓✓ triangulation | ✓ | ✓ | ✓ fast chordal PMFG | **A — fast computation of PMFG** |

## 4a. What the suite found — THE DECISION (supersedes §4)

The full matched-strength suite (trace gate `13` and marker `06` on tmfg/pmfg/disparity vs the
mst family; PMFG-observed readout `25`) **overturned the §4 recommendation**. Planarity is not the
universal principle it looked like — it is *wrong for the trace* and only *modestly better for the
marker*. The evidence, per-scale, matched-strength:

**Trace (cophenetic ρ_sym gate, clears / 16 scales — want α,β CLEAR and θ,low_γ NULL):**

| band | mst@0.20 | mst@0.10 | mst@0.05 | disparity@0.2 | **TMFG** | PMFG (obs ρ) |
|---|---|---|---|---|---|---|
| **alpha** | **12/16** | 8/16 | 15/16 | 11/16 | **0/16 ✗** | ρ 0.11 ✗ |
| **beta**  | **16/16** | 9/16 | 5/16 | 11/16 | 3/16 ✗ | ρ 0.18 ✗ |
| theta (null) | 0/16 ✓ | 0/16 | 0/16 | 1/16 | 3/16 ✗ | — |
| low_γ (null) | 0/16 ✓ | 9/16 ✗ | 1/16 | 8/16 ✗ | 7/16 ✗ | — |

- **Planarity discards the strong edges the trace lives in.** At *identical* density 0.05, α ρ_sym
  is 0.24 for strength-ranked mst but only **0.11 for both PMFG (exact) and TMFG (greedy)** — halved.
  So it is the *planar family*, not TMFG's shortcut (PMFG≈TMFG; §6-iii resolved). It is **not**
  sparsity: mst@0.05 at the same density keeps α fully.
- **mst@0.20 is the unique clean fraction.** α and β have opposite density needs (α peaks sparse
  15/16@0.05; β needs density 16/16@0.20, dies 5/16@0.05). Only 0.20 gives α strong + β maximal +
  θ/low_γ null. Sparser mst leaks low_γ; disparity/planar leak low_γ *and* lose specificity.
- → **Trace backbone = `mst_union_top_fraction` @ 0.20.** Pre-registered outcome "keep mst@0.20"
  (§6): no alternative beats it on the trace.

**Marker (SOZ seeded-diffusion, matched-strength beats-null of 10; fused AUC @ s=10):**

| backbone | fused AUC | R-prec | β beats-null | low_γ | high_γ |
|---|---|---|---|---|---|
| mst@0.20 | 0.816 | 0.409 | 8/10 | 8/10 | 7/10 |
| mst@0.10 | 0.859 | 0.467 | 8/10 | 9/10 | 7/10 |
| mst@0.05 | 0.874 | 0.469 | 8/10 | 8/10 | 7/10 |
| disparity@0.2 | — | — | 8/10 | 9/10 | 6/10 |
| **TMFG** | **0.904** | **0.550** | **10/10** | 9/10 | 4/10 ↓ |

- Fused AUC is **monotone in sparsity within the mst family** (0.742@0.50 → 0.874@0.05): most of
  the marker gain is just "go sparser," which mst@0.05 already captures. **TMFG adds a further +0.03
  AUC / +0.08 R-prec / β 10-of-10** over mst@0.05 *at the same density* — a real planar-specific
  detection edge, at the cost of **losing high_γ (4/10)**.
- → **Marker backbone = TMFG** (user decision 2026-07-15: for *detection* we optimise detection).
  The strength-ranked alternative mst@0.05 (0.874) is recorded so the choice is transparent, not a
  hidden knob.

**Settled pipeline:** trace/cognition → `mst_union_top_fraction(W, 0.20)`; marker/epilepsy →
`tmfg_backbone(W)`. Both via `select_backbone(W, kind)` and the `SA_BACKBONE` env. The two-backbone
split is defensible **because each is matched-strength-validated as its readout's optimum**, and the
magnitude-vs-skeleton dissociation is a reportable methods result — *not* because we shopped filters.

## 4. Recommendation as first written — the PMFG/TMFG family (choice A) — SUPERSEDED by §4a

> Kept as the reasoning trail. The suite (§4a) showed this "planar for everything" recommendation
> was right for the marker but wrong for the trace. Read §4a for the decision.

Why it is principled, not a cherry-pick:
- **Parameter-free by construction.** Planarity fixes |E| = 3(N−2); density = 6(N−2)/[N(N−1)] ≈
  **6/N** → auto-sparsifies as N grows, always in the sparse regime. No knob.
- **A principled enrichment of the tree.** PMFG **provably contains the MST**, then adds exactly the
  3- and 4-cliques (triangles/tetrahedra) that a tree lacks — precisely the cycles audit_175 says
  the diffusion needs. It is "MST + the minimal non-trivial structure," derived, not tuned.
- **Chordal clique-tree = multiscale by construction (TMFG).** TMFG is a triangulation → a
  decomposable graph with a nested clique tree; that hierarchy is inherently multiscale, aligned
  with the LRG reading.
- **Domain-appropriate & established.** PMFG/TMFG were invented for **financial correlation
  matrices** — structurally identical to FC (dense, weighted, heterogeneous) — and PMFG has been
  applied to brain functional networks. Real citations, not invented.
- **Performance is confirmation, not the reason.** On the marker: fused AUC 0.816→0.904, β
  matched-strength 8/10→**10/10** @0.87, α 0.61→0.79, hub patient Pat_15 0.49→0.88 — and the method
  ranking is **monotone in the principle** (dense worst → tree-ish poor → cycle-rich sparse best),
  the signature of a mechanism, not a lucky filter.

## 5. The honest caveat and why we still prefer A

PMFG/TMFG imposes **planarity** (genus-0 embedding), inherited from econophysics and only weakly
motivated for brains. The only method that drops the topological prior while staying principled is
the **disparity filter** (B) — purely statistical, multiscale-by-design — but it carries the α
parameter (violating criterion 3) and can fragment the graph (violating 4). Given your priorities
(no parameters, non-trivial, principled, established), **A wins on the criteria**; B is the
**cross-check** that tells us the A-result is not a planarity artifact. If A and B agree
directionally, the finding is "principled sparse backbone," reviewer-proof.

## 6. Pre-registered adoption rule (anti-cherry-pick — decide the rule BEFORE seeing the suite)

Run the marker **and every arc analysis** on the PMFG/TMFG suite, matched-strength, per-scale.
Then:

- **Adopt PMFG/TMFG globally (single-backbone paper)** IFF: (i) the **marker** improves or holds vs
  mst@0.20 under matched-strength (already indicated), **and** (ii) the **cross-phase α/β trace
  holds-or-improves** under matched-strength, read **per-scale not best-scale** (β stays ≈16/16, α
  ≈12/16 or better; δ/high-γ remain un-claimed collapse artifacts), **and** (iii) **PMFG ≈ TMFG**
  (the win is the planar-filter principle, not TMFG's greedy shortcut), **and** (iv) the **disparity
  cross-check** agrees directionally.
- **Adopt for the marker only (two-backbone paper)** IF the marker passes but the trace **degrades**
  on TMFG — defensible: the directive already says the marker "reads the propagator differently."
- **Keep mst@0.20** IF neither improves, or if PMFG≠TMFG (the effect is a chordal-approximation
  artifact) — then the current pipeline stands and the marker exploration is logged as negative.

Do **not** relax this rule after seeing the numbers. "Better" = the matched-strength, per-scale
criteria above — never a single best-scale scalar or an AUC that only wins in-sample.

**OUTCOME (2026-07-15, rule applied to the suite):**
- (i) marker improves on TMFG ✓ (fused AUC 0.816→0.904, β 8/10→10/10).
- (ii) trace holds-or-improves per-scale on TMFG ✗✗ — **α 12/16→0/16, β 16/16→3/16**, θ/low_γ
  false-positives. The planar family *fails* (ii); PMFG≈TMFG confirms it is planarity, not the greedy
  shortcut, so (iii) is moot for global adoption.
- ⇒ "Adopt PMFG/TMFG globally" is **rejected**. The pre-registered branches give **keep mst@0.20 for
  the trace** (no alternative beats it) and license **marker-only TMFG** (marker passes, trace
  degrades) — the two-backbone paper. The user elected the two-backbone route (2026-07-15), so the
  marker uses TMFG and the trace stays mst@0.20. This is the pre-registered "marker-only" branch, not
  a post-hoc relaxation. The declined alternative (single-family mst@0.05 marker, AUC 0.874) is logged
  in §4a for transparency.

## 7. What the suite run needed — DONE (2026-07-15)

Library (`src/lrg_eegfc/utils/fc/backbone.py` — general names, all in place):
- `pmfg_backbone(W)` — exact planar-maximal filter (networkx `check_planarity`, contains MST,
  3N−6 edges, ~9 s/graph → observed-only, not inside a null).
- `disparity_backbone(W, alpha, ensure_connected=True)` — Serrano per-edge p-value, MST-union to
  stay spanning; fragments natively at α≤0.1, connected/cycle-rich only at α≈0.2.
- `select_backbone(W, kind, frac, disparity_alpha)` — one dispatcher for all pipelines.

Pipelines re-pointed via `SA_BACKBONE` (+ `SA_FRAC`, `SA_DISP_ALPHA`) envs:
- `06_epi_arc.py` (marker) and `13_matched_strength_mst020.py` (trace gate) — **run** across
  tmfg / disparity / pmfg / mst@{0.05,0.10,0.20}; per-backbone OUT dirs under `data/sparsified_arc/`.
- `25_pmfg_observed_trace_readout.py` — PMFG≈TMFG observed-ρ_sym mechanism check (**run**).
- **Not yet re-pointed (only needed to *propagate* the settled marker, not to decide it):**
  `05_enc_inf_arc.py`, `12_sparsification_recovery.py`, `14_controls_ladder_mst020.py`,
  `15_patient_variation.py`. The trace stays mst@0.20, so these downstream trace analyses need **no**
  change. Only the §3 marker artifacts move to TMFG (see §7a).

## 7a. Propagating the settled marker to §3 (follow-on writeup task)

The trace side is unchanged (mst@0.20) → §1 and all trace downstream (localization, laterality,
enc/inf, reinstatement) need nothing. Only the **epilepsy marker (§3)** moves mst@0.20 → **TMFG**:
- Marker data already exist: `data/sparsified_arc/epi_arc_tmfg/` (β 0.869/**10-of-10**, low_γ
  0.843/9, δ 0.830/7, α 0.792/7, θ 0.743/6, high_γ 0.688/4; fused AUC 0.904).
- Re-source `scripts/01_compute/figures_embedded/fig_epi_a_relational_marker.py` from
  `epi_arc_mst020/per_cell.csv` → `epi_arc_tmfg/per_cell.csv`; update `results_sec_3.tex` marker
  AUCs + the fig:epi1(b) caption; keep the distant-seed/detector paragraphs consistent.
- Add the **dissociation** sentence to Methods/§3: the trace uses the strength-preserving
  mst backbone, the marker the planar TMFG backbone, and the split reflects magnitude- vs
  skeleton-carried structure (see §4a).

## 8. Do-not (carry the honesty rules)

- Matched-strength is the only null; per-scale, never best-scale (Villegas-collapse δ/high-γ stay
  un-claimed). Don't re-attribute the Grassmann (dense/whole-task) to any sparsified backbone.
- Don't adopt on the marker win alone — the trace must be checked (its mst@0.20 was chosen because
  **α needs a strong-edge backbone**; TMFG keeps the strongest edges + triangles, so it *should*
  preserve α, but VERIFY — see the trace handoff).
- Don't invent a filter — PMFG/TMFG/disparity are named literature methods; anything else needs a
  citation before use.
