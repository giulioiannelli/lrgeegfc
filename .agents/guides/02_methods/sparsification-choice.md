---
name: sparsification-choice
type: methods-guide
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery → sparsifier-decision pending)
status: decision-pending (recommendation = PMFG/TMFG family; to be settled by the full-suite run)
created: 2026-07-13
scope: PRINCIPLED choice of the graph sparsifier/filter for the LRG diffusion pipeline. Written
  BEFORE committing, to stop the empirical bake-off (pick-the-AUC-winner = cherry-picking the
  method) and replace it with a principle-first decision + a pre-registered adoption rule + a
  de-risk protocol. Read this before touching any backbone in the arc or marker pipelines.
pointers:
  - .agents/reports/2026-07-13_multiscale-marker-exploration-scope.md   # where TMFG beat mst@0.20 (marker)
  - .agents/reports/2026-07-13_TMFG-backbone-for-trace-HANDOFF.md        # does the TRACE improve on TMFG?
  - src/lrg_eegfc/utils/fc/backbone.py                                   # tmfg_backbone / mst_union_top_fraction / percolation
  - memory/villegas2025_sparsification_grounding_2026_07_11.md           # dense→degenerate, sparsify→multiscale
  - memory/audit_175_sparse_backbone_verdict_2026_07_11.md               # cycle-rich d≈0.10–0.20 gives the ladder
  - memory/audit_174_diffusion_vs_rawfc_gate_2026_07_11.md               # fully-connected propagator degenerate at τ_min
---

# Choosing the sparsifier — principle first, performance only to confirm

## Head

We must pick the graph filter **from principle**, then let performance *confirm* it — never pick
the AUC winner (that is cherry-picking the method itself). Derived from our own LRG framework
(Villegas: dense ⇒ degenerate single-peak; sparsify ⇒ multiscale) the filter must be **sparse,
cycle-rich (non-trivial, not a tree), parameter-free, connectivity-preserving, and
literature-established for correlation-like networks**. Exactly one method family clears all five:
the **PMFG/TMFG** planar-maximally-filtered graph (Tumminello 2005 / Massara 2016) — it is
parameter-free (edge budget fixed at 3(N−2) by planarity, so density ≈ 6/N auto-sparsifies with
N), provably contains the MST (a principled *enrichment* of the tree with the cycles the diffusion
needs), and was built for financial correlation matrices — the same object as our FC. Its only
honest weakness is the **planarity prior** (why genus-0 for a brain?); its one serious principled
rival is the **disparity filter** (Serrano 2009, "multiscale backbone" — no topological prior, but
carries an α parameter and can fragment). **Recommendation: adopt the PMFG/TMFG family**, cite
PMFG as the principled parent, use TMFG as its fast chordal computation, and **cross-check against
disparity** so the result is not a planarity artifact. Adoption is **conditional on the
pre-registered rule below** and settled by the full-suite run (marker + all arc analyses on the
TMFG/PMFG suite), not by this document.

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

## 4. Recommendation — the PMFG/TMFG family (choice A)

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

## 7. What the suite run needs (post-compaction plan)

Library first (`src/lrg_eegfc/utils/fc/backbone.py` — general names):
- **Have:** `tmfg_backbone`, `mst_union_top_fraction`, `percolation_backbone`, `maximum_spanning_tree`.
- **Add:** `pmfg_backbone(W)` (exact planar-maximal filter — edge insertion under a planarity test,
  e.g. networkx `check_planarity`; contains MST; 3(N−2) edges) and `disparity_backbone(W, alpha)`
  (Serrano per-edge disparity p-value). Both parameter-free to call except disparity's α (report a
  small α-sweep, not a tuned point).

Re-point each pipeline via the `SA_BACKBONE` env (branch already in `06_epi_arc.py`:
`elif BACKBONE=="tmfg": B = tmfg_backbone(W)` — replicate for `pmfg`, `disparity`):
- `06_epi_arc.py`            — marker (done for tmfg → `data/sparsified_arc/epi_arc_tmfg/`).
- `13_matched_strength_mst020.py` — the τ-resolved trace gate (the decisive one for §1).
- `05_enc_inf_arc.py`        — encoding/inference dissociation.
- `12_sparsification_recovery.py` — backbone sweep (add tmfg/pmfg/disparity as named points).
- `14_controls_ladder_mst020.py`, `15_patient_variation.py` — downstream confirmations.

Run all with `SA_BACKBONE ∈ {tmfg, pmfg, disparity}` vs the `mst020` baseline; env pins
`OMP/OPENBLAS/MKL=1`, python `/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python`.
Compare per-scale matched-strength clears (trace) and matched-strength beats-null + fused AUC
(marker) against the pre-registered rule. **Then, and only then, settle the method.**

## 8. Do-not (carry the honesty rules)

- Matched-strength is the only null; per-scale, never best-scale (Villegas-collapse δ/high-γ stay
  un-claimed). Don't re-attribute the Grassmann (dense/whole-task) to any sparsified backbone.
- Don't adopt on the marker win alone — the trace must be checked (its mst@0.20 was chosen because
  **α needs a strong-edge backbone**; TMFG keeps the strongest edges + triangles, so it *should*
  preserve α, but VERIFY — see the trace handoff).
- Don't invent a filter — PMFG/TMFG/disparity are named literature methods; anything else needs a
  citation before use.
