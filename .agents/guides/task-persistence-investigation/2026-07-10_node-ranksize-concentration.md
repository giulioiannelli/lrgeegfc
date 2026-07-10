---
name: node-ranksize-concentration
era: IMCOH_ABS_COHORT_N10
status: result
kind: analysis
scope: Rank-size (Zipf) decomposition of the per-node trace contribution T_i. Tests whether the cohort trace ρ_sym is carried by a few driver nodes (heavy-tailed) or delocalised. Verdict — REFUTED for the real signal — delocalisation IS the trace; a heavy tail is trace absence.
updated: 2026-07-10
---

# Rank-size (Zipf) view of the per-node trace contribution

**Head.** Decompose the cohort trace `ρ_sym` into per-node "sizes" and rank-size
plot them. The user's hypothesis (Zipf spirit): most nodes inert, a few drivers →
a power law. **Data verdict: the inversion.** A *strong* trace is **delocalised**
(β strong-trace patients: ~88% of nodes co-move, top-10 carry <20%, Gini|T|≈0.32);
the **Zipf-like concentrated** profile (a few idiosyncratic nodes dominate,
Gini|T|≈0.51) is the signature of trace **absence** (resetters, θ/δ null bands).
Cohort `Spearman(ρ_sym, Gini|T|) = −0.50` (sign-independent). A heavy tail is what
*noise* looks like, not the mechanism of the trace.

## Notation / the "size"
Reuses the per-node decomposition of audit_154 (no new measure):
- Per-pair rank concordance `c_p = (rank(ΔD_task)_p − mid)(rank(ΔD_rest)_p − mid)`,
  symmetrized over both split-half arm assignments (`c_sym`); `mean_p c_p ∝ ρ_sym`.
- **Per-node size** `T_i = node_incidence_mean(c_sym)` = mean of `c_p` over pairs
  incident to node `i`. **Exactly `mean_i T_i ∝ ρ_sym`** ⇒ the trace is the *sum of
  node sizes*, so "a few nodes drive it" is a well-posed rank-size question.
- Cached for all 6 bands, n=10: `data/audit/per_node_trace_decomposition_rhosym/per_node.csv`.

## Predicates / concentration metrics
- **frac co-moving** = `mean(T_i > 0)` — fraction of nodes moving *with* the trace.
  (Interpretable but partly definitional: a positive mean ⇒ more positive nodes.)
- **Gini(|T_i|)** — sign-INDEPENDENT inequality of contribution magnitudes. Chosen
  as the headline because the naive `Gini(positive part)` is slaved to `sign(mean
  T)=sign(ρ)` and gives a tautological −0.96; `Gini(|T|)` gives an honest −0.50.
- **PR/N** = `(Σ|T_i|)² / (N·Σ T_i²)` — effective participating fraction (1 = uniform
  / delocalised; →1/N = one node).

## Properties (results)
| band | ρ_sym (med) | Gini\|T\| | PR/N | frac co-moving |
|---|---|---|---|---|
| **β** | 0.198 | 0.433 | 0.614 | 0.684 |
| **α** | 0.101 | 0.484 | 0.525 | 0.620 |
| low_γ | 0.089 | 0.412 | 0.604 | 0.562 |
| high_γ | 0.072 | 0.440 | 0.571 | 0.554 |
| δ | 0.032 | 0.519 | 0.423 | 0.518 |
| θ (null) | −0.038 | 0.531 | 0.459 | 0.368 |

β per-patient: **strong-trace** (ρ>0.3, Pat_08/05/03/02) Gini\|T\| 0.32, PR/N 0.75,
frac 0.88; **reset** (ρ<0, Pat_15/10/14) Gini\|T\| 0.51, PR/N 0.44, frac 0.25.
Cohort `Spearman(ρ_sym, Gini|T|) = −0.50` (p=4e-5). The **log-binned** size density
(log-spaced bins, counts/width/N, pooled per group) is **peaked around the typical
size** — a characteristic scale, not a straight line (OLS slopes −0.43 strong / −0.61
reset, far from a clean power-law −1..−3) ⇒ **not a power law**; the reset pool is
only mildly heavier-tailed.

## Caveats
- **frac co-moving is partly definitional** (tracks `sign(mean T)`); the load-bearing,
  non-tautological claim is on **Gini(|T|) / PR/N** (magnitude spread, sign-free).
- Descriptive/exploratory: no surrogate null on the concentration itself (the
  underlying `ρ_sym` and per-node polarity are MS-tested in audit_154). The rank-size
  is a *shape* statement about an already-verified quantity, not a new hypothesis test.
- Reset patients' low PR/N is partly mechanical (negative mean ⇒ few positive nodes),
  but the *magnitude* concentration (Gini|T| 0.51 vs 0.32) is a genuine structural
  difference (coherent vs idiosyncratic reorganization), not a sign artifact.

## Connection to prior tools
- Quantifies [[dendrogram_persistence_gate_2026_06_04]] "broad rest_post = delocalised
  trace" and [[figure_trace_heat_brain_2026_07_09]] "β delocalised single-contact".
- Reconciles the delocalised-vs-OFC tension: the β trace is delocalised at the *single
  node* level (this analysis) yet localizes to OFC at the *pair / system + gated* level
  ([[localization_audit_plan_2026_05_29]], [[figure_3d_pairglow_ofc_node_not_pair]]).
- Per-node size = audit_154 `T_i` (scope `2026-06-25_per-node-trace-decomposition.md`).

## Open questions
- Is the delocalisation itself a matched-strength-survivable property (would need a
  per-node-concentration MS null), or a generic property of any positive rank
  correlation? Deferred — descriptive for now.

Build `scripts/01_compute/figures_embedded/fig_trace_ranksize_nodes.py`;
figure `data/outputs/figures/trace_ranksize_nodes.pdf`.
