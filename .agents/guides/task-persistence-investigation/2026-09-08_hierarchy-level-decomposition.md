---
name: hierarchy-level-decomposition
type: scope
era: IMCOH_ABS × COHORT_N10 → PAPER_FINALIZATION (Wave 2)
status: draft
created: 2026-09-08
updated: 2026-09-08
pointers:
  - .agents/plans/active/2026-09-08_new-paths-to-objectives.md
  - .agents/preprint/locked/PIPELINE_CONTRACT.md
  - .agents/reports/2026-08-31_lane-s-scale-variability.md
  - scripts/01_compute/paper_final/feasibility/feas_level_decomp.py
  - data/paper_final/feasibility/level_decomposition.csv
---

# Hierarchy-level decomposition of the cross-phase trace

## Head

Bin node pairs by the level of the rest_pre diffusion tree at which they merge, and compute the cross-phase trace statistic separately inside each level. Levels are disjoint pair sets, so per-level statistics are independent by construction, and "which scale carries the trace" becomes a within-patient profile over levels instead of a family of near-identical whole-vector correlations at different τ. Same-day feasibility (no null): the β trace is absent at the top 2–4 modules and present at fine levels, with a within-patient trend in 9/10 patients.

## Five-point preamble

1. **Claim.** The retained task change (rest_pre → task → rest_post) is concentrated at specific levels of the diffusion hierarchy, in a band-dependent way, consistently across patients.
2. **Null.** Three rungs. (a) Ordered sham per level: the same five-phase arc carved from a single resting recording in true temporal order, binned by its own rest-half tree, gives the level profile that temporal drift alone produces. (b) Matched-strength per level: strength-preserving rewiring of each phase's adjacency, re-binned by the rewired reference tree, gives the profile that degree structure alone produces. (c) Backbone-fixed N1 as in the contract. The scale-specificity statistic is the within-patient level trend (or fine-minus-coarse contrast) of the real arc minus that of the sham; the absolute per-level value is never the claim.
3. **Strongest alternative.** Fine-level pairs are the strongest, most reliable pairs; any positive-biased correlation statistic will be larger on them. A level gradient could therefore be a signal-to-noise gradient with no task content.
4. **Does the null control for it.** The sham does, by mechanism: it has the same pair-strength and reliability structure (same recording, same tree construction) and no task, so any gradient it shows is the reliability gradient. Matched-strength controls the degree part but not the reliability part, so it is the weaker rung here. What neither can reject: a task effect that is itself a pure gain change of the strongest pairs (that would still be a task effect, but not a hierarchical one); the module-level localisation (§9) distinguishes the two.
5. **Falsifier and limitations.** The claim is false if the sham reproduces the level gradient, or if the gradient disappears under knob integration over f ∈ [0.07, 0.20] and reference-tree scale, or if it fails LOO. Limitations: seven levels only; the reference tree is built on the full rest_pre while the statistic uses its halves (mild dependence, §6); level membership is patient-specific, so cross-patient aggregation is on k, not on anatomy.

## 1. Notation

- N nodes; P = {(i, j) : i < j} the set of node pairs, |P| = N(N − 1)/2. Condensed index p ∈ P.
- Phases X ∈ {A, B, learn, test, post}, with A and B the contiguous halves of rest_pre. W_X ∈ [0, 1]^{N×N} the dense |ImCoh| adjacency; B_X = backbone(W_X; f) the mst-union backbone at fraction f; L_X = D_X − B_X.
- ρ̂_X(s) = e^{−τ L_X} / Tr, with τ = s / λ_max(L_X). The reference tree T_ref = UPGMA(1 / ρ̂_pre(s_ref)) built on B_pre, the backbone of the **full** rest_pre.
- For pair p, h_p = cophenetic height of p in T_ref; k_p = number of clusters remaining immediately after the merge that joins p, k_p = N − #{merges of T_ref with height ≤ h_p}. Range 1 ≤ k_p ≤ N − 1.
- Levels ℓ = 1..7 with edges in k: {2}, (2, 4], (4, 8], (8, 16], (16, 32], (32, 64], (64, ∞). P_ℓ = {p : k_p ∈ level ℓ}. The P_ℓ partition P.
- Same-shaft mask S ⊂ P; the working pair set is P_ℓ \ S.
- Readout vectors x_X ∈ ℝ^{|P|}: either raw, x_X = W_X[p], or heat, x_X = ρ̂_X(s_read)[p]. Never the cophenetic vector (§5).
- Statistic inside level ℓ: ρ_sym,ℓ = ½[ρ_S(x_test − x_A, x_post − x_B) + ρ_S(x_test − x_B, x_post − x_A)] with ρ_S the Spearman correlation over p ∈ P_ℓ \ S. Same for T_learn,ℓ with x_learn in place of x_test.

## 2. Definitions

- Level profile of a patient and band: r = (ρ_sym,1, …, ρ_sym,7).
- Split-half reliability per level: rel_ℓ = ρ_S(x_A, x_B) over P_ℓ \ S.
- Within-patient trend: t = ρ_S(ℓ, ρ_sym,ℓ) over the levels with |P_ℓ \ S| ≥ 30.
- Coarse-vs-fine contrast: c = mean_{ℓ ≥ 5} ρ_sym,ℓ − mean_{ℓ ≤ 2} ρ_sym,ℓ.
- Cohort test: Wilcoxon signed-rank of t (or c) across patients, one-sided, on the real arc minus the patient's own sham; LOO on the same; BH only across the coordinated family of six bands.
- Band × level interaction: Friedman across bands of c, then post-hoc pairs.

## 3. Properties

- Range of ρ_sym,ℓ: [−1, 1]; positive bias of unknown size inherited from ρ_sym, which is why only profile shape (t, c) and sham-referenced margins are admissible.
- Independence: P_ℓ are disjoint, so the seven statistics share no pairs; their only shared content is through the phase matrices' global structure (e.g. a common gain), which the sham also carries.
- Invariance: t and c are invariant to a level-independent additive bias of ρ_sym; they are not invariant to a reliability gradient, which is the job of the sham.
- Comparability across patients: k is dimensionless and comparable across implants; τ, s and anatomy are not. The construction aggregates on k by design.
- Complexity: one eigendecomposition per phase per backbone fraction, one linkage per reference scale; seconds per cell. Sham cost equals Lane E's sham.
- What it cannot detect: a task effect confined to pairs that change level between phases (level membership is fixed by T_ref); an effect expressed as a re-ordering of the hierarchy rather than as a change of coupling within a fixed hierarchy. Those are the questions Grassmann and tree-distance measures address; this measure is orthogonal to them.

## 4. Caveats and failure modes

- **Reliability gradient** (strongest alternative). Mitigation: sham per level; reporting ρ_sym,ℓ / rel_ℓ is descriptive only.
- **Reference-tree dependence.** T_ref uses the full rest_pre while the statistic uses A and B. Mitigation: repeat with T_ref built on A only (then the statistic is asymmetric in A/B and the role-swap term must be built on B's tree); report both.
- **Level population.** Level 1 (k = 2) and level 7 (k > 64) are the smallest (≈ 170–730 pairs) and level 7 has the highest same-shaft share (23 %). Mitigation: same-shaft mask always on; minimum 30 pairs per level; sensitivity with edges shifted by half an octave.
- **Ties in merge heights.** UPGMA on 1/ρ̂ rarely ties; if it does, k_p uses the last tied merge (side = right in the search). Deterministic.
- **Backbone fraction.** Levels move with f. Mitigation: knob integration over f ∈ {0.07, 0.10, 0.14, 0.20} as in the contract; today's probe shows the β gradient in the same direction at every f.
- **Reference scale s_ref.** The tree changes little between s_ref = 1 and 4.7 but the probe shows weaker p at 4.7 for the heat readout. Report both; do not pick the better.
- **Ultrametric readout.** The cophenetic vector is constant within a merge node and must not be used inside levels; only the raw adjacency or the heat kernel at fixed s.

## 5. Pseudocode

```
for patient in cohort:
  for band in bands:
    load W_X for X in {A, B, learn, test, post}; load W_pre (full rest_pre)
    S := same-shaft pair mask
    for f in fracs:
      B_pre := backbone(W_pre, f); (ev, V) := eig(L(B_pre))
      for s_ref in {1, 4.72}:
        Z := UPGMA(1 / rho_hat(ev, V, s_ref)); h := cophenet(Z)
        k_p := N - count(Z.heights <= h_p) for each p; level_p := bin(k_p)
        for readout in {raw, heat(s_read = 1)}:
          x_X := readout(W_X or backbone(W_X, f)) for each X
          for level in 1..7:
            P_l := {p : level_p == level and not S_p}
            if |P_l| < 30: rho_l := NaN; continue
            rho_l := rho_sym(x_A, x_B, x_test, x_post restricted to P_l)
            rel_l := spearman(x_A, x_B restricted to P_l)
          t := spearman(level, rho_l); c := mean(rho_{5..7}) - mean(rho_{1..2})
          store (patient, band, f, s_ref, readout, rho_l, rel_l, t, c)
    repeat the whole block on each ordered-sham arc of the patient -> t_sham, c_sham
cohort: for each band, readout: Wilcoxon over patients of (t - t_sham), one-sided; LOO; BH over bands
```

## 6. Visualization

One panel per band: per-patient level profiles as thin lines over k (log₂ axis, "modules remaining"), the cohort median as a thick line, the sham median as a dashed line, and the per-level sham spread as a band. No threshold lines; every significant level marked uniformly. Second figure: the fine-level modules of one representative patient per band drawn on the implant with their Desikan-Killiany labels (§9).

## 7. Connection to prior tools

- Whole-vector ρ_sym over s (Lane S): the special case where every level is pooled; its n_eff ≈ 1 is the reason this decomposition exists.
- Per-node trace decomposition (2026-06-25): node-level projection of the same signal; the level decomposition is the pair-level analogue with the hierarchy as coordinate.
- Module-retention landscape (April, archived): asked whether modules exist in task and persist in post; this asks whether coupling within fixed hierarchy levels persists. It uses no partition-agreement scalar (no ARI/NMI/VI), so the April prohibition on partition metrics does not apply.
- Same-shaft `pair_mask` (W0-A): reused unchanged.

## 8. Open questions

- Should levels be defined per patient on k (current) or on N_eff(s) of the reference spectrum, which would tie the level coordinate to the LRG's own resolution measure?
- Is the β negative value at k ∈ (2, 4] (5/10 patients) a reset of top-level coupling or noise? Only the sham says.
- Does the fine-level trace exclude SOZ nodes (the physiology–pathology link in the plan §5)?

## 9. Module-level localisation (for objective 2)

For each patient and band, the levels that carry the trace define a set of modules of T_ref (the clusters present when k modules remain). Each module has a Desikan-Killiany footprint; pool across patients by region with the resolution-preserving cohort statistic (`cohort_localization.py`), whole-grid BH and LOO. This is a different object from the node-score localisation that returned "delocalised" in July and has not been tested.
