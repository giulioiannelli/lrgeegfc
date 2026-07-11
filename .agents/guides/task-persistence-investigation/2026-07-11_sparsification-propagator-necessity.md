---
name: sparsification-propagator-necessity
type: guide
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
status: scope-report-ready-for-implementation
created: 2026-07-11
updated: 2026-07-11
pointers:
  - scripts/01_compute/audit/audit_174_diffusion_vs_rawfc_gate.py
  - scripts/01_compute/audit/audit_150_rho_sym_gate.py
  - .agents/guides/02_methods/lrg-framework-guide.md
  - memory/audit_174_diffusion_vs_rawfc_gate_2026_07_11.md
  - memory/lrg_outlier_case_fully_connected.md
---

# Does the diffusion propagator become load-bearing on a sparse backbone? (scope → audit_175)

## Head

**audit_174 (2026-07-11) found that on the FULLY-CONNECTED FC graph at τ=1/λ_max, the
β cophenetic trace is IDENTICAL with or without diffusion** (raw-FC `D=1/A` UPGMA tree
clears the matched-strength ρ_sym gate at p=.032, same as the diffusion tree). α needs the
diffusion (CLEAR .024 → marginal-fail .053 without); low_γ inverts (raw CLEAR .019 → diffusion
suppresses .080). **This is a pay-attention flag on the LRG-diffusion *framing* (NOT the
empirical trace), and it is the quantification of a known risk: LRG is degenerate on
fully-connected graphs (guide §6; `lrg_outlier_case_fully_connected`).**

**Hypothesis (this scope):** LRG is designed for **sparse / topological** graphs. On a
principled sparse backbone of the FC graph, the propagator does something raw-FC
*structurally cannot* — it assigns a finite distance between contacts with **no direct edge**
(multi-step reachability). So on the backbone the diffusion should (a) diverge materially from
raw-FC, (b) keep or strengthen the band trace, and (c) let C(τ) develop the multi-scale ladder
the LRG framework promises. **If confirmed, this rescues and sharpens the diffusion framing; if
refuted, the framing must be reworked. DO NOT CASCADE audit_174 until this is run.**

---

## 5-point critical preamble

1. **Claim.** On a pre-committed sparse backbone (swept over a density grid), the diffusion
   cophenetic tree diverges materially from every no-diffusion baseline AND the β/α ρ_sym gate
   survives or strengthens — demonstrating the propagator is load-bearing in the graph regime
   LRG is built for.
2. **Null.** Sparsification does not change the picture: diffusion ≈ the no-diffusion baseline
   on the backbone too (β identical), OR sparsification destroys the trace (β fails at all
   reasonable densities). Either means sparsification does not rescue the diffusion framing.
3. **Strongest alternatives the design MUST control for.**
   - *p-hacking via the sparsity parameter* (picking the density where diffusion wins).
     CONTROL: a pre-committed **density grid**; report every quantity **as a function of
     density**; the effect must be robust/monotone across the grid, never a single magic point.
   - *The raw-FC baseline becomes ill-defined on sparse graphs* (`1/A=∞` for removed edges).
     This is itself informative (raw-FC breaks → a multi-step distance is *required*), but it
     is a weak "win". The SHARP test is diffusion vs a **well-defined multi-step baseline that
     is NOT diffusion**: the **graph geodesic (shortest-path) distance** on the backbone. If
     diffusion ≈ geodesic, the message is "you need topology, but shortest-path suffices —
     diffusion is one valid choice, not uniquely necessary"; if diffusion ≠ geodesic and
     separates trace-from-null better, the LRG all-paths integration is specifically valuable.
   - *Sparsification removes the edges that carry the trace.* CONTROL: verify β survives
     sparsification; report β obs ρ_sym vs density (it must not collapse).
   - *Shared-draw / estimator artifact.* CONTROL: identical surrogate draws feed every arm;
     the diffusion arm at density=1.0 must reproduce audit_174 bit-for-bit (anchor).
4. **Cannot.** A backbone that rescues diffusion does NOT retroactively validate the
   fully-connected pipeline (still degenerate) — it MOTIVATES migrating the pipeline to the
   backbone (a large re-validation of β→OFC, encoding/inference, epi). Filter choice
   (disparity vs proportional) is itself a modeling choice → robustness across filters required.
5. **Falsify.** If across the WHOLE density grid the diffusion tree stays ≈ geodesic/raw tree
   AND β clears identically in both, the propagator is inessential even on the backbone → the
   diffusion framing is not rescued. If β dies at all reasonable densities, the trace lives in
   the weak-edge bulk and the backbone approach fails.

---

## Notation

- `A` — per (patient, phase, band) `|ImCoh|` FC matrix, `N×N`, symmetric, non-negative,
  fully connected (all off-diagonal > 0). Load via `workflow.fc.load_fc_matrix` /
  `audit_150.load_phase` (which also serves the split-halves A, B of rest_pre).
- `F_ρ(A)` — a **sparsification filter** at density parameter ρ (see below) → sparse backbone
  `A_s` (a subset of edges kept, weights preserved), then restricted to the **giant component**
  `G_s` (LRG works on the giant component; drop isolated nodes, record how many).
- `d ≡` resulting **edge density** = kept edges / `N(N−1)/2`. The grid is reported in `d`.
- Diffusion arm (LRG): `L_s = D_s − A_s`, `τ = 1/λ_max(L_s)`, `K_s = e^{−τ L_s}`,
  `D^{diff}_{ij} = (1−δ_{ij})/K_{s,ij}` → UPGMA → cophenetic `C^{diff}`.
- Geodesic arm (no diffusion, still topological): edge length `ℓ_{ij} = 1/A_{s,ij}` (or
  `−log(A_{s,ij}/A_max)` as a robustness variant); `D^{geo}_{ij} =` shortest-path length
  (`scipy.sparse.csgraph.shortest_path`) → UPGMA → cophenetic `C^{geo}`.
- Plain-raw arm (control that BREAKS): `D^{raw}_{ij} = 1/A_{s,ij}`, `∞` for non-edges → record
  that UPGMA is undefined / degenerate (this is a *result*, not a bug).
- ρ_sym estimator, matched-strength surrogate: reuse `audit_174` verbatim (which reuses
  `audit_150`), only swapping the tree-construction function.

## Sparsification filters (pre-committed)

- **PRIMARY — disparity filter** (Serrano, Boguñá, Vespignani 2009): reuse
  `lrg_eegfc.visuals.network_templates.disparity_backbone(A, alpha)` (already in the repo; also
  `top_k_backbone(A, keep_frac)`; and `utils/fc/msc/sparsify.py:disparity_filter`). Sweep the
  significance `alpha_ds` over a grid so the resulting density `d` spans e.g. **{0.03, 0.05,
  0.10, 0.20, 0.35, 0.50, 1.0}** (1.0 = no sparsification = audit_174 anchor).
- **SECONDARY — proportional threshold** (robustness): `top_k_backbone(A, keep_frac)` at the
  same target densities. Report both filters; the verdict must not depend on the filter.
- Always symmetrize + re-restrict to the giant component after filtering.

## Predicates (what the run must decide)

- `P_divergence(d)` — cohort-mean `Spearman(C^{diff}, C^{geo})` at density `d`. Expectation:
  falls below the fully-connected 0.83 as `d` shrinks (diffusion departs from shortest-path).
- `P_gate_diff(d)`, `P_gate_geo(d)` — per-band cohort ρ_sym matched-strength gate p-values
  (Wilcoxon one-sided, as audit_150) for each arm, vs density.
- `P_beta_alive(d)` — β obs ρ_sym (cohort median) vs density; must not collapse.
- `P_ladder(d)` — number of interior peaks of `C(τ) = −dS/dlog τ` on the backbone (reuse the
  LRG `entropy`/`entropy_C` machinery). Fully-connected has ONE peak (collapse scale, guide §6);
  a **multi-peak ladder emerging as `d` shrinks** is the strong positive signal that the
  backbone restores genuine LRG multiscale structure.

**Verdict logic.** Diffusion is *rescued* on the backbone iff, robustly across the grid and both
filters: `P_divergence` drops well below fully-connected, `P_gate_diff` keeps β/α CLEAR while
`P_gate_geo` (or raw) does NOT reproduce it, and ideally `P_ladder` > 1 emerges. Diffusion is
*not rescued* iff diffusion ≈ geodesic across the whole grid, or β dies.

## Properties / anchors

- **Anchor:** at `d = 1.0` the diffusion arm MUST reproduce audit_174 (= audit_150): α p=.024,
  β p=.032 CLEAR; δ/θ/low_γ/high_γ fail; and `Spearman(C^{diff}, C^{raw-1/A})` ≈ the 0.83–0.93
  from audit_174. If not, the harness is wrong.
- The null must be applied the SAME way in every arm: **sparsify the matched-strength surrogate
  with the same filter at the same `d`**, then diffuse/geodesic → ρ_sym. (Keeps the null
  pipeline identical to the observed pipeline; the strength-matched surrogate has the same
  strengths so its backbone density is comparable — record it.) Flag the alternative
  (sparse-native degree+strength-preserving null) as an open robustness item.

## Caveats

- Object change: a rescued backbone means MIGRATING the pipeline → all downstream results
  (β→OFC, enc/inf, epi) re-validate on the backbone. Big program; do not promise a quick patch.
- Disconnection: aggressive filters fragment the graph; restrict to the giant component and
  report node loss per `d`. If the giant component is tiny at small `d`, that density is invalid.
- Geodesic on `1/A` weights over-weights weak long paths; also try `−log A` lengths.
- β may live in the weak-edge bulk → could die under sparsification; that is a real (negative)
  outcome, report it honestly, do not tune `d` to keep β.

## Pseudocode

```
GRID = [1.0, 0.50, 0.35, 0.20, 0.10, 0.05, 0.03]          # target densities
for filt in [disparity_backbone, top_k_backbone]:
  for d in GRID:
    for (pat, band) in COHORT × BANDS:
        Ws = {A,B,task_test,rest_post -> load_phase}
        Ws_s = {ph: giant_component(filt(W, d)) for ph in Ws}   # same node set across phases!
        # observed
        Cdiff = {ph: cophenet(UPGMA(1/K(tau_min, Ws_s[ph])))}
        Cgeo  = {ph: cophenet(UPGMA(shortest_path(1/Ws_s[ph])))}
        obs_diff = rho_sym(Cdiff);  obs_geo = rho_sym(Cgeo)
        # matched-strength null: sparsify the SAME surrogate the SAME way
        for r in R:
            Wsurr = {ph: giant_component(filt(shuffle(Ws[ph], rng), d))}   # same rng/order as audit_150
            surr_diff[r] = rho_sym({cophenet(UPGMA(1/K(Wsurr[ph])))})
            surr_geo[r]  = rho_sym({cophenet(UPGMA(shortest_path(1/Wsurr[ph])))})
        record obs/surr for both arms, tree-divergence Spearman(Cdiff,Cgeo), C(tau) peak count
    cohort_gate per band (Wilcoxon obs-surr_p50, one-sided) for diff and geo, at this (filt,d)
```
CAUTION on the null + giant component: filtering a surrogate can change `N` (node loss) per
draw → the four phases may end with different node sets. Enforce a **common kept-node set per
cell** (intersect the four phases' surviving nodes, or filter the observed rest_post backbone
and reuse that edge/node mask for all phases + surrogates) so ρ_sym's per-pair alignment holds.
This is the single most error-prone part — get it right and assert equal `n_pairs` across phases.

## Visualization

- `gate_p(diff) vs gate_p(geo)` per band, as a function of `d` (one line per band per arm).
- tree-divergence `Spearman(C^{diff}, C^{geo})` vs `d` (should fall as `d` shrinks).
- `C(τ)` peak-count vs `d` (ladder emergence).
- β obs ρ_sym vs `d` (survival).

## Connection to prior tools

- Direct sequel to `audit_174` (fully-connected ablation) — reuse its harness + `audit_150`
  ρ_sym gate + numba surrogate.
- `audit_172/173` swept τ on the fully-connected graph (β scale-broad, α mesoscale, coarse=
  collapse). This scope sweeps DENSITY instead of τ — orthogonal axis.
- `audit_166` (pairwise-descriptor-ladder): cophenetic hierarchy beats raw SCALARS — the
  hierarchy is load-bearing; this scope asks whether DIFFUSION is needed to build it.
- guide §6 + `lrg_outlier_case_fully_connected`: predicts diffusion is degenerate fully-connected
  and should revive on a sparse/topological graph — this is the direct test of that prediction.

## Open questions

1. Is there a *principled, data-independent* density (e.g. percolation threshold of the giant
   component, or the disparity `alpha` that maximizes backbone significance) to pre-commit to,
   rather than a grid? Pre-registering one avoids the p-hacking critique entirely.
2. Effective-resistance distance (Laplacian pseudoinverse) as a third multi-step baseline
   between geodesic and diffusion.
3. If diffusion is rescued: the migration plan (which downstream results re-run, in what order).
