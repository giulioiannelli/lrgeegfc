---
name: occult-epi-node-marker
type: scope
era: COHORT_N10
status: built-negative
created: 2026-06-05
updated: 2026-06-05
pointers:
  - .agents/guides/task-persistence-investigation/2026-06-05_epi-node-trace-marker.md
  - scripts/01_compute/audit/audit_79_epi_node_behavior.py
  - scripts/01_compute/audit/audit_80_epi_marker_classifier.py
  - scripts/01_compute/audit/audit_90_occult_node_lno_influence.py
  - scripts/01_compute/audit/audit_77_epi_stratified_cophenetic.py
  - src/lrg_eegfc/utils/surrogate/matched_strength.py
  - src/lrg_eegfc/utils/io/patient.py
---

# Occult (non-labeled) epileptic-node marker from cross-phase LRG trace behaviour

**Head.** This is the *discovery* reframe of the per-node marker question:
not "classify the KNOWN epileptic contacts" (that scope is
[`2026-06-05_epi-node-trace-marker.md`](2026-06-05_epi-node-trace-marker.md),
verdict NEGATIVE — `geometry_dominant` in every band) but "do any
**non-labeled** contacts BEHAVE structurally like the labeled epileptic ones,
in a way not explained by node strength or sEEG-shaft geometry — i.e. candidate
*occult* epileptogenic tissue?". The honest design residualizes every per-node
structural-behaviour score against node strength **and** same-shaft/probe
geometry **first** (the two confounds that sank `audit_79`/`audit_80`), then
asks whether residual behaviour flags non-labeled nodes inside the labeled-epi
range. **Without surgical-outcome ground truth, any "occult candidate" is a
hypothesis, not a validated marker** — and the a-priori-likely outcome, given
the negative known-node classifier, is "no occult signal beyond strength+shaft".

> **RESULT (2026-06-05 — `audit_90` BUILT; verdict: no transferable occult
> marker, with one genuine within-patient nuance).** The new leave-node-out
> (LNO) cohort-trace influence `Δρ_i` and the residualized per-node behaviour
> scores were computed for the three trace-bearing bands (α, β, γ_l). Two
> findings, in order of decisiveness:
>
> 1. **NEW (positive, within-patient):** at **β**, `Δρ_i` IS a real,
> strength-INDEPENDENT epileptic signal — per-patient Spearman(Δρ, strength) ≈
> −0.06 (so not a strength re-read), epi-vs-non-epi cohort AUC **0.634, 8/9
> patients, Wilcoxon p = 0.010**, and it **survives the rank-partial
> nonlinear-strength control** (AUC 0.650, p = 0.0076). Sign: labeled-epi nodes
> have positive median `Δρ` (trace-CARRIERS — removing them weakens the trace),
> non-epi negative (trace-suppressors). This is genuinely new relative to
> `audit_79`/`audit_80`: the row-restricted `rho_split_node` (F1) was at chance
> (β AUC 0.50), because LNO *leverage on the cohort trace* is a different
> statistic from *the trace measured on the node's row*. **Caveat:** the
> within-patient β AUC (0.634) is comparable to — not above — the strength
> baseline AUC (0.652), so even within-patient it is not a dominant axis.
>
> 2. **DECISIVE (negative, cross-patient):** in **leave-one-patient-out**
> transfer — the actual occult/discovery use-case, where an unlabeled patient
> has NO within-patient labels to set a threshold — `Δρ_i` is **at chance**
> (β PR-AUC 0.123 vs prevalence 0.106, **lift ≈1.2×**, ROC 0.562; α/γ_l the
> same). This is the SAME ceiling the `audit_80` classifier hit. The
> within-patient β separation does **not transfer** to held-out patients.
>
> **Verdict: no transferable occult marker** — consistent with `audit_80`
> (`geometry_dominant`). The 154 candidates emitted
> (`occult_candidates.csv`, within-patient top-decile residual-influence
> non-labeled contacts, 122 on epi-free shafts) are **hypothesis-generating
> ONLY**; with no surgical-outcome ground truth they cannot be validated, and
> the cross-patient test says a non-labeled patient could not be flagged from
> this score. The direction is closed without ground truth. Artifacts:
> `data/audit/occult_node_marker/{lno_node_influence,occult_gate_per_band,occult_lopo_transfer,occult_candidates}.csv`
> + `README.md`.

---

## 1. Notation

Inherits the notation of the companion scope
([`2026-06-05_epi-node-trace-marker.md`](2026-06-05_epi-node-trace-marker.md) §1)
and adds the leave-node-out operator.

- Patient `p`, `N_p` contacts, FC-order node index `i ∈ {0,…,N_p−1}`, band
  `b ∈ {δ,θ,α,β,γ_l,γ_h}` (each band an independent experiment).
- Phases `Φ = {rest_pre_A, rest_pre_B, task_test, rest_post}` (split-baseline).
- `W^{φ}_p(b) ∈ [0,1]^{N×N}` symmetric zero-diag FC (`imcoh_abs`).
- `D^{φ} = lrg_ultrametric_condensed(W^{φ})` — condensed LRG cophenetic
  distance at τ=1/λ_max, average linkage (audit_63 primitive).
- `dT = D^{task} − D^{preA}`, `dR = D^{post} − D^{preB}` (locked sign).
- `ρ^{coph}_p(b) = Spearman(dT, dR)` — the patient-level cophenetic ρ_split
  trace (positive = TRACE), the SAME statistic the cohort headline uses.
- **Leave-node-out (LNO) operator:** `W^{φ}_{∖i} = W^{φ}[K_i, K_i]` with
  `K_i = {0,…,N−1}∖{i}` (drop row+col `i`). `ρ^{coph}_{p,∖i}(b)` is `ρ^{coph}`
  recomputed on the four LNO submatrices. (Each LNO rebuild is a fresh
  `N−1`-node LRG: degree recomputed, Laplacian, eigendecomposition, τ=1/λ_max,
  average linkage, cophenetic.)
- `E_p ⊆ {0,…,N_p−1}` labeled epileptic set (`build_epi_masks(p).epi_mask`),
  `y_i = 1[i ∈ E_p]`. `U_p = {i : y_i = 0}` = NON-labeled (candidate-occult
  pool).
- `s_i = mean_φ Σ_j W^{φ}_{ij}` — node strength (4-phase mean), THE dominant
  confound.
- `probe(i)` shaft id, `contact(i)` along-shaft index (`parse_seeg_label`).
  `nse_i = Σ_{j≠i} 1[probe(j)=probe(i) ∧ y_j=1]` (same-shaft labeled-epi
  count); `pos_i = contact(i)/max_contact(probe(i))` (normalized along-shaft
  position). These are the **geometry/leakage** confounds.

---

## 2. Definitions

### 2.1 LNO cohort-trace influence (the new, deferred-F8-cophenetic score)

For node `i`:

```
Δρ_i  =  ρ^{coph}_p(b)  −  ρ^{coph}_{p,∖i}(b)            ∈ [−2, 2]
```

`Δρ_i > 0`  ⇒ removing node `i` *weakens* the patient's trace ⇒ node `i`'s
incident pairs *carry* the trace (a trace-supporting node).
`Δρ_i < 0`  ⇒ removing node `i` *strengthens* the trace ⇒ node `i` is a
*trace-suppressor* (its row contributes noise/anti-trace).

This is the causal-flavoured complement to the row-restriction feature
`rho_split_node` (F1 in the companion scope): F1 asks "is the trace present on
node `i`'s row?"; Δρ_i asks "does the cohort-level trace *depend* on node `i`?".
It is the cophenetic analogue of the eigen-subspace LOO-resect F8 that both
`audit_79` and `audit_80` explicitly **deferred** — and the one genuinely new
per-node statistic relative to those audits.

### 2.2 Residualization against the known confounds (done FIRST)

The hard lesson of `audit_79`/`audit_80`: any per-node epi signal is, by
default, a re-read of strength or shaft geometry. So BEFORE any epi-vs-non-epi
comparison, regress every candidate score on the confound design and keep the
residual. Per (patient, band), OLS within the patient's nodes:

```
score_i  =  γ0 + γ1·s_i + γ2·nse_i + γ3·pos_i + γ4·shaftlen_i + e_i
score⊥_i =  e_i           (the residual — confound-orthogonal behaviour)
```

with `score ∈ {Δρ_i , rho_split_node_i , coherence_node_i}`. The residual
`score⊥` is the part of the node's structural behaviour NOT linearly
attributable to strength or shaft geometry. (Linear residualization is the
*weakest* defensible control; §4 lists what it cannot remove.)

### 2.3 Occult-candidate predicate

A non-labeled node `i ∈ U_p` is an **occult candidate** in band `b` iff its
confound-residual score lies inside the labeled-epi residual range AND above
the non-labeled bulk:

```
occult_i(b)  =  1[ score⊥_i ≥ Q_{0.90}( {score⊥_j : j ∈ U_p} ) ]
                · 1[ score⊥_i ∈ [ min, max ]( {score⊥_j : j ∈ E_p} ) ]
```

i.e. top-decile among non-labeled AND within the labeled-epi spread. This is a
**similarity** statement ("behaves like labeled epi, beyond strength/shaft"),
deliberately NOT a probability or a validated prediction.

### 2.4 Cohort discrimination gate (does the score even separate KNOWN epi?)

Before trusting any occult flag, the score must FIRST separate the *labeled*
epi from non-epi after residualization — otherwise "behaves like epi" is
meaningless. Per patient: `AUC_p = P(score⊥ on E_p > score⊥ on U_p)`
(Mann-Whitney). Cohort: `wilcoxon_z(AUC − 0.5, alternative='greater')` and the
count `n>0.5`. **Gate:** the residualized trace score is a usable occult axis
only if its cohort `|median AUC − 0.5|` exceeds the residualized
strength-baseline AUC. (Mirrors the `audit_79` gate, but on the residual and
including the LNO score.)

---

## 3. Properties

- **Δρ_i range** `[−2, 2]` (difference of two Spearman ρ); in practice
  `|Δρ_i| ≪ 1` since one node out of ~115 perturbs `N−1 ≈ 6.5k` pairs only at
  the margin. The score is a *leverage*, not an effect size.
- **Sum/centering:** `Δρ_i` does not sum to a fixed total (unlike the
  orthonormal mode-load F7); each LNO is an independent rebuild. No
  share-conservation constraint.
- **Residual orthogonality:** by construction `score⊥ ⟂ {s, nse, pos,
  shaftlen}` *linearly*. Any monotone-nonlinear strength dependence survives
  (mitigation: §4 rank fallback).
- **Identifiability of "occult":** the predicate is identifiable only relative
  to the labeled set `E_p`; with `E_p = ∅` (Pat_15) it is undefined — Pat_15
  contributes non-labeled scores but no occult test.
- **What it CANNOT detect:** (a) an occult focus whose cross-phase trace
  behaviour is identical to surrounding healthy tissue (no contrast — the most
  likely real situation); (b) epileptogenicity with no LRG cross-phase
  signature at all (static pathology); (c) any focus in a band with no cohort
  trace (δ/θ/γ_h cophenetic — no trace to be influential on); (d) a focus whose
  only signal is strength or shaft-position (residualized away by design —
  this is intentional, but it means a *genuinely strength-mediated* focus is
  invisible here).
- **Complexity:** `O(N)` LRG rebuilds per (patient, band, phase), each
  `O((N−1)³)` eigh — measured ~1 ms each at N≈115; full α/β/γ_l job ≈ 15 s.
  Cheap; no surrogate ensemble needed for the LNO score (the residualization
  *is* the strength control, far cheaper than per-node matched-strength).

---

## 4. Caveats & failure modes

| caveat | mechanism | mitigation |
|---|---|---|
| **No ground truth** | the cohort has clinical labels, not surgical-outcome maps; an "occult" flag cannot be confirmed or refuted | report flags as hypotheses ONLY; never call them predictions; state the validation gap in the first paragraph |
| **Linear residualization is weak** | strength→behaviour may be nonlinear; OLS leaves curvature in the residual | report a Spearman-rank-partial fallback; if the rank-partial AUC agrees with OLS, the linear control suffices |
| **Label leakage in `nse_i`** | same-shaft-epi count uses labels of *other* nodes | acceptable here (within-patient descriptive, not LOPO prediction); the occult flag is about the residual AFTER removing this geometry, so leakage makes the test *more* conservative (more variance attributed to geometry) |
| **epi spread defines the window** | tiny `\|E_p\|` (Pat_03 has few epi) → wide/unstable `[min,max]` window | require `\|E_p\| ≥ 3` for the occult-window test; report `\|E_p\|` per patient |
| **Pat_15 zero epi** | no labeled set → no occult window | excluded from the occult test; its nodes still scored for completeness |
| **multiple comparisons** | every non-labeled node is a candidate | no per-node p-value claimed; the deliverable is the cohort gate (does the axis separate KNOWN epi at all) + a ranked descriptive list, not a hypothesis test per node |
| **adjacency-to-known-epi confound** | a non-labeled contact one position away from a labeled focus is field/volume continuation, not a novel focus | report each candidate's along-shaft + 3D distance to nearest labeled epi; flag `on_epi_free_shaft` |
| **band with no trace** | δ/θ/γ_h cophenetic null → Δρ_i is leverage on noise | restrict to the trace-bearing bands α/β/γ_l (audit_77) |

---

## 5. Pseudocode

```
# audit_90 — occult-node LNO influence + residualized behaviour (CHEAP, no surrogate)
BANDS_TRACE = {alpha, beta, low_gamma}            # only bands with a cohort trace
for band b in BANDS_TRACE:
  rows = []
  for patient p in COHORT:
    load W^phi for phi in {preA, preB, task, post}; N = N_p
    masks = build_epi_masks(p); epi = masks.epi_mask; chans = masks.channels
    # full-graph patient trace
    for phi: D[phi] = lrg_ultrametric_condensed(W^phi)
    dT = D[task]-D[preA]; dR = D[post]-D[preB]; rho_full = Spearman(dT,dR)
    s_i  = mean_phi(rowsum(W^phi))                # strength (N,)
    (shaft_i, pos_raw_i) = parse_seeg_label(chans[i]) for each i
    nse_i, pos_i, shaftlen_i  from shafts + epi   # geometry confounds
    for node i in 0..N-1:                          # LNO influence
      Ki = all nodes except i
      for phi: Dk[phi] = lrg_ultrametric_condensed(W^phi[Ki,Ki])
      rho_loo = Spearman(Dk[task]-Dk[preA], Dk[post]-Dk[preB])
      dRho_i  = rho_full - rho_loo
      # per-node row-restricted features (reuse companion F1/F2 defn)
      rho_node_i = Spearman(dT[row i], dR[row i]); coh_node_i = signagree(...)
      rows.append(p,b,i,chan,shaft,y_i, dRho_i, rho_node_i, coh_node_i,
                  s_i, nse_i, pos_i, shaftlen_i)
  # residualize each score on the confound design, WITHIN patient
  for score in {dRho, rho_node, coh_node}:
    for p: score_perp = OLS_residual(score ~ s + nse + pos + shaftlen)  on p's nodes
  # cohort discrimination gate (does the residual axis separate KNOWN epi?)
  for score in {dRho_perp, rho_node_perp, coh_node_perp, s (baseline)}:
    for p with |E_p|>=3: AUC_p = MannWhitney(score on E_p vs U_p)
    cohort: med AUC, |AUC-0.5|, n>0.5, wilcoxon_z(AUC-0.5, greater)
  # occult candidates (descriptive, hypothesis-only)
  for p with |E_p|>=3:
    win = [min,max](dRho_perp on E_p)
    for i in U_p with dRho_perp_i >= quantile_0.90(dRho_perp on U_p) and dRho_perp_i in win:
       emit candidate(p,b,chan, dRho_perp_i, pctile_nonepi, dist_to_nearest_epi, on_epi_free_shaft)
  write per_node CSV, gate CSV, candidates CSV, README
```

---

## 6. Visualization spec

(Deferred unless the gate fires; per `feedback_no_pre_registered_acceptance`,
decide post-hoc.) If built: per-band strip/violin of `dRho_perp` for labeled-epi
vs non-labeled, ≥3 bands, with the strength-baseline residual in a companion
row, and the occult candidates over-plotted as points. **Reading rule:** *if
the labeled-epi residual violin sits indistinguishably on top of the
non-labeled violin, there is no occult axis — the candidates are a
strength/geometry re-read.* No in-axes text; PDF, full vector, `use_lrg_style`.

---

## 7. Connection to prior tools

| prior tool | relation |
|---|---|
| companion scope `2026-06-05_epi-node-trace-marker.md` (audit_79/80) | this is its *discovery* reframe; reuses F1/F2 per-node features; the NEGATIVE known-node classifier is the prior that makes "no occult signal" the expected outcome |
| audit_77 epi-stratified cophenetic | supplies which bands have a cohort trace (α/β/γ_l) → the only bands where an influence score is meaningful |
| audit_80 F8 (deferred LOO-resect) | Δρ_i is the cophenetic analogue of the deferred eigen-subspace LOO-resect — the one new statistic here |
| audit_63 `lrg_ultrametric_condensed` | the LNO rebuild reuses it verbatim (no reinvention) |
| `build_epi_masks` / `parse_seeg_label` (io.patient) | epi labels + shaft geometry (canonical, no local copy; avoids the `np.isin(int,str)` no-op) |

Does not replace the cohort trace results, nor the known-node classifier; it
asks the orthogonal discovery question "do unlabeled nodes behave like epi
beyond strength+shaft?".

---

## 8. Implementation plan

- **audit_90** (`scripts/01_compute/audit/audit_90_occult_node_lno_influence.py`):
  LNO cohort-trace influence Δρ_i + per-node F1/F2 + confound residualization +
  cohort discrimination gate + occult-candidate list, bands α/β/γ_l, n=10.
  Reuses `lrg_ultrametric_condensed` (audit_63), `load_phase_fc`,
  `ensure_half_fcs`, `build_epi_masks`, `parse_seeg_label`, `wilcoxon_z`
  (hypothesis). Outputs `data/audit/occult_node_marker/`:
  `lno_node_influence.csv`, `occult_gate_per_band.csv`,
  `occult_candidates.csv`, `README.md`.
- **No new library helper** is required: residualization is a 4-column OLS
  (`numpy.linalg.lstsq`), used once. If a second caller appears, promote a
  `residualize(y, X)` to `utils/metrics/` (general name).
- **No surrogate ensemble** for the LNO score: the within-patient OLS
  residualization against strength IS the strength control, and is far cheaper
  than per-node matched-strength. (The companion F1/F2 surrogate-z numbers
  already exist in `audit_79`; not recomputed.)

---

## 9. Open questions

- **Nonlinear strength control:** if OLS-residual and rank-partial disagree,
  which is canonical? (Default: report both; if they disagree the negative is
  even safer.)
- **Occult window definition:** top-decile-among-nonlabeled ∩ labeled-epi-span
  is one choice; a surrogate-calibrated window is heavier and only worth it if
  the gate fires.
- **Cross-band fusion of candidates:** weight by each band's gate AUC
  (empirical), never assume β — but only if any band's gate fires.
- **Validation:** the only real validation is surgical-outcome ground truth
  (Engel class at resected vs spared candidate sites). Absent that, this
  direction terminates at "hypothesis-generating, gate-negative".
