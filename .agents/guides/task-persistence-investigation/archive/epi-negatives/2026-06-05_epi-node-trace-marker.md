---
name: epi-node-trace-marker
type: scope
era: COHORT_N10
status: built-negative
created: 2026-06-05
updated: 2026-06-05
pointers:
  - .agents/guides/task-persistence-investigation/2026-05-30_asymmetric-pair-trace-regression.md
  - .agents/guides/task-persistence-investigation/2026-05-08_eigenmode_embedding.md
  - .agents/guides/task-persistence-investigation/2026-05-08_epi-eigenmode-localization.md
  - scripts/01_compute/audit/audit_76_pair_displacement_taxonomy.py
  - scripts/01_compute/audit/audit_77_epi_stratified_cophenetic.py
  - scripts/01_compute/audit/audit_78_epi_stratified_grassmann.py
  - src/lrg_eegfc/utils/surrogate/matched_strength.py
  - src/lrg_eegfc/utils/metrics/spectral.py
---

# Per-node epileptic marker from cross-phase LRG trace behaviour

**Head.** Lift our validated per-*pair* cross-phase trace statistics (cophenetic
ρ_split, audit_63/76; Grassmann subspace participation, audit_66/eigenmode) to a
per-*node* feature vector, expressed as **excess over each node's own
matched-strength surrogate**, and ask whether epileptic nodes occupy a
distinct, band-specific region of that feature space — distinct enough to flag
unlabeled epileptogenic contacts. **This marker is a real finding only if a
trace feature beats BOTH a node-strength baseline AND a probe-adjacency
baseline** in leave-one-patient-out discrimination; the project's own evidence
(per-leaf cophenetic ρ is global-position/strength dominated; epi contacts
cluster on sEEG shafts) makes "the marker is just re-reading strength or shaft
proximity" the default hypothesis to defeat, not a footnote.

> **Sequencing (gated).** This scope splits into a *descriptive* step
> (`audit_79`, built now) that characterizes per-band epi-vs-non-epi behaviour,
> and a *classifier* step (`audit_80`, **deferred**) built only if the
> descriptive step shows at least one band where epi nodes separate beyond
> surrogate spread. Direction 2 is gated on Direction 1 (audit_77/78) and on the
> descriptive step. If neither separates, the honest deliverable is
> "no trace-based marker; epi recoverability is a strength/probe artifact."

> **RESULT (2026-06-05 — `audit_80` BUILT; verdict NEGATIVE).** The classifier
> was built (LOPO + within-patient permutation null R=200 + whole-shaft masking
> + discovery) and the honest deliverable above is the outcome: **no
> trace-based epi marker in any band.** Cross-patient LOPO trace PR-AUC sits at
> prevalence (best β 0.124 vs prev 0.095, **perm-p 0.080 — not significant**),
> adds nothing to the baseline (dAUC2 ≤ 0). **Node strength does NOT transfer**
> (PR ≤ prevalence, perm-p ≈ 1.0): the within-patient "epi hyperconnected"
> effect (audit_79 AUC δ 0.76 / β 0.65) is a patient-specific FC-scale
> phenomenon, not a cross-patient threshold. **The only transferable signal is
> label-free electrode geometry** (`probe_geom` = shaft length + along-shaft
> position; PR 0.163 ≈ 1.7× chance, masked ROC 0.685) — real but **trivial**
> (implantation strategy, not a biomarker). The **whole-shaft masking control
> fired correctly**: the label tautology (`n_same_shaft_epi` etc.) collapses to
> masked ROC ≈ 0.38 (below chance). Discovery candidates
> (`epi_marker_candidates.csv`, `epi_marker_band_agreement.csv`) are therefore
> **low-confidence / hypothesis-generating only** — a trace+strength re-read,
> not validated predictions. Artifacts: `data/audit/epi_marker/clf_*.csv` +
> `README_classifier.md`; figure
> `data/preprint/figures/all_bands/fig_epi_marker_verdict.pdf` (preprint_33).
> Implications: the deferred F8 LOO-resect / F3–F5 magnitude features are NOT
> worth building (the rank/sign + mode-load features already at chance set the
> ceiling); the marker direction is closed unless a new question reframes it
> (e.g. *within-patient* augmentation with partial labels, where the tautology
> is legitimately available — a different clinical scenario).

---

## 1. Notation

- Patient `p` with `N_p` contacts; FC-order node index `i ∈ {0, …, N_p−1}`.
- Band `b ∈ {δ, θ, α, β, γ_l, γ_h}`. Each band is an **independent**
  experiment (band-oriented rule); no pooling across bands.
- Phases `Φ = {rest_pre_A, rest_pre_B, task_test, rest_post}` (split-baseline:
  `rest_pre_A`/`rest_pre_B` are disjoint rsPre half-FCs).
- `W^{φ}_p(b) ∈ [0,1]^{N×N}` symmetric zero-diagonal FC (`imcoh_abs`).
- `D^{φ} ∈ R^{N(N−1)/2}` condensed LRG cophenetic distance at τ=1/λ_max
  (`lrg_ultrametric_condensed`, audit_63). Square form `D^{φ}_{ij}`.
- Row of node `i`: `R(i) = {(i,j) : j ≠ i}`, the `N−1` pairs incident to `i`.
- Per-pair displacements (locked sign convention, audit_63):
  `dT_{ij} = D^{task}_{ij} − D^{preA}_{ij}`,
  `dR_{ij} = D^{post}_{ij} − D^{preB}_{ij}`,
  `dBase_{ij} = D^{preA}_{ij} − D^{preB}_{ij}` (split-half noise floor).
- `E_p ⊆ {0,…,N_p−1}` = epileptic node set (`build_epi_masks(p).epi_mask`);
  label `y_i = 1[i ∈ E_p]`.
- `probe(i)` = sEEG shaft of contact `i`; `contact(i)` = integer along-shaft
  index (`parse_seeg_label`).
- Matched-strength surrogate ensemble `r = 1..R` (R=200, 4-cycle ±δ, seed
  20260511): cached Laplacian eigendecompositions
  `(Λ^{φ}_r, V^{φ}_r)` → cophenetic `D^{φ}_r = cophenetic_condensed_from_eigs(·)`.

---

## 2. Definitions

### 2.1 Per-node trace features (candidate signal)

Each is the **row-`R(i)` restriction** of a validated per-pair statistic — no
new statistic is invented, only restricted to a node's incident pairs.

| # | name (`col`) | formula (node `i`, band `b`) | range | reads |
|---|---|---|---|---|
| F1 | `rho_split_node` | `Spearman(dT[R(i)], dR[R(i)])` | [−1,1] | per-node version of validated ρ_split; rank co-movement of task vs recovery on node `i`'s row. **Headline.** |
| F2 | `coherence_node` | `mean( sgn(dT[R(i)]) = sgn(dR[R(i)]) )` | [0,1] | per-node audit_76 coherence; sign-only, magnitude-free. |
| F3 | `m_task_node` | `median|dT[R(i)]| / median|dBase[R(i)]|` | ≥0 | per-node task displacement vs own split-half floor. **magnitude → strength-fragile (see audit_73).** |
| F4 | `m_rest_node` | `median|dR[R(i)]| / median|dBase[R(i)]|` | ≥0 | per-node recovery displacement vs floor. **strength-fragile.** |
| F5 | `slope_TR_node` | `⟨dT[R(i)],dR[R(i)]⟩ / ⟨dT[R(i)],dT[R(i)]⟩` | R | per-leaf `s_TR` (2026-05-30 §3.5); fraction of node-row task imprint recovered. **audit_73 proved this re-reads strength.** |
| F6 | `taxonomy_node` | audit_76 decision tree on (F1, F3, F4) → {TRACE, ANCHOR, RESET, REST-DRIFT, RANDOM, WEAK-TRACE}, one-hot | categorical | per-node lift of the audit_76 cell reading. Uses surrogate-z'd F3/F4 in the moved/not-moved gate so the class is strength-controlled. |
| F7 | `grass_modeload_node` | `Σ_{k ∈ K_b} v^{post}_k(i)²`, `v^{post}_k` = slowest non-trivial Laplacian eigenvectors of the **rest_post** graph; `K_b` = the band's audit_66-significant k-window | [0,|K_b|] | eigenvector mode-loading: how much of the trace-carrying subspace lives on node `i`. (`grassmann_to_coord_subspace`, spectral.py) |
| F8 | `grass_loo_resect_node` | `chordal_full_vs_resect(V_k_full, V_k_resect, retained=all\{i})` on task_test | ≥0 | leave-one-out resect-Δ: how much removing node `i` reorients the trace subspace. Causal-flavoured complement to F7. **O(N) eigh/band → cost; deferred to classifier step.** |
| F9 | `cross_within_reorg_node` | `median|dT[R_cross(i)]| / median|dT[R_within(i)]|`, `R_within(i)`=same-probe pairs, `R_cross(i)`=cross-probe pairs | ≥0 | does node `i` reorganize with distant tissue vs its own shaft. |

### 2.2 Surrogate normalization (strength-independence by construction)

For each trace feature `f` and node `i`, recompute `f` on every surrogate
realization `r` (using `D^{φ}_r` for cophenetic features; `V^{φ}_r` for
Grassmann), giving `f_i^{surr} ∈ R^R`. Then store:
- `z_i = (f_i^{obs} − mean(f_i^{surr})) / std(f_i^{surr})` — **primary**, signed.
- `p_i = (1 + #{r : f_i^{surr} ≥ f_i^{obs}}) / (R + 1)` — Phipson–Smyth upper-tail.

Because the surrogate fixes node `i`'s strength exactly (4-cycle ±δ preserves
marginals to 1e−4), any discriminative power left in `z_i` is **by
construction not explainable by that node's strength**. This is the formal
control that turns "feature beat the strength baseline" into "feature carries
strength-orthogonal trace information" — the exact test ρ_split passed and
`s_TR` failed at the cohort level (audit_73).

### 2.3 Baseline / artifact-detector features (their OWN models too)

| # | name | formula | role |
|---|---|---|---|
| B1–B3 | `strength_{preA,task,post}` | `Σ_j W^{φ}_{ij}` | node strength per phase — THE dominant confound |
| B4 | `strength_mean`, `strength_cv` | mean, CoV of B1–B3 | strength level + cross-phase volatility |
| B5 | `probe_has_other_epi` | `1[∃ j≠i: probe(j)=probe(i) ∧ y_j=1]` (TRAIN labels only) | spatial-autocorrelation shortcut |
| B6 | `n_same_shaft_epi` | `Σ_{j≠i} 1[probe(j)=probe(i) ∧ y_j=1]` (TRAIN labels) | count of same-shaft epi neighbours |
| B7 | `min_contact_dist_to_epi` | `min_{j: y_j=1, probe(j)=probe(i)} |contact(i)−contact(j)|`, ∞ if none | along-shaft distance to nearest epi (TRAIN labels) |
| B8 | `probe_size`, `contact_index_norm` | shaft length; `contact(i)/max_contact(probe)` | shaft geometry (foci cluster mid-shaft) |

No-consensus-scalar rule: features remain a **vector**; the classifier weights
them. We never pre-pool into one "epi score".

---

## 3. Properties

- **F1/F2 are rank/sign based** → invariant to monotone reparametrization of
  distances; mute (not immune to) the magnitude/strength channel that sank
  `s_TR`. Still require surrogate-z (per-leaf ρ is global-position dominated;
  [[dendrogram_persistence_gate_2026_06_04]]).
- **F3/F4/F5 are magnitude based** → guilty-until-proven: included only as
  surrogate-z, and only in an ablation model `trace-mag`, never the headline.
- **F7 ∈ [0, |K_b|]**, sums to |K_b| over all nodes (eigenvectors orthonormal)
  → it is a *share* of subspace mass; a node can score high only at others'
  expense. Cannot detect a node that moves rigidly with the bulk.
- **F9** is defined only for nodes with ≥1 same-probe and ≥1 cross-probe pair;
  single-contact shafts → NaN (rare).
- **What the marker cannot detect:** (a) epileptogenicity with no cross-phase
  signature (purely static pathology); (b) an epi focus whose trace is identical
  to surrounding healthy tissue (no contrast); (c) band where the cohort trace
  itself is null (δ/θ/γ_h on cophenetic) — there is no trace to read.
- **Complexity:** F1–F6, F9 are O(R·N²) per (p,b) from cached eigs; F7 O(R·N·|K_b|);
  F8 O(R_resect·N·eigh(N)) — the only expensive one, deferred.

---

## 4. Caveats & failure modes

| caveat | mechanism | mitigation |
|---|---|---|
| **Global-position contamination** | per-leaf cophenetic ρ tracks a node's far-bulk hierarchy position ≈ strength | surrogate-z (§2.2) is the formal control; Tier-A masking (§6) is the backstop |
| **Probe spatial autocorrelation** | epi contacts cluster on shafts (NMI 0.65 with putative anatomy); a masked epi node is recoverable from probe identity alone | mandatory probe-only baseline (B5–B8) the trace model must beat; **whole-shaft masking** validation (§6) |
| **epi-only underpower** | n_epi∈{6..30}, Pat_15=0; tiny-graph features noisy | epi-subgraph features (F7/F8 on W[epi]) are exploratory; per-node features here use the FULL-graph LRG row R(i), which is well-powered |
| **Pat_13 leverage** | 30 epi nodes (25%) dominate any pooled epi statistic | patient-stratified tests; LOPO; report leave-Pat_13-out |
| **Pat_15 zero epi** | no positives → cannot train/validate | Pat_15 is a pure-prediction target, never a test fold for AUC; all-negative in training |
| **pre_B surrogate cache gap** | canonical seed-20260511 ensemble omits rest_pre_B | audit_77 extends the cache with deterministic per-cell rng; reused here |
| **B5–B7 label leakage** | neighbour-epi features use labels | computed from TRAIN labels only; whole-shaft masking blanks the shaft so no neighbour leaks |

---

## 5. Pseudocode

```
# audit_79 — descriptive per-node behaviour (BUILT THIS PASS: F1, F2, F7 + baselines)
for band b in BANDS:
  for patient p in COHORT:
    load W^phi for phi in {preA, preB, task, post}; N = N_p
    epi = build_epi_masks(p).epi_mask                       # length N
    D_obs[phi] = lrg_ultrametric_condensed(W^phi)           # observed cophenetic
    dT = sq(D_obs[task]) - sq(D_obs[preA]); dR = sq(D_obs[post]) - sq(D_obs[preB])
    for phi: (Lam_r, V_r) = load_or_compute_surrogate_eigs(p,b,phi,W^phi, R, seed)  # cached
    for r in 1..R: D_r[phi] = cophenetic_condensed_from_eigs(Lam_r[phi], V_r[phi])
    for node i in 0..N-1:
      rows_i = pairs incident to i
      F1_obs = Spearman(dT[rows_i], dR[rows_i]); F2_obs = sign_agree(...)
      F7_obs = sum_{k in K_b} V_post_obs[i,k]^2
      for r in 1..R:                                         # surrogate distribution
        dT_r = sq(D_r[task]) - sq(D_r[preA]); dR_r = sq(D_r[post]) - sq(D_r[preB])
        F1_surr[r] = Spearman(dT_r[rows_i], dR_r[rows_i]); F2_surr[r] = ...
        F7_surr[r] = sum_{k in K_b} V_r[post][i,k]^2
      z_F1_i = (F1_obs - mean(F1_surr)) / std(F1_surr)       # strength-independent
      ... z_F2_i, z_F7_i, p_F1_i ...
      B1..B8 from W^phi strengths + probe(i)/contact(i) + epi (train-label-safe)
      emit row(p, b, i, y_i, z_F1, z_F2, z_F7, B1..B8)
  # per band: epi-vs-non-epi separation, patient-stratified
  for feature f in {z_F1, z_F2, z_F7, B1.., probe-features}:
    per patient: wilcoxon_z(f[epi] vs f[nonepi]); cohort = combine over patients
  emit epi_behavior_per_band[b]

# audit_80 — classifier (DEFERRED; specified, not built)
for band b:
  models = {strength_only(B1-4,B8), probe_only(B5-8), baseline(B1-8),
            trace(zF1,zF2,zF6,zF7,zF8,F9), trace+baseline}
  LOPO: for held-out patient q: train on others, predict q's nodes
        record ROC-AUC, PR-AUC per model; permutation null (shaft-respecting)
  decide: dAUC1 = AUC(trace) - max(AUC(strength), AUC(probe))
          dAUC2 = AUC(trace+baseline) - AUC(baseline)
  masking-Tier-A: blank whole epi-containing shafts, retrain, recovery(trace) vs recovery(probe)
  if dAUC2>0 (perm p<.05) and Tier-A trace>probe and surviving coeffs are z-features:
       verdict = "real strength/probe-orthogonal marker"; emit candidates
  else: verdict = "strength/probe re-read; no trace marker"
```

---

## 6. Visualization spec

- **preprint_32 (this pass):** per-band epi-vs-non-epi violin/strip of each
  surrogate-z trace feature (z_F1, z_F2, z_F7), ≥3 bands shown, baseline
  strength-z and probe features in a companion row. Reading rule: *a band where
  the epi violin sits above the non-epi violin AND above the strength-z violin
  is a marker candidate; a band where epi tracks strength is an artifact.*
- **preprint_33 (deferred):** per-band PR/ROC curves for the five models
  side-by-side (the verdict figure) + Tier-A whole-shaft recovery (trace vs
  probe). PR-AUC headline with per-fold prevalence chance line.

---

## 7. Connection to prior tools

| prior tool | relation |
|---|---|
| audit_76 per-pair taxonomy | F1–F6 are its per-pair stats restricted to node `i`'s row; F6 reuses its decision tree |
| audit_73 asymmetric `s_TR` | F5 is its per-leaf form; audit_73's matched-strength failure is the central design constraint (prefer rank/sign + surrogate-z) |
| audit_66 / eigenmode_embedding (2026-05-08) | F7 is the mode-loading mass behind the Grassmann; eigenmode_embedding's per-contact Procrustes motion `T_i^Y` is a sibling per-node Grassmann feature (candidate F10 later) |
| epi-eigenmode-localization (2026-05-08) | static `m^E_k` mode mass on E_p; F7 is its per-node, trace-window analogue — this scope is the cross-phase extension that 2026-05-08 deferred |
| audit_77/78 (Direction 1) | provide the per-band epi-stratified cohort context that GATES whether a per-node marker is viable |

Does not replace the cohort trace results; it asks the orthogonal node-level
question "do epi nodes carry the trace differently?".

---

## 8. Implementation plan

- **audit_79** (`scripts/01_compute/audit/audit_79_epi_node_behavior.py`, this
  pass): F1, F2, F7 surrogate-z per node × band + baselines B1–B8; per-band
  epi-vs-non-epi patient-stratified `wilcoxon_z`; outputs
  `data/audit/epi_marker/{node_features_per_band.csv, epi_behavior_per_band.csv,
  README.md}`. Reuses `cophenetic_condensed_from_eigs`,
  `load_or_compute_surrogate_eigs`, `build_epi_masks`, `parse_seeg_label`,
  `wilcoxon_z`, the audit_77 `_epi_stratify` helpers.
- **audit_80** (deferred): the classifier + LOPO + masking + discovery
  (`epi_marker_candidates.csv`, `epi_marker_band_agreement.csv`).
- **Library promotion:** the per-node-row pair-stat aggregator (used by F1/F2
  and again by audit_80) promotes to
  `src/lrg_eegfc/utils/metrics/tree_metrics.py` as
  `per_node_row_pair_stat(condensed, i, fn)` — a **general** name (no epi /
  manuscript token), per `feedback_library_names_general`.
- **Feature cache:** `data/cache/epi_marker/Pat_NN/{band}_node_features_*.npz`
  (z + p per node, provenance) so audit_80 re-runs fast. Surrogate eigs stay in
  the shared `matched_strength_surrogate_lrg` cache (reused, not duplicated).

---

## 9. Open questions

- **F8 cost vs value:** run LOO-resect only for bands that pass the descriptive
  gate? default R_resect (50)?
- **K_b for F7:** the audit_66-significant k-window per band — for bands with no
  significant window (δ/θ/γ_h cophenetic), use a fixed fallback (k≤5) or skip F7?
- **Masking granularity:** whole-shaft (Tier A) is primary; is the
  probe-residualized single-node Tier-B worth building, or does Tier-A suffice?
- **Discovery threshold:** what epi-probability percentile defines a "candidate"
  — fixed top-q, or surrogate-calibrated? Decide post-hoc with the user
  (no pre-registered acceptance gate, per `feedback_no_pre_registered_acceptance`).
- **Across-band fusion for candidates:** weight by each band's LOPO PR-AUC
  (empirical), never assume β.

---

## 10. Leverage / influence feature family (2026-06-05, user-directed; `audit_81`)

**Motivation.** The audit_80 features ask whether node `i`'s *own* trace
behaviour looks epileptic (an identity feature) — null cross-patient. The user
observes that *including* epi nodes drastically restructures the dendrogram and
shifts the trace in some patients (= the cohort Direction-1 result, pushed to
node level). That is an **influence / leverage** property — the jackknife of the
tree/trace with respect to node `i` — distinct from F1–F9 identity features.

**5-point preamble.** (1) **Claim:** epi nodes have outsized leverage on the LRG
hierarchy; removing one markedly restructures the rest_post/task dendrogram (L1)
and/or shifts ρ_split (L2), and that leverage flags them. (2) **Null:**
matched-strength surrogate leverage (Stage 2). (3) **Strongest alternative:**
leverage = hubness; high-strength nodes perturb the Laplacian spectrum most, and
absolute strength does NOT transfer across patients (audit_80). (4) **Reach:**
matched-strength removes the strength-magnitude channel; it does NOT by itself
fix cross-patient transfer — the failure mode that killed audit_80. The mitigant
is that L1/L2 are **rank/dimensionless**, so the within-patient *percentile* of
leverage is scale-free and may transfer where dimensioned strength could not.
(5) **Falsification:** if leverage separates epi within-patient but is at chance
in LOPO (like strength), it is dead as a discovery tool.

**Definitions** (condensed cophenetic `D^φ` at τ=1/λ_max; `R(i)` = pairs
incident to `i`; pairs not involving `i` align elementwise between the N-triu
subsequence and the (N−1)-triu — verified).

| # | name | formula (node `i`, band `b`) | reads |
|---|---|---|---|
| L1 | `lev_restruct_node` (per phase φ∈{task,post}) | `1 − Spearman(D^φ_full[¬i pairs], D^φ_resect_i)`, `D^φ_resect_i` = cophenetic of the LRG rebuilt on `W^φ` with row/col `i` deleted | how much removing `i` reorganizes the φ dendrogram. **Direct "drastic effect on dendrogram".** |
| L2 | `lev_trace_node` (signed) | `ρ_split_full − ρ_split(¬i pairs)`, `ρ_split = Spearman(dT, dR)` | signed influence on the cross-phase trace: **>0 = `i` inflates the trace, <0 = `i` masks it** (epi predicted negative per Direction-1). |
| L2a | `lev_trace_abs_node` | `|L2|` | magnitude of trace leverage regardless of sign. |
| (F7) | `grass_modeload_node` | reused from audit_79 | spectral-leverage sibling (eigenvector participation); already LOPO-null. |

**Scale-free transforms** (the LOPO-transferable candidates): within each
(patient, band), the **percentile rank** of L1_post, L1_task, L2a, and the
signed L2 → `pct_*`. These are dimensionless by construction.

**Staged plan.**
- **Stage 1 (`audit_81`, build now):** raw L1/L2 + within-patient percentile;
  (a) epi-vs-non-epi within-patient AUC per band (descriptive, audit_79 idiom);
  (b) Spearman(leverage, strength_mean) per patient — the **confound check**;
  (c) fold the leverage feature sets into the audit_80 LOPO + whole-shaft
  masking harness (import, don't duplicate) and read dAUC vs geometry.
  **No surrogate-z yet** (first pass is NOT strength-controlled; the strength
  correlation is reported alongside, per `feedback_no_pre_registered_acceptance`
  — compute, look, evaluate post-hoc).
- **Stage 2 (gated on Stage 1 transfer):** matched-strength surrogate-z of L1/L2
  to separate pattern (interesting) from hubness (trivial). L1 surrogate-z is
  O(R·N·eigh(N)) per (p,b,φ) — expensive; run only for bands that transfer.

**Cross-check.** Per patient, the epi-set total leverage should track the
audit_77 patient-level `trace(exclude_epi) − trace(full)` shift — validates that
node leverage points at epi tissue where the cohort shift is largest.

**Cost.** L2 is eigh-free (pair-drop on cached `D`). L1 is N resect-eigh per
(patient,band,φ) on {task,post} ≈ 1.4e4 eigh of ~120² → a few minutes. Outputs
`data/audit/epi_marker/node_leverage_per_band.csv` +
`leverage_lopo_per_band.csv` + `leverage_behavior_per_band.csv`.

> **RESULT (2026-06-05 — `audit_81` Stage 1 BUILT; per-node leverage does NOT
> make epi discoverable).** Five of six bands null cross-patient. **One weak but
> real signal: θ** lev_pct LOPO PR-AUC 0.137 vs prevalence 0.095, **perm-p
> 0.005** — but it does NOT beat trivial electrode geometry (0.163; dAUC −0.026)
> and barely survives masking (ROC 0.545). **Its direction is opposite to the
> hypothesis:** the within-patient driver is `lev_trace_abs` AUC≈0.34 (p≈0.008) =
> epi nodes have *LOWER* trace-leverage → epi are **trace-INERT**, not
> trace-disruptive (coheres with Direction-1: trace lives in healthy tissue, epi
> are dead weight on it). Per-node leverage is genuinely strength-orthogonal
> (|median strength corr| ≈0.25–0.29, NOT ≈1) — so this is a real "leverage
> doesn't flag epi", not a strength re-read. **Reconciliation with the real
> cohort epi-shift (Direction-1):** single-node LOO underestimates a GROUP
> effect — with 6–30 epi nodes, removing one barely moves the tree; the epi
> influence is **collective / delocalized**, not localizable to individual
> high-leverage contacts (same signature as the spatially-delocalized trace).
> **Stage 2 (matched-strength surrogate-z) NOT warranted** (θ is weak,
> sub-geometry, counter-direction, on a null-trace band). Two independent
> per-node angles — identity (audit_79/80) and influence (audit_81) — both null
> cross-patient → individual epi contacts are not distinguishable by
> trace-derived properties; the discoverable epi structure is set-level, not
> node-level. Artifacts: `data/audit/epi_marker/{node_leverage_per_band,
> leverage_behavior_per_band,leverage_lopo_per_band}.csv` + `README_leverage.md`.

---

## 11. WITHIN-PATIENT propagator-subspace recovery (2026-06-05, user-directed; `audit_85` + cert `audit_86`)

**Reframe.** audit_79/80/81 all asked a CROSS-PATIENT question (global threshold
transfers from 9 patients to the 10th) — doomed by patient-specific FC scale.
The clinically real and answerable question is WITHIN-patient: *given that some
epi contacts are labelled in this patient, can the rest be rediscovered, and can
unlabelled look-alikes be flagged?* That is anomaly / template recovery in a
rich feature subspace, self-normalised within the patient — scale is no longer a
confound because nothing transfers across patients.

**Propagator features** (the user's "ρ(τ) / information diffusion"). From the
LRG density matrix `ρ(τ) = e^{−τL}/Z` (= lrgsglib `compute_laplacian_properties`
`rho`), computed from the cached `eigh(L^post)` over a per-patient τ-grid
`geomspace(1/λ_max, 10/λ_max, 6)`:
- `rho_ii(τ)` return probability (diffusion localisation at node `i`);
- `Hdiff_i(τ) = −Σ_j p_ij log p_ij`, `p_ij = K_ij/Σ_j K_ij` (diffusion reach);
- communication centrality at τ=1/λ_max: `Tcomm_mean_i = mean_j 1/ρ_ij`,
  `Tcomm_min_i` (how central `i` is in the communication/ultrametric).
Plus the existing per-node trace (audit_79: ρ_split, coherence, mode-load) and
leverage (audit_81: L1/L2) features, and `strength_mean` as the trivial
baseline.

**5-point preamble.** (1) **Claim:** in a within-patient propagator+trace
subspace, epi contacts cluster tightly enough that a held-out epi node (label
hidden) re-ranks to the top via similarity to the other epi, and unlabelled
look-alikes are candidates. (2) **Null:** random epi-set permutation (same size)
→ null recovery AUC; AND a strength-only ablation. (3) **Strongest
alternatives:** (a) it is all hubness (strength); (b) epi cluster on shafts → any
smoothly-shaft-varying feature makes neighbours look alike (geometry). (4)
**Reach:** ablation `strength` vs `diffusion` vs `trace` vs `full` quantifies the
strength share; a **shaft control** recomputes recovery with negatives
restricted to epi-free shafts; leave-one-epi-out keeps the held-out label
unused. (5) **Falsification (of the INTERESTING claim):** if `full` recovery AUC
≈ `strength` AUC ≈ shaft-proximity AUC, there is no propagator-specific epi
subspace — the recovery is real but is hubness/geometry, and must be reported as
such (not a novel diffusion biomarker).

**Method.** Per (patient, band): standardise features within patient
(z per feature → scale-free). Leave-one-epi-out: template = mean z-vector of the
OTHER epi; score every node `−‖z_j − template‖`; recovery AUC = P(held-out epi
ranks above non-epi). Pool over epi nodes and patients per band; permutation
null (R=200 random epi sets). Shaft control: negatives = non-epi on epi-free
shafts only. Discovery: template = mean of ALL epi; rank non-epi; candidates with
3D distance to nearest known epi + epi-free-shaft flag + cross-band agreement.
Pat_15 (0 epi) excluded (no template). Honest expectation: strength already gives
~0.76 within-patient, so the headline number is whether `full > strength` and
whether it survives the shaft control.

**Outputs.** `data/audit/epi_marker/{node_propagator_features.csv,
recovery_per_band.csv, propagator_candidates.csv}` + `README_recovery.md`;
figure `preprint_34_epi_recovery.pdf`.

> **RESULT (2026-06-05 — `audit_85` BUILT; ⚠️ FRAMING-CORRECTED by `audit_86`,
> read the CERTIFICATION below — the headline numbers in this banner used an
> UNFAIR strength baseline and DO NOT stand).**
> Within-patient propagator-subspace recovery WORKS, and it is the right framing
> the cross-patient classifiers (audit_80/81) were missing. Leave-one-epi-out
> recovery AUC: full subspace **β 0.715, α 0.629, low-γ 0.637, θ 0.679** (per
> patient up to 0.85, Pat_07/14). **The propagator ρ(τ) subspace ADDS beyond
> node strength in 5/6 bands** (cross-band median full 0.633 vs strength 0.599,
> Δ+0.060; at **β the diffusion features are the PRIMARY driver: 0.720 vs
> strength 0.622**, trace 0.679). δ is the exception (pure hubness, propagator
> hurts). **Three controls pass:** (1) off-shaft recovery ≈ all-negative
> recovery (β 0.730 vs 0.715) → NOT shaft proximity; (2) permutation null (random
> same-size sets) beaten in 6/9 patients at β & θ; (3) strength-orthogonality is
> proven by Pat_08 — epi there are NOT hubs (strength recovery 0.165) yet the
> diffusion subspace recovers them at 0.758. Honest bounds: within-patient only
> (needs ≥3 labelled epi to build the template — the partial-label clinical
> scenario, exactly "rediscover the rest"); AUC ~0.6–0.72 = useful not
> clinical-grade; heterogeneous (Pat_02 strength-better, Pat_13/30-epi hard).
> Discovery `propagator_candidates.csv`: 947 ranked non-epi/patient; 7 contacts
> top-ranked in ≥4 bands, 2 of them far (>20 mm) from any known epi (Pat_06 U9
> @41 mm, Pat_13 W'1 @26 mm) = the genuinely novel candidates (the rest likely
> field continuation). Artifacts `data/audit/epi_marker/{node_propagator_features,
> recovery_per_band,propagator_candidates,propagator_candidates_agreement}.csv`
> + `README_recovery.md`; figure `preprint_34 fig_epi_recovery.pdf`. **This
> supersedes the "no marker" conclusion: there IS a propagator-subspace epi
> marker, it just lives within-patient, not as a cross-patient threshold.**
> Optional next: matched-strength surrogate-z of the diffusion features to
> certify the β diffusion add is pattern not residual hubness.

> **CERTIFICATION (2026-06-05 — `audit_86`; supersedes the banner above).**
> The matched-strength certification + FAIR strength baseline (within-patient
> Mann-Whitney AUC of strength, not the unfair template-distance used in
> audit_85) overturns the β headline and relocates the real signal:
> **β is HUBNESS** — fair strength recovery 0.738 > propagator raw 0.720 >
> strength-removed (MS-z) 0.606; the audit_85 "diffusion 0.720 > strength 0.622"
> was a template-distance artifact (strength's fair AUC is 0.738, not 0.622).
> **The genuine strength-orthogonal signal is at α, not β:** the MS-z
> (strength-removed) propagator subspace recovers α epi at **0.667 > fair
> strength 0.600 (Wilcoxon p=0.049)** — modest, borderline, but real and
> matched-strength-certified, and it COHERES with the independent finding that α
> uniquely recruits epileptic tissue (α epi↔epi trace +0.404,
> `2026-06-05_epilepsy-headline-and-occult-node-marker.md` §1.3). MS-z recovery
> clears chance (Wilcoxon p<0.05) at θ/α/β/low-γ but only EXCEEDS fair strength
> at α; θ comparable (0.646 vs 0.672); β strength-dominated. Per-feature MS-z
> AUCs depart from 0.5 (β Tcomm_min 0.75, Hdiff ~0.71; α Hdiff ~0.62) but are
> NOT individually significant at n=9 — the signal lives in the combined
> subspace. Cross-patient: still null (audit_80/81). **Honest verdict: within-
> patient epi recoverability is hubness-dominated (strongest β); there is a
> modest, certified, strength-orthogonal, α-specific propagator component; NOT a
> cross-patient marker and NOT a β breakthrough.** Artifacts
> `data/audit/epi_marker/{propagator_ms_features,propagator_ms_certify}.csv` +
> `README_ms_certify.md`; report
> `.agents/reports/2026-06-05_propagator-subspace-epi-recovery.md`.

> **DEPTH-CONFOUND CLOSE-OUT (2026-06-05 — `audit_87` within-shaft + depth
> baseline).** The hide-50%-of-epi within-shaft recovery (audit_87) initially
> looked like a real shaft-independent propagator marker (MS-z recovery >chance
> all bands). It is NOT: along-shaft DEPTH (`contact_index_norm`) recovers hidden
> epi within-shaft at 0.72–0.74 every band, beats the propagator (0.60–0.67), and
> the propagator adds nothing to depth (p≥0.33). Same for the signed magnetic
> Laplacian (audit_88). FINAL: within-patient epi recoverability = hubness +
> along-shaft depth (geometry); no dynamics marker; cross-patient null. Report
> `.agents/reports/2026-06-05_propagator-subspace-epi-recovery.md` §10. Always add
> the along-shaft depth baseline to any within-shaft epi-contact test.
