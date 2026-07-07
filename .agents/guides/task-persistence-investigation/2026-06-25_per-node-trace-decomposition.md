---
name: per-node-trace-decomposition
type: scope
era: IMCOH_ABS × COHORT_N10
status: draft
created: 2026-06-25
updated: 2026-06-25
pointers:
  - .agents/reports/2026-06-25_trace-heterogeneity-handoff.md
  - .agents/reports/2026-06-25_cophenetic-gate-presence-vs-consistency.md
  - scripts/01_compute/audit/audit_63_split_baseline_surrogate.py
  - scripts/01_compute/audit/audit_83_localization_matched_strength.py
  - scripts/01_compute/audit/audit_85_wm_decimation_control.py
  - src/lrg_eegfc/utils/metrics/node_localization.py
  - data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
---

> **Head.** Decompose the cophenetic trace `ρ_split` into per-node contributions
> `T_i` (the mean rank-concordance of all pairs touching node *i*), calibrate each
> node against its own matched-strength null to label it **carrier** / **neutral** /
> **anti-trace**, pool nodes across the cohort, and ask the one question the
> localization atlas (audit_83) never asked: **what do the anti-trace nodes share**
> (tissue / anatomy / epilepsy / sampling-depth / hubness)? If the anti nodes share
> an *a-priori* property **P**, the publishable test of the patient heterogeneity is
> whether excluding nodes **by P** (not by their trace score) raises the cohort
> trace **beyond size-matched random-node decimation** and **out-of-sample**. Only
> then is P the explanator of "who traces". This measure is the negative-space
> complement of the β→OFC carrier localization; it is ρ^coph-only (the Grassmann
> per-node analog is `node_localization.grassmann_trace_contributions`, kept
> distinct on purpose).

---

## 0. Why this measure exists (the heterogeneity problem it serves)

The per-band trace verdict is fragile machinery (internal fact, never a paper
claim). The robust, publishable facts are (a) **a trace exists** (β, gray-cortical,
→OFC) and (b) it **appears more in some bands than others**. But there is large
**per-patient heterogeneity**: between-patient variance ≫ between-band variance
(7.8× cophenet). Two *top-down* explanations of that heterogeneity were refuted
(OFC/network coverage; epileptic-gray contamination — failed the mandatory
decimation control). The PI's directive is to go **bottom-up, per node**: find the
nodes that carry vs oppose the trace empirically, then discover what the opposers
share. The scientific bet (constraint: "the patient is not an axis of variance
unless implant-explained") is that a node-level property P separates carriers from
anti-trace nodes, and that non-tracer patients are non-tracers *because* they
sample more P-tissue. This measure is built to find P if it exists, and to test it
non-circularly if it does.

## 1. Five-point critical preamble (mandatory; `feedback_critical_null_preamble`)

1. **Claim.** The per-patient heterogeneity in the cophenetic trace is carried by a
   node-level structure: each patient has trace-**carrier** nodes (pairs move
   concordantly across task→rest) and **anti-trace** nodes (pairs move
   anti-concordantly), and the anti-trace nodes share an a-priori property **P**
   (anatomy / tissue / epilepsy / depth) such that patients sampling more P-tissue
   trace less.
2. **Null.** Per node, the matched-strength 4-cycle ±δ Laplacian surrogate
   (audit_63 ensemble, R=200, cached at seed 20260511) reconstructed exactly as in
   audit_83: per realization recompute the per-pair rank-concordance on
   strength-preserving surrogate graphs, aggregate to the node. A node's polarity is
   real only if `T_i` exceeds (carrier) / falls below (anti) its own surrogate tail.
   For the **closing test**, the null is **size-matched random-node decimation**
   (audit_85 `--mode node`): rebuilding LRG on any submatrix inflates `ρ_split` for
   *any* K-node removal, so "exclude P → trace rises" means nothing unless it beats
   removing K random nodes.
3. **Strongest alternative the null must cover.** (i) High-strength (hub) nodes
   accumulate larger |concordance| and happen to be the carriers — the
   matched-strength surrogate holds each node's strength fixed and randomizes
   topology, so a strength-driven polarity is reproduced under the surrogate and
   does not survive. (ii) The "P explains it" recovery is just the generic node-count
   inflation of rebuilding LRG on fewer nodes — the random-node decimation null
   isolates exactly this. (iii) P is collinear with along-shaft depth or with node
   strength (both known confounds) — the characterization adjusts for both before
   crediting P.
4. **What the null cannot do.** The matched-strength surrogate conditions on the
   observed per-node strengths; if P is *perfectly* collinear with a strength
   pattern that is itself anatomically structured, the surrogate cannot separate
   "anatomy" from "anatomically-structured strength" — that case is reported as a
   companion `strength_localizes` descriptive (as audit_83 does), not hidden. The
   per-node polarity does not re-test trace existence (locked at §5.2/5.4); it
   conditions on the trace and asks where its carriers and opponents sit.
5. **Falsification.** If the anti-trace nodes (pooled across the cohort) are
   indistinguishable from carrier nodes on every a-priori property after adjusting
   for strength and depth — OR if a candidate P separates them but excluding-by-P
   fails to raise the cohort trace beyond random-node decimation (or fails
   out-of-sample) — then the heterogeneity is **not** node-property-explained and
   that is the honest verdict (reported as such; the PI's prior does not override a
   failed closing test). The β→OFC carrier result and the band result stand
   regardless.

## 2. Notation

| symbol | meaning |
|---|---|
| `p ∈ {2,…,15}` | patient (cohort n=10: 02,03,05,06,07,08,10,13,14,15) |
| `b` | band ∈ {δ, θ, α, β, γ_l, γ_h}; surrogate cache exists at R=200 for {α, β, γ_l} |
| `N = N_p` | number of contacts (nodes) for patient p (Pat_10 task rows [53,54,55] dropped at load) |
| `P = N(N−1)/2` | number of unordered node pairs (length of a condensed cophenetic vector) |
| `Â^φ` | imcoh_abs FC matrix, phase φ ∈ {pre_A, pre_B, task_test, rest_post}; pre_A/pre_B = disjoint rsPre halves |
| `D^φ ∈ ℝ^P` | LRG cophenetic condensed distance vector at τ=1/λ_max (average-linkage on `T=1/ρ`), built canonically (`cophenetic_condensed_from_eigs`) so observed and surrogate share construction |
| `Δ_task = D^tt − D^pre_A`, `Δ_rest = D^post − D^pre_B` | per-pair split-baseline displacement vectors (audit_63) |
| `ρ_split = Spearman(Δ_task, Δ_rest)` | the §5.3 cophenetic trace statistic of record; `>0` = trace |
| `c_p` | per-pair rank-concordance contribution (def. §3.1) |
| `T_i` | per-node raw trace score (def. §3.2), un-demeaned; `mean_i T_i ∝ ρ_split` |
| `t_i = T_i − mean_k T_k` | per-node **relative** concentration (audit_83's quantity at node granularity) |
| `{T_i^{(r)}}_{r=1..R}` | matched-strength surrogate per-node scores (def. §3.3) |
| `π_i ∈ {carrier, neutral, anti}` | per-node polarity label (def. §3.4) |
| `X_i` | node property vector: system, DK region, lobe, gray/WM, epi membership, hemisphere, along-shaft depth `d_i`, strength `σ_i`, coords (x,y,z) |
| `S_P` | the set of nodes carrying candidate property P |

## 3. Definitions

### 3.1 Per-pair rank-concordance `c_p` (reuse audit_83.concordance — to be promoted)

For the P pairs, let `r^t_p = rank(Δ_task_p)` and `r^r_p = rank(Δ_rest_p)` (average
ranks over ties), and centre them: `u_p = r^t_p − (P+1)/2`, `v_p = r^r_p − (P+1)/2`.

    c_p = u_p · v_p          (domain ℝ, the per-pair contribution to Spearman)

With no ties, `ρ_split = (Σ_p c_p) / [ P(P²−1)/12 ]`. This is **rank-based on
purpose**: the matched-strength surrogate produces near-disconnected rewired graphs
whose cophenetic magnitudes explode (`T = 1/ρ` on a vanishing ρ), which corrupts the
raw product `Δ_task·Δ_rest`; the rank concordance is immune and is consistent with
the locked `ρ_split = Spearman` measure (audit_73; audit_83 docstring §concordance;
`feedback_no_partition_metrics_use_rho_coph`). **`c_p > 0` = node-pair moves in the
trace direction; `c_p < 0` = anti-trace direction.**

### 3.2 Per-node trace score `T_i` (the new per-node object)

Let `E_i = { p : i is an endpoint of pair p }`, `|E_i| = N−1`. The raw per-node
score is the endpoint-incidence mean of `c_p`:

    T_i = (1 / (N−1)) · Σ_{p ∈ E_i} c_p          (domain ℝ)

**Conservation.** `Σ_i (N−1) T_i = Σ_p (endpoints of p) c_p = 2 Σ_p c_p`, hence
`mean_i T_i = (2/N) · (1/(N−1)) Σ_p c_p ∝ ρ_split` with a positive constant. The
node decomposition is therefore *conservative*: the average node carries the global
trace, and a patient's low `ρ_split` is literally a node-mean of `T_i` dragged down
by negative-`T_i` (anti) nodes. This is the mechanism the heterogeneity claim rests
on. The **relative** companion `t_i = T_i − mean_k T_k` is audit_83's
`unit_means_from_s` quantity at node granularity (demeaned endpoint-incidence mean),
retained for continuity with the β→OFC localization and for the concentration view.

### 3.3 Matched-strength surrogate per-node null

Reuse the audit_63/audit_83 R=200 ensemble (cached eigs at
`data/cache/matched_strength_surrogate_lrg/{pat}/{band}_{phase}_R200_swap20_seed20260511_imcoh_abs.npz`).
For realization r, recompute per-phase surrogate cophenetic condensed distances
`D^{φ,(r)}` (via `cophenetic_condensed_from_eigs`), form
`c_p^{(r)} = concordance(D^{tt,(r)} − D^{preA,(r)}, D^{post,(r)} − D^{preB,(r)})`,
and aggregate to nodes exactly as in §3.2 → `T_i^{(r)}`. Because independent
per-phase rewiring breaks cross-phase coherence (audit_63 README), `{T_i^{(r)}}` is
centred near 0 for every node regardless of its strength.

### 3.4 Polarity classification `π_i`

Per node, with upper/lower surrogate quantiles `Q^+_i = Q_{0.95}({T_i^{(r)}})`,
`Q^−_i = Q_{0.05}({T_i^{(r)}})` and one-sided matched-strength tail probabilities

    p_i^+ = (1 + #{r : T_i^{(r)} ≥ T_i}) / (R+1)     (carrier evidence)
    p_i^− = (1 + #{r : T_i^{(r)} ≤ T_i}) / (R+1)     (anti evidence)

    π_i = carrier   if T_i ≥ Q^+_i        (equivalently p_i^+ ≤ 0.05)
        = anti      if T_i ≤ Q^−_i        (equivalently p_i^− ≤ 0.05)
        = neutral   otherwise

Continuous polarity (for pooled regression / maps):
`π̃_i = −log10(p_i^+) + log10(p_i^−)` (positive = carrier-leaning, negative =
anti-leaning, magnitude = confidence). No BH across nodes (we are not gating a
per-node claim; we pool — see §3.5). The 5% tails are descriptive thresholds, not a
hypothesis gate.

### 3.5 Anti-node characterization (the discovery step)

Pool every node across the cohort, per band → `n ≈ Σ_p N_p` (~1100–1300 nodes). Each
node has `(π_i or π̃_i, X_i)`. Find the property P that separates anti from carrier:

- **Univariate.** Categorical X (system, DK region, lobe, gray/WM, epi, hemisphere):
  carrier/anti contingency, reported as anti-fraction per level with a within-patient
  permutation reference (shuffle `π_i` within patient, preserving each patient's
  carrier/anti counts — this is the honest null because patient is the blocking
  factor). Continuous X (depth `d_i = contact_index/shaft_length`, strength `σ_i`,
  |z| coord): anti vs carrier distribution contrast.
- **Multivariate, within-patient.** Logistic `1[π_i = anti] ~ X_i` with patient fixed
  effects (or a patient random intercept), to find the property that survives
  **adjusting for strength `σ_i` and depth `d_i`** (the two known confounds — depth is
  the recurring epi-confound, `epi_stratified` §10; strength is the universal one).
  The separating property must be (a) a-priori / independent of the trace score, and
  (b) significant after strength+depth adjustment.

The output of §3.5 is a ranked shortlist of candidate P's with effect sizes — *not* a
verdict. The verdict comes only from §3.6.

### 3.6 Non-circular closing test (the only thing that licenses "P explains it")

Selection bias warning: excluding nodes *because* their `T_i` is low trivially raises
`ρ_split` (the conservation identity §3.2 guarantees it). That is **not** evidence.
The evidence is:

1. The anti nodes share a-priori property **P** (from §3.5).
2. For each patient, build `S_P` = nodes with property P, **exclude `S_P`**, rebuild
   LRG on the `(N−|S_P|)`-node submatrix, recompute `ρ_split^{(−P)}`.
3. Compare `Δρ_p = ρ_split^{(−P)} − ρ_split` to the **random-node decimation null**:
   remove `|S_P ∩ p|` *random* nodes, K≥200 times, build `{Δρ_p^{rand}}`
   (audit_85 `--stratify <P> --mode node`). P passes iff `Δρ_p` exceeds the random
   null — cohort-wide, and **specifically in the non-tracers** (the patients whose
   recovery is the whole point).
4. **Out-of-sample.** Leave-patient-out: define P (the property + threshold) on 9
   patients, apply to the 10th, confirm its recovery is predicted. P that only works
   in-sample is overfit.

Only a property surviving (3) **and** (4) is reported as the explanator. The
pair-class variant (same-graph, decimation-exempt; audit_85 `--mode pairclass`) is
the companion when P is naturally a pair property rather than a node property.

## 4. Properties

- **Range / sign.** `T_i, t_i ∈ ℝ`; sign = trace vs anti direction; `mean_i T_i ∝
  ρ_split` (conservation, §3.2). `π̃_i ∈ ℝ`, sign = polarity, magnitude = confidence.
- **Scale invariance.** Rank-based ⇒ invariant to any monotone reparametrization of
  the cophenetic distances; immune to the surrogate magnitude explosion (§3.1).
- **Strength control by construction.** Polarity is defined *relative to a
  strength-preserving surrogate*, so a node is never a "carrier" merely for being a
  hub. (Companion `strength_dev` column reports the residual collinearity.)
- **Decomposition fidelity.** The per-node scores sum back to the global statistic;
  no information is created. The decomposition is a linear re-attribution of the same
  `c_p` to endpoints.
- **What it CANNOT detect (negative properties — part of the contract):**
  - It does **not** measure subspace / eigenmode persistence — that is the Grassmann
    probe (`grassmann_trace_contributions`), kept **distinct** (constraint 3: the two
    probes must not tell the same story). This tool is pairwise-hierarchy only.
  - It does **not** re-test trace existence; it conditions on it.
  - It cannot localize a trace that lives in *which pairs* but is endpoint-balanced
    (a pure pair-structure with no node concentration washes out under endpoint
    aggregation — caught only by the per-pair / pair-class view, §3.6 companion).
  - It cannot, on its own, distinguish anatomy from anatomically-structured strength
    if the two are perfectly collinear (§1.4).
  - The closing test cannot credit P if `|S_P|` is so large that the submatrix is
    degenerate — guard with a `|S_P|/N` ceiling and report it (no silent caps).

## 5. Caveats & failure modes (each with a mitigation)

| failure mode | mitigation |
|---|---|
| Raw-product magnitude explosion under surrogate (`1/ρ` on rewired near-disconnected graphs) | use the **rank** concordance §3.1 (audit_83 fix), never `dD_task*dD_rest` |
| Endpoint double-counting inflates apparent n | nodes are the unit; pooled tests use a within-patient block null (§3.5), never treat the 2P endpoint incidences as independent |
| Demeaning choice changes which nodes look "anti" | report **both** raw `T_i` (the conservation axis, primary for the heterogeneity mechanism) and demeaned `t_i` (audit_83-compatible concentration); classification uses the surrogate tail on raw `T_i` |
| Selection-bias circularity ("exclude low-T_i → trace rises") | the closing test §3.6 excludes by **a-priori P**, never by `T_i`, and must beat random-node decimation + out-of-sample |
| WM / non-anatomical nodes manufacture cohort effects (sampled by every patient) | gray/WM is itself a candidate P; if WM, the decimation control is mandatory (already known: WM "sharpening" is generic node-count — `white_matter_exclusion_2026_06_08`) so the §3.6 test is exactly the right gate |
| Along-shaft depth confound (epi tests recurringly reduce to depth) | depth `d_i` is always in the §3.5 adjustment; P must survive depth |
| Coordinate quirks (mixed scales in x; some huge, some ~0) | use coords only for descriptive maps / hemisphere split, not as a quantitative property until normalized; prefer `system`/`lobe`/`hemisphere` labels |
| Ties in cophenetic heights (ultrametric has only N−1 distinct merge heights → many tied pairs) | average ranks (scipy `rankdata` default); the Spearman normalization with ties is handled by working from the same `concordance` the locked measure uses |
| Surrogate cache only for {α,β,γ_l} | primary analysis on those three; for δ/θ/γ_h either generate the R=200 cache (canonical seed) or use a lighter within-patient sign-flip / phase-pair permutation null, clearly labeled as a weaker calibration (band result needs only the carrier *counts*, which the lighter null supports) |
| Pat_10 dropped task rows | load via the standard loader (drops [53,54,55]); node indexing stays consistent across phases because all phases use the same 113-channel mask |

## 6. Pseudocode

```
# ---- per-patient, per-band: node scores + polarity ----
for each (p, b):
    load Â^{pre_A}, Â^{pre_B}, Â^{tt}, Â^{post}            # halves from imcoh_halves_fc
    D^φ  <- canonical_cophenet(Â^φ)   for φ in {preA,preB,tt,post}   # eigh(D−W) → cophenetic_condensed_from_eigs
    c    <- rank_concordance(D^tt − D^preA, D^post − D^preB)         # length P  (§3.1)
    for node i in 0..N-1:
        T[i] <- mean over pairs p∈E_i of c[p]                        # endpoint-incidence mean (§3.2)
    t <- T − mean(T)                                                 # relative companion

    # surrogate null (cached eigs)
    load surrogate eigs (R=200) for the four phases
    Tsurr <- zeros(R, N)
    for r in 0..R-1:
        cr <- rank_concordance(Dsurr^tt[r]−Dsurr^preA[r], Dsurr^post[r]−Dsurr^preB[r])
        for node i: Tsurr[r,i] <- mean over p∈E_i of cr[p]
    Qplus[i]  <- quantile(Tsurr[:,i], 0.95)
    Qminus[i] <- quantile(Tsurr[:,i], 0.05)
    pplus[i]  <- (1 + #{r: Tsurr[r,i] ≥ T[i]}) / (R+1)
    pminus[i] <- (1 + #{r: Tsurr[r,i] ≤ T[i]}) / (R+1)
    polarity[i] <- carrier if T[i]≥Qplus[i] else anti if T[i]≤Qminus[i] else neutral
    pi_tilde[i] <- -log10(pplus[i]) + log10(pminus[i])

    # node properties (read-only)
    X[i] <- {system, region, lobe, is_wm, is_epi, hemisphere, depth, strength, x,y,z}
    emit row (p, b, i, T[i], t[i], pplus[i], pminus[i], polarity[i], pi_tilde[i], X[i])

# ---- cohort pooling + characterization (§3.5) ----
pool all rows for band b
for each categorical property C:
    anti_fraction[level] vs within-patient-shuffle reference
for each continuous property V in {depth, strength, |z|}:
    contrast distribution(anti) vs distribution(carrier)
fit logistic 1[anti] ~ X  with patient fixed effects   # P must survive strength+depth
rank candidate P by adjusted effect size

# ---- non-circular closing test (§3.6), per candidate P ----
for each patient p:
    rho0      <- rho_split(p, b)                         # full
    rho_excl  <- rho_split(p, b, exclude = S_P ∩ p)      # rebuild LRG on submatrix
    drho[p]   <- rho_excl − rho0
    drho_rand <- { rho_split(p,b, exclude=K random) − rho0 : K=|S_P∩p|, ×200 }   # audit_85 --mode node
    pass[p]   <- drho[p] > quantile(drho_rand, 0.95)
report cohort pass-rate, non-tracer pass-rate, and leave-patient-out prediction
```

## 7. Visualization spec

1. **Per-patient node-polarity brain map** (one row per patient, columns = bands).
   Glass-brain (nilearn, zoomed to electrode bbox per CLAUDE rule), node colour =
   `π̃_i` on a **diverging saturated** map (carrier = warm, anti = cool, neutral =
   grey) — NOT a near-white-interior map (`feedback_no_near_white_cmaps`); use a
   custom warm/cool with a grey midpoint or turbo split. Reading rule: *a non-tracer
   patient should show spatially clustered cool (anti) nodes; the question the figure
   poses is whether those clusters land on one tissue type / system.*
2. **Cohort anti-fraction-by-property bars** (per band): for each property level, the
   fraction of its nodes that are anti vs carrier, with the within-patient-shuffle
   band. Reading rule: *a level whose anti-fraction clears the shuffle band is the
   candidate P.*
3. **The heterogeneity scatter** (the money figure): x = patient `ρ_split`, y =
   fraction of the patient's nodes carrying property P (e.g. WM fraction, or
   anti-system fraction). A negative slope = "patients sampling more P-tissue trace
   less" — the visual statement of the explanator. One panel per band; β is the
   anchor.
4. **Closing-test forest**: per patient, `Δρ_p` (exclude-by-P) vs the random-node
   decimation P5–P95 span, sorted by baseline `ρ_split`, non-tracers highlighted.
   Reading rule: *blue dot right of the grey span = real recovery beyond node-count.*

All PDFs, `use_lrg_style()`, no suptitle, no watermark, ≥3 patients/bands shown.

## 8. Connection to prior tools

| tool | relation |
|---|---|
| **audit_63** (`ρ_split` generator) | parent statistic; this measure is its exact per-node decomposition (`mean_i T_i ∝ ρ_split`). Reuses the same canonical cophenetic + split-baseline construction. |
| **audit_83** (β→OFC localization) | direct ancestor: its `concordance` (rank per-pair) and `unit_means_from_s` (endpoint-incidence mean, demeaned) ARE steps §3.1–§3.2 at unit granularity. This measure keeps node granularity, adds explicit **anti** classification, and characterizes the **negative space** audit_83 left blank (it only described where the trace *concentrates*, →OFC; never where it *opposes*). The carriers here must recapitulate audit_83's OFC result — a built-in sanity check. |
| **audit_85** (decimation / pair-class control) | the closing-test §3.6 null; extend `--stratify` to whatever P emerges (wm / epi exist; a DK-region or hemisphere P needs a new stratifier — small follow-on). |
| **node_localization.py** | the library home. `cophenet_trace_contributions` there is the **raw-product** (pre-audit_83) version — do NOT use it for surrogate-calibrated polarity (magnitude explosion); promote audit_83's rank `concordance` here as `rank_concordance` and add `canonical_cophenet(W)`. Reuse `epi_keep_mask`, `shaft_of`, `eta2`, `perm_p_eta2`, `nmi`, `NON_ANATOMICAL`. |
| **grassmann_trace_contributions** (node_localization) | the per-node Grassmann analog. **Kept distinct** (constraint 3): cophenet = pairwise-hierarchy carriers; Grassmann = subspace-participation carriers. They may be compared (do carriers coincide?) but never merged into one score. |
| **cross_phase_taxonomy** (audit_105–110) | the per-*pair* trace/reset/anchor decomposition. This is its per-*node* sibling; an anti-trace node ≈ a node whose pairs are reset-dominated. Reuse the taxonomy vocabulary, not its per-pair code path. |

This measure **complements** audit_83 (negative space), **does not replace** the
locked β localization, and **does not subsume** the Grassmann per-node probe.

## 9. Implementation plan

- **Library promotions** (`src/lrg_eegfc/utils/metrics/node_localization.py`,
  same commit as second use per `coding-rules.md`):
  - `rank_concordance(dD_task, dD_rest)` — move audit_83's `concordance` verbatim
    (bit-identical; update audit_83 to import it — audit_83 results are LOCKED, so the
    move must be character-for-character, verified by a reproduction check on one
    cell before committing).
  - `canonical_cophenet(W)` — `cophenetic_condensed_from_eigs(*eigh(diag(W.sum)−W))`
    (audit_83's `_canon_cophenet`).
  - `node_incidence_mean(c, iu_i, iu_j, N, keep=None, demean=False)` — the §3.2
    endpoint-incidence mean at node granularity (generalizes `unit_means_from_s`
    with `unit = node`).
- **Script**: `scripts/01_compute/audit/audit_144_per_node_trace_decomposition.py`
  - reuse `surrogate_cache_path`, `cophenetic_condensed_from_eigs`
    (`lrg_eegfc.utils.surrogate`), `load_channel_regions`, `load_epileptic_nodes`,
    `epi_keep_mask`, `shaft_of` (`node_localization`), `load_fc_matrix`
    (`workflow.fc`), the half-FC cache (`imcoh_halves_fc`).
  - CLI: `--bands` (default β α γ_l = the cached trace bands), `--patients`,
    `--R 200`, `--epi-mode {include,exclude}`.
- **Outputs** (`data/audit/per_node_trace_decomposition/`):
  - `per_node.csv` — one row per (patient, band, node): `T_i, t_i, p_plus, p_minus,
    polarity, pi_tilde` + all properties `X_i`.
  - `cohort_characterization_{band}.csv` — anti-fraction per property level +
    within-patient-shuffle reference + logistic adjusted effects.
  - `closing_test_{P}.csv` — per-patient `Δρ` vs random-node decimation (only after a
    candidate P is identified; this re-invokes audit_85 logic).
  - `README.md` (frontmatter; head-first).
- **Figures** (separate `02_preprint`/`figures_embedded` script, only after the
  characterization lands): the four panels of §7.
- **Env**: `/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python` (the
  `conda run` wrapper errors on a g++ hook).

## 10. Open questions (deferred, flagged load-bearing)

- **Q1.** Polarity tail at 5% vs 10% — descriptive threshold; the pooled
  characterization uses the continuous `π̃_i` so the verdict should be
  threshold-robust, but report carrier/anti counts at both.
- **Q2.** Demean before classification? Decision: **no** for classification (raw `T_i`
  is the conservation axis the heterogeneity mechanism needs), **yes** for the
  concentration companion `t_i` (audit_83 continuity). Revisit if the raw scores are
  dominated by a per-patient additive offset the surrogate doesn't absorb.
- **Q3.** Bands beyond {α,β,γ_l}: generate the R=200 surrogate cache (canonical seed)
  or accept the lighter within-patient null? Default: lighter null for the
  carrier-count band comparison; full cache only if a band becomes load-bearing.
- **Q4.** If multiple properties co-separate (e.g. WM ∧ depth), the closing test runs
  per candidate; report which single P, if any, carries the recovery — do not pool
  into a composite (`never pool metrics into a consensus scalar`).
- **Q5.** Whether the anti-node characterization becomes a manuscript control figure
  (likely) or a methods note — deferred to the PI after the result lands.
```
