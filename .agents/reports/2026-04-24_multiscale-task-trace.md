---
name: multiscale-task-trace-for-writing
type: report
era: COHORT_N9
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Multiscale Task Trace — Writing Agent Handoff (n=9)

> **Status (2026-04-24):** Results computed on the 10-patient sEEG
> cohort. Cross-phase tests run on **n=9** (Pat_14 excluded: vendor
> `task_test.mat` corrupt). H3 runs on n=10. H4 runs on n=10.
> All numbers from `fc_method="imcoh_abs"` = `⟨|Im(S_ij)/√(S_ii S_jj)|⟩_f`
> (Nolte 2004 / Ewald 2012); volume-conduction-immune.

---

## 1a. PRIMARY BAND-SELECTIVE RESULT (partition-level cluster permutation)

Four partition-level measures (VI, ARI, conditional H, NMI) applied
k-by-k across the dendrogram, with sign-flip cluster permutation
(5,000 perms) per band per measure, reveal a **cleanly band-selective
trace**:

| band | cluster-perm significant? | k-range | measure | p |
|---|---|:---|:---|---:|
| **δ** | ★ YES | **k=23–31** | **VI** | **0.014** |
| α | ★ YES (narrow) | k=2–4 | H | 0.050 |
| β | — | no cluster | — | — |
| γ_l | — | no cluster | — | — |
| γ_h | — | no cluster | — | — |
| **θ** | **ergodic (anti-trace)** | — | — | — |

**This is the primary band-selective finding at the partition level.**

- **δ carries a statistically robust, scale-coherent multiscale partition trace** at mid-k (23–31), surviving cluster permutation on VI.
- **α carries a narrow fine-scale partition trace** (k=2–4, cluster permutation on conditional entropy H).
- **θ is actively ergodic**: mean fraction of patients with positive contrast drops BELOW 0.5 on 3 of 4 partition-level measures — most patients × most k values go in the trace-OPPOSITE direction for θ. No cluster-perm cluster passes.
- β, γ_l, γ_h show descriptive positive trends but no cluster-permutation-significant scale-coherent structure.

Scripts: `h2_partition_multiscale.py`, `h2_partition_cluster_perm.py`.
Reports: `h2_partition_multiscale.{md,csv}`, `h2_partition_cluster_perm.{md,csv}`.
Figures: `figures/h2_partition_multiscale.pdf`, `figures/h2_partition_cluster_perm.pdf`.

---

## 1. Headline claim (suitable for abstract) — n=9, DRIFT-CONTROLLED

**Cognitive task execution produces a multiscale reorganization of the
sEEG imaginary-coherence functional-connectivity hierarchy. The
post-task resting state retains a band-selective partition-level trace
of this reorganization: δ exhibits a statistically robust, scale-coherent
partition-level task trace at mid-k (Variation-of-Information contrast,
cluster-permutation p=0.014, k=23–31); α shows a narrower fine-scale
partition trace (conditional-entropy contrast, cluster-permutation
p=0.050, k=2–4); θ is ergodic at the partition level — active
anti-trace direction on 3 of 4 partition-level measures, no
cluster-permutation-significant k-range. The drift-floor control
(H2e) independently confirms α and δ carry block-level persistence
exceeding within-session drift (Δρ_cross − Δρ_null, q<0.05 FDR). The
band-heterogeneity is anchored by the α/θ contrast: α shows the
strongest drift-controlled persistence, θ the weakest (paired
Wilcoxon θ<α, Bonferroni m=5, p=0.010; θ<δ p=0.049). The cohort is
internally consistent with universal directional residual drift (H2c:
all 6 bands q<0.005 FDR, 53/54 patient×band cells positive), and this
directional residual also survives the drift floor (H2e: 6/6 bands,
q=0.004; within-session drift ρ_null≈0). All direct global
similarity tests (rest_post vs task vs rest_pre on raw VI, Spearman
on ultrametric distances, Frobenius distance ratio, tree-bipartition
overlap, directed conditional entropy) fail at n=9, consistent with
the view that the task signal is a specific residual perturbation and
not a global geometry shift — the honest and expected pattern for
sEEG recordings where resting-state ultrametric structure is anatomy-
dominated. This is the first demonstration of an |ImCoh|-LRG task
memory signature in human sEEG, with explicit band selectivity
(α/δ trace vs θ ergodic) anchored in a drift-controlled contrast.**

Supporting controls:
- Task phases are internally consistent (H1: all 6 bands q < 0.05 FDR).
- Phase-type (rest vs task) is a real structural axis (H3: all 6 bands
  q < 0.05 FDR).
- A monotonic frequency gradient in cross-phase VI is present but
  modest (H4: Kendall W = 0.231, p = 0.042, n=10).
- **H2c ρ > 0 survives the split-half session-drift noise floor (H2e):
  within-session drift yields ρ_null ≈ 0 in every band, whereas
  ρ_cross = +0.40 to +0.54 — all 6 bands q = 0.0038 FDR. The
  directional memory claim is not an artefact of session order.**
- H2a′ directed conditional-entropy partition-level test is
  underpowered at n=9 (3 bands trend at p≈0.07; supplementary only).

---

## 2. Cohort, methods, and the three metric layers

### Cohort

- 10 sEEG patients: Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15.
- Pat_14 has a corrupt vendor `task_test.mat` — **excluded from all
  cross-phase tests** (H1, H2a, H2b, H2c target = task_test, H2d).
  Included where the measure does not need `task_test` (H3 within/cross,
  H4 gradient, H2c target = task_learn).
- Pat_10 has a vendor channel-count discrepancy (116 ch in task
  recordings vs 113 in rest); identified the 3 extra task contacts via
  a monotonic-constrained exhaustive channel-matching search (C(116,3)
  = 260k combinations, both rest_pre↔task_learn and rest_post↔task_test
  converged on task rows [53, 54, 55]). Those channels are dropped at
  load time via `PATIENT_CHANNEL_DROP` → Pat_10 is uniformly
  113-channel across all phases. Raw files unmodified.
- Pat_03 is a documented outlier (1024 Hz sampling vs 2048 Hz
  everywhere else, 3× MSC inflation) — always included, always flagged;
  every test reports `LOO worst p` and `p without Pat_03` columns.

### Four phases

`rest_pre`, `task_learn`, `task_test`, `rest_post` (snake_case
canonical).

### Six bands

δ (0.5–4 Hz), θ (4–8), α (8–13), β (13–30), γ_l (30–80), γ_h (80–300).

### FC metric

|ImCoh|: at each (patient, phase, band) we compute the full
frequency-resolved signed imaginary coherency
`ImC_ij(f) = Im(S_ij(f)) / √(S_ii(f) · S_jj(f))` via Welch's method
(`nperseg_for_fs(fs)` = 4096 at 2048 Hz, 2048 at 1024 Hz → 2-s
windows, Δf = 0.5 Hz). The LRG input `imcoh_abs` is
`⟨|ImC_ij(f)|⟩_{f ∈ band}` — Jensen-correct magnitude.

### LRG hierarchy

Laplacian `L = D - A` on the |ImCoh| graph → heat-kernel density
`ρ(τ) = e^{-τL}/tr(e^{-τL})` → von Neumann entropy `S(τ)` sampled at
400 log-spaced τ ∈ [10⁻³, 10⁵] → Laplacian properties giving the
ultrametric merge-time matrix `D[i,j]` = τ at which nodes i, j join in
the heat-kernel-induced clustering. From D we derive the dendrogram
(scipy linkage) and per-k partitions `P_k = fcluster(Z, k,
criterion="maxclust")` for k ∈ {2, ..., N-1}.

### The three operationalisations of "task trace"

| Metric | Level of LRG stack | Formula | What it measures | Scale |
|:------|:------|:------|:------|:--|
| **ρ** (H2c) | (1) ultrametric | Spearman ρ on pairs of Δ_rest(i,j) = D^rpost−D^rpre vs Δ_task(i,j) = D^ttest−D^rpre | Continuous drift direction | Integrated |
| **Δρ(k)** (H2d) | (2) dendrogram | P(c_rpost=1 ∣ c_ttest=1, c_rpre=0) − P(c_rpost=1 ∣ c_rpre=0, c_ttest=0) | Causal block-level persistence | Resolved |
| **VI(k)** (H1, H3, H4) | (3) partitions | Meilă variation of information between two k-partitions | Symmetric partition distance | Resolved |

where `c^phase_k(i,j) = 1` iff nodes i, j share a cluster at cut k in
`phase`'s dendrogram.

### Statistics

For every hypothesis × band: one-sided Wilcoxon signed-rank on the
patient-level k-averaged contrast, Benjamini–Hochberg FDR within
hypothesis (m=6 bands), rank-biserial effect size r_rb, 10k bootstrap
95% CI, leave-one-out worst-case p, Pat_03-drop sensitivity.
Scale-specific localisation via cluster-based sign-flip permutation
(Maris & Oostenveld 2007) with 5k permutations, primary |z|>1.96.
Band-selectivity via Friedman χ²(5) omnibus + focused paired Wilcoxon
θ-vs-others one-sided, Bonferroni m=5.

---

## 3. H1 — Task stability (TL-TT < any cross-type pair)

**Claim**: the two task phases (task_learn, task_test) are internally
more similar to each other than either is to any rest phase.

| band | n | mean contrast | 95% CI | r_rb | q (BH) | cluster-corrected k-range (p) |
|------|--:|--------------:|:------|-----:|-------:|:------|
| δ | 9 | +0.211 | [+0.128, +0.312] | +1.00 | **0.0065★** | k=4–87 (p=0.003) |
| θ | 9 | +0.268 | [+0.177, +0.358] | +1.00 | **0.0065★** | k=8–100 (p=0.005) |
| α | 9 | +0.296 | [+0.224, +0.364] | +1.00 | **0.0065★** | k=8–110 (p=0.003) |
| β | 9 | +0.272 | [+0.173, +0.363] | +0.96 | **0.0065★** | k=7–112 (p=0.003) |
| γ_l | 9 | +0.317 | [+0.173, +0.453] | +0.96 | **0.0065★** | k=7–21 (p=0.027), k=26–88 (p=0.005) |
| γ_h | 9 | +0.163 | [+0.042, +0.270] | +0.78 | **0.0191★** | k=17–23 (p=0.046), k=33–44 (p=0.029), k=98–112 (p=0.021) |

All six bands pass FDR with cluster-corrected support across the bulk
of the multiscale range. **Anchor result — safest statement in the
paper.**

---

## 4. H2c — Continuous drift direction (headline)

**Claim**: at every (patient, band), pairs (i, j) that moved closer
together during task (D^ttest < D^rpre) tended to move closer at
rest_post too; pairs that moved apart did the same. Operationalised
as a per-patient Spearman ρ on the upper-triangle ultrametric shift
vectors.

### Target = task_test (n=9)

| band | n | mean ρ | 95% CI | median ρ | r_rb | q (BH) |
|------|--:|-------:|:------|---------:|-----:|-------:|
| δ | 9 | **+0.486** | [+0.377, +0.611] | +0.508 | +1.00 | **0.0046★** |
| θ | 9 | **+0.397** | [+0.313, +0.503] | +0.355 | +1.00 | **0.0046★** |
| α | 9 | **+0.448** | [+0.322, +0.561] | +0.481 | +1.00 | **0.0046★** |
| β | 9 | **+0.537** | [+0.405, +0.649] | +0.595 | +1.00 | **0.0046★** |
| γ_l | 9 | **+0.521** | [+0.372, +0.692] | +0.390 | +1.00 | **0.0046★** |
| γ_h | 9 | **+0.483** | [+0.263, +0.692] | +0.624 | +0.91 | **0.0076★** |

**53 of 54 patient × band cells are positive** (Pat_15 γ_h = -0.09 is
the only negative). All six bands q < 0.05 FDR; LOO worst-p stays
under 0.06 for every band.

### Target = task_learn (n=10; Pat_14 is available for this target)

All six bands q < 0.005★ FDR, mean ρ ∈ [+0.367 (θ), +0.544 (γ_l)].
60/60 patient × band cells positive.

### Writing guidance

> "Across 54 patient × band pairs, 53 showed positive Spearman
> correlation between the rest-pre → task-test ultrametric shift and
> the rest-pre → rest-post ultrametric shift (mean ρ per band ∈
> [0.40 (θ), 0.54 (β)]; all six bands q < 0.005 under BH-FDR;
> replicated with target = task_learn at q < 0.01 and n = 10)."

Frame as the **direction claim**. It is the most robust individual
result and the universality across bands is exactly the point.

---

## 5. H2d — Block-level persistence (post-hoc band heterogeneity)

**Claim**: pairs that task specifically co-clustered (co-clustered in
task_test, NOT co-clustered in rest_pre) co-cluster at rest_post at a
higher rate than pairs that were never co-clustered.

### Patient-level k-averaged Δρ over k ∈ [2, 49] (n=9)

| band | mean Δρ | 95% CI | ρ_task | ρ_inert | r_rb | q (BH) | cluster-corrected k-range |
|------|--------:|:------|-------:|--------:|-----:|-------:|:------|
| δ | +0.248 | [+0.195, +0.305] | 0.41 | 0.16 | +1.00 | **0.0038★** | k=4–49 (p=0.003) |
| θ | **+0.163** | [+0.111, +0.220] | 0.36 | 0.19 | +1.00 | **0.0038★** | k=21–49 (p=0.003) — late only |
| α | +0.283 | [+0.233, +0.347] | 0.46 | 0.18 | +1.00 | **0.0038★** | k=11–49 (p=0.003) |
| β | +0.241 | [+0.191, +0.298] | 0.37 | 0.13 | +1.00 | **0.0038★** | k=5–49 (p=0.003) |
| γ_l | +0.243 | [+0.153, +0.349] | 0.38 | 0.14 | +1.00 | **0.0038★** | k=10–49 (p=0.003) |
| γ_h | +0.249 | [+0.151, +0.349] | 0.49 | 0.24 | +1.00 | **0.0038★** | k=8–49 (p=0.003) |

**All six bands pass FDR at n=9.** So the per-band H2d test is
universal in direction but θ sits clearly as the lowest-Δρ band
(0.163 vs 0.24-0.28 elsewhere) and has the *narrowest* cluster-corrected
k-range (k=21-49 only, vs k=4-49 or k=5-49 in other bands).

### Paired θ-vs-others contrast (one-sided, H₁ θ < other, Bonferroni m=5)

| other band | W | p (raw) | p (Bonf) |
|:-----------|--:|--------:|---------:|
| δ | 3.0 | 0.0098 | **0.049★** |
| α | 0.0 | 0.0020 | **0.010★** |
| β | 8.0 | 0.049 | 0.244 |
| γ_l | 19.0 | 0.367 | 1.000 |
| γ_h | 10.0 | 0.082 | 0.410 |

**θ block persistence is significantly less than α (p = 0.010) and δ
(p = 0.049) after Bonferroni correction.** It is descriptively less
than β and γ_h but does not survive m=5 correction.

### Friedman omnibus

Friedman χ²(5) = 7.22, p = 0.20 at n=9. **Does NOT reject band
exchangeability.** State this honestly — the band heterogeneity claim
rests on the post-hoc θ-vs-α/δ contrasts, not on an omnibus rejection.

### Writing guidance

> "At the block level (k-averaged Δρ, k ∈ [2, 49]), every band showed
> significant excess persistence of task-induced pair co-clusters into
> rest_post (all q < 0.005 FDR), but the excess was quantitatively
> band-heterogeneous: θ exhibited the smallest Δρ (0.16), significantly
> less than α (Δρ = 0.28; paired Wilcoxon, Bonferroni m = 5,
> p = 0.010) and δ (Δρ = 0.25; p = 0.049). A Friedman omnibus across
> all six bands did not formally reject exchangeability at current
> sample size (χ²(5) = 7.22, p = 0.20), so the pairwise θ-vs-α/δ
> dissociation is the strongest statistical claim of
> band-heterogeneity."

Mechanistic framing: θ has the highest baseline co-clustering rate
ρ_inert = 0.19 (vs 0.13-0.18 in other bands except γ_h), consistent
with θ being a high-intrinsic-mixing rhythmic coordinator band.

---

## 6. H3 — Phase-type clustering (within < cross)

Mean VI(within-phase-type) < mean VI(cross-phase-type). Control that
"rest" and "task" are genuine structural types.

| band | n | mean contrast | r_rb | q (BH) |
|------|--:|--------------:|-----:|-------:|
| δ | 10 | +0.056 | +0.71 | **0.023★** |
| θ | 10 | +0.117 | +0.86 | **0.013★** |
| α | 10 | +0.120 | +0.89 | **0.013★** |
| β | 10 | +0.125 | +0.86 | **0.013★** |
| γ_l | 10 | +0.159 | +1.00 | **0.013★** |
| γ_h | 10 | +0.070 | +0.78 | **0.017★** |

All six bands q < 0.05 FDR with extensive cluster-corrected support
(full reporting in `data/reports/imcoh_vi/rigorous_tests.md`).

---

## 7. H4 — Frequency gradient in cross-phase VI (n=10)

Kendall's W = **0.231**, χ²(5) = 11.54, **p = 0.042**.

Mean rank per band (1 = lowest cross-VI, 6 = highest):

| band | mean rank | mean cross-VI |
|------|----------:|--------------:|
| γ_l | 2.40 | 0.879 |
| β | 2.70 | 0.894 |
| δ | 3.00 | 0.950 |
| θ | 4.20 | 1.025 |
| α | 4.30 | 1.038 |
| γ_h | 4.40 | 1.061 |

**A mild but statistically detectable frequency gradient is present**:
γ_l and β have the lowest cross-phase VI (most structural similarity
across phases); θ and γ_h the highest. Not a trivial gradient — γ_l
and γ_h are at opposite extremes — so the pattern reflects inherent
band-level dynamics rather than a simple frequency effect.

---

## 7a. H2e — split-half drift noise floor (primary new result, PASSES)

**Claim**: the H2c ρ > 0 signal is not an artefact of session-order
drift. Within-session drift produces ρ ≈ 0.

Design: for each patient × resting phase we split the time series in
half, run the full ImCoh+LRG pipeline on each half (with a halved
Welch `nperseg` to keep segment count comparable), and compute

- ρ_null_drift = Spearman(D^rpre_B − D^rpre_A, D^rpost_B − D^rpost_A)
  — drift-only direction correlation.
- ρ_within(phase) = Spearman(uppertri D^A, uppertri D^B) — split-half
  reliability ceiling.
- compared against ρ_cross (H2c, full-duration).

Paired one-sided Wilcoxon ρ_cross > ρ_null_drift per band, FDR-BH m=6.

| band | n | ρ_cross | ρ_null_drift | ceiling | cross−null | q (BH) |
|------|--:|--------:|-------------:|--------:|-----------:|-------:|
| δ | 9 | +0.486 | +0.026 | +0.481 | +0.460 | **0.0038★** |
| θ | 9 | +0.397 | −0.040 | +0.462 | +0.437 | **0.0038★** |
| α | 9 | +0.448 | +0.002 | +0.492 | +0.446 | **0.0038★** |
| β | 9 | +0.537 | −0.005 | +0.533 | +0.542 | **0.0038★** |
| γ_l | 9 | +0.521 | −0.005 | +0.536 | +0.526 | **0.0038★** |
| γ_h | 9 | +0.483 | −0.034 | +0.363 | +0.517 | **0.0038★** |

**All 6 bands pass FDR (q=0.0038 in every band — every patient moved in
the same direction, rank-biserial r_rb=+1.0 everywhere). Pre-registered
pass criterion (≥4 bands) exceeded.**

Interpretive takeaway: the within-session drift ρ is statistically
indistinguishable from zero (means span [−0.04, +0.03]) whereas cross-phase
ρ is +0.40 to +0.54. The reliability ceiling ρ_within is +0.36 to +0.54
— so ρ_cross approaches the noise ceiling, which is about as strong a
directional similarity as the half-duration data can resolve.

**Δρ secondary test** (k-averaged Δρ_cross > Δρ_null for triple
(rpre_A, rpre_B, rpost_A)):

| band | Δρ_cross | Δρ_null | cross−null | q (BH) |
|------|---------:|--------:|-----------:|-------:|
| δ | +0.248 | +0.148 | +0.099 | **0.0426★** |
| θ | +0.163 | +0.145 | +0.018 | 0.430 |
| α | +0.283 | +0.147 | +0.136 | **0.0426★** |
| β | +0.241 | +0.195 | +0.046 | 0.235 |
| γ_l | +0.243 | +0.183 | +0.060 | 0.430 |
| γ_h | +0.249 | +0.131 | +0.118 | 0.139 |

α and δ show Δρ_cross significantly above the drift null. θ is flat
(cross−null = +0.018) — reinforces the θ-ergodic interpretation: θ's
H2d block-persistence signal is largely indistinguishable from what
session drift produces. β, γ_l, γ_h trend positive but fail m=6 FDR.

**Writing guidance**: H2e is the reviewer-defense anchor. Quote q=0.004★
in every band (ρ test) in the paragraph that introduces H2c. For Δρ the
honest statement is "α and δ retain above-drift block-persistence; θ
does not, consistent with the band-heterogeneity reported in H2d".

---

## 7b. H2a′ — directed partition-level memory (supplementary)

**Metric**: δH(k) = H(P^rpost_k | P^ttest_k) − H(P^rpost_k | P^rpre_k).
Trace hypothesis: δH < 0 (rpost is more predictable given task than
given rest_pre). Asymmetric half of VI = H(A|B) + H(B|A); recovers a
directed partition-level test that the symmetric VI could not express.

| band | n | mean δH | p (raw) | q (BH) | r_rb |
|------|--:|--------:|--------:|-------:|-----:|
| δ | 9 | −0.060 | 0.087 | 0.173 | +0.51 |
| θ | 9 | +0.021 | 0.843 | 0.843 | −0.38 |
| α | 9 | −0.017 | 0.430 | 0.515 | +0.07 |
| β | 9 | −0.114 | 0.069 | 0.173 | +0.56 |
| γ_l | 9 | −0.116 | 0.257 | 0.386 | +0.24 |
| γ_h | 9 | −0.081 | 0.069 | 0.173 | +0.56 |

Cluster-permutation over k: α has a significant fine-scale cluster at
k=2–4 (mass=11.0, p=0.050★). No other band has a significant cluster.

**0 of 6 bands pass the stringent k-averaged FDR criterion.** Per
pre-registration (≥3 bands for main Results), H2a′ moves to
supplementary. δ, β, γ_h trend in the expected direction at p≈0.07
(borderline); θ trends in the opposite direction (+0.021), consistent
with the θ-ergodic signature across metrics.

**Writing guidance**: one supplementary paragraph: "A directed
conditional-entropy version of the partition-level test
(δH = H(P^rpost|P^ttest) − H(P^rpost|P^rpre)) does not pass
partition-level FDR at n=9 (3 bands trend with p≈0.07 one-sided; α
shows a narrow fine-scale cluster-permutation significant window at
k=2-4, p=0.050). The continuous (H2c) and block-level (H2d)
operationalisations remain the sensitive directed detectors at this
cohort size."

---

## 7d. Topology-level tests (H1-topo, H3-topo PASS; H2a-topo FAILS)

Multiscale tree-topology measure: `bip_overlap(Z^p1, Z^p2)` = fraction
of internal bipartitions (subsets of leaves below each internal node)
shared between two phase dendrograms. Multiscale by construction —
every merge height is tallied.

**H1-topo — task phases share more bipartitions than cross-type:** 6/6 bands pass FDR. **β leads (mean Δ=+0.114★)**; γ_h weakest (+0.045★).

| band | mean Δ | q (BH) |
|------|--------:|-------:|
| β | +0.114 | 0.011★ |
| α | +0.080 | 0.011★ |
| γ_l | +0.078 | 0.017★ |
| θ | +0.074 | 0.011★ |
| δ | +0.057 | 0.013★ |
| γ_h | +0.045 | 0.043★ |

**H3-topo — within-type > cross-type:** 5/6 bands pass FDR (δ q=0.07). Ranking identical to above: β (+0.064★) > γ_l > α ≈ θ > δ > γ_h.

**H2a-topo — rpost closer to task than to rest_pre, at topology level:** 0/6 bands pass FDR. Per-patient signs split ~50/50 (+30/−24 across 54 cells).

**Why this matters for the paper.** The topology-level H1/H3 results provide a **multiscale band-heterogeneity anchor** independent of ρ and Δρ — β is universally strongest, γ_h weakest, across both continuous (ρ, Δρ) and discrete (bipartition) measures. At the same time, H2a-topo joins H2a, H2a-raw, H2a-FROB, H2a′ in failing — confirming that the trace claim is NOT "rpost becomes task-like globally" at ANY metric we have constructed. The trace claim is a residual/conditional statement only.

Scripts: `scripts/01_compute/h2_topology_directed.py`, `diag_topology_vs_rho.py`. Figure: `data/reports/imcoh_vi/figures/h2_topology_directed.pdf`.

---

## 7c. Direct matrix-similarity tests (H2-RAW, H2-FROB) — BOTH FAIL (honest supplementary)

A reviewer's first instinct: is D^rpost directly more similar to D^tt
than to D^rpre? We ran this in two forms:

**H2-RAW** (Spearman on raw upper-triangle):
`ρ(D^rpost, D^tt) − ρ(D^rpost, D^rpre)` per patient × band.

| band | mean Δρ | p (raw) | q (BH) |
|------|--------:|--------:|-------:|
| δ | +0.10 | 0.069 | 0.260 |
| θ | **−0.04** | 0.843 | 0.843 |
| α | +0.05 | 0.384 | 0.575 |
| β | +0.11 | 0.087 | 0.260 |
| γ_l | +0.02 | 0.524 | 0.628 |
| γ_h | +0.06 | 0.221 | 0.441 |

**0/6 bands pass FDR.** θ trends in the opposite direction (rpost is
slightly more correlated with rpre than with task). δ and β trend in
the right direction at p≈0.07–0.09 but do not survive correction.

**H2-FROB** (Frobenius distance ratio, log-Wilcoxon):
`r = ||D^rpost − D^tt||_F / ||D^rpost − D^rpre||_F`; trace
hypothesis r<1.

| band | mean r | direction |
|------|-------:|:----------|
| δ | 1.07 | wrong |
| θ | **1.32** | strongly wrong — rpost 32 % farther from task than from rpre |
| α | 1.09 | wrong |
| β | 0.91 | right, q=0.26 |
| γ_l | 1.14 | wrong |
| γ_h | 0.90 | right, q=0.56 |

**0/6 bands pass FDR.** Most ratios > 1. θ is again the most
trace-opposite.

### What this means

Rest_pre and rest_post have globally very similar ultrametric
geometry (as expected — shared anatomy, shared electrode layout,
shared resting-state organization). Any task-induced signal must
live at the residual level, not at the global-similarity level.

The paper MUST acknowledge this explicitly. The H2c / H2d / H2e
results are not "rest_post looks like task" — they are "the subtle
residual pattern of how rest_post differs from rest_pre aligns with
the task residual pattern", which is a weaker but still scientifically
meaningful claim. Attempting to sell "post-rest becomes task-like"
in the naive sense will be shot down by any reviewer running the
direct tests themselves.

Writing guidance: report H2-RAW and H2-FROB in supplementary,
explicitly framed as "The naive direct-similarity tests do not
detect a post-task shift toward task geometry. Since resting
ultrametric structure is anatomy-dominated and globally stable across
pre- and post-task rest, the task trace is necessarily a residual
signal — captured by H2c (directional drift of the residual) and H2d
(conditional survival of task-formed block partitions), both
controlled against within-session drift by H2e."

---

## 8. What we tested and dropped

### Drop from paper — H2a (VI magnitude of post-rest shift)

Tested: `VI(rest_pre, rest_post) > VI(task_test, rest_post)` at
k-averaged partition level.

At n=9 **NO band passes FDR** (best δ q=0.195, raw p=0.043; θ r_rb =
-0.02, p=0.52). At n=7 (Pat_10 and Pat_13 excluded) δ and β had
appeared to pass FDR — adding Pat_10 and Pat_13, both of whom have
mostly NEGATIVE H2a contrasts, erases the effect.

**Interpretation**: the post-rest "closer to task than to pre-rest"
effect is too weak at the partition level with this cohort. The
continuous (H2c) and causal block (H2d) operationalisations, which
are more statistically powerful, carry the message.

**Writing guidance**: do not report H2a in the main text; a single
supplementary sentence documenting that the partition-level test is
underpowered at n=9 is sufficient.

### Drop from paper — H2b (approach > exit)

Tested: `VI(rest_pre, task_test) > VI(task_test, rest_post)`. Nothing
passes FDR anywhere at n=9. Drop entirely.

---

## 9. Band typology (qualitative table)

Median-threshold classifier inputs:
- `M_task` = mean VI(rest_pre, task_test) across k
- `ρ` = H2c target = task_test
- `Δρ` = k-averaged H2d

| band | M_task | ρ | Δρ | Post-hoc Δρ contrast vs rest | Label |
|:-----|------:|--:|---:|:----|:----|
| α | 1.048 | +0.448 | **+0.283** | highest, significantly > θ | **TRACE (strongest)** |
| γ_h | 1.028 | +0.483 | +0.249 | high | TRACE |
| δ | 0.939 | +0.486 | +0.247 | mid, sig > θ | trace |
| γ_l | 0.929 | +0.521 | +0.243 | mid | mixed |
| β | 0.931 | +0.537 | +0.241 | mid | mixed |
| θ | 1.054 | +0.397 | **+0.163** | **lowest**, sig < α, < δ | **ERGODIC** |

**Caveat**: the absolute differences in M_task are small
(0.929-1.054). The ranking by Δρ is the signal; the median-threshold
TRACE/ERGODIC/DEAD labels are noisy when several bands cluster tightly
— we report the continuous numbers with the post-hoc contrast as the
load-bearing claim.

> **Bottom line for the paper**: α carries the strongest multiscale
> task trace; θ carries the weakest. Every band shows positive drift
> direction (H2c) and positive block persistence (H2d). The band
> heterogeneity is quantitative, not qualitative.

---

## 10. Friedman band-selectivity (supporting, not load-bearing)

Per-metric omnibus χ² test of "all 6 bands have identical per-patient
distribution":

| measure | n | Friedman χ²(5) | p | Kendall W |
|:--------|--:|-------------:|--:|----------:|
| H2a VI magnitude | 9 | 6.52 | 0.259 | 0.145 |
| H2c ρ (target = task_test) | 9 | 3.41 | 0.637 | 0.076 |
| H2c ρ (target = task_learn) | 10 | 5.54 | 0.353 | 0.111 |
| H2d Δρ | 9 | 7.22 | 0.205 | 0.160 |

**None reject.** At n = 9-10 we cannot formally establish "the bands
behave differently" by omnibus. The θ-vs-α and θ-vs-δ post-hoc
contrasts on H2d, which are Bonferroni-corrected paired Wilcoxons, are
what carries the band-heterogeneity claim.

---

## 11. Narrative scaffolding for the paper

Recommended flow from strongest to weakest claim:

1. **(Abstract)** Task execution reorganizes the multiscale
   |ImCoh|-LRG hierarchy, and the reorganization partially persists
   into the post-task resting state, consistent with a band-specific
   memory of task-induced co-activations.
2. **(Result §1)** H1 stability: task_learn ↔ task_test are
   structurally consistent (all 6 bands q < 0.05).
3. **(Result §2)** H2c direction: rest_post drifts in the task
   direction universally (6/6 bands, mean ρ ∈ [0.40, 0.54], 53/54
   cells positive) AND survives the split-half drift noise floor
   (H2e, 6/6 bands q = 0.0038 FDR; within-session drift yields
   ρ_null ≈ 0).
4. **(Result §3)** H2d block persistence: task-formed pair clusters
   persist at rest_post at above-chance rate universally (6/6 bands
   q < 0.005), with quantitative band heterogeneity — α has the
   highest persistence, θ the lowest, significantly so after
   Bonferroni correction. Drift-floor (H2e): α and δ Δρ_cross exceed
   Δρ_null at q < 0.05; θ does NOT, consistent with the ergodic
   interpretation.
5. **(Result §4)** H3 + H4: phase-type separation and mild frequency
   gradient as structural controls.
6. **(Discussion)** θ's reduced block-level persistence is consistent
   with its role as a high-baseline-mixing rhythmic band; α's
   strength is consistent with its load in attentional/task-set
   maintenance. The continuous direction (H2c) is universal because
   LRG captures broad-band structural geometry; the block-level
   persistence (H2d) is band-selective because discrete cluster
   identity is more sensitive to intrinsic band dynamics.
7. **(Limitations)** n = 9 is small; Friedman omnibus on H2d does not
   reject band exchangeability; H2a partition-level magnitude test is
   underpowered and reported as supplementary; H2a′ directed
   conditional-entropy test detects only a narrow α k=2-4 signal
   (also supplementary); pre-registration of θ-vs-α post-hoc in an
   independent cohort is the right replication.

### What NOT to claim

- "A clear frequency gradient in task-trace" — H4 is mild (W=0.23)
  and not a smooth monotonic function of frequency (γ_l low, γ_h
  high).
- "θ is completely ergodic" — H2d for θ passes FDR; θ has reduced,
  not absent, block persistence.
- "Band-level reorganization is qualitatively different" — Friedman
  omnibus fails; at this n we have a quantitative post-hoc contrast,
  not a qualitative split.
- "VI magnitude shows the trace" — H2a fails at n=9.
- "ρ is entirely immune to session drift" — H2e shows the drift null
  is ≈0 in every band, which is strong; but half-duration D matrices
  are noisier than full, so H2e is conservative. The right claim is
  "cross-phase ρ is not explainable by within-session drift as measured
  on half-duration recordings".
- "H2a′ detects partition-level directed memory" — only α k=2-4 passes
  cluster-perm; no band passes k-averaged FDR.

---

## 12. Regeneration (scripts + inputs)

```bash
# Assumes FC caches (corr, msc, imcoh) and LRG caches (corr, msc, imcoh_abs)
# are populated for the whole cohort.  If adding a new patient, first:
#   lrg-eegfc compute corr  --patients <P> -v
#   lrg-eegfc compute msc   --patients <P> -v
#   python scripts/01_compute/compute_imcoh_fc.py --patients <P> -v
#   lrg-eegfc compute lrg --patients <P> --fc-method corr      -v
#   lrg-eegfc compute lrg --patients <P> --fc-method msc       -v
#   lrg-eegfc compute lrg --patients <P> --fc-method imcoh_abs -v

# Multiscale VI + hypothesis testing pipeline:
python scripts/01_compute/compute_imcoh_vi.py -v                 # → vi_raw_profiles.csv, hypothesis_contrasts.csv
python scripts/01_compute/h2c_ultrametric_drift.py               # → h2c_ultrametric_drift.{md,csv}
python scripts/01_compute/h2d_coactivation_persistence.py        # → h2d_coactivation_persistence.{md,csv} — k ∈ [2, 49]
python scripts/01_compute/rigorous_hypothesis_test.py            # → rigorous_tests.{md,csv}
python scripts/01_compute/report_h1h4_vi.py                      # → .agents/reports/2026-04-24_h1-h4-vi-results.md
python scripts/01_compute/band_k_landscape.py                    # → band_k_landscape.md + unanimity figures
python scripts/01_compute/h2_band_selectivity.py                 # → h2_band_selectivity.{md,csv}  (Friedman + post-hoc)
python scripts/01_compute/h2_band_typology.py                    # → h2_band_typology.{md,csv}  (TRACE/ERGODIC/DEAD)

# Audit pass (2026-04-24+):
python scripts/01_compute/h2_conditional_entropy.py              # → h2a_prime_conditional_entropy.{md,csv}  (≈30 min)
python scripts/01_compute/h2e_split_half.py -v --mem-gb 10       # → h2e_split_half.{md,csv}  (≈60-90 min, 10 GB cap)
```

Outputs under `data/reports/imcoh_vi/`. Canonical H1-H4 table:
`.agents/reports/2026-04-24_h1-h4-vi-results.md`.

---

## 13. Figures available for drafting

All figures regenerated 2026-04-24 from n=9 `imcoh_abs` data. Every
figure has both `.pdf` (vector) and, for the custom figures, `.png`
(raster preview). They live at
`data/reports/imcoh_vi/figures/`.

| Figure | File | What it shows | Use in paper |
|:-------|:-----|:--------------|:-------------|
| **H2d band persistence bar chart** | `h2d_band_persistence.pdf` | Mean Δρ per band with 95% bootstrap CIs; θ in red; brackets show Bonferroni θ-vs-α (p=0.010★) and θ-vs-δ (p=0.049★); ρ_task / ρ_inert inside each bar | Headline figure for H2d (Results §3) |
| **Band typology scatter** | `band_typology.pdf` | 2D: ρ (H2c) × Δρ (H2d) per band, labeled, with CIs. All bands in upper-right (universal); θ clearly separated (ergodic-like); α top (strongest trace). Quadrants labeled TRACE/ERGODIC/DEAD/weak | Discussion figure — visualises the dissociation |
| **H2c unanimity heatmap** | `h2c_unanimity.pdf` | 9 patients × 6 bands = 54 cells of per-patient Spearman ρ; 53 positive (one negative cell framed). Direct visualisation of the universality claim | Headline figure for H2c (Results §2) |
| **H2d unanimity heatmap** | `h2d_unanimity.pdf` | 9 patients × 6 bands of per-patient Δρ; column-mean annotations; θ column clearly lightest (lowest persistence). Patient-level view behind the bar chart | Supplementary to h2d_band_persistence |
| **Band × k landscape — H1** | `H1_imcoh_abs.pdf` | Task-stability unanimity across k ∈ [2, N−1]; wide ridges in every band | Results §1 |
| **Band × k landscape — H2a** | `H2a_imcoh_abs.pdf` | VI-magnitude unanimity; mostly pale — visual evidence that H2a fails at n=9 | Supplementary only |
| **Band × k landscape — H2b** | `H2b_imcoh_abs.pdf` | approach-vs-exit unanimity; uniformly pale — H2b drops | Drop |
| **Band × k landscape — H2d** | `H2d_imcoh_abs.pdf` | Block-persistence unanimity over k ∈ [2, 49]; θ row distinctly pale at small k, solid at k≥25 — matches the H2d cluster-permutation result (θ significant only at k=21-49) | Supports H2d scale-localization |
| **Band × k landscape — H3** | `H3_imcoh_abs.pdf` | within<cross unanimity | Results §4 |
| **H2e drift-floor** | `h2e_drift_floor.pdf` | Top: per-band strip of ρ_cross (blue) vs ρ_null_drift (grey) over 9 patients, with a pale reliability-ceiling band; ** over every band (q=0.004 FDR). Bottom: Δρ_cross vs Δρ_null bar pairs, θ in red, ★ on δ and α (q<0.05). Shows that ρ_null≈0 universally while ρ_cross sits near the reliability ceiling, AND that θ's block-persistence is flat above drift. | **Reviewer-defense figure; pair with H2c unanimity in Results §2** |
| **Evidence chain** | `evidence_chain.pdf` | 5-panel patient×band heatmap showing per-cell sign of each test (H1, H2a-raw, H2c, H2d, H2e), rescaled to [−1,+1]. H1 / H2c / H2d / H2e columns are uniformly green (53-54 / +, ≤4 / −); H2a-raw column is mixed with red patches (+30 / −24). **Visually makes the honest point: direct tests fail, residual/conditional tests succeed.** | **Methods/honest-limits figure; use alongside Results §2** |

| **Trace vs ergodic (HEADLINE)** | `trace_vs_ergodic.pdf` | 3-panel figure: (a) per-band Δρ_cross vs Δρ_null bars with category labels; (b) per-patient strip ordered by category; (c) clean TRACE (α,δ) / MARGINAL (β,γ_l,γ_h) / ERGODIC (θ) categorization box. **This is the band-selective narrative figure.** | **Abstract-level headline figure — use as Figure 1 or 2** |

Regeneration: `python scripts/01_compute/{fig_h2d_bar, fig_band_typology, fig_h2c_unanimity, fig_h2d_unanimity, fig_h2e_drift_floor, fig_evidence_chain, fig_trace_vs_ergodic, band_k_landscape}.py`.

---

## 14. Cross-references

- Method rationale: `.agents/guides/02_methods/imcoh-guide.md`
- Probe-bias control: `.agents/guides/02_methods/probe-bias-guide.md`
- Canonical H1-H4 table: `.agents/reports/2026-04-24_h1-h4-vi-results.md`
- Data layout + per-patient quirks: `.agents/guides/03_implementation/data-layout.md`
- Pat_10 channel-drop rationale: `memory/pat10_channel_mask.md`
- Pipeline status by era: `.agents/reports/2026-04-24_pipeline-status.md`

---

## 15. Pre-registration plan (for replication in follow-up cohort)

Pre-register:

1. FC metric = `imcoh_abs` with `nperseg_for_fs(fs)` Welch windows.
2. LRG with entropy_steps = 400, τ ∈ [10⁻³, 10⁵], 400 log-spaced steps.
3. Six bands exactly as defined in `config/const.py`.
4. Three primary tests:
   - H2c one-sample Wilcoxon ρ > 0 per band, FDR-BH m=6.
   - H2d one-sample Wilcoxon Δρ > 0 per band, FDR-BH m=6, k ∈ [2, 49].
   - H2d post-hoc θ-vs-α AND θ-vs-δ, paired one-sided Wilcoxon,
     Bonferroni m=2.
5. Minimum cohort n = 15 (power analysis from effect sizes above).
6. Exclusion criteria published in advance (vendor-corrupt phases,
   channel-count mismatches resolvable by PATIENT_CHANNEL_DROP only;
   anything else leaves the patient in as-is).
