---
name: mrl-vs-cbr-reconciliation
type: post-mortem
era: COHORT_N9
status: current
created: 2026-04-25
updated: 2026-04-25
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-25_module-retention-landscape.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cbr-investigation.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cbr-classifier-failure-analysis.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cohesion-cbr.md
  - data/audit/per_patient_hierarchy_cohesion/
  - data/reports/imcoh_mrl/
---

# MRL ↔ CBR reconciliation — what MRL measured, why it under-detected the trace, and where to go next

**MRL was framed as a discrete, hard-thresholded count of "task-anchored
subtrees absent from rest_pre and present in rest_post". The cohort
result is a null: the (J_pre, J_post) cloud lies on the diagonal at
high Jaccard, the trace zone is essentially empty, no band passes
cohort-wide. The same data analysed by Cohesion-CBR (`audit_12`) shows
abundant per-leaf trace classifications visible on the
`per_patient_hierarchy_cohesion` figures. The discrepancy is not in the
data — it is in the operationalisation. MRL filtered out the limbo
zone where most real signals live, used Jaccard which conflates
cohesion with size match, and produced a scalar field instead of a
per-leaf rendering. Cohesion-CBR replaces hard thresholds with soft
2D affinities on `(J_rpre, J_task)` and projects to per-leaf colour.
This sidecar reconciles the two analyses, names the failure modes,
and frames the path forward: the discrete-module question lives in
the soft-affinity / per-leaf rendering family (CBR), not the
threshold-counting family (MRL); for a continuous claim of trace, the
right object is the ultrametric distance matrix or its differential,
not the dendrogram identity.**

## 1. The arc on the MRL side (2026-04-25 evening)

The chronology so the next session does not relitigate it:

| step | what we did | what we found |
|---:|:--|:--|
| 1 | **MRL forward score** — per (patient, band), pick the rest_post subtree maximising `J*(v, T_test) − J*(v, T_pre)`, size ≥ 5 with score ≥ 0.5. Permutation null with 500 size-matched random subsets. | 60 cells, all with positive median trace score (+0.23 to +0.32 per band). Cohort-wide significance: 8–10/10 patients per band at p<0.01 against the random-subset null. **Looked very strong cohort-wide.** |
| 2 | **Per-band dotplot.** | Every band sat above the null cloud. Cohort medians flat across bands (≈ +0.30). Looked great until we noticed: nothing was band-specific. |
| 3 | **Backward control** — same logic but anchored on `T_pre`, score = `J*(pre_subtree, T_learn) − J*(pre_subtree, T_post)`. | Backward distributions almost identical to forward; in 3 bands (δ, β, low_γ) backward was *stronger*. The "trace" was largely temporal-proximity / anatomy-shared structure between adjacent phases, not a task-induced asymmetric memory. |
| 4 | **Forward − backward asymmetry test** (Wilcoxon paired, FDR over 6 bands). | No band passes q<0.05. Best q_BH = 0.22 (high_γ raw p=0.037). α and high_γ have median Δ +0.14 / +0.21 with 7/10 patients positive — directional but underpowered at n=10. |
| 5 | **TAM (Task-Anchored Modules)** — re-anchored: candidates are leafsets in `T_task_learn` that match `T_task_test` at `J_anchor ≥ 0.7`, then test `J_pre` and `J_post` against rest. 1621 TAMs across 60 cells. | Mean `J_pre ≈ 0.65–0.88`, mean `J_post ≈ 0.65–0.88`, mean Δ ≈ 0. The (J_pre, J_post) cloud sits on the diagonal at high Jaccard. The strict trace zone (`J_pre < 0.5 ∧ J_post ≥ 0.7`) is essentially empty. **Discrete-module trace claim collapses.** |

After step 5, the user asked the right question: the parallel
`per_patient_hierarchy_cohesion` audit (Cohesion-CBR) **does** show
visible task-induced traces on actual dendrograms. So why does MRL
return null when Cohesion-CBR does not?

## 2. What MRL was *thought* to measure

The MRL scope (`2026-04-25_module-retention-landscape.md`) declared a
simple discrete observable:

> Per internal node `v ∈ V(T_task_learn) ∪ V(T_task_test)`,
> `retained(v; J_min) ⇔ J*(v, T_rest_pre) < J_min ∧ J*(v, T_rest_post) ≥ J_min`.
> Cohort field `M̄(b, h_rel) = ⟨retained⟩` per (band, log-spaced height bin).

The intuition: count the modules that "appear" during task and survive
into rest_post. The framing is correct in spirit; what's wrong is the
operationalisation.

## 3. What MRL *actually* measures

Three pathologies, in order of severity (all consistent with the
CBR-classifier failure analysis; this is the MRL-specific
restatement):

### 3.1 The Jaccard limbo zone

`retained` is a logical AND of two hard thresholds: `J_pre < J_min`
AND `J_post ≥ J_min`. The empirical distribution of `(J_pre, J_post)`
for task-anchored modules in our data is concentrated in the
**diagonal band at J ∈ (0.5, 0.85)** — neither side of the threshold
fires cleanly. `retained = 1` requires the leafset to be *both*
genuinely absent in rest_pre and genuinely present in rest_post at
*the same* `J_min`. Almost no modules satisfy this; the modal
behaviour is "moderately matched in both", which gives `retained = 0`.

This is the same limbo zone diagnosed in the CBR-classifier failure
analysis (`§3.R1`) — modes 0.5 < J < 0.75 were absorbing 86 of 86
internal nodes in Pat_15 β rest_post. The Jaccard threshold collapses
the rich middle-ground into "not enough" and produces a 100%-zero
result.

### 3.2 Jaccard conflates cohesion with size match

A real trace candidate has two properties:

- *In rest_pre*, the candidate leaves are **scattered**: their
  smallest containing subtree (LCA) is much larger than `|S|`.
  Tightness `T_pre = |S| / |LCA(S, T_pre)|` is small.
- *In task and rest_post*, the same leaves are **bunched**:
  `T_task ≈ T_post ≈ 1`.

Jaccard `J(S, T_pre) = |S ∩ leaves(LCA)| / |S ∪ leaves(LCA)| =
|S| / |LCA|` for the LCA — same numerator, but `J` penalises by
`|LCA| − |S|` extra leaves. For small dispersion (LCA only slightly
bigger than S) `J` is high. For large dispersion (LCA = root) `J → 0`.
The full range `(0, 1]` is informative but the *identification of
"scattered" leaves* requires looking at `T` (tightness), not `J`. Many
real traces sit in the intermediate region where `J` is moderate
(≈ 0.4–0.7) — the limbo zone.

### 3.3 Scalar-field rendering hides the per-leaf identity

MRL outputs `M̄(b, h_rel)` — a number per (band, height bin). The
modules themselves are aggregated away. To show "what trace looks
like", we had to re-pick a strongest candidate and plot dendrograms
manually (the `fig_trace_multipatient_*.pdf` family). This per-cell
selection privileges high-Jaccard candidates and reproduces the
selection bias that made the forward-only score look good cohort-wide
in the first place.

The CBR family by contrast renders the trace **as colour on each
leaf**, where the colour comes from a soft 2D classification on
`(J_rpre, J_task)` aggregated across all eligible ancestors of that
leaf. The figure *is* the proof: the same leaves carry the same colour
across all four phase dendrograms, and the eye reads the trace
directly from how those colours bunch or scatter in each panel.

## 4. The CBR side — what worked

Cohesion-CBR (`audit_12`, scope at `2026-04-25_cohesion-cbr.md`):

- **Candidate space**: every `v ∈ Internal(T_rest_post)` with
  `size ∈ [5, 60]`. The same anchor structure as legacy CBR.
- **Score**: Jaccard on `(J_rpre, J_tl, J_tt)` → derived
  `J_task := mean(J_tl, J_tt)` (the implementation backed off from
  pure tightness because a few stragglers were collapsing tightness
  to ~0; Jaccard with soft affinities turned out to behave well in
  practice — see the `implementation_note` in that scope).
- **Soft affinities**: each anchor `η` produces four affinities in
  [0, 1] computed by corner-distance on the `(J_rpre, J_task)` unit
  square (no thresholds). Trace at `(0, 1)`; persist at `(1, 1)`;
  reset at `(1, 0)`; rearrange at `(0, 0)`.
- **Per-leaf projection**: leaf colour = `argmax_p mean_{η ∋ ℓ} a_p(η)`.
  Saturation = top1 − top2 (confidence margin).

Result: every patient × band figure under
`data/audit/per_patient_hierarchy_cohesion/` shows a per-leaf colouring
that the eye can read directly. Pat_06 α has clear red (trace) bands
that bunch in the task panels and rest_post but are scattered across
the rest_pre dendrogram — the classical trace signature.

## 5. Reconciliation — why the same data gives null on MRL and signal on CBR

| dimension | MRL | Cohesion-CBR |
|:---|:---|:---|
| Anchor | `V(T_task_learn) ∪ V(T_task_test)` (multi-source) | `V(T_rest_post)` (single anchor) |
| Predicates | hard `J_low / J_high` thresholds on `(J_pre, J_post)` | continuous corner-distance affinities on `(J_rpre, J_task)` |
| Limbo zone | filters out everything with `J ∈ (0.5, 0.7)` | every node contributes; no rejection |
| Output | scalar field `M̄(b, h_rel)` + count | per-leaf colour, four classes, soft confidence |
| Visual proof | requires re-selecting modules ad-hoc | figure is the proof |
| Scientific reading | "did any task-anchored module clearly trace?" — answers no at strict cutoffs | "what does each leaf's behaviour pattern look like across phases?" — shows mixed populations |

**Same data, different question.** MRL asks a binary "is there a
discrete trace event?" and the threshold rejects the limbo cases
where most real signals live. Cohesion-CBR asks "for each leaf,
which pattern dominates in its ancestors' joint behaviour?" and lets
the visualisation render the answer.

## 6. What this means for the scientific claim

Two distinct claims, with two distinct proofs:

1. **Discrete cohort-level claim** — "≥ X/10 patients exhibit at
   least one task-induced subtree that persisted into rest_post and
   was absent from rest_pre, at strict thresholds":
   - **MRL**: null result. Diagonal in (J_pre, J_post). Frac_trace
     ≈ 0 across cohort.
   - **TAM forward−backward asymmetry**: directional positive (5/6
     bands) but no band clears FDR. α median Δ +0.14, high_γ +0.21.
   - **Verdict**: cannot be claimed at n=10 with the strict
     hard-threshold operationalisation. The post-mortem
     (`2026-04-24_post-mortem-scalar-session.md`) says exactly this.

2. **Per-patient visual claim** — "this patient's dendrograms show
   leaves that bunch in task and rest_post and scatter in rest_pre":
   - **Cohesion-CBR**: gives a per-leaf classification + figure.
     The eye verifies the trace pattern on actual dendrograms in
     `data/audit/per_patient_hierarchy_cohesion/`. Per-patient
     proofs exist; cohort census is conditional on the classifier.
   - **MRL multi-patient figure** (`trace_multipatient_alpha.pdf`
     etc.): same kind of visual proof, restricted to one
     pre-selected module per patient.

The cohort question and the visual question are not the same
question. MRL conflated them and got the worst of both: a binary
filter too strict to find cohort-wide signal, plus per-cell-best
selection that quietly biased the visual exemplars toward the
diagonal.

## 7. Where the trace genuinely lives in our data

Confirming the prior project finding (`2026-04-24_multiscale-task-trace.md`),
which the MRL detour did not invalidate:

- **H2c continuous Spearman ρ** on ultrametric distance shifts
  `D^post − D^pre` vs. `D^test − D^pre` is the only test that returns
  cohort-wide significance with the right asymmetry sign. 53/54 cells
  positive, all 6 bands q < 0.005 FDR. **Mean ρ ∈ [0.40, 0.54]**.
- **H2e split-half drift floor**: within-session drift gives `ρ ≈ 0`,
  so the H2c signal is not session-order. Already in the writeup.
- **H2d block-pair persistence**: cohort-wide all 6 bands, with
  band-heterogeneity (α strongest, θ weakest). Pair-level continuous
  asymmetry signal.

The trace exists as a **continuous, residual distortion of the
ultrametric distance matrix**, not as the formation/persistence of
discrete modules. This is the load-bearing scientific claim, supported
by every test that operates on the *continuous* distance landscape;
it is *not* supported by the binary-threshold subtree-identity tests.

## 8. Path forward — two parallel directions

### 8.1 If we want the per-patient discrete proof: stay with Cohesion-CBR

The CBR family already provides the right operational shape. The
implementation in `audit_12` and the figures in
`per_patient_hierarchy_cohesion/` are what visually prove the trace
on individual patients. Open work documented in
`2026-04-25_cbr-investigation.md` and the variant scopes; nothing in
this MRL detour changes the CBR roadmap.

The MRL family (this scope and its variants in this folder) is
**superseded by Cohesion-CBR for the discrete claim**. Mark MRL
status `superseded` next session; keep the file with this reconciliation
as the explainer.

### 8.2 If we want the cohort-level visible proof of continuous persistence: pivot to ultrametric matrices

The honest cohort-wide signal is in the **continuous** trace, not the
discrete one. The dendrogram alone cannot show it because the
dendrogram is a discretisation. The natural visual object is the
**ultrametric distance matrix** itself (`D^{p,b,φ}` ∈ ℝ^{N×N}) and how
it transforms across phases.

Concrete proposals (no code yet — scope reports first per the folder
rule):

- **Distance-shift heatmap.** Per (patient, band), compute
  `Δ_task(i,j) = D^test(i,j) − D^pre(i,j)` and
  `Δ_rest(i,j) = D^post(i,j) − D^pre(i,j)`. Show each as an
  N×N matrix. The H2c claim: the *patterns* of Δ_task and Δ_rest are
  positively correlated. A reader sees this as "the same off-diagonal
  cells light up in both heatmaps". This is the continuous analogue
  of "the same leaves bunch together".
- **Co-shift map.** A single matrix `S(i,j) = sign(Δ_task(i,j))
  · sign(Δ_rest(i,j)) · |Δ_task| |Δ_rest|`. Positive cells are pairs
  that moved in the same direction in task and post — the trace.
- **Per-pair scatter.** `Δ_task(i,j)` vs `Δ_rest(i,j)` for all
  upper-triangle pairs. Cohort-wide positive Spearman ρ → cloud
  along positive diagonal.

These figures show the trace **on the matrix that the dendrogram is
built from**, before the discretisation that loses information. They
do not require thresholds; the visual proof is the geometric
correlation.

A scope file `<YYYY-MM-DD>_continuous-trace-matrix.md` should be
written first, per folder rule. After that, two implementation files
(distance-shift heatmaps, co-shift maps) wired into the
`task-persistence-investigation` workflow.

## 9. Files this scope produced (kept for now)

Under `data/reports/imcoh_mrl/`:

- `trace_sweep.csv`, `trace_null.csv`, `trace_backward_sweep.csv`,
  `trace_backward_null.csv`, `trace_asymmetry.csv`,
  `trace_tam_per_module.csv`, `trace_tam_summary.csv`,
  `mrl_per_node.csv`, `mrl_cohort_step.csv`, `mrl_cohort_smooth.csv`,
  `mrl_bin_edges.txt`.

Under `data/reports/imcoh_mrl/figures/`:

- `trace_step1_Pat_06_alpha.pdf` — single-cell sanity (real, but
  with the MRL caveat above).
- `trace_per_band.pdf` — forward-only cohort dotplot vs. random null.
  **Misleading absent the backward control; keep with caveat.**
- `trace_multipatient_{alpha,beta,delta,theta,low_gamma,high_gamma}.pdf` —
  the per-cell-best dendrogram reprojections. Each is a real visual
  pattern for the chosen patients; as cohort proof they are
  insufficient (per the asymmetry test).
- `trace_forward_vs_backward.pdf` — the killing diagnostic.
- `trace_asymmetry.pdf` — paired Wilcoxon, no band clears FDR.
- `trace_anchor.pdf` — TAM cohort + (J_pre, J_post) cloud.
- (legacy) `mrl_landscape*.pdf`, `mrl_per_patient_strips.pdf` — the
  discarded heatmap iteration.

Scripts: `scripts/01_compute/hypothesis_tests/{trace_sweep,trace_null,
trace_backward,trace_anchor,mrl_landscape}.py`,
`scripts/01_compute/figures_embedded/{fig_trace_step1,fig_trace_per_band,
fig_trace_multipatient,fig_trace_forward_vs_backward,fig_trace_asymmetry,
fig_trace_anchor,fig_mrl_landscape}.py`.

## 10. Summary in one paragraph

We treated a continuous-residual phenomenon (the H2c trace, supported
by the ultrametric distance correlation) with a binary
subtree-identity tool (MRL with hard Jaccard thresholds). The tool
filtered out the limbo zone where most real signal lives, conflated
cohesion with size match, and aggregated to a scalar that hid
per-leaf identity. The cohort signal it produced (forward score) is
mostly temporal-proximity / anatomy-shared structure; the
asymmetry-controlled component is small and doesn't survive FDR. The
parallel CBR family, especially Cohesion-CBR with soft affinities and
per-leaf rendering, is the right tool for the **per-patient discrete
visual claim**; the **cohort continuous claim** belongs to the
ultrametric distance matrix family, which is the next direction worth
scoping. MRL is superseded for both claims.
