---
name: critical-appraisal-of-the-measure
type: report
era: "IMCOH_ABS × COHORT_N10"
status: current
created: 2026-07-06
audience: PI / brutally-honest methodological self-critique of the cophenetic trace ρ_sym
pointers:
  - .agents/reports/2026-07-06_rho-sym-panoramic-and-methodology.md
  - data/audit/rho_sym_gate/
  - data/audit/cross_phase_taxonomy_rhosym/
  - data/audit/per_node_trace_decomposition_rhosym/
---

# Critical appraisal — what ρ_sym actually measures, and where it is weak

## Head (the honest verdict)

The measure is **real, noise-limited, and dual-probe-confirmed — but it is a directional,
magnitude-blind, whole-graph scalar built from a deep transform stack, and the headline
number (ρ≈0.2) reads as "small" only until it is null-referenced (16× a strength-matched
surrogate).** It measures what we want — *does the direction in which the task rearranged the
network's multiscale hierarchy persist into later rest* — but it measures it as a **rank
concordance of a double difference over node-pairs' positions in a global tree**, which is a
specific and lossy operationalization. The biggest honest caveats: (1) ρ=0.2 is a **minority
fraction** of the total cross-phase structure (~a quarter of the "mover" energy; the majority
is rigid anchor); (2) for ~half the cohort the task-induced shift is **at the within-rest
noise floor**, so "no trace" there is partly *unresolvable*, not proven biological null;
(3) the pipeline has many defensible-but-arguable knobs and no behavioral anchor. **The
headline defense is the 16× matched-strength null, not a ceiling ratio** — a crude
noise-ceiling proxy tempted an overclaim (see §3a correction) and does not bear weight. It is
defensible and publishable **if presented as a null-referenced, minority-but-specific effect**
— not as "most of the structure persists," and not as "90% of movers traced" (a global-shift
confound).

---

## 1. Does it measure what we want?

**We want:** does `task_test` leave a band-specific structural trace in `rest_post`?

**It computes:** `ρ_sym = Spearman( D_task − D_preA , D_post − D_preB )` (symmetrized over
A/B), where `D_x` is the per-pair cophenetic (tree-merge-height) vector of the LRG dendrogram
in phase x. So it is the **rank concordance between how the task shifted every pair's position
in the hierarchy and how rest_post shifted it** — a *direction-of-change persistence*.

**Verdict: yes, this is a faithful operationalization of "trace"** — it looks at *changes*
(not the static anchor), so a positive value genuinely means "the task's reorganization
direction survived into rest." But three properties define — and limit — what "trace" then means:

- **It is directional, not quantitative.** Spearman discards magnitude by design (audit_73:
  magnitude weighting reintroduces strength structure and fails matched-strength). So it
  answers "do the same pairs move the same *way*," never "how *much* structure persists." Two
  patients with wildly different amounts of persistent reorganization can score identically.
- **It is a double difference**, `(task−preA)` vs `(post−preB)`. Double-differencing is the
  right way to isolate the task-induced component, but it amplifies noise (see §3) and it is
  what makes the ceiling low.
- **It is symmetric in sign**: negative ρ = reset/rebound (task reorganized, rest_post
  over-corrected). "Trace" and "reset" are the same axis; the sign convention (T>0 = trace)
  is a modeling choice, defensible but not forced.

**One controlled-but-not-eliminated confound:** a monotonic non-task state drift (fatigue,
arousal) that shifts structure steadily through task→post would masquerade as a trace. The
drift-floor null (C2) and matched-strength (C3) are designed against this, and β clears both;
but the measure cannot by itself distinguish "memory trace" from "slow state drift that
happens to be monotonic."

## 2. Nodes or edges? — neither, and this matters

**It is a per-PAIR measure whose per-pair values are GLOBAL tree properties, aggregated to one
whole-graph scalar.**

- The native object is the condensed vector over all `N(N−1)/2` node-pairs.
- But a pair's cophenetic distance is **not** its FC edge — it is the height of the pair's
  lowest common ancestor in the dendrogram, a property of the **entire** hierarchy. Two nodes
  with zero direct connectivity can have a small cophenetic distance if they join early via
  intermediaries. (Confirmed empirically: per-leaf ρ^coph is dominated by global position, not
  local overlap — `perpair_cophenetic_trace_visualization`, `dendrogram_persistence_gate`.)
- So it is **emphatically not an edge measure** (not local FC) and **not natively a node
  measure**. The per-node carrier/anti decomposition (audit_154) is a *re-projection*
  (endpoint incidence), useful but derived — the object lives on pairs-in-a-global-tree.
- **Discreteness caveat:** cophenetic distance is quantized to the `N−1` merge heights, so
  `D_coph` is a coarse ordinal encoding (many ties). The *difference* vectors are far less
  tied (this session verified ~400–560 distinct values, 0% exact zeros — the earlier "98%
  tied" worry was about the raw distances, not the differences), so Spearman is safe; but a
  referee will note the underlying quantization.

**Implication:** the measure's strength (it captures *emergent, multiscale, relational*
structure a raw-FC edge test misses — which is why β localizes to OFC and θ vanishes under
LRG but not raw, audit_146) is the same as its opacity (a single pair's score is not
interpretable locally; only the aggregate is).

## 3. Is ρ≈0.2 justifiable as a "trace"? Does it represent a large fraction of structure?

This is the sharpest question. **Two honest answers that must be stated together:**

### 3a. Against the right yardstick, 0.2 is large — it is near the measurement ceiling.

ρ=0.2 looks "small" only against a mythical ceiling of 1.0. That ceiling is unreachable here:

| quantity (β) | median | reading |
|---|---|---|
| cophenetic-vector reliability (D_preA vs D_preB) | **0.56** | the tree itself is only ~half-reproducible at half-data |
| task-shift SNR (RMS dD_task / RMS dD_base) | **1.16** | the task shift barely exceeds the within-rest noise floor |
| crude achievable ceiling on ρ(dD_task,dD_rest) | **≈0.22** | the *most* ρ could be, given the noise |
| **observed ρ_sym** | **0.20** | see the correction below — do NOT read as "90% of ceiling" |
| matched-strength surrogate median | 0.012 | what strength-matched non-trace gives at this noise |
| **obs / surrogate ratio** | **16.6×** | the null-referenced effect size (THE headline comparison) |

**Correction (2026-07-06, after a deeper check — the crude ceiling does not bear weight, and
"90% traced" is a confound).** The attenuation CONCEPT is valid (a correlation of two noisy
vectors is capped below 1), and the measure is genuinely noise-limited (task-shift SNR ≈ 1.16).
But the crude `1−1/SNR²` ceiling (≈0.22) is too rough to support a "we captured ~90% of the
ceiling" claim, and two findings on re-examination kill the tempting shortcuts:
- **Restricting to movers does NOT raise ρ.** Spearman over the pairs the task actually
  perturbed (`|dD_task| > σ`) ≈ Spearman over all pairs (β cohort median 0.157 vs 0.198). The
  0.2 is not "diluted by irrelevant non-movers" — it is a genuinely modest rank-correlation
  among the movers too. So the ceiling is NOT "about only those that changed"; conditioning on
  changers doesn't move it.
- **"Of those who moved, ~90% traced" is a GLOBAL-SHIFT CONFOUND, not a clean reading.**
  Sign-agreement among movers is 82–99% for strong tracers, BUT that is inflated by a global
  level shift of the whole tree (Pat_03/Pat_06: 99% movers × 99% agreement = near-global
  compression that persists, not pair-specific structural persistence). Spearman ρ deliberately
  strips the global shift (rank-invariant to a constant), which is exactly why it is lower and
  why it is the honest, structure-specific measure.

**So do NOT ceiling-weight and do NOT report "90% traced".** The headline effect size is the
**matched-strength ratio (16×)** — confound-free and already locked; "16× a strength-matched
null" is a stronger and safer sentence than any ceiling ratio. The defensible "fraction that
persisted" is the taxonomy's confound-free energy split (§3b): among movers ~54% of the
perturbation energy persisted (trace) vs ~46% reverted (reset) — a real but modest majority,
NOT 90%. A rigorous ceiling (split-half test-retest of dD_task/dD_rest) remains a worthwhile
SUPPORTING analysis, but it is not the headline and the crude proxy above must not be cited.

### 3b. But it is still a MINORITY of the total cross-phase structure — do not oversell.

The taxonomy (audit_153) gives the direct "fraction of structure" answer, and it is sobering:

| β cross-phase energy | share |
|---|---|
| **anchor** (rigid backbone, unchanged) | **0.63** |
| **trace** (persistence channel) | **0.26** |
| reset (excursion) | 0.11 |

So **~63% of the network's cross-phase structure is rigid scaffolding that never moved, and
only ~26% of the "mover" energy is the persistent trace.** This is not a contradiction with
3a — it is the other half of the truth. The trace is a real, specific, reliably-detectable
*overlay* on a network that is mostly stable. That is exactly what an offline memory trace
*should* look like (you don't rewire the whole brain to consolidate one learned hierarchy),
but it means any sentence implying "the task reorganization broadly persists" is false. The
correct claim is **"a specific, reliably-detectable minority (~a quarter) of the task
reorganization persists, far above chance."**

### 3c. The cohort median hides that where it traces, it traces strongly.

ρ_sym median 0.20 mixes two populations. In the **strong tracers** (Pat_02/03/05/06/08) ρ_sym
is **0.27–0.54** (rank-variance 7–29%; against a null of ~0.01 that is 25–50×) with SNR 2–8×.
In the rest, the task shift is at the noise floor (SNR ~0.8–1.0), ceiling ≈0. So the median
of 0.2 *understates* the effect where it exists. The measure's per-patient spread — which we
treat as signal (`feedback_fluctuations_are_signal`) — means the median is not the right
summary of "how big is the trace when present."

## 4. Is it clearly presentable? Is it "nice"?

**Concept: yes, genuinely elegant.** "Build a principled multiscale hierarchy of the network
in each phase; ask whether the way the task rearranged the hierarchy persists into later rest."
One sentence, no jargon, and it captures a *multiscale, emergent* notion of reorganization that
a flat FC-edge test cannot.

**Implementation: honestly, a deep and knob-laden stack.** FC → graph Laplacian → matrix-
exponential propagator at `τ=1/λmax` → invert to a distance `1/ρ` → UPGMA average-linkage tree
→ cophenetic distances → split-baseline double difference → symmetrized rank correlation →
matched-strength surrogate. Every step is defensible (τ-sweep robustness; UPGMA is standard;
split-baseline decorrelates noise; ρ_sym removes the arbitrary half; matched-strength is the
mandatory null), but collectively it is **six or seven transformations each carrying a
justification**, and a skeptical referee will ask "how much do the specific choices drive the
result?" We have partial answers (τ-sweep robust; dual-probe convergence with Grassmann; raw
baseline comparison) but the pipeline was **not pre-registered**, so implicit selection cannot
be fully excluded by argument alone.

**Presentational risk: the number.** ρ=0.2 will read as "weak" to anyone who skips the
null-referencing and the ceiling. The measure is **only presentable if the null (16×) and the
ceiling (≈0.2 achievable) are shown alongside the point estimate**, and if the taxonomy 26%
is used to state honestly that it is a minority overlay. A bare "r=0.2, p=0.03" invites
rejection.

## 5. What I would actually worry about (steelman the referee)

1. **Small, near-noise-floor effect with many pipeline knobs, no pre-registration.** The
   strongest objection. Defenses: matched-strength 16×, τ-robustness, dual-probe (Grassmann)
   convergence for β, raw-vs-LRG dissociation. But it remains a small effect from a chosen
   pipeline.
2. **No behavioral anchor.** We cannot show the trace predicts memory/inference performance
   (no behavioral data, and none is coming — `task_paradigm_transitive_inference`). So the
   trace is "reliably present," never "behaviorally relevant." This ceilings the claim.
3. **Half the cohort is unresolvable, not proven null.** SNR ~1 for the non-tracers means the
   measure cannot see a trace even if present. "α/β trace in 5–6/10" is honest; "the others
   biologically lack it" partly overreaches — for several it is a detection-limit statement.
   (ρ_sym's "undetermined" tier now says this correctly for the near-zero patients.)
4. **Magnitude-blindness cuts both ways.** It protects against strength confounds but discards
   the very quantity ("how much persisted") a reader most wants.
5. **Cophenetic quantization + global-ness** make individual pairs/nodes uninterpretable; only
   the aggregate scalar and the anatomical re-projection carry meaning.

## 6. What redeems it (why it still stands)

- **Dual-probe convergence.** β clears BOTH the cophenetic trace AND the independent Grassmann
  subspace probe (`per_band_phenomenology_vision`). Two different reductions of the LRG object
  agreeing is the strongest evidence the result is not an artifact of the tree/cophenetic
  choice specifically.
- **It is 16× the mandatory null and noise-limited** (§3a) — the effect is small in absolute
  terms because the *measurement* is noisy (task-shift SNR ≈ 1), not because the *signal* is
  weak where it exists (strong tracers ρ = 0.4–0.55).
- **The LRG layer earns its place**: raw-FC gives a trivial/near-null answer; the multiscale
  transform is what makes β localize to OFC and θ vanish (audit_146) — the structure the
  measure adds is doing real work, not decoration.
- **Every robustness rung passes under the arbitrary-half-free estimator** (this session):
  gate, R=1000 localization + LOO + shaft, inference arc, taxonomy, per-node — nothing rests
  on the coin-flip half anymore.

## Bottom line

ρ_sym measures the right thing (directional persistence of multiscale task reorganization),
measures it near the ceiling the noise allows, and survives every control — **but it is a
lossy, directional, whole-graph scalar reporting a minority overlay (~26%) on a mostly-rigid
network, with a headline number that only reads as strong once null- and ceiling-referenced,
and with ~half the cohort below the detection floor.** Present it as *"a specific,
multiscale, null-referenced minority trace, dual-probe-confirmed in β"* — never
as *"most of the task structure persists."* The measure is honest and defensible at that
claim strength, and overclaims at any stronger one.
