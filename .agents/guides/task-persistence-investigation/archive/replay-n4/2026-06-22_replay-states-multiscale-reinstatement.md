---
name: replay-states-multiscale-reinstatement
era: IMCOH_ABS_COHORT_N10
status: parked_negative
kind: scope-report
measure: replay-states (time-resolved transient reinstatement of the multiscale-multiband task signature)
headline: N4 (.agents/preprint/headlines/04_replay_states.md)
owner_agent: dynamic-FC / replay agent + null-model chat
created: 2026-06-22
---

# Replay states — transient multiscale-multiband reinstatement of the task configuration in rest

**Head.** A mathematically-scoped, pre-code measure for a **conditional** dynamic
question: during rest, does the network *also* **re-enter** the task connectivity
configuration in brief, recurring **bursts**, identifiable as recurring
**multiscale (ρ^coph) + multiband** reorganization patterns? This is the
*time-resolved* version of the N1 trace (same quantity, phase-averaged). **No code
until this scope is reviewed** (project rule). This file is both the scope report
and the execution+verification brief for the replay agent.

> **⚠ PARKED 2026-06-22 (PI) — read first.** The replay narrative has been **pulled
> from the paper**: removed from the headlines/results structure, and the former N4
> headline is archived (`.agents/preprint/headlines/archive/2026-06/`). The
> bursty/transient-replay test came back **cohort-negative** (audit_124–131), and the
> "sustained reinstatement" framing once carried on N1 was dropped (it was the
> persistence trace relabeled). This scope/brief is **kept as the restore point** —
> revisit only if something new surfaces. Original de-risking note follows.
>
> **⚠ DE-RISKED 2026-06-22 (PI).** The replay/consolidation **hook no
> longer rides this analysis.** It is delivered by **N1 (the β trace is a
> *sustained* reinstatement — proven; see headline N1.7).** This file tests only
> the narrower **bonus** question: does reinstatement *also* occur as **transient
> bursts**? Two consequences for execution: **(1)** test it where short windows are
> feasible — the **gamma** band is the make-or-break case (most cycles per short
> window). **(2)** if gamma shows no bursts, **bursty replay is dropped and folds
> into N1 (sustained)** — that is an acceptable, expected outcome, not a failure.
> The paper does not depend on a positive result here. Current read: increasingly
> unlikely.

---

## 0. Mandatory 5-point critical preamble

1. **Claim.** In `rest_post` (and, as a within-subject control, `rest_pre`) there
   exist **transient, recurring time windows** whose multiscale-multiband
   connectivity configuration matches the task-state configuration (the N1
   signature) — i.e. discrete *replay states* — and they are **more frequent
   and/or more intense in `rest_post`** than in `rest_pre`.
2. **Null.** Window-level "task-likeness" shows **no bursty/transient structure
   beyond a stationary process** and **no `rest_post` > `rest_pre` excess**, under
   (a) a **time-shuffle / phase-randomized** null (destroys temporal structure,
   keeps marginal window statistics) and (b) a **matched-strength per-window**
   surrogate (keeps per-node strength, destroys higher-order structure).
3. **Strongest plausible alternative the null MUST control.**
   - (A) **Smeared stationary mean** — the elevated task-likeness is just the N1
     static shift spread uniformly across rest, with *no* transient states
     (elevated mean, flat dynamics). This is the most dangerous confound.
   - (B) **Arousal / slow drift** — non-specific drift toward *any* non-rest
     configuration (not the task specifically).
   - (C) **Estimator artifact** — |ImCoh| instability at short windows manufactures
     spurious configuration excursions.
4. **Does the null control it — by mechanism.**
   - vs (A): the **time-shuffle null** preserves the per-window task-likeness
     *marginal distribution* but destroys *temporal ordering*; a real replay
     process beats it on a **dynamics statistic** (event rate / dwell time /
     burst clustering), not on the mean. The mean elevation alone is conceded to
     N1 and is *not* the N4 claim. **This is the make-or-break test.**
   - vs (B): require **specificity** — task-likeness must exceed likeness to a
     **random alternative configuration** (e.g. a scrambled or unrelated target
     graph) by the same dynamics statistic; drift toward "anything" fails this.
   - vs (C): a **per-window |ImCoh| stability precheck** (below) is a *gating
     precondition*; replay is only claimed at window lengths where the estimator
     variance is acceptable per band. What it **cannot** fully reject: a subtle
     window-length-dependent bias that co-varies with phase — mitigated by using
     identical windowing across phases and the rest_pre control.
5. **Falsification + residual limitations.** Falsified if, at estimator-stable
   window lengths, (i) `rest_post` task-likeness dynamics do not exceed the
   time-shuffle null, OR (ii) the excess is fully explained by the stationary mean
   (no burst structure), OR (iii) it is not task-specific (fails the
   random-alternative control), OR (iv) `rest_post` ≯ `rest_pre`. Residual limits
   regardless of outcome: n=10, sEEG spatial sparsity, low-band windowability,
   no behavioral validation (performance data unavailable, PI 2026-06-22).

---

## 1. Notation

- Phases `x ∈ {pre, learn, test, post}` with FC graphs from `<|ImCoh|>` per band
  `b`. Canonical task target `T = test` (and optionally `learn`).
- Sliding windows `w = 1..W_x` of length `L` (seconds), hop `h`, within phase `x`.
- Per (x, w, b): `A_{x,w,b}` = |ImCoh| graph on the window → LRG → `D(τ_max)` →
  per-pair cophenetic vector `c_{x,w,b} ∈ R^{P}` (P = N·(N−1)/2 pairs).
- **Phase-reference configurations** (from the *full-phase* graphs, i.e. the N1
  objects): `c̄_{pre,b}` (baseline), `c̄_{T,b}` (task target).
- **Multiband stack:** concatenate / aggregate over a band set `B*` (the
  trace-positive bands, e.g. {α, β}, plus an honest all-band variant).

## 2. Predicates / the measure

- **Per-window task-likeness** (multiscale, one band):
  `s_{x,w,b} = sim( c_{x,w,b} − c̄_{pre,b} ,  c̄_{T,b} − c̄_{pre,b} )`,
  where `sim` = Spearman or cosine over pairs (Spearman to match N1's ρ^coph
  convention). Positive `s` = window leans toward the task configuration relative
  to baseline.
- **Multiband task-likeness:** `S_{x,w} = aggregate_b∈B* ( s_{x,w,b} )` (mean or a
  learned combination; report mean for interpretability). The **multiband**
  requirement = a replay state must show coordinated likeness across `B*`, not one
  band.
- **Replay state (event):** a window (or maximal run) with `S_{x,w}` exceeding a
  threshold `θ` calibrated on the null (see §3). Per phase derive **event rate**
  `R_x = (#replay windows)/W_x` and **dwell time** distribution.
- **Consolidation contrast:** `ΔR = R_post − R_pre` (within subject), cohort-tested.

## 3. Nulls (the gate)

- **N_time (primary):** circularly shift / block-permute the window time-series of
  `S_{x,w}` (or phase-randomize the windowed signals before FC) to destroy temporal
  ordering while preserving the marginal `S` distribution. Replay must beat
  `N_time` on **R_x and dwell-time**, not on `mean S`.
- **N_strength:** matched-strength surrogate per window (preserve node strength,
  scramble structure) → recompute `S`; controls "strength fluctuations look like
  reinstatement."
- **N_altconfig (specificity):** replace `c̄_{T,b}` with a random / unrelated
  target; task-likeness dynamics must exceed this.
- Cohort test: one-sided Wilcoxon on `ΔR` (and on the observed-vs-N_time dynamics
  statistic), LOO-max reported (no single-patient-p-driven), per
  `feedback_no_single_patient_p_driven`.

## 4. Properties / sanity

- At `L → full phase`, `S` collapses to the N1 per-phase ρ^coph (continuity check —
  the measure must reduce to N1 in the limit; bit-level anchor like audit_63).
- `s` is sign-aligned with the N1 trace convention (positive = toward task) —
  reuse the locked sign convention (`feedback_td_sign_convention`).
- Multiband aggregation must not let one band dominate (z-score per band within
  phase before aggregating).

## 5. Caveats (first-class)

- **|ImCoh| windowability** is the binding precondition (§6 step 0). Low bands may
  be unusable at replay-relevant `L`; report the per-band stable-`L` and restrict
  `B*` accordingly. Do **not** claim δ/θ replay at a window where the estimator is
  unstable.
- **Burstiness ≠ mean** — never report `mean S` as the result; the result is the
  dynamics statistic vs `N_time`.
- **Heterogeneous implants** — events are within-subject; cohort claims aggregate
  `ΔR`, not pooled windows.
- **Outliers** — N1 anti patients (Pat_15/10) may have no task signature to
  replay; per-patient reporting mandatory.

## 6. Pseudocode

```
# STEP 0 (PRECONDITION) — windowability
for b in bands:
    for L in candidate_lengths:
        est_var[b,L] = variance of |ImCoh| across repeated/ jackknifed sub-windows
    L*[b] = shortest L with est_var below tolerance
B* = bands with a usable L* (expect alpha/beta; check gammas; flag delta/theta)

# STEP 1 — per-window cophenetic configs
for x in {pre, learn, test, post}:
    for w in windows(x, L*, hop):
        for b in B*:
            A = imcoh_graph(x, w, b)            # reuse load/compute infra
            c[x,w,b] = cophenetic_from_lrg(A)   # reuse N1 ρ^coph pipeline

# STEP 2 — references from full-phase graphs (the N1 objects)
cbar_pre[b], cbar_T[b]  = cophenetic_full_phase(pre,b), cophenetic_full_phase(T,b)

# STEP 3 — task-likeness
s[x,w,b] = spearman(c[x,w,b]-cbar_pre[b], cbar_T[b]-cbar_pre[b])
S[x,w]   = mean_b(zscore_within_phase(s[x,w,b]))

# STEP 4 — events + dynamics
theta    = quantile of S under N_time
R[x], dwell[x] = event_rate_and_dwell(S[x,:] > theta)

# STEP 5 — nulls + cohort test
compare (R, dwell) vs N_time, N_strength, N_altconfig
dR = R[post] - R[pre]; wilcoxon(dR > 0), LOO-max
```

## 7. Visualization

- Task-likeness `S` time-course per phase, replay states shaded; `post` vs `pre`.
- Cross-patient event raster + dwell-time distributions vs `N_time`.
- One replay-window dendrogram snapping to the task config vs a baseline window.
- Ablation: detectability under edge-only / single-band / Grassmann-only vs full
  multiscale-multiband (the CORE-payoff panel). X-epi variant; no C4.

## 8. Connection to prior tools

- **Reduces to N1** (`01_trace.md`, ρ^coph, audit_63) at full-phase `L` — same
  object, time-resolved.
- **Complements the consolidation arc** (`audit_103`, N2): the arc shows *what
  content* persists on average; replay asks *whether it is rehearsed in discrete
  states*. `task_learn` vs `task_test` targets connect the two.
- **Reuses** the sliding-window infra (`compute_time_windows.py`,
  `visualize_time_windows.py`) and the N1 cophenetic pipeline; nulls reuse the
  matched-strength machinery (audit_63 family) + a new time-shuffle null.

## 9. Open questions

- Does replay prefer the **encoding (`task_learn`)** or **inference (`task_test`)**
  target configuration?
- Is per-patient replay-state **rate** related to the N1 **trace magnitude**
  (within-subject, no behavior needed)?
- Cross-band **coordination** of replay states (do bands replay together)?
- Relationship to the epileptogenic network (do replay states avoid/recruit the
  SOZ — ties to N1.6 / N3)?
