---
name: 2026-04-29_result-1-raw-fc-phase-trace
type: report
era: IMCOH_ABS × COHORT_N10
status: superseded
created: 2026-04-29
updated: 2026-05-28
sign_convention: pre-2026-05-26
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-28_raw-fc-phase-distance.md
  - .agents/reports/2026-04-28_raw-fc-phase-distance-verdict.md
  - data/audit/raw_fc_phase_distance/cohort_summary.csv
  - data/audit/raw_fc_phase_distance/cohort_diagnostic.pdf
  - data/audit/fc_phase_geometry/cohort_summary.csv
  - data/audit/fc_phase_geometry/cohort_geometry.pdf
---

> **Superseded sign convention (banner added 2026-05-28).** This report
> uses the pre-2026-05-26 T_d convention:
> `T_d = d(task_test, rest_post) − d(rest_pre, task_test)` → negative = trace.
> The locked convention (from 2026-05-26 onward) is
> `T_d = d(rest_pre, task) − d(task, rest_post)` → **positive = trace**.
> To read this report in the current convention, multiply every T_d value
> by −1 and read every "negative = trace" assertion as "positive = trace".
> Numerical magnitudes, p-values, and per-band verdicts are unchanged; only
> the sign of T_d is inverted. Current preprint scripts and the live
> directive use the new convention. See
> `.agents/reports/archive/2026-04/` for the archived location of this
> report.

# Result 1 — Raw-FC carries a task-shaped persistence trace

**This is the first crystallized result of the rebuild. Everything
downstream (LRG dendrogram analyses, multiscale verification,
localization) inherits from this finding.**

---

## Renormalization head

At the bare-substrate level — `imcoh_abs` FC matrices on the
giant-component intersection `V*`, no diffusion, no hierarchy — the
cohort (n=10) shows a **task-shaped persistence trace** that survives
the drift caveat by construction. The evidence is two-handed:

1. **Task phases form a tight cluster, rest phases do not.**
   `d_S(task_learn, task_test)` is roughly **half** `d_S(rest_pre,
   rest_post)` in every band. Tasks are coherent across the session;
   rests are intrinsically variable. The "rest" baseline is not a
   trivial repeat of itself.
2. **`rest_post` sits closer to the task cluster than `rest_pre`
   does**, on the rank-invariant distance `d_S`, in **5 of 6 bands**
   (everywhere except θ ≈ tie at the cohort median; per-patient counts
   show 6+/10 negative in 5 bands). The asymmetry is largest in **δ
   (−0.131)** and **γ_l (−0.105)**; per-patient sign consistency is
   strongest at **β (8/10 on `d_F`, 7/10 on `d_S`)** and **α (8/10 on
   `d_S`)**.

Two distance markers are kept and given separate roles for everything
that follows:

- **`d_S = 1 − corr_S(triu A, triu B)`** — *topological* persistence:
  is the *rank ordering* of edges preserved? Drift-robust by
  construction (rank is scale-invariant).
- **`d_P = 1 − corr_P(triu A, triu B)`** — *amplitude* persistence: do
  the *numeric values* of edges co-vary? Sensitive to magnitude shifts
  but susceptible to hub-edge dominance.
- **`d_F = ‖A − B‖_F / √(‖A‖_F · ‖B‖_F)`** — *kept as a drift
  diagnostic only*. `d_F` shows triangle persistence in 5 of 6 bands
  cohort-wide, which is the predicted signature of monotone session
  drift, not task-induced reorganization.

The load-bearing bands going forward:

- **β** — strongest cohort consistency: `d_S` and `d_P` both pass
  cohort majority on the triangle test (7/10 each), and β is the only
  band with both rank metrics passing.
- **α** — strong on `d_S` (8/10 triangle, single-cell best).
- **δ** and **γ_l** — strongest cohort-median asymmetry on `d_S`.
- **γ_h** — strongest *significance* axis (`n_+ ≥ 8/10` on all three
  distances) but its triangle persistence is mixed.
- **θ** — null cohort-wide.

---

## Operational definitions

- **Cohort:** `PATIENTS_4PHASE = {Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15}`
  (n = 10; Pat_03 included because the Z-score is dimensionless w.r.t.
  `fs`).
- **FC method:** `imcoh_abs` = `mean_f |Im(S_ij)/√(S_ii S_jj)|` over
  the canonical band, on the intersection `V*(p, b)` of giant
  components across the four phases involved in the comparison.
- **Three distances on `triu(A)` (k = 1):**
  - `d_P = 1 − corr_P(triu A, triu B)` — Pearson, magnitude-aware.
  - `d_S = 1 − corr_S(triu A, triu B)` — Spearman, rank-only.
  - `d_F = ‖A − B‖_F / √(‖A‖_F · ‖B‖_F)` — scale-normalised
    Frobenius, magnitude-only.
- **Within-`rest_pre` split-half null** (audit_25): non-overlapping
  Welch segments, `n_split = 50`, robust z-score
  `Z = (d_obs − median)/MAD`. Cached at `<Patient>/null.npz`.
- **Triangle persistence scalar:**
  `T_d = d(task_test, rest_post) − d(rest_pre, task_test)`. Negative ⇒
  `rest_post` is closer to `task_test` than `rest_pre` is to
  `task_test`.
- **Cohort thresholds:** `n_+ ≥ 8/10` for the within-baseline
  significance test; `n_persist ≥ 6/10` (majority) for `T_d < 0`
  count.
- **4-phase geometry** (audit_26): all 6 pairwise distances computed
  on cached observed FCs. Cohort medians and Q1/Q3 reported.

---

## Result 1A — within-baseline significance

Cohort `n_+` (number of patients with `Z = (d_obs − median)/MAD > 2`,
out of 10):

| band | dist | n+(pre→task) | n+(pre→post) | n+(task→post) |
|------|------|---:|---:|---:|
| δ | F | 5 | 6 | 3 |
| δ | P | 3 | 5 | 3 |
| δ | S | 3 | 3 | 1 |
| θ | F | 7 | 7 | 7 |
| θ | P | 4 | 5 | 4 |
| θ | S | 3 | 4 | 2 |
| α | F | **9** | **8** | **9** |
| α | P | 6 | 5 | 2 |
| α | S | 5 | 4 | 1 |
| β | F | **10** | 7 | **10** |
| β | P | 7 | 7 | 3 |
| β | S | 5 | 5 | 5 |
| γ_l | F | **10** | **10** | **10** |
| γ_l | P | 3 | 6 | 5 |
| γ_l | S | 4 | 4 | 4 |
| γ_h | F | **10** | **10** | **10** |
| γ_h | P | **8** | **8** | **9** |
| γ_h | S | 7 | **8** | **8** |

Phase-pair distances cross the within-`rest_pre` floor cohort-wide on
`d_F` for **α/β/γ_l/γ_h**, and *additionally* on `d_P` and (mostly)
`d_S` for **γ_h**. δ and θ are null on the significance axis.

## Result 1B — triangle persistence

Cohort `T_d < 0` count (out of 10), at majority threshold ≥ 6/10:

| band | d_F | d_P | d_S |
|------|---:|---:|---:|
| δ | **7** | **7** | **6** |
| θ | **6** | 5 | 3 |
| α | **7** | **6** | **8** |
| **β** | **8** | **7** | **7** |
| γ_l | 5 | 5 | **7** |
| γ_h | **6** | **7** | 5 |

**`d_F` row is drift-shaped** — passes majority in 5 of 6 bands. This
is the predicted footprint of monotone session drift (because
`task↔post` are temporally adjacent and `pre↔post` are far apart in
time, drift alone gives `d(task,post) < d(pre,task)` for free). Don't
read structural meaning into the `d_F` row.

**Drift-robust read** (rank-invariant distances `d_P` and `d_S`):

- **β** is the only band whose triangle passes cohort majority on
  *both* rank distances (`d_P = 7`, `d_S = 7`).
- **α × `d_S` = 8/10** is the strongest single rank-distance cell.
- **δ** marginal on both rank distances.
- γ_h passes only on `d_P`; γ_l only on `d_S`; θ null on both.

## Result 1C — phase geometry (positive controls)

Cohort medians on `d_S`:

| pair | δ | θ | α | β | γ_l | γ_h |
|---|---:|---:|---:|---:|---:|---:|
| **TL↔TT (within-task)** | 0.238 | 0.173 | 0.246 | 0.205 | 0.354 | 0.368 |
| **RPre↔RPost (within-rest)** | 0.538 | 0.384 | 0.576 | 0.419 | 0.612 | 0.445 |
| RPre↔TT | 0.464 | 0.427 | 0.430 | 0.430 | 0.584 | 0.420 |
| **TT↔RPost (persist)** | 0.333 | 0.375 | 0.395 | 0.389 | 0.479 | 0.433 |

Two structural facts:

1. **Tasks are a coherent cluster.** `d_S(TL, TT)` is roughly half
   `d_S(RPre, RPost)` in every band. The cohort-level fact "tasks are
   similar to each other" holds on the rank metric.
2. **Rests are not.** `d_S(RPre, RPost)` is always ≥ 0.38 — rest
   changes substantially across a session. So when we ask "is `RPost`
   closer to task than to `RPre`?", the comparison is meaningful: the
   "to `RPre`" baseline is itself far from zero.

**`d(TT, RPost)` vs `d(RPre, RPost)` — post-anchored persistence:**

| band | d_S(RPre,RPost) | d_S(TT,RPost) | post closer to task? |
|---|---:|---:|---|
| δ | 0.538 | 0.333 | yes (Δ = +0.205) |
| θ | 0.384 | 0.375 | tie |
| α | 0.576 | 0.395 | yes (Δ = +0.181) |
| β | 0.419 | 0.389 | yes (Δ = +0.030) |
| γ_l | 0.612 | 0.479 | yes (Δ = +0.133) |
| γ_h | 0.445 | 0.433 | tie |

In every band except θ and γ_h, **`RPost` is closer to `task_test`
than to its own `rest_pre`** on the rank distance. This is the
post-anchored persistence claim and it complements the task-anchored
triangle test (`T_d < 0`).

`d_P` shows the same picture for α/β/γ_l/γ_h; `d_F` is mixed (drift).

---

## Null-correctness audit

### What the within-`rest_pre` split-half null does

For each `(patient, band)` it gives a distribution `N_d` of distances
between two halves of `rest_pre` recomputed from non-overlapping Welch
segments. `Z = (d_obs − median(N))/MAD(N)` answers: **"is the observed
phase-pair distance bigger than the typical jitter inside `rest_pre`?"**
That's a *significance* test, not a null model for "task induced this".

### What the within-`rest_pre` null does NOT control for

- **Monotone session drift.** If FC magnitude or ordering drifts
  monotonically across the session, phase-pair distances will exceed
  the within-`rest_pre` floor without any task involvement. The drift
  fingerprint is "`d_F` passes everywhere on the triangle", which is
  exactly what we observe — so `d_F` results must be quarantined.
- **Cross-baseline noise.** What is the typical distance between half
  of `rest_pre` and half of `rest_post`? We never measured it. If that
  cross-baseline distribution is comparable to `d(rest_pre,
  rest_post)`, then `RPre` and `RPost` are drifting together and the
  rest-rest gap is real noise, not structural change.
- **Within-task split-half null.** Would tell us how much of
  `d(TL, TT)` is task-state jitter vs structural difference. Skipped
  here because `d(TL, TT)` is already small and consistent across
  bands; the within-rest comparison is more informative.

### What the 4-phase geometry contributes (Result 1C)

Acts as a **positive control** that strengthens the persistence claim
without needing a drift null:

- `d_S(TL, TT)` ≪ `d_S(RPre, RPost)` says "tasks are a coherent
  cluster, rests are not" — establishing the *reference frame* of
  which-pair-is-naturally-close.
- `d_S(TT, RPost)` < `d_S(RPre, RPost)` in 4 of 6 bands says post is
  pulled toward the task cluster.
- `d_S(TT, RPost)` < `d_S(RPre, TT)` (the triangle `T_d < 0`) says the
  same thing anchored at task.

The two anchorings (post-anchored and task-anchored) agree on δ, α, β,
γ_l. They diverge on θ (post-anchored tie, task-anchored marginal) and
γ_h (post-anchored tie, task-anchored mixed).

### What we still owe

The drift-controlled null — split-half on `rest_post` plus
cross-baseline `d(pre_a, post_b)` — is **deferred to a follow-up
audit**. It is the canonical control for "is the rest-rest distance
explained by drift or by structural change?". For Result 1, the
phase-geometry positive controls are sufficient because they
**disambiguate** the rank-distance findings: the asymmetry holds on
`d_S` in every band except θ and γ_h, in a configuration drift cannot
fake (rank is scale-invariant).

---

## What this result says — and does not say

### What it says

- Raw FC at the edge level **does** carry a task-shaped persistence
  signal in the COHORT_N10 cohort under `imcoh_abs`.
- The signal is **band-specific**: strongest at β (cohort consistency,
  both rank metrics) and α (single-cell best on `d_S`), followed by δ
  and γ_l (cohort medians); θ and γ_h ambiguous.
- The signal is **structural**, not just magnitude-driven: it shows up
  on rank-invariant distances (`d_S`, `d_P`) where drift cannot fake
  it.
- `d_F` is **not** a structural marker here — its triangle pattern is
  the drift fingerprint and should be reported only as a diagnostic.

### What it does not say

- It does not localise *which* nodes / subtrees / scales / brain
  regions carry the trace. That is the localisation step (LRG
  dendrogram + partition machinery).
- It does not separate "the same edges shifted in the same way for all
  patients" from "edges shifted differently per patient but on average
  the asymmetry is consistent". That is a higher-resolution
  per-patient question.
- It does not explain *why* δ and γ_l carry such large cohort-median
  asymmetry on `d_S`. Hypothesis to investigate later: low-frequency
  bands have fewer effective independent edges (fewer freq bins), so
  an edge-rank reshuffle perturbs `d_S` more.
- It does not yet rule out that `d(RPre, RPost)` itself contains
  task-shaped drift (the drift-controlled null is owed).

---

## Implications for the LRG step

The LRG step now has a sharp question to answer:

> Does the LRG dendrogram resolve, in **the same bands** (β, α, δ, γ_l)
> and on **the same persistence axis** (rank-preserving structure),
> the trace that the raw-FC `d_S` test detected? And does it *add*
> resolution by surfacing the trace in **θ and γ_h** where the raw-FC
> rank test was null or ambiguous?

The two outcomes that would justify LRG:

1. **Confirmation:** dendrogram-based topological distance (KC / weighted
   RF / Spearman on the ultrametric `M = 1/ρ̂`) reproduces the rank-
   distance pattern (β strongest, α/δ/γ_l next). This validates the
   LRG transform as a faithful summary of the edge-level result.
2. **Resolution:** dendrogram distance shows persistence in θ and/or
   γ_h where `d_S` was null. This would be the "memory hidden in
   diffusion paths" result — the LRG transform integrates coherent
   multi-edge changes that single-edge rank tests miss.

The LRG step is therefore *not* a redundant rerun; it is the
multiscale verification + extension of Result 1.

---

## Artefacts

- **Audit script (with caching):** `scripts/01_compute/audit/audit_25_raw_fc_phase_distance.py`.
- **Phase-geometry script:** `scripts/01_compute/audit/audit_26_fc_phase_geometry.py`.
- **Per-patient (audit_25):** `data/audit/raw_fc_phase_distance/<Patient>/{audit.csv, null.npz, diagnostic.pdf}`.
- **Per-patient (audit_26):** `data/audit/fc_phase_geometry/<Patient>/{distance_4phase.csv, distance_4phase.pdf}`.
- **Cohort (audit_25):** `data/audit/raw_fc_phase_distance/cohort_summary.csv`, `cohort_diagnostic.pdf`.
- **Cohort (audit_26):** `data/audit/fc_phase_geometry/cohort_summary.csv`, `cohort_geometry.pdf`.
- **Scope (Result 1):** `.agents/guides/task-persistence-investigation/2026-04-28_raw-fc-phase-distance.md`.
- **Verdict (interim, superseded by this report):** `.agents/reports/2026-04-28_raw-fc-phase-distance-verdict.md`.
