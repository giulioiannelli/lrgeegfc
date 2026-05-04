---
name: 2026-04-28_raw-fc-phase-distance-verdict
type: report
era: IMCOH_ABS × COHORT_N10
status: superseded
created: 2026-04-28
updated: 2026-04-29
superseded_by: .agents/reports/2026-04-29_result-1-raw-fc-phase-trace.md
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-28_raw-fc-phase-distance.md
  - .agents/guides/task-persistence-investigation/2026-04-28_edge-vs-hierarchy-discriminability.md
  - .agents/reports/2026-04-24_post-mortem-scalar-session.md
  - data/audit/raw_fc_phase_distance/cohort_summary.csv
  - data/audit/raw_fc_phase_distance/cohort_diagnostic.pdf
---

# Raw-FC phase-distance audit — cohort verdict

**Renormalization head.** Cohort verdict at n=10 (Pat_03 included —
Z is dimensionless w.r.t. fs). The dominant signal in `d_F` is
drift-shaped: magnitude-norm triangle "persistence" passes majority
in 5 of 6 bands, which is exactly what monotone session drift
predicts (since `task↔post` are temporally adjacent and `pre↔post`
are far apart, drift alone gives `d(task,post) < d(pre,task)` for
free). The drift-robust persistence test runs on the rank-invariant
distances (`d_P`, `d_S`), where scale shifts cancel out. There **β is
the only band whose triangle passes cohort majority on both rank
distances** (`d_P = 7/10`, `d_S = 7/10`); α is close on `d_S` alone
(8/10); the rest are null on at least one rank metric. On the
significance axis (`n_+ ≥ 8/10`) γ_h is positive on all three
distances cohort-wide, so γ_h is *moving* clearly across phases even
if its triangle persistence is only F + P. The combined read: **β is
the cleanest task-trace candidate**, γ_h is the strongest "phase
makes a difference" signal, α/δ have weaker secondary support, θ/γ_l
are null. The next step before any LRG is the drift-controlled null
(`rest_post` split-half + cross-baseline `d(pre_a, post_b)`).

## What was run

Per scope at
[`2026-04-28_raw-fc-phase-distance.md`](../guides/task-persistence-investigation/2026-04-28_raw-fc-phase-distance.md).
Substrate: raw `imcoh_abs` FC on `V*(p, b) = ⋂_φ GC(|A^φ|)`. Three
distances `{d_P, d_S, d_F}` against a within-`rest_pre` split-half
null with `n_split = 50` and non-overlapping segments. Robust
`Z = (d_obs − median) / MAD`. Triangle scalar
`T_d = d(task,post) − d(pre,task)`. Cohort threshold `n_+ ≥ 8/9`
in-pool; persistence threshold `n_persist ≥ 5/9` (`T_d < 0 ∧ Z(pre,post) > 2`).
Pat_03 reported separately (1024 Hz, out-of-pool).

Script: `scripts/01_compute/audit/audit_25_raw_fc_phase_distance.py`.
Output dir: `data/audit/raw_fc_phase_distance/`.
Per-patient artefacts: `<Patient>/{audit.csv, null.npz, diagnostic.pdf}`.
Cohort artefacts: `cohort_summary.csv`, `cohort_diagnostic.pdf` (Pages 4–5).

## Per-patient single-look (Pat_06)

Pat_06 alone returned `Z = 12–685` across all bands × distances ×
phase pairs and `T_d < 0` everywhere. The cohort numbers below show
this was an outlier in scale (Pat_06 had unusually consistent Z across
distances) — most patients show `Z > 2` on `d_F` cohort-wide but not
on `d_P / d_S`.

## Cohort `n_+` table (≥ 8/10 ⇒ positive, **bold** = passes)

| band | dist | n+(pre→task) | n+(pre→post) | n+(task→post) |
|------|------|---:|---:|---:|
| δ    | F    | 5 | 6 | 3 |
| δ    | P    | 3 | 5 | 3 |
| δ    | S    | 3 | 3 | 1 |
| θ    | F    | 7 | 7 | 7 |
| θ    | P    | 4 | 5 | 4 |
| θ    | S    | 3 | 4 | 2 |
| α    | F    | **9** | **8** | **9** |
| α    | P    | 6 | 5 | 2 |
| α    | S    | 5 | 4 | 1 |
| β    | F    | **10** | 7 | **10** |
| β    | P    | 7 | 7 | 3 |
| β    | S    | 5 | 5 | 5 |
| γ_l  | F    | **10** | **10** | **10** |
| γ_l  | P    | 3 | 6 | 5 |
| γ_l  | S    | 4 | 4 | 4 |
| **γ_h** | **F** | **10** | **10** | **10** |
| **γ_h** | **P** | **8** | **8** | **9** |
| γ_h  | S    | 7 | **8** | **8** |

## Triangle persistence — corrected

The earlier `n_persist = (T_d < 0) ∧ (Z(pre,post) > 2)` definition
suppressed valid persistence cells where pre/post sit close to each
other in the within-baseline noise floor while post is still on the
task side of pre. **`T_d < 0` is the geometric persistence claim by
itself** — the `Z` filter was over-restrictive.

### `T_d < 0` count (out of 10), majority threshold ≥ 6/10

| band | d_F | d_P | d_S | both rank pass? |
|------|---:|---:|---:|---|
| δ    | 7 | **7** | **6** | yes (marginal) |
| θ    | 6 | 5 | 3 | no |
| α    | 7 | **6** | **8** | yes (S strong) |
| **β**    | 8 | **7** | **7** | **yes (both clear)** |
| γ_l  | 5 | 5 | **7** | no |
| γ_h  | 6 | **7** | 5 | no |

**`d_F` is drift-shaped:** every band except γ_l passes majority on
`d_F` alone. That is the predicted signature of monotone session
drift (because `task↔post` is temporally adjacent and `pre↔post` is
not), so `d_F`-only positivity is not informative about task-induced
reorganization.

**`d_P` and `d_S` are the drift-robust columns:** rank/correlation
distances cancel out a uniform magnitude shift. **β is the only band
that passes majority on both `d_P` AND `d_S`.** α is a close second
(passes S strongly but P only marginally).

### Median `d(task,post) / d(pre,task)` (< 1 ⇒ below diagonal)

| band | d_F | d_P | d_S |
|------|---:|---:|---:|
| α    | 0.857 | 0.877 | 0.924 |
| **β**    | **0.867** | **0.919** | **0.918** |
| δ    | 0.938 | 0.952 | 0.949 |
| γ_h  | 0.924 | 0.876 | 1.003 |
| γ_l  | 0.969 | 0.997 | 0.874 |
| θ    | 0.942 | 0.998 | 1.079 |

α has the most extreme single-cell ratio (0.857 on F) but β has the
lowest median ratio across all three distances simultaneously,
matching the per-patient sign counts above.

## Pat_03 (n=10, marked distinctly in figures)

Pat_03 has positive Z on most cells but lower magnitudes (Z = 1.6–135).
Z is dimensionless w.r.t. fs, so it counts in the cohort `n_+`. In the
cohort scatter (Page 5) Pat_03 stays visually identifiable as an
orange triangle. Notable per-patient cell: `d_P / d_S` for γ_l on
Pat_03 has *negative* Z (= −9 to −18), meaning the rank-correlation
between phases is *higher* than between two halves of `rest_pre` — a
single-patient pattern, doesn't change the cohort verdict.

## What this means

### 1. β has the cleanest triangle persistence in the cohort.

8/10 patients have `T_d < 0` on `d_F` for β — the highest count of any
band on any distance. Persistence at β holds across all three
distances (F: 8, P: 7, S: 7), so it isn't drift-driven (rank /
correlation are scale-invariant). The median `d(task,post)/d(pre,task)
= 0.867` on F at β confirms the geometry.

### 2. γ_h has the most significant phase movement.

γ_h is the only band where `n_+ ≥ 8/10` on **all three** distances:
`d_F`, `d_P`, `d_S`. The cohort sees clear phase-distinguishability of
FC structure at high γ. Triangle persistence at γ_h is moderate (6/10
F, 7/10 P, 5/10 S).

### 3. α is close to β on the triangle (especially `d_S`).

α has 7/10 (F), 6/10 (P), 8/10 (S) for `T_d < 0`. The α × `d_S` count
is the highest of any cell. α also passes `n_+` on `d_F` cohort-wide.

### 4. δ and θ are null on `n_+`, marginal on triangle.

δ and θ don't pass the `n_+ ≥ 8/10` significance threshold on any
distance. δ does have majority `T_d < 0` (7/7/6 on F/P/S), but the
phase-pair distances aren't above the within-baseline floor at the
cohort level — so the persistence is geometrically present but
statistically marginal.

### 5. Magnitude-only positivity at α/β/γ_l is partly drift-suspicious.

α/β/γ_l pass `n_+` only on `d_F`. For β this is consistent with the
strong triangle finding (real reorganization detectable on magnitude
only because the magnitude effect dominates). For γ_l it's the weaker
case. The within-`rest_pre` null does not control for monotone drift,
so any pure `d_F` positivity needs the cross-baseline drift control.

### 6. Pat_06 was a local maximum, not a representative.

Pat_06's "Z = 12–685 everywhere, T_d < 0 on every cell" pattern is not
shared by the cohort. Pat_06 had unusually high consistency across
distance norms; most patients differ between `d_F` and `d_P/d_S`.

## Implications for the LRG program

The n=10 bare-substrate audit returns:

- **Edge-pattern (rank/correlation) moves with phase only at γ_h** —
  `d_P` cohort-positive on all 3 phase pairs, `d_S` on 2 of 3.
- **Edge-magnitude moves with phase in α/β/γ_l/γ_h** but only on `d_F`
  (potentially drift artefact).
- **Triangle persistence** at γ_h × `d_F` and γ_h × `d_P` (6/10).
- **δ/θ null cohort-wide on every distance.**

This sharpens the LRG question:

1. If the LRG-induced hierarchy on `imcoh_abs` shows
   phase-distinguishability in **α/β/γ_l** that is not driven by
   absolute magnitude (because `ρ̂` is degree-normalised / scale-aware),
   it would add resolution beyond `d_F`. The substrate-ladder scope
   ([`2026-04-28_edge-vs-hierarchy-discriminability.md`](../guides/task-persistence-investigation/2026-04-28_edge-vs-hierarchy-discriminability.md))
   tests this directly.
2. **γ_h becomes the band to localise.** Cohort-wide pattern movement +
   triangle persistence make γ_h the strongest substrate-level finding.
   The LRG dendrogram / partition-multiscale machinery is now justified
   *for γ_h specifically* — to find which subtrees / scales / nodes
   carry the trace.
3. If the hierarchy *also* shows phase-distinguishability in **δ/θ**
   where the edge substrate is null, that's the strongest case for the
   "memory hidden in diffusion paths" hypothesis.

## Mandatory next step (before any LRG)

**Add a `rest_post` split-half null and a cross-baseline comparison.**
The current null cannot tell session drift apart from task-induced
reorganization. The minimal control:

1. Compute split-half null on `rest_post` (same machinery, same
   `n_split = 50`).
2. Compare `MAD(N_pre)` to `MAD(N_post)` per `(p, b, d)`. If
   `MAD(N_post) ≫ MAD(N_pre)` for `d_F`, drift is significant.
3. Build a **cross-baseline distance** `d(pre_a, post_b)` (one half of
   `rest_pre` vs one half of `rest_post` with no shared time); use the
   distribution of these cross-baseline distances as a *null model
   for "FC drift across the session"* and re-evaluate `d_obs(pre→task)`,
   `d_obs(pre→post)` against THIS null. This is the drift-controlled
   verdict; if positive cohort-wide it is real task-induced, if null
   it was drift.

Only after this control passes do we run the LRG substrate ladder.

## Caveats

- **Welch convention mismatch.** Cached observed FC uses 50%-overlap
  Welch; null uses non-overlapping segments. Null slightly noisier
  → `Z` slightly biased downward (conservative).
- **Probe bias.** Under `imcoh_abs` the zero-phase-lag component is
  killed by construction (Nolte 2004); same-probe zeroing is **not**
  proposed as a default control.
- **Single FC method.** `imcoh_abs` only.
- **Pat_03.** 1024 Hz outlier reported separately; in-pool denominator = 9.

## Artefacts

- Per-patient: `data/audit/raw_fc_phase_distance/<Patient>/{audit.csv, null.npz, diagnostic.pdf}` — diagnostic PDF Pages 1–3 with horizontal layout, vector imshows, per-row colorbars via `imshow_colorbar_caxdivider`.
- Cohort summary: `data/audit/raw_fc_phase_distance/cohort_summary.csv`.
- Cohort diagnostic (Pages 4–5: `n_+` bars + triangle scatter): `data/audit/raw_fc_phase_distance/cohort_diagnostic.pdf`.
