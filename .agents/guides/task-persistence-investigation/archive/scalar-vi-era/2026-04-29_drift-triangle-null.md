---
name: 2026-04-29_drift-triangle-null
type: scope_report
era: IMCOH_ABS × COHORT_N10
status: scope
created: 2026-04-29
updated: 2026-04-29
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-28_raw-fc-phase-distance.md
  - .agents/reports/2026-04-28_raw-fc-phase-distance-verdict.md
---

# Drift-triangle null — scope

**Renormalization head.** A within-`rest_pre` split-half null only
controls for sampling variance under stationarity. It does **not**
control for monotone session drift across the recording (vigilance,
electrode settling, low-frequency physiological drift). Without a
drift baseline, "persistence" (`T_d < 0` at the FC-edge level) is
indistinguishable from "the brain at session end differs from the
brain at session start regardless of task." The drift-triangle null
builds an empirical drift curve from within-baseline temporal chunks,
extrapolates it to the inter-rest gap, and reports the residual
between observed cross-phase distance and predicted-under-drift.

## Notation

- Phase set `Φ = {rest_pre, task_test, rest_post}`.
- For phase `φ`, time series `X^φ ∈ ℝ^{N × L_φ}` at sampling rate `f_s`.
- Substrate `V*(p, b)` as in
  [`2026-04-28_raw-fc-phase-distance.md`](2026-04-28_raw-fc-phase-distance.md).
- Distances `d_P, d_S, d_F` as in `audit_25_raw_fc_phase_distance.py`.
- `K` = number of consecutive non-overlapping temporal chunks per phase.
- `T_chunk = L_φ / (K · f_s)` = chunk duration (s); require `T_chunk ≥ 30 s`
  for stable `imcoh_abs` estimates.
- Chunk centroid time (relative to phase start): `t_k = (k − 0.5) · T_chunk`,
  `k = 1..K`.
- Chunk FC: `A^φ_k = imcoh_abs(X^φ chunk k, band b)` restricted to `V*`.

## Predicate (per patient × band × distance)

Within-phase pairwise drift distances:
```
D^φ[k, l] = d(A^φ_k, A^φ_l)     k, l ∈ {1..K}, k < l
gap^φ[k, l] = |t_l − t_k|         (seconds, same units across phases)
```
Drift curve: scatter `(gap^φ[k,l], D^φ[k,l])` for `φ ∈ {rest_pre, rest_post}`.

Linear drift fit (per `φ`, per `(b, d)`):
```
D ≈ α^φ + β^φ · gap
```
(Square-root or log alternatives reported as residual diagnostics; default
fit is linear because K=4 gives 6 points which under-supports nonlinear fits.)

Predicted cross-phase distance under drift alone:
```
gap_cross = T_pre/2 + T_inter + T_post/2
d^drift_pred = α + β · gap_cross
```
where `T_inter = L_task_learn/f_s + L_task_test/f_s` (Q5 column
`inter_rest_gap_s`). The fit uses pooled `rest_pre + rest_post` points
under the assumption of a single drift slope across the session;
per-phase fits are reported as a robustness check.

Residual = drift-corrected cross-phase signal:
```
res(φ_A, φ_B) = d_obs(φ_A, φ_B) − d^drift_pred(gap(φ_A, φ_B))
```

A positive residual means the observed cross-phase movement exceeds
what the empirical drift curve predicts at that gap. A residual ≈ 0
means the apparent persistence/movement is consistent with monotone
session drift.

## Properties

1. **Scale-matched.** `D^φ` and `d_obs` use the same distance, the same
   `V*`, and the same FFT-from-segments path → the drift baseline lives
   on the same axis as the observed cross-phase distance.
2. **K underspecification.** `K = 4` gives 6 within-phase pairs; the
   linear fit uses 12 points (pooled across `rest_pre + rest_post`).
   Standard error of the slope is reported per `(band, distance)` and
   the residual confidence interval is propagated.
3. **Non-overlapping segments.** Same convention as the split-half
   null in `audit_25` — preserves chunk independence.
4. **No surrogate.** No phase-shuffled or trial-shuffled null is built
   here; the drift curve IS the null model.

## Verdict logic

For each patient × band × distance:
- `T_d < 0` AND `res(rest_pre, rest_post) ≪ 0`: persistence is robust
  to drift (cross-phase movement smaller than what drift alone would
  produce — the brain at rest_post is *more similar* to rest_pre than
  drift predicts, consistent with anchor / re-equilibration).
- `T_d < 0` AND `res(rest_pre, rest_post) ≈ 0`: apparent persistence
  is drift artefact; downgrade.
- `T_d < 0` AND `res(rest_pre, task_test) ≫ 0` AND `res(rest_pre,
  rest_post) ≪ 0`: the task pushes farther than drift, the post
  returns closer than drift predicts — strongest persistence signature.
- `res(rest_pre, task_test) ≈ 0`: no task-specific signal; cross-phase
  movement is just drift extrapolation.

Cohort threshold: ≥ 5/9 in-pool patients with the strong-signature
pattern in any (band, distance) is the gating criterion for a
drift-controlled persistence verdict.

## Pseudocode

```
for patient in PATIENTS_4PHASE:
    fs = FS_OVERRIDES.get(patient, DEFAULT_SAMPLE_RATE)
    nperseg = nperseg_for_fs(fs)
    for phase in [rest_pre, rest_post]:
        ts = load_timeseries(patient, phase)
        # split into K equal non-overlapping chunks
        chunks = consecutive_chunks(ts, K)
        for k in 1..K:
            fft_k = segment_ffts(chunks[k], fs, nperseg)
            for band:
                A_k_band = imcoh_abs_from_ffts(fft_k, band)
                store A_k_band
    V_star = common_giant_indices(per-phase full-band FCs)  # reuse from audit_25
    for band, distance:
        D_pre = pairwise(A_pre[k] restricted to V*, k=1..K)
        D_post = pairwise(A_post[k] restricted to V*, k=1..K)
        gaps = pairwise centroid gaps
        fit linear D ~ a + b*gap on pooled (D_pre, D_post) ∪ (gaps_pre, gaps_post)
        predict d_drift at gap_cross = T_pre/2 + T_inter + T_post/2
        compute residuals for (pre, task_test) and (pre, rest_post)
        write row to drift_audit.csv
aggregate cohort.
```

## Caveats

- **Linear assumption.** The drift curve may be nonlinear; with K=4,
  6 pairs per phase, only a 2-parameter fit is supportable. Report
  the slope `β` and a goodness-of-fit `R²` per fit; flag low `R²` as
  "drift model unreliable, residual interpretation degraded".
- **Stationarity within chunks.** `imcoh_abs` is a within-chunk Welch
  estimate over `T_chunk ≥ 30 s` — stable in expectation but with
  ~half the segment count of the full-phase observed FC, so
  finite-sample variance inflates `D^φ` slightly. The drift curve is
  therefore biased *upward*; residuals against it are conservative
  for "no persistence beyond drift".
- **Inter-rest gap is the extrapolation distance.** Typical
  `T_inter ≈ 1800 s`, while the within-phase chunk gaps span
  `T_chunk` to `(K-1) · T_chunk` ≈ 100-450 s. The drift fit is
  extrapolated 4-15× past its observed range. Extrapolation error
  dominates the residual error budget for any (band, distance) where
  `R²` is not close to 1.
- **No task-phase chunks.** This audit deliberately fits drift only
  on resting baselines; task-phase FC is excluded because task itself
  is the perturbation under test.

## Connection to prior tools

- Direct extension of `audit_25_raw_fc_phase_distance.py` (same V*,
  same distances, same FFT path).
- Provides the "drift-corrected residual" column requested in
  [`2026-04-28_raw-fc-phase-distance-verdict.md`](../../reports/2026-04-28_raw-fc-phase-distance-verdict.md)
  §"Time-drift confound".

## Open questions

- Does pooling `rest_pre + rest_post` chunks for one slope assume
  drift is monotone across the whole session, or should we fit two
  separate slopes and report both predictions?
- Per-band band-pass filtering within chunks: is the implicit `imcoh_abs`
  band-averaging at the chunk level enough, or do we need to refilter
  the time-series per band before chunking? (Current implementation:
  no refiltering; we band-average inside imcoh_abs from FFTs, same as
  the cached pipeline.)
