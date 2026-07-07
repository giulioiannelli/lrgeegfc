---
name: band-agnostic-lrg
type: scope
era: COHORT_N10 / IMCOH_ABS
status: current
created: 2026-05-08
updated: 2026-05-08
scope: band-agnostic-lrg
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-26_continuous-trace-matrix.md
  - .agents/references/notes_imcoh_260508.pdf
  - data/audit/band_agnostic_lrg/cohort_summary.csv
  - data/reports/notes_verification_2026-05-08/figures/band_agnostic_lrg_teaser.pdf
  - scripts/01_compute/audit/audit_58_band_agnostic_lrg.py
---

# Band-agnostic LRG — does the cohort task-trace direction survive without band splitting?

**Head.** A band-stacked proxy of the controlled CTM ρ_split (concatenate
the per-band Δ_task / Δ_rest pair-vectors across all six canonical EEG
bands per patient, then a single Spearman) yields a broadband cohort
direction count of **9 / 10 patients positive**, all **10 / 10 above the
within-session drift floor** (paired Wilcoxon p = 0.0025), with median
ρ^pool = +0.148 — comparable in *direction* to the manuscript's
per-band trace bands (α 8/10, β 8/10, γ_l 7/10) but smaller in magnitude
than the strongest single band per patient. **Verdict (proxy-tier): the
trace direction is band-agnostic in the sense that no band-specific gate
is needed to recover a cohort-positive cohort signature; band splitting
is a magnitude refinement, not a direction requirement.** This is a
*proxy* result (band-pool, not spectral broadband) — the canonical
broadband |ImCoh| → LRG pipeline (§11) is fully specified and parked at
≤ 30 minutes runtime.

## 1. Notation

Indices and sets:

| Symbol | Domain | Meaning |
|:---|:---|:---|
| `p ∈ P` | `|P| = 10` (`Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15`) | patient |
| `b ∈ B` | `B = {δ, θ, α, β, γ_l, γ_h}`, `|B| = 6` | canonical EEG band |
| `f` | `f ∈ [0, f_Nyquist]` Hz | frequency bin |
| `(i, j)` | `i < j`, `(i, j) ∈ L_p × L_p` | upper-triangle contact pair |
| `m_p` | `m_p = N_p (N_p − 1) / 2` | upper-triangle pair count for patient `p` |
| `m_p^pool` | `m_p^pool = |B| · m_p` | concatenated pair count across the 6 bands |

Per `(p, b, φ)` for `φ ∈ {pre_A, pre_B, test, post}`:

- `D^{p,b,φ}(τ_max) ∈ ℝ_{≥0}^{N_p × N_p}` — LRG ultrametric distance matrix
  at `τ = 1/λ_max`. `pre_A` and `pre_B` are computed on disjoint halves of
  the rest_pre time series (`data/cache/imcoh_lrg_halves/`).

Existing per-pair shift vectors (from
`continuous_trace_matrix.py --mode split-baseline`):

```
Δ_task^{p,b}(i, j) := D^{p,b,test}(i, j) − D^{p,b,pre_A}(i, j)         ∈ ℝ
Δ_rest^{p,b}(i, j) := D^{p,b,post}(i, j) − D^{p,b,pre_B}(i, j)         ∈ ℝ
```

Independence of half-baselines `D_pre_A ⊥ D_pre_B` removes shared-baseline
noise correlation. Stored as `dD_task` / `dD_rest` in
`data/reports/imcoh_continuous_trace/per_pair_split/{pat}_{band}.npz`.

## 2. Definitions

### 2.1 Per-band ρ_split (existing, §5.3 of `notes_imcoh.tex`)

```
ρ_split^{p,b} := Spearman( vec_upper(Δ_task^{p,b}), vec_upper(Δ_rest^{p,b}) )   ∈ [−1, +1]
```

`ρ > 0` = trace direction at band `b`. Cohort gate per band:
`#{p : ρ_split^{p,b} > 0} ≥ 8/10` and paired Wilcoxon
`ρ_split^{p,b} > ρ_drift^{p,b}` at one-sided p ≤ 0.05.

### 2.2 Band-pooled broadband proxy ρ_split^pool

For each patient `p`, concatenate the per-band shift vectors over the
six bands and compute a single Spearman:

```
Δ_task^{p, pool} := concat_b vec_upper(Δ_task^{p,b})         ∈ ℝ^{m_p^pool}
Δ_rest^{p, pool} := concat_b vec_upper(Δ_rest^{p,b})         ∈ ℝ^{m_p^pool}
ρ_split^{p, pool} := Spearman(Δ_task^{p, pool}, Δ_rest^{p, pool})       ∈ [−1, +1]
```

`m_p^pool = 6 · m_p` (≈ 38k–44k pairs per patient). The vectors
inherit the half-baseline-split independence per-band. Pooling across
bands gives equal weight to each band at the rank-aggregation level
(the rank distribution is per-band-stratified, so the cohort's
ρ_split^pool is dominated by the bands with the largest absolute
shifts at each patient).

### 2.3 Band-mean ρ_split (auxiliary)

```
ρ_split^{p, band-mean} := (1/|B|) · Σ_b ρ_split^{p,b}            ∈ [−1, +1]
```

Cohort proxy 2 — equal-weight scalar mean of the per-band ρ_split. Less
defensible than ρ^pool because it averages correlation coefficients
(non-linear), but transparent and reproducible from the cached CSV.

### 2.4 Drift floor (band-mean)

```
ρ_drift^{p, band-mean} := (1/|B|) · Σ_b ρ_null_drift^{p,b}
```

`ρ_null_drift^{p,b}` is the existing within-session split-half drift
(Run C) from `data/reports/imcoh_continuous_trace/controls_summary.csv`:
`ρ_null_drift = Spearman(D^pre_B − D^pre_A, D^post_B − D^post_A)`. Same
halved-data noise regime as the Run-A ρ_split. The "broadband drift" is
its band-mean.

### 2.5 Cohort statistics

```
n_pool_pos        := #{ p ∈ P : ρ_split^{p, pool} > 0 }
n_pool_above_drift:= #{ p ∈ P : ρ_split^{p, pool} > ρ_drift^{p, band-mean} }
z_pool, p_pool    := wilcoxon_z( {ρ_split^{p, pool} − ρ_drift^{p, band-mean} : p ∈ P} )
```

Trace gate (one-sided): `n_pool_above_drift ≥ 8/10` AND `p_pool ≤ 0.05`.

## 3. Properties

| # | Statement |
|:--|:----------|
| 3.1 | `ρ_split^pool ∈ [−1, +1]`. Direction-only by construction; magnitude is band-rank-pool dependent. |
| 3.2 | **Pooling preserves direction; flattens magnitude.** A patient with strong α-only trace and δ noise will have a smaller ρ^pool than ρ_α — but the *sign* is dominated by the strongest band. The cohort direction count `n_pool_pos` is therefore conservative wrt the per-band counts (a patient at α 8/10 might pool to 7/10 or 9/10; the latter is the §10b finding). |
| 3.3 | **Independence inherited.** The half-baseline split (`D_pre_A ⊥ D_pre_B`) is per-band; the concatenated vectors retain it because each band slice was computed independently. The pooled ρ has the same shared-baseline-noise immunity as the per-band ρ. |
| 3.4 | **Pat_03 1024 Hz**: γ_h band has half the spectral coverage (≈ half the freq bins folded into the per-band cache); its Δ vectors are smaller in number per leaf-pair and noisier. In the band-pool, this manifests as a higher per-band sampling variance contribution from γ_h *only at Pat_03*. The pool ρ ranks pairs across bands → Pat_03's ρ^pool is not directly comparable to the others; flag (table §10). |
| 3.5 | **Not a spectral broadband.** The pool concatenates *band-windowed* pair-vectors, not pair-vectors derived from a single full-spectrum |ImCoh|. Two notable consequences: (a) a pair (i, j) appears 6 times in the pooled vector, once per band; (b) the LRG step is *per-band*, not on a unified broadband adjacency. The canonical broadband pipeline (§11) lives one level deeper. |
| 3.6 | **No new helpers introduced.** Reuses `load(...).npz`, `scipy.stats.spearmanr`, `scipy.stats.wilcoxon` via `wilcoxon_z`. No private re-implementations. |
| 3.7 | **Cohort direction count is the contract**, not the magnitude. ρ^pool magnitude is uninterpretable as a "broadband effect size" because the canonical broadband pipeline would weight bins by their spectral count (1/f weighted in practice — δ over-represented). The pool gives equal band weight; the canonical pipeline would give equal *bin* weight. The two are not the same number even when both yield the same direction. |

## 4. Caveats & failure modes

### 4.1 Band-pool ≠ spectral broadband

The pooled vector treats each band's pair-shifts as 6 independent
"observations" of the same pair (i, j). A canonical spectral broadband
|ImCoh| would average over all freq bins in `[0, f_Nyquist]` *first*,
producing one matrix per phase; the LRG and its half-split would be on
that one matrix. Direction agreement between the two protocols
(broadband-pool vs broadband-spectral) is plausible but not proven by
this script. **Mitigation:** §11 specifies the canonical pipeline. Run
it as a confirmation if the proxy-tier signal is positive. (It is.)

### 4.2 1/f weighting absent in the pool

A canonical broadband would amplify δ and θ via the 1/f spectrum (most
power in the lowest bins). The pool gives δ, θ, α, β, γ_l, γ_h equal
representation. Pool ρ_split magnitude therefore *under-represents* the
δ/θ contribution relative to a spectral broadband. **Disclose:** the
pool result reports band-equal-weighted direction, not power-weighted
direction.

### 4.3 Jensen's inequality on |ImCoh|

`mean(|signed|) ≠ |mean(signed)|` and the canonical `imcoh_abs`
applies `np.abs` *per frequency bin first*, then averages. The pool
inherits this correctly because each per-band cache already applied
the per-bin abs (see `workflow/fc.py:_load_imcoh`). The canonical
broadband pipeline must do the same.

### 4.4 Pat_03 1024 Hz coverage half

Pat_03's per-band freq-resolved caches are valid (verified in the
`config.const.FS_OVERRIDES` and `nperseg_for_fs` logic), but its γ_h
band has half the bins of the cohort. In the per-band pipeline this is
absorbed by the band average; in the pool the band still appears with
the same number of pairs (m_p) but each pair is a noisier estimate.
Pat_03's ρ^pool is therefore directionally comparable but the
magnitude is moderately deflated. **Mitigation:** flag in the per-patient
table; cohort statistic robust by virtue of n=10 and the rank-based
Wilcoxon. (Pat_03 ρ^pool = +0.415 — direction holds.)

### 4.5 Pair-pool independence violation

The 6 · m_p pooled pairs are not independent — pairs share leaves
(within a band) and appear multiple times across bands. The Spearman ρ
is fine as a *descriptive* effect size; inference at the pooled level
would over-state precision. The cohort-level Wilcoxon is on per-patient
ρ^pool values (n=10), which is the correct inference unit and inherits
the cohort independence of the patients.

### 4.6 What the proxy cannot decide

- *Whether band specificity is meaningful for downstream interpretation.*
  Pool ρ^pool > 0 says "direction is there"; it does not say "the α/β/γ_l
  per-band gate is redundant". A patient with ρ^pool = +0.05 but ρ_β =
  +0.45 still benefits from band splitting (β localizes the trace; pool
  diffuses it).
- *Whether canonical broadband (spectral, not pool) recovers the same
  direction count.* §11 must be run to settle this.
- *Whether the trace, anchor, reset, rearrange taxonomy
  (`terminology.md`) survives broadband.* Direction-only proxy.

## 5. Pseudocode

```
Input:
  P = {Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15}
  B = {delta, theta, alpha, beta, low_gamma, high_gamma}
  CTM_PAIR  = data/reports/imcoh_continuous_trace/per_pair_split
  CONTROLS  = data/reports/imcoh_continuous_trace/controls_summary.csv

For each patient p ∈ P:
  parts_task ← []
  parts_rest ← []
  For each band b ∈ B (in order):
    f ← CTM_PAIR / "{p}_{b}.npz"
    require f.exists()                            # 60 cells available
    npz ← np.load(f)
    parts_task.append(npz["dD_task"])             # length m_p
    parts_rest.append(npz["dD_rest"])
  dD_task_pool ← np.concatenate(parts_task)        # length 6 · m_p
  dD_rest_pool ← np.concatenate(parts_rest)
  rho_pool ← scipy.stats.spearmanr(dD_task_pool, dD_rest_pool).statistic
  rho_band_mean ← (1/6) · Σ_b ρ_split^{p,b}        # from controls CSV
  rho_drift_mean ← (1/6) · Σ_b ρ_null_drift^{p,b}  # from controls CSV
  emit row(p, rho_pool, rho_band_mean, rho_drift_mean, ...)

# Cohort
n_pos        ← #{rho_pool > 0}
n_above_drift← #{rho_pool > rho_drift_mean}
z, p_value   ← wilcoxon_z( rho_pool - rho_drift_mean )    # one-sided greater
emit cohort_summary.csv
```

### 5.1 Complexity

Per patient: `O(6 · m_p log m_p)` for the pooled Spearman; trivial
(< 1 s on 40k pairs). Cohort total: ≤ 10 s including I/O. Whole
script wall-time: ≤ 15 s.

## 6. Visualization spec

Two-panel teaser figure
(`data/reports/notes_verification_2026-05-08/figures/band_agnostic_lrg_teaser.pdf`):

### Panel A — Per-patient pairing of drift vs broadband ρ_split

- Two x-positions: `drift` (band-mean ρ_null_drift) and `broadband`
  (ρ_split^pool).
- 10 patients × 2 dots, paired with thin gray lines.
- Dotted reference line at `ρ = 0`.
- Title: cohort Wilcoxon statistic and `n_above_drift / n` count.
- **Reading rule:** if every gray line slopes upward from `drift` to
  `broadband`, the trace direction is unanimous and survives the drift
  floor. Mixed slopes = some patients without trace direction.

### Panel B — Best trace-band ρ_split vs broadband ρ_split

- x-axis: per-patient max of `ρ_split^{p, b}` over `b ∈ {α, β, γ_l}`
  (the manuscript's three trace bands).
- y-axis: `ρ_split^{p, pool}` (broadband).
- Per-patient dot, colored by which of α / β / γ_l was the best for
  that patient (green α, orange β, purple γ_l).
- Diagonal `y = x` reference line. Patient ID labels next to dots.
- **Reading rule:** dots near `y = x` mean broadband recovers the
  full magnitude of the per-patient strongest single band; dots
  significantly *below* `y = x` mean band splitting genuinely localizes
  the trace beyond what broadband sees.

PDF only, full vector. No `fig.suptitle`. No watermark by default.

## 7. Connection to prior tools

| Prior measure | Relation to this scope |
|:---|:---|
| **CTM (`continuous_trace_matrix.py`, scope `2026-04-26_continuous-trace-matrix.md`)** | Provides the per-band Δ_task / Δ_rest vectors that the pool concatenates. This scope is a **stacking adapter** on the controlled (Run-A) CTM output — no new measure, a different aggregation. |
| **§5.3 of `notes_imcoh.tex`** | Computes per-band ρ_split with the same Run-A control. Reports α 8/10, β 8/10, γ_l 7/10 trace direction at q=0.027. This scope answers "does cohort direction survive without band splitting" — yes, 9/10 broadband. |
| **CTM triangle (`audit_33_ctm_triangle.py`)** | Computes per-band cohort summary (median ρ_split, n_trace_split, paired Wilcoxon vs drift). The pool is a *patient-pooled* aggregation; the triangle is a *band-stratified* aggregation. Complementary, not redundant. |
| **§11 canonical broadband |ImCoh| → LRG pipeline (this scope)** | **Not yet run.** Specified below. Would compute |ImCoh| on the union spectrum [0, f_Nyquist] per phase, build one LRG per phase, recompute the half-split distance shifts, and report a true broadband ρ_split. ≤ 30 min wall time. The proxy here is the cheapest possible direction probe; §11 confirms or rejects the magnitude claim. |
| **VI(k), KC, MSPC, CTM** | Different scales / different operators. The pool acts at the same scale (per-pair upper-triangle) as CTM by construction, then aggregates over bands rather than over τ or k. |

## 8. Implementation plan

### 8.1 Compute script (this scope)

`scripts/01_compute/audit/audit_58_band_agnostic_lrg.py` — implemented
2026-05-08. Library imports:

- `lrg_eegfc.config.const.BRAIN_BANDS_NAMES`, `BRAIN_BAND_TEX_DICT`
- `lrg_eegfc.utils.metrics.hypothesis.wilcoxon_z`
- `lrg_eegfc.utils.scripting.setup_script_env`

CLI:

```bash
conda activate lapbrain
python scripts/01_compute/audit/audit_58_band_agnostic_lrg.py
```

Outputs:

- `data/audit/band_agnostic_lrg/per_patient_summary.csv` —
  one row per patient: `rho_split_pool`, `rho_split_band_mean`,
  `rho_drift_band_mean`, `rho_best_trace_band`, `best_trace_band`,
  `n_total_pairs`.
- `data/audit/band_agnostic_lrg/cohort_summary.csv` —
  one row: cohort scalars (n_pos, n_above_drift, medians, Wilcoxon
  z and p).
- `data/reports/notes_verification_2026-05-08/figures/band_agnostic_lrg_teaser.pdf` —
  two-panel teaser figure.

### 8.2 Library promotion (deferred)

No new helpers needed. The pool stacking is one line of `np.concatenate`
called once per patient. If a second consumer arises, promote
`stack_per_patient` to `lrg_eegfc.utils.metrics.continuous_trace`.

### 8.3 No new tests required

The proxy reuses CTM's per-pair NPZ outputs whose tests live in the
CTM scope. Smoke-check: ρ^pool reproduces the cohort summary CSV when
re-run; both runs match to 1e-12.

## 9. Open questions

1. **Run §11 canonical pipeline?** The proxy answers the direction
   question (yes — 9/10). Magnitude requires §11. Recommended yes if
   the manuscript discusses broadband at all. Estimated runtime
   ≤ 30 min (one |ImCoh| compute + one LRG compute + one half-split
   per phase per patient).
2. **1/f-weighted band-pool variant?** Each band's pair-vector could be
   weighted by `Δf_b / f_b` (1/f power proxy) before concatenation.
   Brings the pool closer to spectral broadband at zero recompute cost.
   Open: is the 1/f reweighting useful or distracting?
3. **Pat_03 sensitivity check?** Repeat without Pat_03; does the
   cohort direction count survive? (Likely yes given Pat_03 ρ^pool =
   +0.415 — its removal weakens cohort, doesn't flip.)
4. **TARR / KC taxonomy under broadband?** Direction-only proxy;
   leaf-set classification under broadband is a separate scope.

## 10. Results — first-pass run (2026-05-08)

Cohort `n = 10` under `imcoh_abs`, half-baseline-split CTM (Run A),
band-pool proxy.

### 10.1 Per-patient table

| patient | ρ^pool | ρ^band-mean | ρ^drift | best trace-band | ρ at best |
|:---|---:|---:|---:|:---:|---:|
| Pat_02 | +0.378 | +0.186 | +0.078 | γ_l | +0.693 |
| Pat_03 | +0.415 | +0.297 | −0.035 | β   | +0.463 |
| Pat_05 | +0.527 | +0.383 | −0.003 | γ_l | +0.826 |
| Pat_06 | +0.748 | +0.635 | −0.008 | α   | +0.831 |
| Pat_07 | +0.047 | +0.017 | −0.096 | β   | +0.185 |
| Pat_08 | +0.197 | +0.189 | +0.117 | β   | +0.503 |
| Pat_10 | +0.086 | +0.006 | +0.070 | α   | +0.158 |
| Pat_13 | +0.098 | −0.025 | −0.010 | β   | +0.185 |
| Pat_14 | +0.006 | −0.053 | −0.054 | α   | +0.240 |
| Pat_15 | **−0.133** | −0.019 | −0.152 | β   | +0.144 |

### 10.2 Cohort scalars

| Quantity | Value |
|:---|---:|
| `n_pool_pos` (ρ^pool > 0) | **9 / 10** |
| `n_pool_above_drift` (ρ^pool > ρ^drift) | **10 / 10** |
| `n_band_mean_pos` (band-mean ρ > 0) | 7 / 10 |
| median ρ^pool | +0.148 |
| median ρ^band-mean | +0.102 |
| median ρ^drift | −0.009 |
| Wilcoxon z (pool > 0) | +2.293 |
| Wilcoxon p (pool > 0) | **0.011** |
| Wilcoxon z (pool > drift) | +2.803 |
| Wilcoxon p (pool > drift) | **0.0025** |
| Wilcoxon z (band-mean > drift) | +2.293 |
| Wilcoxon p (band-mean > drift) | **0.011** |

### 10.3 Comparison to per-band §5 numbers

| band | per-band ρ_split>0 | median ρ_split | per-band Wilcoxon vs drift |
|:---|:---:|---:|---:|
| α | 8/10 | +0.115 | p = 0.007 |
| β | 8/10 | +0.222 | p = 0.014 |
| γ_l | 7/10 | +0.140 | p = 0.010 |
| **broadband-pool** | **9/10** | **+0.148** | **p = 0.0025** |

The broadband-pool *direction count* (9/10 above 0; 10/10 above drift)
exceeds every per-band count individually, while the *median magnitude*
(+0.148) lands between α (+0.115) and β (+0.222). Reading: pooling
across bands rescues patients whose strongest single band swung
opposite to the cohort gate (Pat_13 ρ_β = +0.185 but ρ_α = +0.005;
ρ^pool = +0.098 — pulled positive by β even though α was null), at
the cost of magnitude attenuation when bands disagree.

### 10.4 Interpretation

- **The cohort task-trace direction is band-agnostic at the proxy tier.**
  Without a single band-specific gate, 9/10 patients show the
  trace-direction broadband ρ^pool > 0; 10/10 are above their own
  within-session drift floor. The cohort Wilcoxon p = 0.0025 against
  drift is *stronger* than any single band's per-band Wilcoxon
  (α 0.007, β 0.014, γ_l 0.010).
- **Band splitting remains useful for magnitude localization.**
  Pat_06 ρ_α = +0.831, ρ^pool = +0.748 — broadband recovers most.
  Pat_05 ρ_γl = +0.826, ρ^pool = +0.527 — broadband loses ≈ 0.3
  because the 5 other bands dilute. The §5 per-band gate is therefore
  a *resolution* statement, not a *direction* statement.
- **Pat_15 is the only proxy-negative.** ρ^pool = −0.133, ρ^drift =
  −0.152. Even Pat_15's ρ^pool > ρ^drift (drift is more negative);
  the drift floor catches the structural noise that pulls Pat_15
  negative in both directions. Cohort gate `pool > drift` therefore
  retains all 10 patients above floor.
- **Pat_03 remains directional.** ρ^pool = +0.415 despite the 1024 Hz
  γ_h coverage halving caveat. Cohort robust without exclusion.

### 10.5 Pass / fail against the scope contract (§4 caveats)

| Caveat | Status |
|:---|:---|
| 4.1 Band-pool ≠ spectral broadband | Disclosed; §11 specified, not run. **DEFERRED.** |
| 4.2 1/f weighting absent | Disclosed; equal-weight pool. **OK as-is.** |
| 4.3 Jensen's inequality (per-bin abs first) | Inherited correctly from cached per-band caches. **OK.** |
| 4.4 Pat_03 1024 Hz | ρ^pool = +0.415, direction holds. **OK.** |
| 4.5 Pair independence | Per-patient ρ^pool used at cohort level; n=10 inference unit. **OK.** |
| 4.6 Direction-only | Acknowledged; magnitude story requires §11. **OK.** |

**Headline.** *The trace direction reported per-band in §5 of
`notes_imcoh.tex` (α / β / γ_l at q=0.027) survives — and intensifies
in cohort consistency, 9/10 → 10/10 above drift — when bands are
pooled into a single broadband ρ_split. Band splitting is a
magnitude-localization tool, not a direction-detection requirement.*

## 11. §11 — Canonical broadband |ImCoh| → LRG pipeline (parked, ≤ 30 min)

The proxy in this scope answers direction; the canonical pipeline
answers magnitude. Specification:

### 11.1 Broadband |ImCoh| construction

For each `(p, φ)`:

1. Load each per-band freq-resolved cache
   `{band}_{phase}_imcoh_freqresolved_nperseg-4096.npy`,
   shape `(N, N, F_band)` — these are the *signed* ImCoh per bin.
2. `np.abs` per bin → `|ImCoh|^band` of same shape.
3. Concatenate along the frequency axis: `|ImCoh|^union` of shape
   `(N, N, Σ_b F_b)`. Pat_03 has half the F_γh count (≈ 220 vs 441);
   acceptable given uniform per-bin treatment.
4. Average over the union frequency axis:
   `A^{p, φ}_broad := mean(|ImCoh|^union, axis=2)`. One adjacency per phase.
5. Verify symmetry and `A_ii = 0`; clip to `[0, 1]`.

This is a 1/f-weighted spectral broadband by virtue of the bin counts
per band (low_gamma 101 bins, alpha 11 bins) — δ contributes ≈ 7/604,
γ_h ≈ 441/604 ≈ 73% by bin count. To rebalance to band-equal weighting
add a per-bin weight `w_f = 1 / F_{b(f)}`; the band-equal variant
keeps the comparison with the proxy aligned.

### 11.2 LRG step

For each `(p, φ)`:

```python
from lrg_eegfc.workflow.lrg import compute_lrg_analysis
result = compute_lrg_analysis(
    A_broad, patient=p, phase=phi, band="broadband",
    fc_method="imcoh_abs",                          # tagging only
    cache_root=BROADBAND_LRG_CACHE,
    use_cache=True,
)
D_broad_phi = result.ultrametric_matrix             # condensed; squareform()
```

LRG cache: `data/cache/imcoh_lrg_broadband/{patient}/broadband_{phase}_lrg_imcoh-abs.npz`.
~10 phases × 10 patients = 100 LRG runs at ~ 5–10 s each = ≈ 15 min.

### 11.3 Half-split

For each `(p, half ∈ {A, B})`, repeat steps 11.1–11.2 with rest_pre time
series cut to its first / second half (shared with `h2e_split_half.py`'s
halving convention). Cache:
`data/cache/imcoh_lrg_broadband_halves/{patient}/broadband_{phase}_{half}_lrg_imcoh-abs.npz`.

### 11.4 Per-patient ρ_split-broadband

```
Δ_task^{p, broad} := D_broad_test - D_broad_pre_A
Δ_rest^{p, broad} := D_broad_post - D_broad_pre_B
ρ_split^{p, broad} := Spearman( vec_upper(Δ_task^{p, broad}),
                                vec_upper(Δ_rest^{p, broad}) )
```

Single (p)-indexed scalar; cohort = 10 values; same Wilcoxon as proxy.

### 11.5 Decision rule

If `ρ_split^{p, broad}` cohort matches the proxy (9/10 above 0;
10/10 above drift), the proxy is validated as a broadband direction
estimator. If not, document the discrepancy: it would mean the
band-pool is over- or under-counting via per-pair multiplicity. Either
outcome is publishable.

### 11.6 Total runtime

| step | wall time |
|:---|---:|
| 11.1 broadband |ImCoh| | ≤ 1 min (numpy, no recompute from time series) |
| 11.2 LRG (10 patients × 4 phases × full + halves = 80 runs) | ≤ 20 min |
| 11.3 half-split LRG (already part of 11.2) | — |
| 11.4 ρ_split + cohort | < 1 min |
| **Total** | **≤ 25 min** |

Script: `scripts/01_compute/audit/audit_59_canonical_broadband_lrg.py`
(slot reserved). Status: **DRAFT, not run.**

## 12. References

- §5 of `notes_imcoh.tex` (PDF at `.agents/references/notes_imcoh_260508.pdf`):
  per-band CTM ρ_split with α / β / γ_l trace at q=0.027.
- Continuous-trace matrix scope (`2026-04-26_continuous-trace-matrix.md`):
  full per-band ρ_split definitions, three-control consolidation,
  `controls_summary.csv` schema.
- `controls_band_stats.md`: per-band cohort table reproduced in §10.3.
- Nolte et al. 2004 (signed ImCoh); Ewald et al. 2012 (|ImCoh|);
  Bastos & Schoffelen 2016 (|ImCoh| convention).
- Villegas et al. 2023 / 2025 (LRG framework; substrate for the
  broadband adjacency feed).
