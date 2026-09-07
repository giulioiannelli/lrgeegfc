---
name: 2026-08-31_scale-local-trace-readouts
type: scope
era: PAPER_FINALIZATION (Wave 0, lane W0-S)
status: current
created: 2026-08-31
updated: 2026-09-03
scope: Pre-registration for lane W0-S. Asks whether the diffusion-scale axis is genuinely flat or whether the incumbent readout (a global Spearman over all contact pairs) cannot see scale. Defines the dilution diagnostic, three candidate scale-local readouts, and — before any number — the three criteria a readout must satisfy to replace the incumbent, plus the criterion under which the dilution hypothesis is declared wrong. Amendment B (2026-09-03) extends the pre-registration to a *contrast between functionals* (T_learn - T_test, and T_infspec alone), which Part A did not cover, and pre-registers the low_gamma band-selectivity call.
pointers:
  - .agents/reports/2026-08-25_w0c-cohort-gate-and-tau.md
  - .agents/preprint/locked/PIPELINE_CONTRACT.md
  - src/lrg_eegfc/utils/fc/heat_multiscale.py
  - src/lrg_eegfc/utils/metrics/cohort_gate.py
  - scripts/01_compute/paper_final/w0s_01_scale_locality_grid.py
  - data/paper_final/lane_e_encinf/grid/cells
  - data/paper_final/lane_e_encinf/sham/cells
  - .agents/reports/2026-08-31_lane-s-scale-variability.md   # OUTCOME: negative, all criteria failed
  - .agents/reports/2026-09-03_lane-s-contrast-scale-structure.md   # OUTCOME of Amendment B: negative, B1/B2/B3 all failed; low_gamma = gate artefact
---

# Scale-local trace readouts — is the scale axis flat, or is the readout blind to it?

## Head

Two W0-C numbers sit badly together: across the scale sweep the hierarchy changes thirty-fold (the tree resolves 118 components at the finest scale and 4 at the coarsest), yet the cross-phase statistic barely moves (a 28-point sweep is worth 1.2–1.9 independent tests, mean cross-scale correlation +0.70 to +0.93). This scope pre-registers the test of the obvious suspect — that `rho_sym`, a Spearman over all ~7 000 contact pairs, is scored on nearly the same information at every scale because the global tree ordering is preserved as the diffusion coarsens — and, critically, states in advance what would count as a better instrument and what would count as the dilution hypothesis being wrong.

**This document is written before any number is computed. Nothing below may be edited in response to a result; corrections are appended with a date.**

---

## 1. Five-point critical preamble

**(1) The claim under test.** The near-constancy of the cross-phase trace along the diffusion-scale axis `s = tau * lambda_max` is a property of the *readout*, not of the phenomenon. A readout that is local in the hierarchy will show scale structure that `rho_sym` averages away.

**(2) The null.** Every gated readout is referenced to the same matched-strength surrogate: the dense FC matrix is 4-cycle shuffled at fixed node strength, then pushed through the *identical* sparsify -> combinatorial Laplacian -> heat kernel -> UPGMA -> stratify -> rank-correlate pipeline. The cohort gate is a one-sided signed-rank test on the per-patient margin `m_k = o_k - median_r s_{k,r}` (`utils.metrics.cohort_gate`). **No readout is ever tested against zero.** W0-B showed all four cross-phase functionals are significantly positive at 16/16 scales on a no-task sham arc with temporal order destroyed, so zero is not the reference and never was.

**(3) The strongest plausible alternatives the null must control for.**

- **A3.1 — "a scale-local readout is just a noisier readout".** Restricting a correlation to a subset of pairs raises its variance. Two distinct consequences, and they need different controls. *For the p-value*: controlled by construction, because the surrogate ensemble passes through the same stratification with the same pair counts, so the margin subtracts the stratum-specific noise floor. *For the scale-structure claim*: **not controlled by the margin at all.** Independent noise decorrelates the columns of the per-patient × per-scale margin matrix, so the effective-number-of-independent-tests statistic `n_eff` **rises automatically with noise**. Any criterion built on `n_eff` alone would be satisfied by simply making the statistic worse. The mitigation is stated as part of criterion (b) below: `n_eff` must be **null-referenced** — compared against the `n_eff` of the same readout computed on held-out surrogate realizations, which carry the identical noise level and no signal.
- **A3.2 — "the incumbent's flatness is a power artifact".** Not available: the incumbent does clear at beta (W0-C: cluster p = 0.0032, worst leave-one-out p = 0.0091), so it is not a powerless statistic. If the flatness were pure noise the profile would be flat *and* null; it is flat *and* significantly positive.
- **A3.3 — "the strata are themselves scale-dependent, so cross-scale decorrelation is induced by the definition, not by the data".** Real, and the reason the surrogate must be stratified with its own trees rather than with the observed graph's strata. Stated again as a caveat in §6.

**(4) What the null cannot do.** Matched-strength conditions on the finished N×N connectivity matrix. Every verdict in this lane is of the form "given this FC estimate": nothing here tests the coherency estimator, the band split, session nonstationarity, drift or artifact epochs. That is W0-B's timeseries ladder and it is a live dependency. Separately, and following W0-C §2.2: a flat result means "no scale structure detectable at n = 10", **never** "scale-invariant". Equivalence cannot be established at this cohort size, and cross-patient inconsistency is indistinguishable from within-patient flatness under a cohort test.

**(5) Falsification, and the limitations that remain either way.** The dilution hypothesis is declared **wrong** under the S1 condition in §5.1. The lane's positive claim is declared **failed** if no construction satisfies (a) + (b) + (c) of §5.2, in which case the flatness is reported as the headline. Limitations that survive any outcome: n = 10 puts an exact signed-rank floor at 1/1024; the merge-height readout matches heights by rank position, not by cluster identity; the stratification is read off the baseline arm's tree, so a phase that reorganises the hierarchy wholesale is scored in the baseline's coordinates.

---

## 2. Notation

| symbol | meaning |
|---|---|
| `N` | contacts in a patient's implant (cohort median 118 on the mst-union backbone) |
| `M = N(N-1)/2` | contact pairs; the length of every condensed vector below |
| `W_ph` | sparsified adjacency of phase `ph`, from `canonical_graph` |
| `L = D - W` | combinatorial Laplacian (never a normalised variant) |
| `s = tau * lambda_max` | dimensionless diffusion scale |
| `rho(s) = e^{-tau L} / Tr` | diffusion density operator at scale `s` |
| `D(s) = (1 - delta_ij) / rho_ij(s)` | communication distance |
| `Z_ph(s)` | UPGMA (average-linkage) linkage of `D_ph(s)` |
| `C_ph(s) in R^M` | condensed **cophenetic** distances of `Z_ph(s)` |
| `h_ph(s) in R^{N-1}` | sorted merge heights of `Z_ph(s)` |
| `A, B` | split halves of `rest_pre`; `t = task_test`, `p = rest_post` |
| `N_eff(s)` | `exp` of the von Neumann entropy of `rho(s)` — resolved components |
| `m(s)` | mean heat-kernel row participation ratio — communication neighbourhood |
| `f` | mst-union edge fraction, `f in {0.07, 0.10, 0.14, 0.20}` |
| `R` | matched-strength realizations per (cell, fraction) |

Phases carry a patient and band index throughout; both are suppressed.

---

## 3. The incumbent, and where the suspicion enters

The locked cross-phase trace is the symmetric split-half Spearman

```
T(s) = 1/2 [ rho_S( C_t(s) - C_A(s), C_p(s) - C_B(s) )
           + rho_S( C_t(s) - C_B(s), C_p(s) - C_A(s) ) ]
```

with `rho_S` Spearman over all `M` pairs. `T` is `rho_sym` in `heat_multiscale.py` and `T_probe` in `cross_phase_functionals`; the two are numerically identical.

Spearman is an inner product of normalised mean-centred ranks, so `T(s)` decomposes **exactly and additively** over pairs:

```
T(s) = sum_{e=1..M} c_e(s),
c_e(s) = 1/2 [ (a_e b_e)/(|a||b|) + (a'_e b'_e)/(|a'||b'|) ]
```

where `a = rank_c(C_t - C_A)`, `b = rank_c(C_p - C_B)` and the primed pair is the other arm (`rank_c` = mean-centred ranks). This identity is the backbone of the whole lane: any partition of the pair set gives shares that sum to the incumbent, so criterion (a) below is satisfied *by construction* for the contribution readout and has to be *earned* for the others.

The suspicion is that `c_e(s)` is dominated, at every `s`, by the same pairs — those whose cophenetic distance is set by the handful of top merges — because UPGMA on `D = 1/rho` puts most pairs at or near the root, and the global ordering of the tree is largely preserved as `s` grows.

---

## 4. Definitions

### 4.1 Tree octave of a pair (the stratification)

For a reference tree `Z(s)` with `N-1` merges taken in increasing height, let merge `l in {1..N-1}` join two clusters. Define

```
k(l) = N - l + 1        (number of clusters present JUST BEFORE merge l; root has k = 2)
```

and for a pair `e = (i,j)` let `l(e)` be the merge at which `i` and `j` first share a cluster. Define the pair's **tree octave**

```
o(e) = floor( log2 k(l(e)) )  in  {1, ..., floor(log2 N)}
```

Reading: `o = 1` (`k in [2,4)`) means the pair only joins when at most three clusters remain — a pair spanning the coarsest split of the implant. `o = 6` (`k in [64,128)`) means the pair joins while more than 64 clusters are still separate — near-neighbours in the hierarchy. Six octaves cover `N ~ 118`.

`o` is defined on the **baseline arm's** tree so that both members of a cross-phase pair are scored on the same index set: arm 1 (`C_t - C_A` vs `C_p - C_B`) stratifies by `Z_A(s)`, arm 2 by `Z_B(s)`. The octave boundaries are fixed powers of two, not quantiles, so no boundary is fitted to the data.

**Amendment, 2026-08-31, before any number was computed.** A UPGMA tree on a communication distance puts most pairs near the root, so the *top* octaves are pair-starved: on a 118-leaf tree octave 6 can hold only a few dozen of the ~6 900 pairs, and a rank correlation on a few dozen pairs is a poor estimator at exactly the fine end of the scale axis where the band-pass is supposed to look. The identical decomposition is therefore **also** run on an **equal-count** stratification — quintiles `q in {1..5}` of the same merge-level ordering, so every stratum holds `M/5 ~ 1 380` pairs whatever the tree's shape — giving `Qcon(s,q)` and `Qloc(s,q)` alongside `Ccon(s,o)` and `Tloc(s,o)`. Its active stratum is `q*(s)` = the quintile containing merge level `l* = N - N_eff(s) + 1`, the same derivation as `o*(s)`. Both stratifications are pre-registered here, before the first result; neither may be selected after the fact, and both are reported whatever they show.

**Why this is the right band-pass.** `N_eff(s)` is the number of components the diffusion resolves at `s`. A pair with `k(l(e)) > N_eff(s)` was already merged below the current resolution; a pair with `k(l(e)) < N_eff(s)` is still unresolved. The pairs being *resolved at* `s` are those with `k ~ N_eff(s)` — i.e. octave `o* (s) = floor(log2 N_eff(s))`. That diagonal of the (scale × octave) surface is the derived scale-local readout; it is a consequence of the definitions, not a tuned window.

### 4.2 Readout R1 — stratified contribution `Ccon(s, o)`

```
Ccon(s, o) = sum_{e : o(e) = o} c_e(s)        (arm-symmetrised as in §3)
```

Range `[-1, 1]` per stratum; `sum_o Ccon(s, o) = T(s)` exactly. This is simultaneously the S1 diagnostic (how the incumbent's value is distributed over tree levels) and a scale-local readout (each octave's own cross-phase evidence).

### 4.3 Readout R2 — within-octave local trace `Tloc(s, o)`

Restrict to `E_o = {e : o(e) = o}`, **re-rank inside the stratum**, and correlate:

```
Tloc(s, o) = 1/2 [ rho_S( (C_t - C_A)|_{E_o^A}, (C_p - C_B)|_{E_o^A} )
                 + rho_S( (C_t - C_B)|_{E_o^B}, (C_p - C_A)|_{E_o^B} ) ]
```

Genuinely different from R1: re-ranking discards the global ordering entirely, so `Tloc` asks whether — *among pairs sitting at the same level of the hierarchy* — the task change predicts the rest change. It does **not** sum to `T(s)`, so criterion (a) must be checked empirically.

### 4.4 Readout R3 — merge-height trace `Thei(s)`

Not pairwise at all. Using the sorted merge heights `h_ph(s) in R^{N-1}` and their logs `H_ph = log h_ph`:

```
Thei(s) = 1/2 [ rho_S( H_t - H_A, H_p - H_B ) + rho_S( H_t - H_B, H_p - H_A ) ]
```

This scores *changes in the merge heights themselves* — where in the tree's own height spectrum the task inflated or deflated the hierarchy, and whether that pattern persists. Logs are used so the differences are relative changes; without them the vector is dominated by the absolute height scale, which grows exponentially in `s`.

### 4.5 Diagnostics (observed only, no gate)

- `xs_coph(s_i, s_j)` = Spearman between `C_A(s_i)` and `C_A(s_j)` — how much the readout's *input* changes with scale.
- `xs_reorg(s_i, s_j)` = the same on the symmetric task reorganisation `u(s) = 1/2[(C_t - C_A) + (C_t - C_B)]`.
- `n_distinct(s)` = `exp` of the Shannon entropy of the pair-count distribution over the `N-1` merge heights — the effective number of distinct cophenetic values a Spearman over `M` pairs actually sees.
- `pairfrac(s, o)` = share of pairs in each octave.
- `N_eff(s)`, `m(s)` — the locked scale units, recomputed here, never sourced.

---

## 5. Pre-registered decision rules

### 5.1 When the dilution hypothesis is wrong

Declared **wrong**, and abandoned rather than pursued, if **both** hold on the observed data:

1. the contribution shares `|Ccon(s,o)| / sum_o |Ccon(s,o)|` are approximately flat across octaves — no octave holds more than 2× the share of the median octave — at the majority of scales; **and**
2. the per-octave local traces are as cross-scale-correlated as the global one: median over octaves of the mean off-diagonal cross-scale correlation of `Tloc(., o)` is within 0.10 of the incumbent's.

If that is the finding it is reported plainly and the band-pass constructions are dropped.

### 5.2 When a readout replaces the incumbent as the scale-resolved instrument

All three, no exceptions, no partial credit:

**(a) Same phenomenon, not a new one.** Aggregated over scales, the readout must agree with the incumbent. Operationally: (i) the cohort Spearman between the readout's scale-aggregated per-patient margin and the incumbent's, computed across the 10 patients, is `>= 0.5` in beta; and (ii) the per-band axis-cluster verdict pattern matches the incumbent's on the two bands the project treats as known — beta clears, theta does not.

**(b) More independent information across scales, above its own noise.** Measured exactly as W0-C did — the eigenvalue participation ratio `n_eff_pr` of the cross-scale correlation matrix of the per-patient margins (`cohort_gate.effective_tests`). Two parts, and both are required:
   - **absolute**: `n_eff_pr` at least 2× the incumbent's in at least half the bands;
   - **null-referenced**: `n_eff_pr` on the observed margins must exceed the distribution of `n_eff_pr` computed from held-out surrogate realizations of the *same readout* (promote realization `r` to "observed", margin it against the remaining `R-1`, form the same `K x n_scales` matrix), at `p < 0.05` in at least half the bands. **Without this second part the criterion is trivially satisfied by any noisier statistic**, per A3.1.

**(c) Calibrated, and margin-referenced.** Held-out-realization false-positive rate at nominal 0.05 must be `<= 0.10` in at least 95% of cells (`cohort_gate.calibrate_from_surrogates`, `refuse_uncalibrated=True`), and the readout is gated only as a margin against the matched-strength null. A readout failing (c) has its p-values withheld, not caveated.

**A readout that shows scale variation but fails (a) or (c) is a worse instrument, not a better one, and is reported as a failed construction.**

### 5.3 Reporting rule

Every construction attempted is reported, including the ones that fail, and in the order they were run. No construction may be added after the first verdict is read without being labelled post-hoc.

---

## 6. Properties, caveats, failure modes

| # | issue | mitigation |
|---|---|---|
| 1 | **Noise inflates `n_eff`.** A stratified statistic uses fewer pairs and is noisier; noise decorrelates scales. | Null-referenced `n_eff` in criterion (b). Non-negotiable. |
| 2 | **Strata are scale-dependent**, so the pair set changes with `s` and could induce decorrelation by itself. | The surrogate is stratified by *its own* trees at each `s`, so the induced component is present in the null and removed by the margin; and `pairfrac(s,o)` is reported so the reader sees the pair budget. |
| 3 | **Octave 1 and 6 have small or degenerate pair counts** at some scales (at coarse `s` almost everything is in octave 1; at fine `s` almost nothing is). | Cells with fewer than 100 pairs in a stratum return `NaN` and are excluded from that cell, not imputed. The count is reported. |
| 4 | **`Thei` matches heights by rank position, not cluster identity.** Two phases' `l`-th merges need not join the same nodes. | Stated as a limitation; `Thei` is a *shape-spectrum* comparison and is never described as a per-cluster result. |
| 5 | **Stratification is read off the baseline arm.** A phase that reorganises the hierarchy globally is scored in the baseline's coordinates. | Deliberate — the alternative (stratifying by the task tree) makes the reference depend on the effect being measured. Both arms are used and symmetrised. |
| 6 | **Cophenetic ties.** UPGMA gives at most `N-1` distinct cophenetic values across `M ~ 7000` pairs, so Spearman is run on a heavily tied vector. | This is measured, not assumed: `n_distinct(s)` is a headline diagnostic of S1. |
| 7 | **Single-fraction numbers are tuned numbers.** | Everything is knob-integrated over `f in {0.07, 0.10, 0.14, 0.20}` via `canonical_graph_ensemble`; the per-patient margin is the median over fractions, and the across-fraction spread is reported. |
| 8 | **Whole-grid BH over (band × scale) overcharges ~20×** at this `n_eff` (W0-C §1.4). | Both families reported: whole-grid BH and the per-band `axis_cluster_gate`. Neither is hidden. |
| 9 | **`s < 1` is the raw-FC limit.** As `tau -> 0`, `K ~ I + tau W`, so `D ~ 1/W` and the cophenetic tree is UPGMA on a monotone transform of the raw FC. | Stated wherever fine-scale results are reported; it is why the fine grid is a *diagnostic* of the axis, not an independent measurement. |

---

## 7. Pseudocode

```
GRID:  S = 32 log-spaced scales spanning [0.02, 180]
       F = (0.07, 0.10, 0.14, 0.20)            # locked plateau fractions
       R = matched-strength realizations per (cell, fraction)

for each patient P, band b:
    Wdense[ph] <- dense imcoh_abs FC for ph in {A, B, task_test, rest_post}

    for each realization index r in {observed} u {1..R}:
        if r is observed: Wr[ph] <- Wdense[ph]
        else:             Wr[ph] <- matched_strength_shuffle(Wdense[ph])   # same draw reused across F
        for each fraction f in F:
            G[ph] <- mst_union_backbone(Wr[ph], f)
            (ev[ph], V[ph]) <- eig(combinatorial_laplacian(G[ph]))
            for each scale s in S:
                for ph: Z[ph] <- UPGMA(1 / heatkernel(ev[ph], V[ph], s))
                        C[ph] <- cophenetic(Z[ph]);  H[ph] <- log(sorted_merge_heights(Z[ph]))
                T      <- symmetric_spearman(C_t - C_A, C_p - C_B ; C_t - C_B, C_p - C_A)
                oA, oB <- tree_octave(Z_A), tree_octave(Z_B)
                for octave o in 1..6:
                    Ccon[s,o] <- sum of per-pair Spearman contributions with octave o
                    Tloc[s,o] <- symmetric_spearman restricted to octave o, re-ranked
                Thei[s]   <- symmetric_spearman(H_t - H_A, H_p - H_B ; H_t - H_B, H_p - H_A)
            record T, Ccon, Tloc, Thei  (+ diagnostics for the observed pass only)

AGGREGATE:
    per patient: margin(readout, s) <- median_f [ obs(f) - median_r surr(f, r) ]
    gate:        cohort_margin_gate / gate_grid over (band x scale), both families
    n_eff:       effective_tests on the K x |S| margin matrix, observed AND held-out-surrogate
    calibrate:   calibrate_from_surrogates per cell, refuse_uncalibrated = True
```

---

## 8. Visualization spec

- **Scale × octave surface** of `Ccon` and of `Tloc`, one panel per band, `s` on a log x-axis, octave on y, diverging colormap centred at zero, with `o*(s) = floor(log2 N_eff(s))` overdrawn as a line. *Reading rule*: if the colour is a horizontal band (constant along `s`) the readout has no scale structure; if it tracks the overdrawn diagonal, the trace is carried by whichever pairs the diffusion is currently resolving.
- **Cross-scale correlation matrix** of the per-patient margins, incumbent vs each candidate, one panel each. *Reading rule*: a matrix that is uniformly near +1 is a statistic seeing the same thing at every scale.
- **`n_eff` bar chart** per band: incumbent, each candidate, and the held-out-surrogate null band behind each bar. *Reading rule*: a bar that does not stand clear of its own null band is noise, not information.

No `fig.suptitle`; PDF only, fully vector; `use_lrg_style()`; band colours from `visuals.styles.band_color`.

---

## 9. Connection to prior tools

| tool | relation |
|---|---|
| `rho_sym` / `T_probe` (`heat_multiscale`) | The incumbent. R1 is its exact additive decomposition; R2 and R3 are alternatives, not refinements. |
| `effective_cluster_count`, `communication_neighbourhood_size` | Supply `N_eff(s)`, which *defines* the octave diagonal `o*(s)`. This lane consumes them, never redefines them. |
| `cohort_gate` (`gate_grid`, `axis_cluster_gate`, `effective_tests`, `calibrate_from_surrogates`) | The locked gate. Used unmodified; this lane writes no gate of its own. |
| `fcluster` / ARI / NMI partition metrics | **Deliberately not used.** `feedback_no_partition_metrics_use_rho_coph` forbids hard-partition cross-phase similarity. The octave stratification indexes pairs by tree level but never compares partitions; every readout stays a continuous rank correlation. |
| `2026-06-22_tau-sensitivity-cophenetic-trace.md` | Earlier tau sweep of the same incumbent statistic. This lane does not re-run it; it asks whether that sweep's flatness was a readout artifact. |
| W0-C §2.2 (characterization) | W0-C tested whether the *shape* of `T(s)` is band-discriminative and found nothing. This lane asks the prior question: whether `T(s)` has a shape to discriminate. |

---

## 10. Implementation plan

| stage | script | output |
|---|---|---|
| master compute | `scripts/01_compute/paper_final/w0s_01_scale_locality_grid.py` | `data/paper_final/lane_s_scale/grid/cells/<patient>__<band>.npz`, `descriptive.csv`, `config.json` |
| S1 dilution diagnostics | `scripts/01_compute/paper_final/w0s_02_dilution_diagnostics.py` | `dilution_*.csv` |
| S3 fine-grid incumbent | folded into the master grid (the 32-point grid extends to `s = 0.02`) | `gate_grid_incumbent.csv` |
| S0 criteria + S4 verdict | `scripts/01_compute/paper_final/w0s_03_scale_local_verdict.py` | `criteria.csv`, `gate_*.csv`, `n_eff_*.csv` |

Library entry points reused unmodified: `workflow.substrate.{canonical_graph_ensemble, canonical_scale_grid, CANONICAL}`, `utils.fc.heat_multiscale.{laplacian_eig, linkage_at_scale, cophenetic_at_scale, effective_cluster_count, communication_neighbourhood_size}`, `utils.metrics.surrogate.matched_strength_shuffle`, `utils.metrics.cohort_gate.*`. New general helpers, if any, are promoted to `utils/metrics/tree.py` (tree octaves) with general names — never `laneS`, `w0s` or `paper`.

---

## 11. Open questions (declared, not resolved)

1. **`R`**: set after a timing probe; stated in `config.json` and in the report. Not tuned to a result.
2. **Octave count**: fixed at powers of two, so it is `floor(log2 N)` and varies with implant size (5–7 across the cohort). Cells are compared within an octave index, not renormalised.
3. **`Thei` on log heights** is a choice; the linear-height variant is not run, and that is recorded as untested rather than claimed equivalent.
4. **The fine end of the grid** (`s < 1`) is the raw-FC limit by construction. Whether the trace there is "the hierarchy" at all is a framing question this lane surfaces but does not settle.

---

## 12. Amendment B (2026-09-03) — scale structure in a *contrast between functionals*

**Written before any contrast statistic was computed.** Part A (§1–§11) tested scale-locality of **one functional at a time**: the incumbent `T` and five scale-local reconstructions of it. It returned a negative. A contrast between two functionals is a **different object** and is not covered by that verdict, because two quantities can each be flat along τ while their difference is not — the flat parts cancel and only the scale-dependent residue survives. This amendment pre-registers the test of that object.

### 12.1 The object

Lane E's five-phase grid (`data/paper_final/lane_e_encinf/grid/cells/*.npz`) stores, per patient × band, `obs (nF, nS, nFunc)` and `surr (nF, R, nS, nFunc)` over `funcs = (T_test, T_learn, T_infspec, T_infspec_pe)` on the locked contract grid `s ∈ [1, 180]`, 16 points, `R = 200`. The surrogate realization index `r` is **shared across functionals** within a cell, so a contrast can be formed with its pairing intact.

The **paired contrast** (primary) forms the difference at the raw level, before any surrogate subtraction, so that the locked `patient_margin` contract stays exact rather than approximated:

    C_k(f, s)      = obs_k[f, s, T_learn] − obs_k[f, s, T_test]
    C_k(f, r, s)   = surr_k[f, r, s, T_learn] − surr_k[f, r, s, T_test]
    M_k(s)         = median_f C_k(f, s) − median_r median_f C_k(f, r, s)

The **unpaired contrast** (sensitivity only) is `margin[T_learn] − margin[T_test]`, which subtracts two independently-drawn surrogate medians and therefore discards the pairing. It is reported beside the paired one; if the two disagree, the paired one stands and the disagreement is reported as a caveat, because the unpaired variant is strictly the noisier estimator of the same quantity.

`T_infspec` alone is carried as a second object under the identical protocol, since the lead cites a slope in it independently of the contrast.

### 12.2 Five-point critical preamble

1. **The claim.** `C(s) = T_learn(s) − T_test(s)` carries scale structure that neither functional carries alone, and specifically **changes sign** along τ; and `T_infspec(s)` has a real monotone scale slope.
2. **The null.** Two, and both are required. (a) **Held-out matched-strength realizations of the same contrast**: promote surrogate draw `r₀` to the observed slot and re-reference to the remaining `R−1`. This is the contrast, at the contrast's own noise level, with no task-order information. (b) **Lane E's ordered sham** (`data/paper_final/lane_e_encinf/sham/cells/*.npz`): a fake five-phase arc carved from one rest recording with block order preserved and no task, which reproduces 89 % of β `T_infspec`. Its own contrast profile is the reference for any scale-structure statistic.
3. **The strongest plausible alternative.** A difference of two noisy quantities is noisier than either of them. Independent noise decorrelates the columns of the per-patient × per-scale margin matrix, so the effective number of independent scales `n_eff` **rises with noise**; and a noisier per-patient profile has a larger `|Spearman(margin, log s)|` by chance. Both of the statistics that would be cited as evidence for the claim move in the claim's direction for a purely mechanical reason. This is the exact trap that caught four of the five Part-A candidates. The second alternative is **within-recording drift**: `rest_pre → rest_post` ordering produces a monotone τ-trend with no task involved.
4. **Does the null control for it — by mechanism.** (a) controls the noise inflation exactly, because the held-out draw *is* the contrast at the contrast's noise level, so the inflation is already inside the null distribution; this is not a vibes argument but an identity of construction. (a) **cannot** control drift: matched-strength surrogates are redrawn per phase and carry no recording-order information. (b) controls drift by construction, because the ordered sham has real within-recording drift and no task. Neither null controls for **low power**: with `K = 10` patients and an axis Part A measured as worth ≈ 1.1 independent tests, a negative means *not detectable at n = 10*, never *proven absent*.
5. **What would falsify the claim, and what remains.** The contrast claim dies if `n_eff(C)` fails to exceed its own held-out noise floor, **or** if the reversal statistic fails to exceed the ordered sham's. It survives only if **both** hold. Remaining regardless of outcome: no timeseries-level null (W0-B's open dependency); `s < 1` is unsampled on Lane E's grid, so this amendment speaks only for `s ∈ [1, 180]`; and the sham exists for four bands (δ, θ, α, β) only, so low_γ and high_γ get null (a) but not null (b).

### 12.3 Pre-registered criteria

- **B1 — more independent cross-scale information.** `n_eff(C)` must exceed the 95th percentile of its own held-out-null `n_eff` (upper-tail `p < 0.05`) in at least one band, after BH over the six bands. The `n_eff` of each raw functional is reported beside it, so that "the contrast decorrelates more than its parts" is visible as a number rather than asserted. *This is Part A's criterion (b), unchanged, on a new object.*
- **B2 — reversal, not merely slope.** A sign flip is the stronger and more falsifiable claim, and it is tested as one. All three must hold in the same band: (i) a **negative** supra-threshold cluster and a **positive** supra-threshold cluster each clear `axis_cluster_gate` (the locked sign-flip cluster-mass test, run on `M` and on `−M`); (ii) the negative cluster lies at **smaller** `s` than the positive one; (iii) the per-patient zero-crossing scales `s*_k` are **more concentrated** (smaller MAD in `log s`) than under null (a) *and* than under null (b). Patients with no crossing are counted and reported, never silently dropped.
- **B3 — sham-referenced.** Every scale-structure statistic is reported as a pair, `p_vs_surrogate` and `p_vs_sham`. A statistic that clears (a) but not (b) is **within-recording drift** and is reported as such, not as scale structure.
- **Calibration.** As in §5.2: cells failing held-out-realization calibration are withheld, not caveated.

### 12.4 Reporting rule (binding)

Every functional and every contrast attempted is reported, including failures: `T_test`, `T_learn`, `T_infspec`, `T_infspec_pe`, the paired `C = T_learn − T_test`, and the unpaired variant. Per-scale, never best-scale. For each object: the number of **distinct values** it takes, its **split-half reliability**, and its scale-variation number are reported together, so that a scale-variation figure can never be read without the noise level that generates it. A flat answer is the expected outcome under the Part-A result and will be the headline if that is what the numbers say; "the contrast is flat too" is a cleaner and more final negative than Part A alone and is to be stated that way.

### 12.5 The low_γ ownership question (declared in advance)

Part A's own grid puts low_γ at margin `+0.096` against β's `+0.105` with a cluster spanning the whole axis and `q = 0.063`. Part A tabulated it and never discussed it. Combined with W0-A's finding that the band-selectivity window is only 0.51 octaves wide, two lanes now independently suggest band-selectivity is weaker than the project claims. The call is pre-registered here so it cannot be made after seeing which answer is convenient:

- low_γ is a **genuine near-miss** — a real threat to the second result — if it passes calibration, its cluster survives leave-one-patient-out, its margin is not carried by a minority of patients, and the **paired** per-patient difference `β − low_γ` fails to clear. That combination says the gate is behaving correctly and the two bands are not separable at `n = 10`.
- low_γ is a **gate artefact** if it fails calibration, or its clearing collapses under LOO, or the paired `β − low_γ` difference clears. That combination says the near-miss is an artefact of reading two marginal `q` values side by side rather than testing their difference.
