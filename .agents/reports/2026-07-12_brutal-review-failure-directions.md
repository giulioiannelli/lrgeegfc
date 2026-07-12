---
name: 2026-07-12_brutal-review-failure-directions
type: report
era: IMCOH_ABS × COHORT_N10 (sparsified-arc audit)
status: current
created: 2026-07-12
updated: 2026-07-12
pointers:
  - .agents/reports/2026-07-11_sparsified-arc-D0-D2.md
  - .agents/guides/task-persistence-investigation/2026-07-11_sparsified-arc-master-scope.md
  - CLAUDE.md
---

# Brutal review — every direction this analysis could be failing

**Head (the juice).** The project's load-bearing claim — *"LRG reveals a
band-specific MULTISCALE cophenetic trace of the task in rest_post"* — is **not
currently defensible**, and the reason is structural, not a tuning problem. The
full `imcoh_abs` matrix is **not multiscale at all** (audit_174: one specific-heat
peak, Wigner-semicircle regime; LRG ≈ raw FC at τ_min). The entire multiscale
framing therefore rests on a **percolation sparsification introduced *after* that
finding, to rescue the framing** — and on that backbone the trace (a) **loses
band-selectivity under matched-strength** and (b) is **SNR-confounded / inconclusive
under the drift null**. Meanwhile we have cycled β→OFC ⇒ α scale-tuned ⇒ high_γ on
backbone ⇒ β survives-detrend, a **new candidate every iteration**, each complicated
by the next control. That instability is itself the strongest evidence that we may
be walking a **garden of forking paths** and fitting noise. This document enumerates,
by pipeline stage, every way the result could be an artifact, the cheapest test that
would falsify it, and a severity. It is written for a fresh agent to investigate one
lane at a time.

**One-line recommendation:** run the **trivial-multiscale check** (§C8/T1) first —
it is the cheapest single experiment that could collapse the whole edifice — and
replace the FC-matrix nulls with a **timeseries-generated null** (§3) before any
further claim.

---

## 0. The meta-risk: garden of forking paths

We have, across the era, selected post-hoc among: **6 bands** × **~16 scales** ×
**≥4 measures** (raw FC, LRG cophenetic, Grassmann, per-node) × **≥4 nulls**
(within-baseline, matched-strength, drift-windowed, drift-detrend) × **≥3 backbones**
× **multiple τ conventions**. Every "result" (β@OFC, α@mesoscale, high_γ@backbone)
emerged from a *different* cell of that grid, and each was later softened by a control.

- The **effective multiple-comparison burden is enormous and uncontrolled.**
- "A different story each day" is not bad luck; it is the **expected signature of an
  underpowered search over a large configuration space** with n=10.
- **No single analysis path was pre-specified.** Until one is, every p-value is
  conditional on a search we did not correct for.

**T0 (meta):** write down ONE primary analysis (one FC, one band-treatment, one
scale rule, one measure, one null) as a frozen pre-registration; run everything else
as explicitly-labelled secondary/exploratory. Severity: **CRITICAL / process.**

---

## 1. Pipeline map (where failure can enter)

```
raw sEEG timeseries
   │  C1  FC estimation (Welch cross-spectrum → imcoh_abs = <|ImCoh|>_f)
   ▼
per-band static FC matrix   ── C2  frequency band splitting (δ..γ_high, band-average)
   │  C3  sparsification (percolation backbone)
   ▼
sparse weighted graph
   │  C4  LRG diffusion (K=e^{-τL}) → cophenetic UPGMA distances D(τ)
   ▼
cophenetic distance vectors per phase (A,B,task,post)
   │  C5  trace measure ρ_sym(s) = split-half Spearman of reorganisations
   ▼
ρ_sym(s) per (patient,band)
   │  C6  nulls (matched-strength, drift)  +  C7 cohort stats (Wilcoxon,LOO)
   ▼
verdict  ── C8  the multiscale claim (is any of this beyond raw FC / beyond generic sparsification?)
```

Each Cn below: **ASSUMPTION → HOW IT FAILS → CHEAPEST TEST → SEVERITY.**

---

## 2. Failure directions by component

### C1 — FC protocol: `imcoh_abs = <|ImCoh|>_f`
**Assumption.** Band-averaged magnitude of imaginary coherency is a valid,
volume-conduction-immune FC that carries task-relevant coupling.
**How it fails.**
- **Static.** Averaged over the *entire* recording → all dynamics discarded; a
  "trace" is a difference of two static snapshots. Task coupling may be transient.
- **`|·|` discards sign** (lead/lag direction) and, band-averaged, conflates
  phase-lag consistency; a broadband `<|ImCoh|>` can be high with inconsistent
  per-bin phase.
- **ImCoh discards zero-lag coupling** — genuine in-phase neural coupling is
  invisible by construction. Some task reorganisation may be exactly there.
- **Strength-dominated.** `|ImCoh|` magnitude tracks amplitude/SNR per contact →
  the whole matrix is strength-driven, which is *why* matched-strength erases so
  much. The measure may be mostly amplitude.
- **Welch/`nperseg`** and the Pat_03 1024 Hz path change magnitudes; cross-patient
  comparability of `|ImCoh|` is assumed, not verified.
**Cheapest test.** Recompute the whole trace with **wPLI / dwPLI** (debiased,
also VC-immune) and with **lagged coherence**; if the β/trace story does not
reproduce, it is an `imcoh_abs`-specific artifact. Also: **time-resolved FC** to
test the static assumption.
**Severity: HIGH.**

### C2 — Frequency splitting (fixed canonical bands)
**Assumption.** δ/θ/α/β/γ_low/γ_high are the right decomposition and
"band-selectivity" (β special) is meaningful.
**How it fails.**
- **Arbitrary boundaries.** β≡13–30 Hz. The entire headline is *conditional on this
  cut*. We have **never tested boundary sensitivity**; a 15–25 or split-β could move
  or dissolve "β-selectivity".
- **γ_high = 80–300 Hz is 220 Hz wide** — averaging `|ImCoh|` across it mixes
  high-γ, ripples/HFO, and EMG/noise regimes; "high_γ trace" may be a binning
  artifact (it *only* appeared on the backbone).
- **Band-averaging is a strong prior**, not a measurement. The real spectral
  structure of coupling is imposed away.
- Feeds the forking-path multiplicity (6 bands).
**Cheapest test.** **Band-edge jitter sweep** (±2–3 Hz on each edge; does β survive?)
+ **frequency-resolved trace** (per-bin ρ_sym: is there a *frequency* where it peaks,
or is it a band-average artifact?) + **narrow γ_high**. Optionally data-driven
factorisation (NMF on the cross-spectrum).
**Severity: HIGH** ("band-selectivity IS the contribution" — if it is a binning
artifact, the contribution evaporates).

### C3 — Sparsification (percolation backbone)
**Assumption.** The percolation backbone (edges ≥ θ*, θ* = weakest max-spanning-tree
edge) is a parameter-free, principled sparsification that reveals *true* multiscale
structure the full matrix hides.
**How it fails.**
- **Post-hoc rescue.** Introduced *after* audit_174 showed the full matrix is
  diffusion-degenerate. It exists to manufacture multiscale-ness. The multiscale
  may be a property of **the sparsification pattern, not the brain**.
- **One edge sets the density.** θ* is fixed by the single weakest bottleneck edge;
  D0 shows density swinging **2–70%** across patients — wildly unstable, patient-set
  by a noisy edge.
- **Regime-dependent.** audit_175: near-tree fails diffusion, cycle-rich moderate
  works. The result *requires* landing in a density window; Pat_14 (2.5%) and Pat_15
  (43%) are outside it.
- **Hard threshold on a noisy static FC** → discards weak-but-real edges; may keep
  **same-probe** edges (CLAUDE.md invariant 5: 2–8× trivially high MSC/coupling on
  one shaft). Backbone community structure could be probe geometry.
- **Method-fragile.** Trace robustness to TMFG / k-NN / threshold / MST-union
  backbones untested.
**Cheapest test (decisive).** **Sparsify a matched-strength surrogate and a random
FC the same way** — do THEY show the same C(τ) multiscale ladder and comparable
ρ_sym? If yes, multiscale-ness is a **generic consequence of sparsifying any dense
matrix** and carries no brain-specific signal. Also: **same-probe fraction** of
backbone edges; **sparsifier-swap** robustness.
**Severity: CRITICAL** (load-bearing, post-hoc, patient-unstable).

### C4 — LRG diffusion + cophenetic
**Assumption.** The heat-kernel cophenetic distance captures multiscale community
structure *beyond* what raw FC shows.
**How it fails.**
- audit_174: at τ_min on the full graph **D=1/K ≈ D=1/A** — LRG adds nothing over
  raw FC. audit_176: **raw −log A reproduces the whole portrait** (β.019, α.0098).
  LRG's value-add is *unproven* except through the (post-hoc) sparsification.
- **UPGMA/average-linkage** cophenetic distance is linkage-dependent and known to
  distort; the "tree" is a modeling artifact of the clustering, not intrinsic.
- **τ / s-grid** is a chosen axis; "multiscale trace" risks being reparametrisation.
**Cheapest test.** At the **cohort trace level**, does **−log A** (raw, no diffusion)
reproduce every LRG verdict? (audit_176 says largely yes on the substrate — confirm
for the cross-phase trace.) If nothing survives that LRG shows and raw-FC+clustering
does not, **the "LRG framework" is decoration** and should be dropped from the
headline.
**Severity: HIGH** (method is oversold if LRG ≈ raw FC).

### C5 — Trace measure ρ_sym (shared-baseline positivity)
**Assumption.** ρ_sym measures task-induced-**and**-persistent reorganisation.
**How it fails.**
- **Shared-baseline bias.** ρ_sym correlates `(D_task − D_pre)` with
  `(D_post − D_pre)`; the shared `−D_pre` term induces **positive correlation even
  with no task** (documented in `node_localization.py`). The symmetric split-half
  form mitigates the *arbitrary-half* problem but **not** the shared-baseline
  positivity. This bias may *be* the "trace".
- It is a **Spearman on cophenetic distances** — heavily tied, structured values.
- A correlation, not an interpretable effect size.
- **This is why the drift null is so competitive:** the drift arc shares the exact
  same `−D_pre`-style structure, so it reproduces the positivity. The "trace" and
  the "drift floor" may both be the shared-baseline artifact.
**Cheapest test.** Quantify the shared-baseline positivity with **fully independent
baselines** (three disjoint rest_pre thirds for the two differences); measure the
null ρ_sym from baseline structure alone. If it equals the observed, the trace is
the bias.
**Severity: HIGH** (may be the whole effect).

### C6 — Nulls (matched-strength, drift): needed? sufficient?
**Assumption.** Matched-strength (strength-preserving FC shuffle) + drift (temporal)
are the right and sufficient controls.
**How it fails.**
- **Both are FC-matrix-level**, downstream of the signal. Neither preserves the
  actual timeseries statistics. They answer "given this matrix, …", not "given this
  brain signal, would this FC arise?"
- **MS (4-cycle ±δ)** preserves node strength and destroys everything else — it may
  be **too weak** (any non-strength structure clears) and tests a narrow hypothesis.
- **Drift** is **SNR-confounded** (this session): full-obs-vs-short-drift says trace,
  SNR-matched says drift; inconclusive. The detrend variant (§current work) is the
  cleaner form.
- **The gold-standard connectivity null — generated from the timeseries — has never
  been run** (§3). *This is the user's central point, and it is correct.*
**Are they needed?** MS is needed IF the claim is "structure beyond node strength"
(it is). Drift is needed because the design is temporally ordered
(pre<task<post). **But neither should be the PRIMARY null** — a timeseries-level
surrogate should, with MS/drift as secondary decompositions of *why* it clears/fails.
**Cheapest test.** Build §3 and compare.
**Severity: CRITICAL** (wrong null ⇒ every verdict uninterpretable).

### C7 — Statistics / cohort
**Assumption.** n=10 Wilcoxon + LOO is adequate.
**How it fails.** n=10, large per-patient heterogeneity (β 6/2/2), LOO borderline in
places; **forking-path multiplicity uncontrolled**; instability across iterations =
overfitting signature.
**Cheapest test.** Pre-register one path (T0); estimate forking multiplicity; simple
power analysis for the β effect size at n=10.
**Severity: HIGH** (the instability is the symptom of this).

### C8 — The multiscale claim itself (foundational)
**Assumption.** There is genuine multiscale structure carrying the trace.
**How it fails.**
- Full matrix: **not multiscale** (single C(τ) peak; Wigner). Multiscale appears
  **only after sparsification** → possibly manufactured.
- On the backbone the trace's **coarse-scale onset coincides with where the MS null
  collapses** — "multiscale trace" may just be "the null goes weak at coarse scales".
- **Villegas grounding (our own memory):** high-⟨κ⟩ dense → Wigner → single peak
  (= our full matrix); sparsify → off-semicircle → multiscale. So multiscale-ness is
  the **generic expected consequence of sparsifying any sufficiently dense weighted
  matrix** — brain or surrogate.
**Cheapest test (T1, do this FIRST).** Sparsify a **matched-strength surrogate** and
a **random** weighted matrix identically; compute C(τ) peaks and the ρ_sym arc. If
the surrogate backbone is *also* multiscale with comparable ρ_sym, **the multiscale
trace is trivial** and the project must pivot (to the epilepsy marker, or to a raw-FC
framing, or to a genuinely dynamic measure).
**Severity: CRITICAL / foundational.**

---

## 3. The null we actually want — generated from the timeseries

The user is right: the strongest null is built from the raw signal, not the FC matrix.

**Primary proposal — multivariate phase-randomisation (Prichard & Theiler 1994).**
Per phase, per patient: FFT every channel; apply an **independent** uniform random
phase rotation to each channel at each frequency (destroys ALL cross-channel phase
structure, i.e. kills coherence, while **exactly preserving each channel's power
spectrum**); inverse-FFT; run the **entire pipeline end-to-end** (imcoh_abs →
backbone → LRG → ρ_sym). Repeat R times → the honest null distribution of the trace.
- Answers: *does the observed FC/trace exceed what channels with identical spectra
  but no coupling produce?*

**Sharper variant — imaginary-part null.** Randomise phases with a **common** rotation
per frequency across channels (preserves the full coherence *magnitude*, destroys the
*imaginary* (lagged) part) → tests specifically whether `imcoh_abs` exceeds what the
zero-lag coherence structure implies. This isolates what ImCoh is *supposed* to add.

**IAAFT** if we also want to preserve the amplitude distribution (guards against
non-Gaussianity driving the result).

**Caveats to state up front:** phase-randomisation assumes (weak) stationarity within
phase and a linear-Gaussian null; it is the right null for "is there coupling beyond
spectrum", not for "is the coupling task-specific" (that still needs the cross-phase
contrast). So the end-to-end design is: **timeseries surrogate → same cross-phase
ρ_sym contrast** → does the *contrast* exceed the surrogate contrast. Cost: R × 3
phases × full pipeline; optimise with the per-window-C-bank trick already in
`08_drift_null_fair.py`.

---

## 4. Priority plan for the next agent (cheapest-falsification first)

| # | Lane | Question it kills-or-confirms | Cost |
|---|------|------------------------------|------|
| **T1** | **Trivial-multiscale** (C8/C3) | Is multiscale-ness generic to sparsifying any dense matrix? Sparsify MS-surrogate + random FC, compare C(τ) & ρ_sym. **If they match → pivot.** | low, ~30 min |
| **T2** | **Timeseries null** (C6/§3) | Does the trace beat a phase-randomised end-to-end surrogate? Replaces MS/drift as primary. | med, ~1–2 h |
| **T3** | **LRG-vs-rawFC necessity** (C4) | Does −log A reproduce every trace verdict? If yes, drop LRG framing. | low |
| **T4** | **Frequency binning** (C2) | Band-edge jitter + per-bin ρ_sym + narrow γ_high. Is β-selectivity a cut artifact? | med |
| **T5** | **Shared-baseline bias** (C5) | Independent-thirds baseline null ρ_sym. Is the trace the bias? | low |
| **T6** | **FC alternative** (C1) | wPLI/dwPLI cross-check of the β/trace story. | med |
| **T0** | **Pre-register one path** (meta) | Freeze one analysis; label the rest exploratory. | process |

Run **T1 first** — it is the single cheapest experiment that could end the project's
central claim. T2 is the null we should have had from the start.

---

## 5. What is currently defensible (honest inventory)

- **DEFENSIBLE — the epilepsy marker (D4).** Band-selective, matches/beats prior AUC
  (best-scale δ.88/γl.84/β.81), sound across C1–C5 controls, and it is a
  **within-phase spatial marker, not a fragile cross-phase trace**. This is the one
  robust thread and a viable paper on its own. It does **not** depend on the
  multiscale/trace edifice.
- **WEAK HINT — β.** β is the *consistently best* band across every trace analysis
  (MS scale-max, drift per-scale, and the 1-patient detrend where β alone survives).
  A persistent hint, **never a clean cohort win**. Worth one clean confirmatory test
  (T2), not a headline.
- **NOT DEFENSIBLE (currently)** — "LRG multiscale band-specific task trace" as a
  result; the multiscale framing; any specific band/scale/null cell as THE finding;
  the cophenetic/Grassmann apparatus as necessary over raw FC.

**Bottom line for the writing:** do not regenerate the §1 preprint trace figures on
either the full-graph or the backbone pipeline until **T1 + T2** return. The epilepsy
section (§3) can proceed independently. If T1 shows trivial multiscale and T2 shows no
timeseries-level trace, the honest paper is **"an epilepsy-localisation marker + a
negative/limits result on the cross-phase trace,"** not a multiscale-trace discovery.

---

## Live thread (do not lose)
- `09_drift_detrend.py` (full-coverage linear-detrend drift control) is **running**
  cohort-wide. 1-patient (Pat_02): **β retained 0.79 (survives), δ/θ/α/high_γ → negative
  (were drift)** — promising for β, awaiting n=10 cohort Wilcoxon. This is the cleanest
  drift control we have and feeds C6; it does **not** by itself rescue C3/C8.
