---
name: intrinsic-scale-reorganization
type: scope
era: IMCOH_ABS × COHORT_N10
status: executed_mostly_negative
created: 2026-06-22
updated: 2026-06-22
status_note: executed first-pass (audit_142) — persistence framing FALSIFIED; see VERDICT banner.
pointers:
  - .agents/guides/02_methods/lrg-framework-guide.md
  - .agents/reports/2026-04-29_critical-post-mortem.md
  - data/audit/entropy_crossphase/README.md
  - .agents/guides/task-persistence-investigation/2026-06-12_cophenetic-consolidation-arc.md
  - .agents/guides/task-persistence-investigation/2026-06-22_tau-sensitivity-cophenetic-trace.md
  - src/lrg_eegfc/workflow/diagnostics.py
  - .agents/preprint/headlines/02_encoding_vs_inference.md
  - .agents/reports/2026-06-22_arc-mesoscale-inference-null.md
---

> **Head.** Every task-trace verdict in this project reads the **cophenetic
> geometry** `D_coph` — a *derived* per-pair object built from the propagator.
> But the Laplacian-Renormalization-Group framework has a **native** observable
> we have never compared across phases: the **von Neumann entropy** `S(τ)` and its
> **specific heat** `C(τ) = −dS/d log τ` — the *thermodynamics of the diffusion
> process*, whose peaks mark the network's **characteristic integration scales**
> (Villegas et al., Nat. Phys. 2023). **Intrinsic-Scale Reorganization (ISR)** asks
> a milestone question the cophenetic/Grassmann probes structurally cannot: does
> learning a transitive-inference hierarchy **reorganize the brain's intrinsic
> multiscale architecture** — shift its characteristic scale `τ*` and/or compress
> its **effective dimensionality** `d_eff = e^{S}` — and does that reorganization
> **persist into post-task rest**? It is **static** (so it sidesteps the ≥10 s
> window floor that foreclosed replay), it is **cache-ready** (`entropy_C`,
> `entropy_1_minus_S` are stored in every LRG result), and it is the one read-out
> **no naive-FC and no spectral-clustering pipeline possesses** — so a positive ISR
> is *simultaneously* the strongest method-superiority demonstration and a genuine
> neuroscience advance: a thermodynamic signature of **schema / cognitive-map
> consolidation**, the milestone N1–N3 currently only gesture at.

---

> ## ⛔ VERDICT (2026-06-22, `audit_142` first-pass) — PERSISTENCE FRAMING FALSIFIED
>
> The primary milestone claim — a **persisting** reorganization of intrinsic scale
> structure (`rest_post − rest_pre`) — is **cohort-NULL** (6–7/10, n.s. for every
> scale observable in every band; falsification criterion §0.5.i triggered). The
> probe also **independently re-derived the 2026-04-29 retirement**
> (`.agents/reports/2026-04-29_critical-post-mortem.md`): on this data the von
> Neumann entropy is **saturated** (`S/log N ≈ 0.994`, within-patient phase std
> ≈ 1e-3) and `C(τ)` is **single-peaked** in all 240 cells — the scalar
> thermodynamic observable is **genuinely low-information here.** ISR-as-persistence
> is a **dead end**; do not pursue it.
>
> **The one signal found is NOT a milestone:** `τ*` (the `C(τ)` peak) shifts
> *coarser* **during `task_test`** (10/10, α/β/γ; β median Δ`log10 τ*` = +0.058,
> Wilcoxon p≈0.0025) and **relaxes back at rest** — i.e. a *transient within-task*
> state, the **opposite** of a persisting trace. It fails three honesty checks:
> (a) **non-persisting** (irrelevant to the paper's thesis); (b) **band-NONSPECIFIC**
> (α/β/γ all fire — the red flag for generic task-engagement / arousal, exactly
> alternative §0.3a, not a cognitive signature); (c) **duration-confounded**
> (`task_test` > `task_learn`, the same confound the cophenetic inference arc had to
> fight) and **likely redundant** with the already-locked mesoscale inference result
> (`audit_103d`, N2.2 — inference favours the coarser/mesoscale). **Marginal hint at
> best; do not inflate into a headline.**
>
> **Methodological keeper (one Methods sentence):** entropy saturation means the
> task-relevant information lives in the **per-pair cophenetic geometry**, NOT the
> global spectrum — which *post-hoc justifies* the project's choice of `D_coph` over
> a global thermodynamic read-out. The §10 **directed-flow companion is untouched by
> this** (different observable). Scope retained per the keep-don't-delete rule;
> superseded as a persistence milestone.

---

## 0. Critical preamble (mandatory 5-point — design before code)

1. **Claim.** Learning the TI structure changes the brain's intrinsic multiscale
   organization — it shifts the characteristic integration scale `τ*` (peak of the
   specific heat `C(τ)`) and/or reduces the effective dimensionality `d_eff` of the
   diffusion geometry (entropy `S`) — and the change **persists into `rest_post`**
   (a consolidation signature, stronger for the *inference* demand and concentrated
   in the OFC/relational system).

2. **Null.** The **matched-strength surrogate** (C3): same per-node strengths,
   scrambled structure (`strength_preserving_shuffle`), recompute `S/C` from the
   surrogate eigenvalues, recompute the cross-phase shift. Plus the within-baseline
   split-half (C1) and drift-floor (C2) nulls for the persistence contrast.

3. **Strongest plausible alternatives the null must defeat.**
   (a) **Global arousal / vigilance drift** `rest_pre → rest_post` that rescales the
   whole spectrum nonspecifically (drowsiness flattens `C(τ)`). (b) A **trivial
   degree/strength change** (the spectrum's scale is set by node strengths).
   (c) **Electrode-count / spectrum-scale** differences (`τ*` depends on `N` and on
   `λ_max`). (d) **Redundancy with `ρ^coph`** — that ISR is just the locked
   cophenetic trace re-described, carrying no independent information.

4. **Does the null actually control them — by mechanism.**
   - (b) **Yes** — matched-strength preserves strengths and kills exactly the
     degree-driven spectrum scale; an ISR shift beyond it is structural, not degree.
   - (c) **Yes, by construction** — ISR is evaluated **within patient** (`N`, `λ_max`
     fixed across that patient's four phases); we *never* compare absolute `τ*`
     across patients, only within-patient cross-phase deltas. Surrogates are built
     per (patient, phase) so the same `N` is held.
   - (d) **Tested, not assumed** — we require ISR to fire where `ρ^coph` is silent
     (e.g. a Grassmann-only band) and/or to survive partialling out the per-patient
     `ρ^coph` magnitude. If ISR ≡ `ρ^coph`, it is not a new result.
   - (a) **NO — matched-strength does not control arousal.** This is the binding
     honesty point. A vigilance drift is a *real* structural change, so a
     strength-null will not remove it. ISR-as-consolidation therefore inherits the
     **whole paper's ceiling: there is no task-free control session**
     ([[feedback_no_opaque_matrix_nulls]] — separating task-specific change from a
     trait/state backbone has no clean null without one). The defense is *exactly*
     the N1/N2 defense, ported: **cognitive specificity within the design** —
     inference (`test−learn`) > encoding, **band-specificity** (it tracks the
     trace bands, not all bands), and **anatomical concentration** (OFC sub-graph).
     A pure arousal effect would be band-nonspecific, phase-monotone, and
     anatomically flat. State this ceiling in the first paragraph of any writeup.

5. **What would falsify it (and what stays unprovable).** Retire ISR as a milestone
   if **any** of: (i) no consistent cohort cross-phase *direction* in `S`/`C`
   (compute+plot says noise); (ii) the persistence shift fails matched-strength;
   (iii) ISR is fully explained by `ρ^coph` (no independent information, alt-d);
   (iv) it appears **equally** in a no-trace band (θ / γ_high) — that would mark it
   as a nonspecific spectrum drift, not a cognitive reorganization. **Permanently
   unprovable without new data:** that the shift is *caused by* consolidation rather
   than by any post-task state change (no control session, no behavior). The claim
   ceiling is "inference-specific, OFC-concentrated, matched-strength-surviving
   reorganization of the brain's intrinsic scale structure that persists offline" —
   never "consolidation" as mechanism, never a behavioral correlation.

---

## 1. Why this is the overlooked milestone (positioning, read this once)

N1 says *a* trace persists; N2 says it carries encoding vs inference content; N3
says the same operator flags epileptogenic tissue. All true — but all **descriptive
of displacement**: "something moved (in per-pair geometry) and stuck." A referee
asks *what it means for the brain's organization*. ISR answers at the level
neuroscience cares about: **integration vs segregation, effective dimensionality,
and characteristic scale** — the language of cognitive-map / schema formation,
neural-manifold compression, and complexity. It reframes the headline from
"connectivity reorganizes" to **"learning a relational structure reshapes the
brain's intrinsic multiscale architecture, and the reshaping consolidates
offline."** Crucially the observable (`S`, `C`) is the RG framework's *defining*
quantity and has **no analogue** in thresholded-FC graph theory or
spectral-clustering — so it converts the method-superiority argument from "our
probe dissociates from Grassmann" (internal) into "we measure a thermodynamic order
parameter of cognition that standard pipelines cannot even define" (external,
NN-shaped). This is the fusion of *method power* and *neuroscience milestone* the
PI asked for, and replay's failure (a *dynamic* milestone foreclosed by the α/β
timescale floor) does not touch it — ISR is **static**.

---

## 2. Notation

- `p ∈ P` patient (n=10); `b ∈ B` band; `φ ∈ Φ = {pre, learn, test, post}` phase.
- `L̂_{p,b,φ}` symmetric-normalized Laplacian of the `|ImCoh|` graph; eigenvalues
  `0 = λ_1 ≤ … ≤ λ_{N_p}`; `N_p` = electrode count (fixed across `φ` for a patient).
- Propagator `ρ̂(τ) = e^{−τL̂}/Z`, `Z = Tr e^{−τL̂} = Σ_i e^{−τλ_i}`.
- Boltzmann weights `π_i(τ) = e^{−τλ_i}/Z`.
- `S_φ(τ) = −Σ_i π_i log π_i` — von Neumann entropy. Range `[0, log N_p]`.
- `C_φ(τ) = −dS_φ/d log τ` — specific heat (a.k.a. spectral susceptibility).
- Diffusion-time grid `T = {τ_k}` log-spaced on `[1/λ_max, 1/λ_2]` padded
  `PAD_FACTOR` decades each side (the cached `entropy_tau`).

All three (`entropy_tau`, `entropy_1_minus_S ≡ 1 − S/log N`, `entropy_C`) are
**already cached** per (p, b, φ) in the LRG result (`load_lrg_result`), and
recomputable from any eigenvalues via `compute_entropy_curve` — including
matched-strength surrogate eigenvalues (`load_or_compute_surrogate_eigs`).

---

## 3. Definitions (the scalars ISR reads off the curves)

For each (p, b, φ) define three **within-patient-comparable** scalars:

- **Characteristic scale** `τ*_φ = argmax_{τ ∈ T_core} C_φ(τ)`, reported as
  `log10 τ*`. `T_core` excludes the padded grid extremes (peaks at the pad edge are
  artifacts, §6). Multi-peak networks: keep the prominence-ranked peak set
  `{τ*^{(1)} > τ*^{(2)} …}` via `compute_susceptibility_diagnostics`; the primary
  scalar tracks the dominant peak with within-patient peak-matching across `φ`.
- **Peak sharpness** `Cmax_φ = max_{τ ∈ T_core} C_φ(τ)` — how cleanly scales
  separate (a sharp peak = a well-defined mesoscale; a flat `C` = scale-free /
  drowsy).
- **Effective dimensionality** `deff_φ = exp(S_φ(τ_ref))` at the canonical
  `τ_ref = 1/λ_max` (finest scale; the project's locked τ). `deff` = effective
  number of active diffusion modes; a *drop* = compression onto a lower-dimensional,
  more integrated geometry. (Companion: area `∫ (1−S/log N) d log τ`, a
  scale-integrated complexity.)

**Cross-phase contrasts** (for any scalar `X ∈ {log10 τ*, Cmax, deff}`), using the
**locked sign convention `T > 0 = trace`** ([[feedback_td_sign_convention]]):

- **Persistence / trace** (task as anchor):
  `T^X_p = |X_pre − X_test| − |X_post − X_test|`
  (`> 0` ⇔ `rest_post` is closer to the task value than `rest_pre` was — the
  scalar analogue of the cophenetic triangle).
- **Encoding shift** `e^X_p = X_learn − X_pre`.
- **Inference-specific shift** `f^X_p = X_test − X_learn`.
- **Directional persistence** `Δ^X_p = X_post − X_pre` (sign-readable: did the brain
  move *up* or *down* in scale / dimensionality and stay there?).

**Cohort aggregation** (per b, X, contrast): median, **sign count `/10`** (lead with
this per the regularity-over-null rule [[feedback_regularity_over_bh_null]]),
`wilcoxon_z`, and the **matched-strength upper-tail** `p = mean(s_finite ≥ obs)`.
Inference-specificity is the contrast that matters most: `f^X` significant while a
no-trace band's `f^X` is null.

---

## 4. What ISR is (and is not) — properties

- **Global order parameter, not a localizer.** `S/C` are *spectral invariants*
  (basis-free, relabel-invariant) — a PLUS for cross-phase comparison: **no
  eigenvector alignment / Procrustes** is needed (unlike Grassmann, which needs
  subspace matching and caused an anatomy retraction). But the price is that ISR by
  itself says *that* the architecture's scale structure moved, **not which pairs**.
  Localization is recovered two ways: (i) pair ISR with the locked `ρ^coph`
  OFC-concentration (ISR = "the scale shifts", cophenet = "and it concentrates in
  OFC"); (ii) the **sub-graph extension** — compute `S/C` on a-priori system-induced
  subgraphs (OFC, cingulate, …) for a per-system characteristic scale (§11).
- **Within-patient only.** `τ*`, `deff` depend on `N_p` and the spectrum scale;
  cross-patient *absolute* values are meaningless. Only within-patient cross-phase
  contrasts enter the cohort statistic. (This is the identifiability guard.)
- **Detects:** changes in integration scale, dimensionality, and complexity of the
  diffusion geometry. **Cannot detect:** which contacts/pairs carry it (needs §11);
  changes that preserve the spectrum but rotate modes (that is exactly Grassmann's
  job — ISR and Grassmann are complementary; a `deff`-flat / subspace-rotated change
  shows on Grassmann, not ISR).
- **Complexity:** `O(|T|·N)` per curve, all cached — the observed pass is minutes;
  the null is `R=200 × 4 phases` surrogate eigendecompositions (reuse the cached
  surrogate-eig store from audit_66+).

---

## 5. Caveats & failure modes (each with a mitigation)

- **Arousal/vigilance confound** (the headline risk, §0.3a). *Mitigation:*
  inference>encoding specificity + band-specificity + OFC sub-graph concentration;
  honest no-control-session ceiling stated up front. A pure arousal drift is
  band-flat, phase-monotone, anatomically flat — the opposite of the predicted
  signature.
- **`C(τ)` is a derivative** (`−dS/d log τ`) → noise amplification. *Mitigation:*
  light smoothing before peak-finding; verify `τ*` stability under the smoothing
  kernel; prefer `deff` (an entropy *level*, not a derivative) as the robust primary
  if `C` peaks are jittery.
- **Multi-peak / peak-swap across phases.** A patient may have 2–3 `C` peaks that
  exchange dominance. *Mitigation:* within-patient peak-matching by `log10 τ`
  proximity, report the full peak set, and run the contrast on `deff` (peak-free) as
  a cross-check.
- **Pad-edge peaks.** `compute_entropy_curve` pads the grid; a spurious `argmax` can
  sit at the pad boundary. *Mitigation:* restrict `argmax` to `T_core` (inside
  `[1/λ_max, 1/λ_2]`), drop pad decades.
- **Finite-N spectrum noise** (small implants, e.g. low-coverage patients).
  *Mitigation:* within-patient contrasts (own baseline), per-patient reporting,
  bootstrap CI; never a single-patient "strong" tag ([[feedback_no_single_patient_p_driven]]).
- **Redundancy with `ρ^coph`** (§0.3d). *Mitigation:* the explicit independence
  test — ISR must fire in a band where `ρ^coph` is silent (γ_low/δ are
  Grassmann-only candidates) or survive partialling `ρ^coph`. If it cannot, report
  it honestly as "a global re-description of the cophenetic trace," not a new axis.
- **Never as `optimal_threshold`.** `τ*` here is a *characteristic scale of the
  entropy*, **not** `lrg.optimal_threshold` and **never** used as the diffusion τ
  for any downstream cophenetic computation (never/always list).

---

## 6. Pseudocode (algorithmic)

```
# ---- observed ----
for p in patients:
  for b in bands:
    for φ in {pre, learn, test, post}:
      R = load_lrg_result(p, φ, b, "imcoh_abs")
      τ, S, C = R.entropy_tau, 1 - R.entropy_1_minus_S*logN, R.entropy_C   # or recompute from R.eigenvalues
      Ccore = restrict(C, τ ∈ [1/λmax, 1/λ2])         # drop pad edges
      Csm   = smooth(Ccore)
      X[p,b,φ].tau_star = log10( argmax_τ Csm )
      X[p,b,φ].Cmax     = max Csm
      X[p,b,φ].deff     = exp( S at τ=1/λmax )
    # within-patient cross-phase contrasts
    for Xname in {tau_star, Cmax, deff}:
      T[p,b,Xname]  = |X.pre−X.test| − |X.post−X.test|     # trace (T>0)
      e[p,b,Xname]  = X.learn − X.pre                       # encoding
      f[p,b,Xname]  = X.test  − X.learn                     # inference-specific
      d[p,b,Xname]  = X.post  − X.pre                       # directional persistence

# ---- cohort (compute + PLOT first; evaluate signal/noise post-hoc) ----
for b, Xname, contrast in product:
  vals = [ contrast[p,b,Xname] for p in patients ]
  report median(vals), sign_count(vals), wilcoxon_z(vals)

# ---- matched-strength null (only AFTER observed signal exists) ----
for p, b:
  for r in 1..R=200:
    for φ: eig_s = load_or_compute_surrogate_eigs(p, b, φ, W, R, ...)[r]
           Xs[φ] = scalars( compute_entropy_curve(eig_s) )
    Ts[r] = contrast(Xs)                                    # surrogate cross-phase shift
  p_MS[p,b] = mean( Ts ≥ T_obs[p,b] )                       # upper tail
# cohort null: combine per-patient via the project's surrogate-p convention
```

---

## 7. Visualization spec

- **Money figure — phase-overlaid thermodynamics.** Per band: cohort-mean `S(τ)`
  (top) and `C(τ)` (bottom) vs `log10 τ`, one line per phase
  (pre / learn / test / post), z-scored *within patient* on a common `log10 τ` grid
  before averaging. *Reading rule:* **a rightward/leftward shift of the `C(τ)` peak
  pre→post = the characteristic integration scale moved and stuck**; **a downward
  shift of `S` = dimensionality compression**. β / α first; θ / γ_high as the flat
  negative controls.
- **Per-patient `τ*` slopegraph** across the four phases (10 thin lines + cohort
  median) — shows consolidation as a within-subject trajectory (pre → learn → test →
  **post stays displaced**), no behavior needed.
- **Band × observable cohort-trace forest** (`T^X` with matched-strength band),
  mirroring the `ρ^coph` verdict forest so ISR sits beside the locked probes.
- Honor plotting rules: PDF-only vector, `use_lrg_style()`, no `suptitle`, X-epi
  variant where localized, `plt.close(fig)`.

---

## 8. Connection to prior tools (what it complements, subsumes, does not replace)

| Prior tool | Reads | ISR vs it |
|---|---|---|
| `ρ^coph` (N1/N2) | per-pair cophenetic geometry of `D(τ_max)` | **Complementary axis.** Same `ρ̂(τ)`, but ISR reads the *global spectral thermodynamics* (all modes, all scales), not pairwise distance. Localizable cophenet says *where*; ISR says *what kind* (scale/dimensionality). Independence is **tested**, not assumed. |
| Grassmann `d_G(k)` | leading-k subspace rotation | **Complementary.** Grassmann needs alignment and reads a *finite* leading subspace; ISR is basis-free and reads the *whole* spectrum's entropy. A pure mode-rotation shows on Grassmann not ISR; a scale/dimensionality change shows on ISR not necessarily Grassmann. |
| Mesoscale inference bump (`audit_103d`) | `T_infspec·e` peaking at τ≈2.6–6.8 | **ISR is its global/order-parameter form.** That result already says inference favors a mesoscale; ISR predicts a *corresponding `C(τ)` characteristic-scale shift* — the same physics read as a thermodynamic peak. Strong prior; direct tie-in. |
| τ-sensitivity scope (`2026-06-22_tau-sensitivity…`) | `ρ_split` vs τ | ISR uses the *whole* `S(τ)`/`C(τ)` curve as the object, not a τ-swept distance; the τ-sweep showed the trace is fine-scale + collapse-artifact-prone past Fiedler — ISR's `T_core` restriction respects that. |

**Never-list tension (declared):** the rule "*never start a new scalar hypothesis
test when the signal is visible in existing VI(k) / partition-multiscale / H2d
artifacts*." **Resolution:** ISR is a **different observable** (spectral
thermodynamics) that is **not present** in those artifacts (they read partitions /
distances), and it is run as **compute+plot+adjudicate** (no pre-registered gate,
[[feedback_no_pre_registered_acceptance]]), not a scalar FDR gate. It is exempt, but
the exemption is *conditional on the independence test* (§0.3d): if ISR turns out to
be `ρ^coph` re-described, it collapses back into the existing trace and the rule
would have been right.

---

## 9. Implementation plan

- **First-pass probe — `audit_142_entropy_crossphase.py` (RUNNING 2026-06-22).**
  Observed-only: extract `τ*`, `Cmax`, `deff` per (p, b, φ); cohort cross-phase
  contrasts (sign count + Wilcoxon); phase-overlay figures; novelty grep. Output
  `data/audit/entropy_crossphase/`. *Goal: is there ANY cohort signal?* If null →
  stop (falsified, §0.5.i).
- **Full version (only if first pass shows signal):**
  1. **Matched-strength null** (C3) via `load_or_compute_surrogate_eigs` +
     `compute_entropy_curve` — the mandatory referee.
  2. **Independence vs `ρ^coph`** (§0.3d): does ISR fire in a Grassmann-only band /
     survive partialling `ρ^coph`?
  3. **Inference-specificity** (`f^X` vs encoding) + **band-specificity** (vs θ/γ_h
     negative controls) — the arousal defense.
  4. **OFC sub-graph ISR** (§11) — per-system characteristic scale for localization.
- **Library reuse (no new private helpers):** `load_lrg_result`,
  `compute_entropy_curve`, `compute_susceptibility_diagnostics` (diagnostics.py),
  `load_or_compute_surrogate_eigs` (surrogate/matched_strength.py), `wilcoxon_z`,
  `bh_fdr` (metrics/hypothesis.py), `use_lrg_style`. A within-patient
  **peak-matching** helper (`track_C_peaks`) is the only candidate new function;
  promote to `lrg_eegfc.utils.metrics` only on a 2nd caller (coding-rules).
- **Preprint home if positive:** this becomes the *milestone framing of N1* (turns
  "a trace persists" into "the brain's intrinsic architecture is reorganized and
  consolidated") **or** a 4th headline — PI decision, deferred until the null lands.

---

## 10. Secondary, higher-risk companion (flagged, not primary): directed-flow installation

The `|ImCoh|` substrate is the *magnitude* of a **signed** imaginary coherency —
i.e. it discards a **direction-of-flow** sign already in hand (the signed-ImCoh
Hodge decomposition into gradient / curl / harmonic flow was the surviving
epileptogenic curl marker, [[epi_seedless_signed_flow_curl_2026_06_22]]). Transitive
inference is intrinsically **directed** (A→B→C…). **High-risk milestone:** does
inference install a **gradient (potential) flow hierarchy** — a directed ordering of
information flow — that persists into `rest_post`, i.e. a *directed* cognitive-map
trace? *Why risky:* mapping a flow potential to "the learned order" has **no clean
identification without item/stimulus labels** (which we lack), so the claim ceiling
is "a directed-flow hierarchy reorganizes and persists," not "it encodes the trained
order." Park behind ISR; the Hodge machinery exists (epi curl), the cognition
application does not. One-line verdict path: gradient-component triangle `T` per band
under the magnitude-matched control that already retracted the replay directional
candidate (audit_140) — if it cannot beat `|A|`, it dies the same way.

---

## 11. Open questions (load-bearing, deferred)

- **Primary scalar:** `deff` (robust, derivative-free) vs `τ*` (more interpretable,
  noisier)? Decide post-first-pass from which carries cleaner cohort signal.
- **Localization:** per-system sub-graph `S/C` — does the OFC sub-graph carry the
  shift (ties ISR to the locked N1.3 anatomy)? Power floor at ≤5-patient systems
  (same as audit_112) — may only support the demeaned cohort form.
- **Encoding vs inference:** which drives the architectural shift? Prediction from
  audit_103d: inference favors the mesoscale ⇒ `f^{τ*}` should be the live cell.
- **Redundancy threshold:** how much `ρ^coph`-independent variance must ISR carry to
  count as a new axis vs a re-description? Set *after* seeing the partial-correlation.
- **Directed-flow companion (§10):** scope its own report only if ISR lands and the
  PI wants the directed extension.
