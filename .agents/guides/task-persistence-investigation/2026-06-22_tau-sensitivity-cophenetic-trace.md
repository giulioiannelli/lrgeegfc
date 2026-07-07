---
name: tau-sensitivity-cophenetic-trace
type: scope
era: IMCOH_ABS × COHORT_N10
status: active
created: 2026-06-22
updated: 2026-06-22
pointers:
  - .agents/guides/02_methods/lrg-framework-guide.md
  - scripts/02_preprint/preprint_06_beta_tau_sweep_raw_D.py
  - scripts/01_compute/audit/audit_85_epi_propagator_recovery.py
  - src/lrg_eegfc/utils/metrics/cross_phase.py
  - memory/lrg_tau_choice.md
---

> **Head.** Every cophenetic task-trace verdict in this project — β→OFC,
> the per-band `T_d`/`ρ_split` cohort trace, the whole §5.3 edifice — was
> computed at a **single** diffusion time `τ = 1/λ_max`, the *finest*
> (fastest) scale. That value is borrowed from the LRG papers, where it is
> "the right place to start" for **sparse topological networks with discrete
> spectra**. Our graphs are **fully connected, weight-heterogeneous, with a
> near-continuous spectrum** — the regime where the paper's scale-selection
> machinery is explicitly stated not to apply (guide §5.2, §6). At `τ = 1/λ_max`
> the heat kernel has barely diffused, so `D_coph(1/λ_max)` is approximately a
> monotone reweighting of the **raw FC** — the multi-step / mesoscale structure
> that is LRG's actual value-add only develops at **larger τ**. The SOZ marker
> already demonstrated, *in this exact dataset*, that τ matters and that the
> informative regime is **slower diffusion** (δ distant-SOZ lives at
> `τ ≈ 4–10/λ_max`, not at `1/λ_max`; audit_85/101/102). This scope defines a
> τ-resolved diagnostic of the cohort trace to settle whether `τ = 1/λ_max` is
> **robust**, **special**, or **suboptimal** — discharging guide open-question
> #1, which to date is flagged on every cohort claim as "sensitivity to τ not
> exhaustively tested".

---

## 0. Five-point critical preamble (mandatory; written before any code)

**(1) The claim under interrogation.** Not a new positive result. The object
is the *robustness* of the existing cophenetic trace verdicts to the choice of
diffusion time τ. Concretely: "the per-band cohort trace (`T_d^coph`,
`ρ_split^coph`) and its band ranking are a property of the LRG geometry, not an
artifact of evaluating it at the single finest scale `τ = 1/λ_max`."

**(2) The reference / null.** The reference is the trace statistic's value **at
`τ = 1/λ_max`** (the current canonical, which the diagnostic must reproduce
bit-exactly — correctness anchor). The diagnostic measures the **observed**
trace as a function of τ over the meaningful spectral window. The question is
binary at first order: is the τ-response **flat** (τ-irrelevant → current
verdicts robust) or **structured** (τ matters → verdicts are scale-specific)?
This is an **observed-statistic** sweep. **No matched-strength surrogate null is
run here** — that is mandatory only for a *new positive claim* at some τ ≠
1/λ_max, and is deferred to the follow-up that such a claim would trigger.

**(3) Strongest alternative the diagnostic must control for.** That any apparent
τ-dependence is a **degeneracy artifact, not signal**. As `τ → ∞` the heat
kernel collapses to the stationary uniform mode (`e^{−τL̂} → (1/N)·𝟙𝟙ᵀ`), so
`ρ(τ) → (1/N)` everywhere, `D_coph(τ) → const`, and `ρ^coph` is computed on a
near-zero-variance vector → numerically unstable. A trace that "moves" purely
because the geometry is dissolving into noise is **not** evidence that τ
matters scientifically. The mirror alternative: at the fine end the result must
**equal** the locked cache, or the sweep machinery is simply wrong.

**(4) Whether the diagnostic controls it — by mechanism.**
- *Collapse:* the per-phase **self-similarity** `G_self(τ) = ρ^coph(D_coph(τ),
  D_coph(1/λ_max))` and the per-phase cophenetic-distance variance are tracked
  at every τ; cells whose `D_coph(τ)` variance falls below a floor are
  NaN-masked, and the collapse onset is reported explicitly (not silently
  truncated). The meaningful window is bounded **at the coarse end by the
  Fiedler time `1/λ_gap`** (beyond which only the bipartition mode survives);
  we sweep to `K·/λ_gap` with `K = 4` precisely to *observe* the collapse, then
  exclude it from any verdict.
- *Anchor:* the α≡τλ_max = 1 point reconstructs `D_coph` from cached eigenpairs
  and is **verified equal to the cached `ultrametric_matrix`** (max abs diff
  2.8 × 10⁻¹⁴, Spearman 1.000000 on Pat_02/β; asserted per-cell in code).
- *Cross-patient comparability:* each (patient, band) is swept on its **own**
  `[1/λ_max, K/λ_gap]` grid and cohort aggregation is by **tau-index**
  (idx 0 = finest for everyone, idx max = `K/λ_gap` for everyone), so the
  cohort compares like spectral positions, not absolute times. The absolute
  `α = τ·λ_max` is also recorded for the absolute-axis view.
- *What it cannot do:* it **cannot** certify that a larger trace at some τ ≠
  1/λ_max is *significant* — only a matched-strength null at that τ can. It
  shows only whether the observed signal is τ-flat or τ-structured, and where.

**(5) What would falsify / what remains.**
- *Robustness verdict (caveat discharged):* if the observed cohort trace and the
  band ranking are essentially **flat** across the finest→Fiedler window for
  every band, then `τ = 1/λ_max` is representative and the standing "sensitivity
  to τ not exhaustively tested" caveat is downgraded.
- *τ-matters verdict (re-analysis warranted):* if the β (or any) trace
  **systematically strengthens toward a particular τ**, or a band that is null
  at `1/λ_max` becomes positive at coarser τ (the SOZ-δ analogue), then
  `τ = 1/λ_max` is **not** representative and a τ-resolved re-analysis — *with*
  matched-strength nulls and per-node re-localization — is the next phase.
- *Remaining limitation regardless of outcome:* this diagnostic speaks only to
  the **cohort scalar** trace. The per-node localization (β→OFC) τ-sensitivity
  is a separate, heavier follow-up (rebuild the per-pair trace incidence +
  matched-strength localization at the τ this diagnostic flags).

---

## 1. Notation

Per (patient `p`, band `b`, phase `s ∈ {pre, task, post}` and baseline halves
`{preA, preB}`), the LRG cache stores the Laplacian eigenpairs
`(λ_i, v_i)_{i=0..N−1}` of `L̂ = D̂ − Â` on the giant component (verified
present: keys `eigenvalues`, `eigenvectors`). For any `τ > 0` the heat-kernel
density and communication distance reconstruct **exactly** from the eigenpairs,
with no expm and no FC reload:

```
ρ(τ)      = V diag(e^{−τλ}) Vᵀ / Σ_i e^{−τλ_i}          (V = [v_0 … v_{N−1}])
D(τ)_ij   = 1 / ρ(τ)_ij ,  max-symmetrised, zero diagonal   (Villegas 2025 Eq. 1)
D_coph(τ) = cophenet( UPGMA( squareform D(τ) ) )            (the load-bearing object)
```

`D_coph(1/λ_max)` ≡ the cached `ultrametric_matrix` (anchor, verified). Let
`u_s(τ) = triu_{k=1}( D_coph^s(τ) )` be the per-pair cophenetic vector.

---

## 2. Measures (functions of τ)

For each (p, b) over a τ-grid `τ ∈ geomspace(1/λ_max, K/λ_gap, M)` (`K=4`,
`M=30`; per-phase geometric-mean spectrum defines the grid, all phases
evaluated at the same τ):

**Primary — cohort trace, cophenetic (load-bearing):**
- Locked triangle: `T_d^coph(τ) = ρ^coph(u_task, u_post) − ρ^coph(u_pre, u_task)`
  (i.e. `d(pre,task) − d(task,post)` with `d = 1 − ρ^coph`; `T_d > 0 = TRACE`,
  locked sign convention `feedback_td_sign_convention`). Full phases.
- Split-baseline (validated): `ρ_split^coph(τ) = Spearman(u_task − u_preA,
  u_post − u_preB)` via `cross_phase` semantics. Requires halves cache.

**Secondary — pre-cophenetic (raw `D(τ)`, no UPGMA):** identical formulas on
`triu D(τ)` instead of `triu D_coph(τ)`. Connects to / extends preprint_06 and
isolates *what the UPGMA projection adds* to the τ-response.

**Geometry diagnostics (collapse + value-add):**
- `G_self^s(τ) = ρ^coph(u_s(τ), u_s(1/λ_max))` — how fast the cophenetic
  geometry drifts from the canonical finest scale (1.0 at α=1 by construction).
- `var[u_s(τ)]` — collapse monitor; NaN-mask `ρ^coph` when below floor.
- Fiedler marker `α_F = λ_max/λ_gap` per cell (where `τ = 1/λ_gap` falls on the
  α axis).

**Cohort aggregation** (per band, per tau-idx): median trace, `n_trace =
#{T_d>0}` (or `#{ρ_split>0}`), one-sided `wilcoxon(..., alternative='greater')`
p, median `α`, fraction of cells collapsed. **No acceptance gate is
pre-registered** (`feedback_no_pre_registered_acceptance`): compute, plot,
evaluate signal/noise post-hoc with the PI.

---

## 3. Properties

- **Exact at the anchor.** α=1 ⇒ `D_coph(τ)` = cache ⇒ `T_d^coph(1/λ_max)`,
  `ρ_split^coph(1/λ_max)` = the locked per-band verdicts (asserted in code).
- **Cheap.** Reconstruction is `O(N² M)` per cell from cached eigenpairs; UPGMA
  on ~100 nodes is sub-ms. Full 10×6×30×~4-phase sweep ≈ seconds–minute.
- **Monotone reference frame.** α<1 probes *finer-than-current* (toward the raw
  1-step limit); α∈(1, α_F) probes the genuine LRG mesoscale; α>α_F probes the
  post-Fiedler collapse (diagnostic only).

## 4. Caveats

- Observed-statistic only; **a positive claim at τ≠1/λ_max needs its own
  matched-strength null** (CLAUDE.md mandatory-null rule).
- The cophenetic projection is a nonlinear UPGMA map; `T_d^coph` need not move
  monotonically with `T_d^raw`. That divergence is itself informative, not a
  bug.
- Collapse regime ρ^coph values are unstable and must be excluded from verdicts.

## 5. Pseudocode

```
for (p, b):
    load eigenpairs for pre, task, post (+ preA, preB if halves present)
    λmax_g, λgap_g = geomean over phases
    grid = geomspace(1/λmax_g, K/λgap_g, M)
    for τ in grid:
        D_s     = reconstruct 1/ρ(τ) from eigenpairs           # per phase
        Dcoph_s = UPGMA-cophenetic(D_s)   [lrgsglib.core]       # per phase
        assert α≈1 ⇒ Dcoph_s == cached ultrametric (once)
        record T_d^{coph,raw}(τ), ρ_split^{coph,raw}(τ),
               G_self^s(τ), var[u_s(τ)], α=τ·λmax_g, α_F
aggregate cohort per (band, tau_idx);  write per-patient + cohort CSV;  plot
```

## 6. Connection to prior tools

- **guide §5.1/§5.3/§10.1 + `memory/lrg_tau_choice.md`** — names the three
  canonical τ (1/λ_max, 1/λ_gap, τ*) and flags τ-choice as open question #1.
  This diagnostic operationalises it for the *trace* (not for C(τ), which is
  uninformative in our continuous-spectrum case — guide §6, L2 removed).
- **preprint_06_beta_tau_sweep_raw_D.py** — swept τ on **raw** `D(τ)` with the
  old within-baseline ρ_split and CSV-only output. This scope generalises it to
  the **cophenetic** (load-bearing) object, the **locked** `T_d`, all bands,
  collapse handling, and a figure.
- **audit_85/101/102 (SOZ)** — the in-project precedent that τ matters; their
  grid `geomspace(1/λmax, 10/λmax, 6)` is a *narrower, absolute-α* window. We
  extend to the Fiedler-relative window and apply it to the trace, not the
  node marker.
- **`utils.metrics.cross_phase`** — canonical per-pair trace decomposition;
  reused verbatim, only `D_coph` is swapped for `D_coph(τ)`.

## 7. Open questions (post-diagnostic)

1. If τ matters: which τ? 1/λ_gap (Fiedler), an intermediate `α*`, or a band-
   specific value? Does the SOZ-style slow-diffusion regime light up a trace
   that is null at the finest scale?
2. Does the UPGMA cophenetic projection *amplify* or *attenuate* the τ-response
   vs raw `D(τ)`?
3. Localization τ-sensitivity: does β→OFC survive / sharpen / move at the τ
   this diagnostic flags? (Separate matched-strength localization follow-up.)
