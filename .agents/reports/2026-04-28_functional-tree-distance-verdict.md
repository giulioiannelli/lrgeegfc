---
name: 2026-04-28_functional-tree-distance-verdict
type: handoff
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-28
updated: 2026-04-28
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-28_functional-tree-distance.md
  - .agents/reports/2026-04-28_residual-subspace-diagnostic.md
  - .agents/reports/2026-04-28_psi-tau-scan-verdict.md
  - data/cache/functional_tree_distance/summary_n10_imcoh_abs.csv
  - data/cache/functional_tree_distance_null/summary_n10_imcoh_abs.csv
  - data/outputs/figures/2026-04-28_functional_tree_distance_null.pdf
---

# Functional tree distance — n=10 verdict (negative)

**Renormalization head.** The functional-tree-distance framework I
proposed earlier today (τ-resolved L²-distance on the LRG flow:
rank-based `Δ_S` on `D = 1/ρ` and height-based `Δ_P` on `ρ` itself,
integrated over `[1/λ_max, τ\*]` against a within-phase split-half
noise floor) is the **third measure today** that fails the cohort
task-trace gate. At n=10 with the within-baseline noise floor, only
**Pat_06** clears the rank gate (3/6 bands), and only **5/10** cells
across the *entire cohort* reach the strict gate; the height view is
slightly more permissive (10/60 cells trace) but still nowhere near
the CLAUDE.md cohort threshold (≥ 8/10 patients per band). **Pat_15**
is a clean anti-trace dissenter at 4/6 bands. No band reaches even
3/10 trace patients in the rank view; the best band (`low_gamma`,
height view) is 3/10. **The method is shit. It cost a half-day of
session time and produced no cohort signal.** Adds to the eight
previous task-trace operationalisations that have failed at n=10.
H2c α/β under continuous controls remains the only path forward.

## What I proposed and built

A τ-indexed L² distance between LRG flows:

- For each `(p, b, φ)`: build `L^φ`, compute `ρ(τ) = e^{−τL}/Z(τ)` and
  `D(τ) = 1/ρ(τ)` on a log-spaced τ-grid `[1/λ_max, τ\* = argmax C(τ)]`.
- Two parallel views per phase pair:
  - **`δ_S(τ)` (rank, dendrogram-topology):** `1 − Spearman(triu(D_A), triu(D_B))`.
  - **`δ_P(τ)` (height, scale-magnitude):** `1 − Pearson(triu(ρ_A), triu(ρ_B))`.
  - Integrate as `Δ_X = (log τ_max − log τ_min)^{−1} ∫ δ_X(τ) d(log τ)`.
- **Within-baseline noise floor:** split each phase's timeseries at
  `T//2`, compute `|ImCoh|` per half, evaluate the same `Δ_S, Δ_P` on
  the half-Laplacians at the same anchor τ-grid. The 4 within-phase Δ
  values per `(p, b)` form the floor.
- **Verdict per cell:** `C = Δ(test, post) − Δ(pre, post)`.
  - `trace` ⇔ `C < 0` AND `|C| > within-phase median`.
  - `anti-trace` ⇔ `C > 0` AND `|C| > within-phase median`.
  - `noise` ⇔ `|C| ≤ within-phase median`.

Scope report: `.agents/guides/task-persistence-investigation/2026-04-28_functional-tree-distance.md`.

## Cost

- 144 obs cells (4 patients × 6 bands × 6 phase pairs) on the
  primitive set, plus 216 obs cells across the 6 added cohort patients
  (cached from existing imcoh / imcoh_lrg). Obs compute: ~2 min.
- 60 null cells (10 patients × 6 bands), each requiring Welch on
  full timeseries split into 2 halves across all 4 phases. Null
  compute: ~30 min for the n=10 pass.
- Library: `src/lrg_eegfc/utils/metrics/functional_tree_distance.py`
  (~330 lines, 5 exports, reuses lrgsglib + scipy primitives, no
  private re-implementations).
- Compute scripts:
  `scripts/01_compute/diagnostics/diag_functional_tree_distance.py` and
  `scripts/01_compute/diagnostics/diag_functional_tree_distance_null.py`.
- Figure script:
  `scripts/01_compute/figures_embedded/fig_functional_tree_distance_null.py`
  — 7-page PDF, contrast-vs-noise scatter + per-band 10-patient panels.

## What the data says (n=10, within-baseline gate)

### Rank view `Δ_S` (strict gate `|C_S| > within_S_med`)

```
band         Pat_02 Pat_03 Pat_05 Pat_06 Pat_07 Pat_08 Pat_10 Pat_13 Pat_14 Pat_15  (trace/anti)
δ            noise  noise  noise  TRACE  noise  noise  noise  noise  noise  noise   1/0
θ            noise  noise  noise  TRACE  noise  noise  noise  noise  noise  ANTI    1/1
α            noise  noise  noise  noise  noise  noise  noise  noise  noise  ANTI    0/1
β            noise  noise  noise  noise  noise  noise  noise  noise  noise  ANTI    0/1
γ_l          noise  noise  noise  noise  noise  noise  noise  noise  noise  ANTI    0/1
γ_h          noise  noise  noise  TRACE  noise  noise  noise  noise  noise  noise   1/0
```

Cohort total: **3 trace cells (all Pat_06), 4 anti-trace cells (all
Pat_15), 53 noise cells.**

### Height view `Δ_P` (more permissive)

```
band         Pat_02 Pat_03 Pat_05 Pat_06 Pat_07 Pat_08 Pat_10 Pat_13 Pat_14 Pat_15  (trace/anti)
δ            noise  TRACE  noise  TRACE  noise  noise  noise  noise  noise  noise   2/0
θ            noise  noise  noise  TRACE  TRACE  noise  noise  noise  noise  ANTI    2/1
α            noise  noise  noise  noise  noise  noise  noise  noise  noise  ANTI    0/1
β            TRACE  noise  noise  TRACE  noise  noise  noise  noise  noise  ANTI    2/1
γ_l          TRACE  noise  TRACE  TRACE  noise  noise  noise  noise  noise  ANTI    3/1
γ_h          noise  noise  noise  TRACE  ANTI   noise  noise  noise  noise  noise   1/1
```

Cohort total: **10 trace cells, 5 anti-trace cells, 45 noise cells.**
Pat_06 dominates (5/6 bands trace under height); Pat_15 is the clean
anti-trace dissenter (4/6 bands).

### Best-band cohort threshold

CLAUDE.md cohort threshold = ≥ 8/10 patients per band.

| view | best-band | count | passes? |
|:---|:---|:---:|:---:|
| rank `Δ_S` | δ, θ, γ_h | 1/10 | no |
| height `Δ_P` | γ_l | 3/10 | no |

## Why this failed (and where I spent the day mistakenly)

The framework is mathematically clean — Villegas 2025 Eq. 1 is the
right primitive (`D_ij(τ) = 1/ρ_ij(τ)` is the genuine
ultrametric-precursor distance), the two-view rank/height split is a
real degree of freedom, and the within-baseline split-half null is the
correct gate for "is the contrast above noise". None of those choices
are wrong individually.

What's wrong is that **the cohort-task-trace claim — "post is closer
to test than to pre on the τ-axis" — is empirically not there at the
geometric whole-tree level under our continuous-spectrum, fully-
connected, weight-heterogeneous FC graphs.** This is the same outcome
as:

- **MRL** (`2026-04-25_module-retention-landscape.md`) — discrete
  subtree retention via Jaccard: cohort null. Same data.
- **Residual subspace** (`2026-04-28_residual-subspace-diagnostic.md`) —
  alignment of leading eigenvectors of `D^task − D^pre` and
  `D^post − D^pre`: below phase-label-permutation null. Same data.
- **Ψ τ-scan** (`2026-04-28_psi-tau-scan-verdict.md`) — 84% deep /
  16% trivial-coarse / 0.7% mid: useless as partition selector. Same
  data.
- **k-cuts / h_rel / τ̃ / mm** (across the April scalar session) —
  artefact-driven or dispersion-driven; closed in
  `2026-04-24_post-mortem-scalar-session.md`.
- **FTD with within-baseline floor** (this report) — 5/60 strict-gate
  trace cells, no band ≥ 3/10. **Same data, same negative outcome.**

The session-2 obs-only heatmap on n=4 (Pat_02/03/06/13) showed
5/6 bands in the cohort-median trace direction — that *looked* like
signal because the n=4 set was biased toward Pat_03 and Pat_06, and
because no noise gate was applied. Once n=10 expands the
patient-pool and the noise floor gates the contrast, almost
everything collapses to noise. That collapse is the third reading
of "**we keep proposing geometric primitives that look promising at
small-n unfiltered, then collapse under proper controls**". It's
the post-mortem-scalar-session pattern, repeated on an L²
function-space primitive.

## What this is *not* — to keep the framing honest

- **Not in tension with H2c α/β survival.** H2c measures
  *directional alignment* of pair-wise distance shifts across the
  upper triangle: per-patient Spearman of `D^test − D^pre` against
  `D^post − D^pre`. FTD measures *closeness ranking* of trees
  themselves: is the post-tree closer to the test-tree than the
  pre-tree under a τ-resolved L² distance? These are different
  geometric statements; one can hold while the other doesn't. H2c
  α/β under continuous controls remains the only path with positive
  cohort evidence at n=10.
- **Not "Pat_06 is the answer".** A single-patient case study is not
  a cohort claim. Pat_06 deserves separate scrutiny (implant, lesion
  location, task engagement) but does not rescue the cohort claim.
- **Not a methodology problem.** The FTD framework is the right
  framework — Villegas 2025's communication distance evaluated on a
  τ-grid, with rank and height views, with a within-phase noise
  floor. There is no obvious cleaner formulation. The negative
  result is about the *data*, not the *measure*.

## What I'd do tomorrow

1. **Stop adding measures.** This is the eighth task-trace
   operationalisation to fail at n=10. The probability that the
   ninth one succeeds where eight have failed is low; the user
   already flagged this in the post-mortem-scalar-session and the
   `feedback_dont_rerun_scalar_tests.md` memory.
2. **Commit to H2c α/β under continuous controls** as the cohort
   deliverable. Surfacing memo:
   `.agents/reports/2026-04-24_h1-h4-vi-results.md` already has the
   numbers; the controls are in place; the writeup is the deliverable.
3. **Pat_06 case study as a per-patient companion** (not a cohort
   claim). FTD figures + clinical context. Defensible if framed as
   "we observe robust FTD trace in a single patient; the cohort does
   not show this".
4. **Pat_15 as a dissenter case** — 4/6 bands anti-trace under both
   views. Worth understanding (anatomy? clinical?) but not as a
   counterargument to H2c α/β (different geometric statement).

## Status

`current` (verdict). Companion scope:
`.agents/guides/task-persistence-investigation/2026-04-28_functional-tree-distance.md`
status `draft` — to be moved to `superseded` after this report is
read by the user, with this report as the successor pointer.

## Cross-references

- `.agents/reports/2026-04-28_psi-tau-scan-verdict.md` — Ψ-as-selector died this morning.
- `.agents/reports/2026-04-28_residual-subspace-diagnostic.md` — residual-subspace died this morning.
- `.agents/reports/2026-04-28_raw-fc-phase-distance-verdict.md` — sibling raw-FC verdict.
- `.agents/reports/2026-04-24_post-mortem-scalar-session.md` — what the project already concluded about scalar gates in April.
- `.agents/diary/2026-04-28.md` — full session log.
- `data/cache/functional_tree_distance/summary_n10_imcoh_abs.csv` — observed Δs.
- `data/cache/functional_tree_distance_null/summary_n10_imcoh_abs.csv` — within-phase nulls.
- `data/outputs/figures/2026-04-28_functional_tree_distance_null.pdf` — 7-page n=10 PDF.
