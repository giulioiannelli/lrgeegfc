---
name: linearity-and-higher-order-structure
type: guide
era: COHORT_N10
status: current
created: 2026-06-16
updated: 2026-06-16
pointers:
  - .agents/guides/02_methods/lrg-framework-guide.md
  - .agents/guides/02_methods/imcoh-guide.md
---

# What our pipeline does (and doesn't) capture: second-order ImCoh, linear diffusion, higher-order *graph* structure

## Renormalization head

`|ImCoh|` is a **second-order** (linear-dependence) connectivity measure — it is
built from the cross-spectrum and is *lagged*, not *nonlinear* (the two get
conflated). LRG is **linear dynamics** (the diffusion `ρ = e^{−τL}` solves a linear
equation, diagonalized by the Laplacian eigenbasis) but the **map from the pairwise
matrix `W` to the cophenetic geometry is nonlinear**, so it manufactures
**higher-order *graph* features** — multi-step paths, triangles, mesoscale
communities — out of the pairwise edges. Crucially those are higher-order **in the
network**, not **in the signal statistics**: by the data-processing inequality a
deterministic function of `W` cannot carry more information about the signals than
`W` (second-order) already does. So the pipeline *reveals/reorganizes* structure
latent in pairwise coupling; it does **not** recover nonlinear neural coupling that
`|ImCoh|` discarded. **Manuscript wording:** say "multiscale/mesoscale structure
emergent from pairwise connectivity", never "higher-order interactions" or
"nonlinear analysis".

This note exists because three distinct conflations keep recurring and each one, if
it leaks into the manuscript, invites a correct reviewer objection.

## 1. `|ImCoh|` is second-order — and "lagged" ≠ "nonlinear"

`ImCoh_jk(f) = Im S_jk(f) / √(S_jj S_kk)` is a function of the **cross-spectrum**
`S_jk(f) = ⟨X_j(f) X_k(f)*⟩` — a **second-order** statistic (products of *pairs* of
signals). It captures linear/Gaussian dependence, frequency-resolved.

- **What makes it special is *lag*, not nonlinearity.** Nolte (2004): keeping the
  *imaginary* part isolates **non-zero-phase-lag** (time-delayed) coupling and
  discards the zero-lag part dominated by **volume conduction** (instantaneous
  source spread). A *purely linear* system with a transmission delay produces strong
  `ImCoh`. **Lagged ≠ nonlinear.**
- **What it does not capture:** higher-order / nonlinear statistical coupling
  (phase–amplitude, cross-frequency, non-Gaussian dependence). Those need different
  tools (bispectrum, transfer entropy, mutual information).
- This is not a weakness peculiar to us — essentially all mainstream EEG/MEG FC
  measures are second-order (coherence, PLV, wPLI, PLI). `|ImCoh|` is the
  volume-conduction-immune member of that family.

## 2. LRG is linear *dynamics* but a nonlinear *map* of `W`

- **Linear dynamics (the user's correct point).** `ρ(τ) = e^{−τL}` solves the linear
  diffusion equation `∂_τ ρ = −L ρ`, `L = D − W`, diagonalized by the Laplacian
  eigenbasis `L = V Λ Vᵀ`. The dynamics/physics is linear. (We use a single
  `τ = 1/λ_max`; we do not iterate the full RG flow.)
- **Nonlinear map.** The composition `W ↦ L ↦ e^{−τL} ↦ (1/ρ) ↦ average-linkage ↦
  cophenetic` is a **nonlinear, non-smooth function of the entries of `W`**: the
  matrix exponential is `I − τL + τ²L²/2 − …` and hierarchical clustering is
  combinatorial (sorting/merging). "Nonlinear function of the input", **not**
  "nonlinear dynamics" — two different senses.

## 3. The nonlinear map *does* build higher-order **graph** features

This is the legitimate value of LRG and the user's sharpening is correct:

```
L² = (D − W)²  contains  W²,   (W²)_ik = Σ_j W_ij W_jk     # products of two edges = 2-step paths i→j→k
```

`e^{−τL}` sums such edge-products to all orders, so the cophenetic distance for a
pair `(i,j)` is a nonlinear function of the **entire** pairwise matrix, not just
`W_ij`. Two contacts with **no direct edge** can come out "close" if many indirect
paths connect them. Hence LRG exposes **multi-step / relational / mesoscale**
structure (paths, triangles, communities, multiscale geometry) that is invisible in
the raw matrix. *This is why one runs LRG instead of reading `W` directly.*

## 4. …but **not** higher-order **signal** statistics

The features are a deterministic function `σ(W)`. The Markov chain
`signals → W → σ(W)` gives, by the **data-processing inequality**,
`I(signals; σ(W)) ≤ I(signals; W)`. `W` is second-order, so `σ(W)` carries **only
second-order** information about the signals. LRG **reveals/reorganizes** structure
already latent in the pairwise coupling — like community detection or PCA exposing
structure without injecting new information — it does **not add** information about
nonlinear neural coupling that `|ImCoh|` discarded.

So both statements hold simultaneously, with no contradiction:

| | emergent from pairwise edges | recovered |
|:--|:--:|:--:|
| higher-order **network/graph** features (paths, communities, mesoscale) | **yes** | ✓ revealed |
| higher-order **statistical** coupling of the signals (3rd+ order) | — | **no** (DPI) |

## 5. What we can / cannot claim — and exact wording

- **Can claim:** multiscale reorganization of (lagged, second-order) connectivity;
  mesoscale/community structure emergent from pairwise coupling; relational,
  multi-step graph geometry.
- **Cannot claim:** that the pipeline captures nonlinear neural coupling, or
  "higher-order interactions" in the network-neuroscience sense (genuine multi-way
  statistical dependence — hypergraphs, simplicial complexes from multivariate
  information). It does not measure those; it builds higher-order *graph* features
  from *pairwise* data.
- **Bounded honest caveat:** coupling that is *purely* nonlinear, leaving **no**
  second-order footprint, is invisible to `|ImCoh|`, and nothing downstream recovers
  it. Pure-nonlinear-only coupling is rare (most coupling has a second-order shadow),
  so this is a standard, bounded limitation, not a hole.
- **Wording to use:** "multiscale / mesoscale structure emergent from pairwise
  connectivity." **Avoid:** "nonlinear analysis", "higher-order interactions",
  "captures nonlinear coupling" — each invites a correct reviewer objection.

## Connection

- LRG primitive, formulas, `e^{−τL}` / eigenbasis / ultrametric:
  [`lrg-framework-guide.md`](lrg-framework-guide.md).
- `|ImCoh|` definition, volume-conduction immunity, taxonomy:
  [`imcoh-guide.md`](imcoh-guide.md).
- Why a physical phase-randomized null does not help an `|ImCoh|` pipeline
  (multivariate phase randomization preserves the cross-spectrum → no-op; univariate
  destroys connectivity) and why matrix-level nulls exist:
  `scripts/archive/2026-06_opaque-matrix-nulls/POSTMORTEM.md`.
