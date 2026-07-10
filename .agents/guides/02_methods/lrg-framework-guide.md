---
name: lrg-framework-guide
type: guide
era: CROSS_ERA
status: current
created: 2026-04-29
updated: 2026-04-29
pointers:
  - .agents/guides/02_methods/imcoh-guide.md
  - .agents/guides/02_methods/probe-bias-guide.md
  - .agents/guides/02_methods/linearity-and-higher-order-structure.md
  - .agents/guides/02_methods/multiscale-notions-propagator-vs-grassmann.md
  - .agents/guides/task-persistence-investigation/archive/scalar-vi-era/2026-04-29_decision-rules.md
  - .agents/guides/task-persistence-investigation/archive/scalar-vi-era/2026-04-29_measure-correctness-audit.md
---

> **What this pipeline does/doesn't capture (read before writing "nonlinear" /
> "higher-order" in the manuscript):** `|ImCoh|` is second-order (lagged ≠
> nonlinear); LRG is *linear* diffusion (`e^{−τL}`, Laplacian eigenbasis) but a
> *nonlinear map* of `W` that builds higher-order **graph** features (multi-step
> paths, communities) — **not** higher-order signal statistics. See
> [`linearity-and-higher-order-structure.md`](linearity-and-higher-order-structure.md).

# LRG framework — single source of truth

**The Laplacian Renormalization Group (LRG) primitive applied to our
fully-connected weighted ImCoh-FC graphs. This guide is the agent
reference: it states the canonical formulas (Villegas 2023 *Nat Phys*
+ Villegas 2025 *PRR*), maps each symbol to the exact lrgsglib /
workflow function that implements it, calls out where our pipeline
deviates from the paper (UPGMA shortcut at fixed τ, normalisation
rescaling, no τ-sweep), and explicitly enumerates which LRG-paper
arguments do NOT apply to our outlier case (fully connected,
continuous spectrum, weight-heterogeneity-driven). When in doubt,
re-read this — DO NOT freelance LRG semantics.**

---

## 0. References

- Villegas, Gili, Caldarelli, Gabrielli (2023). *Laplacian
  renormalization group for heterogeneous networks*. **Nat Phys** 19,
  445–450. PDF: `/home/giulio/Documents/papers/lrgeegfc/villegas2023laplacian.pdf`.
  Defines: L̂, K̂(τ), ρ̂(τ), S(τ), C(τ); Kadanoff supernodes via
  `ζ̂_ij = Θ(ρ'_ij − 1)` where `ρ'_ij = ρ_ij / min(ρ_ii, ρ_jj)`.
- Villegas, Gabrielli, Poggialini, Gili (2025). *Multi-scale Laplacian
  community detection in heterogeneous networks*. **PRR** 7, 013065.
  PDF: `/home/giulio/Documents/papers/lrgeegfc/villegas2025multiscale.pdf`.
  Defines: communication distance `D_ij(τ) = (1−δ_ij)/K_ij(τ)` (Eq. 1),
  proves ultrametricity, prescribes UPGMA dendrogram, defines Ψ(n;τ)
  partition stability (Eq. 2) and local Ψ_L(τ); introduces
  metastable nodes via Sankey across τ.

---

## 1. Notation

- Graph `G = (V, E)` with `|V| = N` nodes; weighted, symmetric,
  non-negative adjacency `Â ∈ ℝ_{≥0}^{N×N}`. For our case
  `Â_ij = ⟨|ImCoh_ij(f)|⟩_{f ∈ band}` (see `imcoh-guide.md`).
- Degree matrix `D̂_ii = Σ_j Â_ij`.
- Fluid Laplacian `L̂ = D̂ − Â`. Symmetric PSD.
- Eigendecomposition `L̂ = Σ_i λ_i |λ_i⟩⟨λ_i|`, with
  `0 = λ_0 < λ_1 ≤ … ≤ λ_{N−1} = λ_max`.
- λ_gap (= λ_1, "Fiedler") = smallest *positive* eigenvalue.
- Diffusion time `τ ∈ ℝ_{>0}`.

---

## 2. Core formulas (with codebase mapping)

| # | Quantity | Formula | Codebase |
|---|---|---|---|
| 2.1 | **Network propagator** | `K̂(τ) = e^{−τL̂}` | `expm(-tau * L)` at `lrgsglib/utils/lrg/spectral.py:121` |
| 2.2 | **Partition function** | `Z(τ) = Tr(K̂(τ)) = Σ_i e^{−τλ_i}` | `den = np.trace(num)` at `spectral.py:122` |
| 2.3 | **Density operator** | `ρ̂(τ) = K̂(τ) / Z(τ)` (trace 1, PSD) | `rho = num / den` at `spectral.py:123` |
| 2.4 | **Default τ** | `τ' = 1/λ_max` (finest meaningful resolution) | `tau = 1/max(spectrum)` at `spectral.py:120` |
| 2.5 | **Communication distance** | `D_ij(τ) = (1 − δ_ij) / K_ij(τ) = (1 − δ_ij) / ρ_ij(τ)` (Villegas 2025 Eq. 1). Properties: ≥ 0; symmetric; **strong triangle inequality** ⇒ ULTRAMETRIC | `Trho = 1/rho`, then `np.maximum(Trho, Trho.T)`, diag set to 0 — at `spectral.py:124-126`. Stored as `LRGResult.ultrametric_matrix` (after UPGMA cophenetic projection) |
| 2.6 | **Normalized von Neumann entropy** | `S(τ) = −(1/log N) Σ_i μ_i(τ) log μ_i(τ)` where `μ_i(τ) = e^{−τλ_i}/Z(τ)` are eigenvalues of ρ̂. Range `[0, 1]` | `S[i] = -np.nansum(rho * np.log(rho)) / np.log(N)` at `infocomm.py:150` |
| 2.7 | **`1−S(τ)`** (stored field) | Monotonically increasing with τ; 0 at τ→0 (uniform), 1 at τ→∞ (single dominant mode) | `1 - S` first return of `entropy()` at `infocomm.py:160`; stored as `entropy_1_minus_S` |
| 2.8 | **Specific heat** | `C(τ) = −dS/dlog τ` | `dS = log(N) · np.diff(1-S) / np.diff(log(t))` at `infocomm.py:158`; stored as `entropy_C`. **Length is `len(τ) − 1`** (one fewer than `entropy_tau`) |
| 2.9 | **τ-grid** | 400 log-spaced points over `[10^{-3}, 10^5]` for `compute_lrg_analysis` defaults | `t = np.logspace(t1, t2, steps)` at `infocomm.py:137`; `entropy_t1=-3, entropy_t2=5, entropy_steps=400` in `workflow/lrg.py` |

---

## 3. Dendrogram construction

### 3.1 What the paper prescribes (Villegas 2025 §I)

1. Pick a τ.
2. Compute `D̄_ij(τ) = 1/K_ij(τ)` (or equivalently `1/ρ_ij(τ)`).
3. Order all node-pairs by `|ρ̂(τ)|` from highest to lowest.
4. Aggregate sequentially → produces a hierarchical tree whose merge
   distances are in `D̄` units.
5. Use **UPGMA / average-linkage** to balance outlier sensitivity vs
   cluster compactness (the paper's explicit choice).

### 3.2 What lrgsglib + workflow.lrg implement

Sequence in `workflow/lrg.py:339-353`:

```python
sm1, spec, *_, tau = entropy(giant, steps=400, t1=-3, t2=5)
_, _, _, Trho, _ = compute_laplacian_properties(giant)         # uses default τ = 1/λ_max
dists = squareform(Trho)                                       # condensed 1/ρ matrix
linkage_matrix, labels, _ = compute_normalized_linkage(dists, giant, method="average")
threshold, *_ = compute_optimal_threshold(linkage_matrix)
ultrametric_square = extract_ultrametric_matrix(linkage_matrix, n_nodes)
```

The dendrogram is **built once at τ = 1/λ_max** via UPGMA on the condensed `Trho`. This is the **paper's recipe at one canonical τ**, NOT a τ-scan.

### 3.3 The `compute_normalized_linkage` rescaling (visualization-only)

`clustering.py:55-68`:

```python
linkage_matrix = linkage(dists, method)              # initial UPGMA
max_distance = linkage_matrix[-1, 2]                 # root height
tmax = max_distance + 0.01 * max_distance            # = 1.01 · root
linkage_matrix = linkage(dists / tmax, method)       # re-link on rescaled dists
```

After rescaling, **every dendrogram has root at exactly `1/1.01 ≈ 0.9901`**.
This is **visualization-comparability only**: relative cophenetic
distances and partition assignments at any cut are unchanged. The
underlying merge ORDER and ratios are preserved; only the absolute
units change. **Effect on results: none.** Effect on cross-patient
comparability of fractional depth `h_rel = h(k)/dmax`: enables it
(otherwise dmax would vary cohort-wide).

### 3.4 Cophenetic ultrametric matrix

`extract_ultrametric_matrix` (`infocomm.py:22-44`):

```python
cophenetic_dists = cophenet(linkage_matrix)       # scipy
ultrametric_matrix = squareform(cophenetic_dists)
```

So **`LRGResult.ultrametric_matrix[i, j]` is the cophenetic
distance** = the merge height of the LCA of `i, j` in the (rescaled)
dendrogram. NOT a τ value. NOT a `1/ρ` value. It's the
*dendrogram-derived* ultrametric (which by UPGMA construction
satisfies the strong triangle inequality, and which numerically
*approximates* the original `1/ρ` distances after the average-linkage
projection).

---

## 4. Ψ partition stability (`optimal_threshold`)

### 4.1 Definition (Villegas 2025 Eq. 2)

For the dendrogram at fixed τ, with merge thresholds Δ ordered from
root downward (Δ_1 = root, Δ_n = nth merge from root):

```
Ψ(n; τ) = N · [log10 Δ_n(τ) − log10 Δ_{n+1}(τ)]
N        = [log10 Δ_1(τ) − log10 Δ_{n_max}(τ)]^{-1}      (normalisation)
```

Ψ ∈ [0, 1]. **The optimal cut is `n* = argmax_n Ψ(n; τ)`** — the
partition that survives the largest log-merge-distance gap.

### 4.2 Codebase implementation (`compute_optimal_threshold`)

`clustering.py:73-115` implements EXACTLY the Ψ formula:

```python
dendro_thresholds = linkage_matrix[:, 2]
D_values = dendro_thresholds[::-1]                            # root → leaves
N = 1 / (np.log10(D_values[0]) - np.log10(D_values[-1]))
sigma_i = N * (np.log10(D_values[i]) - np.log10(D_values[i+1]))   # = Ψ(i; τ)
optimal_branch_index = np.argmax(sigma)
optimal_threshold    = D_values[optimal_branch_index + 1]
FlatClusteringTh     = 0.9 * optimal_threshold                # what LRGResult stores
```

**`LRGResult.optimal_threshold` IS the Ψ-selected cut height** scaled
by `0.9`. Use it as `fcluster(Z, FlatClusteringTh, criterion="distance")`.

### 4.3 What `optimal_threshold` is NOT (anti-hallucination)

- **NOT a diffusion time τ.** Different space entirely (cophenetic
  merge distance in normalized `1/ρ` units, not time).
- **NOT a single canonical cluster count.** It's the height *below the
  most stable Ψ branch* — which corresponds to a cluster count
  determined by the tree's particular merge ordering.
- **NOT scale-invariant.** Specific to the patient/band/phase tree at
  τ = 1/λ_max.
- **NOT a partition quality score.** It's a position in the
  log-merge-distance axis. Use `Ψ(n; τ)` itself for quality.

The CLAUDE.md never-list rule "Never use `lrg.optimal_threshold` as a
diffusion time τ" enforces this anti-hallucination guard.

---

## 5. τ choice — open question for our case

### 5.1 Three canonical τ values (Villegas 2025 §I, §III)

| τ | Meaning | Codebase use |
|---|---|---|
| **τ' = 1/λ_max** | Finest resolution; no info integrated out yet | Default in `compute_laplacian_properties`; `audit_17_taumin_four_phase.py` |
| **τ = 1/λ_gap = 1/λ_1 (Fiedler)** | Defocuses microscopic detail; partition reflects bipartition modes | Not currently used; possible follow-up |
| **τ\* = argmax_τ C(τ)** at long times | Coarsest meaningful scale; beyond τ\* nothing new | `audit_18_taustar_four_phase.py` (re-runs the four-phase view at this τ) |

### 5.2 Why our case is hard

Villegas's BA / RR / real-PPI examples have **multiple discrete C(τ) peaks** that identify a *ladder* of characteristic mesoscale scales. Our ImCoh FC graphs have **continuous spectra** (close to Wigner-semicircle limit; see Villegas 2025 §V, Fig. 7d) — so:

- **C(τ) HAS a single, well-defined interior peak τ\*** — NOT smooth/monotonic/peakless.
  Verified empirically 2026-07-09: unimodal in every (patient, band) cell (30/30 +
  full cohort), prominence ≈ 1.0–1.4, matching the pipeline's own cached `entropy_C`
  to ~1e-15. τ\* sits at an intermediate scale of order **10/λ_max** (α ≡ τ\*·λ_max
  ≈ 8–15), just coarser than the Fiedler time 1/λ_gap (which is an *interior* scale
  of the window, at α ≈ 3, ~4× finer than τ\*).
- What is **absent** is a *family of multiple* gap-induced peaks. The single τ\* marks
  the **collapse scale** (where the propagator relaxes onto the trivial λ_0 mode and
  inter-contact structure is lost), not a nested mesoscale community ladder — so `Ψ`
  isolates no privileged interior cut, and spectrum-based *mesoscale-community*
  identification fails. Identification of the single collapse scale τ\* does **not** fail.
- τ\* ≈ the empirical **collapse onset** of audit_121 (α ≈ 10–15) — two independent
  diagnostics (spectral C-peak; cross-phase placebo/self-similarity) locating the
  same boundary. τ\* is therefore a valid coarse-τ **ceiling**: read at/below it, not past.

Consequently the LRG pipeline commits to a single τ = 1/λ_max for the **cross-phase
trace** (its multiscale content coming from the **dendrogram's own UPGMA hierarchy**,
not a τ-scan — and coarsening the *cross-phase* comparison toward τ\* manufactures a
collapse artifact, audit_121). A **within-phase** read-out (the epi SOZ marker) is
immune to that cross-phase collapse and is legitimately read near the ceiling
(τ_5 = 10/λ_max ≈ τ\*); see `overleaf/methods.tex` §`ssec:methods_epi`.

### 5.3 Open methodological question

Whether the right τ for the task-trace question is (a) 1/λ_max (current), (b) some intermediate τ tuned to the FC graph's effective spectral dimension, or (c) a τ-sweep producing a Sankey of partition-evolution → all UNANSWERED. Flag every cohort claim with "computed at τ = 1/λ_max; sensitivity to τ not exhaustively tested".

---

## 6. Our outlier case — what does NOT translate from the paper

| Paper assumes | Our case has | Implication |
|---|---|---|
| Topological structure (sparse adjacency, BA / RR / ER / real biological) | Fully connected weighted graph (every pair has non-zero edge weight) | Coarse-graining "by edge presence" (Kadanoff supernode rule via `ρ'_ij ≥ 1`) is degenerate — at any τ > 0 the threshold rule produces either everything-merged or nothing-merged |
| *Multiple* discrete C(τ) peaks → a mesoscale community ladder | Continuous spectrum (semicircle-like for high spectral dimension d_S) → a **single** well-defined susceptibility peak τ\* (≈10/λ_max, the collapse scale), no multi-peak ladder | Cannot read a *nested mesoscale ladder* off C(τ), and `1−S` / `C` **values** are not a cross-phase discriminator — **L2 (entropy-curve task-trace rung) permanently removed 2026-04-29** (comparing C-curves across phases is theatre). BUT the peak **location** τ\* is a real, well-defined scale marker (verified 2026-07-09) and a valid coarse-τ *ceiling*. Do **not** write "C(τ) is smooth / has no peak" — that is false; write "no *multiple* gap-induced peaks". |
| Mesoscale community detection at τ' < τ < τ\* | Single τ = 1/λ_max + dendrogram hierarchy | "Multiscale" in our work is **dendrogram-multiscale** (cuts at integer-k or h_rel), not τ-multiscale |
| Static network analysis | Cross-phase comparison (rest_pre, task_test, rest_post per patient) | Need per-patient within-baseline null (halves cache) — no equivalent in the paper |
| No anatomy bias | sEEG: contacts on the same probe carry trivially high coupling | Need probe-bias control (cross-probe / same-probe split per `probe-bias-guide.md`) |
| Scale-invariance / fractal / power-law degree | Degree distribution is uniform-ish across nodes (fully connected) | Cannot use "constant-C plateau" or "scale-invariant" framing |

**Practical rule:** cite Villegas 2023+2025 for the LRG primitive (`ρ(τ)`, communication distance, ultrametric, UPGMA dendrogram, Ψ optimal cut) — these are valid for any positive-weight graph. **Do NOT cite them as authority for "mesoscale communities at τ\*", "informational phase transitions", "scale-invariant networks", or "metastable bridge nodes detected by τ-sweep" — those arguments do not go through for us.**

---

## 7. Verification snippets (run these to confirm formulas vs codebase)

All snippets assume `lapbrain` env activated and ROOT = `/home/giulio/Documents/research/neural_networks/lrgeegfc`.

### 7.1 Check that `Trho = 1/ρ` element-wise

```python
import numpy as np
from scipy.linalg import expm
from lrgsglib.utils.lrg.spectral import compute_laplacian_properties
import networkx as nx

# Tiny graph
G = nx.path_graph(4)
spectrum, L, rho, Trho, tau = compute_laplacian_properties(G)
expected = expm(-tau * L) / np.trace(expm(-tau * L))
assert np.allclose(rho, expected), "ρ definition mismatch"
# Trho = 1/rho symmetrised, diag zero
T_expected = np.maximum(1.0/rho, (1.0/rho).T)
np.fill_diagonal(T_expected, 0)
assert np.allclose(Trho, T_expected), "Trho ≠ 1/ρ"
# τ default = 1/λ_max
assert abs(tau - 1.0/max(spectrum)) < 1e-12, "default τ ≠ 1/λ_max"
print("OK: ρ, Trho, default τ all match definitions")
```

### 7.2 Check that `optimal_threshold` is the Ψ-selected cut

```python
import numpy as np
from lrgsglib.utils.lrg.clustering import compute_optimal_threshold

# Build a fake linkage matrix with a clear Ψ-optimal gap
Z = np.array([
    [0, 1, 0.10, 2],
    [2, 3, 0.12, 2],
    [4, 5, 0.50, 4],   # large log-gap above 0.12 → Ψ peaks here
    [6, 7, 0.55, 6],
], dtype=float)  # placeholder; in practice from scipy.linkage

threshold, opt, sigma, idx = compute_optimal_threshold(Z, scaling_factor=0.9)
# Verify Ψ formula by hand
D = Z[:, 2][::-1]
N = 1 / (np.log10(D[0]) - np.log10(D[-1]))
sigma_hand = [N * (np.log10(D[i]) - np.log10(D[i+1])) for i in range(len(D)-1)]
assert np.allclose(sigma, sigma_hand), "Ψ formula mismatch"
assert idx == int(np.argmax(sigma_hand)), "argmax mismatch"
assert opt == D[idx+1], "optimal_threshold ≠ D_{argmax+1}"
assert threshold == 0.9 * opt, "FlatClusteringTh ≠ 0.9 · opt"
print("OK: compute_optimal_threshold implements Ψ exactly")
```

### 7.3 Check that `entropy_C` is `−dS/dlogτ`

```python
import numpy as np
import networkx as nx
from lrgsglib.utils.lrg.infocomm import entropy

G = nx.barabasi_albert_graph(50, 2, seed=0)
sm1, dS, varL, t = entropy(G, steps=300, t1=-2, t2=4)
# sm1 = 1 - S, length len(t)
# dS  = log(N) · diff(1-S) / diff(log t), length len(t)-1
N = G.number_of_nodes()
dS_hand = np.log(N) * np.diff(sm1) / np.diff(np.log(t))
assert np.allclose(dS, dS_hand), "C(τ) formula mismatch"
print(f"OK: entropy_C[i] = -dS_unnorm/dlogτ at midpoint of (t[i], t[i+1])  (len = {len(dS)})")
```

### 7.4 Check that the cophenetic ultrametric matches the linkage's `1/ρ` projection

```python
import numpy as np
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import squareform
import networkx as nx
from lrgsglib.core import (
    compute_laplacian_properties, compute_normalized_linkage,
    extract_ultrametric_matrix, get_giant_component,
)

G = nx.barabasi_albert_graph(50, 2, seed=0)
giant = get_giant_component(G)
_, _, _, Trho, _ = compute_laplacian_properties(giant)
dists = squareform(Trho)
Z, labels, tmax = compute_normalized_linkage(dists, giant, method="average")
ultra = extract_ultrametric_matrix(Z, giant.number_of_nodes())
# Sanity: cophenetic distances all ≤ root height, ≥ 0
assert ultra.min() >= 0.0
assert ultra.max() <= 1.0 / 1.01 + 1e-9, "root height exceeds 1/1.01 (rescaling broken)"
# Strong triangle inequality on a few random triples
rng = np.random.default_rng(0)
N = ultra.shape[0]
for _ in range(100):
    i, j, k = rng.choice(N, size=3, replace=False)
    assert ultra[i, j] <= max(ultra[i, k], ultra[j, k]) + 1e-12, "ultrametric violated"
print("OK: cophenetic ultrametric is in [0, 1/1.01] and satisfies strong triangle inequality")
```

### 7.5 Check that LRG cache fields are consistent with what `LRGResult` expects

```python
import numpy as np
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE

r = load_lrg_result("Pat_02", "rest_pre", "alpha", "imcoh_abs", IMCOH_LRG_CACHE)
# Shapes
assert r.linkage_matrix.shape == (r.n_nodes - 1, 4), "Z shape wrong"
assert len(r.entropy_tau) == 400, "τ-grid length wrong (default 400)"
# entropy_C is dS = diff(1-S), so length is steps-1 = 399
# (lrgsglib stores it in entropy_C anyway; the alignment quirk for plotting is steps-1)
assert len(r.entropy_C) == len(r.entropy_tau) - 1, \
    "entropy_C length should be entropy_tau length − 1 (np.diff)"
# Ultrametric is condensed; full triangular = N(N-1)/2
N = r.n_nodes
assert r.ultrametric_matrix.shape == (N * (N - 1) // 2,), "ultrametric not in condensed form"
# optimal_threshold is in (0, root_height ≈ 0.9901)
assert 0 < r.optimal_threshold < 1.0
print(f"OK: LRGResult for {r.patient}/{r.band}/{r.phase} has consistent shapes")
```

---

## 8. The dendrogram is our load-bearing object

Per §6 above, our case sits outside the LRG-paper sweet spot. The dendrogram at τ = 1/λ_max is what we have:

- **Subtree leafsets** at native heights (or fcluster cuts) → L4 measures (trace-modules, MRL, CBR variants).
- **Flat partitions at integer-k or h_rel** → L5 measures (Δ_VI, Δ_H, Δ_NMI).
- **Per-leaf cophenetic vector** (row of `ultrametric_matrix`) → L6 measures (CNP per-leaf neighbourhood).
- **Per-(leaf, scale) cluster-mate set** → L7 measures (MSPC).
- **Whole-tree distance** between phase-pairs → L3 measures (Kendall–Colijn λ-sweep, Matching-Cluster).
- **Whole-tree communication distance matrix** (the cached `ultrametric_matrix`) → L1 measures (continuous H2c).

These rungs are codified in the rebuild plan (`/home/giulio/.claude/plans/i-think-there-are-binary-puppy.md` §2) and the decision-rules scope (`.agents/guides/task-persistence-investigation/archive/scalar-vi-era/2026-04-29_decision-rules.md`).

---

## 9. Things to NOT freelance (anti-hallucination checklist)

- `Trho` is **`1/ρ_ij(τ)`** (element-wise inverse of the heat-kernel density at the chosen τ). It is the *communication distance* per Villegas 2025 Eq. 1, NOT a τ-related quantity.
- `LRGResult.ultrametric_matrix` is the *cophenetic ultrametric* derived from UPGMA on `Trho` AFTER the `1.01·max` rescaling. Its values are in `[0, 1/1.01 ≈ 0.9901]` by construction.
- `LRGResult.optimal_threshold` is `0.9 · Δ_n*` where `n* = argmax_n Ψ(n; τ=1/λ_max)`. It's a cophenetic cut height. **Never a τ.**
- `LRGResult.entropy_1_minus_S` is `1 − S(τ)` where `S` is normalised von Neumann entropy of ρ̂(τ). Ranges 0 (small τ) → 1 (large τ).
- `LRGResult.entropy_C` is `−dS/dlog τ` returned by `np.diff` so its **length is `entropy_tau − 1`**. Plot it at `(τ[i] + τ[i+1])/2` if you need alignment.
- The dendrogram is **built once at `τ = 1/λ_max`**, not via a τ-scan. Multiscale comes from the dendrogram cuts, not from τ evolution.
- **Ψ is NOT a τ-scan stability metric.** It's a fixed-τ index over consecutive merge-distance log-gaps in a single dendrogram.
- **C(τ) has a single well-defined peak τ\* (the collapse scale, ≈10/λ_max), not a family of mesoscale peaks** (verified 2026-07-09). Don't read a nested community ladder off it, and don't call C(τ) "smooth/peakless" — the correct statement is "no *multiple* gap-induced peaks". τ\* is a valid coarse-τ *ceiling*; the cross-phase trace is still read at τ_min for a cross-phase-collapse reason (audit_121), **not** because C(τ) has no peak.
- **Don't claim "scale-invariant" or "informational phase transition" findings** in our dataset — those are paper-only arguments for fractal/sparse networks.

---

## 10. Open questions (research-level)

These are NOT settled by the framework guide; they are flagged for explicit investigation:

1. **τ choice for the task-trace question.** DISCHARGED for the cohort trace
   (audit_121, 2026-06-22): τ = 1/λ_max vindicated; the trace is τ-robust
   finest→Fiedler and coarse-τ "gains" are a placebo-confirmed collapse artifact.
   The C(τ) peak τ\* (≈10/λ_max; verified single/unimodal 2026-07-09) coincides with
   that collapse onset and is the natural coarse *ceiling* — used by the epi
   within-phase marker (τ_5 = 10/λ_max ≈ τ\*), NOT by the cross-phase trace. The
   trace/marker split is "scale follows the question": trace at the floor τ_min,
   marker near the ceiling τ\*. Per-node (β→OFC) τ-sensitivity remains a lighter
   open follow-up.
2. **Logarithmic communication distance.** `−log K_ij(τ)` is a metric (not strictly ultrametric except when `K` is itself an ultrametric kernel). It's information-theoretically natural ("nats"). Could change dendrogram shape; may or may not change cohort verdicts. Worth a follow-up audit.
3. **Random-walk Laplacian L_RW = D^{-1} L̂.** Used in the paper §V for community detection; differs from L̂ in its treatment of degree-heterogeneity. May give different dendrograms on the same FC matrix. Not currently in the codebase.
4. **Effective spectral dimension d_S.** For a fully-connected graph the spectral dimension is poorly defined; characterising d_S for our FC matrices would tell us *how close* we are to the Wigner-semicircle limit and therefore how much the paper's "scale-invariant" framing applies (probably very little).
5. **What does it mean to "find structure" in a fully connected weighted graph?** The paper offers **no equivalent question** — they study graphs where structure = topology. We need our own positive definition. Current working answer: "patterns of edge-weight clustering that produce reproducible UPGMA hierarchy, with cross-phase reproducibility tested against within-baseline shuffling".

---

## 11. Cross-references

- Implementation correctness audit (Phase 0): `.agents/guides/task-persistence-investigation/archive/scalar-vi-era/2026-04-29_measure-correctness-audit.md`.
- Decision rules for cohort verdicts: `.agents/guides/task-persistence-investigation/archive/scalar-vi-era/2026-04-29_decision-rules.md`.
- Our outlier case in detail: `memory/lrg_outlier_case_fully_connected.md`.
- τ choice options: `memory/lrg_tau_choice.md`.
- ImCoh derivation: `.agents/guides/02_methods/imcoh-guide.md`.
- Probe-bias control: `.agents/guides/02_methods/probe-bias-guide.md`.
- H2 metric family: `.agents/guides/02_methods/h2-metrics.md`.
- Rebuild plan (the geometric ladder): `/home/giulio/.claude/plans/i-think-there-are-binary-puppy.md` §2 + §3.
