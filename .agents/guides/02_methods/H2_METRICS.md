# Multiscale task-trace metrics: VI, ρ, Δρ

This guide explains the three quantitative metrics used to detect
*task-induced memory in the LRG multiscale structure*: the symmetric
partition distance VI(k), the continuous drift correlation ρ, and the
conditional block-persistence Δρ(k). Each operates on a different level
of the LRG output, each answers a different scientific question, and
each has specific things it can and cannot prove.

This document is deliberately critical — for every metric we state
what it measures, what it does NOT measure, and what the honest
interpretational limits are.

---

## 0. What the LRG produces

For each `(patient, phase, band)` we have an FC adjacency matrix
$A = |\mathrm{ImCoh}|$ and its combinatorial graph Laplacian
$L = D - A$. The LRG pipeline (via `lrgsglib.compute_laplacian_properties`)
produces the following objects at a **principled, information-maximizing
diffusion time** $\tau$:

### (0.a) Choice of τ (Villegas-canonical for fine-grained structure)

$$\tau = \frac{1}{\lambda_{\max}(L)}$$

This is not an arbitrary convenience. In the LRG framework, $\rho(\tau)
\to I/N$ as $\tau \to 0$ (no structural information) and $\rho(\tau)$
loses fine-grained structure as $\tau \to \tau^\star$ (specific-heat peak)
and beyond (coarse-graining takes over). **$\tau = 1/\lambda_{\max}$ is
the finest-resolution diffusion time at which the density matrix carries
meaningful pairwise structural information**, i.e. the τ at which
information retention is maximal before any coarsening has occurred. Each
phase has its own Laplacian and therefore its own $\lambda_{\max}$; the
per-phase choice of τ is the information-maximizing choice adapted to
that phase's graph.

### (0.b) Villegas density matrix and pairwise distance

$$\rho_{ij}(\tau) = \frac{[e^{-\tau L}]_{ij}}{\mathrm{tr}(e^{-\tau L})}, \qquad T_\rho[i,j] = \frac{1}{\rho_{ij}(\tau)}$$

This is the canonical Villegas $1/\rho_{ij}$ information-theoretic
distance (symmetric, zero diagonal by convention). In the implementation
the small asymmetry introduced by floating-point evaluation is removed
via $T_\rho \leftarrow \max(T_\rho, T_\rho^{\!T})$.

### (0.c) Hierarchical clustering on $T_\rho$ (UPGMA)

$$Z = \mathrm{linkage}\bigl(\mathrm{squareform}(T_\rho / t_{\max}),\;\text{method="average"}\bigr)$$

where $t_{\max} = \max(Z[:, 2]) \cdot 1.01$ normalizes the linkage
heights into $\approx [0, 1]$. UPGMA (average linkage) on the Villegas
distance gives the dendrogram.

### (0.d) Cophenetic ultrametric matrix (what H2c actually uses)

$$D_{ij} = \mathrm{cophenet}(Z)_{ij}$$

i.e. the height in $Z$ at which nodes $i$ and $j$ are first merged. This
is by construction an ultrametric matrix (satisfies the strong triangle
inequality). It is stored in every LRG `.npz` cache as the
`ultrametric_matrix` field (condensed form).

### (0.e) Partitions

$$\mathcal{P}_k = \mathrm{fcluster}(Z,\;k,\;\text{criterion="maxclust"})$$

for any $k \in \{2, \ldots, N-1\}$. Each $k$-partition is a cut of $Z$
at the height that produces exactly $k$ clusters.

### Summary of the stack

```
A = |ImCoh|
     │
     ↓  L = D − A,   τ = 1/λ_max(L)
     │
(0.b) T_ρ[i,j] = 1 / ρ_ij(τ)           ← Villegas inverse-density distance
     │
     ↓  UPGMA average-linkage clustering
     │
(0.c) Z = dendrogram (linkage matrix)
     │
     ├──→ (0.d) D[i,j] = cophenetic(Z)   ← ultrametric distance matrix (used by H2c)
     │
     └──→ (0.e) 𝒫_k = fcluster(Z, k)     ← partition at scale k (used by H2d, VI)
```

**The multiscale character comes from the dendrogram hierarchy** — the
cophenetic heights in $Z$ span from fine (near-leaf merges) to coarse
(near-root merges). The τ of the underlying density matrix is fixed at
the information-maximizing value per graph; it is not swept.

Each phase `rpre | tlearn | ttest | rpost` produces its own
$(T_\rho, Z, D, \mathcal{P}_k)$ at its own $\tau = 1/\lambda_{\max}$.
The cross-phase metrics below compare these objects.

---

## 1. VI(k) — symmetric partition distance

### Construction

Given two partitions `𝒫_a` and `𝒫_b` of the same N nodes at the same
scale k:

$$p_a(c) = \frac{|\{i : \mathcal{P}_a(i) = c\}|}{N}, \qquad p_b(c') = \frac{|\{i : \mathcal{P}_b(i) = c'\}|}{N}$$

$$p_{ab}(c, c') = \frac{|\{i : \mathcal{P}_a(i) = c \wedge \mathcal{P}_b(i) = c'\}|}{N}$$

Define Shannon entropy and mutual information of the partitions:

$$H(\mathcal{P}_a) = -\sum_c p_a(c) \log p_a(c)$$

$$I(\mathcal{P}_a; \mathcal{P}_b) = \sum_{c, c'} p_{ab}(c, c') \log \frac{p_{ab}(c, c')}{p_a(c)\, p_b(c')}$$

Then Meilă's (2007) **variation of information**:

$$\boxed{\;\mathrm{VI}(\mathcal{P}_a, \mathcal{P}_b) = H(\mathcal{P}_a) + H(\mathcal{P}_b) - 2\, I(\mathcal{P}_a; \mathcal{P}_b)\;}$$

Equivalently, using conditional entropy $H(X|Y) = H(X) - I(X;Y)$:

$$\mathrm{VI}(\mathcal{P}_a, \mathcal{P}_b) = H(\mathcal{P}_a \mid \mathcal{P}_b) + H(\mathcal{P}_b \mid \mathcal{P}_a)$$

Range: $[0, \log N]$. A proper metric (non-negative, symmetric, satisfies
the triangle inequality).

### What VI(k) measures

- $\mathrm{VI} = 0$ ⟺ $\mathcal{P}_a$ and $\mathcal{P}_b$ are identical partitions.
- $\mathrm{VI}$ small ⟺ most pairs of nodes that co-cluster in one
  also co-cluster in the other (at that k).
- $\mathrm{VI}$ large ⟺ the two partitions disagree about block
  membership.

### What VI(k) does NOT measure

1. **Direction** — $\mathrm{VI}(\mathcal{P}_a, \mathcal{P}_b) =
   \mathrm{VI}(\mathcal{P}_b, \mathcal{P}_a)$. VI cannot answer "did
   A become more like B than B became like A"; it's symmetric.
2. **Magnitude of ultrametric change** — Two very different ultrametric
   matrices can yield the same partition at a particular k if the
   changes happen *within* clusters rather than across cluster
   boundaries. VI at one k is a lossy projection of (1).
3. **Specific pair identity** — VI aggregates over the whole partition.
   It does not isolate pairs that specifically moved.

### Where VI(k) is the right tool (H1, H3, H4)

- **H1 (task stability)**: "are task_learn and task_test similar to
  each other at every k?" — symmetric, partition-level. Exactly what
  VI answers.
- **H3 (phase-type clustering)**: "is mean VI(within-type) < mean
  VI(cross-type)?" — symmetric pair-type comparison, averaged.
- **H4 (frequency gradient)**: "do some bands show systematically
  smaller cross-phase VI than others?" — again symmetric.

### Where VI(k) was the wrong tool (H2a, currently failing)

**H2a** as originally formulated:

$$\mathrm{VI}(\mathcal{P}^{rpre}_k, \mathcal{P}^{rpost}_k) \stackrel{?}{>} \mathrm{VI}(\mathcal{P}^{ttest}_k, \mathcal{P}^{rpost}_k)$$

— "is rpost's partition closer to ttest's than to rpre's?" This fails
FDR at n=9 for every band. Why:

1. **Symmetric distance loses direction.** The trace claim is directed
   ("rpost inherited structure from ttest") but VI is symmetric.
2. **Difference of differences.** Subtracting two symmetric distances
   stacks variance; patient heterogeneity overwhelms the effect.
3. **Partition-level coarsening at small k.** At k=2 or k=3 most
   pairs are co-clustered everywhere, and VI is near zero in all
   directions; signal lives at larger k but is variable there.

A better VI-family test would be the **directed** conditional-entropy
version:

$$H(\mathcal{P}^{rpost}_k \mid \mathcal{P}^{ttest}_k) \stackrel{?}{<} H(\mathcal{P}^{rpost}_k \mid \mathcal{P}^{rpre}_k)$$

"Given task, rpost is more predictable than given rest_pre." Directed,
partition-level, cleaner. This is **H2a′** — see §1a below.

---

## 1a. H2a′ — directed conditional entropy (partition-level memory)

### Construction

For each (patient, band, k):

$$\delta H_k \;=\; H\bigl(\mathcal{P}^{rpost}_k \mid \mathcal{P}^{ttest}_k\bigr) \;-\; H\bigl(\mathcal{P}^{rpost}_k \mid \mathcal{P}^{rpre}_k\bigr)$$

where $H(Y \mid X) = H(Y) - I(X; Y)$ is the standard conditional-entropy
decomposition (nats). **Trace hypothesis**: $\delta H_k < 0$.

Implemented via `conditional_entropy(labels_y, labels_x)` in
`src/lrg_eegfc/utils/metrics/vi.py`, which reuses the same joint-histogram
machinery as `compute_vi` — so $\mathrm{VI}(A,B) = H(A|B) + H(B|A)$ holds
to machine precision (verified).

### What δH_k measures (and what it does not)

- $\delta H_k < 0$ ⟺ the task partition is a strictly better predictor of
  rpost's clustering than rest_pre is. This IS a directed,
  partition-level statement of memory.
- **Lossy in the same way VI is at a given k**: conditional entropy
  aggregates all rpost↔x co-occurrence; it does not isolate specific
  pair identities, nor the magnitude of merge-height change within
  stable partitions.
- **Vulnerable to partition-granularity coarsening at small k**: at
  k = 2 most items are in one big cluster, so $H(\cdot \mid \cdot) \to 0$
  and $\delta H_k \to 0$ regardless of structural memory.
- **Same time-drift concern as VI** — a trace-shaped session drift would
  also reduce $H(rpost \mid ttest)$ more than $H(rpost \mid rpre)$. H2e
  (§2a) is the cross-check.

### Why H2a′ is the right complement to ρ and Δρ

ρ (H2c) uses continuous ultrametric distances, Δρ (H2d) uses the
co-cluster indicator. Both live downstream of the partition. δH_k acts
directly on partition labels as categorical variables, asymmetrically —
so it closes the "directed, partition-level" cell of the metric matrix
that VI left blank. A pass here, combined with H2c and H2d, gives
**three-level convergence** (continuous ultrametric + conditional
block + categorical partition) as the memory claim.

### Pre-registered test (matches the script)

Per band, k-averaged $\delta H$ per patient; one-sample Wilcoxon on
$-\delta H$ (alternative greater); FDR-BH m=6; rank-biserial effect
size; 10k bootstrap CI; LOO worst-p; Pat_03-drop sensitivity.
Scale-localization via sign-flip cluster-permutation over k (reusing
the H2d implementation). Pre-reg pass: ≥3 of 6 bands q<0.05.

**Current status (n=9, 2026-04-24+)**: see
`data/reports/imcoh_vi/h2a_prime_conditional_entropy.md`. 0 of 6 bands
pass at the stringent k-averaged FDR level; δ, β, γ_h all trend in the
expected direction with p ≈ 0.07–0.09, and α shows a narrow
cluster-permutation-significant window at k=2-4 (p=0.050). θ trends in
the opposite direction, consistent with the θ-ergodic interpretation
from H2d. **Conclusion**: H2a′ moves to supplementary as an honest
limitation — partition-level directed memory is at the edge of
detectability at n=9; ρ (continuous) and Δρ (conditional) remain the
sensitive detectors.

---

## 2. ρ — continuous ultrametric drift correlation (H2c)

### Construction

For each patient × band with all four phases:

1. Retrieve the cophenetic ultrametric matrices $D^{rpre}, D^{rpost}, D^{ttest}$
   (each one computed from its phase's UPGMA dendrogram on
   $T_\rho^{phase} = 1/\rho^{phase}(\tau = 1/\lambda_{\max}(L^{phase}))$;
   see §0.d).
2. For each upper-triangle pair $(i,j)$ with $i<j$, compute two change
   vectors:

   $$\Delta_{\text{rest}}(i,j) = D^{rpost}_{ij} - D^{rpre}_{ij}, \qquad \Delta_{\text{task}}(i,j) = D^{ttest}_{ij} - D^{rpre}_{ij}$$

   Both are $\mathbb{R}$-valued, one number per pair. Length of each
   vector: $N(N-1)/2$.

3. Compute the **Spearman rank correlation**:

   $$\boxed{\;\rho = \frac{\mathrm{cov}(\mathrm{rank}\,\Delta_{\text{rest}},\;\mathrm{rank}\,\Delta_{\text{task}})}{\sigma_{\mathrm{rank}\,\Delta_{\text{rest}}}\,\sigma_{\mathrm{rank}\,\Delta_{\text{task}}}}\;}$$

   One scalar ρ per (patient, band). Range: $[-1, +1]$.

### What ρ measures

- $\rho > 0$: pairs that became *closer* in task (negative $\Delta_{\text{task}}$)
  also tended to become closer in rest_post (negative $\Delta_{\text{rest}}$);
  and pairs that became *farther* in task became farther in rest_post.
- $\rho \approx 0$: the rest-post pair-wise shift is unrelated to the
  task shift.
- $\rho < 0$: rest_post drifted in the opposite direction from task
  (unusual; would falsify the trace claim in reverse).

The test H2c is a **one-sample Wilcoxon signed-rank** on the per-patient
ρ values per band, with H₀: ρ = 0 vs H₁: ρ > 0. FDR-BH across the 6
bands, with 10k bootstrap CI on the patient-mean.

### Why ρ is statistically powerful

- **Each ρ aggregates across** $N(N-1)/2$ pairs. Pat_02 at 116 channels
  = 6670 pairs feeding one ρ value. Per-patient noise is suppressed.
- **Spearman is rank-based**: robust to nonlinear monotonic
  transformations of ultrametric distance, unaffected by outlier
  magnitudes.
- **Integrated across scales**: ρ uses the full ultrametric matrix
  simultaneously, without choosing a single k. Escape the k-choice
  problem entirely.

### What ρ does NOT measure

1. **Magnitude of the task effect.** ρ is *scale-free*: pairs that
   barely moved but moved "in the right direction" count as much in
   rank terms as pairs that moved dramatically. ρ > 0 does not imply
   "large memory"; it implies "directionally consistent memory."
2. **Specific pair identity.** ρ averages over all pairs. It does not
   isolate which pairs carry the memory.
3. **Causal claim.** ρ is a correlation on change vectors. It says
   "these two perturbations are similar in direction" — it does not
   prove the rest_post perturbation was *caused* by task. See §5.

### Why ρ is the right tool (H2c)

"Did the overall structural reshaping at rest_post resemble the
overall reshaping during task?" — continuous, global, directional.
This is the most general test of memory-as-directional-inheritance.

### The time-drift vulnerability (critical)

sEEG phases are recorded sequentially:
`rest_pre → task_learn → task_test → rest_post`. Any property that
drifts monotonically across the session (electrode adherence,
cortical excitability, fatigue, temperature, regression to the mean)
would affect later phases more than earlier, and would produce:

- $\Delta_{\text{task}}$ partially shifted from rpre in some direction
- $\Delta_{\text{rest}}$ further shifted in the same direction

Under pure time drift with zero task-specific memory, **ρ > 0 would
be observed anyway**. This is a serious confounder for the
"task-induced memory" interpretation.

Partial defenses already in the paper:
- **H1 stability** shows task_learn and task_test look task-like to
  each other, not just "late-in-session" similar.
- **Replication with target = task_learn** (earlier than task_test)
  also gives ρ > 0 — weakens a pure monotonic-drift story.

Needed but not yet run:
- **Split-half noise floor (H2e)**: quantify how much of the ρ magnitude
  is above what is observed between two halves of the same resting
  recording.
- **Condition on task-specific signatures** that can't drift with time
  (e.g. spatial locus of the change).

---

## 2a. H2e — split-half noise floor for ρ (and Δρ)

### Construction

For each patient × band × resting phase $\phi \in \{rpre, rpost\}$:

1. Split the raw time series $X_\phi$ in half along time into $X_\phi^A$
   and $X_\phi^B$.
2. Compute $|\mathrm{ImCoh}|$ on each half with $\mathrm{nperseg} =
   \mathrm{nperseg\_for\_fs}(fs) / 2$ (1-s windows instead of 2-s) so
   that the spectral sample count per half stays comparable to the
   full-duration run.
3. Run LRG on each half → cophenetic matrices $D_\phi^A, D_\phi^B$ and
   linkage matrices $Z_\phi^A, Z_\phi^B$.

### Three quantities per patient × band

- **Drift-only ρ** (null that would be observed under pure session
  drift, zero task memory):

  $$\rho_{\text{null,drift}} = \mathrm{Spearman}\bigl(D_{rpre}^B - D_{rpre}^A,\;\; D_{rpost}^B - D_{rpost}^A\bigr)$$

  Each $\Delta^{split}_\phi = D_\phi^B - D_\phi^A$ captures whatever
  drift-consistent change happened during $\phi$. If only drift is
  present, those two vectors align regardless of task.

- **Reliability ceiling** (noise ceiling — how well the ultrametric
  matrix reproduces itself across halves of the same recording):

  $$\rho_{\text{within},\phi} = \mathrm{Spearman}\bigl(\mathrm{uppertri}(D_\phi^A),\;\; \mathrm{uppertri}(D_\phi^B)\bigr)$$

  The true cross-phase $\rho$ should be below this and ideally above
  $\rho_{\text{null,drift}}$.

- **Δρ drift-null** (triple-analog): substitute (rpre_A, rpre_B,
  rpost_A) into the H2d triple. Pairs "newly co-clustered in rpre_B
  vs rpre_A" play the role of "task-induced"; their persistence into
  rpost_A is measured as $\rho_{\text{task}}^{null}$. Then
  $\Delta\rho^{null}(k) = \rho_{\text{task}}^{null}(k) -
  \rho_{\text{inert}}^{null}(k)$.

### Pre-registered test

Paired Wilcoxon per band across patients, one-sided "cross > null":

- ρ primary: $\rho^{cross}_{pat,band} > \rho^{null\,drift}_{pat,band}$,
  FDR-BH m=6. Pass = ≥4 of 6 bands q<0.05.
- Δρ secondary: k-averaged $\Delta\rho^{cross}_{pat,band} >
  \Delta\rho^{null}_{pat,band}$, FDR-BH m=6.

### What H2e proves (and what it does not)

- **Pass**: ρ > 0 is genuinely above what any within-session drift
  process yields; the task is contributing a signal beyond generic
  session non-stationarity. H2c survives this cross-check.
- **Partial**: drift is part of the ρ signal but not all of it;
  H2d's $\rho_{\text{inert}}$-corrected metric becomes the load-bearing
  directed claim.
- **Fail**: ρ magnitude is not distinguishable from a drift null;
  H2c drops as a standalone result; H2d + H2a′ carry whatever
  remains.

**Honest caveat**: half-duration D matrices are noisier than
full-duration (half the integration time per Welch CSD), so the null
distribution is modestly inflated. H2e is therefore a conservative
floor, biased AGAINST the trace claim — passing it is a strong result;
failing it needs careful interpretation (could be true null, could be
the half-duration noise penalty eating the signal). This is documented
in the writeup explicitly.

**Current status (n=9, 2026-04-24+)**: see
`data/reports/imcoh_vi/h2e_split_half.md`.

---

## 3. Δρ(k) — conditional block persistence (H2d)

### Construction

Define the co-cluster indicator at scale k:

$$c^{phase}_k(i,j) = \mathbb{1}\bigl[\mathcal{P}^{phase}_k(i) = \mathcal{P}^{phase}_k(j)\bigr]$$

This is a binary $\{0,1\}$ scalar per (pair, phase, k).

Partition all $N(N-1)/2$ pairs into four conditional sets based on
$(c^{rpre}_k, c^{ttest}_k)$:

|  | $c^{rpre}_k = 0$ | $c^{rpre}_k = 1$ |
|---|:---:|:---:|
| $c^{ttest}_k = 1$ | **task-induced** | always-together |
| $c^{ttest}_k = 0$ | **never-together** | task-dissolved |

The **task-induced** set = pairs that were *not* co-clustered at
rest_pre but *were* co-clustered during task_test. These are the
candidate "task-memory-carrying" pairs.

The **never-together** set = pairs that were *not* co-clustered at
either rest_pre or task_test. These are the baseline "inert" pairs —
no task event brought them together.

Now measure the rate at which each set has $c^{rpost}_k = 1$:

$$\rho_{\text{task}}(k) = \frac{\bigl|\{(i,j):\, c^{ttest}_k=1,\, c^{rpre}_k=0,\, c^{rpost}_k=1\}\bigr|}{\bigl|\{(i,j):\, c^{ttest}_k=1,\, c^{rpre}_k=0\}\bigr|}$$

$$\rho_{\text{inert}}(k) = \frac{\bigl|\{(i,j):\, c^{rpre}_k=0,\, c^{ttest}_k=0,\, c^{rpost}_k=1\}\bigr|}{\bigl|\{(i,j):\, c^{rpre}_k=0,\, c^{ttest}_k=0\}\bigr|}$$

And the **task-trace excess**:

$$\boxed{\;\Delta\rho(k) = \rho_{\text{task}}(k) - \rho_{\text{inert}}(k)\;}$$

One scalar per (patient, band, k). We sweep k ∈ [2, 49] and both
report the k-averaged Δρ per patient and perform cluster-based
permutation over k to localize the scale range.

### What Δρ(k) measures

- $\rho_{\text{task}}(k)$ answers "of all the pairs that task specifically
  brought together (that weren't already together pre-task), what
  fraction stayed together at rest_post?"
- $\rho_{\text{inert}}(k)$ answers "of all the pairs that nobody ever
  clustered, what fraction ended up clustered at rest_post anyway?"
  (This is the **baseline rate** for spontaneous clustering at rest_post.)
- $\Delta\rho(k) > 0$ ⟺ task-induced pairs persist *at a higher rate
  than chance baseline*. This IS the block-level memory claim.

### Why Δρ is stronger than plain "fraction persisting"

Without the $\rho_{\text{inert}}$ baseline, a reviewer can say:
"maybe rest_post just spontaneously has more clusters than rest_pre;
every random pair is equally likely to be together." The $\rho_{\text{inert}}$
measures exactly that spontaneous rate. $\Delta\rho > 0$ rejects
"spontaneous clustering" as a sufficient explanation.

This is conceptually a **causal** metric: we conditioned on pairs that
were specifically changed by task. The question is whether that
change persists *above the baseline that applies to all pairs*.

### What Δρ(k) does NOT measure

1. **Global structural similarity** — Δρ counts pair-wise agreements,
   not holistic dendrogram shape. Two dendrograms with very different
   merge-height structures can give similar Δρ if they happen to
   agree on the specific pairs conditioned on.
2. **Magnitude of task-induced change** — Δρ doesn't know how much
   the ultrametric distance changed, only whether a pair crossed the
   co-cluster threshold.
3. **Immunity to time drift** (though stronger than ρ here). If rest_post
   simply has higher baseline clustering rates at all pairs (uniform
   drift), both ρ_task and ρ_inert inflate roughly equally — so
   Δρ is somewhat protected. However, if the drift is *pair-specific*
   and happens to correlate with task-induced pairs for non-task
   reasons, Δρ can still be biased.

### Why Δρ is the right tool (H2d)

"Do task-induced pair co-activations specifically persist at rest_post
at a higher rate than pairs that were never co-active?" — directed,
causal, scale-resolved. The closest thing we have to a clean memory
claim at the block level.

### Scale-range caveat

The k range is a judgment call:
- k = 2 or 3: most pairs are in the same cluster (almost all
  $c^{phase}_k = 1$) — Δρ is near zero and noisy.
- k near N−1: most pairs are in different clusters — induced set
  size shrinks, variance blows up.
- Practical range: k ∈ [2, 49] sweeps the "meaningful community"
  regime. Results shown to be robust in that window via cluster-based
  permutation.

---

## 4. How the three metrics relate

All three are functions of the same underlying LRG output but operate
at different levels of abstraction and ask different directional /
conditional questions.

| Property | VI(k) | ρ | Δρ(k) |
|:---------|:-----:|:-:|:-----:|
| Operates on | partitions | signed Δ vectors on ultrametric distances | partitions (triple) |
| LRG stack level | (3) | (1) | (2) → binary from (3) |
| Arg-swap symmetric? | yes (VI(A,B)=VI(B,A)) | yes (Spearman(X,Y)=Spearman(Y,X)) | **no** (rpre/tt/rpost play distinct roles) |
| Sign-aware on inputs? | no (partitions are categorical) | **yes** (each Δ is signed after−before) | yes (co-cluster change) |
| Scale-resolved? | yes | no (integrated across all cophenetic heights) | yes |
| Captures specific pair ID? | no | no | yes |
| Captures magnitude? | partially | no (rank-only) | no (binary) |

**Note on "direction":** ρ is often described informally as "directional"
because ρ > 0 is a meaningful claim about the *signs* of the Δ vectors
(pairs moved the same way in both rest-drift and task-drift). But Spearman
itself is symmetric in its two arguments — ρ(Δ_rest, Δ_task) = ρ(Δ_task,
Δ_rest). The directional content lives in the SIGN of each Δ
(after−before), not in the order of arguments to Spearman. Δρ, by
contrast, is genuinely asymmetric at the role level: swapping the roles
of ttest and rpost changes the value.

### Logical relationships

- If two ultrametric matrices $D^a$ and $D^b$ are close in $L_\infty$
  (all merge heights similar), their dendrograms are close, their
  partitions at every k are close, and hence $\mathrm{VI}(\mathcal{P}^a_k,
  \mathcal{P}^b_k) \approx 0$ and $\rho \to +1$ and Δρ is large.
- **The converse fails**: VI ≈ 0 at one k does NOT imply $D^a = D^b$.
  Different dendrograms can produce the same partition at a specific k.
- **ρ and Δρ do NOT need to agree perfectly**. ρ captures directional
  drift of the full geometry; Δρ captures specific pair agreement
  conditional on task. A band could have ρ > 0 (directional drift
  preserved) but Δρ ≈ 0 (no specific block memory above baseline).
  This is exactly what θ appears to show in our n=9 data: "the
  geometry drifted roughly the right way, but no specific task-induced
  block preferentially persisted."

### Why we need all three

Each metric isolates a different aspect of "structural memory":

- **VI** → "does the task phase structurally look different from
  rest, and does rest_post partially look like task at the partition
  level?" (H1/H3 yes; H2a no.)
- **ρ** → "did the overall geometry drift in the same direction?"
  (H2c: yes, universally.)
- **Δρ** → "did specific task-induced pair clusters persist above
  baseline?" (H2d: yes in 5 bands, θ is the weakest.)

Dropping any one loses information. Dropping VI loses the cluster-range
localization of H1 and H3; dropping ρ loses the most powerful single-
number universality result; dropping Δρ loses the causal conditional
claim.

---

## 5. How they prove (or fail to prove) the thesis

**Thesis**: Cognitive task execution leaves a multiscale trace in the
post-task resting-state functional-connectivity hierarchy.

This decomposes into three sub-claims, each with a proper test:

### Sub-claim A — "Task induces a real structural change"

Required to rule out the null "task and rest are indistinguishable."

Test: VI(task_learn, task_test) < VI(anything-cross-type) at the same k.
That's H1. Passes in every band (all q < 0.02 FDR at n=9) with
cluster-corrected support across most of k.

**Also required implicitly**: mean VI(rpre, ttest) is non-zero — task
actually moves the structure. Covered by the H3 within < cross result
(which requires VI(within-type) < VI(cross-type), hence non-zero
cross-VI).

### Sub-claim B — "Rest_post is NOT a random resting state: it
carries directional information from task"

This is the directional-memory claim.

Test: ρ(Δ_rest, Δ_task) > 0. That's H2c. Passes in every band at
n=9 with 53/54 patient×band cells positive.

**Critical vulnerability**: time drift (§2, time-drift subsection).

Partial defenses in place (H1 task stability, task_learn replication).
Not-yet-run defense: H2e split-half noise floor.

### Sub-claim C — "Rest_post preserves specific task-induced block
structures, above baseline spontaneous clustering"

This is the block-level-memory claim.

Test: Δρ(k) > 0 for k in the meaningful community range. That's H2d.
Passes in every band at n=9 (all q < 0.005 FDR, k ∈ [2, 49]), with
band-heterogeneity: α strongest (Δρ = 0.28), θ weakest (Δρ = 0.16,
significantly < α and δ by Bonferroni post-hoc).

The band heterogeneity partially defends against time-drift: a
uniform time drift would inflate all bands equally, but θ's Δρ is
visibly smaller than α's.

### Sub-claim D — "The trace is quantitatively band-specific"

This is the band-heterogeneity claim.

Tests:
- **Descriptive**: per-band FDR tables show all 6 bands pass H2c and
  H2d, but Δρ magnitude varies (α > γ_h ≈ δ ≈ γ_l ≈ β > θ).
- **Omnibus**: Friedman χ²(5) on per-patient Δρ vectors. *Does not
  reject* exchangeability at n=9 (p = 0.20).
- **Post-hoc**: paired Wilcoxon θ-vs-α on Δρ with Bonferroni m=5:
  *p = 0.010 ★* (significant). θ-vs-δ: p = 0.049 ★.

**So the band-heterogeneity claim is not carried by an omnibus
rejection** — that would require more patients. It is carried by the
Bonferroni-corrected post-hoc dissociation θ vs α, which is the
strongest statistically defensible version of "different bands carry
different trace magnitudes" at current cohort size.

---

## 6. Honest limits of the framework

What the current metric suite **cannot** establish with current data:

1. **Clean causal attribution of ρ > 0 to task** (as opposed to time
   drift) — requires H2e split-half noise floor.
2. **Partition-level directional memory** at the level of conditional
   entropy — requires implementing $H(\mathcal{P}^{rpost}_k \mid
   \mathcal{P}^{ttest}_k) < H(\mathcal{P}^{rpost}_k \mid
   \mathcal{P}^{rpre}_k)$. This would corroborate or contradict ρ/Δρ
   at the partition level and is a ~30-minute implementation.
3. **Band-exchangeability rejection at omnibus level** — Friedman
   doesn't reject at n=9. Requires a larger cohort for confident
   "bands behave systematically differently" claim.
4. **Spatial specificity of the trace** — we haven't mapped which
   brain regions / electrodes drive ρ and Δρ. If the trace is
   specifically in task-engaged regions, that's strong support; if
   it's diffuse, the interpretation weakens.
5. **Decay dynamics** — does the trace attenuate over the rest_post
   recording? Currently rest_post is treated as a single block.

What the current metric suite **does** establish, conservatively:

- Task state is structurally distinct from rest state (H1, H3).
- Rest_post is structurally closer to task state than to rest_pre
  *in directional drift* (H2c) and *in specific-pair block persistence*
  (H2d), at rates exceeding spontaneous-clustering baselines, in every
  frequency band.
- Quantitatively, the block-level trace is heterogeneous across bands,
  with α carrying the strongest trace and θ the weakest (significant
  α-vs-θ post-hoc dissociation; other band dissociations descriptive
  only).

---

## 7. Quick reference: what number supports what claim

| Claim | Primary evidence | File |
|:------|:-----------------|:-----|
| Task is structurally distinct | H1 VI(TL,TT) < cross, all bands q<0.02 | `data/reports/imcoh_vi/rigorous_tests.md` |
| Phase-type is a real axis | H3 within<cross, all bands q<0.03 | `data/reports/imcoh_vi/rigorous_tests.md` |
| Directional memory in rest_post | H2c ρ>0, all bands q<0.005, 53/54 cells + | `data/reports/imcoh_vi/h2c_ultrametric_drift.md` |
| Block-level memory above baseline | H2d Δρ>0, all bands q<0.005 | `data/reports/imcoh_vi/h2d_coactivation_persistence.md` |
| Band heterogeneity (α > θ block memory) | post-hoc θ vs α, Bonf p=0.010 | `data/reports/imcoh_vi/h2_band_selectivity.md` |
| Frequency gradient | H4 Kendall W=0.23, p=0.04, n=10 | `data/reports/imcoh_vi/rigorous_tests.md` |
| Ergodic θ interpretation | Δρ_θ = 0.16 vs Δρ_α = 0.28, ρ_inert_θ highest | `data/reports/imcoh_vi/h2d_coactivation_persistence.md` |

---

## 8. Provenance — what's from the original LRG literature, what's ours

It matters for both honest attribution and for anticipating reviewer
questions to be clear about which pieces of this framework come from
the established LRG literature and which are our contribution.

### From the original LRG framework (Villegas–Gabrielli et al.)

Primary references: Villegas, Gabrielli, Gili, Caldarelli (2023,
*Nature Physics*), *"Laplacian renormalization group for
heterogeneous networks"*; and related work by Poggialini, Villegas,
Muñoz, Gabrielli.

- The **heat-kernel diffusion** on the graph Laplacian:
  $\rho(\tau) = e^{-\tau L} / \mathrm{tr}(e^{-\tau L})$.
- The **von Neumann entropy** $S(\tau) = -\mathrm{tr}(\rho \log \rho)$
  and **specific heat** $C(\tau) = -dS/d\log\tau$.
- The **characteristic scale** $\tau^\star$ identified at the $C(\tau)$
  peak as the diffusion time at which coarse-graining proceeds fastest;
  beyond $\tau^\star$, fine-grained information has been integrated out.
- The **information-maximizing choice** $\tau = 1/\lambda_{\max}(L)$ —
  the finest-resolution diffusion time at which the density matrix
  carries meaningful pairwise structural information (below it
  $\rho \approx I/N$, at or past $\tau^\star$ coarsening has erased
  fine detail). This is the τ we evaluate $\rho$ at.
- The **information-theoretic pairwise distance** $T_\rho[i,j] =
  1/\rho_{ij}(\tau)$ (Villegas inverse-density distance).

The `lrgsglib` library used by this project implements the above
directly.

### What we added (original to this project)

- **LRG as a cross-state comparison backbone.** Applying the LRG's
  hierarchical output to compare *the same physical FC network
  across distinct brain states* (rest_pre → task → rest_post) rather
  than analyzing one network in isolation.
- **ρ = Spearman($\Delta_\text{rest}$, $\Delta_\text{task}$) on the
  ultrametric shift vectors** (H2c). Related conceptually to
  Representational Similarity Analysis / Mantel tests, but applied
  to LRG ultrametric matrices across brain states is new in this
  project.
- **Δρ — conditional block persistence** (H2d). The four-way pair
  conditioning on $(c^{rpre}_k, c^{ttest}_k)$ and the $\rho_\text{inert}$
  baseline are our formulation. Related to "partition stability" and
  "co-activation persistence" concepts in community-detection
  literature but not a drop-in published metric.
- **The ρ vs Δρ dissociation** as a continuous-direction vs
  discrete-block-memory framing.
- **Band-selectivity interpretation** — α as strongest-trace band,
  θ as ergodic-baseline band — inferred from the ρ / Δρ / Friedman
  post-hoc pattern.

### τ and k are two different parameters — do not conflate them

The pipeline has **two distinct parameters** that live in different
spaces:

1. **τ — network resolution (diffusion time).** Enters the density matrix
   $\rho(\tau) = e^{-\tau L}/\mathrm{tr}(e^{-\tau L})$ and therefore the
   pairwise distance $T_\rho = 1/\rho_{ij}(\tau)$. Units: diffusion time.
   Set to $\tau = 1/\lambda_{\max}(L)$ (Villegas-canonical,
   information-maximizing choice — finest resolution at which the density
   matrix carries meaningful pairwise structural information; below it
   $\rho \approx I/N$, past $\tau^\star$ coarsening has already erased
   fine information). One τ per graph (per phase per patient per band).
   *Fixed by the Villegas rationale — we do not sweep it.*
2. **k — hierarchy cut (partition scale).** Indexes partitions of the
   UPGMA dendrogram $Z$ produced by clustering $T_\rho$:
   $\mathcal{P}_k = \mathrm{fcluster}(Z, k, \text{"maxclust"})$. Each $k$
   corresponds to a specific merge height $D_k$ in $Z$ (the height that
   yields exactly $k$ clusters). Units: number of clusters (or equivalently
   cophenetic distance along $Z$). *Swept in H2d, integrated by using the
   full cophenetic distance matrix in H2c.*

An "optimal" $k^\star$ for a natural flat clustering is identified via the
**partition stability index** Ψ (Poggialini / Villegas — gaps on the
dendrogram, not C(τ) peaks on the diffusion). Let $D_i$ be the linkage
merge distances in descending order ($D_0$ = root, $D_{N-1}$ = first
merge). Define

$$\sigma_i = \frac{\log_{10} D_i - \log_{10} D_{i+1}}{\log_{10} D_0 - \log_{10} D_{N-1}}, \qquad k^\star = 1 + \arg\max_i \sigma_i$$

Ψ peaks where the log-gap between consecutive merges is largest — i.e.
where the dendrogram has its most stable (widest) plateau. The cached
`optimal_threshold` field in every LRG `.npz` is $D_{k^\star} \cdot 0.9$;
it is a **property of the hierarchy**, not of the underlying diffusion,
and has **nothing to do with C(τ) peaks**. Cached for diagnostic /
visualization use only.

In our H2 pipeline:
- $\tau = 1/\lambda_{\max}(L)$ is fixed per graph; not swept.
- $k$ IS swept: H2d integrates over $k \in [2, 49]$; H2c uses the full
  cophenetic matrix of all merge heights simultaneously (all $k$ at
  once). We do **not** use `optimal_threshold` (= the Ψ-picked $k^\star$)
  for any statistical test.

The multiscale character of H2c and H2d therefore comes entirely from
the hierarchy on $T_\rho$ — from cophenetic distances and k-cuts on the
UPGMA tree — not from any τ sweep. We evaluate the density matrix at
the Villegas-canonical finest-resolution τ and let the dendrogram
carry the multiscale axis.

---

## 9. Cross-references

- FC metric rationale: `.agents/guides/02_methods/IMCOH_GUIDE.md`
- Probe-bias caveat: `.agents/guides/02_methods/PROBE_BIAS_GUIDE.md`
- Writing-agent handoff (cohort-level results): `.agents/reports/MULTISCALE_TASK_TRACE_FOR_WRITING.md`
- Era index: `.agents/reports/PIPELINE_STATUS.md`
- Implementation scripts:
  - VI: `scripts/01_compute/compute_imcoh_vi.py`
  - ρ (H2c): `scripts/01_compute/h2c_ultrametric_drift.py`
  - Δρ (H2d): `scripts/01_compute/h2d_coactivation_persistence.py`
  - Band typology: `scripts/01_compute/h2_band_typology.py`
  - Friedman: `scripts/01_compute/h2_band_selectivity.py`
