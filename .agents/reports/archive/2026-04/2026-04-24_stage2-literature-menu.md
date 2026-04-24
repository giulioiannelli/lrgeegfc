---
name: stage2-literature-menu
type: report
era: COHORT_N9
status: dead
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Stage 2 Literature Menu — Scalar Tree-Similarity Measures

**Date:** 2026-04-24
**Scope:** Scalar dissimilarity for rooted binary trees with **aligned
leaves** (same channel identity across every tree). N ≈ 84–116 leaves per
tree; one UPGMA dendrogram per (patient, phase, band) from an LRG Laplacian.
Goal: one number per pair that jointly uses **merge heights AND branching
topology**.

**Stage 1 tried (skip):** cophenetic Pearson/Spearman, Baker's gamma, L1/L2
on ultrametric vectors, weighted/mean ARI, weighted/mean VI (k ∈ 2..20),
top-k merge overlap (k = 5, 10, 20) — all cophenetic or flat-cut, none blend
topology with heights.

**Data-specific constraint:** every tree has the same labelled leaf set, so
methods designed for *labelled* trees on a common leaf set are the natural
fit.

---

## 1. Kendall–Colijn metric (KC)

**Paper.** Kendall M. & Colijn C. (2016), "Mapping Phylogenetic Trees to
Reveal Distinct Patterns of Evolution," *Mol. Biol. Evol.* 33(10): 2735–2743.
DOI: [10.1093/molbev/msw124](https://doi.org/10.1093/molbev/msw124).
Preprint: [arXiv:1507.05211](https://arxiv.org/abs/1507.05211).

**Construction.** For each of the C(N, 2) ordered leaf pairs build two
vectors:
- **m(T)** — number of *edges* from root to MRCA of the pair (leaf
  self-entries = 1). Pure topology.
- **M(T)** — *summed branch lengths* from root to that MRCA. Pure heights.

The KC vector is v_λ(T) = (1 − λ) m(T) + λ M(T), λ ∈ [0, 1]. Distance:
d_KC(T₁, T₂; λ) = ‖v_λ(T₁) − v_λ(T₂)‖₂. λ = 0 topology only, λ = 1 heights
only, λ ∈ (0, 1) genuine blend. Report either the curve d_KC(λ) or a
cohort-fitted λ*.

**What it captures vs the 14.** The m-vector is MRCA *depth in edges*
per pair, genuinely new — not a cophenetic quantity (that would be M
alone, close to Stage-1 L1/L2). The joint (m, M) formulation penalizes
trees that agree on heights but disagree on which leaves coalesce first
(and vice versa). None of the Stage-1 scalars has this split.

**Python implementation.** No maintained Python package. Two routes:
(a) port `treespace::treeVec` / `TreeDist::KendallColijn` from R — the
algorithm is ~30 LOC (one post-order traversal to tag MRCA depth per
pair); (b) `rpy2` wrapper around `treespace::treeVec`.

**Complexity.** Vector construction O(N²) per tree with LCA pre-processing;
distance O(N²). Milliseconds at N ≈ 100.

**Caveats.** Unbounded — rescale before cross-cohort comparison.
Root-sensitive (fine here: UPGMA has an unambiguous root). λ must be
swept; the 2016 supplement recommends reporting the full curve.

---

## 2. Weighted Robinson–Foulds / Branch-Score Distance (wRF, BSD)

**Papers.** Robinson D.F. & Foulds L.R. (1981), "Comparison of phylogenetic
trees," *Math. Biosci.* 53: 131–147. Weighted version: Kuhner M.K. &
Felsenstein J. (1994), "A simulation comparison of phylogeny algorithms …,"
*Mol. Biol. Evol.* 11(3): 459–468 (branch-score distance).

**Construction.** Every internal edge induces a cluster (the leaves below
it). For trees T₁, T₂ take the multiset of (cluster, branch-length) pairs.
d_wRF = Σ_c |w₁(c) − w₂(c)|, setting w(c) = 0 if c is absent. BSD is the L2
variant. Clusters unique to one tree contribute their own length; shared
clusters contribute length differences — topology + heights in one number.

**What it captures vs the 14.** Weights each internal edge by its length,
simultaneously scoring topology (absent clusters) and heights (matched
clusters with different lengths). Stage-1 ARI/VI are flat-cut partition
scores with no length weighting.

**Python implementation.** Yes, first-party:
`dendropy.calculate.treecompare.weighted_robinson_foulds_distance` and the
L2 `euclidean_distance` (BSD). `ete3.Tree.robinson_foulds` is topology-only
— do **not** use for wRF.

**Complexity.** O(N log N) via split-hash tables.

**Caveats.** Unbounded; normalize by Σw₁ + Σw₂. Notoriously noisy
between close trees (a relocated subtree flips many edges). Use as
baseline only, pair with MC.

---

## 3. Matching Cluster (MC) distance — Bogdanowicz & Giaro

**Paper.** Bogdanowicz D. & Giaro K. (2013), "On a matching distance between
rooted phylogenetic trees," *Int. J. Appl. Math. Comput. Sci.* 23(3):
669–684. DOI:
[10.2478/amcs-2013-0050](https://doi.org/10.2478/amcs-2013-0050).

**Construction.** For rooted binary T₁, T₂ on the same leaf set, build a
bipartite graph between the internal clusters of T₁ and T₂ with edge cost
= Hamming / symmetric-difference distance between the leaf subsets
(optionally branch-length-weighted). Solve min-weight perfect matching
(Hungarian). d_MC = matching cost. A relocated subtree costs only O(subtree
size), not "every affected edge" — the fix for wRF's flip problem.

**What it captures vs the 14.** Topology via which clusters match, AND
(weighted variant) heights via matched-pair cost. Crucially, grades
"almost-matching" clusters — exactly the regime where close dendrograms
differ by a few leaves migrating between sibling communities. Stage-1
ARI/VI see only one k at a time; MC fuses all internal partitions in
one scalar.

**Python implementation.** Partial. The canonical implementation is Java
(TreeCmp: github.com/TreeCmp/TreeCmp, Bogdanowicz, Giaro & Wróbel 2012).
In Python, `scipy.optimize.linear_sum_assignment` on a hand-built cost
matrix reproduces MC in ~40 LOC. No mature standalone package.

**Complexity.** O(N³) for Hungarian on N − 1 clusters. Milliseconds at
N ≈ 100.

**Caveats.** Branch-length-weighted MC is still a metric (2013 Thm. 4).
Unbounded — normalize by cohort diameter. MC is rooted (ours);
unrooted variant is Matching Split (inapplicable). Empirically agrees
with RF on identical topologies, degrades gracefully under subtree
relocation (2013 Fig. 3).

---

## 4. Interleaving distance on labelled merge trees

**Papers.** Morozov D., Beketayev K., Weber G. (2013), "Interleaving
Distance between Merge Trees," *TopoInVis 2013*. Labelled variant:
Gasparovic E. et al. (2019), [arXiv:1908.00063](https://arxiv.org/abs/1908.00063);
Yan L. et al. (2022), *IEEE TVCG*.

**Construction.** A UPGMA dendrogram **is** a merge tree (leaves at height
0, internal nodes at merge heights). d_I(T₁, T₂) = inf{ε : there exist
ε-maps T₁ → T₂[+ε] and T₂ → T₁[+ε] commuting up to a 2ε-shift} — the
smallest ε such that lifting all merge heights by ε makes the trees
equivalent as labelled rooted topologies. In the *labelled* variant each
leaf carries a fixed identity (our setting).

**What it captures vs the 14.** d_I = 0 iff labelled merge trees are
isomorphic with identical heights. d_I grows with *both* height
perturbation and topology change — a bona fide metric on the joint
(topology, height) object. Stage-1 scalars collapse these axes.

**Python implementation.** Partial.
- `tdavislab/MergeTreeMetric` — Python, labelled interleaving +
  bottleneck + Wasserstein; depends on TTK / ParaView (adaptor needed
  for scipy dendrograms).
- `trneedham/Decorated-Merge-Trees` — pure Python, related decorated
  merge trees.
- `gudhi.bottleneck_distance` on the PD of merge heights is a cheap
  *lower bound* on interleaving but **drops leaf labels** — two trees
  with identical heights and disjoint topologies score 0. Use only as
  a sanity check.

**Complexity.** Labelled interleaving is O(N³) or better under mild
structural assumptions; bottleneck on PDs is O(N² log N).

**Caveats.** "Interleaving distance" names subtly different objects in
different papers — verify the package enforces leaf labels.
Infrastructure overhead is the main practical cost.

---

## 5. Path-difference distance (Steel & Penny) — **subsumed, skip**

**Paper.** Steel M.A. & Penny D. (1993), "Distributions of tree comparison
metrics — some new results," *Syst. Biol.* 42(2): 126–141.

**Construction.** For each leaf pair (i, j) record p(i, j) = edges
(topological) or summed branch length (patristic) on the i–j path;
d_path = L2 of P-vector differences.

**Verdict: subsumed by KC, skip as a separate method.** Patristic
version ≈ Stage-1 L2-on-ultrametric; topological version = KC at λ = 0.
KC with intermediate λ strictly generalizes path-difference. Python:
`dendropy.calculate.treecompare.path_distance` (patristic).

---

## 6. Billera–Holmes–Vogtmann (BHV) geodesic — **not recommended here**

**Papers.** Billera L.J., Holmes S.P., Vogtmann K. (2001), "Geometry of the
space of phylogenetic trees," *Adv. Appl. Math.* 27: 733–767. Polynomial
algorithm: Owen M. & Provan J.S. (2011), "A Fast Algorithm for Computing
Geodesic Distances in Tree Space," *IEEE/ACM TCBB* 8(1): 2–13.

**Construction.** Embed each labelled tree as a point in a CAT(0)
piecewise-Euclidean orthant complex (coordinates = branch lengths); d_BHV
= length of the unique geodesic. Theoretically the gold standard for a
joint topology + length metric.

**Python.** No native. `Pathtrees` (Python) wraps the Java GTP binary.
Usable but dependency-heavy.

**Complexity.** O(N⁴) worst case; fine for N ≈ 100.

**Caveats.** For ultrametric trees BHV is dominated by height differences
and behaves like a smoothed wRF; rarely applied to dendrograms. Skip
unless KC + MC disagree structurally — the Java dependency is not worth
it for a pilot.

---

## Recommendation — run KC and MC first

**Start Stage 3 with the Kendall–Colijn metric and the Matching Cluster
distance.** Rationale:

- **KC directly fixes the Stage-1 gap.** Explicit topology (m) and
  height (M) components blended by a knob λ; reporting the cohort
  statistic at λ ∈ {0, 0.5, 1} cleanly separates whether
  cross-phase / cross-patient discriminability is driven by topology,
  heights, or their combination — a story the 14 Stage-1 scalars cannot
  tell. Short R-to-Python port, O(N²), easy to normalize.
- **MC complements KC** because it is not a pair-path quantity. It
  operates on induced clusters via bipartite matching and degrades
  gracefully under subtree relocation — the failure mode that makes wRF
  noisy between close dendrograms. The branch-length-weighted variant is
  a genuine topology+height scalar and pairs with KC as a robustness
  check.

Keep **wRF/BSD** as a 5-minute DendroPy baseline (expect it to be noisy).
Hold **labelled interleaving** in reserve for Stage 4 — promote it if
KC + MC motivate a merge-tree-metric analysis; the GUDHI bottleneck is a
cheap lower bound but loses labels. **Skip** path-difference (subsumed
by KC) and BHV (dependency cost > payoff).
