---
type: report
status: current
era: IMCOH_ABS / COHORT_N10
date: 2026-05-08
audience: writing-agent (§5.6 polish)
sources:
  - scripts/01_compute/audit/audit_47_kc_trace_network_view.py
  - scripts/01_compute/audit/audit_48_kc_reset_module_view.py
  - scripts/01_compute/audit/audit_49_kc_rearrangement_module_view.py
  - scripts/01_compute/audit/audit_50_kc_anchor_module_view.py
  - scripts/01_compute/audit/audit_51b_taxonomy_pat07_beta_composite.py
  - scripts/01_compute/audit/audit_52_section5_class_showcases.py
  - data/audit/kc_trace_network_view/trace_subtrees_summary.csv
  - data/audit/kc_reset_module_view/reset_subtrees_summary.csv
  - data/audit/kc_rearrangement_module_view/rearrangement_subtrees_summary.csv
  - data/audit/kc_anchor_module_view/anchor_subtrees_summary.csv
---

# §5.6 verification — taxonomy definitions + cohort-ranking phrasings

**Head.** Definitional gaps A1-A5 close cleanly from the four `audit_4{7,8,9}` /
`audit_50` scripts and `audit_51b` / `audit_52` (which define strict-intersection
shading); the matching rule is **greedy argmax Jaccard with a 1.5× size band**,
not bipartite. The writing agent's strict-intersection reading for trace / reset
/ anchor needs **one fix**: the strict set is `leaves_tt ∩ leaves_post` for
trace (and the analogous symmetric pair for reset), **not** the three-way
construction over rspre-fragments. On B/C: **B2 holds** (Pat_07 β is the only
cell with `min(n_trace, n_anchor) = 6`). **B1 overclaims** — Pat_10 low-γ ties
Pat_13 at 5 reset clades and is *beaten* on reset-leaf count by Pat_13 (18 vs
29 / 32). **C overclaims locally** — Pat_03 high-γ has 12 rearrange clades but
is **tied** with Pat_02 (12) and Pat_06 (12); the global "most abundant
rearrange class" claim is fine (rearrange 560 clades vs trace 107 + reset 83 +
anchor 144 across 60 cells).

---

## A. Definitional gaps in the matching rule

All four module audits share the same scaffolding (same `all_subtrees`,
same FRAG_FACTOR, same MIN_SIZE, same MAX_OVERLAP_FRAC, same MAX_*_TOTAL).
Numbers below cite `audit_47_kc_trace_network_view.py` lines unless
otherwise noted, but they hold verbatim for `audit_48` / `audit_49` /
`audit_50` (modulo the per-class Jaccard gates).

### A1. Subtree pool — option (i) with a size window

**Answer: option (i), all internal nodes, with `MIN_SIZE = 3` and
`max_size = max(MIN_SIZE+1, n // 2)`.**

`all_subtrees(Z, n, min_size, max_size)` (audit_47:93–104) walks all
N − 1 internal nodes of the linkage `Z`, building each subtree's leaf
set, then keeps every internal node whose leaf-set size lies in
`[min_size, max_size]`. There is **no cut** on the dendrogram. The
caller passes `min_size=3` (`MIN_SIZE`, audit_47:80) and
`max_size = max(MIN_SIZE+1, n // 2)` (audit_47:124–125), which for our
typical `n ≈ 110–120` cohort gives `max_size ≈ 55–60`.

So for a phase tree with `n ≈ 115`, the candidate pool is roughly all
internal nodes with leaf-set size between 3 and 57 — about
`O(n) − ε` candidates per phase.

### A2. Closest counterpart — option (i), greedy + size band

**Answer: option (i), greedy `argmax Jaccard` over the destination
pool, restricted to a size band `[s/1.5, s × 1.5]`.**

For each source subtree of size `s`, the destination is sought by
plain `argmax_T J(L_S, L_T)` (audit_47:140–151). There is **no
bipartite matching**, **no removal of used subtrees**, and **no
two-pass refinement** — same source can in principle be matched to
the same destination as another source (the deduplication kicks in
only at selection time, via `MAX_OVERLAP_FRAC`).

Two side conditions on the destination search:

1. **Size band.** `size_lo_post = max(min_size, int(s/1.5))`,
   `size_hi_post = min(max_size, int(s*1.5))` (audit_47:135–136).
   Destinations outside this band are skipped before the Jaccard test.
2. **Best-so-far loop.** Only destinations that **strictly improve**
   the running best `j_post` are stored (audit_47:149–151), so ties
   resolve to the **first encountered** subtree (i.e. by the
   linkage-row order).

**Anchor (audit_50)** does the same per source `L_pre` but searches
**two destination pools simultaneously** (`L_tt` and `L_post`) under
the same 1.5× size band, and selects the `(L_tt, L_post)` pair that
maximizes the geometric mean `(J_pre,tt · J_pre,post · J_tt,post)^(1/3)`
(audit_50:137–172). Each constituent Jaccard must individually be
≥ 0.6.

**Reset (audit_48)** is the symmetric inverse of trace: source = rPre
subtrees; destination = rPost subtrees with the same 1.5× size band.

**Rearrange (audit_49)** is the only audit with no destination
search — every gate is on the rPost source itself (rPre fragment, tt
fragment, plus containment guards).

**Selection-time deduplication (all four audits).** After the
candidate list is built, candidates are sorted by `(-size, -score)`
(audit_47:201) and selected greedily, accepting a candidate only
if its **union leaf set** overlaps the cumulative `used` set by less
than `MAX_OVERLAP_FRAC = 0.30` (audit_47:202–211). Selection stops
at `MAX_TRACE_TOTAL = 12` (or `MAX_RESET_TOTAL`,
`MAX_REARRANGE_TOTAL`, `MAX_ANCHOR_TOTAL`, all 12).

### A3. Jaccard gray zone [0.5, 0.6) — both gates strict; rejection (not "fragmented")

**Answer: a candidate whose key Jaccard lands in `[0.5, 0.6)` is
rejected entirely from the trace-eligible set, on either side. The
gray zone is a deliberate buffer that excludes ambiguous matches; it
does not get classified as "fragmented".**

Trace gates (audit_47:153, 172):

```
if best_j_post < J_POST_MIN (= 0.6):  continue   # not coherent enough → reject
if j_pre        >= J_PRE_MAX (= 0.5): continue   # not fragmented enough → reject
```

So the buffer is asymmetric and the gray zone is between the two
strict thresholds:

| `j_post` (tt-rPost match)  | label                  | outcome          |
|----------------------------|------------------------|------------------|
| `≥ 0.6`                    | coherent               | continue         |
| `[0.0, 0.6)`               | not-coherent-enough    | **reject**       |

| `j_pre` (best rPre match)  | label                  | outcome          |
|----------------------------|------------------------|------------------|
| `< 0.5`                    | fragmented             | continue         |
| `[0.5, 1.0]`               | not-fragmented-enough  | **reject**       |

A candidate with `j_pre ∈ [0.5, 0.6)` is **NOT** declared fragmented.
It is **rejected** — the leaves are co-clustered too tightly in rPre
to be called fragmented but not tightly enough to be called an
anchor; the trace pipeline has no class for them and they fall
through to "diffuse" downstream.

The §5.6 prose phrasing should therefore be: "coherent if Jaccard
≥ 0.6, fragmented if Jaccard < 0.5, **rejected** if in between" — not
"fragmented if below 0.5" (which silently invents a fragmented-side
buffer of 0.1).

### A4. FRAG_FACTOR = 2 — smallest-containing-subtree containment guard

**Answer: a source subtree of size `s` is rejected if the smallest
disrupted-phase subtree fully containing all of `L_S` has size
`< FRAG_FACTOR × s = 2 × s`.**

Code (audit_47:179–188):

```python
smallest_containing = float("inf")
for st_pre in sts_pre:
    if L_tt <= st_pre["leaves"]:
        if st_pre["size"] < smallest_containing:
            smallest_containing = st_pre["size"]
if smallest_containing < FRAG_FACTOR * s:
    continue
```

So if `|L_tt| = 4` and the smallest rPre subtree fully containing all
4 leaves has size 5 or 6 or 7, the candidate is rejected
(`5 < 2 × 4`, `6 < 8`, `7 < 8`). The candidate survives only if the
smallest containing rPre subtree is at least size 8 — i.e. the leaves
are co-located inside a substantially larger rPre block. The intuition
is that even when leaf-set Jaccard is low (< 0.5), if the leaves are
all packed inside a tight rPre subtree only slightly bigger than `s`,
they are anchored, not dispersed.

The same containment guard appears in `audit_48` (Gate 3, "tt
containment fragmentation", line 175), `audit_49` Gates 2 and 4 (rPre
containment AND tt containment), and is **absent** from `audit_50`
(anchor) since anchor's coherent-everywhere definition has nothing to
fragment against. `FRAG_FACTOR = 2.0` is the same constant in all
three.

### A5. Strict-intersection vs union — code-level definitions

**Answer to the writing agent's reading: needs ONE FIX.** The strict
intersection for trace is **only over the matched coherent pair**
`leaves_tt ∩ leaves_post`, not a three-way construction across
rspre-fragments. The same correction applies to reset symmetrically.

#### A5.1 Strict-intersection (used for D̂ block reordering, dendrogram leaf-bar shading, and module-membership counts)

From `audit_51b:102–111` and `audit_52:166–175` (verbatim, identical):

```python
def strict_module_leaves(cls, module):
    if cls == "trace":
        return module["leaves_tt"] & module["leaves_post"]
    if cls == "reset":
        return module["leaves_pre"] & module["leaves_post"]
    if cls == "anchor":
        return module["leaves_pre"] & module["leaves_tt"] & module["leaves_post"]
    if cls == "diffuse":
        return module["leaves"]
    return module.get("leaves_post_residual", module["leaves_post"])  # rearrange
```

Per class:

| class      | strict-intersection leaf set                                                        |
|------------|--------------------------------------------------------------------------------------|
| trace      | `L_tt ∩ L_post`                                                                      |
| reset      | `L_pre ∩ L_post`                                                                     |
| anchor     | `L_pre ∩ L_tt ∩ L_post`                                                              |
| rearrange  | `L_post  −  ( strict_trace_leaves  ∪  strict_reset_leaves  ∪  strict_anchor_leaves )`, kept only if residual size ≥ 3 |
| diffuse    | `range(n)  −  ( all classified strict leaves )`                                      |

The rearrange rule is built in `assemble_class_modules`
(audit_52:188–198 / audit_51b:213–222): for each rearrangement
candidate `m`, `residual = m["leaves_post"] − non_rearrange_leaves`,
keep iff `len(residual) ≥ 3`. Diffuse is whatever leaves remain
unclassified after the previous four classes pass.

**Where the writing agent's reading is wrong.** The agent wrote
trace as "leaves in (rspre fragments) AND in matched rspre-taskt
subtree AND in matched taskt-rspost subtree". There is **no matched
rspre-taskt subtree** in the trace definition — rPre's role is to
fragment, not to provide a matched leaf set. The strict-trace set is
just `L_tt ∩ L_post`. Likewise reset is `L_pre ∩ L_post` only, and
anchor is the three-way intersection.

Rearrange and diffuse readings are correct in spirit (rPost-anchored
residual; everything else); the formal rule keeps a clade-size floor
of 3 on the residual (residuals smaller than that drop the candidate,
and those leaves end up in `diffuse` instead of `rearrange`).

#### A5.2 Union convention (network panels)

From `audit_47:241–246`, `audit_48:245–246`, `audit_49:197–198`,
`audit_50:230–231`:

| class      | network-panel leaf set                          |
|------------|-------------------------------------------------|
| trace      | `L_tt ∪ L_post`                                 |
| reset      | `L_pre ∪ L_post`                                |
| anchor     | `L_pre ∪ L_tt ∪ L_post`                         |
| rearrange  | `L_post`  (after residual subtraction in audit_52) |
| diffuse    | the strict diffuse leaf set (no further union)  |

So the writing agent's reading for the network panels is right: the
network shows the **union** of leaves across the matched phase pairs,
not the strict intersection. This is the per-class audit's canonical
visual convention; `audit_52` (the §5.6 single-class showcase script)
delegates to each per-class `draw_network` so this convention is
inherited unchanged.

**One subtlety on rearrange networks.** In `audit_52`, the
"`leaves_post`" that the network panel reads is overwritten with
`leaves_post_residual` (audit_52:382–386), so the rearrange network
panel shows the **residual** rPost leaves (after subtracting strict
trace ∪ reset ∪ anchor leaves), not the raw rearrangement-candidate
leaves. This is why the §5.6 rearrange network is visually less
populated than the audit_49 rearrange network for the same cell.

---

## B. Backing numbers for two cohort-ranking phrasings

### B1. "Cohort's reset specialist" — Pat_10 low-γ — **OVERCLAIMS**

**Pat_10 ties Pat_13 at 5 reset clades and is beaten by Pat_13 on
reset-leaf count.** "Specialist" should be retired or weakened to
"tied for cohort top by reset-clade count".

Per-patient reset (λ=0) at low_gamma, full table:

| patient | n_clades | sizes_pre        | n_leaves_pre | sizes_post       | n_leaves_post |
|---------|----------|------------------|--------------|------------------|---------------|
| Pat_02  | 2        | [9, 49]          | 58           | [8, 55]          | 63            |
| Pat_03  | 3        | [3, 4, 6]        | 13           | [4, 4, 4]        | 12            |
| Pat_05  | 1        | [3]              | 3            | [3]              | 3             |
| Pat_06  | 1        | [3]              | 3            | [3]              | 3             |
| Pat_07  | 1        | [3]              | 3            | [3]              | 3             |
| Pat_08  | 3        | [4, 4, 8]        | 16           | [3, 4, 11]       | 18            |
| **Pat_10** | **5** | [3, 3, 3, 4, 5]  | **18**       | [3, 3, 3, 4, 5]  | **18**        |
| **Pat_13** | **5** | [3, 4, 4, 5, 13] | **29**       | [3, 4, 4, 6, 15] | **32**        |
| Pat_14  | 2        | [3, 3]           | 6            | [3, 3]           | 6             |
| Pat_15  | 3        | [3, 4, 5]        | 12           | [4, 4, 6]        | 14            |

Sorted descending by clade count:

| rank | patient | n_clades | n_leaves_pre |
|------|---------|----------|--------------|
| 1    | Pat_13  | 5        | 29           |
| 1    | Pat_10  | 5        | 18           |
| 3    | Pat_03  | 3        | 13           |
| 3    | Pat_08  | 3        | 16           |
| 3    | Pat_15  | 3        | 12           |

The §5.6 prose claim was "Pat_10 carries five reset clades, 17 reset
leaves of 113 total" → on cohort `_summary.csv` the leaf count is 18,
not 17 (minor); but the bigger issue is that Pat_13 has the same
clade count and **substantially more reset leaves** (29 / 32 vs Pat_10's
18 / 18), with one large 13-leaf reset clade that Pat_10 lacks. So
Pat_10 is *not* the cohort top in any simple sense — it ties on
clade count and loses on leaf count.

**Recommended phrasing fix.** Either drop the "specialist" framing
("Pat_10 carries the cohort's largest equal-rank reset profile at low-γ")
or pivot the showcase to Pat_13 low-γ (which is the unambiguous
reset cohort top). Pat_10 low-γ remains a fine reset showcase
visually; the prose just shouldn't claim it's the standalone winner.

### B2. "Cohort's joint-richest trace-and-anchor profile" — Pat_07 β — **HOLDS**

**Pat_07 β is uniquely the only cell with both `n_trace ≥ 6` and
`n_anchor ≥ 6`, and is the cohort top under `min(n_trace, n_anchor)`.**

Per-patient (n_trace, n_anchor) at β, λ=0:

| patient | n_trace | n_anchor | min(t,a) | sum(t,a) |
|---------|---------|----------|----------|----------|
| **Pat_07** | **6** | **6**  | **6**    | **12**   |
| Pat_02  | 4       | 9        | 4        | 13       |
| Pat_03  | 3       | 3        | 3        | 6        |
| Pat_13  | 2       | 7        | 2        | 9        |
| Pat_06  | 2       | 3        | 2        | 5        |
| Pat_14  | 3       | 2        | 2        | 5        |
| Pat_10  | 1       | 6        | 1        | 7        |
| Pat_05  | 1       | 5        | 1        | 6        |
| Pat_15  | 1       | 1        | 1        | 2        |
| Pat_08  | 0       | 1        | 0        | 1        |

Pat_07 is uniquely cohort top by `min(n_trace, n_anchor) = 6`. Pat_02
has higher `sum` (13) but the imbalance (4 trace vs 9 anchor) is
exactly what "joint-richest" rules out — Pat_02 is anchor-heavy.
The phrasing "six multiscale trace modules and six anchor modules
of comparable count, expressing the cohort's joint-richest
trace-and-anchor profile" is supported.

### Note on a prose typo

Earlier the agent's brief mentioned "five reset clades **(17 reset leaves
of 113 total)**" — the cached summary records **18** reset leaves at
Pat_10 low-γ on each side (rPre and rPost). 17 might come from a
strict-intersection count (`leaves_pre ∩ leaves_post`) rather than the
union or one-sided sum; if §5.6 is reporting the strict-intersection
count, the number is plausible but the cited table column should be
named accordingly.

---

## C. Pat_03 high-γ as rearrange showcase — **OVERCLAIMS LOCALLY, GLOBAL CLAIM HOLDS**

### C.1 Pat_03 high-γ: rearrange counts

Pat_03 high-γ has **12 rearrangement clades, 70 rearrangement leaves**:

```
sizes_post = [4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 20]
```

### C.2 Per-patient rearrangement clade count at high-γ — Pat_03 is **tied**, not unique top

| patient | n_clades | n_leaves | sizes_post (largest few)                |
|---------|----------|----------|------------------------------------------|
| **Pat_02** | **12** | 95     | [..., 7, 7, 21, 26]                      |
| **Pat_03** | **12** | 70     | [..., 5, 5, 6, 20]                       |
| **Pat_06** | **12** | 61     | [..., 5, 6, 10, 12]                      |
| Pat_05  | 10       | 69       | [..., 4, 5, 5, 35]                       |
| Pat_14  | 10       | 99       | [..., 7, 8, 8, 9, 52]                    |
| Pat_10  | 9        | 85       | [..., 3, 4, 6, 7, 53]                    |
| Pat_13  | 9        | 89       | [..., 3, 4, 5, 20, 45]                   |
| Pat_15  | 7        | 25       | [..., 3, 3, 4, 6]                        |
| Pat_07  | 4        | 20       | [3, 3, 5, 9]                             |
| Pat_08  | 4        | 12       | [3, 3, 3, 3]                             |

Pat_03 is **tied for cohort top** at 12 rearrangement clades along
with Pat_02 and Pat_06, and is ranked **3rd** by leaf count (70, vs
Pat_02 95 and Pat_14 99). The §5.6 prose claim that Pat_03 high-γ
has "the cleanest rearrange-dominant high-γ visual" is plausible
*visually* (large 20-leaf clade dominates without the noise of
Pat_14's 52-leaf super-clade or Pat_02's 21+26 pair), but the
cohort-rank framing should not present Pat_03 as the standalone
winner. **Recommended phrasing:** "Pat_03 high-γ — one of the cohort's
top-ranked rearrange cells (tied at 12 clades with Pat_02 and Pat_06)
— gives the cleanest single-clade dominant visual".

### C.3 Rearrange is the most abundant cross-phase class — **HOLDS**

Per-patient rearrange-clade total summed across the 6 bands (after
the audit_49 `MIN_SIZE = 3` and `MAX_OVERLAP_FRAC` filters):

| patient | rearrange total clades | rearrange total leaves |
|---------|------------------------|------------------------|
| Pat_06  | 70                     | 394                    |
| Pat_03  | 65                     | 476                    |
| Pat_15  | 65                     | 312                    |
| Pat_02  | 57                     | 504                    |
| Pat_14  | 56                     | 511                    |
| Pat_05  | 54                     | 487                    |
| Pat_13  | 52                     | 298                    |
| Pat_08  | 50                     | 411                    |
| Pat_10  | 47                     | 435                    |
| Pat_07  | 44                     | 395                    |

Cohort grand totals (λ=0 for trace / reset / anchor; raw audit_49
output for rearrange) over 60 (10 patients × 6 bands) cells:

| class      | total clades | per-cell mean |
|------------|--------------|---------------|
| trace      | 107          | 1.78          |
| reset      | 83           | 1.38          |
| anchor     | 144          | 2.40          |
| rearrange  | **560**      | **9.33**      |

Rearrange is **3.9× the next-most-abundant class (anchor)** and
**5.2× trace** at the per-cell mean. The §5.6 prose claim that
rearrangement is "the most abundant cross-phase module class across
the cohort" is well supported.

(Note: this is the audit_49 raw output. After the §5.6 residual
subtraction in `audit_52` / `audit_51b`, some rearrangement clades
collapse below the size-3 floor and migrate to diffuse; the residual
counts are smaller but rearrangement is still the dominant class.
The grand totals above use the raw audit_49 counts because that's
what the per-class summary CSVs hold.)

---

## D. Pat_13 β composite — shade-tracking convention

**Confirmed: shade is keyed to strict-intersection leaf set identity;
the same module receives the same shade in every column; in columns
where the module's matched subtree at that phase is fragmented or
absent, the shaded leaves still appear at their actual dendrogram
positions (just no longer connected by colored internal links).**

`audit_51b:114–115` defines `class_leaves_for_phase` as a thin alias
to `strict_module_leaves` — phase-agnostic. The leaf-bar shading
loop (`audit_51b:184–191`) draws each leaf at its actual horizontal
position with the module's `_shade`, where `_shade` is assigned
once per module by `shades_for_class(cls, n)` (audit_51b:75–81 /
audit_52:178–184) on the largest-strict-intersection-first
ordering. So:

- Module identity is fixed across columns (the same dict object is
  reused).
- The strict-intersection leaf set is fixed across columns
  (phase-agnostic by definition).
- The shade is keyed to module identity, hence stable across columns.
- The internal-link color (`make_link_color_func`,
  audit_51b:128–155) requires `subtree_leaves ⊆ strict_module_leaves`
  to inherit the module shade; in phases where the strict set is
  fragmented across the dendrogram, no large internal subtree
  satisfies this containment and the colored region appears as
  isolated leaf bars + small internal links, exactly as the writing
  agent describes.

The proposed caption sentence is correct as written; one suggested
tightening: replace "the same module receives the same shade across
all three phase columns" with "the same strict-intersection leaf set
receives the same shade across all three phase columns" — this makes
clear what "same module" means operationally (the strict set, not the
matched-subtree leaf set, which differs per phase).

---

## Notes for prose polish

- A3: the "fragmented if Jaccard < 0.5" phrasing in §5.6 is too loose
  — the actual rule is **strict-rejection** in the gray zone, not
  "fragmented".
- A5: trace strict-intersection is **two-way** (`L_tt ∩ L_post`), not
  three-way over "rspre-fragments". Same structure for reset.
- B1: drop "cohort's reset specialist" or weaken to "tied-top"; if a
  unique winner is wanted, pivot to **Pat_13 low-γ** (5 clades,
  29 / 32 leaves).
- B2: holds as written — Pat_07 β is the only cell with
  `min(n_trace, n_anchor) = 6`.
- C: drop the "cleanest rearrange-dominant high-γ visual" implication
  that Pat_03 leads on count; keep the visual claim, but caveat the
  cohort-rank position (tied at 12 with Pat_02 and Pat_06; 3rd by
  leaf count).
- C3 global "most abundant" claim holds (560 vs 144 anchor / 107
  trace / 83 reset across 60 cells).
