---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
scope: the bare label "m=48" is used for two unrelated multiple-comparison families; rename or scope-tag every occurrence
---

# Bonferroni / BH `m=48` overload

**Head.** The label `m=48` appears in §5.3 and §5.5 referring to **two unrelated** multiple-comparison families. §5.3 `m=48` is the joint LRG-probe family of (KC × 2 axes + CTM + 3 matrix distances + VI(k) + Grassmann) × 6 bands used after the within-probe per-band BH at m=6. §5.5 `m=48` is the (band, region) anatomy cell space passing the strict-eligibility filter (n_contacts ≥ 5 AND n_pat_trace ≥ 2 AND enrichment > 1). These are different probabilistic spaces, different statistics (continuous-statistic Wilcoxon vs count-statistic hypergeometric), and different rejection thresholds (BH-q vs Bonferroni α/m). The manuscript uses the bare label `m=48` for both. Rename to **`m=48^{(LRG)}`** and **`m=48^{(anat)}`** (or spell out scope on every mention) to prevent silent conflation.

## The two families, exactly

### §5.3 `m=48` — LRG joint-BH probe family

**Composition.** 48 cells assembled from five LRG-probe sub-families × 6 bands:

| sub-family | cells | construction |
|---|---:|---|
| KC at λ=0 (topology) and λ=1 (heights) | 12 | 6 bands × 2 KC parameter values |
| CTM | 6 | 6 bands × 1 split-vs-drift Wilcoxon |
| Matrix distances on `D(τ)` | 18 | 6 bands × {`d_S`, `d_P`, `d_F`} |
| VI(k) summary | 6 | 6 bands × 1 best-k summary scalar |
| Grassmann summary | 6 | 6 bands × 1 best-k summary scalar |
| **TOTAL** | **48** | (matches §5.3 prose breakdown) |

**Statistic.** Each cell is a **continuous-statistic** one-sided Wilcoxon signed-rank test on a per-patient probe scalar (`T_KC`, `T_CTM`, `T_VI`, `T_d_*`, `d_G`).

**Rejection rule.** **BH-FDR** (Benjamini-Hochberg) at q = 0.05.

**Outcome at m=48.** Smallest joint q = **0.164** at four cells (CTM α p=0.0068, CTM low-γ p=0.0098, CTM β p=0.0137, Grassmann low-γ k=13 p=0.0137). **Zero cells survive q ≤ 0.05.**

**Cross-reference.** Verified in `.agents/reports/2026-05-08_section5_joint_bh.md` (Family B, with alternative families A at m=84 and C at m=36 also reported). The manuscript currently quotes "the smallest joint q ≤ 0.05" outcome for this family but the headline q value cited in the §5.3 prose ("smallest joint at q ≤ 0.05; the smallest CTM trace bands") is at the m=48 family level.

### §5.5 `m=48` — anatomy Bonferroni cell space

**Composition.** 48 (band, region) cells passing the strict-eligibility filter applied to the cohort cortical contact pool (`N_total = 866`):

- Filter: `n_contacts(r) ≥ 5` AND `n_pat_trace(r) ≥ 2` AND `enrichment(r) > 1`.
- Each band yields ≈ 8 strict-eligible regions (the audit_50 `tail(8)` cap on the bar chart for visualization, but actually the strict-eligible regions per band).
- Sum over 6 bands = 48 cells.
- Trace-band breakdown: α has 5 strict-eligible cells, β has 5, γ_l has 7 (= 17 cells in the trace bands; the rest of the m=48 budget is consumed by δ / θ / γ_h non-trace bands).

**Statistic.** Each cell is a **count-statistic** one-sided hypergeometric tail test:
`p_hyper(r) = P(X ≥ n_trace(r) | N_total, K, n_contacts(r))`,
where `K` is the cohort baseline trace-leaf count for the band and `n_trace(r)` is the observed trace-leaf count in region r.

**Rejection rule.** **Bonferroni** at α/m = 0.05 / 48 ≈ **1.04 × 10^−3**.

**Outcome at m=48.** **Exactly one** cell survives: γ_l × ctx-lh-fusiform (13/38, 3.95× cohort baseline, p_hyper = 5.2 × 10^−6).

**Cross-reference.** Verified in `.agents/reports/2026-05-08_section5_5_verify.md` (full-cohort + Pat_03 dropout + pro-cohort sensitivity).

## Why these are not the same family

1. **Different probabilistic objects.** §5.3 cells test continuous probe scalars (Wilcoxon on per-patient `T_*`). §5.5 cells test count statistics (hypergeometric on per-region `n_trace`).
2. **Different correction policies.** §5.3 uses BH-FDR (false discovery rate control). §5.5 uses Bonferroni (familywise error rate control). The two policies have different operational meanings and different conservativeness.
3. **Different per-cell semantics.** §5.3 cells say "is this probe / band combination directionally consistent across patients?" §5.5 cells say "does this anatomical region carry an over-representation of trace-leaves in this band?"
4. **Different units of analysis.** §5.3 sums over (probe, band) combinations and aggregates over patients within each cell. §5.5 sums over (band, region) and aggregates over leaves within each cell.

Mixing the two would (i) double-count the β / low-γ cells that appear in both families at the cohort level (the §5.3 β cell and the §5.5 β anatomy cells both reference the β trace direction), and (ii) merge a continuous-statistic family with a count-statistic family without a coherent joint distribution.

## Occurrences of "m=48" / "m = 48" in the manuscript prose

The following are the explicit occurrences of the label as observed in `notes_imcoh_260508.pdf`. Line refs are page-relative — writing agent should grep the source `.tex` to enumerate the exact `.tex` line numbers.

### §5.3 occurrences

- §5.3 closing paragraph (`p.33` of PDF):
  > "Under joint BH across the **forty-eight** (probe, band) cells of the five LRG probe families assembled in this section — KC at `λ ∈ {0, 1}`, CTM, the three matrix distances on `D(τ)`, VI(k), and the Grassmann probe at one summary scalar per band — no single cell survives at `q ≤ 0.05`; the smallest joint `q` is 0.164, reached by the three CTM trace bands."

  This is the §5.3 LRG joint-BH family.

### §5.5 occurrences

- §5.5 anatomy paragraph (`p.37` of PDF):
  > "The multiple-comparison correction is Bonferroni across the `m = 48` (band, region) cells passing the strict-eligibility filter, with corrected threshold `α/m ≈ 1.04 × 10^−3`."

  This is the §5.5 anatomy Bonferroni family.

- §5.5 anatomy headline (`p.37` of PDF):
  > "the only (band, region) cell to clear Bonferroni correction across the `m = 48` cohort-eligible cells (`p_hyper = 5.2 × 10^−6`)".

  Same §5.5 family.

### §1 introductory preview

- §1 closing paragraph (`p.3` of PDF):
  > "an anatomical headline at low-γ ctx-lh-fusiform that survives Bonferroni correction across the (band, region) cell space"

  This refers to the §5.5 family. The label `m=48` is not numeric here; the label "Bonferroni correction across the (band, region) cell space" is unambiguous about the family.

### §6.1 synthesis

- §6.1 (`p.48` of PDF):
  > "The single-region anatomical address that distinguishes low-γ from β sits at left fusiform cortex, the only (band, region) cell to clear Bonferroni correction across the `m = 48` cohort-eligible cells".

  This is the §5.5 family.

## Recommendation — writing agent

Apply one of the following two patterns globally:

### Pattern A — explicit superscript (recommended)

Rename in prose:
- `m=48^{(LRG)}` for the §5.3 joint-BH family.
- `m=48^{(anat)}` for the §5.5 anatomy Bonferroni family.

Apply on first use in each subsection; subsequent uses in the same subsection can drop the superscript. Cross-section back-references (e.g. §6.1 referring to §5.5 anatomy) should always include the superscript or full-phrase scope.

### Pattern B — full phrase scope (alternative)

Always spell out scope:
- "Bonferroni correction across the `m = 48` (band, region) **anatomy** cells" — for §5.5.
- "BH-FDR across the `m = 48` (probe, band) **LRG-probe** cells" — for §5.3.

Pattern B is more verbose but reads cleanly without a typeset superscript convention.

## Why this is non-trivial

A reader who skims the §6.1 synthesis and sees "Bonferroni `m = 48`" without scope-tagging is likely to associate it with whichever m=48 they last encountered (§5.5 if they're reading sequentially; §5.3 if they're reading as a citation chain from a per-pair claim). The two families have **opposite outcomes** at q/α = 0.05 (§5.3 has zero survivors, §5.5 has one survivor). Conflation here flips the headline.

The §5.3 family also has an alternative-family choice (`m = 84` or `m = 36`, see `2026-05-08_section5_joint_bh.md`); the writing agent should check that the §5.3 prose breakdown (`12 + 6 + 18 + 6 + 6 = 48`) is the chosen family and the cited smallest q is the one corresponding to that family (`q_min = 0.164` at m=48). The current PDF prose appears consistent with Family B at m=48 → q_min = 0.164.

## Action — writing agent

1. Decide between Pattern A (superscript) and Pattern B (full-phrase scope).
2. Apply globally.
3. Cross-check that the §5.3 cited `q_min` matches the family-arithmetic of `2026-05-08_section5_joint_bh.md` (Family B → q_min = 0.164).
4. Cross-check that the §5.5 cited `p_hyper = 5.2 × 10^−6` and `α/m ≈ 1.04 × 10^−3` reproduce against `2026-05-08_section5_5_verify.md` (full-cohort 13/38 cell, 3.95× enrichment).
5. The §6.1 back-reference to §5.5 anatomy must explicitly say "anatomy m=48" (or `m=48^{(anat)}`), since the §6.1 paragraph immediately precedes a discussion of "every layer" of §5 which includes §5.3 (also m=48).

## Action — coding agent

No coding changes required — the m=48 family arithmetic is already in audit_round2_section5_v2.py (LRG joint) and audit_50_anatomy_55_figure.py (anatomy Bonferroni). What's needed is a one-line comment header in each script noting "this is the LRG-probe m=48 family" or "this is the anatomy m=48 family" so future readers don't conflate them.
