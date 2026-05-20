---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
section: 6
---

# Section 6 — fix list

**Head.** §6 nails the band-resolved reading correctly: β multifacet (every layer
fires), low-γ multiscale-with-heights-axis-only (KC λ=1 + matrix d_F + VI(k)
mid-coarse + Grassmann broad-positive-diagonal), α per-pair-only with mixed
parsopercularis disclosure, δ edge-only-and-LRG-washout, θ drift-only across both
substrate and LRG, high-γ borderline. Three issues need writing-level action:
(i) some inherited "almost universal in β" phrasing in earlier drafts is too
strong for 7-8/10 — current §6.1 prose uses measured language and should be
preserved against any future tightening; (ii) the "deferred controls" paragraph
(LRG-layer drift-triangle null and symmetric cross-baseline at the diffusion
layer) is one paragraph in §6.3 but is load-bearing — these gate the persistence
language; (iii) §6.3 Outlook's per-patient anatomical mapping of trace vs anchor
vs reset should explicitly reference the anchor anatomy-baseline question
(§5.6 caveat) and the multiscale taxonomy question (τ-sweep).

## Confirmed (no action)

- **§6.1 β headline**: "**The β band is the cohort's multifacet imprint. It registers at every geometric level the analysis tests: the per-pair correlation on D(τ) (ρ_split = +0.222, controlled p = 0.014, q = 0.027), the tree topology and tree heights of the dendrogram (within-probe BH-q = 0.006 on both KC axes, 10/10 patients below their split-baseline null), the matrix Pearson and Frobenius distances on D(τ) (p = 0.042 on each), the partition-cut view across the dendrogram continuum (sustained 7-9/10 cohort agreement on VI(k) for k ∈ [17, 32]), and the leading-eigenmode subspace (positive principal-angle diagonal on the Grassmann probe sustained from k ≈ 12 out to k = 79, above a fine-k cold zone)**". All numbers reproduce against the §5 verifier reports. ✓
- **§6.1 low-γ headline**: per-pair p = 0.010, q = 0.027 (controlled); 9/10 patients above their drift floor (NOT 7/10 positive — this distinction is correctly drawn here per `_ctm_lowgamma_counts.md`); tree heights only (KC λ=1, 8/10, p=0.042, within-probe BH q = 0.055, "directional companion"); matrix Frobenius (p=0.053); intermediate-coarse VI(k) at k ∈ [45, 56]; broad principal-angle diagonal. **Single-region anatomical headline**: γ_l ctx-lh-fusiform Bonferroni-survived (p_hyper = 5.2×10⁻⁶ under m=48, "with the disclosed double dependence on Pat_02"). ✓
- **§6.1 α headline**: "**The α band is the geometrically narrowest of the three controlled bands. It registers at the per-pair correlation alone (controlled p = 0.007, q = 0.027) and is silent on the global tree distance at both KC axes (6/10 patients on each, neither approaching significance). It surfaces in the multiscale companions at fine partition cuts (VI(k) peak 9/10 at k=3) and in the spectral subspace diagonal above k ≈ 10, and admits one suggestive anatomical cell at left parsopercularis with two-patient base aggregating one anti-aligned and one pro-aligned CTM contribution**". ✓ Matches §5.5 Pat_07/Pat_14 mixed-CTM disclosure.
- **§6.1 δ / θ / high-γ negative bands**:
  - δ: edge-only, propagator washout — consistent with `_ctm_lowgamma_counts.md` δ ρ̃_split = +0.031, p = 0.246.
  - θ: substrate drift-only AND LRG drift-only.
  - high-γ: borderline at both layers; the reduced edge-weight dynamic range in this band is a documented substrate caveat (§4.6).
  - "**The Pat_07 and Pat_15 cross-probe-aggregate anti-alignment is preserved by every LRG probe and is read as a scientifically meaningful structural finding rather than as cohort noise, consistent with the implant-geometry interpretation introduced at §4.4.**" ✓
- **§6.2 diffusion-geometry enrichment statement**: "**The α band is the cleanest illustration of the enrichment: its substrate trace is the strongest of the cohort at the directional level (8/10 patients on d_S) but does not register at the dendrogram tree distance under controls; only the per-pair geometry of D(τ) crystallizes the cohort signal into a controlled claim. The α memory imprint is invisible to a direct edge-rank reading and to a cluster-assignment reading; it lives in the diffusion geometry that reads how information flows.**" ✓ Matches `feedback_lrg_step_is_confirmation_not_resolution`.
- **§6.2 Ψ-irrelevance methodological cornerstone**: "**The standard log-gap detector Ψ(n; τ) was designed to find topological scale boundaries on networks where the propagator vanishes across edge-absence boundaries; on the dense weight-heterogeneous |ImCoh| substrate no such gaps exist, the merge-height profile decays linearly rather than log-linearly, and Ψ pins to dendrogram boundaries as a generic property of UPGMA on this substrate (89/180 cells at strict argmax-at-an-extremum). The reframing that follows — reading D(τ) directly rather than through a privileged partition cut, and using the dendrogram as a hierarchical encoding rather than a partition selector — is not a workaround but the principled adaptation of LRG to dense weighted graphs, and it is the layer at which the band-resolved memory imprint surfaces.**" ✓ Matches `_psi_argmax_pinning.md`.

## Numerical corrections (action: writing agent)

- §6.1 β: "**within-probe BH-q = 0.006 on both KC axes**" — confirm m=12 within-probe BH-q. Per `_kc_beta_pvalue_check.md` and `_kc_crossprobe_check.md`, q = 0.006 on both. ✓
- §6.1 β: "**10/10 patients below their split-baseline null**" — verify against `_kc_beta_pvalue_check.md`. The §5.2 verifier shows the paired Wilcoxon vs within-baseline null gives p = 0.001 on both axes; need to confirm the 10/10 patient count specifically (the prose is cited in the §5.2 manuscript text and reproduces in `result_2_lrg_beta_trace.md` per memory). If only 9/10 fall below null with one tied, the count should be 9/10. **Verify against** `data/audit/section5_v2_kc_controls/per_patient_table.csv`.
- §6.1 low-γ: "**KC λ=1, 8/10**" with within-probe q = 0.055 — matches §5.2 cross-probe restriction reading. ✓
- §6.1 low-γ: "**(VI(k)) k ∈ [45, 56]**" — matches `_section5_4_verify.md` (12-cell contiguous span). ✓
- §6.1 closing single-region claim at γ_l: "**p_hyper = 5.2×10⁻⁶**" + "**only Bonferroni-survived (band, region) cell across the m=48 cohort-eligible cells**" — matches `_section5_5_verify.md`. ✓

## Framing rewrites (action: writing agent)

### Issue 1 — "almost universal in β" phrasing — DON'T LET IT CREEP IN

- The user's brief note flagged that an earlier draft framing of "almost universal in β" is **too strong for 7-8/10 patients**. The current §6.1 prose uses **measured language**: "registers at every geometric level the analysis tests"; "10/10 patients below their split-baseline null" (the only "10/10" claim is at the within-baseline null layer, which is a paired test against the patient's own drift floor — not a population-universal claim).

  **Action**: keep §6.1 as written. If a future writing-pass tightens the language to "almost universal in β" or similar, push back: 7-8/10 across the cohort is the right framing for the count-level claim; "registers at every layer" is the right framing for the multifacet claim; "10/10 below their split-baseline null" is the right framing for the within-baseline-null layer. Do not collapse these three into a single "universal" framing.

- The §6.1 closing: "**the post-task rest is not statistically equivalent to the rsPre and the four-phase paradigm is non-ergodic at the rest level in those bands**" — measured. ✓ Keep.

### Issue 2 — deferred controls should be a dedicated paragraph

- §6.3 Outlook currently mentions in a single paragraph: "**The cohort claim of §5 is established under the within-baseline null and drift-floor controls applied to the per-pair correlation; two cohort-wide drift insulations remain owed before the persistence language can be elevated from conditional to established. The first is a cohort-wide drift-triangle null at the LRG layer, mirroring Control 1 of §4.6, currently run only at Pat_06 where R² ≤ 0.05 in nearly every (band, distance) cell. The second is a symmetric cross-baseline check at the diffusion layer, mirroring Control 3 of §4.6, comparing D(τ) computed on (taskTest-late, rsPost-early) against (rsPre-late, rsPost-early) on matched half-segment noise budgets.**"

  This is **load-bearing for the persistence claim** but currently buried in a single paragraph among the other Outlook items. Recommend **promoting it to a dedicated short subsection** §6.3 ("Deferred controls") between §6.1/§6.2 and the Outlook's per-patient anatomical / cross-probe-aggregate / τ-sweep items.

  **Suggested rewrite**:

  > "## 6.3.1 Deferred controls (load-bearing for elevation from conditional to established).
  >
  > Two LRG-layer controls are owed before the persistence language can be elevated from conditional to established. Both have direct counterparts in the substrate audit of §4.6 and have been argued for in the §5.3 framing of the controlled per-pair correlation as the load-bearing claim under within-baseline null + drift-floor.
  >
  > **(i) LRG-layer drift-triangle null** (cohort-wide). Currently run on Pat_06 only, where R² ≤ 0.05 in nearly every (band, distance) cell at the LRG ultrametric level. The cohort-wide version mirrors §4.6 Control 1: split each rest phase into four consecutive chunks, compute all within-rest D(τ) distances at fixed τ = 1/λ_max, fit the linear within-rest drift fit against temporal gap, and report cohort-median R² per (band, distance) cell.
  >
  > **(ii) Symmetric cross-baseline at D(τ)** (cohort-wide). Mirrors §4.6 Control 3: compute D(τ) on (taskTest-late, rsPost-early) versus (rsPre-late, rsPost-early), both quantities using the same half-segment noise budget so the comparison is noise-symmetric. The directional-pass count at the LRG layer is the test-of-task-specificity.
  >
  > Until these two controls complete, the §5.3 cohort claim stands as a **controlled per-pair signal under the documented null architecture** and **NOT as a fully drift-insulated persistence statement**; the manuscript will be updated accordingly when the controls land."

  This makes the gating role explicit (per `feedback_renormalization_style` honest-head-first principle).

### Issue 3 — Outlook anchor anatomy + multiscale taxonomy

- §6.3 Outlook closes with three open directions:
  1. Per-patient anatomical mapping of trace-class versus anchor-class versus reset-class spatial distributions at the per-patient level, **deferred from the structural-illustration role taken in §5.6**.
  2. Cross-probe-aggregate anti-alignment of Pat_07 and Pat_15 merits an implant-geometry analysis.
  3. Multiscale τ-sweep at coarser τ.

  **Action**: in (1), explicitly cross-reference the §5.6 caveat that **the anchor class has no anatomy-only baseline**: "**...per-patient anatomical mapping ... deferred from the structural-illustration role taken in §5.6, including the anchor-class hypergeometric enrichment against the cohort cortical-contact pool that has no analogue in the §5.5 trace-class anatomy**."

  In (3), connect explicitly to the §5.6 τ-fixed caveat: "**The diffusion resolution is fixed at τ' = 1/λ_max throughout §5 for the reasons developed in §5.1, but a controlled τ-sweep that reads D(τ) at coarser scales tests whether the band-resolved imprint sharpens or dissolves as the propagator integrates over longer pathways, and is the natural extension of the per-pair reading to a fully multiscale diffusion-geometry claim. The same τ-sweep applied to the §5.6 taxonomy would test whether the trace-class leaves at τ_max are also trace-class at τ ≫ τ_max, addressing the multiscale class-assignment open question of §5.6.**"

  These are minor textual extensions but they tighten the §5 → §6.3 outlook hand-off.

## Caveats to add (action: writing agent)

- Promote the deferred-controls paragraph to a dedicated §6.3.1 subsection (see Issue 2). This is the most important §6 action item.
- Add the §5.6 caveat cross-references in §6.3 Outlook (see Issue 3).

## Figure actions (coding agent)

- None at §6 level. §6 is text-only.

## Deferred / questions

- The two cohort-wide drift insulation controls (LRG drift-triangle null at all 10 patients × 6 bands × 3 distances; symmetric cross-baseline on D(τ)) should be queued. Both are well-defined extensions of §4.6 Controls 1 and 3 lifted to the LRG layer. Estimate: same compute scale as the §4.6 substrate audit (already cached), so should be on the order of one workpackage script + one figure.
- Should the `\textsc{Notes\textsubscript{MSC}}` follow-up note these two controls explicitly when it is updated? Probably yes — they form the "bridge to NOTES_MSC" piece of the §5 → §6 hand-off.
