---
name: writing-directive-raw-copula-figure-swap
era: IMCOH_ABS_COHORT_N10
status: current
kind: directive
scope: Swap the raw |ImCoh| per-pair joint-density figure in Results §ssec:raw_fc from the empirical-KDE render to the Gaussian-copula render, with the body sentence and caption edited to match and a short formal justification of the analytic representation. Single subsection, single figure. Changes no locked verdict; the trace reading (non-specific across all bands, then exposed by the nulls) is unchanged.
created: 2026-06-03
supersedes: the raw joint-density figure recipe + caption in directives/writing_directive_2026-06-03_raw-substrate-results-figures.md §2–§3 (that file assumed the empirical KDE). All other content of that directive (null triangle, raw→LRG magnitude-collapse) stands.
verified_against: data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv (per-band obs_median_rho, surr_median_rho_p95); figure data/preprint/figures/all_bands/fig_bands_joint_density_copula_rawfc.pdf rendered + inspected this session.
---

# Raw subsection: swap empirical-KDE joint density → Gaussian-copula joint density

**Head.** The raw `|ImCoh|` per-pair joint-density figure moves from the
empirical-KDE render to the **Gaussian-copula** render
(`fig_bands_joint_density_copula_rawfc.pdf`). The scientific reading is
identical — a positive diagonal in **every** band, so the substrate reports a
trace everywhere, which the nulls then expose. Only the *representation*
changes: each panel now shows the one-parameter copula implied by the cohort
rank-correlation instead of a pooled-rank histogram. This directive gives the
three exact edits and the short justification to insert.

> **Do not** mention the cophenetic-substrate copula here — that subsection is
> not written yet. Keep this confined to §`ssec:raw_fc`. Leave the third
> paragraph ("Read across the six bands…") and the null-triangle figure
> untouched; both remain correct verbatim.

---

## 1. Why the analytic copula (the justification to convey, briefly)

The per-pair trace test is **defined on the rank correlation** `\rhorawsplit`.
Rendering both the observed dependence and its null as the one-parameter
Gaussian copula fixed to that correlation makes each panel a direct,
like-for-like map of *observed-minus-null in the very statistic that defines
the trace* — rather than the difference of two nonparametric density estimates,
which would also differ in the incidental sampling structure of a finite,
discrete rank cloud. Keep this to **two sentences** (drafted in §2). Do **not**
expand into copula theory, do **not** claim "maximum entropy", do **not**
apologise for it being a model — state it as the natural coordinate for the
measure.

---

## 2. The three edits (find-and-replace against the current LaTeX)

### Edit A — figure include (one line)

**Find:**
```
\includegraphics[width=\linewidth]{figures/fig_bands_joint_density_empirical_rawfc.pdf}
```
**Replace:**
```
\includegraphics[width=\linewidth]{figures/fig_bands_joint_density_copula_rawfc.pdf}
```

### Edit B — second paragraph, the rendering sentence

**Find** (the clause beginning "rank them, and pool…"):
```
--- rank them, and pool the rank pairs over the cohort; their joint density is what each panel displays.
```
**Replace with** (analytic-render description + the two-sentence justification):
```
--- and rank them; their cohort dependence is summarized by the median per-patient rank correlation \(\rhorawsplit\). Each panel renders this dependence analytically, as the bivariate Gaussian copula whose single parameter equals \(\rhorawsplit\), set against the copula of the level reachable under matched-strength surrogates. We display the implied copula rather than a pooled-rank histogram because the trace test is itself defined on the rank correlation: rendering observation and null in this one-parameter form makes each panel a like-for-like map of observed minus null in the statistic that defines the trace, without the sampling structure a finite rank cloud would impose.
```
The sentence that follows it in the source — "A trace, that is a positive
\(\rhorawsplit\), corresponds to probability mass gathering along the main
diagonal: …" — stays **verbatim**: it is exactly correct for a positive-ρ
copula (diagonal concentration) and its ρ=0 isotropic limit.

### Edit C — the caption (full replacement)

**Replace the entire `\caption{…}` body with:**
```
\textbf{Per-pair trace at the raw \(\abs{\ImCoh}\) substrate, resolved by frequency band.} Each panel renders the cohort per-pair dependence as a bivariate Gaussian copula on the within-patient rank coordinates \(u=\operatorname{rank}\Delta^{\mathrm{raw}}_{\mathrm{task}}(i,j)\) (horizontal) and \(v=\operatorname{rank}\Delta^{\mathrm{raw}}_{\mathrm{rest}}(i,j)\) (vertical), its single parameter fixed to the cohort-median \(\rhorawsplit\) (top row \(\delta,\theta,\alpha\); bottom row \(\beta,\gamma_{\mathrm{l}},\gamma_{\mathrm{h}}\); \(n=10\)). The colour field is the signed enrichment of this copula over the copula fixed to the 95th-percentile \(\rhorawsplit\) reachable under matched-strength surrogates: green mass on the main diagonal (\(u\simeq v\)) marks the trace direction (\(\rhorawsplit>0\), a task-induced displacement preserved into post-task rest), red mass on the anti-diagonal marks reversal, and cream marks an observed correlation at or below the surrogate level. The white solid and dashed guides mark the main diagonal (\(u=v\), the trace axis) and the anti-diagonal (\(u=1-v\), the reversal axis); the marginal histograms give the empirical conditional split of the pooled ranks --- \(u\) for pairs in the upper versus lower half of \(v\) (top) and symmetrically for \(v\) (right).
```

Notes for the caption:
- The heatmap (signed excess) is the copula (analytic); the white solid/dashed
  lines are fixed reference axes (main = trace, anti = reversal), **not**
  contours; the **marginal histograms are empirical** (the pooled within-patient
  ranks). Keep that distinction as written — it is honest and it shows the data
  directly. (HDR density contours were removed: on an analytic copula they are a
  deterministic function of ρ and trace z_main, not the plotted signed excess.)
- The floor is the **95th-percentile** matched-strength level (more
  conservative than the mean-surrogate "noise floor" wording of the old
  empirical caption — update accordingly, do not carry the old phrase over).

---

## 3. Keep unchanged

- The first paragraph (substrate definition, `\rhorawsplit` as the hallmark).
- The third paragraph ("Read across the six bands… drift null… cross-probe…
  matched-strength surrogate sets the level…") — still correct verbatim; the
  copula figure still shows a diagonal in every band.
- The null-triangle figure `fig:rawfc_null_triangle` and its prose.
- All notation (`\rhorawsplit`, `\Adjacency`, `\Band`, `u`, `v`,
  `\Delta^{\mathrm{raw}}_{\cdot}`).

---

## 4. Verified figure facts (so the prose stays accurate)

Render = signed excess of `copula(cohort-median ρ_obs)` over
`copula(surr_median_rho_p95)`, per band, from
`raw_fc_matched_strength/cohort_summary_all_bands.csv`:

| band | obs_median_rho | surr p95 | obs − p95 | diagonal in figure |
|------|----------------|----------|-----------|--------------------|
| δ    | 0.111 | 0.032 | +0.079 | green |
| θ    | 0.117 | 0.013 | +0.104 | green |
| α    | 0.158 | 0.020 | +0.138 | green |
| β    | 0.258 | 0.037 | +0.221 | green (tightest) |
| γ_l  | 0.126 | 0.048 | +0.078 | green |
| γ_h  | 0.111 | 0.046 | +0.065 | green |

All six positive → **every band lit** (the non-specificity the paragraph
relies on). β has the largest correlation, hence the tightest diagonal; the
prose need not rank the bands, but if it does, β is the strongest, α second.
No band is anti or null at this substrate — that selectivity is the
cophenetic story, told later.
