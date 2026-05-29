---
name: methods-directive-joint-density-symbol
era: IMCOH_ABS_COHORT_N10
status: pending
kind: methods-agent-directive
date: 2026-05-28
target: §results / Fig 1 caption — colorbar label of the joint-rank density panel
priority: low (notation lock; no number changes)
source_of_truth:
  - scripts/02_preprint/preprint_17_band_joint_density_single.py
  - scripts/02_preprint/preprint_09_bands_joint_density.py
  - scripts/02_preprint/preprint_11_bands_joint_density_3d_grid.py
  - scripts/02_preprint/preprint_12_bands_joint_density_3d_stacked.py
applied_in: []
---

# Methods directive — joint-density colorbar symbol `Δc⋆`

## Head

The Fig 1 (and §results 3D-grid / 3D-stacked figures) **joint-rank density
colorbar** label is locked at the compact symbol

    Δc⋆(u, v; b)

read aloud as "delta-c-star", where

    Δc⋆(u, v; b) = sgn(c_obs(u, v; b) - 1) · [ |c_obs(u, v; b) - 1|
                                              - |c_surr^p95(u, v; b) - 1| ]_+

and:

* `c_obs(u, v; b)` is the joint density on the rank grid for band `b`,
  taken either as the analytic Gaussian copula at the cohort
  `ρ_obs(b)` (the model) or the per-patient KDE on rank pairs
  averaged across the cohort (the data — selectable via
  `preprint_17_band_joint_density_single.py --source {copula,empirical}`);
* `c_surr^p95(u, v; b)` is the matched-strength noise-floor density
  evaluated at the surrogate Spearman p95 (always analytic Gaussian
  copula at `ρ_surr^p95(b)`, regardless of `c_obs` source);
* `[·]_+ = max(·, 0)` clips below-noise-floor deviations to zero;
* `sgn(·)` carries the trace / anti direction (positive → diagonal
  enrichment / trace; negative → anti-diagonal enrichment / anti).

The previous verbose label
**"signed enrichment over matched-strength noise floor"**
is retired in figure axes — the symbol carries the meaning, the full
expansion lives only in the manuscript caption and this directive.

## Why

* The user finds long descriptive labels in figures unacceptable (see
  [[feedback-no-text-in-figures]]).  Compact symbol + caption is the
  norm in this preprint.
* The quantity was previously named `signed_excess` only as a python
  variable (`scripts/02_preprint/preprint_09_*`, `_11_*`, `_12_*`,
  `_17_*`); no symbolic notation existed.  Without a locked symbol
  every figure ended up with a different long label.
* `Δc⋆` parallels the project's existing `Δ_task` / `Δ_rest` /
  `Δ_D_coph` notation (capital delta for an observed minus reference
  quantity); the star marks the matched-strength clipping.  Does NOT
  collide with `ρ_split^coph` (Spearman correlation, §5.3) or
  `T_d` (triangle scalar, [[feedback-td-sign-convention]]).

## How to apply

Required wherever the joint-rank density figure (Fig 1 panel B candidate,
or the §results 6-band atlas) is referenced or captioned.

### In Fig 1 caption (LaTeX)

    \emph{Δc⋆} is the matched-strength noise-floor-corrected joint-rank
    density excess on the unit-square rank grid, defined as
    \[
        \Delta c^{\!\star}(u, v; b) =
            \mathrm{sgn}\bigl(c_{\rm obs}(u, v; b) - 1\bigr) \cdot
            \bigl[ |c_{\rm obs}(u, v; b) - 1|
                  - |c_{\rm surr}^{p95}(u, v; b) - 1| \bigr]_+ ,
    \]
    where $c_{\rm obs}$ is the cohort joint density (analytic Gaussian
    copula at $\rho_{\rm obs}(b)$ in panel/source = ``copula''; per-
    patient KDE cohort mean in panel/source = ``empirical''),
    $c_{\rm surr}^{p95}$ is the analytic Gaussian copula at the
    matched-strength surrogate Spearman p95, and $[\cdot]_+ =
    \max(\cdot, 0)$.  Positive $\Delta c^{\!\star}$ marks
    diagonal-trace enrichment; negative marks anti-diagonal
    enrichment.

### In the 6-band atlas (preprint_09)

Replace the `cb.set_label(r"signed enrichment over matched-strength
noise floor", ...)` line with `cb.set_label(r"$\Delta c^{\!\star}$",
...)`.  Same applies to `preprint_11` and `preprint_12`.

### In the companion .md (when Fig 1 lands)

Cite Δc⋆ in the figure's `.md` notation block.  No need to redefine —
point at this directive.

## Open

* If a reviewer asks why `c_surr^p95` is always the analytic copula
  even when `c_obs` is empirical KDE: because we lack per-pair surrogate
  rank pairs cached, only the cohort-level `surr_median_rho_p95`
  scalar.  The matched-strength surrogate at this Spearman level IS
  Gaussian-copula-shaped on the rank grid by construction (Spearman
  of paired matched-strength surrogates is the only DoF the surrogate
  leaves on the unit square).  This is consistent across sources.
