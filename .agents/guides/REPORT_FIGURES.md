# Report Figure Registry

Master index of all figures in the technical report, how to generate them, and where they live.

## Output Locations

- **Fast figures (cache-only):** `data/figures/report_fc_section/Pat_02/`
- **Precomputed data:** `data/fc_fig_cache/`
- **Generation scripts:** `scripts/gen_fc_figures_fast.py`, `scripts/precompute_fc_figures.py`
- **Notebook:** `ipynb/06_presentation_figures/07_fc_section_figures.ipynb`

---

## Section: Signal Processing and FC Estimation (`sec:fc`)

### Quick regenerate

```bash
# Fast figures (2,3,4,6) — ~10 seconds
python scripts/gen_fc_figures_fast.py

# Precompute slow data for figs 1,5 — ~3 minutes
python scripts/precompute_fc_figures.py

# Then run notebook for figs 1,5 (or full notebook)
jupyter nbconvert --execute ipynb/06_presentation_figures/07_fc_section_figures.ipynb
```

### Figure Index

| Label | File | Script/Cell | Data Source | Speed |
|-------|------|-------------|-------------|-------|
| `fig:msc_example` | `fig_msc_example.pdf` | notebook cell `fig1-plot` | `fc_fig_cache/*_coh_*.npz` | Needs precompute |
| `fig:adjacency_heatmaps` | `fig_adjacency_heatmaps.pdf` | `gen_fc_figures_fast.py` | `msc_cache/` (dense) | Fast |
| `fig:weight_distributions` | `fig_weight_distributions.pdf` | `gen_fc_figures_fast.py` | `msc_cache/` (dense) | Fast |
| `fig:network_properties` | `fig_network_properties.pdf` | `gen_fc_figures_fast.py` | `msc_cache/` (dense) | Fast |
| `fig:surrogates` | `fig_surrogates.pdf` | notebook cell `fig5-plot` | `fc_fig_cache/*_null_*.npz` | Needs precompute |
| `fig:surrogate_comparison` | `fig_surrogate_comparison.pdf` | `gen_fc_figures_fast.py` | `msc_cache/` (dense + soft) | Fast |

### Parameters Used

- **Patient:** Pat_02 (representative)
- **Phase:** rsPre (default), all 4 phases for grid/metric plots
- **Bands:** All 6 canonical bands
- **nperseg:** 1024
- **Surrogates:** 100 (comparison), 200 (null distribution)
- **Colormap:** magma (matrices), per-band colors for plots
- **vmax:** 0.6 (heatmaps shared scale)

### Metrics in Fig 4 (network_properties)

- (a) Node strength distribution (boxplot per band, rsPre)
- (b) Weighted clustering coefficient (all phases x bands)
- (c) Strength assortativity (Pearson corr of endpoint strengths)
- (d) Weight concentration (fraction of total weight in top 5% edges)

---

## Section: Laplacian Diffusion States and Communication Ultrametrics (`sec:ultrametrics`)

### Quick regenerate

```bash
# All 4 figures — ~3 seconds
python scripts/gen_ultrametric_figures.py
```

### Figure Index

| Label | File | Script | Data Source | Speed |
|-------|------|--------|-------------|-------|
| `fig:lrg_overview` | `fig_lrg_overview.pdf` | `gen_ultrametric_figures.py` | None (schematic) | Fast |
| `fig:rho_tau` | `fig_rho_tau.pdf` | `gen_ultrametric_figures.py` | `msc_cache/` (beta, rsPre) | Fast |
| `fig:D_tau` | `fig_D_tau.pdf` | `gen_ultrametric_figures.py` | `msc_cache/` (beta, rsPre) | Fast |
| `fig:spectrum_tau_window` | `fig_spectrum_tau_window.pdf` | `gen_ultrametric_figures.py` | `msc_cache/` (beta, rsPre) | Fast |

### Output Directory

`data/figures/report_ultrametric_section/Pat_02/`

### Parameters Used

- **Patient:** Pat_02, **Phase:** rsPre, **Band:** beta
- **Laplacian:** Combinatorial (L = D - A), NOT normalized
- **τ values:** τ' = 1/λ_max = 0.056, τ_mid = 0.392 (geometric mean), 1/λ_2 = 2.73
- **Spectrum:** λ_2 = 0.3658, λ_max = 17.82 (N = 117 channels)
- **Colormaps:** inferno (ρ), cividis (D, log10 scale)
- **Leaf ordering:** average-linkage dendrogram from D(τ_mid)

---

## Section: Multiscale Community Detection (`sec:mslcd`)

### Quick regenerate

```bash
# All 8 figures (6 PDF + 2 interactive HTML) — ~10 seconds
python scripts/gen_mslcd_figures.py
```

### Figure Index

| Label | File | Script | Data Source | Speed |
|-------|------|--------|-------------|-------|
| `fig:dendrogram_n2` | `fig_dendrogram_n2.pdf` | `gen_mslcd_figures.py` | `lrg_cache/` (beta, rsPre) | Fast |
| `fig:dendrogram_n8` | `fig_dendrogram_n8.pdf` | `gen_mslcd_figures.py` | `lrg_cache/` (beta, rsPre) | Fast |
| `fig:dendrogram_n16` | `fig_dendrogram_n16.pdf` | `gen_mslcd_figures.py` | `lrg_cache/` (beta, rsPre) | Fast |
| `fig:dendrogram_n31` | `fig_dendrogram_n31.pdf` | `gen_mslcd_figures.py` | `lrg_cache/` (beta, rsPre) | Fast |
| `fig:psi_profile` | `fig_psi_profile.pdf` | `gen_mslcd_figures.py` | `msc_cache/` (beta, rsPre) | Fast |
| `fig:entropy_susceptibility` | `fig_entropy_susceptibility.pdf` | `gen_mslcd_figures.py` | `msc_cache/` (beta, rsPre) | Fast |
| `fig:nstar_tau` | `fig_nstar_tau.pdf` | `gen_mslcd_figures.py` | `msc_cache/` (beta, rsPre) | Fast |
| `fig:metastable_nodes` | `fig_metastable_nodes.pdf` | `gen_mslcd_figures.py` | `msc_cache/` + electrode coords (nilearn) | Fast |
| `fig:community_alluvial` | `fig_community_alluvial.pdf` | `gen_mslcd_figures.py` | `msc_cache/` (beta, rsPre) | Fast |
| _(interactive)_ | `fig_sankey_community_beta.html` | `gen_mslcd_figures.py` | `msc_cache/` via lrgsglib (Plotly Sankey) | Fast |
| _(interactive)_ | `fig_brain_connectome_beta.html` | `gen_mslcd_figures.py` | `msc_cache/` + lrg_cache + fsaverage mesh (Plotly 3D) | Fast |

### Output Directory

`data/figures/report_mslcd_section/Pat_02/`

### Parameters Used

- **Patient:** Pat_02, **Phase:** rsPre, **Band:** beta
- **Laplacian:** Combinatorial (L = D - A)
- **Spectrum:** λ_2 = 0.3658, λ_max = 17.82 (N = 117)
- **τ window:** [0.0561, 2.7334]
- **Characteristic scales:** 1 peak at τ* = 1.0596
- **n*(τ):** 2 (global PSI max at all scales — simple 2-community structure)
- **Prominent PSI peaks (log-scale half-rule):** n* = 2, 8, 16, 31 (4 peaks, LOG_DROP_MAX=0.5 decades, MIN_N_SPACING=5)
- **Alluvial/Sankey pairs:** [(0.1,15), (0.3,12), (0.5,8), (0.8,6), (1.0,5), (2.0,3), (5.0,2)] — shared by alluvial PDF + Sankey HTML
- **Metastable nodes:** 101/117 (μ > 0), max μ = 0.500, mean μ = 0.178
- **Dendrogram:** Cached LRG linkage (lrgsglib default τ), optimal_threshold = 0.3644, split into M separate PDFs (one per prominent PSI peak)

### Interactive Figure Parameters

- **Sankey (fig 7):** tau/cluster pairs = [(0.1,15), (0.3,12), (0.5,8), (0.8,6), (1.0,5), (2.0,3), (5.0,2)], via `compute_clustering_across_tau()` + `create_sankey_diagram()` from `visuals/metastable.py`
- **Brain connectome (fig 8):** Plotly 3D + fsaverage pial mesh (opacity=0.08), ALL 6786 edges with alpha ramp (0.02→0.9, power=1.5), power-law width (w^2.5, 0.3→20px, 20 buckets), Greys colormap, LRG cluster colors (fcluster maxclust=31 from finest PSI peak), tab20 colormap, node size=8, colorbar via dummy scatter, MNI-transformed coords

### Key Finding

Pat_02 beta band has a simple dominant 2-community split (n* = 2 at all τ). Log-scale PSI analysis reveals 4 prominent hierarchical levels: n* = 2, 8, 16, 31. Metastability analysis (101/117 nodes with μ > 0) shows fine-grained community dynamics across 7 shared τ/n scales.

---

## Section: [NEXT SECTION PLACEHOLDER]

_To be filled as figures are generated._

---

## Notes

- **nperseg=4096** is the new default (was 1024). See `memory/nperseg_migration.md` for full rationale. df=0.5 Hz gives 7 bins in delta (was 2). Beta-band MSC matrices differ by < 1.4% (r=0.9999).
- Dense MSC matrices are fully connected → path-based graph metrics (global efficiency, degree assortativity) are trivial. Use weighted/strength-based alternatives.
- Soft sparsification has negligible effect on Pat_02 data (r=1.000000, max diff ~2.5e-3). This supports the report's conclusion that MSC estimates are stable for these long recordings.
- **50 Hz power line harmonics** visible in high_gamma at nperseg=4096: 300 Hz (median MSC=0.556, 10.8x neighbors), 100 Hz (4.5x). Band-average inflation < 1%. No pipeline change needed — note in caption.
- Channel pair changed from (21, 24) to (9, 34) after nperseg migration (95th percentile of beta coherence shifted slightly).
