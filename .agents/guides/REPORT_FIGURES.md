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

## Section: Cross-Condition MSLCD Diagnostics (`sec:cross_condition`, Subsection 5.2)

This subsection presents a systematic comparison of the LRG multiscale community detection diagnostics across all patients, experimental phases, and frequency bands. It includes cross-condition heatmaps and boxplots (Part B), edge-weight threshold analysis (Part C), and a summary LaTeX table (Part D).

### Quick regenerate

```bash
# Step 1: Compute master diagnostics CSV + coarsening trajectories (~5 min)
lrg-eegfc compute diagnostics -v

# Step 2: Compute threshold analysis for Pat_02/rsPre (~2 min)
lrg-eegfc compute threshold-analysis -v

# Step 3: Generate all 9 figures + LaTeX table (~15 sec)
lrg-eegfc plot cross --figure-type all -v
```

### Figure Index

| Label | File | CLI command | Data Source | Speed |
|-------|------|-------------|-------------|-------|
| `fig:coarsening_beta` | `fig_B1_coarsening_beta.pdf` | `plot cross --figure-type coarsening-beta` | `mslcd_coarsening_trajectories.npz` | Fast |
| `fig:coarsening_patient` | `fig_B2_coarsening_Pat_02.pdf` | `plot cross --figure-type coarsening-patient` | `mslcd_coarsening_trajectories.npz` | Fast |
| `fig:metastability_heatmaps` | `fig_B3_metastability_heatmaps.pdf` | `plot cross --figure-type metastability-heatmaps` | `mslcd_diagnostics_master.csv` | Fast |
| `fig:metastability_boxplots` | `fig_B4_metastability_boxplots.pdf` | `plot cross --figure-type metastability-boxplots` | `mslcd_diagnostics_master.csv` | Fast |
| `fig:partition_richness` | `fig_B5_partition_richness.pdf` | `plot cross --figure-type partition-richness` | `mslcd_diagnostics_master.csv` | Fast |
| `fig:community_balance` | `fig_B6_community_balance.pdf` | `plot cross --figure-type community-balance` | `mslcd_diagnostics_master.csv` | Fast |
| `fig:threshold_susceptibility` | `fig_C1_threshold_susceptibility.pdf` | `plot cross --figure-type threshold-susceptibility` | `threshold_analysis_Pat02.npz` | Fast |
| `fig:threshold_connectivity` | `fig_C2_threshold_connectivity.pdf` | `plot cross --figure-type threshold-connectivity` | `threshold_analysis_Pat02.npz` | Fast |
| `tab:mslcd_diagnostics` | `mslcd_diagnostics_Pat02.tex` | `plot cross --figure-type latex-table` | `mslcd_diagnostics_master.csv` | Fast |

### Output Directories

- **B-figures:** `data/figures/report_mslcd_section/cross_condition/`
- **C-figures:** `data/figures/report_mslcd_section/threshold_analysis/`
- **Tables:** `data/tables/`

### Data Files Produced

| File | Content | Size |
|------|---------|------|
| `data/tables/mslcd_diagnostics_master.csv` | 114 rows x 30 columns — all scalar diagnostics for every (patient, phase, band) triplet | 27 KB |
| `data/tables/mslcd_coarsening_trajectories.npz` | Per-triplet arrays: `{patient}__{phase}__{band}__tau`, `__nmax`, `__eigenvalues` | 229 KB |
| `data/tables/threshold_analysis_Pat02.npz` | 6 bands x 9 thresholds: eigenvalues, S(tau), C(tau), connectivity metrics (all LCC-based) | 367 KB |
| `data/tables/mslcd_diagnostics_Pat02.tex` | Booktabs LaTeX table with 6 diagnostic columns | 2 KB |

### Parameters Used

- **Patients:** Pat_02, Pat_03, Pat_05, Pat_08 (Pat_06, Pat_07 excluded — missing phases)
- **Phases:** rsPre, taskLearn, taskTest, rsPost
- **Bands:** delta, theta, alpha, beta, low_gamma, high_gamma
- **FC method:** MSC (dense, sparsify=none, nperseg=4096)
- **Metastability:** n_tau=20 log-spaced tau in [1/lambda_max, 1/lambda_2], monotonic non-increasing n_max(tau) constraint, normalization by coarsening events only
- **Sensible peaks:** mean-floor + log-drop chain (LOG_DROP_MAX=0.5 decades, MIN_N_SPACING=5)
- **Threshold analysis:** Patient=Pat_02, Phase=rsPre, percentiles=[0, 50, 70, 80, 85, 90, 95, 97, 99]
- **Community balance:** H_size = normalized size entropy, H / ln(n*), range [0,1]

### Figure Descriptions

#### Part B: Cross-Condition Figures

**B1 — Coarsening trajectories (beta band).** Step-plots of n_max(tau) vs log10(tau) for the beta band across all patients. Each subplot is one patient, phase colours distinguish rsPre (blue), taskLearn (orange), taskTest (green), rsPost (red). Common legend below grid. No titles. Shows how hierarchical resolution degrades with diffusion time — steeper curves indicate stronger multiscale structure.

**B2 — Coarsening trajectories (Pat_02, all bands).** Same layout as B1 but fixed patient = Pat_02, one subplot per frequency band. Band names annotated in LaTeX (e.g., $\delta$, $\theta$). Reveals band-dependent coarsening speed: high_gamma collapses quickly (few effective channels), while delta/theta sustain high n_max over a wider tau range.

**B3 — Metastability heatmaps.** Two rows: mean metastability mu_bar (top) and max metastability mu_max (bottom). One column per patient. Each cell is a band x phase value. Square cells with magma colormap and luminance-adaptive text annotations (white on dark, black on light). Shared colorbar per row. Highlights which band/phase combinations exhibit the most community reassignment under the LRG coarsening flow.

**B4 — Metastability boxplots.** One subplot per band, boxplots of mu_bar grouped by phase. Individual patient values overlaid as shaped markers (circle=Pat_02, square=Pat_03, diamond=Pat_05, down-triangle=Pat_08). Combined legend outside the grid shows both phase colours and patient marker shapes. Reveals inter-patient variability within each band/phase cell.

**B5 — Partition richness heatmap.** N_sens (number of sensible PSI peaks) as band x phase heatmap, one panel per patient. YlGnBu colormap, integer annotations, square cells. High N_sens indicates rich multiscale structure; low N_sens (e.g., alpha taskLearn/taskTest for Pat_02 = 1) indicates a single dominant partition scale.

**B6 — Community balance heatmap.** H_size (normalized community size entropy) as band x phase heatmap, one panel per patient. RdYlGn colormap (red=unbalanced, green=balanced), range [0,1], square cells. H_size near 1 means communities are roughly equal-sized; near 0 means one giant community dominates.

#### Part C: Threshold Analysis (Pat_02, rsPre)

**C1 — Entropic susceptibility at thresholds.** 2x3 grid (one subplot per band). Curves show C(tau) = -dS/d(log10 tau) at 9 edge-weight threshold percentiles (0% to 99%), coloured by plasma colormap. Log-scale x-axis (tau), extending to tau = 10^3. Shared colorbar indicates threshold percentile. Reveals how progressive edge removal exposes hidden multi-peak structure in C(tau): at low thresholds the susceptibility has one dominant peak, but at high thresholds (90-99%) secondary peaks emerge or the dominant peak splits as weak edges are pruned.

**C2 — Connectivity metrics vs threshold.** 2x2 grid: (a) connected components, (b) LCC fraction |C_1|/N, (c) LCC density (edges_in_lcc / max_possible), (d) transitivity (global clustering coefficient). X-axis: "Edges remaining (%)" on inverted log scale (100% at left, 1% at right). One curve per band. Band-dependent disconnection thresholds visible: delta network remains connected until ~85% removal, while high_gamma fragments already at ~50%. LCC density actually increases at extreme thresholds (surviving core is dense clique-like). Transitivity reveals band-specific clustering topology under progressive pruning.

#### Part D: LaTeX Table

**Tab — MSLCD diagnostics for Pat_02.** Booktabs table with multirow band grouping. 6 diagnostic columns: N_sens (partition richness), n_max(tau') (communities at finest scale), n_max(1/lambda_2) (communities at coarsest scale), mu_bar (mean metastability), mu_max (max metastability), H_size (community balance). Ready for direct inclusion in Overleaf.

### Code Architecture

| Module | Function | Purpose |
|--------|----------|---------|
| `workflow.diagnostics` | `compute_triplet_diagnostics()` | Full LRG pipeline for one (patient, phase, band) triplet — spectral properties, entropic susceptibility, partition richness, coarsening + metastability, community balance |
| `workflow.diagnostics` | `compute_threshold_analysis()` | Progressive edge-weight pruning with LCC extraction, eigenvalue decomposition, entropy/susceptibility curves, connectivity metrics |
| `workflow.diagnostics` | `threshold_matrix()` | Zero out edges below given weight percentile |
| `workflow.diagnostics` | `sensible_psi_peaks()` | Mean-floor + log-drop chain peak detection for PSI |
| `visuals.cross_condition` | `plot_coarsening_beta()` | B1 figure |
| `visuals.cross_condition` | `plot_coarsening_patient()` | B2 figure |
| `visuals.cross_condition` | `plot_metastability_heatmaps()` | B3 figure |
| `visuals.cross_condition` | `plot_metastability_boxplots()` | B4 figure |
| `visuals.cross_condition` | `plot_partition_richness()` | B5 figure |
| `visuals.cross_condition` | `plot_community_balance()` | B6 figure |
| `visuals.cross_condition` | `plot_threshold_susceptibility()` | C1 figure |
| `visuals.cross_condition` | `plot_threshold_connectivity()` | C2 figure |
| `visuals.cross_condition` | `generate_latex_table()` | LaTeX table |
| `visuals.cross_condition` | `generate_cross_condition_figure()` | Master dispatcher |
| `cli.compute` | `diagnostics` command | CLI entry: compute master CSV + trajectories NPZ |
| `cli.compute` | `threshold-analysis` command | CLI entry: compute threshold analysis NPZ |
| `cli.plot` | `cross` command | CLI entry: generate all figures from cached data |

### Implementation Notes

- **LCC extraction is critical for threshold analysis.** After removing low-weight edges the graph can disconnect. Computing the Laplacian on the full disconnected graph injects spurious zero eigenvalues (appearing as ~10^-16) that corrupt the entropy/susceptibility curves and spectral analysis. All eigenvalues and entropy curves are computed on the Largest Connected Component (LCC) subgraph.
- **Connectivity metrics are non-trivial.** Edge fraction and mean degree were replaced with LCC density and transitivity because the former are trivially determined by the threshold percentile (edge_frac = 1 - pct/100 by definition). LCC density and transitivity reveal actual topological changes: density increases at extreme thresholds as the surviving core becomes clique-like, and transitivity shows band-specific clustering decay.
- **Data-driven metastability.** The coarsening trajectory uses a monotonic non-increasing constraint on n_max(tau) to prevent flickering from marginal PSI peaks. Metastability mu_i is normalized by the number of coarsening events (steps where n_max actually changes), not total steps, so increasing n_tau doesn't dilute the signal.
- **C2 eigenvalue spectrum figure was dropped** as uninformative for the cross-condition narrative. The three-panel layout (sorted eigenvalues, lambda_2 vs threshold, spectral gap ratio) was messy and didn't add insight beyond what C1 and C2-connectivity already show.
- **Square heatmap cells** via `aspect="equal"` in imshow, with luminance-adaptive text colour (white if cell luminance < 0.5, black otherwise) for readability on both dark and light cells.

### Key Findings (from generated data)

- **Band-dependent coarsening:** High_gamma networks collapse to 1-2 communities quickly (n_max(tau') = 6-19), while delta/theta networks sustain 40-57 communities over a wide tau range. This reflects the sparser effective connectivity in high-frequency bands.
- **Phase sensitivity of metastability:** Task phases (taskLearn, taskTest) show markedly different metastability from resting phases in delta/theta (e.g., delta rsPre mu_bar=0.436 vs taskLearn mu_bar=0.051). Alpha taskLearn/taskTest collapse to a single partition (N_sens=1, n_start=6).
- **Threshold-revealed multi-scale structure:** Progressive edge removal in C1 exposes secondary C(tau) peaks that are hidden in the full-weight graph, particularly in delta and theta bands where a single dominant peak at low thresholds splits into multiple peaks at 90-95% removal.
- **Band-dependent disconnection:** Delta/beta networks remain connected until 85-90% edge removal; high_gamma fragments already at ~50% removal, consistent with sparser high-frequency FC.
- **LCC densification:** At extreme thresholds (>95%), LCC density increases despite fewer total edges — the surviving core forms a dense subgraph, suggesting a hub-like backbone structure in the FC network.

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
