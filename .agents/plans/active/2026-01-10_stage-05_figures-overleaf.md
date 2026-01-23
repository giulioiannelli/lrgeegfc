# Plan (2026-01-10): Stage 05 - Figures and Overleaf bundle

## Goals
- Produce the full figure set A-T with a consistent naming scheme.
- Split single-case (notebooks) vs batch (scripts) generation.
- Package figures/tables for Overleaf.

## Current status
- `data/figures/` exists with multiple subfolders from prior runs.
- Script wrappers exist in `src/visualize_*.py` but need alignment to the
  figure list below.
- No Overleaf bundle folder yet.
- Visualization functions currently live under `visuals/` (kept at top level).
- Time-window scripts added: `src/compute_time_windows.py`,
  `src/visualize_time_windows.py`.
- Time-window guide added: `.agents/guides/TIME_WINDOW_GUIDE.md`.
- Stage 05 notebooks created for legacy figure prototypes:
  - `ipynb/05_figures/01_timeseries_emd_singlepat.ipynb`
  - `ipynb/05_figures/02_corr_networks_per_band_singlepat.ipynb`
  - `ipynb/05_figures/03_lrg_full_panel_singlepat.ipynb`
  - `ipynb/05_figures/04_sankey_cluster_transitions.ipynb`
  - `ipynb/05_figures/05_band_pca_exploration.ipynb`
  - `ipynb/05_figures/06_poster_layouts.ipynb`
  - `ipynb/05_figures/07_single_patient_experiments.ipynb`

## Output locations
- Figures: `data/figures/<category>/Pat_XX/`
- Bundle: `outputs/overleaf/` (copy or link from cached outputs)

## Figure map

A) Time series + EMD + band overlays
- Type: single-case notebook
- Output: `data/figures/timeseries/Pat_XX/`
- Implementation target: `ipynb/05_figures/01_timeseries_emd_singlepat.ipynb`

B) MSC network + FC matrix (single case)
- Type: single-case notebook
- Output: `data/figures/msc/Pat_XX/`
- Implementation target: `ipynb/05_figures/02_msc_network_singlepat.ipynb`

C) All phases, all bands for one patient (grid)
- Type: batch script (per patient)
- Output: `data/figures/msc/Pat_XX/phase_band_grid.png`
- Implementation target: new script `src/visualize_msc_grid.py`
- Cache requirement: reads `data/msc_cache/Pat_XX/` only.

D) Same band/phase across all patients
- Type: batch script
- Output: `data/figures/msc/all_patients_<band>_<phase>.png`
- Implementation target: extend `src/visualize_msc.py` or new script
- Cache requirement: reads `data/msc_cache/Pat_XX/` only.

E) MSC vs correlation matrix comparison
- Type: notebook (single) + batch script (all)
- Output: `data/figures/comparison/`
- Implementation target: reuse `src/visualize_comparison.py`
- Cache requirement: reads `data/msc_cache/` and `data/corr_cache/`.

F) Network analysis comparison (MSC unvalidated vs validated vs correlation)
- Type: notebook (single) + batch script (all)
- Output: `data/figures/comparison/` and `results/`
- Implementation target: new script `src/compare_network_variants.py`
- Cache requirement: reads cached MSC/corr matrices only.

G) Marchenko-Pastur noise analysis + spectrum
- Type: notebook (single) + batch script (all)
- Output: `data/figures/cleaning/`
- Implementation target: extend `src/clean_correlation_matrices.py` with plots
- Cache requirement: reads `data/corr_cache/` or writes `data/corr_cache/` if missing.

H) Unvalidated MSC vs validated MSC (matrix + network)
- Type: notebook (single) + batch script (all)
- Output: `data/figures/msc_validation/`
- Implementation target: new script `src/visualize_msc_validation.py`
- Cache requirement: reads `data/msc_cache/` (sparsify none/soft).

I) LRG analysis pipeline summary
- Type: notebook (single) + batch script (all)
- Output: `data/figures/lrg/`
- Implementation target: extend `src/visualize_lrg.py`
- Cache requirement: reads `data/lrg_cache/`.

L) Density matrix video across tau
- Type: batch script (per patient/band/phase)
- Output: `data/figures/lrg_video/`
- Implementation target: new script `src/visualize_lrg_video.py`
- Cache requirement: reads `data/lrg_cache/` only.

M) Specific heat + entropy + PSI + PSI-cut network
- Type: notebook (single) + batch script (all)
- Output: `data/figures/lrg/`
- Implementation target: `visuals/lrg.py` + new notebook in `ipynb/03_lrg/`
- Cache requirement: reads `data/lrg_cache/` only.

N) Dendrograms across phases (per patient/band)
- Type: batch script
- Output: `data/figures/lrg/`
- Implementation target: extend `src/visualize_lrg.py` with phase grids
- Cache requirement: reads `data/lrg_cache/` only.

O) Reorganization metric distance per method
- Type: batch script
- Output: `data/figures/reorganization/`
- Implementation target: new script `src/visualize_reorganization_metrics.py`
- Cache requirement: reads `results/reorganization/`.

P) Cross-metric correlation matrix
- Type: batch script
- Output: `data/figures/reorganization/`
- Implementation target: new script `src/visualize_metric_correlation.py`
- Cache requirement: reads `results/reorganization/`.

Q) Metastable nodes (Sankey)
- Type: notebook (single) + batch script (all)
- Output: `data/figures/metastable/`
- Implementation target: `visuals/metastable.py` + `src/visualize_metastable.py`
- Cache requirement: reads `data/lrg_cache/` only.

R) Cross-subject statistical summary
- Type: batch script
- Output: `results/` and `data/figures/summary/`
- Implementation target: new script `src/aggregate_subject_stats.py`
- Cache requirement: reads `results/reorganization/` and writes summaries.

S) Time-window analysis + entropy animation
- Type: notebook (single) + batch script (all)
- Output: `data/figures/time_windows/` and video folder
- Note: define new, consistent window length/overlap defaults for reproducibility.
- Defaults:
  - Window length: max(10 s, 10 cycles at lowest band frequency)
  - Overlap: 25%
  - Minimum segments: 20 (reduce window length if needed)
- Compute target: new script `src/compute_time_windows.py` (writes window caches)
- Visualization target: new script `src/visualize_time_windows.py`
- Documentation target: `.agents/guides/TIME_WINDOW_GUIDE.md`
- Single-case notebooks added:
  - `ipynb/01_preprocessing/02_time_window_split_singlepat.ipynb`
  - `ipynb/01_preprocessing/03_time_window_corr_singlepat.ipynb`
  - `ipynb/05_figures/08_specific_heat_animation_singlepat.ipynb`
  - `ipynb/05_figures/09_time_window_panels_singlepat.ipynb`
- Validation approach:
  - Dev mode: unvalidated MSC (`sparsify="none"`).
  - Final run: per-window surrogate validation (`sparsify="soft"`).
  - Cache per-window MSC and per-window LRG outputs to avoid recomputation.
- Cache layout (run directories):
  - FC windows: `data/*_cache_windows/Pat_XX/<phase>/<band>/<run_id>/win-0000.npy`
  - LRG windows: `data/lrg_cache_windows/Pat_XX/<phase>/<band>/<run_id>/win-0000_lrg.npz`

T) Band-dependent reorganization summary
- Type: batch script
- Output: `data/figures/summary/` and `results/`
- Implementation target: reuse outputs from Stage 04 and plot summaries

U) Spatial embedding of nodes with cluster coloring (single case)
- Type: notebook (single patient with implant coordinates)
- Use a neuroscience visualization library (e.g. MNE or Nilearn) to plot
  electrode positions and color nodes by LRG clusters.
- Output: `data/figures/spatial/Pat_XX/`
- See `2026-01-10_stage-05u_spatial-embedding.md` for full tasks.

## Status checklist (2026-01-10)
| Figure | Type | Target | Status | Notes |
| --- | --- | --- | --- | --- |
| A | notebook | `ipynb/05_figures/01_timeseries_emd_singlepat.ipynb` | present | Cache-first single-case. |
| B | notebook | `ipynb/05_figures/02_msc_network_singlepat.ipynb` | present | MSC matrix + network single-case. |
| C | script | `src/visualize_msc_grid.py` | present | Phase x band grid from MSC cache. |
| D | script | `src/visualize_msc_all_patients.py` | present | All patients for same band/phase. |
| E | notebook/script | `ipynb/02_fc_msc/05_msc_vs_corr_singlepat.ipynb`, `src/visualize_comparison.py` | present | Uses cached corr + MSC. |
| F | script | `src/compare_network_variants.py` | present | Compare dense MSC vs validated vs corr. |
| G | script | `src/visualize_corr_cleaning.py` | present | MP spectrum + cleaned matrix plots. |
| H | notebook/script | `ipynb/02_fc_msc/04_msc_dense_vs_validated_singlepat.ipynb`, `src/visualize_msc_validation.py` | present | Notebook + batch script. |
| I | notebook/script | `ipynb/03_lrg/01_lrg_singlepat.ipynb`, `src/visualize_lrg.py` | present | LRG summary panels. |
| L | script | `src/visualize_lrg_video.py` | present | Density frames across tau (optional GIF). |
| M | notebook/script | `ipynb/03_lrg/02_lrg_entropy_psi_singlepat.ipynb`, `src/visualize_lrg.py` | present | Full panel includes entropy/PSI network. |
| N | script | `src/visualize_lrg_phase_grid.py` | present | Dendrogram grids across phases. |
| O | script | `src/visualize_reorganization_metrics.py` | present | Heatmaps from `results/reorganization`. |
| P | script | `src/visualize_metric_correlation.py` | present | Cross-metric correlation. |
| Q | notebook/script | `ipynb/05_figures/04_sankey_cluster_transitions.ipynb`, `src/visualize_metastable.py` | present | Metastable Sankey. |
| R | script | `src/aggregate_subject_stats.py` | present | Cross-subject summary tables + heatmaps. |
| S | notebook/script | `ipynb/05_figures/08_specific_heat_animation_singlepat.ipynb`, `src/visualize_time_windows.py` | present | Time-window analysis. |
| T | script | `src/visualize_reorganization_summary.py` | present | Band-dependent summary. |
| U | notebook | `ipynb/05_figures/06_spatial_embedding_singlepat.ipynb` | present | Spatial embedding (matplotlib fallback). |

Bundle status:
- `src/build_overleaf_bundle.py` present (writes manifest + copies).
- `outputs/overleaf/manifest.md` generated by script.

## Deliverables
- A figure generation matrix that tags each item as notebook vs script.
- Batch scripts for grid and cross-patient figures.
- `outputs/overleaf/` folder containing final figures/tables.
- Overleaf manifest: `outputs/overleaf/manifest.md` listing A-U.
- Bundle builder script: `src/build_overleaf_bundle.py` (collects/copies outputs).

## Compute vs visualize
- All figure scripts read cached outputs; they do not recompute FC or LRG.
- Single-case notebooks should call `load_*` and compute only if cache is missing.

## Exit criteria
- All A-T figures can be regenerated from cached outputs.
- Overleaf bundle contains the final figure set and summary tables.
