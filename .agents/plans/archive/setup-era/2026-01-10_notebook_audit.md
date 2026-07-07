---
name: notebook-audit
type: plan
era: IMCOH_SQ
status: dead
created: 2026-01-10
updated: 2026-04-24
pointers: []
---

# Notebook Migration Audit (2026-01-10)

Audit scope: `ipynb/90_archive/.old_reviewed` (32 notebooks). This is a
lightweight mapping to current notebooks/scripts; **partial** means the target
exists but the migration needs verification.

| Legacy notebook | Status | Target(s) |
| --- | --- | --- |
| `DSTCMP_all_distance_measures.ipynb` | partial | `ipynb/04_reorganization/01_distance_metrics_singlepat.ipynb`, `ipynb/04_reorganization/02_distance_metrics_crosspat.ipynb` |
| `DSTCMP_permutation_robust.ipynb` | partial | Stage 04 notebooks + `src/compute_reorganization_metrics.py` |
| `FIGMNTGN01.ipynb` | partial | `ipynb/03_lrg/01_lrg_singlepat.ipynb`, `src/lrg_eegfc/visuals/lrg.py` |
| `FIGMNTGN02.ipynb` | partial | `ipynb/03_lrg/02_lrg_entropy_psi_singlepat.ipynb`, `src/lrg_eegfc/visuals/lrg.py` |
| `FIGMNTGN03.ipynb` | partial | `ipynb/05_figures/04_sankey_cluster_transitions.ipynb`, `src/visualize_metastable.py` |
| `FIGMNTGN04.ipynb` | migrated | `ipynb/05_figures/04_sankey_cluster_transitions.ipynb` |
| `MLB-F01.ipynb` | drop | superseded |
| `NEW_distance_of_distances.ipynb` | partial | Stage 04 notebooks + `src/compute_reorganization_metrics.py` |
| `NEW_standard_disXpat.ipynb` | partial | Stage 04 cross-patient notebook |
| `NEW_tree_measures_pat_comparison.ipynb` | partial | Stage 04 cross-patient notebook |
| `NEW_ultrametric_quantile_rmse_pat_comparison.ipynb` | partial | Stage 04 cross-patient notebook |
| `NEW_ultrametric_rank_correlation_pat_comparison.ipynb` | partial | Stage 04 cross-patient notebook |
| `NEW_ultrametric_scaled_distance_pat_comparison.ipynb` | partial | Stage 04 cross-patient notebook |
| `band_pca.ipynb` | migrated | `ipynb/05_figures/05_band_pca_exploration.ipynb` |
| `distance_of_distances.ipynb` | partial | Stage 04 notebooks |
| `pat02_corrnet_bands.ipynb` | migrated | `ipynb/05_figures/02_corr_networks_per_band_singlepat.ipynb` |
| `pat02_corrnet_emd.ipynb` | migrated | `ipynb/05_figures/01_timeseries_emd_singlepat.ipynb` |
| `poster01.ipynb` | migrated | `ipynb/05_figures/06_poster_layouts.ipynb` |
| `poster02.ipynb` | migrated | `ipynb/05_figures/06_poster_layouts.ipynb` |
| `single_patient_experiments.ipynb` | migrated | `ipynb/05_figures/07_single_patient_experiments.ipynb` |
| `specific_heat_animation.ipynb` | migrated | `ipynb/05_figures/08_specific_heat_animation_singlepat.ipynb` |
| `TEST_distance_of_distances.ipynb` | drop | superseded |
| `TEST_distance_of_distances_2.ipynb` | drop | superseded |
| `TEST_improve_speed.ipynb` | drop | superseded |
| `TEST_interactive_single_patient.ipynb` | drop | superseded |
| `TEST_per_patient_analysis.ipynb` | drop | superseded |
| `TEST_per_patient_time_windows.ipynb` | partial | `ipynb/01_preprocessing/02_time_window_split_singlepat.ipynb`, `ipynb/01_preprocessing/03_time_window_corr_singlepat.ipynb` |
| `tests.ipynb` | partial | Time-window notebooks (Stage S) |
| `tmp.ipynb` | drop | scratch |
| `UTILS-COHERENCE_NETWORKS.ipynb` | partial | `ipynb/02_fc_msc/03_coherence_fc.ipynb`, `ipynb/02_fc_msc/02_build_msc_networks_singlepat.ipynb` |
| `UTILS-FC_COMPARISON_PIPELINE.ipynb` | partial | `ipynb/02_fc_msc/05_msc_vs_corr_singlepat.ipynb`, `src/visualize_comparison.py` |
| `UTILS-FILEREADER.ipynb` | drop | superseded |

Notes:
- `.old` is empty; `.old_reviewed` contains all legacy notebooks listed above.
- "partial" indicates a likely replacement exists, but the extraction needs
  a notebook-by-notebook verification pass.
