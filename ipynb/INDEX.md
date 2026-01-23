# Notebook Index

Legacy migration audit: `.agents/plans/active/2026-01-10_notebook_audit.md`.

All notebooks are listed here with a stable ID. Use this index to track consolidation and migration into reusable modules.

| ID | Path | Status | Summary |
| --- | --- | --- | --- |
| NB-001 | `ipynb/01_preprocessing/time_windows_analysis.ipynb` | active | Time window sensitivity checks for SEEG data. |
| NB-002 | `ipynb/02_fc_msc/01_build_corr_networks_singlepat.ipynb` | active | Dev notebook to build correlation FC matrices for a single patient. |
| NB-003 | `ipynb/02_fc_msc/02_build_msc_networks_singlepat.ipynb` | active | Dev notebook to build MSC FC matrices for a single patient. |
| NB-004 | `ipynb/02_fc_msc/03_coherence_fc.ipynb` | active | MSC/Coherence functional connectivity tutorial. |
| NB-005 | `ipynb/04_reorganization/distance_measures_comparison.ipynb` | active | Compare distance measures for reorganization analyses. |
| NB-041 | `ipynb/04_reorganization/01_distance_metrics_singlepat.ipynb` | active | Single-patient reorganization metrics (phase distance matrices). |
| NB-042 | `ipynb/04_reorganization/02_distance_metrics_crosspat.ipynb` | active | Cross-patient metric ranking and consistency checks. |
| NB-039 | `ipynb/03_lrg/01_lrg_singlepat.ipynb` | active | Single-patient LRG full-panel check (cached inputs only). |
| NB-040 | `ipynb/03_lrg/02_lrg_entropy_psi_singlepat.ipynb` | active | Single-patient LRG entropy + dendrogram checks (cached inputs only). |
| NB-038 | `ipynb/00_intake/01_sanity_imports.ipynb` | active | Sanity check for imports, data loading, and cached outputs. |
| NB-043 | `ipynb/05_figures/01_timeseries_emd_singlepat.ipynb` | active | Figure A: timeseries + EMD example (single patient). |
| NB-044 | `ipynb/05_figures/02_corr_networks_per_band_singlepat.ipynb` | active | Figure B/C: correlation matrices across bands (single patient). |
| NB-045 | `ipynb/05_figures/03_lrg_full_panel_singlepat.ipynb` | active | Figures M/N: LRG full panel (single patient). |
| NB-046 | `ipynb/05_figures/04_sankey_cluster_transitions.ipynb` | active | Figure Q: Sankey for cluster transitions across tau. |
| NB-047 | `ipynb/05_figures/05_band_pca_exploration.ipynb` | active | Optional: band PCA exploration (single patient). |
| NB-048 | `ipynb/05_figures/06_poster_layouts.ipynb` | active | Poster layout scaffold combining figures. |
| NB-049 | `ipynb/05_figures/07_single_patient_experiments.ipynb` | active | Single-patient sandbox for cached FC checks. |
| NB-050 | `ipynb/02_fc_msc/04_msc_dense_vs_validated_singlepat.ipynb` | active | Dense vs validated MSC comparison (single patient). |
| NB-051 | `ipynb/02_fc_msc/05_msc_vs_corr_singlepat.ipynb` | active | MSC vs correlation comparison (single patient). |
| NB-052 | `ipynb/00_intake/02_data_reader_singlepat.ipynb` | active | Data reader walkthrough (single patient). |
| NB-053 | `ipynb/01_preprocessing/02_time_window_split_singlepat.ipynb` | active | Time-window split setup (single patient). |
| NB-054 | `ipynb/01_preprocessing/03_time_window_corr_singlepat.ipynb` | active | Correlation matrices across time windows (single patient). |
| NB-055 | `ipynb/05_figures/08_specific_heat_animation_singlepat.ipynb` | active | Specific heat animation across time windows. |
| NB-056 | `ipynb/05_figures/09_time_window_panels_singlepat.ipynb` | active | Time-window 4-panel summary (corr/dendrogram/heat/network). |
| NB-057 | `ipynb/03_lrg/03_corr_threshold_jumps_singlepat.ipynb` | active | Threshold jumps + dendrogram (single patient). |
| NB-058 | `ipynb/03_lrg/04_interactive_dendrogram_singlepat.ipynb` | active | Interactive dendrogram threshold slider. |
| NB-059 | `ipynb/dev/01_speed_smoke.ipynb` | active | Dev speed smoke test (corr + MSC). |
| NB-060 | `ipynb/dev/02_graph_sandbox.ipynb` | active | Graph sandbox for thresholding + entropy sanity checks. |
| NB-061 | `ipynb/05_figures/02_msc_network_singlepat.ipynb` | active | Figure B: MSC network + FC matrix (single patient). |
| NB-062 | `ipynb/05_figures/06_spatial_embedding_singlepat.ipynb` | active | Figure U: spatial embedding with cluster coloring (single patient). |
| NB-006 | `ipynb/90_archive/.old_reviewed/DSTCMP_all_distance_measures.ipynb` | merge | Distance-of-distances (all measures) -> Stage-04 merge. |
| NB-007 | `ipynb/90_archive/.old_reviewed/DSTCMP_permutation_robust.ipynb` | merge | Permutation-robust ultrametric distances -> Stage-04 merge. |
| NB-008 | `ipynb/90_archive/.old_reviewed/FIGMNTGN01.ipynb` | merge | Figure prototype (dendrogram/entropy/PSI) -> Stage-03/05. |
| NB-009 | `ipynb/90_archive/.old_reviewed/FIGMNTGN02.ipynb` | merge | Figure prototype (dendrogram/entropy) -> Stage-03/05. |
| NB-010 | `ipynb/90_archive/.old_reviewed/FIGMNTGN03.ipynb` | merge | Figure prototype (sankey/metastable) -> Stage-05. |
| NB-011 | `ipynb/90_archive/.old_reviewed/FIGMNTGN04.ipynb` | keep | Sankey diagram for cluster transitions -> Stage-05. |
| NB-012 | `ipynb/90_archive/.old_reviewed/MLB-F01.ipynb` | drop | Data loader prototype (superseded). |
| NB-013 | `ipynb/90_archive/.old_reviewed/NEW_distance_of_distances.ipynb` | keep | Comprehensive distance-of-distances pipeline -> Stage-04. |
| NB-014 | `ipynb/90_archive/.old_reviewed/NEW_standard_disXpat.ipynb` | keep | Ultrametric matrix distances across patients -> Stage-04. |
| NB-015 | `ipynb/90_archive/.old_reviewed/NEW_tree_measures_pat_comparison.ipynb` | keep | Tree-based distance comparisons -> Stage-04. |
| NB-016 | `ipynb/90_archive/.old_reviewed/NEW_ultrametric_quantile_rmse_pat_comparison.ipynb` | keep | Ultrametric quantile RMSE -> Stage-04. |
| NB-017 | `ipynb/90_archive/.old_reviewed/NEW_ultrametric_rank_correlation_pat_comparison.ipynb` | keep | Ultrametric rank-correlation distances -> Stage-04. |
| NB-018 | `ipynb/90_archive/.old_reviewed/NEW_ultrametric_scaled_distance_pat_comparison.ipynb` | keep | Scaled ultrametric distances -> Stage-04. |
| NB-019 | `ipynb/90_archive/.old_reviewed/TEST_distance_of_distances.ipynb` | drop | Test notebook (duplicate distance-of-distances). |
| NB-020 | `ipynb/90_archive/.old_reviewed/TEST_distance_of_distances_2.ipynb` | drop | Test notebook (duplicate distance-of-distances). |
| NB-021 | `ipynb/90_archive/.old_reviewed/TEST_improve_speed.ipynb` | drop | Speed test notebook. |
| NB-022 | `ipynb/90_archive/.old_reviewed/TEST_interactive_single_patient.ipynb` | drop | Interactive test notebook. |
| NB-023 | `ipynb/90_archive/.old_reviewed/TEST_per_patient_analysis.ipynb` | drop | Test per-patient analysis. |
| NB-024 | `ipynb/90_archive/.old_reviewed/TEST_per_patient_time_windows.ipynb` | merge | Time-window experiments -> Stage S. |
| NB-025 | `ipynb/90_archive/.old_reviewed/UTILS-COHERENCE_NETWORKS.ipynb` | merge | MSC/Coherence pipeline notes -> Stage-02. |
| NB-026 | `ipynb/90_archive/.old_reviewed/UTILS-FC_COMPARISON_PIPELINE.ipynb` | merge | Correlation vs MSC pipeline -> Stage-02/05. |
| NB-027 | `ipynb/90_archive/.old_reviewed/UTILS-FILEREADER.ipynb` | drop | Data reader prototype (superseded). |
| NB-028 | `ipynb/90_archive/.old_reviewed/band_pca.ipynb` | merge | Band PCA exploration -> optional Stage-05. |
| NB-029 | `ipynb/90_archive/.old_reviewed/distance_of_distances.ipynb` | merge | Early distance-of-distances version -> Stage-04. |
| NB-030 | `ipynb/90_archive/.old_reviewed/pat02_corrnet_bands.ipynb` | merge | Correlation networks across bands -> Stage-05. |
| NB-031 | `ipynb/90_archive/.old_reviewed/pat02_corrnet_emd.ipynb` | merge | EMD + correlation network (figure A) -> Stage-05. |
| NB-032 | `ipynb/90_archive/.old_reviewed/poster01.ipynb` | merge | Poster figure prototype -> Stage-05. |
| NB-033 | `ipynb/90_archive/.old_reviewed/poster02.ipynb` | merge | Poster figure prototype -> Stage-05. |
| NB-034 | `ipynb/90_archive/.old_reviewed/single_patient_experiments.ipynb` | merge | Single-patient experiments -> Stage-02/05. |
| NB-035 | `ipynb/90_archive/.old_reviewed/specific_heat_animation.ipynb` | keep | Specific heat animation -> Stage S. |
| NB-036 | `ipynb/90_archive/.old_reviewed/tests.ipynb` | merge | Time-window split tests -> Stage S. |
| NB-037 | `ipynb/90_archive/.old_reviewed/tmp.ipynb` | drop | Scratch notebook. |
