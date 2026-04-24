# Function Map - lrg_eegfc Package

Quick reference for all public functions organized by module.

---

## 1. VISUALIZATION (`lrg_eegfc.visuals`)

### Correlation (`visuals.correlation`)
| Function | Purpose |
|----------|---------|
| `plot_correlation_heatmap(corr_matrix, output_path, ...)` | Single correlation matrix heatmap |
| `plot_correlation_summary(patient, phase, band, ...)` | 4-panel: raw, abs, thresholded, network |
| `plot_correlation_and_network(patient, phase, band, ...)` | Side-by-side matrix + network graph |
| `plot_marchenko_pastur_comparison(patient, phase, band, ...)` | 3-panel: original, cleaned, eigenvalue distribution |
| `plot_percolation_curves(patient, phase, band, ...)` | P_inf and E_inf vs threshold |

### MSC (`visuals.msc`)
| Function | Purpose |
|----------|---------|
| `plot_msc_heatmap(msc_matrix, output_path, ...)` | Single MSC matrix heatmap |
| `plot_msc_and_network(patient, phase, band, ...)` | MSC matrix + network graph |
| `plot_msc_comparison_dense_vs_validated(patient, phase, band, ...)` | 3-panel: dense, validated, difference |
| `plot_msc_summary(patient, phase, band, ...)` | 3-panel: heatmap, network, percolation |

### LRG (`visuals.lrg`)
| Function | Purpose |
|----------|---------|
| `plot_lrg_entropy_curves(patient, phase, band, fc_method, ...)` | Entropy (1-S, C) vs tau |
| `plot_lrg_dendrogram(patient, phase, band, fc_method, ...)` | Hierarchical dendrogram with threshold |
| `plot_ultrametric_heatmap(patient, phase, band, fc_method, ...)` | Ultrametric distance matrix |
| `plot_lrg_full_panel(patient, phase, band, fc_method, ...)` | 5-panel comprehensive analysis |
| `compute_partition_stability_index(linkage_matrix)` | Compute PSI (Ψ) from clustering |

### Spatial 3D (`visuals.spatial`)
| Function | Purpose |
|----------|---------|
| `load_spatial_metadata(patient, dataset_root)` | Load electrode 3D coordinates |
| `plot_spatial_network_3d(patient, phase, band, ...)` | Interactive Plotly 3D network |
| `plot_spatial_network_3d_mpl(patient, phase, band, ...)` | Static Matplotlib 3D network |
| `plot_spatial_clusters_comparison(patient, phase_a, phase_b, ...)` | Side-by-side 3D phase comparison |
| `view_brain_connectome(patient, phase, band, ...)` | Interactive nilearn glass brain |
| `plot_brain_connectome(patient, phase, band, ...)` | Static nilearn multi-view |

### Reorganization (`visuals.reorganization`)
| Function | Purpose |
|----------|---------|
| `plot_phase_reorganization(patient, band, fc_method, ...)` | Multi-phase network/ultrametric/dendrogram |
| `plot_reorganization_distance_matrix(patient, band, fc_method, ...)` | Pairwise phase distance matrix |
| `compute_cluster_membership_correlation(linkage_1, linkage_2, ...)` | Correlation between hierarchies |

### Plotting Helpers (`visuals.plotting`)
| Function | Purpose |
|----------|---------|
| `plot_correlation_matrix(matrix, output_path, ...)` | Basic matrix plot with coolwarm |
| `plot_entropy(graph, output_path, ...)` | Entropy curves from LRG graph |
| `prepare_dendrogram(graph)` | Compute linkage matrix |
| `plot_dendrogram(linkage_matrix, labels, ...)` | Dendrogram with log scale |
| `plot_graph(graph, dendro, channel_map, ...)` | Network with dendrogram colors |

---

## 2. WORKFLOWS (`lrg_eegfc.workflow`)

### Correlation (`workflow.corr`)
| Function | Returns | Purpose |
|----------|---------|---------|
| `get_corr_cache_path(patient, phase, band, ...)` | `Path` | Cache file path |
| `load_corr_matrix(patient, phase, band, ...)` | `ndarray \| None` | Load from cache |
| `compute_corr_matrix(patient, phase, band, ...)` | `CorrResult` | Compute and cache |
| `compute_corr_for_patient(patient, ...)` | `dict` | All phases/bands for patient |

### MSC (`workflow.msc`)
| Function | Returns | Purpose |
|----------|---------|---------|
| `get_msc_cache_path(patient, phase, band, ...)` | `Path` | Cache file path |
| `load_msc_matrix(patient, phase, band, ...)` | `ndarray \| None` | Load from cache |
| `compute_msc_matrix(patient, phase, band, ...)` | `MSCResult` | Compute and cache |
| `compute_msc_for_patient(patient, ...)` | `dict` | All phases/bands for patient |

### LRG (`workflow.lrg`)
| Function | Returns | Purpose |
|----------|---------|---------|
| `get_lrg_cache_path(patient, phase, band, ...)` | `Path` | Cache file path |
| `load_lrg_result(patient, phase, band, fc_method, ...)` | `LRGResult \| None` | Load from cache |
| `compute_lrg_analysis(fc_matrix, ...)` | `LRGResult` | Compute LRG from FC matrix |
| `compute_lrg_for_patient(patient, fc_method, ...)` | `dict` | All phases/bands for patient |

### Cleaning (`workflow.cleaning`)
| Function | Returns | Purpose |
|----------|---------|---------|
| `clean_correlation_matrix_full(corr_matrix, n_channels)` | `ndarray` | Marchenko-Pastur cleaning |
| `load_cleaned_corr_matrix(patient, phase, band, ...)` | `ndarray \| None` | Load cleaned from cache |
| `get_cleaned_corr_cache_path(patient, phase, band, ...)` | `Path` | Cache file path |
| `compute_cleaned_corr_for_patient(patient, ...)` | `dict` | All cleaned matrices |

### Time Windows (`workflow.time_windows`)
| Function | Purpose |
|----------|---------|
| `build_window_run_id(window_sec, overlap, ...)` | Unique run ID for params |
| `get_window_cache_dir(run_id, cache_root, fc_method)` | Cache directory for windowed FC |
| `suggest_window_sec(phase_duration, cycles)` | Suggest window size |
| `compute_window_params(n_samples, sample_rate, ...)` | Window params (step, n_windows) |
| `generate_window_indices(n_samples, window_size, ...)` | Window start/stop indices |

---

## 3. DATA LOADING (`lrg_eegfc.utils.io`)

### Loaders (`utils.io.loaders`)
| Function | Purpose |
|----------|---------|
| `load_mat_pat_data(patient, phase, root_path)` | Load patient phase data from .mat |
| `load_data_dict(mat_path)` | Load raw .mat file as dict |

### Patient Data (`utils.io.patient`)
| Function | Purpose |
|----------|---------|
| `load_timeseries(patient, phase, root_path)` | Load (channels × time) array |
| `load_patient_metadata(patient, root_path)` | Load metadata dataframe |
| `load_patient_dataset(patient, root_path, phases)` | Load all recordings for patient |
| `load_dataset(root_path, patients, phases)` | Load entire dataset |

### Robust Loaders (`utils.io.patient_robust`)
| Function | Purpose |
|----------|---------|
| `load_timeseries_robust(patient, phase, root_path)` | Load with scipy/h5py fallback |
| `load_patient_metadata_robust(patient, root_path)` | Robust metadata loading |
| `load_patient_dataset_robust(patient, root_path, phases)` | Robust full patient load |

### Inspection (`utils.io.inspect`)
| Function | Purpose |
|----------|---------|
| `inspect_mat_file(mat_path)` | Inspect .mat file structure |
| `inspect_patient(patient, root_path)` | Inspect all patient recordings |
| `inspect_all_patients(root_path, patients)` | Inspect entire dataset |
| `generate_report(results)` | Generate text report |
| `save_report(results, output_path)` | Save report to file |

---

## 4. FC UTILITIES (`lrg_eegfc.utils.fc`)

### Correlation Base (`utils.fc.corr.base`)
| Function | Purpose |
|----------|---------|
| `apply_threshold_filter(matrix, threshold)` | Zero values below threshold |
| `build_corr_network(data, filter_type, zero_diagonal)` | Build correlation network |
| `clean_correlation_matrix(X, rowvar)` | Clean NaN/inf values |

### Correlation Bands (`utils.fc.corr.bands`)
| Function | Purpose |
|----------|---------|
| `build_corrmat_perband(patient, phase, band, ...)` | Correlation matrix for band |
| `build_band_correlation_matrices(patient, phase, ...)` | Matrices for multiple bands |

### Correlation Thresholds (`utils.fc.corr.thresholds`)
| Function | Purpose |
|----------|---------|
| `find_exact_detachment_threshold(corr_mat)` | First node detachment threshold |
| `find_threshold_jumps(G, return_stats)` | Percolation threshold jumps |

### MSC Core (`utils.fc.msc.msc`)
| Function | Purpose |
|----------|---------|
| `compute_msc_welch(data, nperseg, ...)` | MSC via Welch's method |
| `band_average_msc(msc_full, band_freqs)` | Average MSC over band |

### MSC Sparsification (`utils.fc.msc.sparsify`)
| Function | Purpose |
|----------|---------|
| `soft_sparsify_surrogate(msc_matrix, surrogates, threshold)` | Surrogate-based sparsification |

### MSC Surrogates (`utils.fc.msc.surrogates`)
| Function | Purpose |
|----------|---------|
| `circular_shift_surrogates(data, n_surrogates, seed)` | Generate surrogate data |
| `surrogate_msc_null(data, n_surrogates, ...)` | MSC null distribution |

### MSC Pipeline (`utils.fc.msc`)
| Function | Purpose |
|----------|---------|
| `coherence_fc_pipeline(data, sample_rate, bands, ...)` | Full MSC pipeline |

---

## 5. METRICS (`lrg_eegfc.utils.metrics`)

### Ultrametric Comparison (`utils.metrics.compare`)
| Function | Purpose |
|----------|---------|
| `compare_ultrametric_matrices(u1, u2)` | Compare with multiple distance metrics |
| `compare_fc_methods(patient, phase, band, ...)` | Compare corr vs MSC |
| `compare_phases(patient, band, fc_method, ...)` | Compare across phases |
| `batch_compare_fc_methods(patients, bands, ...)` | Batch FC method comparison |
| `batch_compare_phases(patients, bands, fc_method, ...)` | Batch phase comparison |

### Comparison Metrics (`utils.metrics.comparison`)
| Function | Purpose |
|----------|---------|
| `compute_phase_distance_matrix(patient, band, ...)` | Pairwise phase distances |
| `compute_cross_patient_consistency(patients, band, ...)` | Cross-patient consistency |
| `rank_distance_measures(distance_dict, ...)` | Rank by value |

### Reorganization (`utils.metrics.reorganization`)
| Function | Purpose |
|----------|---------|
| `compute_metric_matrix(phase_results, metric_key)` | Distance matrix for metric |
| `compute_cluster_labels(linkage, threshold)` | Extract cluster labels |
| `cluster_swap_coefficient(labels_a, labels_b)` | Cluster swap coefficient |
| `compute_cluster_swap_matrix(phase_results)` | Swap matrix across phases |
| `compute_ari_matrix(phase_results)` | Adjusted Rand index matrix |

---

## 6. LRG UTILITIES (`lrg_eegfc.utils.lrg`)

### Hierarchical (`utils.lrg.hierarchical`)
| Function | Purpose |
|----------|---------|
| `fcluster_with_outliers(linkage, threshold)` | Cluster with outlier detection |
| `get_dendrogram_consistent_clusters(linkage, ...)` | Dendrogram-consistent clusters |
| `compute_optimal_clusters_auto(linkage, ...)` | Auto-determine cluster count |
| `get_outlier_nodes(clusters)` | Extract outlier nodes |
| `compute_cluster_statistics(clusters)` | Basic cluster statistics |

---

## 7. CONFIGURATION (`lrg_eegfc.config`)

### Constants (`config.const`)
| Constant | Value | Purpose |
|----------|-------|---------|
| `sEEG_DATAPATH` | `Path("data/stereoeeg_patients")` | Default data path |
| `PHASE_LABELS` | `('rest_pre', 'task_learn', 'task_test', 'rest_post')` | Recording phases |
| `BRAIN_BANDS` | `{'delta': (1,4), 'theta': (4,8), ...}` | Frequency bands |
| `BRAIN_BANDS_NAMES` | `['delta', 'theta', 'alpha', ...]` | Band names list |
| `BRAIN_BAND_TEX_DICT` | `{'delta': r'$\delta$', ...}` | LaTeX labels |
| `DEFAULT_N_SURROGATES` | `200` | MSC surrogates default |
| `DEFAULT_SAMPLE_RATE` | `2048.0` | Sample rate (Hz) |

| Function | Purpose |
|----------|---------|
| `list_patients(dataset_root)` | List available patient IDs |

---

## 8. RESULT DATACLASSES

| Class | Module | Fields |
|-------|--------|--------|
| `CorrResult` | `workflow.corr` | adjacency_matrix, graph, patient, phase, band, n_channels, mean_corr, filter_type, zero_diagonal, filter_order |
| `MSCResult` | `workflow.msc` | adjacency_matrix, graph, patient, phase, band, n_channels, mean_msc, sparsify, n_surrogates, nperseg |
| `LRGResult` | `workflow.lrg` | linkage_matrix, ultrametric_matrix, optimal_threshold, entropy_tau, entropy_1_minus_S, entropy_C, n_nodes, patient, phase, band, fc_method |
| `CleanedCorrResult` | `workflow.cleaning` | adjacency_matrix, original_matrix, graph, metadata |
| `PatientRecording` | `utils.io.patient` | timeseries, parameters, channel_labels |

---

## Quick Import Patterns

```python
# Standard notebook header
from lrg_eegfc.notebook import *
move_to_root(pathname="lrgeegfc")

# Visualization imports
from lrg_eegfc.visuals import (
    plot_correlation_and_network,
    plot_lrg_full_panel,
    plot_spatial_network_3d,
    plot_phase_reorganization,
)

# Workflow imports
from lrg_eegfc.workflow import (
    compute_corr_matrix,
    compute_msc_matrix,
    compute_lrg_analysis,
    load_lrg_result,
)

# Configuration
from lrg_eegfc.config import BRAIN_BANDS, PHASE_LABELS, DEFAULT_SAMPLE_RATE
```
