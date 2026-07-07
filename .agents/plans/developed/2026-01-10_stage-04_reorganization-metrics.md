---
name: stage-04-reorganization-metrics
type: plan
era: IMCOH_SQ
status: superseded
created: 2026-01-10
updated: 2026-04-24
pointers: []
---

# Plan (2026-01-10): Stage 04 - Reorganization metrics

## Goals
- Quantify structural reorganization across phases for each band/patient.
- Keep all existing metrics and compare them against each other.

## Inputs
- Cached LRG results from Stage 03
- `visuals/reorganization.py` (post-refactor: `lrg_eegfc.utils.metrics.reorganization` + `lrg_eegfc.visuals.reorganization`),
  `compare.py` (post-refactor: `lrg_eegfc.utils.metrics.compare`), `utils/distances/`
- Legacy notebooks for distance metrics

## Current status
- Phase reorganization plots exist in `visuals/reorganization.py`.
- Ultrametric comparisons exist in `compare.py`.
- New metric helpers added under `lrg_eegfc.utils.metrics.reorganization`.
- Single-patient + cross-patient notebooks created under `ipynb/04_reorganization/`.
- No consolidated per-metric distance tables yet (script added, not run).

## Tasks
1) Metric inventory and retention
- List every reorganization metric already implemented.
- Keep all of them; no deletions without replacements.

1b) Cross-subject aggregation strategy
- Aggregate per-subject metrics without node-to-node alignment (subjects have
  different electrode placements).
- Use subject-level summaries (per band/phase) for group averages.
- Record any spatial analyses separately (see Stage 05U).
- Primary summary: unweighted mean ± SEM across patients.
- Robust check: median ± IQR across patients (supplementary).
- Optional sensitivity: weighted by channel count (report separately).

2) Cluster correlation coefficient (new)
- Define a node swap coefficient per phase transition:
  - For each consecutive phase pair, compute the fraction of nodes that change
    cluster label after dendrogram cut.
  - Define correlation-style score: 1 - (swapped_nodes / total_nodes).
- Implement as a reusable function and include in summary tables.
- Implementation target: add helper to `visuals/reorganization.py` or
  `utils/clustering/hierarchical.py` and expose via `__init__.py`.

3) ARI/NMI agreement (new)
- Compute ARI (and optionally NMI) for cluster assignments across phases.
- Report alongside the swap coefficient for complementary views.
- Implementation target: `visuals/reorganization.py` (or new module in Stage 06).

4) Distance matrices per metric
- For each patient/band, compute phase-to-phase distance matrices for all
  metrics (ultrametric rank correlation, cophenetic distance, etc.).
- Output as CSV/NPZ under `results/reorganization/Pat_XX/`.
- Implementation target: new script `src/compute_reorganization_metrics.py`
  that writes per-patient CSVs + NPZ matrices.
  - Script added; run to populate `results/reorganization/`.

5) Cross-metric comparisons
- Build a correlation matrix comparing the metrics themselves.
- Summarize consistency between metrics across patients.
- Output as CSV + heatmap under `results/reorganization/` and `data/figures/reorganization/`.
- Implementation target: `src/visualize_metric_correlation.py`

## Compute vs visualize
- Compute: `src/compute_reorganization_metrics.py` writes CSV/NPZ matrices.
- Visualize: `src/visualize_metric_correlation.py` and plotting notebooks only
  read cached metrics and do not recompute LRG.

## Deliverables
- Per-metric distance matrices (per patient/band).
- Cluster correlation coefficient outputs.
- Cross-metric correlation matrix.
- Cross-subject summary tables under `results/reorganization/summary/`.
 - Updated notebooks:
   - `ipynb/04_reorganization/01_distance_metrics_singlepat.ipynb`
   - `ipynb/04_reorganization/02_distance_metrics_crosspat.ipynb`

## Exit criteria
- All metrics run from cached LRG results.
- Outputs are ready for figures O/P/R.
