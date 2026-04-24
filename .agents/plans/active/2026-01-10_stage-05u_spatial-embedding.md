---
name: stage-05u-spatial-embedding
type: plan
era: IMCOH_SQ
status: dead
created: 2026-01-10
updated: 2026-04-24
pointers: []
---

# Plan (2026-01-10): Stage 05U - Spatial embedding and cluster coloring

## Goals
- Visualize electrode nodes in anatomical space for a single patient.
- Color nodes by LRG cluster membership for interpretability.
- Keep this as a notebook-only step (no batch run yet).

## Inputs
- `data/stereoeeg_patients/Pat_XX/Implant_pat_XX.csv`
- `data/stereoeeg_patients/Pat_XX/channel_labels.csv`
- Cached LRG results (Stage 03)

## Tasks
1) Patient selection
- Choose a patient with implant coordinates (from Stage 01 report).
- Confirm channel labels align between FC matrices and implant metadata.

2) Spatial plotting notebook
- Create `ipynb/05_figures/06_spatial_embedding_singlepat.ipynb`.
- Use a neuroscience library (prefer MNE; fallback Nilearn) to render
  electrode positions in 3D and color by cluster assignment.
- Add a second view with clusters from different phases to show reorganization.

3) Cluster mapping
- Extract cluster labels from LRG dendrogram cut (use optimal threshold).
- Map labels to channel names, then to implant coordinates.

4) Outputs
- Save static figures under `data/figures/spatial/Pat_XX/`.
- Optional: save interactive HTML (if supported by the library).

## Compute vs visualize
- This stage is visualization-only.
- Use cached LRG outputs and implant metadata; do not recompute FC or LRG.

## Deliverables
- Notebook showing spatial embedding with cluster coloring.
- A small helper function in `src/lrg_eegfc/visuals/` if reusable.

## Exit criteria
- At least one patient has a spatially embedded network plot
  with phase-specific cluster coloring.
