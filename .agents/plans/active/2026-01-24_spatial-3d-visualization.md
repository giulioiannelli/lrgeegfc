# Plan: 3D Spatial Brain Network Visualization

## Goal
Create visualization tools to display SEEG networks with nodes positioned at their real 3D anatomical coordinates (from electrode implant data), showing network edges/connections in spatial context.

---

## Key Design Decisions

| Aspect | Decision | Rationale |
|--------|----------|-----------|
| **Primary library** | Plotly 3D | Already in project, interactive, excellent for 3D networks |
| **Fallback** | Matplotlib 3D | Static PNG for publications |
| **Coordinates** | Native patient space (mm) | No MNI transform needed; data in micrometers, divide by 1000 |
| **Edge rendering** | Straight lines + opacity by weight | Performance, clarity |
| **Node coloring** | LRG clusters or atlas regions | Integrate with existing analysis |

---

## Data Structure (from Implant_pat_*.csv)

```
label,x,y,z,Desikan-Killany
A1,-20528,-53975,-41970," ctx-lh-fusiform,80.0,..."
```
- Coordinates in micrometers (divide by 1000 for mm)
- Desikan-Killiany atlas with probability scores
- Existing loader: `load_patient_metadata()` in `src/lrg_eegfc/utils/io/patient.py`

---

## Files to Create/Modify

### 1. New Module: `src/lrg_eegfc/visuals/spatial.py`

**Functions:**
- `prepare_spatial_coordinates(metadata, scale="mm", center=True)` - Extract/transform coords
- `build_edge_traces(coords, adjacency, threshold, max_edges)` - Create Plotly edge lines
- `plot_spatial_network_3d(patient, phase, band, ...)` - Main interactive 3D plot
- `plot_spatial_network_3d_mpl(...)` - Matplotlib fallback for static PNG
- `plot_spatial_clusters_comparison(patient, phase_a, phase_b, ...)` - Side-by-side 3D comparison

### 2. Update: `src/lrg_eegfc/visuals/__init__.py`
Add new exports for spatial functions.

### 3. New Notebooks

| Notebook | Purpose |
|----------|---------|
| `ipynb/05_figures/10_spatial_network_3d.ipynb` | Demo: 3D network with cluster coloring, edge filtering |
| `ipynb/05_figures/11_spatial_cluster_comparison.ipynb` | Compare phases in 3D space |

---

## Implementation Steps

1. Create `spatial.py` module with coordinate utilities
2. Implement `build_edge_traces()` for efficient edge rendering
3. Implement `plot_spatial_network_3d()` with Plotly
4. Implement `plot_spatial_network_3d_mpl()` matplotlib fallback
5. Implement `plot_spatial_clusters_comparison()` for phase comparison
6. Update `visuals/__init__.py` with new exports
7. Create demo notebook with interactive examples
8. Create comparison notebook for pre/post analysis

---

## Key Implementation Details

**Edge trace efficiency** (single trace with None separators):
```python
x_edges, y_edges, z_edges = [], [], []
for i, j in edges:
    x_edges.extend([coords[i,0], coords[j,0], None])
    ...
```

**Node hover info**: electrode label + cluster ID + atlas region

**Edge filtering**: threshold by weight + max_edges limit for performance

---

## Verification

1. Run demo notebook with Pat_02/rsPre/beta
2. Verify interactive rotation/zoom works
3. Confirm cluster colors match LRG dendrogram
4. Test edge threshold slider
5. Export HTML (interactive) and PNG (static)
6. Compare pre vs post phases in 3D

---

## Dependencies

No new packages needed. Uses existing:
- `plotly.graph_objects` (already in project)
- `matplotlib` with `mpl_toolkits.mplot3d`
- `pandas`, `numpy`, `networkx`
