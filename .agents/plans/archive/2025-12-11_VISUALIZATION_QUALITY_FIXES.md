---
name: visualization-quality-fixes
type: plan
era: MSC
status: superseded
created: 2025-12-11
updated: 2026-04-24
pointers: []
---

# Visualization Quality Fixes

## Summary

Fixed critical quality issues in MSC and LRG visualizations based on user feedback and correct implementations from FIGMNTGN notebooks.

## Issues Fixed

### 1. MSC Visualization Quality (CRITICAL)

**Problems Identified:**
1. Matrix not displayed with 1:1 aspect ratio
2. Network appeared random despite having weighted structure

**Root Cause:**
- Matrix used `aspect="auto"` instead of `aspect="equal"`
- Spring layout used `k=1.0` which is too large for fully-connected weighted networks

**Fixes Applied:**

**File:** `src/lrg_eegfc/visuals/msc.py`

**Line 181-183** - Fixed matrix aspect ratio:
```python
# BEFORE:
im = ax[0].imshow(msc_matrix, cmap="viridis", vmin=0, vmax=1,
                  interpolation="none", aspect="auto")

# AFTER:
im = ax[0].imshow(msc_matrix, cmap="viridis", vmin=0, vmax=1,
                  interpolation="none", aspect="equal")
```

**Line 206-207** - Fixed network layout for weighted networks:
```python
# BEFORE:
pos = nx.spring_layout(G, seed=42, k=1.0, iterations=50)

# AFTER (k=0.1 reveals structure in fully connected weighted networks):
pos = nx.spring_layout(G, seed=42, k=0.1, iterations=50)
```

**Result:**
- Matrix now displays with proper 1:1 aspect ratio
- Network layout reveals clear community structure

---

### 2. LRG Visualization Quality (CRITICAL OVERHAUL)

**Problems Identified:**
1. Specific heat multiplied by N, causing values to exceed axis bounds
2. Dendrogram settings incorrect (missing log scale, wrong tmin/tmax)
3. Ultrametric matrix plot uninterpretable - should be PSI plot instead
4. Network layout not optimized for weighted networks

**Root Cause:**
- Implementation deviated from correct FIGMNTGN notebook examples
- Missing PSI (Partition Stability Index) computation function
- Incorrect thermodynamic observable scaling

**Fixes Applied:**

**File:** `src/lrg_eegfc/visuals/lrg.py`

#### Added PSI Computation Function (lines 54-107):
```python
def compute_partition_stability_index(linkage_matrix: np.ndarray):
    """Compute Partition Stability Index from hierarchical clustering.

    PSI(n) = N * (log10(Δ_n) - log10(Δ_{n+1}))

    Quantifies stability of partitions with n communities.
    """
    N = linkage_matrix.shape[0] + 1
    deltas = linkage_matrix[:, 2]

    psi_values = []
    n_communities = []

    for i in range(len(deltas) - 1):
        n = N - i - 1
        delta_n = deltas[i]
        delta_n_plus_1 = deltas[i + 1]

        if delta_n <= 0: delta_n = 1e-10
        if delta_n_plus_1 <= 0: delta_n_plus_1 = 1e-10

        psi = N * (np.log10(delta_n) - np.log10(delta_n_plus_1))
        psi_values.append(psi)
        n_communities.append(n)

    # Normalize to [-1, 1]
    psi_array = np.array(psi_values)
    if np.max(np.abs(psi_array)) > 0:
        psi_normalized = -psi_array / np.max(np.abs(psi_array))
    else:
        psi_normalized = psi_array

    return psi_normalized, np.array(n_communities)
```

#### Replaced plot_lrg_full_panel (lines 446-760):

**Change 1 - Panel Layout:** Changed from 4-panel (2×2) to 5-panel layout:
```python
# BEFORE: 2x2 grid
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

# AFTER: Custom 5-panel layout (from FIGMNTGN01)
gs = gridspec.GridSpec(4, 6, figure=fig, hspace=0.3, wspace=0.4)
ax_matrix = fig.add_subplot(gs[0:2, 0:2])   # (a) FC Matrix
ax_entropy = fig.add_subplot(gs[0:2, 2:4])  # (b) Entropy/Heat
ax_dendro = fig.add_subplot(gs[0:4, 4:6])   # (c) Dendrogram
ax_network = fig.add_subplot(gs[2:4, 0:3])  # (d) Network
ax_psi = fig.add_subplot(gs[2:4, 3:4])      # (e) PSI
```

**Change 2 - Specific Heat:** Removed N multiplication (lines 601-606):
```python
# BEFORE (WRONG - caused values to exceed axis):
line_heat = ax_entropy.plot(
    entropy_tau[1:], N * entropy_C, "-",
    label="C (Specific Heat)", color="blue"
)

# AFTER (CORRECT - direct use):
line_heat = ax_entropy.plot(
    entropy_tau[1:], entropy_C, "-",
    label="C (Specific Heat)", color="blue", linewidth=2
)
```

**Change 3 - Dendrogram Settings:** Correct settings from FIGMNTGN03 (lines 655-665):
```python
# BEFORE (missing log scale and correct limits):
ax_dendro.axhline(optimal_th, color="blue", linestyle="--", lw=2)
ax_dendro.set_yscale("log")

# AFTER (CORRECT - from FIGMNTGN03):
tmin = linkage_matrix[:, 2][0] * 0.8
tmax = linkage_matrix[:, 2][-1] * 1.01
ax_dendro.set_xscale("log")
ax_dendro.axvline(
    optimal_threshold, color="b", linestyle="--",
    linewidth=2, label=r"$\mathcal{D}_{\rm th}$"
)
ax_dendro.set_xlim(tmin, tmax)
ax_dendro.set_xlabel(r"$\mathcal{D}/\mathcal{D}_{\max}$", fontsize=11)
```

**Change 4 - Network Layout:** k=0.1 for weighted networks (line 680):
```python
# BEFORE:
pos = nx.spring_layout(G, seed=42, k=1.0, iterations=50)

# AFTER:
pos = nx.spring_layout(G, seed=43, scale=1, k=0.1, iterations=50)
```

**Change 5 - PSI Plot:** Replaced ultrametric matrix (lines 710-740):
```python
# BEFORE: Ultrametric heatmap (uninterpretable)
ultrametric_square = squareform(result.ultrametric_matrix)
im = ax.imshow(ultrametric_square, cmap="viridis")

# AFTER: PSI Plot (interpretable stability metric)
psi_values, n_communities = compute_partition_stability_index(linkage_matrix)
ax_psi.plot(n_communities, psi_values, "-o", color="green", linewidth=2)

# Mark optimal partition from PSI
max_psi_idx = np.argmax(-psi_values)
optimal_n_communities = n_communities[max_psi_idx]
ax_psi.axvline(optimal_n_communities, ls="--", c="red",
               label=f"PSI optimal: n={optimal_n_communities}")
ax_psi.axvline(n_clusters, ls=":", c="blue",
               label=f"Used: n={n_clusters}")
```

**Result:**
- Specific heat now displays correctly within axis bounds
- Dendrogram uses correct log scale and limits from FIGMNTGN03
- PSI plot provides interpretable partition stability information
- Network layout reveals community structure

---

## References

### Correct Implementations
- **FIGMNTGN01.ipynb**: 5-panel layout, PSI calculation
- **FIGMNTGN03.ipynb**: Dendrogram settings (tmin, tmax, log scale, k=0.1)
- **LRG_VISUALIZATION_INSTRUCTIONS.md**: Comprehensive requirements

### Key Parameters
- **MSC network layout:** `k=0.1` for fully-connected weighted networks
- **LRG network layout:** `k=0.1` for weighted networks
- **Dendrogram limits:** `tmin = linkage[:, 2][0] * 0.8`, `tmax = linkage[:, 2][-1] * 1.01`
- **Specific heat:** NO multiplication by N or log(N)
- **Matrix aspect:** `aspect="equal"` for 1:1 ratio

---

## Testing

### Test Commands
```bash
# Test MSC visualization (single)
python src/visualize_msc.py --patient Pat_03 --phase rsPre --band beta --plot-type summary --verbose

# Test LRG visualization (single)
python src/visualize_lrg.py --patient Pat_03 --phase rsPre --band beta --fc-method corr --plot-type full --verbose

# Test MSC batch (all bands, all phases)
python src/visualize_msc.py --patient Pat_03 --batch --plot-type summary --verbose

# Test LRG batch (all bands, all phases)
python src/visualize_lrg.py --patient Pat_03 --fc-method corr --batch --plot-type full --verbose
```

### Verification Results
✅ MSC: Matrix 1:1 aspect ratio, network shows clear structure
✅ LRG: 5 panels with correct specific heat, dendrogram, and PSI plot
✅ Batch mode: All 24 visualizations (6 bands × 4 phases) generated successfully

---

## Impact

### Before Fixes
- MSC: Distorted matrix, random-looking network
- LRG: Specific heat out of bounds, missing PSI analysis, wrong dendrogram scale

### After Fixes
- MSC: Professional publication-quality visualizations
- LRG: Complete hierarchical analysis with interpretable stability metrics
- Both: Network layouts reveal actual community structure

---

## Files Modified

1. `src/lrg_eegfc/visuals/msc.py` (lines 181-183, 206-207)
2. `src/lrg_eegfc/visuals/lrg.py` (complete overhaul of plot_lrg_full_panel, added compute_partition_stability_index)

## Backup Files

- `src/lrg_eegfc/visuals/lrg_backup.py` - Original implementation
- `src/lrg_eegfc/visuals/lrg_revised.py` - Development version with corrections

---

## Next Steps

1. Generate visualizations for all patients using corrected implementations
2. Verify PSI analysis provides meaningful partition stability insights
3. Document interpretation guidelines for PSI plots
