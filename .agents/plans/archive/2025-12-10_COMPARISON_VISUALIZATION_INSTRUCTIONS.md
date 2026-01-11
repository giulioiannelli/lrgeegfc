# Instructions: Implement Comparison Visualizations for FC Methods

**Date:** 2025-12-10
**Task:** Implement comparison visualizations between Correlation-based FC and MSC-based FC
**Goal:** Create minimal, summarizing figures for comparing the two functional connectivity methods

---

## Overview

You will implement a comparison visualization function that compares two functional connectivity methods:
1. **Correlation-based FC**: Pearson correlation on bandpass-filtered time series
2. **MSC-based FC**: Magnitude-squared coherence using Welch's method

The visualization should be **minimal and summarizing** - one comprehensive figure per patient/phase/band combination.

---

## Reference Implementation

### Existing Notebook Analysis

File: `ipynb/UTILS-FC_COMPARISON_PIPELINE.ipynb`

**Current visualizations:**
1. **6x4 grid (24 panels)** showing all bands with:
   - Column 1: Coherence heatmaps
   - Column 2: Correlation heatmaps
   - Column 3: Scatter plots (Coherence vs Correlation)
   - Column 4: Distribution histograms

2. **Distribution comparison** for single band (dense vs sparse)

**Key observations from notebook:**
- Both matrices are symmetric, range [0,1] for coherence, [-1,1] for correlation
- Use absolute correlation for fair comparison: `np.abs(corr_matrices[band])`
- Compute agreement using Pearson correlation of upper-triangle values
- Common threshold for sparsity analysis: 0.01

---

## Implementation Requirements

### File Structure

Create new file: `src/lrg_eegfc/visuals/compare.py`

```python
"""Comparison visualizations for functional connectivity methods."""

from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

__all__ = ["plot_fc_comparison"]
```

### Function Signature

```python
def plot_fc_comparison(
    patient: str,
    phase: str,
    band: str,
    cache_root_corr: Path = Path("data/corr_cache"),
    cache_root_msc: Path = Path("data/msc_cache"),
    output_path: Optional[Path] = None,
    nperseg: int = 1024,
    figsize: tuple = (18, 6),
) -> Path:
    """Create comparison visualization between Correlation and MSC FC methods.

    Generates a 3-panel horizontal layout comparing two FC methods:
    - Panel 0: Correlation FC matrix (absolute value)
    - Panel 1: MSC FC matrix
    - Panel 2: Method agreement analysis (scatter + distributions)

    Parameters
    ----------
    patient : str
        Patient identifier (e.g., "Pat_02")
    phase : str
        Recording phase (e.g., "rsPre", "taskLearn")
    band : str
        Frequency band (e.g., "beta", "alpha")
    cache_root_corr : Path
        Root directory for correlation cache
    cache_root_msc : Path
        Root directory for MSC cache
    output_path : Path, optional
        Output file path. If None, uses default location
    nperseg : int
        Window length for MSC (must match cached files)
    figsize : tuple
        Figure size (width, height)

    Returns
    -------
    Path
        Path to saved figure

    Raises
    ------
    FileNotFoundError
        If cache files not found for either method

    Notes
    -----
    - Uses absolute correlation for fair comparison
    - MSC values are already in [0,1]
    - Agreement computed using Pearson correlation of edge weights
    - Only upper triangle values used (excluding diagonal)
    """
```

---

## Implementation Details

### Step 1: Load Both FC Matrices

```python
from lrg_eegfc.workflow_corr import load_corr_matrix
from lrg_eegfc.workflow_msc import load_msc_matrix

# Load correlation matrix
corr_matrix = load_corr_matrix(
    patient, phase, band,
    cache_root=cache_root_corr,
    use_cleaned=False  # Use raw correlation for fair comparison
)

if corr_matrix is None:
    raise FileNotFoundError(
        f"Correlation matrix not found: {patient} {phase} {band}\n"
        f"Run: python src/compute_corr_matrices.py --patient {patient}"
    )

# Load MSC matrix
msc_matrix = load_msc_matrix(
    patient, phase, band,
    cache_root=cache_root_msc,
    sparsify="none",  # Use dense MSC for fair comparison
    n_surrogates=0,
    nperseg=nperseg
)

if msc_matrix is None:
    raise FileNotFoundError(
        f"MSC matrix not found: {patient} {phase} {band}\n"
        f"Run: python src/compute_msc_matrices.py --patient {patient}"
    )

# Take absolute correlation for comparison
corr_matrix_abs = np.abs(corr_matrix)
```

### Step 2: Extract Upper Triangle Values

```python
# Get upper triangle indices (exclude diagonal)
N = corr_matrix_abs.shape[0]
triu_idx = np.triu_indices(N, k=1)

# Extract values for comparison
corr_values = corr_matrix_abs[triu_idx]
msc_values = msc_matrix[triu_idx]

# Compute agreement (Pearson correlation)
agreement = np.corrcoef(corr_values, msc_values)[0, 1]

# Compute statistics
corr_stats = {
    'mean': float(corr_values.mean()),
    'std': float(corr_values.std()),
    'median': float(np.median(corr_values)),
}

msc_stats = {
    'mean': float(msc_values.mean()),
    'std': float(msc_values.std()),
    'median': float(np.median(msc_values)),
}
```

### Step 3: Create 3-Panel Figure

```python
fig, axes = plt.subplots(1, 3, figsize=figsize)
fig.suptitle(f"{patient} {phase} {band.upper()} - FC Method Comparison",
             fontsize=14, fontweight='bold')
```

### Step 4: Panel 0 - Correlation Matrix

```python
ax0 = axes[0]

# Plot absolute correlation matrix
im0 = ax0.imshow(corr_matrix_abs, cmap='viridis', vmin=0, vmax=1,
                 origin='upper', aspect='auto')

ax0.set_title('Correlation FC\n(Absolute)', fontsize=12, fontweight='bold')
ax0.set_xlabel('Channel', fontsize=10)
ax0.set_ylabel('Channel', fontsize=10)

# Colorbar
cbar0 = plt.colorbar(im0, ax=ax0, fraction=0.046, pad=0.04)
cbar0.set_label('|Correlation|', fontsize=10)

# Add statistics text
stats_text = f"Mean: {corr_stats['mean']:.3f}\nStd: {corr_stats['std']:.3f}"
ax0.text(0.02, 0.98, stats_text, transform=ax0.transAxes,
         fontsize=9, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
```

### Step 5: Panel 1 - MSC Matrix

```python
ax1 = axes[1]

# Plot MSC matrix
im1 = ax1.imshow(msc_matrix, cmap='viridis', vmin=0, vmax=1,
                 origin='upper', aspect='auto')

ax1.set_title('MSC FC\n(Magnitude-Squared Coherence)', fontsize=12, fontweight='bold')
ax1.set_xlabel('Channel', fontsize=10)
ax1.set_ylabel('Channel', fontsize=10)

# Colorbar
cbar1 = plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
cbar1.set_label('MSC', fontsize=10)

# Add statistics text
stats_text = f"Mean: {msc_stats['mean']:.3f}\nStd: {msc_stats['std']:.3f}"
ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes,
         fontsize=9, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
```

### Step 6: Panel 2 - Method Agreement Analysis

This panel combines scatter plot and distribution comparison:

```python
ax2 = axes[2]

# Create twin axis for distributions
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax2)
ax2_top = divider.append_axes("top", size="30%", pad=0.1, sharex=ax2)

# --- Top subplot: Distributions ---
ax2_top.hist(corr_values, bins=40, alpha=0.6, label='Correlation',
             color='steelblue', density=True, edgecolor='none')
ax2_top.hist(msc_values, bins=40, alpha=0.6, label='MSC',
             color='coral', density=True, edgecolor='none')
ax2_top.set_ylabel('Density', fontsize=9)
ax2_top.legend(fontsize=9, loc='upper right')
ax2_top.grid(alpha=0.3)
ax2_top.tick_params(labelbottom=False)

# --- Main subplot: Scatter plot ---
# Subsample for visualization if too many points
n_points = len(corr_values)
if n_points > 5000:
    sample_idx = np.random.choice(n_points, 5000, replace=False)
    corr_plot = corr_values[sample_idx]
    msc_plot = msc_values[sample_idx]
else:
    corr_plot = corr_values
    msc_plot = msc_values

ax2.scatter(corr_plot, msc_plot, alpha=0.3, s=2, color='steelblue')
ax2.plot([0, 1], [0, 1], 'r--', linewidth=1.5, alpha=0.7, label='Identity')

ax2.set_xlabel('Correlation FC', fontsize=10)
ax2.set_ylabel('MSC FC', fontsize=10)
ax2.set_title('Method Agreement', fontsize=12, fontweight='bold')
ax2.grid(alpha=0.3)
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)

# Add agreement metric
agreement_text = (
    f'Pearson r = {agreement:.3f}\n'
    f'N edges = {n_points}'
)
ax2.text(0.05, 0.95, agreement_text, transform=ax2.transAxes,
         fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
```

### Step 7: Save Figure

```python
if output_path is None:
    output_dir = Path("data/figures/comparison") / patient
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{band}_{phase}_fc_comparison.png"

plt.tight_layout()
plt.savefig(output_path, dpi=150, bbox_inches='tight')
plt.close(fig)

return output_path
```

---

## CLI Integration

### Create Script: `src/visualize_comparison.py`

```python
#!/usr/bin/env python3
"""CLI script for FC method comparison visualizations.

Examples
--------
# Single comparison
python src/visualize_comparison.py --patient Pat_02 --phase rsPre --band beta --verbose

# Batch mode: all bands for a phase
python src/visualize_comparison.py --patient Pat_02 --phase rsPre --batch --verbose
"""

import argparse
import sys
from pathlib import Path

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.visuals.compare import plot_fc_comparison


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare Correlation vs MSC functional connectivity methods"
    )

    parser.add_argument("--patient", required=True, help="Patient ID")
    parser.add_argument("--phase", help="Recording phase")
    parser.add_argument("--band", help="Frequency band")

    parser.add_argument(
        "--nperseg",
        type=int,
        default=1024,
        help="MSC window length (must match cache). Default: 1024"
    )

    parser.add_argument(
        "--cache-root-corr",
        type=Path,
        default=Path("data/corr_cache"),
        help="Correlation cache directory"
    )

    parser.add_argument(
        "--cache-root-msc",
        type=Path,
        default=Path("data/msc_cache"),
        help="MSC cache directory"
    )

    parser.add_argument(
        "--batch",
        action="store_true",
        help="Generate for all bands (requires --phase)"
    )

    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    # Validation
    if args.batch and not args.phase:
        parser.error("--batch requires --phase")

    if not args.batch and (not args.phase or not args.band):
        parser.error("--phase and --band required unless using --batch")

    # Determine combinations
    if args.batch:
        bands = list(BRAIN_BANDS.keys())
        phases = [args.phase]
    else:
        bands = [args.band]
        phases = [args.phase]

    # Process
    for band in bands:
        for phase in phases:
            if args.verbose:
                print(f"Comparing FC methods: {args.patient} {phase} {band}")

            try:
                output_path = plot_fc_comparison(
                    patient=args.patient,
                    phase=phase,
                    band=band,
                    cache_root_corr=args.cache_root_corr,
                    cache_root_msc=args.cache_root_msc,
                    nperseg=args.nperseg,
                )

                if args.verbose:
                    print(f"  ✓ Saved: {output_path}")

            except FileNotFoundError as e:
                print(f"  ERROR: {e}")
            except Exception as e:
                print(f"  ERROR: Failed for {args.patient} {phase} {band}: {e}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
```

---

## Testing Instructions

### Test 1: Single Comparison

```bash
# Ensure both cache files exist
python src/compute_corr_matrices.py --patient Pat_02 --verbose
python src/compute_msc_matrices.py --patient Pat_02 --verbose

# Generate comparison
python src/visualize_comparison.py --patient Pat_02 --phase rsPre --band beta --verbose
```

Expected output: `data/figures/comparison/Pat_02/beta_rsPre_fc_comparison.png`

### Test 2: Batch Mode

```bash
# All bands for one phase
python src/visualize_comparison.py --patient Pat_02 --phase rsPre --batch --verbose
```

Expected: 6 comparison figures (one per band)

---

## Expected Output Format

### Figure Layout

```
┌──────────────────┬──────────────────┬──────────────────┐
│  Correlation FC  │     MSC FC       │ Method Agreement │
│  (Absolute)      │  (Coherence)     │                  │
│                  │                  │ ┌──────────────┐ │
│  [Heatmap]       │  [Heatmap]       │ │Distributions │ │
│  viridis         │  viridis         │ └──────────────┘ │
│  [0, 1]          │  [0, 1]          │                  │
│                  │                  │  [Scatter Plot]  │
│  Stats box       │  Stats box       │  + Identity line │
│                  │                  │  Agreement box   │
└──────────────────┴──────────────────┴──────────────────┘

Title: Pat_02 rsPre BETA - FC Method Comparison
```

### Panel Descriptions

**Panel 0 (Correlation FC):**
- Heatmap of absolute correlation matrix
- Colorbar: |Correlation| [0, 1]
- Stats box: Mean, Std

**Panel 1 (MSC FC):**
- Heatmap of MSC matrix
- Colorbar: MSC [0, 1]
- Stats box: Mean, Std

**Panel 2 (Method Agreement):**
- Top subplot: Overlaid histograms (Correlation vs MSC distributions)
- Main subplot: Scatter plot (Correlation vs MSC edge weights)
- Red dashed identity line
- Agreement box: Pearson r, N edges

---

## Integration with Pipeline

### Update `scripts/run_full_analysis.sh`

Add after MSC visualization (Step 5b):

```bash
# -------------------------------------------------------------------------
# Step 5c: Generate FC method comparison visualizations
# -------------------------------------------------------------------------
echo -e "${BLUE}[5c/8]${NC} Generating FC comparison visualizations..."
if [ -f "src/visualize_comparison.py" ]; then
  if python src/visualize_comparison.py --patient "$PATIENT" --batch --verbose 2>/dev/null; then
    echo -e "${GREEN}✓${NC} FC comparison visualizations generated"
  else
    echo -e "${YELLOW}⚠${NC} FC comparison skipped or failed for $PATIENT"
  fi
else
  echo -e "${YELLOW}⚠${NC} Comparison visualization script not available yet"
fi
```

---

## Additional Requirements

### Export Function

Add to `src/lrg_eegfc/__init__.py`:

```python
from .visuals.compare import plot_fc_comparison
```

### Error Handling

The function MUST check:
1. Both cache files exist
2. Matrices have same shape
3. Matrices are symmetric
4. Valid band name

### Performance Considerations

- Subsample scatter plot if > 5000 points (already included above)
- Use `density=True` for histograms (fair comparison)
- Use `edgecolor='none'` for histogram performance

---

## Validation Checklist

Before marking as complete:

- [ ] Function created in `src/lrg_eegfc/visuals/compare.py`
- [ ] CLI script created: `src/visualize_comparison.py`
- [ ] Exports added to `__init__.py`
- [ ] Test with Pat_02 rsPre beta (single)
- [ ] Test with Pat_02 rsPre --batch (all bands)
- [ ] Verify output in `data/figures/comparison/Pat_02/`
- [ ] Verify figure has 3 panels as specified
- [ ] Verify agreement metric is displayed
- [ ] Verify distributions overlay correctly
- [ ] Integration with pipeline script (optional)

---

## Key Design Decisions

1. **Use absolute correlation**: Fair comparison since MSC is always positive
2. **Dense matrices only**: Compare raw methods, not sparsification strategies
3. **Upper triangle only**: Avoid double-counting symmetric edges
4. **Subsample scatter**: Performance optimization for large networks
5. **Combined panel**: Distributions + scatter in one panel saves space
6. **Horizontal layout**: 3 panels side-by-side (not 2x2 grid)

---

## Notes for Implementation

- Follow existing code style in `src/lrg_eegfc/visuals/corr.py` and `msc.py`
- Use same color scheme as other visualizations (viridis for heatmaps)
- Keep consistent with correlation/MSC summary layouts
- Function should be ~200-250 lines including docstring
- CLI should follow same pattern as `visualize_correlation.py` and `visualize_msc.py`

---

## Example Expected Statistics

For Pat_02 rsPre beta:

```
Correlation FC - Mean: ~0.15, Std: ~0.12
MSC FC - Mean: ~0.06, Std: ~0.08
Agreement (Pearson r): ~0.45-0.65 (typical range)
```

These values will vary by patient, phase, and band.

---

**End of Instructions**

When implementing this, you should be able to run:
```bash
python src/visualize_comparison.py --patient Pat_02 --phase rsPre --band beta --verbose
```

And get a publication-quality 3-panel comparison figure.
