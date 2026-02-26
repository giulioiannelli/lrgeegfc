# CLI Reference: lrg-eegfc

The `lrg-eegfc` command-line tool provides unified access to the full analysis
pipeline: computing functional connectivity matrices, running LRG analysis,
generating figures, and managing cache files.

## Installation

```bash
pip install -e .          # installs the lrg-eegfc entry point
lrg-eegfc --help          # verify
```

## Architecture

The CLI uses **Click** with a `LazyGroup` pattern so that `lrg-eegfc --help`
responds instantly without importing numpy/scipy/matplotlib. Each subcommand
module is imported only when invoked.

```
src/lrg_eegfc/cli/
├── __init__.py      # exports `app`
├── _app.py          # LazyGroup root + 7 lazy subcommands
├── _common.py       # shared option factories, resolvers, CliReporter
├── compute.py       # 6 compute commands
├── plot.py          # 16 plot commands (4 are TODO stubs)
├── show.py          # 4 show commands (query cached results)
├── data.py          # 3 data inspection commands
├── cache.py         # 4 cache management commands
├── config_cmd.py    # 2 config display commands
├── bundle.py        # 1 publication bundle command
└── _legacy.py       # old cli.py (kept for lrg-eegfc-corr compat)
```

---

## Command Hierarchy

```
lrg-eegfc [--verbose/-v] [--quiet/-q]
│
├── compute                         # FC + analysis computation (cached)
│   ├── corr                        # Correlation matrices
│   ├── msc                         # MSC matrices (6 sparsification methods)
│   ├── lrg                         # LRG ultrametric analysis
│   ├── clean                       # Marchenko-Pastur spectral cleaning
│   ├── time-windows                # Sliding-window FC
│   └── reorganization              # Phase distance metrics
│
├── plot                            # Visualization (reads from cache)
│   ├── corr                        # Correlation heatmap/network/percolation
│   ├── msc                         # MSC heatmap/network/comparison
│   ├── lrg                         # Entropy, dendrogram, ultrametric, 5-panel
│   ├── reorganization              # Cross-phase comparison panels
│   ├── metastable                  # Sankey cluster evolution (HTML)
│   ├── cleaning                    # MP eigenvalue diagnostics
│   ├── comparison                  # Corr vs MSC side-by-side
│   ├── msc-grid                    # Band x phase MSC grid
│   ├── msc-all-patients            # Cross-patient MSC overview
│   ├── lrg-phase-grid              # LRG across phases
│   ├── lrg-video                   # (TODO) Ultrametric threshold animation
│   ├── msc-validation              # Dense vs validated MSC
│   ├── time-windows                # (TODO) Time-window FC panels
│   ├── reorg-metrics               # Reorganization distance matrices
│   ├── reorg-summary               # (TODO) Band-level metric summary
│   └── metric-correlation          # (TODO) Cross-metric agreement
│
├── show                            # Query cached results (no figures)
│   ├── corr                        # Correlation matrix stats
│   ├── msc                         # MSC matrix stats
│   ├── lrg                         # LRG analysis summary
│   └── cleaned                     # MP cleaning metadata
│
├── data                            # Inspection & statistics
│   ├── inspect                     # Patient data report
│   ├── stats                       # Aggregate cross-patient stats
│   └── compare                     # FC method comparison (tabular)
│
├── cache                           # Cache management
│   ├── list                        # List cached files with sizes
│   ├── verify                      # Check cache integrity
│   ├── clean                       # Remove stale/dev caches
│   └── info                        # Print naming conventions
│
├── config                          # Configuration
│   ├── show                        # Display bands, phases, defaults
│   └── paths                       # Show/verify data paths
│
└── bundle                          # Publication utilities
    └── overleaf                    # Collect figures into Overleaf bundle
```

---

## Global Options

| Flag | Description |
|------|-------------|
| `-v`, `--verbose` | Enable verbose output (works at root or command level) |
| `-q`, `--quiet` | Suppress non-error output |
| `--version` | Show package version |

---

## Common Option Patterns

Most commands share these option decorators (defined in `_common.py`):

| Option | Description | Default |
|--------|-------------|---------|
| `--patients` / `-p` | Patient IDs (multiple). Omit for all. | auto-discover |
| `--patient` | Single patient ID (plot commands) | required |
| `--bands` | Band names (multiple). Omit for all. | all bands |
| `--band` | Single band name | - |
| `--phases` | Phase names (multiple). Omit for all. | all phases |
| `--phase` | Single phase name | - |
| `--fc-method` | `corr` or `msc` | `msc` |
| `--cache-root` | Cache directory | varies by command |
| `--overwrite` | Recompute even if cached | false |
| `--output-dir` | Figure output directory | varies by command |
| `--dpi` | Figure resolution | 150 |
| `--filter-time` | Dev mode: limit to N samples | disabled |

---

## Compute Commands

### `lrg-eegfc compute corr`

Compute correlation-based FC matrices for all band/phase combinations.

```bash
# All patients, all bands, all phases
lrg-eegfc compute corr -v

# Single patient, specific band
lrg-eegfc compute corr --patients Pat_02 --band alpha --phase rsPre -v

# Custom filter type
lrg-eegfc compute corr --patients Pat_02 --filter-type pos --zero-diagonal -v
```

Delegates to: `workflow.corr.compute_corr_for_patient()`

### `lrg-eegfc compute msc`

Compute magnitude-squared coherence FC matrices with optional sparsification.

```bash
# Dense MSC (no sparsification)
lrg-eegfc compute msc --patients Pat_02 --band beta --phase rsPre -v

# Surrogate-based soft thresholding (auto-sets n_surrogates=200)
lrg-eegfc compute msc --patients Pat_02 --sparsify soft -v

# FDR sparsification with custom q-value
lrg-eegfc compute msc --patients Pat_02 --sparsify fdr --fdr-q 0.01 -v

# ECM (Exponential Configuration Model)
lrg-eegfc compute msc --patients Pat_02 --sparsify ecm --ecm-alpha 0.05 -v

# Dev mode (first 10000 samples only, uses _dev cache)
lrg-eegfc compute msc --patients Pat_02 --filter-time 10000 -v
```

Sparsification methods: `none`, `soft`, `fdr`, `disparity`, `hybrid`, `ecm`

Delegates to: `workflow.msc.compute_msc_for_patient()`

### `lrg-eegfc compute lrg`

Run LRG ultrametric analysis on cached FC matrices.

```bash
# Using MSC FC (default)
lrg-eegfc compute lrg --patients Pat_02 --band alpha --phase rsPre -v

# Using correlation FC
lrg-eegfc compute lrg --patients Pat_02 --fc-method corr -v

# Custom entropy range
lrg-eegfc compute lrg --patients Pat_02 --entropy-steps 600 --entropy-t1 -4 --entropy-t2 6 -v
```

Delegates to: `workflow.lrg.compute_lrg_for_patient()`

### `lrg-eegfc compute clean`

Apply Marchenko-Pastur spectral cleaning to correlation matrices.

```bash
lrg-eegfc compute clean --patients Pat_02 --band alpha --phase rsPre -v
```

Delegates to: `workflow.cleaning.compute_cleaned_corr_for_patient()`

### `lrg-eegfc compute time-windows`

Compute sliding-window FC matrices.

```bash
lrg-eegfc compute time-windows --patients Pat_02 --band alpha --phase rsPre \
    --window-sec 10 --overlap 0.5 -v
```

### `lrg-eegfc compute reorganization`

Compute reorganization metrics from cached LRG results.

```bash
lrg-eegfc compute reorganization --patients Pat_02 --fc-method msc -v
```

---

## Plot Commands

All plot commands read from cache and never recompute. Existing figures are
skipped unless `--overwrite` is used (via cache) or you delete the output.

### `lrg-eegfc plot corr`

```bash
# Summary plot (combined figure)
lrg-eegfc plot corr --patient Pat_02 --band alpha --phase rsPre -v

# Network visualization (cleaned matrices)
lrg-eegfc plot corr --patient Pat_02 --plot-type network --cleaned -v

# All plot types at once
lrg-eegfc plot corr --patient Pat_02 --plot-type all -v
```

Plot types: `summary`, `network`, `mp`, `percolation`, `all`

### `lrg-eegfc plot msc`

```bash
# MSC summary
lrg-eegfc plot msc --patient Pat_02 --band alpha --phase rsPre -v

# Network with soft sparsification
lrg-eegfc plot msc --patient Pat_02 --plot-type network --sparsify soft -v
```

Plot types: `network`, `comparison`, `summary`, `all`

### `lrg-eegfc plot lrg`

```bash
# Full 5-panel figure (entropy + dendrogram + ultrametric + FC + network)
lrg-eegfc plot lrg --patient Pat_02 --band alpha --phase rsPre --fc-method msc -v

# Just entropy curves
lrg-eegfc plot lrg --patient Pat_02 --band alpha --phase rsPre \
    --fc-method msc --plot-type entropy -v

# All LRG plot types
lrg-eegfc plot lrg --patient Pat_02 --fc-method msc --plot-type all -v
```

Plot types: `entropy`, `dendrogram`, `ultrametric`, `full`, `all`

### `lrg-eegfc plot reorganization`

```bash
lrg-eegfc plot reorganization --patient Pat_02 --fc-method msc -v
```

### `lrg-eegfc plot metastable`

Generates interactive Sankey diagrams (HTML output).

```bash
lrg-eegfc plot metastable --patient Pat_02 --band alpha --phase rsPre --fc-method msc -v
```

### `lrg-eegfc plot comparison`

Side-by-side correlation vs MSC comparison.

```bash
lrg-eegfc plot comparison --patient Pat_02 --band alpha --phase rsPre -v
```

### `lrg-eegfc plot msc-grid`

Band x phase grid of MSC matrices for one patient.

```bash
lrg-eegfc plot msc-grid --patient Pat_02 -v
```

### `lrg-eegfc plot msc-all-patients`

MSC matrices across all patients for given band/phase.

```bash
lrg-eegfc plot msc-all-patients --band alpha --phase rsPre -v
```

### `lrg-eegfc plot lrg-phase-grid`

LRG full panels across all phases.

```bash
lrg-eegfc plot lrg-phase-grid --patient Pat_02 --fc-method msc -v
```

### `lrg-eegfc plot msc-validation`

Dense vs surrogate-validated MSC comparison.

```bash
lrg-eegfc plot msc-validation --patient Pat_02 --band alpha --phase rsPre -v
```

### `lrg-eegfc plot reorg-metrics`

Reorganization distance matrices.

```bash
lrg-eegfc plot reorg-metrics --patient Pat_02 --fc-method msc -v
```

---

## Show Commands

Query cached results and print summary statistics to the terminal.
No figures are generated. Omit `--band` or `--phase` to loop over all values.

### `lrg-eegfc show corr`

```bash
# Single band/phase
lrg-eegfc show corr --patient Pat_02 --band alpha --phase rsPre

# All bands for one phase
lrg-eegfc show corr --patient Pat_02 --phase rsPre
```

Prints: shape, mean, std, min, max, median, density, edge count.

### `lrg-eegfc show msc`

```bash
lrg-eegfc show msc --patient Pat_02 --band alpha --phase rsPre
lrg-eegfc show msc --patient Pat_02 --band alpha --phase rsPre --sparsify soft --n-surrogates 200
```

Prints: same matrix stats as `show corr`, plus sparsification method.

### `lrg-eegfc show lrg`

```bash
# One band
lrg-eegfc show lrg --patient Pat_02 --band alpha --phase rsPre --fc-method msc

# All bands at once
lrg-eegfc show lrg --patient Pat_02 --phase rsPre --fc-method msc
```

Prints: n_nodes, optimal_threshold, entropy tau range, entropy steps,
1-S range, C range, C peak (value and tau location), ultrametric distance
stats (mean, std, range), linkage merges, max merge distance.

### `lrg-eegfc show cleaned`

```bash
lrg-eegfc show cleaned --patient Pat_02 --band alpha --phase rsPre
```

Prints: matrix stats + MP cleaning metadata (threshold, lambda_min/max,
signal component count, signal eigenvalue range).

---

## Data Commands

### `lrg-eegfc data inspect`

```bash
# Single patient
lrg-eegfc data inspect --patients Pat_02 -v

# All patients
lrg-eegfc data inspect -v
```

### `lrg-eegfc data stats`

```bash
lrg-eegfc data stats -v
```

### `lrg-eegfc data compare`

Compare MSC vs correlation FC methods (tabular output).

```bash
lrg-eegfc data compare --patients Pat_02 --band alpha --phase rsPre -v
```

---

## Cache Commands

### `lrg-eegfc cache list`

```bash
# All cache types
lrg-eegfc cache list

# Only MSC cache for a specific patient
lrg-eegfc cache list -t msc --patients Pat_02
```

### `lrg-eegfc cache verify`

```bash
lrg-eegfc cache verify
lrg-eegfc cache verify -t corr --patients Pat_02
```

### `lrg-eegfc cache info`

Print cache file naming conventions.

```bash
lrg-eegfc cache info
```

### `lrg-eegfc cache clean`

```bash
# Remove only dev caches
lrg-eegfc cache clean --dev-only

# Dry run (show what would be removed)
lrg-eegfc cache clean --dry-run
```

---

## Config Commands

### `lrg-eegfc config show`

Display frequency bands, phase labels, and default parameters.

```bash
lrg-eegfc config show
```

### `lrg-eegfc config paths`

Show and verify data/cache directory paths.

```bash
lrg-eegfc config paths
```

---

## Bundle Commands

### `lrg-eegfc bundle overleaf`

Collect figures into an Overleaf-ready bundle.

```bash
lrg-eegfc bundle overleaf --dry-run
lrg-eegfc bundle overleaf --patient Pat_02
```

---

## Typical Workflows

### Full pipeline for one patient

```bash
# 1. Inspect data
lrg-eegfc data inspect --patients Pat_02 -v

# 2. Compute FC matrices
lrg-eegfc compute corr --patients Pat_02 -v
lrg-eegfc compute msc --patients Pat_02 --sparsify none -v

# 3. Run LRG analysis
lrg-eegfc compute lrg --patients Pat_02 --fc-method msc -v

# 4. Check results
lrg-eegfc show lrg --patient Pat_02 --phase rsPre --fc-method msc
lrg-eegfc show msc --patient Pat_02 --phase rsPre

# 5. Generate figures
lrg-eegfc plot lrg --patient Pat_02 --fc-method msc --plot-type full -v
lrg-eegfc plot corr --patient Pat_02 --plot-type all -v
lrg-eegfc plot msc --patient Pat_02 --plot-type all -v

# 6. Phase reorganization
lrg-eegfc plot reorganization --patient Pat_02 --fc-method msc -v
```

### Batch processing all patients

```bash
# Compute everything (omitting --patients processes all)
lrg-eegfc compute corr -v
lrg-eegfc compute msc -v
lrg-eegfc compute lrg --fc-method msc -v

# Verify caches
lrg-eegfc cache list
lrg-eegfc cache verify
```

### Quick dev iteration (limited data)

```bash
# Use --filter-time to limit samples (auto-switches to _dev cache)
lrg-eegfc compute msc --patients Pat_02 --band alpha --phase rsPre \
    --filter-time 10000 -v
```

---

## Key Implementation Notes

- **Verbose flag**: `-v` works at both root and command level.
  `lrg-eegfc -v compute corr` and `lrg-eegfc compute corr -v` are equivalent.
- **Auto-discovery**: Omitting `--patients`, `--bands`, or `--phases` defaults
  to processing all available values.
- **Skip existing**: Plot commands skip existing output files. Delete the file
  to regenerate.
- **Dev mode**: `--filter-time N` limits input to N samples and automatically
  switches to `*_dev` cache roots to avoid polluting production caches.
- **MSC surrogates**: Methods requiring surrogates (`soft`, `fdr`, `hybrid`)
  auto-set `n_surrogates=200` if not specified.

## Source Files

| File | Purpose |
|------|---------|
| `src/lrg_eegfc/cli/__init__.py` | Package entry, exports `app` |
| `src/lrg_eegfc/cli/_app.py` | LazyGroup root, 7 lazy subcommands |
| `src/lrg_eegfc/cli/_common.py` | Shared options, resolvers, CliReporter |
| `src/lrg_eegfc/cli/compute.py` | 6 compute commands |
| `src/lrg_eegfc/cli/plot.py` | 16 plot commands |
| `src/lrg_eegfc/cli/show.py` | 4 show commands (query cached results) |
| `src/lrg_eegfc/cli/data.py` | 3 data commands |
| `src/lrg_eegfc/cli/cache.py` | 4 cache commands |
| `src/lrg_eegfc/cli/config_cmd.py` | 2 config commands |
| `src/lrg_eegfc/cli/bundle.py` | 1 bundle command |
| `src/lrg_eegfc/cli/_legacy.py` | Old cli.py (lrg-eegfc-corr compat) |
