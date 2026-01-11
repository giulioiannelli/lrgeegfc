# Complete Analysis Command List

Run these commands in sequence to perform the full analysis pipeline for all patients (Pat_02, Pat_03, Pat_05, Pat_07, Pat_08).

**Important Notes:**
- Pat_05, Pat_07, Pat_08 will use default 2048 Hz (fs missing in data)
- Missing phases will be automatically skipped with warnings
- Each step builds on the previous (FC → LRG → Comparison)

---

## Step 1: Data Inspection (Already Done)

```bash
# Already completed - review reports:
cat PATIENT_DATA_REPORT.txt
cat DATA_STATUS.md
```

**Expected:** Report shows Pat_02 and Pat_03 complete, Pat_05/07/08 missing fs

---

## Step 2: Compute FC Matrices

### 2A. Correlation-Based FC (All Patients)

```bash
# Compute correlation matrices (absolute values, zero diagonal)
python src/compute_corr_matrices.py --verbose
```

**Expected Output:**
- Pat_02: 24 matrices (6 bands × 4 phases)
- Pat_03: 24 matrices (6 bands × 4 phases)
- Pat_05: 24 matrices (6 bands × 4 phases)
- Pat_07: ~18 matrices (taskTest missing/corrupted)
- Pat_08: 24 matrices (6 bands × 4 phases)

**Cache Location:** `data/corr_cache/{patient}/{band}_{phase}_corr_ftype-abs_zdiag-True.npy`

**Estimated Time:** ~5-10 minutes

---

### 2B. MSC-Based FC (All Patients, Dense)

```bash
# Compute dense MSC matrices (no validation yet)
python src/compute_msc_matrices.py --verbose
```

**Expected Output:**
- Same matrix counts as correlation
- Pat_07 taskTest will fail/skip

**Cache Location:** `data/msc_cache/{patient}/{band}_{phase}_msc_sparsify-none_nperseg-256.npy`

**Estimated Time:** ~10-15 minutes (MSC is slower than correlation)

---

### 2C. Clean Correlation Matrices (Optional but Recommended)

```bash
# Clean correlation matrices using Marchenko-Pastur + percolation thresholding
python src/clean_correlation_matrices.py --verbose
```

**What this does:**
1. **Marchenko-Pastur spectral filtering** - Removes noise eigenvalues from correlation matrices
2. **Percolation thresholding** - Thresholds network just before first node detachment

**Expected Output:**
- Same matrix counts as Step 2A
- Cleaned matrices with `_cleaned.npy` suffix
- Metadata files with eigenvalue and threshold information

**Cache Location:**
- Matrices: `data/corr_cache/{patient}/{band}_{phase}_corr_cleaned.npy`
- Metadata: `data/corr_cache/{patient}/{band}_{phase}_corr_cleaned_meta.npz`

**Estimated Time:** ~10-15 minutes (similar to regular correlation)

**Why clean?**
- Removes spurious correlations due to finite-size effects
- Preserves network connectivity (avoids over-thresholding)
- More robust for downstream LRG analysis
- MSC matrices don't need cleaning (already frequency-specific)

---

## Step 3: Compute LRG Analysis

### 3A. LRG from Correlation Matrices

```bash
# Compute ultrametric distances from correlation FC
python src/compute_lrg_analysis.py --fc-method corr --verbose
```

**Expected Output:**
- ~114 LRG analyses (same as correlation matrices that succeeded)
- Each contains: ultrametric matrix, linkage, entropy curves

**Cache Location:** `data/lrg_cache/{patient}/{band}_{phase}_lrg_corr.npz`

**Estimated Time:** ~15-20 minutes

---

### 3B. LRG from MSC Matrices

```bash
# Compute ultrametric distances from MSC FC
python src/compute_lrg_analysis.py --fc-method msc --verbose
```

**Expected Output:**
- ~114 LRG analyses (same as MSC matrices that succeeded)

**Cache Location:** `data/lrg_cache/{patient}/{band}_{phase}_lrg_msc.npz`

**Estimated Time:** ~15-20 minutes

---

## Step 4: Compare Networks

### 4A. MSC vs Correlation Comparison

```bash
# Create results directory
mkdir -p results

# Compare MSC vs Correlation across all patients/phases/bands
python src/compare_fc_methods.py --mode fc-methods --output results/msc_vs_corr.csv --verbose
```

**Expected Output:**
- CSV with ~114 rows (one per successful patient/phase/band)
- Columns: patient, phase, band, matrix_distance, multiscale_distance, etc.

**Output:** `results/msc_vs_corr.csv`

**Estimated Time:** <1 minute (reads from cache)

---

### 4B. Memory Effects Analysis (Correlation)

```bash
# Analyze pre vs post changes using correlation FC
python src/compare_fc_methods.py --mode memory --fc-method corr --output results/memory_corr.csv --verbose
```

**Expected Output:**
- Comparisons for: rsPre vs rsPost, taskLearn vs taskTest
- Per patient/band
- Shows network changes due to task/memory

**Output:** `results/memory_corr.csv`

**Estimated Time:** <1 minute

---

### 4C. Memory Effects Analysis (MSC)

```bash
# Analyze pre vs post changes using MSC FC
python src/compare_fc_methods.py --mode memory --fc-method msc --output results/memory_msc.csv --verbose
```

**Expected Output:**
- Same as 4B but for MSC-based networks

**Output:** `results/memory_msc.csv`

**Estimated Time:** <1 minute

---

## Complete Pipeline (All Steps)

Run everything in one go:

```bash
# Step 2: FC matrices
echo "=== Step 2A: Correlation FC ==="
python src/compute_corr_matrices.py --verbose

echo "=== Step 2B: MSC FC ==="
python src/compute_msc_matrices.py --verbose

# Step 3: LRG analysis
echo "=== Step 3A: LRG from Correlation ==="
python src/compute_lrg_analysis.py --fc-method corr --verbose

echo "=== Step 3B: LRG from MSC ==="
python src/compute_lrg_analysis.py --fc-method msc --verbose

# Step 4: Comparisons
echo "=== Step 4A: MSC vs Correlation ==="
mkdir -p results
python src/compare_fc_methods.py --mode fc-methods --output results/msc_vs_corr.csv --verbose

echo "=== Step 4B: Memory Effects (Correlation) ==="
python src/compare_fc_methods.py --mode memory --fc-method corr --output results/memory_corr.csv --verbose

echo "=== Step 4C: Memory Effects (MSC) ==="
python src/compare_fc_methods.py --mode memory --fc-method msc --output results/memory_msc.csv --verbose

echo "=== PIPELINE COMPLETE ==="
```

---

## Expected Warnings (Normal)

These warnings are expected and can be ignored:

```
WARNING: No sampling rate found in data, using provided value: 2048 Hz
  → Pat_05, Pat_07, Pat_08 (missing fs parameter)

WARNING: No cached FC matrix for Pat_07 taskTest msc
  → taskTest.mat corrupted for Pat_07

WARNING: Failed to compute Pat_07 taskTest beta: No data variable found
  → taskTest missing for Pat_07
```

---

## Verify Results

After running, check:

```bash
# Count cached matrices
echo "Correlation matrices:"
find data/corr_cache -name "*.npy" | wc -l

echo "MSC matrices:"
find data/msc_cache -name "*.npy" | wc -l

echo "LRG analyses (corr):"
find data/lrg_cache -name "*_lrg_corr.npz" | wc -l

echo "LRG analyses (msc):"
find data/lrg_cache -name "*_lrg_msc.npz" | wc -l

# Check comparison results
echo "MSC vs Corr comparisons:"
wc -l results/msc_vs_corr.csv

echo "Memory effects (corr):"
wc -l results/memory_corr.csv

echo "Memory effects (msc):"
wc -l results/memory_msc.csv
```

**Expected Counts:**
- Correlation matrices: ~114 (depends on missing phases)
- MSC matrices: ~114
- LRG (corr): ~114
- LRG (msc): ~114
- Comparisons: ~114 rows in CSV

---

## Custom Analysis Examples

### Analyze Specific Band

```bash
# Only beta band
python src/compute_corr_matrices.py --patients Pat_02 --verbose
python src/compute_lrg_analysis.py --fc-method corr --patients Pat_02 --verbose
```

### Different Correlation Filter

```bash
# Positive correlations only
python src/compute_corr_matrices.py --filter-type pos --verbose

# Then run LRG on these
python src/compute_lrg_analysis.py --fc-method corr --verbose
```

### Validated MSC with Surrogates

```bash
# MSC with surrogate-based validation (much slower!)
python src/compute_msc_matrices.py --sparsify soft --n-surrogates 200 --patients Pat_02 --verbose
```

---

## File Sizes (Approximate)

- Correlation matrix (.npy): ~100 KB each
- MSC matrix (.npy): ~100 KB each
- LRG analysis (.npz): ~50 KB each
- Comparison CSV: ~50 KB total

**Total disk usage:** ~50-100 MB for all patients

---

## Troubleshooting

### Out of Memory Errors

If you get memory errors, process patients one at a time:

```bash
for patient in Pat_02 Pat_03 Pat_05 Pat_07 Pat_08; do
    echo "Processing $patient..."
    python src/compute_corr_matrices.py --patients $patient --verbose
    python src/compute_msc_matrices.py --patients $patient --verbose
    python src/compute_lrg_analysis.py --fc-method corr --patients $patient --verbose
    python src/compute_lrg_analysis.py --fc-method msc --patients $patient --verbose
done
```

### Long Runtime

If computation is too slow, reduce to 2 patients for testing:

```bash
# Quick test with Pat_02 and Pat_03 only
python src/compute_corr_matrices.py --patients Pat_02 Pat_03 --verbose
python src/compute_msc_matrices.py --patients Pat_02 Pat_03 --verbose
python src/compute_lrg_analysis.py --fc-method corr --patients Pat_02 Pat_03 --verbose
python src/compute_lrg_analysis.py --fc-method msc --patients Pat_02 Pat_03 --verbose
python src/compare_fc_methods.py --patients Pat_02 Pat_03 --output results/test.csv
```

---

## Next Steps After Analysis

1. **Examine comparison results:**
   ```bash
   python -c "import pandas as pd; df = pd.read_csv('results/msc_vs_corr.csv'); print(df.groupby('band')['matrix_distance'].describe())"
   ```

2. **Statistical analysis:**
   - Test if MSC vs Corr differences are significant
   - ANOVA across bands/patients
   - Post-hoc tests

3. **Visualization:**
   - Plot entropy curves
   - Dendrograms
   - Heatmaps of distance matrices

4. **Memory effects:**
   - Analyze pre-post changes
   - Correlation with behavioral metrics

5. **Validation pipeline:**
   - Implement full surrogate-based sparsification
   - Compare validated vs dense MSC
