---
name: data-status-report
type: plan
era: MSC
status: superseded
created: 2025-12-09
updated: 2026-04-24
pointers: []
---

# Patient Data Status Report

Generated: 2025-12-09

## Summary

- **Total patients**: 6 (Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08)
- **Fully complete**: 2 (Pat_02, Pat_03)
- **Patients with issues**: 4 (Pat_05, Pat_06, Pat_07, Pat_08)
- **Total issues**: 19

## Patient-by-Patient Status

### Pat_02 ✓ COMPLETE
- **Phases**: rsPre, taskLearn, taskTest, rsPost (all present)
- **Sampling rate**: 2048 Hz (all phases)
- **Channel labels**: ✓ channel_labels.csv, ChannelNames.mat
- **Issues**: None

### Pat_03 ✓ COMPLETE
- **Phases**: rsPre, taskLearn, taskTest, rsPost (all present)
- **Sampling rate**: 1024 Hz (all phases) - **Note: Different from Pat_02!**
- **Channel labels**: ✓ channel_labels.csv
- **Issues**: None (data shapes auto-corrected by transpose detection)

### Pat_05 ⚠ INCOMPLETE
- **Phases**: rsPre, taskLearn, taskTest, rsPost (all present)
- **Sampling rate**: ✗ NOT FOUND in any phase
- **Channel labels**: ✓ channel_labels.csv
- **Issues**:
  - Missing sampling rate (fs) in all 4 phases
  - **Action needed**: Add 'fs' field to Parameters struct or as standalone variable

### Pat_06 ✗ INCOMPLETE
- **Phases**: rsPre ✓, rsPost ✓, taskLearn ✗, taskTest ✗
- **Sampling rate**: ✗ NOT FOUND in rsPre, rsPost
- **Channel labels**: ✗ NO channel label files
- **Issues**:
  - Missing taskLearn.mat
  - Missing taskTest.mat
  - Missing sampling rate in 2 phases
  - No channel label files
  - **Action needed**: Provide missing .mat files, add fs, provide channel_labels.csv or ChannelNames.mat

### Pat_07 ✗ INCOMPLETE
- **Phases**: rsPre ✓, taskLearn ✓, rsPost ✓, taskTest ✗ (corrupted)
- **Sampling rate**: ✗ NOT FOUND in rsPre, taskLearn, rsPost
- **Channel labels**: ✗ NO channel label files
- **Issues**:
  - taskTest.mat exists but contains no data variable (corrupted/empty)
  - Missing sampling rate in 3 phases
  - No channel label files
  - **Action needed**: Fix taskTest.mat, add fs, provide channel_labels.csv or ChannelNames.mat

### Pat_08 ⚠ INCOMPLETE
- **Phases**: rsPre, taskLearn, taskTest, rsPost (all present)
- **Sampling rate**: ✗ NOT FOUND in any phase
- **Channel labels**: ✗ NO channel label files
- **Issues**:
  - Missing sampling rate in all 4 phases
  - No channel label files
  - **Action needed**: Add fs, provide channel_labels.csv or ChannelNames.mat

## Critical Divergences

### Sampling Rates
- **Pat_02**: 2048 Hz
- **Pat_03**: 1024 Hz ⚠ **DIFFERENT!**
- **Pat_05-08**: Unknown (missing)

This divergence is handled automatically by the pipeline, which now extracts sampling rates from the data files.

### Data Variable Names
All files use 'Data' as the variable name (consistent).

### File Formats
All files are MATLAB v5-v7 format (scipy.io.loadmat compatible).

## Actions Required from Collaborators

### High Priority (Blocking Analysis)

1. **Pat_06**:
   - Provide taskLearn.mat
   - Provide taskTest.mat
   - Add 'fs' parameter to rsPre.mat and rsPost.mat
   - Provide channel_labels.csv or ChannelNames.mat

2. **Pat_07**:
   - Fix taskTest.mat (currently corrupted/empty)
   - Add 'fs' parameter to rsPre.mat, taskLearn.mat, rsPost.mat
   - Provide channel_labels.csv or ChannelNames.mat

### Medium Priority (Analysis Can Proceed with Warnings)

3. **Pat_05**:
   - Add 'fs' parameter to all .mat files (rsPre, taskLearn, taskTest, rsPost)
   - Currently defaults to 2048 Hz with warning

4. **Pat_08**:
   - Add 'fs' parameter to all .mat files (rsPre, taskLearn, taskTest, rsPost)
   - Provide channel_labels.csv or ChannelNames.mat
   - Currently defaults to 2048 Hz with warning

## How to Add Missing Sampling Rate (fs)

The pipeline expects 'fs' (sampling rate in Hz) in one of these formats:

**Option 1: In Parameters struct**
```matlab
Parameters.fs = 2048;  % or 1024, or actual sampling rate
```

**Option 2: As standalone variable**
```matlab
fs = 2048;  % or 1024, or actual sampling rate
save('phase.mat', 'Data', 'fs');
```

## Current Pipeline Capabilities

The analysis pipeline now includes:

1. **Robust data loading**: Automatically handles different .mat formats (v5-v7, v7.3)
2. **Auto-transpose detection**: Corrects data shape to (n_channels, n_samples)
3. **Auto-sampling rate extraction**: Reads fs from data files
4. **Multiple data variable name support**: Searches 'Data', 'data', 'EEG', etc.
5. **Comprehensive error reporting**: Detailed logs of missing data

## Recommendations

### For Immediate Use
- Run analyses on **Pat_02 and Pat_03** only (fully complete)
- Both have all phases, metadata, and sampling rates

### For Complete Dataset
- Request collaborators complete Pat_05-08 data
- Share the detailed report: `PATIENT_DATA_REPORT.txt`

## Technical Notes

- Default sampling rate is 2048 Hz (used when fs not found in data)
- Pipeline warns when using default sampling rate
- Different sampling rates are automatically accommodated
- Cache files are patient-specific, so different fs values don't cause conflicts
