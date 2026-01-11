# Project Status

## ✅ Implemented Features

### Core FC Pipelines
- **Correlation-based FC** (`src/lrg_eegfc/workflow_corr.py`) - Pearson correlation matrices
- **MSC-based FC** (`src/lrg_eegfc/workflow_msc.py`) - Magnitude-squared coherence
- **LRG Analysis** (`src/lrg_eegfc/workflow_lrg.py`) - Hierarchical network analysis

### Visualizations
- **MSC Summary** (`src/lrg_eegfc/visuals/msc.py`) - Matrix (1:1 aspect), network (k=0.1), percolation curves
- **LRG Full Panel** (`src/lrg_eegfc/visuals/lrg.py`) - 5 panels: Matrix, Thermodynamics, Dendrogram, Network, PSI
- **FC Comparison** (`src/lrg_eegfc/visuals/comparison.py`) - Correlation vs MSC comparison
- **Phase Reorganization** (`src/lrg_eegfc/visuals/reorganization.py`) - Network structure across phases
- **Metastable Nodes** (`src/lrg_eegfc/visuals/metastable.py`) - Sankey diagrams showing cluster evolution across tau values

### CLI Scripts
- `src/compute_corr_matrices.py` - Batch correlation computation
- `src/compute_msc_matrices.py` - Batch MSC computation
- `src/compute_lrg_analysis.py` - Batch LRG analysis
- `src/visualize_msc.py` - MSC visualizations
- `src/visualize_lrg.py` - LRG visualizations
- `src/visualize_comparison.py` - FC method comparison
- `src/visualize_phase_reorganization.py` - Phase reorganization analysis
- `src/visualize_metastable.py` - Metastable nodes Sankey diagrams

### Pipeline
- `scripts/run_full_analysis.sh` - Complete analysis pipeline (12 steps)

## 📁 Documentation Organization

### Active Development (`.agents/plans/dev/`)
- `ToDo.md` - Current tasks
- `NEXT_SESSION_START_HERE.md` - Session continuity
- `DATA_STATUS.md` - Data availability
- `PATIENT_DATA_REPORT.txt` - Detailed patient data inspection report
- `ANALYSIS_COMMANDS.md` - Quick reference commands
- `VALIDATION_PIPELINE.md` - Validation procedures

### Completed (`.agents/plans/past/`)
- `COHERENCE_FC_README.md` - MSC pipeline documentation (complete)
- `VISUALIZATION_FIXES.md` - Batch mode fixes
- `VISUALIZATION_QUALITY_FIXES.md` - MSC/LRG quality improvements
- `SESSION_MSC_NPERSEG_FIX.md` - MSC nperseg parameter fix
- `CHECKPOINT_VISUALIZATION_IMPLEMENTATION.md` - Implementation checkpoints
- `PIPELINE_GUIDE.md` - Pipeline documentation
- `LRG_VISUALIZATION_INSTRUCTIONS.md` - LRG viz specs
- `COMPARISON_VISUALIZATION_INSTRUCTIONS.md` - Comparison viz specs
- `VISUALIZATION_PLAN.md` - Original visualization plan
- `VISUALIZATION_PLAN_DETAILED.md` - Detailed specifications

### Root Documentation
- `README.md` - Project overview
- `CLAUDE.md` - AI assistant instructions

## 🎯 Next Steps

1. **Testing & Validation**
   - Test all visualizations on all patients
   - Validate batch processing
   - Check figure quality across all bands and phases

2. **Pipeline Integration**
   - Consider adding metastable visualization to full analysis pipeline
   - Optimize batch processing performance

3. **Documentation Updates**
   - Update README with metastable visualization features
   - Document usage examples
   - Add notebooks demonstrating new features
