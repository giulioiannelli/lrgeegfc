# Archive — 2026-04 failed-scalar-session

**Scripts from the April 2026 gauntlet that tried to prove the task trace
with a single per-(patient, band) scalar test. All failed the pre-reg
gate (n=9 cannot clear q<0.05 FDR m=6 at rb ≈ 0.6 effect sizes).**

Preserved for historical record. **Do not re-run.** Do not cite their
outputs. See `.agents/reports/2026-04-24_post-mortem-scalar-session.md`.

## Contents

- `stage0b_dmax_stability.py` — dmax(phase) comparability diagnostic
  (result: ≈ 0.99 constant, h_rel cuts are comparable across phases —
  this finding survives; the script is historical).
- `stage1_14metrics_rerun.py` — the 14 MSC-era metrics re-run under
  `imcoh_abs`. 0/14 pass pre-reg.
- `stage3_tree_distance.py` — KC / MC / wRF implementations. 0/8 pass.
- `stage4_kc_exploration.py` — KC λ-sweep per band. Per-band signatures
  differ but 6–7/9 ceiling stays.
- `diag_pi_modules.py` — abandoned π-weighted module diagnostic.
- `gallery_cbr_modules.py` — abandoned cluster-birth-retention gallery.
- `h2_cluster_birth_retention.py` — abandoned H2 sub-hypothesis.

## Imports

All scripts import from `lrg_eegfc.utils.metrics.hypothesis`
(post-3.8 elevation). They would run from this archive location if
invoked — but per the post-mortem, **don't**.
