# Archive — 2026-04 imcoh-sq-era

**Scripts from the brief `|ImCoh|²` era (2026-04-10 → 2026-04-15).
Retained for the transition record; numerically superseded by the
`imcoh_abs` era.**

See `.agents/era-map.md` and `.agents/reports/2026-04-15_imcoh-reset.md`
for the transition context.

## Contents

- `h2d_imcoh_sq_comparison.py` — comparison between `|ImCoh|²` and
  `|ImCoh|` LRG outputs, used during the reset diagnosis.
- `verify_imcoh_reset.py` — one-shot validation script that proved the
  ordering bug in `compute_lrg_analysis`.

## Status

Don't re-run. Their FC cache (`data/cache/imcoh_sq/`) is superseded by
`data/cache/imcoh/` (frequency-resolved, canonical). The reset report
summarizes everything these scripts found.
