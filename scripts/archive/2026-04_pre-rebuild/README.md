---
name: 2026-04_pre-rebuild
type: archive-readme
era: PRE-REBUILD (MSC + MRL + dead branches)
status: archived
created: 2026-04
updated: 2026-05-28
---

# 2026-04 pre-rebuild archive

Pre-reset (April 2026) artifacts from before the IMCOH_ABS era was
established. README added 2026-05-28 to document the three sub-eras
this folder bundles. No live code imports anything here; every script
was retired before the n=10 cohort lock and the 2026-05-26 T_d
sign convention.

## Subfolders

- `msc-era/` — MSC-only compute scripts (4 files). Superseded by the
  IMCOH_ABS era 2026-04-15. Knowledge of the workflow shape survives
  in `src/lrg_eegfc/`; the per-band MSC numbers don't (volume-conduction
  bias retired).
- `mrl-era/` — Module-Retention Landscape hypothesis machinery
  (8 trace / figure scripts). The MRL multiscale measure was reframed
  into the cohesion-CBR / Grassmann-cluster-extent family by 2026-05.
  Scope reports for both directions live under
  `.agents/guides/task-persistence-investigation/`.
- `dead-branches/` — Pre-reset H2a-family audit scripts (6 files) that
  produced inconclusive verdicts under the |ImCoh|² compute. Reset
  invalidated the numerics; the scripts are kept as a record of what
  was tried.

## Recovery

Same as other archive folders: `git log --follow` traces history;
`git mv` restores. None of the scripts here will run cleanly against
the current cache layout — the cache paths and FC method routing
changed at the 2026-04-15 reset.
