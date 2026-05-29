---
name: start-here
type: guide
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-24
updated: 2026-05-28
pointers:
  - .agents/era-map.md
  - .agents/preprint/WRITING_GUIDE.md
  - .agents/reports/archive/2026-05/2026-05-05_result-2-lrg-beta-trace.md
  - .agents/reports/2026-05-07_epileptic-n10-revisit.md
  - data/reports/section_5_lrg_trace/README.md
---

# Start Here

Entry point for agents. Points to the current state of the scientific
pipeline and the central numerical results.

## Current scientific state (2026-05-28)

- **Era**: `IMCOH_ABS × COHORT_N10`.
  - FC carrier: `imcoh_abs` = `⟨|ImCoh(f)|⟩_f` (Ewald 2012 magnitude-averaged
    Nolte 2004 imaginary coherency). Volume-conduction immune.
    MSC kept accessible for diagnostics; do not use it for current
    hypotheses.
  - Cohort: 10 sEEG patients, locked 2026-04-25 after Pat_14's vendor
    `task_test.mat` replacement: `Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15`.
    All cross-phase tests now run at n=10. The earlier n=9 framing
    (Pat_14 excluded) is retired.
- **Locked T_d sign convention (2026-05-26)**: at every layer (raw FC,
  D_coph, KC, Grassmann)
  `T_d := d(rest_pre, task) − d(task, rest_post)` — **positive = TRACE**,
  negative = ANTI-TRACE, zero = NO TRACE. Wilcoxon trace-direction tests
  use `alternative='greater'`. Per-patient counts use `(T > 0).sum()`.
  Surrogate upper-tail p = `mean(s ≥ obs)`.
- **Preprint probe for raw FC and D_coph (locked 2026-05-26)**:
  `ρ_split = Spearman(Δ_task, Δ_rest)` (audit_33 / audit_63 cross-phase
  rank correlation, under C3 matched-strength gate). NOT the
  `T_d^(d_S)` triangle scalar (that family stays in the codebase as a
  diagnostic; the manuscript reports `ρ_split` + paired Wilcoxon p + LOO).

## Central numerical results

- **Result 2 (Section 5 headline)** — LRG β-band trace at n=10.
  `.agents/reports/archive/2026-05/2026-05-05_result-2-lrg-beta-trace.md`. Survives
  matched-strength surrogate (audit_66), epi-zone exclusion (audit_67),
  cluster-extent permutation (audit_70). LOO robust (Pat_02-driven but
  not single-patient-leveraged).
- **Grassmann verdicts (audit_70, locked 2026-05-19 pm)**:
  β strong `p_mass=0.005` LOO=0.005 (Pat_02), γ_l strong LOO=0.040 (Pat_05),
  δ strong LOO=0.055 (Pat_08) — **δ is single-patient-leveraged, flag in
  manuscript text**.
- **Epileptic-zone results (n=10)**: δ cross-probe 1.55× is known-biology
  confirmation, not novel discovery; β V1 retired. See
  `.agents/reports/2026-05-07_epileptic-n10-revisit.md`.
- **Section-5 measure index**:
  `data/reports/section_5_lrg_trace/README.md`.

## Where to begin

1. **If you are writing manuscript text**: the central routing document is
   `.agents/preprint/WRITING_GUIDE.md`. It defines the 5-subfolder
   structure (`locked/`, `bands/`, `methods/`, `directives/`, `responses/`)
   and the update-existing-not-create-new default. Read
   `.agents/preprint/HANDOFF_INDEX.md` for the live status of every
   verdict + the EVALUATION_PROTOCOL.
2. **If you are running new analysis**: read
   `.agents/guides/01_project/agent-playbook.md` for the session workflow,
   `.agents/guides/03_implementation/data-layout.md` for per-patient
   quirks (Pat_03 1024 Hz handled at config layer; Pat_10 task rows
   53/54/55 dropped at load; Pat_14 vendor-replaced 2026-04-25). Then
   `.agents/era-map.md` to era-tag every artefact before citing or
   re-running.
3. **If you are adding a method**: every new multiscale task-trace measure
   lands under `.agents/guides/task-persistence-investigation/` as a
   mathematically rigorous scope report (notation → predicates →
   properties → caveats → pseudocode → visualization → connection-to-prior-tools
   → open questions) BEFORE any code. See that folder's README.md for the
   required structure.

## Always-applicable rules

- Matched-strength surrogate null is **mandatory** before any FC-derived
  cohort claim. Within-baseline / split-half / drift nulls are
  diagnostics, not verification.
- Brutal scientific honesty: any measure not yet matched-strength tested
  is marked "unverified" in writeups and chat.
- Default posture: critical questioning of every methodology. If a
  control hasn't been run, say so on the first line.
- Library-first: before writing a helper in a script, grep
  `src/lrg_eegfc/`. Helpers with ≥2 callers must be promoted.

## Quick references

- `.agents/era-map.md` — era landmarks with what-invalidated-what
- `.agents/preprint/WRITING_GUIDE.md` — preprint routing source of truth
- `.agents/preprint/HANDOFF_INDEX.md` — live verdict status
- `.agents/preprint/EVALUATION_PROTOCOL.md` — gating protocol (no
  cross-band BH-FDR; per-band independent controls)
- `.agents/guides/INDEX.md` — full guide index
- `.agents/guides/03_implementation/cli-reference.md` — `lrg-eegfc` CLI
- `.agents/guides/03_implementation/caching-guide.md` — cache layout
- `data/reports/section_5_lrg_trace/README.md` — Section-5 measure index

## Notebook header (canonical)

```python
from lrg_eegfc.notebook import *
move_to_rootf(pathname="lrgeegfc")
```

## Key invariants

- Cache-first: never recompute in visualization scripts or notebooks.
- Cache naming is parameter-sensitive (MSC: `sparsify-*`, `nperseg-*`,
  `n_surrogates-*`; ImCoh: `nperseg-*` only, no surrogates).
- `fs` defaults to 2048 Hz; Pat_03 (1024 Hz) is handled at config layer
  via `FS_OVERRIDES` + `nperseg_for_fs(fs)` in `src/lrg_eegfc/config/const.py`.
  Pat_03 is a full cohort member at analysis layer (no outlier flag,
  no dropout robustness check).
- Pat_10 is canonically 113-channel in every phase via load-time drop
  (`PATIENT_CHANNEL_DROP` in `config/const.py`). Raw files unmodified.
- TRACE / ANCHOR / RESET / EMERGENT taxonomy is the disambiguation
  vocabulary; bare "persistence" is acceptable in unambiguous
  task-trace contexts but never in cross-phase taxonomy tables or
  mixed-band paragraphs. See `.agents/guides/01_project/terminology.md`.
