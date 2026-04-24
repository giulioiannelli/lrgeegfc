---
name: surface-multiscale-trace
type: plan
era: COHORT_N9
status: draft
created: 2026-04-25
updated: 2026-04-25
pointers:
  - .agents/reports/2026-04-24_multiscale-task-trace.md
  - .agents/reports/2026-04-24_h1-h4-vi-results.md
  - .agents/reports/2026-04-24_pipeline-status.md
  - .agents/reports/2026-04-24_post-mortem-scalar-session.md
  - .agents/guides/04_rules/never-always-list.md
---

# Surface the multiscale task-trace from existing |ImCoh| × n=9 evidence

**Head: the signal is already in the cache. Read the three canonical
reports + the two hypothesis-test outputs, synthesize a band-by-scale
map, write the paper. No new scalar hypothesis tests.**

## Context

What we're trying to show: `task_test` leaves a band-specific, multiscale
structural trace in `rest_post` LRG dendrograms, cohort-wide (≥ 7/9) at
n=9 under `imcoh_abs`.

What we already have (per `.agents/reports/2026-04-24_multiscale-task-trace.md`):

- **H1 (VI(k)):** θ unanimous cohort-wide — strongest band.
- **H2a:** demoted (single-patient fluctuations; not the primary claim).
- **H2c (cophenetic drift):** passes all 6 bands × 9 patients — too
  coarse to be band-specific on its own.
- **H2d (coactivation persistence):** θ has the least block memory —
  genuine band-specific cohort-wide result.
- **Partition-multiscale cluster-perm:** δ k=23–31 Δ_VI p=0.014;
  α k=2–4 Δ_H p=0.050.
- **Band typology:** per-band signature of `h2_band_typology.py`.

What we do NOT need: new scalar gate at q < 0.05 FDR m=6. Why:
`.agents/reports/2026-04-24_post-mortem-scalar-session.md` — unreachable
for effect sizes of rb ≈ 0.5–0.7 at n=9.

## Sequence

### Stage A — read existing evidence (no compute)

Read in order and take notes:

1. `.agents/reports/2026-04-24_h1-h4-vi-results.md` — VI(k) unanimity tables.
2. `.agents/reports/2026-04-24_multiscale-task-trace.md` — writing handoff (729 lines, has the regen recipes + band typology).
3. `.agents/reports/2026-04-24_pipeline-status.md` — era index; confirms every artifact cited is IMCOH_ABS / COHORT_N9 current.
4. `scripts/01_compute/hypothesis_tests/h2_partition_multiscale.py` — cached outputs under `data/reports/imcoh_vi/`.
5. `scripts/01_compute/hypothesis_tests/h2d_coactivation_persistence.py` — cached outputs.
6. `scripts/01_compute/hypothesis_tests/h2c_ultrametric_drift.py` — for the cohort-wide drift signal across all bands.

Deliverable: one note per report in the diary (`.agents/diary/`), each
with renormalization head + 3–5 bullets of "what's in it".

### Stage B — build the band × scale map (minimal new work)

Using already-cached data only, produce a single synthesis figure:
the **band × scale map** showing where each of the known effects lives.

- x-axis: scale (k or h_rel).
- y-axis: band (δ θ α β γ_l γ_h).
- colored cells: the evidence type (VI unanimity, cluster-perm p, H2d
  persistence, etc.) — aggregated from cached outputs only.
- No new statistical test.

Output:
- `data/outputs/figures/multiscale/band_scale_map.{pdf,png}` + sidecar
  `.md` with the legend and the cited sources.
- Generator script: `scripts/01_compute/figures_embedded/fig_band_scale_map.py`
  (new, but minimal — aggregates cached artifacts).

### Stage C — descriptive writing

Write a synthesis section (for paper or internal handoff) with the
band × scale map + renormalized prose:

- Head: cohort-wide multiscale trace is present; it's band-specific;
  here's the landscape.
- Body per band: cite the specific k-range or h-range where effect
  holds + patient count + source artifact.
- Caveats: Pat_03 (1024 Hz) flagged; H2a demoted; n=9 constrains the
  claim to descriptive + cohort-pattern rather than single-scalar test.

Deliverable: `.agents/reports/YYYY-MM-DD_multiscale-trace-synthesis.md`
(renormalized, short head, body expands).

## Pass criterion (descriptive, not gated)

- The band × scale map shows a *visibly* heterogeneous landscape —
  some bands populate the map, others don't.
- Each cited effect is traceable to a cached file + frontmatter era.
- The prose describes the landscape without a single-gate pass/fail
  frame.

**No FDR q < 0.05 m=6 gate.** See post-mortem.

## Reused utilities (do NOT reinvent)

| Helper | Location |
|:-------|:---------|
| Wilcoxon / FDR / rank-biserial / bootstrap / cluster-stats | `lrg_eegfc.utils.metrics.hypothesis` |
| Tree helpers | `lrg_eegfc.utils.metrics.tree` (`tree_internal_nodes`, `jaccard_leafsets`, `fcluster_at_h_rel`, `h_log_grid`) |
| KC / MC / wRF (descriptive, if useful) | `lrg_eegfc.utils.metrics.tree_distance` |
| LRG loader | `lrg_eegfc.workflow.lrg.load_lrg_result` |
| FC loader | `lrg_eegfc.workflow.fc.load_fc_matrix` |
| Band / patient constants | `lrg_eegfc.config.const` |
| Report paths | `lrg_eegfc.config.paths` |

## Guardrails

- **Cohort:** COHORT_N9 (Pat_02, 03, 05, 06, 07, 08, 10, 13, 15).
- **FC:** `imcoh_abs` only.
- Pat_14 excluded. Pat_03 flagged (1024 Hz). Pat_10 task rows [53,54,55] dropped.
- No new FC or LRG cache computation needed.
- No pooling of metrics into consensus.
- No new scalar gate.

## Out of scope (deferred)

- Replication on held-out cohort (separate plan, for writeup §15).
- New FC methods.
- New LRG variants.
- Any Stage 5 scalar gauntlet.

## Verification

```bash
# Stage A — read-only; no commands, but confirm artifacts exist:
ls .agents/reports/2026-04-24_{h1-h4-vi-results,multiscale-task-trace,pipeline-status}.md
ls scripts/01_compute/hypothesis_tests/{h2_partition_multiscale,h2d_coactivation_persistence,h2c_ultrametric_drift}.py
ls data/reports/imcoh_vi/ data/reports/imcoh/ 2>/dev/null

# Stage B — new synthesis figure
python scripts/01_compute/figures_embedded/fig_band_scale_map.py -v

# Stage C — synthesis report exists
test -f .agents/reports/$(date +%Y-%m-%d)_multiscale-trace-synthesis.md
```

## Critical files

- `/home/giulio/Documents/research/neural_networks/lrgeegfc/.agents/reports/2026-04-24_multiscale-task-trace.md`
- `/home/giulio/Documents/research/neural_networks/lrgeegfc/.agents/reports/2026-04-24_h1-h4-vi-results.md`
- `/home/giulio/Documents/research/neural_networks/lrgeegfc/.agents/reports/2026-04-24_pipeline-status.md`
- `/home/giulio/Documents/research/neural_networks/lrgeegfc/scripts/01_compute/hypothesis_tests/h2_partition_multiscale.py`
- `/home/giulio/Documents/research/neural_networks/lrgeegfc/scripts/01_compute/hypothesis_tests/h2d_coactivation_persistence.py`
- `/home/giulio/Documents/research/neural_networks/lrgeegfc/scripts/01_compute/hypothesis_tests/h2c_ultrametric_drift.py`
- `/home/giulio/Documents/research/neural_networks/lrgeegfc/src/lrg_eegfc/utils/metrics/hypothesis.py`
- `/home/giulio/Documents/research/neural_networks/lrgeegfc/src/lrg_eegfc/config/paths.py`
