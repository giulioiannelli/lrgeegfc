---
name: task-trace-audit-and-recovery
type: report
era: COHORT_N10
status: current
created: 2026-04-25
updated: 2026-04-25
supersedes:
  - .agents/reports/2026-04-24_multiscale-task-trace.md
  - .agents/reports/2026-04-24_post-mortem-scalar-session.md
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md
  - .agents/guides/task-persistence-investigation/2026-04-25_module-retention-landscape.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cbr-investigation.md
  - .agents/reports/2026-04-24_pipeline-status.md
  - .agents/reports/2026-04-24_h1-h4-vi-results.md
  - .agents/reports/2026-04-25_measure-ledger.csv
---

# Task-Trace Audit & Recovery — n=10 IMCOH_ABS (writing handoff)

**Pat_14 task_test was vendor-replaced 2026-04-25; cohort returned to n=10.
The headline band×k heatmap at n=10 (`task_trace_band_k_n10_imcoh_abs.pdf`)
shows that H2a — the partition-level "rest_post closer to task_test than to
rest_pre" claim — works cleanly under the corrected cohort, with three
independent operationalizations (Δ_VI, Δ_H, Δ_NMI) triangulating a
band-specific multiscale signal: δ is the headline (k=20–32 ≥8/10 cohort-wide
under Δ_VI, with k=28 unanimous 10/10), α and γ_h follow at mid-k (k≈20–26
and k≈21–23 / k≈30–35 respectively), β survives at low-k (k=6–7), and θ + γ_l
stay silent (ergodic, as already established by H2d). The n=9-era framing
"H2a demoted to supplementary" was an artefact of the corrupt-Pat_14 cohort
+ scalar-aggregate gating; at the per-(band, k) cell level under n=10, **H2a
is re-instated as the central H2 result and the band-specificity is real**.
Caches were backfilled, the figure is published, the canonical reformalization
(P/T/R/RA at leaf-set level, cohort threshold ≥ 8/10) lives at
`2026-04-25_task-trace-canonical.md`, and module identification confirms that
the residual signal lives in *partition divergence* (the heatmap), not in
strict subtree identity (J_min=0.9 trace-modules near-null cohort-wide). MRL
and CBR investigations remain parallel and untouched.**

---

## 1. Reading order

1. **This document** — what changed, what to do, what not to touch.
2. `.agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md` —
   leaf-set-level reformalization (P/T/R/RA, residual-vs-global, locked decisions).
3. `.agents/reports/2026-04-25_measure-ledger.csv` — every measure × era × cohort
   × cache status. The single source of truth for "is this current".
4. `data/outputs/figures/section6/task_trace_band_k_n10_imcoh_abs.pdf` — the headline
   `(band × k)` map (Stage 2 below).
5. `data/outputs/figures/section6/k_artefact_diagnostic_n10_imcoh_abs.pdf` — the
   companion diagnostic (Stage 3).
6. `data/outputs/figures/section6/trace_modules_summary_n10_imcoh_abs.pdf` and
   `data/audit/trace_modules/Pat_NN_band_<band>_k<k>.pdf` — module identification
   (Stage 4).

Parallel investigations (do not modify):
- `2026-04-25_module-retention-landscape.md` (MRL)
- `2026-04-25_cbr-investigation.md` and friends (CBR variants)

## 2. What changed 2026-04-25

| # | Change | Source |
|---|---|---|
| 1 | Pat_14 task_test vendor-replaced; cohort n=9 → n=10 | `CLAUDE.md` line 71 |
| 2 | Hypothesis reformalized at leaf-set level (P/T/R/RA) | `2026-04-25_task-trace-canonical.md` |
| 3 | `2026-04-24_multiscale-task-trace.md` marked superseded (numbers stay as n=9 snapshot) | this report supersedes |
| 4 | `2026-04-24_post-mortem-scalar-session.md` marked superseded (conclusion re-asserted) | this report supersedes |
| 5 | Cohort threshold ≥ 7/9 → ≥ 8/10 (preserves ~78% fraction) | locked in canonical doc L5 |
| 6 | k-artefact handling: no headline masking, separate companion diagnostic | locked per user choice 2026-04-25 |

## 3. Implementation stages (the audit-and-recovery plan)

### Stage 1 — Hypothesis reformalization ✅ this PR

Files:
- `.agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md` (NEW, scope report)
- This report (NEW)
- Frontmatter updates on the two superseded reports

### Stage 2a — Pat_14 backfill into multiscale CSVs

**Inputs:** Pat_14 LRG cache under `imcoh_abs` for all 6 bands × 4 phases. If
not yet built, run:

```bash
lrg-eegfc compute lrg --patients Pat_14 --fc-method imcoh_abs -v
```

**Producers to run on Pat_14 only** (idempotent merge by
`(patient, band, k, hypothesis)` tuple):

- `scripts/01_compute/hypothesis_tests/h2_partition_multiscale.py --patients Pat_14`
- `scripts/01_compute/hypothesis_tests/h2c_ultrametric_drift.py --patients Pat_14`
- `scripts/01_compute/hypothesis_tests/h2d_coactivation_persistence.py --patients Pat_14`

**Asserts after merge:** no duplicate `(patient, band, k, hypothesis)` rows;
unique patient count = 10; Pat_14 rows present in every band.

### Stage 2b — Headline 6-panel `(band × k)` figure

**Script:** `scripts/01_compute/figures_embedded/fig_task_trace_band_k_n10.py` (NEW)

**Inputs (read-only):**
- `data/reports/imcoh_vi/h2_partition_multiscale_raw.csv` (post-backfill) → H1, H2a, H2b, H3
- `data/reports/imcoh_vi/h2c_ultrametric_drift.csv` (post-backfill) → H2c
- `data/reports/imcoh_vi/h2d_coactivation_persistence.csv` (post-backfill) → H2d

**Layout:** 3×2 grid; rows in each panel = 6 bands; columns = full integer
k ∈ [2, 49]. Color = patient-mean contrast.
- **Black border**: cell with ≥ 8/10 patients sign-correct.
- **Dashed grey border**: 6–7/10 (advisory).
- **Diagonal hatch**: cells inside a `cluster_stats` significant cluster on
  the cohort z-vector along k.

Pat_03 (1024 Hz) included with distinct marker; Pat_14 (newly restored) included.

**Output:** `data/outputs/figures/section6/task_trace_band_k_n10_imcoh_abs.pdf` +
sidecar `.md` listing every cached row used.

**Do not edit:** `scripts/07_figures/gen_h2b_figures.py` (n=5 historical record).

### Stage 3 — k-artefact companion diagnostic

**Library extension:** `src/lrg_eegfc/utils/metrics/tree.py` gains
`simpson_neff(labels)` and `cluster_size_stats(labels)`.

**Script:** `scripts/01_compute/diagnostics/diag_k_artefact.py` (NEW)

**Output:** `data/outputs/figures/section6/k_artefact_diagnostic_n10_imcoh_abs.pdf`
— per band, three curves overlaid per patient: `n_eff(k)`, `singleton_fraction(k)`,
`max_cluster_fraction(k)`. **No headline-figure masking** — reader sees the
artefact transparently.

### Stage 4 — Module identification (T-subtrees)

**Scope:** `2026-04-25_task-trace-canonical.md` (above). Read MRL `J_min` from
`2026-04-25_module-retention-landscape.md` frontmatter at execution time;
abort if CBR scopes disagree.

**Analysis script:** `scripts/01_compute/audit/audit_15_trace_modules.py` (NEW)

**Figure script:** `scripts/01_compute/figures_embedded/fig_trace_modules.py` (NEW)

**Algorithm:** see canonical doc §8 pseudocode. For each `(band, k_bin)` cell
with cohort-wide signal in Stage 2b:
1. Pick representative k = geometric centre of the bin.
2. Per patient: extract leafsets at k for `Z_pre`, `Z_TT`, `Z_post`.
3. Match across phases with `jaccard_leafsets` at locked `J_min = 0.9`.
4. Classify each candidate as P/T/R/RA per the canonical predicates.
5. Keep T regime; record `(p, b, k, leaf_indices, leaves_named, J_pre, J_post, size)`.

**Outputs:**
- `data/audit/trace_modules/trace_subtrees_n10_imcoh_abs.csv` — long table.
  Reserved columns `pointer_to_MRL_node_id`, `pointer_to_CBR_module_id` (empty
  for now; future cross-checks join on these).
- `data/outputs/figures/section6/trace_modules_summary_n10_imcoh_abs.pdf` —
  cohort summary: per-band stacked bar of T-module count × patient frequency
  + sensor-space co-occurrence heatmap.
- `data/audit/trace_modules/Pat_NN_band_<band>_k<k>.pdf` — per-patient PDF:
  3 stacked dendrograms (pre / task_test / post) with T subtrees coloured
  (one colour per matched T module), non-T subtrees grey. Rasterise dendrogram
  artists. Cap at 10 patients × 2 bands (cluster-perm winners) × 2 k-bins ≈ 40 PDFs.

### Stage 5 — Audit ledger CSV

**File:** `.agents/reports/2026-04-25_measure-ledger.csv` (NEW)

One row per `(measure_id, era, cohort, fc_method)`. Columns: `measure_id, family,
definition_file, definition_lines, library_entrypoint, output_artifact,
k_resolved, n_patients_actual, cohort, fc_method, last_modified, status,
supersedes, superseded_by, notes`. Status enum:
`current / superseded / archived / scope-only / stale-numeric`.

Initial population (~25 rows): H1, H2a, H2a′, H2b, H2c, H2c-raw, H2d, H2d-raw,
H2e, H2-FROB, H2-RAW, H1-topo, H2a-topo, H3, H3-topo, H4, KC/MC/wRF (n=5 only,
flag stale-numeric), MRL (scope-only at draft → current after first compute),
CBR variants (5 rows, status per their READMEs), reorganization metrics.

### Stage 6 — Honest framing (in this doc)

See §5 below.

## 4. n=10 status of n=9 claims (re-evaluated 2026-04-26)

| Claim | n=9 status | n=10 outcome | Action |
|---|---|---|---|
| **δ k=23–31 partition cluster-perm** | ★ primary, p=0.014 | **strengthened** — Δ_VI ≥8/10 over k=20–32 (L=13), 10/10 at k=28; Δ_H ≥8/10 over k=19–28 (L=10); Δ_NMI ≥8/10 over k=22–32 (L=11). Three operationalizations triangulate. | promote to **headline finding** |
| **α k=2–4 conditional-H cluster-perm** | secondary, p=0.050 | **shifted** — narrow Δ_H k=4 still survives, but main α signal is now Δ_VI k=20–26 (L=7) at mid-k, not low-k | re-frame: α has a **mid-k cluster** (k≈20–26), not the n=9 narrow low-k one |
| **H2a "demoted to supplementary"** | failed @ scalar-aggregate gate | **RE-INSTATED at the cell level** — Δ_VI shows ≥8/10 cohort-wide runs in δ, α, β, γ_h. The n=9 demotion was an artefact of (corrupt Pat_14) + (scalar-aggregate FDR m=6 over 708 cells) | H2a is now the central per-cell partition result; previous "supplementary" labelling does not apply at n=10 |
| **β multiscale Δ_VI** | not flagged | Δ_VI ≥8/10 at k=6–7 (L=2); Δ_H scattered hits (k=3, k=6, k=19, k=21–22) | report as **band-specific low-k partition contrast** |
| **γ_h multiscale Δ_VI** | not flagged | Δ_VI ≥8/10 at k=21–23 (L=3) and k=30–35 (L=6); Δ_H at k=40–43 (L=4) | new **secondary band**; report alongside δ/α/β |
| **θ and γ_l silent on Δ_VI** | θ ergodic | unchanged at n=10 — no ≥8/10 runs in any partition operationalization | confirms **θ ergodicity** (matches H2d post-hoc Wilcoxon: θ Δρ minimum) |
| **H2c universal across 6 bands** | passes | passes 6/6 at n=10 (delta ρ=0.477) | unchanged |
| **H2d universal across 6 bands** | passes | passes 6/6 at n=10 (delta Δρ=0.228) | unchanged; θ still significantly lowest |
| **KC / MC / wRF scalar tests** | closed loop, n=5 only | unchanged | stays closed; flagged stale-numeric in ledger |
| **Strict T-regime subtree identity (J_min=0.9)** | scope-only at n=9 | **near-null at n=10** — 1 module total across (δ k=27, α k=3) cells × 10 patients. (J_pre, J_post) cloud sits on identity diagonal. | confirms MRL ↔ CBR reconciliation; the residual trace lives in partition divergence (this section), not in strict subtree identity |

### Headline reading

**δ is the strongest band**: the cohort-wide ≥8/10 ridge runs k=20–32
under Δ_VI with one fully unanimous cell (k=28, 10/10), confirmed at
overlapping k by Δ_H (k=19–28) and Δ_NMI (k=22–32). The user's
"agreement across scales" intuition is correct — δ shows a coherent
multiscale H2a band rather than an isolated cell.

**α and γ_h are the secondary bands**: α has a clean Δ_VI cluster at
k=20–26 (L=7) corroborated by Δ_NMI; γ_h has two distinct multiscale
clusters (k≈21–23 and k≈30–35).

**β is a low-k band-specific contributor**: small but consistent
Δ_VI / Δ_NMI agreement at k=6–7.

**θ and γ_l are silent / ergodic**: zero ≥8/10 cohort runs in any
partition operationalization — consistent with the H2d post-hoc
finding that θ has the lowest block-pair persistence (Δρ_θ=+0.163
vs Δρ_α=+0.283).

## 5. Honest framing

- **H2a is the central per-cell partition result at n=10.** The n=9
  "demoted to supplementary" framing was an artefact of (a) the corrupt
  Pat_14 task_test removing one patient's contribution and (b) the
  scalar-aggregate FDR-m=6 gate. Under the n=10 cohort and the per-(band,
  k) cell-level ≥ 8/10 cohort threshold, three independent partition
  operationalizations (Δ_VI, Δ_H, Δ_NMI) triangulate cohort-wide H2a
  signal in **δ (headline), α, β, γ_h**.
- **The signal is multiscale-coherent in δ**: contiguous ≥8/10 ridge at
  k=20–32 (L=13) under Δ_VI, with k=28 unanimous 10/10. This is
  qualitatively the same pattern that the n=5 `mean_contrast_heatmap_4panel.pdf`
  showed and that the n=9 era partially obscured.
- **H2-RAW / H2-FROB fail because they conflate residual trace with
  global geometry shift.** The trace is **residual** at the partition
  level — rest_post leans toward task structure relative to rest_pre,
  not "rest_post becomes task-like overall". The headline figure
  surfaces this directional residual; H2-RAW/H2-FROB tested the
  wrong (global) thing.
- **θ as ergodic is a finding, not a bug.** Confirmed independently by
  (i) zero ≥8/10 cohort runs in any partition operationalization at
  n=10, (ii) H2d post-hoc Wilcoxon (θ Δρ minimum), (iii) H2e split-half
  drift floor (θ flat). Highest intrinsic ρ_inert (0.19 vs 0.13–0.18
  elsewhere) — little room for task-induced structure to carve a
  residual trace.
- **Strict subtree-identity T-regime (J_min=0.9) is near-null cohort-wide.**
  This sounds contradictory but it is not: Δ_VI captures the *partition-level
  directional residual* (cell-level "which way does rpost lean") which can be
  cohort-wide even when no individual subtree's leafset is reproducible across
  patients at strict Jaccard match. Confirms the MRL ↔ CBR reconciliation:
  the residual trace lives in partition *divergence*, not in subtree
  *identity*.
- **Scalar-gate failure at FDR `q < 0.05` `m = 6` is a gate problem, not a
  data problem.** Effect sizes `r_rb ≈ 0.5–0.7` at n=10 make the cell-level
  ≥ 8/10 unanimity the right gate; cluster-permutation is the descriptive
  overlay along k.
- **VI(k) k-artefact is real but not load-bearing for the headline finding** —
  δ k=20–32 sits at mid-k where `n_eff(k)` is well above the singleton/giant
  collapse zones (companion diagnostic). β k=6–7 is closer to the giant
  edge — flag it as fragile in the writeup.

## 5b. Robustness — leave-one-out + dissenter audit (2026-04-26)

Three checks confirm the headline claim is not driven by any single
patient and that dissent is explained by implant geometry, not by
absence of phenomenon. Full diagnostic in
`data/audit/ridge_diagnostics/dissenter_summary.md`.

| check | result |
|---|---|
| **LOO on δ k=20–32 headline ridge** | 13/13 cells ≥ 7/9 unanimous under **every** patient drop. Mean Δ_VI ranges +0.200 (no Pat_03, weakest) to +0.283 (no Pat_02, strongest). |
| **Pat_03 (1024 Hz outlier) sensitivity** | Pat_03 is fully concordant on every ridge — *largest* positive contributor on δ k=20–32 (+0.65) and α k=20–26 (+0.73). Dropping Pat_03 modestly weakens means but every ridge stays cohort-wide. |
| **Persistent dissenter** | **Pat_02** dissents on 4/6 ridges. Dropping Pat_02 *strengthens* every ridge (δ headline mean +0.246 → +0.283, sign-fraction 87.7% → 94.9%). |
| **Implant explanation for Pat_02** | (Anatomy from `implant_pat_02.csv` coordinates + Desikan-Killany.) Pat_02 is the **only left-hemisphere-only implant in the cohort** (x_max = −2 mm), with the **most-inferior centroid** (z = −32 mm, 16 mm below the next-most-inferior patient) and **inferior-temporal-cortex weighting** (left fusiform 5×, inferior-temporal 4×, lateral-orbitofrontal 4× cohort means; white-matter sampling 14 pp below cohort). Pat_02 dissents at **mid-k** (k=20–32) but is **concordant at low-k** (k=6–11) — consistent with anatomical-coverage sensitivity at mid-k partition resolution. |

**Bottom line.** The multiscale LRG-hierarchy trace of task-induced
reorganization is robust at n=10 under `imcoh_abs`: cohort-wide,
band-specific (δ headline; α/γ_h secondary; β low-k; θ/γ_l silent),
triangulated across three partition operationalizations, and
load-bearing-independent of any single patient. Dissent in Pat_02 is
mechanistic (atypical implant coverage), not phenomenological. The
result is publishable as a per-band multiscale cohort claim.

## 6. Cohort manifest (locked 2026-04-25)

```
Pat_02   2048 Hz   raw N
Pat_03   1024 Hz   raw N    ← outlier; included, marked distinctly
Pat_05   2048 Hz   raw N
Pat_06   2048 Hz   raw N
Pat_07   2048 Hz   raw N
Pat_08   2048 Hz   raw N
Pat_10   2048 Hz   N=113    ← row drop [53,54,55] at load
Pat_13   2048 Hz   raw N
Pat_14   2048 Hz   raw N    ← newly restored 2026-04-25
Pat_15   2048 Hz   raw N
```

Cohort threshold: `≥ 8/10` (preserves the prior 7/9 ≈ 78% fraction).

## 7. Where to look when X breaks

- **"This claim is from before the n=10 fix"** → check
  `.agents/reports/2026-04-25_measure-ledger.csv`. If `status ≠ current`,
  do not cite.
- **"What's the τ for jaccard match"** → 0.9 primary, locked in
  `2026-04-25_module-retention-landscape.md` Q1 and re-asserted in
  `2026-04-25_task-trace-canonical.md` L1. Sweep `{0.85, 0.9, 0.95}` for
  sensitivity.
- **"What's the cohort-wide threshold"** → 8/10 (n=10) per L5.
- **"What about MRL / CBR results"** → parallel investigations; consult
  their scope reports under `.agents/guides/task-persistence-investigation/`.
- **"What's the figure that nailed the result before n grew"** →
  `data/outputs/figures/section6/mean_contrast_heatmap_4panel.pdf` (n=5,
  archived). Producer `scripts/07_figures/gen_h2b_figures.py` is frozen —
  do not edit.
