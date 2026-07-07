---
name: agentic-folder-reorg
type: plan
era: IMCOH_ABS × COHORT_N10
status: active
created: 2026-07-04
scope: execution checklist for the `.agents/` folder reorganization — analysis (6-agent swarm) + plan approved 2026-07-04. Living doc; check off as stages land.
pointers:
  - .agents/START_HERE.md
  - .agents/preprint/WRITING_GUIDE.md
---

# `.agents/` reorganization — execution checklist

**Head.** The `.agents/` folder has no junk (zero stubs, all frontmatter present)
but ~148 of 408 files are superseded and sitting un-archived next to live work,
and three stacked generations of entry-point docs route agents into retired
eras. This plan archives the sediment **by research era**, collapses navigation
to **one entry point + one write-router**, and reframes the preprint folder
around **manuscript sections** — mechanically, without touching any scientific
verdict.

## Progress (RESUME POINT — updated 2026-07-05)

**Done:**
- **Stage 1** ✅ — `plans/` (active 13→2; 8→`developed/`, 4→`archive/setup-era/`;
  INDEX rebuilt; `msc-method-guide` citer fixed) + `writing-bundles/` (7 dup PDFs
  →`archive/superseded-renders/`; bundle README bannered `historical`).
- **Stage 2a** ✅ — `task-persistence-investigation/`: 43 retired scope-docs →
  `archive/{scalar-vi-era 23, kc-era 10, epi-negatives 4, replay-n4 4,
  opaque-null-era 2}`; 10 live + README remain; README rewritten (live/archived
  index); live citers repointed (`CLAUDE.md`/`AGENTS.md` MRL bullet, `guides/INDEX.md`,
  `lrg-framework-guide.md`, `headlines/verification/README.md`, epilepsy plan).

**Tally so far:** 54 files archived, 0 broken links (all via `git mv`/`mv`;
renames preserve history).

**▶ RESUME HERE → Stage 2b (`reports/`)**, then Stage 2c (preprint), Stage 3
(routing). Content-freshness pass stays separate/deferred.

**Citer policy reminder:** live→archived = UPDATE; archived↔archived + diary =
LEAVE. Reports being archived that cite already-moved TPI docs = leave
(archived↔archived after 2b). The ~5 KEEP reports citing moved epi/replay
scope-docs get their links fixed during 2b.

## Decisions locked (PI, 2026-07-04)

1. **Spine = manuscript sections** (Intro / Methods / Results / Discussion).
   Applied as a *navigational* map in the entry docs; physical preprint folders
   stay put for now (we are mid Results-writing — a physical `sections/` move is
   a later, separate step).
2. **Archive = by research era** (not by date). Era vocabulary below. Existing
   dated archive folders get **converted** to era buckets too.
3. **Routing = one entry point** (`START_HERE.md`) **+ one write-router table**
   (the "where does X go?" table, sole authority) **+ an ask-escape** ("if none
   fit, stop and ask the user").
4. **Scope = mechanical reorg + routing + index/pointer fixes ONLY.** Do **not**
   touch `preprint/locked/` ledgers or the γ_low→PFC verdict cascade — those are
   PI-owned and on hold in a separate chat.

## Era vocabulary (archive bucket names)

Applied per-parent as `<parent>/archive/<era>/` (keeps the existing per-parent
convention; just renames date-buckets to era-buckets). Shared vocabulary so the
archaeology is legible across the whole tree:

| Bucket | What belongs |
|---|---|
| `msc-era/` | MSC / coherence FC method, pre-ImCoh |
| `imcoh-sq-era/` | `\|ImCoh\|²` substrate + April ImCoh reset/recovery |
| `scalar-vi-era/` | April VI(k) / H1–H4 / scalar-session / CBR-MRL / functional-tree-distance / E1-pivot / eigenmode-subspace dead-ends |
| `kc-era/` | May KC / CTM / RF tree-distance + old "Section 5" manuscript numbering |
| `replay-n4/` | the pulled replay headline thread (former N4) |
| `completed/` | spent handoffs & executed directives of the **current** era (not method-retired — the one non-era bucket) |
| `setup-era/` | pre-research project/infrastructure/refactor plans (2025-11 → 2026-01) — **plans/ only** |
| `superseded-renders/` | duplicate/old figure renders — **writing-bundles/ only** (figures don't map to research eras) |

## Citer-handling policy (from the citer-web recon, 2026-07-04)

The retired files are cited almost entirely by *other retired files* and by
*diary* entries. Policy:

- **LIVE file → archived file:** UPDATE the path. (Entry docs, INDEXes, KEEP
  reports, headlines, `CLAUDE.md`/`AGENTS.md`.)
- **Archived file → archived file:** LEAVE as historical record (both end up in
  the archive; low value to rewrite, and it preserves the trail).
- **Diary → anything:** LEAVE. Diary is append-only history; never rewrite it.

---

## Stage 1 — low-coupling territories — ✅ DONE 2026-07-04

### 1a · `plans/` — `active/` holds only live plans
- `git mv` → `plans/developed/` (decisions still embodied in the shipped package):
  `2026-01-10_stage-00_notebook-consolidation.md`, `stage-01_data-qc.md`,
  `stage-02_fc-msc-corr.md`, `stage-03_lrg.md`, `stage-04_reorganization-metrics.md`,
  `stage-05_figures-overleaf.md`, `stage-05u_spatial-embedding.md`,
  `2026-01-24_spatial-3d-visualization.md` (8 files).
- `git mv` → `plans/archive/setup-era/` (point-in-time infra notes):
  `2026-01-10_data_inventory.md`, `2026-01-10_fc_params.md`,
  `2026-01-10_lrg_params.md`, `2026-01-10_notebook_audit.md` (4 files).
- KEEP in `active/`: `2026-05-08_lrg-epilepsy-research-directions.md` (⚠ verify:
  its "Direction E" may already be the shipped Grassmann probe → possibly retire)
  + this reorg plan.
- Citer fixes: rebuild `plans/INDEX.md`; repoint `msc-method-guide.md:73`
  (`fc_params` TODO).

### 1b · `writing-bundles/raw-fc/` — drop the duplicate render batch
- `git mv` → `writing-bundles/raw-fc/archive/superseded-renders/` (7 pre-`imcoh_abs`
  duplicate PDFs, byte-identical / metadata-only-diff vs the `_imcoh_abs` set,
  zero external citers): `fig1_trace_scatter_dS.pdf`, `fig3_dS_dP_convergence.pdf`,
  `fig4_structural_vs_drift.pdf`, `fig5_distance_class_examples.pdf`,
  `figS1_z_inflation.pdf`, `figS2_Td_swarm_all_distances.pdf`,
  `figS3_significance_counts.pdf`.
- KEEP: the 9 `_imcoh_abs` PDFs.
- Note on the bundle README: flag the "raw FC as headline Result 1" framing as
  superseded by the 2026-06-18 Laplacian-only lock (framing note only).

---

## Stage 2 — the big by-era archive (PENDING — review the classification below)

### `guides/task-persistence-investigation/` — ✅ DONE 2026-07-04 (43 archived; 10 live + README)
- → `archive/scalar-vi-era/` (24): the April CBR/MRL family (`2026-04-25_cbr-*`,
  `cohesion-cbr`, `containment-cbr`, `consensus-subtree`, `task-anchored-cbr`,
  `module-retention-landscape`, `mrl-vs-cbr-reconciliation`, `cophenetic-neighbourhood`,
  `trace-modules`, `multiscale-partition-coherence`, `cbr-classifier-failure-analysis`,
  `cbr-investigation`), the April substrate/rebuild-ladder diagnostics
  (`functional-tree-distance`, `residual-subspace-trace`, `edge-vs-hierarchy-discriminability`,
  `raw-fc-phase-distance`, `drift-triangle-null`, `cohort-coverage-matrix`,
  `decision-rules`, `measure-correctness-audit`, `e1-spectral-subspace-alignment`),
  **+ the 3 PI-directed ambiguous files → `scalar-vi-era/`**
  (`2026-04-27_perspective-multidirectional-multiscale`, `2026-05-08_eigenmode_embedding`,
  `2026-05-08_epi-grassmann-embedding`).
- → `archive/kc-era/` (11): `2026-05-04_lrg-trace-investigation-plan`,
  `2026-05-06_rf-k-multiscale-measure`, the KC family (`2026-05-07_kc-anchor-modules`,
  `kc-rearrangement-modules`, `kc-reset-modules`, `kc-leaf-topological-memory`,
  `2026-05-08_taxonomy_mutual_exclusivity`, `2026-05-11_kc-partition-merge-node-trace`,
  `2026-05-08_epi-cross-phase-rigidity`, `2026-05-08_band_agnostic_lrg`) + the epi
  per-node negatives (`2026-06-05_epi-node-trace-marker`, `occult-epi-node-marker`,
  `signed-magnetic-laplacian-epi`, `epi-eigenmode-localization` — negative marker
  attempts; distinct from the shipped relational-propagator R3 marker).
- → `archive/replay-n4/` (4): `2026-06-22_replay-states-multiscale-reinstatement`,
  `replay-states-v3-scale-target-pair`, `replay-sustained-and-timescale`,
  `2026-06-22_intrinsic-scale-reorganization`. (Join the existing
  `archive/2026-06/` CRC-null pair → convert those to `scalar-vi-era` or keep as
  `opaque-null-era`? — flag.)
- LIVE-citer fixes: `CLAUDE.md`/`AGENTS.md:580` (MRL bullet → drop/repoint);
  `headlines/verification/README.md:46` (replay scope restore-point → new path);
  rewrite the TPI `README.md` index (live-vs-archived table).
- Rewrite each moved file's frontmatter with `archived: 2026-07-04` + one-clause
  `reason:` (several currently rely on a *later* doc to explain why they're dead).

### `reports/` (37 keep; 45 files + the 14-file `full-verification/` archive)

**KEEP — do NOT archive (37; cited by the live spine or uncontradicted):**
`2026-05-07_epileptic-n10-revisit`, `2026-05-11_grassmann-matched-strength-verification`,
`2026-05-11_implant-geometry-and-kc-null-verification`, `2026-05-14_grassmann-epi-exclusion-sensitivity`,
`2026-05-15_errata-corrige-matched-strength-null`, `2026-05-15_grassmann-regating-no-7of10-filter`,
`2026-05-28_matched-strength-null-defense`, `2026-05-30_pair-trace-measures-methodology-audit`,
`2026-06-01_trace-concordance-vs-blind-fc`, `2026-06-05_eigenmode-localization`,
`2026-06-05_epilepsy-headline-and-occult-node-marker`, `2026-06-05_propagator-subspace-epi-recovery`,
`2026-06-05_results-assessment-and-nature-strategy`, `2026-06-05_rhocoph-grassmann-dissociation`,
`2026-06-08_epi-propagator-soz-marker`, `2026-06-08_white-matter-exclusion`,
`2026-06-11_epi-soz-marker-literature-positioning`, `2026-06-12_decimation-controls-handoff`,
`2026-06-12_propagator-distant-soz-marker`, `2026-06-12_trace-impact-and-leverage`,
`2026-06-16_epi-distant-marker-walkthrough`, `2026-06-18_inference-mark-handoff`,
`2026-06-18_inference-mark-localization`, `2026-06-22_arc-mesoscale-inference-null`,
`2026-06-22_learning-phase-own-trace`, `2026-06-22_lowgamma-focal-cingulate-encoding`,
`2026-06-22_replay-states-verification`, `2026-06-22_task-learn-length-quality-audit`,
`2026-06-22_tau-sensitivity-of-the-trace`, `2026-06-22_tau-sweep-arc-result`,
`2026-06-25_multiphase-snr-reliability`, `2026-06-25_per-node-trace-anatomy-and-heterogeneity`,
`2026-06-25_raw-vs-multiscale-trace`, `2026-06-25_snr-band-taxonomy-handoff`,
`2026-06-26_per-band-phenomenology-vision`, `2026-07-04_results-draft-review-handoff`,
`2026-07-04_results-readiness-map`. Everything else in `reports/` top-level archives by era.

**KEEP-report → moved-scope citer fixes owed in 2b** (live→archived, must UPDATE):
`epileptic-n10-revisit`→`archive/kc-era/2026-05-08_epi-cross-phase-rigidity`;
`propagator-subspace-epi-recovery`→`archive/epi-negatives/2026-06-05_epi-node-trace-marker`;
`epilepsy-headline-and-occult-node-marker`→`archive/epi-negatives/{2026-06-05_occult-epi-node-marker,epi-node-trace-marker}`;
`eigenmode-localization`→`archive/epi-negatives/2026-05-08_epi-eigenmode-localization`;
`replay-states-verification`→`archive/replay-n4/2026-06-22_replay-states-v3-scale-target-pair`
(all under `guides/task-persistence-investigation/`).

⚠ **Classification is by content, not date** — a 2026-05-05 file can be KC-era
*or* early-current-LRG. Proposed mapping (review before executing):
- → `archive/imcoh-sq-era/`: `2026-04-15_imcoh-recovery-plan`, `imcoh-reset`,
  `2026-04-18_epileptic-imcoh-final`.
- → `archive/scalar-vi-era/`: the April VI/H1–H4/scalar/E1 cluster (14 files) +
  `2026-04-25_measure-ledger.csv`.
- → `archive/kc-era/`: the May 5–9 KC/CTM/§5 cluster (18) + `full-verification/`
  (14) + 3 May transitional (`matched-strength-surrogate-batch`,
  `full-null-verification-evidence-review`, `pre-preprint-cleanup`).
- → `archive/completed/`: the 7 June session-handoffs whose job is done
  (`wm-exclusion-continuation`, `cross-phase-taxonomy`, `consolidation-arc-handoff`,
  `n2-flagship-handoff`, `cophenetic-gate-presence-vs-consistency`,
  `per-band-consistency-taxonomy-lock`, `trace-heterogeneity-handoff`).
- **Convert existing dated archives**: `reports/archive/2026-04/` (14) → split
  imcoh-sq-era / scalar-vi-era; `reports/archive/2026-05/` (result-2, section-5-critical,
  trace-minus-epi) → kc-era **except** ⚠ `2026-05-05_result-2-lrg-beta-trace.md`
  is the LRG-β lineage (current science, old writeup) + a **live landmark** cited
  by `START_HERE.md` + `bands/01_beta.md` — classify as early-current, update its
  2 live citers.
- Set `status: superseded` + `superseded_by:` on the two known-wrong-headline
  files before archiving (`2026-06-12_cross-phase-taxonomy`,
  `2026-06-25_per-band-consistency-taxonomy-lock`).
- Add `reports/README.md` (newest-first, "current handoff = 2026-07-04_results-readiness-map",
  era-bucket legend).

### `preprint/` (38 keep; ~32 archive)
- **New** `preprint/directives/archive/{scalar-vi-era,kc-era,completed}/` +
  `preprint/responses/archive/completed/` + `preprint/methods/archive/` +
  `preprint/headlines/verification/archive/`.
- → `directives/archive/completed/` (spent one-shot writing directives, ~17):
  the `2026-05-18…2026-06-04` beta/alpha/anatomy/methods/substrate directives
  (executed against the retired per-band Results architecture).
  ⚠ **Extract-first**: lift the reusable guards out of the two `2026-06-03`
  substrate-figure directives into the relevant headline before archiving.
- → archive: the superseded `2026-06-25_headline-restructure-broadband-snr`
  directive (reversed within 48h), the misfiled `investigation_directive_2026-06-08`,
  the 6 executed `responses/`, `methods/methods_section_review_2026-05-19`,
  `verification/verify_stage2_arc_duration` (RESOLVED), both `established_results/`.
- ⚠ **Extract-first**: lift still-live methods items (B1 sign, B3 citation, τ_min)
  out of `METHODS_AUDIT_ISSUES.md` before archiving it.
- Fix ~11 `status:` frontmatter drift cases so the freshness signal is trustworthy.

---

## Stage 3 — routing consolidation (PENDING)

- **One entry point:** make `START_HERE.md` canonical; add the write-router table
  (below) + forward pointers (headlines/README, 2026-06-26 vision, 2026-07-04
  readiness-map); demote `README.md`, `preprint/HANDOFF_INDEX.md`, and the
  `CLAUDE.md`/`AGENTS.md` "Start points" list to redirect stubs.
- **Manuscript-section map** in `START_HERE` + `preprint/README`: Methods →
  `preprint/methods/` + `locked/*CONTROLS`; Results → `preprint/headlines/`
  (R1/R2/R3) + `directives/results_paragraphs/*.tex`, underpinned by
  `locked/*LEDGER` + `bands/` (provenance); Intro/Discussion → new stubs.
- **Write-router table** (sole authority, lives in `START_HERE`):

  | Output | Destination |
  |---|---|
  | analysis result / session handoff | `reports/YYYY-MM-DD_<topic>.md` + bump `reports/README.md` current-handoff line |
  | new measure (before code) | `guides/task-persistence-investigation/YYYY-MM-DD_<measure>.md` |
  | manuscript verdict change | dated entry in `preprint/locked/VERDICT_LEDGER.md` + matching `bands/NN_<band>.md` |
  | Results narrative / headline | `preprint/headlines/NN_*.md` |
  | Results prose (LaTeX) | `preprint/directives/results_paragraphs/RX.Y_<slug>.tex` |
  | Methods content | `preprint/methods/` |
  | writing directive / cascade note | `preprint/directives/…` / `preprint/responses/…` |
  | diary / figure | `diary/YYYY-MM-DD.md` / `data/preprint/figures/…` |
  | superseded material | `<parent>/archive/<era>/` |
  | **none of these fit** | **stop and ask the user before creating** |

- **Index rebuilds:** `guides/INDEX.md` (add `04_rules/`, the 3 canonical methods
  docs, `network_templates/`), `plans/INDEX.md` (done in Stage 1a).
- **Pointer fixes (7):** `era-map.md`→`2026-04-22_migration-plan-dryrun.md`;
  `README.md`→`guides/02_methods/imcoh-guide.md`; 4 dead `guides/INDEX.md`→reports
  links; TPI README frontmatter pointer.
- **Downgrade status** on stale-subject guides (`h2-metrics`, `time-window-guide`
  → superseded/historical); fix `cli-reference.md` (`--fc-method` default
  `imcoh_abs`, 5 choices; `cli/plot/` package).

---

## Content-freshness pass — SEPARATE, flagged (NOT part of the mechanical reorg)

These are content/correctness, some PI-owned — recommend a distinct pass:
- **γ_low→PFC overclaim** cascade (`headlines/01_trace`, `RESULTS_OUTLINE`,
  `bands/00_cohort`, `03_gammalow`) — **PI-owned, on hold.** Untouched here.
- `EVALUATION_PROTOCOL.md` still enforces the *retired* multi-region anatomy
  anti-patterns (would reject correct β→OFC prose).
- `cli-reference.md` factual (pre-ImCoh CLI) — folded into Stage 3.
- Commit the untracked load-bearing `references/{acronyms,math_commands}.tex`
  + `results-section-full-draft.tex`.
- `patient-warnings.md` missing Pat_15 (canonical dropout exception), Pat_02/06/08.
- Diary gap 2026-06-22 → 2026-07-04 (backfill from git log + reports).

## Provenance
- Analysis: 6-agent read-only swarm, 2026-07-04 (reports / preprint / TPI /
  reference-guides+plans / small-subtrees+top-level / routing-integrity).
- Citer-web recon: `rg` citer maps, 2026-07-04 (this session).
