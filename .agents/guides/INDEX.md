---
name: index
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Agent Guides Index

Quick-reference index of all agent guides. Start with the category most relevant to your task.

## 01_project/ — Project understanding & workflow

| Guide | What it covers |
|-------|---------------|
| [agent-playbook.md](01_project/agent-playbook.md) | Session workflow, task patterns, how to approach work |
| [agent-structure-guide.md](01_project/agent-structure-guide.md) | Repository structure and module organization |
| [agent-tasks.md](01_project/agent-tasks.md) | Common task patterns and templates |

## 02_methods/ — Scientific methodology

| Guide | What it covers |
|-------|---------------|
| [msc-method-guide.md](02_methods/msc-method-guide.md) | Magnitude-squared coherence: computation, surrogates, sparsification |
| [imcoh-guide.md](02_methods/imcoh-guide.md) | ImCoh: volume-conduction-immune FC, why it replaces MSC for community analysis |
| [probe-bias-guide.md](02_methods/probe-bias-guide.md) | **CRITICAL**: same-probe MSC bias, enrichment quantification, debiasing rules |
| [h2-metrics.md](02_methods/h2-metrics.md) | **The VI / ρ / Δρ framework**: formulas, what each measures, what each can't prove, how they decompose the task-trace thesis |
| [time-window-guide.md](02_methods/time-window-guide.md) | Sliding-window FC analysis methodology |

## 03_implementation/ — Code & tooling reference

| Guide | What it covers |
|-------|---------------|
| [data-layout.md](03_implementation/data-layout.md) | **START HERE** for anything under `data/raw/`: canonical per-patient layout, `lrg-eegfc data normalize`, vendor → canonical rules, per-patient quirks (Pat_03 @ 1024 Hz, Pat_15 `implant_CM.xlsx`, Pat_10 vendor-filename drift, ...) |
| [cli-reference.md](03_implementation/cli-reference.md) | `lrg-eegfc` command reference and examples |
| [caching-guide.md](03_implementation/caching-guide.md) | Cache structure, paths, load-first pattern |
| [function-map.md](03_implementation/function-map.md) | Complete function lookup table by module |
| [figure-patterns.md](03_implementation/figure-patterns.md) | Figure generation templates and styling (legacy reference; `05_plotting/` is canonical) |

## 05_plotting/ — Plotting style guide (READ BEFORE PLOTTING)

| Guide | What it covers |
|-------|---------------|
| [README.md](05_plotting/README.md) | **Entry point.** Six rules + file map. Read first. |
| [legends.md](05_plotting/legends.md) | Figure-level legend placement; never per-axis for shared content |
| [colorbars.md](05_plotting/colorbars.md) | `imshow_colorbar_caxdivider` is the default; per-row colorbars on multi-axis grids |
| [multi-axis-figures.md](05_plotting/multi-axis-figures.md) | Grid layouts, shared axes, figure-level decoration |
| [colormaps-and-styles.md](05_plotting/colormaps-and-styles.md) | Colormap conventions + project-wide phase-pair colour code |
| [library-helpers.md](05_plotting/library-helpers.md) | Index of reusable helpers in `lrgsglib.plotlib` and `lrg_eegfc.visuals`; open coding tasks (helpers worth promoting) |
| [output-and-rasterization.md](05_plotting/output-and-rasterization.md) | PDF only, full vector — **never rasterise**; watermark/provenance footer is opt-in (`watermark=True`), off by default |
| [captions.md](05_plotting/captions.md) | How to write a `<figure>.md` caption *when asked* — short, plain language, three blocks |
| [per-figure-style-sheets/](05_plotting/per-figure-style-sheets/) | One sheet per figure family (heatmap-grid, scatter-with-identity, swarm-or-strip, cohort-conventions) |
| [fc_templates/](05_plotting/fc_templates/README.md) | **Canonical FC adjacency-matrix templates** (single, row-per-phase, mosaics). Each template = `.md` style sheet + `.py` script producing a vector PDF. Read this folder before any FC heatmap. |

Skill entry point: `/plotguide`. Read this folder before any new figure code.

## task-persistence-investigation/ — Active research question

| Guide | What it covers |
|-------|---------------|
| [README.md](task-persistence-investigation/README.md) | **Canonical home for every multiscale measure built to investigate task-induced reorganization in `rest_post`.** All new task-trace tooling MUST land here as a mathematically rigorous scope report (notation → predicates → properties → caveats → pseudocode → visualization → prior-tool connection → open questions) BEFORE any code. |
| [2026-04-25_module-retention-landscape.md](task-persistence-investigation/2026-04-25_module-retention-landscape.md) | **MRL** — `M̄(b, ξ) ∈ [0, 1]` cohort field counting `task_test` subtrees absent from `rest_pre` and present in `rest_post` (Jaccard match threshold τ). Multiscale, band-specific, no scalar gate. Status: `draft`. |

## ../reports/ — Investigation results

| Report | What it covers |
|--------|---------------|
| [2026-04-15_imcoh-process.md](../reports/2026-04-15_imcoh-process.md) | Full discovery narrative: problem → failed approaches → ImCoh solution |
| [2026-04-24_imcoh-results-for-writing.md](../reports/2026-04-24_imcoh-results-for-writing.md) | MSC vs ImCoh comparison for paper writing |
| [2026-04-15_imcoh-verification.md](../reports/2026-04-15_imcoh-verification.md) | Verified cell counts and flags |
| [2026-02-26_report-figures.md](../reports/2026-02-26_report-figures.md) | Publication figure inventory |
