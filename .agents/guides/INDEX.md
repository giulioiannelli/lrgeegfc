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
| [AGENT_PLAYBOOK.md](01_project/AGENT_PLAYBOOK.md) | Session workflow, task patterns, how to approach work |
| [AGENT_STRUCTURE_GUIDE.md](01_project/AGENT_STRUCTURE_GUIDE.md) | Repository structure and module organization |
| [AGENT_TASKS.md](01_project/AGENT_TASKS.md) | Common task patterns and templates |

## 02_methods/ — Scientific methodology

| Guide | What it covers |
|-------|---------------|
| [MSC_METHOD_GUIDE.md](02_methods/MSC_METHOD_GUIDE.md) | Magnitude-squared coherence: computation, surrogates, sparsification |
| [IMCOH_GUIDE.md](02_methods/IMCOH_GUIDE.md) | ImCoh: volume-conduction-immune FC, why it replaces MSC for community analysis |
| [PROBE_BIAS_GUIDE.md](02_methods/PROBE_BIAS_GUIDE.md) | **CRITICAL**: same-probe MSC bias, enrichment quantification, debiasing rules |
| [H2_METRICS.md](02_methods/H2_METRICS.md) | **The VI / ρ / Δρ framework**: formulas, what each measures, what each can't prove, how they decompose the task-trace thesis |
| [TIME_WINDOW_GUIDE.md](02_methods/TIME_WINDOW_GUIDE.md) | Sliding-window FC analysis methodology |

## 03_implementation/ — Code & tooling reference

| Guide | What it covers |
|-------|---------------|
| [DATA_LAYOUT.md](03_implementation/DATA_LAYOUT.md) | **START HERE** for anything under `data/raw/`: canonical per-patient layout, `lrg-eegfc data normalize`, vendor → canonical rules, per-patient quirks (Pat_03 @ 1024 Hz, Pat_15 `implant_CM.xlsx`, Pat_10 vendor-filename drift, ...) |
| [CLI_REFERENCE.md](03_implementation/CLI_REFERENCE.md) | `lrg-eegfc` command reference and examples |
| [CACHING_GUIDE.md](03_implementation/CACHING_GUIDE.md) | Cache structure, paths, load-first pattern |
| [FUNCTION_MAP.md](03_implementation/FUNCTION_MAP.md) | Complete function lookup table by module |
| [FIGURE_PATTERNS.md](03_implementation/FIGURE_PATTERNS.md) | Figure generation templates and styling |

## ../reports/ — Investigation results

| Report | What it covers |
|--------|---------------|
| [IMCOH_PROCESS_REPORT.md](../reports/IMCOH_PROCESS_REPORT.md) | Full discovery narrative: problem → failed approaches → ImCoh solution |
| [IMCOH_RESULTS_FOR_WRITING.md](../reports/IMCOH_RESULTS_FOR_WRITING.md) | MSC vs ImCoh comparison for paper writing |
| [IMCOH_VERIFICATION_RESULTS.md](../reports/IMCOH_VERIFICATION_RESULTS.md) | Verified cell counts and flags |
| [REPORT_FIGURES.md](../reports/REPORT_FIGURES.md) | Publication figure inventory |
