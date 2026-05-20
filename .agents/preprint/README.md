---
name: preprint-folder-readme
era: IMCOH_ABS_COHORT_N10
status: current
kind: index
scope: per-band frozen result reports for the preprint manuscript
---

# `.agents/preprint/` — per-band frozen result reports

This folder gathers the **frozen per-band result reports** that the preprint manuscript draws from. Each report is the single source of truth for the band's numbers, methodology, cache provenance, controls, and open issues. Audit-script outputs (`data/audit/<probe>/cohort_summary.csv`) are the raw substrate; these reports translate them into preprint-ready narrative with every number cross-checked against the CSVs at the time of writing.

## Folder layout (5 subfolders + 4 top-level entry points; reorganized 2026-05-19)

```
.agents/preprint/
├── README.md               # this file
├── WRITING_GUIDE.md        # routing rules for new files (read before creating)
├── HANDOFF_INDEX.md        # LaTeX writing-agent entry point
├── EVALUATION_PROTOCOL.md  # protocol Claude runs on writing-agent LaTeX output
│
├── locked/                 # locked sources of truth (VERDICT_LEDGER, CONTROLS, ANATOMY_*)
├── bands/                  # the 7 fixed per-band briefs (00_cohort + 01_beta … 06_delta)
├── methods/                # long-lived methods companions (cophenet directive, Grassmann methodology, methods review)
├── directives/             # dated writing-agent directives (writing_directive_YYYY-MM-DD_*.md)
├── responses/              # dated cascade summaries, preflight notes, feedback replies
├── established_results/    # frozen methodology Q&A from the n=10 era
└── tables/                 # rendered tables (.md + .tex)
```

See [`WRITING_GUIDE.md`](WRITING_GUIDE.md) for the decision tree on where
new files go (and when to *not* create a new file but update an existing
one — the default).

## Source of truth (read before any band brief)

### Trace detection (locked 2026-05-18)
- [`locked/CONTROLS.md`](locked/CONTROLS.md) — locked 4-control battery + C5 epi-X sensitivity layer; per-probe gating rules. Verdict vocabulary: `strong trace` / `weak trace` / `no trace` per probe, combined into a per-band coverage tag (`both probes` / `only D_coph` / `only Grassmann` / `none`).
- [`locked/VERDICT_LEDGER.md`](locked/VERDICT_LEDGER.md) — per-band locked verdicts + transcribed CSV evidence + decision log + anti-revisitation clause. **No brief, figure, or LaTeX paragraph may introduce a verdict not in this ledger.** Verdicts are revised only via a dated revision entry citing a new audit.

### Anatomy localization (locked 2026-05-19)
- [`locked/ANATOMY_CONTROLS.md`](locked/ANATOMY_CONTROLS.md) — locked 4-control anatomy battery (A1 hypergeometric, A2 sampling-corrected [deferred], A3 matched-strength surrogate **mandatory**, A4 implant-geometry regression [deferred]). Verdict vocabulary: `strong localized` / `weak localized` / `not localized` per (band, probe).
- [`locked/ANATOMY_LEDGER.md`](locked/ANATOMY_LEDGER.md) — per-(band, probe) **locked** anatomy verdicts under A1+A3 (cophenet probes) and A3 alone (Grassmann probes). Lab labels A1/A3 throughout (manuscript uses compressed A1/A2; mapping note in ledger). All four trace-positive bands locked as **strong localized**. **KC-era anatomy artifacts** (`lrg_localization_anatomy/`, the "Hippocampus + left fusiform" memory claim) are retired — Hippocampus survives at β Grassmann; **left fusiform appears nowhere under the locked cluster-extent paradigm `S(b)`** (audit_72 --cluster-extent rerun 2026-05-19 pm). Five audits landed under `S(b)`: β cophenet, β Grassmann (`S(β)` — region set identical to retired `K*(β)` run), α cophenet (full + epi-X reproduces identically), γ_l Grassmann (`S(γ_l)` — 7-region occipito-temporal + frontal + medial-OFC network), δ Grassmann (`S(δ)` full + `S^epiX(δ)` — **fully disjoint** at the anatomy level, **0 shared regions**).

### Methods directive
- [`methods/methods_revision_2026-05-18_cophenet.md`](methods/methods_revision_2026-05-18_cophenet.md) — binding methods directive (KC retired, VI(k) retired, τ-sweep retired, `D_coph` adopted as canonical per-pair object).
- [`methods/methods_grassmann_cluster_extent.md`](methods/methods_grassmann_cluster_extent.md) — locked Grassmann probe methodology (subspace `U_k`, chordal distance, per-`k` Wilcoxon, cluster-extent permutation null with co-primary `LR` and cluster-mass statistics, disjunctive verdict gate). **Defines the k-independent band-level scalar `T_G^*(band) = Σ_{k ∈ C*} (−log10 p_k)`** that the manuscript reports for the Grassmann probe.

## Naming convention

`NN_<band>.md` in strength-of-results order:

| File | Band | Trace verdict | Anatomy verdict |
|---|---|---|---|
| `bands/01_beta.md` | β (13–30 Hz) | strong trace, both probes | strong localized, both probes |
| `bands/02_alpha.md` | α (8–13 Hz) | strong trace, only D_coph | strong localized (epi-X identical) |
| `bands/03_gammalow.md` | γ_low (30–80 Hz) | weak trace, only Grassmann | strong localized (temporal-cortex) |
| `bands/04_theta.md` | θ (4–8 Hz) | no trace | n/a |
| `bands/05_gammah.md` | γ_high (80–300 Hz) | no trace | n/a |
| `bands/06_delta.md` | δ (0.53–4 Hz) | weak trace, only Grassmann | strong localized (full ≠ epi-X) |

The ordering reflects the strength of the LRG-layer cohort-controlled trace under both nulls (within-baseline + matched-strength), not raw substrate counts.

## What each band report contains

A complete per-band report has the following structure:

1. **Head** — 1–2 sentences. Renormalization style.
2. **Scientific claim** — the structural-memory question being tested at this band.
3. **Cohort + substrate** — patient list, FC method, τ choice, cache paths.
4. **Probe-by-probe results** — the two locked probes:
   - **Per-pair multiscale**: `ρ_split^coph` on the cophenetic ultrametric `D_coph = cophenet(UPGMA(D(τ_max)))` at `τ_max = 1/λ_max`. Primary controls: C1 split + C2 drift + C3 matched-strength + C4 cross-probe. Sensitivity (C5): epi-X if available.
   - **Subspace**: Grassmann chordal distance `d_G(k)` on the leading-`k` Laplacian eigenmodes `U_k = span{φ_2, …, φ_{k+1}}`. Primary control: C3 matched-strength gated by cluster-extent permutation (audit_70). Sensitivity (C5): epi-X.

   The raw FC substrate (`|ImCoh|`) and raw `D(τ_max)` are reported as a **sensitivity layer** in the headline three-layer cohort table — they are *not* independent probes. Matrix-distance `T_d^(d_S/d_P)` on `D_coph` is reported as **descriptive only** (no matched-strength companion).

   For each probe: critical preamble (5-point) → method → results table → C5 epi-X sensitivity → cache + script provenance.
5. **Multiscale structure** of the trace — how it distributes across `k` (Grassmann subspace) and across cophenetic merge heights (`D_coph` per-pair).
6. **Anatomical distribution** — per-(band, probe) DK region enrichment under A1+A3 (cophenet probes) or A3 alone (Grassmann probes); locked verdicts in `locked/ANATOMY_LEDGER.md`.
7. **Patient-by-patient reading** — per-patient sign and magnitude.
8. **Sensitivity panel** — C5 epi-X (where run) + descriptive Pat_03 dropout + three-layer cohort table.
9. **Interpretation** — neuroscientific framing in plain language.
10. **Source-of-truth references** — `.agents/preprint/locked/VERDICT_LEDGER.md`, `.agents/preprint/locked/CONTROLS.md`, `data/audit/<probe>/cohort_summary.csv`.

## Rules anchoring per-band numbers

- Every numerical claim must come from a `data/audit/<probe>/cohort_summary.csv` row. Hand-quoted numbers from older reports are not allowed; the CSV is authoritative.
- `ρ_split^coph` (per-pair on `D_coph`) and Grassmann `d_G(k)` (subspace on `U_k`) are the **two matched-strength-controlled probes**. The per-band verdict is the coverage tag from `locked/VERDICT_LEDGER.md` (one of seven: `strong trace, both probes` / `strong trace, only D_coph` / `strong trace, only Grassmann` / `weak trace, both probes` / `weak trace, only D_coph` / `weak trace, only Grassmann` / `no trace`).
- Verdicts are **locked** in `locked/VERDICT_LEDGER.md`. Briefs document them; they do not re-derive them. To change a verdict, run a new audit and add a dated revision entry to the ledger.
- KC `T_KC(λ)` and VI(k) partition cuts are **retired** (2026-05-18). They are not reported in any per-band brief. Archived artefacts live at `data/audit/archive/2026_05_18/`.
- Use the **trace / anchor / reset / emergent** taxonomy; never bare "persistence" (cf. `.agents/guides/01_project/terminology.md`).

## `scripts/02_preprint/` layout (current)

Per-band figure regeneration lives here. Existing β scripts will be parameterized over `--band` to produce per-band F1–F4 figures:

- `preprint_07_beta_rho_split_figure.py` → F1: 4-control panel for `ρ_split^coph`.
- `preprint_08_beta_grassmann_figure.py` → F2: Grassmann `T_G(k)` cohort cells.
- `preprint_09_beta_grassmann_heatmap.py` → companion heatmap.

Outputs land at `data/preprint/figures/<band>/`. The audit scripts in `scripts/01_compute/audit/` (including `audit_70_grassmann_cluster_extent.py` introduced 2026-05-19) remain the source for any number cited in the ledger.

## Companion guides

- [`bands/00_cohort.md`](bands/00_cohort.md) — cross-band synthesis (the headline three-layer table + verdict matrix + per-band paragraphs + anatomy synthesis).
- [`HANDOFF_INDEX.md`](HANDOFF_INDEX.md) — single entry point for the writing agent producing LaTeX (anti-pattern checklist + per-band cheatsheet + CSV reference index).
- [`EVALUATION_PROTOCOL.md`](EVALUATION_PROTOCOL.md) — protocol Claude runs when the user pastes writing-agent LaTeX output back into a session (number-level CSV cross-check + verdict consistency check + anti-pattern scan).
- [`locked/CONTROLS.md`](locked/CONTROLS.md) + [`locked/VERDICT_LEDGER.md`](locked/VERDICT_LEDGER.md) — trace lockdown.
- [`locked/ANATOMY_CONTROLS.md`](locked/ANATOMY_CONTROLS.md) + [`locked/ANATOMY_LEDGER.md`](locked/ANATOMY_LEDGER.md) — anatomy lockdown (locked 2026-05-19).
- [`methods/methods_revision_2026-05-18_cophenet.md`](methods/methods_revision_2026-05-18_cophenet.md) — binding methods directive.
- [`methods/methods_grassmann_cluster_extent.md`](methods/methods_grassmann_cluster_extent.md) — Grassmann probe methodology + `T_G^*` scalar (locked 2026-05-19).
- `.agents/guides/04_rules/never-always-list.md` — coding + plotting + verification rules.
- `.agents/guides/01_project/terminology.md` — trace / anchor / reset / emergent taxonomy.
- `.agents/guides/05_plotting/` — figure conventions (read before any preprint figure).
- `data/audit/grassmann_cluster_extent/README.md` — cluster-extent permutation null (audit_70).
- `.agents/reports/2026-05-15_errata-corrige-matched-strength-null.md` — matched-strength source-of-truth.
