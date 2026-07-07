---
name: trace-heterogeneity-handoff
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-06-25
updated: 2026-06-25
pointers:
  - .agents/reports/2026-06-25_cophenetic-gate-presence-vs-consistency.md
  - data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
  - data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv
  - data/audit/localization_atlas/per_patient_decomp.csv
  - data/audit/epi_stratified/decimation_control_per_patient.csv
  - scripts/01_compute/audit/audit_63_split_baseline_surrogate.py
  - scripts/01_compute/audit/audit_83_localization_matched_strength.py
  - scripts/01_compute/audit/audit_85_wm_decimation_control.py
  - scripts/02_preprint/preprint_23_bands_perpair_rho_null_forest.py
  - .agents/preprint/locked/VERDICT_LEDGER.md
  - .agents/preprint/locked/CONTROLS.md
---

> **Head — resume token: "continue the trace-heterogeneity / per-node investigation".**
> We are making sense of the cophenetic (and Grassmann) task-trace results. The
> per-band trace/no-trace verdict is fragile (an internal machinery issue, NOT a
> paper claim) and there is large per-patient heterogeneity in who shows the
> trace. **Two top-down explanations for the heterogeneity were REFUTED this
> session** (OFC/network coverage; epileptic-gray contamination — fails the
> mandatory decimation control). **The PI does not accept "it's the patient" as a
> result and does not accept the negative result.** The agreed next move is a
> **bottom-up, per-NODE** analysis: decompose ρ^coph into per-pair contributions,
> establish per patient which **nodes carry the trace vs are anti-trace**, pool
> nodes across the cohort, and find what the **non-carrier / anti-trace nodes have
> in common** (anatomy / tissue / epilepsy / sampling / depth). That common
> property — if it survives a random-node control — is the publishable explanation
> of the heterogeneity. **Do this BEFORE any further writeup or cascade.**

---

## 0. The resume contract

The PI will `/compact` and then say "continue the …". When resuming:
1. Read this file + `2026-06-25_cophenetic-gate-presence-vs-consistency.md`.
2. Per the task-persistence rule (`.agents/guides/task-persistence-investigation/README.md`),
   write a **mathematically rigorous scope report for the per-node decomposition
   tool BEFORE writing code** (notation → predicates → properties → caveats →
   pseudocode → connection-to-prior-tools → open-Qs). The §5 spec below is the seed.
3. Then implement + run, honoring every constraint in §4.

Python for all compute: `/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python`
(the `conda run -n lapbrain` wrapper errors on a g++ hook — call the interpreter directly).

## 1. Origin — the PI's complaint that opened this

Starting figure: `data/preprint/figures/all_bands/patients_forest_plot/fig_bands_perpair_rho_null_forest_coph_sort-patient.pdf`.
Per-patient ρ^coph forest, grouped by band. The PI's points:
- γ_l "looks as traced as β" yet the headline says γ_l has no trace — **unsellable**.
- Pat_06 traces in all bands and strongest — must be explained.
- The p-value gate is the wrong lens; extract the **per-band trace message** as a
  meaningful, presentable thing.

## 2. Established this session (keep)

**Metric.** `ρ_split = Spearman(Δ_task, Δ_rest)`, `Δ_task = D^tt − D^pre_A`,
`Δ_rest = D^post − D^pre_B`, on **disjoint** rsPre halves (audit_63). It is a
baseline-referenced double-difference, NOT a raw similarity. `>0 = trace`.

**The cophenetic gate is one-axis (presence vs consistency).** The locked verdict
is `C3 paired one-sided Wilcoxon p<0.05 else "no trace"` (`CONTROLS.md:228`). That
Wilcoxon tests cohort *direction-consistency*; it cannot tell "absent" from
"split". Per-band ground truth (`per_patient_per_band.csv`, n=10):

| band | pos-clear / neg-clear | gate Wilcoxon p | presence binom p | honest read |
|---|---|---|---|---|
| β | 7 / 1 | 0.005 ✓ | 4e-7 ✓ | cohort-wide |
| α | 5 / 0 | 0.002 ✓ | 6e-5 ✓ | cohort-wide |
| γ_l | 5 / 1 | 0.116 ✗ | 6e-5 ✓ | **present but split** (NOT "no trace") |
| δ | 4 / 2 | 0.278 ✗ | 1e-3 ✓ | present but split |
| γ_h | 4 / 4 | 0.246 ✗ | 1e-3 ✓ | bidirectional |
| θ | 2 / 2 | 0.722 ✗ | 0.086 ✗ | absent (only one) |

"presence binom" = `binomtest(pos_clear, 10, 0.05, 'greater')` = more patients
clear own matched-strength null than the 5% chance baseline. **Only θ fails it.**
No BH across bands is applied (retired 2026-05-20).

**Integration of the two probes (per-patient).**
- Patient axis ≫ band axis: between-patient/between-band variance ratio **7.8×
  (cophenet), 2.4× (Grassmann)**.
- A latent "trace propensity" PC1 explains **39%** of the 12 cells (6 bands × 2
  probes); order: Pat_06 > 03 > 08 > 05 > 02 > 07 > 13 > 10 > 14 > 15.
- **Nesting**: on cophenet, the clearer-set of δ, θ, γ_l, γ_h is a **subset of
  β's** clearers; γ_l = core {02,03,05,06,08} ⊂ β's {…,07,13}. **α is the lone
  exception** (adds Pat_10, Pat_14, otherwise non/anti). This is why "γ_l tracers
  are also β tracers".
- **Probes are COMPLEMENTARY, not redundant** (PI requirement — they must not tell
  the same story): cophenet declares a trace in **α + β** (pairwise hierarchy
  persists); Grassmann in **β + γ_l** (+δ weak) (leading subspace persists); they
  share only **β**. α = ρ^coph-only; γ_l = Grassmann-only. Per-patient cross-probe
  agreement is high only where both fire (β +0.54, γ_l +0.58), ~0 elsewhere.

**Figure updated.** `preprint_23` now has `--verdict two-axis` (default): band
letter green=cohort-wide / amber=present-but-split / grey=absent, with a
presence/consistency glyph (full / half / open disc). `--verdict binary` keeps the
legacy colouring. PDFs regenerated (sort-patient, sort-rho). Counts + both p's go
to stdout/caption (no-in-axes-text rule).

## 3. REFUTED this session (do NOT resurrect without new evidence)

- **OFC / β-network coverage does NOT explain who traces.** PC1 ~ n_OFC = +0.16.
  Pat_14 has 15 OFC contacts and is the 9th tracer; Pat_03/08 have 0 OFC and rank
  2nd/3rd. (`beta_trace_anatomy_coverage/per_patient.csv` also showed Pat_10, a
  β-anti patient, samples the β-net the MOST.)
- **Epileptic-gray contamination does NOT explain it — fails the decimation
  control.** The correlation looked strong (n_epi_gray ~ β z = −0.63; total n_epi
  = −0.09, i.e. specifically gray), BUT excluding epi-gray contacts and recomputing
  does NOT raise the non-tracers' trace beyond size-matched random node removal
  (`epi_stratified/decimation_control_per_patient.csv`, β cophenet): the
  most-contaminated non-tracer **Pat_13 p_dec=0.255** (no recovery); **Pat_14
  p_dec=0.785, Pat_15 p_dec=0.845** (epi-exclusion is *worse* than random);
  Pat_10 p_dec=0.120 (trend, trace stays ≈0); only **Pat_06 p_dec=0.025** but it is
  already a strong tracer. Cohort p_dec(β)=0.39. Pair-class (decimation-exempt) is
  mixed/patient-specific (for Pat_06 epi pairs carry MORE trace). **The −0.63 is
  correlational, not causal.**
- **Therefore the per-patient SCALAR + implant-covariate route is exhausted.** The
  PI's instruction is to go per-NODE instead (§5).

## 4. HARD CONSTRAINTS (the PI corrected me repeatedly — obey these)

1. **"It's the patient" is not a result.** Per-patient variance is publishable
   ONLY if explained by implant / epileptic placement / sampling. Otherwise it is
   noise the cohort statistic averages over.
2. **The BAND result stays.** "A trace appears more in some bands than others" is a
   real result — frame it better, do NOT dissolve it into a patient story.
3. **ρ^coph and Grassmann must stay DISTINCT.** If they tell the same story one is
   trivial. Use the complementarity (α pairwise / γ_l subspace / β both).
4. **Machinery-dependence (the α↔γ_l probe swap, window/aggregation fragility) is
   an INTERNAL fact. NEVER surface it in the paper.** Use it only to justify
   reading structure instead of the per-band bit.
5. **Prioritize the OLD results** when unifying: (a) the trace (β, gray-cortical,
   →OFC); (b) the learning/consolidation dissociation (task_learn vs task_test —
   encoding vs inference, N2 / `headlines/02_encoding_vs_inference.md`). New work
   serves these.
6. **Decimation control is MANDATORY** for any "remove subset X → trace changes"
   claim (`feedback_decimation_control_for_subset_exclusion`). Rebuilding LRG on a
   submatrix inflates ρ_split for ANY K-node removal. Pair-class / same-graph
   analyses are exempt but need their own random-pair null.
7. **Brutal honesty, no confidence laundering.** A control that kills a hypothesis
   is reported as such (as in §3).
8. **PI does not accept the negative result of §3** — believes a real explanator
   for the heterogeneity exists and the per-node route will find it.

## 5. THE NEXT TASK — per-node trace decomposition (PI's directive, verbatim intent)

> "see all the rho^coph_ij, then establish the nodes carrying the trace per each
> patient, and then figure out if the nodes that do not carry trace or are
> completely anti-trace have something in common." + "Recomputing traces when the
> explanation is removed should raise trace."

**Spec (seed for the scope report):**
1. **Per-pair contribution.** For each (patient, band) get the pairwise
   displacement vectors `Δ_task_ij` and `Δ_rest_ij` (cophenetic distance changes).
   Define a per-pair trace contribution `c_ij` = concordance of (Δ_task_ij,
   Δ_rest_ij) — e.g. the rank-product term entering the Spearman, or signed
   concordance `sign(Δ_task_ij)·sign(Δ_rest_ij)` weighted by magnitude. `ρ_split`
   is (proportional to) the mean of `c_ij`.
2. **Per-node trace score.** `t_i = aggregate_j c_ij` (mean over node i's pairs).
   Classify nodes per patient: **carrier** (t_i strongly +), **neutral**,
   **anti-trace** (t_i strongly −). Calibrate the threshold against a within-patient
   null (e.g. matched-strength surrogate per-node, or sign-flip).
3. **Pool across cohort** (~1170 nodes) → far more power than n=10 patients.
4. **Characterize the non-carrier / anti-trace nodes** by INDEPENDENT properties:
   DK region / a-priori system, lobe, gray vs WM, epileptic (SOZ) membership,
   hemisphere, **along-shaft depth** (contact_index_norm — a known confound, always
   include), node strength / hubness, spatial coordinates. Which property separates
   carriers from anti-trace nodes?
5. **Close the loop NON-circularly.** Removing nodes you selected *for being
   anti-trace* trivially raises the trace (selection bias) — that is NOT evidence.
   The evidence is: the anti-trace nodes share an a-priori property **P**; then test
   that excluding nodes **by P** (not by their trace score) raises the cohort trace
   **beyond random-node decimation** (`audit_85 --stratify <P> --mode node`) and
   ideally out-of-sample (leave-patient-out). Only then is P the explanation.
6. **Connect to established structure.** The β trace is already gray-cortical-
   dominant (gray↔gray pair-count hotspot, core-spared) and →OFC (audit_83). The
   per-node carriers should recapitulate that; the question is what the **anti**
   nodes are (the negative space audit_83 did not characterize).

**Reuse, don't reinvent:** `audit_83` already computes a per-node demeaned trace
incidence ("obs_trace") + matched-strength for the β localization — its per-node
machinery is the starting point; extend it to (a) all bands, (b) explicit
anti-trace classification, (c) anti-node characterization. The per-pair ρ^coph
visual tools exist: `preprint_24` (movers-density), `preprint_25` (consistency
t-map), `preprint_27` (3D). See memory `perpair_cophenetic_trace_visualization`.

## 6. Data & tooling map

| what | path |
|---|---|
| cophenet per-pt×band (obs_rho, obs_z, surr percentiles) | `data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv` |
| ρ_split generator (per-pair Δ recomputable here) | `scripts/01_compute/audit/audit_63_split_baseline_surrogate.py` |
| per-NODE β trace incidence + localization | `scripts/01_compute/audit/audit_83_localization_matched_strength.py`; `data/audit/localization_atlas/per_patient_decomp.csv` (β only) |
| Grassmann per-pt×band×k (obs_z) + windows | `data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv`; `data/audit/grassmann_cluster_extent/` |
| epi decimation + pair-class controls | `data/audit/epi_stratified/{decimation_control,pairclass_decimation}_per_patient.csv`; `audit_85_wm_decimation_control.py --stratify epi --mode {node,pairclass}` |
| region / lobe / system / hemisphere / x,y,z per contact | `lrg_eegfc.utils.io.regions.load_channel_regions(pat)`; `_normalise_label`; systems in `regions.py` (OFC, PFC, MTL, cingulate, insula, lateral_temporal, sensorimotor, parietal, occipital) |
| epileptic nodes (strings → normalise → intersect) | `lrg_eegfc.utils.io.patient.load_epileptic_nodes(pat)` |
| forest figure (two-axis) | `scripts/02_preprint/preprint_23_bands_perpair_rho_null_forest.py` |
| LRG cophenetic distance D_coph cache | `IMCOH_LRG_CACHE` (`config.paths`); halves at `data/cache/imcoh_lrg_halves/`, `imcoh_halves_fc/` |

Cohort n=10: Pat_02,03,05,06,07,08,10,13,14,15. Drop Pat_10 task rows [53,54,55].

## 7. Artifacts produced this session

- **Report** `.agents/reports/2026-06-25_cophenetic-gate-presence-vs-consistency.md`
  — presence-vs-consistency gate diagnosis + two-axis fix + cross-probe table.
  **NEEDS REVISION**: its §6 Pat_06 "implant coverage" explanation is REFUTED (§3
  here); update once the per-node analysis lands. Its gate/cascade content (§1-§5,
  §8-§9) stands.
- **Figure** `preprint_23` two-axis verdict; PDFs regenerated.
- **No locked-ledger edits made.** The "execute full cascade" the PI approved
  earlier is **on hold** — superseded by "fix the storytale first" then this
  per-node redirect. Do not cascade until the heterogeneity is explained and the
  framing is ratified.

## 8. Open decisions (deferred, pending the per-node result)

- **Band framing** (PI deferred until the heterogeneity control resolves): β =
  complete trace (both probes); α / γ_l = partial (one probe each). Lock after §5.
- **γ_l verdict language**: "present but split", "subgroup trace", or
  "heterogeneous" — naming unresolved.
- Whether the per-node anti-trace characterization becomes a manuscript control
  (likely yes) or a methods note.
