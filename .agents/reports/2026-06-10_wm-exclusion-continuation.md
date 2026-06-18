---
name: wm-exclusion-continuation
type: plan
era: IMCOH_ABS_COHORT_N10
status: current
created: 2026-06-10
updated: 2026-06-10
scope: post-compaction continuation plan for the white-matter-exclusion control (C6) — the two honesty-flag fixes + the full file map to resume work
pointers:
  - .agents/reports/2026-06-08_white-matter-exclusion.md
  - .agents/preprint/responses/2026-06-08_wm-exclusion-cascade.md
  - .agents/preprint/locked/CONTROLS.md
  - memory/white_matter_exclusion_2026_06_08.md
---

# WM-exclusion (C6) — post-compaction continuation plan

> **✅ 2026-06-12 — BOTH FIXES DONE + RE-CASCADED.** `audit_85` (Fix A,
> random-node decimation) and `audit_86` (Fix B, C3 cluster-extent gate on the WM
> submatrix) are built, run, and cascaded. **Key outcome:** Fix A produced a
> genuine correction — the cophenetic "sharpening" is **generic node-count, NOT
> WM-specific** (cohort `p_dec`: β 0.145, α 0.305, γ_l 0.185), so the "WM removal
> sharpens / gray-resident" framing was **withdrawn** for the cophenetic
> substrate; the C6 claim is now **survival/robustness**. Raw β/γ_l ARE WM-specific
> (`p_dec` 0.000/0.030). Fix B: β/γ_l Grassmann **re-pass the locked C3 gate**
> (`cluster_p_mass^wmX`=0.005, strong, LOO-robust), δ fails (0.105). No verdict
> flips. Cascaded into CONTROLS §C6, 4 briefs, VERDICT_LEDGER (2026-06-12
> amendment), responses (Resolved section), report (§4/§5), memory. **Remaining
> open:** PTD<0 companion (§5); NOT committed. The §2/§3 below are the original
> as-built spec, retained for provenance.

**Head.** The white-matter-exclusion analysis is **done and wrapped into the
preprint as control C6** (secondary sensitivity layer, no verdict changes). Two
honesty flags remain and require **new compute** (not text edits): (A) the
~40 %-node-count confound on the `exclude_wm` `weaken`/`emerge` calls, and (B)
C6-Grassmann was judged on the lighter per-`k` Wilcoxon basis, not the locked C3
cluster-extent mass gate. This doc specifies both fixes precisely (so they're
mechanical to build after compaction) and lists every file needed to resume.

**Resume command after compaction:** read this file, then build **audit_85** (Fix
A) and **audit_86** (Fix B) per §2, run them, then re-cascade per §3.

---

## 1. Current state (DONE)

- **Analysis:** `audit_83` (WM-stratified cophenetic + raw ρ_split), `audit_84`
  (WM-stratified Grassmann T_G(k)), `audit_83b` (figures). WM = dominant
  Desikan-Killiany tissue `Wm` (`load_channel_regions(pat)["region"]=="Wm"`),
  30–57 % of montage (median 117→74 nodes under `exclude_wm`). n=10, 6 bands,
  matched-strength R=200 regenerated per submatrix.
- **Result (headline):** WM removal does NOT weaken the trace, it sharpens it.
  Cophenet α persist+strengthen, β preserved/strengthened (p=0.053 marginal but
  LO-P15 0.037 / gray_gray 0.032 / raw-β 0.024), γ_l+θ emerge; raw substrate
  clears at δ/α/β/γ_l once WM gone; Grassmann β/γ_l/γ_h per-`k` counts unchanged
  (β 40→40); only δ Grassmann weakens (23→7, = the epi-biology δ channel).
  Trace is gray-matter-RESIDENT (defensive result).
- **Report:** `.agents/reports/2026-06-08_white-matter-exclusion.md` (canonical).
- **Preprint wrap (C6):** `locked/CONTROLS.md` §C6 + head; `locked/VERDICT_LEDGER.md`
  2026-06-08 entry; `bands/{01_beta,02_alpha,03_gammalow,06_delta}.md` frontmatter
  `verdict_layers.{rho_split_coph_wm_excluded, grassmann_wm_excluded}` + revision
  lines; `responses/2026-06-08_wm-exclusion-cascade.md`. No verdict changes.
- **NOT committed.** (Date note: the cascade text is stamped 2026-06-08 — when the
  compute ran — but was written 2026-06-10; cosmetic, optional to relabel.)

---

## 2. The two fixes (FUTURE MOVES — new compute)

### Fix A — node-count confound → `audit_85_wm_decimation_control.py`

**Why.** `exclude_wm` drops ~40 % of nodes, so full-vs-exclude mixes biology with
power. A `weaken` is ambiguous (WM-carried vs noisier Spearman); an `emerge`
could be a cleaner gray-only renormalization regardless of which nodes go.

**Design (size-matched random-node decimation; cheap, NO nested matched-strength).**
For each (patient, band): let `K = #WM nodes`. For `r in 1..R_dec` (R_dec=200),
draw a uniformly-random node-keep mask of size `N−K` (drop K random nodes, the
removed set need NOT be WM), submatrix all 4 phases, recompute `D_coph`
(`lrg_ultrametric_condensed`) and raw (`raw_condensed`), and `ρ_split` (cophenetic
+ raw) on each random subgraph → a null distribution of the trace statistic under
random K-node loss. Compare the **observed `exclude_wm` ρ_split** to that
distribution.
- One-sided percentile p_dec = fraction of random subgraphs with
  `ρ_split_random ≥ ρ_split_exclude_wm`.
- **Reading:** `exclude_wm` ≈ random-decimation median → the strengthen/emerge is
  a GENERIC node-reduction effect (or pure power), NOT WM-specific. `exclude_wm`
  ABOVE the random distribution → WM removal SPECIFICALLY helps (WM was diluting).
  BELOW → WM was carrying.
- Cohort: median p_dec per (band, substrate) + n patients where `exclude_wm`
  beats its own random-decimation median; report the headline bands α/β/γ_l/θ.
- **Reuse (no forks):** `load_phase_fc`, `lrg_ultrametric_condensed` (audit_63);
  `raw_condensed`, the `ρ_split` Spearman (audit_83); `wm_mask_for`,
  `cell_rng`-style deterministic RNG (`_wm_stratify`). Per-cell RNG seeded by
  (pat,band) for reproducibility.
- **Output:** `data/audit/wm_stratified/decimation_control_per_patient.csv`
  (pat×band×substrate: obs_exclude_wm, dec_median, dec_p5/p50/p95, p_dec) +
  `decimation_control_cohort.csv` (band×substrate: median p_dec, n_above_dec_median,
  verdict {WM_specific | generic_nodecount | WM_carried}). README.
- **Cost:** R_dec=200 × 4 phases × ~60 cells eigh of ~74-node matrices ≈ minutes.
- **OPTIONAL heavier follow-up (only if a borderline cell needs it):**
  matched-strength-nested decimation (each random subgraph gets its own R=50
  matched-strength null) to compare the strength-corrected z, not just raw ρ.
  Skip unless Fix A's first pass leaves β/γ_l genuinely ambiguous.

### Fix B — C6-Grassmann on the real C3 gate → `audit_86_wm_grassmann_cluster_extent.py`

**Why.** C6-Grassmann reported per-`k` Wilcoxon significant-cell counts (β 40→40)
— the audit_66 basis, NOT the locked C3/C5 Grassmann gate (audit_70 cluster-extent
mass `cluster_p_mass`). So "β Grassmann survives WM-X" is currently a sensitivity
read, not a re-pass of the locked gate.

**Design.** Run the **audit_70 cluster-extent permutation mass test** on the
`exclude_wm` (and optionally `wm_only`) per-`k` T_G curves → `cluster_p_mass^wmX(b)`;
compare to the C3 gate (`< 0.05`). The cluster-extent null is a phantom-surrogate
permutation (each of R surrogates treated as obs vs the remaining R−1), so it
needs the **per-surrogate** T_G(k) stack, which audit_84 did NOT persist (only
summary percentiles). Regenerate it cheaply from the cached `wmX` eigvecs:
- Load `data/cache/matched_strength_surrogate_wm_excluded_lrg/Pat_NN/{band}_{phase}_wmX_R200_swap20_seed20260608_imcoh_abs.npz` (3 phases: rest_pre_A, task_test, rest_post).
- Per surrogate r: `topk_basis(evecs[r], kmax)` → `t_g_at_all_k` (audit_67/84) → surr T_G(k).
- Obs T_G(k): `laplacian_eig(Wk)` per phase → `topk_basis` → `t_g_at_all_k`.
- Per-k cohort Wilcoxon p (`audit_70.wilcoxon_per_k_greater`) → `audit_70.cluster_mass(p, ALPHA_K)` for obs and each phantom-surrogate → `cluster_p_mass^wmX = (1+#{null≥obs})/(R+1)`. LOO max diagnostic per audit_70.
- **Reuse (no forks):** `audit_70.{cluster_mass, wilcoxon_per_k_greater}`;
  `audit_67.{laplacian_eig, topk_basis, t_g_at_all_k, K_GRID, load_phase_fc_full}`;
  `_wm_stratify.{wm_mask_for, node_mask_for_config, surr_eig_path}`;
  `load_or_compute_eigs_at_path` (caches warm).
- **Output:** `data/audit/wm_stratified/grassmann_cluster_extent_wmX.csv`
  (band: obs_cluster_mass, cluster_p_mass^wmX, LOO max, verdict vs C3 gate).
- **Expectation:** β/γ_l likely hold `cluster_p_mass^wmX < 0.05` (per-`k` counts
  were unchanged) — but this must be MEASURED, not assumed. δ likely fails (it
  weakened on per-`k` already).
- **Cost:** caches warm → minutes.

---

## 3. Re-cascade after fixes A+B (replace caveats with results)

Once audit_85/86 land, update — same files as the original C6 cascade:
- **`locked/CONTROLS.md` §C6:** replace the "node-count caveat" bullet with the
  decimation-control verdict (per band: WM-specific vs generic); replace the
  "C6-Grassmann uses per-`k` Wilcoxon, not C3 cluster-extent" caveat with
  `cluster_p_mass^wmX` per band (now judged on the real gate). Keep C6 SECONDARY.
- **`bands/{01,02,03,06}.md`:** update `verdict_layers.{rho_split_coph_wm_excluded,
  grassmann_wm_excluded}` — swap the caveat phrasing for the decimation p_dec and
  `cluster_p_mass^wmX`. Add a 2026-06-10 revision line each.
- **`locked/VERDICT_LEDGER.md`:** dated 2026-06-10 amendment to the 2026-06-08
  entry (flags resolved; still no verdict change unless a fix surprises us — if
  Fix A shows an `emerge` is generic-node-count, DEMOTE that emerge to "not
  WM-specific" in C6 text, still no locked-verdict change since C6 is secondary).
- **`responses/2026-06-08_wm-exclusion-cascade.md`:** update §"Honest caveats" →
  "Resolved" with the two results.
- **Report `2026-06-08_white-matter-exclusion.md` §4-5:** replace the caveat prose
  with the decimation + cluster-extent outcomes.
- **Memory `white_matter_exclusion_2026_06_08.md` + MEMORY.md line:** note flags
  resolved.

---

## 4. File map — everything to resume

**Build (scripts):**
- `scripts/01_compute/audit/_wm_stratify.py` — WM mask + config vocab + cache paths + cohort_verdict (reused from `_epi_stratify`).
- `scripts/01_compute/audit/audit_83_wm_stratified_cophenetic.py` — cophenetic + raw (defines `raw_condensed`, `_surr_stacks`, `_stats_row`).
- `scripts/01_compute/audit/audit_84_wm_stratified_grassmann.py` — Grassmann (re-emit full from audit_66 + fresh exclude_wm/wm_only).
- `scripts/01_compute/audit/audit_83b_wm_stratified_figures.py` — figures.
- **Reused engines:** `audit_63_split_baseline_surrogate.py` (`load_phase_fc`, `lrg_ultrametric_condensed`, `ensure_half_fcs`); `audit_67_grassmann_epi_exclusion.py` (`laplacian_eig`, `topk_basis`, `t_g_at_all_k`, `K_GRID`, `load_phase_fc_full`); `audit_70_grassmann_cluster_extent.py` (`cluster_mass`, `wilcoxon_per_k_greater` — Fix B); `src/lrg_eegfc/utils/surrogate/matched_strength.py` (`load_or_compute_eigs_at_path`, `cophenetic_condensed_from_eigs`, `adjacency_from_laplacian_eigs`).

**Data:**
- `data/audit/wm_stratified/{cophenetic_raw_per_patient,cophenetic_raw_cohort,grassmann_per_patient_per_k,grassmann_cohort}.csv` + `README.md`, `README_grassmann.md`.
- `data/audit/wm_stratified/figures/fig_wm_{verdict_matrix_coph,verdict_matrix_raw,exclude_forest,grassmann_kspan}.pdf`.
- Surrogate caches: `data/cache/matched_strength_surrogate_wm_{excluded,only}_lrg/Pat_NN/{band}_{phase}_wm{X,ONLY}_R200_swap20_seed20260608_imcoh_abs.npz` (warm — Fix B reads `wmX` eigvecs).
- Canonical full-graph cache (full/pair-class): `data/cache/matched_strength_surrogate_lrg/...seed20260511...` (warm).

**Preprint (the C6 wrap):**
- `.agents/preprint/locked/CONTROLS.md` §C6 (+ head, "## The sensitivity layers").
- `.agents/preprint/locked/VERDICT_LEDGER.md` — 2026-06-08 revision entry.
- `.agents/preprint/bands/01_beta.md`, `02_alpha.md`, `03_gammalow.md`, `06_delta.md` — `verdict_layers.*_wm_excluded` + revision lines.
- `.agents/preprint/responses/2026-06-08_wm-exclusion-cascade.md`.
- Routing rules (if more cascade needed): `.agents/preprint/WRITING_GUIDE.md`.

**Report + memory:**
- `.agents/reports/2026-06-08_white-matter-exclusion.md`.
- `memory/white_matter_exclusion_2026_06_08.md` + MEMORY.md index line.

**Env:** run with `/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python`
(direct binary — `conda run` swallows tracebacks and emits a benign cross-compiler
warning). Filter stderr with `grep -vE "cross-compiler|x86_64-conda"`.

---

## 5. Other open items (lower priority)

- **PTD<0 robustness companion** — alternate WM definition (proximal-tissue-density
  < 0, Mercier 2017) vs the atlas-dominant `Wm` used. Conclusions that survive a
  40 % cut are unlikely to flip under a borderline reclassification, but untested.
- **Commit** — nothing committed yet (branch `audit/cohort-n10-diagnostic`).
- **Date stamps** — cascade text says 2026-06-08 (compute date); actually written
  2026-06-10. Cosmetic; relabel only if doing a commit pass.
