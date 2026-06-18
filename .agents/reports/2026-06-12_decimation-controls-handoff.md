---
name: decimation-controls-handoff
type: plan
era: IMCOH_ABS_COHORT_N10
status: current
created: 2026-06-12
scope: post-compaction handoff — node-count decimation control on C5 epi-X + C6 WM-X (done), plus the one OPEN item (pair-level random-pair-subset control) and the PI-review decisions
supersedes: .agents/reports/2026-06-10_wm-exclusion-continuation.md (WM-only; this extends it to epi-X + pair-class)
pointers:
  - data/audit/wm_stratified/decimation_control_cohort.csv
  - data/audit/epi_stratified/decimation_control_cohort.csv
  - memory/feedback_decimation_control_for_subset_exclusion.md
  - .agents/preprint/locked/VERDICT_LEDGER.md   # 2026-06-12 amendments (WM + epi)
---

# Decimation-control session — handoff (2026-06-12)

**Resume command:** read this file, then act on §3 (the open item is the
pair-level random-pair-subset control; the rest are decisions for the user).

## 1. What was established (the science)

The node-count confound: rebuilding the LRG on a node submatrix inflates cophenetic
ρ_split for **any** ~K-node removal, so "exclude subset X → trace strengthens/
emerges" is confounded with node count. Control = `audit_85_wm_decimation_control.py
--stratify {wm,epi}` (drop K random nodes ×200, compare obs `exclude_X` to the
random-decimation distribution; cohort `p_dec`).

- **Cophenetic strengthening/emergence is GENERIC NODE-COUNT for BOTH exclusions**
  (none `p_dec < 0.05`): WM-X α 0.305 / β 0.145 / γ_l 0.185; epi-X α 0.135 /
  β 0.390 / γ_l 0.130 (epi drops only ~9.5% of nodes, yet random removal
  reproduces it). The "strengthens" / "gray-resident" / "epileptic zone masks it"
  framings are **withdrawn for the cophenetic substrate**.
- **Genuinely tissue-specific:** raw |ImCoh| WM β/γ_l (`p_dec` 0.000/0.030,
  `WM_specific`). epi has NO tissue-specific trace cell (raw δ/θ are `epi_carried`
  = epi channel).
- **Grassmann (audit_86, WM-X on the locked C3 cluster-extent gate):** β/γ_l
  RE-PASS (`cluster_p_mass^wmX`=0.005, strong, LOO-robust); δ FAILS (0.105);
  γ_h emerges but LOO-fragile (0.050). Survival, not node-count.
- **Consequence for C5 (PRIMARY lens):** the 2026-06-05 PI elevation rested on the
  (now node-count) subgraph strengthening → re-anchored on the PAIR-CLASS
  decomposition in the records; flagged ⚠ for PI review. γ_l lost its only cophenet
  evidence (verdict stays "strong, only Grassmann").

## 2. What was cascaded (all done, consistent)

- Scripts: `audit_85_wm_decimation_control.py` (generalized `--stratify {wm,epi}`),
  `audit_86_wm_grassmann_cluster_extent.py`.
- Data: `data/audit/{wm,epi}_stratified/decimation_control_{cohort,per_patient}.csv`
  + `README_decimation.md`; `data/audit/wm_stratified/grassmann_cluster_extent_wmX.csv`.
- Preprint: `CONTROLS.md` §C5 (node-count bullet) + §C6 (rewritten); briefs
  `01_beta`/`02_alpha`/`03_gammalow` (epi entries + 2026-06-12 lines), `06_delta`
  (WM); `VERDICT_LEDGER.md` (two 2026-06-12 amendments: WM, then epi);
  `responses/2026-06-08_wm-exclusion-cascade.md` (Resolved); report
  `2026-06-08_white-matter-exclusion.md` (§4/§5).
- Rule locked: `.agents/guides/04_rules/never-always-list.md` +
  `memory/feedback_decimation_control_for_subset_exclusion.md`.
- Memory: `white_matter_exclusion_2026_06_08`, `epi_stratified_trace_2026_06_05`,
  MEMORY.md lines updated.

## 3. OPEN (pick up here)

1. **Pair-class specificity — DONE 2026-06-12 (`audit_85 --mode pairclass`).**
   Built the random-same-size-PAIR-subset null: draw M = #class-pairs random pairs
   from the FULL-graph cophenetic/raw distances R=1000× and locate the observed
   class ρ_split. `p_pair<0.05`=**concentrated** (genuine localization),
   `>0.95`=**depleted** (spared/anti class), else=**generic_paircount**
   (underpowered, indistinguishable from a random M-pair subset). Orthogonal to
   matched-strength. Outputs `data/audit/{wm,epi}_stratified/
   pairclass_decimation_{cohort,per_patient}.csv` + `README_pairclass.md`.
   **⚠ Runs at the CURRENT n=10 (Pat_15 now has 17 epi nodes from the 2026-06-11
   SOZ labelling), so it SUPERSEDES the Jun-5 epi pair-class CSV (n=9):** observed
   epi_epi α +0.423 (was cascaded +0.404), epi_epi β +0.142 (was +0.240).
   **RESOLVED 2026-06-16: audit_77 epi-stratified REGENERATED at n=10** (Jun-5 n=9
   CSVs backed up to `cophenetic_{cohort,per_patient}_n9_2026-06-05.bak.csv`). NO
   sensitivity-flag flips on any load-bearing cell — α cross/epi_epi/nonepi `persist`
   (epi_epi α STRENGTHENS, p 0.006→0.003), β cross/nonepi `persist`, **β epi_epi
   stays `weaken`** (obs +0.240→+0.142). Only 2 immaterial verdict-bucket relabels
   on γ null cells. Observed ρ_split matches pair-decimation `cohort_obs_class`
   bit-for-bit. **The two nulls now AGREE on the epi core:** α `persist`+`concentrated`,
   β `weaken`+`depleted` → "α recruits / β spares the core" holds on both.
   **Verdicts that decide the re-anchor (focus on small/specific classes — big
   classes like nonepi_nonepi @78% give near-deterministic verdicts that just
   track the graph average):**
   - **α-exception CONFIRMED:** epi_epi α cophenetic obs +0.423 vs null +0.129,
     p_pair<0.001 → **concentrated**. α genuinely recruits the diseased core,
     not a pair-count fluctuation. (Raw not concentrated → cophenetic/hierarchical,
     consistent with α being cophenet-only.) cross α also concentrated.
   - **"β spares the core" CORROBORATED (reverses the earlier "soften to
     underpowered" guess):** epi_epi β **depleted** on raw (obs −0.153 vs +0.235,
     p_pair=1.0), trends depleted on cophenetic (p=0.945). Core carries LESS trace
     than random pairs. β trace is gray-concentrated (gray_gray concentrated both
     substrates) and epi-core/WM-depleted.
   - **WM↔WM β: matched-strength-survival ≠ pair-concentration.** wm_wm β is
     `separated` under matched-strength (not a strength artifact) yet **depleted**
     vs random pairs both substrates (p_pair=1.0). The two nulls are orthogonal:
     β trace is gray-DOMINANT, WM-WM a below-average (strength-clean) contributor.
     Confirms "gray-dominant-not-exclusive" + adds the unconfounded "WM-WM depleted".
2. **PI decisions (user's call — NOW INFORMED by item 1, flagged ⚠ in C5 + briefs):**
   (a) the pair-class CAN carry a re-anchored C5 lens, but the headline is a **band
   dissociation**, not one "epi-X strengthens" lens: **α RECRUITS the core**
   (epi_epi concentrated), **β SPARES it** (epi_epi depleted). Decide: keep C5
   *primary* re-anchored on this α-recruits/β-spares pair-class dissociation, or
   demote to *secondary* like C6. (b) "β spares the diseased core" → **CORROBORATE**
   (raw depleted p=1.0; cophenet p=0.945; matched-strength epi_epi β `weaken`), do
   NOT soften to "underpowered". (c) ✅ DONE — audit_77 regenerated at n=10, no flag
   flips (see item 1). **(a) still open: re-anchor framing held for PI discussion.**
3. **Commit:** nothing committed; branch `audit/cohort-n10-diagnostic`. Three new
   scripts (audit_85 +node+pairclass modes, audit_86) + 2 data dirs + the cascade +
   memory are all uncommitted.
4. **Naming collision (flagged):** `audit_85`/`audit_86` collide with the
   2026-06-05 `audit_85_epi_propagator_recovery.py` /
   `audit_86_epi_propagator_matched_strength.py` (consistent with the project's
   per-track suffix-disambiguated convention — audit_83/84/90/91/92/99 all double).
   Renumber to globally-unique IDs only if the user wants.

## 4. Env

Run with `/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python` (direct
binary — `conda run` swallows tracebacks). Filter stderr with
`grep -vE "cross-compiler|x86_64-conda"`.
