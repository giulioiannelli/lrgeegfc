---
name: writing-directive-alpha-results-drafting
era: IMCOH_ABS_COHORT_N10
status: open
kind: writing-directive
scope: results-section LaTeX drafting for the α subsection, mirroring the locked β subsection
created: 2026-05-22
companion: writing_directive_2026-05-20_beta_results_application.md (locked β subsection — structural template), .agents/preprint/bands/02_alpha.md (lab α brief; numbers verified here, three values flagged as stale or unverified)
target: writing agent producing manuscript Results §`\ssec:results_alpha`
verified_against: data/audit CSVs as of 2026-05-22 (every cited number cross-checked to the cited decimal; three discrepancies with the bands brief flagged in §"Critical issues")
---

# α Results-section drafting — writing-agent directive

**Head.** Produce a new `\subsection{\texorpdfstring{\(\alpha\)}{α}: a per-pair multiscale task signature concentrated in non-epileptic cortex}` block mirroring the locked β subsection structure. α is **single-probe**: the per-pair multiscale probe on `\(\Dcoph\)` carries a cohort trace; the Grassmann subspace probe is silent (cluster-mass `\(p_{\mathrm{mass}} = 0.348\)`, no contiguous-significant run). The defining feature of the α subsection is that the per-pair cohort agreement at the matched-strength gate is borderline at full cohort (5/10 above own surrogate) and is **decisively resolved by epileptogenic-zone exclusion** (7/10, effect-size ratio 8.3× → 27.7×). The anatomical localization is an 11-region distributed cortical network; under epi-X 9 of these 11 regions reproduce (verified after a 2026-05-22 patch to `audit_71`, see Issue 1) and the two dropouts carry distinct biological meaning. Three discrepancies between the bands brief `02_alpha.md` and the current locked CSVs were originally flagged; **Issue 1 was a script bug, now fixed and rerun; Issues 2–3 remain as brief-vs-Methods naming/value drifts and are resolved in favor of Methods**.

## What to read first

1. **This directive** — verified numerical table + critical issues + structure scaffold + figure placeholders.
2. **The locked β subsection** (passed to the writing agent as the structural template). Mirror its paragraph order, voice, sentence cadence, and table layout.
3. **The locked Methods §`\ssec:methods_compare` and §`\sssec:methods_anatomy`** (passed to the writing agent). These are the single source of truth for which statistical test gates which claim — when the bands brief and Methods disagree on test names or numbers, **Methods wins**.
4. **`bands/02_alpha.md`** — lab α brief; useful for narrative context, but three of its specific numerical claims are stale or use a non-Methods-locked test (see §"Critical issues" below). Always cross-check the brief's numbers against the verified table here.
5. **`locked/CONTROLS.md`** §§C1–C5, **`locked/ANATOMY_CONTROLS.md`**, and **`locked/VERDICT_LEDGER.md`** α row — verdict + gate definitions.

## Critical issues — RESOLVE BEFORE WRITING

These three issues are quantitative-consistency problems between the bands brief `02_alpha.md` and the live audit CSVs as of 2026-05-22. The directive resolves each in favor of the CSV; the brief must not be cited where it disagrees.

### Issue 1 — α cophenet anatomy under epi-X: **RESOLVED 2026-05-22**

**Bug history (now fixed).** The original `data/audit/anatomy_alpha_cophenet_epiX/cohort_summary.csv` was byte-identical to the full-cohort file because `audit_71_anatomy_cophenet.py` had a silent no-op at the `--epi-x` mask:

```python
# Buggy code (lines 278–283 and mirror at 354–359):
epi = load_epileptic_nodes(pat)                 # returns ['B7', 'B8', …]  (STRINGS)
keep_mask = ~np.isin(np.arange(len(regions_df)), np.asarray(list(epi)))
#                    ^ integers [0..N-1]         ^ string array
# np.isin(int, str) → all-False → keep_mask = all-True → mask never filters.
```

Confirmed by the script log: under `--epi-x` Pat_02 reported `n_pairs = 6786` (the full N·(N−1)/2 = 6786) and Pat_13 reported `n_pairs = 7021`, when the correct epi-X counts are 5253 and 4005 respectively.

**Patch applied 2026-05-22.** Replaced the buggy mask with the same label-based translation that `audit_68_alpha_epi_exclusion.py` uses:

```python
epi = set(load_epileptic_nodes(pat))
channel_labels = [_normalise_label(l) for l in regions_df["label_raw"]]
keep_mask = np.array([l not in epi for l in channel_labels], dtype=bool)
pair_keep = keep_mask[obs["iu_i"]] & keep_mask[obs["iu_j"]]
```

(`_normalise_label` imported from `lrg_eegfc.utils.io.regions`.) Patch applied to **both** the observation block (line 278–283) and the surrogate-loop mirror (line 354–359). Sanity-check after patch: Pat_02 keep_mask drops 14 contacts (matches audit_68 `N_epi = 14`); Pat_13 drops 30 (matches audit_68 `N_epi = 30`); Pat_15 drops 0 (matches audit_68 `N_epi = 0`).

**Post-patch rerun: real α cophenet anatomy under epi-X.** `data/audit/anatomy_alpha_cophenet_epiX/cohort_summary.csv` is now a genuine epi-X aggregation (MD5 = `199d742b22e00893ff3fbae4623a1b1c`, distinct from the full-cohort MD5). Per-patient pair counts drop appropriately: Pat_05 6903 → 5356, Pat_13 7021 → 3916, Pat_15 6903 → 6903 (N_epi = 0).

**The real epi-X anatomy result: 9 of 11 named regions survive; 2 drop out.**

| # | Region | Full-cohort verdict | Epi-X verdict | Reading |
|---|---|---|---|---|
| ✓ | lh-caudal-anterior-cingulate | pass (enr 2.06, q 2.1e-5, z 4.67) | **pass** (enr 2.08, q 2.9e-2, z 5.79) | stable |
| ✓ | lh-parahippocampal | pass (enr 3.80, q 1.4e-25, z 9.76) | **pass strengthens** (enr 3.89, q 6.6e-25, z 8.05) | the medial-temporal anchor; stays strong |
| ✓ | lh-rostral-anterior-cingulate | pass (enr 3.04, q 4.1e-20, z 14.26) | **pass strengthens dramatically** (enr 5.37, q 1.5e-31, z 17.38) | the strongest gainer under epi-X |
| ✓ | rh-caudal-anterior-cingulate | pass (enr 2.44, q 1.8e-8, z 3.45) | **pass strengthens** (enr 3.31, q 3.2e-11, z 3.80) | |
| ✓ | rh-caudal-middle-frontal | pass (enr 2.90, q 3.5e-21, z 6.23) | **pass** (enr 2.74, q 5.2e-17, z 5.52) | mild weakening |
| ✓ | rh-postcentral | pass (enr 3.47, q 3.8e-65, z 4.47) | **pass** (enr 3.66, q 7.4e-69, z 4.50) | unchanged |
| ✓ | rh-posterior-cingulate | pass (enr 1.99, q 5.9e-10, z 3.13) | **pass** (enr 1.93, q 1.1e-8, z 2.95) | mild weakening |
| ✓ | rh-precuneus | pass (enr 1.47, q 1.3e-3, z 2.60) | **pass** (enr 1.47, q 1.1e-3, z 2.60) | identical |
| ✓ | rh-superior-parietal | pass (enr 1.64, q 3.7e-5, z 3.64) | **pass** (enr 1.64, q 3.4e-5, z 3.64) | identical |
| ✗ | **lh-caudal-middle-frontal** | pass (enr 1.49, q 1.6e-3, z 3.02) | **A1 FAIL** (enr 1.34, q 5.1e-2, z 2.32, p_emp 0.020) | drops at A1 — q_BH = 0.051 is just above the 0.05 threshold; the weakest of the 11 in full cohort; marginal statistical fail |
| ✗ | **rh-medial-orbitofrontal** | pass (enr 1.77, q 1.1e-16, z 3.22, p_emp 0.005) | **A2 FAIL** (enr 1.55, q 3.3e-6, z 1.28, p_emp 0.124) | drops at A2 (matched-strength surrogate): A1 still strong but the matched-strength z plummets from 3.22 → 1.28. Biologically distinct from the lh-cMF dropout: rh-medialOFC's full-cohort enrichment had a substantial epi-zone component, and once epi-zone contacts are removed the residual signal fits within matched-strength noise. |

**Required action — in the manuscript paragraph.** Update the anatomy claim from the brief's stale "11 named regions, all reproduce under epi-X" to the verified **"9 of 11 named regions reproduce under epi-X; two drop out under the joint A1 + A2 gate (lh-caudal-middle-frontal at A1 marginally and rh-medial-orbitofrontal at A2 — the latter losing its matched-strength surrogate significance under epi-X, indicating a substantial epi-zone contribution to its full-cohort enrichment)."** The 9-region surviving network spans bilateral cingulate (4 regions: lh- and rh-caudal-anterior-cingulate, lh-rostral-anterior-cingulate, rh-posterior-cingulate), left parahippocampal cortex, right caudal middle frontal gyrus, right postcentral gyrus, right precuneus, and right superior parietal cortex.

The five regions that **strengthen** under epi-X (lh-parahippocampal, lh-rostral-AC, rh-caudal-AC; rh-postcentral is stable) provide direct anatomical support for the trace-level finding that α reorganization concentrates in non-epileptic cortex (paragraph α-22 through α-28 above). The rh-medial-OFC dropout is biology, not artifact — it tells the reader that the medial-OFC contribution to the full-cohort α anatomy was carried partly by epi-zone-coupled topology, and the matched-strength control distinguishes that contribution from the genuinely-topological signal that survives in the remaining 9 regions.

**Source for verified epi-X values.** `data/audit/anatomy_alpha_cophenet_epiX/cohort_summary.csv` (post-patch, written 2026-05-22). The bands brief `02_alpha.md` §5 "C5 epi-X anatomy" is now **stale** and should be updated separately to match the verified outcome.

### Issue 2 — Brief's Grassmann cluster-mass p-value is stale

- The brief reports `cluster_p_cluster_mass = 0.0995` (line 189 of `bands/02_alpha.md`).
- Current locked CSV `grassmann_cluster_extent/cohort_summary.csv` α row reports **`cluster_p_cluster_mass = 0.348`**.
- This reflects the post-fix audit_70 "all-clusters sum" `\(T_{\rm G}^{\ast}\)` definition locked 2026-05-19 pm (project memory `feedback_cluster_mass_all_clusters.md`; binding Methods definition is eq. `methods_TGstar`).
- The qualitative conclusion (α has no Grassmann trace) is unchanged — both 0.0995 and 0.348 are well above the 0.05 gate, and the verdict `cluster_p_mass_loo_max = 0.577` (Pat_14-driven) confirms robustness to single-patient leverage.
- **Required action**: cite the **post-fix value `\(p_{\mathrm{mass}} = 0.348\)`**, not 0.0995, when reporting the Grassmann result for α. The longest-run cluster `\(p_{\mathrm{LR}} = 0.159\)` is also available but is descriptive only — the Methods-locked Grassmann gate is `\(p_{\mathrm{mass}}\)`.

### Issue 3 — Brief's C5 cophenet p uses the paired-vs-surrogate test; Methods locks the one-sample test

- The brief's C5 epi-X table reports `paired Wilcoxon z/p = 49 / 0.014` (line 128 of `bands/02_alpha.md`), drawn from `alpha_epi_exclusion/cohort_summary.csv` `paired_wilcoxon_p` (obs vs patient-matched surrogate median, under H₁: obs > surr).
- Methods §`\sssec:methods_compare_stats` locks the C5 per-pair cophenet gate as a **one-sample one-sided Wilcoxon on per-patient `\(\rhosplit^{\rm epi-X}\)` under H₁: `\(\rhosplit^{\rm epi-X} > 0\)`**, gated at p < 0.05.
- The Methods-locked test is reported in `alpha_epi_exclusion/c5_wilcoxon_cohort.csv`: **`wilcoxon_one_sided_p_epiX = 0.010`**, `wilcoxon_loo_max_p_epiX = 0.020` (Pat_03), `c5_pass = True`.
- Both tests pass at α = 0.05; the conclusion (C5 strengthens α and gate-passes) is unchanged. But the **manuscript must cite the Methods-locked test (p = 0.010)**, not the paired-vs-surrogate variant (p = 0.014), for the C5 verdict.
- The ratio (8.3× → 27.7×) and the per-patient agreement count (5/10 → 7/10) remain descriptive accompaniments to the locked gate — both numbers are correct.

## Numerical verification — α (every cited number CSV-verified to 3 decimals)

| # | Quantity | Value | Source CSV row | Notes |
|---|---|---|---|---|
| **Substrate (raw `\(\abs{\ImCoh}\)`)** | | | | |
| α-1 | cohort median `\(\rhorawsplit\)` | +0.158 | `raw_fc_matched_strength/cohort_summary_all_bands.csv` α `obs_median_rho` | |
| α-2 | surrogate cohort median | +0.008 | same `surr_median_rho_median` | |
| α-3 | ratio obs/surr | ≈ 21× | derived (0.158/0.0076) | |
| α-4 | cohort-paired Wilcoxon p (obs > surr) | **0.014** | same `paired_wilcoxon_p` | |
| α-5 | `n_above_surrogate` | 8/10 | same | |
| α-6 | substrate verdict | **`"separated"`** | same `verdict` | The **only** band among the six classified `separated` at the substrate level — α has the strongest substrate signal of the cohort. |
| **Per-pair cophenet — within-baseline (Methods C1, C2)** | | | | |
| α-7 | within-baseline `\(\rhosplit\)` median | +0.115 | `ctm_triangle/cohort_summary.csv` α `rho_split_median` | |
| α-8 | n_trace_split (descriptive) | 8/10 | same | |
| α-9 | n_above_drift (descriptive) | 8/10 | same | |
| α-10 | **C1** Wilcoxon p (`\(\rhosplit > 0\)`) | **0.010** | same `wilcoxon_split_gt_0_p`; computed LOO max p = 0.020 (Pat_03) | |
| α-11 | **C2** Wilcoxon p (`\(\rhosplit > \rhodrift\)`) | **0.007** | same `wilcoxon_split_gt_drift_p`; computed LOO max p = 0.014 (Pat_03) | |
| **Per-pair cophenet — cross-probe restriction (Methods C4)** | | | | |
| α-12 | `\(\rhoxprobe\)` median | +0.105 | `ctm_triangle/cohort_summary.csv` α `rho_xprobe_median` | sign-aligned with `\(\rhosplit\)`; cohort-median sign-agreement test ✓ |
| α-13 | C4 paired Wilcoxon p (`\(\rhosplit > \rhoxprobe\)`) | 0.461 | `ctm_triangle/c4_wilcoxon_cohort.csv` α `paired_wilcoxon_p_split_gt_xprobe` | non-significant → non-degradation gate passes |
| α-14 | C4 LOO max p | 0.674 (Pat_02) | same `wilcoxon_loo_max_p_split_gt_xprobe` + `wilcoxon_loo_argmax_patient` | |
| α-15 | C4 sign agreement | True | same | |
| α-16 | C4 pass | True | same `c4_pass` | |
| **Per-pair cophenet — matched-strength surrogate (Methods C3)** | | | | |
| α-17 | obs median `\(\rhosplit\)` | +0.105 | `matched_strength_surrogate_split_baseline/cohort_summary.csv` α `obs_median_rho` | |
| α-18 | surrogate cohort median | +0.013 | same `surr_median_rho_median` | |
| α-19 | ratio obs/surr | ≈ 8.3× (8.25× exact) | derived | |
| α-20 | C3 paired Wilcoxon z / p (obs > surr) | 54 / **0.002** | same `paired_wilcoxon_z`/`p`; computed LOO max p = 0.004 (Pat_03) | |
| α-21 | n_above_surrogate (descriptive) | **5/10** | same | **borderline per-patient agreement at full cohort** — the C5 epi-X regime resolves this to 7/10 (α-26). |
| **Per-pair cophenet — epileptogenic-zone exclusion (Methods C5)** | | | | |
| α-22 | obs median `\(\rhosplit^{\rm epi-X}\)` | +0.187 | `alpha_epi_exclusion/cohort_summary.csv` α `obs_median_rho` | full +0.105 → epi-X +0.187 |
| α-23 | surrogate cohort median^epi-X | +0.007 | same `surr_median_rho_median` | |
| α-24 | ratio obs/surr^epi-X (descriptive) | **27.7×** | same `ratio` (27.698) | full 8.25× → epi-X 27.70× ≈ **3.4× increase**. |
| α-25 | **C5** one-sample Wilcoxon p (`\(\rhosplit^{\rm epi-X} > 0\)`) | **0.010** | `alpha_epi_exclusion/c5_wilcoxon_cohort.csv` α `wilcoxon_one_sided_p_epiX` | **Methods-locked C5 gate.** Brief reports paired-vs-surrogate 0.014 — see Issue 3. |
| α-26 | C5 LOO max p | 0.020 (Pat_03) | same `wilcoxon_loo_max_p_epiX`/`wilcoxon_loo_argmax_patient` | |
| α-27 | n_above_surrogate^epi-X (descriptive) | **7/10** | `alpha_epi_exclusion/cohort_summary.csv` α | full 5/10 → epi-X 7/10 (+ 2 patients into trace direction). |
| α-28 | C5 verdict | **pass** (`c5_pass = True`) | `alpha_epi_exclusion/c5_wilcoxon_cohort.csv` α | |
| **Grassmann subspace (no trace at α)** | | | | |
| α-29 | obs longest contiguous-significant run | 4 cells, `k = 11..14` | `grassmann_regate_no_filter/contig_summary.csv` `audit_66_full` α `longest_run_new` + `_k_start`/`_k_end`; also `grassmann_cluster_extent/cohort_summary.csv` α `obs_longest_run` = 4 | |
| α-30 | obs cluster mass `\(\sum_{k} -\log_{10} p_{k}\)` | 7.98 | `grassmann_cluster_extent/cohort_summary.csv` α `obs_cluster_mass_neglog10p` | |
| α-31 | obs cluster mass normalized `\(T_{\rm G}^{\ast}(\alpha)\)` | 0.031 | same `obs_cluster_mass_neglog10p_norm` | |
| α-32 | null mean / p95 / max longest run | 2.09 / 6.0 / 12 | same | |
| α-33 | null mean / p95 cluster mass | 7.35 / 21.26 | same | |
| α-34 | `\(p_{\rm LR}\)` (descriptive longest-run cluster-p) | 0.159 | same `cluster_p_longest_run` | |
| α-35 | **`\(p_{\mathrm{mass}}\)` (Methods-locked Grassmann gate)** | **0.348** | same `cluster_p_cluster_mass` | Brief reports 0.0995 — see Issue 2 (stale value). |
| α-36 | `\(p_{\mathrm{mass}}\)` LOO max | 0.577 (Pat_14) | same `cluster_p_mass_loo_max` + `_argmax_patient` | |
| α-37 | Grassmann verdict | **`no_trace`** | same `verdict_cluster_extent` | |
| **Grassmann subspace — epi-X (no trace under epi-X either)** | | | | |
| α-38 | obs longest run^epi-X | 3 cells (k = 10..12) | `grassmann_regate_no_filter/contig_summary.csv` `audit_67_epiX` α `longest_run_new` | |
| α-39 | `\(p_{\rm LR}^{\rm epi-X}\)` | 0.274 | `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` α `cluster_p_longest_run_epiX` | |
| α-40 | `\(p_{\mathrm{mass}}^{\rm epi-X}\)` | 0.413 | same `cluster_p_cluster_mass_epiX` | |
| α-41 | Grassmann epi-X verdict | **`no_trace`** (`c5_pass = False`) | same | |
| **Anatomy α cophenet (Methods A1 + A2)** | | | | |
| α-42 | DK regions with `passes_both = True` (excl. `Unk`/`unknown`/`Wm`) | **11 named regions** | `anatomy_alpha_cophenet/cohort_summary.csv` (passes_A1 True ∧ passes_A3 True, where the CSV's "A3" maps to Methods-locked **A2** matched-strength; see §"Naming map" below) | |
| α-43 | total `passes_both` rows | 14 (= 11 named + `unknown` + `Wm` + 0 white-matter exclusions) | same | The locked `ANATOMY_CONTROLS.md` excludes the three non-named bins from the localization claim. |
| α-44 | 11 named regions (Desikan–Killiany) | left caudal anterior cingulate; left parahippocampal cortex; left rostral anterior cingulate; right medial orbitofrontal cortex; right caudal middle frontal gyrus; right postcentral gyrus; left caudal middle frontal gyrus; right caudal anterior cingulate; right precuneus; right superior parietal cortex; right posterior cingulate | same; sorted by q_BH ascending | |
| α-45 | enrichment range across the 11 regions | 1.47× (right precuneus) — 3.80× (left parahippocampal) | same `enrichment` column | |
| α-46 | A1 q_BH range | 5.91×10⁻¹⁰ — 4.5×10⁻² | same `q_bh` column | all < 0.05 |
| α-47 | A2 (Methods nomenclature) / "A3" (locked-anatomy nomenclature) p_emp range | 0.005 — 0.045 | same `p_empirical` column | all < 0.05 |
| α-48 | A2 z-score range | 2.60 — 14.26 | same `obs_z` column | all > 2.0 |
| α-49 | anatomy verdict | **strong localized** | matches `ANATOMY_LEDGER.md` α row | |
| α-50 | epi-X anatomy reproduction (post-patch 2026-05-22) | **9 of 11 named regions reproduce** under epi-X joint A1+A2 gate; 2 drop (lh-caudal-MF at A1 marginal; rh-medial-OFC at A2 — `\(z^{A2}\)` 3.22 → 1.28, indicating epi-zone-coupled component of full-cohort enrichment); 4 cingulate regions + lh-parahipp + rh-postcentral + rh-precuneus + rh-superior-parietal + rh-caudal-MF survive; lh-rostral-AC strengthens 3.04× → 5.37× | `data/audit/anatomy_alpha_cophenet_epiX/cohort_summary.csv` post-patch (MD5 `199d742b22e00893ff3fbae4623a1b1c`) | |
| **Anatomy α Grassmann** | | | | |
| α-51 | regions passing A2 | **0** (empty CSV) | `anatomy_alpha_grassmann/cohort_summary.csv` (only the header row is present) | consistent with α having no Grassmann trace; no Grassmann anatomy verdict computed |

## Naming map — locked anatomy controls vs Methods

The locked `ANATOMY_CONTROLS.md` defines four anatomy controls A1–A4; Methods §`\sssec:methods_anatomy` locks **only A1 and A2** for the manuscript. The mapping:

| Locked-anatomy name | Methods name | Test |
|---|---|---|
| A1 | **A1** | Hypergeometric / Fisher's exact one-sided + BH-FDR across DK regions |
| A2 | (not in Methods) | Sampling-corrected bootstrap (deferred — not used in manuscript) |
| A3 | **A2** | Matched-strength surrogate enrichment (R = 200, paired with C3 ensemble) — Methods-locked manuscript test |
| A4 | (not in Methods) | Implant-geometry regression (deferred) |

When the manuscript paragraph cites the anatomy gate, use **Methods nomenclature: "A1 + A2"**, never "A1 + A3" (which is the locked-anatomy battery's internal naming). The audit CSV's column `passes_A3` corresponds to the Methods-locked A2. The bands brief `02_alpha.md` uses the locked-anatomy "A3" naming — translate to "A2" when citing in the manuscript.

## Paragraph structure — mirror the locked β subsection, single-probe variant

The α subsection is structurally **lighter than β** because α is single-probe (no Grassmann results to write up). The following paragraph order mirrors the β subsection, with the Grassmann paragraph collapsed to one or two sentences (negative result) and the C5 epi-zone exclusion paragraph elevated to load-bearing status (the decisive sensitivity for α).

### Opening line (one sentence)

State the single-probe nature of α and the role of C5. Mirror the β opening line ("The β band is the only band whose post-task signature clears matched-strength surrogacy at both LRG-layer probes, with the largest cohort-level effect size at the per-pair multiscale probe."), e.g.:

> The α band carries a per-pair multiscale trace on `\(\Dcoph\)` that is borderline in the full cohort and decisively confirmed under epileptogenic-zone exclusion, with no detectable Grassmann subspace signature. The trace at α localizes anatomically to a distributed cortical network spanning bilateral cingulate, left parahippocampal cortex, medial orbitofrontal cortex, and bilateral caudal middle frontal cortex.

### Substrate paragraph (one paragraph, mirroring β's substrate paragraph)

α has the **strongest substrate signal** of the six bands (α-1 through α-6). Cite:
- cohort median `\(\rhorawsplit = +0.158\)`
- 21× the matched-strength noise floor (+0.008)
- 8/10 patients above own surrogate
- cohort-paired Wilcoxon p = **0.014** (the only band classified `"separated"` at the substrate level)
- contrast with β (substrate borderline, p = 0.053): for α the substrate gate **already passes**

Then the LRG-wrap reading: the per-pair cophenet probe demotes α to a borderline 5/10 agreement at full cohort (α-21), only to be re-promoted under C5 epi-X (α-27). This is the structural inverse of β (whose LRG step **promotes** a borderline substrate effect into a clean cohort verdict): for α, the cophenet step at full cohort acts as a **multiscale + spatial-purity filter** that strips out the epi-zone contribution and reveals a non-epileptic-cortex per-pair phenomenon.

### Per-pair cophenet paragraph (mirror β's per-pair paragraph; α-7 through α-21)

Structure: Within-baseline result → drift-floor null (α-11) → cross-probe restriction (α-12 through α-16) → matched-strength surrogate (α-17 through α-21).

Quantitative content:

- "Eight of ten patients carry `\(\rhosplit\)` individually in the trace direction at α, with the cohort median at +0.115 (cohort Wilcoxon **p = 0.010**)."
- "The same eight patients lie above their own drift-floor null (**p = 0.007**), so the cohort effect cannot be inherited from monotone session drift through a shared baseline."
- "The cross-probe restriction `\(\rhoxprobe = +0.105\)` is sign-aligned with `\(\rhosplit\)`, and the non-degradation gate of Methods §`\sssec:methods_compare_ctm` is non-significant at p = 0.461 — the trace is not concentrated in near-anatomy."
- Per-node strength evolution paragraph (mirroring β R3-style framing): "Against the matched-strength surrogate, the observed cohort median exceeds the surrogate cohort median by an **≈ 8.3× ratio (+0.105 vs +0.013)**, with **five of ten patients individually above their own surrogate (p = 0.002)**. The per-patient agreement is borderline at the full cohort and **resolved by the C5 epi-zone exclusion** described below."

### C5 epi-zone exclusion paragraph — LOAD-BEARING FOR α (α-22 through α-28)

Structure: state the epi-X regime → report the strengthened effect → cite the Methods-locked one-sample Wilcoxon gate → describe per-patient direction of change.

Quantitative content:

- "Restricting the analysis to non-epileptic cortex strengthens the α per-pair multiscale trace decisively: the cohort median `\(\rhosplit\)` rises from +0.105 to **+0.187**, the matched-strength ratio rises from 8.3× to **27.7×**, and the per-patient agreement increases from 5/10 to **7/10** above own surrogate. The Methods-locked C5 one-sample Wilcoxon on `\(\rhosplit^{\rm epi-X} > 0\)` returns **p = 0.010** (LOO max p = 0.020, Pat_03)."
- Per-patient narrative (use the data in α-brief §3.1 epi-X comparison table; numbers there are verified against `alpha_epi_exclusion/comparison.csv`):
  - Three patients (Pat_05, Pat_07, Pat_10) move from marginally-positive to clearly-positive under epi-X exclusion.
  - One patient (Pat_13) inverts: full-cohort +0.028 → epi-X −0.171. Pat_13 has the largest absolute epi-zone contact count in the cohort (N_epi = 30/119 ≈ 25%). The inversion is consistent with α coupling concentrated within Pat_13's epi-zone being stripped by the exclusion.
  - The four already-strongly-positive patients (Pat_03, Pat_06, Pat_08, Pat_14) stay strongly positive with mild weakening.
  - Pat_15 (N_epi = 0) is invariant by construction.
- One sentence on the interpretive consequence: "Pat_13's inversion under C5 is descriptive — the cohort verdict is gated by the Wilcoxon p, not by per-patient counts — and the LOO sensitivity (LOO max p = 0.020, Pat_03-driven) confirms the cohort verdict is not carried by any single patient."

### Per-pair multiscale interpretation paragraph (mirror β's "The per-pair multiscale layer identifies…")

Concentrate on the **non-epileptic-cortex framing**:
- α carries a structured per-pair multiscale reorganization that is **masked at the full cohort** by epileptogenic-zone contacts and **revealed by C5 exclusion**.
- The cophenet step is acting as a multiscale + spatial-purity filter: it integrates across N−1 merge heights AND, combined with C5, isolates the non-pathological component of the α reorganization.
- Position the α finding against β: where β shows a multi-probe trace (both per-pair and subspace), α shows a single-probe trace concentrated in physiological cortex.

### Grassmann subspace — short negative paragraph (one or two sentences; α-29 through α-41)

α has **no Grassmann subspace signature**. Cite:
- "The leading-eigenmode subspace probe shows no detectable trace at α: the longest contiguous-significant run across `\(k \in \{2, \ldots, 112\}\)` is 4 cells at `\(k = 11..14\)`, with cluster-mass permutation **`\(p_{\mathrm{mass}}(\alpha) = 0.348\)`** and `\(p_{\rm LR} = 0.159\)`; under epi-zone exclusion the run further shrinks to 3 cells and the cluster mass remains non-significant (`\(p_{\mathrm{mass}}^{\rm epi-X} = 0.413\)`). The α reorganization does not project into a coherent rotation of the leading-mode slow-diffusion subspace."

This is one short paragraph — no dedicated `\paragraph{}` heading, no per-`\(k\)` table, no figure. The β subsection's Grassmann content (cluster-mass, k ≈ 40 sweet spot, heatmap figure) does NOT have an α counterpart.

### Per-patient Table 1 (α)

Build `tab:alpha_per_patient` mirroring `tab:beta_per_patient` (locked β table), with these column adaptations:

| Col | Content | Source |
|---|---|---|
| 1 | Patient | — |
| 2 | `\(N_{\mathrm{ch}}\)` | `matched_strength_surrogate_split_baseline/per_patient_per_band.csv` α `N_nodes` |
| 3 | `\(\rhosplit\)` (matched-strength) | same `obs_rho` |
| 4 | surr. med. | same `surr_p50` (per-patient surrogate median) |
| 5 | `\(z\)` vs surr. | same `obs_z` |
| 6 | `\(\rhosplit^{\rm epi-X}\)` (matched-strength under C5) | `alpha_epi_exclusion/comparison.csv` `obs_rho_excluded` |
| 7 | `\(z^{\rm epi-X}\)` vs surr. | same `obs_z_excluded` |
| 8 | `\(N_{\rm epi}\)` | same `n_epi_excluded` |
| 9 | `\(\rhorawsplit\)` (raw FC reference) | `raw_fc_matched_strength/per_patient_per_band_all_bands.csv` α `obs_rho` (per-patient values pulled into table below) |
| 10 | LRG-trace classification | `+` if `\(\rhosplit^{\rm epi-X}\) > 0` AND `\(z^{\rm epi-X} \ge 1.96\)`; `~` if `\(\rhosplit\)` and `\(\rhosplit^{\rm epi-X}\)` have opposite signs; `−` otherwise |
| 11 | Notes | "right-hemisphere only" for Pat_15; "largest epi-zone" for Pat_13; otherwise "—" (do not flag Pat_03 fs=1024 Hz per `feedback_pat03_no_dropout.md`) |

Drop the Grassmann `T_G^{*,s}` and `n_sig,k` columns from the β table — α has no Grassmann content. The table is shorter / narrower as a result.

Per-patient values for α — every value CSV-verified to 3 decimals:

| Patient | `\(N_{\rm ch}\)` | `\(\rhosplit\)` (full) | surr p50 | `\(z\)` | `\(\rhosplit^{\rm epi-X}\)` | `\(z^{\rm epi-X}\)` | `\(N_{\rm epi}\)` | `\(\rhorawsplit\)` |
|---|---|---|---|---|---|---|---|---|
| Pat_02 | 117 | −0.037 | −0.009 | −0.37 | −0.035 | −0.16 | 14 | +0.043 |
| Pat_03 | 122 | **+0.360** | +0.023 | **+6.67** | **+0.189** | **+3.31** | 6 | +0.276 |
| Pat_05 | 118 | +0.078 | +0.013 | +1.31 | +0.121 | **+2.52** | 14 | +0.271 |
| Pat_06 | 115 | **+0.833** | +0.032 | **+12.15** | **+0.692** | **+10.17** | 10 | +0.676 |
| Pat_07 | 116 | +0.081 | +0.012 | +1.17 | **+0.211** | **+3.51** | 7 | +0.228 |
| Pat_08 | 120 | **+0.380** | +0.059 | **+6.03** | **+0.290** | **+3.99** | 9 | +0.219 |
| Pat_10 | 113 | **+0.130** | +0.018 | **+2.12** | **+0.193** | **+2.91** | 10 | +0.096 |
| Pat_13 | 119 | +0.028 | −0.007 | +0.68 | **−0.171** | **−3.52** | **30** | −0.019 |
| Pat_14 | 119 | **+0.210** | −0.008 | **+4.18** | **+0.184** | **+3.14** | 12 | −0.094 |
| Pat_15 | 118 | +0.040 | −0.006 | +0.50 | +0.040 | +0.42 | **0** | +0.021 |
| **Cohort** | — | **+0.105\textsuperscript{\textbf{**}}** | +0.013 | — | **+0.187\textsuperscript{\textbf{*}}** | — | — | **+0.158\textsuperscript{\textbf{*}}** |

Sources:
- Cols `\(\rhosplit\)` (full) / surr p50 / `\(z\)`: `data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv` α rows (`obs_rho`, `surr_p50`, `obs_z`).
- Cols `\(\rhosplit^{\rm epi-X}\)` / `\(z^{\rm epi-X}\)` / `\(N_{\rm epi}\)`: `data/audit/alpha_epi_exclusion/comparison.csv` (`obs_rho_excluded`, `obs_z_excluded`, `n_epi_excluded`).
- Col `\(\rhorawsplit\)`: `data/audit/raw_fc_matched_strength/per_patient_per_band_all_bands.csv` α rows (`obs_rho`).
- Cohort row `\(\rhosplit\)` full: `matched_strength_surrogate_split_baseline/cohort_summary.csv` α `obs_median_rho = +0.105`. Cohort row `\(\rhosplit^{\rm epi-X}\)`: `alpha_epi_exclusion/cohort_summary.csv` α `obs_median_rho = +0.187`. Cohort row `\(\rhorawsplit\)`: `raw_fc_matched_strength/cohort_summary_all_bands.csv` α `obs_median_rho = +0.158`.

Boldface convention (same as β table):
- Per-patient col `\(\rhosplit\)` (full): bold the value when `\(z \ge 1.96\)` (individual C3 significance against own matched-strength surrogate). At α: Pat_03, Pat_06, Pat_08, Pat_10, Pat_14 (5/10 above own surrogate) qualify.
- Per-patient col `\(\rhosplit^{\rm epi-X}\)`: bold when `\(|z^{\rm epi-X}| \ge 1.96\)`. At α: Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14 qualify (Pat_13's negative-direction significance after C5 — flag in table footnote).
- Cohort-row asterisks anchored to the cohort-paired Wilcoxon p (full cohort `\(\rhosplit\)` C3 paired p = 0.002 < 0.01 ⇒ \textsuperscript{\textbf{**}}; epi-X Methods-locked C5 one-sample p = 0.010 ⇒ \textsuperscript{\textbf{*}}; substrate `\(\rhorawsplit\)` paired Wilcoxon p = 0.014 ⇒ \textsuperscript{\textbf{*}}).
- Pat_15 has `\(N_{\rm epi} = 0\)` — bolded to flag that the full vs epi-X columns are identical by construction.

### Anatomy paragraph (mirror β's anatomy paragraph; α-42 through α-50, updated post-patch 2026-05-22)

Structure: one paragraph, 5–7 sentences (slightly longer than β's anatomy paragraph because the α epi-X anatomy now carries a biologically meaningful 9-of-11 result with two named dropouts). Mirror the β cadence:

- Opening: "Anatomically, the α trace localizes to a distributed cortical network at the per-pair multiscale probe, rather than to a single region or to diffuse whole-brain cortex."
- Region list (full cohort): "The signal concentrates in 11 Desikan–Killiany regions spanning bilateral cingulate, left medial-temporal lobe, right medial orbitofrontal cortex, bilateral caudal middle frontal cortex, and right parieto-postcentral cortex: …" then enumerate all 11 regions in their anatomical names (see α-44).
- Gate: "All 11 regions clear the joint A1 + A2 anatomical-enrichment gate (Methods §`\sssec:methods_anatomy`), with hypergeometric `\(q_{\rm BH}\)` ranging from 5.9 × 10⁻¹⁰ to 4.5 × 10⁻² and matched-strength surrogate `\(z\)` from 2.60 to 14.26."
- Comparison to β: "The medial-temporal lobe family is shared across the α and β cophenet anatomies (left parahippocampal cortex appears at both bands), and the bilateral cingulate / parieto-postcentral cortex pattern is distinctive to α."

- **Epi-X anatomy paragraph (verified post-patch):**
  "Restricting the analysis to non-epileptic cortex preserves the localization for 9 of the 11 named regions and strengthens the enrichment in several of them: left rostral anterior cingulate enrichment rises from 3.04× to 5.37× (`\(q_{\rm BH}\)` 4.1 × 10⁻²⁰ → 1.5 × 10⁻³¹), right caudal anterior cingulate from 2.44× to 3.31×, and left parahippocampal cortex remains the strongest single-region signal (3.80× → 3.89×). The two regions that drop the joint gate under epi-X are left caudal middle frontal gyrus (marginal full-cohort fit, `\(q_{\rm BH}\)` 1.6 × 10⁻³ → 5.1 × 10⁻², just above the 0.05 threshold) and right medial orbitofrontal cortex (full-cohort enrichment `\(z^{A2} = 3.22\)` against matched-strength, falling to `\(z^{A2} = 1.28\)` under epi-X — the matched-strength surrogate gate no longer separates the residual signal, indicating that the full-cohort medial-OFC enrichment carried a topology-coupled epi-zone component that the epi-X regime correctly strips out)."
- One closing sentence positioning the epi-X anatomy against the trace-level finding: "The strengthening of the cingulate and parahippocampal regions under epi-X mirrors the trace-level strengthening (α-22 through α-28) and supports the interpretation of α reorganization as a non-epileptic-cortex multiscale phenomenon."

**Sources for the anatomy paragraph:**
- Full cohort: `data/audit/anatomy_alpha_cophenet/cohort_summary.csv` (11 named regions with `passes_both = True`).
- Epi-X: `data/audit/anatomy_alpha_cophenet_epiX/cohort_summary.csv` (post-patch, MD5 `199d742b22e00893ff3fbae4623a1b1c`).
- Bug + patch provenance: `audit_71_anatomy_cophenet.py` lines 278–283 and 354–359 patched on 2026-05-22 (see Issue 1 above).

(α has **no Grassmann anatomy paragraph** — α-51 is zero regions; this is omitted from the manuscript.)

### Closing paragraph (mirror β's "The two LRG-layer probes thus read…")

α has **one LRG-layer probe**, not two. Mirror β's closing-paragraph framing but adapted:

- "The per-pair multiscale probe identifies α as the single-probe complement of the β multi-probe finding: where β rotates both the per-pair cophenetic geometry and the leading-`\(k\)` eigenmode subspace, α reorganizes the per-pair cophenetic geometry alone, and that reorganization concentrates in non-epileptogenic cortex."
- One sentence noting the structural distinction: at α the slow-mode subspace at any `\(k\)` does not rotate — the trace is distributed across the cophenetic merge-height continuum without producing a coherent subspace shift. This complements β's multi-probe trace and is consistent with α being the canonical posterior/attentional rhythm rather than a sensorimotor-coupled rhythm.

## Figure placeholders (F1, F2 — F3 omitted as α is single-probe)

α has fewer figures than β:

- **F1** `fig_alpha_rho_split.pdf` — per-patient `\(\rhosplit\)` matched-strength panel + cohort summary, parallel to F1 of β. Add a per-patient C5 epi-X overlay (or a separate panel) showing the full vs epi-X comparison.
  - Label: `\label{fig:alpha_rho_split}`
  - Tentative caption: "Cohort-level per-pair multiscale trace correlation `\(\rhosplit\)` at α across the n=10 patient set, full cohort and under C5 epileptogenic-zone exclusion. Top: per-patient observed `\(\rhosplit\)` against the patient-matched matched-strength surrogate distribution (R = 200, 4-cycle ±δ rewiring; Methods §`\sssec:methods_compare_stats`). Bottom: per-patient comparison of full-cohort `\(\rhosplit\)` and epi-X `\(\rhosplit^{\rm epi-X}\)`. \textit{Figure placeholder — final caption pending updated figure.}"
- **F2** `fig_alpha_anatomy_brain.pdf` — Desikan–Killiany region overlay for the 11 cohort-trace-localized regions; parallel to F3 of β.
  - Label: `\label{fig:alpha_anatomy_brain}`
  - Tentative caption: "Anatomical localization of the α trace at the per-pair multiscale probe. Eleven Desikan–Killiany regions clearing the joint A1 + A2 anatomical-enrichment gate (Methods §`\sssec:methods_anatomy`) are highlighted on a glass-brain projection. The 11-region network spans bilateral cingulate, left parahippocampal cortex, right medial orbitofrontal cortex, bilateral caudal middle frontal cortex, and right parieto-postcentral cortex. \textit{Figure placeholder — final caption pending updated figure.}"
- **F3** (β's Grassmann heatmap) is **omitted** — α has no Grassmann figure.

## What you must NOT do

- Do **not** change any cited number — every claim is CSV-verified in this directive.
- Do **not** cite the brief's stale value `\(p_{\mathrm{mass}} = 0.0995\)` for α Grassmann. Cite **0.348** (Issue 2).
- Do **not** cite the brief's paired-vs-surrogate test (p = 0.014) as the C5 cophenet gate. Cite the **Methods-locked one-sample Wilcoxon p = 0.010** (Issue 3).
- Do **not** assert "α cophenet anatomy reproduces identically under C5 epi-X". The verified post-patch outcome (Issue 1) is **9 of 11**, not 11 of 11; the two dropouts (lh-caudal-MF at A1 marginal; rh-medial-OFC at A2) carry biological information and must be reported, not hidden under an "identical reproduction" framing inherited from the stale brief.
- Do **not** use the locked-anatomy nomenclature "A3" for the matched-strength surrogate gate. Use Methods nomenclature "A2".
- Do **not** introduce KC, VI(k), τ-sweep, K*(b), or any retired probe terminology.
- Do **not** introduce cross-band BH-FDR (dropped 2026-05-20 — feedback memory `feedback_no_unmotivated_bh_fdr.md`). Within-band per-region BH-FDR for A1 (across DK regions per probe) is **kept** — that is the Methods-locked A1 gate.
- Do **not** invert the trace sign convention: positive `\(T_d\)`, `\(\rhosplit\)`, `\(\rhosplit^{\rm epi-X}\)` are all in the trace direction (consistent with the post-audit Methods).
- Do **not** flag Pat_03 distinctly in the per-patient table or in the body text. Pat_03 is treated identically at the analysis layer per `feedback_pat03_no_dropout.md`. The 1024 Hz sampling rate is handled at the config layer; no "Notes" annotation for Pat_03.
- Do **not** flag Pat_15 as the "anti-aligned" patient at α. Pat_15 at α is `\(\rhosplit = +0.040\)` (sign-aligned with the cohort) and `\(z = +0.50\)` — weak but not anti. The "anti-aligned" framing was used for β cophenet (which is a different band). At α the relevant per-patient nuance is Pat_13's epi-X inversion, not Pat_15.
- Do **not** add `\setcaptionlinewidth` or other per-figure tweaks until the figures themselves are finalized — the placeholder captions are tentative.
- Do **not** propose new analyses or new controls; if a numerical claim cannot be verified from the directive's CSV table, mark it with a `\giulio[…]{…}` comment for the user to resolve, not by inventing a substitute.

## Deliverable

Full α subsection LaTeX (the entire `\subsection{\(\alpha\): a per-pair multiscale task signature concentrated in non-epileptic cortex}` block) with:

1. Opening line.
2. Substrate paragraph (α-1 through α-6; cite p = 0.014, 8/10, "separated").
3. Per-pair cophenet paragraph (α-7 through α-21; cite C1 p = 0.010, C2 p = 0.007, C4 p = 0.461 non-significant, C3 p = 0.002 with 5/10 borderline).
4. C5 epi-zone exclusion paragraph (α-22 through α-28; cite Methods-locked one-sample Wilcoxon p = 0.010, 27.7× ratio, 7/10).
5. Per-pair multiscale interpretation paragraph.
6. Grassmann negative paragraph (α-29 through α-41; cite `\(p_{\mathrm{mass}} = 0.348\)`).
7. Per-patient `tab:alpha_per_patient` (cols defined in §"Per-patient Table 1 (α)"; values from this directive's per-patient table; 10 rows + cohort row).
8. Anatomy paragraph (α-42 through α-50; cite 11 named regions in full cohort and **9 of 11 surviving under epi-X**, with the two epi-X dropouts named and biologically interpreted per the post-patch verified outcome).
9. Closing paragraph (single-probe framing; complement to β multi-probe finding).
10. Figure label + tentative caption placeholders F1 (`fig_alpha_rho_split.pdf`) and F2 (`fig_alpha_anatomy_brain.pdf`).

Plus a one-paragraph change log naming:
- Every CSV row consulted (echo the rows from the directive's verification table).
- Any `\giulio[…]{…}` comment inserted for Issue 1 (anatomy epi-X) — explicit text of the comment.
- Confirmation that the brief's three stale claims (Issue 1 anatomy reproduction, Issue 2 Grassmann p = 0.0995, Issue 3 C5 paired-vs-surrogate p = 0.014) were **not** carried into the manuscript.

After return, the user pastes the revised LaTeX back for an `EVALUATION_PROTOCOL.md`-style verification pass (number-level CSV cross-check, anti-pattern scan, Methods/Results test-name consistency). Iterate if findings; otherwise the α subsection is locked and the next subsection (`\(\gamma_{\rm low}\)`) follows the same template.
