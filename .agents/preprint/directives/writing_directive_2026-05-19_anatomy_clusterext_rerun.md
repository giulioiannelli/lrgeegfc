---
name: writing-directive-anatomy-clusterext-rerun
era: IMCOH_ABS_COHORT_N10
status: directive
kind: writing-agent-directive
date: 2026-05-19 (pm)
target: writing-agent — cascade the cluster-extent anatomy rerun into LaTeX, briefs, and the methods text
priority: high (changes per-band Grassmann anatomy region sets; the β results LaTeX paragraph and bands/01_beta.md anatomy section must be updated before draft of α / γ_l / δ paragraphs)
source_of_truth:
  - .agents/preprint/locked/ANATOMY_LEDGER.md      # locked verdicts + 2026-05-19 (pm) revision note
  - .agents/preprint/locked/ANATOMY_CONTROLS.md    # locked A1-A4 battery (updated to S(b) aggregation)
  - .agents/preprint/directives/writing_directive_2026-05-19_methods_anatomy.md  # methods text spec
  - data/audit/anatomy_<band>_grassmann_clusterext/cohort_summary.csv  # NEW source CSVs
  - data/audit/anatomy_delta_grassmann_epiX_clusterext/cohort_summary.csv
companion_to:
  - writing_directive_2026-05-19_methods_anatomy.md       # methods §sssec:methods_anatomy
  - writing_directive_2026-05-19_beta_post_methods_revision.md  # β results LaTeX
---

# Writing directive — cascade the cluster-extent anatomy rerun (2026-05-19 pm)

**Manuscript ↔ lab label mapping (locked 2026-05-19 pm).** In the manuscript
methods, the anatomy controls are **A1 = hypergeometric** and **A2 =
matched-strength surrogate**. In the lab artifacts (`ANATOMY_CONTROLS.md`,
`ANATOMY_LEDGER.md`, audit CSV column names) the labels **A1 / A3 stay**, with
**lab A2 (sampling-corrected bootstrap)** and **lab A4 (implant-geometry
regression)** reserved for the sensitivity supplement and **not disclosed in
the manuscript methods**. Mapping: manuscript A1 ↔ lab A1; manuscript A2 ↔
lab A3. The diff tables and prose below use lab labels (z values come from
the matched-strength surrogate test = lab A3 = manuscript A2). Writing agent:
when generating LaTeX, use manuscript labels (A1 / A2).

**Head.** The Grassmann anatomy methodology has been corrected from the
retired longest-contiguous-significant `k`-window (`K*(b)`) to the locked
all-clusters significance-thresholded set `S(b) = {k : p_k(b) < α_k}` — the
**support of the cluster-mass statistic `T_G^*`**. Audit_72 has been re-run
for the four affected (band, probe) cells and the cohort_summary CSVs at
`data/audit/anatomy_<band>_grassmann{_epiX}_clusterext/` are now the **single
source of truth** for any Grassmann anatomy citation in the manuscript. **β
Grassmann anatomy is unchanged (same 7 regions)**; **γ_l, δ-full, and
δ-epi-X anatomy shifts substantially** — the writing agent must rewrite the
corresponding paragraphs from the new CSVs and retract the
locked-but-now-superseded region lists.

---

## What changed methodologically (and why)

### What was wrong

The previous audit_72 ran the Grassmann anatomy with per-node participation
averaged across a **contiguous-significant `k`-window** chosen from
`data/audit/grassmann_regate_no_filter/contig_summary.csv` (longest contiguous
run of `k` cells with cohort-paired Wilcoxon `p < 0.05`):

- β: `K*(β) = [27, 55]`, 29 cells
- γ_l: `K*(γ_l) = [12, 23]`, 13 cells
- δ-full: `K*(δ) = [57, 63]`, 7 cells
- δ-epi-X: `K*^epiX(δ) = [33, 39]`, 7 cells

Under the locked **all-clusters cluster-mass `T_G^*`** paradigm
(`feedback_cluster_mass_all_clusters.md`, `methods_grassmann_TG_star.md`),
there is **no canonical "significant window"** — the verdict gate is
`cluster_p_mass < 0.05` on the empirical-null permutation of
`T_G^*(b) = Σ_{k : p_k(b) < α_k} (−log10 p_k(b))`. The support of this sum is
the **set of all `k` cells clearing the per-`k` threshold**, irrespective of
contiguity:

```
S(b) = { k ∈ {2, ..., 112} : p_k(b) < α_k = 0.05 }
```

`S(b)` is generally a **superset** of `K*(b)` (it includes the longest
contiguous run plus any scattered cells outside it). Aggregating
participation over `K*(b)` was effectively reading a *subset* of the
significant subspaces — fine as a sensitivity, wrong as the canonical anatomy
measure.

### What is locked now

- **Trace-flagged unit at the Grassmann probe** = top decile (cohort-pooled)
  of per-node participation `part(i; b)`.
- **Per-node participation** aggregated **unweighted** across `S(b)`:
  ```
  part(i; b) = (1 / |S(b)|) · Σ_{k ∈ S(b)} ||U_k^T e_i||²_2
  ```
- **Significance-weighted aggregation is NOT used** — `T_G^*` already weights
  by `−log10 p_k`, so re-weighting in participation would double-count
  the same evidence; unweighted reads the cleaner question "across
  significance-thresholded subspaces, where do the heavy contacts sit".
- **`S(b)` source** for full-cohort runs:
  `data/audit/grassmann_cluster_extent/per_k_obs_p.csv`, filter
  `band=<b>, obs_p_one_sided_less < 0.05`.
- **`S^epiX(b)` source** for epi-X runs: derived from per-k cohort-paired
  Wilcoxon on `data/audit/grassmann_epi_exclusion/per_patient_per_band_per_k.csv`
  `obs_T_G` column, `alternative='less'`, `p < 0.05`.

### S(b) sizes (locked)

| Band | `|K*(b)|` (retired) | `|S(b)|` (locked) | k-cells in `S(b)` outside `K*(b)` |
|---|---|---|---|
| β | 29 (contiguous [27, 55]) | 40 (non-contiguous, span [21, 90]) | k = 21, 22, 57, 58, 59, 60, 61, 63, 73, 89, 90 |
| γ_l | 13 (contiguous [12, 23]) | 41 (non-contiguous) | many; see CSV |
| δ-full | 7 (contiguous [57, 63]) | 23 (non-contiguous) | many; see CSV |
| δ-epi-X | 7 (contiguous [33, 39]) | 22 (non-contiguous, [2, 34–52, 87, 88]) | k = 2, 34, 40–52, 87, 88 |

---

## Per-band diff — locked ledger vs S(b) rerun

### β Grassmann — verdict UNCHANGED

| Region | Locked z, p_emp | S(β) z, p_emp | Change |
|---|---|---|---|
| Hip | >>2, 0.005 | >>2, 0.005 | same |
| ctx-lh-insula | >>2, 0.005 | >>2, 0.005 | same |
| ctx-lh-middletemporal | 5.67, 0.005 | 5.27, 0.005 | minor |
| ctx-lh-lateralorbitofrontal | >>2, 0.005 | >>2, 0.005 | same |
| ctx-rh-medialorbitofrontal | >>2, 0.005 | >>2, 0.005 | same |
| ctx-lh-superiortemporal | 28.28, 0.005 | 28.28, 0.005 | same |
| ctx-rh-rostralmiddlefrontal | 3.58, 0.030 | 5.21, 0.010 | **strengthens** |

**Action for writing agent**: no change needed to the β Grassmann region list
in bands/01_beta.md or the β results LaTeX. **Update z and p_emp values for
ctx-lh-middletemporal (z 5.67 → 5.27) and ctx-rh-rostralmiddlefrontal
(z 3.58 → 5.21, p_emp 0.030 → 0.010)** from
`data/audit/anatomy_beta_grassmann_clusterext/cohort_summary.csv`.

### γ_l Grassmann — verdict SHIFTED (3 retained / 3 dropped / 4 added)

| Region | Locked | S(γ_l) | Note |
|---|---|---|---|
| ctx-lh-inferiortemporal | ✓ z=3.05 p=0.005 | — | **dropped** |
| ctx-lh-fusiform | ✓ z>>2 p=0.005 | — | **dropped** — KC-era "left fusiform at γ_l" claim retracts |
| ctx-lh-middletemporal | ✓ z=5.00 p=0.005 | ✓ z=3.78 p=0.005 | retained |
| ctx-lh-superiortemporal | ✓ z=11.06 p=0.005 | ✓ z=7.16 p=0.005 | retained |
| ctx-rh-paracentral | ✓ z>>2 p=0.005 | — | **dropped** |
| ctx-rh-parstriangularis | ✓ z=10.13 p=0.005 | ✓ z=9.01 p=0.005 | retained |
| ctx-lh-lateraloccipital | — | ✓ z>>2 p=0.005 | **new** |
| ctx-lh-rostralmiddlefrontal | — | ✓ z=16.30 p=0.005 | **new** |
| ctx-rh-medialorbitofrontal | — | ✓ z=3.36 p=0.005 | **new** |
| ctx-lh-cuneus | — | ✓ z=5.25 p=0.040 | **new** |

**Narrative change**: the "temporal-cortex-dominant" reading of γ_l is wrong
under the locked methodology. The new network is **occipito-temporal +
medial-OFC + frontal**: superior + middle temporal cortex retained, fusiform
+ inferior temporal + paracentral retracted, lateral occipital + cuneus +
rostral middle frontal + medial OFC added.

**Action for writing agent**: when drafting the γ_l anatomy paragraph,
rewrite the network description from scratch using
`data/audit/anatomy_low_gamma_grassmann_clusterext/cohort_summary.csv`.
Retract the "left fusiform" claim.

### δ-full Grassmann — verdict SHIFTED (1 retained / 5 dropped / 3 added)

| Region | Locked | S(δ) | Note |
|---|---|---|---|
| ctx-lh-fusiform | ✓ z=12.62 p=0.005 | — | **dropped** — KC-era fusiform retracts further |
| ctx-lh-inferiorparietal | ✓ z=5.07 p=0.005 | ✓ z>>2 p=0.005 | retained |
| ctx-rh-medialorbitofrontal | ✓ z=4.08 p=0.005 | — | **dropped** |
| ctx-rh-caudalanteriorcingulate | ✓ z=2.84 p=0.005 | — | **dropped** |
| Amy | ✓ z=2.33 p=0.030 | — | **dropped** — δ "anchor anatomy" (Amy + cingulate) retracts |
| ctx-lh-bankssts | ✓ z=2.42 p=0.045 | — | **dropped** |
| ctx-lh-inferiortemporal | — | ✓ z>>2 p=0.005 | **new** |
| ctx-rh-parstriangularis | — | ✓ z>>2 p=0.005 | **new** |
| ctx-lh-superiortemporal | — | ✓ z=14.11 p=0.010 | **new** |

**Narrative change**: the locked "δ anchor anatomy" framing (Amy + cingulate
+ fusiform overlap with low-frequency epi-near synchronization) is **not
supported** under S(δ). The new full-cohort δ Grassmann network is
**parietal + temporal + frontal** — distinct from the C4 cross-probe
known-biology layer of `memory/epileptic_imcoh_universal.md` (which is a
SEPARATE descriptive observation about δ cross-probe 1.55× ratio at the
substrate, not about Grassmann anatomy). The two layers should now be
discussed independently: (i) δ Grassmann anatomy = parietal-temporal-frontal,
(ii) δ substrate cross-probe = known epi-near low-frequency synchronization.

**Action for writing agent**: rewrite the δ-full Grassmann anatomy paragraph
using
`data/audit/anatomy_delta_grassmann_clusterext/cohort_summary.csv`. Move any
"anchor anatomy" framing to the δ substrate / C4 cross-probe discussion,
which the Grassmann probe is not testing.

### δ-epi-X Grassmann — verdict SHIFTED (2 retained / 4 dropped / 1 added)

| Region | Locked | S^epiX(δ) | Note |
|---|---|---|---|
| ctx-lh-fusiform | ✓ z>>2 p=0.005 | — | **dropped** |
| ctx-lh-inferiorparietal | ✓ z>>2 p=0.005 | — | **dropped** |
| ctx-lh-superiorparietal | ✓ z>>2 p=0.005 | ✓ z>>2 p=0.005 | retained |
| ctx-rh-postcentral | ✓ z=3.03 p=0.005 | — | **dropped** |
| ctx-rh-rostralmiddlefrontal | ✓ z=3.52 p=0.005 | ✓ z=2.71 p=0.005 | retained |
| ctx-lh-inferiortemporal | ✓ z=5.69 p=0.035 | — | **dropped** |
| ctx-lh-superiorfrontal | — | ✓ z=14.11 p=0.010 | **new** |

**Narrative change**: the parietal-dominant reading is preserved (superior
parietal kept) but the network is more compact. The
"parietal-postcentral-fusiform" combination of the locked ledger is replaced
by a tighter "superior-parietal + rostral-middle-frontal +
superior-frontal" set.

**Action for writing agent**: rewrite the δ-epi-X Grassmann anatomy
paragraph using
`data/audit/anatomy_delta_grassmann_epiX_clusterext/cohort_summary.csv`.

### δ full-vs-epiX dissociation — STRENGTHENS

Locked: 2/6 shared regions (fusiform, inferiorparietal). New under S(b):
**0/3 shared regions**. The "mixture of two distinct phenomena" claim
strengthens because the dissociation is now complete.

**Action for writing agent**: keep the dissociation claim; quantify
correctly as **0 shared regions** between the full-cohort δ Grassmann
network (parietal + temporal + frontal) and the C5 epi-X δ Grassmann
network (superior parietal + rostral middle frontal + superior frontal).

---

## What unchanged

- **β Grassmann region set**: identical 7 regions; only z and p_emp updates needed.
- **β cophenet anatomy**: untouched (no `k`-aggregation; top-decile per-pair `|Δρ_split^coph|`). All 7 regions remain locked as in `ANATOMY_LEDGER.md` β cophenet subsection.
- **α cophenet anatomy**: untouched. All 11 regions remain locked.
- **Methods §sssec:methods_anatomy**: the equation Eq.~\eqref{eq:methods_anatomy_node_unit} and the surrounding paragraph have been updated in `writing_directive_2026-05-19_methods_anatomy.md` (see that directive for the drop-in LaTeX replacement).
- **Verdict tags (`strong localized` / `weak localized` / `not localized`)** at every (band, probe) cell — all four trace-positive cells remain **strong localized** under the new region sets.
- **Cluster-mass `T_G^*` verdict gate**: untouched. `cluster_p_mass < 0.05` remains the trace-side gate; only the *anatomy-side aggregation set* changed.

---

## Files to update (writing-agent action list)

In strict order:

### 1. `bands/01_beta.md` anatomy section (cascading from locked ledger)

- Update β Grassmann region table to use the cluster-extent CSV values
  (z + p_emp changes for ctx-lh-middletemporal and ctx-rh-rostralmiddlefrontal).
- Source CSV row references must point to
  `data/audit/anatomy_beta_grassmann_clusterext/cohort_summary.csv` (not the
  retired `anatomy_beta_grassmann/cohort_summary.csv`).
- Region list itself is unchanged.

### 2. β results LaTeX anatomy paragraph

If the writing agent has already drafted the β anatomy paragraph (per
`writing_directive_2026-05-19_beta_post_methods_revision.md`), update only
the cited z and p_emp values for ctx-lh-middletemporal and
ctx-rh-rostralmiddlefrontal. The region list is unchanged. Add a parenthetical
("Methods §sssec:methods_anatomy aggregation over `S(β)`") if a single
citation of the methodology is desired.

### 3. Future α anatomy paragraph

Untouched — α cophenet anatomy is the load-bearing localization at α (D_coph
verdict); no Grassmann anatomy was audited for α (α has no Grassmann trace).
Use `data/audit/anatomy_alpha_cophenet/cohort_summary.csv` and the C5 epi-X
companion `anatomy_alpha_cophenet_epiX/cohort_summary.csv` as before.

### 4. Future γ_l anatomy paragraph

**Rewrite from scratch** using
`data/audit/anatomy_low_gamma_grassmann_clusterext/cohort_summary.csv`. The
new network (occipito-temporal + medial-OFC + frontal, 7 regions) replaces
the locked-but-superseded temporal-cortex-dominant 6-region list. Do **not**
cite "left fusiform at γ_l" — that claim is retracted under S(γ_l).

### 5. Future δ anatomy paragraphs

**Rewrite from scratch** in two subsections, one per ensemble:

- Full cohort: 4 regions from
  `data/audit/anatomy_delta_grassmann_clusterext/cohort_summary.csv`. Move
  any "anchor anatomy" framing to the δ substrate / C4 cross-probe layer.
- C5 epi-X: 3 regions from
  `data/audit/anatomy_delta_grassmann_epiX_clusterext/cohort_summary.csv`.
- Quantify dissociation as **0 shared regions** (was 2/6 under the retired
  windows).

### 6. Cross-band anatomy comparison table (if present in 00_cohort.md)

The cross-band Region-overlap table in `ANATOMY_LEDGER.md` (under "Cross-band
anatomy comparison") needs cascading updates:

- "Hippocampus / parahippocampal" row — unchanged (β / α only via cophenet).
- "Fusiform" row — **drop left from γ_l Grassmann**; **drop left from δ Grassmann (both windows)**. Fusiform now appears **nowhere** in the cluster-extent anatomy.
- "Insula (any hemisphere)" — unchanged.
- "Cingulate (any subregion)" — drop caudal ACC (rh) from δ Grassmann full.
- "Medial / lateral OFC" — add **medial (rh)** to γ_l Grassmann; remove medial (rh) from δ Grassmann full.
- "Inferior parietal" — keep left at δ Grassmann full; remove left from δ Grassmann epi-X.
- "Superior parietal / postcentral" — keep superior (lh) at δ Grassmann epi-X; remove postcentral (rh).
- "Temporal cortex (mid / sup / inf)" — γ_l keeps middle + superior, drops inferior; δ full **adds** superior + inferior (was none).
- "Amygdala" — **remove** δ Grassmann full entry. Amy now appears **nowhere** in the cluster-extent anatomy.
- "Bankssts" — **remove** δ Grassmann full entry.
- "Pars (orbitalis / triangularis)" — keep triangularis (rh) at γ_l; **add** triangularis (rh) at δ full.
- "Rostral / caudal middle frontal" — γ_l **adds** rostral (lh); δ full **drops** all entries; δ epi-X keeps rostral (rh).
- "Superior frontal" — **add** (lh) at δ epi-X. Unchanged at β.
- "Paracentral" — **drop** (rh) from γ_l (was the only entry).

---

## Anti-pattern reminders for the writing agent

- **Do not cite the retired `K*(b)` windows** anywhere — never write
  "across the contiguous-significant window `k = 27..55`" or analogous
  phrases. Use `S(b)` and reference `sssec:methods_anatomy`.
- **Do not cite the retired `anatomy_<band>_grassmann/cohort_summary.csv`
  CSVs**. The new source-of-truth CSVs all have the `_clusterext` suffix.
- **Do not carry the KC-era "Hippocampus + left fusiform at β" claim**. The
  cluster-extent rerun confirms: Hippocampus is at β Grassmann (Hip);
  left fusiform is **NOT** at β under either probe, NOT at γ_l under
  S(γ_l), NOT at δ under S(δ) or S^epiX(δ).
- **Do not lump "δ anchor anatomy" into the δ Grassmann claim**. Under
  S(δ), δ Grassmann anatomy is parietal-temporal-frontal. The δ C4
  cross-probe 1.55× known biology is a substrate-layer observation and
  goes in the substrate / C4 discussion only.
- **Do not double-count weighting**. Aggregation across `S(b)` is
  unweighted; `T_G^*` already weights by `−log10 p_k`. Do not add a
  second weighting in the participation aggregation.

---

## Cross-reference to source-of-truth files

| Reference | Path |
|---|---|
| Methods text spec | `.agents/preprint/directives/writing_directive_2026-05-19_methods_anatomy.md` |
| Anatomy battery lock | `.agents/preprint/locked/ANATOMY_CONTROLS.md` |
| Anatomy verdict ledger | `.agents/preprint/locked/ANATOMY_LEDGER.md` (see 2026-05-19 pm revision entry) |
| β cluster-extent CSV | `data/audit/anatomy_beta_grassmann_clusterext/cohort_summary.csv` |
| γ_l cluster-extent CSV | `data/audit/anatomy_low_gamma_grassmann_clusterext/cohort_summary.csv` |
| δ-full cluster-extent CSV | `data/audit/anatomy_delta_grassmann_clusterext/cohort_summary.csv` |
| δ-epi-X cluster-extent CSV | `data/audit/anatomy_delta_grassmann_epiX_clusterext/cohort_summary.csv` |
| S(b) source | `data/audit/grassmann_cluster_extent/per_k_obs_p.csv` |
| S^epiX(δ) source | `data/audit/grassmann_epi_exclusion/per_patient_per_band_per_k.csv` (derived per-k Wilcoxon) |

---

## Revision history

- **2026-05-19 (pm)** — Initial directive. Cluster-extent anatomy rerun
  complete for β, γ_l, δ-full, δ-epi-X. β verdict unchanged; γ_l, δ-full,
  δ-epi-X verdicts shifted. β cophenet, α cophenet untouched.
