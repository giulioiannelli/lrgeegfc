---
date: 2026-05-08
era: COHORT_N10 / IMCOH_ABS
status: current
type: verification-report
scope: §5.5 anatomy — p95 calibration definition, α ctx-lh-parsopercularis sensitivity, Pat_07/Pat_15 per-band cortical counts
inputs:
  - scripts/01_compute/audit/audit_39_per_leaf_sigma_aggregate.py (per-leaf rho)
  - scripts/01_compute/audit/audit_39b_per_leaf_rho_demeaned.py (demeaning)
  - scripts/01_compute/audit/audit_42_per_leaf_rho_null.py (within-baseline null + p95)
  - scripts/01_compute/audit/audit_49_anatomy_sensitivity_55.py (Bonferroni m=48 + dropouts)
  - data/audit/lrg_localization_anatomy/per_trace_leaf.csv
  - data/audit/per_leaf_rho_null/calibrated_trace_leaves.csv
  - data/audit/per_leaf_rho_null/cohort_calibration_summary.csv
  - data/raw/stereoeeg_patients/Pat_NN/implant_pat_NN.csv
---

# §5.5 anatomy — verification

**Head.** **The §5.5 prose has a wrong patient attribution at α ctx-lh-parsopercularis.** The cell does **not** rest on Pat_02 + Pat_13; it rests on **Pat_07 (3 trace-leaves) and Pat_14 (2 trace-leaves)**. Pat_02 contributes zero α parsopercularis trace-leaves. This is consequential because **Pat_07 is one of the cross-probe anti-aligned patients of §5.3/§5.4**, so the α parsopercularis cell is primarily driven by an anti-aligned patient — the inverse of the γ_l fusiform pattern (where the load-bearing patient Pat_02 is one of the strongest cohort-pro-trace patients).

The sensitivity battery is finished here:
- **Pat_07 dropout collapses the cell** (5/26 → 2/20, p_hyper 0.0014 → 0.12, **fails filter** since only Pat_14 remains).
- Pat_14 dropout also collapses the filter (3/16 from Pat_07 alone, p=0.019, fails n_pat_trace ≥ 2).
- Pro-cohort (drop Pat_07 + Pat_15) collapses to 2 leaves, single-patient (Pat_14), p=0.13.
- **Pat_02 dropout paradoxically strengthens the cell** to Bonferroni significance (5/20, 7.35×, p=3.4e-4, < 0.00104) because Pat_02's α trace-leaves lie in *other* regions and removing them sharpens the parsopercularis enrichment.
- Pat_03 dropout is essentially a no-op (5/26, 5.11×, p=0.0020).

**Verdict:** the α ctx-lh-parsopercularis cell is **two-patient, Pat_07-dominant, Pat_14-secondary**, and is fragile under any single-contributor dropout. It should be downgraded in §5.5 prose to "directional indication, two-patient base (Pat_07, Pat_14), Pat_07-dominant, fragile under dropout — Pat_07 is also the cross-probe anti-aligned patient at β/low-γ" with the dropout numbers shown.

The other two asks (A and C) reproduce cleanly. Per-band cortical baseline rates match the §5.5 prose to two decimals (α 3.46%, β 5.31%, γ_l 8.66%). Pat_07 / Pat_15 per-band cortical totals confirm the prose's "9 total β trace-leaves, mostly Wm/Unk" framing for both.

---

## A. p95 calibration — methods paragraph

**A.1 Per-leaf score s_ℓ.** For each (patient p, band b, leaf ℓ), the score is the **demeaned per-leaf Spearman correlation** `rho_demeaned(p, b, ℓ)`:

```
Δ_task[ℓ, j] = D^test[ℓ, j] − D^pre_A[ℓ, j]
Δ_rest[ℓ, j] = D^post[ℓ, j] − D^pre_B[ℓ, j]
rho_ℓ        = Spearman( Δ_task[ℓ, ·], Δ_rest[ℓ, ·] )    (over j ≠ ℓ)
rho_demeaned(ℓ) = rho_ℓ − mean_ℓ' rho_ℓ'                  (mean over leaves of patient/band)
```

i.e. for each leaf ℓ, the per-leaf Spearman is computed over the N−1 communication distances containing ℓ in the ultrametric (`D^test`, `D^pre_A`, `D^pre_B`, `D^post` are the post-LRG ultrametric matrices on `task_test`, `rest_pre_A`, `rest_pre_B`, `rest_post`). The score for leaf ℓ is then *demeaned* by subtracting the patient/band leaf-mean ρ to remove patient-global rotation effects (this is the audit_39b correction that fixes Pat_06's apparent localization, which was a global rank-rotation rather than a localized signal). The score family is in `data/audit/per_leaf_rho_demeaned/leaf_rho_demeaned.csv`.

**A.2 Null distribution.** The threshold is calibrated against a **within-baseline split-half null** (option *i* in the verification request — a control-anchored permutation null, not the empirical patient distribution of *s_ℓ*). The null replaces the (task, post) contrast with two within-rest contrasts on the same patient:

```
Δ_drift_a[ℓ, j] = D^pre_B[ℓ, j]  − D^pre_A[ℓ, j]    (within-RPre split-half noise)
Δ_drift_b[ℓ, j] = D^post_B[ℓ, j] − D^post_A[ℓ, j]   (within-RPost split-half noise)
rho_null_ℓ      = Spearman( Δ_drift_a[ℓ, ·], Δ_drift_b[ℓ, ·] )
rho_null_demeaned(ℓ) = rho_null_ℓ − mean_ℓ' rho_null_ℓ'
```

This is the **per-leaf analogue** of the cohort-level `rho_split_drift` floor used in the §5.3 CTM controlled test. The null measures how much per-leaf rank correlation arises from within-rest split-half noise alone, with no task. The null distribution is over leaves of the patient (same demeaning, same sample size), built from the same halves cache as the §5.3 CTM null.

**A.3 Calibration scope: per-(patient, band).** A separate threshold τ_95(p, b) is computed for every (patient, band) cell — **band-specific within patient**, and patient-specific within band. There is no pooling across bands. The 95th percentile is computed over the patient/band leaves of `rho_null_demeaned`.

**A.4 One-sided upper tail.** A leaf ℓ is flagged as a calibrated trace-leaf if `rho_demeaned(ℓ) ≥ τ_95(p, b)`. This is the upper tail only — no two-tailed flag, no negative-tail flag.

**A.5 Cohort-aggregate trace-leaf yields per band (cortical, after Wm/Unk drop).** Numerator = K (band-specific cortical trace-leaf count, pooled across n=10 patients); denominator = N_total = 866 cortical contacts pooled across cohort:

| band | K (cortical trace-leaves) | N_total | base_rate | §5.5 prose |
|---|---:|---:|---:|---:|
| δ | 56 | 866 | 6.47% | — |
| θ | 40 | 866 | 4.62% | — |
| α | **30** | **866** | **3.46%** | 3.5% ✓ |
| β | **46** | **866** | **5.31%** | 5.3% ✓ |
| low-γ | **75** | **866** | **8.66%** | 8.7% ✓ |
| high-γ | 64 | 866 | 7.39% | — |

The three §5.5-prose baseline rates reproduce to the second decimal. The numerators and denominators above are the missing K / N_total to insert into the methods paragraph.

---

## B. α ctx-lh-parsopercularis sensitivity battery

### B.1 — Patient breakdown of the 5 trace-leaves

The trace-leaves at α ctx-lh-parsopercularis (full cohort, p95-calibrated):

| patient | leaf_id | channel | rho_demeaned | τ_95(patient, α) |
|---|---:|---|---:|---:|
| **Pat_07** | 81 | Y2 | 0.2786 | 0.2657 |
| **Pat_07** | 82 | Y3 | 0.4024 | 0.2657 |
| **Pat_07** | 84 | Y5 | 0.4000 | 0.2657 |
| **Pat_14** | 58 | R'3 | 0.2925 | 0.2722 |
| **Pat_14** | 60 | R'5 | 0.3608 | 0.2722 |

**Patient totals:** Pat_07 = 3, Pat_14 = 2, all other patients = 0.

**The §5.5 prose's "two-patient base (Pat_02 and Pat_13)" is wrong.** The two patients are **Pat_07 and Pat_14**. Pat_07 carries **3/5 = 60%** of the trace-leaves and is the load-bearing patient — directly parallel to Pat_02 / γ_l fusiform but with a more-fragile 3:2 split.

### B.2–B.4 + B.5–B.6 — Sensitivity battery

`enrich = trace_rate / cohort_base_rate`. p_hyper = `hypergeom.sf(n_trace − 1, N_total, K, n_contacts)`. Bonferroni α/m = 0.00104 (m = 48). Filter = (n_contacts ≥ 5) ∧ (n_pat_trace ≥ 2) ∧ (enrichment > 1).

| scenario | n_trace | n_contacts | n_pat | base_rate | enrich | p_hyper | Bonf? | filter? |
|---|---:|---:|---:|---:|---:|---:|:-:|:-:|
| **full n=10** | 5 | **26** | 2 | 3.46% | **5.55×** | **0.0014** | no | yes |
| drop Pat_02 (n=9) | 5 | 20 | 2 | 3.40% | **7.35×** | **3.4e-4** | **YES** | yes |
| drop Pat_03 (n=9) | 5 | 26 | 2 | 3.77% | 5.11× | 0.0020 | no | yes |
| **drop Pat_07** (n=9) | 2 | 20 | **1** | 3.04% | 3.29× | **0.121** | no | **NO** |
| drop Pat_14 (n=9) | 3 | 16 | **1** | 3.72% | 5.04× | **0.019** | no | **NO** |
| pro-cohort (n=8, drop Pat_07+Pat_15) | 2 | 20 | **1** | 3.13% | 3.20× | **0.127** | no | **NO** |

**Reading:**

1. **Full cohort (verification of §5.5 numbers).** 5/26, 5.55×, p_hyper = 0.0014. The §5.5 prose figures (5/26, 5.55×, 1.4e-3) reproduce exactly. The cell narrowly misses Bonferroni m=48 (threshold 0.00104). ✓

2. **Pat_02 dropout sharpens to Bonferroni-pass.** The cell promotes to **5/20, 7.35×, p=3.4e-4** under Pat_02 dropout — well below the 0.00104 threshold. This is because Pat_02 carries 4 α trace-leaves in *other* cortical regions (none in parsopercularis) plus 6 α parsopercularis *contacts that are not trace-leaves*, so Pat_02 dilutes both the cohort K and the parsopercularis denominator. Removing Pat_02 sharpens the cell rather than collapsing it. (This is the **opposite direction** from the γ_l fusiform Pat_02 dropout, which collapsed the fusiform cell from p=5e-6 to p=0.027.) Pat_02 is *not* the load-bearing patient at α parsopercularis; it is a dilution noise contributor.

3. **Pat_03 dropout is a near-no-op.** Pat_03 contributes one α cortical trace-leaf (not at parsopercularis) and 2 contacts in Pat_03's right-hemisphere implant don't include parsopercularis. Drop changes K=30→29, N=866→770, enrichment 5.55→5.11, p 0.0014→0.0020. Headline survives.

4. **Pat_07 dropout collapses the cell.** Pat_07 carries 3/5 trace-leaves and 6 of the 26 cohort parsopercularis contacts (Y2..Y7). After drop: 2/20, 3.29×, p=0.121, **fails the n_pat_trace ≥ 2 filter** (only Pat_14 remains). The cell does not survive.

5. **Pat_14 dropout collapses the filter.** Pat_14 carries 2/5 trace-leaves and ~10 of the cohort parsopercularis contacts (R'1..R'10 left hemisphere). After drop: 3/16, 5.04×, p=0.019, **fails n_pat_trace ≥ 2** (only Pat_07 remains). The 3.4× enrichment with single-patient backing is suggestive but does not survive any cohort-level filter.

6. **Pro-cohort restriction (drop Pat_07 + Pat_15) collapses the cell.** This drop is dominated by Pat_07 dropout — the additional Pat_15 drop is a no-op at α parsopercularis (Pat_15 contributes zero α parsopercularis trace-leaves). 2/20, 3.20×, p=0.127, single-patient (Pat_14), fails filter.

### B.7 — Verdict and prose recommendation

The α ctx-lh-parsopercularis cell is:
- **Real at full cohort** (5.55×, p=0.0014, narrowly below Bonferroni threshold of 0.00104).
- **Two-patient (Pat_07, Pat_14), Pat_07-dominant** (3:2 split).
- **Fragile under any single-contributor dropout** — fails the n_pat_trace ≥ 2 filter under either Pat_07 drop or Pat_14 drop.
- **Strengthened under Pat_02 dropout** to Bonferroni-pass (3.4e-4, **0.0014×Pat_02 dilution becomes a Bonferroni-clean 3.4e-4**), because Pat_02's α implant has 6 parsopercularis contacts contributing zero trace-leaves and 4 α trace-leaves in other regions.
- **Anti-aligned-patient driven**: Pat_07, the largest contributor (3 leaves), is one of the two cross-probe anti-aligned patients of §3.5/§5.3/§5.4. This inverts the γ_l fusiform Pat_02 pattern, where the load-bearing patient is a strong cohort-pro-trace contributor.

The §5.5 prose should be edited to:
- Replace "two-patient base (Pat_02 and Pat_13)" → "two-patient base (Pat_07 and Pat_14), Pat_07-dominant".
- Disclose the Pat_07 anti-aligned status: "the cell is driven by Pat_07 (the cross-probe anti-aligned patient of §3.5), and is fragile under Pat_07 dropout (collapses to single-patient at p=0.12)".
- Optionally note the inverse Pat_02 effect: "removing Pat_02 sharpens the cell to Bonferroni significance (3.4e-4), opposite the γ_l fusiform pattern; Pat_02 here is a dilution contributor, not a load-bearing patient."

This downgrades the cell to "directional indication, two-patient, fragile, anti-aligned-patient-driven", parallel to the γ_l fusiform single-patient disclosure but in the opposite cohort-direction.

---

## C. Pat_07 / Pat_15 cortical trace-leaf counts (figure-caption check)

Per-band counts (after Wm/Unk drop):

### Pat_07

| band | total | cortical | Wm/Unk | cortical regions |
|---|---:|---:|---:|---|
| α | 14 | 5 | 9 | parsopercularis 3, rostralmiddlefrontal 1, superiorfrontal 1 |
| **β** | **9** | **1** | **8** | **parsopercularis 1** |
| low-γ | 10 | 1 | 9 | rostralmiddlefrontal 1 |

### Pat_15

| band | total | cortical | Wm/Unk | cortical regions |
|---|---:|---:|---:|---|
| α | 3 | 1 | 2 | rh-paracentral 1 |
| **β** | **9** | **3** | **6** | **rh-isthmuscingulate 1, rh-paracentral 1, rh-posteriorcingulate 1** |
| low-γ | 4 | 1 | 3 | rh-precentral 1 |

### Verification of the §5.5 prose claims

| §5.5 body claim | verdict |
|---|---|
| Pat_07 has 9 total β trace-leaves of which most are Wm/Unk | ✓ (9 total, 8 Wm/Unk, 1 cortical) |
| Pat_15 has 9 total β trace-leaves of which most are Wm/Unk | ✓ (9 total, 6 Wm/Unk, 3 cortical) |
| Pat_15 contributes "one [trace-leaf] per region" cortically | ✓ (β: 1 in each of isthmuscingulate, paracentral, posteriorcingulate. α: 1 region. low-γ: 1 region.) |
| Pat_07 contributes "single ctx-lh-parsopercularis contact" at β | ✓ (β cortical = 1 leaf at parsopercularis, channel Y3 leaf 82) |

**No corrections needed** for the figure-caption / body Pat_07 / Pat_15 claims at β. The prose's "sparse and regionally uncoordinated" framing matches the data (Pat_15 spreads 3 cortical β leaves over 3 distinct right-hemisphere cortical regions; Pat_07's single cortical β leaf sits at parsopercularis — the same region that drives the α cell).

**Side observation worth flagging in §5.5 prose (cosmetic):** Pat_07's *single* cortical β trace-leaf is at the *same region* (ctx-lh-parsopercularis) that anchors the α cohort enrichment cell. So Pat_07's β cortical contribution and α cohort contribution co-localize at parsopercularis. This is a cross-band Pat_07-internal consistency that supports the "Pat_07 has parsopercularis-localized lagged-coupling reorganization across α and β" reading, separate from the cohort-direction interpretation.

---

## Sufficient numbers for §5.5 prose corrections

- **Methods paragraph (Section A):** insert the 5-bullet definition above; cohort baselines reproduce as α 30/866=3.46%, β 46/866=5.31%, γ_l 75/866=8.66%.
- **α parsopercularis attribution (Section B):** Pat_07=3 + Pat_14=2 (NOT Pat_02 + Pat_13). Pat_07 is the load-bearing anti-aligned patient.
- **α parsopercularis dropout numbers:** drop Pat_07 → 2/20, p=0.121, single-patient (filter fail); drop Pat_14 → 3/16, p=0.019, single-patient (filter fail); pro-cohort → 2/20, p=0.127, single-patient (filter fail); drop Pat_02 → 5/20, 7.35×, p=3.4e-4, Bonferroni-pass; drop Pat_03 → 5/26, 5.11×, p=0.0020, headline survives.
- **Verdict:** downgrade α parsopercularis to "two-patient, Pat_07-dominant, anti-aligned-patient driven, fragile under any contributor dropout". Pat_02 dropout *strengthens* the cell rather than collapsing it — opposite of γ_l fusiform.
- **Pat_07 / Pat_15 figure-caption claims (Section C):** all verified, no corrections.
