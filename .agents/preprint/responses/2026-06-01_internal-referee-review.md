---
name: internal-referee-review-2026-06-01
era: IMCOH_ABS_COHORT_N10
status: open
kind: response
scope: referee-grade review of the locked two-probe preprint (cophenet rho_split + Grassmann T_G*); every load-bearing statistic recomputed from on-disk CSVs and each alleged defect adversarially refuted before inclusion
companion: 2026-06-01_review-issue-tracker.md, locked/VERDICT_LEDGER.md, locked/ANATOMY_LEDGER.md, METHODS_AUDIT_ISSUES.md
date: 2026-06-01
---

# Internal referee review — LRG functional-connectivity preprint (2026-06-01)

**Head.** The trace science is sound and survives independent recomputation: the
**β trace is real on both LRG probes, α on cophenet alone**, every load-bearing
p-value reproduces to machine precision, and the controls and outlier handling are
honest. The paper is **not** ready to draw LaTeX from, for one reason — the
**anatomy ("where") claim has been correctly retracted at the locked source-of-truth
(`ANATOMY_LEDGER`, 2026-05-30) but that retraction has not been cascaded into the five
band/cohort briefs that feed the prose**, and two headline framings ("band resolution
*not* amplification", "two *independent* probes") overstate what the numbers support.
Recommendation: **major revision — documentation cascade + framing, no new experiments
required** to defend the core finding. The contribution shifts honestly from *where*
the trace lives to *that* a band-specific multiscale trace persists and *at which
scales* — a thinner but defensible neuroscience claim.

This review recomputed all 17 checks in the `lapbrain` env against the canonical CSVs
and ran an adversarial refute pass on every alleged defect. The refute pass earned its
keep: it dissolved the one apparent blocker, falsified one tempting downgrade, and
correctly bounded the scariest bug. Net after refute: **0 surviving blockers, 4 major,
10 minor.** The companion issue tracker is the actionable list; this document is the
argument.

---

## 1. Summary verdict

| Layer | Status | One-line |
|---|---|---|
| **β trace (both probes)** | ✅ holds | cophenet p=0.0049 (7/10), Grassmann p_mass at floor with robust LOO; survives band-Bonferroni |
| **α trace (cophenet only)** | ✅ holds | cophenet p=0.0020 (5/10); Grassmann null — verdict correctly tagged "only D_coph" |
| **γ_l / δ Grassmann** | ⚠️ weak/fragile but honestly tagged | both single-probe; δ correctly demoted to "weak" (LOO Pat_08); γ_l rests on one probe |
| **Anatomy / localization** | ⛔ retracted at ledger, **not cascaded** | β/α cophenet DK lists not cohort-supported (max region 5/10); only Hippocampus@βGrassmann survives, and that is an *anchor* |
| **Headline framing** | ⚠️ overstated | "band resolution NOT amplification" and "two independent probes" both fail on the numbers |
| **Controls & outliers** | ✅ honest | n=10 everywhere, Pat_15 inclusion is conservative, δ/Pat_08 correctly handled |
| **Statistical faithfulness** | ✅ exact | every load-bearing number reproduced from CSV |

**Bottom line for the science:** nothing in this review touches a trace *verdict*. The β
finding is the headline and it is defensible. The required work is editorial-but-mandatory:
propagate the anatomy retraction and reword two slogans before any results prose is drawn.

---

## 2. Major points (must fix before drafting LaTeX)

These do not invalidate a result; they are load-bearing *documentation/framing* defects
that would propagate a false claim into the manuscript if the briefs are used as-is.

### M1 — The 2026-05-30 anatomy retraction is locked in `ANATOMY_LEDGER` but **not cascaded** into the 5 briefs that feed the prose. (issue B4)
`ANATOMY_LEDGER.md` (2026-05-30 revision) correctly retracts **β cophenet 0/7, α cophenet
0/11, γ_l Grassmann 0/7, δ Grassmann full 0/4 + epi-X 0/3**, and softens **β Grassmann to
Hippocampus only** — and explicitly lists its own "Cascade TODO (not yet done)" for
`bands/00_cohort.md`, `01_beta.md`, `02_alpha.md`, `03_gammalow.md`, `06_delta.md`,
`HANDOFF_INDEX.md`, `VERDICT_LEDGER.md`. **None of those edits have been made.** All five
briefs still assert the retracted anatomy as live fact (e.g. `00_cohort.md:29` "anatomically
localized to a band-specific distributed cortical network"; `01_beta.md:45` the 7+7-region
network; `02_alpha.md` the 11-region list + "9 of 11 reproduce under epi-X"). The decisive
fact behind the retraction, independently reproduced here to ~1e-16: **no DK region is
sampled by more than 5/10 patients**, and the "locked ∩ localized ∩ n≥5" count is **0/7 (β
cophenet), 0/11 (α cophenet), 1/7 (β Grassmann = Hippocampus)**. Do **not** change any trace
verdict — only the anatomy claims. The full file→line→edit map is in the tracker (B4).

### M2 — `audit_72` epi-X mask is a no-op bug; it taints δ Grassmann **anatomy** (and the bug is uncommitted). (issue B2)
`audit_72_anatomy_grassmann.py` (HEAD, commit 4871de1) masks epileptic nodes with
`~np.isin(np.arange(N), <string labels>)` — an int-vs-string comparison that excludes **zero**
nodes. Proof: `anatomy_delta_grassmann_epiX/per_patient_trace_flags.csv` shows `n_nodes_kept =
full N` for every patient (Pat_02 = 117, not 103; Pat_13 = 119, not 89), versus the correct
counts in `grassmann_epi_exclusion/epi_counts.csv`. So the locked **δ Grassmann epi-X "fully
disjoint networks"** claim (`bands/06_delta.md §5`, `00_cohort.md:143`) rests on a mask that
excluded nothing — it differs from the full-data network *only by its k-window*. This is the
**same `np.isin(int,str)` bug class fixed once in `audit_71` (commit e97d249) that recurred**.
Two scoping points the adversarial pass nailed and that must be stated in the errata so the
fix is not over-applied:
- **The C5 *trace* verdict is CLEAN.** The locked C5 cluster-mass numbers (δ p_mass^epiX=0.005,
  T_G*=43.99; β; γ_l) come from `audit_67_grassmann_epi_exclusion.py`, whose mask is the
  correct label-based `keep = [lbl not in epi_set for lbl in labels]` (`epi_counts.csv` proves
  real exclusions). The bug is confined to the per-region *anatomy* script.
- **The fix is in the working tree but uncommitted** (HEAD still ships the bug). Commit it,
  regenerate the two tainted anatomy dirs, then retract the δ "disjoint networks" prose.

### M3 — Grassmann anatomy is **phase-averaged participation (an anchor)**, not a cross-phase trace. (issue B3)
`audit_72`'s `participation_phase_avg()` averages per-node participation over (rest_pre, task,
rest_post) — verified bit-for-bit equal to `(pre+task+post)/3` on Pat_02 β. It is therefore
"where the leading-k Laplacian eigenmodes sit *on average*" — an **anchor** quantity — and is
nearly uncorrelated (node-level Pearson ρ ≈ 0.16) with a genuine differencing trace field. The
manuscript attaches this anatomy to the *cross-phase* T_G\* trace claim, implying the trace is
localized. Even the surviving "Hippocampus @ β Grassmann" must be framed as anchor-flavored
("where modes concentrate"), **not** trace localization. If genuine trace localization is
wanted, a differencing quantity (per-node participation change, or contribution to the
cross-phase d_G) must be built and A1/A3 re-run — do not reuse `participation_phase_avg`.

### M4 — "Band resolution NOT amplification" is overstated; the cophenet step also **amplifies within band**, strongest in the primary band. (issue D1)
The three-layer table reproduces exactly, but the slogan fails on the project's own
matched-strength gate (p<0.05): for **β**, raw FC p=0.0527 (**fails** the gate) → cophenet
p=0.0049 (**passes**) — a cross-gate *promotion* and >10× within-band strengthening; α
strengthens 0.0137 → 0.0020 (already passing). At the gate, raw FC detects only **{δ, α}**, so
"raw FC already detects everything" conflates the descriptive 6–8/10 `n_above` count with
gate-level detection. Only **δ** is a true across-band demotion (0.042 → 0.278); θ and γ_h were
never gate-detected, so "demoting" them is consistency, not resolution. Defensible reframe:
**"band-selective sharpening — across-band demotion of δ + within-band strengthening of β/α"**,
which is in fact the *stronger* methodological argument. Retire the "NOT amplification"
anti-pattern in `HANDOFF_INDEX.md`.

---

## 3. Minor points

- **m1 (A1)** — α "strong" is legitimate against the gate of record (C3 paired Wilcoxon p=0.00195,
  W=54/55), and the median-vs-p95 distinction *is* disclosed in the methods revision. Two cosmetic
  fixes: the CSV `verdict` column still says `intermediate` (stale pre-2026-05-19 8/10 scheme), and
  the column labeled `paired_wilcoxon_z` is actually the Wilcoxon **W** statistic, not a z-score.
  Surface α's 5/10-vs-β's-7/10 cohort-composition asymmetry in the ledger, not only in the methods table.
- **m2 (A2)** — β passes all four primary controls, but **C2 (split>drift, p=0.0137) is the marginal
  gate** — it clears 0.05 by a ~3× narrower margin than C1/C3 (both 0.00488). State this in limitations.
- **m3 (A3)** — `cluster_p_mass = 0.004975 = 1/(R+1)` is the empirical floor; β/δ/γ_l all sit at it
  (0/200 surrogates reach observed mass). The strong-vs-weak split rests entirely on the LOO column —
  **and that column is itself floor-quantized** (γ_l 8/201 vs δ 11/201), so "γ_l strong vs δ weak" hangs
  on a ~single-surrogate Monte-Carlo margin. State "p ≤ 0.005 at R=200" and **raise R (≥2000)** to both
  sharpen the bound and stabilize the LOO knife-edge.
- **m4 (A4)** — `grassmann_cluster_extent/cohort_summary.csv` `verdict_cluster_extent` still labels δ
  `strong` (pre-Decision-12); the locked ledger says `weak`. Regenerate the column or add a stale-warning.
- **m5 (A5)** — every three-layer cell reproduces exactly (all apparent mismatches are correct 2-sig-fig
  rounding); see M4 for the framing.
- **m6 (A6)** — multiple comparisons are handled honestly: cophenet β/α survive band-Bonferroni
  (0.05/6 = 0.0083) *and* Holm; Grassmann has within-band FWE (cluster-extent) and the cross-band BH(m=6)
  sentence was correctly deleted from the manuscript. Keep the explicit "no across-band correction"
  statement so three simultaneous 0.005 Grassmann bands aren't read as a corrected family.
- **m7 (C3)** — concrete instance of the anatomy defect: the β cophenet region `ctx-rh-insula` is supplied
  by **Pat_10 alone (8 contacts; 9/10 patients have 0)** — and Pat_10 is the **most β-anti** patient
  (obs_rho = −0.0906). An anti-trace patient certifies a "trace anatomy" region, a direction-blindness
  artifact of the absolute-value top-decile pipeline. Retracted at the ledger; still live in `01_beta.md:404`.
- **m8 (C4)** — one live retired-probe citation: `01_beta.md:528` calls γ_l's anatomy a "left-fusiform
  anatomical anchor", contradicting the locked retraction (fusiform appears at no band under S(b)). Strike it.
- **m9 (D2)** — "two **independent** probes" overstates: cophenet and Grassmann are two deterministic
  functions of **one eigendecomposition of the same Laplacian L̂ at τ_max**, so they share FC-estimation
  noise and are not statistically independent. The genuine support is the **double dissociation** (β both;
  α cophenet-only; γ_l/δ Grassmann-only; 3/6 bands disagree) → use **"complementary"**, not "independent".
  Surgical word swap at `01_beta.md:263,386,528`; leave the legitimate other uses (rsPre_A/B half-baselines).
- **m10 (D3)** — the residual message is coherent and verified, but the two cited source docs
  (`00_cohort.md`, `methods_neurophysiological_interpretation_2026-05-26.md`) predate and contradict the
  retraction (see M1). The neurophys companion's *body* is mechanism-level (no DK names) and salvageable.

---

## 4. Statistical faithfulness — recompute ledger

Every load-bearing number was independently re-derived from the on-disk CSVs in `lapbrain`.
**Result: 100% reproduce.** No fabricated or drifted statistics were found in the trace battery.

| Claim | Ledger | Recomputed | Match |
|---|---|---|---|
| β cophenet C3 | p=0.00488, 7/10, 24× | W=52/55 → p=0.0048828; 7/10 > p95; obs 0.2206 vs surr 0.0093 | ✅ |
| α cophenet C3 | p=0.00195, 5/10, 8.3× | W=54/55 → p=0.0019531; 5/10 > p95; obs 0.1054 vs surr 0.0128 | ✅ |
| β controls C1/C2/C4 | all pass | C1 0.00488, C2 0.0137 (marginal), C4 pass (non-degradation) | ✅ |
| Grassmann β/γ_l/δ p_mass | 0.005 (strong/strong/weak) | exactly 1/(R+1)=0.0049751; 0/200 nulls ≥ obs for all three | ✅ (floored) |
| Grassmann LOO | β 0.005, γ_l 0.040, δ 0.055 | β 0.00498 (Pat_02), γ_l 0.0398 (Pat_05), δ 0.0547 (Pat_08) | ✅ |
| Three-layer table | (6 bands × 3 layers) | every cell reproduces (rounding-consistent) | ✅ |
| band-Bonferroni | β/α survive | 0.00488 & 0.00195 < 0.0083; Holm agrees | ✅ |
| Anatomy localization | β 0/7, α 0/11, βGr 1/7 | reproduced to ~1e-16; max region n=5/10 | ✅ |

The desyncs are in *labels*, not *values*: δ Grassmann `verdict` column (m4), α `verdict`
column + mislabeled `z`/`W` (m1). The arithmetic is trustworthy throughout.

---

## 5. Anatomy & control integrity

The anatomy story is the one place the manuscript currently *over-claims in its consumer
docs*, and it is the spine of M1–M3. Three independent problems compound:

1. **Undersampling (B1, reproduced):** at this implant coverage no DK region reaches >5/10
   patients; the locked DK lists rest on 1–4 patients each. Cohort localization is simply not
   resolvable at n=10 — neither "distributed network" nor "diffuse brain-wide" is supportable.
2. **Direction-blindness (C3):** the absolute-value top-decile flagging pools trace and
   anti-trace pairs, which is how a β-*anti* patient (Pat_10) ends up certifying β "trace
   anatomy" (insula).
3. **Wrong quantity (B3):** the Grassmann anatomy is phase-averaged (anchor), not cross-phase
   (trace), so even the well-sampled survivor (Hippocampus) localizes the wrong thing.

The 2026-05-30 audit and the locked `ANATOMY_LEDGER` already encode (1)–(3) correctly; the
manuscript fails only at propagation (M1). The **trace controls (C1–C5) are intact** — the
epi-X bug (M2) is confined to the anatomy script and does **not** reach the C5 trace gates.

---

## 6. Outliers — handled honestly

- **Pat_15** (right-hemisphere-only, 0 epi, β-anti): included at **n=10 everywhere**; the n=9
  drop-Pat_15 narrative is vestigial (KC-era) and correctly confined to a supplementary
  sensitivity. Pat_15 is a *strong anti* contributor (drop-Pat_15 LOO mass jumps to β 184 / γ_l
  151 vs ~66 full), so **including her makes the n=10 β verdict conservative** — a point to make
  explicitly in the manuscript's favor.
- **Pat_08** (height-only / T:R=0 KC outlier): the **δ** verdict correctly demotes to "weak"
  because dropping Pat_08 crosses the LOO gate (0.0547). The tempting further downgrade to "no
  trace / single-patient-leveraged" is **wrong** — all 10 δ LOO masses exceed the null p95, so δ
  is a genuine (weak) cohort signal, not a Pat_08 artifact.
- **Pat_10** (β-anti): the only patient implanted in rh-insula, and the sole certifier of that
  retracted region (m7/C3).
- **Pat_03** (1024 Hz): correctly absorbed at the config layer; no analysis-level distinction. ✅

---

## 7. Neurophysiological message — is it worth a paper?

With anatomy retracted, the contribution moves from **WHERE** to **THAT + AT WHICH SCALES**.
Strongest defensible sentence (CSV-supported):

> *Under matched-strength surrogacy a band-specific multiscale structural trace of the task
> persists into rest_post: β carries it at every rung (edge FC, the dendrogram merge-order
> D_coph, and the principal communication subspace at coarse scales k=27..55), α only at the
> per-pair hierarchical level, and γ_l only at the global subspace level, while θ and γ_h are
> null — but at n=10 with this implant coverage the trace cannot be pinned to specific cortical
> regions.*

**Is that worth a neuroscience paper?** Conditionally yes — but the contribution is now thinner
than a regional-localization headline, and its novelty rests entirely on the two LRG probes being
read as **genuinely multiscale, complementary objects** rather than re-statements of edge-level
FC. The single most compelling, genuinely novel and CSV-supported residual finding is the
**band × mechanism dissociation** (the 2×2 probe map: β both / α cophenet-only / γ_l·δ
Grassmann-only / θ·γ_h null). Lead with that. The biggest remaining threat to the contribution is
that the non-β cells are each fragile (δ weak/Pat_08-leveraged; γ_l single-probe; all Grassmann
epi-X "secondary reassurances" need M2 cleared before they can be cited) — so the paper is
effectively a **β paper with a supporting band-dissociation map**, and should be scoped as such.

---

## 8. Not covered here — recommended further checks (from the completeness critic)

These are outside the 17-check scope but a referee will raise them; ranked by leverage:

1. **`imcoh_abs` substrate validity (foundation under *both* probes).** No check examined whether
   `<|ImCoh|>_f = mean(|signed|)` is a sound substrate for a combinatorial Laplacian — the
   Jensen order-of-operations (`|mean|` vs `mean|·|`), whether abs-rectification preserves Nolte
   volume-conduction immunity after band-averaging, and sign/lag destruction. This is the bedrock.
2. **τ_max = 1/λ_max single-scale choice.** Both probes live entirely at the finest scale; λ_max
   is dominated by the largest-degree node (probe geometry can inflate it). `memory/lrg_tau_choice`
   flags this as the open methods question. A referee will ask whether the β trace survives at other τ.
3. **Figure scripts (`scripts/02_preprint/`, 18 scripts) — unaudited.** `preprint_10_bands_brain_anatomy`
   / `preprint_11_beta_anatomy_brain` likely render the **now-retracted** DK region-lists; also check
   the never/always plot rules (no `axhline` threshold lines, RdBu_r polarity, n/10 reference lines).
4. **α C5 epi-X possibly bug-adjacent.** The α "8.3× → 27.7× strengthening under epi-exclusion"
   that `00_cohort.md` cites should be re-checked: confirm the `audit_72` int-vs-str mask bug (M2)
   did **not** also corrupt the α cophenet epi-X path (`alpha_epi_exclusion/`).
5. **γ_l "strong, only Grassmann" not independently stress-tested** (single probe; C5 epi-X LOO
   fragile at 0.159, Pat_05). **γ_h** no-trace sits on the 0.05 boundary (0.0597 mass) with
   inconsistent quotes (0.055/0.060) across docs — one permutation could flip it.
6. **`METHODS_AUDIT_ISSUES.md` open ledger** (separate Q1–Q4/N1–N4/etc.) carries pre-existing
   desyncs (e.g. Q2: `00_cohort.md §3` cluster-extent table holds pre-fix masses; Q3/Q4: γ_l/δ
   "weak" vs "strong" drift) — reconcile in the same pass.
7. **raw-FC substrate (Result-1)** robustness and the d_S/d_P ρ≈0.85–0.95 correlation were not
   re-audited; the "cophenet sharpens a borderline β substrate" narrative depends on that layer.

A separate **literature-positioning memo** (β post-task persistence / hippocampal involvement /
consolidation / iEEG band-connectivity priors) is queued as the next deliverable.

---

## 9. Recommendation

**Major revision.** No experiment is required to defend the core β (and secondary α) trace
finding — every load-bearing statistic is faithful and the controls/outliers are honest. Before
any results prose is drawn, the authors must: **(M1)** cascade the anatomy retraction into the
five briefs + `HANDOFF_INDEX` + `VERDICT_LEDGER` anatomy refs; **(M2)** commit the epi-X mask
fix, regenerate the two δ anatomy dirs, and errata-note that the C5 trace gate was never tainted;
**(M3)** reframe all Grassmann anatomy as anchor-flavored phase-averaged participation; **(M4)**
replace "band resolution NOT amplification" with "band-selective sharpening". The minors are
low-cost label/wording fixes. Done, the paper is an honest, defensible **β multiscale-trace
paper with a band×mechanism dissociation map** — worth publishing, provided §8.1–8.2
(substrate + τ) are addressed in a methods/limitations pass.
