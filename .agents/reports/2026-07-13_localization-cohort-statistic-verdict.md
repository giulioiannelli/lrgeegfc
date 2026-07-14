---
name: localization-cohort-statistic-verdict
kind: verdict
era: IMCOH_ABS_COHORT_N10 (mst@0.20 recovery)
status: current
created: 2026-07-13
scope: Resolves whether the mst@0.20 "localization is dead / β delocalized" verdict was a
  testing artifact. Answer — PART artifact (the Wilcoxon cohort test was DOA at n=5 and hid
  a real low_γ localization) but the β localization is GENUINELY delocalized, now on a valid
  powered statistic. The only robust localization is low_γ (disease band).
pointers:
  - .agents/guides/task-persistence-investigation/2026-07-13_localization-cohort-statistic-audit.md  # scope + method
  - scripts/01_compute/sparsified_arc/23_localization_statistic_bakeoff_mst020.py                     # C1 bake-off
  - scripts/01_compute/sparsified_arc/24_localization_rawcoph_bakeoff_mst020.py                       # C4 raw-vs-coph
  - src/lrg_eegfc/utils/metrics/cohort_localization.py                                                # 4 cohort statistics
  - data/sparsified_arc/localization_bakeoff_mst020/                                                  # C1 outputs
---

# Localization cohort-statistic verdict (mst@0.20)

## Head

The user suspected the scale-dependent framework's "localization washed out" was a testing
artifact. **It is — partly — but fixing it does not bring β→OFC back.** The mst@0.20
localization scripts (16/17) used a per-patient Wilcoxon signed-rank across the 5–8 patients
that sample each region; its one-sided floor at n=5 is p=0.031, so a region sampled in 5
patients (OFC, insula, MTL) **cannot clear BH regardless of effect size** — a genuinely
broken test. But a four-statistic bake-off on the *identical* backbone + matched-strength
null (holding everything fixed, varying only the cohort statistic) shows that the only band
beating the null under a *powered, consistency-favouring, whole-grid-BH, LOO-robust* test is
**`low_gamma`** — but its concentration is **DIFFUSE** (8/9 systems positive, top-ranked
region is a 3-patient coverage artifact), NOT a focal anatomical home; the one *specific*
low_γ result is **low_γ → SOZ** (disease marker). **β (trace, encoding, inference → OFC,
cingulate, anywhere) is genuinely delocalized** — dead under all four statistics, at the 2%
surrogate-beating noise floor. So the delocalization verdict for β stands, now for the right
reason, and no band gives a focal *cognitive* localization on the backbone.

## Method (what was actually run)

`cohort_localization.py` computes four cohort statistics on identical per-patient
(obs, matched-strength-surrogate-ensemble) cells:
- **(a) pooled** — cohort-median-of-obs vs R surrogate-medians (floor 1/(R+1), the whole-graph-era statistic).
- **(b) wilcoxon** — per-patient (obs−surr) signed-rank (floor 1/2^K, the scripts-16/17 statistic; DOA at K=5).
- **(c) stouffer** — fixed-effect combine of per-patient z (powered, but omnibus / 1-patient-drivable).
- **(d) sign-consistency** — cohort median-of-z with a matched-strength permutation null (powered AND
  median-robust → consistency-favouring). Reported with `frac_pos`.

Gate reported (no single committed gate — user decision): sign-consistency **whole-grid BH**
(family = all systems × scales × targets within a band × grouping) **+ LOO-worst < 0.05**.
Substrate + null IDENTICAL to script 17 (mst@0.20, matched-strength R=200). Three granularities
(system, supersystem=limbic K=10, soz). Verified: independent recompute == library.

## Findings

### 1. The wash-out was PART test-artifact
Under the powered statistics, localization is NOT empty — the Wilcoxon's "nothing anywhere"
was under-powered. Recovered cells (sign-consistency whole-grid q<0.05, LOO-robust) exist.

### 2. …but everything recovered is `low_gamma`, nothing is β/α/δ
Floor audit: **low_γ = 66/432 (15%)** system cells beat the *entire* matched-strength null;
**β/α/δ = 8–10/432 (~2%)** = pure noise floor. β localization is at the noise level in every
region and target. **β → OFC / cingulate / any system: DEAD under all four statistics + LOO.**
The pre-registered β→OFC, encoding→OFC, inference→cingulate all fail. (The whole-graph β→OFC was
itself LOO-fragile per audit_162 — nothing lost.)

### 3. The `low_gamma` signal is DIFFUSE, not a focal anatomical home
low_γ beats matched-strength but its per-system concentration profile is **broad, not
focal**: 8/9 systems have positive Θ and clear (trace: insula +3.7, MTL +2.8, occipital
+2.8, PFC +2.6, sensorimotor, OFC, cingulate, parietal all +). Two red flags against a
"focal home" reading: (a) **`occipital` (K=3, 3 patients only) sits at/near the top for all
three targets** — a coverage artifact; (b) the top systems form no coherent story (insula /
occipital / MTL / PFC / parietal). So low_γ is NOT delocalized the way β is (β = 2% noise
floor), but it does **not localize to a specific region** either — it is a **diffuse** low_γ
cross-phase concordance field. **The one specific, interpretable low_γ concentration is
`low_γ → SOZ`** (seizure-onset zone; clears trace/enc/inf q≈0.04, while β→SOZ does NOT) —
the *disease* marker (binary tissue split, matches the epi line AUC .82), not an anatomical
system. Correction of an earlier over-generous "low_γ → PFC/limbic" phrasing.

### 4. Statistic disagreement is diagnostic (validates the 4-way design)
Stouffer (omnibus) fires for β-PFC (STOU q≈0.002) while the median-based pooled + sign-
consistency + LOO do not → those β "hits" are single-patient-driven artifacts. The trustworthy
(median) statistics agree, and only in low_γ.

### 5. C3 — duration (NON-DESTRUCTIVE, no truncation)
Every inference/encoding localization cell is duration-clean: Spearman(test/learn ratio,
per-patient loc) is ns (|ρ|≤0.5, p≥0.25). Truncation null retired-invalid; ratio regression used.

### 6. C4 — raw-vs-cophenetic value-add
The cophenetic localization value-add is REAL but exists **only in low_γ**: 16 cells where coph
concentrates OFC/PFC/cingulate/sensorimotor/limbic while raw edges are placeless; + 26 raw-only
cells (raw localizes elsewhere too). **Zero β value-add** (β delocalized in both representations).

### 7. C5 — backbone-density robustness (frac 0.10 / 0.20 / 0.50) — COMPLETE
The pruning hypothesis is **ruled out** and the band dissociation is **density-invariant**.
LOO-robust survivor cells per band:

| frac | δ | α | β | low_γ |
|------|---|---|---|-------|
| 0.10 | 0 | 0 | 0 | 32 |
| 0.20 | 0 | 0 | 0 | 49 |
| 0.50 | 0 | 0 | 0 | 56 |

**β = exactly 0 at every density** (β→OFC does not reappear even at frac=0.50, closest to the
dense graph where it originally showed) → β delocalization is not a frac=0.20 pruning artifact.
low_γ is present at all densities and strengthens with density. Caveat: at frac=0.50 low_γ
survivors span nearly all systems (localization broadens toward a distributed/global low_γ
concordance); the focal reading (PFC/limbic/MTL/sensorimotor) is sharpest at frac=0.10–0.20.

### 8. β trace is SOZ-INDEPENDENT (cognitive, not epileptic) — 3 ways
The β trace must be a cognitive result, not one riding on epileptic tissue. Confirmed three
independent ways:
- **Pair-class** (script 19): β trace on non-SOZ pairs (`nonepi_nonepi`) clears MS **16/16
  scales** (q<.05 15/16, ~5400 pairs); in the SOZ core (`epi_epi`) only 1/16 (~55 pairs) —
  β **spares** the seizure zone.
- **Node-exclusion** (script 25): physically drop SOZ contacts → rebuild mst@0.20 → recompute
  under MS: **β 14/16 scales clear** (best p=.002, ρ=.324) vs 16/16 full — essentially intact.
- **Decimation control** (10 seeds, size-matched random non-SOZ drop): SOZ-excluded β ρ sits
  *within* (or *above*, at s≈23) the random-drop band at every scale — SOZ removal is no more
  disruptive than equivalent random node removal.
**Mirror image: low_γ trace does NOT survive SOZ removal (node-excl 0/16; nonepi_nonepi 0/16).**
The low_γ cross-phase structure IS the SOZ tissue. Clean dissociation: **β = SOZ-independent
cognitive trace; low_γ = SOZ/disease.** Caveat: low_γ→SOZ may be an *anchor* (stably
distinctive SOZ low_γ FC) rather than a task-*trace* — MS controls strength, not tissue
stability; flag for a later trace-vs-anchor check.

### 9. Scale-matched parcellation ladder (C7) — β localizes COARSELY (left hemisphere); α placeless
The fixed 9-system atlas was a granularity mismatch: the LRG diffusion at scale s coarse-grains
space, so a coarse-scale trace must be tested at a coarse parcellation. Re-run over a coarse→fine
ladder (hemisphere / lobe / supersystem / system; `SA_GRANLADDER=1`, data
`localization_bakeoff_mst020_glad/`):
- **β: LEFT hemisphere ENRICHED** (z=+1.38, pooled p=.005, sign-consistency p=.005, **LOO=.005**,
  5/7 patients, peak meso scale s≈11–16); **RIGHT hemisphere DEPLETED** (negative z all scales).
  This is the laterality (ρ=.685, pillar 6 / audit_148) surfacing in the localization frame.
  **β is NOT placeless — it has a COARSE (hemisphere) home (left), just no fine one.** Marginal
  under strict whole-grid BH (q=.119 — only 2 hemisphere units) but robust under pooled + LOO, and
  left-lateralization is an a-priori hypothesis.
- **α: NO localization at ANY granularity** (only the K=3 occipital coverage artifact). Genuinely
  placeless at every resolution — the "α localizes at its meso scale" prediction FAILS.
- **low_γ: FRONTAL / PFC** (lobe=frontal K=8 q=.030; system=PFC K=7 q=.033; LOO-robust) — the
  cleanest low_γ home; given low_γ→SOZ + vanishes-without-SOZ, this frontal/PFC is very likely the
  seizure zone (OPEN: confirm SOZ is frontal-dominant in this cohort).
- **The clean scale↔granularity DIAGONAL is NOT confirmed** — β's hemisphere signal sits at *meso*
  scale (s≈11), not the coarsest; a crude K(s) proxy hints s≈11≈K2 (hemisphere-matched) but is
  unreliable, so not leaned on. The ladder is a good *lens* (revealed β's coarse home) but scale and
  localization-granularity do not align in a simple diagonal.

**Corrected per-band picture (supersedes "all delocalized"):**
| band | localization | reading |
|---|---|---|
| **β** | **left hemisphere** (coarse), meso scale; no fine home | cognitive, SOZ-independent, left-lateralized |
| **α** | none at any resolution | genuinely placeless |
| **low_γ** | frontal/PFC ≈ SOZ | disease band, not a cognitive trace |
| δ | none clean | — |

## Implication for the paper

- β→OFC / β→cingulate / encoding→OFC / inference→cingulate stay **retired** (no fine anatomical
  home; broken-statistic verdicts both ways). The paper must NOT reinstate a fine β localization.
- **REFRAME (not "placeless"):** β localizes **coarsely — to the LEFT hemisphere** (=the laterality),
  and is **SOZ-independent** (holds on non-SOZ tissue, survives SOZ node removal). β = a
  left-lateralized cognitive trace with a coarse home, no fine one.
- **α** is genuinely placeless at every resolution.
- **low_γ** is the disease band: diffuse but **frontal/PFC ≈ SOZ**-dominant; its trace **vanishes
  without the SOZ**. This is the anatomical signal the broken Wilcoxon hid (the "localization
  value-add" is low_γ-only, NOT a β cognitive value-add — soften audit_171 framing in directives).
- Never say "localization entirely dead / everything delocalized" — state **per band**.

## NEXT (post-compact) — CRYSTALLIZE into current results (the agreed next step)

Fold this verdict into the live results so agents read one consistent story:
1. **β localization**: everywhere it says "delocalized / placeless / β→OFC" → "**coarsely
   LEFT-LATERALIZED (hemisphere), SOZ-independent; no fine home**". Targets: settled-doc R1.8,
   `sec1` directive, talk slide 16, the R1/R2 cascade tables, and the 4 MEMORY figure-lines flagged
   `⚠ at risk (β→OFC)`.
2. **low_γ**: NEW — add the frontal/PFC≈SOZ disease localization to the epi / band-taxonomy lines;
   note it's what the DOA Wilcoxon hid. Carry the **anchor-vs-trace caveat** (low_γ→SOZ may be stable
   SOZ tissue, not a task-trace).
3. **α**: placeless at every resolution (reinforces "α no home").
4. **"localization value-add" claim** (audit_171 + directives): it's **low_γ-only**, not a β
   cognitive value-add — soften/scope it.
5. **Method LOCK to record**: cohort localization MUST use a resolution-preserving + consistency +
   whole-grid-BH + LOO statistic (`utils/metrics/cohort_localization.py`), NEVER a per-patient
   Wilcoxon at K=5 (DOA); and test COARSE parcellations too (fixed atlas = granularity mismatch).
6. **Open / not-yet-done**: confirm SOZ-is-frontal for low_γ; build the scale×granularity ladder
   FIGURE (script `fig_localization_statistic_bakeoff.py` exists for the system bake-off); **COMMIT**
   everything (all UNCOMMITTED, new branch off `audit/cohort-n10-diagnostic`).

## Reproduce

```
PY=/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 $PY \
  scripts/01_compute/sparsified_arc/23_localization_statistic_bakeoff_mst020.py --stage ab   # C1 (~7min)
$PY scripts/01_compute/sparsified_arc/24_localization_rawcoph_bakeoff_mst020.py --stage ab   # C4 (~2min)
# C5: SA_FRAC=0.50 / 0.10 python 23_...py --stage ab
```
