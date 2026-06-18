---
name: rhocoph-grassmann-dissociation
type: report
era: IMCOH_ABS
status: current
created: 2026-06-05
updated: 2026-06-05
pointers:
  - .agents/preprint/locked/VERDICT_LEDGER.md
  - .agents/reports/2026-06-01_trace-concordance-vs-blind-fc.md
  - .agents/reports/2026-05-30_pair-trace-measures-methodology-audit.md
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
  - data/audit/grassmann_cluster_extent/cohort_summary.csv
  - data/audit/epi_stratified/cophenetic_cohort.csv
  - scripts/02_preprint/preprint_34_bands_coph_grassmann_dissociation_map.py
---

# Where the per-pair cophenetic trace and the Grassmann subspace trace dissociate

**Head.** Across the six bands the two LRG trace probes decompose into four clean
cells, and the scientifically interesting output is *where they disagree*. **β**
is the only band with **both** a per-pair cophenetic trace (the novel backbone
probe: do pairs keep their relative communication-distance position task→rest?)
**and** a Grassmann subspace trace (do the dominant collective modes rotate?).
**α is the flagship dissociation — COPH-ONLY**: a robust per-pair trace
(p = 0.0020) with **no** trace in the dominant collective modes
(Grassmann p_mass = 0.0995, gate p = 0.348). The task reorganizes fine
relational structure that persists into rest *without* reorganizing the
dominant collective pattern. **δ and γ_l (on the full graph) are
GRASSMANN-ONLY** (global-mode trace, no per-pair trace), δ LOO-fragile. **θ
and γ_h are NEITHER.** γ_l's classification is **conditional on the epileptic
confound**: on the full graph it is Grassmann-only, but its cophenetic trace
**emerges** under epileptic-node exclusion (full p = 0.116 → exclude_epi
p = 0.024, 8/10) while the Grassmann signal simultaneously *contracts* under
the same exclusion — so γ_l is global-only at the whole brain and
**conditionally both** once epileptic contacts are removed.

This is an aggregation of locked matched-strength outputs (R = 200, 4-cycle
±δ, split-baseline). **No surrogate was recomputed.** Every number is
verified against `VERDICT_LEDGER.md`; the one mismatch found (δ cophenet
ratio) is flagged in §5.

---

## 1. The 6-band dissociation table

x-probe = **ρ^coph** (per-pair cophenetic, the backbone): one-sided paired
Wilcoxon over matched-strength surrogates, `ρ_split^coph = Spearman(ΔD_coph^task,
ΔD_coph^rest)`. y-probe = **Grassmann** `T_G*` (leading-subspace
span{φ_2..φ_{k+1}} chordal rotation): cluster-mass permutation gate over the
same R = 200 surrogates.

| Band | ρ^coph median | ρ^coph ratio | n/10 | ρ^coph p | ρ^coph verdict | Grassmann p_mass | LR (k-cells) | Grassmann LOO-max | Grassmann verdict | **Dissociation cell** |
|---|---|---|---|---|---|---|---|---|---|---|
| **β** (13–30) | +0.221 | 23.7× | 7/10 | **0.0049** | **trace** | **0.005** | 29 | 0.005 | **strong trace** | **BOTH** |
| **α** (8–13) | +0.105 | 8.25× | 5/10 | **0.0020** | **trace** | 0.0995 (gate 0.348) | 4 | 0.577 | no trace | **COPH-ONLY** (flagship) |
| **γ_l** (30–80) | +0.083 | 30.2× | 5/10 | 0.116 | no trace | **0.005** | 13 | 0.040 | **strong trace** | **GRASSMANN-ONLY** (full) → conditional BOTH (epi-X) |
| **δ** (0.53–4) | +0.008 | 0.98× (see §5) | 4/10 | 0.278 | no trace | **0.005** | 7 | **0.055** | **weak trace** (LOO-fragile) | **GRASSMANN-ONLY** |
| **γ_h** (80–300) | +0.0003 | 0.082× | 4/10 | 0.246 | no trace | 0.060 | 9 | 0.164 | no trace | **NEITHER** |
| **θ** (4–8) | −0.040 | −10.2× (anti) | 2/10 | 0.722 | no trace | 0.159 | 5 | 0.826 | no trace | **NEITHER** |

Notes on the gate columns:
- ρ^coph verdict = `p < 0.05` AND `median > 0`. The `n/10` (`n_above_surrogate`)
  and the ratio are **descriptive**; the Wilcoxon p is the gate (locked rule).
- Grassmann verdict follows the locked Decision-12 tier: strong needs
  `p_mass < 0.01` AND `LR_p < 0.05` AND `LOO-max < 0.05`. δ clears the cohort
  gate at the empirical floor (`p_mass = 0.005`) but fails the LOO precondition
  (LOO-max = 0.055, Pat_08) → **weak**, not strong. This is marked in the figure
  with a crimson marker edge.
- Grassmann `p_mass = 0.005` is the empirical-null floor `1/(R+1)` for R = 200;
  β, γ_l, δ all hit it, so on the figure their y-coordinates coincide at the
  ceiling −log10(0.005) = 2.305. The LOO-max column is what discriminates them.

### γ_l — conditional dissociation (the epileptic confound)

| γ_l view | ρ^coph p | ρ^coph median | n/10 | Grassmann |
|---|---|---|---|---|
| full graph | 0.116 (no trace) | +0.083 | 5/10 | strong (p_mass 0.005, LR 13) |
| exclude_epi | **0.024 (trace EMERGES)** | +0.174 | **8/10** | **contracts**: p_mass^epi-X 0.030, mass 66→33, LR 13→10, LOO-max^epi-X 0.159 (fragile) |

The two probes move in **opposite directions** under epileptic-node exclusion:
the per-pair cophenetic trace *emerges* (the epileptic subgraph was masking a
non-epileptic-cortex per-pair trace), while the Grassmann subspace trace
*weakens* (the full-graph collective-mode rotation was partly carried by
epileptic contacts). So γ_l is **global-only on the whole brain, conditionally
both once epileptic contacts are removed** — characterized explicitly rather
than collapsed to one label.

---

## 2. Neurophysiological frame (the interpretive spine)

- **ρ^coph** asks whether each **pair** keeps its relative ultrametric
  (communication-distance) position across phases. Positive = the fine, local
  *relational* fabric of the hierarchy is reorganized by the task and the
  reorganization persists. This is the per-pair, multiscale-aware backbone the
  paper introduces.
- **Grassmann** asks whether the **dominant collective modes** (the leading
  eigenvector subspace) rotate. Positive = coarse, global reorganization of the
  network's principal collective pattern.
- **α = COPH-ONLY** is therefore a genuine, interpretable finding, not a
  discrepancy to hide: *the task reorganizes fine relational structure that
  persists into rest, without reorganizing the dominant collective pattern.*
  Locally (per-pair) there is a trace; the global/dominant modes show none.
- **β = BOTH**: the reorganization reaches all the way up — both the fine
  relational fabric and the dominant collective modes carry the persistent
  imprint. This is why β is the manuscript backbone.
- **δ / γ_l (full) = GRASSMANN-ONLY**: the dominant modes rotate but the
  per-pair relational ranks do not move coherently — a coarse-only
  reorganization.

---

## 3. The per-patient picture — "nested yet dissociating"

Source: `2026-06-01_trace-concordance-vs-blind-fc.md` §2–3 (per-patient
concordance cross-tab, integrity-gated against `beta_per_patient.tex`).

Per-patient classification (B = both probes positive, C = cophenet only,
G = Grassmann only, · = neither; `z ≥ 1.96` for ρ^coph, `n_sig_k ≥ 20` for
the per-patient Grassmann flag — both verbatim manuscript rules):

```
band        02 03 05 06 07 08 10 13 14 15   both onlyC onlyG none
delta        B  B  G  B  B  G  G  G  .  .     4    0     4    2
theta        G  G  B  B  G  G  G  G  G  .     2    0     7    1
alpha        G  B  G  B  G  B  B  G  B  .     5    0     4    1
beta         B  B  B  B  B  B  G  B  G  .     7    0     2    1
low_gamma    B  B  B  B  .  B  G  G  G  .     5    0     3    2
high_gamma   G  B  B  B  .  B  G  G  G  .     4    0     4    2
```

### What "nested yet dissociating" means precisely

Two facts that look contradictory but are not:

1. **Nested at the per-patient level: `onlyC = 0` in *every* band.** There is
   no patient who is cophenet-positive but Grassmann-negative. The per-patient
   Grassmann flag is a strict **superset** of the per-patient cophenet flag:
   whenever a patient's per-pair merge-height ranks move coherently, that
   patient's leading subspace has *always already* rotated. So at the level of
   "is there *any* signal in this patient", the subspace probe is the more
   sensitive of the two and the cophenet probe never fires alone.

2. **Dissociating at the band/cohort level.** The *cohort gates* are different
   objects and disagree per band: ρ^coph fires at α (and β), Grassmann fires at
   δ, γ_l (and β). The two co-fire on **β alone**. The dissociation lives in
   *what each cohort gate aggregates*, not in a per-patient orthogonality.

The reconciliation: the **per-patient Grassmann flag is permissive** (`n_sig_k
≥ 20` fires for 8–9/10 patients in *every* band, including θ and α where the
cohort Grassmann cluster-extent gate says no trace). What is band-selective is
the **cohort-level cluster extent** (count + contiguity of significant k-cells
over matched-strength), *not* the per-patient flag. So "nested" is a statement
about the permissive per-patient flags (G ⊇ C), and "dissociating" is a
statement about the demanding cohort gates (different bands clear cophenet vs
Grassmann). Both are true simultaneously; they describe two different layers.

### β single-patient breakdown — Pat_10 and Pat_14 carry β purely in subspace

At β: Grassmann fires for **9/10** but ρ^coph for **7/10**. The two
patients who split:

| patient | ρ^coph (z) | Grassmann T_G* (n_sig_k) | class |
|---|---|---|---|
| Pat_10 | −0.091 (z −1.70) ✗ | 0.612 (71 cells) ✓ | only Grassmann |
| Pat_14 | −0.049 (z −1.00) ✗ | 0.695 (78 cells) ✓ | only Grassmann |

Pat_10 and Pat_14 carry the β reorganization **entirely in the leading-
eigenmode subspace rotation** (z ≈ +5 on T_G* at audit_66, ~100% of the
k-window above their own surrogate) while being flat/anti on the per-pair
merge-height distance. A real geometric dissociation, not noise: "where is the
β trace in Pat_10?" has an answer (in the slow-mode subspace) under Grassmann,
and is a genuine dissent under cophenet. The lone complete β dissenter is
**Pat_15** (right-hemisphere-only implant, the pre-registered biology
exception).

---

## 4. The figure

`scripts/02_preprint/preprint_34_bands_coph_grassmann_dissociation_map.py`
→ `data/preprint/figures/all_bands/fig_bands_coph_grassmann_dissociation_map.pdf`

A 2-D dissociation plane, one point per band:
- **x** = signed cophenetic effect `sign(median) · (−log10 p_Wilcoxon)`.
- **y** = Grassmann effect `−log10 p_mass`.
- Four quadrants tinted + labeled: neither / local-only (coph) / global-only
  (Grassmann) / both. p = 0.05 guides on both axes (vertical at +1.301 for
  cophenet, horizontal at 1.301 for Grassmann); a x = 0 dotted line separates
  trace from anti-trace direction.
- Points coloured by cell; **crimson marker edge = Grassmann LOO-fragile** (δ).
- **γ_l plotted twice**: full-graph circle (global-only) → epi-excluded diamond
  (crosses into both), joined by a dashed arrow — the conditional emergence.

PDF-only, full vector (no raster XObject, verified), `use_lrg_style()`,
transparent, no suptitle, no in-axes numeric text (the numeric table prints to
stdout), single-letter panel tag. 23 KB, 1 page.

---

## 5. Honest caveats

**(a) The "nested yet dissociating" story is partly an artefact of differing
per-measure gates.** The cophenet probe is gated by a **paired one-sided
Wilcoxon on the per-patient ρ_split** (a per-patient scalar), whereas Grassmann
is gated by a **cluster-extent / cluster-mass FWE permutation across the k-grid**
(a curve-level statistic). These are not the same test on the same footing:
- The per-patient Grassmann flag (`n_sig_k ≥ 20`) used in §3 is **far more
  permissive** than the cohort cluster-mass gate used in §1, which is why
  `onlyC = 0` (every cophenet-positive patient trivially clears the permissive
  Grassmann per-patient bar). The strict superset relation G ⊇ C is therefore
  *partly mechanical* — a consequence of the permissive per-patient Grassmann
  rule, not proof that the subspace is intrinsically the more sensitive
  geometry. The cohort-level dissociation (different bands clear cophenet vs
  Grassmann) is the more defensible statement; the per-patient nesting should
  be read as "the permissive subspace flag never undershoots the cophenet flag",
  not as a deep ordering of the two geometries.
- A fully symmetric comparison would gate both probes with the same statistic
  (e.g. both at cohort level, or both per-patient at matched thresholds). That
  has not been run; the table mixes a per-patient-scalar Wilcoxon (cophenet)
  with a curve-FWE permutation (Grassmann). The cell assignments are correct
  under the *locked* gates, but the cross-probe sensitivity comparison is not
  apples-to-apples and should not be over-read.

**(b) δ ratio mismatch vs the ledger (flagged, immaterial to verdict).** The
ledger δ cophenet line reports "ratio 1.69×"; the cohort_summary.csv gives
ratio-of-cohort-medians = obs 0.00765 / surr 0.00782 = **0.98×**. For all other
bands the cohort-CSV ratio matches the ledger exactly (α 8.25×, β 23.7×, γ_l
30.2×, γ_h 0.082×). The δ discrepancy is almost certainly ledger 1.69× =
median-of-per-patient-ratios vs CSV 0.98× = ratio-of-medians; the two diverge
most for δ because its values are near zero. **Immaterial**: δ cophenet is "no
trace" at p = 0.278 under either ratio, and the figure x-axis uses signed
−log10(p), not the ratio. Flagged per the brutal-honesty rule.

**(c) The empirical-null floor flattens the Grassmann y-axis.** β, γ_l, δ all
sit at `p_mass = 0.005 = 1/(R+1)`, so the cohort p alone cannot rank them; the
y-coordinate is a censored value, not a graded effect size. The LOO-max column
(§1) is the discriminator and is reported alongside.

**(d) Every band is a cohort split; none is universal.** No band has a clean
null in all 10 patients and none has a trace in all 10. The same patients recur
as dissenters (Pat_10/14/15 at β cophenet). The band-selectivity rests on the
cohort gate crossing threshold, which is itself a function of cohort
composition — the matched-strength null validates the cohort median against a
strength-matched null but does not make the effect present in every patient.

**(e) Third-common-cause confound is unfalsifiable** in this four-phase design
for every measure here (a shared task→rest drift that both probes pick up).
This aggregation compares the two probes under a shared null; it does not
certify the underlying trace beyond what the locked §5.3 controls establish.

---

## Provenance

| Quantity | Source |
|---|---|
| ρ^coph cohort gate (all bands) | `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` |
| ρ^coph epi-X (γ_l emergence) | `data/audit/epi_stratified/cophenetic_cohort.csv` (audit_77) |
| Grassmann cluster-mass gate + LOO | `data/audit/grassmann_cluster_extent/cohort_summary.csv` (audit_70) |
| Grassmann epi-X contraction (γ_l, δ) | `VERDICT_LEDGER.md` C5 rows + `data/audit/epi_stratified/README_grassmann.md` (audit_78) |
| Per-patient concordance cross-tab + β breakdown | `.agents/reports/2026-06-01_trace-concordance-vs-blind-fc.md` §2–3 (audit_75, integrity-gated to `beta_per_patient.tex`) |
| Grassmann > cophenet patient-coherence at β | `.agents/reports/2026-05-30_pair-trace-measures-methodology-audit.md` §"Grassmann vs cophenet" |
| Locked verdicts cross-check | `.agents/preprint/locked/VERDICT_LEDGER.md` |

Reproduce the figure:
`python scripts/02_preprint/preprint_34_bands_coph_grassmann_dissociation_map.py`
(surfacing only; reads the three cohort CSVs, recomputes nothing).
