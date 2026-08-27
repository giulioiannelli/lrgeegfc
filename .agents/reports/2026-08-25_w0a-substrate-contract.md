---
name: 2026-08-25_w0a-substrate-contract
kind: report
era: PAPER_FINALIZATION (Wave 0, lane W0-A)
status: in-progress
created: 2026-08-25
scope: Evidence report for the W0-A substrate contract. Decides the two choices that underpin every result in the paper -- the |ImCoh| magnitude transform and the sparsified backbone -- under a decision rule pre-registered before any number was computed. Reports the full sparsification stability surface, the parameter-free-filter verdict, and an explicit list of what remains unsettled.
pointers:
  - scripts/01_compute/paper_final/w0a_00_preregistration.md
  - .agents/preprint/locked/PIPELINE_CONTRACT.md
  - src/lrg_eegfc/workflow/substrate.py
  - data/paper_final/w0a_substrate/
  - .agents/plans/active/2026-08-25_paper-finalization-master-plan.md
---

# W0-A — the substrate contract

## Head

PLACEHOLDER-HEAD

---

## 0. What this lane was asked to settle, and the firewall around it

Two choices underpin every multiscale result in the paper and neither was decided by evidence. `imcoh_abs = ⟨|Im C|⟩_f` was chosen on 2026-04-15 by argument (a genuine Jensen order-of-operations fix plus a citation to Ewald 2012 / Bastos & Schoffelen 2016), with no head-to-head comparison against the equally citable `imcoh_sq = ⟨(Im C)²⟩_f`. The backbone fraction `f = 0.20` was chosen because, in the methods guide's own words, "only 0.20 gives α strong + β maximal + θ/low_γ null" — selecting the knob so the answer comes out right, and failing that same guide's own criterion 3 ("parameter-free"), which lists `mst@f` in its method table as "reject (free parameter f)".

The decision rule for both was written and frozen in `scripts/01_compute/paper_final/w0a_00_preregistration.md` **before any number was computed**, with the mandatory 5-point critical preamble. It selects on hypothesis-independent properties only: estimator reliability, distance from a lag-destroying null floor, weight heterogeneity, and invariance of the verdict under the knob. Trace gates and marker AUCs under each option are computed and recorded but are **not admissible** at any step of the rule. The rule was not relaxed after the numbers arrived; one dated addendum (R1-bis) was added after A1 and before any A2 trace number, and it is flagged as such in place.

Everything here is conditioned on **matched-strength** being the null — a null injected at the finished-FC-matrix stage, which therefore cannot test the coherency estimator, the band split, session nonstationarity, or drift. That limitation is lane W0-B's to remove and is not removed here.

---

## 1. A1 — the transform, head to head

### 1.1 What was computed

Every number in this section is recomputed fresh from the raw timeseries, not sourced from any prior report or CSV. For each of the 10 patients: one Welch pass per `rest_pre` half (at `nperseg_for_fs(fs) // 2`, the incumbent split-half setting) and one at full duration and full `nperseg`, giving the complex band coherency `C(f)` per band. From the same object both transforms are derived per frequency bin and then band-averaged — the order the canonical loader uses.

**Reproduction check.** The `abs` arm reproduces the existing `imcoh_halves_fc` cache to `max|diff| ≤ 3.0e-8` for every patient (float32 round-trip of the incumbent files), so this is the same estimator, not a re-implementation that drifted.

### 1.2 Step 1 — the rank-equivalence gate: the transform is NOT a no-op, but its effect is bounded

The brief's prior was that squaring is a monotone map on non-negative weights, so any rank-based backbone should be invariant to it. **That is false here, and the reason is instructive**: the transform is applied per frequency bin and the band average is taken afterwards, so `mean(|x|)` and `mean(x²)` are not related by a monotone map — Jensen's inequality again, the same mechanism as the 2026-04-15 reset.

| band | Spearman(abs, sq), median | mst-union edge-set Jaccard @ f=0.05 | @ f=0.10 | @ f=0.20 |
|---|---|---|---|---|
| delta | 0.987 | 0.823 | 0.847 | 0.862 |
| theta | 0.987 | 0.782 | 0.823 | 0.846 |
| alpha | 0.986 | 0.784 | 0.808 | 0.834 |
| beta | 0.990 | 0.828 | 0.854 | 0.882 |
| low_gamma | 0.987 | 0.896 | 0.903 | 0.881 |
| high_gamma | 0.990 | 0.813 | 0.867 | 0.891 |

Across all 60 (patient, band) cells the edge-rank Spearman spans `[0.882, 0.998]` and the backbone edge-set Jaccard at `f = 0.20` spans `[0.510, 0.987]`. **So the transform moves between 1 % and 49 % of the retained backbone edges depending on the cell, with a cohort median of 12–17 %.** The gate does not fire; the choice is not vacuous. But the effect is bounded, and that bound — a rank correlation never below 0.88 and a typical backbone overlap of ~0.86 — is itself the cleanest thing to quote about the transform question, in either direction.

### 1.3 Step 2 — split-half reliability: the rule selects `imcoh_sq`, by a small margin

Spearman is the primary statistic because it is invariant under any strictly monotone map, so a difference can only come from the band-averaging order and not from units. Pearson is reported as secondary and was not used to decide.

| band | Spearman rel. abs | sq | Δ(abs−sq) | winner | [Pearson abs / sq] |
|---|---|---|---|---|---|
| delta | 0.5950 | 0.5971 | −0.0021 | sq | 0.820 / 0.875 |
| theta | 0.6585 | 0.6652 | −0.0067 | sq | 0.795 / 0.821 |
| alpha | 0.6248 | 0.6399 | −0.0151 | sq | 0.776 / 0.816 |
| beta | 0.7313 | 0.7402 | −0.0089 | sq | 0.845 / 0.889 |
| low_gamma | 0.5481 | 0.5587 | −0.0106 | sq | 0.848 / 0.864 |
| high_gamma | 0.6297 | 0.5852 | **+0.0445** | **abs** | 0.701 / 0.634 |

Cohort over the 60 cells: median Δ = −0.0078, Wilcoxon two-sided p = 5.0e-6, band wins abs 1 / sq 5. The pre-registered decisiveness condition (≥ 4/6 band wins **and** p < 0.05) is met, so **R1 stops at step 2 and selects `imcoh_sq`**.

Two honest qualifications, neither of which changes the rule's output. First, the margin is ~1.2 % relative on reliabilities of 0.55–0.74 — systematic, but small. Second, `high_gamma` is the one band where `abs` wins, and it wins by five times the median margin.

**Post-hoc diagnostic (not part of the rule).** The strongest alternative reading of a reliability win is that it is inherited from coherence *magnitude* — anatomy and volume conduction, which are highly reproducible across halves and have nothing to do with lag estimation. It is not: `Spearman(observed edges, band-averaged |C|)` differs between the two transforms by ≤ 0.005 in five of six bands, and the one appreciable difference (high_gamma, +0.055) is in the band where `abs` wins reliability. The diagnostic can only weaken the step-2 conclusion and does not.

### 1.4 Step 3 — the circular-shift floor, and a polarity warning that matters for lane W0-B

A circular time shift of channel `c` by `t_c` maps `C_ij(f) → C_ij(f)·e^{−2πif(t_i−t_j)}`: it preserves `|C_ij(f)|` — hence every bit of volume conduction and of the coupling magnitude — exactly, and randomises only the phase, i.e. the lag. Thirty-two draws per (patient, band), frequency-domain phase-ramp form.

**The observed `|ImCoh|` sits far BELOW its lag-randomised floor, in 120 of 120 cells.** Median `obs/floor` is 0.13–0.32 for `abs` and 0.02–0.11 for `sq`; median per-edge `z` runs from −3.3 (delta) to −96 (high_gamma); the fraction of edges above their own floor is 0.3–5 %. This is expected from the mechanism, not a pathology: under uniform random phase `⟨|Im C|⟩ → (2/π)·⟨|C|⟩`, and the true phase distribution is concentrated near zero lag (which is what volume conduction does), so destroying the lag structure *inflates* `|ImCoh|` rather than deflating it.

The consequence is a design warning, and it is the most transferable thing in this section: **"distance above the circular-shift floor" is the wrong polarity for a per-edge comparison and must never be reported as a null-clearing statistic.** The floor comparison has to be framed as a *structure* test — does the surrogate reproduce the FC *pattern*? — rather than a magnitude test. Framed that way it is informative and reassuring: `Spearman(observed edges, floor edges)` has a cohort median of −0.13 to +0.09 in five bands (high_gamma −0.26 / −0.31), so the observed edge ordering is essentially **uncorrelated** with the magnitude-driven floor ordering. Neither transform is separated by this, so step 3 does not arbitrate — and it would not have been reached anyway, since step 2 was decisive.

### 1.5 Step 4 — weight heterogeneity (recorded; step 2 had already decided)

| band | Gini abs / sq | CV abs / sq | participation ratio abs / sq |
|---|---|---|---|
| delta | 0.345 / 0.606 | 0.70 / 1.47 | 0.670 / 0.317 |
| theta | 0.292 / 0.516 | 0.55 / 1.24 | 0.767 / 0.394 |
| alpha | 0.301 / 0.531 | 0.60 / 1.29 | 0.734 / 0.374 |
| beta | 0.271 / 0.514 | 0.57 / 1.49 | 0.756 / 0.310 |
| low_gamma | 0.216 / 0.457 | 0.53 / 1.47 | 0.783 / 0.316 |
| high_gamma | 0.122 / 0.255 | 0.26 / 0.67 | 0.937 / 0.692 |

`imcoh_sq` roughly **doubles** the edge-weight Gini in every band and drops the participation ratio from ~0.75 to ~0.32 — i.e. the effective fraction of edges carrying the weight falls by more than half. On a dense graph, weight heterogeneity is precisely what drives the propagator into the Villegas single-peak regime, so this is a real cost of the transform the rule selected. The rule is lexicographic and step 4 is a final tie-break, so it does **not** override step 2 — but it is why R1-bis exists.

### 1.6 R1-bis — the dated addendum

Written after A1 and **before any A2 trace number**: a 1.2 % reliability margin is not a sufficient basis on which to regenerate every artifact in the project, when the same switch moves 10–22 % of the backbone and doubles the weight heterogeneity. So A2 is run under **both** transforms. If the per-band, per-scale verdict is invariant across the transform inside the plateau, the contract records the transform as **immaterial within the plateau** and retains `imcoh_abs` for continuity, citing the demonstrated invariance — not the reliability comparison — as the justification. If the verdict is not transform-invariant, `imcoh_sq` is adopted exactly as R1 dictates and every downstream artifact is regenerated on it.

---

## 2. A2 — the sparsification stability surface

PLACEHOLDER-A2

---

## 3. A3 — killing the free parameter

PLACEHOLDER-A3

---

## 4. A4 — one substrate or two?

PLACEHOLDER-A4

---

## 5. The library deliverable

`src/lrg_eegfc/workflow/substrate.py` holds the substrate in one place so the three science lanes cannot silently diverge. `CANONICAL` is a frozen `Substrate` dataclass (transform, backbone, frac, disparity_alpha, plateau, plateau_fracs); `canonical_graph(patient, phase, band)` returns the analysed adjacency; `canonical_eig`, `canonical_phase_graphs`, `canonical_phase_eigs`, `canonical_structure` and `canonical_scale_grid` cover the rest. Per-call overrides exist for robustness sweeps but are explicit arguments, and `dense=True` is the visible escape hatch for a raw-FC baseline rather than a separate code path.

**The phase set is data, never a hardcoded tuple.** `CANONICAL_PHASES` is the five-phase set with `task_learn` first-class, and `canonical_phase_eigs` accepts any phase set and feeds `cross_phase_functionals_over_scales`, which reads the semantic roles (`baseline_a`, `baseline_b`, `encode`, `probe`, `follow`) off whichever phases are present. Dropping `task_learn` returns the standard trace alone through the *same* call, so a four-phase caller and a five-phase caller share one code path.

`canonical_graph_ensemble(patient, phase, band)` yields `(frac, graph)` across the declared invariance plateau — the knob-integrated readout, for use whenever the backbone is not parameter-free, so no single fraction is ever reported as the pipeline setting.

Acceptance-checked on real FC (`scripts/01_compute/paper_final/w0a_verify_canonical_graph.py`), all passing: `canonical_graph` reproduces the incumbent inline pipeline **bit-for-bit** (`max|diff| = 0`); the five-phase path yields all four functionals; the four-phase path yields exactly `{T_probe}`, equal to the five-phase `T_probe` to 0; split-half phases resolve under both spellings; overrides and the plateau ensemble behave.

Supporting library work, all with general names and no manuscript-local tokens:

- `utils/fc/heat_multiscale.py` — `CROSS_PHASE_ROLES`, `cross_phase_functionals`, `cross_phase_functionals_over_scales`. Verified on real FC to reproduce `rho_sym_over_scales` **and** the `05_enc_inf_arc.py` reference implementation of all four functionals to 1e-16. All Spearmans now come from one rank correlation matrix per scale instead of ~10 pairwise `spearmanr` calls, which is what makes carrying the functional axis nearly free.
- `utils/fc/split_half.py` — `imcoh_split_half_adjacencies`, `band_transform_signed`. Generalises the script-local `compute_imcoh_abs_halves` over the transform, so the abs/sq comparison runs through one code path.
- `utils/fc/backbone.py` — `top_fraction_threshold` (the plain global threshold: the mechanism contrast for density sweeps), `backbone_structure` (structural covariates), and the `"thresh"` branch of `select_backbone`.
- `utils/metrics/graph_descriptors.py` — `gini_coefficient`, `weight_heterogeneity`.
- `config/paths.py` — `IMCOH_HALVES_CACHE`, the canonical home for split-half FC under both transforms (the legacy `imcoh_halves_fc` holds the `abs` arm only).

---

## 6. What is NOT settled

PLACEHOLDER-OPEN

---

## 7. Compute provenance

PLACEHOLDER-PROVENANCE
