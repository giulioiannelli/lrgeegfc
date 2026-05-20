---
name: 2026-05-07_epileptic-n10-revisit
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-05-07
updated: 2026-05-07
supersedes:
  - .agents/reports/2026-04-18_epileptic-imcoh-final.md
pointers:
  - data/audit/epileptic_n10_revisit/figures/
  - data/audit/epileptic_n10_revisit/M1_raw_imcoh_enrichment.csv
  - data/audit/epileptic_n10_revisit/M2_lrg_enrichment.csv
  - data/audit/epileptic_n10_revisit/M2b_tau_sweep.csv
  - data/audit/epileptic_n10_revisit/M3_mrca_height.csv
  - data/audit/epileptic_n10_revisit/M3_kc_epi_subtree.csv
  - data/audit/epileptic_n10_revisit/M3_kc_epi_subtree_triangle.csv
  - scripts/01_compute/audit/audit_48_epileptic_n10_compute.py
  - scripts/01_compute/audit/audit_48b_epileptic_n10_figures.py
---

# Epileptic-node revisit at n=10 — figure-first report

## Renormalization head (read first)

The n=5 V1 finding (cross-probe |ImCoh| epi-epi/non-non = 1.45× at β,
70/120 BH-FDR) **does not survive the cohort expansion**. With Pat_06,
10, 13, 14 added (Pat_15 has no annotation; n=9 effective), the β
cross-probe enrichment drops to **1.16×** (cohort-Wilcoxon p=0.013, log
ratio one-sided), losing ~half its effect size. The replacement
headline is **δ-band cross-probe enrichment 1.55×** (p<10⁻⁴), which is
broadband-stable and the strongest band signature in the n=9 cohort.
The LRG step at τ_max **passes through** rather than amplifies: ρ̂
ratios match raw |ImCoh| ratios within ±0.05 cross-probe, so heat-
diffusion does not concentrate epi-flow beyond what direct edges
already give. Three new structural findings emerge from the LRG layer
that the raw FC layer cannot see: (i) epi-epi pairs cluster tighter
(lower MRCA merge heights) than non-non in 5/6 bands, (ii) the
epileptic sub-dendrogram shows a **trace** pattern in α-topology (8/9
patients, p=0.006) and γ_l-heights (7/9, p=0.037) — distinct from the
cohort-wide β trace at the global network — and (iii) same-probe |ImCoh|
hyperconnectivity at β / γ_l (2.0× cohort-wide) is real focal
hypersynchrony under ImCoh (Nolte 2004 kills the volume-conduction
confound), and the LRG step dilutes it (Δ ≈ -0.2 cohort-wide), meaning
the focal-edge story does **not** propagate through the heat kernel.

**Verdict in plain words.** At n=10 the strongest single number for
"do epileptic contacts share more network structure than chance?" is
the δ-band cross-probe |ImCoh| ratio (1.55×, all 6 phase × all-pair
cells stable). The epileptic signature in the LRG layer is **not** an
amplification of edge-magnitude — it is a **hierarchical-position**
signature: epi nodes form tight clusters in the LRG dendrogram and the
α / γ_l epi-sub-dendrogram reorganizes during task in a way that
persists into post-task rest. The β finding from n=5 is no longer
defensible at the global cross-probe level; the |ImCoh|² → |ImCoh|
reset and the cohort expansion together account for the loss.

---

## Cohort and eligibility

n=10 cohort with 9 epi-annotated patients (Pat_15 has zero
red-coloured contacts in the implant XLSX). The previous **≥3-probe
eligibility rule** has been **dropped**: it was an MSC-era spatial-
confound mitigation, but under |ImCoh| the volume-conduction
confound is killed by construction (Nolte 2004 phase-lag
suppression — see `probe_bias_critical.md`). Same-probe pairs are
now reported as a *diagnostic* alongside cross-probe, not as an
exclusion criterion. Two patients (Pat_06, Pat_07) have only 2
distinct probes and would have been excluded under the old rule;
they are now included and contribute to the per-patient table.

| Patient | n_epi (channel-matched) | n_probes | new in n=10? |
|:--|:-:|:-:|:--|
| Pat_02 | 14 | 4 (B,G,L,T) | (n=5 baseline) |
| Pat_03 |  6 | 2 (L,O) | (n=5 baseline; was borderline) |
| Pat_05 | 14 | 5 (L',P',X,Y,Y') | (n=5 baseline; cleanest) |
| Pat_06 | 10 | 2 (J,O) | **NEW** |
| Pat_07 |  7 | 2 (A,K) | (n=5 baseline; was borderline) |
| Pat_08 |  9 | 2 (T,W) | (n=5 baseline; was inconclusive) |
| Pat_10 | 10 | 5 | **NEW**, strong implant coverage |
| Pat_13 | 30 | 6 | **NEW**, **largest in cohort** |
| Pat_14 | 12 | 4 | **NEW** (vendor task_test fixed 2026-04-25) |
| Pat_15 |  0 | 0 | dropped — no annotation |

Note: counts above are *channel-matched* (the loader cross-references
red-coloured XLSX labels against the actual recorded channels in
`channel_labels.csv`). The raw XLSX scan I ran in the previous
session reported higher numbers (e.g. Pat_02 25 epi / 6 probes); the
delta is XLSX entries that don't appear in the recorded channel set
(probably reference electrodes or non-recorded contacts). The
matched count is canonical and matches the n=5-era memories.

---

## Figure deliverables

All six PDFs live under
`data/audit/epileptic_n10_revisit/figures/`. Each is fully vector,
no PNG sibling, no rasterised artists, no provenance footer (per
the canonical plotting rules).

### `fig_01_amplification_scatter.pdf`

Per-(patient, phase) scatter, R_imcoh_cp on x and R_ρ_cp at τ_max
on y, one panel per band, dots colored by phase. Diagonal = LRG
passes raw enrichment through. Reads in one glance:

- All six bands cluster tightly along the diagonal — LRG **does
  not amplify or dilute** the raw cross-probe enrichment at τ_max.
- δ panel has the largest spread along the diagonal (R values
  span 0.7–2.2) — broadband epi enrichment with high inter-
  patient variability.
- γ_h panel has all dots clustered near (1, 1) — null finding,
  same conclusion at both layers.

### `fig_02_enrichment_heatmap.pdf`

3 × 3 grid: rows = pair-class (all / cross-probe / same-probe);
cols = raw R_imcoh, LRG R_ρ(τ_max), Δ = R_ρ − R_imcoh. Each
cell is cohort-median over n=9 patients. Stars: BH-FDR q < 0.05
on a per-panel Wilcoxon (one-sided log-ratio > 0).

Headline reads:
- **All-pairs δ-Pre 1.80\*** is the strongest single cell — the
  δ-band epi enrichment is concentrated at rest_pre.
- **Cross-probe δ all phases 1.44–1.83** — δ is the only cross-
  probe-stable enrichment band (p < 10⁻⁴ pooled).
- **Same-probe β-Pre 2.45\*, γ_l-Pre 2.07\*** — focal
  hypersynchrony at rest_pre. Persists through task and post-rest
  but weakens.
- **Δ panel (right)**: cells are mostly blue (negative) at
  same-probe entries (Δ -0.20 to -0.35) — LRG dilutes the focal
  hypersynchrony story but preserves the cross-probe story (Δ ≈
  0 there).

### `fig_03_tau_sweep.pdf`

Cohort R_ρ_cp(τ) median + IQR shading, one panel per band, three
phase lines. x-axis = τ × λ_max (cohort median; per-patient τ grids
were aggregated by τ-index = relative position within each
patient's [τ_min, τ\*] sweep). Reads:

- All bands: enrichment **decays monotonically** as τ increases,
  converging to ~1 by τ × λ_max ≈ 4. No τ regime where the LRG
  step amplifies above the τ_min value. Multi-hop heat diffusion
  does NOT concentrate epi flow.
- δ: highest peak (~1.5–1.7), longest decay tail.
- α: peak ~1.4 at task_test; rest_pre and rest_post slightly lower.
- β: rest_pre rises slightly then decays, peaks ~1.35 at intermediate τ.
- γ_l: weakest peak (~1.1), narrow decay.
- γ_h: flat ≈ 1.0 throughout — null.

This rules out the "multi-hop information shortcut" interpretation:
the heat-diffusion ρ̂ does not preferentially route through epi tissue
at any scale tested, beyond what raw |ImCoh| says at τ = 1/λ_max.

### `fig_04_mrca_violins.pdf`

Cohort-pooled violins of MRCA merge height (= dendrogram height at
which a pair first joins) per pair-class (epi-epi red, non-non
green, epi-non blue), one panel per band, three sub-violins per
phase. Visual answer to "where in the LRG hierarchy do epi pairs
first meet?".

- δ, θ, α, β, γ_l: epi-epi (red) violins are visibly **lower** than
  non-non (green). Cohort-median MRCA ratio epi/non in task_test:
  δ 0.87, θ 0.94, α 0.97, β 0.96, γ_l 0.90, γ_h 0.99.
- Most pronounced at δ and γ_l. β shifts are real but small.
- γ_h: all three classes pile up at MRCA ≈ 1.0 — no discrimination.
- Pattern is consistent across all three phases (rest_pre / task /
  rest_post), so it is a **band signature**, not a phase-dependent
  effect.

Interpretation: epileptic contacts form clusters that merge **earlier**
in the dendrogram than random non-epi pairs. They are tight
sub-networks at low diffusion times.

### `fig_05_kc_epi_triangle.pdf`

Box-and-dot plot of T_KC^epi = d_KC(task, rsPost; epi-epi pairs only)
− d_KC(rsPre, task; epi-epi pairs only). Three panels: λ = 0
(topology only), 0.5 (balanced), 1.0 (heights only). Negative =
the epi-induced sub-dendrogram is closer to its task shape after
task than before — i.e., **trace at the epi sub-network level**.

- **λ = 0 (topology), α**: 8/9 patients in trace direction, Wilcoxon
  one-sided p = 0.006 (only band/λ cell with a star).
- **λ = 1.0 (heights), γ_l**: 7/9 patients in trace direction, p = 0.037.
- λ = 0.5 (balanced): no significant cells.
- β: at most 5/9 trace; no signal at the epi sub-network level —
  the cohort-wide β trace from Result 2 does **not** localize to
  the epi sub-network.

This is a **novel** finding distinct from the cohort-wide Result 2: the
epileptic sub-dendrogram has its own trace signature in α (topology)
and γ_l (heights), in different bands than the global β trace.

### `fig_06_dendrograms_beta.pdf`

9 patients × 3 phases at β. Dendrograms with epi leaves coloured
red (and merges entirely above-epi-leaves coloured red); non-epi
grey. Visual support for fig_04: in most patients the red leaves
form one or two compact subtrees rather than being scattered across
the dendrogram. β dendrograms specifically — chosen because β is the
load-bearing global Section-5 band and it is interesting that β is
*absent* from the epi-sub-network trace (fig 5), so this figure
documents what the β epi-subtree actually looks like (compact, but
not phase-reorganizing).

---

## Cohort-Wilcoxon table (pooled phases, n=27 per cell, log-ratio one-sided > 0)

### Cross-probe (the canonical conservative test)

| band | R_imcoh_median | p | R_rho_median | p | LRG vs raw |
|:--|:-:|:-:|:-:|:-:|:--|
| **δ** | **1.551** | **<10⁻⁴** | **1.527** | **<10⁻⁴** | preserved |
| α | 1.299 | 0.006 | 1.297 | 0.003 | preserved |
| θ | 1.266 | 0.031 | 1.210 | 0.008 | preserved |
| β | 1.162 | 0.013 | 1.177 | 0.013 | preserved |
| γ_l | 1.033 | 0.089 | 1.045 | 0.067 | null |
| γ_h | 0.989 | 0.524 | 0.990 | 0.607 | null |

### Same-probe (focal hypersynchrony diagnostic)

| band | R_imcoh_median | p | R_rho_median | p | LRG vs raw |
|:--|:-:|:-:|:-:|:-:|:--|
| **β** | **2.033** | **<10⁻⁴** | **1.819** | **<10⁻⁴** | -10% (diluted) |
| **γ_l** | **1.969** | **<10⁻⁴** | **1.719** | **<10⁻⁴** | -13% (diluted) |
| δ | 1.498 | <10⁻⁴ | 1.431 | <10⁻⁴ | -5% |
| α | 1.308 | 0.016 | 1.155 | 0.008 | -12% |
| γ_h | 1.079 | 0.022 | 1.073 | 0.033 | -1% |
| θ | 1.095 | 0.214 | 1.037 | 0.105 | null |

### KC epi-subtree triangle (n=9 patients per cell, one-sided p < 0)

| band | λ=0 (topology) | λ=0.5 (balanced) | λ=1 (heights) |
|:--|:--|:--|:--|
| δ | +0.77 (3/9) | +0.03 (4/9) | +0.58 (2/9) |
| θ | +0.42 (4/9) | +0.01 (4/9) | +0.07 (3/9) |
| **α** | **−0.52 (8/9, p=0.006)** | +0.06 (4/9) | −0.18 (7/9, p=0.21) |
| β | −0.06 (5/9) | −0.10 (5/9) | +0.15 (4/9) |
| γ_l | −0.08 (5/9) | −0.02 (5/9) | **−0.18 (7/9, p=0.037)** |
| γ_h | +0.04 (4/9) | −0.12 (5/9) | −0.09 (5/9) |

---

## What changed vs the n=5 era

The previous canonical writeup
(`.agents/reports/2026-04-18_epileptic-imcoh-final.md`, now
superseded) reported:

> "Cross-probe edge ratio: 1.45× mean, 70/120 BH-FDR (58%), 50/120
> Bonferroni. Pooled across 4/5 patients with epi contacts on ≥3
> probes."

The 70/120 BH-FDR figure was 5 patients × 4 phases × 6 bands = 120
cells. With the n=5 cohort restricted to ≥3-probe eligibility (n=4
for V1: Pat_02, 03, 05, 07; Pat_08 inconclusive), the cohort
Wilcoxon at β cross-probe was implicitly significant.

At n=9 (no probe filter, all 9 patients with epi annotations), the
β cross-probe ratio drops to 1.16× cohort-median (vs 1.39× at n=5).
The drop has two components:

1. **Cohort dilution**: Pat_06, Pat_07 (only 2 probes each) and
   Pat_08 (was inconclusive) now contribute their per-patient
   cross-probe ratios, which are below 1.0 at β for several phase
   cells. Pat_13 (largest in cohort, 30 epi) has β cross-probe
   ratios near 1.0 in all phases — its sheer mass pulls the cohort
   median toward 1.

2. **|ImCoh|² → |ImCoh| reset (April 2026)**: the n=5 V1 number
   1.45× was on `imcoh_sq` mislabelled as `imcoh`; recomputing on
   `imcoh_abs` (the canonical metric per Ewald 2012 / Bastos &
   Schoffelen 2016) shifts ratios by roughly √. This was already
   acknowledged in the April reset memo (see
   `imcoh_taxonomy.md`), and the n=5 numbers were updated then to
   1.45× on `imcoh_abs`. The further drop to 1.16× is purely the
   cohort expansion.

The **δ-band cross-probe finding** (1.55× cohort median, p<10⁻⁴) is
new and stronger than the n=5-era β finding ever was. It is also
robust: 9/9 patients have R_cp > 1 at δ in at least 2/3 phases.

---

## Per-patient table (cross-probe, β-band, all 3 phases — for reference)

| Patient | rest_pre | task_test | rest_post |
|:--|:-:|:-:|:-:|
| Pat_02 | 1.40 | 1.16 | 0.91 |
| Pat_03 | 1.39 | 1.41 | 1.36 |
| Pat_05 | 1.31 | 0.86 | 0.71 |
| Pat_06 | 1.27 | 1.30 | 0.99 |
| Pat_07 | 0.83 | 1.18 | 0.71 |
| Pat_08 | 0.93 | 1.18 | 1.04 |
| Pat_10 | 1.06 | 1.15 | 1.55 |
| Pat_13 | 1.16 | 1.11 | 1.20 |
| Pat_14 | 0.85 | 0.78 | 0.83 |

Pat_05 (cleanest in n=5), Pat_07, and Pat_14 contribute the
β cross-probe ratios < 1 (anti-direction) that pull the cohort
median down. Pat_03 is the only patient strongly above 1 in all
3 phases; it is also the 1024 Hz outlier.

---

## What the LRG step adds beyond the raw FC layer (the user's question)

You asked: *"do LRG measures tell us something about epileptic node
connectivity patterns — information flow encoded by edge
heterogeneity rather than topological modular structure"*?

Three findings answer this directly.

**1. ρ̂ at τ_max is approximately a pass-through.**
fig_01 + fig_02 right column: R_ρ ≈ R_imcoh cross-probe
within ±0.05 cohort-median. The heat kernel at the finest scale
does not concentrate or dilute epi-flow. The τ-sweep (fig_03)
confirms this: increasing τ monotonically erodes the enrichment
toward 1.0 — there is no multi-hop regime where epi tissue acts
as a network shortcut.

*Implication*: the user's hypothesis "edge heterogeneity propagates
into preferential information flow at multi-step scales" is **not
supported** at n=9. Epi enrichment is a local-edge phenomenon at the
cross-probe scale.

**2. Hierarchical position is altered.**
fig_04 (MRCA depth violins): epi-epi pairs systematically merge at
LOWER heights than non-non pairs in 5/6 bands. fig_06: visual
confirmation that epi leaves form compact subtrees in most patients.

*Implication*: epi nodes form **tight clusters in the LRG dendrogram**.
This is a multiscale-position signature, not a heat-flow signature.
It is independent of the raw FC enrichment and is visible across all
phases.

**3. The epi sub-network reorganizes in α (topology) and γ_l
(heights), not β.**
fig_05 (KC epi-subtree triangle): T_KC^epi shows trace direction at
α/λ=0 (8/9 patients, p=0.006) and γ_l/λ=1 (7/9, p=0.037). β is null
at the epi-subtree level, despite being the load-bearing global band
in Result 2.

*Implication*: the **epileptic sub-dendrogram has its own task-trace
signature, in different bands than the global cohort trace**. This
is genuinely orthogonal to Result 2 and is invisible at the raw FC
layer (which has no notion of subtree).

---

## What this is NOT (defensive framing)

- **NOT** a stronger finding than Result 2 (β global trace at the
  cohort network level). It is a *different kind of* finding: a
  localized sub-network reorganization within the epi-only induced
  subtree.
- **NOT** clinically actionable at n=9. Effect sizes are modest
  (δ cross-probe 1.55×, α epi-subtree T_KC ≈ -0.5 in normalized KC
  units), and per-patient heterogeneity is high (Pat_15 has zero
  annotation; Pat_06/07 only 2 probes).
- **NOT** a multi-hop information-shortcut story. The τ-sweep
  rules this out cleanly. Whatever edge-heterogeneity-encoded
  information flow the user was asking about, the heat-diffusion
  ρ̂ does not amplify it beyond the raw |ImCoh| enrichment.
- **NOT** independent of Result 2. Both use the same imcoh_abs
  matrices and the same LRG cache. The KC epi-subtree triangle is
  derived from the same dendrograms as Result 2's global KC
  triangle, just restricted to the epi-epi pair subset.
- **NOT** Bonferroni-surviving under the joint corrections of
  Section 5. The α/p=0.006 finding survives BH-FDR within its own
  panel (m=18) but would not survive joint Bonferroni against the
  Section-5 controls (m=42–48). For an epilepsy-specific report it
  is a per-measure-significant finding; for a manuscript claim it
  needs explicit caveat.

**Direction A follow-up (2026-05-08): band-resolved cross-phase
phenomenology.** The induced-subtree analysis at
[`.agents/reports/2026-05-08_direction-a-induced-subtree.md`](2026-05-08_direction-a-induced-subtree.md)
reveals a structured per-band signature on the epi-induced subtree
(scope at
[`.agents/guides/task-persistence-investigation/2026-05-08_epi-cross-phase-rigidity.md`](../guides/task-persistence-investigation/2026-05-08_epi-cross-phase-rigidity.md)).
Three regularities pass the IQR > 20% × |median| floor:

| (band, λ) | direction | %IQR | reading |
|---|---|---|---|
| **δ / λ=1** | **RESET (7/9)** | **75%** | strongest cohort effect — epi heights at slow rhythm shift in task and revert |
| **β / λ=1** | **RESET (6/9)** | **43%** | epi heights reset at the same band where the global cohort network TRACES (Result 2) — band-specific functional dissociation |
| **α / λ=0, 0.5** | **TRACE (5/9, 5/9)** | **55%, 58%** | α-topology of the epi sub-network reorganises and persists |

The β-dissociation is the new theoretical headline: the global β
trace (Result 2, 10/10 patients) is carried by the *non-epi* part
of the network — within the epi sub-network at β/λ=1, the heights
revert. Direct testable consequence: subtracting the epi nodes from
the global tree should sharpen the Section 5 β trace.

The fig_05 α-topology cell (8/9 p=0.006 under pair-mask) attenuates
to 5/9 trace under the induced subtree (still trace direction at
55% of IQR). The γ_l/λ=1 cell replicates identically (7/9 in both
operationalisations, small amplitude). The cross-probe δ-enrichment
headline (1.55×, edge magnitude) is unaffected — it is the
*magnitude* of the coupling and lives at a different layer than
the cross-phase signature.

---

## Bottleneck and downstream

**Bottleneck**: more patients with denser epi annotation. The
cohort effect-size dilution from n=5 (β 1.39×) to n=9 (β 1.16×)
shows that the original V1 was not robust to additional patients.
δ (1.55× at n=9) is the new defensible cross-probe finding but
also needs more patients before any clinical claim can be made.

**Downstream tasks**:

1. (Optional, if the user wants) **anatomical localization** of the
   tight-cluster epi nodes (fig_04 → which Desikan-Killany regions
   do they live in?). Mirrors task #28 from Result 2 but for the
   epi subset specifically.
2. (Optional) **per-band epi MRCA-shift quantification** — turn
   fig_04 into a per-patient scalar (epi MRCA median / non MRCA
   median) and test whether this scalar is correlated with clinical
   variables (seizure-onset zone overlap, post-surgery outcome,
   etc.) once those metadata are available.
3. (Optional) **probe-bias decomposition figure** — split same-probe
   from cross-probe contributions per patient and per band, to show
   the focal-vs-distributed structure of each patient's epileptic
   network.

None are blocking. The current report is self-contained as a
figure-first answer to the cohort revisit.

---

## Reproducibility

```bash
cd /home/giulio/Documents/research/neural_networks/lrgeegfc
conda activate lapbrain
python scripts/01_compute/audit/audit_48_epileptic_n10_compute.py
python scripts/01_compute/audit/audit_48b_epileptic_n10_figures.py
```

Both scripts run in < 5 minutes on the n=9 cohort. The compute step
emits 5 CSVs into `data/audit/epileptic_n10_revisit/`; the figures
step reads only those CSVs and emits 6 PDFs into
`data/audit/epileptic_n10_revisit/figures/`.

---

## TL;DR

The n=5 cross-probe β finding (1.39× → 1.16×) does not survive the
n=10 cohort expansion. The replacement headline at n=9 is **δ-band
cross-probe enrichment 1.55×** (p<10⁻⁴), broadband-stable. The LRG
layer does not amplify this — ρ̂ at τ_max passes through, τ-sweep
shows monotonic decay to 1. What the LRG layer adds is structural:
**epi-epi pairs cluster tighter in the dendrogram (lower MRCA
heights)** in 5/6 bands, and the **epileptic sub-dendrogram shows a
trace pattern in α-topology (8/9 patients, p=0.006) and γ_l-heights
(7/9, p=0.037)** — a band signature distinct from the cohort-wide β
trace at the global network level. The β finding from the n=5 era is
retired; the n=9 picture is a δ-broadband edge-magnitude story plus
a localized α/γ_l hierarchical reorganization story.

---

## Audience verdict — for an epilepsy reader

**The δ-band cross-probe finding is a confirmation, not a discovery.**
Slow / delta hypersynchrony in and around the epileptogenic zone is
already a well-established marker (Bartolomei, Wendling, Spencer et
al., two decades of literature). A 1.5× cross-probe ratio is real but
unremarkable in effect size, and the band (δ) is exactly the band the
literature already flags. We have not added a new phenomenon.

What this report contributes is methodologically tighter rather than
phenomenologically new:

- **Volume-conduction-immune by construction.** |ImCoh| (Nolte 2004
  phase-lag suppression) kills the zero-lag confound that affects
  scalp EEG / MEG studies; we further restrict to *cross-probe* pairs
  in sEEG, eliminating the same-electrode-shaft confound. So the
  same δ phenomenon is shown with two layers of confound control on
  top of what the field usually has.
- **Direct intracranial cohort, n=9.** Most prior literature is scalp
  EEG / MEG; sEEG cohorts are typically 5–10. We confirm the slow-band
  hypersynchrony marker at the upper end of that range with a
  reproducible pipeline.
- **A useful sanity check for the methods paper.** If the LRG framework
  is published as a tool, this finding is the "method recovers known
  epilepsy biology" check that reviewers ask for. It belongs in a
  short section, not a headline.

What this is **not**:

- Not a clinical biomarker. Effect size 1.5×, per-patient
  heterogeneity high.
- Not a multi-hop information-shortcut story. The LRG ρ̂ at τ_max
  passes the raw |ImCoh| ratio through unchanged; the τ-sweep
  monotonically decays to 1. No τ regime where epi-tissue acts as a
  network shortcut.
- Not an LRG-specific result. The dendrogram MRCA tightness (fig_04)
  and the dendrogram-based view (fig_06) are consistent with the same
  δ-coupling story viewed through a clustering lens — they confirm
  but do not add an independent signal.
- Not Bonferroni-surviving against Section-5 joint corrections. The
  α/topology epi-subtree trace (p=0.006) and γ_l/heights (p=0.037)
  survive within-panel BH-FDR but would not survive joint correction
  against the Section-5 controls (m=42–48).

**One-line for an epilepsy talk or paper paragraph**:

> *"We confirm, with intracranial sEEG recordings and a
> volume-conduction-immune connectivity measure, that epileptic
> contacts on different probes are about 1.5× more strongly coupled
> in the delta band than non-epileptic contacts (n=9, p<10⁻⁴). The
> renormalised heat-kernel communicability passes this raw enrichment
> through unchanged; LRG does not reveal an additional multi-hop
> information-shortcut beyond what local edge magnitudes already say."*

---

## Forward look — what would be a real finding?

**The LRG framework has communication-path information that current
epilepsy literature has not exploited. Four directions could yield
something genuinely new (full scope and decision rules in
[`.agents/plans/active/2026-05-08_lrg-epilepsy-research-directions.md`](../plans/active/2026-05-08_lrg-epilepsy-research-directions.md)).**

Ranked by expected novelty × feasibility:

1. **Cross-phase rigidity of the epileptic sub-network**
   ("are epi nodes anchors while the rest of the brain traces?").
   Restrict the Section-5 trace / anchor / reset / emergent taxonomy
   to epi-induced leaf sets. Tests an intuitive but rarely-tested
   hypothesis (epileptic networks are pathologically rigid across
   cognitive states). **Strongest** because it leverages our unique
   3-phase setup (rest_pre / task / rest_post) which is rare in the
   epilepsy literature.

2. **Virtual resection in LRG ρ̂ at the task-emergent scale**.
   Compute global and pair-wise communicability before vs after
   removing the epi set, at τ where the cohort dendrogram structure
   stabilizes. Predicts which patients have epi-as-critical-relay vs
   epi-as-redundant. Connects directly to surgical-outcome literature
   (Kini 2019, Sinha 2017, Jirsa 2017) but uses path-integrated ρ̂
   instead of raw FC.

3. **Eigenmode localization on the epi set** (inverse participation
   ratio of L̂ eigenmodes restricted to epi nodes). Borrows from
   solid-state physics (Anderson localisation). If specific modes are
   sharply localised on the epi set at characteristic eigenvalues,
   that is a fundamentally new descriptor of the epileptic network's
   communication structure. Does not require cross-phase data.

4. **ρ̂-leakage → propagation-zone candidates**. For each non-epi
   node, what fraction of its total communicability flows into the
   epi set? Identifies non-epi nodes most coupled to the epileptic
   network — candidate Propagation Zone (PZ) nodes (Bartolomei
   sEEG-PZ literature). Per-patient figure with anatomical
   localisation. Most clinically interpretable.

The δ cross-probe finding survives any of these as the "known-biology
sanity check"; the four directions are layered on top, not
replacements for it.
