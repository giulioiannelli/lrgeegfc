---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
scope: ranked list of analyses that the manuscript or referee critique flags as missing/deferred
---

# New analyses owed — ranked list

**Head.** The manuscript explicitly defers several analyses to "future work" and another set is implicit-but-owed (referee-defense surface). They divide into **CRITICAL** (the analysis is required for an internally consistent reading of §5.6 / §5.5), **HIGH** (controls that mirror §4 controls but at the LRG layer), and **MEDIUM** (additional probes that strengthen but do not gate the cohort claim). Each item lists what it tests, what data + script it needs, and the corresponding §-reference. Ten items total.

## Ranking

### 1 — CRITICAL — Anchor anatomy baseline (same-probe vs cross-probe fraction at β taxonomy)

**What it tests.** Whether anchor-class modules (the 144 anchor clades at β across the cohort, per §5.6) sit predominantly on same-probe leaf sets (anatomically tight, recording-montage-driven) versus cross-probe leaf sets (anatomically dispersed, neural-driven). If anchors are dominated by same-probe leaf sets, "anchor" is partly a re-detection of the recording-shaft anatomy that §2 already documented for MSC and that §2.2.3 confirms is suppressed under |ImCoh|. If anchors are dominated by cross-probe pairs, anchors are a genuine phase-invariant neural backbone.

**Why CRITICAL.** §5.6 currently treats trace and anchor symmetrically as two of four module-level classes, but the manuscript does not disclose whether the anchor class is anatomy-driven. Without that disclosure, a reader might over-interpret the "144 anchor clades at β" as a rich phase-invariant backbone; if half of them are actually same-probe, the count is inflated.

**Data + script.**
- Inputs: per-cell anchor-class strict-intersection leaf sets from `audit_50_kc_anchor_module_view.py`.
- Per-leaf probe label from `extract_probe_labels` (`src/lrg_eegfc/utils/probe.py`).
- For each anchor clade, compute the fraction of its leaf set on the same probe (modal probe over leaves) vs cross-probe.
- Report: cohort-aggregate same-probe fraction per band, with per-patient breakdowns.

**Section.** Add a one-paragraph control to §5.6 between Tab 4 and the prose count summary. Or move to a §5.6 supplementary control panel.

**Cost.** ≤ 1 hour (data is already cached in audit_50 outputs).

### 2 — CRITICAL — Mutual-exclusivity rule for §5.6 taxonomy

**What it tests.** A leaf can sit in *both* a trace strict-intersection set AND an anchor strict-intersection set (the matching procedure operates per source-clade, no leaf-level uniqueness constraint). The §5.6 prose Tab 4 lists the strict-intersection sets but does not state a mutual-exclusivity rule. **Pat_13 β** in Fig. 34 is the explicit example: trace 2 + reset 4 + rearrange 3/10 + anchor 7 + diffuse 36/119 leaves — these add to 52 + 36 = 88 < 119; the difference (31 leaves) is leaves that participate in multiple class strict-intersections.

**Why CRITICAL.** Without a tie-breaker, a single leaf can be reported in multiple class counts simultaneously, inflating cohort-aggregate counts and making the "560 rearrange clades against 144 anchor, 107 trace, 83 reset" headline ambiguous (are these clade counts or leaf-membership counts? The clade count is fine; the leaf-membership count over classes is over-counted).

**Proposed rule.** A leaf belongs to the **largest** class it satisfies, with priority `trace > anchor > reset > rearrange > diffuse`. Or alternatively `anchor > trace > reset > rearrange > diffuse` (if the goal is to demote trace claims when the leaf is also a phase-invariant anchor).

**Data + script.**
- Inputs: per-cell strict-intersection leaf sets from each of `audit_15_trace_modules.py`, `audit_48_kc_reset_module_view.py`, `audit_49_kc_rearrangement_module_view.py`, `audit_50_kc_anchor_module_view.py`.
- For each (patient, band), apply the priority rule and re-tabulate per-class counts.
- Report: cohort-aggregate counts under the priority rule vs the raw (overlap-allowing) counts.

**Section.** §5.6 prose: add a one-line rule statement before Tab 4. Either keep the raw clade-count headline (560/144/107/83) and clarify in prose that clade overlap is allowed, or apply the priority rule and report adjusted counts.

**Cost.** 2–3 hours (need to thread the priority rule through the 5-class union).

### 3 — HIGH — Multiscale taxonomy (τ-sweep, per-leaf class vector across scales)

**What it tests.** §5.6 fixes `τ' = 1/λ_max` (the finest scale resolved by the propagator). At coarser τ, the dendrogram structure changes and leaves may switch class. A leaf that is "trace" at fine τ may be "rearrange" at coarse τ. This dependence is the multiscale extension of the per-pair claim of §5.3 to the module-level taxonomy of §5.6.

**Section reference.** §6.3 outlook: "Third, the diffusion resolution is fixed at `τ' = 1/λ_max` throughout §5 ... a controlled τ-sweep that reads `D(τ)` at coarser scales tests whether the band-resolved imprint sharpens or dissolves as the propagator integrates over longer pathways". The manuscript explicitly defers this.

**Data + script.**
- Inputs: LRG NPZ caches from `data/cache/lrg_imcoh/Pat_*/<band>_<phase>_lrg_imcoh-abs.npz` for τ-sweep at e.g. `τ ∈ {1/λ_max, 0.5/λ_max, 0.1/λ_max, 1/λ_gap, τ*}`.
- Generator: new script `scripts/01_compute/audit/audit_57_tau_sweep_taxonomy.py` (slot 57 unused; 56 = trace_minus_epi).
- For each τ, recompute the §5.6 taxonomy class assignment per leaf.
- Report: cohort-aggregate class transition matrix between fine and coarse τ for each band.

**Cost.** 1–2 days (the LRG dendrogram at coarse τ is fast; the taxonomy per-cell pass is the bottleneck).

**Connection to scope reports.** Should land as a scope report under `.agents/guides/task-persistence-investigation/` first (e.g. `2026-05-09_module-taxonomy-tau-sweep.md`) before the audit script.

### 4 — HIGH — Independent matched-strength surrogate Laplacian null

**What it tests.** §5.3 uses a **split-half within-baseline null** (`ρ_drift`): rsPre is split into two halves and the LRG is built on each half independently. This null shares baseline correlation structure with the cohort target. An **independent surrogate null** would draw a Laplacian with the same edge-strength distribution but randomized topology; this tests whether the trace direction is specific to the observed network structure or recoverable from any matched-strength random network.

**Section reference.** §5.3 mentions only the split-half null; the independent surrogate is not discussed. §6.3 mentions "two cohort-wide drift insulations remain owed" — this is one of them.

**Data + script.**
- Inputs: per-(patient, band, phase) edge-weight distribution; resample to obtain a permuted upper-triangular vector; re-form the symmetric matrix.
- Generator: new script `scripts/01_compute/audit/audit_58_matched_strength_null.py`.
- For R = 100 surrogate Laplacians, compute the §5.3 `ρ_split` distribution under the null.
- Compare observed `ρ_split` against this null at each (patient, band) cell.

**Section.** Add to §5.3 control battery (between drift-floor and cross-probe restriction).

**Cost.** 1–2 days (LRG eigendecomposition is the bottleneck; cache the surrogate spectrum once).

### 5 — HIGH — LRG drift-triangle null (cohort-wide)

**What it tests.** §4.6 Control 1 fits a linear within-rest drift over time gap and reports cohort-median `R²`. This is a per-band scalar diagnostic. The LRG-layer mirror would compute the analogous drift-triangle on `D(τ)` (per-pair distances within rsPre as a function of half-segment time gap) and report the same `R² < 0.05` floor. Currently this is run only at Pat_06 (per the §6.3 outlook reference).

**Section reference.** §6.3: "first, a cohort-wide drift-triangle null at the LRG layer, mirroring Control 1 of §4.6, currently run only at Pat_06 where `R² ≤ 0.05` in nearly every (band, distance) cell".

**Data + script.**
- Inputs: per-patient rsPre split into 4 chunks; build LRG per chunk; compute `D(τ)` per chunk.
- For each pair of chunks, compute distance vs time gap; fit linear drift; report `R²`.
- Generator: extend `scripts/01_compute/audit/audit_28_drift_triangle_null.py` to LRG layer (currently substrate-only) — or new script `scripts/01_compute/audit/audit_59_lrg_drift_triangle.py`.

**Cost.** 2 days (per-chunk LRG is fast; per-pair distance is the bottleneck; cohort × band × half-pairs = 10 × 6 × 6 = 360 cells).

### 6 — HIGH — Symmetric cross-baseline check at LRG layer

**What it tests.** §4.6 Control 3 compares `d_S(taskTest.late, rsPost.early)` against `d_S(rsPre.late, rsPost.early)` on matched half-segment noise budgets. The LRG-layer mirror computes the analogous comparison on `D(τ)` per-pair distances at the same matched half-segments.

**Section reference.** §6.3: "the second is a symmetric cross-baseline check at the diffusion layer, mirroring Control 3 of §4.6, comparing `D(τ)` computed on (taskTest-late, rsPost-early) against (rsPre-late, rsPost-early) on matched half-segment noise budgets".

**Data + script.**
- Inputs: per-patient phase-half segments; build LRG per half; compute `D(τ)` per half.
- Compare task-late→post-early distance against rsPre-late→post-early distance.
- Generator: extend `scripts/01_compute/audit/audit_31_cross_baseline.py` to LRG layer — or new script `scripts/01_compute/audit/audit_60_lrg_cross_baseline.py`.

**Cost.** 2 days (similar profile to drift-triangle).

### 7 — HIGH — Same-probe vs cross-probe disclosure for §4 within-baseline null

**What it tests.** §4 within-rsPre split-half null (`d_S` between two halves of rsPre) is computed over **all** edge pairs. The cross-probe restriction is only applied at the LRG layer in §5.3. A reader cannot tell whether the §4 substrate trace is sharpened or attenuated by cross-probe restriction (even though §2.2.2 says |ImCoh| same-probe enrichment is below 1.5×).

**Section reference.** §4.6 caveats and the §5.3 cross-probe restriction together imply a parallel disclosure at substrate level is reasonable.

**Data + script.**
- Inputs: per-(patient, band) edge-weight vectors with cross-probe mask.
- Compute §4 distances (`d_S`, `d_P`, `d_F`) on cross-probe-only edges.
- Compare `T_d` cohort-median values to full-cohort values.
- Generator: extend `scripts/01_compute/audit/audit_25_raw_fc_phase_distance.py` with a `--cross-probe-only` flag.

**Cost.** ≤ 1 day (substrate-level analysis is fast).

### 8 — MEDIUM — Band-agnostic LRG read (broadband |ImCoh|)

**What it tests.** Whether the band-resolved trace at α / β / low-γ also surfaces in a broadband |ImCoh| matrix (averaged over the full 0.5–300 Hz frequency range). If yes, the band-resolved view is sharpening a pre-existing broadband signal. If no, the trace lives genuinely in band-specific dynamics.

**Section reference.** Not currently in the manuscript; would be a §5 supplementary control or a §6.3 forward-look entry.

**Data + script.**
- Inputs: per-(patient, phase) freq-resolved ImCoh cache; band-average over the full 0.5–300 Hz range.
- Build LRG; compute `D(τ)`; run §5.3 `ρ_split` test.
- Generator: new script `scripts/01_compute/audit/audit_61_broadband_lrg.py`.
- Scope report first under `.agents/guides/task-persistence-investigation/2026-05-09_broadband-lrg-control.md`.

**Cost.** 1–2 days (one extra LRG build per (patient, phase) cell = 40 cells).

### 9 — MEDIUM — Eigenmode embedding view (per-contact motion across phases)

**What it tests.** Per-contact motion in the leading-eigenmode subspace across phases. For a given leading-`k` eigenmode subspace `U_k^Φ`, project each contact `i` to `(U_k^Φ)^T e_i`. Compare the per-contact projection across phases to identify which contacts move most in the eigenmode embedding from rsPre → taskTest → rsPost. A trace-direction contact would have its rsPost projection sit closer to taskTest than to rsPre.

**Section reference.** §5.4 Grassmann measures the chordal angle between subspaces (cohort-aggregate); the per-contact motion within the subspace is not measured.

**Data + script.**
- Inputs: per-(patient, band, phase) leading-`k` eigenvectors from `IMCOH_LRG_CACHE` NPZs (with the eigenvector cache extension from `2026-04-29_eigenvector-direct-pivot-plan.md`).
- For each contact, compute per-phase projection vector; compute pairwise distance.
- Generator: new script `scripts/01_compute/audit/audit_62_eigenmode_embedding.py`.
- Scope report first under `.agents/guides/task-persistence-investigation/2026-05-09_eigenmode-embedding.md`.

**Cost.** 2 days (eigenvector cache extension is the gating factor; once that lands, the analysis is straightforward).

### 10 — MEDIUM — Epi-contact fraction in §5.5 region pools

**What it tests.** For each of the §5.5 anatomy regions (γ_l ctx-lh-fusiform Bonferroni-survived, β Hippocampus / fusiform / superior-temporal uncorrected), compute the fraction of the cohort contact pool labeled "epileptic" (per `load_epileptic_nodes` from `utils.io.patient`). If the fusiform Bonferroni-survived cell has > 30% epi contacts, the trace-leaf signal at γ_l fusiform may be partly driven by epileptic dynamics.

**Why MEDIUM, not CRITICAL.** The §5.5 anatomy is independent of the cohort §5.3 claim (per §5.5 prose: "by construction, the localization layer and the controlled cohort claim are mutually consistent rather than independent"). An epi-contact contribution to the γ_l fusiform cell would not invalidate §5.3 but would contextualize §5.5.

**Section reference.** §5.5 anatomy paragraph; would add a one-paragraph epi-fraction disclosure.

**Data + script.**
- Inputs: per-region cohort contact pool from `audit_50_anatomy_55_figure.py`; per-contact epi labels from `load_epileptic_nodes()`.
- Compute epi fraction per region in {γ_l fusiform, β Hippocampus, β fusiform, β superior-temporal, α parsopercularis}.
- Generator: extend `scripts/01_compute/audit/audit_50_anatomy_55_figure.py` with an epi-fraction column — or a new helper script.

**IMPORTANT — NOT a hypergeometric epi × β-trace test.** Per `feedback_epilepsy_not_trace_locked.md` and `16_epilepsy_pointer.md`, the epi track is parallel and not folded into the manuscript via TARR. This analysis is a **descriptive disclosure** ("X% of contacts in region R are epi") not a hypothesis test ("epi contacts trace more").

**Cost.** 0.5 days (data is cached).

## Summary table

| # | Priority | Owed deliverable | Section |
|---:|---|---|---|
| 1 | CRITICAL | Anchor anatomy baseline (same-probe vs cross-probe) | §5.6 |
| 2 | CRITICAL | Mutual-exclusivity rule for §5.6 taxonomy | §5.6 |
| 3 | HIGH | Multiscale taxonomy (τ-sweep, per-leaf class) | §6.3 |
| 4 | HIGH | Independent matched-strength surrogate Laplacian null | §5.3 |
| 5 | HIGH | LRG drift-triangle null (cohort-wide) | §6.3 |
| 6 | HIGH | Symmetric cross-baseline check at LRG layer | §6.3 |
| 7 | HIGH | Same-probe vs cross-probe disclosure for §4 within-null | §4.6 |
| 8 | MEDIUM | Band-agnostic LRG read (broadband \|ImCoh\|) | §5 supp |
| 9 | MEDIUM | Eigenmode embedding view (per-contact motion) | §5.4 |
| 10 | MEDIUM | Epi-contact fraction in §5.5 region pools | §5.5 |

## Action — next session

Items 1 and 2 (CRITICAL) should be folded into the §5.6 verification cycle in the same revision pass. Items 3–7 (HIGH) should each spawn a scope report under `.agents/guides/task-persistence-investigation/` before any code lands. Items 8–10 (MEDIUM) can be batched as a single supplementary control supplement once the CRITICAL + HIGH list is closed.

**Order of operations for the next revision cycle:**

1. CRITICAL items 1 + 2 → §5.6 prose update + cohort table re-cite. ≤ 1 day.
2. HIGH item 7 (cross-probe disclosure for §4) → §4.6 caveats expansion. ≤ 1 day.
3. HIGH items 5 + 6 (drift-triangle + cross-baseline at LRG) → §5.3 control battery expansion. 4 days.
4. HIGH item 4 (matched-strength null) → §5.3 control battery expansion. 2 days.
5. HIGH item 3 (multiscale taxonomy τ-sweep) → §6.3 forward-look preview + scope report. 2 days.
6. MEDIUM items 8–10 → batch as one supplementary section.

Total estimated time-to-completion: ~12 working days for the closure of the CRITICAL + HIGH list.
