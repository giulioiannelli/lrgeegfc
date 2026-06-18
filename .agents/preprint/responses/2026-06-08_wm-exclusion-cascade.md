---
name: wm-exclusion-cascade
era: IMCOH_ABS_COHORT_N10
status: current
kind: response
scope: cascade summary for the C6 white-matter-exclusion sensitivity layer — what changed in the locked artifacts + briefs, and where the numbers live
---

# C6 white-matter-exclusion — cascade summary (2026-06-08)

**Head.** A new sensitivity layer — **C6, white-matter exclusion** — was added to
the control battery and cascaded into the locked ledgers and four band briefs. It
is the **second sensitivity layer** alongside C5 (epi-zone exclusion), **secondary
and mechanistic, not a verdict gate**. It **confirms every trace verdict and
changes none**: removing white matter from the montage does not weaken the trace;
the trace **survives** a gray-only montage and is not a volume-conducted
white-matter recording artifact.

> **⚠ 2026-06-12 amendment (two honesty flags resolved — read this).** The
> original 2026-06-08 framing below said WM removal *sharpens* the trace
> (gray-resident). The **`audit_85` random-node-decimation control** shows the
> cophenetic "sharpening/emergence" is **NOT WM-specific** — removing any ~40 % of
> nodes does the same (cophenet cohort `p_dec`: β 0.145, α 0.305, γ_l 0.185, all
> `generic_nodecount`). So the "sharpening" claim is **withdrawn for the
> cophenetic substrate**; the load-bearing C6 claim is **survival/robustness**
> (the trace still clears matched-strength on the gray-only graph; WM is neither
> carrying nor diluting it). The **raw |ImCoh|** substrate IS genuinely WM-specific
> (β/γ_l `p_dec` 0.000/0.030). Separately, **`audit_86`** re-runs the locked C3
> cluster-extent gate: **β/γ_l Grassmann RE-PASS** (`cluster_p_mass^wmX` 0.005,
> strong, LOO-robust), **δ FAILS** (0.105). Read the per-band table below with
> these corrections; the "Resolved" section at the end has the details. No verdict
> flips.

Canonical analysis report (numbers, methods, limitations):
`.agents/reports/2026-06-08_white-matter-exclusion.md`. Audits: `audit_83`
(WM-stratified cophenetic + raw), `audit_84` (WM-stratified Grassmann). Caches:
`data/audit/wm_stratified/`. Figures: `data/audit/wm_stratified/figures/`.

## What C6 is

Recompute each trace probe with **dominant-white-matter contacts removed**
(`load_channel_regions(pat)["region"] == "Wm"`, atlas argmax). WM is 30–57 % of
every montage (median 117→74 nodes under `exclude_wm`). Because ImCoh is pairwise,
the gray-only FC equals the gray-only submatrix of the cached FC to numerical
precision, so `exclude_wm` is the exact "drop WM channels before computing FC"
test; only the node-coupled LRG step changes. The matched-strength null (R=200) is
regenerated on each submatrix. Full definition: `locked/CONTROLS.md` §C6.

## What changed (no verdict flips)

| File | Change |
|---|---|
| `locked/CONTROLS.md` | New §C6 block; head updated to "two sensitivity layers"; "## The sensitivity layers" |
| `locked/VERDICT_LEDGER.md` | Dated 2026-06-08 revision-history entry (no verdict / verdict-layer changes) |
| `bands/01_beta.md` | frontmatter `verdict_layers.{rho_split_coph_wm_excluded, grassmann_wm_excluded}` + revision_history line |
| `bands/02_alpha.md` | same two frontmatter entries + revision_history line |
| `bands/03_gammalow.md` | same + revision_history line |
| `bands/06_delta.md` | same + revision_history line |

No `bands/00_cohort.md`, anatomy ledger, or prose/table edits in this cascade —
the structured `verdict_layers` entries + this summary are the authoritative
record; the writing agent can lift prose/table rows from them via a directive if
desired.

## The result, per band

Reading corrected per the 2026-06-12 amendment (decimation `p_dec` from audit_85,
Grassmann gate `cluster_p_mass^wmX` from audit_86):

| band | cophenet under exclude_wm | Grassmann on the C3 gate (audit_86) | reading (corrected 2026-06-12) |
|---|---|---|---|
| α | survives (obs +0.105→+0.183, p=0.014, LO-P15 0.027) | no_trace→no_trace (`cp_mass^wmX`=0.488; α cophenet-only) | survives; strengthening NOT WM-specific (`p_dec`=0.305) |
| β | survives (obs +0.221→+0.324; p=0.053 marginal, LO-P15 0.037, gray_gray 0.032, raw-β 0.024) | **RE-PASSES gate** (`cp_mass^wmX`=0.005, strong, LOO 0.005) | survives; strengthening NOT WM-specific (`p_dec`=0.145); Grassmann re-passes locked gate |
| γ_l | cophenet "emerges" (p=0.019) — but NOT WM-specific | **RE-PASSES gate** (`cp_mass^wmX`=0.005, strong, LOO 0.005) | cophenet emergence is node-count (`p_dec`=0.185), NOT a WM-unmasked trace; **raw** γ_l IS WM-specific (`p_dec`=0.030); Grassmann re-passes gate |
| θ | "emerges" (p=0.005) | absent (`cp_mass^wmX`=0.413) | node-count effect (`p_dec`=0.225), not WM |
| δ | absent (p=0.500); raw-δ "emerges" p=0.019 but NOT WM-specific (`p_dec`=0.355) | **FAILS gate** (`cp_mass^wmX`=0.105, no_trace) | δ Grassmann WM-dependent (confirmed on gate) — the cross-probe epi-biology channel (Decision 5), not the task trace |
| γ_h | absent | emerges on gate (`cp_mass^wmX`=0.005) but LOO-fragile (0.050) | no task trace at full graph; gate emergence LOO-fragile |

The raw |ImCoh| substrate (historically weak: audit_74 best p≈0.057) clears
matched-strength at δ/α/β/γ_l once WM is removed — and for **β/γ_l this IS
WM-specific** (decimation `p_dec` 0.000/0.030): its full-graph weakness genuinely
was WM strength-structure the null was matching. (For δ/α the raw clearing is
within the generic-node-count band, `p_dec` 0.355/0.300.)

## Resolved (2026-06-12 — the two honesty flags now have compute behind them)

1. **Node-count confound → RESOLVED by `audit_85` (random-node decimation).** For
   each (patient, band) we drop K = #WM nodes *at random* R=200 times and rebuild
   ρ_split, then ask where the observed `exclude_wm` ρ_split sits in that
   distribution (cohort-level `p_dec`). **Verdict (`decimation_control_cohort.csv`):**
   - **Cophenetic: NOT WM-specific** — no band crosses `p_dec < 0.05`
     (β 0.145, α 0.305, γ_l 0.185, θ 0.225, all `generic_nodecount`). Removing any
     ~40 % of nodes raises the cophenet cohort ρ_split the same amount. **The
     "WM removal sharpens the trace / gray-resident-because-it-sharpens" framing is
     withdrawn for the cophenetic substrate.** WM is neither carrying nor diluting
     it; the trace simply **survives** WM removal (still clears matched-strength on
     the gray-only graph — that test is untouched).
   - **Raw |ImCoh|: genuinely WM-specific** — β `p_dec`=0.000, γ_l `p_dec`=0.030
     (`WM_specific`). The raw substrate's matched-strength failure really was
     white-matter strength-structure; gray-only raw FC is a real tissue improvement.
2. **C6-Grassmann gate → RESOLVED by `audit_86` (locked C3 cluster-extent on the
   WM submatrix).** Re-running the audit_70 cluster-extent **mass** gate on the
   `exclude_wm` per-`k` T_G curves (`grassmann_cluster_extent_wmX.csv`):
   **β and γ_l RE-PASS** (`cluster_p_mass^wmX` = 0.005 both, strong, LOO-robust) —
   the per-`k` "40→40 / 41→40" read is now backed by the actual locked gate; **δ
   FAILS** (0.105, no_trace — confirming the 23→7 per-`k` weakening on the gate;
   WM-dependent epi channel); γ_h emerges (0.005) but LOO-fragile (0.050); α/θ
   absent both.
3. **Single WM definition** (atlas-dominant `Wm`; PTD<0 still not run — open).
4. **Not elevated to a primary interpretive lens.** C6 stays secondary; the
   resolution does not change that (and reinforces it — the cophenetic
   "strengthening" is not a tissue result).

## Pointers

- `.agents/reports/2026-06-08_white-matter-exclusion.md` — full report.
- `locked/CONTROLS.md` §C6 — locked control definition.
- `locked/VERDICT_LEDGER.md` — 2026-06-08 revision entry.
- `data/audit/wm_stratified/{cophenetic_raw_cohort,grassmann_cohort}.csv` — numbers.
- `data/audit/wm_stratified/figures/fig_wm_{verdict_matrix_coph,verdict_matrix_raw,exclude_forest,grassmann_kspan}.pdf`.
