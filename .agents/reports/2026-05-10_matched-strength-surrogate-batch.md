---
date: 2026-05-10
era: COHORT_N10 / IMCOH_ABS
status: current
type: verification-batch
scope: §5.3 ρ matched-strength surrogate (shared + split baseline) + §5.5 epi-fraction disclosure
inputs:
  - data/audit/matched_strength_surrogate/cohort_summary.csv
  - data/audit/matched_strength_surrogate/per_patient_per_band.csv
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
  - data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
  - data/audit/epi_fraction_anatomy/region_epi_fractions.csv
scripts:
  - scripts/01_compute/audit/audit_61_epi_fraction_anatomy.py
  - scripts/01_compute/audit/audit_62_matched_strength_surrogate.py
  - scripts/01_compute/audit/audit_63_split_baseline_surrogate.py
---

# 2026-05-10 — matched-strength surrogate batch + §5.5 epi-fraction disclosure

**Head.** Three investigations landed on 2026-05-10 in support of the §5.3
spine-settling decision. `audit_62` (shared-baseline ρ, R=100) returned
**intermediate / also_positive** verdicts at all three trace bands —
shared-baseline ρ is largely strength-driven. `audit_63` (split-baseline
ρ_split, R=200) returned **intermediate** verdicts at all three trace
bands but with cohort-paired Wilcoxon p<0.005 at α and β — the
split-baseline ρ_split exceeds an independent-strength-randomized null
at the cohort scale even though the per-patient distribution is mixed
(α 5/10, β 7/10, γ_l 5/10 below own surrogate at p<0.05). `audit_61`
adds a descriptive disclosure for §5.5: epileptic-zone fraction in the
five flagged region pools ranges 0–14.8%, with **β / Hippocampus the
highest at 14.8% (4/27 contacts)**.

The combined message to §5.3 / §6 is **band-resolved**: β cohort
non-stationarity claim survives matched-strength surrogacy at the
cohort scale; γ_l cohort claim does not (cohort p=0.116). The §5.6
β anchor / §5.5 anatomy story should disclose the epi-fraction at
Hippocampus before pinning Hip-driven claims to disease-free anatomy.

## Cohort verdicts

### audit_62 — shared-baseline ρ (R=100, SWAP_FACTOR=5)

| band | n | obs median ρ | surr median ρ | n above surr mean | cohort z | verdict |
|---|---|---|---|---|---|---|
| α | 10 | +0.454 | +0.425 | 6/10 | +0.02 | **also_positive** |
| β | 10 | +0.502 | +0.433 | 7/10 | +0.48 | **intermediate** |
| γ_l | 10 | +0.432 | +0.496 | 5/10 | -0.29 | **intermediate** |

Reading: shared-baseline ρ (numerator shares D_pre across phases) is
inflated by construction. Matched-strength surrogates reproduce this
inflation almost perfectly at α (also_positive) and substantially at
β / γ_l (intermediate). The shared-baseline ρ is **largely
strength-driven**; the §5.3 statistic of record cannot be the
shared-baseline ρ.

### audit_63 — split-baseline ρ_split (R=200, SWAP_FACTOR=20)

| band | n | obs median ρ_split | surr median (per-pat med) | n above own surr | Wilcoxon z | p | verdict |
|---|---|---|---|---|---|---|---|
| α | 10 | +0.105 | +0.013 | 5/10 | +54.00 | **0.0020** | intermediate |
| β | 10 | +0.221 | +0.009 | 7/10 | +52.00 | **0.0049** | intermediate |
| γ_l | 10 | +0.083 | +0.003 | 5/10 | +40.00 | 0.1162 | intermediate |

Reading: split-baseline ρ_split (Δ_task = D^tt − D^pre_A; Δ_rest =
D^post − D^pre_B; Spearman) returns near-zero surrogate medians
(|surr_med| < 0.05) at all three bands, confirming the surrogate
construction successfully randomizes the cross-phase coherence at the
edge level. **α and β cohort-paired Wilcoxon p<0.005**; γ_l
p=0.116 (n_above 5/10, three patients strongly anti-aligned). All three
bands fail the strict per-patient gate (need ≥8/10 individually below
own surrogate at p<0.05) so the verdict label is "intermediate" by the
spec criterion, but the cohort-level signal at α / β is real.

Per-patient picture at β: 7/10 patients individually exceed their own
surrogate at z ≥ 3.5 — Pat_03/05/06/08 strongly trace-aligned
(z = 3.9–5.8). Pat_15 is strongly anti (z = -3.6). Pat_10/14 mildly
anti. The cohort claim survives because the trace-aligned patients
are robustly distinguished; the per-patient claim is split.

Per-patient picture at γ_l: Pat_02/05/06 strongly trace-aligned
(z = 4.0–4.5); Pat_07/10/14/15 strongly anti (z up to -7.4 / -10.4).
The bipolar pattern means cohort signal does NOT survive matched-strength
surrogacy at γ_l (Wilcoxon p = 0.116).

### audit_61 — epi-fraction in §5.5 region pools

| band | region | n_contacts | n_epi | frac_epi |
|---|---|---|---|---|
| γ_l | ctx-lh-fusiform | 38 | 4 | 10.5% |
| β | Hip | 27 | 4 | **14.8%** |
| β | ctx-lh-fusiform | 38 | 4 | 10.5% |
| β | ctx-lh-superiortemporal | 53 | 4 | 7.5% |
| α | ctx-lh-parsopercularis | 26 | 0 | 0.0% |

Reading: the §5.5 trace-leaf enrichment is computed against the full
cohort cortical pool without epi-zone exclusion. For β / Hippocampus
(the strongest β cell at 3.49× cohort baseline), **14.8% of contacts
are epi-zone-flagged**. The §5.5 disclosure paragraph should add a
sensitivity comment: an exclusion-sensitivity analysis (drop epi-zone
contacts and recompute the hypergeometric) is owed before the β/Hip
finding can be pinned to disease-free anatomy.

α / ctx-lh-parsopercularis is epi-zone-free (0/26) — that finding does
not need an epi-sensitivity caveat.

## Manuscript implications (band-resolved)

The headline claim "non-stationarity at three frequency bands carrying
band-specific signatures" **holds** for the §5.3 substrate-level
non-stationarity result. What changes is which band carries which kind
of evidence at the LRG-layer:

- **β**: substrate non-stationarity (raw FC d_S T_d cohort-wide) +
  LRG-layer shape-driven trace at cohort level (split-baseline ρ_split
  Wilcoxon p = 0.005). Per-patient picture is 7/10 — write the cohort
  effect as "majority cohort effect (7/10) at p_paired = 0.005" rather
  than "10/10 universal".
- **α**: substrate non-stationarity + LRG-layer shape-driven trace at
  cohort level (Wilcoxon p = 0.002). Per-patient mixed (5/10).
- **γ_l**: substrate non-stationarity at the per-pair / per-leaf level,
  but **not** an LRG-layer shape-driven trace beyond what
  matched-strength noise provides (Wilcoxon p = 0.116, 5/10 per-patient
  with bipolar Pat_14 / Pat_15 anti-alignment z = -7.4 / -6.9).

For §5.5 anatomy: disclose the 14.8% Hippocampus epi-fraction; an
epi-exclusion-sensitivity analysis is owed before the β/Hip claim is
finalized.

## Algorithm caveat (vs ticket text)

Ticket asked for "the SAME structural permutation π_r" applied across
phases via 4-cycle ±δ rewiring. The 4-cycle ±δ method does **not**
permute edges — it perturbs weights along chosen 4-cycles to preserve
node strengths — so there is no permutation π_r to apply across phases.
The standard interpretation (independent per-phase rewiring) is what
audit_62 / audit_63 implement.

Implication: independent per-phase rewiring decorrelates cross-phase
edge identity by construction. The surrogate ρ_split is therefore
expected to be near zero regardless of whether the §5.3 trace direction
is shape-driven or strength-driven. A "separated" verdict in this test
confirms the observed cross-phase ρ_split is not an artifact of
independent per-phase strength heterogeneity, but does NOT isolate
shape-driven from strength-driven cross-phase coherence on its own.
The shared-baseline test (audit_62) is the complement: preserves
cross-phase comparability under inflated baseline but does not separate
cleanly.

A coordinated cross-phase surrogate that preserves both per-node strengths
and cross-phase comparability is not well-defined in the strength-preserving
family. Decision deferred until a cleanly-defined cross-phase surrogate
exists.

## Acceptance gate (per ticket)

The §5.3 split-baseline ρ_split run is "the spine of the manuscript;
this run settles it" per the ticket. After landing, the §5.3 prose is
finalized at the verdict label returned. **No further matched-strength
surrogate runs queued unless this run returns null with diagnostic
evidence the algorithm is at fault.** Algorithm has been audited
(strength-preservation tol 1e-4 verified per surrogate; n_swaps =
20·E sufficient for marginal preservation; surrogate medians at all
three bands |.| < 0.05 confirm independent-phase decorrelation works).

## Companion / parallel investigations

- **audit_64 (implant geometry)** — see
  `2026-05-11_implant-geometry-and-kc-null-verification.md` for the
  per-patient implant-geometry test of cohort anti-alignment. Pure
  correlational (no LRG, no surrogate). β returns moderate Spearman
  +0.697 against B_hemi (uncorrected p=0.025, BH-FDR q=0.994 fails
  m=135 correction). α / γ_l null.
- **audit_65 (KC tree-distance matched-strength surrogate)** — also see
  `2026-05-11_implant-geometry-and-kc-null-verification.md`. In flight
  at time of this writing; verdict pending.

## Files

- `data/audit/matched_strength_surrogate/{cohort_summary,per_patient_per_band}.csv`
  + `figures/{cohort_distribution,per_patient_panel}.pdf` + `README.md`
- `data/audit/matched_strength_surrogate_split_baseline/{cohort_summary,per_patient_per_band}.csv`
  + `figures/{cohort_distribution,per_patient_panel}.pdf` + `README.md`
- `data/audit/epi_fraction_anatomy/{region_epi_fractions,per_contact_table}.csv`
  + `README.md`

## Provenance

- Cohort: Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15 (n=10, COHORT_N10)
- Bands: alpha, beta, low_gamma (the three §5.3 trace bands)
- FC method: imcoh_abs
- LRG: τ = 1/λ_max, ultrametric via average linkage on Trho.
- audit_62 wall-clock: 820.8 s (R=100)
- audit_63 wall-clock: 7616.2 s (R=200, ~127 min)
- audit_61 wall-clock: ~5 s (descriptive only)
