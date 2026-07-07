---
name: 2026-07-06_rho-sym-gate-marginality-flag
type: report
era: IMCOH_ABS_COHORT_N10
status: current
created: 2026-07-06
updated: 2026-07-06
pointers:
  - memory/feedback_rho_sym_gate_marginal_not_clean_bands.md
  - .agents/preprint/bands/00_cohort.md
  - .agents/preprint/locked/VERDICT_LEDGER.md
audience: headline-agent
---

# Flag: the ρ_sym cophenet gate is marginal, and its per-band p-values compress vs ρ_split

**Head.** Good news first: the preprint taxonomy is **already graded and
dual-probe** — there is NO clean "α/β trace vs rest no-trace" binary to retract
(`bands/00_cohort.md:29`: β both probes / α cophenet-only / γ_l+δ Grassmann-only
/ γ_h+θ neither). The real item to flag for the ρ_sym number migration: the
per-band **cophenet cohort-gate p-values** cited in the briefs are `ρ_split`,
and under `ρ_sym` they **compress** — γ_l/γ_h/δ move from clearly-failing
(0.12–0.28) to a **marginal p=0.080** (one patient off the p<0.05 line). No
verdict flips, but "γ clearly fails cophenet" should become "**marginal**".

## The compression (cohort gate p, `data/audit/rho_sym_gate/cohort_summary.csv`)

| band | ρ_split gate p (in briefs) | ρ_sym gate p (canonical) | note |
|---|---|---|---|
| α | 0.0029 | **0.024** | strong band, less extreme |
| β | 0.0049 | **0.032** | strong band, less extreme |
| γ_l | 0.116 | **0.080** | now marginal (was clearly failing) |
| γ_h | 0.246 | **0.080** | now marginal |
| δ | 0.216 | **0.080** | now marginal |
| θ | 0.722 | **0.784** | ~unchanged, absent/anti |

ρ_sym pulls both ends toward the middle (symmetrizing the arbitrary A/B half
regresses the ρ_split extremes). β/α still clear p<0.05; γ_l/γ_h/δ still fail
the 0.05 gate — **the classification does not flip** — but the γ/δ cophenet
margin is now thin. At n=10 the signed-rank p is discrete: γ/δ sit at `W+=42`,
and `W+=45` (p=0.042) passes → ~1 patient from crossing. β itself is only 7/10
positive (same count as δ/γ_h) and clears on **magnitude**, not consistency.

## What this does and does NOT change

- **Does not change** any verdict tag. β = trace both probes; α = cophenet-only;
  γ_l = Grassmann-only; δ = weak Grassmann-only; γ_h/θ = no trace. All hold.
- **Does not touch** the Grassmann probe (independent of the ρ_split→ρ_sym
  cophenet estimator) or the β→OFC localization.
- **Does change wording**: anywhere a brief says γ_l/γ_h/δ cophenet is "clearly
  no trace / p≈0.12–0.25 / indistinguishable from surrogate", the ρ_sym number
  is **p=0.080 (marginal)**. Frame as "marginal on the cophenet gate, verdict
  carried by the Grassmann probe" (γ_l/δ) or "marginal on both" (γ_h).

## Specific locations carrying ρ_split cophenet numbers (revisit on migration)

- `bands/00_cohort.md` — lines ~29, 40–43, 106, 108–109, 112 (per-band cophenet
  p-values + the summary taxonomy table).
- `bands/03_gammalow.md` — C3 cophenet `p = 0.116`, `n_above 5/10` (several spots).
- `bands/05_gammah.md` — C3 cophenet `p = 0.246`, `4/10` (several spots).
- `bands/04_theta.md` — cophenet anti framing (θ ρ_sym p=0.784, still absent — safe).
- `locked/VERDICT_LEDGER.md` — already has the 2026-07-06 α/β ρ_sym swap; γ/δ C3
  rows still cite ρ_split.

## Recommended reading (already on the headline-agent TODO per memory)

The single ρ_sym cohort gate cannot license band separation on its own — present
graded tiers (**clear** α/β · **marginal** γ_l/γ_h/δ · **absent** θ) + the
per-patient spread. β's flagship status rests on **magnitude + multi-probe
convergence** (cophenet gate AND Grassmann AND OFC), not the gate. Figure:
`scripts/01_compute/figures_embedded/fig_rho_sym_band_map.py`.
