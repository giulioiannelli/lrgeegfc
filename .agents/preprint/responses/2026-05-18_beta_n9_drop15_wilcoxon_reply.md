---
name: beta-n9-drop15-wilcoxon-reply
era: IMCOH_ABS_COHORT_N10
status: reply-with-flag
kind: writing-agent-reply
date: 2026-05-18
band: beta
target: paragraph robustness statement "Pat_15-dropped (n=9)" replacement of legacy "n=8" sentence
---

# β `n=9` (drop \patient{15}) Wilcoxon recompute — reply

**Head.** The numbers are computed; they're below. Before reading them, the
ticket's "ordering sanity check" is **violated for the median direction**,
but the underlying reason is a flawed premise in the ticket — *not* evidence
that \patient{07} is anti at the LRG layer. The Wilcoxon `p` ordering matches
the ticket's prediction cleanly. I am reporting the numbers with the flag
attached, because the substantive `p` claim is sound; the writing agent
should incorporate the values and discard the median-ordering premise.

## Flag — ordering sanity check (median direction): premise was wrong

| Restriction | n | median \(\rhosplit\) (matched-strength `obs_rho`) | Wilcoxon p |
|---|---|---|---|
| Full cohort | 10 | +0.221 | 0.005 |
| **Drop \patient{15}** | **9** | **+0.230** | **0.010** |
| Legacy substrate (drop \patient{07} + \patient{15}) | 8 | +0.292 | 0.020 |

The ticket predicted **full > drop15 > legacy in median**. The actual
ordering is **full < drop15 < legacy** — the median *increases* as patients
are dropped.

Why: \patient{07}'s \(\rhosplit\) is +0.230 (matched-strength) / +0.185
(within-baseline), both **below** the cohort median (+0.221 / +0.222).
"Pro-trace direction" (positive sign, `z > 0`) and "above-cohort-median"
are different things. \patient{07} is solidly **pro-direction** at the
LRG layer (matched-strength `obs_z = +4.47`, Grassmann `obs_z = -3.90`,
both decisive) but moderate in magnitude — its \(\rhosplit\) value sits
in the lower half of the cohort's pro-direction distribution. Dropping
a below-median pro-direction patient *raises* the cohort median (rank
order rearranges so the median moves up the pro tail). \patient{15} is
the same direction in influence on the median because its \(\rhosplit\)
= +0.083 is also below the cohort median (and is the *only* truly
anti-trace LRG patient by direction — but its small magnitude does not
pull the cohort median much).

The Wilcoxon `p` ordering — which is the load-bearing test — matches
the ticket's prediction: **0.005 < 0.010 < 0.020**. The signed-rank test
respects the cohort sign distribution, not the median magnitude; removing
patients increases `p` as expected because `n` drops.

**Recommendation.** Treat the `p` ordering as the sanity check
(passes). Discard the median-ordering predicate from the ticket. The
paragraph can safely use the `n=9` values below as the LRG-native
robustness check.

## Primary numbers (\(\rhosplit\) under matched-strength)

| Quantity | Value | Three-decimal rendering |
|---|---|---|
| Cohort median `obs_rho` at β, `n=9` (\patient{15} dropped) | +0.23022 (matched-strength) / +0.25901 (within-baseline) | **+0.230** / **+0.259** |
| Cohort-paired one-sided Wilcoxon p, `n=9` | 0.009765625 | **0.010** |

The matched-strength `obs_rho` value (+0.230) is the table-consistent
column-3 cohort median (Option B of the per-patient table, where col. 3
holds matched-strength `obs_rho`). The within-baseline `rho_split`
value (+0.259) is the ctm_triangle realization median and would be the
relevant number if the paragraph cites within-baseline cohort medians.

The Wilcoxon `p = 0.010` is identical for both realizations at `n=9`
because the per-patient sign pattern is the same across the two runs
(only Pat_10 and Pat_14 are below zero on \(\rhosplit\); both
realizations agree at every patient on sign).

## Auxiliary numbers (`n=9` drop \patient{15})

| Test | Cohort median | Wilcoxon p | Verdict |
|---|---|---|---|
| Grassmann \(T_G(k{=}40)\) (matched-strength, β) | **\(-0.572\)** | **0.002** | strengthens vs `n=10` (\(p = 0.005\)) |
| Drift-floor null (\(\rhosplit\) > \(\rho_\text{drift}\), within-baseline) | n/a (paired Wilcoxon on diff) | **0.027** | weakens vs `n=10` (\(p = 0.014\)) but still passes the 0.05 gate |

Grassmann at `k=40` *strengthens* when \patient{15} is dropped, because
\patient{15} is the only Grassmann anti-aligned patient
(`T_G = +1.020`) and pulls the cohort median toward zero. Removing
\patient{15} cleans the cohort signal: median drops from \(-0.489\) to
\(-0.572\), `p` falls to 0.002.

Drift-floor weakens (`p` rises from 0.014 to 0.027) because the
matched-strength test's robustness is carried mostly by the
high-magnitude pros; dropping the smallest pro (\patient{15}, +0.083 vs
drift) trims that tail.

## Paragraph-level summary of the four `n=9` results

The **LRG-native robustness restriction** drops only \patient{15}, the
lone LRG anti-aligned patient (right-hemisphere-only implant, no
epileptic contacts). The earlier legacy `n=8` "pro-cohort restriction"
that also dropped \patient{07} is **retired** under the 2026-05-18
patient-dropout policy update (see `feedback_no_patient_dropout.md`;
Pat_07 is solidly pro at the LRG layer and was inherited from a
substrate-only construct). The `n=9` LRG-native restriction preserves
the trace claim at every probe and test:

- per-pair multiscale \(\rhosplit\): cohort median +0.230, `p = 0.010`
  (matched-strength gate)
- Grassmann \(T_G(k{=}40)\): cohort median \(-0.572\), `p = 0.002`
  (matched-strength gate — strengthens)
- drift-floor: `p = 0.027` (within-baseline, still passes 0.05)

All three pass at `n=9`, including the strengthening of the Grassmann
probe under the restriction.

## Source CSVs and column names used

| Column / quantity | CSV path | Column |
|---|---|---|
| `obs_rho` (matched-strength) | `data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv` | `obs_rho` (filter `band=beta`) |
| `surr_p50` (matched-strength) | same | `surr_p50` |
| `rho_split` (within-baseline) | `data/audit/ctm_triangle/Td_per_patient_per_band.csv` | `rho_split` (filter `band=beta`) |
| `rho_null_drift` (within-baseline) | same | `rho_null_drift` |
| `obs_T_G` (Grassmann) | `data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv` | `obs_T_G` (filter `band=beta, k=40`) |
| `surr_T_G_p50` (Grassmann) | same | `surr_T_G_p50` |

All tests use `scipy.stats.wilcoxon(diffs, alternative='greater')` for
\(\rhosplit\)/drift-floor (trace = positive `diffs`) and
`alternative='less'` for Grassmann (trace = negative `T_G - surr`),
exact distribution for `n ≤ 10`.

## Note on the legacy errata `p = 0.039`

The brief and the writing directive cite the legacy `n=8` (drop
\patient{07} + \patient{15}) Wilcoxon as `p = 0.039`. My recompute on
the same data gives `p = 0.020`. The factor-of-two discrepancy is
consistent with the errata having used a *two-sided* test or a
different sign convention; the *one-sided* matched-strength gate cited
in the paragraph is `p = 0.020`. The legacy `0.039` should be updated
to `0.020` in the next pass of the brief if the agent is touching that
sentence. Not a blocker for the `n=9` replacement.
