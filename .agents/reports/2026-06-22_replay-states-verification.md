---
name: replay-states-verification
era: IMCOH_ABS_COHORT_N10
status: CLOSED — transient negative; sustained folded into N1.7 as headline
kind: report
scope: N4 (replay states) executed + verified — verdict + full method/result ledger
headline: .agents/preprint/headlines/archive/2026-06/04_replay_states.md (archived) → folded into 01_trace.md §N1.7
scope_doc: .agents/guides/task-persistence-investigation/2026-06-22_replay-states-v3-scale-target-pair.md
created: 2026-06-22
closed: 2026-06-23
---

# N4 replay states — verification ledger (CLOSED 2026-06-23)

> **⛔ INVESTIGATION CLOSED 2026-06-23 (PI).** Burst/transient replay state = exhaustive
> negative (5 propagator representations × all scales × all targets × all bands × all
> timescales down to the 0.2 s high-γ coherence floor; audit_124–141). The **sustained
> reinstatement** is the kept result and is now the **headline** — folded into N1 as
> **§N1.7** (`01_trace.md`): the β trace is a *state the brain holds*, dynamically,
> scale-invariantly, across 5 representations. N4 headline archived. No further replay
> compute planned (single ~100 ms ripple is below the coherence floor of ANY connectivity
> method → out of reach of this data, not a gap to fill).

**Head.** N4's distinctive claim — the task config RE-ENTERS rest as **discrete,
transient, recurring replay STATES** — is a rigorous cohort NEGATIVE across every
in-framework method and every measurable timescale/band. What IS strongly true is
**sustained reinstatement**: rest_post DWELLS in the task connectivity configuration
(window-level, **10/10, p=0.001**) — i.e. N1 expressed dynamically as occupancy, a
*state the brain holds*, not a flashing replay process. **DONE: folded into N1.7 as
the headline; transient-replay framing dropped.**

## Window pipeline anchor (audit_124)
Windowed |ImCoh| recipe == canonical `imcoh_abs` EXACTLY (Spearman=Pearson=1.000),
so the time-resolved ρ^coph reduces to N1 in the full-window limit. Lib promoted:
`lrg_eegfc.utils.fc.coherence.windowed` (segment_ffts, imcoh_abs_cube, band_abs_average).

## Result ledger (all n=10 unless noted)
| audit | test | timescale/band | replay STATES? | sustained level? |
|---|---|---|---|---|
| 124 | continuity anchor + windowability | — | (gate) | recipe==N1 ✓ |
| 124b | downstream reliability | 20s α/β | per-window s ~constant (sd~0.01–0.05) | — |
| 125 | template burstiness (AC1 vs time-shuffle) | 20s α/β | **NO** (cohort p=0.50; pre≥post) | — |
| 126 | k-means state discovery | 20s α/β | **NO** (lift~1.2, p=0.31–0.54) | occupancy↑ |
| 126b | shift-vs-states decomposition | 20s α/β | **NO** (occ=100% uniform shift; spread −) | — |
| 127 | four-phase contrasts (incl inference test−learn) | 20s α/β | **NO**; a_infer/g_IE null/neg (no inference-pref replay) | g_TP 10/10 p=0.001 |
| 128 | short-window precondition | 0.2–20s | (gate) | **α/β unmeasurable <10s**; only γ faster |
| 129 | sustained/stabilized (Dir A) | 20s α/β | — | **A1 10/10 p=0.001**; A2 tight null; A2b settle p=0.08 |
| 130 | fast replay, γ combined | 2s γ | **NO** (pilot 3/3 AC1 n.s.; shift heterog.) | shift big in Pat_05 |
| 131 | fast replay, γ **low vs high** | 2s γ | **NO** (low γ all state-neg incl anti-burst; high γ excess-flash 8/10 p=0.065 SUB-THRESH) | shift low γ 7/10 p=0.065, high γ 6/10 p=0.19 |
| 132 | localized OFC (Dir B1) | 20s α/β | burstiness null (2/20 cells) | **OFC #1 β system both epi modes; epi-excl 4/5 p=0.094** (underpowered, K=5; α null) |
| 133 | **τ-resolved (v3 R1)** — burstiness PER diffusion scale, time-shuffle + node-perm placebo | 20s α/β, α=0.5→19/λmax | **NO at EVERY scale** (burst best α p=0.16/β p=0.19; specific null; clean coverage frac_collapsed=0 to α≈13) | **shift POSITIVE at ALL τ 8–10/10 p≤0.007 = sustained hold is SCALE-INVARIANT** (anchor α=1 ✓ 10/10) |
| 134 | **dual-target trajectory (v3 R2)** — test-vs-learn drift over rest | 20s α/β | **NO bursts** (0/10 indiv both bands) | α slow drift test→learn 8/10 two-sided p≈0.06 (footnote, α not β, drift≠state); β flat |

## v3 closure (audit_133/134): the three collapsed axes are now open
PI hypothesis "the transient hides in the τ value" → **falsified, cleanly**, across the
WHOLE non-collapsed scale axis (finest→Fiedler). Combined with R2 (target axis: no
transient toward learn or test) and B1 (per-pair/OFC subtree burstiness null), **all three
axes prior tests collapsed are now opened and the transient is absent in every one.** The
SEQUENCE axis is logically foreclosed too: sequential replay needs distinguishable
transient states to order, but k-means (126) found none separable. NEW positive: the
sustained hold is **scale-invariant** (R1 shift 8–10/10 at every τ). α drift = footnote.

## Propagator-PORTFOLIO closure (audit_136–139): PI "don't stop until found within a
## Laplacian-propagator approach" → 5 DISTINCT propagator representations, all converge
| P# | representation | transient STATE? | sustained hold? |
|---|---|---|---|
| P-std (133) | combinatorial L, cophenetic, all τ | ✗ clean | ✓ scale-invariant |
| P1 (136) | **magnetic L_H=D̄−iA, DIRECTIONAL** | ✗ (placebo caught coarse-τ pseudo-burst, specific_p≈0.5) | ✓ + **directional-flow candidate** (shift α 10/10 p=0.001) |
| P5 (137) | **raw propagator 1/ρ̂(τ)**, no dendrogram | ✗ (quantization NOT the culprit) | ✓ strongest (β +1.2, 10/10 every τ) |
| P4 (138) | **subspace / Grassmann** (k modes) | ✗ (no k* target-specific) | ✓ (shift α/β k=16 10/10/9-10) |
| P3 (139) | normalized L_sym=I−D^{-1/2}WD^{-1/2} | ✗ (best burst α p=0.053, not specific) | ✓ 8–10/10 every τ |

**Net: replay-STATE hypothesis EXHAUSTED across the propagator space.** The obstacle is NOT
the representation — it is an UPSTREAM measurement floor: the α/β trace is intrinsically slow
(≥10 s window to estimate |ImCoh|), longer than any replay event, so windowing averages out
fast dynamics BEFORE any propagator sees them; per-window task-likeness is near-constant
(audit_124b sd~0.01–0.05). No propagator variant is downstream-fixable for that.
**NEW POSITIVE from the chase:** sustained reinstatement is multi-representation robust +
scale-invariant (strengthens N1). The directional-flow candidate (P1) was tested with a
magnitude-matched control (audit_140: L_H=D̄−iA vs L_mag=D̄−|A|, same D̄) and **RETRACTED** —
`dir_reinstated` significant (α=1 10/10 p=0.001) but so is `mag_reinstated` (9/10 p=0.014)
and the `directional_excess` (dir−mag) is NULL (med +0.013, 7/10, p=0.138 best; negative at
α=2 & β). So the magnetic reinstatement is magnitude-inherited, NOT directional-specific.

## FORECLOSURE: why no propagator approach can find the transient on α/β
The transient is blocked UPSTREAM of every propagator: α/β need ≥10 s windows to estimate
|ImCoh| (audit_128), longer than any replay event, so windowing averages the fast dynamics
away BEFORE any operator (cophenetic/raw/magnetic/subspace/normalized/Hodge) sees the data;
per-window task-likeness is near-constant (sd~0.01–0.05). This is a measurement floor, not a
modeling choice — it forecloses a transient α/β replay STATE for ALL propagator-based methods,
Hodge curl included (curl is computed from the same windowed signed ImCoh → same floor). Only
γ is fast-resolvable and γ showed no states (131) and is not a trace band. Last genuinely-
distinct lead not yet RUN: P2 Hodge curl propagator (best epi prior, but downstream of the
same floor; expensive Hodge-L1 build) — run-on-request, predicted to confirm.

## Sub-second high-γ — the ripple regime (audit_141, PI pushback 2026-06-23, CORRECT)
The ≥10 s floor is α/β-ONLY; high γ is feasible to 0.2 s (audit_128 n_bins=3@0.2s/14@0.5s) =
the actual ripple-replay timescale (80–150 Hz, ~50–150 ms) — audit_130/131 stopped at 2 s out
of habit. audit_141 swept high+low γ at L∈{0.2,0.5,1.0,2.0}s (3000+ windows @0.2s), burstiness
+ time-shuffle + node-perm placebo: **NO replay state at any sub-second window.** Clustered
burst null everywhere (best p=0.31, target-non-specific); the 2 s high-γ excess-flash hint
(8/10 p=0.065) did NOT sharpen toward the ripple scale — non-monotonic (0.2s null, 1.0s 7/10
p=0.080 sub-thresh, 2.0s null) = noise shape, not a timescale-locked signal. REAL RESIDUAL
(honest): a single ~100 ms ripple is below the COHERENCE floor itself (coherence needs ≥several
segments × several cycles ⇒ ~0.2 s min even for high γ) → single-EVENT replay unmeasurable by
ANY connectivity method, not just LRG; but the whole RESOLVABLE fast range (0.2–1 s) is clean
negative. ⇒ foreclosure now EARNED at the ripple timescale. Honest floor = "0.2 s for high-γ
connectivity," not "20 s"; tested, no state.

## Key conclusions
1. **No discrete replay states** at 20s (α/β) by burstiness, state-clustering, or
   shift-decomposition; rest_post is *less* variable than rest_pre (anti-states).
2. **No preferential inference (test-over-learn) replay** (audit_127 a_infer/g_IE).
3. **Sustained reinstatement is real & strong** (g_TP / A1: 10/10, p=0.001, LOO 0.002)
   = N1 window-level. "Tighter attractor" only a weak trend (A2b p=0.08).
4. **The β/α trace is intrinsically SLOW**: imaginary coherence at those freqs needs
   ~10s; trace bands cannot be resolved <10s, so fast replay of the trace is not
   *measurable* (audit_128). Only γ is fast-measurable, and γ@2s shows a (heterogeneous)
   stationary shift, **no burstiness** (audit_130; per-band low/high in audit_131).
5. **γ (n=10, audit_131) gives NO replay states either.** Low γ is flatly
   state-negative (AC1 p=0.95, even *anti*-bursty; excess-flash p=0.90); high γ shows the
   ONLY whisper — isolated excess-flash 8/10 same-sign but p=0.065 (sub-threshold, dies
   under correction), and its per-patient AC1 "hits" (Pat_06/07/15) look like slow
   vigilance drift (big AC1, ~zero shift) not discrete replay. Both γ bands carry a weak
   heterogeneous stationary shift (the N1 level), low γ 7/10 p=0.065.
6. **B1 localization (audit_132): the sustained state has the trace's anatomy but is
   underpowered.** OFC is the #1 β system for the post>pre dwelling in BOTH epi modes
   (epi-excl conc=+0.027, 4/5, p=0.094, LOO 0.188; epi-incl 3/5, p=0.22); α null at OFC
   (band-specific). Not significant — OFC implanted in only K=5, so p<0.05 needs 5/5.
   Within-OFC burstiness null (2/20). Matched-strength moot (no significant claim).

## OPEN / RESUME (all compute DONE — only editorial remains)
- ALL planned tests executed (audit_124–132). γ (131) and B1 OFC (132) are final;
  CSVs in `data/audit/replay_states/` (gamma_perband_*, localized_ofc_*).
- The ONLY sub-threshold flicker is high-γ excess-flash (8/10, p=0.065). A single
  matched-strength referee at high-γ excess-flash could formally kill-or-confirm it, but
  it will not survive correction; recommended NOT to chase (KC-style trap). PI's call.
- **Pending editorial (PI decision):** reframe headline 04 N4 "replay states" →
  "sustained reinstatement" and fold into N1; mark frontmatter executed/not-supported.
  NOT done unilaterally — PI is invested in N4 as the "central idea."

## Provenance (CSVs under data/audit/replay_states/ and data/audit/replay_windowability/)
burstiness_summary.csv, cohort_verdict.csv, state_decomposition_cohort.csv,
shift_vs_states.csv, fourphase_cohort.csv, sustained_reinstatement_cohort.csv,
short_window.csv, shortwindow_gamma_per_patient.csv, gamma_perband_*.csv.
Scripts: audit_124,124b,125,126,126b,127,128,129,130,131 in scripts/01_compute/audit/.
