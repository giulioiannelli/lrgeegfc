---
name: replay-sustained-and-timescale
era: IMCOH_ABS_COHORT_N10
status: scope_executing
kind: scope-report
measure: replay states v2 — (A) sustained/stabilized reinstatement; (B) localized OFC + multi-timescale transient
headline: N4 (.agents/preprint/headlines/04_replay_states.md)
supersedes_framing: transient-bursts-only (audit_125/126/127 = NEGATIVE for discrete flashes)
created: 2026-06-22
---

# Replay states v2 — sustained reinstatement + multi-timescale transient (in-framework)

**Head.** Transient-flash replay (scoped N4) is a rigorous cohort NEGATIVE at 20 s
windows across 4 methods (audit_125/126/126b/127): rest_post is task-like via N1's
STATIONARY uniform shift, not bursts; no preferential inference replay. PI redirect
(2026-06-22): pursue two in-framework directions instead. Stay in ρ^coph on |ImCoh|
windows (recipe == N1, anchored audit_124).

## Direction A — Sustained, stabilized reinstatement (REFRAME of the trace)

**Claim.** "Replay state" = rest_post DWELLS in the task connectivity state: it sits
closer to the task config (level) AND more tightly/stably there (lower window-to-
window dispersion = a tighter attractor) than rest_pre. This is the trace (N1)
deepened into a state statement, NOT transient flashes.

**Tests (cohort Wilcoxon one-sided, LOO-max; from existing 20 s cophenetic cache):**
- A1 LEVEL: per-window task-similarity rho^coph(c_w, c_test) higher in post than pre.
  (Have it: g_TP shift 10/10, p=0.001.)
- A2 STABILIZATION: within-phase window dispersion LOWER in post than pre. Dispersion =
  median pairwise (1 - Spearman) among a phase's window cophenetic vectors (how tightly
  the windows cluster). post < pre ⇒ tighter task-state attractor.
- A3 ATTRACTOR: rest_post windows are closer (rho^coph) to the task centroid than
  rest_pre windows are (mean window→task-centroid distance, post < pre).
- Caveat: A1 ≈ N1; A2/A3 are the genuinely new "it's a STATE not just a shift" content.
  Honest framing: sustained reinstatement, an attribute OF the trace.

## Direction B — Localized (OFC) + MULTI-TIMESCALE transient

**Claim.** Even if whole-brain config is stationary at 20 s, a TRANSIENT flash may live
(i) in the trace's anatomical carrier (OFC system, N1.3), and/or (ii) at SHORTER
timescales / faster bands than 20 s probed.

**B0 — short-window precondition (GATE, audit_128).** Sweep window L ∈ {0.2, 0.5, 1, 2,
5, 10, 20} s, choosing nperseg per L to keep ~12–16 Welch segments, and report per-band
split-half |ImCoh| reliability + frequency-resolution feasibility (n_bins in band ≥ 2).
EXPECTED physics: α/β unmeasurable below ~10 s (band only a few Hz wide vs coarse df);
only γ_low/γ_high survive at ≤2 s. Confirm empirically; restrict B to feasible (band, L).

**B1 — localized OFC task-likeness (audit_129).** Restrict the per-window cophenetic
task-likeness to OFC-system pairs (anatomical carrier of the β trace, N1.3). Test
burstiness (AC1 vs time-shuffle) + post-vs-pre, at the feasible window lengths.

**B2 — fast-band transient at short L (audit_130, only for feasible bands).** For bands
that survive B0 at short L (expected γ), recompute per-window cophenetic + task-likeness
at the short L and run the burstiness make-or-break. This is where fast (ripple-like)
replay would appear if it exists.

**Nulls/controls (unchanged):** time-shuffle (burstiness), shift-vs-states decomposition
(audit_126b logic), rest_pre within-subject control, specificity (shuffled target).
Matched-strength only if a POSITIVE survives (referee for a claim, not a null result).

## Execution order
1. B0 short-window precondition (audit_128) — establishes the timescale floor per band.
2. A1–A3 sustained reinstatement (audit_129) — from existing cache, fast.
3. B1 localized OFC (needs OFC node sets) + B2 fast-band short-L (gated by B0).

## Provenance / prior
- Negative transient verdict: data/audit/replay_states/{burstiness_summary,cohort_verdict,
  state_decomposition_cohort,shift_vs_states,fourphase_cohort}.csv (audit_125/126/126b/127).
- Windowed |ImCoh| lib: lrg_eegfc.utils.fc.coherence.windowed (segment_ffts, imcoh_abs_cube,
  band_abs_average); cophenetic: utils.surrogate.matched_strength.cophenetic_condensed_from_adjacency.
- OFC localization: data/audit/localization_atlas/ (N1.3, audit_83/92).
