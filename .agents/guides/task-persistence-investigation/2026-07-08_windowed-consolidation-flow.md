---
name: windowed-consolidation-flow
type: scope
era: IMCOH_ABS
status: current
created: 2026-07-08
updated: 2026-07-08
pointers:
  - .agents/guides/task-persistence-investigation/2026-06-12_cophenetic-consolidation-arc.md
  - scripts/01_compute/audit/audit_152_consolidation_arc_rhosym.py
  - src/lrg_eegfc/utils/fc/coherence/windowed.py
---

# Windowed consolidation flow — the encoding/inference plane as a dynamical trajectory

**Head.** The four-phase consolidation arc (audit_152: encoding echo `T_learn`,
inference-specific persistence `T_infspec_pe`) is currently four static points.
This scope turns it into a **dynamical trajectory**: slide a 30 s window across the
whole recording, compute the LRG cophenetic state per window, and project each
window onto a **fixed, interpretable plane** whose axes are the two things the task
did — an *encoding axis* (`D_taskLearn − D_restPre`) and an *inference-specific
axis* (`D_taskTest − D_taskLearn`, orthogonalized to encoding). Watch the state
flow: baseline hovers at the origin, the task drives it out along encoding then up
into inference, and after the task it either relaxes back (reset) or stays
displaced (trace). The claim to be *visualized* (not re-tested) is that **only β
fails to relax on the inference axis** — its loop stays open — while α/δ slide back
to the encoding axis and the remaining bands barely leave the origin. This is a
faithful re-projection of an already matched-strength-verified effect, not new
evidence.

This is a **visualization measure**, not a new hypothesis test. The verdict of
record stays audit_152 (matched-strength Wilcoxon on the ρ_sym functionals). No new
null is introduced; the matched-strength surrogate landing (already computed in
`arc_null_per_patient.csv`) is shown as the reference blob.

## 5-point critical preamble

1. **Claim (visual).** The offline persistence of task reorganization is
   band-dissociated: encoding echo is broadband (δ/α/β), inference-specific
   persistence is β-only. Restated dynamically: every band's windowed state travels
   out along encoding and up along inference *during* the task; at rest_post only
   β's inference coordinate stays elevated above its rest_pre baseline scatter.
2. **Null.** The apparent "flow" is a projection artefact — any dense set of noisy
   window states, projected onto axes built from the same data, will drift
   plausibly. Equivalently: the rest_post displacement is non-specific (global
   strength/arousal drift), not aligned with the inference axis.
3. **Strongest alternative the null must cover.** (a) A **data-driven** embedding
   (PCA/UMAP of the ~7000-dim cophenetic window vectors) would be dominated by
   non-specific temporal drift, and a ρ≈0.1–0.2 trace would be invisible or, worse,
   confounded with drift that *looks* like structure. (b) Even with fixed axes, the
   rest_post cloud could sit displaced from rest_pre purely by matched-strength.
4. **Does the design cover it?** (a) is defeated by construction: the axes are
   **not** data-driven — they are fixed to the encoding/inference task contrasts, so
   the plane cannot manufacture structure the task contrasts do not contain; the
   rest_pre window cloud is drawn as the empirical no-displacement reference. (b) is
   **not** re-tested here — it was already settled by audit_152's matched-strength
   Wilcoxon (β `T_infspec_pe` p=0.0098). The matched-strength surrogate landing is
   overlaid so the reader sees the observed rest_post cloud sitting *outside* the
   null blob. **What this figure cannot do:** it cannot itself establish
   significance (overlapping windows are autocorrelated; within-phase scatter is not
   an independent sample). It shows the *shape* of a verified effect.
5. **Falsification (of the visual claim).** If, at window scale, β's rest_post
   inference coordinate is indistinguishable from its rest_pre baseline scatter, or
   if α/δ retain as much inference height as β, the dynamical reframing fails and we
   revert to the static forest (fig_arc_a). Cross-check: the cohort-median windowed
   endpoint must reproduce the sign and band-ordering of the audit_152 phase-level
   `T_learn` / `T_infspec_pe`.

## Notation

- `N` contacts (patient-specific); `P = N(N−1)/2` contact-pairs.
- Phases in recording order: `rest_pre (R0) → task_learn (TL) → task_test (TT) → rest_post (RP)`.
- Window `w`: a 30 s slice wholly inside one phase; `φ(w)` its phase label,
  `t(w)` its centre time.
- `Wᵦ(w) ∈ ℝ^{N×N}`: windowed `|ImCoh|` FC in band `b`
  (`windowed_imcoh_abs`; per-bin `|Im coherency|` then band-average).
- `D(w) ∈ ℝ^{P}`: condensed LRG cophenetic distance vector of `Wᵦ(w)`
  (`lrg_ultrametric_condensed`, τ = 1/λ_max — the locked choice).
- `D̄_φ = mean_{w: φ(w)=φ} D(w)`: window-scale phase mean state.

Everything is defined **at window scale** (window means, not full-phase FC) so the
displacement subtractions are scale-consistent; Spearman ranks make absolute scale
irrelevant regardless.

## Axes (the fixed plane), per patient per band

    e  = D̄_TL − D̄_R0                      (encoding direction)
    f₀ = D̄_TT − D̄_TL                      (inference-specific contrast)
    Δ(w) = D(w) − D̄_R0                     (displacement of window w from baseline)

Projection of any window onto the plane (rank world, matching audit_152):

    x(w) = ρ_S( Δ(w), e )                   encoding-axis coordinate
    y(w) = ρ_S( Δ(w), f₀ | e )              inference-specific coordinate (partial)

where `ρ_S(a,b|c)` is the first-order partial Spearman (the exact `_partial_rho`
audit_152 uses). `y` conditions on `e`, so it measures inference-specific alignment
*over and above* encoding — the same construction as `T_infspec_pe`.

## Properties (sanity anchors the compute asserts)

- **P1 — baseline at origin.** For `w ∈ R0`, `Δ(w) ≈ noise ⊥ e` ⇒ `x,y ≈ 0`; the
  rest_pre cloud defines the empirical origin scatter.
- **P2 — encoding pole.** `mean_{w∈TL} Δ(w) = e` ⇒ `mean x(TL) ≈ ρ_S(e,e)=1`,
  `mean y(TL) ≈ 0`. task_learn sits at the +x pole by construction.
- **P3 — inference lift during task.** `mean_{w∈TT} Δ(w) ≈ e+f₀`; since it contains
  `f₀`, `mean y(TT) > 0` for **every** band (all bands do the inference online).
- **P4 — endpoint = the trace.** `mean_{w∈RP} (x,y) ≈ (T_learn, T_infspec_pe)` at
  window scale — the audit_152 numbers. This is the faithfulness cross-check
  (sign + band-ordering must match; exact equality is not expected — windowed &
  unsplit vs full-phase ρ_sym).
- **P5 — full-window limit.** By audit_124, one window = whole phase reproduces the
  cached phase FC to ρ≈1, so the construction degenerates to the static arc.

## Caveats (stated, not narrated away)

- The plane is **hypothesis-shaped**: axes are built from the task contrasts. It
  can only redistribute what those contrasts contain; it cannot certify
  significance. The gate remains audit_152.
- Windows **overlap** (25%) ⇒ autocorrelated ⇒ the trajectory is smooth but points
  are not independent; no per-window inferential statistic is computed.
- Window FC uses ~29 Welch segments (30 s, nperseg = `nperseg_for_fs(fs)`); delta
  (0.53 Hz low edge) is the binding band — 30 s = 16 delta cycles, the minimum for a
  stable band-|ImCoh|. Shorter windows would invalidate delta; longer windows blur
  the relaxation timescale. W is a documented parameter; a W-sweep is a sensitivity,
  not a result.
- Relaxation timescale (rest_post `y(t)` decay) is **descriptive/exploratory** —
  reported as a curve against the rest_pre baseline band, not a fitted τ with a
  p-value (n windows small, autocorrelated).

## Pseudocode

    for pat in COHORT:
        # window states (heavy step; reuse windowed.py, all 6 bands per FFT cube)
        for phase in (R0, TL, TT, RP):
            X = load_timeseries(pat, phase); mask channels as phase FC does
            for w in slide(X, W=30s, step=0.75W):        # window inside phase
                freqs, ff = segment_ffts(X[:,w], fs, nperseg_for_fs(fs), fmax=305)
                cube = imcoh_abs_cube(ff)                # one (N,N,F) cube / window
                for b in BANDS:
                    Wb   = band_abs_average(cube, freqs, *BRAIN_BANDS[b])
                    D[b].append( lrg_ultrametric_condensed(Wb) )   # cophenetic
                record (pat, b, phase, w_idx, t_center)
        for b in BANDS:
            R0m, TLm, TTm = mean D over each phase's windows
            e  = TLm - R0m;  f0 = TTm - TLm
            for every window w:  x = rho(D(w)-R0m, e);  y = partial_rho(D(w)-R0m, f0, e)
            assert sign/order( mean_{RP}(x,y) ) matches audit_152 (T_learn, T_infspec_pe)  # P4
    write windowed_flow_per_window.csv  (pat, band, phase, w_idx, t_center, x, y)

Cost: ~1400 windows × (one CSD cube + 6 eigh(N)); time 1 patient, extrapolate,
print live `[i/N] ETA`, launch. Estimate ~10–20 min single-thread; parallel over
patients if needed.

## Visualization (fig_arc_h — 2 exemplars + cohort overlay)

- **Exemplar panels** (Pat_08 = strong β-tracer; Pat_10 = canonical reset): the
  windowed flow in the plane. β bold + time-coloured (viridis along `t`); α and a
  null band (θ or high_γ) thin for contrast. Phase clouds shaded; arrows mark the
  mean pre→learn→test→post path. Pat_08's β loop stays open on +y; Pat_10's
  collapses to the origin.
- **Cohort overlay**: β rest_post endpoints for all 10 patients (arrows from
  origin) + matched-strength null blob (from `arc_null_per_patient.csv`); an inset
  fan of across-band median endpoints (α/β/δ stretch along +x, only β lifts +y) —
  the disentangling summary.
- Honesty overlays: rest_pre window cloud (empirical origin scatter) and the
  matched-strength null landing, so "displaced beyond noise/null" is shown, not
  asserted. use_lrg_style; PDF-only vector; no near-white cmaps.

## Connection to prior tools

- **audit_152 / consolidation_arc_rhosym** — endpoint is its `(T_learn,
  T_infspec_pe)`; verdict of record; supplies the null landing.
- **windowed.py (audit_124 replay backend)** — the window FC engine; continuity to
  phase FC already verified there.
- **fig_arc_a** — the static forest this dynamizes; stays as the quantitative
  companion.

## Outcome (2026-07-08 run — audit_162, full cohort)

Ran R=all, W=30 s, 3 projection variants (A window-mean raw / B window-mean
z-scored / C phase-anchor z-scored). **Key finding: the windowed view is per-patient
FAITHFUL but does NOT isolate β across bands in the raw cohort median.**

- Faithfulness (variant B): `corr(x_win, T_learn)` = .73–.82 (α/β/δ),
  `corr(y_win, T_infspec_pe)` = **.68 for β**. So the dynamical flow tracks the
  verified per-patient trace — it is real, not a projection artefact.
- BUT the cohort-median inference coordinate does NOT separate β: β y≈.10 sits
  mid-pack among δ (.10), low-γ (.13), high-γ (.12). **Window noise inflates the
  non-β bands' inference coordinate and washes out the β-isolation that is clean at
  phase scale** (phase β `T_infspec_pe`=.091 vs next δ=.036). Encoding stays broadband
  at window scale (all bands x≈.16–.25 ✓). So the windowed medians show HALF the
  disentangling (encoding broadband) but not the β-only half.
- Variant A fails faithfulness (β `corr(y,Tinf)`=.07 — glued rest_post); phase-anchor
  raw (no z) GLUES all windows (scale-dominated `D(w)−D0`); z-scoring fixes the glue.
  **Chosen: variant B** (window-mean z-scored) — best-balanced faithfulness + flow.

**Decision (user, 2026-07-08): TWO figures, each honest about its job.**
- `fig_arc_h1_gated_loop` (phase-scale, gated) — the DISENTANGLING. online panel (all
  bands rise) + offline panel (β-only inference, green gate; α/δ/β encoding rings).
  rest_post = verified ρ_sym exactly (audit_163 cross-check machine precision).
- `fig_arc_h2_windowed_flow` (windowed, variant B) — the DYNAMICS. Pat_08 retains vs
  Pat_10 relaxes-back + cohort β displacement gated by verified own-null (6/10).
  Visualization only; autocorrelated windows ⇒ no per-window significance claimed.

## Open questions

- Does the β relaxation curve show a *monotone* non-return, or a transient dip then
  re-elevation (replay-like reinstatement)? If the latter, it connects to the N4
  sustained-reinstatement fold (`replay_states_verification_2026_06_22`). (fig_arc_h2
  relaxation-timescale panel — not yet built; add if asked.)
- The unsplit single-baseline trajectory is inflated + NON-specific (every band inf
  ≈0.19) — the shared rest_pre baseline induces the very correlation ρ_sym removes.
  ALWAYS draw the trajectory split-half-consistent (moving point on OPPOSITE half from
  the axis); audit_163 does this and rest_post == audit_152 to 1e-16.
- W-sensitivity: does the β open-loop survive W ∈ {20, 30, 45} s? Sensitivity panel
  if asked.

## State-space attractor embedding (audit_164) — 2026-07-09, SUPERSEDES h1/h2

**Why.** User rejected `fig_arc_h1/h2` (2026-07-09) as flat 2-D axis panels / 1-D
curves where the huge `task_test` inference peak buries the small `rest_post` residue.
Fix = a dynamical-systems / state-space portrait: amplitude of the task excursion is a
transient; the result is the **limit set** (where `rest_post` settles).

**Embedding.** Per patient, per band: windowed condensed LRG cophenetic states
`S_ph(w) ∈ R^P` (audit_162 backend), z-scored per window (unit-variance PATTERN), stacked
`Z ∈ R^{W×P}`, PCA → top-8 scores saved (`{pat}_{band}_pca.npz`; W×8, evr, phase,
t_center, w_idx). ~16 min/cohort, 6 bands. Two readouts built from the scores (no recompute):
- **raw-PCA displacement axis** `a1 = mean(task)−mean(pre)` — TOTAL structural persistence.
  Dramatic two-lobe separation (Pat_06 d'=6.37, 0% overlap) but NOT the verified quantity.
- **task-contrast axes** `a1=encoding (learn−pre)`, `a2=inference (test−learn)⊥a1`,
  `a3=residual top-PC`. rest_pre≈origin (osc σ≈7–11), task_learn jumps +40 along encoding,
  task_test drifts +9 along inference, rest_post STAYS (tracer) or RETURNS (resetter).

**Drift control (the load-bearing validation).** rest_pre is the earliest recording;
a slow drift would carry EVERY `rest_post` outward. Discriminator = resetters, who share
the time-order: `post/task` along the task axis = **0.99 / 1.00 (tracers Pat_06/08)** vs
**0.01 / −0.05 (resetters Pat_10/13)**. Resetters snap back ⇒ the separation is task-driven,
not drift; the tracer-vs-resetter contrast IS the drift control.

**Honest division of labour (brutal-honesty gate).** The windowed per-patient separation
is a NOISY proxy for the verified trace — cohort corr(gate p, windowed sep) ≈ −0.37, with
real exceptions (Pat_07 verified-trace but windowed-returns; Pat_02 windowed-stays but
unverified). So NO windowed axis reproduces the phase-scale matched-strength ρ_sym gate.
Therefore: the two attractors **ILLUSTRATE** two clean, verified exemplars only; the cohort
statistic of record stays the verified `T_infspec_pe` gate (audit_152/163), 6/10 clearing.
Never color a cohort element by the gate on a windowed axis (it scatters — incoherent).

**Figure.** `fig_arc_attractor_3d.py` → results_section2/. Two live time-ordered windowed
trajectories (phases concatenated = "one long recording"), colored ONLY by time (phase
emerges from motion+position; plasma, rest_pre→rest_post), within-phase Hann-smoothed flow
+ faint raw-window swarm + faint ghost-basin averages; 1:1:1 cube (encoding leap large,
inference drift honestly small). Cohort strip = verified ρ_sym trace per patient vs own
null, exemplars pinned in. Fully vector.

Build: `audit_164_windowed_attractor_embedding.py`, `fig_arc_attractor_3d.py`.
