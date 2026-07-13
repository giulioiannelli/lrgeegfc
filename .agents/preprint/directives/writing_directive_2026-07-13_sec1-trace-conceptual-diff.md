---
name: writing_directive_2026-07-13_sec1-trace-conceptual-diff
kind: writing-directive
era: IMCOH_ABS_COHORT_N10
status: current
created: 2026-07-13
scope: conceptual diff for the agent rewriting results_sec_1.tex (§"A held, multiscale trace of learning and reasoning"). Story-level — tell what changed and why, keep the arc, weave in τ. Read the OVERVIEW first.
pointers:
  - .agents/preprint/directives/writing_directive_2026-07-13_sparsified-recovery-OVERVIEW.md
  - .agents/preprint/overleaf/results_sec_1.tex
  - .agents/reports/2026-07-12_mst020-recovery-arc.md
---

# §1 (the trace) — conceptual diff

## Head

Your section is the one that changes most, but the **narrative arc is intact**:
reasoning leaves a held multiscale trace; it lives in \(\beta\) and \(\alpha\); in
\(\beta\) it concentrates in orbitofrontal cortex; it rides healthy cortical
coupling; it is held rather than replayed; and over the seizure zone \(\beta\) and
\(\alpha\) diverge. **Every one of those beats stays.** What changes is the
*machinery underneath the first three paragraphs*: the graph is now sparsified, and
\(\tau\) is swept rather than fixed — and that upgrade lets you make the
band-specificity claim far more strongly, because you can now show a whole **ladder
of simpler methods failing** where the multiscale hierarchy succeeds. Your job is to
re-tell paragraphs 1–3 and the "multiscale geometry" paragraph on the new footing,
and to convert the "spread is not noise" paragraph from a *control* into an
*explanation*. Don't touch the OFC / tissue / reinstatement / SOZ-divergence beats
beyond swapping "dense graph at \(\tau=1/\lambda_{\max}\)" language where it appears.

## The issue (what broke, and what the fix buys)

The trace was read from a **dense, unthresholded** connectivity graph at a **single**
diffusion time \(\tau=1/\lambda_{\max}\). We later found that on the fully-connected
graph the propagator is **degenerate at that scale** — the heat-kernel geometry is
essentially identical to the raw connectivity, so the "diffusion hierarchy" was not
yet earning its keep (audit_174; the fully-connected weight-heterogeneous graph is a
Villegas outlier where a single characteristic scale dominates). The fix is a
**minimal sparsification** to a backbone (maximum spanning tree + strongest 20 % of
edges, `mst@0.20`). On that backbone the propagator is non-degenerate, \(\tau\)
becomes a real **multiscale scanner**, and — crucially — **the previously published
gate reproduces exactly at the fine scale** and then *sharpens and differentiates by
band* as \(\tau\) grows. So this is a strengthening, not a retraction.

## What STAYS (do not weaken these)

- \(\beta\) and \(\alpha\) are the two cohort trace bands; \(\theta\) is a clean
  negative; \(\delta\)/low-\(\gamma\)/high-\(\gamma\) do not clear the cohort gate.
- The **β→OFC** two-sided localization (OFC enriched, prefrontal/sensorimotor
  depleted; \(q=0.010\)–\(0.015\); LOO-robust in 3/4; OFC sampled 5/10, expressed
  4/5). Unchanged.
- **Healthy-cortex carrier** (gray–gray/gray–white carry β, not white–white or SOZ;
  decimation control). Unchanged.
- **Held-not-replayed** reinstatement (10/10, \(p=0.001\); no bursts/sequences).
  Unchanged.
- **SOZ divergence** (β steers around the SOZ, α recruits it, \(p=0.005\)). Unchanged.
- The **Grassmann table** (Table `tab:probe_band`): keep it. It is a complementary
  spectral read on the same operator, matched-strength-verified on the dense scheme
  (audit_66); it is **not** re-derived on `mst@0.20`. Do not re-attribute it to the
  sparsified pipeline — if anything, footnote that it reads the *global modes* at a
  scale while the cophenetic distance reads *across* scales.

## What CHANGES (re-tell these)

**(A) Paragraph 1–2 framing: "dense graph, \(\tau=1/\lambda_{\max}\)" → "sparsified
backbone, \(\tau\) swept across scales."** Introduce the backbone in one clause (MST +
strongest 20 %, needed because the fully-connected propagator is degenerate at the
fine scale) and make \(\tau\) an *axis*: the trace is read as a function of scale. The
two carrier verdicts become **scale-resolved**:
- \(\beta\): a **scale-invariant** trace — clears the matched-strength gate at **all
  16 scales**, from \(p=0.014\) at \(\tau_{\min}\) (the old operating point) deepening
  to \(p=0.001\) at the mesoscale. Median \(\rho_{\mathrm{sym}}\approx+0.20\)–0.29.
- \(\alpha\): a **mesoscale** trace — clears **12/16** scales, peaking near \(s\approx5\)
  (\(p=0.007\)), and its \(\tau_{\min}\) value \(p=0.024\) **is exactly the published
  full-graph gate**. Say this explicitly: the old number is the fine-scale slice of
  the new curve.
- Keep the consistency-not-amplitude point (low-γ carries the largest single tracers
  yet fails) and the "\(\alpha\) and \(\beta\) part company over *where* it lands" hinge
  into the OFC paragraph.
- Add the **robustness sentence**: \(\beta\) is sparsification-invariant (holds on
  dense, every `mst@f`, percolation, TMFG); \(\alpha\) requires a backbone that keeps
  the *strongest* edges (`mst@0.20` keeps it; percolation at the same density kills
  it). The trace is a property of the network, not of the thresholding choice.

**(B) The "spread is not measurement-reliability" paragraph → an EXPLANATION of the
spread.** It currently only rules out reliability (\(R^2=0.06\)). Keep that, but now
*explain* the between-patient spread: the \(\beta\) trace is **left-lateralized** —
per-patient \(\beta\) strength tracks the left-hemisphere contact fraction (Spearman
\(\rho=+0.685\), \(p=0.029\)); left-dominant implants trace, the right-only implant
(Pat_15) runs null. \(\alpha\) shows **no** laterality (\(\rho=0.16\), ns), matching its
"trace without an address" status. This converts a defensive control into a positive
finding: *who traces is set by which cortex the electrodes sample* — the heterogeneity
is structure, not noise.

**(C) The "band-specificity is a property of the multiscale geometry" paragraph →
promote to a full CONTROLS LADDER.** This is the section's biggest upgrade and the
cleanest evidence for the whole thesis. Currently it contrasts raw edges vs
cophenetic. Now walk a ladder of read-outs on the *same* sparsified graph and show the
band-selective trace is **unique to the multiscale cophenetic distance**:
- **raw FC edges** register the reorganization but **non-selectively** — they fire in
  four bands (\(\delta\)/\(\alpha\)/\(\beta\)/low-\(\gamma\)), unable to separate a band
  that traces from one that does not (keep the existing "θ ≈ +0.12, indistinguishable"
  point — it now sits inside this ladder).
- **node strength**: degenerate — the strength-matched surrogate preserves it exactly,
  so \(p=1\) by construction.
- **clustering coefficient**: catches \(\alpha\) only, misses \(\beta\).
- **graph geodesic distance**: fires the **wrong** bands (\(\delta\)/low-\(\gamma\)),
  misses both carriers.
- **spectral resistance distance** (effective resistance = Laplacian pseudo-inverse):
  the strongest foil — a *global, multi-step* Laplacian read-out that is nonetheless
  **multiscale-blind**, and it finds **nothing in any band**. This is the key sentence
  for the "not even a spectral method can do it" claim: the trace needs the *nested,
  scale-resolved* hierarchy specifically, not merely a Laplacian.
- **cophenetic distance** (ours): the **only** read-out that is band-selective —
  \(\alpha\)/\(\beta\) and nothing else, sharpening from the fine scale (β \(p=0.014\),
  α \(p=0.019\)) to the mesoscale (β \(p=0.001\), α \(p=0.007\)).
Land the paragraph on the mechanism: **band-specificity lives in how the trace depends
on scale (\(\tau\)); it is invisible to any single-scale, edge-level, or spectral
summary and emerges only when the connectivity is read as a nested multiscale
geometry.** That *is* the section's thesis, now demonstrated rather than asserted.

## τ — how to weave it (this section owns the concept)

This is where \(\tau\) is introduced for the whole paper, so define it here as the
**diffusion depth / scale knob** and let the later sections lean on it. The clean
arc: at \(\tau_{\min}\) you recover the published gate; sweeping \(\tau\) up resolves
the bands by scale signature (β flat, α mesoscale, non-tracers silent at every \(\tau\));
and the reason the hierarchy beats the ladder of simpler methods is precisely that it
carries this scale dependence. Keep it physical and brief — one or two sentences of
"let heat diffuse for a time \(\tau\); short \(\tau\) sees fine splits, long \(\tau\)
sees coarse communities; the trace's band-identity is written in that dependence."

## Numbers you may cite (fresh, mst@0.20)

- Trace gate (matched-strength, per-scale): β 16/16, \(p_{\tau_{\min}}=0.014\),
  \(p_{\mathrm{meso}}=0.001\); α 12/16, \(p_{\tau_{\min}}=0.024\),
  \(p_{\mathrm{meso}}=0.007\); θ 0/16, low-γ 0/16; δ/high-γ coarse-only artifacts (not
  claimed).
- Controls ladder (T\_test gate \(p\)): see the OVERVIEW table — raw FC
  δ.024/α.019/β.024/low-γ.032; strength all 1.0; clustering α.019; geodesic
  δ.042/low-γ.014; resistance all ≥.12 (dead); cophenetic-meso α.007/β.001.
- Laterality: β trace vs left-contact-fraction Spearman \(\rho=+0.685\), \(p=0.029\);
  α \(\rho=0.16\) ns.
- Robustness: β best-scale \(\rho\) ≈ 0.20–0.30 across dense/mst-any/perc/TMFG; α best
  \(\rho\) ≈ 0.20 on `mst@0.20` vs ≈ 0.08 on percolation at the same density.
- Multiscale onset (if you mention thermodynamics): single specific-heat peak in 100 %
  of cells when dense; a second peak in ~38 % of cells at `mst@0.20`; ~0 above density
  0.5. **State as a fraction, not "the curve has N peaks."**

## Do-not

- Do not claim the C(τ) curve is multi-peaked in general (it is single-peaked in the
  median; the ladder is a minority-of-cells phenomenon — lean on the τ-signature).
- Do not re-attribute the Grassmann result to the sparsified backbone.
- Do not report best-scale p-values as the verdict; always fine + mesoscale + N/16.
