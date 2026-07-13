---
name: writing_directive_2026-07-13_sec2-encinf-conceptual-diff
kind: writing-directive
era: IMCOH_ABS_COHORT_N10
status: current
created: 2026-07-13
scope: conceptual diff for the agent rewriting results_sec_2.tex (§"Offline abstraction of a learned structure"). Encoding HOLDS on the recovered scheme; the "inference in β alone" headline WEAKENS honestly and reframes to localization. Read the OVERVIEW first.
pointers:
  - .agents/preprint/directives/writing_directive_2026-07-13_sparsified-recovery-OVERVIEW.md
  - .agents/preprint/overleaf/results_sec_2.tex
  - .agents/reports/2026-07-12_mst020-recovery-arc.md
  - raw_vs_coph_inference_localization_2026_07_10   # audit_171: inference→cingulate is coph-only; raw FC is placeless
---

# §2 (encoding vs inference) — conceptual diff

> ⚠ **STOP — do NOT rewrite this section around "encoding is washed out / null in the
> multiscale read."** That premise is a **bug artifact** and is FALSE. The shared-baseline
> error in `T_learn` made encoding look null; with it fixed, **encoding is recovered, and
> β encoding is the single most multiscale-exclusive result in the whole section** (only
> the cophenetic hierarchy sees it; raw FC, clustering, geodesic, and the spectral read
> all miss it). If your draft says encoding is strength-explained or invisible to the
> hierarchy, discard it and use the band-split below.

## Head

Two things happened when we re-derived this section on the sparsified `mst@0.20`
backbone. **Good news, and the load-bearing half:** the *encoding* decomposition
**holds** — encoding persists in \(\alpha\) and \(\beta\), exactly as the paper says.
(An earlier internal re-run made encoding look null; that was a **script bug** — a
shared-baseline error in the encoding functional, now fixed — not biology.) **The
honest correction:** the flagship "inference-specific persists in \(\beta\) **alone**"
claim does **not** survive as stated. On the recovered scheme the inference-specific
component clears at the mesoscale in \(\delta\), \(\alpha\), **and** \(\beta\) — not
\(\beta\) alone — and, separately, the controls ladder shows the \(\beta\)
inference-specific component is **also readable from raw connectivity**, so it is not
something only the multiscale hierarchy can see. The section's real, defensible
contribution shifts from *"only \(\beta\), only the hierarchy"* to **localization**:
the hierarchy's unique value is that it places encoding in OFC and inference in the
cingulate, where raw connectivity is anatomically **placeless**. Rewrite the two
opening paragraphs around that; keep the anatomy paragraphs, which are the point.

## What STAYS (these are the section's spine now)

- **Encoding persists — but the two bands are NOT the same kind of result; split them.**
  - **β encoding = the multiscale-exclusive result.** Scale-invariant (16/16 scales,
    \(p_{\mathrm{meso}}=0.003\), \(p_{\min}=0.001\)) AND **cophenetic-unique**: on the
    controls ladder only the hierarchy clears it (coph \(p=0.003\)); raw FC is a near-miss
    (\(p=0.053\)), clustering (\(0.216\)), geodesic (\(0.116\)) and the spectral resistance
    (\(0.138\)) all fail. So the premises the brain was shown leave a trace that *only the
    multiscale read recovers* — this is a strong sentence, use it.
  - **α encoding = the "stored in a basic metric" result.** It is readable straight from
    raw FC (\(p=0.032\)) and clustering (\(p=0.010\)) — it does **not** need the hierarchy;
    and in the cophenetic read it is marginal and run-fragile (3/16 scales; \(p=0.042\) at
    \(\tau_{\min}\), \(0.032\) at the mesoscale in the enc/inf run but \(0.053\) in the
    controls-ladder run — it flips across \(p=0.05\)). Do **not** sell α encoding as
    multiscale; frame it as "α keeps what was shown, and a simple network read already
    sees it."
  - This maps to the paper's "α keeps encoding, β keeps encoding + inference" but sharpens
    it: β's encoding is the hierarchy-only one, α's is the basic-metric one.
- **Encoding anchors in OFC**, and learning sets the anchor (β; \(q=0.010\)/\(0.040\)
  with/without SOZ; same orbitofrontal hotspot as the test-phase trace). Unchanged and
  important — this is the anatomical backbone of the section.
- **Inference localizes to the cingulate**, away from the OFC encoding anchor
  (\(q=0.030\)/\(0.040\)). Unchanged — and now *elevated* to the headline (below).
- **low-γ cingulate memory trace** (region-by-region read; \(\rho_{\mathrm{sym}}=+0.25\),
  \(q=0.035\), 8/8 sampled; +0.39 without SOZ) and the **cingulate multiplexing** beat
  (memory at low-γ, inference at β in the same cortex). Unchanged.
- **Duration control**: β inference is uncorrelated with the test/learn length ratio
  (\(\rho=+0.10\), \(p=0.78\)) while the length-tracking bands (α, high-γ) carry no
  inference-specific trace. Keep it — it still does its job for the β inference
  component that survives.

## What CHANGES (rewrite the two opening paragraphs)

**(A) Drop "\(\beta\) alone" for the inference-specific component.** On the recovered
scheme, \(T_{\mathrm{inf}\cdot e}\) (inference-specific, encoding partialled out) clears
the matched-strength gate at the mesoscale in **three** bands: \(\delta\) (10/16,
\(p_{\mathrm{meso}}=0.019\)), \(\alpha\) (7/16, \(p=0.024\)), and \(\beta\) (7/16,
\(p=0.005\)). \(\beta\) is still the *strongest and cleanest* (deepest \(p\), the band
that also localizes), but it is **not exclusive**. Rewrite the claim as "the
inference-specific component persists most strongly and specifically in \(\beta\)"
rather than "in \(\beta\) alone", and drop the "every other band null, next-smallest
\(p=0.25\)" sentence — that was the old dense-scheme result.

**(B) Move the value-add from *detection* to *localization*.** The controls ladder
(§sec 1 / OVERVIEW) shows the \(\beta\) inference-specific component is **readable from
raw FC edges too** (\(p=0.010\)), the same order as the cophenetic read (\(p=0.010\)).
So the hierarchy is **not** what makes inference *detectable*. What the hierarchy —
and only the hierarchy — supplies is **where** it lives: read through the cophenetic
geometry the inference component **localizes to the cingulate**, whereas the same
component read through raw connectivity is **anatomically placeless** (audit_171: raw
FC gives no surviving system). Make that the paragraph's punchline: *the multiscale
read-out earns its place here by localizing a component that simpler measures can
detect but cannot place.* This is more honest and, frankly, a better story — it ties
§2 to the same "the hierarchy tells you **where**" thread as the β→OFC result in §1.

**(C) Reframe the "offline abstraction" interpretation accordingly.** The persistent
\(\beta\) component is still an offline abstraction of the inferred order, held not
replayed, anchored where reasoning-relevant cortex sits (cingulate for the
inference-specific part, OFC for the encoding part). Just stop leaning on
band-exclusivity as the evidence; lean on the **encoding/inference anatomical
dissociation** (OFC vs cingulate) plus the duration control.

## τ — how to weave it

Lighter here than in §1, but present: encoding \(\beta\) is **scale-invariant** (holds
across \(\tau\)), α encoding is **fine-and-mesoscale**, and the inference-specific
components surface at the **mesoscale**. One sentence noting that the decomposition was
read across scales (not at a single \(\tau\)) and that β encoding holds at every scale
while the inference-specific parts are mesoscale is enough. Do not best-scale.

## Numbers you may cite (fresh, mst@0.20, T_learn bug fixed)

- Encoding \(T_{\mathrm{learn}}\): β 16/16, \(p_{\tau_{\min}}=0.032\),
  \(p_{\mathrm{meso}}=0.003\); α 3/16, \(p_{\tau_{\min}}=0.042\),
  \(p_{\mathrm{meso}}=0.032\) (α-meso is fragile: 0.053 in the controls-ladder run);
  δ/θ/low-γ 0/16.
- Controls ladder, ENCODING cell (`14`, gate p): **β cophenetic-unique** — coph_τmin
  0.024 / coph_meso **0.003**; raw FC 0.053 (near-miss, fails), clustering 0.216,
  geodesic 0.116, resistance 0.138, strength 1.0. **α basic-metric** — raw FC **0.032**,
  clustering **0.010**, coph_τmin 0.032, coph_meso 0.053. ⇒ β encoding needs the
  hierarchy; α encoding does not.
- Inference-specific \(T_{\mathrm{inf}\cdot e}\): β 7/16 (\(p_{\mathrm{meso}}=0.005\)),
  α 7/16 (0.024), δ 10/16 (0.019); low-γ/θ 0/16; high-γ 1/16. **Not β-only.**
- Controls ladder, inference cell: readable from raw FC (β \(p=0.010\), α 0.024,
  low-γ 0.024) **and** cophenetic-meso (β \(p=0.010\)) → detection is not
  multiscale-exclusive; localization is (cingulate survives only in the cophenetic
  read; raw FC placeless — audit_171).
- Anatomy (unchanged): encoding→OFC \(q=0.010\)/\(0.040\); inference→cingulate
  \(q=0.030\)/\(0.040\); low-γ cingulate memory \(q=0.035\); duration β \(\rho=+0.10\),
  \(p=0.78\).

## Do-not

- Do not write "\(\beta\) alone" / "only the hierarchy detects inference" — both are
  now false on the recovered scheme.
- Do not present the encoding recovery as a *change* to the reader — the paper already
  claimed α/β encoding; the bug was internal. Just make sure the section's encoding
  numbers match §3-OVERVIEW.
- Do not best-scale; fine + mesoscale + N/16.
