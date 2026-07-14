---
name: multiscale-narrative-audit-2026-07-13
kind: talk-directive
era: IMCOH_ABS_COHORT_N10 (mst@0.20 recovery)
status: current
created: 2026-07-13
scope: THE narrative audit for the talk "The multiscale nature of inference". The 8 intended talk points plus the offline-abstraction climax (the R2 cognitive interpretation), each graded against the settled mst@0.20 evidence — what to claim with numbers, what NOT to claim, and how to reframe so the method's necessity stays front-and-center. Hand this to the presentation orchestration agent as the single source; every number is fresh from this session's mst@0.20 runs. Do NOT let the agent reintroduce localization, "β-only inference", "encoding single-scale", or "only multiscale detects it".
pointers:
  - .agents/preprint/established_results/2026-07-13_settled-three-results.md   # the three results, locked
  - .agents/reports/2026-07-13_CHECKPOINT-localization-scale.md
  - data/sparsified_arc/figures/fig_scale_signatures.pdf                        # talk fig 1
  - data/sparsified_arc/figures/fig_method_ladder.pdf                           # talk fig 2
---

# The multiscale nature of inference — narrative audit (2026-07-13)

## Spine (the one line the whole talk defends)

**The inferred order persists as a mesoscale-emergent, distributed structure that
single-scale reads can *detect* but cannot *characterize* — and the same diffusion
operator frames both cognition and epilepsy.** The title holds if "multiscale nature"
means **mesoscale-emergence + scale-characterization**, NOT "only multiscale can see it".

Verdict grades below: 🟢 defend as-is · 🟡 defend with precise wording · 🔴 cannot claim / reframe.

---

## ★★ THE fundamental only-multiscale headline (hardened 2026-07-13) 🟢

**"Reasoning writes a *scale-invariant* trace in β."** The β trace is significant at all
16 diffusion scales with **no significant scale-dependence** (Friedman p = 0.17) and **no
monotonic trend** (p = 0.56) — statistically consistent with scale-invariance / no
characteristic scale. **α, by contrast, is significantly scale-dependent** (Friedman
p = 0.004, only 12/16 scales) and β encoding *strengthens* toward coarse scales (trend
p = 0.010). So β is *uniquely* the scale-free band.

**Why this is the one:** scale-invariance is a property that **cannot be stated off a
scale axis** — a pairwise or single-scale method literally cannot express "present at
every scale with no characteristic scale." It is robust (Friedman + trend both ns, 16/16),
and fundamental (a scale-free, self-similar memory of the inferred structure). This is the
cognitive feature only multiscale can frame; build the central slide on it.

**Honesty guardrail (say exactly this):** "**no characteristic scale / consistent with
scale-invariance**", NOT "proven flat" — with n = 10 Friedman is fail-to-reject (limited
power). That wording is bulletproof.

### The two NON-CIRCULAR wins (added 2026-07-13 PM — these answer "is scale-invariance just self-referential?")

**(A) Band-selectivity — the discrimination win.** The read-out × band table (trace,
matched-strength) shows the cophenetic hierarchy is the **only** read-out that fires in
*exactly* the cognitive bands and rejects the rest:

| read-out | fires in | verdict |
|---|---|---|
| **cophenetic (meso / τmin)** | **α (.007/.019), β (.001/.014) only** | **selective — cognitive bands only** |
| raw pairwise FC | δ, α, β, low-γ (4/6) | non-selective — conflates cognition + disease |
| clustering | α only (misses β) | incomplete |
| geodesic | δ, low-γ (WRONG bands) | anti-selective |
| resistance / strength | none | dead |

Raw FC **passes** the matched-strength null (it robustly traces) — it is NOT fragile to the
null. What it cannot do is **discriminate**: it lights up δ and low-γ, which are the *disease*
bands (R3). Only the hierarchy isolates α/β. **Say "raw FC cannot separate the cognitive trace
from the disease/global structure", NOT "raw FC fails a control".** Non-circular: it
cross-validates against R3's band dissociation. (`14_controls_ladder`, `22_raw_dense_selectivity`)

**(B) Coarse higher-order β component.** Residualizing each phase's cophenetic geometry on
its raw edges (coph⊥raw = the part of the hierarchy NOT in the pairwise edges), the β trace
**still persists** and *strengthens toward coarse scales*: coph⊥raw β is ns at τ_min (.16) →
p=.032 at meso → p=.032 at coarse, effect +0.15→+0.24 (`21_coph_beyond_raw`). So the β
reorganization carries persistent **higher-order (multi-path) structure the edges do not
contain**, and it is a *coarse-scale* phenomenon. Honest caveats: BH-marginal for β (q=.064;
the clean BH survivor is δ at coarse, q=.039), and **raw FC is complementary, NOT subsumed**
(raw⊥coph traces in all bands) — the hierarchy *adds* a coarse higher-order view, it does not
defeat pairwise FC. Non-circular (it is about information beyond pairwise, not "we have a
scale axis"). Serves point 4.

**What did NOT harden — do not claim it.** "β encoding is cophenetic-only / no single-scale
method sees it" is dead: raw FC sits on the boundary (gate p = 0.053, unchanged at R=1000)
and clears in 4/10 leave-one-out drops. The honest, still-pro-multiscale fallback is a
**robustness** claim: cophenetic detects β encoding in **10/10** LOO subsets (p = 0.003);
the best single-scale competitor is a boundary near-miss that flips with one patient. Use
"robust reader", never "exclusive".

---

## 1 — Title: "The multiscale nature of inference" 🟡

**Claim.** Transitive inference leaves a trace that is uniquely multiscale-framed.

**Defend with.** Inference-specific persistence is **mesoscale-emergent**: it does NOT
fire at the finest diffusion scale in any band (cohort gate p at s₁: β 0.097, α 0.278,
δ 0.188 — all ns) and clears only at the mesoscale (β 0.005, α 0.005, δ 0.001).
Encoding, by contrast, fires already at s₁ (β 0.032*, α 0.042*). ⇒ the *inferred* order
is literally invisible at the local scale; you must climb the hierarchy to see it.

**Cannot claim.** "Only the hierarchy detects inference." Raw pairwise FC detects β
inference (gate p = 0.010). The title is about *where on the scale axis* inference lives
(mesoscale) and its *characterization*, not exclusive detection.

**Narration.** "Inference doesn't live at the level of pairs — it's a mesoscale
reorganization of the relational hierarchy. That is the multiscale nature: the inferred
order is only legible after coarse-graining."

---

## 2 — Three-fold multiscale (time · space · topology) 🟡 (reframe the "space" leg)

**Claim.** Multiscale in time, space, topology; topology proxies space (mm electrodes).

**Defend with.** Two genuinely swept axes: **frequency bands** (δ→γ, clean band
dissociation) and **diffusion scale** (16 scales, s = τλ_max).

**Cannot claim / caveat.** (a) In LRG, **diffusion-time and topology are the SAME axis** —
the heat kernel sets s = τλ_max, so sweeping τ *is* sweeping topological scale. That is
elegance, not three independent axes. (b) The **space** leg is undercut: localization is
dead (see below), so the trace cannot be placed in physical space.

**Narration.** Present the three as **frequency × diffusion-time(≡topology) ×
consolidation timeline (pre→task→post)**, and meet the space mandate through topology
*honestly*: "because contacts are mm-scale, graph organization proxies spatial
organization — and it shows the trace is a **distributed, whole-network object**, not a
focal one." Space enters as a scale statement (network-scale), which is true, rather than
a localization claim, which is not.

---

## 3 — The trace hides properties blind to single-scale methods 🟢 (strongest point)

**Defend with three hard results.**
- **Depth grading:** β trace scale-invariant (16/16 scales, p 0.014→0.001) vs α single
  mesoscale peak (12/16). A single-scale probe sees both at the mesoscale and cannot tell
  they differ in depth.
- **Mesoscale-emergence** of inference (point 1) — invisible at the fine scale.
- **Method ladder (fig 2):** the cophenetic hierarchy detects β encoding **robustly**
  (p = 0.003, 10/10 LOO) while clustering (0.216), geodesic (0.116), and spectral
  resistance (0.138) fail robustly; raw FC is a boundary near-miss (0.053, LOO-fragile 4/10)
  — say "robust reader", NOT "only the hierarchy sees it" (see ★★ headline caveat).

Lead the methods-justification with the **scale-invariance** result (★★) — it is the clean,
robust, only-multiscale claim; the ladder is supporting robustness, not exclusivity.

---

## 4 — Multiscale "is" a higher-order feature 🟡 (graph higher-order, NOT statistical)

**Defend with.** The propagator e^{−τL} sums walks of **all lengths**, so cophenetic
affinity reflects *all paths* between two contacts, not their single shared edge. |ImCoh|
is a **second-order (pairwise)** statistic; the LRG step builds **higher-order graph
features** on top of it. Evidence the working feature is beyond pairwise: the whole ladder
of pairwise/local/spectral scalars fails where cophenetic works (point 3 / fig 2).

**Must not say.** "higher-order interactions" / "nonlinear". It is higher-order in the
**path/graph** sense only. Use "multi-step / all-paths graph structure".

---

## 5 — Multiscale enriches both spectral-clustering and raw-pairwise views 🟢

**Defend with the method ladder (fig 2), a per-read-out table for β:**
- **Raw pairwise FC:** sees β trace (0.024) and β inference (0.010) — but misses β
  encoding (0.053) and has no scale axis to grade depth.
- **Spectral clustering / effective resistance:** dead in all bands (0.116 trace, 0.138
  encoding, 0.812 inference).
- **Clustering coefficient:** catches α encoding (0.010) but misses β encoding (0.216).
- **Cophenetic hierarchy:** catches all of the above *and* supplies the scale-signature.

**Narration.** "Each single-scale view catches a fragment; the hierarchy catches the union
and adds the scale axis none of them have." Frame as **enrichment** (adds an axis), not
replacement — honest, since multiscale does not dominate raw FC on detection, only on
characterization.

---

## 6 — ρ^coph_sym as the measure; two orthogonal axes 🟢 (methodological clarity)

**Defend with definition (no null needed).** ρ^coph_sym = Spearman on **cophenetic**
distance vectors (relational, not positional), symmetrized over the rest_pre split-halves
(removes the arbitrary-half artifact). **One UPGMA hierarchy per scale τ.** Keep the two
axes un-confounded:
- **Scale axis (τ):** *which* hierarchy — the curve ρ_sym(τ) (fig 1).
- **Cross-phase axis:** *how similar* two phases' hierarchies are at that τ.
They are orthogonal (ρ_sym is computable at every τ). This is the "don't confound the two
hierarchies" slide, and it is what licenses "one number per (band, scale)".

---

## 7 — Encoding / inference dissociation 🟡 (a SCALE dissociation; the old ones died)

**Defend with.**
- **Scale-shape:** encoding fine-legible + β scale-invariant (fires s₁, all 16 scales);
  inference mesoscale-emergent (ns at s₁, emerges meso). Holds for α *and* β.
- **Construction:** inference-specific **survives partialling out encoding**
  (T_infspec_pe) and still persists — separable components.

**Cannot claim (all died on the backbone).** encoding→OFC / inference→cingulate
(delocalized); "inference β-only" (it is δ/α/β); and "encoding single-scale, inference
multiscale" (β encoding is the *most* scale-broad thing, 16/16). The **band** leg is also
weak now (α and β carry both).

**Narration.** "Premises and inferred order dissociate by the **scale** at which they
persist: premises are encoded locally and at all scales; the inferred order is a
mesoscale-emergent structure." This is the honest dissociation and it directly serves the
title.

---

## ★ Offline abstraction — rest keeps the order *inferred*, not just the pairs *shown* 🟡 (the R2 climax / cognitive interpretation)

**Claim.** What persists into post-task rest is not merely the premise pairs the subject
was *shown* (encoding), but the full transitive order the subject had to *infer* — an
**offline abstraction** of the latent structure, held during rest. This is the cognitive
reading of the encoding/inference dissociation (point 7) and the payoff of the title.

**Defend with.**
- **A persistent component beyond the shown pairs.** Inference-specific persistence
  (T_infspec_pe) survives partialling out encoding and still clears: β 7/16
  (p_meso 0.005), α 7/16 (0.005). So rest holds something the premises alone do not
  explain — the inferred order.
- **The abstraction lives at a coarser scale than the surface pairs.** Inference-specific
  is mesoscale-emergent (ns at s₁: β 0.097, α 0.278; emerges at the mesoscale), exactly
  what an *abstraction* of relational structure should look like — coarser than the
  literal, fine-scale premises.
- **Not a recording-length artifact.** Duration control ns (β inference vs test/learn
  ratio: s₁ ρ +0.03 p 0.93; peak ρ +0.14 p 0.70).
- **Held, not replayed.** Sustained-reinstatement companion: rest_post dwells in the task
  configuration *more tightly* than rest_pre (10/10, p 0.001, LOO 0.002) — an offline
  *state*, not a bursty replayed sequence.

**Cannot claim.**
- **No cortical address.** The old "premises in OFC / inferred order in cingulate"
  dissociation is delocalized on the backbone — the abstraction is **distributed**, not
  placed. Frame the *where* as network-scale, per the localization reframe below.
- **Not β-only** (δ/α/β), and δ inference-specific is CUT (partial-corr artifact).
- **"Held not replayed" is a whole-graph temporal companion** (τ_min, dense windows), not
  a mst@0.20/scale result — present it as corroboration, not a backbone finding.

**Narration.** "After the task, rest doesn't just echo the pairs the subject saw — it
holds the order the subject *worked out*. The inferred transitive structure persists as a
sustained, offline abstraction: coarser than the surface pairs (mesoscale), distributed
across the network, and held as a state rather than replayed as a sequence." This is the
sentence the title is built to earn.

---

## 8 — Epilepsy and cognition in one multiscale framework 🟢 (strong closer)

**Defend with.**
- **Same operator** e^{−τL}: cross-phase trace (cognition) and seed-affinity community
  (epilepsy) from one propagator.
- **Band dissociation:** cognition α/β, epilepsy δ/low-γ/β, **β the bridge** — mutually
  corroborated: β *spares* the SOZ (tissue epi_epi clears only 1/16 scales) while
  **low-γ concentrates in the SOZ** (the one anatomical signal that survives correction,
  q = 0.010).
- **Numbers:** seed-affinity AUC δ 0.83 / low-γ 0.82 / β 0.745; beats matched-strength
  fake-SOZ null 7–8/10; robust to sparsification (dense→mst@0.20 reproduces).
- **Multiscale property of the marker:** τ lifts **ranking** (multi>single AUC 7–8/10) but
  not per-contact precision.

**Narration.** "The same coarse-graining that reveals the inferred order reads the
epileptogenic network — cognition and disease separate by band inside one operator."

---

## The elephant: localization is dead — turn it into a result

On mst@0.20, **nothing localizes** — trace, encoding, and inference are all delocalized
(0/16 BH-clearing cells, every band); the β trace is **tissue-distributed** (gray q 0.015,
white q 0.015, cross q 0.015; 14–16/16 scales). This kills "the hierarchy tells you
*where*". Do NOT hide it — state it:

> "Reasoning rewrites the network at the **whole-network scale** — a distributed
> reorganization (~88% of nodes co-move), not a hotspot."

Consistent across three independent views (localization null + tissue-distributed +
rank-size co-movement). It *replaces* the anatomy leg cleanly and is itself a
spatial-scale statement (points 2, 3).

---

## The two talk figures (built, cached)

- **fig_scale_signatures.pdf** — cohort median ρ_sym vs scale for trace/encoding/inference,
  bands β/α/δ; filled marker = cohort-significant; vertical line at s=1 (local). Carries
  points 1, 3, 7 (encoding fires at the local scale; inference silent locally → mesoscale).
- **fig_method_ladder.pdf** — β gate p per read-out (cophenetic meso/τmin, geodesic,
  clustering, raw FC, resistance, strength) × component (trace/encoding/inference), red box
  = p<.05. Carries points 4, 5 (cophenetic catches β encoding where raw/spectral/clustering
  fail).

---

## Locked numbers the agent may cite (fresh, mst@0.20)

| item | number |
|---|---|
| β trace | 16/16 scales, p 0.014→0.001 |
| **β scale-invariance** | **Friedman p 0.17 (ns), trend p 0.56 (ns), CV 0.25 → no characteristic scale** |
| **α scale-dependence** | **Friedman p 0.004 (sig), 12/16 scales → single-scale/mesoscale** |
| **β encoding robustness** | cophenetic 10/10 LOO (p 0.003); raw FC boundary 0.053, LOO-fragile 4/10 (R=1000-confirmed) |
| α trace | 12/16 scales, mesoscale peak |
| encoding β / α | β 16/16 (p_meso 0.003); α 3/16 (fine, fragile) |
| inference-specific | β 7/16, α 7/16, δ 10/16 (mesoscale); ns at s₁ (β .097/α .278/δ .188) |
| δ inference | CUT (partial-corr artifact — δ null for trace & encoding) |
| localization | delocalized, 0/16 BH-clearing, every target/band |
| tissue (β) | gray/white/cross all q 0.015, 14–16/16 scales; epi_epi spared (1/16) |
| laterality (β) | ρ(left-frac) >0 all 16 scales, sig 6/16, peak ρ 0.74 p 0.014; α flat |
| duration (β inf) | ns: s₁ ρ +0.03 p 0.93; peak ρ +0.14 p 0.70 |
| epilepsy AUC | δ 0.83 / low-γ 0.82 / β 0.745; beats null 7–8/10 |
| low-γ → SOZ | q 0.010 (the surviving anatomical signal) |
| null | matched-strength ONLY (drift retired); read per-scale, never best-scale |

**Do-not list for the agent:** no localization / "where it lives" for cognition; no
"β-only inference"; no "encoding single-scale, inference multiscale"; no "only multiscale
detects it"; no drift; no best-scale; no calling any band "false".
