---
name: talk-slide-15-encoding-vs-inference
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 15
status: draft
updated: 2026-07-13
canva: NEW — merges old 17 (decomposition setup) + old 18 (both persist). Dissociation is now by
  SCALE, not anatomy; the "where" (anatomical localization) is CUT — there is no localization slide (dropped 2026-07-14; the cognitive trace is delocalized).
---

# Slide 15 — Encoding vs inference: the order it never saw, held offline

1. TITLE
Held offline — not just the pairs it saw, but the order it inferred

2. MAIN CONCEPT
- **Head (the juice):** the held reorganization is not one thing. In the task the patient was *shown*
  adjacent pairs (**encoding**) and had to work out an order they were *never shown* (**inference**).
  The four-phase design separates them, and rest keeps **both**. The interesting one — the inferred
  order — is **not β-only**, and it is **mesoscale-emergent**: it is invisible at the finest scale and
  appears as the network coarse-grains, exactly the signature of an abstraction rather than a replayed
  detail.
- **The decomposition.** Four phases → three hierarchy-difference quantities: encoding
  `Δ_e = D(learn) − D(pre)` (reorganization while the pairs are on screen), inference-specific
  `Δ_f = D(test) − D(learn)` with encoding partialled out (the structure never shown), persistence
  `Δ_p = D(post) − D(pre)` (what lasts into rest). Encoding and inference are **co-equal** persisting
  results — encoding is not a nuisance we scrub away.
- **Encoding persists.** β at **every scale** (16/16, p_meso ≈ 0.003, scale-invariant); α too (3/16;
  τ_min 0.042, meso 0.032), and there **ordinary tools already register it** (raw FC 0.032, clustering
  0.010). So encoding shows the same split the whole trace does — β genuinely multiscale, α resolvable
  by simple measures.
- **Inference-specific persists — and the correction to the old story: NOT β-only.** β at 7/16 scales
  (p = 0.005) **and** α at 7/16 (p = 0.024). The order the brain reasoned out — encoding controlled away
  — survives into rest. *(The earlier δ inference-specific effect is a partial-correlation artifact — δ is
  null for both encoding and the whole trace — and is CUT.)*
- **The dissociation is by SCALE, not anatomy.** Inference is **mesoscale-emergent**: at the finest scale
  it is not significant (s₁ p ≈ 0.097) and it becomes significant only as the network coarse-grains
  (mesoscale p ≈ 0.005). That within-hierarchy climb — both scales on the *same* backbone — is clean, and
  it is the point: a re-grouping that only exists once you zoom out is the *shape of an abstraction*, not a
  memorized pair.
- **Detection vs form.** Raw FC already **detects** β inference (p = 0.010) — the cophenetic distance is
  built from the same edges, so detection is tied by construction. What the hierarchy adds is *form*: it
  shows the inferred order as a mesoscale coarse re-grouping. (*Where* it lives is the next slide.)
- **Held, not replayed; distributed; not a length artifact.** The reinstatement is a sustained state
  across post-task rest (10/10), not transient bursts; it is **distributed** — no hotspot, ~88% of nodes
  co-move; and the inference signal does **not** scale with test/learn duration (ρ = +0.10, p = 0.78).

3. ON-SLIDE TEXT

persistence: memory or inference?
encoding (e), inference (f), rest persistence (p)

  $\Delta_e(A) = D_{\mathrm{learn}}-D_{\mathrm{pre},A}$   (encoding)
  $\Delta_p(A) = D_{\mathrm{post}}-D_{\mathrm{pre},A}$   (persistence)
  $\Delta_f = D_{\mathrm{test}}-D_{\mathrm{learn}}$   (inference)     [$A,B$: rest$_\mathrm{pre}$ split-halves]

  $T_{\mathrm{enc}}=\rho_S\big(\Delta_e(A),\,\Delta_p(B)\big)$   (encoding; sym. $A\!\leftrightarrow\!B$)
  $T_{\mathrm{inf}}=\rho_S\big(\Delta_f,\,\Delta_p(B)\mid \Delta_e(A)\big)$   (inference; sym. $A\!\leftrightarrow\!B$)
  $\displaystyle \rho_S(x,y\mid z)=\frac{\rho_{xy}-\rho_{xz}\,\rho_{yz}}{\sqrt{(1-\rho_{xz}^{2})(1-\rho_{yz}^{2})}}$   (partial Spearman; $\rho_{xy}\!=\!\rho_S(x,y)$)

encoding persists — β every scale, α too
inference persists — β AND α, only at the mesoscale
mesoscale-emergent = the shape of an abstraction

4. SPEECH
The trace holds — but holding what? In the task the patient did two things: they were shown adjacent
pairs, and from those they had to work out an order they were never shown. Encoding versus inference.
The four phases let us separate them cleanly, and rest keeps both. Let me show you how the measure is
built, because it's the same idea as before, just resolved into two pieces. In every phase we read the
hierarchy as a cophenetic distance — how far apart two contacts sit in that phase's tree. Encoding is how
much the tree moved from rest to learning, while the pairs were on screen. Inference is the *further* move
from learning to test — the structure the pairs never showed. And persistence is what rest holds onto,
the shift from before the task to after. Then, across every pair of contacts, we ask by rank correlation:
did the task's reshaping predict what rest kept? For inference we partial encoding out, so what's left is
only what the inferred order added beyond the pairs it saw. And one honest detail that makes the whole
thing trustworthy: we never measure a task-shift and a rest-shift against the *same* baseline — we split
the pre-task rest into two independent halves and cross-pair them, then average. That symmetrization is
what stops a shared baseline from manufacturing a correlation; drop it and the number more than doubles.
Encoding persists strongly — in beta
at every scale, and in alpha too, where ordinary tools already see it. Now the inferred order — and here
is the correction to the older version of this talk. It is not beta alone: beta and alpha both, seven of
sixteen scales each. And it has a beautiful property — it is mesoscale-emergent. At the finest scale it
isn't there; it only becomes significant as the network coarse-grains. Think about what that means: a
re-grouping that exists only when you zoom out is the signature of an abstraction, not a replayed detail.
One honesty point: raw connectivity already detects the beta inference signal, because the hierarchy is
built from the same edges — so detection is tied. What the multiscale read adds is the form, that this
persistence lives as a coarse re-grouping. And it's held, not replayed — a sustained state across the ten
minutes of rest, spread across the network with no single hotspot, and it doesn't scale with how long the
phases ran. Rest keeps not just what the brain saw, but the order it built.

Careful: inference is α AND β, NEVER "β-only"; the δ inference-specific effect is a partial-correlation
artifact, CUT (never mention it). Do NOT say "only multiscale sees inference" / "invisible to simple
methods" — raw detects β inference at .010; the value-add is FORM (mesoscale re-grouping), and LOCALIZATION
is CUT — no localization slide, NOT detection. Do NOT name OFC / cingulate on this slide — the "where" (β coarse
left-hemisphere, α placeless, low-γ→SOZ) is DROPPED; the old "encoding→OFC / inference→cingulate different
cortex" is RETIRED (delocalized under a valid statistic). Encoding is a co-equal persisting result, not a
nuisance partialled away. Null = matched-strength only (no drift, no "second null"). Read per-scale — never
a single +0.091 / p=.0098 scalar. "inference / abstraction" is licensed by design + the mesoscale-emergence,
never a brain–behaviour correlation (no behavioural data).

5. FIGURES  (produced 2026-07-13; talk PNGs in data/outputs/figures/talk/)
- **[BUILT] Trajectory taxonomy in the enc/inf space** — `slide15_encoding_inference_gallery.png`
  (src `data/preprint/figures/new_results_sec2/dynamics_gallery_s5.6.pdf`; script
  `scripts/01_compute/sparsified_arc/17_dynamics_gallery_mst020.py`). Four 3-D portraits whose axes ARE the
  decomposition: encoding x = state·(learn−pre), inference y = state·(test−learn)⊥x, residual z. Use as the
  OPENING geometric visual — "here is the space encoding and inference define, and states move through it."
  ⚠ HONESTY: the four archetype labels (anchor/reset/reorganize/TRACE) are of the RAW windowed trajectory,
  which is ρ≈0 with the matched-strength gate; MS collapses the 3-way, only TRACE-like holding survives.
  Caption it as the SPACE + range of motion, NOT four validated cross-phase categories; the "trace" panel is
  state-space-held, not gate-cleared.
- **[BUILT] Encoding vs inference branch-origin network** — `slide15_encoding_inference_chord_Pat05_beta.png`
  (src `data/preprint/figures/new_results_sec1/_drafts/fig_hierarchy_chord_network_Pat_05_beta_rest_post_mst020_DRAFT.pdf`;
  script `.../new_results_sec1/fig_hierarchy_chord_network_mst020.py`, mst@0.20 backbone, s=5.6). Bundled
  |ImCoh| network in the centre; the held rest_post hierarchy drawn OUTSIDE, each branch tinted amber
  (encoding-origin) ↔ teal (inference-origin). Shows the co-equal claim directly: the held tree carries BOTH
  encoding- and inference-dominated branches (inf-merges 0.37). Illustrative single-patient tint, not a gated
  per-node statistic; not a cohort claim.
- **[BUILT — the scientific payoff]** Mesoscale-emergence curve — `slide15_encinf_scale_payoff.png`
  (src `data/preprint/figures/new_results_sec2/fig_encinf_scale_payoff.pdf`; script
  `scripts/01_compute/sparsified_arc/fig_encinf_scale_payoff.py`). Two panels (β, α): cohort concordance
  with persistence vs diffusion scale s, encoding `T_learn` (amber) + inference-specific `T_infspec_pe` (teal),
  over the matched-strength null floor; filled marker = gate_p<0.05 at that scale, open = ns, "ns" flags the
  finest scale. Carries the headline directly: inference-specific is ns at s₁ (β .097 / α .28) → significant
  at the mesoscale, in BOTH bands (7/16 each); encoding is scale-broad in β (16/16). Numbers on-figure.
- ⚠ RETIRE / REBUILD: any figure showing the β inference wedge glowing ALONE (contradicts α+β 7/16 each), or
  encoding→OFC / inference→cingulate "different cortex" (that localization is retired; the localization slide was dropped 2026-07-14).

6. REFERENCES
- Behrens 2018; Wilson 2014; Schuck 2016 (cognitive map / transitive inference — carried from slide 2). The
  four-phase decomposition and the mesoscale-emergence of inference-specific persistence are ours.

7. CANVA STATUS
NEW slide (merges old 17 setup + old 18 persistence). Recast so the encoding/inference dissociation is by
SCALE (inference mesoscale-emergent = offline abstraction), the "different cortex / OFC-vs-cingulate"
localization is CUT (dropped 2026-07-14 — no localization slide; cognitive trace delocalized), inference is α AND β (δ cut),
and detection-vs-form is stated. Missing on deck: decomposition schematic + mesoscale-emergence curve,
compressed on-slide text, presenter notes.
