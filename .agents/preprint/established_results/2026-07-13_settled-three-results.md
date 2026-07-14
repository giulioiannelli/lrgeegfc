---
name: settled-three-results-2026-07-13
kind: verdict
era: IMCOH_ABS_COHORT_N10 (mst@0.20 recovery)
status: current
created: 2026-07-13
scope: THE locked, verified, non-contradictory statement of the three major results for the preprint + talk. Every number re-verified fresh this session on the mst@0.20 backbone. Supersedes the "cophenetic-unique / multiscale-exclusive" framing (an overclaim). Single source of truth — write the presentation from this.
pointers:
  - .agents/reports/2026-07-12_mst020-recovery-arc.md
  - .agents/preprint/directives/writing_directive_2026-07-13_sparsified-recovery-OVERVIEW.md
---

# The three results — settled (2026-07-13)

## Head

Three results, one operator. (1) Reasoning leaves a **held, genuinely multiscale trace**
in β (scale-invariant) plus a single-scale one in α; β localizes to orbitofrontal cortex.
(2) What persists is not just the pairs shown (**encoding**) but the order **inferred**;
encoding and inference sit in different cortex (OFC vs cingulate). (3) The **same
diffusion operator** localizes epileptogenic tissue as a strength-independent co-diffusing
community. The bands dissociate cleanly — cognition in α/β, epilepsy in δ/low-γ/β, β the
bridge — and the whole thing runs on one sparsified diffusion propagator with **one null**
(matched-strength).

---

## 0. The unifying frame (say this once, keep it consistent everywhere)

**Pipeline.** |ImCoh| FC → **sparsify to a backbone** (MST + strongest 20%, `mst@0.20`) →
diffusion \(e^{-\tau L}\) **swept over τ** (scale \(s=\tau\lambda_{\max}\)) → cophenetic
hierarchy → cross-phase \(\rho_{\mathrm{sym}}\). **One null: matched-strength.** (Drift
retired 2026-07-12: a directional task makes a drift null degenerate with the trace.)
Read **per-scale**, never best-scale.

**The frame is information-DEPTH, not exclusive detection.** Simple pairwise FC gives a
**scalar verdict per band** — *"did the reorganization lean back toward the task?"* — and
several bands say yes. That scalar is the ceiling of its information. The multiscale read
gives a **structured object** (a hierarchy per phase, a ρ_sym(τ) curve per band), so it
*characterizes* what the scalar cannot even pose:

1. **Scale-signature — does the firing hold multiscale nature?** ρ_sym(τ) is a curve, so it
   grades each firing band: **β scale-invariant (deep)**, **α single mesoscale peak
   (single-scale)**, **δ fires-but-shallow** (fragmented, fails fine+meso, rides the coarse
   collapse). A scalar has no scale axis to grade this. **HARDENED (2026-07-13 PM): the β
   trace is statistically SCALE-INVARIANT** — no significant scale-dependence (Friedman
   p=0.17) and no monotonic trend (p=0.56), significant at all 16 scales; **α is
   significantly scale-DEPENDENT** (Friedman p=0.004, 12/16). This is THE fundamental
   only-multiscale claim (scale-invariance cannot be stated off a scale axis). Say "no
   characteristic scale / consistent with scale-invariance" (n=10 → Friedman is
   fail-to-reject).
2. **Where it lives — localization is DEAD (full sweep done 2026-07-13).** Per-system
   decomposition of the **trace, encoding, AND inference** over all 16 scales (`17_localization_arc`)
   clears BH in **0** cells, every band — **nothing localizes** (β→OFC, encoding→OFC,
   inference→cingulate all delocalized; the whole-graph anatomy was a degenerate-graph
   artifact). The trace is **tissue-distributed** (β gray/white/cross all q≈0.015, 14–16/16
   scales; `19_tissue_pairclass`). The **only** surviving anatomical signal is **low-γ→SOZ**
   (a disease signal, q=0.010). ⇒ Do NOT claim "the hierarchy tells you where" for cognition;
   frame *where* as **network-scale distribution** (~88% co-move). Lean on ① and ③.
3. **Higher-order multi-step structure (non-circular, tested).** The propagator sums walks
   of all lengths, so cophenetic affinity reflects **all** paths between two contacts, not
   their single shared edge. Tested: residualizing each phase's cophenetic geometry on its
   raw edges (coph⊥raw), the **β trace still persists and strengthens toward coarse scales**
   (coph⊥raw β: ns@τ_min → p=.032@meso/coarse, effect +0.15→+0.24; `21_coph_beyond_raw`) —
   persistent higher-order structure the edges do NOT contain. Honest: BH-marginal for β
   (q=.064; δ@coarse clean q=.039), and **raw FC is COMPLEMENTARY not subsumed** (raw⊥coph
   traces in all bands) — the hierarchy *adds* a coarse higher-order view, it does not defeat
   pairwise FC.
4. **Band-selectivity — the discrimination win (non-circular).** The cophenetic is the ONLY
   read-out that fires in *exactly* the cognitive bands: **α (.007), β (.001) only** (2/6),
   while raw FC fires in δ/α/β/low-γ (4/6, conflating cognition with the δ/low-γ **disease**
   bands of R3), geodesic fires in the WRONG bands (δ/low-γ), clustering misses β, resistance
   /strength are dead (`14_controls_ladder`). Raw FC **passes** the null (robustly traces) —
   it is not fragile to the null; it **cannot discriminate**. Say "raw FC cannot separate the
   cognitive trace from disease/global structure", not "raw FC fails a control".

**Weakness of simple:** it is **information-poor** — it can register that a band moved, but
is silent about the *scale-structure*, the *location*, and the *higher-order pathways*.
**No band is "false"** (the cophenetic result is not ground truth): bands are *characterized*
by scale-shape and anatomy — β (deep-multiscale), α (single-scale), δ (fires-but-shallow),
θ/low-γ (silent).

**❌ What we do NOT claim (guardrails — violating these creates the contradictions):**
- **Not** "only LRG/multiscale sees the trace" / "multiscale-exclusive." Raw FC (β trace
  p=.024) and Grassmann both detect β. The edge is selectivity + scale-shape, not exclusive
  detection.
- **Not** "inference in β alone" (that was the dense single-τ result; on `mst@0.20` it is
  δ/α/β).
- **Not** the δ inference-specific effect — it is a partial-correlation artifact (cut it).
- **Not** "drift survivor / survives both nulls" (drift retired).
- **Never** rank bands by the ×-null ratio.

Status tags below: ✅ = verified on `mst@0.20` this session; ⚠️ = carried from the
whole-graph ρ_sym pipeline, not yet re-derived on `mst@0.20` (downstream of the same
cophenetic distances, so expected to hold); ❌ = do not claim.

---

## RESULT 1 — A held, multiscale trace of reasoning

**Headline.** Transitive-inference reorganizes the resting network, and the change
*persists* into post-task rest — in β and α. β carries a **scale-invariant (genuinely
multiscale)** trace that concentrates in **orbitofrontal cortex**; α carries a
**single-scale** trace with no anatomical home.

| # | subpoint | number | status |
|---|---|---|---|
| 1.1 | Post-task rest holds the task hierarchy | 10/10 patients rest_post closer to task than rest_pre | ⚠️ (robust descriptive) |
| 1.2 | **β trace clears at every scale — scale-invariant** | 16/16 scales, p .014 (τ_min) → .001 (meso) | ✅ |
| 1.3 | **α trace is mesoscale-tuned — single-scale** | 12/16 scales, p .024 → .007; effect peaks at s≈5 | ✅ |
| 1.4 | Clean nulls | θ 0/16, low-γ 0/16 | ✅ |
| 1.5 | Not-claimed bands are artifacts | δ 8/16 patchy (fails τ_min .216 & meso .116); high-γ 4/16 coarse-only | ✅ |
| 1.6 | Multiscale value = **selectivity** | raw FC fires 4 bands; cophenetic fires only α/β; spectral resistance dead in all | ✅ |
| 1.7 | Multiscale value = **scale-shape** | β flat across τ, α peaked — only the sweep shows this | ✅ |
| 1.8 | **β trace is DELOCALIZED on mst@0.20 — no OFC (or any) home at any scale** | OFC +91k(p.50)@s1 → −322k(p.59)@meso; nothing clears BH at s=1/2.83/5.65 (`16_localization`) | ✅ (delocalized) — supersedes the whole-graph "β→OFC" (degenerate-graph artifact) |
| 1.9 | β carried by **healthy cortex** (gray–gray), not the SOZ | gray–gray clears; white–white / SOZ do not | ⚠️ |
| 1.10 | **Held, not replayed** — sustained proximity, no bursts/sequences | reinstatement 10/10, p=.001 (LOO .002) | ⚠️ |
| 1.11 | Over the SOZ, **β spares, α recruits** | α SOZ–SOZ +0.41, p=.005; β generic there | ⚠️ |
| 1.12 | **Inter-patient spread is explained by implant laterality** | β trace ∝ left-contact fraction, ρ=+0.685, p=.029; α not lateralized (ρ=0.16 ns) | ✅ |

**One-line for the talk:** *"Two bands hold the trace; only β is genuinely multiscale
(present at every scale), and only β has an address — orbitofrontal cortex."*

---

## RESULT 2 — Offline abstraction: encoding vs inference

**Headline.** The persisting reorganization contains not only the premises the patient was
*shown* (**encoding**) but the order the patient *inferred* (**inference-specific**).
Encoding and inference occupy **different cortex** — OFC vs cingulate.

| # | subpoint | number | status |
|---|---|---|---|
| 2.1 | **Encoding persists** | β 16/16 (p_meso .003); α 3/16 (τ_min .042, meso .032) | ✅ |
| 2.2 | β encoding is scale-invariant + hierarchy-cleanest; **α encoding is visible to simple tools** | α: raw FC .032, clustering .010 | ✅ |
| 2.3 | **Inference-specific persists — NOT β-only** | β 7/16 (.005), α 7/16 (.024) | ✅ |
| 2.4 | **δ inference-specific is CUT** — partial-correlation artifact | δ 10/16 but δ null for encoding & trace | ❌ |
| 2.5 | Inference multiscale value = **localization, not detection** | raw FC detects β inference (.010) but is placeless | ✅/⚠️ |
| 2.6 | **Encoding → OFC** (learning sets the anchor) | q = .010 / .040 (with/without SOZ) | ⚠️ |
| 2.7 | **Inference → cingulate** (away from OFC) | q = .030 / .040 | ⚠️ |
| 2.8 | **low-γ encoding → cingulate** — a focal memory trace whole-brain averaging hides | ρ_sym +0.25, q=.035, 8/8 sampled | ⚠️ |
| 2.9 | Inference is **not** a recording-length effect | β inference vs test/learn duration ρ=+0.10, p=.78 | ⚠️ |

**Honesty note (do not headline the inversion):** for β it happens that *encoding* is the
scale-invariant/hierarchy-cleanest component while *inference* is raw-visible — but the
β-encoding "hierarchy-only" rests on a raw-FC near-miss (.053 vs .05) and an untested
Grassmann, so present it as **scale-invariance + simple-tool-visibility of α**, not as
β-encoding exclusivity.

**One-line for the talk:** *"Rest keeps not just the pairs seen but the order inferred —
and the two live in different cortex: premises in OFC, inferred order in the cingulate."*

---

## RESULT 3 — The propagator as an epileptogenic marker

**Headline.** The **same diffusion operator**, read within a single recording as seed
affinity, localizes epileptogenic tissue: the seizure-onset contacts form a
**strength-independent co-diffusing community**, not a contact-by-contact label.

| # | subpoint | number | status |
|---|---|---|---|
| 3.1 | SOZ = strength-independent co-diffusing group (relational, not per-contact) | above node-strength baseline | ✅ |
| 3.2 | Seed-affinity **AUC**, strongest in δ | δ .83, low-γ .82, β .745, α .61 | ✅ |
| 3.3 | Beats matched-strength fake-SOZ null | δ 7/10, low-γ 8/10, β 8/10 | ✅ |
| 3.4 | **Band dissociation** — epilepsy = δ/low-γ/β; trace = α/β; β = the bridge | δ marks epi but carries **no** cognitive trace | ✅ (consistent w/ R1) |
| 3.5 | τ (multiscale) lifts **ranking, not precision** | multi>single AUC 7–8/10; single-band prec@5 = .40 | ✅ |
| 3.6 | Precision comes from **band fusion**, not τ | fused prec@5 = 60% (~7× base rate) | ⚠️ |
| 3.7 | **Two populations** — 8 co-diffusing community + 2 right-hemi hub | community AUC up to .99; hubs fail (.49/.57), strength rescues | ⚠️ |
| 3.8 | Calibrated leave-one-patient-out detector | median AUC .87, 9/10 above chance | ⚠️ |
| 3.9 | **Robust to sparsification** (dense → mst@0.20 reproduce) | δ .80→.83, low-γ .74→.82, β .69→.745 | ✅ |
| 3.10 | Scope: SOZ **labels**, not surgical outcome; marker not detector | — | ✅ |

**One-line for the talk:** *"The same operator that carries the cognitive trace reads the
seizure zone — a strength-independent diffusion community, strongest in δ; a triage marker,
not a diagnostic."*

---

## Consistency audit (checked, no contradictions)

- **α across R1/R2.** R1: α is a real trace (12/16). R2: α encoding is simple-tool-visible.
  Consistent — α is the *single-scale (mesoscale-peaked)* carrier, real but resolvable by
  ordinary measures; β is the *scale-invariant* one. State α this way in both.
- **δ across R1/R2/R3.** δ is a trace-null (R1), an inference artifact to cut (R2), and the
  strongest epilepsy band (R3). **Not a contradiction — it is the band dissociation**
  (cognition ≠ epilepsy bands), a selling point. β is the only band in both sets.
- **β everywhere.** Genuine multiscale trace (R1) + encoding & inference persistence (R2) +
  the one cognitive band that also marks epilepsy (R3). Coherent — β is the bridge.
- **Multiscale framing.** Selectivity + scale-characterization + localization, NOT exclusive
  detection — identical in all three. (This is the fix for the earlier "cophenetic-unique"
  overclaim; the directives + diary must be aligned to it.)
- **Null.** Matched-strength only, everywhere. Drift retired everywhere.

## Genuinely open (bounded — everything else is settled)

1. **Localization + tissue + SOZ-divergence + reinstatement are whole-graph (⚠️).** Re-derive
   on `mst@0.20` for full self-consistency; expected to hold (same cophenetic distances). Part
   of task #13.
2. **Grassmann not re-run on `mst@0.20`.** Keep as a whole-graph companion or drop from the
   talk; do not attribute it to the sparsified backbone.
3. **α mesoscale encoding is fragile** (.032 vs .053 across runs). Present α encoding as
   "τ_min-solid, simple-tool-visible"; don't lean on the α-meso cell.
