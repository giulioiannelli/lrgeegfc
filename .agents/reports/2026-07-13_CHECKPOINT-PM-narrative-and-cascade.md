---
name: checkpoint-2026-07-13-PM-narrative-and-cascade
kind: verdict
era: IMCOH_ABS_COHORT_N10 (mst@0.20 recovery)
status: current
created: 2026-07-13
scope: THE resume anchor after compaction (supersedes the AM localization checkpoint). ALL missing results finalized on mst@0.20, scale-swept; the multiscale-value case is settled into 4 pillars; the talk narrative doc is fully folded. ONE thing left — the cascade (task #29) into the 4 directives + 4 headlines + the R1/R2 tables of the settled doc. Every number here is fresh this session; do not recompute.
pointers:
  - .agents/talk/2026-07-13_multiscale-narrative-audit.md          # talk source of truth — FULLY folded
  - .agents/preprint/established_results/2026-07-13_settled-three-results.md  # §0 folded; R1/R2 tables NOT yet
  - .agents/reports/2026-07-13_CHECKPOINT-localization-scale.md    # AM checkpoint (superseded by this)
  - data/sparsified_arc/figures/{fig_scale_signatures,fig_method_ladder}.pdf
---

# CHECKPOINT 2026-07-13 PM — narrative settled, cascade pending (RESUME HERE)

## Head

Every missing result is finalized on `mst@0.20`, scale-swept, matched-strength (scripts
16–22). The talk is titled **"The multiscale nature of inference"**; its narrative doc is
fully folded and defensible. The multiscale-value case rests on **4 pillars** (band-
selectivity, higher-order β, epilepsy, scale-invariance). **Localization is dead** (→
distribution). **Raw FC is complementary/robust, NOT fragile** — its weakness is
*discrimination*, not null-failure. **ONE task remains: the cascade** (#29) — propagate the
locked numbers into the 4 directives + 4 headlines + the R1/R2 tables of the settled doc,
which still carry dead claims (β→OFC, inference-β-only, localization value-add).

## FINALIZED RESULTS (this session — locked numbers, do not recompute)

**Localization — DEAD (`17_localization_arc_mst020.py`; data `localization_arc_mst020/`).**
Trace + encoding + inference, per-system + SOZ, all 16 scales, matched-strength R=200.
**0 BH-clearing cells, every target, every band.** Encoding→OFC best p=.156 (0/16);
inference→cingulate β p=.23. Only survivor: **low-γ→SOZ q=.010** (disease). β→OFC was a
degenerate-graph artifact. ⇒ trace is **distributed** (~88% co-move), not placed.

**Tissue pair-class (`19_tissue_pairclass_mst020.py`; `tissue_pairclass_mst020/`).**
β trace **tissue-DISTRIBUTED**: gray_gray q.015 (14/16), cross_wm q.015 (16/16), wm_wm q.015
(16/16), nonepi q.015 (15/16) — NOT gray-exclusive. β **spares SOZ** (epi_epi clears only
1/16). **α does NOT recruit SOZ** on backbone (epi_epi +0.004, p.19 — softens old "α
recruits +0.41"). low-γ leans epi_epi (p.024). Only β clears BH in tissue.

**Laterality (`18_laterality_duration_mst020.py`; `laterality_duration/`).** β trace ∝
left-contact fraction **positive at all 16 scales**, sig 6/16, peak ρ=+0.74 p=.014 @s≈16;
**α flat 0/16**. Scale-robust, not a single-scale coincidence.

**Duration (`18`).** β inference-specific vs test/learn ratio **ns everywhere**: s₁ ρ+.03
p.93; mesoscale ρ+.48 p.16 (ns trend, honest caveat); peak ρ+.14 p.70. Not a length effect.

**Scale-invariance hardening (`20_hardening_checks.py`; `hardening/`).**
- β trace: Friedman p=**0.17 (ns)**, trend p=0.56 (ns), CV 0.25, 16/16 → **SCALE-INVARIANT**
  (say "no characteristic scale / consistent with scale-invariance"; n=10 → fail-to-reject).
- α trace: Friedman p=**0.004 (sig)**, 12/16 → scale-DEPENDENT (single/mesoscale).
- β encoding: Friedman .0001, rising trend p=.01 → spans all scales but **coarse-strengthening**
  (NOT flat).
- **Leg-2 "encoding cophenetic-only" did NOT harden**: raw FC .053 (R=1000-confirmed), LOO
  **4/10 fragile**; cophenetic 10/10 LOO robust. ⇒ "robust reader", NEVER "exclusive".

**Incremental coph vs raw (`21_coph_beyond_raw.py`; `coph_beyond_raw_s*`).** NOT subsumption
— **complementary**. coph⊥raw β traces at meso/coarse (p.032, effect +0.15→+0.24, BH-marginal
q.064; δ@coarse clean q.039); raw⊥coph traces in **all** bands → raw FC carries persistent
structure the hierarchy lacks. Higher-order β = a **coarse-scale** phenomenon.

**Band-selectivity (`14_controls_ladder` + `22_raw_dense_band_selectivity.py`).** Trace,
read-out × band: **cophenetic fires ONLY α(.007)/β(.001) = 2/6**; raw FC fires δ/α/β/low-γ
(4/6, dense too — θ .053 borderline big effect +.165); geodesic WRONG bands (δ/low-γ);
clustering misses β; resistance/strength dead. **Raw FC passes the null (robust) but cannot
DISCRIMINATE** — conflates cognition with the δ/low-γ disease bands (ties R3). NB: `14` was
edited (confound-fix) to read all descriptors on backbone B, not dense.

**Companions (decisions, not ported).** Reinstatement (held-not-replayed) = combinatorial but
windowed-temporal, orthogonal to backbone/scale → temporal companion at τ_min, NOT ported.
Grassmann → DROP from mst@0.20 story (whole-graph companion only).

## MULTISCALE VALUE — the 4 pillars (honest grades)

1. **Band-selectivity** (non-circular): coph 2/6 (α/β only) vs raw 4/6; cross-validates R3.
2. **Higher-order β** (non-circular): coph⊥raw β persists, coarse-strengthening (`21`).
3. **Epilepsy** (non-circular, external ground truth): community beats strength on SOZ AUC.
4. **Scale-invariance** (partly self-referential — acknowledge): β Friedman .17 vs α .004.

**Through-line (say this):** raw FC is a **robust, complementary** detector — the hierarchy
earns its place by **discriminating + characterizing + adding a higher-order view**, NOT by
seeing what raw FC is blind to. "raw FC is insufficient/fragile" = OVERCLAIM, banned.

## R2 reframe (encoding vs inference)

Dissociation is now by **SCALE**, not anatomy (delocalized) and not band-exclusivity (both
α/β): **encoding fine-legible + scale-invariant/coarse-strengthening; inference
mesoscale-emergent** (ns at s₁: β .097/α .278/δ .188 → meso .005). Offline-abstraction =
R2 climax: the inferred order is held as a **sustained, distributed, mesoscale** abstraction.
δ inference-specific CUT (partial-corr artifact).

## DOCS STATE — what's folded vs pending

- **Talk narrative** (`.agents/talk/2026-07-13_multiscale-narrative-audit.md`) — **FULLY
  FOLDED**: 8 points + offline-abstraction climax + ★★ fundamental headline (scale-invariance)
  + the two non-circular wins + locked-numbers table + do-not list. Source of truth for the
  presentation agent.
- **Settled doc** (`…/established_results/2026-07-13_settled-three-results.md`) — **§0 FOLDED**
  (point 1 scale-invariance, point 2 localization-dead-all-components, point 3 higher-order/
  coph⊥raw, point 4 band-selectivity). **R1/R2 tables NOT yet updated** (see cascade below).
- **4 directives + 4 headlines** — **NOT cascaded** (still carry dead claims).

## ⬜ THE PENDING CASCADE (task #29 — the only work left; actionable)

**(1) Settled doc R1/R2 tables** (`2026-07-13_settled-three-results.md`):
- R1 headline: "concentrates in OFC" → delocalized/distributed.
- 1.6: sharpen to "coph 2/6 α/β vs raw 4/6" (band-selectivity numbers).
- 1.8: already delocalized ✅.
- 1.9: "gray-gray not SOZ" → **tissue-DISTRIBUTED (gray+white+cross all q.015)**; status ✅.
- 1.11: "α recruits SOZ +0.41" → **α does NOT recruit on backbone (epi_epi +.004 ns)**; ✅.
- 1.12: laterality → **scale-robust (all 16 scales, sig 6/16, peak ρ.74)**; ✅.
- 1.2: add Friedman scale-invariance.
- One-line: drop "only β has an address — OFC".
- R2 headline: "different cortex OFC vs cingulate" → **delocalized; scale-dissociation**.
- 2.5: "value = localization" → **DEAD; value = scale-dissociation + higher-order**.
- 2.6 encoding→OFC, 2.7 inference→cingulate: → **DELOCALIZED (❌)**.
- Add scale-dissociation row (encoding fine-legible; inference mesoscale-emergent).
- 2.9: duration → mst@0.20 numbers (ns, s₁ ρ.03 p.93).
- Update "Genuinely open" section: localization/tissue/laterality/duration now DONE.

**(2) 4 directives** (`.agents/preprint/directives/writing_directive_2026-07-13_*.md`):
- **sec2-encinf** — WORST offender: rewrite "value-add = localization / places encoding in
  OFC, inference in cingulate" → **scale-dissociation + offline abstraction + distribution**.
- **OVERVIEW** — kill "value-add = localization", β→OFC/encoding→OFC/inference→cingulate rows;
  add scale-invariance + band-selectivity + coph⊥raw + complementary-not-fragile.
- **sec1-trace** — β→OFC (q.010-.015) → delocalized/distributed + scale-invariance + laterality.
- **sec3-epi** — mostly fine (epi localization survives, low-γ→SOZ); add β-spares-SOZ.

**(3) 4 headlines** (`.agents/preprint/headlines/`):
- **01_trace.md** — scope + body "β→OFC" → delocalized/distributed.
- **02_encoding_vs_inference.md** — "inference β-only… anchors OFC" → δ/α/β, scale-dissociation,
  delocalized.
- **00_CORE_methodology.md** — "localizable" claims → distribution + band-selectivity.
- **03_epileptogenic_markers.md** — mostly fine; low-γ→SOZ is the surviving localization.

**Guardrails for the cascade:** matched-strength only (drift retired); read per-scale never
best-scale; no band "false"; "robust reader" not "exclusive"; "raw FC cannot discriminate"
not "raw FC fails/blind"; localization DEAD for cognition (→ distribution); low-γ→SOZ is the
ONE surviving anatomical signal.

## New scripts + data (this session)

`17_localization_arc_mst020` · `18_laterality_duration_mst020` · `19_tissue_pairclass_mst020`
· `20_hardening_checks` · `21_coph_beyond_raw` · `22_raw_dense_band_selectivity` ·
`fig_multiscale_narrative` (all in `scripts/01_compute/sparsified_arc/`). Data under
`data/sparsified_arc/{localization_arc_mst020, tissue_pairclass_mst020, laterality_duration,
hardening, coph_beyond_raw_s*, raw_dense_selectivity, figures, controls_ladder_beta_R1000}/`.
Env: `/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python`, pin OMP/OPENBLAS/MKL=1.

## Git / state

- All new scripts + `14` edit (—R/—tag + confound-fix) + docs edits **UNCOMMITTED**. `data/`
  gitignored. Branch `audit/cohort-n10-diagnostic`.
- Tasks: #23–#28 completed; **#29 (cascade) + #13 (regenerate preprint figures) pending**.

## Recommended next action

Run the cascade (task #29) per the actionable list above — settled-doc tables first (master),
then the 4 directives (sec2 first — it's the worst), then the 4 headlines. Then #13 (figures).
Nothing else needs computing; the science is settled.
