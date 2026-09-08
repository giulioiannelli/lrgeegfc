---
name: 2026-07-23_old-vs-current-headlines-panoramic
kind: report
era: IMCOH_ABS_COHORT_N10 (mst@0.20 recovery / τ-sweep)
status: current
created: 2026-07-23
scope: The point-of-the-situation the user asked for — a paragraph-by-paragraph diff of the CURRENT Overleaf Results (results_sec_1/2/3.tex, written on the dense graph at a single τ) against what the τ-sweep + sparsification + localization-statistic passes have since established. What STANDS, what SHARPENED, what WEAKENED, what is DEAD, and the exact source for each. Written to make the Overleaf rewrite mechanical. Reconciles the old headlines with scripts 12–27 in scripts/01_compute/sparsified_arc/ and the 2026-07-13 localization verdict.
pointers:
  - .agents/preprint/overleaf/results_sec_1.tex           # §1 trace — biggest change
  - .agents/preprint/overleaf/results_sec_2.tex           # §2 enc/inf — β-alone + cingulate at risk
  - .agents/preprint/overleaf/results_sec_3.tex           # §3 marker — robust, backbone-contested
  - .agents/reports/2026-07-13_localization-cohort-statistic-verdict.md   # the verdict that killed β→OFC
  - .agents/preprint/directives/writing_directive_2026-07-13_sparsified-recovery-OVERVIEW.md  # STALE §5 (says localization stands)
  - scripts/01_compute/sparsified_arc/                    # 12–27: the analyses behind every "new" cell below
---

# Old headlines vs current results — the point of the situation

## Head

The τ-sweep pass was **not additive**. It kept and sharpened the spine of §1 (a held β/α trace, now scale-resolved and better-controlled) and left §3 (the epilepsy marker) essentially intact — but it **retracted the entire anatomical-localization layer** that both §1 and §2 are built around (β→OFC, encoding→OFC, inference→cingulate all fall to a 2% noise floor), it **killed the "inference in β *alone*" headline** (δ/α/β now), and it **downgraded "only the hierarchy sees the trace"** to selectivity, not exclusive detection. Your instinct that "localization seems gone" is correct and important: the old localization was measured on the dense, degenerate graph with a cohort statistic that was dead-on-arrival, and when it is redone properly on the backbone the anatomical homes vanish. What replaces them is much thinner — β only **coarsely left-lateralized**, low_γ sitting on the **seizure zone**, α **placeless** at every anatomical resolution. The rewrite is therefore a *demotion of the "where" story* and a *promotion of the "how it scales" story*, not a cosmetic edit.

## ⚠ Correction after user review (2026-07-23) — two rows below are retracted/unsupported

Two things in the table below are now known to be wrong or unsupported and supersede the corresponding rows. Both came out of quantifying claims the user (correctly) refused to take on faith.

- **Scale labels ("α mesoscale") are RETRACTED — unquantified and false.** From `26_diffusion_spatial_reach_mst020` (cohort span ≈ 85–90 mm, pitch ≈ 3.5 mm): α's "peak" at s≈5.65 → physical reach ℓ ≈ **86 mm ≈ the whole implant span** (global, not meso); even τ_min (s=1) → ℓ ≈ 65 mm ≈ 3/4 of the span; reach saturates to full span by **s≈2**; the only local regime (ℓ<20 mm) is s<1, where the hierarchy degenerates to raw edges. So the entire trace-gate range is physically near-global, and "α clears 12/16 scales" is a scale-*coverage* statement, NOT evidence of a mesoscale community size. Never label a scale micro/meso/macro without ℓ(s) or a cluster count. Structural corollary: the hierarchy is non-trivial only for s>1, exactly where reach is already global → "multiscale hierarchy" and "spatially local" are mutually exclusive here, which is the mechanistic reason localization washes out (supersedes/explains §1.4–§1.5, §2.3–§2.4).

- **The hierarchy does NOT subsume raw FC — the existential caveat (supersedes §1.3, §2.2 framing).** The non-circular residualization test `21_coph_beyond_raw.py` (run at s=1 / 5.6 / 30) returns, on BH-corrected q: `coph|raw` (the hierarchy-only residual) **never clears** for α or β at any scale (α p=0.053–0.161; β q=0.064–0.105), while `raw|coph` (the raw-only residual) **clears strongly** (α p≤0.009, 8/10; β q=0.032–0.053). For α the persistent trace is entirely in the raw edges — the hierarchy adds nothing; for β they are complementary/symmetric (obs 0.155 vs 0.156), never subsumption. By the script's own criterion this is the *weaker claim at best, the reverse for α*. **The strong thesis "the multiscale hierarchy reveals persistent structure raw FC cannot see" is UNSUPPORTED.** The only honest surviving contributions are (i) selectivity (the hierarchy rejects raw's δ/θ/low_γ false positives) and (ii) scale-characterization — both modest. P0 of the re-investigation (τ-resolved beyond-raw across all 16 scales) decides whether even a single scale rescues subsumption; until then the paper must not claim hierarchy-exclusive detection. (Caveat: rank-linear residualization removes only monotone raw dependence — but this is the strong-form negative, not borderline.)

## What changed at the method level (the frame for everything below)

- The old Results read one dense, unthresholded diffusion propagator at a single time τ = 1/λmax. On the fully-connected graph that propagator is **degenerate at τ_min** (heat-kernel geometry ≈ raw connectivity — audit_174), so "diffusion" was buying nothing.
- The fix: sparsify to a **backbone** (maximum spanning tree + strongest 20% of edges, `mst@0.20`) and **sweep τ** across scales (dimensionless s ≡ τλmax, from s=1 at τ_min up to s≈180). On the backbone the propagator is non-degenerate and τ becomes a real multiscale scanner.
- **The old numbers survive at the fine scale.** At s=1 the published α gate reproduces *exactly* (p=0.024) and β reproduces closely; sweeping τ up is an *addition* that resolves the bands by their scale signature.
- **Only null now is matched-strength.** Drift / timeshift / within-baseline nulls are retired (the task is directional and a real trace is monotonic, so those nulls are degenerate with the alternative).
- **Read per-scale, never best-scale.** Reporting the single most-significant scale falsely lights up δ and high_γ at coarse τ where the surrogate median collapses — artifacts, not traces.

## The master diff table

Status legend: ✅ STANDS · ▲ SHARPENED · ⚠ WEAKENED/AT RISK · ✖ DEAD · ↔ KEEP-AS-COMPANION.

| # | Old headline (as written in the .tex) | Status | Current statement | Source |
|---|---|---|---|---|
| §1.1 | Reasoning leaves a held β/α trace; θ null | ▲ | β **scale-invariant 16/16** (p_meso 0.001–0.006, obs_med 0.20–0.32; Friedman p=0.17); α **mesoscale 12/16** (peak s≈5, p 0.006–0.04); θ 0/16; low_γ 0/16; τ_min recovers the old gate | `13` / `ms_mst020` |
| §1.2 | β median ρcoph +0.20, 6/10, p=0.032; α +0.10, p=0.024 | ✅ | Recovered at s=1 (α p=0.024 exactly); now a *slice* of a per-scale curve | `13` |
| §1.3 | Band-specificity invisible in raw edges — "emerges only in the hierarchy" | ⚠ | Raw FC **does** trace α/β (and δ/low_γ false positives): it is **selectivity + scale-characterization**, NOT exclusive detection. Apples ladder: only cophenetic-meso is band-selective (β p=0.001); strength p=1; spectral resistance dead | `14`,`22` |
| §1.4 | **β trace concentrates in OFC** (q=.010–.015); PFC/sensorimotor depleted; "anatomically precise" | ✖ | **DEAD.** β at the 2% noise floor in every system, under all 4 cohort statistics + LOO, **density-invariant** (0 survivors at frac .10/.20/.50). Replaced by a **coarse LEFT-hemisphere** tilt only (z=+1.38, LOO p=.005, 5/7 patients, meso s≈11–16; strict-BH q=.119 marginal) | `16`,`17`,`23`; verdict §2,§7,§9 |
| §1.5 | α "a trace without an address" | ⚠ | Placeless at every *anatomical* granularity — **but the intrinsic-community test at α's own scale (s≈5) has not been run** (see §"your τ↔space instinct") | verdict §9 |
| §1.6 | β carried by **healthy cortex**; white–white (+0.10) and SOZ (+0.08) do not clear | ⚠ | **Complicated on the backbone.** β gray_gray clears (q↓0.015) but **wm_wm (q↓0.016) and cross_wm (q↓0.015) also clear** → β is *broad*, not gray-exclusive. The "white matter doesn't carry it" half no longer holds cleanly | `19` |
| §1.7 | SOZ divergence: β steers around SOZ, **α recruits SOZ** (SOZ–SOZ +0.41, p=0.005) | ⚠ | β-spares-SOZ **holds** (SOZ-independent, 3 ways — below); but **α-recruits-SOZ does NOT reproduce** on mst@0.20 pair-class (α epi_epi null, q=0.237). The α half of this beat is at risk | `19`,`25` |
| §1.8 | β SOZ-independent / healthy-core carrier (decimation control) | ✅ | Confirmed **3 ways**: pair-class non-SOZ 16/16 vs epi_epi 1/16 scales; node-exclusion (drop SOZ, rebuild) β 14/16 (p=0.002); size-matched decimation. Mirror: low_γ trace **dies** without SOZ | `19`,`25`; verdict §8 |
| §1.9 | Held-not-replayed reinstatement (10/10, p=0.001) | ✅ | Unchanged (inherits the static trace's control) | as published |
| §1.10 | Grassmann table (β both, α coph-only, low_γ Grassmann-only) | ↔ | Keep as a **methods companion** — dense/whole-task, matched-strength-verified (audit_66), **NOT re-derived on the backbone**. Do not re-attribute to mst@0.20 | audit_66 |
| §2.1 | Encoding persists in α and β (p=0.014 / 0.032) | ▲ | Holds. β encoding **scale-invariant 16/16** (p_meso 0.003); α marginal 3/16. (An earlier "encoding washed out" was a T_learn **bug**, now fixed) | `05`,`14`,`20` |
| §2.2 | **Inference-specific in β ALONE** (p=0.010; next band p=0.25) | ✖ | **DEAD.** Clears in **δ (10/16, p 0.019), α (7/16, 0.024), β (7/16, 0.005)** — β strongest but not exclusive. And raw FC detects β inference too (p=0.010) | `05`,`14` |
| §2.3 | **Encoding anchors in OFC, learning sets the anchor** (q=.010/.040) | ✖ | **DEAD.** No system clears at any scale | `17` (system, all 16 scales) |
| §2.4 | **Inference concentrates in the cingulate** (q=.030/.040) | ✖ | **DEAD.** Inference localizes to no system (min q≈0.13) | `17`,`23` |
| §2.5 | low_γ **cingulate** memory trace (ρ_sym +0.25, q=.035, 8/8) | ⚠ | Direction survives but the **home moved**: low_γ encoding concentrates on **PFC/occipital/limbic + SOZ** (q≈0.020–0.040); cingulate appears only inside the low_γ coph-value set, not as the primary home. Needs precise re-derivation before any cingulate claim | `17`,`23`,`24` |
| §2.6 | Cingulate multiplexing (memory@low_γ, inference@β, same cortex) | ✖ | Rests on §2.4 + §2.5, both retired/moved → the multiplexing beat collapses | — |
| §2.7 | Duration control: β inference ρ=+0.10, p=0.78 | ✅ | Holds for the β inference component (not length-driven, all p>0.16) | `18` |
| §3.1 | Marker: δ≳low_γ≳β, distant-seed discovery, calibrated LOO, two-populations, prec@5 60% | ✅ | Reproduces — but **on the TMFG (planar) backbone** (§3 tex is already two-backbone). TMFG detector AUC 0.909, prec@5 0.60; cohort recovery (Pat_15 0.49→0.88) is **TMFG-specific** | `06`,`27` |
| §3.2 | (implicit) multiscale τ sharpens the marker | ⚠ | Honest: multiscale does **not** beat fixed single scale out-of-sample; prec@5 is a **band-fusion** effect, not a τ effect. τ helps ranking/AUC only | `16b`,`27` |

## What was lost (be honest about this in the rewrite)

- **The whole "where" spine.** β→OFC (§1), encoding→OFC and inference→cingulate (§2) were three of the paper's most concrete, quotable results and they are gone. §1's Fig-trace2a (the OFC pial map) and §2's Fig-enc_anatomy (a,c panels) no longer have a result behind them. This is the single largest structural loss.
- **"Inference in β alone."** The §2 opening headline — arguably the sharpest cognitive claim in the paper — does not survive; it is now a three-band effect with β merely strongest.
- **"Only the hierarchy detects the trace."** Softened everywhere to selectivity + scale-characterization + (residual) localization. Raw FC detects α/β; Grassmann detects β. The exclusivity framing must come out.
- **The clean tissue story.** "β carried by healthy cortex, not white matter" is now "β is broad" (wm_wm clears). And "α recruits the SOZ" does not reproduce. Both were vivid; both need softening.
- **A fine anatomical address for any cognitive band.** Nothing focal survives for β or α. The only focal anatomy left is **low_γ → seizure zone**, which is the *disease* band (with an anchor-vs-trace caveat: low_γ→SOZ may be stable SOZ tissue, not a task trace).

## What was gained (the new load-bearing story)

- **τ as a scale axis, and band scale-signatures.** β is scale-invariant (16/16, Friedman p=0.17); α is mesoscale-tuned (12/16, peak s≈5). This *dissociation* is a genuinely multiscale result no single-scale tool can produce — this is the new headline of §1.
- **The apples controls ladder.** On the *same* backbone, only the cophenetic mesoscale read-out is band-selective; raw edges fire non-selectively, strength is degenerate (p=1), graph geodesic fires the wrong bands, and the **spectral resistance distance is dead everywhere** (the strongest foil — a Laplacian read-out that is multiscale-blind — finds nothing). This *demonstrates* the thesis instead of asserting it.
- **The spread is explained, not just controlled.** β trace strength tracks the left-hemisphere contact fraction (ρ=0.68, p=0.029; per-scale peak ρ=0.74 at s≈16); who traces is set by which cortex the electrodes sample. α shows no laterality.
- **β SOZ-independence, three ways** — a clean "cognitive, not epileptic" result.
- **Sparsification-robustness** — β holds on dense / every mst@f / percolation / TMFG; α needs only "keep the strong edges."
- **The τ→physical-space mapping** (script 26) — a quantitative answer to your "τ is linked to spatial scale" instinct; see below.

## Your τ↔space instinct — quantified, and it cuts both ways

You are right that τ maps to physical spatial scale. Script `26_diffusion_spatial_reach_mst020.py` computes it directly: ℓ(s) = the largest physical distance (mm) between two contacts that *communicate* (heat ≥ its equilibrium share 1/N). The findings:

- **s = 1 (τ_min) is the raw↔communication boundary.** For s<1, e^{−τL} ≈ I − τL reads single raw edges; multi-hop communication only opens for s>1. So the fine operating point is exactly where "diffusion" begins.
- **The finest connected functional length is ~15 mm, NOT the 3.5 mm electrode pitch.** The |ImCoh| backbone is *functional*, not spatial (median edge ~35 mm; only ~2–4% of edges are same-shaft ≤5 mm), so the diffusion never resolves the electrode pitch — it is spatially delocalized from the start (~5× pitch at the finest scale).
- **ℓ(s) grows to ~85–92 mm ≈ the whole implant span** (cohort median span ≈ 85 mm, up to 109 mm) by s ≳ 15.

So the mapping is real — and it is precisely *why* anatomical localization washes out: by the mesoscale where β and α live, the operator already communicates across most of the implant. **Being scale-localized (α at s≈5) does not imply being anatomically localized** — at s≈5 the physical reach is already several centimetres, and those few-cm co-diffusing neighbourhoods are patient-specific, so they don't line up with any named region. That is the honest reconciliation of "α fires at a precise scale" with "α is anatomically placeless."

**The one test this leaves un-run** (your original α challenge): the localization suite tested *imposed anatomical parcels* (system / hemisphere / lobe / SOZ), never the diffusion's **own intrinsic communities** at α's characteristic scale. Whether the α trace concentrates in identifiable few-cm co-diffusing communities — even if those communities have idiosyncratic anatomy — is a distinct, potentially positive result the current tests are structurally blind to. It would turn "α placeless" into "α has a community home, idiosyncratic anatomy." Still worth running; keep α marked "placeless (anatomical); community test pending" in the interim.

## Open decisions that block a clean rewrite (flags, not questions)

1. **Backbone posture for §3.** The Overleaf §3 (updated 2026-07-15) and the detector code (`27`) are on the **two-backbone** story: marker = planar **TMFG**, trace = strength-preserving **mst@0.20**, carrying the strong numbers (β AUC 0.87, fused 0.91, both right-hemisphere hub patients recovered). A **2026-07-16 decision reverted the *talk* to mst@0.20-only** (marker β 0.745, Pat_15 back at chance, cohort-recovery retired) and explicitly flagged that **§3 + the verdict + the sparsifier-settled memory are still on TMFG, pending a scope call.** So the paper's backbone is *undecided*: keep TMFG (best §3 numbers, two backbones to justify in Methods) or collapse to single mst@0.20 (cleaner, weaker §3). This must be settled before §3 is finalized.
2. **The α intrinsic-community analysis** (above) — run before finalizing α's status, or ship "placeless (anatomical); community test pending."
3. **low_γ → cingulate vs frontal/SOZ** (§2.5) — needs a precise re-derivation on the valid statistic before any low_γ anatomy is written; and the anchor-vs-trace caveat must be stated.

## Rewrite skeleton (how the three sections should land)

- **§1 — reframe from "where" to "how it scales."** Lead paragraphs: held β/α trace, now read across τ (β scale-invariant, α mesoscale). Promote the **controls ladder** to the section's centrepiece (band-selectivity is a multiscale-geometry property, demonstrated). Convert the spread paragraph into the **laterality explanation**. Replace the OFC paragraph with the honest localization: **no fine home; β coarsely left-lateralized; α placeless; low_γ on the SOZ (disease).** Keep healthy-core carrier but soften to "β broad, SOZ-independent" (drop gray-exclusive and α-recruits-SOZ). Keep held-not-replayed and the Grassmann companion.
- **§2 — reframe from "β alone / anatomy" to "encoding vs inference as scale-differentiated components."** Encoding: β scale-invariant, α basic-metric-visible. Inference-specific: three bands, β strongest — drop "β alone." **Remove the OFC and cingulate localization paragraphs and Fig-enc_anatomy panels** (or replace with the low_γ→SOZ/frontal result once re-derived). Keep the duration control.
- **§3 — smallest edit, once the backbone is settled.** Marker reproduces; add the honest τ contribution (ranking not precision; prec@5 is band-fusion). Freeze the TMFG-vs-mst@0.20 wording to whatever decision #1 lands on.

## Note on the directives themselves

The three 2026-07-13 conceptual-diff directives and the OVERVIEW §5 are **stale on localization** — they were written the morning of 2026-07-13 and still say "β→OFC / encoding→OFC / inference→cingulate stand as written" and "localization is the headline value-add for inference." The **2026-07-13 evening verdict overturned all of that.** When the cascade propagates, OVERVIEW §5 and the sec1/sec2 "What STAYS" localization bullets must be corrected to the per-band verdict above before any agent rewrites from them.
