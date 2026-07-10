---
name: writing_directive_2026-06-25_headline-restructure-broadband-snr
kind: writing-directive
era: IMCOH_ABS × COHORT_N10
status: active
created: 2026-06-25
scope: directive for a parallel headline-rewriting agent — reframe R1 (the trace) from "β-specific" to "broadband, θ-exempt, SNR-gated", reconcile with R2-as-flagship. SAFE half only; band-taxonomy verdicts are PENDING live verification (Step 1 + Q2) and must stay placeholders.
pointers:
  - .agents/reports/2026-06-25_snr-band-taxonomy-handoff.md          # full source of truth
  - .agents/reports/2026-06-25_per-node-trace-anatomy-and-heterogeneity.md  # the SNR/detectability result
  - .agents/preprint/headlines/RESULTS_OUTLINE.md                    # the structure to edit
  - .agents/preprint/headlines/02_encoding_vs_inference.md           # R2 / N2 home (flagship, locked)
  - .agents/preprint/headlines/01_trace.md                           # R1 home (the file to reframe)
  - .agents/preprint/locked/VERDICT_LEDGER.md                        # DO NOT EDIT — cascaded separately
---

> **Head.** The Results skeleton is sound and **R2 (offline abstraction of a
> learned structure = the encoding/inference dissociation) stays the flagship** —
> do not move it. The one substantive rewrite is **R1's band story**: the offline
> trace is **NOT β-specific**. Under a measurement-detectability (SNR) lens it is
> **broadband, with θ the single genuine exception**; the apparent β-exclusivity
> was a detectability artifact. β remains the headline band **because it carries
> the cognitive content (R2), not because it is the only band with a bare trace.**
> Two sub-claims — *which* bands are cohort-consistent vs patient-specific, and
> whether the band story is *multiscale-specific* (R1.2) — are under live
> verification THIS session and must remain marked placeholders. Do not edit the
> locked ledgers.

---

## 1. LOCKED — safe to rewrite now (high confidence, will not move)

- **R2 stays flagship.** Encoding/inference dissociation; β carries the
  inference-specific consolidation component, duration-controlled, → OFC.
  Independent of this session's SNR work. Leave R2.1–R2.7 as is
  (they are already solid/flagged correctly). *(The τ≈2.6 "mesoscale-favouring"
  scale claim was retired 2026-07-09, PI call — read at τ=1/λmax only.)*
- **β trace → OFC** (R1.3). Matched-strength q=0.009, bilateral, hotspot-not-
  container. Locked. Keep.
- **The trace is broadband, θ-exempt.** Among high-SNR patients the trace clears
  its matched-strength null in most bands; **θ is the one band that does not, even
  at high SNR.** This is established from cached all-6-band data
  (`per_patient_per_band.csv`), not pending.
- **Per-patient heterogeneity = measurement detectability (SNR), not anatomy /
  implant / epilepsy.** Apparent "non-tracer" patients are SNR-limited, not
  counter-evidence — this DEFENDS the cohort claim. (Anatomy/implant/epi
  explanators were all tested and null; the explanator is task-vs-baseline SNR.)
- **R3 (epileptogenic-tissue markers)** unaffected by this session. Leave as is.
- **Narrow guardrail (Q2):** don't argue method-superiority from *magnitude* — the
  cophenetic `ρ_split` is SMALLER than the raw per-edge analogue in 5/6 bands (the
  hierarchy filters, it doesn't amplify; the bigger raw number is the trivial blob,
  see §1b). Argue superiority from **refinement / separability / localization**,
  not bigger numbers. See `.agents/reports/2026-06-25_raw-vs-multiscale-trace.md`.

## 1b. RESOLVED this session (Q2) — multiscale REFINES the trivial raw trace into structure

**The framing is NOT a magnitude contest (PI-corrected 2026-06-25).** The point of
the raw-vs-multiscale comparison: **raw FC shows a trace, but a blunt,
undifferentiated one that cannot be separated from noise; the multiscale filter
transforms that blob into refined, interpretable structure.**

- **Raw is the trivial, uninterpretable baseline.** raw_ρ is positive and floored
  at ≈ +0.11 in EVERY band — including **θ, where there is no real trace**. Raw
  cannot tell genuine persistence (β) from noise (θ); its *larger* magnitude is
  exactly its triviality — an undifferentiated positive overlap everywhere. So
  raw's bigger number is NOT a point against the hierarchy; it is the noise-blob
  we are filtering.
- **The multiscale filter refines it.** Pass the SAME construction through the LRG
  hierarchy and the blob does NOT survive as-is: **θ collapses to ~0** (correctly
  "no trace"), **β sharpens to +0.22 and stands out**, and the surviving trace
  acquires structure raw never had — the cross-band consistency taxonomy AND an
  anatomical home (β→OFC). The magnitude SHRINKS because the filter strips the
  trivial component raw cannot separate. (Do NOT frame this as "the hierarchy has
  bigger numbers" — it has smaller ones, and that is the point.)
- **β-specific null-clearing confirms it:** at β the multiscale trace clears
  matched-strength (p=0.005, 7/10) where raw falls just short (p=0.053).
- **Complementary evidence (R1.2's existing point):** the same progressive
  refinement appears comparing **spectral clustering vs the multiscale measure** —
  the raw/spectral structure is reorganized into a more refined one by the full
  multiscale read-out. (If wanted in THIS figure, the spectral-clustering
  intermediate is a small follow-up — flag, not yet run here.)
- So R1.2 becomes: **raw edges give a real-but-trivial trace, inseparable from
  noise (positive everywhere, uninterpretable); the multiscale read-out is what
  turns it into differentiated, localized, interpretable structure.** (Downstream
  OFC / encoding-inference claims raw cannot express by construction — out of
  scope for this edge-level foil, unaffected.)

## 2. PENDING — do NOT finalize; leave as clearly-marked placeholders

- **The 3-tier consistency taxonomy — NOW LOCKED (2026-06-25).** Step 1 finished;
  use this table as fact. Locked via audit_83 (cohort-median per-system trace ×
  matched-strength × both epi modes × shaft-collapse, all 6 bands) + audit_83b
  synthesis. Report: `.agents/reports/2026-06-25_per-band-consistency-taxonomy-lock.md`.

  | tier | bands | locked basis |
  |---|---|---|
  | **consistent** | **β→OFC** (alone) | cohort trace + shaft-robust + **LOO-robust** home; the ONLY band that clears every condition; + sensorimotor *depleted* (q=0.025) |
  | **consistent-weak** | **γ_l→PFC** | cohort trace + PFC lean, shaft-robust, but **FAILS LOO** at the β bar (loses BH-q significance in 6/10 patient drops; audit_83c) — a lean, not a hard home |
  | **patient-specific** | **α, δ, γ_h** | trace present but no cohort-consistent home — α has a clear cohort trace yet no surviving system; γ_h leans parietal but on a net-null band |
  | **absent** | **θ** | no net cohort trace (2/10, median ρ=−0.04) |

  Net: **β is the ONLY band with a fully-locked, LOO-robust, anatomically-localized
  cohort trace** — which sharpens β's specialness for the headline. γ_l is a weak
  second; the others have a trace but no robust home.

  TWO refinements to honour in prose: (i) θ is NOT structureless — it has a real,
  shaft-robust MTL concentration on a *net-null* band (a point FOR the method, NOT
  a trace — keep θ "the absent band"); (ii) γ_h is patient-specific by weak net
  trace, not by failing localization (its parietal lean is real but subset-only).
  One open compute: γ_l→PFC still owes a leave-one-patient-out before it is as
  hard as β. The *consistency* axis is ORTHOGONAL to the broadband/SNR
  *detectability* axis — β is special on both.
- **R1.2 "invisible to a connection-by-connection comparison" (multiscale beats
  raw).** This is being verified right now (Q2 = raw per-edge ρ vs cophenetic
  ρ_split on the same patients). If raw reproduces the band story, R1.2's
  method-superiority claim narrows or is partly retracted. **Mark R1.2's
  raw-comparison sentence `[under verification]` until Q2 lands.**
- **RESOLVED — "heterogeneity is detectability" is REVISED by Q1 (do not write the
  old strong version).** The full multi-phase reliability analysis
  (`.agents/reports/2026-06-25_multiphase-snr-reliability.md`) shows the
  heterogeneity is **NOT purely detectability**: after correcting `ρ_split` for the
  reliability of all four phase-distance vectors, reliability explains only **~6%**
  of the trace variance (β R²≈0.001), and the cross-patient "who-traces" axis is
  **orthogonal to reliability** (ρ=+0.08, flat). The clean falsifier: **Pat_15 is
  the HIGHEST-reliability patient yet the LOWEST tracer** — a genuine non-tracer,
  not a noisy one. (The earlier "+0.78 SNR↔ρ, all detectability" was inflated — it
  used a rest_pre-only noise model that was both partly tautological and blind to
  rest_post, the actual reliability bottleneck in 38/60 cells; corrected β SNR↔ρ
  is +0.49.) **So write:** within a band, recording reliability *modulates* whether
  the trace is detectable (detectability is one real factor); but a **residual
  biological "who-traces" axis survives** and is NOT noise. **That residual is
  UNEXPLAINED — do not give it a mechanism.** The natural guess (implant coverage
  of OFC) was TESTED and REFUTED (audit_147, 2026-06-25): Spearman(β ρ_split, OFC
  fraction) = −0.10 (flat), partial-controlling-reliability −0.02; the two
  strongest tracers (Pat_08, Pat_03) have ZERO OFC electrodes, the two with the
  most OFC trace middling/negative. So write the heterogeneity as **a real spread
  that is neither measurement noise nor a known covariate (not reliability, not
  coverage, not implant size)**; the cohort trace stands on its majority (β 7/10,
  p=0.005) regardless ([[feedback_fluctuations_are_signal]] — show the spread,
  don't force a single number OR a single explanation). This is CONSISTENT with the
  locked **"OFC = hotspot, not container"**: patients without OFC still show a
  strong brain-wide trace, so OFC is where the trace is densest WHEN sampled, not
  where it is required. Frame the cohort claim as "robust at the cohort level
  (β 7/10), with a real but unattributed between-patient spread," NOT "universal,
  detectability-limited" and NOT "coverage-gated."

## 3. The concrete R1 edits (this is the actual job)

Current R1 (`RESULTS_OUTLINE.md` / `01_trace.md`) frames the trace as β-specific.
Reframe:

- **R1 title.** Current: *"A held, multiscale β-band trace of learning and
  reasoning."* The trace is broadband, so "β-band" as the *defining* property is
  now wrong. Propose a title where **β is the cognitively-special band (forward
  ref to R2), not the only band with a trace** — e.g. *"A held, multiscale trace
  of learning — broadband in form, β-band in cognitive content."* (workshop it).
- **R1.1.** Current: *"…strongest and fully confirmed in β, while θ shows nothing
  — a selective signature, not generic drift."* → Rewrite to: the reorganization
  persists **across bands** (it is not a β-only effect); **θ is the single
  exception**; the per-patient differences are a **detectability (SNR) axis**, not
  biology — high-SNR recordings reveal the trace broadly, low-SNR ones cannot
  resolve it (so "non-tracers" are SNR-limited). The band specialness is about
  *content* (R2), introduced here and paid off in R2.
- **R1.2.** Keep the "not single links / survives strength-matched control"
  spine, but **flag the raw-comparison clause `[under verification]`** (Q2).
- **R1.3 (OFC), R1.4 (β rides healthy cortex / spares core), R1.5 (held non-
  ergodic state), R1.6 (α recruits core bridge).** Unchanged by this session.
- **Add (provisional) a sub-result for the consistency taxonomy** — the band
  narration recovered as *consistency*, not *presence* — but mark it
  `[provisional — pending per-band localization lock]` per §2.

Net effect: R1 stops claiming β-exclusivity of the bare trace, reframes the band
axis as (i) broadband presence with θ-exempt, (ii) detectability-graded per
patient, (iii) — provisionally — consistency-graded per band; and β's specialness
migrates to *cognitive content*, which is R2's flagship job.

## 4. Coordination boundary (hard rules)

- **DO** edit the draft headline files: `RESULTS_OUTLINE.md`, `01_trace.md`, and
  R1-related headline drafts. Propose new titles/sub-result phrasings.
- **DO NOT** edit `.agents/preprint/locked/VERDICT_LEDGER.md` or
  `.agents/preprint/locked/ANATOMY_LEDGER.md`. Those verdicts are cascaded by the
  main (compute) session AFTER Step 1 + Q1 + Q2 resolve; parallel edits there will
  collide. If you believe a ledger verdict needs changing, write the proposed
  change as a NOTE in the headline draft, not in the ledger.
- **DO NOT** invent tier assignments or assert "multiscale beats raw" as settled —
  see §2. Placeholders only for those two.

## 5. Source of truth

Read `.agents/reports/2026-06-25_snr-band-taxonomy-handoff.md` first (the full
decided/pending picture and the key numbers), then
`2026-06-25_per-node-trace-anatomy-and-heterogeneity.md` (the SNR/detectability
result the R1 reframe rests on). Everything in §1 is drawn from those; everything
in §2 is what the main session is actively locking.
