---
name: results-readiness-map
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-04
supersedes_workflow_anchor: 2026-07-04_results-draft-review-handoff.md
pointers:
  - .agents/preprint/directives/writing_directive_2026-06-24_results-section-full-draft.tex   # THE artifact under review
  - .agents/preprint/headlines/{01_trace,02_encoding_vs_inference,03_epileptogenic_markers}.md # claim sources
  - .agents/preprint/locked/VERDICT_LEDGER.md                                                  # do NOT contradict
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv                    # C3 referee (verified 2026-07-04)
  - data/audit/grassmann_cluster_extent/cohort_summary.csv                                     # Grassmann (verified 2026-07-04)
  - data/audit/localization_atlas/per_band_taxonomy_verdict.csv                                # taxonomy (verified 2026-07-04)
---

> **Head.** This is the working spine for the Nature-Neuroscience Results review:
> every paragraph of the draft `.tex`, its load-bearing numbers (re-verified
> against the cached CSVs today), the figure that must support it, and the single
> sharpest critique a referee will make. We review **paragraph by paragraph** — PI
> brings the LaTeX rewrite, agent fixes/critiques/backs-with-numbers-and-lit,
> **reviewer-mode always on, brutal honesty, no sycophancy, no overclaim.** §3 lists
> the scientific tensions that need a PI framing decision *before* they bite. §5 is
> the running log of the PI's LaTeX commands/preferences — fill as we go.

## 0. Workflow (locked this session, 2026-07-04)

1. **Gather + arrange** (this file). Done for R1–R3 claim-set + R1 anchor numbers.
2. **Review the `.tex` paragraph by paragraph.** For each paragraph: (a) PI supplies
   their LaTeX rewrite with their macro/style preferences; (b) agent fixes prose +
   integrates missing pieces; (c) **re-verify every number against the cited CSV**;
   (d) gather + weave the supporting literature; (e) read the code that produced the
   result if the claim is contestable; (f) design/refine the figure together.
3. **Reviewer mode, always.** Every paragraph gets the hardest fair critique — the
   goal is a real high-impact paper, so weight every word, never overclaim, frame
   each claim to exactly what the null supports, support with numbers + figures.
4. **Capture PI LaTeX preferences** in §5 the moment they appear (memory too).

## 1. The three headlines

- **R1 — A held, multiscale trace of learning and reasoning.** The task reorganises
  resting connectivity and it does **not** wash out; clearest in β; consolidates to a
  consistent cortical home in two rhythms (β→OFC, γ_low→PFC); carried by healthy
  cortex; **held** as a sustained non-ergodic state, not bursty.
- **R2 — Offline abstraction of a learned structure (FLAGSHIP).** The persistent
  trace carries what was **inferred**, not just what was **seen** (β-only,
  inference-specific persistence); mesoscale-favouring; encoding anchored in OFC set
  by learning; α = memory-only, β = memory+inference; low-γ hides a focal cingulate
  memory trace.
- **R3 — Propagator-inspired markers of epileptogenic tissue.** Same operator;
  seizure tissue is relational not contact-by-contact; seed-based 6-band calibrated
  P(SOZ) detector (AUC≈0.81); two-population split; clinical-label not outcome.

## 2. Per-paragraph spine

Legend: ✓ = number re-verified against CSV today; ○ = to re-verify at that paragraph.
"Flag" = the sharpest fair referee critique.

### R1 — held multiscale trace (`.tex` lines 123–258)

| ¶ | Claim | Load-bearing numbers | Figure | Referee flag |
|---|---|---|---|---|
| R1.1 | Reorganisation persists into rest, cohort-level | β +0.22, 7/10, p=0.005 ✓; α p=0.002 ✓; θ null ✓ | Fig 2 forest/gate | "clearest in β" must not read as "most significant" — **α's p (0.002) is smaller than β's (0.005)**; β leads on 7/10 + localizes + both probes, not on p. Word it so. |
| R1.1 | Spread is biology, not reliability | Q1 residual axis survives reliability correction (audit_145) | — (Methods) | SNR is a **control** → one clause only; Pat_15 most-reliable-lowest-tracer is the proof, keep in Methods. |
| R1.2 | β→OFC consolidation | R=1000, BH q≈0.009–0.013, LOO 10/10 top, low-strength/non-hub, bilateral ○ | Fig 3 OFC system brain | OFC **n=5 coverage**; verdict flipped 3× historically → say "concentration above baseline", never "lives in OFC"; washes out at lobe (p≈0.29)/hemi. |
| R1.2 | γ_low→PFC (2nd rhythm) | band gate p=0.116, 5/10 (FAILS C3) ✓; **PFC LOO fails q≤0.05 in 6/10 drops ✓ (carrier_loo.csv)** | Fig 3 | **THE tension (§3.1) — worse than a borderline gate. DOWNGRADED by 2026-06-26 vision doc; draft's "sig every fold, top 8/10" is a misread of carrier_p/is_top, missing carrier_q. Not a consistent home.** |
| R1.3 | Independent read-out confirms β | Grassmann β mass-p=0.005, LOO 0.005 ✓ | Fig 1 band×probe | Grassmann = **triangle baseline**, cophenetic = split-baseline → **not baseline-matched** (§3.4); frame as "confirms at β", the convergence band. |
| R1.3 | Genuinely multiscale | β cophenetic p=0.005 vs matched raw-edge p=0.053 ✓ | Fig 1 3-layer strip | **Do NOT say "invisible to edge-wise"** — raw magnitude is larger, Q2 forbids (§3.3). Allowed: clears MS where raw falls just short. |
| R1.4 | β rides healthy cortex, spares core | gray↔gray hotspot p<0.001 ○; epi-core p_pair≈0.945 ○ | Fig 3 pair-class | p_pair≈0.945 is **non-significant, not significant depletion** → "core spared" is directional only; state both withdrawals (node-count, interface). |
| R1.5 | Held state, not bursty | sustained 10/10 p=0.001, LOO 0.002 ○; 5 read-outs +8–10/10 p≤0.007 ○; transient = clean cohort negative ○ | Fig 2 held-state timecourse | "**held, not tighter**" — tighter-attractor is only p=0.08. Solid otherwise. |
| R1.6 | Bridge: β spares, α recruits core | α epi↔epi clears C5 ○ | — (bridge) | C5 is a **secondary** observation; α has **no** anatomical home → keep as bridge, don't over-weight. |

Table 1 (band×probe): δ Grassmann shown ✓ but ledger says **weak + LOO-fragile** (§3.2).

### R2 — offline abstraction / flagship (`.tex` lines 261–420)

| ¶ | Claim | Load-bearing numbers | Figure | Referee flag |
|---|---|---|---|---|
| R2.1 | Trace carries what was inferred | inference-specific β p≈0.007, LOO-robust, clears MS ○; arc reproduces N1 to ~1e-16 ○ | Fig 4 4-phase + arc | **The paper's central claim.** No behaviour → "abstraction" is interpretation of a *persistence* result, never tied to inference success. State ceiling in-paragraph. |
| R2.2 | Mesoscale of multi-step integration | τ≈2.6 p=0.014, τ≈6.8 p=0.010 ○; obs−null gap +0.11→+0.15 ○; controls null ○ | Fig 4 mesoscale | "favouring, not exclusive" — sig across scales. Mesoscale↔multi-step-paths is interpretation; label it. |
| R2.3/2.6 | OFC anchor set by learning | T_learn α+0.22/β+0.25 p=0.014 ○; learn-β→OFC R=1000 q=0.010 ○ | Fig 4 double-dissoc | Solid. Smooth the **OFC seam** with R1.2 (same hotspot, overall vs encoding component). |
| R2.4 | Memory α / inference β; duration-clean | β length-ratio ρ≈+0.25 n.s.; α ρ≈+0.53 (suspect) ○; truncation null retired ○ | Fig 4 band×content | One methods sentence: whole-brain concordance corruptible by contiguous truncation; per-system demeaned is not. |
| R2.5 | Low-γ hides focal cingulate memory | ρ≈+0.39, 8/8, q≈0.035; epi-excl ρ≈+0.51 7/8 ○ | Fig 4 within-region | Cleanest "averaging hides it." Power floor: only well-sampled systems testable — cingulate is; don't generalise the absolute test to OFC. |
| R2.6 | Inference leans cingulate | **HINT** — fails duration, p≈0.12–0.18, q≈0.42–0.76 ○ | Fig 4 (caption) | Directional hint only; never an established location. |
| R2.7 | Cingulate multiplexes by rhythm | memory low-γ (solid) + inference β (hint) | — | Suggestive; the inference half is a hint → say so. |

### R3 — epileptogenic markers (`.tex` lines 423–511)

| ¶ | Claim | Load-bearing numbers | Figure | Referee flag |
|---|---|---|---|---|
| R3.1 | Relational, not contact-by-contact | δ off-shaft LOSO, label-shuffle p≈0 ○ | Fig 5 community schematic | Node-by-node reduces to hubness+depth; all-contacts version withdrawn (proximity tautology) → honest metric = off-shaft. |
| R3.2 | 6-band calibrated P(SOZ) | AUC≈0.81, prec@5≈60% (~7×), 9/10 ○ | Fig 5 calibration/strip | **Label-leak history (0.97→0.81)** — keep the fix + shuffle-null; it is the credibility. Interpretable logistic over GBM at n=10. |
| R3.3 | Seed-free = open frontier | (under investigation) | — | Report no verdict. |
| R3.4 | Two populations | strong-community AUC 0.90–0.99; hub (Pat_10/15) marginal ○ | Fig 5 per-patient strip | Cohort mean **can lie** → report the split, never the mean alone. |
| R3.5 | Honest scope | clinical-label not outcome; ~41% SOZ are WM ○ | — | Triage/shortlist, not standalone localizer; precision ceiling structural. |

## 3. Critical scientific tensions (need a PI framing call)

**3.1 — γ_low→PFC is NOT a consistent home; the draft overclaims it. [most severe]**
The draft (R1.2, Table 1, abstract-adjacent) presents γ_low→PFC as the second
consistent localized rhythm, co-headline with β→OFC. **Two independent failures say it
is not** (both verified against CSV today):
- **Cohort gate FAILS C3:** γ_low cophenetic `p=0.116, 5/10, median +0.083`
  (`matched_strength_..._split_baseline/cohort_summary.csv`) — the per-pair trace does
  not clear the mandatory matched-strength referee.
- **PFC localization FAILS LOO:** in `carrier_loo.csv`, PFC clears `q≤0.05` in only
  5/10 single-patient drops; **dropping Pat_05/06/07/10/14/15 pushes q>0.05** (worst
  −Pat_07 q=0.20). β→OFC, by contrast, clears `q≤0.05` in every fold, top-ranked in
  every fold.
- **Superseded:** `2026-06-26_per-band-phenomenology-vision.md` (frontmatter
  `supersedes` the 2026-06-25 taxonomy lock) already downgraded γ_low from *consistent*
  → **"strong 5-patient subset, PFC home single-patient-fragile, not cohort-locked."**
  The biology-led R1 (rewritten 2026-07-01) did not incorporate this — its
  "significant in every LOO fold, top-ranked in 8/10" was read off `carrier_p` /
  `carrier_is_top`, missing the BH-corrected `carrier_q`.

So the honest current state is: **β→OFC is the single cohort-consistent, LOO-robust
localized home. α is a real cohort trace (p=0.002) with NO anatomical home
(delocalised). γ_low neither clears the cohort gate nor holds its PFC localization
under LOO.** This collapses R1.2's "two consistent rhythms" frame. Options:
(a) drop γ_low→PFC from the R1 headline; R1.2 becomes β→OFC (the one locked home) with
α-has-no-address as the deliberate contrast — cleaner, fully defensible, but "band-
specific" now rests on β alone (+ the α/θ contrasts + the R2 β-only inference);
(b) keep γ_low but honestly demoted to "a strong-tracer subset shows a PFC lean, not
cohort-locked" (an aside, not a co-home); (c) other.

**→ RESOLVED (PI, 2026-07-04): option (b) — demote γ_low→PFC to a one-line honest
aside.** β→OFC is the single locked cortical home; α = cohort trace, no address;
γ_low = a one-line aside, not a co-home. Apply everywhere:
- **R1.2:** lead on β→OFC; α-has-no-address as the contrast; γ_low → one aside line.
  Candidate: *"In low-γ a prefrontal concentration appears within a strong-tracing
  subgroup, but it clears neither the cohort trace gate nor leave-one-out, so we do not
  treat it as a consistent cortical home."*
- **Table 1:** low-γ cophenetic cell `(✓)` → `--` (fails C3 gate); low-γ stays
  Grassmann-`✓` (Grassmann-only, per ledger). Caption: delete the "borderline gate…
  consistent localisation (low-γ→PFC)" sentence (now false); the parenthesised-check
  device is no longer needed.
- **Caption/abstract:** abstract already localizes β-only (no γ_low claim) — no change
  needed there. Discussion "two rhythms" language, if any, → "β" only.

**3.2 — Table 1 gives δ a full Grassmann ✓, but δ is "weak + LOO-fragile."**
Ledger: δ Grassmann cohort gate clears only at the R=200 floor (mass-p=0.005) and
**LOO fails** (0.055, Pat_08); Decision-12 tags it **weak**. The draft's full checkmark
over-represents it — and δ sits on the un-baseline-matched triangle arm (§3.4), so it's
doubly soft. Fix: parenthesise/footnote δ as weak+LOO-fragile, or drop the cell. Not a
verdict change (ledger stays) — a presentation-honesty fix.

**3.3 — Abstract + Discussion say "invisible to edge-wise": Q2 forbids it.**
`.tex` abstract (l.66–67) and Discussion (l.528) claim the reorganisation is "invisible
to edge-wise" comparison. Q2/audit_146: **raw magnitude is larger**; the ultrametric
*compresses*. Allowed phrasing only: "at β the hierarchy clears its strength-matched
control where the matched raw-edge comparison falls just short (p=0.053)" and "the band
taxonomy is a hierarchy property." Queued fix (determined by the lock, not a PI fork).

**3.4 — Grassmann baseline asymmetry (parked).** Grassmann = triangle
`T_G=d(pre_A,task)−d(task,post)` (single rest_pre); cophenetic = split-baseline ρ_split.
**Not baseline-matched** → the δ Grassmann-only cell is the exposed one; β convergence is
safe (β clears the rigorous cophenetic independently). Fix = one Methods caveat
paragraph (or a split-baseline Grassmann compute follow-up). Not on R1's critical path.

**3.5 — Spectral-superiority (N1.2b) fate — PI has not green-lit.** Abstract currently
leads on "invisible to spectral-clustering/PCA, which a multiscale read-out recovers."
Honest limit (audit_143): **β is a tie** (spectral recovers β); only the **α** component
is spectral-invisible + selectivity (cophenetic fires on α,β; textbook embedding
over-detects into null-θ). Keep as an abstract/Results highlight, or demote to Methods?
→ PI decision (affects abstract + Discussion methodological-lesson ¶).

**3.6 — The multiscale / higher-order-propagator framing (PI: fundamental; develop it). [plan]**
Placement decided 2026-07-04 — PI wants this surfaced as (a)+(b)+(c):
- **(b) R1.1 — DONE.** One foreshadow clause added: the trace is "carried not by the connectivity's
  strongest links but by the higher-order, multiscale organisation the diffusion hierarchy exposes,
  as developed below" (Q2-safe: strength-independence = C3; higher-order = propagator).
- **(a) R1.3 — TODO when we draft it.** Develop fully as a *property of the finding*: the
  \(|\mathrm{ImCoh}|\) graph is essentially complete/densely-weighted (no combinatorial structure —
  the structure is in the weights); the propagator \(e^{-\tau\hat L}\) retains the trace in its
  higher-order, multi-step terms (\(\hat L^2,\hat L^3,\dots\)) that the hierarchy reads out; at β the
  hierarchy clears matched-strength (p=0.005) where the raw weighted edges give only a diffuse,
  non-specific signal (p=0.053). **The \(\hat L^n\)-expansion mechanism stays in Methods.**
- **(c) Discussion — TODO.** Headline methodological point: cognition emerges from the higher-order
  propagation structure of a densely-weighted network, which single-scale / edge reads cannot resolve.
- **Q2 GUARD (everywhere):** never "invisible to the edges" — raw magnitude is LARGER but diffuse;
  say "emergent / higher-order / refined by the hierarchy."

**Citations:** R1.1 refs swapped to **electrophysiology substrate** (PI: Tambini/Liu2022 read fMRI):
**liu2019human** (Cell, MEG — offline reorg of learned structure) + **duan2025awake** (Prog. Neurobiol.,
iEEG — 2025, our substrate); both DOI-verified 2026-07-04. tambini2010enhanced + liu2022decoding kept in `references.bib`
(verified) for possible Discussion/R2 use. **Rule locked:** one authoritative + one recent, both
substrate-appropriate (MEG/EEG/iEEG for this paper), DOI-verified by opening the source.

## 4. Figure inventory (5 figures — PDF vector, use_lrg_style, no PNG, no C4 fig)

- **Fig 1 (method + band resolution).** Diffusion hierarchy + two probes; 3-layer
  band-resolution strip (raw FC → raw D(τ) → cophenetic, β resolved only at cophenetic);
  band×probe dissociation (Table 1). Source: raw_vs_multiscale (audit_146), verdict CSVs.
- **Fig 2 (the trace).** Cohort β forest/gate per band + X-epi variant; held-state
  timecourse (post vs pre rest) + transient-negative panel. Source: matched_strength +
  replay_states. **Show per-patient spread** (locked principle — the fluctuation is signal).
- **Fig 3 (anatomy + tissue).** OFC system-scale on bilateral brain ("hotspot not
  container" in caption); gray↔gray vs WM↔WM vs seizure-core pair-class. Source:
  localization_atlas + epi_stratified.
- **Fig 4 (flagship: abstraction).** 4-phase schematic (encoding vs inference);
  inference-specific β persistence + mesoscale; encoding→OFC vs inference→cingulate
  double-dissoc (inference "duration-downgraded" caption); low-γ within-region cingulate
  panel. Source: consolidation_arc + inference_localization.
- **Fig 5 (clinic).** Community schematic (seeds → distant SOZ, different shaft);
  per-patient detector strip (two-population); calibration + precision@k. Source:
  epi_propagator_detector.

## 5. PI LaTeX preferences (running log — locked 2026-07-04, append as more appear)

- **Output location:** each paragraph is its own numbered `.tex` file in
  `.agents/preprint/directives/results_paragraphs/` (e.g. `R1.1_reasoning-imprint.tex`).
  **Never paste the paragraph `.tex` body into chat** — write it to the file; discuss in chat.
- **No manual line-wrapping** (locked 2026-07-04): write each paragraph as a single unwrapped
  line and let the editor soft-wrap. Do not insert hard newlines mid-paragraph.
- **American English** (locked 2026-07-04): organization, reorganization, behavior, analyze,
  characterize, favor, center, inquiry — never the British `-isation`/`-our`/`-re` variants.
- **Register:** result-forward, high-level. **Minimise "we did X / we compared / we asked"**
  procedural first-person — it reads pedestrian. State the finding, not the action.
- **Math:** inline math with `\(...\)` (not `$...$`). Use the macros in
  `.agents/references/math_commands.tex`. Key ones: `\rhocoph` (= \(\rho^{\rm coph}\) —
  **this is the measure; "ρ_split" is not used**), `\lowgamma` (\(\gamma_{\rm l}\)),
  `\highgamma` (\(\gamma_{\rm h}\)), `\imcoh`, `\Dcoph`, `\pmassemp`. If a frequently-used
  symbol is missing, **propose adding it to `math_commands.tex`** rather than inlining.
- **Acronyms / phases:** use `.agents/references/acronyms.tex`. Phases via **`\acrshort`**
  (short form everywhere — locked 2026-07-04): `\acrshort{rspre}` (rsPre),
  `\acrshort{taskl}` (taskLearn), `\acrshort{taskt}` (taskTest), `\acrshort{rspost}`
  (rsPost). Patients in text: `\patient{15}` → "Pat15".
- **Bands: keep the normal Greek symbols** (locked 2026-07-04) — `\(\alpha\)`, `\(\beta\)`,
  `\(\delta\)`, `\(\theta\)`, `\lowgamma`, `\highgamma`. Do **not** introduce new band macros.
- **Citations:** per claim, prefer **one authoritative older + one recent (field-active)**
  reference, **substrate-appropriate** (MEG/EEG/iEEG for this paper — not fMRI-only). **Never
  invent.** Give the DOI and **open each source to verify** before use; accumulate verified
  entries in `results_paragraphs/references.bib`.
- **Cite-key format (locked 2026-07-04):** `<firstauthorsurname><year><first non-trivial title
  word>`, all lowercase — e.g. `liu2019human`, `duan2025awake`, `tambini2010enhanced`,
  `liu2022decoding`. Skip leading articles/prepositions (the, on, a, an, of, in, …).
- **Avoid defensive hedging** — e.g. "we make no claim that every patient expresses it" is
  banned; state the positive (cohort-level property, established by the group test).
- **Band-selectivity framing:** foreground the pronounced between-patient heterogeneity, then
  β as the band whose signal is strong enough to rise above it and clear the controls.

## 6. Verified-numbers ledger (re-checked against CSV, 2026-07-04)

Cophenetic C3 (`matched_strength_surrogate_split_baseline/cohort_summary.csv`), fmt
band: obs_median_rho / paired_wilcoxon_p / n_above_surrogate:
- β: +0.2206 / 0.00488 / 7 ✓ (clears)
- α: +0.1054 / 0.00195 / 5 ✓ (clears)
- low_γ: +0.0830 / 0.1162 / 5 ✗ (**fails gate — §3.1**)
- δ: +0.0076 / 0.2783 / 4 ✗
- θ: −0.0399 / 0.7217 / 2 ✗ (anti)
- high_γ: +0.0003 / 0.2461 / 4 ✗

Grassmann (`grassmann_cluster_extent/cohort_summary.csv`), band: cluster_p_mass / LOO_max:
- β: 0.005 / 0.005 (strong) ✓
- low_γ: 0.005 / 0.040 (strong) ✓
- δ: 0.005 / **0.055** (weak, LOO-fragile — §3.2) ✓
- high_γ: 0.060 (no trace) ✓
- θ: 0.159 (no trace) ✓
- α: 0.348 (no trace) ✓

Taxonomy (`per_band_taxonomy_verdict.csv`): β consistent→OFC; low_γ tagged consistent→PFC
in this CSV (gate p=0.116) **but SUPERSEDED** — see below; α patient-specific (cohort);
δ patient-specific (subset); θ absent; high_γ patient-specific (subset). θ & high_γ carry
`loc_without_cohort_trace=True` (MTL / parietal spatial structure on a net-null band —
**not** a trace, do not report as one).

Carrier LOO (`carrier_loo.csv`, `carrier_q` = BH-corrected, the locked bar):
- **β→OFC:** q≤0.05 in **all 11 rows** (FULL + 10 drops), top-ranked every fold. Locked. ✓
- **low_γ→PFC:** q≤0.05 only in FULL, drop-Pat_02/03/08/13 (5/10 drops); **q>0.05 in
  drop-Pat_05/06/07/10/14/15 (6/10).** Not LOO-robust. ✓ → §3.1.

**Supersession:** `2026-06-26_per-band-phenomenology-vision.md` (status current;
`supersedes` the 2026-06-25 taxonomy lock) is the LATEST synthesis. It downgrades
γ_low: *consistent → strong-subset, PFC single-patient-fragile*; keeps **β the only band
strong on BOTH axes (cohort gate + LOO-stable home).** Also corrects heterogeneity
framing: spread is **residual biology, unexplained** (split-half reliability explains ~0%
of β; Pat_15 most-reliable-lowest-tracer; OFC coverage refuted ρ=−0.10) — do **not** call
it a "detectability/SNR axis." The biology-led R1 (2026-07-01) predates incorporating this.

**Caveat — do NOT claim without the number:** the *cophenetic cohort-gate* LOO (does the β
paired-Wilcoxon survive dropping each patient?) is NOT in hand. β→OFC *localization* is LOO-robust
(carrier_loo.csv), but that is a different test. Do not assert single-patient robustness of the
*trace gate* until its LOO is verified. (Removed a vague "does not rest on any single patient" line
from R1.1 for this reason, 2026-07-04.)

**Deferred to per-paragraph verification:** all R2 arc numbers (p≈0.007 inference,
T_learn α+0.22/β+0.25 p=0.014, mesoscale τ nulls, OFC R=1000 q, PFC LOO, low-γ cingulate
8/8), all R3 detector numbers (AUC≈0.81, prec@5≈60%, 9/10), held-state 10/10 p=0.001.
