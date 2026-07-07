---
name: writing-directive-coph-substrate-results-figures
era: IMCOH_ABS_COHORT_N10
status: current
kind: directive
scope: Results-section figures, caption specs and narrative skeleton for the LRG cophenetic-distance per-pair trace subsection ("LRG-based cophenetic distance for the disentanglement of multiscale per-pair communication task-induced trace"). All-band cohort overview only (β primary, α secondary; per-band deep-dives land later). Pairs with the raw-substrate directive of the same date. Changes no locked verdict.
created: 2026-06-03
companion: directives/writing_directive_2026-06-03_raw-substrate-results-figures.md (the raw twin — its §0 magnitude-collapse caveat and Fig C caption are SHARED with this file); directives/writing_directive_2026-06-01_per-pair-probe-battery-restructure.md (methods restructure + canonical numbers); methods/methods_displacement_taxonomy_2026-06-03.md (parallel descriptive decomposition — distinct, do not merge)
verified_against: data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv (coph GATE); data/reports/imcoh_continuous_trace/{controls_summary.csv,controls_band_stats.md} (coph drift + cross-probe); data/audit/raw_fc_matched_strength/{cohort_summary_all_bands.csv,drift_floor_band_stats.md,drift_floor_per_patient.csv,cross_probe_per_patient.csv} (raw GATE + drift + cross-probe). Every number below was reproduced this session by running preprint_18 --both and preprint_21 --both.
---

# Cophenetic-substrate Results subsection — figures, captions, narrative skeleton

**Head.** The second Results subsection takes the same two cohort views built
for the raw substrate — a per-band joint-density portrait and a per-band
"null triangle" — and recomputes them on the **LRG cophenetic distance** `D`.
Where the raw substrate lit a positive diagonal in *every* band and failed its
drift control everywhere, the renormalised representation is **band-selective**:
the per-pair trace clears the full control battery (matched-strength gate +
within-session drift + cross-probe) at **β (primary)** and **α (secondary)**
only. The renormalization step does not merely lift significance — it
**relocates** it from the wrong bands the raw gate fired on (δ) to the
physiologically coherent α/β. That relocation is the "disentanglement" in the
subsection title. **No locked verdict changes**; this is figure + narrative for
an already-locked battery.

---

## 0. The arc the subsection must carry — and the four honesty guards

The subsection completes the "naive fails → renormalization succeeds" arc the
raw subsection opened. It should make exactly four moves, in order, and respect
four guards that are easy to get wrong.

**The arc**
1. The raw substrate closed unable to separate task-trace from session drift
   (and its matched-strength gate even fired on a non-trace band, δ, while
   *missing* β — see §3). 
2. Recompute the identical two views on the cophenetic distance `D`.
3. Joint density (Fig A-coph): the diagonal ridge survives the strength floor
   **only at α/β**; the other four bands fall to/below the floor (cream).
4. Null triangle (Fig C-coph): the same control battery now clears at α/β —
   the trace pulls clear of the drift floor and survives the cross-probe
   restriction — and the inversions that sank the raw substrate **collapse in
   magnitude** at the trace bands. Headline: **β primary, α secondary.**

**GUARD 1 — the verdict control is the matched-strength GATE, and it is α/β.**
"Mainly β, secondarily α" is sourced to the gate (`paired_wilcoxon_p` of the
matched-strength split-baseline surrogate): **β p=0.005** (obs +0.221, 23.7× the
surrogate median, 8/10 above floor), **α p=0.002** (obs +0.105, 8.3×, 8/10). β is
primary on magnitude *and* it is the only band carried by a second, independent
LRG probe (Grassmann); **α is cophenetic-distance–specific** (Grassmann does not
carry it). Say so: α is the secondary, cophenet-only trace.

**GUARD 2 — γ_l is NOT a third trace band. Do not let the drift figure say it
is.** γ_l clears the *drift* control (p=0.011) but **fails the gate (p=0.116)**
and sits **below its strength floor** in the joint density (obs +0.083 <
surr-p95 +0.109; it is cream in Fig A-coph). Its high obs-to-surrogate-median
ratio is misleading — the γ_l surrogate distribution is wide, so the obs does
not clear it. γ_l is "drift-cleared, gate-failed": mention as a marginal at
most, never as a headline trace alongside α/β. The clean read is **two trace
bands**.

**GUARD 3 — the raw→coph drift improvement is an inversion-MAGNITUDE collapse,
not a change in inversion count. Never "the red lines turn gray."** This is the
SHARED §0 of the raw twin directive
(`writing_directive_2026-06-03_raw-substrate-results-figures.md`); it governs
Fig C, which both subsections cite. Verified: the inverting patients largely
persist in *number* raw→coph (β stays 2/10; δ 3/10, θ 4/10, γ_h 5/10 unchanged);
what collapses is the summed inversion excess, band-selectively — α 0.69→0.04,
β 1.05→0.35, γ_l 0.48→0.13 — while δ/θ/γ_h barely move. The dissenting patients
still dissent, only mildly. Mechanism: the communication-distance representation
**strips within-session drift's leverage** specifically at the trace bands
(consistent with [[feedback_lrg_step_is_confirmation_not_resolution]] — LRG is
structural enrichment, not resolution of θ/γ).

**GUARD 4 — this is the all-band cohort overview. No anatomy, no per-patient,
no per-band deep-dives here.** Per the user's framing, per-band subsections land
later. In particular: do **not** add any anatomical localization in this (substrate)
subsection — anatomy lands in its own subsection. (⚠ the "spatially delocalized" claim
here is SUPERSEDED 2026-06-10: the multi-region DK lists stay retracted, but β concentrates
in the OFC *system*; α/γ_l/δ have no localization — see `ANATOMY_LEDGER.md` / `HANDOFF_INDEX.md`;
[[localization_audit_plan_2026_05_29]] / [[preprint_referee_review_2026_06_01]]).
Keep the claim at "a cohort-wide, band-selective per-pair trace at α/β."

---

## 1. Resolves the writing agent's two flagged concerns

**(BLOCKING) Verify the drift-null numbers.** Resolved — they are not
provisional. The cohort drift triangle has been run cohort-wide for both layers
and the numbers below were reproduced this session:

- Raw drift: `data/audit/raw_fc_matched_strength/drift_floor_band_stats.md`
  (status: current) + per-patient `drift_floor_per_patient.csv`. One-sided
  paired Wilcoxon `ρ_split^raw > drift`: δ 0.070 · θ 0.142 · α 0.142 · β 0.101 ·
  **γ_l 0.057** · γ_h 0.399 — no band < 0.05. Inversions (10 − n>floor): δ 3,
  θ 4, α 4, **β 2**, γ_l 3, γ_h 5 → "2–5/10". These are the exact numbers in the
  raw caption.
- Coph drift: `data/reports/imcoh_continuous_trace/controls_band_stats.md` +
  `controls_summary.csv`. Paired Wilcoxon `ρ_split^coph > ρ_drift^coph`:
  δ 0.222 · θ 0.254 · **α 0.0083 · β 0.0142 · γ_l 0.0109** · γ_h 0.399.
  *(Source caveat to carry: that file's header notes N=9 for the drift column —
  Pat_14 is not in the H2e run — while its paired table lists n=10. Immaterial
  to the verdict but flag it if a reviewer counts patients.)*

Point the writing agent at those two files; the section's drift claim is safe to
make load-bearing.

**(STRUCTURAL) "The paper holds the naive substrate to a stricter test (drift)
than the result it advances (gate)."** Resolved by symmetry, not by lowering the
bar. Both subsections run the **same three-control battery at both layers** (gate
+ drift + cross-probe; Fig B raw, Fig C coph are the same figure on the two
substrates). The honest framing is therefore:

> Same battery, applied identically to both representations. The raw substrate
> fails it two ways at once — it misses the drift floor in all six bands **and**
> its matched-strength gate fires on the wrong bands (δ yes, β no). The
> cophenetic representation passes the same battery cleanly at α/β: above the
> strength gate, above the drift floor, unchanged by the cross-probe restriction.

This is resolution (a) the writing agent proposed (run the drift triangle at the
coph layer and show it passes where the substrate failed) — and it is already
computed. The gate remains the verdict control everywhere; the drift triangle is
the supplementary "naive-fails / method-succeeds" visual. For α/β the two point
the same way, so there is no double standard.

---

## 2. Figure inventory (both under `data/preprint/figures/all_bands/`)

| fig | file | script | what it shows |
|---|---|---|---|
| **A-coph** | `fig_bands_joint_density_empirical_coph.pdf` | `preprint_18_bands_joint_density_rawfc.py --layer coph` | 2×3 cohort joint density of within-patient rank pairs on `D`. **Diagonal ridge clears the matched-strength floor only at β and α; δ/θ/γ_l/γ_h are cream.** Clean visual swap with the raw Fig A. |
| **C-coph** | `fig_bands_null_triangle_coph.pdf` | `preprint_21_bands_null_triangle.py --layer coph` | 1-row null triangle on `D`: boxes ρ_split/ρ_drift/ρ_xprobe + gray matched-strength floor; per-patient connectors. **Drift cleared at α/β/γ_l; inversions persist in number but collapse in magnitude at α/β/γ_l.** (Already inventoried in the raw twin directive — shared figure.) |

**Engineering note (clean swap).** `preprint_18` now takes `--layer {raw,coph}`
(mirroring `preprint_21`). Both layers run the *identical* panel recipe — same
`TRACE_CMAP`, same `SHARED_VABS=0.025`, same nested HDR contours `(0.10…0.99)`,
same conditional-split marginals — so the raw and cophenetic grids differ only in
substrate and the `Δ` superscript. The coph branch reads the canonical per-pair
cophenetic shift vectors from
`data/reports/imcoh_continuous_trace/per_pair_split/{Pat}_{band}.npz`
(`dD_task`, `dD_rest`) and the matched-strength floor from
`matched_strength_surrogate_split_baseline/cohort_summary.csv`
(`surr_median_rho_p95`).

**`preprint_09_bands_joint_density.py` is superseded for this purpose** — it
carries on-figure ρ text and significance-coloured borders (forbidden by the
no-text-in-figures convention) and is not a pixel-clean swap with the raw
figure. It is retained as history; the writing agent should reference
`preprint_18 --layer coph`. The `Δ` axis superscript on Fig A-coph is currently
`^{\mathrm{coph}}`; **align it with whatever symbol the methods section uses for
the cophenetic per-pair shift** — the figure regenerates in seconds if the
symbol changes.

---

## 3. Caption specs (the writing agent writes the final text)

**Fig A-coph — what the caption must contain** (draft below; tighten to house
style):
- It is the cophenetic twin of the raw joint-density figure — "as Fig.
  \ref{fig:rawfc_joint_density}, but on the LRG cophenetic distance `D`."
- Axes: `u = rank Δ_task(i,j)`, `v = rank Δ_rest(i,j)` on `D`, cohort-pooled,
  n=10, band order δ θ α / β γ_l γ_h.
- Colour field, contours, marginals identical to the raw figure (one sentence,
  do not re-explain the encoding in full — cross-reference).
- **The payload sentence:** where the raw substrate showed a green diagonal in
  every band, the renormalised representation keeps a diagonal above the
  matched-strength floor **only at β (strongest) and α**; δ/θ/γ_l/γ_h fall to or
  below the floor (cream). State that this band-selectivity is the property the
  edge-local substrate lacked.

> **Per-pair trace at the LRG cophenetic substrate, resolved by frequency band.**
> As Fig.~\ref{fig:rawfc_joint_density}, but the two within-patient rank-shift
> vectors are computed on the LRG cophenetic distance $D$ rather than the raw
> $|\mathrm{ImCoh}|$ adjacency: $u=\operatorname{rank}\Delta_{\mathrm{task}}(i,j)$,
> $v=\operatorname{rank}\Delta_{\mathrm{rest}}(i,j)$, cohort-pooled (top row
> $\delta,\theta,\alpha$; bottom row $\beta,\gamma_{\mathrm{l}},\gamma_{\mathrm{h}}$;
> $n=10$). Colour scale, highest-density-region contours and conditional-split
> marginals are identical to Fig.~\ref{fig:rawfc_joint_density}. Where the raw
> substrate showed a green (trace-direction) diagonal ridge in every band, the
> renormalised representation is band-selective: a diagonal clearing the
> matched-strength noise floor survives only at $\beta$ (strongest) and $\alpha$,
> while $\delta$, $\theta$, $\gamma_{\mathrm{l}}$ and $\gamma_{\mathrm{h}}$ fall to
> or below the floor (cream). The communication-geometry representation isolates
> the per-pair trace to the two bands the edge-local substrate could not
> separate.

**Fig C-coph — caption already approved.** Reproduce the "Fig C (cophenetic null
triangle)" caption verbatim from
`writing_directive_2026-06-03_raw-substrate-results-figures.md` §2. It already
encodes the magnitude-collapse framing (Guard 3) and the α/β drift-clearance.
Do not rewrite it from scratch.

---

## 4. Verified numbers (reproduced this session via the figure scripts)

Cophenetic layer, all three controls, with sources. **Bold = clears that
control.**

| band | ρ_obs (median) | obs / surr-median | n>0 | joint-density floor (obs vs p95) | **gate p (MS)** | drift p | inv | Σ-excess raw→coph |
|---|--:|--:|:--:|:--|--:|--:|:--:|:--|
| δ | +0.008 | ~0.8× | 6/10 | below (0.008 < 0.140) | 0.278 | 0.222 | 3/10 | 0.73 → 0.66 |
| θ | −0.040 | <0 | 3/10 | below | 0.722 | 0.254 | 4/10 | 1.05 → 0.58 |
| **α** | **+0.105** | **8.3×** | 8/10 | **above** (0.105 > 0.047) | **0.002** | **0.008** | 2/10 | 0.69 → **0.04** |
| **β** | **+0.221** | **23.7×** | 8/10 | **above** (0.221 > 0.076) | **0.005** | **0.014** | 2/10 | 1.05 → **0.35** |
| γ_l | +0.083 | 30×* | 7/10 | **below** (0.083 < 0.109) | 0.116 | 0.011 | 1/10 | 0.48 → 0.13 |
| γ_h | +0.000 | ~0× | 4/10 | below | 0.246 | 0.399 | 5/10 | 1.59 → 1.00 |

\* γ_l's large obs/median ratio is an artefact of a *wide* surrogate
distribution (high p95) — it is exactly why γ_l fails the gate despite the
ratio. Do not quote the 30× for γ_l as evidence of a trace (Guard 2).

**Sources for every number (so the writing agent can cite):**
- gate p, ρ_obs, surr-p95, surr-median, n-above:
  `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`
  (= audit_63; this is the headline obs/gate source per GAP-2).
- drift p, cross-probe medians, n>0:
  `data/reports/imcoh_continuous_trace/{controls_summary.csv,controls_band_stats.md}`.
- inversion count + Σ-excess: `preprint_21 --layer coph` stdout (recomputed from
  the per-patient `controls_summary.csv` triple).
- raw comparison column: `data/audit/raw_fc_matched_strength/…` (see §1).

**Headline-friendly framing of the two ×-ratios:** β is 23.7× and α 8.3× the
strength-surrogate median; β is also 8/10 above the per-patient p95 floor and α
8/10. (The "23×" that appears elsewhere in the briefs is this β obs/median ratio,
**not** an obs/p95 ratio — keep the two reference points distinct.)

---

## 5. Narrative skeleton (the agent writes the prose; this is the spine)

Three paragraphs, all-band cohort level, picking up exactly where the raw
subsection closed.

- **P1 — the representation that provides the separation (→ Fig A-coph).** Open
  on the raw subsection's unresolved tension ("…provides no further means to
  separate the two"). Introduce the cophenetic distance `D` as the renormalised
  communication-geometry representation, and recompute the identical per-pair
  rank statistic on it. The joint density (Fig A-coph) is now band-selective:
  the diagonal clears the matched-strength floor only at β and α; the other four
  bands return to the floor. The substrate's ubiquity has become a selection.

- **P2 — the controls confirm the selection, same battery (→ Fig C-coph).** The
  same three-control battery the raw substrate failed now clears at α/β: above
  the matched-strength gate (β p=0.005, 23.7×; α p=0.002, 8.3×), above the
  within-session drift floor (β p=0.014, α p=0.008), and unchanged by the
  cross-probe restriction (ρ_split ≈ ρ_xprobe). The mechanism behind the drift
  clearance is the **inversion-magnitude collapse** (Guard 3 — never "turns
  gray"): the dissenting patients persist in number but their inversions shrink
  band-selectively at the trace bands, because the communication-distance
  representation strips within-session drift's leverage there. β is primary (and
  the only band a second LRG probe corroborates); α is the secondary,
  cophenetic-specific trace. (γ_l clears drift but not the gate — Guard 2 — so it
  is a marginal, not a third band.)

- **P3 — what the renormalization buys ("disentanglement").** The edge-local
  substrate could not localise: its gate fired on δ and missed β. The cophenetic
  representation relocates the signal to the physiologically coherent α/β and
  removes it from δ/θ/γ. The renormalization step is not a significance boost on
  the same bands — it is a **disentanglement** of a genuine multiscale per-pair
  trace (α/β) from system-wide drift that the per-edge comparison conflated.
  Close by forward-pointing to the per-band subsections (β, then α) without
  making any per-band/anatomical claim here (Guard 4).

---

## 6. Minor concerns + wording fixes (from the writing agent)

- **`\includegraphics` prefix consistency.** The raw Fig 2 line is missing the
  `figures/` prefix that Fig 1 has (`graphicspath` covers it). Make the new
  cophenetic figures' `\includegraphics` use the same `figures/…` prefix as Fig 1
  for consistency.
- **"three network-based layers" vs "no network notion."** The opening sentence
  calls the substrate one of three "network-based" layers while the same
  paragraph (correctly) says the substrate has no global/network notion. Reframe
  the three layers by *increasing structural abstraction* rather than as
  uniformly "network-based": the **edge-local** substrate (per-pair, no
  integration) → the **communication-geometry** cophenetic layer (this
  subsection) → (the later layer). The substrate is the pre-network/edge-local
  rung, not a network layer; the cophenetic distance is where the network
  geometry first enters.

---

## 7. Placement + scope

- Fig A-coph + Fig C-coph + the three-paragraph skeleton → **the LRG cophenetic
  Results subsection** ("LRG-based cophenetic distance for the disentanglement
  of multiscale per-pair communication task-induced trace").
- **All-band cohort overview only.** Per-band subsections (β primary, α
  secondary) land later; defer Grassmann, per-patient and *all* anatomy to them
  — and remember anatomy localization is retracted (Guard 4).
- **No verdict change.** Per-band trace verdicts stand as locked
  (`locked/VERDICT_LEDGER.md`); matched-strength is the GATE, drift/cross-probe
  are supplementary controls (`locked/CONTROLS.md`). This file is figure +
  narrative for that already-locked battery.
- **Distinct from** `methods/methods_displacement_taxonomy_2026-06-03.md` (a
  per-pair *character* decomposition of ρ_split^coph) — that is a different
  thread; do not merge the two.
- **Sibling of** the raw twin directive of the same date; the two subsections are
  one arc (raw fails → coph disentangles) and share Fig C and its §0 caveat.
