---
name: writing-directive-raw-substrate-results-figures
era: IMCOH_ABS_COHORT_N10
status: current
kind: directive
scope: Results-section figures, captions and prose for the raw |ImCoh| substrate per-pair trace (joint density + null triangle) and the cophenetic null-triangle twin, plus the VERIFIED raw→LRG inversion-magnitude-collapse finding. Pairs with the 2026-06-01 per-pair-probe methods-restructure directive. Changes no locked verdict.
created: 2026-06-03
companion: directives/writing_directive_2026-06-01_per-pair-probe-battery-restructure.md (methods restructure + GAP-closure canonical numbers); methods/methods_displacement_taxonomy_2026-06-03.md (parallel descriptive decomposition — distinct)
verified_against: data/audit/raw_fc_matched_strength/{cross_probe_per_patient.csv,drift_floor_per_patient.csv,per_patient_per_band_all_bands.csv,cohort_summary_all_bands.csv}; data/reports/imcoh_continuous_trace/controls_summary.csv; data/audit/matched_strength_surrogate_split_baseline/{per_patient_per_band.csv,cohort_summary.csv}. All per-band counts/p-values recomputed this session.
---

# Raw-substrate Results figures + the raw→LRG magnitude-collapse finding

**Head.** The Results section opens at the raw `|ImCoh|` substrate with two
cohort figures — a per-band joint-density portrait and a per-band "null
triangle" — that establish a per-pair trace signal which is **non-specific**
(every band) and **not separable from within-session drift** at the cohort
level. The LRG cophenetic twin of the null triangle then shows the trace
clearing its drift null at α/β/γ_l. This file records the figures, the
approved captions and prose, and — most importantly — the **verified mechanism
of the raw→LRG difference**, which is easy to mis-state.

## 0. CRITICAL caveat for the writing — do NOT mis-state this

The raw→cophenetic improvement is an **inversion-magnitude collapse at the
trace bands, NOT a change in inversion count, and NOT "the red lines turn
gray."** Verified this session:

- The patients who *invert* (their drift null exceeds their observed signal,
  `ρ_drift > ρ_split`) **persist in number** from raw to cophenetic — β stays
  2/10, δ/θ/γ_h are identical, only α (4→2) and γ_l (3→1) lose a couple.
- What changes is the **magnitude** of those inversions, **band-selectively**:
  the summed inversion excess collapses at α/β/γ_l and barely moves at
  δ/θ/γ_h (table in §4). The dissenting patients still dissent — only mildly —
  so the rank-based paired test clears.
- Mechanistic framing (consistent with "LRG = structural enrichment", see
  `feedback_lrg_step_is_confirmation_not_resolution`): the communication-
  distance representation **strips the within-session drift's leverage**
  specifically at the trace bands. Never write "the inversions disappear" or
  "the cohort becomes unanimous".

## 1. Figure inventory (all under `data/preprint/figures/all_bands/`)

| fig | file | script | what it shows |
|---|---|---|---|
| **A** | `fig_bands_joint_density_empirical_rawfc.pdf` | `scripts/02_preprint/preprint_18_bands_joint_density_rawfc.py` | 2×3 cohort joint density of within-patient rank pairs (rank Δ_task^raw, rank Δ_rest^raw), signed enrichment over the matched-strength floor + conditional-split marginals. **Every band shows a positive diagonal ridge → non-specific.** |
| **B** | `fig_bands_null_triangle_rawfc.pdf` | `preprint_21_bands_null_triangle.py --layer raw` | 1-row null triangle: boxes ρ_split^raw / ρ_drift^raw / ρ_xprobe^raw + gray matched-strength floor; per-patient connectors red (magnitude-encoded) on inversion. **Signal > strength floor (median) and ≈ cross-probe, but fails drift in all 6 bands.** |
| **C** | `fig_bands_null_triangle_coph.pdf` | `preprint_21_bands_null_triangle.py --layer coph` | Same on LRG cophenetic D, identical axes/encoding. **Inversions persist in number but collapse in magnitude at α/β/γ_l → drift test clears there.** |

Fallbacks (kept, not for publication): `fig_bands_controls_whisker_rawfc.pdf`
(`preprint_19`, three-distribution whisker) and
`fig_bands_signal_vs_noise_rawfc.pdf` (`preprint_20`, signal-vs-zone). They
were judged less immediate; retained only as alternates.

Connector encoding (B & C, shared scale `MARGIN_REF=0.6`): red when
`ρ_drift > ρ_split`, thickness + saturation ∝ `ρ_drift − ρ_split`; thin grey
otherwise. Per-patient x-offset is shared by a patient's connector and its
three markers (alignment fix). Boxes `showfliers=False` (outliers are the
per-patient circles; default fliers double-plotted them).

## 2. Approved captions

**Fig A (raw joint density)** — see the locked caption text in
§ "Fig A caption" below (already drafted/approved 2026-06-01; reproduced for
the writing agent). [Joint-density caption: cohort-pooled within-patient rank
pairs; green diagonal = trace direction; a diagonal ridge in every band
illustrates the substrate registers coupling throughout the spectrum but does
not single out a band-specific reorganization.]

**Fig B (raw null triangle).**
> **The raw per-pair trace against its nulls, by band.** Each panel (one row:
> δ θ α then β γ_l γ_h) summarises the cohort (n=10) over three boxes: the
> observed coupling ρ_split^raw, the within-session *drift null* ρ_drift^raw
> (the identical statistic with no task, aligning two rest-only difference
> maps, so pure session drift alone yields a positive value), and the
> *cross-probe* restriction ρ_xprobe^raw (ρ_split^raw recomputed on contact
> pairs spanning different probes). The gray band is the *matched-strength*
> surrogate floor (cohort p5–p95 of the strength-preserving null). Each
> patient contributes one value per box, joined by a connector that is red
> when its drift null exceeds its signal (an *inversion*), thickness +
> saturation scaling with the excess, thin grey otherwise. In every band the
> signal sits above the strength floor in the cohort median and is essentially
> unchanged by the cross-probe restriction (ρ_split^raw ≈ ρ_xprobe^raw); the
> binding control is the drift null, which it fails in all six bands — thick
> red inversions throughout (2–5 of 10 patients per band) and the one-sided
> paired test ρ_split^raw > ρ_drift^raw clears no band (best γ_l, p=0.057).

**Fig C (cophenetic null triangle).**
> **The same per-pair trace and nulls on the renormalised communication
> geometry.** As Fig. B, but every quantity is computed on the LRG cophenetic
> distance D rather than the raw |ImCoh| adjacency; axes, boxes and connector
> encoding identical. The inverting patients largely persist in *number*, but
> their inversions *collapse in magnitude* at α, β and γ_l (summed inversion
> excess 0.69→0.04, 1.05→0.35, 0.48→0.13 relative to Fig. B): the signal and
> cross-probe boxes pull clear of the drift box and the paired test
> ρ_split > ρ_drift clears at α (p=0.008), β (p=0.014) and γ_l (p=0.011). The
> non-trace bands δ, θ, γ_h retain large inversions and do not separate. The
> change is the within-session drift losing its leverage specifically at the
> trace bands — a property of the communication-distance representation, not a
> reduction in the number of dissenting patients.

## 3. Approved Results prose (raw substrate subsection)

**Paragraph 1 — after the substrate-intro paragraph, introducing Fig A**
(joint density): the cohort hallmark of ρ_split^raw > 0 reads as rank mass on
the diagonal; a diagonal concentration is present in *every* band, so the raw
per-pair statistic registers a coherent displacement everywhere and cannot, on
its own, separate a band-specific reorganization from system-wide drift.
(Full LaTeX already delivered 2026-06-01; reproduce verbatim.)

**Paragraph 2 — introducing Fig B (raw null triangle), raw-only, no coph:**
> This ubiquity is deceptive: reading the substrate critically requires putting
> ρ_split^raw on the same footing as its nulls (Fig. B). We compare it against
> a within-session *drift null* ρ_drift^raw, which applies the identical
> correlation to two task-free rest-only difference maps and so measures the
> alignment slow session drift produces on its own, and a *cross-probe*
> restriction across different probes; a matched-strength surrogate sets the
> level reachable from electrode-strength structure alone. Across every band
> the raw signal sits above the matched-strength floor in the cohort median
> and is essentially unchanged by the cross-probe restriction
> (ρ_split^raw ≈ ρ_xprobe^raw), so neither the strength structure nor
> same-probe short-range coupling is what carries it; the binding control is
> the within-session drift null, which it fails in all six bands. In each band
> two to five of the ten patients *invert* — their drift null exceeds their
> observed coupling (red connectors) — and the paired test
> ρ_split^raw > ρ_drift^raw reaches significance in none (best γ_l, p=0.057).
> At the substrate, then, the positive diagonal of Fig. A is not
> distinguishable from within-session drift: present in most patients but
> undercut by a high-leverage minority in whom drift alone reorganises
> connectivity at least as strongly as the task, and the raw per-pair
> comparison provides no further means to separate the two.

**Paragraph 3 — the cophenetic resolution (Fig C)** lives in the LRG
subsection, NOT the raw one (user directive 2026-06-03: the raw paragraph
speaks only about raw). To be drafted there using the §0 magnitude-collapse
framing — never "turns gray".

## 4. Verified numbers (recomputed 2026-06-03)

Per band, raw vs cophenetic: inversion count, summed inversion excess
(Σ over inverters of ρ_drift − ρ_split), drift-null paired-Wilcoxon p,
matched-strength gate p.

| band | raw inv | raw Σexcess | raw drift p | raw gate p | coph inv | coph Σexcess | coph drift p | coph gate p |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| δ | 3/10 | +0.73 | 0.070 | 0.042 | 3/10 | +0.66 | 0.222 | 0.278 |
| θ | 4/10 | +1.05 | 0.142 | 0.138 | 4/10 | +0.58 | 0.254 | 0.722 |
| α | 4/10 | +0.69 | 0.142 | **0.014** | 2/10 | **+0.04** | **0.008** | **0.002** |
| β | 2/10 | +1.05 | 0.101 | 0.053 | **2/10** | **+0.35** | **0.014** | **0.005** |
| γ_l | 3/10 | +0.48 | 0.057 | 0.053 | 1/10 | **+0.13** | **0.011** | 0.116 |
| γ_h | 5/10 | +1.59 | 0.399 | 0.188 | 5/10 | +1.00 | 0.399 | 0.246 |

Reading: raw clears the drift null in **no** band; cophenetic clears α/β/γ_l.
The discriminator is Σexcess (magnitude), not inv count — β is 2/10 in both
yet Σexcess falls 1.05→0.35. δ/θ/γ_h Σexcess barely moves → those bands stay
unseparated at both layers.

NB cophenetic ρ_split source: the per-patient triple (split/drift/xprobe) in
Fig C is from `imcoh_continuous_trace/controls_summary.csv` (internally
consistent with its own drift/xprobe). The **headline obs number** elsewhere
remains audit_63 (`matched_strength_surrogate_split_baseline`) per GAP-2; the
~0.01 α difference is immaterial to the figure.

## 5. Placement + scope

- Figs A, B + Paragraphs 1, 2 → **raw substrate Results subsection**
  (`sssec:methods_compare_rawfc` companion).
- Fig C + Paragraph 3 → **LRG cophenetic Results subsection**.
- **No verdict change.** Per-band trace verdicts stand as locked
  (`locked/VERDICT_LEDGER.md`); this is figure + narrative for an
  already-established control battery (drift = supplementary control;
  matched-strength = GATE; see `locked/CONTROLS.md`).
- Distinct from `methods/methods_displacement_taxonomy_2026-06-03.md` (that is
  a per-pair *character* decomposition of ρ_split^coph; this is the
  raw-vs-LRG control-separation figure set).
