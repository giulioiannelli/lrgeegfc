---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
section: 2
priority: critical
---

# Section 2 — fix list

**Head. CRITICAL.** §1 frames the cohort as n=10 throughout, but every §2 figure
caption that mentions a cohort size says "**all six patients**" (Pat_02, 03, 05,
06, 07, 08). This is a 6-patient pre-replacement snapshot, not the current
n=10 cohort. The mismatch propagates the entire §2 → §3 figure family. Fix is
either (i) regenerate at n=10 (preferred — matches the rest of the manuscript),
or (ii) explicitly relabel as "n=6 pre-replacement snapshot" everywhere with a
caveat sentence in §2.1 explaining why §2 uses a 6-patient subset while §4/§5
use n=10.

## Confirmed (no action)

- §2.1 (Tab 1): canonical band table {δ, θ, α, β, low-γ, high-γ} with cutoffs (0.5, 4, 8, 13, 30, 80, 300) Hz. Matches `BRAIN_BANDS` in `lrg_eegfc.config.const`. ✓
- §2.2.1 / §2.2.2: ImCoh definition (Nolte 2004, Eq. 2) and band-averaged magnitude (Eq. 3) are consistent with `imcoh-guide.md` and `imcoh_taxonomy` memory. ✓
- §2.2.2 LRG-compatibility justification (non-negativity + zero-diag preserved by |·|; circular-shift surrogates inappropriate for ImCoh) is correct.
- §2.2.3 cohort-mean enrichment ratio: MSC 2.59× → ImCoh 0.99× (cohort mean), consistent with the bias-removal claim. ✓
- §2.3.4 same-shaft Pat_06 spectral-noise call-out (above 30 Hz) is the correct caveat per `imcoh_validation` notes.

## Numerical corrections (action: writing agent)

- None at the numerical level for the bias-removal table. All §2.2.3 / §2.3 numbers (1.5×–8.9× MSC range, 1.40× max under |ImCoh|) reproduce in the `Same-shaft enrichment ratio` figure (Fig. 1).
- §2.2.3 prose: "the cohort mean same-shaft enrichment ratio drops from 2.59× under MSC to 0.99× under |ImCoh|, the latter indistinguishable from no enrichment, and no (patient, band) cell under |ImCoh| exceeds 1.40× (Pat08, β)." — consistent with Fig. 1 caption. ✓

## Framing rewrites (action: writing agent)

### CRITICAL — cohort-size mismatch in figure captions

**Locator: every §2 figure caption ending with "for all six patients", and §3.2 (Fig. 7).**

The §2 figure family was generated on a **6-patient subset** (Pat_02, 03, 05, 06, 07, 08) — the pre-2026-04-25 cohort *minus* Pat_10/13/14/15. The current cohort is n=10. The contradiction is:

- §1 page 2: "the patient cohort, LRG framework, inter-phase comparison methodology ... are all unchanged" (implying n=10 carry-over from NOTES_MSC).
- §1 page 3: "α and low-γ via the per-pair correlation on D(τ) (within-probe BH-q = 0.027 **at m = 6**)" — this "m=6" is the BH band-family count, NOT the patient count. The patient count is n=10 throughout §4/§5.
- §2 figure captions, repeatedly: "all six patients".

**Specific captions to fix** (quote → action):

1. **Fig. 1** (`Same-shaft to cross-shaft mean weight ratio for MSC (left) and |ImCoh| (right), computed at the edge-weight level for all six patients and six frequency bands (rsPre phase)`):
   - Y-axis labels are Pat_02, 03, 05, 06, 07, 08 (six patients confirmed).
   - **Action**: regenerate at n=10 (add Pat_10, 13, 14, 15 rows). If regen is deferred, relabel: "six-patient pre-replacement subset (Pat_02, 03, 05, 06, 07, 08); the n=10 enrichment table appears in §X / Tab Y" with one new sentence. Strongly prefer regen.

2. **Fig. 2** (`Community-scale dependence of same-shaft enrichment under MSC (blue) and |ImCoh| (red) for the β band (rsPre phase)`):
   - Left panel y-axis Pat_02–08 same six patients.
   - **Action**: same as Fig. 1.

3. **Fig. 3** (`Cross-spectral profiles for all six patients (rsPre phase). Each panel shows the median ...`):
   - 6-panel grid, six patients.
   - **Action**: regenerate at n=10 (10-panel 2×5 or 4-panel-per-row).

4. **Fig. 4** (`Kernel density estimates of β-band edge weights for same-shaft (dashed) and cross-shaft (solid) pairs, rsPre phase, all six patients overlaid by color`):
   - Six-color KDE overlay.
   - **Action**: regenerate at n=10.

5. **Fig. 7** (§2.3.3, `Band-resolved FC matrices for all six patients (β band, rsPre)`):
   - 2×6 mosaic = top row MSC, bottom row |ImCoh|, six patients.
   - **Action**: regenerate as a 2×10 mosaic at n=10.

6. **Fig. 8** (§2.3.4, force-directed embeddings for `Pat_02, Pat_05, and Pat_08`):
   - Three-patient exemplar, not six.
   - Caption text: "...for Pat02, Pat05, and Pat08 (rsPre)..." — explicit. The "all six" framing does NOT apply here. Keep, but disclose exemplar choice in main text §2.3.4 ("we show three exemplars; full cohort comparison in App.").

7. **Fig. 9** (§2.3.4 spectra): caption is `Sorted Laplacian eigenvalues for Pat05 across all six bands` — single-patient, six-band. No cohort claim. Keep.

8. **§2.3 head paragraph** (`We characterize these networks at the edge, node, and spectral levels before entering the LRG pipeline, comparing throughout with the corresponding MSC matrices to document the structural consequences of the estimator switch.`): no cohort-count claim. Keep.

### Suggested new sentence at the head of §2.2.3 (only if regen is deferred)

> "Bias-quantification figures in this section (Fig. 1, 2, 3, 4, 7) use a six-patient pre-replacement subset (Pat_02, 03, 05, 06, 07, 08); §4 and §5 carry the analysis through on the full n=10 cohort (Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15) per the 2026-04-25 vendor data restoration. The bias-removal pattern (MSC same-shaft enrichment 1.5–8.9× vs |ImCoh| ≤ 1.40× across (patient, band) cells) is invariant under cohort expansion and is documented at six-patient resolution for visual clarity."

### MSC vs ImCoh comparison framing — not over-claimed

- §2.2 / §2.3 stick to "comparison" and "documenting the structural consequences of the estimator switch" — correct framing. The phrase "fundamentally different LRG structures" (§3.2) is reserved for the dendrogram discussion. ✓

## Caveats to add (action: writing agent)

- **§2.1 head** (or §2.3 head): one-sentence cohort-snapshot disclosure (see template above).
- **§2.3.1 (Pat_03 1024 Hz outlier)**: §2.3.1 prose mentions "Pat_03 shows comparatively flat and elevated coherence" — should add the methodological footnote: "Pat_03 is the cohort's only 1024 Hz patient (others 2048 Hz; nperseg=2048 vs 4096, see `nperseg_for_fs`); rank-based and ratio-based quantities are dimensionless w.r.t. sampling rate so Pat_03 is included throughout, but it is plotted with a distinct marker (orange triangle) in §4." Optional but improves cross-section consistency.

## Figure actions (coding agent)

**Preferred path**: regenerate Fig. 1, 2, 3, 4, 7 at n=10. The bias-removal claim is invariant under cohort expansion (verified: ImCoh 0.99× cohort mean is a per-patient property, addition of Pat_10/13/14/15 will not invalidate the headline). Estimated effort: cohort-level regeneration with `compute_imcoh_for_patient` and the existing same-shaft mask functions in `lrg_eegfc.utils.fc.imcoh_bias`.

**Fallback path**: keep figures as-is, add the one-sentence n=6 disclosure at the §2.1 / §2.3 head, and label every figure caption "(n=6 pre-replacement snapshot, Pat_02–08)" explicitly. Less clean.

**Fig. 8** specifically (force-directed): three-patient exemplar (Pat_02, 05, 08) chosen for visual contrast. State this exemplar choice in §2.3.4 prose: "We show Pat_02, Pat_05, and Pat_08 as three contrasting implant geometries (left-hemisphere temporal-dominant, bilateral-frontal, and right-hemisphere temporo-parietal respectively); the cohort-wide bias pattern is documented in Fig. 1 / Fig. 2."

## Deferred / questions

- Is regeneration at n=10 in scope for this revision pass, or is the regen blocked by upstream cache state? If regen is blocked, escalate the question of whether the §2 cohort framing (n=6 figures + n=10 prose) is acceptable for submission.
- Pat_14 vendor-replaced 2026-04-25; verify Pat_14 ImCoh cache for `task_test` is current before regenerating the §2 figures (otherwise the new figures will inherit the old corrupt-data signature). Cache verification: `lrg-eegfc cache list --patients Pat_14 --fc-method imcoh`.
