---
name: talk-slide-08-imaginary-coherence
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 8
status: draft
updated: 2026-07-16
canva: page 9 · ~30%
---

# Slide 8 — Imaginary Coherence for FC estimation

1. TITLE
Imaginary Coherence for FC estimation

2. MAIN CONCEPT
- Coherence measures how consistently two signals share a phase relationship at a frequency. Build it from the cross-spectrum: Welch-estimate the complex cross-spectral density S_ij(f) = ⟨X_i(f) X_j*(f)⟩ (average of windowed-FFT cross-products over Hann segments); its diagonal S_ii(f) is the power spectrum.
- Coherency C_ij(f) = S_ij / √(S_ii S_jj) is the normalized cross-spectrum (|C| ≤ 1); |C_ij|² is the usual (magnitude-squared) coherence.
- The trap: volume conduction — and our common reference — make one generator appear in many contacts at the same instant (zero phase lag). Zero-phase mixing is a *purely real* contribution to S_ij, so ordinary coherence counts the leakage as connectivity.
- Nolte's fix (2004): keep the IMAGINARY part. ImCoh_ij(f) = Im C_ij(f) ∈ [−1, 1]. Zero-phase leakage has no imaginary part → it cancels exactly; what survives is genuine, time-lagged coupling. (Signed ImCoh is skew-symmetric; its sign encodes lead/lag.)
- Bastos & Schoffelen (2016) confirmation: their tutorial review shows side-by-side that ordinary coherence lights up with field-spread artefacts while imaginary coherence stays clean (the right-hand figure) — the independent confirmation that ImCoh removes the volume-conduction bias. In our sEEG this bites hardest for contacts on the *same electrode shaft*, where leakage is worst and a raw measure like MSC is fooled; ImCoh stays robust there.
- Our band FC: magnitude per frequency bin, then average across the band — W_ij = ⟨|ImCoh_ij(f)|⟩_{f∈b} (Ewald 2012 / Bastos-Schoffelen 2016). Non-negative, as the diffusion Laplacian (slide 10) needs; keeps strength, drops lead/lag direction.

3. ON-SLIDE TEXT
(matches the rendered deck — 3 bullets, each carrying ONE formula, in this order; the volume-conduction beat is now carried by the right-hand figure, NOT a bullet. ⟨·⟩ = average over Welch Hann segments; S_ii is the i=j case = the power spectrum.)

- coherence = phase consistency at a frequency
    C_{ij}(f)=\dfrac{S_{ij}(f)}{\sqrt{S_{ii}(f)\,S_{jj}(f)}}
- built from the cross-spectrum (Welch)
    S_{ij}(f)=\big\langle X_i(f)\,X_j^{*}(f)\big\rangle
- our FC = band-averaged magnitude of ImCoh
    \mathrm{ImCoh}_{ij}(f)=\operatorname{Im}\,C_{ij}(f)\in[-1,1]
    \;\longrightarrow\; W_{ij}=\big\langle\,\big|\mathrm{ImCoh}_{ij}(f)\big|\,\big\rangle_{f\in b}

Figure caption (right, on-slide): same-shaft volume conduction effect robust (e.g. vs MSC)

4. SPEECH
How do we measure connectivity between two contacts? We start from coherence — at a given frequency, do two signals hold a consistent phase relationship? We build it from the cross-spectrum, the Welch method, and normalize.

The trap, especially in intracranial recordings: volume conduction. One source bleeds into many contacts at the same instant — zero lag — and zero-lag mixing is purely *real*. Ordinary coherence counts that leakage as a connection that isn't there.

Nolte's fix, in 2004: keep only the *imaginary* part. Zero-phase leakage has no imaginary component, so it cancels — what's left is genuine, time-lagged coupling. And Bastos and Schoffelen confirm it directly [gesture]: ordinary coherence lights up with field-spread artefacts, imaginary coherence stays clean — robust even for contacts on the same shaft, where leakage is worst.

So our connectivity, in one line: per band, take the magnitude of imaginary coherence and average across the band. That non-negative number, W-i-j, is the edge weight — it keeps the strength, drops the direction, and being non-negative is exactly what the diffusion step needs.

Careful (do NOT say / keep honest):
- The three beats map to the slide: (1) what coherence / ImCoh IS, in plain terms; (2) Nolte 2004 *introduced* the fix — keep the imaginary part to null zero-phase leakage; (3) Bastos & Schoffelen 2016 is the *confirmation* that it removes volume-conduction bias (the side-by-side figure). Credit Nolte for the measure, Bastos for the robustness demonstration — don't swap them.
- ImCoh's price (say it if asked, don't hide it): rejecting zero-lag mixing also discards any *genuine* zero-lag coupling. That's the accepted trade for volume-conduction immunity — do NOT claim ImCoh captures all true coupling; it captures lag-bearing coupling and rejects instantaneous leakage.
- The "vs MSC" on the figure is a same-shaft leakage-robustness *motivation*, not a results claim — no AUC / no numbers here.
- W_ij is non-negative BY the magnitude — required for the diffusion Laplacian (slide 10). It keeps strength, drops lead/lag direction. Do not call it "directed."

5. FIGURES
- Volume-conduction / field-spread figure — ✅ EXTRACTED: `data/outputs/figures/talk/_external/bastos2016_fig6_field_spread_imcoh.png` = Bastos & Schoffelen 2016 **Fig. 6** "Effects of field spread on the estimation of connectivity" (coherence vs imag-coh, panels A–D). The "why ImCoh" motivator. 🟢 CC BY 4.0 → on-slide credit: *Bastos & Schoffelen, Front. Syst. Neurosci. 9:175 (2016), Fig. 6 — CC BY 4.0*.
- ImCoh spectral distribution, ONE patient — data/outputs/figures/section2/fig_C/fig_C2_spectral_distribution_imcoh_rsPre.pdf (imcoh only; drop the msc sibling). Lane F: confirm single-patient rendering, else produce one.
- (optional "real ImCoh" exemplar, ✅ EXTRACTED) `…/_external/nolte2004_fig12_imcoh_beta_headmap.png` = Nolte 2004 **Fig. 12**, β head-map (the clean pick); alt `…/_external/nolte2004_fig11_significant_imcoh.png` = Fig. 11 significance maps. 🟡 © Elsevier → talk-only, credit *Nolte et al., Clin. Neurophysiol. 115 (2004), Fig. 12/11*.
- The four formulas render as on-slide LaTeX (compiled in Canva), not Lane-F PDFs.
- ~~MSC vs ImCoh dendrograms (fig_H1)~~ — REMOVED from slide 8 (user 07-11: drop the MSC comparison here). Still lives on slide 11 (pipeline).
- ~~ImCoh band-average (imcoh_band_average.pdf)~~ — REMOVED earlier (wrong kind of figure).

6. REFERENCES
- Nolte, G., et al. (2004), "Identifying true brain interaction from EEG data using the imaginary part of coherency", Clinical Neurophysiology 115(10), 2292–2307. DOI 10.1016/j.clinph.2004.04.029. ✅ verified. — source of the ImCoh formula Im(S_ij)/√(S_ii·S_jj) + the volume-conduction-immunity argument. (The paper's finger-movement 20 Hz result IS the Fig. 12 β head-map extracted for this slide.)
- Ewald, A., Marzetti, L., Zappasodi, F., Meinecke, F. C., Nolte, G. (2012), "Estimating true brain connectivity from EEG/MEG data invariant to linear and static transformations in sensor space", NeuroImage 60(1), 476–488. DOI 10.1016/j.neuroimage.2011.11.084. ✅ verified. — the <|ImCoh|> magnitude convention we use.
- Bastos, A. M. & Schoffelen, J.-M. (2016), "A tutorial review of functional connectivity analysis methods and their interpretational pitfalls", Frontiers in Systems Neuroscience 9, 175. — source of the volume-conduction figure for this slide + review. [Lane R: fetch the figure + pin DOI.]

7. CANVA STATUS
Page 9 · ~30%. Have: title, Canva auto-filler bullets, the references. Missing: the ImCoh figures, compressed on-slide text, presenter notes.
