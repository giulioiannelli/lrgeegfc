---
name: talk-slide-08-imaginary-coherence
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 9
status: draft
updated: 2026-07-11
canva: page 9 · ~30%
---

# Slide 8 — Imaginary coherence for FC

1. TITLE
Imaginary coherence for FC

2. MAIN CONCEPT
- Coherence measures how consistently two signals share a phase relationship at a frequency. Build it from the cross-spectrum: Welch-estimate the complex cross-spectral density S_ij(f) = ⟨X_i(f) X_j*(f)⟩ (average of windowed-FFT cross-products over Hann segments); its diagonal S_ii(f) is the power spectrum.
- Coherency C_ij(f) = S_ij / √(S_ii S_jj) is the normalized cross-spectrum (|C| ≤ 1); |C_ij|² is the usual (magnitude-squared) coherence.
- The trap: volume conduction — and our common reference — make one generator appear in many contacts at the same instant (zero phase lag). Zero-phase mixing is a *purely real* contribution to S_ij, so ordinary coherence counts the leakage as connectivity.
- Nolte's fix: keep the IMAGINARY part. ImCoh_ij(f) = Im(S_ij)/√(S_ii S_jj) ∈ [−1, 1]. Zero-phase leakage has no imaginary part → it cancels exactly; what survives is genuine, time-lagged coupling. (Signed ImCoh is skew-symmetric; its sign encodes lead/lag.)
- Our band FC: magnitude per frequency bin, then average across the band — A_ij = ⟨|ImCoh_ij(f)|⟩_{f∈b} (Ewald 2012 / Bastos-Schoffelen 2016). Non-negative, as the diffusion Laplacian (slide 7) needs; keeps strength, drops lead/lag direction.

3. ON-SLIDE TEXT
Bullet list (on-slide):
- coherence = phase consistency at a frequency
- built from the cross-spectrum S_ij(f) (Welch)
- volume conduction = zero-phase → purely real → fake coupling
- keep the imaginary part → only genuine, lagged coupling
- our FC = band-averaged magnitude ⟨|ImCoh|⟩

Formulas — compile in Canva's LaTeX tool (each line introduces ONE new quantity; ⟨·⟩ = average over Welch Hann segments; S_ii is the i=j case = power spectrum):

(1) cross-spectrum
S_{ij}(f)=\big\langle X_i(f)\,X_j^{*}(f)\big\rangle

(2) coherency
C_{ij}(f)=\dfrac{S_{ij}(f)}{\sqrt{S_{ii}(f)\,S_{jj}(f)}}

(3) imaginary coherence (Nolte 2004)
\mathrm{ImCoh}_{ij}(f)=\operatorname{Im}\,C_{ij}(f)\in[-1,1]

(4) our band FC
A_{ij}=\dfrac{1}{|b|}\sum_{f\in b}\big|\mathrm{ImCoh}_{ij}(f)\big|

4. SPEECH
Now the connectivity itself. We start from coherence — how consistently two signals hold the same phase relationship at a given frequency. To build it we take the cross-spectrum: cut each channel into overlapping windowed segments, Fourier-transform them, and average the cross-products — that's the cross-spectral density, S-i-j; its diagonal is just the power spectrum. Normalize it and you get coherency; its squared magnitude is the coherence people usually quote. But here's the trap: volume conduction — and our common reference — make one generator show up in many contacts at the very same instant, at zero phase lag. Zero-phase mixing is purely real, and ordinary coherence swallows it as connectivity. Nolte's trick is to keep only the imaginary part: zero-phase leakage has no imaginary component, so it cancels exactly. What's left — the imaginary coherence — is genuine, time-lagged coupling. For each band we take its magnitude and average across frequencies, and that non-negative number is our connectivity strength — non-negativity being exactly what the diffusion Laplacian needs.

5. FIGURES
- Volume-conduction / field-spread figure — ✅ EXTRACTED: `data/outputs/figures/talk/_external/bastos2016_fig6_field_spread_imcoh.png` = Bastos & Schoffelen 2016 **Fig. 6** "Effects of field spread on the estimation of connectivity" (coherence vs imag-coh, panels A–D). The "why ImCoh" motivator. 🟢 CC BY 4.0 → on-slide credit: *Bastos & Schoffelen, Front. Syst. Neurosci. 9:175 (2016), Fig. 6 — CC BY 4.0*.
- ImCoh spectral distribution, ONE patient — data/outputs/figures/section2/fig_C/fig_C2_spectral_distribution_imcoh_rsPre.pdf (imcoh only; drop the msc sibling). Lane F: confirm single-patient rendering, else produce one.
- (optional "real ImCoh" exemplar, ✅ EXTRACTED) `…/_external/nolte2004_fig12_imcoh_beta_headmap.png` = Nolte 2004 **Fig. 12**, β head-map (the clean pick); alt `…/_external/nolte2004_fig11_significant_imcoh.png` = Fig. 11 significance maps. 🟡 © Elsevier → talk-only, credit *Nolte et al., Clin. Neurophysiol. 115 (2004), Fig. 12/11*.
- The four formulas render as on-slide LaTeX (compiled in Canva), not Lane-F PDFs.
- ~~MSC vs ImCoh dendrograms (fig_H1)~~ — REMOVED from slide 9 (user 07-11: drop the MSC comparison here). Still lives on slide 11 (pipeline).
- ~~ImCoh band-average (imcoh_band_average.pdf)~~ — REMOVED earlier (wrong kind of figure).

6. REFERENCES
- Nolte, G., et al. (2004), "Identifying true brain interaction from EEG data using the imaginary part of coherency", Clinical Neurophysiology 115(10), 2292–2307. DOI 10.1016/j.clinph.2004.04.029. ✅ verified. — source of the ImCoh formula Im(S_ij)/√(S_ii·S_jj) + the volume-conduction-immunity argument. (The paper's finger-movement 20 Hz result IS the Fig. 12 β head-map extracted for this slide.)
- Ewald, A., Marzetti, L., Zappasodi, F., Meinecke, F. C., Nolte, G. (2012), "Estimating true brain connectivity from EEG/MEG data invariant to linear and static transformations in sensor space", NeuroImage 60(1), 476–488. DOI 10.1016/j.neuroimage.2011.11.084. ✅ verified. — the <|ImCoh|> magnitude convention we use.
- Bastos, A. M. & Schoffelen, J.-M. (2016), "A tutorial review of functional connectivity analysis methods and their interpretational pitfalls", Frontiers in Systems Neuroscience 9, 175. — source of the volume-conduction figure for this slide + review. [Lane R: fetch the figure + pin DOI.]

7. CANVA STATUS
Page 9 · ~30%. Have: title, Canva auto-filler bullets, the references. Missing: the ImCoh figures, compressed on-slide text, presenter notes.
