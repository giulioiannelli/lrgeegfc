---
name: talk-slide-09-frequency-bands
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 10
status: draft
updated: 2026-07-11
canva: page 10 · ~35%
---

# Slide 9 — The temporal scale: frequency bands

1. TITLE
The temporal scale — frequency bands

2. MAIN CONCEPT
- Along the temporal axis of the multiscale brain we split the signal into frequency bands — different rhythms index different temporal scales of dynamics.
- Our six bands (the actual edges we use, from config — NOT textbook): δ 0.53–4 · θ 4–8 · α 8–13 · β 13–30 · γ_low 30–80 · γ_high 80–300 Hz.
- γ_high runs to 300 Hz — intracranial-specific: scalp EEG can't resolve high-frequency activity that far, but our clean depth contacts can (ties back to slide 3's resolution pro). Flag it as OUR definition, not the textbook ~60–150 Hz. (δ's 0.53 Hz low-cut is the acquisition band-pass edge, not a round 0.5.)
- Band-specificity is a built-in control: a genuine oscillatory effect is band-selective, whereas a broadband change — overall arousal, SNR drift, a non-oscillatory artefact — hits every band alike. (SEPARATE from the matched-strength null of slide 14, which handles the per-band magnitude confound; don't conflate the two controls.)
- We compute a separate FC per band — each band is one choice of b in slide 9's band-average A_ij = ⟨|ImCoh|⟩_{f∈b}.

3. ON-SLIDE TEXT
Bullets:
- rhythms = temporal scales (slow → fast)
- band-selective → a specific rhythm; broadband → a global shift
- one FC per band (6 bands)

Band table (colour slow=red → fast=blue per the band palette; γ_high = our intracranial definition):
δ 0.53–4 Hz · θ 4–8 · α 8–13 · β 13–30 · γ_low 30–80 · γ_high 80–300

Compile-ready LaTeX (band table, if you want it typeset):
\begin{array}{ll}\delta & 0.53\text{–}4\ \mathrm{Hz}\\ \theta & 4\text{–}8\\ \alpha & 8\text{–}13\\ \beta & 13\text{–}30\\ \gamma_{\mathrm{l}} & 30\text{–}80\\ \gamma_{\mathrm{h}} & 80\text{–}300\end{array}

4. SPEECH
The temporal axis of the multiscale brain is its rhythms, so we split the signal into six frequency bands — from delta below four hertz, up through beta, and gamma all the way to three hundred hertz. That top end is unusual: scalp EEG can't see high-frequency activity that far, but our clean intracranial contacts can, so we keep it. Splitting into bands isn't bookkeeping — it's a control. A genuine oscillatory effect should be selective for a band; a broadband change — overall arousal, a drift in signal quality — hits every band the same. So band-specificity, when we see it later, is itself a piece of evidence that the effect is real. We build one functional network per band.

5. FIGURES
- Per-band FC matrices, one patient — data/outputs/figures/fc_templates/row_per_band/Pat_05_rest_pre_imcoh_abs_chnames_log_shared.pdf

6. REFERENCES  (all verified M 07-11)
- LEAD — rhythms across scales (the temporal-axis anchor): Buzsáki, G. & Draguhn, A. (2004), "Neuronal oscillations in cortical networks", Science 304(5679), 1926–1929. DOI 10.1126/science.1099745. ✅ (Rhythms span ~5 orders of magnitude in frequency — the literal "temporal scales" claim.)
- Band-specific large-scale interactions: Siegel, M., Donner, T. H., Engel, A. K. (2012), "Spectral fingerprints of large-scale neuronal interactions", Nature Reviews Neuroscience 13(2), 121–134. DOI 10.1038/nrn3137. ✅
- Optional intracranial-specific support (our modality): Kucyi, A. et al. (2018), "Intracranial electrophysiology reveals reproducible intrinsic functional connectivity within human brain networks", J. Neurosci. 38(17), 4230–4242. ✅ title/locators verified (jneurosci.org/content/38/17/4230). DOI not asserted — pin from that page if you keep this ref.
- ✂ DROPPED (M 07-11): Betzel et al. (2019), Nat. Biomed. Eng. 3(11), 902–916 = "Structural, geometric and genetic factors predict interregional brain connectivity … electrocorticography" — tangential (predicts connectivity from structure/geometry/genetics, not "bands = temporal scale").

7. CANVA STATUS
Page 10 · ~45%. Content (M 07-11): added the REAL band table (δ 0.53–4 … γ_high 80–300 Hz, actual config edges) + compile-ready LaTeX array; refined the "control" claim (band-selective vs broadband; kept separate from the matched-strength null). Figure ready (per-band FC row, Pat_05 — F 07-11). Missing on the deck: paste band table + bullets, place the figure, presenter notes. Refs → Lane R (verify Betzel/Kucyi relevance; pin Buzsáki 2004 + Siegel 2012).
