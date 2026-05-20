---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
section: 5
---

# Section 5 — fix list

**Head.** §5 is the load-bearing analytical section and its numbers reproduce
across all six verifier reports. The ψ-irrelevance, KC β/low-γ, CTM trace bands,
VI(k) / Grassmann shapes, anatomy γ_l fusiform Bonferroni cell, and §5.6 cohort
totals all check out. Five issues need writing-level action: (i) §5.1 lacks an
independent matched-strength surrogate confirming Ψ-failure is generic to dense
weighted graphs; (ii) §5.4 has a Grassmann-formula inconsistency (two equations
that differ by √2) — only `d_G^2 = sum sin²θ_i` matches the code; (iii) §5.5 should
disclose what fraction of Hippocampus / fusiform / superior-temporal contact pools
are epileptic contacts; (iv) §5.6 classes are NOT mutually exclusive at fixed
scale — same leaf can be in a trace L1 AND an anchor L2 because the matching
pipeline allows independent acceptance; (v) τ is fixed at 1/λ_max — no multiscale
class assignment, owed in Outlook.

## §5.1 — Ψ-irrelevance + ρ̂(τ)→D(τ) reframing

### Confirmed (no action)

- **89/180 cells (49%) strict argmax-at-extremum** under the top-down convention. 21 top-cut, 68 bottom-cut. Bottom-pinning dominates 3:1 (`_psi_argmax_pinning.md`). ✓
- **94% within 5 indices of a boundary**, only **2/180 (1%) at d > 10** (cohort-level deep-interior outliers). ✓
- Per-band gradient: δ 27% → high-γ 70% strict-pinned, with β = 43% / α = 50% / θ = 63% / low-γ = 43%. ✓
- The "merge-height profile decays linearly rather than log-linearly" mechanism (Fig. 21 prose) matches the dense-weight-heterogeneous outlier-case framing of `lrg_outlier_case_fully_connected.md`.
- τ' = 1/λ_max as the finest-resolved scale is correct; the prose's reading "coarser scales τ ≫ τ' wash out the heterogeneity that carries the trace" is consistent with `lrg_tau_choice.md`. ✓

### Numerical corrections (action: writing agent)

- §5.1 cites "**89 sitting within ten merges of an extremum**" — verify against `_psi_argmax_pinning.md` Tab. The 89 figure refers to **strict d=0** count (argmax exactly at extremum), not "within 10". Within-10 counts are 178/180. Recommend rephrasing: "of the 91 strictly-interior cells, only 2 have argmax more than ten merges from the nearer boundary, with 89 sitting within ten merges of an extremum." — current phrasing is correct as written but easy to misread (89 vs 91 vs 178). Tighten by re-indexing to "the 91 strictly-interior cells of which 89 sit within ten of an extremum and only 2 reach d > 12".
- §5.1 last paragraph: "169/180 (94%) within d ≤ 5; only 2 cells (1%) beyond d = 12" — matches `_psi_argmax_pinning.md`. ✓

### Framing rewrites (action: writing agent)

- **§5.1 prose currently reasons from one example cell (Pat_02 β rsPre, Fig. 20) to the 180-cell cohort statistic (Fig. 21)**. Keep — this is the right rhetorical shape. But a one-sentence forward to a **matched-strength surrogate** would close the "is this generic to dense weighted graphs or specific to |ImCoh|?" question:

  > "An independent matched-strength surrogate (preserving each node's strength sequence under permuted edge identities, see Maslov-Sneppen 2002) reproduces the same Ψ-pinning pattern on randomized |ImCoh| matrices, confirming that the failure mode is generic to dense weight-heterogeneous graphs rather than specific to the |ImCoh| cohort substrate."

  If this surrogate has not been run, flag as **deferred**. The §5.1 reading does not depend on it (the empirical 89/180 cohort statistic and the merge-height-linearity mechanism are sufficient), but it would harden the §5.1 → §5.2 transition.

- **§5.1 closing — `ρ̂(τ)→D(τ)` reframing language**: "the natural LRG object is therefore not the single Ψ-selected partition of ρ̂(τ) but the multiscale structure encoded by ρ̂(τ) and its derived dendrogram — read as a continuous hierarchy rather than collapsed onto one privileged cut." This is the correct reframing per `feedback_no_l2_entropy_curve` (no S(τ) / C(τ) curve as a trace metric — the L2 entropy/susceptibility rung is permanently retired) and `feedback_lrg_step_is_confirmation_not_resolution` (LRG step is structural enrichment, not pure confirmation). ✓

### Caveats to add (action: writing agent)

- **§5.1 head**: insert a one-line acknowledgment that the failure mode is hypothesized-generic but not independently verified by surrogate: "The mechanism (linearly-decaying merge sequence on dense weight-heterogeneous graphs amplifying the smallest-magnitude tail) is generic in principle; an independent matched-strength surrogate confirming this on randomized |ImCoh| matrices is deferred to §5.6 / future work." ← only if the surrogate has not been run; otherwise replace with the suggested forward sentence above.

### Figure actions (coding agent)

- **Fig. 20** (Ψ profile on Pat_02 β rsPre): keep. This is the canonical exemplar.
- **Fig. 21** (cohort distribution of Ψ argmax position folded to d): caption says "lower bound at d ≥ 12 with 2 cells (1%)" — verify caption text matches `_psi_argmax_pinning.md` Tab. ✓

## §5.2 — KC tree distance

### Confirmed (no action)

- **β λ=0**: 7/10 imprint, p = 0.0186 (W=7), within-probe BH q = 0.006. λ=0 mean |T| = 1.02 over positives vs 5.61 over negatives (5.5× louder). Pat_03 anchors at |T|=19.73 but signal is not Pat_03-only (Pat_06 5.09, Pat_08 5.07, Pat_13 3.16, Pat_07 2.43, Pat_10 2.45 all aligned). Per `_kc_beta_pvalue_check.md`. ✓
- **β λ=1**: 8/10 imprint, p = 0.0420 (W=10), within-probe BH q = 0.006. ✓
- **γ_l λ=1**: 8/10 imprint, p = 0.0420, within-probe BH q = 0.146; carried as a "directional companion" because it does NOT pass within-probe q ≤ 0.05. ✓
- **Cross-probe**:
  - β λ=0: 7/10, p = 0.0244, within-probe BH q = 0.146 (cross-probe table).
  - β λ=1: 8/10, p = 0.0420, q_BH 0.168.
  - low-γ λ=1: 9/10, p = 0.0244, q_BH 0.146.
  - α λ=0/λ=1: 6/10 / 6/10, p > 0.46 — null control intact.
  - Per `_kc_crossprobe_check.md`. ✓
- **Pat_03 dropout** at β: λ=0 and λ=1 both preserve directional sign at p = 0.002 on n=9 paired Wilcoxon vs the within-baseline null. ✓
- **Joint-BH at m=12** (KC family): both β cells survive at q = 0.006 (within-probe correction). ✓
- **Per-leaf decomposition (Fig. 23)**: per-patient counts |P| = 37 / 31 / 38 (Pat_07 / Pat_05 / Pat_10) at fixed-predicate β; |P| = 1 at high-γ Pat_14 (sanity-check trivially). f* values 0.32 / 0.26 / 0.34 / 0.01 respectively. Matches `kc_trace_module_visualization.md` and 2026-05-07 audit_47.

### Numerical corrections (action: writing agent)

- §5.2 prose says "Both β cells survive at p_BH ≤ 0.006" — verify the q-value is 0.006 (not 0.005 or 0.0058). Per `_kc_beta_pvalue_check.md` and `_kc_crossprobe_check.md`, q_BH = 0.006 at within-probe family-of-12. ✓
- §5.2 prose claims **"all ten patients fall below their own null"** for the within-baseline-null β paired Wilcoxon. The §5.2 verifier `_kc_crossprobe_check.md` does not explicitly tabulate this, but the resulting paired Wilcoxon p = 0.001 / 0.001 at λ=0/λ=1 implies sign-uniform alignment of the 10 patients with their respective nulls — verify in `data/audit/section5_v2_kc_controls/per_patient_table.csv`. If 9/10 patients fall below their own null but Pat_X does not, this should be a one-line caveat.
- §5.2 cross-probe restriction: "**at λ = 0, the same-probe pairs contribute a coherent cohort-mean imprint of +0.205 per patient** (cohort standard deviation 0.195) in the trace direction" — matches `_kc_crossprobe_check.md` β λ=0 mean(xpr−full) = +0.205, std = 0.195. ✓
- §5.2 low-γ: "**under the cross-probe restriction the directional count moves to 9/10 at p = 0.024, but the shift is driven by a single patient (Pat_07) whose T_KC sign-flips by less than 0.2 at the noise floor of this cell**" — exactly matches `_kc_crossprobe_check.md` (Pat_07: +0.087 → −0.084, both within noise of zero). ✓ This caveat is correctly stated in the §5.2 prose.

### Framing rewrites (action: writing agent)

- The §5.2 framing of low-γ as "trace at the dendrogram level is therefore a shift in the ultrametric heights without a corresponding rearrangement of which leaf-pairs share their most-recent common ancestor" is the right reading per `_kc_beta_pvalue_check.md` (λ=1 captures heights, λ=0 captures topology, low-γ registers only at λ=1). ✓
- The "β: only band whose trace registers in both tree topology AND tree heights under the within-baseline-null control with all ten patients falling below the null on both axes" is the load-bearing β-trace claim. Verified at q = 0.006 within-probe BH(m=12), per `result_2_lrg_beta_trace.md`. ✓
- The figure-23 reading ("Pat_07/Pat_05/Pat_10 reproduce the same β topology trace signature on different cortical territories, confirming the β topology trace is reproducible across patients and not a single-patient artifact") is well-supported by the |P| = 37 / 31 / 38 counts and the 2026-05-07 audit_47 visualization. ✓

### Caveats to add (action: writing agent)

- **§5.2 explicit framing of "louder patients with quieter dissenters"** at λ=0: per `_kc_beta_pvalue_check.md`, the smaller-count-at-λ=0-but-smaller-p-than-λ=1 paradox resolves into a magnitude effect (negative-direction patients are 5.5× louder than positive-direction dissenters at λ=0; only 1.6× louder at λ=1). The §5.2 prose ("the cleaner-by-magnitude λ=0 cell carries cohort mean |T_KC| of 5.61 over its 7 imprint-direction patients against 1.02 over the 3 dissenters, a 5.5× ratio; at λ=1 the corresponding ratio is only 1.6× (2.87 versus 1.83)") matches the verifier exactly. ✓ Keep as written.

### Figure actions (coding agent)

- **Fig. 22** (KC λ=0/λ=1 cohort-median imprint direction, per band): both panels match the verifier numbers. ✓
- **Fig. 23** (per-leaf decomposition for Pat_07/Pat_05/Pat_10/Pat_14): linear-y, per-row labels (focal i*, f*, |P|), red bundles up to MRCA depth. Matches `kc_trace_module_visualization.md`. ✓ Keep.

## §5.3 — Controlled trace via per-pair correlation on D(τ) (CTM)

### Confirmed (no action)

- **α**: ρ̃_split = +0.115, p = 0.007, q_BH(m=6) = 0.027, **8/10 positive AND 8/10 above drift** (intersection 8/10).
- **β**: ρ̃_split = +0.222, p = 0.014, q_BH(m=6) = 0.027, **8/10 positive ≈ 8/10 above drift** (intersection 7/10 — one-patient swap inside the same headline count).
- **low-γ**: ρ̃_split = +0.140, p = 0.010, q_BH(m=6) = 0.027, **7/10 positive < 9/10 above drift** (Pat_07 and Pat_13 sit slightly negative on ρ_split but exceed their own drift baselines; the controlled criterion is the load-bearing one).
- δ: 6/10 / 7/10, ρ̃_split = +0.031, p = 0.246. Borderline.
- θ: 3/10 / 6/10, ρ̃_split = −0.049, p = 0.278. Drift-only.
- high-γ: 4/10 / 5/10, ρ̃_split = −0.014, p = 0.423. Borderline.
- Per `_ctm_lowgamma_counts.md`. ✓
- **Cross-probe ratios**: α 1.9×, β 1.25×, low-γ 1.21× (same-probe / cross-probe ρ_split cohort-median ratios). Per §5.3 prose. ✓
- **Pat_03 dropout** (n=9): α p=0.014, β p=0.027, low-γ p=0.014. α and low-γ remain BH-surviving; β remains directional and uncorrected-significant. ✓
- **Pro-cohort restriction (n=8, drop Pat_07 + Pat_15)**: p = 0.027 at α, 0.039 at β, 0.027 at low-γ — sensitivity check, all preserve signal. ✓
- **Joint BH at m=48**: smallest q = 0.164 at CTM α, CTM low-γ, CTM β, Grassmann low-γ k=13 (per `_section5_joint_bh.md` family B). 0 cells survive q ≤ 0.05. The cited m=48 family in §5.3 closing matches the verifier's family B (KC ×2 + Grassmann ×1 + d_rank ×3 + CTM ×1 + VI(k) ×1 = 8 per band × 6 bands).

### Numerical corrections (action: writing agent)

- **§5.3 prose at low-γ**: the verifier `_ctm_lowgamma_counts.md` flags the apparent contradiction "7/10 ... 9/10" as **NOT a typo** — both numbers are real (7/10 positive on ρ_split + 9/10 positive on ρ_split > ρ_drift). Current §5.3 prose text "**nine of ten patients above their drift floor and seven of ten with ρ_split positive in absolute terms, with two patients (Pat_07 and Pat_13) slightly negative on ρ_split but more negative on ρ_drift and therefore passing the controlled criterion while failing the absolute one**" — reproduces both counts and the mechanism precisely. ✓ Keep.
- **§5.3 closing joint-BH count**: prose says "the forty-eight (probe, band) cells of the five LRG probe families assembled in this section". The §5.3 verifier `_section5_joint_bh.md` confirms the breakdown (KC at λ ∈ {0, 1} = 12 cells; CTM = 6 cells; matrix distances on D(τ) at d_S/d_P/d_F = 18 cells; VI(k) = 6 cells; Grassmann = 6 cells; total = 48). ✓ Smallest joint q = 0.164 = "no single cell survives at q ≤ 0.05; the smallest joint q is 0.164, reached by the three CTM trace bands."

  **Verify**: §5.3 prose says "**the smallest joint q ≤ 0.05; the smallest joint q is 0.164**". Per `_section5_joint_bh.md` family B, the four cells achieving q = 0.164 are CTM α, CTM low-γ, CTM β, AND Grassmann low-γ k=13. So either (i) §5.3 should mention all four cells reach the joint smallest q, or (ii) restrict the "three CTM trace bands" claim to CTM-only at q = 0.164. The current prose says "reached by the three CTM trace bands" — this is **incomplete** because Grassmann low-γ k=13 ties with them; correct to "reached by the three CTM trace bands and the Grassmann low-γ k=13 cell". Minor.

### Framing rewrites (action: writing agent)

- **§5.3 head**: clean, includes the controlled criterion definition (ρ_split > ρ_drift) and the BH(m=6) survival statement. Keep.
- **§5.3 cross-probe restriction**: "the cross-probe restriction discards real signal rather than removing bias" at α/β/low-γ. Correctly framed per `feedback_renormalization_style` and per the §5.3 verifier ratio readings. The β 1.25× and low-γ 1.21× ratios indicate that same-probe pairs at these bands carry near-baseline signal contributions, not a multiplicative bias as under MSC. ✓

### Caveats to add (action: writing agent)

- None at §5.3 level beyond what is already there. The within-probe BH(m=6) → joint BH(m=48) hierarchy is the right honesty hand-off.

### Figure actions (coding agent)

- **Fig. 24** (CTM cohort ρ_split per band): three boxplots per band (split / drift / xprobe), per-patient lines, p-value asterisks. Matches §5.3 prose. ✓
- **Fig. 25** (per-pair scatter at three exemplar (patient, band) cells: Pat_06 β + Pat_05 α + Pat_08 θ): same-probe vs cross-probe blue/grey, diagonal y=x reference. ✓ Per `_ctm_per_pair_scatter.md`. Caveat: the manuscript displays "horizontal/vertical banding" and attributes it to "pairs sharing a common contact, which carry correlated communication distances by construction on D(τ)". Verify this banding interpretation reproduces in the cached scatter PDF. ✓

## §5.4 — VI(k) and Grassmann

### Confirmed (no action)

- **VI(k) longest contiguous trace span**: β k=17–32 (16 cells, n_trace ≥ 7). One-cell dip at k=33; resumes at k=34–35. Peak at k=22 with **9/10**. Per `_section5_4_verify.md`. ✓
- **VI(k) low-γ**: contiguous span k=45–56 (12 cells), peak 8/10 at k=48. ✓
- **VI(k) α**: peak 9/10 at k=3, with non-contiguous extension at k=13–15, 18–21, 23–24, 28, 34–37. ✓
- **Pat_03 dropout** at α k=3, β k=22, low-γ k=48: all three preserve majority trace under Pat_03 drop. β k=22 promotes to 9/9 unanimous under Pat_15 drop. ✓
- **Mask sensitivity** (singleton fraction 0.4 / 0.5 / 0.6): β k=17–32 contiguous span unchanged; low-γ k=45–56 unchanged. ✓
- **Grassmann diagonal positive at i ≈ k**: β sustained from k ≈ 12 to k = 80 (81% of k ≥ 12 cells at n_trace ≥ 7); low-γ sustained k ≈ 10 to k = 80; α emerging at k ≈ 10 and sustained. ✓
- **Grassmann uncorrected p < 0.05 cells**: ~35 cells over the full {2..80} k grid; dominated by β at k ≥ 40 and low-γ at k = 10–28. **The strongest single cells are low-γ k=13 and k=14, both at p=0.014**. ✓ Per `_section5_4_verify.md`.

### Numerical corrections (action: writing agent)

- **§5.4 Grassmann formula inconsistency** (`_section5_4_verify.md` Section A). The §5.4 prose currently states two equations:
  - `d_G(U, U') = ||U U^T − U' U'^T||_F` (Frobenius identity, **option a**).
  - `d_G²(U, U') = sum sin²(θ_i)` (matching code, **option c**).

  These differ by a factor of √2 (`||P − P'||_F² = 2 · sum sin²θ`). The code computes `d_G = sqrt(sum sin² θ)` (option c).

  **Action**: replace the Frobenius identity by either of the equivalents:
  - `d_G(U, U') = sqrt(sum sin²θ_i)` (matches code), or
  - `d_G(U, U') = (1/√2) · ||U U^T − U' U'^T||_F` (equivalent to the code).

  As-written, the Frobenius and the sum-of-sin² lines are mutually inconsistent.

- **§5.4 β VI(k) span claim** (current prose: "**a sixteen-cell block at n_trace ≥ 7/10 between k=17 and k=32, with a one-cell dip at k=33 and resumption at k=34–35**"). Per `_section5_4_verify.md`, this is correct. ✓ Keep.

- **§5.4 low-γ VI(k) cells**: prose says "**a twelve-cell block at k=45–56 reaching n_trace = 8/10 at k=48, with continuation through k=77**" — matches `_section5_4_verify.md` (k=45–56 contiguous, with extensions at 58, 60–64, 67–69, 72–77). ✓

- **§5.4 Grassmann strongest cell**: the prose (in §5.4 closing on the chordal scalar) currently says **"the one-sided cohort Wilcoxon on the chordal scalar reaches p = 0.014 at k=13 (8/10) and k=14 (9/10)"**. ✓ Matches `_section5_4_verify.md`.

- **§5.4 ~35 (band, k) cells uncorrected p<0.05** statement reproduces. The key qualifier is that "these counts are not reported as controlled cohort claims (no split-baseline null is applied at this layer)". ✓ Honest framing.

### Framing rewrites (action: writing agent)

- **CRITICAL**: the manuscript caveat at the end of §5.4 — "**the one caveat that does not already follow from the descriptive role announced in the opening of this subsection is the meaning of k itself: cutting the dendrogram at the same k on networks of different size samples the tree at different relative depths, and the leading k Laplacian eigenmodes capture different fractions of the diffusive response on networks of different size and edge-weight density**" — should be **lifted into the head of §5.4**, not buried at the end. Per `feedback_renormalization_style`: VI(k) at fixed k across patients of different N (range 113–122 here) is apple-to-oranges; the caveat is the load-bearing methodological caveat for the whole section.

  **Suggested head insertion** (after the introductory two paragraphs on the two probes' definitions):

  > "Both VI(k) and Grassmann d_G(k) are evaluated at the **same nominal k across patients with different contact counts N ∈ [113, 122]**. Cutting the dendrogram at k = k₀ samples each patient's tree at a slightly different relative depth (k₀/N), and the leading k Laplacian eigenmodes capture different fractions of the diffusive response on networks of different size. The cohort-aggregate counts reported below are therefore patient-count claims at fixed k, not per-patient claims at the same partition resolution; this is the §5.4 layer-specific caveat that motivates the per-patient evaluation in §5.6."

- **§5.4 Grassmann uncorrected p<0.05 framing**: current prose ("**roughly 35 (band, k) cells reach uncorrected one-sided p<0.05 on the chordal scalar, dominated by β at k ≥ 40 and low-γ at k ≈ 10–28**") matches the verifier. The framing as a directional consistency statement ("not reported as controlled cohort claims") is honest. ✓

### Caveats to add (action: writing agent)

- See "lift the apple-to-oranges k caveat into the §5.4 head" above.

### Figure actions (coding agent)

- **Fig. 26** (cohort VI(k) heatmap + per-patient panels): top panel is the cohort n_trace(b, k)/10 heatmap; bottom is 10-patient panel grid with band on vertical axis. Matches `_section5_4_verify.md` reading. ✓
- **Fig. 27** (Grassmann mode-resolved heatmap per band): six-panel band layout, each with mode-index × cutoff-k heatmap of Δθ_i(k) and a strip below showing chordal-scalar n_trace/k at the dashed 7/10 cohort line. Light-shaded background highlights α/β/low-γ trace bands. ✓ Matches verifier.

## §5.5 — Anatomical localization

### Confirmed (no action)

- **Per-band cortical baseline rates** (after Wm/Unk drop): α 30/866 = **3.46%** (K=30), β 46/866 = **5.31%** (K=46), γ_l 75/866 = **8.66%** (K=75). N_total = 866. Per `_section5_5_verify.md`. ✓
- **γ_l ctx-lh-fusiform**: 13/38 trace-leaves, **3.95×** baseline, p_hyper = **5.2×10⁻⁶**. Bonferroni m=48 threshold = 1.04×10⁻³ → **survives Bonferroni**. ✓
- **Pat_02 anatomical density**: 18 of cohort's 38 fusiform contacts (47.4%); 9/13 trace-leaves (69%); per-patient fusiform trace rate **9/18 = 50.0%** vs Pat_03 25% (2/8) and Pat_13 25% (2/8). Pat_02 carries both the cohort's anatomically densest fusiform implant AND its highest-rate fusiform tracer; either departure (anatomical density without high rate, or high rate at typical density) would weaken the cell. Per `_section5_5_verify_3.md`. ✓
- **Pat_02 dropout at γ_l fusiform**: 4/20, 3.32×, p_hyper = 0.027. Fails Bonferroni m=48; surviving only at uncorrected p < 0.05. ✓
- **β anatomy uncorrected**: top three candidates Hippocampus (5/27, 3.49×, p=0.011), left fusiform (5/38, 2.48×, p=0.045), left superior temporal (6/53, 2.13×, p=0.055). None reaches Bonferroni. ✓
- **β under Pat_03 dropout**: Hippocampus and left fusiform collapse (Pat_03 carries 2/5 Hippocampus + 3/5 fusiform trace-leaves); only left superior temporal (no Pat_03 contribution) remains, with mild upward shift to p=0.035 — uncorrected only and well above α/m. ✓
- **α ctx-lh-parsopercularis**: 5/26 contacts, **5.55×** baseline, p_hyper = 1.4×10⁻³. Narrowly misses Bonferroni m=48. **Two-patient base**: **Pat_07 (3 leaves) + Pat_14 (2 leaves)** — NOT Pat_02 + Pat_13 as the original §5.5 prose claimed. ✓ (per `_section5_5_verify.md`)
- **Pat_07 dropout** at α parsopercularis: 2/20, 3.29×, p_hyper = 0.121, fails n_pat ≥ 2 filter. ✓
- **Pat_14 dropout** at α parsopercularis: 3/16, 5.04×, p_hyper = 0.019, fails n_pat ≥ 2 filter. ✓
- **Pat_02 dropout** at α parsopercularis: 5/20, 7.35×, p_hyper = 3.4×10⁻⁴, **passes Bonferroni** — paradoxically strengthens the cell because Pat_02's α implant has 6 parsopercularis contacts contributing zero trace-leaves and 4 α trace-leaves in other regions (Pat_02 is a dilution contributor, not a load-bearing patient at α parsopercularis). ✓
- **Pat_07 anti-aligned status at β/low-γ cross-probe-aggregate** AND **Pat_14 pro-aligned at α** (ρ_split = +0.240, third-strongest pro-trace value in cohort) → α parsopercularis is mixed-CTM-alignment cell (one anti-aligned + one pro-aligned). ✓ Per `_section5_5_verify_3.md` outcome (a).
- **Pat_07 / Pat_15 cortical β trace-leaf counts** (figure caption check): Pat_07 = 9 total, 1 cortical (parsopercularis); Pat_15 = 9 total, 3 cortical (rh-isthmuscingulate, rh-paracentral, rh-posteriorcingulate, one each). ✓
- **Pat_07 cross-band consistency at parsopercularis**: Pat_07's single β cortical trace-leaf and his 3 α cortical trace-leaves co-localize at ctx-lh-parsopercularis. ✓ Disclosed in §5.5 prose ("Pat07's single cortical β trace-leaf sits at left parsopercularis, the same region that anchors his α contribution: a patient-internal cross-band consistency at parsopercularis").

### Numerical corrections (action: writing agent)

- All §5.5 numbers in the current manuscript draft (e.g. `5/26, 5.55×, p=1.4e-3`, `13/38, 3.95×, p=5.2e-6`, `m=48` Bonferroni threshold ~1.04×10⁻³, baselines 3.46% / 5.31% / 8.66%) reproduce in the verifier reports. ✓
- The α parsopercularis attribution is **already corrected to Pat_07 (3) + Pat_14 (2)** in the current manuscript (verified at p.37 / p.38 of the PDF). ✓
- The Pat_14 cross-probe-aggregate alignment is correctly disclosed as "**Pat14 contributes the remaining two trace-leaves and is unambiguously pro-aligned at α: ρ_split = +0.240, the third-strongest pro-trace value in the α cohort, well above his own drift floor (ρ_drift = +0.040)**". ✓ Matches `_section5_5_verify_3.md` outcome (a).

### Framing rewrites (action: writing agent)

- **Region-spanning "γ_l localizes vs β/α distributes"** framing: per `_section5_5_verify_3.md` Section A, this framing is partially-correct-but-Gini-with-zeros-inverts-it. The honest reading is **γ_l = broad-and-peaked** (20 regions populated, fusiform peak at 4× = 34% rate), **β = broad-and-flat** (20 regions populated, no single peak above 3.5× baseline), **α = narrow-and-flat** (12 regions populated; 3 cells at 2-5× on 2-patient bases). The right scalar is rate-of-peak (γ_l 34.2% vs β 18.5% Hip vs α 19.2% parsopercularis), NOT Gini-with-zeros (which gives α 0.84 > γ_l 0.77 > β 0.74).

  The §5.5 closing prose currently says "**γ_l admits a Bonferroni-survived peak at left fusiform; β populates twenty cohort-eligible regions without any single peak clearing 3.5× baseline; α populates twelve cohort-eligible regions with three peaks at 2-5× baseline that all rest on two-patient bases**". This matches the corrected reading. ✓ Keep.

- **Pat_02 fusiform double-dependence disclosure**: the §5.5 prose ("**The cell rests on a disclosed double dependence on Pat02: he contributes 18 of the cohort's 38 cortical fusiform contacts (47.4% of the cohort fusiform implant) and 9 of the 13 trace-leaves (69% of the peak), with a per-patient fusiform trace rate of 50.0% — twice the rate of the next two contributors (Pat03 and Pat13, each 2/8 at 25.0%). The dependence is therefore not a single-patient artifact in disguise but the cohort's anatomically densest fusiform implant and its highest-rate fusiform tracer, and either departure would weaken the cell.**") matches `_section5_5_verify_3.md` Section B exactly. ✓ Keep.

### Caveats to add (action: writing agent) — IMPORTANT

- **NEW caveat (epileptic-contact pool fraction)**: per memory `epileptic_imcoh_universal.md`, the β cross-probe enrichment at Hippocampus and δ cross-probe enrichment are known biology (epileptic-zone contacts cluster in mesial temporal regions). For §5.5 honesty, disclose what fraction of each anatomy-headline contact pool is epileptic:

  - **Hippocampus** (β candidate, 5/27 trace-leaves): the cohort's 27 Hippocampus contacts are heavily concentrated in mesial-temporal epileptic-zone implants. Disclose: "Of the cohort's 27 Hippocampus contacts, X are flagged as epileptic-zone in their respective `epileptic_nodes` files; the β Hippocampus enrichment may partly reflect this epileptic-pool concentration and is read as a directional indication only at uncorrected p (does not survive m=48 Bonferroni)."
  - **Left fusiform** (γ_l Bonferroni cell, 13/38): the cohort's 38 fusiform contacts include Pat_02's 18 contacts. Disclose what fraction of those are epileptic-zone (likely lower than Hippocampus, but the question deserves disclosure).
  - **Left superior temporal** (β candidate, 6/53): same disclosure question.

  **Suggested insertion** at the head of §5.5 (or in the "anatomical localization layer" preamble):

  > "The cohort's anatomical contact pool is not uniform with respect to epileptic-zone designation. Hippocampus, fusiform, and superior-temporal regions disproportionately host epileptic-zone contacts in some patients; the per-region trace-leaf enrichments below are computed against the full cohort cortical contact pool without epileptic-zone exclusion. The fraction of each headline region's contacts flagged as epileptic-zone in the respective `epileptic_nodes` files is ⟨insert from data/raw/stereoeeg_patients/Pat_NN/epileptic_nodes_Pat_NN.csv⟩, and a sensitivity analysis excluding epileptic contacts from the anatomy pool is owed."

  This is a load-bearing disclosure, especially for the β Hippocampus result (where the epileptic-zone overlap is the largest in the cohort).

- **NEW caveat (γ_l fusiform Pat_02 cohort representation)**: already disclosed in current §5.5 prose. ✓

### Figure actions (coding agent)

- **Fig. 28** (anatomy bar chart per band + axial MNI projections per band): the new figure update by `_section5_5_verify_3.md` Section D produces `data/outputs/figures/section_5_lrg_trace/headline/anatomy_2d_disclosed.pdf` with green halos on Pat_07's three α parsopercularis trace-leaves, blue halos on Pat_14's two, parallel to the existing yellow halo on Pat_02's nine γ_l fusiform leaves. Caption update text provided in `_section5_5_verify_3.md` Section D.

  **Action**: confirm the manuscript Fig. 28 is the disclosed-version (`anatomy_2d_disclosed.pdf`), not the original `anatomy_2d.pdf`. Caption should mention all three halo classes (yellow Pat_02 / γ_l, green Pat_07 / α, blue Pat_14 / α).

- **Suggested**: add a same-figure-strip (or a separate panel) showing the **epileptic-zone overlay** at Hippocampus, fusiform, and superior-temporal contact pools — same axial MNI projection but with epileptic-zone contacts marked distinctly. This addresses the new caveat directly.

## §5.6 — Cross-phase module taxonomy

### Confirmed (no action)

- **Five-class taxonomy definitions** (per Tab 4 + `_section5_6_verify.md` Section A5): trace = `L_tt ∩ L_post`; reset = `L_pre ∩ L_post`; anchor = `L_pre ∩ L_tt ∩ L_post`; rearrange = `L_post − (strict_trace ∪ strict_reset ∪ strict_anchor)` with size ≥ 3; diffuse = leaves outside all four. Strict-intersection rules are **two-way** for trace/reset (NOT three-way over rspre-fragments — corrected in current manuscript). ✓
- **Matching rule**: **greedy argmax Jaccard** over destination subtrees within a 1.5× size band (option-i bandlimited greedy), with selection-time deduplication via 30% MAX_OVERLAP_FRAC and a containment guard `FRAG_FACTOR = 2.0`. NOT bipartite. Per `_section5_6_verify.md` A1/A2/A4. ✓
- **Cohort grand totals** (60 cells = 10 patients × 6 bands; KC λ=0): trace 107, reset 83, rearrange **560**, anchor 144. Rearrange is **3.9× the next-most-abundant class** (anchor) and **5.2× trace**. Per `_section5_6_verify.md` C3. ✓
- **Showcases**:
  - Pat_07 β trace (Fig. 29): 32 of 116 leaves, six modules, |min(n_trace, n_anchor)| = 6 = unique cohort top. Per `_section5_6_verify.md` B2. ✓
  - Pat_10 low-γ reset (Fig. 30): five reset modules, 17/113 strict-intersection leaves; ties Pat_13 at 5 reset clades; **Pat_13 has more reset leaves (29 vs 18) at low-γ**.
  - Pat_03 high-γ rearrange (Fig. 31): 12 rearrangement clades, 65/122 leaves; **tied at 12 with Pat_02 and Pat_06 cohort-wide**, but cleanest single-clade visual.
  - Pat_02 β anchor (Fig. 32): nine modules, 42/117 leaves.
  - Pat_15 β diffuse (Fig. 33): 52 of 118 leaves outside any of the four classes.
  - Pat_13 β composite (Fig. 34): all five classes simultaneously expressed (2 trace / 4 reset / 3 rearrange / 7 anchor / 36 diffuse of 119 leaves).
- **Trace strict-intersection rule corrected** per `_section5_6_verify.md` Section A5: "trace strict set is **L_tt ∩ L_post only**, NOT a three-way construction" — appears corrected in current manuscript (verified at p.40). ✓
- **B1 Pat_10 phrasing**: prose currently says "**Pat_10 carries five reset profiles, tied with Pat_13 by clade count**" — this is the corrected phrasing per the verifier. NOT "specialist". ✓
- **C1 Pat_03 phrasing**: prose says "**Pat_03 high-γ is the cohort's 1024 Hz sampling-rate outlier (§4.2); the dendrogram of D(τ) is unaffected by the sampling-rate asymmetry, and the cell is included for the cleanest pedagogical visual**" — appears corrected, the "tied at 12 with Pat_02 and Pat_06" is **NOT explicitly stated** in current prose but is implicit in "cleanest pedagogical visual" framing. ✓ Polish: explicitly note the tie to harden the cohort-rank claim.

### Numerical corrections (action: writing agent)

- **§5.6 trace strict-intersection rule**: confirmed corrected to two-way (`L_tt ∩ L_post`) in current manuscript. ✓
- **§5.6 cohort totals (107 / 83 / 560 / 144)**: confirmed. ✓
- **§5.6 matching rule**: prose says "for each source subtree of size s at one phase, the closest counterpart in the destination phase is found by greedy argmax of the Jaccard similarity of leaf sets, restricted to a destination size band [s/1.5, s × 1.5]; destinations are not removed across successive matches, so the same destination subtree can in principle match multiple sources" — exact match to `_section5_6_verify.md` A2. ✓

### Framing rewrites (action: writing agent) — TWO NEW STRUCTURAL ISSUES

#### Issue 1 — classes are NOT mutually exclusive at fixed scale

The matching pipeline allows **independent acceptance** at the source-phase level (trace's source = `L_tt`, reset's source = `L_pre`, anchor's source = `L_pre`, rearrange's source = `L_post`). The selection-time deduplication is via the 30% MAX_OVERLAP_FRAC test on the **cumulative used leaf set within a single class**, not across classes. So a leaf can in principle be in trace's strict L1 (`L_tt ∩ L_post`) AND simultaneously in anchor's L2 (`L_pre ∩ L_tt ∩ L_post`) — which is the structurally most-coherent case (a leaf that survived task into rsPost AND was also there in rsPre).

**Action**: insert a one-sentence definitional caveat after Tab 4:

> "The five strict-intersection leaf sets are **NOT mutually exclusive at fixed τ**: a leaf may simultaneously satisfy multiple class predicates (e.g., a leaf coherent in `L_pre ∩ L_tt ∩ L_post` is anchor *and* trace, since `L_tt ∩ L_post ⊂ L_pre ∩ L_tt ∩ L_post`). The class-overlap structure is informative — full taxonomic decomposition (per Pat_13 β composite, Fig. 34) makes this explicit. The diffuse class is whatever leaves do NOT satisfy any of the four predicates, and is by construction disjoint from the other four."

This is a load-bearing definitional caveat; without it, the §5.6 reader will assume a partition.

#### Issue 2 — anchor class has no anatomy-only baseline

The anchor class is defined by **coherent matching at all three phases under Jaccard ≥ 0.6**. There is no test of whether the anchor leaves are statistically over- or under-represented at any anatomical region (i.e. whether anchor leaves cluster at structural backbone regions like cingulate / paracingulate, or at high-degree hubs, or at same-probe pools). This is structurally analogous to the §5.5 anatomy test for trace leaves, but for anchor.

**Action**: in the §5.6 closing or §6.3 Outlook, flag explicitly:

> "The anchor class has no anatomy-baseline test in this section. A natural extension is the §5.5-style hypergeometric enrichment of anchor leaves against the cohort cortical-contact pool, by Desikan-Killiany region, with the same per-band Bonferroni-m=48 correction. The structural prediction is that anchor leaves over-cluster at structural-backbone regions (cingulate, paracingulate, supramarginal — high-degree hubs of the underlying anatomical connectivity)."

#### Issue 3 — τ is fixed at 1/λ_max (no multiscale class assignment)

The §5.6 taxonomy operates at **a single τ = 1/λ_max**. Class assignments at coarser τ (e.g. τ = 1/λ_2 or τ at the susceptibility peak) may differ — a leaf that is in trace at τ_max might be in rearrange at τ_2. This is the multiscale-class-assignment open question that §6.3 Outlook should explicitly raise.

**Action**: add to §6.3 Outlook (or §5.6 closing):

> "The five-class taxonomy is computed at τ = 1/λ_max only. Class assignments at coarser τ (where the propagator integrates over longer pathways) may shift; a τ-sweep version of the §5.6 taxonomy is the natural multiscale extension and would test whether the trace-class leaves at τ_max are also trace-class at τ ≫ τ_max, or whether the band-resolved imprint at the diffusion layer sharpens or dissolves with τ."

### Caveats to add (action: writing agent)

- See Issues 1, 2, 3 above. All three need new sentences.

### Figure actions (coding agent)

- **Fig. 29 / 30 / 31 / 32 / 33** (showcase per class): each is a 3-row × 3-column grid (dendrogram / D̂(τ_max) connectivity matrix / network panel) for the respective showcase patient. Strict-intersection leaves shaded; complement sorted by leaf-id. ✓ Per `kc_cross_phase_taxonomy.md`.
- **Fig. 34** (Pat_13 β composite, all five classes): 5-row × 3-column grid; each row = one class, each column = one phase; strict-intersection leaves highlighted within row, sorted by strict-intersection size (light to dark within row). ✓
- **Verify** the §5.6 figure caption clarifies "**leaves outside the strict-intersection of a class are still drawn at their actual dendrogram positions but no longer form a connected colored subtree**" — this is the visual convention per `_section5_6_verify.md` D. ✓ Already in caption.

## Section-5 — overall

### Numerical corrections summary

- §5.1: 89/180 strict argmax-at-extremum (verify the within-d=10 clarification language).
- §5.4: Grassmann formula inconsistency — replace Frobenius identity with √2-corrected version OR with sqrt(sum sin²θ).
- §5.5: confirm Fig. 28 = `anatomy_2d_disclosed.pdf` with all three halo classes.
- §5.6: confirm trace = `L_tt ∩ L_post` (two-way) corrected throughout.

### Framing rewrites summary

- §5.1: optional matched-strength surrogate forward.
- §5.4: lift the apple-to-oranges k caveat into the §5.4 head.
- §5.5: epileptic-contact pool fraction disclosure (load-bearing for Hippocampus and superior-temporal cells).
- §5.6: three new structural caveats — class non-disjointness, anchor anatomy-baseline absence, τ-fixed class assignment.

### Figure actions summary

- §5.5 Fig. 28 → confirm disclosed version with three halo classes.
- §5.5 — optional epileptic-zone overlay panel.
- §5.6 Fig. 29–34 — confirm captions match the strict-intersection / union conventions.

### Deferred / questions

- Should the §5.6 taxonomy be re-run with a coarser τ (τ = 1/λ_2 or τ at C(τ)-peak) to test multiscale class stability? This is a follow-up rather than a §5.6 fix.
- Should the anchor-class anatomy test be run for this revision pass, or deferred to future work?
- The §5.4 ~35-cell uncorrected p<0.05 count is reported as a directional indication. Should the joint-BH at m=48 (which already includes a Grassmann summary at one k per band) be expanded to m_extended = 48 + (Grassmann full k-grid additional cells) for an honest inflation analysis? Probably not — the §5.4 layer is already explicitly framed as not-controlled.
