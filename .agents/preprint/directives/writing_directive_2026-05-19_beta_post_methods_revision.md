---
name: writing-directive-beta-post-methods-revision
era: IMCOH_ABS_COHORT_N10
status: directive
kind: writing-agent-directive
date: 2026-05-19
target: manuscript β subsection (\ssec:results_beta) — full alignment with the revised Methods section
band: beta
source_of_truth:
  - .agents/preprint/bands/01_beta.md
  - .agents/preprint/locked/VERDICT_LEDGER.md (Decision 7: cluster_mass co-primary)
  - .agents/preprint/locked/CONTROLS.md
  - .agents/preprint/locked/ANATOMY_LEDGER.md
  - .agents/preprint/methods/methods_revision_2026-05-18_cophenet.md
  - writing_directive_2026-05-18_beta_paragraph.md (initial fixes)
  - writing_directive_2026-05-19_beta_red_paragraphs.md (Grassmann red + anatomy red)
methods_changes_this_directive_responds_to:
  - Grassmann band-level gate is now \(p_{\rm mass}(b) < 0.05\) on the cluster mass \(T_{\rm G}^{\ast}\) ALONE.
  - The longest contiguous-significant run \(L_{\rm obs}\) is dropped from the results paragraph and from the per-patient table (Methods retains the definition as a descriptive companion only; the results section does NOT cite it). [Updated 2026-05-19, second pass: user directive — drop contiguous-run entirely from the main paragraph and the main table.]
  - \(T_{\rm G}^{\ast}\) is the canonical and only band-level Grassmann scalar reported in the paragraph.
  - The phantom-surrogate permutation construction is explicit in Methods; the paragraph cites verdict + p-value only.
  - The single systematic robustness regime is epi-X. Patient-dropout sensitivities (Pat_15) are NOT in the methods battery; if cited in the paragraph they must be framed as supplementary sensitivity, not as a primary control.
---

# Writing directive — β paragraph, post-Methods-revision alignment (2026-05-19, second pass)

## Head

The Methods section has been revised, and on user directive (2026-05-19, second pass) the **longest contiguous-significant run statistic is now dropped from the results paragraph and from the per-patient main table entirely**. The β paragraph reports only the **cluster mass `T_G^*`** and its cohort gate `p_mass`. The Methods retains the definition of `L_obs` as a descriptive companion (Eq. `methods_Lobs`) — that is methods-side context — but the results paragraph does **not** cite `L_obs` or `p_LR`. The current β results paragraph still says `"both criteria"` and reports both p-values — that framing is now superseded. This directive enumerates the exact rewrites, notation updates, table-footnote updates, and number-fidelity checks needed to bring the β paragraph into agreement with the revised Methods + the cluster-mass-only directive. After this pass the β subsection is methods-aligned.

Three downstream items also need touching: (i) the table footnote that still says `T̄_G` is "placeholder pending cluster-mass-weighted band integral" — this directive resolves the per-patient Grassmann column choice cleanly; (ii) a residual `\patient{03} dropout` reference in the per-pair stability sentence that should not be there under the locked patient-dropout policy; (iii) one numerical discrepancy in `data/audit/grassmann_cluster_extent/cohort_summary.csv` between the stored `obs_cluster_mass_neglog10p` and the Methods equation, flagged here so the agent does not inline the at-risk value.

---

## Authority and constraint (unchanged from prior directives)

- Plain English first; numerical evidence appended at the close of each load-bearing claim.
- One p-value per load-bearing claim. Counts (7/10, 8/10, 9/10) carry the per-patient cohort picture as descriptive context.
- Methods stay in the Methods section. Do **not** re-state cluster-mass construction, phantom-surrogate permutation, the matched-strength algorithm, or the cophenet-vs-raw-D rationale. The paragraph cites the **verdict** and the **gate p-value**; methodology lives one section away.
- Patient-dropout policy: Pat_03 is **never** a sensitivity row. The only LRG-native sensitivity is `n=9` (drop Pat_15). Patient-dropout sensitivities are **outside the Methods control battery** — they are reported as **supplementary** robustness, not as a primary control. The paragraph may cite them; do not frame them as part of the locked CONTROLS battery.
- "Persistence" / "trace" interchangeable in the β subsection.
- "Inter-shaft" / "cross-probe" synonyms.

---

## Mandatory text rewrites

### Fix R1 — Grassmann paragraph: cluster-mass-only, drop longest-run statistic entirely

**Current text (the load-bearing two-sentence block):**

> "The contiguous block is large enough — and the aggregate strength of the band's significant cells high enough — to clear the empirical-null distributions built from the matched-strength surrogate ensemble on both criteria: longest-contiguous-run permutation \(p = 0.005\) and aggregate cluster-mass permutation \(p = 0.005\) (both at the \(1/(R+1)\) floor of the \(R = 200\) ensemble)."

**Problem.** Under the revised Methods + the 2026-05-19 user directive, only the **cluster mass** `T_G^*` is reported in the results paragraph. The longest-run statistic `L_obs` is **not** cited in the paragraph at all — neither as a co-gate nor as a descriptive companion. The current sentence reports both, in a "both criteria" framing that suggests a disjunctive gate, which is wrong on both counts (the gate is mass-only, and `L_obs` is methods-side context only).

**Required edit.** Drop the longest-run statistic from the paragraph entirely. Replace the two-sentence block with a single sentence that reports the cluster-mass verdict and nothing more.

**Suggested wording (replacement for the two-sentence block):**

> "The aggregate strength of the band's significant cells, summed as the cluster mass \(T_{\rm G}^{\ast}\) of Methods Eq.~\eqref{eq:methods_TGstar}, clears the empirical-null distribution built from the matched-strength surrogate ensemble at the band level (cluster-mass permutation \(p_{\rm mass} = 0.005\), at the \(1/(R+1)\) floor of the \(R = 200\) ensemble)."

This wording carries through the Methods discipline:
- `T_G^*` is named once and pointed to Eq.~\eqref{eq:methods_TGstar}.
- The cluster-mass `p_mass = 0.005` is identified as the **only** band-level statistic the paragraph reports.
- No mention of `L_obs`, `p_LR`, "longest contiguous run", "single long contiguous window", "both criteria", or any related framing. The Methods section is where the cluster mass is defined and where `L_obs` is introduced as a companion; the results section reports only the mass.

The earlier sentence that introduces the contiguous-significant-band geometry — `"the significant cells concentrating into a single contiguous band of intermediate mode counts (k ∈ [27, 55], 29 cells)"` — is also affected. Two cases:

**Case A (recommended).** Drop the geometric description entirely. The paragraph then says only "the cohort signal clears the matched-strength surrogate across roughly a third of the resolvable subspace dimensions" → "[Fix R2 keeps a single-`k` anchor if the agent wants one]" → "the aggregate strength of the band's significant cells clears the empirical null (\(p_{\rm mass} = 0.005\))". Clean and methods-aligned.

**Case B (if the agent wants to preserve some across-`k` shape information without using `L_obs`).** Replace "a single contiguous band of intermediate mode counts (k ∈ [27, 55], 29 cells)" with a more agnostic phrasing like "an intermediate-mode-count regime" or "across the intermediate \(k\) range". This gives the reader a sense that the significant cells cluster in `k`-range without naming the contiguous-run statistic. Avoid "29 cells" and "k ∈ [27, 55]" in the paragraph (those are the `L_obs` numbers) — those can go in a supplementary figure caption with the per-`k` cohort Wilcoxon profile.

**Recommendation: Case A.** Simpler; less risk of the reader inferring an undocumented contiguous-run claim.

Notation reminder:
- `T_G(k)` per-k statistic (Methods Eq.~\eqref{eq:methods_T_G}); `T_G^*` band-level cluster mass (Eq.~\eqref{eq:methods_TGstar}). Adopt consistently.
- `p_mass` lowercase, matches Methods. Do not write `P_mass`, `pₘₐₛₛ`, or `cluster-mass p-value`.

### Fix R2 — Drop "approximately three times" claim or anchor it explicitly

**Current text:**

> "with the strongest effect near \(k \approx 40\) — where the observed cohort median is approximately three times the surrogate cohort median in magnitude, with nine of ten patients individually below their own surrogate."

**Problem.** This sentence anchors the reader at a **single representative `k` cell**, which under the revised Methods is no longer the load-bearing claim (the load-bearing claim is the band-level cluster mass). The "approximately three times" ratio (actual `|obs|/|surr|` = 2.94 at `k = 40`) is descriptive only; the Methods does not single out any one `k`.

**Required edit.** Keep the sentence as **descriptive context** — it gives the reader a concrete per-k anchor inside the contiguous window. Reframe it to make clear this is a representative descriptive number, not a load-bearing statistic. Suggested:

> "The strongest per-\(k\) effect inside this window sits near \(k \approx 40\), where the cohort-median observed \(T_{\rm G}\) is approximately three times the surrogate cohort median in magnitude and nine of ten patients sit individually below their own surrogate at this \(k\) (descriptive only; the band-level verdict is on the cluster mass \(T_{\rm G}^{\ast}\) above, not on any single \(k\)-cell)."

The parenthetical is short and removes the risk that a reviewer reads the `k = 40` sentence as the load-bearing claim.

Alternative: drop the `k = 40` sentence entirely and replace it with a one-clause description of the across-`k` profile: *"with `T_G(k)` cohort-negative across the full contiguous-significant window."* The agent's call — both are defensible.

### Fix R3 — Remove the Pat_03 dropout reference from the per-pair stability sentence

**Current text (per-pair stability sentence):**

> "The cohort effect is stable across patient-subset restrictions: dropping the \qty{1024}{\hertz} \patient{03}, or dropping \patient{15} — the only \gls{lrg} anti-aligned patient (right-hemisphere-only implant, no epileptic contacts) — leaves the cohort median in the same narrow band (\(+0.221\) to \(+0.230\) across all restrictions; \patient{03} dropout \(p = 0.027\); \patient{15} dropout \(p = 0.010\)), so the per-pair trace is carried by the cohort as a whole rather than by any one patient."

**Problem.** Under the locked patient-dropout policy (`feedback_no_patient_dropout.md`, retired Pat_03 dropout am 2026-05-18), Pat_03 is **not a dropout sensitivity**. Pat_03 is a full cohort member treated identically at the analysis layer; the 1024 Hz sampling rate is absorbed at the config layer (`nperseg_for_fs(fs)`, `FS_OVERRIDES`) and does not propagate to analysis-level treatment. Including Pat_03 as a robustness row reintroduces the retired outlier framing.

**Required edit.** Remove the `\patient{03}` reference from this sentence. The narrow band statement (`+0.221` to `+0.230`) stays — it should now reflect the LRG-native `n = 9` restriction only (drop Pat_15). Suggested rewording:

> "The cohort effect is stable across the LRG-native robustness restriction (\(n = 9\) cohort excluding \patient{15} — the only \gls{lrg} anti-aligned patient, right-hemisphere-only implant, no epileptic contacts): the cohort median sits in the same narrow band (\(+0.221\) to \(+0.230\); \(p = 0.010\) at \(n = 9\)), so the per-pair trace is carried by the cohort as a whole rather than by any one patient."

If the agent wants two restriction comparisons in the sentence, the second sensitivity check that **is** in the methods is the **cross-probe (inter-shaft) restriction** — same statistic on the cross-probe subset of pairs, cohort median `+0.223`, `n_trace_xprobe = 8/10`. That comparison is already covered earlier in the same paragraph ("Restricting the correlation to inter-shaft pairs ... leaves the cohort median essentially unchanged at +0.223"), so doubling it here is redundant; just drop Pat_03 and stay with the Pat_15 sensitivity.

### Fix R4 — Frame Pat_15 dropout as a supplementary sensitivity, not a control

**Current text (Grassmann paragraph):**

> "Dropping \patient{15} sharpens rather than weakens the subspace claim — the cohort median at \(k = 40\) deepens from \(-0.489\) to \(-0.572\) (\(p = 0.002\)) when the single \gls{lrg} anti-aligned patient is excluded."

**Problem (minor).** The wording is fine prose-wise. The only issue is that under the revised Methods the **single systematic robustness regime is epi-X**, not patient dropout. Patient-dropout sensitivities are **not** part of the methods-defined control battery; they are post-hoc cohort-composition sensitivity. The paragraph as written can read as if Pat_15 dropout is a primary control like epi-X — it is not.

**Required edit.** Soften the framing to "single-patient sensitivity" rather than "control". Suggested:

> "As a single-patient sensitivity, dropping \patient{15} — the only \gls{lrg} anti-aligned patient — sharpens rather than weakens the subspace claim: the cohort median at \(k = 40\) deepens from \(-0.489\) to \(-0.572\) (\(p = 0.002\), \(n = 9\)) when this patient is excluded."

The point is the cohort-paired Wilcoxon `p = 0.002` is **on the per-k test at `k = 40`**, not on the cluster mass `T_G^*`. The sentence should make that explicit; otherwise the reader could read this as `p_mass = 0.002` (it isn't — the per-patient single-cell Wilcoxon at `k = 40` and the cohort-level cluster-mass permutation are different statistics).

### Fix R5 — Remove "more than an order of magnitude" claim from per-pair paragraph

**Current text (per-pair paragraph close):**

> "the observed cohort median exceeds the surrogate cohort median by more than an order of magnitude (\(+0.221\) versus \(+0.009\), an \(\approx\!23\times\) ratio)"

**Problem (minor).** "More than an order of magnitude" + "≈ 23×" are redundant. "Order of magnitude" is also imprecise — `23×` is just under 1.4 decades, on the boundary of how a careful physical-scientist reader interprets the phrase.

**Required edit.** Drop "more than an order of magnitude"; keep the numerical anchor:

> "the observed cohort median exceeds the surrogate cohort median by an \(\approx 23\times\) ratio (\(+0.221\) vs \(+0.009\))"

Optional alternative: "by more than twenty-fold (`+0.221` vs `+0.009`, a 23× ratio)" if the agent prefers prose.

### Fix R6 — LOCK the sign convention: POSITIVE = trace, everywhere in the table and text

**Locked convention (user directive 2026-05-19, third pass).** All trace-direction numerical entries in the main table and in the results paragraph are **POSITIVE**. The negative-`T_*` phase-triangle convention used in Methods Eqs. `methods_triangle` / `methods_T_G` is **methods-side notation only**; the **results section reports correlation / mass quantities directly**, where the sign of trace is naturally positive. This eliminates the sign whiplash between `ρ_split > 0 = trace` (correlation, paragraph) and `T_ctm < 0 = trace` (sign-flipped scalar, table).

**Concrete consequences:**

1. **Per-pair multiscale column (col. 3)** — rename from `T_ctm` to **`ρ_split^coph`** (the correlation directly). Sign convention: **positive = trace**. Cohort median **+0.221** (matched-strength) or **+0.222** (within-baseline) — both positive, both with `**` asterisks. Per-patient values **flip sign** vs the current `T_ctm` column: where Pat_02 currently reads `−0.507`, it should read **`+0.507`**, etc. (The current table caption already cites `+0.222` for `ρ_split` in the text — the table column was the only sign-flipped artefact; remove it.)

2. **Grassmann per-patient column (col. 6)** — rename to **`T_G^{*, p}`** (per-patient cluster mass; positive by construction, monotone with trace strength). Definition: per-patient analogue of the band-level `T_G^*`, computed as `Σ_{k : p_k^p < α} −log10 max(p_k^p, 1/(R+1))` where `p_k^p` is the per-patient empirical one-sided lower p-value at cutoff `k` against that patient's R = 200 matched-strength surrogate distribution (regularized at the empirical-null floor `1/(R+1) ≈ 0.005`, max contribution per cell ≈ 2.303). Sign convention: **non-negative; trace strength is monotone in mass**.

3. **Substrate `T_d^(d_S)` sign column (col. 8)** — keep as a `+/−` sign column, but **flip the displayed sign** so that **`+` = trace direction** (currently the column shows `T_d^(d_S)` sign where `−` = trace per the methods phase-triangle convention; flip to display `sign(−T_d^(d_S))` or just relabel as "trace dir."). Caption clarifies the flip. Pat_02 currently shows `−`; under the locked convention should show `+`. Pat_07/Pat_10/Pat_15 currently show `+` (anti); under the locked convention should show `−`.

4. **Cohort-row entries in cols. 3 and 6** — all positive, asterisks indicate the band-level Wilcoxon gate (per-pair) or the cluster-mass permutation gate (Grassmann).

5. **Paragraph text** — already reports `ρ_split^coph = +0.221` (positive). No change needed there. The Grassmann paragraph should report `p_mass = 0.005` only (no inline `T_G^*` value per the open-data flag).

**One inline sentence near the start of the results section** to fix the convention for the reader (suggested, agent's placement):

> "All trace-strength magnitudes are reported with a positive sign throughout: \(\rhosplit > 0\) (per-pair correlation), \(T_{\rm G}^{\ast} > 0\) (cluster-mass aggregate), and trace direction at the substrate marked `+` (Methods Eq.~\eqref{eq:methods_triangle} sign-flipped from \(T_d^{(d_S)}\) for cross-probe consistency)."

The phrase "sign-flipped from `T_d^(d_S)`" is the only place this clarification is needed. After that one sentence the reader carries forward the unified convention.

---

## Table updates (locked positive=trace convention; normalized `T_G^{*, p}` is the per-patient Grassmann column)

### Per-patient Grassmann column — `T_G^{*, p}` (normalized per-patient cluster mass; percentage of saturation)

**Definition (user directive 2026-05-19, fourth pass: normalize to [0, 1]).** Per-patient cluster mass divided by the maximum possible mass over the `k`-grid:

```
T_G^{*, p}(b) = ( Σ_{k : p_k^{p}(b) < α_k} −log10 max(p_k^{p}(b), 1/(R+1)) )
                       / [ K · log10(R+1) ]
                       ∈ [0, 1]
```

with:
- `p_k^{p}(b) = mean_{r ∈ {1, …, R}} 𝟙[T_G^{surr,r}(k; p) ≤ T_G^{obs}(k; p)]` the per-patient one-sided lower empirical p at cutoff `k` against that patient's R = 200 matched-strength surrogate distribution;
- `α_k = 0.05` the per-cell cluster-forming threshold;
- `K = 111` the size of the `k`-grid (`k ∈ {2, …, 112}`);
- `R = 200` the matched-strength surrogate ensemble size;
- regularization at the empirical-null floor `1/(R+1) ≈ 0.00498` so single-cell contributions are bounded by `log10(R+1) ≈ 2.303`.

**The denominator `K · log10(R+1) ≈ 255.65`** is the maximum cluster mass a patient could achieve if every `k`-cell hit the floor. `T_G^{*, p}(b) = 1` ⇔ saturated (every `k`-cell at the empirical floor); `= 0` ⇔ no significant cells.

**The definition is now in Methods** (per the parallel methods-agent directive `methods_directive_2026-05-19_TG_normalization.md`). The table caption references Methods Eq.~\eqref{eq:methods_TGstar_perpatient}; do **not** re-state the definition in the caption.

**Why this normalization (and not the raw mass):**
- **Readable**: values in `[0, 1]`, naturally rendered as percentages. No surprise units.
- **Bounded by construction**: the denominator is the principled empirical-null ceiling, not an arbitrary scale.
- **Verdict gate unchanged**: simultaneous normalization of `T_G^{*, p}` and the null distribution `mass_normalized^{null}` leaves the permutation p-value `p_mass` calibrated (see methods directive).
- **Sign convention**: non-negative by construction; consistent with the locked positive=trace convention throughout the table.

**Per-patient `T_G^{*, p}` values at β (computed; no recompute needed):**

| Patient | `T_G^{*, p}` | percentage | `n_sig, k` |
|---|---|---|---|
| Pat_02 | **0.962** | **96.2 %** | 108 |
| Pat_03 | **0.859** | **85.9 %** | 96 |
| Pat_05 | **0.766** | **76.6 %** | 87 |
| Pat_06 | **0.858** | **85.8 %** | 96 |
| Pat_07 | **0.530** | **53.0 %** | 62 |
| Pat_08 | **0.520** | **52.0 %** | 62 |
| Pat_10 | **0.612** | **61.2 %** | 71 |
| Pat_13 | **0.520** | **52.0 %** | 61 |
| Pat_14 | **0.695** | **69.5 %** | 78 |
| Pat_15 | 0.006 | 0.6 % | 1 |
| **Cohort median** | 0.653 | 65.3 % | 74.5 |
| **Cohort mean** | 0.633 | 63.3 % | 78.2 |

Bold values in the per-patient rows above indicate per-patient significance at the Grassmann probe under the locked rule (`n_sig_k ≥ 20`); Pat_15 (`n_sig_k = 1`) is not bold. Cohort summary rows are not bolded here (they are not per-patient cells).

Pat_02 at 96% essentially saturated (all 111 k-cells at empirical floor); Pat_15 at 0.6% essentially null (1 single marginal cell at `k = 43` — a noise-floor artefact, not a trace signal); cohort median 65% of saturation. The ranking and the LRG anti-aligned identification of Pat_15 are preserved from the raw mass; the normalized scale just makes the table instantly readable.

### Per-pair column — flip to `ρ_split^coph` (positive = trace; remove `T_ctm` sign-flip)

The current table column 3 reports `T_ctm = −ρ_split^coph` (negative = trace, cohort median `−0.221`). Under the locked positive-trace convention this becomes **`ρ_split^coph`** directly. Per-patient values flip sign (Pat_02 `−0.507` → `+0.507`, etc.).

**Per-patient `ρ_split^coph` values at β (matched-strength split-baseline realization, sign-flipped from current table):**

| Patient | `ρ_split^coph` (obs) | surr. median | `z` vs surr. |
|---|---|---|---|
| Pat_02 | **+0.507** | +0.096 | +6.27 |
| Pat_03 | **+0.373** | +0.021 | +9.16 |
| Pat_05 | **+0.491** | +0.051 | +7.32 |
| Pat_06 | **+0.211** | +0.018 | +3.71 |
| Pat_07 | **+0.230** | −0.001 | +4.47 |
| Pat_08 | **+0.502** | +0.040 | +8.76 |
| Pat_10 | −0.091 | −0.000 | −1.70 |
| Pat_13 | **+0.208** | +0.001 | +3.51 |
| Pat_14 | −0.049 | −0.003 | −1.00 |
| Pat_15 | +0.083 | −0.014 | +1.27 |
| **Cohort median** | +0.221 | +0.009 | (7/10 > own surrogate; cohort Wilcoxon `p = 0.005`) |

Bold values in the per-patient rows above indicate per-patient significance at the per-pair probe under the locked rule (`z ≥ 1.96` against own matched-strength surrogate); Pat_10, Pat_14, Pat_15 are not bold (z below threshold). Cohort summary row is not bolded here (not a per-patient cell).

Trace direction = **positive**. Pat_10 and Pat_14 are anti at this probe (`ρ_split^coph < 0`); the other 8/10 are in trace direction.

### Substrate column (col. 8) — kept as comparison reference; bold-significance gates LRG-trace classification

Locked decision (user 2026-05-19, fifth pass): col. 8 stays in the table, but it is **only a comparison reference** to the substrate analysis's sign-of-trace. The **LRG-trace classification gate** is read off cols. 3 and 6 via **bold-significance**.

**Bold-significance rule (per-patient):**
- Col. 3 (`ρ_split^coph`) per-patient value is rendered **bold** iff `z ≥ 1.96` against the patient's own matched-strength surrogate (one-sided `p < 0.05`).
- Col. 6 (`T_G^{*, p}`) per-patient value is rendered **bold** iff `n_sig_k ≥ 20` (well above the null expectation of ≈ 5.5 cells; the patient's normalized cluster mass clears the empirical-null floor with substantive amplitude).
- A patient is **LRG-trace** iff ≥ 1 of cols. 3 or 6 is bold. Col. 8 substrate sign does **not** participate in this gate.

**Col. 8 substrate sign-of-trace as comparison reference:**
- Computed as `sign(−T_d^(d_S))` (Methods Eq.~\eqref{eq:methods_triangle} sign-flipped for cross-probe positive-trace consistency).
- Display rule: `+` = substrate-trace, `−` = substrate-anti, `0` if `|T_d^(d_S)| < 0.01` (below substrate noise budget at β).
- Per-patient cohort agreement count: `7/10` (Pat_02, 03, 05, 06, 08, 13, 14 in trace direction; Pat_10, Pat_15 anti; Pat_07 near zero). **No asterisk** — col. 8 does not gate any verdict.

**Resulting table column count**: 8 data columns (Patient, N_ch, `ρ_split^coph`, surr. med., z vs surr., `T_G^{*, p}`, n_sig_k, `d_S`) + Notes = 9 columns total.

### Table caption — rewrite to reflect locked convention (definitions in Methods, caption slim, bold-significance + LRG-trace gate stated)

**Replacement caption (locked convention: positive = trace; normalized cluster mass `T_G^{*, p}` defined in Methods; bold-significance and LRG-trace gate explicit):**

> "Per-patient evidence at \(\beta\) for the two \gls{lrg}-layer probes on the \(n = 10\) cohort. Trace direction is **positive** across all numerical columns. Cols.\ 3--5: per-pair multiscale correlation \(\rhosplit\) on the cophenetic communication distance \(\Dcoph\) against R\,{=}\,200 matched-strength surrogates per patient (matched-strength split-baseline realization; cohort median \(+0.221\), of which the within-baseline realization gives \(+0.222\)). Cols.\ 6--7: per-patient normalized Grassmann cluster mass \(T_G^{\ast, p} \in [0, 1]\) (Methods Eq.~\eqref{eq:methods_TGstar_perpatient}) and the count \(n_{\rm sig, k}\) of \(k\)-cells passing the per-cell threshold \(\alpha_k = 0.05\) at the patient level. **Bold values** in cols.\ 3 and 6 indicate per-patient significance against the matched-strength surrogate (\(z \ge 1.96\) on \(\rhosplit\); \(n_{\rm sig, k} \ge 20\) on \(T_G^{\ast, p}\), well above the null expectation of \(\approx 5.5\) cells). Trace classification at the \gls{lrg} layer follows the rule: a patient is LRG-trace iff at least one of cols.\ 3 or 6 is bold. Col.\ 8: substrate analysis's sign-of-trace as comparison reference, computed as sign of \(-T_d^{(d_S)}\) (Methods Eq.~\eqref{eq:methods_triangle}, sign-flipped for cross-probe consistency; \(+\) = trace, \(-\) = anti, \(0\) for \(|T_d^{(d_S)}| < 0.01\)). Col.\ 8 does **not** gate LRG-trace classification."

The full definition of `T_G^{*, p}` lives in the Methods section (`methods_directive_2026-05-19_TG_normalization.md` instructs the methods agent to add it). The caption is slim by design — it references the Methods equation rather than re-stating it.

### Table footnote (a) — slim, references Methods

**Replace with:**

> "Normalized per-patient cluster mass; see Methods Eq.~\eqref{eq:methods_TGstar_perpatient}. Source data: `data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv` over the full \(k \in \{2, \ldots, 112\}\) grid. Cohort row reports the cohort median."

### Table footnote (b) — cohort row asterisks (no change in meaning)

Already correct: `*p<0.05`, `**p<0.01`, `***p<0.001`, one-sided cohort-paired Wilcoxon signed-rank test in the trace direction. Under positive=trace, the test is one-sided alternative='greater' (instead of 'less'), but the inferential content is identical.

### Table footnote (c) — Grassmann cohort row asterisk

**Replace with:**

> "Asterisks on the cohort row of cols.\ 6--7 are anchored to the band-level cluster-mass permutation gate \(p_{\rm mass}\) on \(T_{\rm G}^{\ast}\) (\(R = 200\) ensemble, \(1/(R+1)\) floor; see Methods Eq.~\eqref{eq:methods_TGstar} and Eq.~\eqref{eq:methods_cluster_p}). At \(\beta\), \(p_{\rm mass} = 0.005\)."

No mention of `p_LR` or `L_obs` in the table or its footnotes.

### Cohort row entries (normalized)

**Locked-convention cohort row:**

```
\multicolumn{1}{l}{\textbf{Cohort}} & \multicolumn{1}{c}{\textemdash}
  & \multicolumn{1}{c}{$+0.221^{**}$} & \multicolumn{1}{c}{$+0.009$} & \multicolumn{1}{c}{$7/10^{**}$\tnote{b}}
  & \multicolumn{1}{c}{$0.653^{**}$\tnote{c}} & \multicolumn{1}{c}{$9/10^{**}$\tnote{c}}
  & \multicolumn{1}{c}{$7/10$} & --- \\
```

- Col. 3: `+0.221**` (matched-strength cohort median of `ρ_split^coph`).
- Col. 4: `+0.009` (matched-strength surrogate cohort median).
- Col. 5: `7/10**` (patients individually above own surrogate at `ρ_split^coph`).
- Col. 6: `0.653**` (cohort median of normalized per-patient cluster mass `T_G^{*, p}`; render as 0.653 or 65.3% — pick one format and stay consistent through the row).
- Col. 7: `9/10**` (patients with `T_G^{*, p} > 0.01` — i.e., at least one significant `k` cell at the patient level beyond the floor; Pat_15 has only 1 cell at the floor → marginal, hence 9/10 strictly above floor).
- Col. 8: `7/10` (patients with substrate trace direction, no asterisk).

### Table column header renames

- `T_ctm` → `\rhosplit` or `\rho_{\rm split}^{\rm coph}`.
- `\bar{T}_G\tnote{a}` → `T_G^{*,p}\tnote{a}` (or `T_G^{*,p}` with footnote indicating it is the normalized form, agent's call — but the symbol used in the table must match the symbol the methods agent ends up locking).
- `surr.\ med.` (col. 4): per-patient surrogate median of `ρ_split^coph`. No rename.
- `z` vs surr. (col. 5): per-patient z of `ρ_split^coph` vs surrogate. No rename.
- `z` vs surr.\tnote{a} (col. 7): per-patient `n_sig_k`. **Rename** to `n_{\rm sig, k}` since the column under the new definition is a count, not a z-score.
- `d_S` (col. 8): rename to "trace dir." or keep `d_S` with caption flip. Agent's call.

### Per-patient row entries (locked convention; bold-significance; Notes = biology + "mixed" flag)

**Header order (locked, 9 columns):**

| Patient | N_ch | `ρ_split^coph` | surr. med. | z vs surr. | `T_G^{*,p}` | n_sig, k | `d_S` | Notes |
|---|---|---|---|---|---|---|---|---|
| Pat_02 | 117 | **+0.507** | +0.096 | +6.27 | **0.962** | 108 | + | --- |
| Pat_03 | 122 | **+0.373** | +0.021 | +9.16 | **0.859** | 96 | + | 1024 Hz |
| Pat_05 | 118 | **+0.491** | +0.051 | +7.32 | **0.766** | 87 | + | --- |
| Pat_06 | 115 | **+0.211** | +0.018 | +3.71 | **0.858** | 96 | + | --- |
| Pat_07 | 116 | **+0.230** | −0.001 | +4.47 | **0.530** | 62 | 0 | --- |
| Pat_08 | 120 | **+0.502** | +0.040 | +8.76 | **0.520** | 62 | + | --- |
| Pat_10 | 113 | −0.091 | −0.000 | −1.70 | **0.612** | 71 | − | mixed |
| Pat_13 | 119 | **+0.208** | +0.001 | +3.51 | **0.520** | 61 | + | --- |
| Pat_14 | 119 | −0.049 | −0.003 | −1.00 | **0.695** | 78 | + | mixed |
| Pat_15 | 118 | +0.083 | −0.014 | +1.27 | 0.006 | 1 | − | R-hemi only |

**Reading the table at a glance.** Patients with both LRG cells bold (Pat_02, 03, 05, 06, 07, 08, 13) are LRG-trace at both probes. Patients with only the Grassmann cell bold (Pat_10, Pat_14) are **mixed** — anti at `ρ_split^coph`, trace at `T_G^{*, p}` — flagged in Notes. Pat_15 has **neither cell bold** — the only patient with no LRG-significance, the lone LRG-anti-aligned anchor (right-hemisphere-only implant). Col. 8 substrate is comparison only: Pat_07's substrate `0` (near-zero) does **not** demote Pat_07 from LRG-trace status; Pat_10's substrate `−` does **not** demote Pat_10 from LRG-trace (carried by Grassmann); etc.

LaTeX implementation of bold: `{\bfseries VALUE}` inside each `S` column cell. With `\sisetup{detect-weight=true}` (already set), siunitx preserves the decimal alignment under bold. See `tables/beta_per_patient.tex` for the canonical implementation.

---

## Number fidelity check (do this before sending the draft back)

The agent must verify that the following claims survive byte-for-byte. If any number is changed from this directive, flag it and stop.

| Claim | Number | Source row in CSV |
|---|---|---|
| Cluster-mass cohort gate, β | `p_mass = 0.005` | `grassmann_cluster_extent/cohort_summary.csv`, `band=beta`, `cluster_p_cluster_mass = 0.00498` (round to 0.005) |
| Longest-run companion p, β | **DO NOT CITE** (methods-only context per 2026-05-19 user directive) | — |
| Longest contiguous-significant run length, β | **DO NOT CITE** (methods-only context per 2026-05-19 user directive) | — |
| Per-k effect-size anchor at `k = 40` | observed cohort median ≈ 3× surrogate cohort median in magnitude (specifically, `−0.489` vs `−0.166` → 2.94× ratio); 9/10 patients below own surrogate at this `k` | `grassmann_matched_strength_surrogate/cohort_summary.csv`, `band=beta, k=40` |
| Per-pair `ρ_split^coph` cohort median, matched-strength realization | `+0.221` (or `+0.222` for within-baseline realization — pick one and stick to it) | `matched_strength_surrogate_split_baseline/cohort_summary.csv`, `band=beta` and `ctm_triangle/cohort_summary.csv` |
| Per-pair surrogate cohort median | `+0.009` (matched-strength) | same |
| Per-pair effect-size ratio | `≈ 23×` (do not write "more than an order of magnitude") | derived |
| Per-pair cohort Wilcoxon p, matched-strength | `0.005` | same |
| Per-pair n_above surrogate | `7/10` | same |
| Per-pair drift-floor null Wilcoxon p | `0.014` | `ctm_triangle/cohort_summary.csv` |
| Cross-probe cohort median | `+0.223` | same (`rho_xprobe_median`) |
<!-- Per-pair BH-FDR q = 0.027 across 6 bands row RETIRED 2026-05-20: cross-band BH-FDR dropped per writing-agent feedback (no verdict use; per-band controls already gate). See feedback_no_unmotivated_bh_fdr.md + METHODS_AUDIT_ISSUES.md §B2. -->
| LRG-native robustness restriction (drop Pat_15) | `n = 9, ρ_split^coph cohort median = +0.230, p = 0.010` | `responses/2026-05-18_beta_n9_drop15_wilcoxon_reply.md` |
| Pat_15-dropped Grassmann at `k = 40` | cohort median `−0.572`, `p = 0.002`, `n = 9`, single-patient sensitivity (NOT a methods control) | same reply md |
| Epi-X broadened window | `k ∈ [21, 56]`, 36 cells; 29 persist + 7 emerge | `grassmann_epi_exclusion/cohort_summary.csv` + `sensitivity.csv` |
| β cophenet anatomy regions (7) | left isthmus cingulate, right rostral anterior cingulate, left parahippocampal, left entorhinal, right insula, right postcentral, left superior frontal | `anatomy_beta_cophenet/cohort_summary.csv` A1 q_BH<0.05 AND A3 p_emp<0.05 AND z_obs>2 |
| β Grassmann anatomy regions (7) | Hippocampus, left middle temporal, left superior temporal, left lateral orbitofrontal, right medial orbitofrontal, left insula, right rostral middle frontal | `anatomy_beta_grassmann/cohort_summary.csv` A3 p_emp<0.05 AND z_obs>2 |

---

## ⚠ Open data question — cluster mass numerical value (`T_G^*` itself)

**Inline value of `T_G^*(β)` is at risk.** The current implementation of the audit_70 cluster-extent script (`scripts/01_compute/audit/audit_70_grassmann_cluster_extent.py`, lines 154–167) defines `cluster_mass()` as the sum of `−log10 p_k` over **all** cells with `p_k < α`, matching the Methods Eq.~\eqref{eq:methods_TGstar} literally:

```python
def cluster_mass(p_values, alpha):
    sig = p_values < alpha
    return float(np.sum(-np.log10(np.clip(p_values[sig], 1e-300, 1.0))))
```

Recomputing on `data/audit/grassmann_cluster_extent/per_k_obs_p.csv` for β yields **`T_G^*(β) = 69.76` over 40 significant cells across `k ∈ [21, 22] ∪ [27, 55] ∪ [57, 61] ∪ {63, 73} ∪ [89, 90]`**.

The stored value in `cohort_summary.csv` is **`obs_cluster_mass_neglog10p = 52.97`**, which corresponds to the sum over **only the longest contiguous-significant cluster** `k ∈ [27, 55]` (29 cells, mass 52.97). This is **not what the current Methods equation says**.

Hypothesis: the cohort_summary CSV was generated by an earlier version of the script that summed over the longest contiguous cluster only, before the script was refactored to its current sum-over-all-sig-cells form (script mtime 12:08, CSV mtime 11:12 on 2026-05-19).

**Implication for the paragraph.** **Do not inline `T_G^*(β) = 52.97` or `T_G^*(β) = 69.76`.** Either is at risk of disagreeing with the Methods or the CSV. The cohort verdict `p_mass = 0.005` is robust either way (both the 52.97 and 69.76 observed values are at the `1/(R+1)` floor of the R=200 null distribution — null max is around 11.24, both observations sit well outside). The paragraph should cite `p_mass = 0.005` and skip the raw `T_G^*` value. The Methods text already explains what `T_G^*` aggregates without needing a specific number.

**Action item flagged for the user (not for the writing agent):** decide which of the two definitions stands, then re-run audit_70 if needed so the cohort_summary CSV matches the Methods equation. Until that is reconciled, the paragraph should **not** inline a numerical `T_G^*` value.

---

## Section structure (preserve, with one tweak)

The current four/five-paragraph layout is sound:

1. **Lead claim + β-uniqueness + diffuseness flag** — keep as-is, with R6 (trace-direction recap sentence) optionally added.
2. **`ρ_split^coph` results** — apply R3 (remove Pat_03 dropout reference) and R5 (drop "more than an order of magnitude"). Stability sentence rewrites to the LRG-native `n = 9` restriction only.
3. **Interpretation of the per-pair finding** — unchanged.
4. **Grassmann results + epi-X** — apply R1 (cluster-mass-only; drop the entire longest-run / contiguous-run sentence and references), R2 (k≈40 descriptive framing), R4 (Pat_15 as sensitivity, not control). Notation update `T_G^*` and `p_mass` per Methods. No `L_obs` / `p_LR` anywhere in the paragraph.
5. **Anatomy** — unchanged from prior directive (already written, methods-aligned).
6. **Two-probe synthesis closing** — unchanged.

After R1–R6 and the table updates, the β subsection is methods-aligned and ready for layout / figure pairing.

---

## Source-of-truth references

- `.agents/preprint/bands/01_beta.md` — full β brief; every number traces here.
- `.agents/preprint/locked/VERDICT_LEDGER.md` Decision 7 (2026-05-19) — cluster_mass co-primary; supersedes Decision 6 (longest-run primary). Methods now uses mass-alone gate; LR is descriptive.
- `.agents/preprint/locked/CONTROLS.md` lines 44–77 — current trace verdict definitions referencing `cluster_p_longest_run` and `cluster_p_cluster_mass` as co-primary. NOTE: under the new Methods the verdict gate is cluster_mass alone; CONTROLS.md may need a follow-up edit to align with the Methods (separate item — flag for the user, not for the writing agent).
- `.agents/preprint/locked/ANATOMY_LEDGER.md` — β anatomy verdict, locked 2026-05-19.
- `data/audit/grassmann_cluster_extent/cohort_summary.csv` — `p_mass`, `p_LR` values.
- `data/audit/grassmann_cluster_extent/per_k_obs_p.csv` — per-`k` observed Wilcoxon `p` values used in cluster mass and longest-run computation.
- `data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv` — per-patient per-`k` `T_G` values for the table column.
- `data/audit/ctm_triangle/cohort_summary.csv` — `ρ_split^coph` within-baseline, drift floor, cross-probe.
- `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` — `ρ_split^coph` matched-strength contrast.
- `data/audit/grassmann_epi_exclusion/cohort_summary.csv` + `sensitivity.csv` — epi-X 29 persist + 7 emerge → 36 cells at `k ∈ [21, 56]`.
- `responses/2026-05-18_beta_n9_drop15_wilcoxon_reply.md` — Pat_15-dropped sensitivity numbers.

## Revision history

- **2026-05-19** (this directive) — Bring the β paragraph in alignment with the revised Methods section. Six rewrites (R1–R6), table footnote (a)/(c) updates, table column rename `\bar{T}_G` → `T_G(k=40)`, cohort row Grassmann cells filled with `−0.489**` / `9/10**`. Flag `T_G^*` numerical-value discrepancy between Methods equation (sum over all sig cells, 69.76) and cohort_summary CSV (sum over longest contiguous cluster, 52.97); cite `p_mass = 0.005` and skip the raw `T_G^*` value until reconciled.
- **2026-05-19** (prior, `writing_directive_2026-05-19_beta_red_paragraphs.md`) — Un-red the Grassmann subspace + anatomy red paragraphs. Built on the 2026-05-18 directive.
- **2026-05-18** (initial, `writing_directive_2026-05-18_beta_paragraph.md`) — Pat_07 attribution, audit_70 cluster-extent gate, β-uniqueness sentence, p-values at end of claims, number fidelity. Per-patient table requested; anatomy deferred.

## What the agent should produce next

1. The revised β subsection LaTeX with R1–R6 applied.
2. The updated table `beta_per_patient.tex` with `T_G(k=40)` column filled per the per-patient values above, and footnote (a)+(c) updated.
3. **Do not inline `T_G^*` numerical value** — cite `p_mass = 0.005` only.
4. Flag back to the user any number that does not match the fidelity table above.

That closes the β subsection methodologically. After the agent's pass, the next move is figures.
