---
name: writing-directive-per-pair-probe-battery-restructure
era: IMCOH_ABS_COHORT_N10
status: open
kind: writing-directive
scope: Methods restructure of the two per-pair trace subsections — substrate ρ_split^raw on A(B) (§sssec:methods_compare_rawfc) and cophenetic ρ_split on D_coph (§sssec:methods_compare_ctm) — into ONE shared generic-pipeline subsection + two thin object-specific instantiations, with a uniform control battery (matched-strength = GATE; drift-floor + cross-probe = CONTROLS) applied to BOTH probes.
created: 2026-06-01
companion: current LaTeX §sssec:methods_compare_rawfc + §sssec:methods_compare_ctm; compute homes audit_67 (raw obs+matched-strength), audit_74 (raw drift-floor, new 2026-06-01), audit_63 / matched_strength_surrogate_split_baseline (cophenetic matched-strength), continuous_trace_controls + h2e (cophenetic drift-floor + cross-probe)
target: methods writing agent
verified_against: data/audit CSVs as of 2026-06-01 pm (raw obs+MS = raw_fc_matched_strength/cohort_summary_all_bands.csv; raw drift = raw_fc_matched_strength/drift_floor_band_stats.md; raw cross-probe = raw_fc_matched_strength/cross_probe_band_stats.md [audit_74, GAP-1 closed]; cophenetic MS = matched_strength_surrogate_split_baseline/cohort_summary.csv [GAP-2 canonical]; cophenetic drift+xprobe = imcoh_continuous_trace/controls_band_stats.md). Equation↔code matches confirmed for ρ_split^raw (audit_67), both drift floors (h2e / audit_74), and the raw cross-probe restriction (ρ_full == published obs_rho to 1e-16).
---

# Per-pair trace probes — Methods restructure directive

**Status 2026-06-01 (pm): all three data gaps are now CLOSED** — see §5. The raw cross-probe control (GAP-1) is computed, GAP-2 reconciled (canonical = audit_63), GAP-3 resolved (report the matched-strength gate p, ignore the stacked-threshold label). Both probes now have a complete, single-source battery; the Results table (§4) can be typeset in full. The Methods restructure (§1–§3) was already independent of the gaps.

**Sequencing.** The Methods restructure (§1–§3, §6) is purely structural/definitional and was always **independent of the data gaps** in §5. Write the cross-probe control **definition** generically (object-agnostic, valid for both `A` and `D_coph`). The gaps that once blocked the Results per-band **numbers** (§4) are now **all closed** (see §5) — the raw cross-probe number exists, the cophenetic obs is reconciled to a single source, and the gate p-values are tabulated. So both the Methods prose **and** the Results table can be written in full now.

## 0. The decision (context you are writing toward)

Both per-pair probes carry the **full control battery** and each is reported as its *fully-controlled* cohort verdict:

- **`ρ_split^raw`** — per-pair shift correlation on the raw `|ImCoh|` adjacency `A(B)`.
- **`ρ_split`** — the same statistic on the LRG cophenetic communication distance `D_coph`.

The scientific point is the **contrast**, not either probe alone: at the **LRG cophenetic** scale the battery yields a **clean trace** (α/β clear the matched-strength *gate* and the drift/cross-probe *controls*); at the **raw substrate** scale it is a **mess** (numerically *larger* raw correlations that mostly fail the gate and clear no drift floor). **Raw-weak → LRG-clear is the intended message** — it is the empirical justification for the renormalization step, not a weakness to hide.

## 1. What is structurally wrong with the two current subsections

1. **Triplicated machinery.** The split-baseline shift definition, the drift-floor null, and the cross-probe restriction are written out *in full for `D_coph`* in §ctm, and the shift definition is written *again for `A`* in §rawfc. Adding the raw battery the way it is currently written would put each formula on the page **2–3 times**.
2. **Asymmetric controls.** §rawfc describes **no controls** (it frames `ρ_split^raw` as a "descriptive baseline"); §ctm describes drift + cross-probe. Under the decision in §0 this asymmetry is wrong — both probes get the full battery. (Keep the "conservative / blind to communication pathways" language for raw as *interpretive flavor*, but raw now also carries the battery.)
3. **Control hierarchy is mislabeled.** §ctm's cohort-claim sentence reads as if **drift + cross-probe are the gate**. They are **not**. Per the locked project rule, **matched-strength is the GATE** (the mandatory minimum null); **drift-floor and cross-probe are supplementary CONTROLS**. Matched-strength currently has **no home in either subsection** even though it is the primary test for both probes.
4. **Imprecise null wording.** "the drift floor is on the *same halved-data noise budget* as `ρ_split`" — not strictly true: `ρ_split` uses **full** task/rspost against half rspre, the drift floor uses **half** rspost. (Wording fix — your call on phrasing; flagged, not prescribed.)

## 2. Target structure — define once, instantiate twice

Replace the two subsections with **one shared subsection followed by two short instantiations.**

### 2a. Shared subsection — "Per-pair split-baseline trace probes"
Define everything **once** on a *generic* symmetric per-pair object `X^φ` (upper triangle, `N(N−1)/2` entries):

- **Half-split baseline.** Split rspre into two equal-duration halves → `X^{rspre_A}`, `X^{rspre_B}` (independent estimators on disjoint samples).
- **Shift vectors.** `Δ_task^X = X^taskt − X^{rspre_A}`, `Δ_rest^X = X^rspost − X^{rspre_B}` (task and rspost are single full-duration observations; the *independent* half-baselines remove by construction the shared-baseline term that would correlate the two maps trivially).
- **Trace scalar.** `ρ^X = ρ_S(Δ_task^X, Δ_rest^X)` over all `N(N−1)/2` pairs; `ρ^X > 0` is the trace direction. State the **Spearman rationale once** (rank/topology, orthogonal to the edge-strength backbone, robust to the heavy-tailed `|ImCoh|` shift distribution; the `D_coph`-specific tie argument is a one-liner deferred to the cophenetic instantiation).
- **Control battery (state the hierarchy explicitly):**
  - **GATE — matched-strength surrogate.** R=200 4-cycle ±δ strength-preserving rewiring of each phase matrix, `ρ^X` recomputed on surrogates, one-sided upper tail. *Why:* strips the degree/strength structure (the magnitude backbone that inflates any correlation of levels), so a surviving `ρ^X` is not explained by strength alone. **Mandatory minimum null — the cohort claim rests on this.**
  - **CONTROL — drift-floor null.** `ρ_drift^X = ρ_S(X^{rspre_B} − X^{rspre_A}, X^{rspost_B} − X^{rspost_A})` — rest-only, no task; one-sided paired Wilcoxon `ρ^X > ρ_drift^X`. *Why:* monotone within-session drift makes both displacements share a direction and would fake a positive `ρ^X`; the floor measures that no-task component.
  - **CONTROL — cross-probe restriction.** `ρ_xprobe^X` = `ρ^X` recomputed on the cross-probe pair subset only; non-degradation gate. *Why:* same-probe pairs carry residual short-range coupling. (Note for accuracy: under `imcoh_abs` the zero-lag component is killed by construction, so this is a **robustness control, not a de-biasing requirement** — see `probe_bias_critical`.)
- **Verdict logic, stated once:** the cohort claim requires the **matched-strength GATE**; `ρ_drift^X` and `ρ_xprobe^X` are reported **alongside as supplementary controls**.

### 2b. Substrate instantiation — `X = A(B)`
Two–three sentences: "instantiate the operator with `X` the raw `|ImCoh|` adjacency `A(B)`." Keep the conservative / "blind to the communication pathways the LRG propagator unfolds" interpretive framing. **No formula is restated** — `ρ_split^raw`, its drift floor, its cross-probe, its matched-strength surrogate are simply the generic objects with `X = A`.

### 2c. Cophenetic instantiation — `X = D_coph`
"instantiate the operator with `X = D_coph`." **Keep in full** the cophenetic-choice justification that is genuinely object-specific (τ_min = 1/λ_max single-scale propagator; continuous spectrum forbidding the standard τ-sweep; linkage hierarchy / `N−1` merge heights as the multiscale carrier; cophenetic image as the per-pair scale at which a pair coalesces). Add the one-line tie/monotone-rescaling note that motivates Spearman specifically on `D_coph`. Again **no battery formula is restated**.

## 3. Clarifying "the raw drift floor comes for free" (the point that was unclear)

`ρ_drift^raw` is **not a new equation**. It is the generic `ρ_drift^X` of §2a with `X` set to the raw adjacency `A`. Once the generic drift-floor (and cross-probe, and matched-strength) are written **once** in the shared subsection, the substrate subsection introduces **zero** new formulas — it only names its object (`X = A`). Same for the cophenetic subsection (`X = D_coph`). That is the entire anti-repetition mechanism.

## 4. Reporting target (Results, for orientation — not your section to write, but write the Methods toward it)

One results block per probe, 6 bands, columns: obs `ρ` (median) · **matched-strength GATE** (paired-Wilcoxon p) · drift-floor control (median floor + paired p) · cross-probe control (non-degradation). The narrative the Methods must support:

- **`ρ_split` (LRG cophenetic) = clean.** Matched-strength GATE clears **α (p=0.002), β (p=0.005)**; drift control clears α/β/low-γ (p=0.008/0.014/0.011); cross-probe non-degraded (cross ≈ split every band). Coherent trace at α/β.
- **`ρ_split^raw` (substrate) = mess.** Raw correlations are *larger* (β median +0.258 vs cophenetic +0.221) yet the battery is incoherent: the matched-strength GATE clears **α (p=0.014) and δ (p=0.042)** — *not* the LRG trace set — with the headline band β only **marginal (p=0.053; low-γ 0.053)**; the **drift control clears no band** (best low-γ p=0.057). The **cross-probe control is, by contrast, a clean pass** — raw cross-probe ≈ full at every trace band (β +0.258, α +0.159, low-γ +0.122; non-degradation p ≥ 0.17 except null-band γ_h p=0.037), so the raw weakness is **not probe geometry**. Net: the GATE fires on the wrong bands, the drift floor clears nothing, the headline band β fails the gate → the "mess" that motivates the LRG step. **The LRG↔raw contrast lives in the matched-strength GATE and the drift floor — both probes pass the cross-probe control**, so cross-probe is genuinely supplementary, not the discriminating axis.

## 5. Data inventory + GAPS — ALL CLOSED 2026-06-01

| probe | obs ρ + matched-strength (GATE) | drift-floor (CONTROL) | cross-probe (CONTROL) |
|---|---|---|---|
| `ρ_split^raw` | ✓ `raw_fc_matched_strength/cohort_summary_all_bands.csv` (audit_67) | ✓ `raw_fc_matched_strength/drift_floor_band_stats.md` (audit_74) | ✓ `raw_fc_matched_strength/cross_probe_band_stats.md` (audit_74) |
| `ρ_split` (LRG) | ✓ `matched_strength_surrogate_split_baseline/cohort_summary.csv` (audit_63) | ✓ `imcoh_continuous_trace/controls_band_stats.md` | ✓ `imcoh_continuous_trace/controls_summary.csv` (Run B) |

- **GAP-1 — CLOSED.** `ρ_xprobe^raw` built by extending `audit_74` with the Run-B `build_probe_mask` cross-probe restriction on the SAME `ρ_split^raw` per-pair shift vectors. Recomputed `ρ_full` matches audit_67's published `obs_rho` to `max|Δ| = 9.7e-17` over all 60 cells (anti-hallucination guard), so the cross-probe split is provably the only change. **Result: non-degradation passes for every trace band** — raw cross-probe ≈ full (β +0.258 vs +0.258, α +0.159 vs +0.158, low-γ +0.122 vs +0.126; non-degradation paired-Wilcoxon p ≥ 0.17 for δ/θ/α/β/low-γ; only null-band γ_h shows mild degradation p=0.037). The raw probe's weakness is the GATE + drift floor, **not** probe geometry.
- **GAP-2 — CLOSED. Canonical cophenetic `ρ_split` obs = audit_63 (`matched_strength_surrogate_split_baseline/cohort_summary.csv`).** Reason: the matched-strength GATE p-values are computed against *that* obs, so the headline obs and its gate must come from one computation. The continuous_trace Run A value (α +0.115) is the same statistic on an independently-estimated half-baseline used only for the drift/cross-probe controls; it agrees with audit_63 to 0.001 at β (+0.222 vs +0.221) and differs by 0.010 at α (+0.115 vs +0.105) — half-baseline estimation noise that flips no verdict. **Use audit_63 for every obs-ρ cell; never mix the two.**
- **GAP-3 — CLOSED. Report the matched-strength paired-Wilcoxon p as the gate; ignore the CSV `verdict` label** (its "separated/intermediate/also_positive" is the forbidden stacked-threshold compound — `cohort z>2 AND |surr median|<0.05 AND n_above≥8` — banned by `feedback_no_hardcoded_test_thresholds`). The single-source gate p-values are tabulated below.

### Canonical numbers (single source — typeset from these)

**Matched-strength GATE** — cohort paired-Wilcoxon p (one-sided, obs > per-patient surrogate median):

| band | `ρ_split` (LRG, audit_63) obs / p | `ρ_split^raw` (audit_67) obs / p |
|---|---|---|
| δ | +0.008 / 0.278 | +0.111 / 0.042 |
| θ | −0.040 / 0.722 | +0.117 / 0.138 |
| α | +0.105 / **0.002** | +0.158 / **0.014** |
| β | +0.221 / **0.005** | +0.258 / 0.053 |
| low-γ | +0.083 / 0.116 | +0.126 / 0.053 |
| high-γ | +0.000 / 0.246 | +0.111 / 0.188 |

GATE clears (p<0.05): **LRG = {α, β}** (the coherent trace set); **raw = {α, δ}** with β/low-γ marginal at 0.053 — different bands, headline β fails.

**Drift-floor CONTROL** — paired-Wilcoxon p (obs > within-session drift floor):

| band | LRG (continuous_trace Run C) | raw (audit_74) |
|---|---|---|
| α | **0.008** | 0.142 |
| β | **0.014** | 0.101 |
| low-γ | **0.011** | 0.057 (best raw band) |
| δ / θ / high-γ | 0.22 / 0.25 / 0.40 | 0.070 / 0.142 / 0.399 |

Drift clears: **LRG = {α, β, low-γ}; raw = none.**

## 6. Division of labor

- **Methods agent:** the restructure in §2 (shared subsection + two instantiations), the control-hierarchy fix in §1.3 (matched-strength = gate; drift + xprobe = controls), and the wording fix in §1.4.
- **Compute side (not the methods agent): DONE 2026-06-01.** GAP-1 (`ρ_xprobe^raw`, audit_74) computed; GAP-2/3 reconciled into the §5 canonical-numbers block. Nothing further blocks the Results table.
- **Sign convention everywhere:** positive = trace (`feedback_td_sign_convention`).
