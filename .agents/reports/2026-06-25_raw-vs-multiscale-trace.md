---
name: raw-vs-multiscale-trace
type: report
era: "IMCOH_ABS × COHORT_N10"
status: complete
created: 2026-06-25
---

# Raw per-edge vs multiscale cophenetic trace — head-to-head

**HEAD.** The multiscale cophenetic trace does **not** beat the raw per-edge
comparison on raw *magnitude* — raw `rho` is in fact LARGER than cophenetic
`rho_split` in 5 of 6 bands (the ultrametric step compresses, it does not
amplify). The multiscale tool wins on the two axes that actually carry the
manuscript's story: (i) **band-taxonomy** — cophenetic-space lets θ fall to ~0
(absent) and β rise to +0.22 (consistent), a 0.26-wide separation, whereas raw
sits every band on a generic positive floor (~0.11, range only 0.15) and never
expresses the absent-θ / consistent-β contrast; and (ii) **null-clearing at β**
— cophenetic clears its matched-strength null at the flagship band (paired
Wilcoxon p=0.005, 7/10) while raw does **not** (p=0.053, 6/10). SNR-gating is
present in BOTH spaces and is not a clean multiscale win. **Verdict: the
multiscale tool is non-redundant with raw on the taxonomy and on β-specific
null-clearing — exactly the two places the narrative leans on it — but it is
redundant-or-worse on bare trace magnitude. The raw comparison belongs in the
figure as the foil that proves the cophenetic taxonomy is a multiscale result,
not as evidence the trace is "stronger" in the hierarchy.**

**FRAMING NOTE (PI, 2026-06-25) — read axis 1 correctly.** "Raw larger" is NOT a
defeat; it is the POINT. Raw FC gives a blunt trace that is positive in *every*
band including θ (no real trace there) — so raw cannot separate genuine persistence
from noise, and its larger magnitude is exactly that triviality (an undifferentiated
overlap everywhere). The multiscale filter does not preserve this blob: it strips
the trivial component (magnitude shrinks), θ collapses to ~0, β stands out and
localizes. So the comparison shows **raw = trivial, noise-inseparable trace;
multiscale = the filter that refines it into differentiated, interpretable
structure.** Judge method-superiority by refinement/separability (axes 2a/2b),
never by magnitude (axis 1) — the magnitude drop is the filter working.

Compute: `scripts/01_compute/audit/audit_146_raw_vs_multiscale.py`.
Data: `data/audit/raw_vs_multiscale/`.

---

## Construction (identical except the distance space)

For each (patient, band), from the SAME inputs — `task_test` FC, `rest_post`
FC, and the two cached `rest_pre` half-FCs (`A`, `B`) — symmetrized, zero-diag,
clipped to [0,1]:

- **MULTISCALE** `rho_split` = `Spearman(D_tt − D_preA, D_post − D_preB)` over
  the condensed **cophenetic** distance vectors `D_x = canonical_cophenet(A_x)`
  (LRG ultrametric, τ=1/λmax, average-linkage). Cached: audit_63
  (`matched_strength_surrogate_split_baseline`), col `obs_rho`.
- **RAW** `raw_rho` = `Spearman(triu(A_tt) − triu(A_preA),
  triu(A_post) − triu(A_preB))` over the **raw edges**. Cached: audit_67
  (`raw_fc_matched_strength`), col `obs_rho`, an intentional near-clone of
  audit_63 with `cophenet` replaced by `squareform`.

The ONLY difference between the two is cophenetic-distance-space vs raw-edge-
space. Both CSVs are node-matched (60 cells, 0 `N_nodes` mismatch). I
recomputed 5 cells (Pat_06 β, Pat_02 α, Pat_08 δ, Pat_05 γ_l, Pat_03 γ_h)
from the half-FCs and the cached `obs_rho` reproduces bit-exactly (max |diff|
= 8.3e-17), so the head-to-head is built on certified inputs.
(`recompute_certificate.csv`.)

Both measures additionally carry their own R=200 4-cycle ±δ matched-strength
surrogate null (audit_63, audit_67), so "clears a null" is comparable across
the two spaces.

---

## Axis 1 — MAGNITUDE: does the hierarchy persist MORE than its edges?

Paired per (patient, band): `delta = rho_split − raw_rho`. **No.** Raw is the
larger of the two in most bands; cophenetic is never significantly larger.

| band | median multiscale | median raw | median Δ (ms−raw) | n(ms>raw)/10 | p(Δ>0, 1-sided) | p(Δ, 2-sided) |
|------|------:|------:|------:|:---:|------:|------:|
| delta | 0.008 | 0.111 | **−0.054** | 4 | 0.83 | 0.33 |
| theta | −0.040 | 0.117 | **−0.123** | 4 | 0.83 | 0.33 |
| alpha | 0.105 | 0.158 | −0.040 | 7 | 0.19 | 0.39 |
| beta | 0.221 | 0.258 | −0.034 | 6 | 0.40 | 0.80 |
| low_gamma | 0.083 | 0.126 | −0.038 | 3 | 0.78 | 0.44 |
| high_gamma | 0.000 | 0.111 | **−0.093** | 3 | 0.81 | 0.39 |

High-SNR-5 only (n=5; small-sample Wilcoxon, read the medians not the p): same
direction — median Δ negative or ≈0 in 5/6 bands (θ −0.30, β −0.04, γ_l −0.04,
γ_h −0.00; only α +0.08), none significant.

**Read.** Cophenetic distances are quantized to N−1 merge heights and lose the
fine magnitude that a per-edge difference retains; the average-linkage
ultrametric is a compression, so `rho_split ≤ raw_rho` in magnitude is the
expected behaviour. The "rho_split is just raw_rho in a compressed coordinate
system" alternative is **confirmed on the magnitude axis** — and it predicts
exactly this sign. The multiscale claim therefore cannot be "the trace is
stronger in the hierarchy". It must be one of the structural axes below.

## Axis 2 — STORY: is the 3-tier consistency taxonomy multiscale-specific?

The manuscript band picture is `consistent {β, γ_l} | patient-specific
{δ, α, γ_h} | absent {θ}`. **The cophenetic taxonomy is multiscale-specific;
the raw measure flattens it onto a positive floor.**

Per-band cohort median across the 6 bands (`band_taxonomy_raw_vs_multiscale.csv`):

| measure | δ | θ | α | β | γ_l | γ_h | range | min band | max band |
|---------|--:|--:|--:|--:|--:|--:|--:|:--:|:--:|
| **multiscale** | 0.008 | **−0.040** | 0.105 | **0.221** | 0.083 | 0.000 | **0.261** | θ | β |
| **raw** | 0.111 | 0.117 | 0.158 | 0.258 | 0.126 | 0.111 | 0.147 | γ_h | β |

- The multiscale measure lets **θ go NEGATIVE (−0.04)** and **γ_h to ~0** — it
  expresses "absent". It rises to **+0.22 at β** — it expresses "consistent".
  The absent→consistent span is 0.26.
- The raw measure floors every band at **≈ +0.11** (its MINIMUM, at γ_h, is
  still +0.11; θ is +0.12, indistinguishable from the others). Raw range is
  0.15, ALL positive. Raw does **not** reproduce the absent-θ contrast — there
  is a generic positive cross-phase edge-overlap in every band that the
  hierarchy filters out and the raw edges do not.

**Null-clearing per band** (paired Wilcoxon obs-vs-surrogate-median; from each
measure's own matched-strength cohort summary):

| band | MS p | MS n>surr | RAW p | RAW n>surr |
|------|----:|:---:|----:|:---:|
| delta | 0.278 | 4/10 | 0.042 | 7/10 |
| theta | 0.722 | 2/10 | 0.138 | 7/10 |
| alpha | 0.0020 | 5/10 | 0.0137 | 8/10 |
| **beta** | **0.0049** | **7/10** | **0.053** | 6/10 |
| low_gamma | 0.116 | 5/10 | 0.053 | 7/10 |
| high_gamma | 0.246 | 4/10 | 0.188 | 7/10 |

- At the **flagship β band, multiscale clears its matched-strength null
  (p=0.005, 7/10) and raw does NOT (p=0.053, 6/10)** — the cleanest single
  win: the β trace is a multiscale-specific, null-beating effect.
- HONEST counterweight: raw is *not* null-dead. Raw clears at δ (p=0.042) and
  is "separated" at α (p=0.014, 8/10) where multiscale also clears (α p=0.002).
  Raw carries a generic positive bias that survives matched-strength in the
  easy bands. So "multiscale clears a null and raw never does" is FALSE; the
  precise true statement is "**at β specifically, multiscale clears and raw
  falls just short, and only the multiscale measure reproduces the
  absent-θ/consistent-β taxonomy structure**".

## Axis 3 — SNR GATING: raw property or multiscale-specific?

High-SNR-5 (Pat_06,05,02,03,08) vs low-SNR-5 (Pat_15,14,07,13,10), Mann-Whitney
one-sided high>low (`snr_gating.csv`):

| band | MS hi−lo | MS p(hi>lo) | RAW hi−lo | RAW p(hi>lo) |
|------|------:|----:|------:|----:|
| delta | 0.357 | 0.075 | 0.310 | 0.111 |
| theta | 0.077 | 0.500 | 0.430 | 0.075 |
| alpha | 0.280 | 0.210 | 0.250 | **0.028** |
| beta | 0.408 | **0.008** | 0.553 | **0.004** |
| low_gamma | 0.650 | **0.004** | 0.519 | **0.004** |
| high_gamma | 0.650 | 0.048 | 0.531 | 0.075 |

**SNR-gating is present in BOTH spaces and is NOT a clean multiscale win.**
Detectability (high-SNR patients carry more trace than low-SNR) shows up in raw
edges too — at β and γ_l it is significant in BOTH (raw even marginally
sharper at β: hi−lo 0.55 vs 0.41). At α, raw separates (p=0.028) while
multiscale does not (p=0.21); at γ_h/δ both are marginal. The SNR→trace
relationship is therefore a **general property of the cross-phase FC signal**,
inherited by the hierarchy, not manufactured by it. The honest statement for
the detectability framing: SNR-gating is real in both, so it argues the
heterogeneity is detectability-driven *regardless of measure*; it does **not**
argue the hierarchy is what makes detectability visible.

---

## Verdict — on which axes does multiscale beat raw?

| axis | winner | basis |
|------|--------|-------|
| **1. Magnitude** | **RAW** (or tie) | raw_rho ≥ rho_split in 5/6 bands; median Δ ≤ 0; ultrametric compresses, never amplifies. The "compressed coordinate" alternative is confirmed here. |
| **2a. Band taxonomy** | **MULTISCALE** | cophenetic spans [−0.04 θ … +0.22 β] range 0.26 and expresses absent-θ; raw floors every band at ≈+0.11, range 0.15, never reproduces the taxonomy. |
| **2b. β null-clearing** | **MULTISCALE** | β clears matched-strength in cophenetic (p=0.005) but not raw (p=0.053) — but raw clears at δ/α, so it is a β-specific win, not a blanket one. |
| **3. SNR-gating** | **TIE / raw if anything** | present and significant in both at β,γ_l; raw separates α where multiscale does not. Detectability is measure-agnostic. |

**One-line overall verdict:** *The multiscale tool does NOT make the trace
larger — raw edges persist more in magnitude — but it is the ONLY measure that
expresses the absent-θ / consistent-β consistency taxonomy and the ONLY one
that clears the matched-strength null at the flagship β band; the raw per-edge
comparison is the foil that proves the taxonomy (not the magnitude) is the
multiscale result, and it belongs in the narrative as exactly that.*

### What this licenses / forbids in the manuscript
- ALLOWED: "the consistency taxonomy is a property of the LRG hierarchy, not of
  raw edges — raw FC overlap is positive and floored in every band including θ,
  so the absent-θ / consistent-β structure only emerges in cophenetic space"
  (axis 2a, strongest claim).
- ALLOWED: "at β the cophenetic trace clears matched-strength where the
  raw-edge analogue does not" (axis 2b, β-specific).
- FORBIDDEN: "the trace is stronger / better resolved in the hierarchy than in
  raw edges" (axis 1 — false; raw is larger). The earlier "resolution not
  amplification" intuition is correct as *compression*, not as a magnitude win.
- FORBIDDEN: "SNR-gating reveals the hierarchy's detectability advantage" (axis
  3 — gating is present in raw too; it is measure-agnostic).
- NOTE: the *downstream* multiscale claims (OFC localization, encoding/inference
  dissociation) are graph-structure claims raw edges cannot express by
  construction and are out of scope for this edge-level foil — they stand on
  their own audits (audit_83/103/110/112), not on this comparison.

## Caveats
- This contrasts cophenetic *trace magnitude / null-clearing / band pattern*
  against the per-edge analogue. It is the hardest like-for-like foil for the
  trace, and NOT a foil for community/path/hierarchy claims (those have no
  edge-level analogue). Absence of a magnitude win does not bear on them.
- Both nulls are the same R=200 matched-strength surrogate (independent per
  phase); the audit_63 coordinated-cross-phase-null question is unchanged here.
- The hiSNR-5 paired tests are n=5 (scipy small-sample warning); read the
  medians/effect direction, not the p-values, for that split.

## Outputs
- `data/audit/raw_vs_multiscale/per_patient_per_band.csv` — 60 cells:
  `rho_split_multiscale`, `raw_rho`, `delta`, `snr_tier`, per-measure null flags.
- `data/audit/raw_vs_multiscale/cohort_summary.csv` — axis-1 paired stats
  (full10 + hiSNR5).
- `data/audit/raw_vs_multiscale/band_taxonomy_raw_vs_multiscale.csv` — axis 2.
- `data/audit/raw_vs_multiscale/snr_gating.csv` — axis 3.
- `data/audit/raw_vs_multiscale/recompute_certificate.csv` — bit-exact certify.
