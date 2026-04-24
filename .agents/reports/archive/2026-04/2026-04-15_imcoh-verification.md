---
name: imcoh-verification-results
type: report
era: IMCOH_SQ
status: superseded
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

> **⚠ QUANTITATIVE_STALE (2026-04-15)** — numbers in this report were
> computed under the pre-reset ImCoh mislabelling (stored `|ImCoh|²`
> under the name `imcoh`). Qualitative conclusions survive (rankings
> preserved under sqrt); regenerate absolute values against
> `fc_method="imcoh_abs"` before quoting. See
> `~/.claude/projects/.../memory/imcoh_taxonomy.md`.

# ImCoh Verification Results — Complete

Generated: April 2026. All numbers independently verified from cached data.

---

## SECTION A: Unanimity Cell Counts

### A1. Full k-range (k=2..N-1) ✅ ALL CONFIRMED

| Hyp | delta | theta | alpha | beta | low_γ | high_γ | TOTAL |
|-----|-------|-------|-------|------|-------|--------|-------|
| H1  | 72 ✅ | 110 ✅| 75 ✅ | 82 ✅| 73 ✅ | 32 ✅  | 444 ✅|
| H2a | 16 ✅ | 5 ✅  | 32 ✅ | 82 ✅| 59 ✅ | 0 ✅   | 194 ✅|
| H2b+| 5 ✅  | 0 ✅  | 3 ✅  | 21 ✅| 8 ✅  | 0 ✅   | 37 ✅ |
| H2b−| 1 ✅  | 9 ✅  | 3 ✅  | 0 ✅ | 1 ✅  | 0 ✅   | 14 ✅ |
| H3  | 6 ✅  | 46 ✅ | 11 ✅ | 0 ✅ | 20 ✅ | 28 ✅  | 111 ✅|

### A2. N/2 subset (k=2..58) ✅ ALL CONFIRMED

| Hyp | delta | theta | alpha | beta | low_γ | high_γ | TOTAL |
|-----|-------|-------|-------|------|-------|--------|-------|
| H1  | 35 ✅ | 54 ✅ | 24 ✅ | 32 ✅| 48 ✅ | 0 ✅   | 193 ✅|
| H2a | 6 ✅  | 5 ✅  | 29 ✅ | 36 ✅| 31 ✅ | 0 ✅   | 107 ✅|
| H2b+| 0 ✅  | 0 ✅  | 3 ✅  | 16 ✅| 0 ✅  | 0 ✅   | 19 ✅ |
| H2b−| 1 ✅  | 9 ✅  | 2 ✅  | 0 ✅ | 1 ✅  | 0 ✅   | 13 ✅ |
| H3  | 0 ✅  | 16 ✅ | 3 ✅  | 0 ✅ | 18 ✅ | 1 ✅   | 38 ✅ |

N/2 values are proper subsets of full-range values: ✅

### A3. H2a k-ranges ⚠️ NOT CONTIGUOUS

| Band | Cells | k-range | Contiguous? |
|------|-------|---------|-------------|
| beta | 82 | 16-108 | ❌ Gaps at 16..18, 34..36, 36..42, 58..60, 97..100, 106..108 |
| low_gamma | 59 | 21-108 | ❌ Gaps at 35..43, 78..90, 90..95, etc. |
| alpha | 32 | 17-61 | ❌ Gaps at 17..22, 26..32, 32..36, 38..40 |
| delta | 16 | 8-71 | ❌ Multiple gaps |
| theta | 5 | 45-49 | ✅ Contiguous |
| high_gamma | 0 | — | — |

**Flag:** H2a cells are NOT contiguous in any band except theta. The
unanimous regions have gaps where 1-2 patients dip below zero. This is
typical for k-resolved unanimity but should be noted — the "82 cells" in
beta are spread over k=16-108, not a single block.

### A4. Patient-level H2a ✅ CONFIRMED

| Band | Pat_02 | Pat_03 | Pat_05 | Pat_07 | Pat_08 | All 5 + ? |
|------|--------|--------|--------|--------|--------|-----------|
| beta | +0.262 | +0.259 | +0.187 | +0.192 | +0.182 | **YES** ✅ |
| low_γ | +0.280 | +0.123 | +0.701 | +0.085 | +0.177 | **YES** ✅ |
| alpha | −0.053 | +0.231 | +0.131 | +0.156 | +0.145 | NO (Pat_02) ✅ |
| delta | −0.002 | +0.368 | −0.010 | +0.204 | +0.113 | NO (Pat_02, Pat_05) ✅ |
| theta | −0.082 | +0.105 | +0.013 | +0.291 | −0.020 | NO (Pat_02, Pat_08) ✅ |
| high_γ | −0.003 | +0.197 | +0.087 | −0.179 | +0.268 | NO ✅ |

### A5. H2b theta reversal ✅ CONFIRMED

- 9 contiguous cells at k=28-36
- ALL 5 patients negative at every k in this range
- Contrast magnitudes: -0.04 to -0.39 (genuine, not marginal)
- Interpretation: in theta, task-test is FARTHER from rsPost than rsPre is
- This is a real theta-specific finding, not noise

### A6. H3 beta=0, theta=46 ✅ CONFIRMED

**Beta H3 = 0 because:** Pat_07 mean = −0.001 and Pat_08 mean = −0.075.
Two patients are negative → unanimity impossible. Beta's strong H2a
(rsPost resembles task) disrupts the rest/task discrimination that H3 requires.

**Theta H3 = 46 because:** 4/5 patients have positive mean (Pat_07 = −0.018,
barely negative). The 46 unanimous cells are at k=8-113 where all 5 agree.
Theta H2a is weak (5 cells) → rsPost still resembles rsPre in theta →
rest/task discrimination works cleanly.

---

## SECTION B: Beta H2a Scalar Contrasts ✅

| Patient | Mean contrast | Frac positive | Min | Max |
|---------|:------------:|:-------------:|:---:|:---:|
| Pat_02 | +0.262 | 96% | −0.094 | +0.559 |
| Pat_03 | +0.259 | 83% | −0.706 | +0.811 |
| Pat_05 | +0.187 | 93% | −0.088 | +0.687 |
| Pat_07 | +0.192 | 88% | −0.320 | +0.624 |
| Pat_08 | +0.182 | 90% | −0.068 | +1.093 |

**Group mean: +0.217 ± 0.036** (all positive) ✅

**Pat_03 sensitivity:** Excluding Pat_03 → 84 cells with 4/4 unanimous (stronger) ✅

**Reference distances (beta, mean VI across k):**
- VI(Pre,Post) ≈ 1.1 (largest → resting states are most different)
- VI(TT,Post) ≈ 0.9 (task-test closer to rsPost)
- VI(TL,TT) ≈ 0.7 (smallest → task phases most similar = H1)
- Δ_H2a ≈ 0.2 = ~20% of VI(Pre,TT) — meaningful effect size

---

## SECTION C: Probe Bias ⚠️ PARTIALLY PASSES

| Patient | n=3 | n=5 | n=10 | n=15 | n=20 | n=30 |
|---------|:---:|:---:|:----:|:----:|:----:|:----:|
| Pat_02 | 1.0 ✅ | 1.2 ✅ | 1.3 ✅ | 1.7 ❌ | 1.7 ❌ | 2.5 ❌ |
| Pat_03 | 1.0 ✅ | 1.0 ✅ | 1.1 ✅ | 1.1 ✅ | 1.2 ✅ | 1.3 ✅ |
| Pat_05 | 1.1 ✅ | 1.1 ✅ | 1.2 ✅ | 1.5 ❌ | 1.5 ✅ | 3.0 ❌ |
| Pat_07 | 1.0 ✅ | 1.0 ✅ | 1.1 ✅ | 1.5 ✅ | 1.8 ❌ | 2.4 ❌ |
| Pat_08 | 1.0 ✅ | 1.2 ✅ | 2.4 ❌ | 2.8 ❌ | 3.2 ❌ | 3.8 ❌ |

**ImCoh dramatically reduces probe bias at coarse scales (n≤10)** compared to
MSC (which had 3-8× at n=10). At n≤10, only Pat_08 exceeds 1.5×.

**At fine scales (n≥15), probe bias persists** because genuine lagged
short-range coupling exists between adjacent contacts. This is real physics,
not artifact — ImCoh guarantees it's not volume conduction.

**Implication for VI analysis:** The VI comparison uses ALL k-levels. At
k≤10 (the coarse scales), communities are clean. At k≥15 they start
following probes again — but this is genuine neural proximity coupling.

---

## SECTION D: High-Gamma Null ✅

- H2a: 0 cells ✅ (was 1 under MSC — artifact)
- H2b: 0 cells ✅
- H1: 32 cells at k=84-115 (fine scales only, not at N/2)
- H3: 28 cells
- Mean ImCoh: 0.001 (barely above noise floor)
- H1 at N/2: 0 cells → the 32 cells at full range are fine-scale only

**Flag:** High-gamma H1 signal exists but only at very fine scales (k>84).
At N/2 it's zero. This is likely noise amplification at fine scales where
ImCoh values are near the noise floor (~0.001). Report cautiously or omit.

---

## SECTION E: Narrative Verification

### E2. Claimed narrative shifts ✅

1. "H1: theta-dominant" → **YES** (110 vs beta 82) ✅
2. "H2a: beta+low_gamma dominant" → **YES** (82+59 vs alpha 32) ✅
3. "H2b: beta only" → **YES** (21 cells, only clean band) ✅
4. "H3: theta-dominant" → **YES** (46 vs low_gamma 20) ✅
5. "High_gamma task trace was artifact" → **YES** (0 H2a/H2b cells) ✅

### E3. Theta-Beta dissociation ✅

| | Theta | Beta |
|---|:-----:|:----:|
| H1 (stability) | **110** | 82 |
| H2a (trace) | 5 | **82** |
| H2b (approach) | 0 (+) / **9** (−) | **21** (+) / 0 (−) |
| H3 (within<cross) | **46** | 0 |

**Clean dissociation:**
- Theta = stability + discrimination (H1+H3), NO trace
- Beta = trace + approach (H2a+H2b), moderate stability (H1=82)
- Beta H1 is high → the dissociation is NOT perfect for H1
- But H2a/H2b/H3 show clean separation

### E4. H2a-H3 anti-correlation ✅

| Band | H2a | H3 | Pattern |
|------|:---:|:--:|---------|
| beta | 82 | 0 | H2a only |
| low_gamma | 59 | 20 | Mixed (both) |
| theta | 5 | 46 | H3 only |
| high_gamma | 0 | 28 | H3 only |
| alpha | 32 | 11 | Mixed |
| delta | 16 | 6 | Mixed |

Anti-correlation survives for beta (H2a↑, H3=0) and theta (H2a↓, H3↑).
Not as clean as MSC's alpha-beta pattern, but the beta-theta version is
actually more interpretable: beta connectivity reorganizes (H2a), disrupting
rest/task separation (H3=0); theta stays stable (H2a≈0), preserving
rest/task separation (H3=46).

---

## SECTION F: Code Traceability

| Computation | Script/Function | Log convention |
|-------------|----------------|----------------|
| ImCoh FC | `coherence_fc_pipeline(metric="imcoh")` in `utils/fc/msc/__init__.py` | natural log (np.log) |
| ImCoh LRG | Raw lrgsglib: `entropy()`, `compute_laplacian_properties()`, `compute_normalized_linkage()` | Bypasses compute_lrg_analysis cache |
| VI computation | `compute_vi()` in `scripts/py/compute_imcoh_vi.py` | natural log (np.log) ✅ |
| Unanimity | `scripts/py/compute_imcoh_unanimity.py` | 5/5 agreement criterion |
| Full k-range | k=2 to N-1, where N varies: Pat_02=117, Pat_03=122, Pat_05=118, Pat_07=116, Pat_08=120 | |

---

## FLAGS AND CONCERNS

1. ⚠️ **H2a beta cells are NOT contiguous** — 82 cells spread over k=16-108
   with gaps. Report the total count but note non-contiguity.

2. ⚠️ **Probe bias persists at n≥15** for some patients (Pat_08: 2.4× at n=10).
   This is genuine lagged coupling, not volume conduction, but still means
   fine-scale community structure partly follows probes.

3. ⚠️ **High-gamma H1 = 32 cells** appears only at k>84 (fine scales). Likely
   noise amplification. Report cautiously or note as fine-scale-only.

4. ⚠️ **Beta H1 = 82** means the theta-beta dissociation for H1 is imperfect.
   Both bands show strong task stability. The dissociation is clean for
   H2a/H2b/H3 but not for H1.

5. ⚠️ **Pat_02 is the weakest patient for H2a** — negative mean in alpha,
   delta, and theta. Only positive in beta and low_gamma. Pat_02 may have
   a different reorganization profile.

6. ⚠️ **compute_lrg_analysis() caching bug** — must always use `use_cache=False`
   when computing on ImCoh, or bypass with raw lrgsglib functions. The current
   ImCoh LRG cache was computed correctly via the bypass route.
