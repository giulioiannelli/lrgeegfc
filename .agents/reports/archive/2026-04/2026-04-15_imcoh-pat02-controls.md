---
name: imcoh-pat02-and-controls
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

# Pat_02 Investigation + Three Controls

## Step 3: Pat_02

### Finding: Pat_02 is NOT globally weak — it's H2a-specific

Pat_02 ranks **4th or 5th (strongest)** in H1, H2b, and H3 across most bands.
It's specifically weak in **H2a only**, and only in low-frequency bands
(delta, theta, alpha). In beta H2a, Pat_02 is **rank 5/5 (strongest!)** at +0.262.

| Hypothesis | Pat_02 typical rank | Interpretation |
|-----------|:------------------:|---------------|
| H1 | 2-4/5 | Normal to strong |
| **H2a** | **1-2/5 in delta/theta/alpha** | Weak trace in slow bands |
| H2a beta | **5/5 (strongest)** | Strong trace in beta |
| H2b | 4-5/5 | Strong approach |
| H3 | 4-5/5 | **Strongest** in most bands |

**Conclusion**: Pat_02 has a different frequency profile for task trace —
strong in beta/low_gamma, absent in alpha/delta/theta. This is not a
global deficiency but a band-specific pattern. Since the paper's headline
claim is about **beta** H2a (where Pat_02 is the strongest patient),
Pat_02's alpha weakness doesn't undermine the main result.

### No clinical metadata found

Pat_02's data folder contains only standard files (timeseries .mat,
implant .xlsx/.csv, channel labels). No clinical notes, task performance
data, or epilepsy metadata beyond electrode positions.

### ImCoh magnitude: Pat_02 is not an outlier

Pat_02's mean ImCoh values are comparable to other patients across all
bands. No spectral or connectivity anomaly explains the H2a pattern.
LRG n_nodes (117) and threshold (0.88) are within normal range.

---

## Step 4a: Theta as H2a negative control ✅

**Theta shows NO H2a effect.** The fraction of patients with positive
theta H2a contrast fluctuates around 0.60 ± 0.18 across k-levels
(0.5 = pure chance). Only 27% of k-levels have >60% patients positive.

Compare: beta H2a has **88%** of k-levels with 5/5 patients positive.

**Theta H2b goes OPPOSITE:** mean fraction positive = 0.49 (below 0.5),
with 9 k-levels of 5/5 unanimous NEGATIVE. The theta task-approach effect
is reversed — task-test is FARTHER from rsPost than rsPre is, the opposite
of what happens in beta.

**Control verdict**: ✅ Clean negative control. Theta shows the task does
NOT leave a trace in theta-band connectivity, and the approach effect is
reversed. This strengthens the claim that beta's H2a is band-specific,
not a generic artifact.

---

## Step 4b: High-gamma as null control ⚠️ PARTIAL

**Not as clean as expected.** High-gamma shows:
- H1: mean +0.22, 83% positive — there IS some signal
- H2a: mean +0.08, 68% positive — weak but not zero
- H2b: mean +0.05, 60% positive — weak
- H3: mean +0.10, 77% positive — some signal

High-gamma is NOT a clean null. The mean ImCoh is ~0.001 (near noise floor),
but the VI contrasts aren't zero — they're just variable and don't reach
unanimity (0/0 H2a/H2b cells).

**Control verdict**: ⚠️ Partial null. High-gamma has high variance
(std 0.25-0.33) and no unanimity, but the mean contrasts are positive.
Use cautiously — better framed as "high variance prevents unanimity"
rather than "no effect."

---

## Step 4c: H2a-H3 anti-correlation ⚠️ WEAK

**Overall anti-correlation is not significant.** Spearman ρ = −0.37 (p=0.47)
across 6 bands — the right direction but only 6 data points.

The beta-theta contrast IS clean:
- Beta: H2a = +0.217 (strongest), H3 = +0.064 (weakest among positive)
- Theta: H2a = +0.061 (near zero), H3 = +0.145 (strongest)

But across all 6 bands, the anti-correlation is driven by this single
beta-theta pair. Low_gamma breaks it (H2a = +0.272 AND H3 = +0.103,
both high).

**At k=90**, the anti-correlation is very strong (ρ = −0.94, p = 0.005),
suggesting it emerges at fine scales. At coarser scales it's absent.

**Control verdict**: ⚠️ The beta-theta specific dissociation is real
(beta=trace, theta=discrimination). But the global H2a-H3 anti-correlation
across all bands is too weak for a general claim. Frame as a beta-theta
specific dissociation, not a universal principle.

---

## Summary for paper narrative

| Control | Verdict | Use in paper |
|---------|---------|-------------|
| Theta as H2a control | ✅ Clean | "Task trace is absent in theta (fraction positive ≈ chance)" |
| Theta H2b reversal | ✅ Clean | "Task approach reverses in theta (9 unanimous negative cells)" |
| High-gamma null | ⚠️ Partial | "High-gamma shows no unanimous effects" (avoid "null") |
| H2a-H3 anti-correlation | ⚠️ Weak globally | "Beta and theta show complementary roles" (not "anti-correlated") |
| Pat_02 weakness | ✅ Explained | "Band-specific: weak in alpha H2a, strongest in beta H2a" |
