# H1-H4 on VI(k): post-ImCoh-reset authoritative results

**FC method:** `imcoh_abs` = <|ImCoh|>_f (Nolte-2004 signed ImCoh with abs-then-band-average per frequency bin; Jensen's inequality ordering respected).
**LRG cache:** `data/cache/imcoh_lrg/` (representative mtime unknown).  
**Patients (N=10):** Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15. Per-patient caveats: Pat_03 recorded at 1024 Hz (all others 2048 Hz) — kept as a documented outlier / negative control. Pat_10 resting phases ship 113 channels vs 116 in task phases, so cross-phase VI (rest↔task pairs) is skipped for Pat_10. Pat_13 is missing rest_pre (corrupt .mat) and Pat_14 is missing task_test (corrupt .mat) — any contrast requiring those phases is n/a for the affected patient.
**Bands:** delta, theta, alpha, beta, low_gamma, high_gamma.  
**k range:** 2 to N/2 per dendrogram (typically k ≤ 50).  
**Git:** `f51e0aa`  
**Report generated:** 2026-04-24 01:12  
**Source CSVs:** `data/reports/imcoh_vi/vi_raw_profiles.csv`, `data/reports/imcoh_vi/hypothesis_contrasts.csv`.

Signed contrast convention: **positive = hypothesis supported**. Unanimity at a (band, k) cell means every patient that contributes data at that cell has the same sign.  The number of contributing patients per cell is not constant across the cohort — the per-patient columns below show n/a where a patient lacks the phase(s) required by the contrast.

---

## H2 — Task trace (central hypothesis)

H2 captures whether a cognitive task leaves a persistent trace on the multiscale brain organization.  It is split into two complementary contrasts, both evaluated at every dendrogram cut k:

- **H2a (task trace / reset asymmetry):** `VI(rest_pre, rest_post) > VI(task_test, rest_post)` — post-task rest is closer to task_test than to pre-task rest.
- **H2b (approach vs exit asymmetry):** `VI(rest_pre, task_test) > VI(task_test, rest_post)` — entering the task reorganizes more than exiting it.

### H2a

> **⚠ H2a does not pass at n=9.** The unanimity table below shows
> 2/708 (band, k) cells unanimous positive — essentially zero signal
> above patient-level heterogeneity. The parametric Wilcoxon rigorous
> test (see `data/reports/imcoh_vi/rigorous_tests.md`) corroborates:
> no band passes FDR-BH at n=9, best case is δ q=0.195, θ has mean
> contrast ≈ 0 (r_rb = −0.02). δ and β had passed at n=5/n=7 era, but
> Pat_10 and Pat_13 carry strongly negative H2a contrasts that erase
> the effect at n=9. **The multiscale task-trace claim is carried by
> H2c (continuous ultrametric drift direction) and H2d (causal
> block-level persistence), not H2a.** H2a remains here for historical
> continuity of the VI-based framework but should be demoted to
> supplementary reporting in the paper. The load-bearing handoff
> document is `.agents/reports/MULTISCALE_TASK_TRACE_FOR_WRITING.md`.

| band | k cells | unan. + | unan. − | frac+ | frac− | mean | sd |  Pat_02 | Pat_03 | Pat_05 | Pat_06 | Pat_07 | Pat_08 | Pat_10 | Pat_13 | Pat_14 | Pat_15 |
|------|--------:|--------:|--------:|------:|------:|-----:|---:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| $\delta$ | 118 | 1 | 0 |  +0.8% |  +0.0% | +0.0888 | 0.2164 | +67.0% | +82.5% | +47.4% | +100.0% | +95.6% | +39.0% | +36.0% | +25.6% |   n/a | +43.1% |
| $\theta$ | 118 | 0 | 0 |  +0.0% |  +0.0% | +0.0111 | 0.2056 | +47.0% | +79.2% | +52.6% | +92.9% | +97.4% | +32.2% |  +8.1% |  +5.1% |   n/a | +10.3% |
| $\alpha$ | 118 | 1 | 0 |  +0.8% |  +0.0% | +0.0545 | 0.2467 | +26.1% | +93.3% | +81.0% | +89.4% | +78.9% | +61.0% | +16.2% | +19.7% |   n/a | +32.8% |
| $\beta$ | 118 | 0 | 0 |  +0.0% |  +0.0% | +0.0625 | 0.3003 | +88.7% | +64.2% | +81.9% | +95.6% | +82.5% | +77.1% |  +1.8% | +18.8% |   n/a | +21.6% |
| $\gamma_{\mathrm{l}}$ | 118 | 0 | 0 |  +0.0% |  +0.0% | +0.0818 | 0.3781 | +83.5% | +71.7% | +94.8% | +95.6% | +51.8% | +87.3% | +12.6% |  +7.7% |   n/a |  +9.5% |
| $\gamma_{\mathrm{h}}$ | 118 | 0 | 0 |  +0.0% |  +0.0% | +0.0919 | 0.3005 | +70.4% | +70.0% | +66.4% | +94.7% | +30.7% | +87.3% | +13.5% | +63.2% |   n/a | +55.2% |

**Overall H2a:** 2/708 (band, k) cells unanimous positive; 0/708 unanimous negative; 706 split.

### H2b

| band | k cells | unan. + | unan. − | frac+ | frac− | mean | sd |  Pat_02 | Pat_03 | Pat_05 | Pat_06 | Pat_07 | Pat_08 | Pat_10 | Pat_13 | Pat_14 | Pat_15 |
|------|--------:|--------:|--------:|------:|------:|-----:|---:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| $\delta$ | 118 | 0 | 0 |  +0.0% |  +0.0% | +0.0383 | 0.2396 | +96.5% | +67.5% | +69.8% | +97.3% | +18.4% | +33.1% | +36.9% | +67.5% |   n/a | +17.2% |
| $\theta$ | 118 | 0 | 0 |  +0.0% |  +0.0% | +0.0116 | 0.1911 | +80.0% | +67.5% | +28.4% | +86.7% | +64.9% |  +5.1% | +62.2% | +47.0% |   n/a | +11.2% |
| $\alpha$ | 118 | 0 | 0 |  +0.0% |  +0.0% | +0.0274 | 0.2233 | +80.9% | +84.2% | +65.5% | +100.0% |  +7.9% | +28.8% | +53.2% | +94.0% |   n/a |  +6.9% |
| $\beta$ | 118 | 0 | 0 |  +0.0% |  +0.0% | +0.0812 | 0.2475 | +94.8% | +75.8% | +90.5% | +92.0% | +22.8% | +65.3% | +29.7% | +69.2% |   n/a | +14.7% |
| $\gamma_{\mathrm{l}}$ | 118 | 1 | 0 |  +0.8% |  +0.0% | +0.1287 | 0.3614 | +87.0% | +84.2% | +93.1% | +91.2% | +13.2% | +66.9% | +49.5% | +74.4% |   n/a |  +2.6% |
| $\gamma_{\mathrm{h}}$ | 118 | 2 | 0 |  +1.7% |  +0.0% | +0.0723 | 0.3297 | +62.6% | +69.2% | +87.1% | +90.3% |  +7.9% | +87.3% | +51.4% | +53.0% |   n/a | +22.4% |

**Overall H2b:** 3/708 (band, k) cells unanimous positive; 0/708 unanimous negative; 705 split.

---

## H1 — Task stability

**Statement:** `VI(task_learn, task_test) < mean(VI over other phase pairs)`. Positive contrast = the task-learn/task-test hierarchy is closer than the average phase pair at that (band, k).

| band | k cells | unan. + | unan. − | frac+ | frac− | mean | sd |  Pat_02 | Pat_03 | Pat_05 | Pat_06 | Pat_07 | Pat_08 | Pat_10 | Pat_13 | Pat_14 | Pat_15 |
|------|--------:|--------:|--------:|------:|------:|-----:|---:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| $\delta$ | 118 | 44 | 0 | +37.3% |  +0.0% | +0.2113 | 0.2379 | +89.6% | +100.0% | +76.7% | +100.0% | +79.8% | +85.6% | +95.5% | +70.9% |   n/a | +99.1% |
| $\theta$ | 118 | 64 | 0 | +54.2% |  +0.0% | +0.2684 | 0.2548 | +99.1% | +96.7% | +100.0% | +100.0% | +99.1% | +94.9% | +75.7% | +63.2% |   n/a | +99.1% |
| $\alpha$ | 118 | 81 | 0 | +68.6% |  +0.0% | +0.2972 | 0.2459 | +94.8% | +92.5% | +78.4% | +100.0% | +93.0% | +99.2% | +96.4% | +94.9% |   n/a | +99.1% |
| $\beta$ | 118 | 35 | 0 | +29.7% |  +0.0% | +0.2722 | 0.2690 | +99.1% | +99.2% | +97.4% | +92.0% | +98.2% | +50.0% | +85.6% | +90.6% |   n/a | +92.2% |
| $\gamma_{\mathrm{l}}$ | 118 | 38 | 0 | +32.2% |  +0.0% | +0.3172 | 0.3395 | +95.7% | +100.0% | +96.6% | +99.1% | +91.2% | +63.6% | +89.2% | +45.3% |   n/a | +98.3% |
| $\gamma_{\mathrm{h}}$ | 118 | 2 | 0 |  +1.7% |  +0.0% | +0.1639 | 0.2594 | +87.8% | +95.0% | +100.0% | +92.0% | +83.3% | +22.0% | +23.4% | +72.6% |   n/a | +98.3% |

**Overall H1:** 264/708 unanimous positive, 0/708 unanimous negative.

---

## H3 — Within-modality similarity exceeds cross-modality

**Statement:** `mean_VI(within) < mean_VI(cross)`, with within = `{(rest_pre,rest_post), (task_learn,task_test)}` and cross = the four rest-to-task pairs.

| band | k cells | unan. + | unan. − | frac+ | frac− | mean | sd |  Pat_02 | Pat_03 | Pat_05 | Pat_06 | Pat_07 | Pat_08 | Pat_10 | Pat_13 | Pat_14 | Pat_15 |
|------|--------:|--------:|--------:|------:|------:|-----:|---:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| $\delta$ | 118 | 0 | 0 |  +0.0% |  +0.0% | +0.0560 | 0.1102 | +86.1% | +95.8% | +45.7% | +90.3% |  +8.8% | +77.1% | +93.7% | +96.6% | +35.0% | +84.5% |
| $\theta$ | 118 | 6 | 0 |  +5.1% |  +0.0% | +0.1170 | 0.1499 | +99.1% | +98.3% | +96.6% | +92.9% | +32.5% | +91.5% | +90.1% | +94.0% | +24.8% | +100.0% |
| $\alpha$ | 118 | 22 | 1 | +18.6% |  +0.8% | +0.1191 | 0.1519 | +99.1% | +85.0% | +56.9% | +96.5% | +39.5% | +89.8% | +97.3% | +97.4% | +47.0% | +99.1% |
| $\beta$ | 118 | 0 | 0 |  +0.0% |  +0.0% | +0.1247 | 0.1735 | +100.0% | +99.2% | +92.2% | +74.3% | +56.1% |  +5.9% | +98.2% | +95.7% | +82.9% | +95.7% |
| $\gamma_{\mathrm{l}}$ | 118 | 30 | 0 | +25.4% |  +0.0% | +0.1588 | 0.1667 | +98.3% | +100.0% | +56.9% | +97.3% | +79.8% | +66.1% | +95.5% | +93.2% | +85.5% | +93.1% |
| $\gamma_{\mathrm{h}}$ | 118 | 13 | 0 | +11.0% |  +0.0% | +0.0698 | 0.1304 | +96.5% | +81.7% | +99.1% | +76.1% | +66.7% | +31.4% | +99.1% | +72.6% | +70.9% | +78.4% |

**Overall H3:** 71/708 unanimous positive, 1/708 unanimous negative.

---

## H4 — Frequency gradient of reorganization

**Statement:** across patients, is there a consistent band ordering by reorganization strength (mean cross-pair VI)?  Tested via Kendall's coefficient of concordance W (over patients as judges, bands as items) and pairwise Spearman ρ of per-patient band rankings.

### Global (k-averaged) concordance

- **Kendall's W = 0.231** (n_judges = 10, n_items = 6)
- Friedman χ²(5) = 11.54, p = 0.0416
- Interpretation: **weak agreement** (W < 0.3 weak, < 0.5 moderate, < 0.7 good, ≥ 0.7 strong).

**Consensus ranking (lower mean_rank = more reorganized):**

| band | mean rank | median rank | mean cross-VI |
|------|----------:|------------:|--------------:|
| $\gamma_{\mathrm{h}}$ | 2.60 | 1.5 | 1.0611 |
| $\alpha$ | 2.70 | 2.0 | 1.0384 |
| $\theta$ | 2.80 | 3.0 | 1.0246 |
| $\delta$ | 4.00 | 4.5 | 0.9497 |
| $\beta$ | 4.30 | 4.0 | 0.8938 |
| $\gamma_{\mathrm{l}}$ | 4.60 | 4.5 | 0.8786 |

**Per-patient top-3 most-reorganized bands:**

- **Pat_02:** $\gamma_{\mathrm{h}}$ > $\alpha$ > $\theta$
- **Pat_03:** $\gamma_{\mathrm{h}}$ > $\delta$ > $\beta$
- **Pat_05:** $\theta$ > $\gamma_{\mathrm{h}}$ > $\alpha$
- **Pat_06:** $\gamma_{\mathrm{h}}$ > $\gamma_{\mathrm{l}}$ > $\theta$
- **Pat_07:** $\theta$ > $\alpha$ > $\delta$
- **Pat_08:** $\alpha$ > $\theta$ > $\beta$
- **Pat_10:** $\gamma_{\mathrm{h}}$ > $\alpha$ > $\delta$
- **Pat_13:** $\alpha$ > $\gamma_{\mathrm{h}}$ > $\theta$
- **Pat_14:** $\gamma_{\mathrm{h}}$ > $\delta$ > $\alpha$
- **Pat_15:** $\theta$ > $\alpha$ > $\beta$

### k-resolved concordance

Kendall's W and pairwise Spearman ρ at each dendrogram cut k:

| k | W | χ² | p | mean ρ | min ρ | max ρ |
|---|------:|------:|------:|-------:|------:|------:|
| 2 | 0.181 | 9.06 | 0.107 | +0.092 | -0.886 | +0.986 |
| 7 | 0.293 | 14.63 | 0.012 | +0.214 | -0.714 | +0.943 |
| 12 | 0.144 | 7.20 | 0.206 | +0.049 | -0.771 | +0.829 |
| 17 | 0.089 | 4.46 | 0.486 | -0.012 | -0.886 | +0.943 |
| 22 | 0.085 | 4.23 | 0.517 | -0.017 | -0.943 | +0.943 |
| 27 | 0.126 | 6.29 | 0.279 | +0.029 | -0.829 | +0.886 |
| 32 | 0.150 | 7.49 | 0.187 | +0.055 | -0.943 | +0.829 |
| 37 | 0.257 | 12.86 | 0.025 | +0.175 | -0.943 | +0.829 |
| 42 | 0.277 | 13.83 | 0.017 | +0.196 | -0.886 | +0.943 |
| 47 | 0.309 | 15.43 | 0.009 | +0.232 | -0.771 | +0.886 |
| 52 | 0.269 | 13.43 | 0.020 | +0.187 | -0.886 | +0.943 |
| 57 | 0.384 | 19.20 | 0.002 | +0.316 | -0.600 | +0.943 |
| 62 | 0.402 | 20.11 | 0.001 | +0.336 | -0.486 | +0.943 |
| 67 | 0.371 | 18.57 | 0.002 | +0.302 | -0.600 | +0.886 |
| 72 | 0.377 | 18.86 | 0.002 | +0.308 | -0.771 | +1.000 |
| 77 | 0.406 | 20.29 | 0.001 | +0.340 | -0.771 | +0.943 |
| 82 | 0.488 | 24.40 | 0.000 | +0.431 | -0.771 | +1.000 |
| 87 | 0.472 | 23.60 | 0.000 | +0.413 | -0.600 | +0.943 |
| 92 | 0.371 | 18.57 | 0.002 | +0.302 | -0.543 | +0.943 |
| 97 | 0.327 | 16.34 | 0.006 | +0.252 | -0.657 | +0.943 |
| 102 | 0.322 | 16.11 | 0.007 | +0.247 | -0.714 | +0.943 |
| 107 | 0.264 | 13.20 | 0.022 | +0.182 | -0.829 | +0.943 |
| 112 | 0.087 | 4.33 | 0.503 | -0.015 | -0.928 | +0.943 |
| 117 | 0.070 | 2.10 | 0.836 | -0.116 | -0.886 | +0.841 |

---

## Notes

- **Pat_03** was recorded at 1024 Hz (all others at 2048 Hz); its MSC-era FC is ~3× denser than the rest of the cohort and unstable. The ImCoh pipeline uses `nperseg_for_fs(fs)` (2048 for Pat_03, 4096 for the others) so the 2-second segmentation is preserved. Pat_03 remains in the tables as a negative control: deviations from the cohort are expected there and do not invalidate a group-level result.
- **Pat_06** was completed with task_learn + task_test on 2026-04-22 and now contributes to all 4-phase analyses.
- **Pat_10** resting recordings carry 113 channels vs 116 in the task phases and in `channel_labels.csv` (3 channels dropped by the vendor at resting). Cross-phase VI is skipped for Pat_10 on any pair that mixes rest and task (n_nodes mismatch).
- **Pat_13** rest_pre.mat is corrupted (neither scipy v5/7 nor valid HDF5); Pat_14 task_test.mat opens but carries no `Data` key. Both are pending vendor re-supply.
- **Same-probe bias** (CLAUDE.md invariant 4) is *partially* mitigated by ImCoh relative to MSC but not fully — any community-level interpretation should still verify with probe-debiased FC. VI at coarse k inherits probe geometry; the k-resolved tables above expose this directly (compare low-k vs high-k behaviour).

## Reproducing

```bash
conda activate lapbrain
python scripts/01_compute/compute_imcoh_vi.py --dry-run       # coverage check
python scripts/01_compute/compute_imcoh_vi.py                 # compute CSVs
python scripts/01_compute/report_h1h4_vi.py                   # regenerate this markdown
```

For figures on the same data, the multiscale scripts now default to `--fc-method imcoh_abs` (output under the corresponding `.../<fc_method>/` subfolder):

```bash
python scripts/05_multiscale/multiscale_all_hypotheses.py
python scripts/05_multiscale/continuous_multiscale_h2.py        # H2 central
python scripts/05_multiscale/definitive_multiscale_h2.py        # H2 central
python scripts/05_multiscale/continuous_all_pairs.py
python scripts/05_multiscale/analyze_h4_gradient.py
```
