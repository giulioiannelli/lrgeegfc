---
name: established-results-folder-readme
era: IMCOH_ABS_COHORT_N10
status: current
kind: index
scope: frozen verdict results with full computation provenance for preprint
---

# `.agents/preprint/established_results/` — frozen-verdict provenance

Every numerical claim that appears in the preprint manuscript must have one file in this folder. The file freezes:

1. **The claim** — the exact sentence + the band/cohort/probe it applies to.
2. **The number(s)** — to full precision; rounded forms used in the preprint are derived and explicitly noted.
3. **The full computation chain** — every script, cache, parameter, seed, and FC method needed to reproduce the number from raw time series, in order.
4. **The status** — `established` / `provisional` / `open_question` / `withdrawn`.
5. **The verification trail** — what diagnostic was run, on what date, by what means.

If a claim cannot be filed here cleanly, it cannot go in the preprint.

## Status taxonomy

- **`established`** — number reproduces deterministically from the documented chain. Verified against at least two independent reads (e.g., CSV column + hand-recomputation).
- **`provisional`** — number is consistent with the documented chain but has not been fully re-verified, or an owed control is pending. Time-stamped and tracked.
- **`open_question`** — the number's underlying methodology has an ambiguity that must be resolved before the claim can be filed as established. Files prefixed `00_open_methodology_question_*` document these.
- **`withdrawn`** — was claimed, is no longer claimed. Kept for archaeology; never deleted.

## File naming

```
<band>_<probe>_<short_description>.md
```

Examples:
- `beta_rho_split_within_baseline.md` — the +0.222 / 8/10 / p=0.005 claim
- `beta_rho_split_matched_strength.md` — the 23.7× / 7/10 / p=0.005 claim
- `beta_grassmann_window.md` — the k=27..55 / 29 contiguous cells claim
- `beta_grassmann_epi_X.md` — the 29/29 retention claim
- `00_open_methodology_question_<name>.md` — methodology questions blocking establishment of a claim

## Required frontmatter

```yaml
---
name: <slug>
era: IMCOH_ABS_COHORT_N10
status: established | provisional | open_question | withdrawn
kind: verdict | open_question
band: <band>
probe: <substrate | kc | rho_split | grassmann | anatomy | other>
claim_sentence: "<short verbatim claim as it appears or will appear in the preprint>"
number_canonical: "<exact value with units, no rounding>"
number_preprint_form: "<value as cited in the preprint, with rounding>"
computation_chain:
  - script: <relative path>
    inputs: [<cache path>, <CSV path>]
    parameters: {seed: N, nperseg: N, ...}
    output: <cache/CSV path>
verified_by:
  - {date: YYYY-MM-DD, method: <CSV row | hand-recompute | rerun | other>, by: <agent | user>}
---
```

## The discipline

- **No discrepancy gets accepted as a 'small difference'.** If two pipelines produce different numbers for the same claimed statistic, that goes here as an `open_question` until reconciled.
- **No claim graduates to `established` without an independent verification** beyond the original computing script (hand-recompute, alternative pipeline, etc.).
- **Owed controls block establishment.** A claim with an owed control sits as `provisional` until the control runs.
- **Withdrawn claims stay filed.** Never delete; mark `withdrawn` with a one-line reason and a pointer to the replacement.

## Cross-references

- `.agents/preprint/README.md` — top-level preprint folder index
- `.agents/preprint/bands/01_beta.md` — β band preprint result report (this folder's files supersede any inline numbers in that report if they ever drift)
- `.agents/guides/04_rules/never-always-list.md` — rule list including `no_hardcoded_test_thresholds`, `matched_strength_mandatory`, `brutal_honesty_no_sycophancy`
