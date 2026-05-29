---
name: preprint-folder-audit-2026-05-26
era: IMCOH_ABS_COHORT_N10
status: open
kind: response
scope: preprint-folder inconsistency audit; per-item working document, user decides each fix
companion: METHODS_AUDIT_ISSUES.md, locked/VERDICT_LEDGER.md, locked/ANATOMY_LEDGER.md
date: 2026-05-26
---

# Preprint-folder audit — itemized inconsistencies (2026-05-26)

**Head.** Audit of all 35 files in `.agents/preprint/` (~10k lines) executed
2026-05-26 in a session. Findings organised by severity; one row per
issue with the file, line refs, quoted evidence, proposed fix, and a
checkbox for status. The user decides each item individually after
verification. **Do not apply any fix without explicit user approval.**

## How to use this document

- Read each item with the file open at the cited line.
- Verify the quoted text matches the file (rg / Read).
- Decide: `[ ]` open · `[~]` partial / in progress · `[x]` done · `[-]` rejected (will not fix).
- Update the checkbox in this file once a decision is locked.
- Add a one-line note after the proposed fix if the resolution differs
  from what is proposed.
- Cross-references to `METHODS_AUDIT_ISSUES.md` are noted where the
  item duplicates an entry there (Q2/Q3/Q4/V1 etc).

The numbering is stable across iterations — `fix item 5` should still
refer to the same row even after other items are resolved.

---

## ⚠ CRITICAL — locked verdict cascade under refined Decision 12

### ▮ Decision 12 (user-locked 2026-05-26) — LOO + extent robustness as preconditions for strong tier

The 2026-05-19 pm Decision-8 mechanical rule (`strong ⇔ cluster_p_mass < 0.01`)
is refined to a conjunctive **three-condition gate**:

> **Strong-tier verdict requires all three (full data):**
> 1. `cluster_p_mass < 0.01` (locked from Decision 8 — cluster-mass-null clearance)
> 2. `cluster_p_longest_run < 0.05` (cluster-extent-null clearance; weak threshold suffices since mass is the gate)
> 3. `LOO max p_mass < 0.05` (cohort verdict robust to single-patient removal)
>
> C5 epi-X resolutions are **secondary mechanistic observations** and never promote
> a band to strong if conditions (1)–(3) fail at full data.

Rationale: `cluster_p_mass = 0.005 = 1/(R+1)` is the empirical floor for `R = 200`
and cannot discriminate marginal-clear from decisive-clear signal. The LR null and
LOO sensitivity supply the missing discrimination.

**Per-band verdict under refined rule** (from `data/audit/grassmann_cluster_extent/cohort_summary.csv` 2026-05-26):

| Band | `mass_p` < 0.01 | `LR_p` < 0.05 | LOO < 0.05 | **Verdict** |
|---|---|---|---|---|
| β | ✓ 0.005 | ✓ 0.005 | ✓ 0.005 (Pat_02) | **strong** |
| γ_l | ✓ 0.005 | ✓ 0.015 | ✓ 0.040 (Pat_05) | **strong** ↑ |
| **δ** | ✓ 0.005 | ✓ 0.025 | **✗ 0.055 (Pat_08)** | **weak** ← LOO binds |
| γ_h | ✗ 0.060 | — | — | no trace |
| θ | ✗ 0.159 | — | — | no trace |
| α | ✗ 0.348 | — | — | no trace |

**Locked verdict tags** (the four-column row used throughout the folder):

| Band | Grassmann tag |
|---|---|
| β | strong trace, both probes |
| γ_l | **strong trace, only Grassmann** ↑ (Decision-8 promotion confirmed under Decision-12 LOO+extent) |
| **δ** | **weak trace, only Grassmann** (cohort gate passes at floor; full-data LOO Pat_08-fragile; C5 epi-X strengthens, reported as secondary mechanistic observation) |
| γ_h | no trace |

**Lock action**: Decision 12 needs to be added to `locked/VERDICT_LEDGER.md`
decisions list as a new entry (after Decision 11). Separate follow-up.

---

## ⚠ CRITICAL — locked verdict cascade per Decision 12

Original cascade (pre-Decision-12, audit drafted 2026-05-26 morning) flagged
γ_l + δ "weak → strong" sweep. **Under Decision 12 the cascade splits:**

- **γ_l**: "weak" → **"strong"** everywhere it appears.
- **δ**: stays **"weak"** everywhere (the current "weak" text is correct;
  Decision 8's promotion to "strong" is now retracted).

So roughly half the items below become **no-ops** for δ, half stand for γ_l.
Item 10 (the leftover "weak δ" paragraph in `VERDICT_LEDGER.md`) **inverts** —
that paragraph turns out to be the honest residue; what needs rewriting is the
"strong δ" header + narrative at lines 240-321 of the same file.

### Item 1 — `README.md:57` — γ_l coverage tag stale  · `[x]` (2026-05-26: variant A — single-word swap `weak`→`strong`)

```
| `bands/03_gammalow.md` | γ_low (30–80 Hz) | weak trace, only Grassmann | strong localized (temporal-cortex) |
```
Fix: change `weak trace, only Grassmann` → `strong trace, only Grassmann`.
(Anatomy descriptor also needs updating — see item 13.)

### Item 2 — `README.md:60` — δ coverage tag  · `[x]` (2026-05-26: NO-OP under Decision 12 — current "weak" is correct)

```
| `bands/06_delta.md` | δ (0.53–4 Hz) | weak trace, only Grassmann | strong localized (full ≠ epi-X) |
```

Under Decision 12 (locked 2026-05-26), δ stays at "weak trace, only Grassmann"
because full-data LOO max `p_mass = 0.055 (Pat_08)` fails the robustness
precondition. The current README text is **already correct**; no edit required.
(δ anatomy descriptor "(full ≠ epi-X)" qualitatively still correct; precise
counts revisited under item 16.)

### Item 3 — `HANDOFF_INDEX.md` — γ_l + δ "weak" in 4 places  · `[x]` (2026-05-26: 3a + 3b + 3c + 3d all applied under Decision 12)

- Line 30: "γ_l (`bands/03_gammalow.md`) → … (negative reference) → `bands/05_gammah.md` → `bands/06_delta.md`. Read all six" — context line, OK.
- Line 39: "γ_l + δ carry weak Grassmann-only traces with informative anatomy."
- Line 89: "**γ_l** (`bands/03_gammalow.md`): **weak trace, only Grassmann**. Cophenet no trace (C3 p=0.116). Grassmann cluster-extent permutation cluster_p_LR=0.015 …"
- Line 92: "**δ** (`bands/06_delta.md`): **weak trace, only Grassmann**. Cophenet no trace (C3 p=0.278 …"

Fix: 3 substitutions of "weak" → "strong"; add Pat_05 LOO caveat at γ_l
(`p_mass^epi-X = 0.030, LOO max 0.159, Pat_05`) and Pat_08 LOO + C5-resolution
caveat at δ (`p_mass=0.005 full, LOO max 0.055 Pat_08; p_mass^epi-X=0.005 LOO max 0.005 fully robust`).

### Item 4 — `bands/00_cohort.md:39-40` — §1 verdict matrix stale γ_l + δ rows  · `[x]` (2026-05-26: γ_l flipped to strong; δ row no-op under Decision 12)

```
| γ_l | 30–80 | no trace | **weak** | strong localized, occipito-temporal + frontal + medial-OFC (7 named DK regions under `S(γ_l)`) | **weak trace, only Grassmann** |
| δ | 0.53–4 | no trace | **weak** | strong localized, full ≠ epi-X (4+3 regions, **0 shared** under `S(δ)`) | **weak trace, only Grassmann** |
```

Fix: column 4 (`Grassmann trace`) `**weak**` → `**strong**`; column 6
(`Coverage tag`) `**weak trace, only Grassmann**` → `**strong trace, only Grassmann** ↑`
matching VERDICT_LEDGER convention with the upward arrow.

### Item 5 — `bands/00_cohort.md:86-99` — §4 per-band paragraphs stale  · `[x]` (2026-05-26: γ_l flipped to strong; δ stays weak with LOO Pat_08 framing; γ_h Decision 6 → 6+8)

- Line 92: `### γ_l — weak trace, only Grassmann; occipito-temporal + frontal anatomy`
- Line 92 body: "Grassmann cluster-extent permutation **passes at weak** (cluster_p_LR = 0.015 …"
- Line 95: `### δ — weak trace, only Grassmann; full and epi-X read DISJOINT NETWORKS`
- Line 95 body: "Grassmann cluster-extent **passes at weak** (cluster_p_LR = 0.025 …"

Fix: flip "weak" → "strong"; change the cited gating statistic from
`cluster_p_LR = 0.015 / 0.025` to `cluster_p_mass = 0.005`; add LOO
caveats. Section headers and body text both need editing.

### Item 6 — `bands/03_gammalow.md` — whole-file stale "weak" verdict  · `[x]` (2026-05-26: full Decision-12 sweep applied; γ_l = strong; gate cite → cluster_p_mass=0.005; values post-fix all-clusters; C5 reframed as secondary observation)

The brief is anchored to the retired weak-trace framing:
- Frontmatter line 14: `verdict_tag: "weak trace, only Grassmann"`
- Frontmatter line 19: `grassmann: weak_trace (cluster-extent permutation p=0.0149, audit_70 …)`
- Headline (line 41): "Verdict from `locked/VERDICT_LEDGER.md`: **`weak trace, only Grassmann`**"
- §3.2 heading (line 123): "**weak trace at γ_l (LOAD-BEARING)**"
- §3.2 results table caption (line 154): "`cluster_p_longest_run` | **0.0149** | same — **`weak` gate**"
- §3.2 reading (line 174): "γ_l Grassmann passes the matched-strength cluster-extent gate at `weak` strength"
- §3.3 patient-by-patient: "γ_l Grassmann subspace trace per-patient (descriptive)"
- §9 line 292: "**Why γ_l Grassmann is `weak` rather than `strong`.** Three reasons …"
- Verdict line 332: "Verdict ready for writing-agent handoff: **`weak trace, only Grassmann`**"

Fix (substantial): full pass to flip every "weak" → "strong" with respect
to the Grassmann verdict; replace the gate citation from `cluster_p_LR =
0.0149` to `cluster_p_mass = 0.005` (the gate per Decision 8); update the
"Why weak" paragraph (§9) to "Why strong with LOO caveat on Pat_05 epi-X";
update §11 figure plan to cite audit_72 c5_wilcoxon CSV not the retired
audit_67 retention framing. Roughly 15-20 edits.

### Item 7 — `bands/06_delta.md` — gate citation + LOO reframing under Decision 12 (verdict stays "weak")  · `[x]` (2026-05-26: full Decision-12 sweep applied; verdict tag stays "weak"; gate cite → cluster_p_mass=0.005; LOO Pat_08=0.055 cited as the demoter; values post-fix all-clusters; C5 epi-X reframed as secondary observation strengthening trace)

- Frontmatter line 14: `verdict_tag: "weak trace, only Grassmann"`
- Frontmatter line 19: `grassmann: weak_trace (cluster-extent permutation p=0.0249, audit_70 …)`
- Frontmatter line 38: "revision_history: 2026-05-19: … δ Grassmann promoted no trace → **weak trace** via cluster-extent permutation"
- Head line 45: "Verdict from `locked/VERDICT_LEDGER.md`: **`weak trace, only Grassmann`** (promoted from no trace under the previous 8-cell hardcoded threshold via cluster-extent permutation, Decision 6)"
- §3.2 heading (line 123): "**weak trace (cluster-extent promotion)**"
- §3.2 reading (line 165): "δ Grassmann passes cluster-extent at `weak` strength (cluster_p = 0.0249)"
- §8 row line 298: "C3 matched-strength … ✓ cluster_p = 0.0249 (**promoted from no trace** under cluster-extent revision)"
- §9 line 318: "**Why δ is `weak trace` rather than `strong trace` or `no trace`.** Three reasons …"
- Verdict line 358: "Verdict ready for writing-agent handoff: **`weak trace, only Grassmann`**"

Fix (substantial): same pattern as γ_l. Replace gate citation
`cluster_p_LR = 0.0249` → `cluster_p_mass = 0.005`. Add the Pat_08 LOO
flag at full data + C5 epi-X resolution (`p_mass^epi-X = 0.005`, mass
38 → 44, LOO max under epi-X = 0.005 fully robust). §11 figure plan
needs the C5 Wilcoxon-on-epi-X framing (Decision 10), not the
retention/window-shift framing.

### Item 8 — `locked/ANATOMY_LEDGER.md:25-26` — top table γ_l/δ trace verdicts stale  · `[x]` (2026-05-26: γ_l flipped to strong ↑; δ stays weak with Decision-12 LOO Pat_08 caveat inline)

```
| γ_l | weak trace, only Grassmann | Grassmann | **strong localized, only Grassmann** (6 named DK regions, temporal-cortex-dominant) |
| δ | weak trace, only Grassmann | Grassmann (full + C5 epi-X) | **strong localized, only Grassmann; full and epi-X read DIFFERENT NETWORKS** (6 + 6 named regions, only 2 shared) |
```

Fix: column 2 (`Trace verdict`) flip both to `strong trace, only Grassmann`.
This also fixes the trace-side parity within the locked artifact.
Anatomy column also needs updating — see items 13-14.

### Item 9 — `locked/ANATOMY_CONTROLS.md:140-141` — object-to-control mapping stale  · `[x]` (2026-05-26: γ_l flipped to strong ↑; δ stays weak with Decision-12 LOO Pat_08 caveat inline)

```
| γ_l | weak trace, only Grassmann | Grassmann (top participation over `S(γ_l)`) | yes |
| δ | weak trace, only Grassmann | Grassmann (top participation over `S(δ)` full + `S(δ)` epi-X) | yes for both windows |
```

Fix: column 2 flip both to `strong trace, only Grassmann`.

---

## ⚠ CRITICAL — internal self-contradiction in `VERDICT_LEDGER.md`

### Item 10 — `VERDICT_LEDGER.md` δ section + locked verdict table + decisions list  · `[x]` (2026-05-28: full rewrite under Decision 12 — δ verdict flipped strong → weak; orphan residue paragraph removed; locked verdict table row 35 updated; Decision 12 added to decisions list)

After the §"δ" table (line 320 end), this paragraph appears mid-flow
(begins with a sentence fragment):

```
real subspace signature. Cophenet C3 fails so no per-pair trace. The δ C4
anchor-anatomy reading remains a *separate, descriptive* known-biology
observation, not part of the trace verdict. → **weak trace, only Grassmann**
(Grassmann), with descriptive anchor-anatomy note alongside.
```

This concludes with `weak trace, only Grassmann` — directly contradicting
the same file's headline locked table (line 35) and the §"δ (0.53–4 Hz)
— **strong trace, only Grassmann**" header (line 240). The same locked
artifact says both weak and strong for δ.

Fix: delete these 5 lines. The "Verdict reasoning" block at lines 315-321
already states the correct verdict.

### Item 11 — Decision 6/7 ordering inversion  · `[x]` (2026-05-28: applied (B) head note; clarifies file-order vs chronological + flags 8/9/10 living in CONTROLS.md and 11 in revision history)

`VERDICT_LEDGER.md` lists decisions in this file order: 1, 2, 3, 4, 5,
**7** (line 454), **6** (line 479), 8, 9, 10, 11. Not chronological,
not strictly dependency-ordered (Decision 7 disjunctive gate is retired
by Decision 8 mass-only).

Fix (optional): either renumber to file order, or add a head note "Decisions
ordered by dependency, not chronology". Low priority — accuracy is intact,
only readability suffers.

---

## ⚠ HIGH — `T_G*` cluster-extent numbers stale + C1 normalization not propagated

The CSV `data/audit/grassmann_cluster_extent/cohort_summary.csv` is
post-2026-05-19-pm all-clusters fix (β 69.76, γ_l 66.14, δ 38.07, γ_h
31.67, θ 12.97, α 7.98 with `cluster_p_mass` 0.005 for β/γ_l/δ).

### ▮ Lock note 0a — C1 normalization to `[0,1]` (referenced by all items in this section)

Per `METHODS_AUDIT_ISSUES.md:116-129` (Locked decision C1) + user decisions
2026-05-26 on cohort-vs-per-patient `n_k`:

> **Cohort scalar:**
> $T_G^*(b) = \frac{1}{n_k^{\text{cohort}} \log_{10}(R+1)} \sum_{k:p_k(b)<\alpha_k} (-\log_{10} p_k(b)) \in [0,1]$
>
> **Per-patient scalar:**
> $T_G^{*,s}(b) = \frac{1}{n_k^{s} \log_{10}(R+1)} \sum_{k:p_k^s(b)<\alpha_k} (-\log_{10} \max(p_k^s(b), 1/(R+1))) \in [0,1]$

**`n_k` convention** (user-locked 2026-05-26):
- **Cohort**: `n_k^cohort` = **biggest k-grid common to all patients in the cohort
  condition**. For full-data n=10 cohort this is `n_k^cohort = 111` (k ∈ [2, 112],
  set by Pat_10 with `N = 113` contacts being the smallest).
- **Per-patient**: `n_k^s` = patient s's own available k-grid (= `N_s − 1` or the
  effective epi-X-excluded range under C5). Rationale: per-patient `T_G^{*,s} = 1`
  then means "this patient saturates their own significance budget across the
  full available subspace", which is comparable across patients with different `N`.
- `R = 200`, log10(R+1) = log10(201) ≈ 2.3032.
- **Full-data cohort denominator**: `111 × 2.3032 ≈ 255.65`.
- **C5 epi-X cohort denominator**: computed from
  `data/audit/grassmann_epi_exclusion/epi_counts.csv` 2026-05-26.
  `min(N_reduced) = 89` (Pat_13, 30 epi contacts), so cohort-common k-grid
  is `k ∈ [2, 88]`, giving **`n_k^epi-X-cohort = 87`** and denominator
  `87 × 2.3032 ≈ 200.38`. (Caveat: verify audit script
  `scripts/01_compute/audit/audit_72_*.py` actually aggregates the cohort
  over this common grid; if it uses per-patient grids and then aggregates,
  the cohort raw mass values 89.04 / 32.75 / 43.99 require a different
  normalization convention. Flagged for compute-side verification.)

**Locked normalized cohort values** (full-data, `n_k^cohort = 111`):

| Band | Raw `T_G*` | **Normalized `T_G*` (locked citation)** | `p_mass` |
|---|---|---|---|
| β | 69.76 | **0.273** | 0.005 |
| γ_l | 66.14 | **0.259** | 0.005 |
| δ | 38.07 | **0.149** | 0.005 |
| γ_h | 31.67 | **0.124** | 0.060 |
| θ | 12.97 | **0.0507** | 0.159 |
| α | 7.98 | **0.0312** | 0.348 |

**C5 epi-X cohort normalized values** (using `n_k^epi-X-cohort = 87`, denominator `200.38`):

| Band | Raw `T_G*^epi-X` | **Normalized `T_G*^epi-X`** | `p_mass^epi-X` | vs full |
|---|---|---|---|---|
| β | 89.04 | **0.444** | 0.005 | strengthens (full 0.273 → epi-X 0.444) |
| γ_l | 32.75 | **0.163** | 0.030 | contracts (full 0.259 → epi-X 0.163) |
| δ | 43.99 | **0.220** | 0.005 | strengthens (full 0.149 → epi-X 0.220) |

(My earlier tentative computation 0.348 / 0.128 / 0.172 used `n_k = 111`,
which is the wrong denominator for the C5 condition. Corrected to 0.444 /
0.163 / 0.220 above using user-locked "biggest common across the relevant
cohort condition" rule with Pat_13 as the limiting patient at 89 non-epi
contacts. Verification of cohort-common-grid convention in audit script
still flagged above.)

**Citation policy decisions (user-locked 2026-05-26):**
- **Q1 — locked ledgers (`VERDICT_LEDGER.md`)**: **dual format** — raw value
  with the normalized value in parentheses immediately after. Example:
  `mass = 69.76 (T_G* = 0.273)`. Be explicit and consistent (not "same" or
  ambiguous shorthand).
- **Q2 — historical passages** (e.g., `methods_section_review_2026-05-19.md:54-55`
  documenting the pre-fix vs post-fix formula change): **leave raw** with a
  one-line forward pointer to the C1 normalization.
- **Q3 — `n_k` convention**: cohort uses biggest common `n_k^cohort`; per-patient
  uses own `n_k^s` (encoded above).

**Implications for items below:**
- Items 12, 13, 14, 14b, 14c, 14d, 14e all become **dual fix** per Q1 in locked
  ledgers; normalized-only in briefs + manuscript-facing files.
- C5 epi-X normalized values cannot be locked into items 14b et al. until
  `n_k^epi-X-cohort` is determined.

### Item 12 — `bands/00_cohort.md:67-82` — §3 cluster-extent table stale  · `[x]` (2026-05-28: full §3 table rewrite — post-fix all-clusters values, C1-normalized T_G* column added, Decision-8 mass-only gate + Decision-12 LOO precondition cited, verdicts updated)

This is `METHODS_AUDIT_ISSUES.md` Q2 (confirmed not yet done).

Stale row values (`obs_LR`, `obs_mass`, `cluster_p_mass`):
- β: 29 / **52.97** / 0.005 → should be 29 / **69.76** / 0.005
- δ: 7 / **12.78** / **0.025** → 7 / **38.07** / **0.005**
- γ_l: 13 / **19.17** / **0.035** → 13 / **66.14** / **0.005**
- γ_h: 9 / 17 / **15.75** / **0.065** → 9 / 19 / **31.67** / **0.060**
- θ: 5 / 17 / **7.11** / **0.144** → 5 / 17 / **12.97** / **0.159**
- α: 4 / 12 / 7.98 / **0.099** → 4 / 12 / 7.98 / **0.348**

Header at line 69 also stale: "Disjunctive gate as of 2026-05-19: pass
if EITHER cluster_p_LR<0.05 OR cluster_p_mass<0.05." Decision 8 retired
this (mass-only). Replace.

Verdict column at lines 73-78 also needs flipping (γ_l/δ to "strong")
to match items 4-5.

### Item 13 — `methods/methods_grassmann_cluster_extent.md` Head + §5b + §6 + §9 + §10 + §11  · `[x]` (2026-05-28: full Decision-8/12 + C1 normalization sweep — head + §5b + §6 + §9 checklist + §10 refs + §11 revision log all updated; per-patient T_G*,s formula with n_k^s added; gate citations switched disjunctive → mass-only + LOO precondition; γ_h removed from "contributes" band list)

Multiple stale citations in the same file:
- Line 26 (Head): "β at `T_G^* = 52.97`, `cluster_p_mass = 0.005` is the **load-bearing band**" → β = **69.76**, vocabulary swap (see Low §)
- Line 230 (§5b): "CONTROLS.md C3 records the rule as the **disjunctive gate** `min(cluster_p_LR, cluster_p_mass) < α`" → mass-only per Decision 8 (consistent with §5d in same file, lines 267-281)
- Line 596 (§9 checklist): "give the β number first **(52.97**, `cluster_p_mass = 0.005`)" → 69.76
- Line 605 (§9 step 9): "Do not cite the LR-only verdict column from the CSV (`verdict_cluster_extent`); cite the **disjunctive verdict** from CONTROLS.md C3 (computed from `min(cluster_p_LR, cluster_p_mass)`)" → cite mass-only verdict
- Line 626 (§10): "Locked gate: `.agents/preprint/locked/CONTROLS.md` §C3 Grassmann + Decision 7 disjunctive verdict" → Decision 8 mass-only
- Line 629 (§10): "Per-band briefs (where Grassmann is load-bearing or contributes): `bands/01_beta.md`, `bands/03_gammalow.md`, `bands/05_gammah.md`, `bands/06_delta.md`" — `05_gammah.md` is `no trace` so doesn't "contribute" (false advertising), and the wording uses "load-bearing"

§6 table (lines 420-430) is current — uses post-fix values + mass-only
gate. So this file is mid-revision: §6 was updated but §1/§5b/§9/§10
were not.

Fix: 5-6 edits to align head + checklist with §5d + §6.

### Item 14 — `directives/methods_directive_2026-05-19_TG_normalization.md` pre-fix table  · `[x]` (2026-05-28: option (B) — frontmatter `status: applied` + applied_date + applied_in trace + notes_on_application with post-fix values; body kept verbatim as historical record)

- §"Numerical anchors" (line 56-63): "Raw: `T_G^*(β) = 52.97`. Normalized: `T_G^*(β) = 52.97 / 255.65 = 0.2073` ≈ **20.7%**"
- §"Open item" (line 64): "the standing discrepancy between the cohort_summary CSV (sums over the longest contiguous cluster only) and the methods companion equation (sums over all sig cells). Methods agent should pick one and lock it" — RESOLVED 2026-05-19 pm (CSV is now all-clusters)
- Per-band table (lines 70-78): all pre-fix raw values (β 52.97 / γ_l 19.17 / δ 12.78 / α 7.98 / γ_h 15.75 / θ 7.11) and pre-fix `p_mass` (γ_l 0.035, δ 0.025, γ_h 0.065, θ 0.144, α 0.099) — should use post-fix values

Fix: replace raw values; annotate the "Open item" as resolved; recompute
normalized percentages (β raw 69.76 / 255.65 ≈ 27.3%, γ_l 66.14 / 255.65
≈ 25.9%, etc.).

Note: this directive may be functionally superseded — the methods agent
has presumably already applied the normalization. If so, this file should
be marked `status: applied` or moved to a directives/archive. Open
question for the user.

### Item 14b — `locked/VERDICT_LEDGER.md` 6 lines cite raw `T_G*` (post-C1 review)  · `[x]` (2026-05-28: dual format `raw (normalized)` applied to β/γ_l/δ full-data + C5 epi-X tables + revision narratives — 8 line edits)

Found during the 2026-05-26 verification grep. Raw values appear at:
- Line 101 β table: "**Resilient all-clusters mass `T_G*`** | **69.76** | same"
- Line 105 β C5: "C5 epi-X `T_G*^epi-X` | **89.04** (strengthens vs full 69.76)"
- Line 109 β narrative: "β trace strengthens under epi-X (mass 69.76 → 89.04, LR 29 → 36)"
- Line 206 γ_l table: "**Resilient all-clusters mass `T_G*`** | **66.14** | same"
- Line 210 γ_l C5: "C5 epi-X `T_G*^epi-X` | 32.75 (vs full 66.14, contracts)"
- Line 233 γ_l revision narrative: "corrected formula gives `T_G* = 66.14`"
- Line 247 δ revision narrative: "corrected resilient formula gives `T_G* = 38.07`"
- Line 296 δ table: "**Resilient all-clusters mass `T_G*`** | **38.07** | same"
- Line 300 δ C5: "C5 epi-X `T_G*^epi-X` | **43.99** (strengthens vs full 38.07)"

Fix options:
- **(a)** Replace with normalized values (β 0.273 / 0.349 epi-X; γ_l 0.259
  / 0.128 epi-X; δ 0.149 / 0.172 epi-X) — matches C1 brief/manuscript citation.
- **(b)** Keep raw values in locked ledger for forensic clarity; add a
  parenthetical normalized value next to each (e.g., "mass = 69.76 (T_G* = 0.273)").
- **(c)** Migrate the in-table column header from "raw mass" to "T_G* (normalized)"
  with values replaced.

C5 epi-X normalized values (computed): β 89.04/255.65 = 0.348; γ_l 32.75/255.65 = 0.128;
δ 43.99/255.65 = 0.172. (User decision needed on whether to also recompute via
the precise `n_k^epi-X · log10(R+1)` if n_k differs under epi-X exclusion.)

### Item 14c — `bands/02_alpha.md:186` raw α `T_G* = 7.98` cited  · `[x]` (2026-05-28: dual raw + normalized 0.0312; also updated null mean/p95 to post-fix all-clusters values 7.35/21.26; updated cluster_p_mass = 0.348 + Decision-8 gate label)

α brief Grassmann results table cites raw `7.98` not normalized `0.0312`.
Since α is a no-trace verdict it's borderline whether normalization matters
here — but consistency with the panel pushes for normalized citation.

Fix: line 186 swap raw `7.98` → normalized `0.0312` (or dual: `0.0312 (raw 7.98)`).

### Item 14d — `methods/methods_neurophysiological_interpretation_2026-05-26.md` 7 raw citations  · `[x]` (2026-05-28: cohort summary table + §4.3 γ_l + §4.4 δ reframed under Decision 12; T_G* normalized + raw dual format; γ_l strong; δ weak with LOO Pat_08 binding; C5 explicitly secondary)

This methods file was authored 2026-05-26 (this same audit session) but used
raw T_G* values. Stale per C1 from creation.

Lines with raw values:
- Line 73: "β mass=69.76 / p=0.005 / LOOmax=0.005 — decisive, both C5 gates pass"
- Line 74: "γ_l mass=66.14 / p=0.005 / LOOmax=0.040"
- Line 75: "δ mass=38.07 / p=0.005 / LOOmax=0.055"
- Line 76: "α mass=7.98 / p=0.348"
- Line 77: "γ_h mass=31.67 / p=0.060"
- Line 161: "Trace at Grassmann (mass=66.14, p=0.005 …"
- Line 167: "Grassmann strong (mass=38.07, p=0.005) …"

Fix: 7 line-edits swapping `mass=<raw>` → `T_G*=<normalized>` per the locked table.

### Item 14e — `methods/methods_section_review_2026-05-19.md:54-55` raw values  · `[x]` (2026-05-28: per Q2 — raw values kept verbatim as historical record of pre-fix→post-fix formula change; forward pointer added documenting C1 normalization + Decision-12 δ demotion retraction)

Body narrative: "in the manuscript draft (β 52.97, γ_l 19.17, δ 12.78) are
obsolete; re-run values are β 69.76, γ_l 66.14, δ 38.07."

This was correct *at the time* (pre-C1) — narrative is documenting the
formula-revision history. Two readings:
- **Historical record** — leave raw, add a footnote "per C1 lock 2026-05-20
  these are subsequently normalized to β 0.273 / γ_l 0.259 / δ 0.149".
- **Live citation** — replace with normalized.

Since this file is `status: open_for_writing_agent` and items M1 etc. are
RESOLVED, recommend marking it `mostly_superseded` (per item 22) and
treating this passage as historical — leave raw with a one-line forward
pointer.

---

## ⚠ HIGH — `ANATOMY_LEDGER.md` top table vs revision history disagree

### Item 15 — γ_l anatomy: top table 6 regions temporal-cortex; revision history 7 regions occipito-temporal+frontal+OFC  · `[x]` (2026-05-28: top-table descriptor updated + §3 γ_l subsection replaced with current cluster-extent 7-region anatomy + KC-era left fusiform retraction noted)

Top table line 25: "(6 named DK regions, temporal-cortex-dominant)"

Revision history lines 203-204:
> γ_l Grassmann: |S(γ_l)| = 41 cells, was K*=[12,23] 13 contiguous. A3-passing region set **shifted to 7 regions**, only 3/6 overlap with the locked ledger. Retained: ctx-lh-middletemporal, ctx-lh-superiortemporal, ctx-rh-parstriangularis. **Dropped**: ctx-lh-inferiortemporal, ctx-lh-fusiform, ctx-rh-paracentral. **Added**: ctx-lh-lateraloccipital, ctx-lh-rostralmiddlefrontal, ctx-rh-medialorbitofrontal, ctx-lh-cuneus. The "temporal-cortex dominant" narrative weakens — the new network is **occipito-temporal + frontal + medial-OFC**.

The §3 "Per-(band, probe) verdict detail" subsection at lines 87-101
still uses the retired K*(γ_l)=[12,23] 6-region anatomy.

Fix: either (a) replace the top-table cell and §3 detail subsection with
the cluster-extent 7-region anatomy (the right move under "locked file =
single source of truth"); (b) demote the §3 details to a "Retired K*
anatomy (archived)" sub-section and elevate the revision-history numbers
to the top. Option (a) recommended.

### Item 16 — δ anatomy: top table 6+6 regions 2 shared; revision history 4+3 regions 0 shared  · `[x]` (2026-05-28: top-table descriptor → 4+3 / 0 shared; §3 δ Grassmann subsection replaced with S(δ) + S^epiX(δ) regions; cross-band table + Pattern bullets regenerated under cluster-extent paradigm; KC-era retractions documented for fusiform, Amy, bankssts, caudal ACC, paracentral)

Top table line 26: "(6 + 6 named regions, only 2 shared)"

Revision history lines 205-206:
> δ-full Grassmann: |S(δ)| = 23 cells, was K*=[57,63] 7 contiguous. A3-passing region set **shrinks to 4 regions** … The "anchor anatomy" interpretation (Amy + cingulate + fusiform) is not supported under S(δ) — replaced by a parietal-temporal-frontal network.
> δ-epi-X Grassmann: |S^epiX(δ)| = 22 cells … A3-passing region set **shrinks to 3 regions** … δ full-vs-epiX dissociation strengthens under S(b): the two networks now share **0/3 named regions** (locked ledger had 2/6 shared via fusiform + inferiorparietal). The "mixture of two distinct phenomena" claim becomes stronger.

§3 "Per-(band, probe) verdict detail" subsection at lines 103-131 still
narrates K*(δ)=[57,63] 6-region full anatomy + K*^epiX=[33,39] 6-region
epi-X anatomy.

Fix: replace top-table cell ("4 + 3 named regions, **0 shared**"); update
§3 detail subsection or demote to archive.

§"Cross-band anatomy comparison" table at lines 137-153 also uses retired
K*(b) data — has "left fusiform" at γ_l Grassmann + δ Grassmann (both
full and epi-X) + "Amy" at δ Grassmann full. Under cluster-extent paradigm,
left fusiform appears nowhere and Amy retracts. This table must be
regenerated under S(b).

### Item 17 — `README.md:57` γ_l anatomy descriptor "temporal-cortex" stale  · `[x]` (2026-05-28: one-line swap to "occipito-temporal + frontal + medial-OFC")

```
| `bands/03_gammalow.md` | γ_low (30–80 Hz) | weak trace, only Grassmann | strong localized (temporal-cortex) |
```

Fix: descriptor → "occipito-temporal + frontal + medial-OFC". (Same row
also has the verdict-flip item 1.)

### Item 18 — `HANDOFF_INDEX.md:87` β cophenet anatomy missing 7th region  · `[x]` (2026-05-28: added missing rostralanteriorcingulate; isthmus cingulate now spelled out; vocab swap "Load-bearing" → "Primary finding" in both files; explicit "7 named DK regions" count added to both cophenet + Grassmann descriptors)

Line 87: "Anatomy: 7+7 named DK regions across the two probes (cingulate
+ parahippocampal + entorhinal + insula + postcentral + superior frontal
on cophenet; Hippocampus + temporal + orbitofrontal + insula + rostral
middle frontal on Grassmann)."

Cophenet list has 6 region names; the locked ANATOMY_LEDGER (lines
36-43) lists 7: isthmuscingulate + superiorfrontal + insula +
parahippocampal + postcentral + entorhinal + **rostralanteriorcingulate**.
Missing: `rostralanteriorcingulate`.

`bands/00_cohort.md:87` β paragraph has the same 6-region cophenet list
— also missing `rostralanteriorcingulate`.

Fix: add `rostralanteriorcingulate` to both lists.

---

## Medium — folder index missing files / stale pointers

### Item 19 — `README.md` doesn't list new methods files  · `[x]` (2026-05-28: added pointers to METHODS_AUDIT_ISSUES.md + methods_neurophysiological_interpretation + methods_section_review with mostly_superseded tag)

- `methods/methods_neurophysiological_interpretation_2026-05-26.md` —
  created 2026-05-26, not in README
- `METHODS_AUDIT_ISSUES.md` — created 2026-05-20, not in README

Fix: add one-line pointers under §"Companion guides" or §"Methods" section.

### Item 20 — `WRITING_GUIDE.md:48-51` folder-layout block missing methods file  · `[x]` (2026-05-28: added neurophys interpretation file; tagged methods_section_review as mostly_superseded)

Lists 3 methods/ files; folder now has 4 (added neurophysiological_interpretation
2026-05-26).

Fix: add the file to the listing.

### Item 21 — `tables/beta_per_patient.md:9` stale pre-reorg preflight path  · `[x]` (2026-05-28: path corrected to `responses/2026-05-18_beta_per_patient_table_preflight.md`)

```
preflight: .agents/preprint/2026-05-18_beta_per_patient_table_preflight.md
```

Actual location after 2026-05-19 reorg:
`.agents/preprint/responses/2026-05-18_beta_per_patient_table_preflight.md`

Fix: update the path.

### Item 22 — `methods/methods_section_review_2026-05-19.md` status field stale  · `[x]` (2026-05-28: status → mostly_superseded + status_updated + superseded_by METHODS_AUDIT_ISSUES.md; head note added with forward pointer)

Frontmatter: `status: open_for_writing_agent`. But M1/M4/m6/m9 are all
marked RESOLVED/SUPERSEDED in the body, and `METHODS_AUDIT_ISSUES.md`
frontmatter says it "supersedes: methods/methods_section_review_2026-05-19.md
(partially — see resolution map §1)".

Fix: change status to `mostly_superseded` (or similar) + add head note
pointing at `METHODS_AUDIT_ISSUES.md` for live items.

### Item 23 — `methods/methods_revision_2026-05-18_cophenet.md` checklist contradicts own body  · `[x]` (2026-05-28: checklist item 469-472 rewritten to match 2026-05-20 broadening — drop all cross-band BH-FDR at LRG probes per `feedback_no_unmotivated_bh_fdr.md`; keep BH-FDR only for anatomy A1 hypergeometric)

Body lines 281-294 (locked 2026-05-18 + 2026-05-20 supersedure): "**no
cross-band BH-FDR is applied at any LRG probe**".

Checklist line 469-471: "Statistical inference subsection updates family
sizes to `m = 12` (matrix distances on `D_coph`, six bands × two distances)
and `m = 6` (`ρ_split^coph`, six bands). Drops KC, VI, taxonomy family-size
references."

Same file says both "no cross-band BH-FDR" and "BH-FDR at m=6 / m=12".

Fix: rewrite the checklist item to match the 2026-05-20 broadening
("Statistical inference subsection: drop all cross-band BH-FDR statements
across LRG probes; keep BH only for the anatomy A1 hypergeometric across
DK regions per (band, probe)").

### Item 24 — `directives/methods_directive_2026-05-19_TG_normalization.md` filename uses non-canonical prefix  · `[x]` (2026-05-28: not renamed to avoid breaking inbound references; `filename_note:` added to frontmatter explaining that `methods_directive_` prefix is intentional for methods-side directives — distinct from writing_directive_ writing-agent tasks)

File is `methods_directive_2026-05-19_TG_normalization.md`. `WRITING_GUIDE.md:96`
specifies directive filename: `writing_directive_YYYY-MM-DD_<short-topic>.md`.

Fix: rename to `writing_directive_2026-05-19_TG_normalization.md` and
update any references in other files, OR amend WRITING_GUIDE to allow
`methods_directive_` for methods-side directives.

Note: this file may also be functionally superseded (see item 14 — methods
agent has presumably normalized). If superseded, archive rather than
rename.

### Item 25 — `established_results/` folder is essentially stub  · `[x]` (2026-05-28: option (a) — README rewritten as "historical methodology Q&A" + aspirational example list dropped; documents that bands/ briefs + locked/ ledgers absorbed the frozen-claim role in practice; only `00_open_methodology_question_lrg_D_convention.md` (withdrawn) remains as documented historical artifact)

Folder contents:
- `README.md` (describes a frozen-claims-with-provenance system)
- `00_open_methodology_question_lrg_D_convention.md` (`status: withdrawn`)

The README §"Examples" lists `beta_rho_split_within_baseline.md`,
`beta_rho_split_matched_strength.md`, `beta_grassmann_window.md`,
`beta_grassmann_epi_X.md` — none of these files exist; they are
aspirational. The actual frozen-claims work happens in `bands/`.

Fix options:
- (a) Demote folder to "Historical methodology Q&A" + delete the
  aspirational example list from README, leaving only the withdrawn file
- (b) Populate properly (one frozen claim per CSV cited in manuscript) —
  large effort, unclear gain
- (c) Fold the withdrawn file into a `methods/` or `responses/` archive,
  remove the folder

User decision needed on direction.

---

## Medium — claims without quantitative support / outdated control language

### Item 26 — `bands/00_cohort.md:99` γ_h Decision-citation incomplete  · `[x]` (2026-05-28: already applied in item 5 — verified — "Decision 6 demotion" → "Decisions 6 + 8 demotion")

"Demoted to 'no trace' in the 2026-05-19 cluster-extent revision (**Decision 6**)"

Decision 6 (8-cell threshold retired in favor of cluster-extent
permutation). γ_h verdict also satisfies Decision 8 (mass-only,
`p_mass = 0.060`). Both decisions support the no-trace verdict.

Fix: cite "Decision 6 + 8" or "the cluster-extent revision" (band-neutral),
not Decision 6 alone.

### Item 27 — `bands/00_cohort.md:142` §6 cites retired disjunctive gate  · `[x]` (2026-05-28: disjunctive LR ∨ mass → mass-only Decision 8 + Decision-12 LOO precondition)

"Cluster-extent permutation gates Grassmann at p<0.05 on the **disjunctive
LR ∨ mass statistic**."

Decision 8 retired disjunctive in favor of mass-only.

Fix: → "on the cluster-mass statistic (mass-only gate per Decision 8)".

### Item 28 — `bands/03_gammalow.md` §3.2 + §8 sensitivity panel uses retention framing  · `[x]` (2026-05-28: already applied in item 6 — verified — retention rule retired, C5 reframed as secondary observation per Decision 10)

Lines 156-169 narrate C5 as "10-cell longest run at k=19..28 vs 13 cells
at k=12..23" + "Cell-count retention: 10/13 ≈ 77% (below the 80% guide).
Within-window overlap: 5/13 ≈ 42%".

Decision 10 retired the ≳80% retention rule. C5 gate is now
`cluster_p_mass^epi-X < 0.05` (Wilcoxon-based).

Fix: reframe the §3.2 robustness paragraph and the §8 sensitivity panel
to cite the Wilcoxon-on-epi-X gate (`cluster_p_mass^epi-X = 0.030 < 0.05`,
LOO max p_mass^epi-X = 0.159 Pat_05). Retention numbers should be
demoted to descriptive co-statistics. Cross-reference `c5_wilcoxon_cohort.csv`.

### Item 29 — `bands/06_delta.md` §3.2 + §8 + §11 — same retention-framing issue  · `[x]` (2026-05-28: already applied in item 7 — verified — retention rule retired, C5 reframed as secondary mechanistic observation per Decision 10; §11 figure plan updated audit_67 → audit_72 c5_wilcoxon)

Lines 152-163 narrate C5 as "Within-window overlap: 0% — the full-cohort
window k=57..63 and the epi-X window k=33..39 do not intersect at any k.
This is the most extreme window shift in the panel".

§11 figure plan (line 354): "F4 (audit_67 epi-X overlay showing the
window shift k=57..63 → k=33..39 with 0% overlap)"

Decision 10 C5 gate is Wilcoxon-based. Under that gate δ Grassmann
strengthens (mass 38 → 44, `p_mass^epi-X = 0.005`, LOO max = 0.005 fully
robust — `methods/methods_grassmann_cluster_extent.md` §6 LOO table). The
"0% overlap" is descriptive and should be reframed accordingly.

Fix: reframe §3.2/§8/§11 to cite Decision-10 gate + the `p_mass^epi-X =
0.005` strengthening + the LOO resolution of Pat_08 full-data leverage.
Retain "0% overlap" as descriptive but not as the gating statistic.

### Item 30 — `methods/methods_grassmann_cluster_extent.md` §7 C5 retention rule retired  · `[x]` (2026-05-28: §7 rewritten — 80% retention rule retired, locked Decision-10 Wilcoxon-on-epi-X gate cited with `cluster_p_mass^epi-X < 0.05`, audit_67 → audit_72 c5_wilcoxon path; explicit "secondary mechanistic observation, never verdict-promoter" framing)

Lines 531-544 describe C5 as: "≳ 80% retention (cell count or `T_G^*`
retention) supports a strong verdict; < 80% retention or a substantial
`k`-window shift downgrades to weak; Window emergence under epi-X is
interpreted as the epi zone masking a wider trace."

Decision 10 retired this in favor of Wilcoxon-on-epi-X.

Fix: replace §7 content with the locked C5 Wilcoxon gate (from
CONTROLS.md §C5) — `cluster_p_mass^epi-X < 0.05` for Grassmann; one-sample
Wilcoxon `H_1: rho^epi-X > 0` for cophenet α.

### Item 31 — `bands/04_theta.md:245` stub figure reference  · `[x]` (2026-05-28: F_cohort_2 citation annotated as "not yet produced, deferred"; θ panel scoped as standalone descriptive supplement rather than part of an unrealised cross-band suite)

Line 245 (step 1 of §11): "F_cohort_2 cross-band verdict matrix"

`bands/00_cohort.md:157` says: "F_cohort_1_three_layer.pdf,
F_cohort_2_verdict_matrix.pdf, F_cohort_3_grassmann_strip.pdf — see Phase
C figures sub-task; **not yet produced (figure scripts deferred until ready)**"

Stub — no scripts exist, no figures generated.

Fix: either commit to producing F_cohort_*, or remove the citations. If
deferred, status note ("figure deferred to Phase C") should accompany
each citation.

### Item 32 — `bands/01_beta.md` audit-path wildcards  · `[x]` (2026-05-28: `audit_*_ctm_triangle.py` resolved to `audit_33_ctm_triangle.py`; `audit_*_lrg_localization_anatomy + implant_geometry` resolved to the named retired KC-era scripts; data/audit/lrg_phase_distance paths verified present)

Multiple lines reference `audit_*_ctm_triangle.py` (wildcard, not a
specific filename):
- Line 158 cache provenance
- Line 268 within-baseline script
- Line 357 audit data row

Other paths that should be ls-checked:
- `data/audit/lrg_phase_distance/cohort_summary.csv` (lines 130, 152, 153, 155)
- `data/audit/lrg_phase_distance/drift_R2.csv` (line 132)
- `scripts/02_preprint/preprint_03_beta_matched_strength_raw_D.py` (line 273)
- `scripts/02_preprint/preprint_05_allbands_matched_strength_raw_D.py` (line 274)

Fix: resolve the wildcard to actual script names; `ls` each cited path
once and either correct or remove if not present.

Commands to verify:
```bash
ls scripts/01_compute/audit/audit_*_ctm_triangle*.py
ls data/audit/lrg_phase_distance/
ls scripts/02_preprint/preprint_0[35]_*.py
```

---

## Low — vocabulary violations per `feedback_avoid_load_bearing_facet.md` (2026-05-26)

The rule was saved this session; pre-existing files predate it.
Cleanup-grade, not blocker-grade.

### Item 33 — "load-bearing" usage sweep  · `[x]` (2026-05-28: all live preprint files swept — bands/00-06, HANDOFF_INDEX, WRITING_GUIDE, EVALUATION_PROTOCOL, METHODS_AUDIT_ISSUES, methods_grassmann_cluster_extent, methods_revision_cophenet, tables/beta_per_patient. Substitutions: primary / central / principal / decisive. Skipped: directives/, responses/, methods_section_review_2026-05-19 (status: mostly_superseded) — kept verbatim as historical record)

Approximate counts (grep-verifiable):
- `bands/00_cohort.md` (1×)
- `bands/01_beta.md` (~6×)
- `bands/02_alpha.md` (~5×, several **LOAD-BEARING** uppercase in headings)
- `bands/03_gammalow.md` (~3×, including a §3.2 heading)
- `HANDOFF_INDEX.md` (~3×)
- `methods/methods_revision_2026-05-18_cophenet.md` (~4×)
- `methods/methods_grassmann_cluster_extent.md` (~6×, including an emphasized `**load-bearing**` at line 60)
- `methods/methods_section_review_2026-05-19.md` (~3×)
- `tables/beta_per_patient.md` (1× implicit)
- `directives/methods_directive_2026-05-19_TG_normalization.md` (1×, in `priority` frontmatter line)

Substitution menu (from feedback memory): primary, central, headline,
principal, decisive, definitive, "the result that carries the X".

Fix: opt-in sweep when other items in the same file are being touched.
Not a separate pass.

### Item 34 — "facet"/"facets" usage sweep  · `[x]` (2026-05-28: all live preprint files swept — bands/01_beta + bands/02_alpha + methods_revision_cophenet + tables/beta_per_patient. Substitutions: aspect / aspects. Skipped: directives/ (historical task records) — kept verbatim)

- `bands/01_beta.md` §3.3 ("two probes read **distinct facets**")
- `methods/methods_revision_2026-05-18_cophenet.md` ("structurally distinct facets")
- `methods/methods_grassmann_cluster_extent.md` §2 (uses "complementary aspects" — clean)
- `tables/beta_per_patient.md` ("different facets of the same geometry")

Substitution menu: aspect, side, dimension, strand, angle, layer.

Fix: same opt-in pattern as item 33.

Verification:
```bash
rg -i 'load.bearing|facet' .agents/preprint/
```

---

## Recommended fix order (highest leverage first)

Each phase is a coherent unit; do one phase fully before the next.

1. **Phase A — verdict cascade (items 1-9)** ~1-2 hours
   The γ_l + δ "weak → strong" sweep. Q3/Q4/V1 in METHODS_AUDIT_ISSUES.
   Biggest correctness gain. Touch 7 files.

2. **Phase B — internal contradiction cleanup (items 10-11)** ~5 min
   Delete stale δ paragraph in VERDICT_LEDGER; optionally renumber decisions.

3. **Phase C — T_G\* number cascade (items 12-14)** ~30 min
   `bands/00_cohort.md` §3 table (Q2); `methods/methods_grassmann_cluster_extent.md`
   head/§5b/§9/§10; `directives/methods_directive_2026-05-19_TG_normalization.md`
   raw values + open-item.

4. **Phase D — anatomy ledger contradictions (items 15-18)** ~30 min
   Top table γ_l/δ anatomy; cross-band table; β cophenet 7th region missing
   in HANDOFF + 00_cohort.

5. **Phase E — index + status drift (items 19-25)** ~30 min
   README/WRITING_GUIDE missing files; preflight path; methods_section_review
   status; methods_revision checklist contradiction; directive filename
   convention; established_results stub.

6. **Phase F — control-language cascade (items 26-32)** ~1-2 hours
   Decision 6/8 citations; disjunctive→mass-only language sweep; C5
   retention→Wilcoxon language sweep; figure stubs; audit-path wildcards.

7. **Phase G — vocabulary sweep (items 33-34)** opportunistic
   Apply during other phases; don't dedicate a pass.

---

## Verification checklist per phase

Before marking a phase done, verify:

- No file in `.agents/preprint/` says "γ_l … weak trace" or "δ … weak
  trace" (rg -i "weak.{0,5}trace.{0,5}only.{0,5}grassmann")
- `obs_mass` values in any table match
  `data/audit/grassmann_cluster_extent/cohort_summary.csv` exactly
- Anatomy region counts in any table match
  `data/audit/anatomy_<band>_<probe>_clusterext/cohort_summary.csv`
- Every "disjunctive gate" reference is either replaced with "mass-only
  gate" or annotated as historical/superseded
- Every "≳80% retention" reference is either replaced with the Wilcoxon
  gate or annotated as historical
- All cited file paths resolve (no wildcards, no broken pointers)

---

## Items NOT in this audit (out of scope)

- LaTeX manuscript bugs (B1-B4, G1-G5 in METHODS_AUDIT_ISSUES) — outside
  `.agents/preprint/`
- Methods-section narrative-level edits — owned by the writing agent
- Figure regeneration — covered by figure scripts in `scripts/02_preprint/`
- CSV regeneration — only `T_G^*` normalized column (C1 in
  METHODS_AUDIT_ISSUES) is pending compute work

---

## Revision history

- **2026-05-26** — Initial audit. 34 items + 7 fix-order phases logged.
  All items have `[ ]` status — user decides each.
- **2026-05-26 to 2026-05-28** — All 34 items + 4 added sub-items (14b–14e)
  closed across 7 phases:
  - **Phase A** (items 1–9): γ_l + δ verdict cascade under Decision 12.
    γ_l promoted to **strong trace, only Grassmann ↑** (LOO 0.040 Pat_05
    passes Decision-12 precondition); δ stays **weak** (LOO 0.055 Pat_08
    fails precondition). New **Decision 12** locked: strong-tier verdict
    requires `cluster_p_mass < 0.01` AND `cluster_p_longest_run < 0.05`
    AND `LOO max p_mass < 0.05` at full data. C5 epi-X explicitly demoted
    to secondary mechanistic observation per Decision 10.
  - **Phase B** (items 10–11): `VERDICT_LEDGER.md` δ section fully
    rewritten; Decision 12 added to decisions list; head note on
    file-order vs chronological.
  - **Phase C** (items 12–14e): T_G* number cascade + C1 normalization
    fully propagated. Cohort uses `n_k^cohort = 111` (full-data) or
    `n_k^epi-X-cohort = 87` (C5 epi-X, Pat_13 limit). Per-patient uses
    `n_k^s`. Dual `raw (normalized)` in locked ledgers; normalized-only
    in briefs + methods + neurophys interpretation.
  - **Phase D** (items 15–18): anatomy ledger updated under locked `S(b)`
    cluster-extent paradigm. γ_l = 7 regions (occipito-temporal + frontal
    + medial-OFC); δ = 4+3 regions fully disjoint; β cophenet 7th region
    (`rostralanteriorcingulate`) added to HANDOFF + cohort brief; KC-era
    retractions documented for left fusiform / Amy / bankssts / caudal
    ACC / paracentral.
  - **Phase E** (items 19–25): index / status drift fixed. README +
    WRITING_GUIDE list `methods_neurophysiological_interpretation` +
    `METHODS_AUDIT_ISSUES`. `tables/beta_per_patient` preflight path
    corrected. `methods_section_review` → `status: mostly_superseded`.
    `methods_revision_cophenet` checklist rewritten to match
    `feedback_no_unmotivated_bh_fdr.md`. Methods directive frontmatter
    annotated. `established_results/` README rewritten as historical Q&A.
  - **Phase F** (items 26–32): control-language cascade. Decision 6 →
    Decisions 6+8 citations; disjunctive → mass-only gate language;
    retention rule retired (Decision 10 Wilcoxon-on-epi-X); audit-path
    wildcards resolved (`audit_*_ctm_triangle.py` → `audit_33_*`);
    `F_cohort_2` stub annotated as deferred.
  - **Phase G** (items 33–34): vocabulary sweep. "load-bearing" → primary
    / central / principal / decisive; "facet" → aspect. All live preprint
    files swept; directives + responses + `methods_section_review`
    skipped as historical record.

**Final state**: 34 of 34 items closed. The preprint folder is internally
consistent under the locked Decision-8 + Decision-12 gate, the C1
normalization, and the `S(b)` cluster-extent anatomy paradigm. No known
open contradictions or stale-value cascades remain.
