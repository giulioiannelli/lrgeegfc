---
name: review-issue-tracker-2026-06-01
era: IMCOH_ABS_COHORT_N10
status: open
kind: response
scope: ranked, actionable defect list backing the 2026-06-01 internal referee review; each item has location, severity, refute status, and a concrete fix; per-item working document, user decides each fix
companion: 2026-06-01_internal-referee-review.md, locked/ANATOMY_LEDGER.md, locked/VERDICT_LEDGER.md, METHODS_AUDIT_ISSUES.md
date: 2026-06-01
---

# Review issue tracker — preprint referee pass (2026-06-01)

**Head.** 14 surviving defects after adversarial refute: **0 blockers, 4 major, 10 minor.**
No defect invalidates a trace verdict; the 4 majors are documentation-cascade + framing items
that **block drawing LaTeX** because the briefs feed the prose. Every number below was recomputed
from on-disk CSVs in `lapbrain`. "Refute" = an independent skeptic agent tried to overturn the
finding; the adjusted severity is post-refute.

Ranking key: **P0** = must fix before any results prose is drawn; **P1** = fix in the same pass;
**P2** = low-cost label/wording; **P3** = backlog / referee-anticipation (not in the 17-check scope).

---

## P0 — block LaTeX drafting (4 major)

### B4 · Cascade the anatomy retraction into the 5 briefs  · **major**
- **What:** `ANATOMY_LEDGER` (2026-05-30) retracts β cophenet 0/7, α cophenet 0/11, γ_l Grassmann
  0/7, δ Grassmann full 0/4 + epi-X 0/3, and softens β Grassmann → Hippocampus-only; its own
  "Cascade TODO" is undone. All 5 consumer docs still assert the retracted anatomy as live fact.
- **Location → edit:**
  - `bands/00_cohort.md` — Head `:29` ("anatomically localized to a band-specific distributed
    cortical network"); verdict matrix `:37-40` (4 "strong localized N-region" cells); §4 `:87,:92-96`;
    §5 cross-band recurrence table `:104-137`. → Reframe to "localization not cohort-resolvable at
    n=10; only Hippocampus@βGrassmann (anchor)"; retract the 4 region-list cells + §5 table.
  - `bands/01_beta.md` — frontmatter `:22-23`, Head `:45`, §5 `:390-451` (7+7 region tables). →
    Retract β cophenet 7-region list; soften β Grassmann to Hippocampus-only as **anchor**.
  - `bands/02_alpha.md` — frontmatter `:21-22`, §5 `:221-276` (11-region list + "9 of 11 reproduce
    under epi-X"). → Retract list as **diffuse**; strike the epi-X-reproduction narrative as moot.
  - `bands/03_gammalow.md` — frontmatter `:21`, §5 `:204-226` (7-region). → Retract list (A3-only,
    phase-avg/anchor).
  - `bands/06_delta.md` — frontmatter `:23-24`, §5 `:188-237` ("FULLY DISJOINT networks", 4+3 regions).
    → Retract both lists + the disjoint framing; flag the no-op epi-X mask (see B2).
  - Also propagate to `HANDOFF_INDEX.md` + `VERDICT_LEDGER.md` anatomy references.
- **Do NOT touch** any trace verdict (β strong/both, α strong/cophenet, γ_l strong/Grassmann, δ
  weak/Grassmann). Add a `revision_history` frontmatter entry dated 2026-06-01 to each brief.
- **Refute:** not refuted — it is the ledger's own outstanding TODO. (confidence high)

### B2 · `audit_72` epi-X mask no-op bug — taints δ Grassmann **anatomy**; fix uncommitted  · **major**
- **What:** `audit_72_anatomy_grassmann.py` HEAD `~:234` does `~np.isin(np.arange(N), <string
  labels>)` (int-vs-str) → excludes 0 nodes. `anatomy_delta_grassmann_epiX/per_patient_trace_flags.csv`
  shows `n_nodes_kept = full N` (Pat_02 117 not 103; Pat_13 119 not 89) vs `grassmann_epi_exclusion/
  epi_counts.csv`. Same bug class as the fixed `audit_71` (commit e97d249) — recurred.
- **Scope (verified):** taints **only** the δ Grassmann epi-X "disjoint networks" *anatomy*
  (`bands/06_delta.md §5`, `00_cohort.md:143`). The **C5 trace verdict is CLEAN** — it comes from
  `audit_67` whose label-based mask actually excludes nodes (`epi_counts.csv` proves it). Do not
  re-touch any `VERDICT_LEDGER` C5 cluster-mass number.
- **Fix:** (1) commit the working-tree label-based fix (HEAD still ships the bug); (2) regenerate
  `anatomy_delta_grassmann_epiX/` and `_clusterext/` and verify `n_nodes_kept = N_reduced`;
  (3) retract the δ "disjoint networks" prose (folds into B4); (4) errata-note that the bug was
  confined to per-region anatomy and never reached the C5 trace gate.
- **Refute:** not refuted; scoping confirmed (`git show HEAD` + grep that C5 reader does no masking).

### B3 · Grassmann anatomy is phase-averaged participation (anchor), not a trace  · **major**
- **What:** `participation_phase_avg()` averages over (rest_pre, task, rest_post) — verified
  `== (pre+task+post)/3` bit-for-bit; node-level ρ ≈ 0.16 vs a differencing trace field. So Grassmann
  "anatomy" answers "where leading modes sit on average" (anchor), not "where the trace localizes".
- **Location → edit:** every Grassmann-anatomy sentence in the briefs + LaTeX. Frame Hippocampus@β
  strictly as anchor-flavored. Retract γ_l/δ Grassmann anatomy lists (folds into B4).
- **Fix (if trace localization is wanted later):** build a differencing per-node quantity (task −
  mean(rest), or contribution to cross-phase d_G) and re-run A1/A3; never reuse `participation_phase_avg`.
- **Refute:** not refuted; mechanism reproduced. Already internalized in `ANATOMY_LEDGER`, not in briefs.

### D1 · "Band resolution NOT amplification" overstated — β is promoted across the gate  · **major**
- **What:** at the matched-strength gate (p<0.05): β raw FC 0.0527 (**fails**) → cophenet 0.0049
  (**passes**) = cross-gate promotion + >10× strengthening; α 0.0137 → 0.0020. Gate-level raw-FC
  detections are only **{δ, α}**, so "raw FC detects everything" conflates `n_above` (6–8/10) with
  gate detection. Only δ is a true demotion; θ/γ_h were never detected.
- **Location → edit:** `bands/00_cohort.md §2 :63`; `HANDOFF_INDEX.md` anti-pattern `:63`.
  → "band-selective sharpening — across-band demotion of δ + within-band strengthening of β/α";
  retire the "NOT amplification" anti-pattern; separate `n_above` count from gate detection in the
  three-layer narration.
- **Refute:** not refuted; every cell recomputed. Severity capped at major (framing, not a wrong result).

---

## P1 / P2 — same-pass fixes (10 minor)

### A3 · Grassmann `cluster_p_mass` floored at 1/(R+1); LOO knife-edge  · **minor** *(was major; refuted-down)*
- β/δ/γ_l all at 0.0049751 (0/200 nulls ≥ obs). Strong-vs-weak rests on the LOO column, which is
  itself floor-quantized (γ_l 8/201 vs δ 11/201) → ~single-surrogate margin. **Fix:** quote "p ≤ 0.005
  at R=200"; **raise R ≥ 2000** to sharpen the bound and stabilize the LOO. *(Refute: manuscript already
  documents the LOO-driven split; the floored-LOO fragility is the genuine residue.)*
- **Location:** `methods_grassmann_cluster_extent.md §5b :451-453`.

### A4 · δ Grassmann CSV `verdict` column stale (`strong` vs ledger `weak`)  · **minor**
- `grassmann_cluster_extent/cohort_summary.csv` δ `verdict_cluster_extent = strong` (pre-Decision-12)
  vs ledger `weak`. **Fix:** regenerate the column under Decision-12 (`cluster_p_mass_loo_max < 0.05`)
  or add a stale-warning header. One-line gate change in `audit_70`.

### A1 · α "intermediate" CSV label + mislabeled `z`/`W`  · **minor**
- α verdict CSV column says `intermediate` (stale 8/10 scheme); `paired_wilcoxon_z` is actually the
  Wilcoxon **W**, not a z-score. The α "strong" verdict itself is correct (gate of record). **Fix:**
  relabel CSV verdict column to the post-2026-05-19 C3 gate; rename `paired_wilcoxon_z` → `_W`; mirror
  α's 5/10 vs β's 7/10 in the ledger row. **Location:** `matched_strength_surrogate_split_baseline/{cohort_summary,README}`,
  `VERDICT_LEDGER.md` α row.

### A2 · β C2 (split>drift) is the marginal control  · **minor**
- C2 p=0.0137 clears 0.05 by ~3× narrower margin than C1/C3 (0.00488). **Fix:** state in limitations
  that β's strong status is most sensitive to the split-vs-drift contrast. No verdict change.

### A5 · three-layer table — framing only  · **minor** → folds into D1 (numbers all verify exactly).

### A6 · multiple comparisons — honest, keep the disclosure  · **minor**
- cophenet β/α survive band-Bonferroni + Holm; Grassmann within-band FWE, no across-band correction
  (the cross-band BH m=6 sentence was correctly deleted). **Fix:** keep the explicit "no across-band
  correction" sentence; optionally note the 3 Grassmann floor-bands would also pass Bonferroni×6.

### C3 · `ctx-rh-insula` certified by Pat_10 (β-anti) alone  · **minor** *(was major; refuted-down)*
- Pat_10 = sole sampler (8 contacts; 9/10 have 0) and most β-anti (obs_rho −0.0906). Direction-blind
  |Δρ| top-decile artifact. **Fix:** strike from `bands/01_beta.md:404` + 2026-05-19 directives (folds
  into B4); for future anatomy use signed `s_ij = dD_task·dD_rest` + n≥5 gate + not-cohort-anti check.
  *(Refute: real, but already retracted at the ledger → residual is stale downstream text.)*

### C4 · live retired-probe citation: "left-fusiform anatomical anchor"  · **minor**
- `bands/01_beta.md:528` cites a γ_l "left-fusiform anatomical anchor" — retracted at every band under
  S(b). **Fix:** replace with the S(γ_l) anatomy ("occipito-temporal + frontal + medial-OFC") or just
  "at narrow k=12..23". (KC/VI/τ-sweep are otherwise correctly framed as retired everywhere.)

### D2 · "two independent probes" → "complementary"  · **minor** *(was major; refuted-down)*
- Both probes are functions of one eigendecomposition of L̂ at τ_max → not statistically independent;
  the double dissociation (3/6 bands disagree) supports **complementary**. **Fix:** word-swap at
  `01_beta.md:263,386,528`; add one sentence acknowledging the shared operator; leave the legitimate
  "independent rsPre half-baselines" / "per-probe gating" uses. *(Refute: science untouched; word choice.)*

### D3 · residual message coherent, source docs stale  · **minor** *(was major; refuted-down)*
- The THAT+WHICH-SCALES message is verified; `00_cohort.md` + neurophys companion still oversell anatomy
  (folds into B4/M1). **Fix:** strike the anatomy clause/tables in `00_cohort.md`; add a retraction
  cross-ref to the neurophys companion §4/§7. *(Refute correctly removed the overreach "audit_72 invalidates
  all Grassmann epi-X" — the C5 trace path via audit_67 is clean.)*

---

## Resolved / not a defect (recorded for the audit trail)

- **B1 · anatomy retraction itself — OK (refuted).** The retraction is sound, reproduced to ~1e-16,
  **already applied** in `ANATOMY_LEDGER`; the one README fix it proposed was already in place. (The
  *outstanding* work is the cascade — tracked as B4.)
- **C1 · Pat_15 / n=9 — OK.** n=10 everywhere; n=9 vestigial & supplementary; Pat_15 inclusion is
  conservative (drop-LOO mass β 184 / γ_l 151). State the conservativeness as a point in the paper's favor.
- **C2 · δ "weak" — OK.** Honest; the "no trace / single-patient" downgrade is **falsified** (all 10 δ
  LOO masses exceed null p95). Keep "cohort gate clears, LOO fails at Pat_08 → weak"; do not let the C5
  epi-X strengthening creep toward re-promotion.

---

## P3 — referee-anticipation backlog (outside the 17-check scope)

Ranked by leverage; from the completeness critic. Not blockers for the β finding, but a strong
reviewer will raise them.

1. **`imcoh_abs` substrate validity** — Jensen order-of-operations (`|mean|` vs `mean|·|`), whether
   abs-rectification preserves Nolte volume-conduction immunity after band-averaging, sign/lag loss.
   *Foundation under both probes.* Audit against `imcoh-guide`.
2. **τ_max = 1/λ_max single-scale choice** — both probes live at the finest scale; λ_max is
   degree-dominated. Sweep τ ∈ {1/λ_max, 1/λ_2, intermediate} for the β probes (`memory/lrg_tau_choice`).
3. **Figure scripts (`scripts/02_preprint/`, 18 files)** — `preprint_10/_11` likely render retracted DK
   region-lists; static-scan all for stale hardcoded numbers + never/always plot-rule violations
   (`axhline` threshold lines, RdBu_r polarity, n/10 reference lines).
4. **α C5 epi-X possibly bug-adjacent** — confirm the `audit_72` int-vs-str mask bug (B2) did NOT also
   corrupt the α cophenet epi-X "8.3×→27.7×" numbers cited in `00_cohort.md` (`alpha_epi_exclusion/`).
5. **γ_l "strong, only Grassmann"** single-probe, C5-LOO fragile (0.159, Pat_05); **γ_h** no-trace on the
   0.05 boundary (0.0597 vs 0.055/0.060 quotes) — could flip on one permutation. Stress-test both.
6. **`METHODS_AUDIT_ISSUES.md` open ledger** — reconcile its Q2/Q3/Q4 desyncs (pre-fix masses in
   `00_cohort.md §3`; γ_l/δ weak-vs-strong drift) against the current `VERDICT_LEDGER` in the same pass.
7. **raw-FC substrate (Result-1)** robustness + d_S/d_P ρ≈0.85–0.95 correlation — the "cophenet sharpens
   a borderline β substrate" narrative is defined relative to this layer.
8. **γ_l/δ Grassmann anatomy localization** — the 2026-05-30 README marks this cell "Deferred" yet the
   per-region CSVs exist (Jun 1) and retraction items 4/5 cite them; reconcile the README self-contradiction.

---

## Tally

| Severity | Count | IDs |
|---|---|---|
| Blocker (surviving) | 0 | — |
| Major (P0) | 4 | B4, B2, B3, D1 |
| Minor (P1/P2) | 10 | A1, A2, A3, A4, A5, A6, C3, C4, D2, D3 |
| Resolved / not a defect | 3 | B1, C1, C2 |
| Backlog (P3, out of scope) | 8 | substrate, τ, figures, α-C5, γ_l/γ_h, methods-ledger, raw-FC, γ_l/δ-anatomy |
