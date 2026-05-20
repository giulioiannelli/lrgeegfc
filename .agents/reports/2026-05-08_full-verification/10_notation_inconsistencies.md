---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
scope: small notation drift across notes_imcoh.tex; collapse to single canonical symbol per object
---

# Notation inconsistencies — fix list

**Head.** Several physical objects appear under two or three different symbols across the manuscript: the merge-height profile (`Δ_n(τ)` vs `t_n`), the communication distance (`\distm` vs `D` vs `D(τ)`), the propagator (`\propag(τ)` vs inline `e^{-τL}`), and the symbol `n_trace` (count of patients vs count of leaves). None of these change the physics; they create reader friction. Pick one symbol per object and apply globally.

## 1. Merge-height profile

**Two names for the same object.**

- §3.1.5 (eq. 6, "Ψ profile" definition): `Δ_n(τ)` — the linkage merge height at index n.
  > `\Psi(n;\tau) = \mathcal{N}[\log_{10}\Delta_n(\tau) - \log_{10}\Delta_{n+1}(\tau)]`
- §5.1 prose ("merge-height profile t_n decays smoothly..."): `t_n` — same object.
- Fig. 20 caption: `t_n` (right panel y-axis label is "merge height t_n").

**Recommendation.** Pick `t_n` (more readable; consistent with Fig. 20 axis label). If `Δ_n(τ)` stays in eq. 6 for parallelism with the original Villegas 2025 notation, add a one-line bridge after the equation: "Hereafter we write `t_n := Δ_n(τ')` for the merge-height profile at the analysis scale `τ' = 1/λ_max`." Then use `t_n` everywhere in §5.1 prose and figure captions.

## 2. Communication distance

**Three notations across §3, §4, §5.**

- §3.1 eq. 5 introduces `\distm` (e.g. `\distm_{ij}(\tau) = (1 - δ_{ij})/K_{ij}(\tau)`).
- §4 prose treats the substrate edge weight directly without defining a separate `D` — but the §4 figures use `D` informally in axis labels ("D = D_taskTest − D_rsPre,A").
- §5 head prose: "the communication distance `D(τ)`" — collapse-form without the `\distm` macro.
- Fig. 25 axis labels: `D = D_taskTest − D_rsPre,A` — uses `D` not `\distm`.

**Recommendation.** Define one canonical symbol `\mathcal{D}(\tau)` (or `\distm(\tau)` if the macro is preferred) at first introduction in §3.1 and use it in every subsequent section / figure / caption. The `D` in Fig. 25 axis labels and §5 prose should be `\mathcal{D}` for consistency.

## 3. Propagator

**Two interchangeable forms.**

- §3.1 introduces `\propag(\tau)` as the propagator (matrix-valued operator).
- Several inline formulas write `e^{-\tau L}` directly without invoking the macro.
- §3.1 eq. 4: `\hat{\rho}(\tau) = e^{-\tau L} / \mathrm{Tr}[e^{-\tau L}]` — uses inline form, not the macro.
- §5.1 head: `K(\tau) = e^{-\tau L}` — uses K plus inline, no `\propag`.

**Recommendation.** Adopt one of:
- (a) use `K(\tau)` everywhere as the kernel, retire `\propag` macro;
- (b) use `\propag(\tau)` everywhere as the macro, drop `K(\tau)` and `e^{-\tau L}` inline forms.

Either choice is fine, but pick one and edit the inline `e^{-\tau L}` instances to match. Recommend (a) because `K(\tau)` is already used in eq. 5 and the §5 prose.

## 4. Phase variable

**Three notations.**

- `\Phi` (capital phi, generic phase variable in §5.3): `\rho_S(\Delta_{\mathrm{task}}, \Delta_{\mathrm{rest}})` and "for each phase Φ".
- `\acrshort{rspre}` / `\acrshort{rspost}` etc. (acronym macros, most uses in body prose).
- `\txtacr{rspre}` (text-acronym variant, appears in some figure captions).

**Recommendation.** Standard pattern is:
- `\Phi` only when needed as a generic / dummy variable (e.g. "for each phase Φ ∈ {rsPre, taskLearn, taskTest, rsPost}").
- `\acrshort{rspre}` (or just `rsPre` typeset with the canonical short form) for explicit phase references in body prose.
- The `\txtacr{...}` variant in figure captions should be replaced with the same `\acrshort{...}` form to match body prose. Check captions of Fig. 12 / 13 / 14 for the `\txtacr` artifact.

## 5. Greek typography in tables / captions

**Mostly fine in body, drift in captions / tables.**

- Body prose: `\beta`, `\alpha`, `\gamma_l`, etc. — math mode, consistent.
- Some figure captions use plain "beta" or "Beta" inline. Check Fig. 1 caption ("β band (rsPre phase)") — typeset is OK, but the legend renders "beta" as a plain word in some panels.
- Table 1 (band names): typeset uses Greek capitals ("Delta", "Theta", "Alpha", "Beta", "Low gamma", "High gamma") in the **Band** column and `delta`, `theta`, etc. typewriter in the **Label** column. This is intentional (label column is the literal config string). Keep, but prose references should use `\delta` / `\theta` / `\alpha` / `\beta` / `\gamma_l` / `\gamma_h` consistently.

**Recommendation.** Audit figure-caption Greek letters and tabular headers for plain-text vs `$\beta$` consistency. Body prose is fine.

## 6. `n_trace` overloaded

**Two distinct meanings, same notation.**

- §4 (Tab 2, Fig 14, Fig 15): `n_trace^{(d_S)}` = number of patients with `T_d^{(d_S)} < 0` (cohort-aggregate count of trace-direction patients per band, max 10).
- §5.5 prose (anatomy): "trace-leaf count `n_trace(r)`" = count of trace-leaves in region r within the cohort cortical contact pool (typically 5–13 for the surviving cells).

These are unrelated: one counts patients, one counts leaves. The notation `n_trace` is the same.

**Recommendation.** Disambiguate via subscript:
- `n_trace^{patient}` or `n_pat,trace` for the per-band patient count of §4 and §5.3.
- `n_trace^{leaf}(r)` or `n_leaf,trace(r)` for the per-region leaf count of §5.5.

The §4 use is most natural as `n_pat,trace` because Tab 2 already uses that pattern visually with the "Indication" column reading "trace" / "drift-only" / "borderline" alongside the patient count. The §5.5 use is most natural as `n_trace(r)` (with the explicit region argument disambiguating it).

Alternatively keep `n_trace` everywhere but add the explicit context on first use of each subsection: §5.5 already does "trace-leaf count" prose immediately around the symbol; §4 does not, but the table column header could read "n_pat with T_d < 0 / 10" instead of bare `n_trace^{(d_S)}` to make the patient interpretation explicit.

## 7. Trace — covered separately

The word "trace" carries three senses (direction / module / leaf). See `11_terminology_disambiguation.md` for the full rewrite list. Not duplicated here.

## Action priority

1. **High** — `n_trace` overload (§4 ↔ §5.5): the same symbol means "patient count" in §4 and "leaf count" in §5.5. Disambiguate.
2. **Medium** — merge-height profile `Δ_n(τ)` vs `t_n`: pick one and unify.
3. **Medium** — communication distance `\distm` vs `D` vs `D(τ)`: pick one (`\mathcal{D}(\tau)` recommended) and unify.
4. **Low** — propagator `\propag(τ)` vs inline `e^{-\tau L}`: pick one.
5. **Low** — `\acrshort` vs `\txtacr` in figure captions: audit and unify.
6. **Cosmetic** — Greek letters in figure captions / tabular headers: audit.

## Action — writing agent

Apply globally. Use `find` + grep on the source `.tex` to enumerate every site of each symbol before unifying:
```
grep -n "Delta_n\|t_n\|\\\\distm\|\\\\propag\|n_trace\|n_{\\\\rm trace}" notes_imcoh.tex
```

After the unification pass, re-check the §5 figure-caption variants of `\Phi`, `\Phi'`, and the `T_KC(\lambda)` notation block in §5.2 (the `\widetilde{T}_{KC}(\lambda)` variant in Fig. 22 caption uses both forms; that one is intentional and should stay — the tilde indicates cohort median).
