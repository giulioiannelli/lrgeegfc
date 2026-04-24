---
name: coding-rules
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-24
pointers:
  - CLAUDE.md
  - .agents/guides/03_implementation/function-map.md
  - .agents/guides/04_rules/never-always-list.md
  - src/lrg_eegfc/
---

# Coding rules

**Library-first. No reinvention. Scripts call; library computes.**

General-purpose code lives in `src/lrg_eegfc/`. Dataset-specific
experiments, figures, and reports live in `scripts/`. When a helper is
used by ≥2 scripts, it is promoted to the library in the same commit
as its second use.

## Hard rules

1. **Library-first check** — before writing any helper in a script,
   `rg` `src/lrg_eegfc/` for the same work. If it exists, use it.
2. **≥2 callers → library** — any helper used in 2+ scripts is migrated
   to `lrg_eegfc.utils.*` immediately.
3. **No private-copy forks** — don't redefine `load_fc_matrix`,
   `wilcoxon_z`, `bh_fdr`, plot helpers. Import them.
4. **FC-method-agnostic by default** — new helpers take `fc_method` as
   a kwarg; never hardcode `"msc"` / `"imcoh_abs"`.
5. **Config-driven constants** — bands, phases, colors, DPI, nperseg —
   from `lrg_eegfc.config`. Never hardcode.
6. **Frontmatter on every new `.agents/` .md** — schema in
   `frontmatter-schema.md`.

## Tests

- No coverage target. Add a smoke test every time you edit a library
  module. Test path mirrors module: `tests/test_<module>.py`.
- Scripts don't need tests. Library code does.
- If the module you're touching has no test, add a minimal one in the
  same commit.

## File size

- Library modules under ~400 lines. Split by concern.
- Scripts focused: one question, one driver, optional local helpers.
  Past ~500 lines, split the script.

## Comments and docstrings

Apply the renormalization principle:
- **Docstring** = head. One line saying what the function does and why.
- **Body comments** only where the WHY is non-obvious (hidden
  constraint, subtle invariant, workaround for a specific bug).
- Never explain WHAT the code does if the identifiers make it clear.
- Never narrate "added for the Y flow" or "used by X" — those belong in
  the PR description, not in code that will outlive them.

## Duplication triage

When finding duplication, don't silently fix it in a drive-by. Flag it:
1. Name the duplicate locations.
2. Propose a library home.
3. Migrate + update callers in one commit with a clear message.

## Renormalized commit messages

- **Subject** = head. `<type>: <juice>` in ≤ 72 chars.
- **Body** = expansion. Why, not what. One paragraph is plenty.
- Trailer: `Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>`.

## Never / always — see CLAUDE.md

Hard enforcement rules live in CLAUDE.md so every agent session loads
them. This file holds the *rationale*; CLAUDE.md holds the list.
