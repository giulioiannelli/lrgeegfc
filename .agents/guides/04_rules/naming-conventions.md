---
name: naming-conventions
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-24
pointers:
  - .agents/guides/04_rules/frontmatter-schema.md
---

# Naming conventions

**kebab-case. Date prefix where chronology matters. Frontmatter for
semantics.**

## Files

| Type | Pattern | Example |
|------|---------|---------|
| report (dated) | `YYYY-MM-DD_<slug>.md` | `2026-04-15_imcoh-reset.md` |
| plan (dated) | `YYYY-MM-DD_<slug>.md` | `2026-04-25_surface-multiscale-trace.md` |
| diary | `YYYY-MM-DD.md` | `2026-04-24.md` |
| guide | `<slug>.md` | `imcoh-guide.md` |
| rule | `<slug>.md` | `coding-rules.md` |
| post-mortem | `YYYY-MM-DD_post-mortem-<slug>.md` | `2026-04-24_post-mortem-scalar-session.md` |

Slugs are **kebab-case** (lowercase, hyphens). No `CamelCase`,
`snake_case`, or `UPPER_SNAKE` for `.md` files.

## Exceptions

These filenames stay as-is because tooling expects them:
- `README.md`, `START_HERE.md`, `MANIFEST.md`
- `CLAUDE.md`, `AGENTS.md`

## Folders

- Kebab-case when descriptive: `04_rules/`, `figures_embedded/`.
- Numeric prefix (`01_`, `02_`) allowed when ordering matters.
- Date prefix (`2026-04/`) is used for archive buckets only.

## Python

- Snake_case modules, functions, variables.
- PascalCase classes.
- Follow the existing `lrg_eegfc/` package conventions.

## Archive buckets

Dead files move to `<parent>/archive/YYYY-MM/`. Era tags live in
frontmatter, not filenames. The YYYY-MM bucket + frontmatter together
communicate lineage.

## Cross-references

Relative paths in markdown:
`[VI results](../reports/2026-04-24_h1-h4-vi-results.md)`.
After any rename, `rg` for orphaned references and fix.
