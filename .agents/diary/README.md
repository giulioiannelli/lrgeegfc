---
name: diary-readme
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-24
pointers:
  - .agents/guides/04_rules/naming-conventions.md
  - .agents/guides/01_project/renormalization-style.md
---

# Diary — how to use it

**One file per day: `.agents/diary/YYYY-MM-DD.md`. Each session appends
a stamped renorm-style block. Short head on top, body below.**

## Daily file template

```markdown
---
name: diary-YYYY-MM-DD
type: diary
era: IMCOH_ABS          # or whichever applies to the work
status: current
created: YYYY-MM-DD
updated: YYYY-MM-DD
---

# YYYY-MM-DD

## Session HH:MM–HH:MM — <one-line juice>

**Context:** <what we were trying to do, ≤ 1 sentence>
**Done:** <result in one line>
**Issues:** <what went sideways or blocked>
**Decisions:** <what was locked in>
**Files touched:** <short list or pointer to commit SHA>
**Next:** <one pointer for the next session>

## Session HH:MM–HH:MM — <next session juice>
...
```

## Rules

- **Append, don't overwrite.** Each session adds a new `##` block.
- **Head = one-line juice.** If a reader skips the body, the juice must
  carry the day.
- **Don't duplicate** what's already in reports or commit messages —
  point to them instead.
- **Bad day → write it.** Post-mortems live in `.agents/reports/`, but
  the diary records day-by-day reality, warts and all.

## Backfilling

Backfill is allowed for historical days — mark the block with
`backfilled: YYYY-MM-DD` on a line under the session header. Source
from `git log` + report update times. Don't invent details you don't
have.

## The `/diary` skill (once installed)

Runs the `Append` flow: reads today's file, prints the template with
date + current era pre-filled, and writes the new block.
