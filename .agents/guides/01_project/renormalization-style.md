---
name: renormalization-style
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-24
pointers:
  - CLAUDE.md
  - .agents/guides/04_rules/coding-rules.md
---

# Renormalization style

**Head first. Expand on demand.**

Every piece of written output — reports, commits, diary entries, skill
replies, PRs, user-facing chat — leads with a 1–2 sentence *juice tag*
that a tired reader grasps in 5 seconds. Technical detail lives below,
only as deep as asked.

## Why

Signals we had were buried inside 400-line reports. Readers bounced off.
One month got spent rerunning scalar tests whose answer was already
sitting in the VI(k) / partition-multiscale outputs. The remedy is
structural: force the *point* to the top, keep the *evidence* below,
and trust that whoever needs depth will scroll or ask.

## The pattern

```
<HEAD — one or two sentences. What the reader should walk away with.>

<optional: headline number / date / pointer>

## Body (only when needed)

<technical expansion, evidence, caveats>
```

## Examples

**Bad — buried lead**
> This report compiles results from the H1–H4 hypothesis tests run
> under the |ImCoh| framework across 9 patients. The cohort was locked
> on 2026-04-22 [30 more lines before any finding].

**Good — renormalized**
> **VI(k) unanimity passes cohort-wide for H1 (θ) at n=9. H2a demoted.**
>
> Post-reset (|ImCoh|, 9 patients, locked 2026-04-22). Full tables below.

## Where it applies

- **Reports** — first paragraph = head; sections below = body.
- **Commit messages** — subject = head; body = expansion.
- **Diary entries** — first bold line = head; blocks below = body.
- **Plans** — Context section = head; phases = body.
- **Memory files** — description field = head; content = body.
- **User-facing chat** — the first sentence = head; tool calls expand.

## Zoom-in trigger

If the user says "expand", "technical details", or something diverges
from expectation — drop one level of abstraction and expose mechanism.
Do not volunteer depth that wasn't asked for.

## Anti-pattern

A crisp head over a hand-waved body is worse than an honest long-form
draft. Head-first is a *summary contract*: the head must be accurate,
the body must actually support it.
