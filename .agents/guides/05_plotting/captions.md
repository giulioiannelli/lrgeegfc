---
name: 05_plotting / captions
type: style_sheet
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
updated: 2026-04-29
---

# Captions — short, plain, human

## When to write one

Only when explicitly asked, or when the destination folder already
accumulates `<figure>.md` sidecars (e.g. `data/audit/...` audit folders).
Captions are NOT a default by-product of every figure.

## How to write one

Plain English. Talk to the reader, not at them. Keep it short.

### Three blocks, in order

1. **One short sentence** — what is this figure?
   Read like a colleague pointing at it: "Per-band scatter showing
   whether rest_post sits closer to task than rest_pre does, one
   point per patient."
2. **Layout** — bullet list, plain words.
   "Rows are distances. Columns are bands. Each panel: x-axis is
   pre→task distance, y-axis is task→post distance."
3. **Read** — one short paragraph, the take-away.
   "Below the dashed identity line means persistence. β has the
   tightest cluster below; θ scatters around the line."

That's it. Three blocks.

## Style rules

- **Sentence first, paragraph second.** Never lead with a heading.
- **Plain language.** "Spread of points" beats "interquartile range
  of the distribution"; "persistence" beats "T_d < 0 in cohort
  median". When you must use jargon, define it once on first use.
- **No copy of the methods section.** The reader knows what `d_S` is
  by the time they reach the caption. Don't restate the formula —
  point at the file that has it.
- **Specific over generic.** "α has 8/10 patients below the line"
  beats "most bands show the effect".
- **No bullet salad.** Two or three bullets max. Combine related
  points into a sentence.
- **No emoji.** No apologetic hedges ("note that…", "it's worth
  pointing out…"). Just say it.
- **No suptitle or repeat of the figure title.** The title is on the
  figure (or in the folder context). The caption is the *read*.
- **Length.** Aim for half a screen. If it's longer than the figure
  takes to look at, it's too long.

## Anti-patterns

- ❌ Multi-section caption with `## Methods`, `## Results`, `## Discussion`.
  The caption is not the paper.
- ❌ "This figure shows…" as the opening clause. Cut it. Just say
  what it is.
- ❌ Repeating every panel's content as a separate bullet. The
  reader can see them. Layout = how to navigate, not what's there.
- ❌ A "Read" block that just restates layout. Tell me what to look at.
- ❌ Future-tense plans inside a caption ("we will run the cohort
  sweep next…"). That belongs in the report, not on the figure.

## Tiny example

> **Per-band persistence scatter.** One point per patient. Below the
> dashed line means rest_post is closer to task than rest_pre is.
>
> - Rows: nothing (single row).
> - Columns: 6 bands, δ θ α β low_γ high_γ.
> - Blue circles = 9 in-pool patients; orange triangle = Pat_03.
>
> α has the tightest cluster below the line (8/10). β and low_γ are
> next. θ straddles the line — no persistence. high_γ goes the wrong
> way for the cohort median.

That is enough. Stop there.
