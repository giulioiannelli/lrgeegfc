---
name: preprint-writing-guide
era: IMCOH_ABS_COHORT_N10
status: current
kind: guide
scope: where new files go inside `.agents/preprint/` — routing rules for agents adding to the preprint folder
---

# `.agents/preprint/` — writing guide for agents

**Head.** Don't add new files at the top level. The top level is reserved
for four fixed entry points (`README.md`, `WRITING_GUIDE.md`,
`HANDOFF_INDEX.md`, `EVALUATION_PROTOCOL.md`); everything else lives in
one of five subfolders by **kind**. When in doubt, **update an existing
file** rather than creating a new one — the locked artifacts
(`locked/VERDICT_LEDGER.md`, `locked/CONTROLS.md`,
`locked/ANATOMY_LEDGER.md`, `locked/ANATOMY_CONTROLS.md`) and the seven
per-band briefs (`bands/00_cohort.md` + `bands/01_beta.md` …
`bands/06_delta.md`) are intentionally fixed; revisions go into them, not
beside them.

---

## Folder layout (current)

```
.agents/preprint/
├── README.md               # comprehensive index (read first)
├── WRITING_GUIDE.md        # this file — routing rules
├── HANDOFF_INDEX.md        # LaTeX writing-agent entry point
├── EVALUATION_PROTOCOL.md  # LaTeX-output verification protocol
│
├── locked/                 # SOURCES OF TRUTH — append only via dated revision
│   ├── VERDICT_LEDGER.md   #   per-band locked trace verdicts + decision log
│   ├── CONTROLS.md         #   C1–C5 control battery (trace)
│   ├── ANATOMY_LEDGER.md   #   per-(band, probe) locked anatomy verdicts
│   └── ANATOMY_CONTROLS.md #   A1–A4 control battery (anatomy)
│
├── bands/                  # per-band briefs + cohort synthesis (fixed set)
│   ├── 00_cohort.md
│   ├── 01_beta.md
│   ├── 02_alpha.md
│   ├── 03_gammalow.md
│   ├── 04_theta.md
│   ├── 05_gammah.md
│   └── 06_delta.md
│
├── methods/                # long-lived methods companions (one per probe family)
│   ├── methods_revision_2026-05-18_cophenet.md
│   ├── methods_grassmann_cluster_extent.md
│   └── methods_section_review_2026-05-19.md
│
├── directives/             # writing-agent directives (dated, accumulate)
│   └── writing_directive_YYYY-MM-DD_<topic>.md
│
├── responses/              # cascade summaries, preflight notes, replies (dated)
│   └── YYYY-MM-DD_<topic>.md
│
├── established_results/    # historical methodology Q&A (frozen)
└── tables/                 # rendered tables (.md + .tex)
```

---

## Where do I put …?

Decision tree for adding a new file. **Default action is "update an existing
file"** — only create a new file when the decision tree explicitly says so.

### A trace verdict change → don't create a file
- Add a **dated revision entry** at the bottom of `locked/VERDICT_LEDGER.md`.
- Cite the new audit + CSV row.
- Update the corresponding `bands/NN_<band>.md` brief in the same edit.

### A new control or control redefinition → don't create a file
- Update `locked/CONTROLS.md` directly (numbered `Decision N` block).
- If the control replaces an older one, mark the older one as superseded
  in-place, don't delete it.

### An anatomy verdict change → don't create a file
- Add a dated revision entry to `locked/ANATOMY_LEDGER.md`.
- Update `locked/ANATOMY_CONTROLS.md` only if the battery itself changes.

### Per-band findings (any band) → don't create a file
- Update the corresponding `bands/NN_<band>.md`.
- For cross-band synthesis updates, edit `bands/00_cohort.md`.
- Never create `bands/02_alpha_v2.md` or `bands/01_beta_addendum.md` —
  the per-band brief is the single home for that band.

### A new methods directive (probe-level) → `methods/` (new file)
- Create `methods/methods_<topic>_YYYY-MM-DD.md` only when a fundamentally
  new probe or methods convention is being introduced.
- For methods *corrections* to an existing companion, edit in-place.

### A writing-agent directive → `directives/` (new file)
- Filename: `writing_directive_YYYY-MM-DD_<short-topic>.md`.
- One file per directive sent to the LaTeX-producing writing agent.
- These accumulate over time — they are not consolidated.

### A response to writing-agent feedback / cascade summary / preflight note → `responses/` (new file)
- Filename: `YYYY-MM-DD_<short-topic>.md` (or
  `feedback_response_YYYY-MM-DD_<topic>.md` for explicit feedback replies).
- These document what changed *in response to* writing-agent feedback.
- Locked artifacts (verdicts, controls, methods) updated as a result of
  the response live in their canonical home, not in the response file —
  the response file points to them.

### A figure → not here
- Per-band figures: `data/preprint/figures/<band>/`.
- Cohort figures: `data/preprint/figures/cohort/`.
- The brief in `bands/` references the figure paths; it does not embed them.

### Anything else → ask before creating
- If the file doesn't fit any of the above categories, prefer updating an
  existing file. If nothing fits, ask the user explicitly before creating
  a new top-level category.

---

## Frontmatter (mandatory)

Every new `.md` file under `.agents/preprint/` begins with YAML
frontmatter:

```yaml
---
name: <short-kebab-case-slug>
era: IMCOH_ABS_COHORT_N10
status: <current | superseded | draft>
kind: <verdict | controls | brief | methods | directive | response | guide>
scope: <one line on what this file covers>
---
```

`era` is always `IMCOH_ABS_COHORT_N10` for files created in this era.
`kind` is one of the seven values above (drives the routing rule that put
the file in its subfolder).

---

## Cross-reference conventions

Inside `.agents/preprint/`, link with **relative paths**:

- From a subfolder to top level: `[X](../X.md)`.
- From a subfolder to another subfolder: `[X](../target/X.md)`.
- From a subfolder to same subfolder: `[X](X.md)`.
- From top level to a subfolder: `[X](subfolder/X.md)`.

For inline-text references (backtick form), display the subfolder prefix
when it helps navigation (e.g., `` `locked/CONTROLS.md` `` from a brief),
omit the prefix for same-subfolder references.

From **outside** `.agents/preprint/` (e.g., `.agents/reports/`,
`scripts/`), always use the full path: `.agents/preprint/locked/CONTROLS.md`.

---

## Anti-patterns (never do)

- Don't create a new file at the top level — only the four fixed entry
  points live there.
- Don't create `<existing-file>_v2.md`, `<existing-file>_addendum.md`,
  `<existing-file>_update_YYYY-MM-DD.md`. Edit in place.
- Don't create a new per-band brief outside the fixed 7-file set
  (`bands/00_cohort.md` + `bands/0[1-6]_<band>.md`).
- Don't create a parallel locked file (e.g., `locked/VERDICT_LEDGER_v2.md`,
  `locked/CONTROLS_NEW.md`). The locked files are append-only.
- Don't put dated responses or directives in `locked/` or `bands/`. They
  go in `responses/` or `directives/`.
- Don't move files between subfolders without updating cross-references
  (use `git mv` if tracked).

---

## Reorganization history

- **2026-05-19** — flat-folder → 5-subfolder layout
  (`locked/`, `bands/`, `methods/`, `directives/`, `responses/`).
  Pre-reorg top level had ~24 files mixed by kind; reorg consolidated
  by kind with the four fixed entry points retained at top level.
  All cross-references rewritten via `/tmp/update_preprint_refs.py` +
  `/tmp/update_preprint_backticks.py` (one-shot, not committed).

## Source-of-truth references

- [`README.md`](README.md) — comprehensive index (what's here, what's load-bearing).
- [`HANDOFF_INDEX.md`](HANDOFF_INDEX.md) — LaTeX writing-agent entry.
- [`EVALUATION_PROTOCOL.md`](EVALUATION_PROTOCOL.md) — verification of writing-agent output.
- [`locked/VERDICT_LEDGER.md`](locked/VERDICT_LEDGER.md) + [`locked/CONTROLS.md`](locked/CONTROLS.md) — trace lockdown.
- [`locked/ANATOMY_LEDGER.md`](locked/ANATOMY_LEDGER.md) + [`locked/ANATOMY_CONTROLS.md`](locked/ANATOMY_CONTROLS.md) — anatomy lockdown.
