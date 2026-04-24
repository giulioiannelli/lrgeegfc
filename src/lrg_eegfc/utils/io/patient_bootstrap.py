"""Patient directory layout normalization.

Turns a mix of legacy / vendor / new-raw patient directories into the single
canonical layout agreed for this project:

    pat_NN/
    ├── resting/
    │   ├── rest_pre.mat
    │   └── rest_post.mat
    ├── task/                          (optional; absent for pat_06)
    │   ├── task_learn.mat
    │   └── task_test.mat
    ├── implant/                       (optional; absent for pat_06)
    │   └── implant_pat_NN.xlsx
    ├── channel_labels.csv
    ├── implant_pat_NN.csv
    └── provenance.md

The module computes a *plan* (list of MigrationAction) first. Nothing is
written until ``apply_plan`` is called.

Planning invariants (important for correctness of apply):

- We never rename vendor subdirs directly. We ``mkdir`` the canonical
  target subdirs, move files into them one by one, and ``rmdir`` the
  (now-empty) vendor subdirs at the end. This avoids path-invalidation
  bugs where a subdir rename would make subsequent file paths stale.
- Every ``src`` in the plan is a path that exists *before* any action
  runs, and every ``dst`` is a path that exists *after* its action runs.
  Thanks to the invariant above, a straightforward top-to-bottom apply
  always works.
"""

from __future__ import annotations

import csv
import hashlib
import logging
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Canonical target specification
# ---------------------------------------------------------------------------

#: Canonical snake_case phase keys in their new form. The order matches
#: ``PHASE_LABELS`` post-refactor (rest_pre, task_learn, task_test, rest_post).
CANONICAL_PHASES: Tuple[str, str, str, str] = (
    "rest_pre",
    "task_learn",
    "task_test",
    "rest_post",
)

#: Mapping from canonical phase -> (subdir, target_basename).
CANONICAL_PHASE_LAYOUT = {
    "rest_pre":   ("resting", "rest_pre.mat"),
    "rest_post":  ("resting", "rest_post.mat"),
    "task_learn": ("task",    "task_learn.mat"),
    "task_test":  ("task",    "task_test.mat"),
}

#: Vendor / legacy filename candidates for each canonical phase, searched in
#: this order. Paths are relative to the patient root. First hit wins.
PHASE_SOURCE_CANDIDATES = {
    "rest_pre": [
        # canonical-subdirs vendor shape (Pat_05/06/07/08/13/14/15)
        "0.53to300Hz_resting/PRE_STIM_resting_Pre_Data.mat",
        # Pat_10 anomaly (no underscore)
        "0.53to300Hz_resting/PRE_STIM_restingPre_Data.mat",
        # Legacy-flat Pat_02 / Pat_03
        "Resting_PreTask_Data.mat",
        "RestingPreTask_Data.mat",
    ],
    "rest_post": [
        "0.53to300Hz_resting/PRE_STIM_resting_Post_Data.mat",
        "0.53to300Hz_resting/PRE_STIM_restingPost_Data.mat",
        "Resting_PostTask_Data.mat",
        "RestingPostTask_Data.mat",
    ],
    "task_learn": [
        "0.53to300Hz_task/PRE_STIM_task_learning_Data.mat",
        "Task_learning_Data_TimeSeries.mat",
    ],
    "task_test": [
        "0.53to300Hz_task/PRE_STIM_task_test_Data.mat",
        "Task_test_Data_TimeSeries.mat",
    ],
}

#: Vendor subdirs that become empty after migration and must be removed.
OBSOLETE_VENDOR_SUBDIRS = ("0.53to300Hz_resting", "0.53to300Hz_task", "Implant_locations")

#: Candidate locations for duplicate channel-label files that must be
#: verified-and-deleted (the canonical one lives at the patient root).
CHANNEL_LABELS_DUP_CANDIDATES = [
    "0.53to300Hz_task/channel_labels.txt",
    "0.53to300Hz_resting/channel_labels.txt",
    "0.53to300Hz_task/channel_labels.csv",
    "0.53to300Hz_resting/channel_labels.csv",
    "task/channel_labels.txt",
    "resting/channel_labels.txt",
    "task/channel_labels.csv",
    "resting/channel_labels.csv",
    "channel_labels.txt",
]

#: Vendor per-phase INFO.txt files → canonical targets next to their .mat.
#: The keys are patterns relative to the patient root.
INFO_TXT_LAYOUT = {
    "0.53to300Hz_resting/resting_Pre_INFO.txt":  ("resting", "rest_pre.info.txt"),
    "0.53to300Hz_resting/resting_Post_INFO.txt": ("resting", "rest_post.info.txt"),
    "0.53to300Hz_task/task_learning_INFO.txt":   ("task",    "task_learn.info.txt"),
    "0.53to300Hz_task/task_test_INFO.txt":       ("task",    "task_test.info.txt"),
}


# ---------------------------------------------------------------------------
# Plan / action data
# ---------------------------------------------------------------------------

@dataclass
class MigrationAction:
    """One step in the migration plan.

    ``kind`` ∈ {mkdir, rename, generate_csv, generate_labels, verify_labels,
    write_provenance, delete, rmdir, noop}.
    """

    kind: str
    src: Optional[Path]
    dst: Optional[Path]
    notes: str = ""


@dataclass
class MigrationPlan:
    patient: str
    patient_dir: Path
    actions: List[MigrationAction] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def add(self, kind: str, src: Optional[Path], dst: Optional[Path], notes: str = "") -> None:
        self.actions.append(MigrationAction(kind=kind, src=src, dst=dst, notes=notes))

    def warn(self, msg: str) -> None:
        self.warnings.append(msg)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _patnum(patient: str) -> int:
    try:
        return int(patient.split("_")[-1])
    except ValueError:
        return 0


def _first_existing(root: Path, candidates: Iterable[str]) -> Optional[Path]:
    for c in candidates:
        p = root / c
        if p.exists():
            return p
    return None


def _locate_implant_xlsx(patient_dir: Path, n: int, nn: str,
                         plan: MigrationPlan) -> Optional[Path]:
    """Find the implant xlsx among many vendor naming variants.

    Search order:
    1. Exact canonical names at patient root or under Implant_locations/.
    2. Lowercase / mixed-case variants of ``implant_pat_{n}.xlsx``.
    3. Fallback: the single .xlsx in ``Implant_locations/`` if unique,
       with a warning explaining we had to guess.
    """
    named_candidates = [
        patient_dir / f"Implant_pat_{n}.xlsx",
        patient_dir / f"Implant_pat_{nn}.xlsx",
        patient_dir / f"implant_pat_{n}.xlsx",
        patient_dir / f"implant_pat_{nn}.xlsx",
        patient_dir / "Implant_locations" / f"Implant_pat_{n}.xlsx",
        patient_dir / "Implant_locations" / f"Implant_pat_{nn}.xlsx",
        patient_dir / "Implant_locations" / f"implant_pat_{n}.xlsx",
        patient_dir / "Implant_locations" / f"implant_pat_{nn}.xlsx",
        patient_dir / "implant" / f"implant_pat_{nn}.xlsx",
    ]
    for cand in named_candidates:
        if cand.exists():
            return cand

    # Last-ditch: unique xlsx in Implant_locations/ with an unexpected stem.
    impl_dir = patient_dir / "Implant_locations"
    if impl_dir.is_dir():
        xlsxs = sorted(impl_dir.glob("*.xlsx"))
        if len(xlsxs) == 1:
            plan.warn(
                f"non-canonical implant xlsx stem: found "
                f"{xlsxs[0].name} (expected implant_pat_{nn}.xlsx); "
                f"renaming it during migration"
            )
            return xlsxs[0]
        if len(xlsxs) > 1:
            plan.warn(
                f"multiple xlsx in Implant_locations/: {[p.name for p in xlsxs]} "
                f"(cannot disambiguate; add --force or rename by hand)"
            )
    return None


# ---------------------------------------------------------------------------
# Plan construction
# ---------------------------------------------------------------------------

def plan_migration(patient_dir: Path) -> MigrationPlan:
    """Build a migration plan for one patient directory.

    Algorithm (no in-place dir renames):

    1. For each canonical phase, locate the vendor .mat and emit
       ``mkdir resting/`` / ``mkdir task/`` (deduped) + ``rename`` to
       the canonical basename inside the canonical subdir.
    2. Locate the xlsx anywhere; emit ``mkdir implant/`` + ``rename`` to
       ``implant/implant_pat_NN.xlsx``.
    3. Generate ``implant_pat_NN.csv`` from the xlsx; delete mixed-case
       predecessors (``Implant_pat_NN.csv`` / ``Implant_pat_N.csv``).
    4. Consolidate channel labels into a single ``channel_labels.csv``;
       verify and delete all ``.txt`` copies.
    5. Delete stale phase symlinks (``rsPre.mat`` etc.) and duplicate
       ``channel_labels.{txt,mat}`` at root.
    6. ``rmdir`` any now-empty vendor subdirs (``0.53to300Hz_*``,
       ``Implant_locations``).
    7. Write ``provenance.md`` last.
    """
    patient = patient_dir.name
    plan = MigrationPlan(patient=patient, patient_dir=patient_dir)
    n = _patnum(patient)
    nn = f"{n:02d}"

    resting_dir = patient_dir / "resting"
    task_dir = patient_dir / "task"
    implant_dir = patient_dir / "implant"
    mkdirs_scheduled: set[Path] = set()

    def schedule_mkdir(target: Path, note: str) -> None:
        if target in mkdirs_scheduled or target.exists():
            return
        plan.add("mkdir", None, target, note)
        mkdirs_scheduled.add(target)

    # 1. Phase .mats → canonical subdir / canonical basename.
    for phase, (subdir, target_name) in CANONICAL_PHASE_LAYOUT.items():
        dst_subdir = patient_dir / subdir
        dst = dst_subdir / target_name
        # Idempotent: if already in canonical location, no action.
        if dst.exists():
            continue
        src = _first_existing(patient_dir, PHASE_SOURCE_CANDIDATES[phase])
        if src is None:
            has_vendor_task = (patient_dir / "0.53to300Hz_task").exists() \
                              or (patient_dir / "task").exists()
            if phase in ("task_learn", "task_test") and not has_vendor_task:
                plan.add("noop", None, None, f"{phase}: no source (patient has no task data)")
                continue
            plan.warn(f"could not locate source .mat for phase {phase}")
            continue

        schedule_mkdir(dst_subdir, f"create {subdir}/ for canonical phase files")
        if src != dst:
            plan.add("rename", src, dst, f"{phase}: {src.name}")

    # 2. Implant xlsx: rename to implant/implant_pat_NN.xlsx.
    xlsx_src = _locate_implant_xlsx(patient_dir, n, nn, plan)
    xlsx_dst: Optional[Path] = None
    if xlsx_src is not None:
        xlsx_dst = implant_dir / f"implant_pat_{nn}.xlsx"
        schedule_mkdir(implant_dir, "create implant/")
        if xlsx_src != xlsx_dst:
            plan.add("rename", xlsx_src, xlsx_dst, f"xlsx: {xlsx_src.name}")

        # 3. Generate implant CSV from the xlsx (always done).
        csv_dst = patient_dir / f"implant_pat_{nn}.csv"
        plan.add("generate_csv", xlsx_dst, csv_dst,
                 "build implant_pat_NN.csv from xlsx")

        # Delete mixed-case predecessors at patient root (dedup candidates).
        legacy_csvs = {
            patient_dir / f"Implant_pat_{n}.csv",
            patient_dir / f"Implant_pat_{nn}.csv",
        }
        for old_csv in sorted(legacy_csvs):
            if old_csv.exists() and old_csv != csv_dst:
                plan.add("delete", old_csv, None, f"legacy implant csv: {old_csv.name}")
    else:
        plan.warn("no implant xlsx found; implant_pat_NN.csv will not be generated")

    # 4. channel_labels.csv consolidation — single source of truth at root.
    root_csv = patient_dir / "channel_labels.csv"
    dup_sources: List[Path] = []
    for c in CHANNEL_LABELS_DUP_CANDIDATES:
        p = patient_dir / c
        if p.exists() and p != root_csv:
            dup_sources.append(p)

    if root_csv.exists() and dup_sources:
        # Verify each duplicate matches the canonical root file.
        for dup in dup_sources:
            plan.add("verify_labels", dup, root_csv,
                     f"verify {dup.relative_to(patient_dir)} == channel_labels.csv")
    elif not root_csv.exists() and dup_sources:
        # Generate canonical from first duplicate; verify others against first.
        plan.add("generate_labels", dup_sources[0], root_csv,
                 f"channel_labels.csv from {dup_sources[0].relative_to(patient_dir)}")
        for extra in dup_sources[1:]:
            plan.add("verify_labels", extra, dup_sources[0],
                     f"verify {extra.relative_to(patient_dir)} == "
                     f"{dup_sources[0].relative_to(patient_dir)}")
    elif not root_csv.exists() and not dup_sources:
        plan.warn("no channel_labels source found")

    # Delete all duplicate copies.
    for p in dup_sources:
        plan.add("delete", p, None,
                 f"duplicate labels file: {p.relative_to(patient_dir)}")

    # 5. Delete legacy extras at patient root: channel_labels.mat and
    #    obsolete phase symlinks.
    extra_mat = patient_dir / "channel_labels.mat"
    if extra_mat.exists() or extra_mat.is_symlink():
        plan.add("delete", extra_mat, None, "duplicate channel_labels.mat at patient root")

    for name in ("rsPre.mat", "rsPost.mat", "taskLearn.mat", "taskTest.mat"):
        p = patient_dir / name
        if p.exists() or p.is_symlink():
            plan.add("delete", p, None, f"legacy phase symlink: {name}")

    # 6. Move vendor *_INFO.txt files next to their canonical .mat, renamed
    #    to the canonical phase stem (rest_pre.info.txt, task_learn.info.txt...).
    for vendor_rel, (subdir, target_name) in INFO_TXT_LAYOUT.items():
        src = patient_dir / vendor_rel
        if not src.exists():
            continue
        dst = patient_dir / subdir / target_name
        if dst.exists():
            # Already migrated; skip.
            continue
        schedule_mkdir(dst.parent, f"create {subdir}/ (for info.txt)")
        plan.add("rename", src, dst,
                 f"info: {src.relative_to(patient_dir)} -> {dst.relative_to(patient_dir)}")

    # 7. rmdir now-empty vendor subdirs.
    for vendor_name in OBSOLETE_VENDOR_SUBDIRS:
        vendor = patient_dir / vendor_name
        if vendor.is_dir():
            plan.add("rmdir", vendor, None, f"remove empty vendor subdir: {vendor_name}/")

    # 8. Provenance.
    plan.add("write_provenance", None, patient_dir / "provenance.md",
             "patient provenance record")

    return plan


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def render_plan(plan: MigrationPlan) -> str:
    """Return a human-readable report for one patient's plan."""
    lines = [f"## {plan.patient}", f"_{plan.patient_dir}_", ""]
    if plan.warnings:
        lines.append("**Warnings:**")
        for w in plan.warnings:
            lines.append(f"- ⚠ {w}")
        lines.append("")
    if not plan.actions:
        lines.append("_no actions planned_")
        return "\n".join(lines)

    lines.append("| # | kind | src | dst | notes |")
    lines.append("|---|------|-----|-----|-------|")
    for i, a in enumerate(plan.actions, start=1):
        src = a.src.relative_to(plan.patient_dir) if a.src else ""
        dst = a.dst.relative_to(plan.patient_dir) if a.dst else ""
        lines.append(f"| {i} | {a.kind} | `{src}` | `{dst}` | {a.notes} |")
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def build_provenance(plan: MigrationPlan) -> str:
    """Return markdown content for provenance.md.

    Records original vendor filenames, their canonical targets, and sha256
    of each migrated .mat. The hash is computed on the *destination* file,
    whose content is byte-identical to the pre-rename source (``rename``
    does not modify contents); this keeps the hash correct even when
    ``build_provenance`` runs after the rename has executed.

    For idempotent re-runs where no rename action fired (the file was
    already in canonical location), we still enumerate every canonical
    .mat target and record its current hash. The vendor filename column
    is filled from the plan's rename actions when present, or from a
    pre-existing ``provenance.md`` table when available, so information
    is never lost across re-runs.
    """
    nn = f"{_patnum(plan.patient):02d}"

    # 1. Build a {canonical_rel: vendor_rel} map from rename actions,
    #    including both .mat and .xlsx renames.
    vendor_map: dict[str, str] = {}
    for action in plan.actions:
        if action.kind != "rename" or action.src is None or action.dst is None:
            continue
        if action.dst.suffix not in (".mat", ".xlsx"):
            continue
        rel_dst = str(action.dst.relative_to(plan.patient_dir))
        rel_src = str(action.src.relative_to(plan.patient_dir))
        vendor_map[rel_dst] = rel_src

    # 2. Recover any previously-recorded vendor names so we don't lose them
    #    on idempotent re-runs.
    old_prov = plan.patient_dir / "provenance.md"
    if old_prov.exists():
        for line in old_prov.read_text().splitlines():
            # Rows look like: | `resting/rest_pre.mat` | `...` | `...` |
            if line.startswith("| `") and line.count("|") >= 4:
                parts = [c.strip().strip("`") for c in line.strip("|").split("|")]
                if len(parts) >= 2 and parts[0].endswith(".mat"):
                    if parts[0] not in vendor_map or vendor_map[parts[0]].startswith("resting/") \
                            or vendor_map[parts[0]].startswith("task/"):
                        if parts[1] and parts[1] not in ("vendor filename", "(unknown)"):
                            vendor_map.setdefault(parts[0], parts[1])

    # 3. Discover every canonical .mat + the implant xlsx currently on disk.
    canonical_files: list[Path] = []
    for subdir in ("resting", "task"):
        d = plan.patient_dir / subdir
        if d.is_dir():
            canonical_files.extend(sorted(d.glob("*.mat")))
    implant_dir = plan.patient_dir / "implant"
    if implant_dir.is_dir():
        canonical_files.extend(sorted(implant_dir.glob("*.xlsx")))

    lines = [
        f"# Provenance — {plan.patient}",
        "",
        f"- Import date: {date.today().isoformat()}",
        f"- Canonical id: `pat_{nn}`",
        "",
        "## Source → canonical mapping",
        "",
        "| canonical | vendor filename | sha256 |",
        "|-----------|-----------------|--------|",
    ]
    for f in canonical_files:
        rel_dst = str(f.relative_to(plan.patient_dir))
        vendor = vendor_map.get(rel_dst, "(unknown — pre-existing canonical file)")
        try:
            digest = _sha256(f)
        except OSError as exc:
            digest = f"(sha256-error: {exc})"
        lines.append(f"| `{rel_dst}` | `{vendor}` | `{digest}` |")

    lines += [
        "",
        "## Notes",
        "",
        "- The xlsx under `implant/` is the source of truth for epileptic",
        "  contact labels (red-font cells). Never overwrite it.",
        "- `implant_pat_NN.csv` is regenerated from the xlsx at normalization",
        "  time. Edit the xlsx, not the CSV.",
        "- Phase files live under `resting/` and `task/`. The loader maps",
        "  phase key → subdir automatically; no top-level symlinks exist.",
        "",
        "## Adding a new patient",
        "",
        "1. Drop the vendor directory (whatever its shape) under",
        "   `data/raw/stereoeeg_patients/Pat_NN/`.",
        "2. Run `lrg-eegfc data normalize --patients Pat_NN` to preview.",
        "3. Re-run with `--apply` once the plan looks right.",
        "4. This file is (re)generated by the normalizer.",
    ]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Labels helpers
# ---------------------------------------------------------------------------

def _read_label_lines(path: Path) -> List[str]:
    """Read a channel-labels source file (.txt or .csv) into a list of raw labels.

    Label contents are returned **verbatim** (reference suffix ``,G2`` and
    whitespace preserved). The reader only:
    - Skips an optional header row ``label`` (case-insensitive).
    - Strips surrounding double-quotes added by CSV quoting.
    - Rejoins a label that CSV parsing split on an interior comma such as
      ``"A 1,G2"`` (because the source file may or may not be properly
      quoted); the interior comma is put back exactly.

    Empty lines are dropped.
    """
    labels: List[str] = []
    with path.open(newline="") as f:
        reader = csv.reader(f)
        first = True
        for row in reader:
            if not row:
                continue
            # Rejoin any interior commas CSV split apart, preserving
            # whatever separator characters were in the source.
            cell = ",".join(c for c in row if c is not None).strip()
            if not cell:
                continue
            if first:
                first = False
                if cell.lower() == "label":
                    continue
            labels.append(cell.strip('"'))
    return labels


def _label_comparison_key(raw: str) -> str:
    """Return a comparison-only key for a label.

    Used **solely** when comparing two label sources that may differ in:
    - reference-suffix annotation (``A 1`` vs ``A 1,G2``),
    - whitespace (``A 1`` vs ``A1``),
    - character encoding of the prime marker on contralateral probes
      (``P' 12``, ``Pì 12``, ``P� 12`` — per the canonical parser in
      ``lrg_eegfc.utils.io.patient.parse_seeg_label``, ``ì`` is a
      known mis-encoding of ``'`` and ``U+FFFD`` appears when the ``ì``
      itself got decoded into the replacement glyph).

    This key exists for comparison only. It never touches files on disk.
    """
    # Delegate to the canonical parser so any future quirk fixes (e.g.,
    # new vendor encodings) live in one place.
    from .patient import parse_seeg_label

    s = raw.strip().strip('"')
    if "," in s:
        s = s.split(",", 1)[0]
    # Normalise the known U+FFFD replacement to the ì that parse_seeg_label
    # already maps onto an apostrophe. Any other U+FFFD is a genuine data
    # corruption that should fail the comparison loudly.
    s = s.replace("�", "ì")
    probe, contact = parse_seeg_label(s)
    if probe is not None and contact is not None:
        return f"{probe}{contact}"
    # Fallback: just collapse whitespace — better than nothing.
    return "".join(s.split())


def write_labels_csv(source_txt: Path, dst_csv: Path) -> None:
    labels = _read_label_lines(source_txt)
    dst_csv.parent.mkdir(parents=True, exist_ok=True)
    with dst_csv.open("w", newline="") as f:
        w = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
        w.writerow(["label"])
        for lab in labels:
            w.writerow([lab])


def verify_labels_equal(a: Path, b: Path) -> Tuple[bool, str]:
    """Compare two label files for contact-identity equality.

    Byte equality is preferred; if the files differ only by the optional
    ``,G2`` reference suffix or whitespace, the comparison still passes
    (contact sets match). In that case no data is altered on disk.

    Returns ``(ok, msg)`` where ``msg`` is empty on success, otherwise a
    short diff for logging.
    """
    la = _read_label_lines(a)
    lb = _read_label_lines(b)
    if la == lb:
        return True, ""
    # Fall back to contact-identity comparison (ignores annotation only).
    ka = [_label_comparison_key(x) for x in la]
    kb = [_label_comparison_key(x) for x in lb]
    if ka == kb:
        return True, ""
    diff_a = [x for k, x in zip(ka, la) if k not in set(kb)]
    diff_b = [x for k, x in zip(kb, lb) if k not in set(ka)]
    return False, (
        f"labels differ between {a.name} and {b.name}: "
        f"only-in-{a.name}={diff_a[:5]}{'...' if len(diff_a) > 5 else ''}, "
        f"only-in-{b.name}={diff_b[:5]}{'...' if len(diff_b) > 5 else ''}"
    )


# ---------------------------------------------------------------------------
# Xlsx -> CSV helper
# ---------------------------------------------------------------------------

def write_implant_csv_from_xlsx(xlsx_path: Path, dst_csv: Path) -> None:
    """Convert the first meaningful sheet of an implant xlsx into a CSV.

    Mirrors the existing Pat_05 / Pat_07 / Pat_08 CSV schema: values are
    written verbatim. All rows are padded to the sheet's max column count
    so pandas can parse the result as a fixed-width table (legacy CSVs
    use ``Unnamed: 5/6/...`` placeholder headers — we preserve that
    shape by padding with empty strings).
    """
    try:
        import openpyxl  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "openpyxl is required to convert implant xlsx to csv"
        ) from exc

    wb = openpyxl.load_workbook(xlsx_path, data_only=True)

    # Prefer the sheet used by load_epileptic_nodes.
    ws = wb.active
    for name in ("all_leads", "ALL_LEADS"):
        if name in wb.sheetnames:
            ws = wb[name]
            break

    # First pass: collect rows, stop at first fully-empty row, and
    # determine the max column count for padding.
    rows: list[list] = []
    max_cols = 0
    for raw_row in ws.iter_rows(values_only=True):
        row = list(raw_row)
        if all(v is None or (isinstance(v, str) and not v.strip()) for v in row):
            break
        # Trim trailing Nones only to find the actual content width of this row.
        last = len(row)
        while last > 0 and (row[last - 1] is None or row[last - 1] == ""):
            last -= 1
        row = row[:last]
        if not row:
            continue
        rows.append(row)
        if last > max_cols:
            max_cols = last

    # Pad header row with Unnamed: N placeholders (legacy convention) so
    # every row has exactly max_cols fields.
    if rows:
        header = rows[0]
        if len(header) < max_cols:
            header = header + [f"Unnamed: {i}" for i in range(len(header), max_cols)]
            rows[0] = header

    dst_csv.parent.mkdir(parents=True, exist_ok=True)
    with dst_csv.open("w", newline="") as f:
        writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
        for row in rows:
            padded = row + [""] * (max_cols - len(row))
            writer.writerow([("" if v is None else v) for v in padded])


# ---------------------------------------------------------------------------
# Apply
# ---------------------------------------------------------------------------

def apply_plan(plan: MigrationPlan) -> List[str]:
    """Execute all actions. Returns a list of result-line strings."""
    results: List[str] = []
    for a in plan.actions:
        try:
            if a.kind == "mkdir":
                assert a.dst is not None
                a.dst.mkdir(parents=True, exist_ok=True)
                results.append(f"mkdir {a.dst}")
            elif a.kind == "rename":
                assert a.src is not None and a.dst is not None
                a.dst.parent.mkdir(parents=True, exist_ok=True)
                a.src.rename(a.dst)
                results.append(f"rename {a.src} -> {a.dst}")
            elif a.kind == "generate_csv":
                assert a.src is not None and a.dst is not None
                write_implant_csv_from_xlsx(a.src, a.dst)
                results.append(f"generated {a.dst}")
            elif a.kind == "generate_labels":
                assert a.src is not None and a.dst is not None
                write_labels_csv(a.src, a.dst)
                results.append(f"generated {a.dst}")
            elif a.kind == "verify_labels":
                assert a.src is not None and a.dst is not None
                ok, msg = verify_labels_equal(a.src, a.dst)
                if not ok:
                    raise RuntimeError(msg)
                results.append(f"verified {a.src.name} == {a.dst.name}")
            elif a.kind == "delete":
                assert a.src is not None
                if a.src.is_symlink() or a.src.is_file():
                    a.src.unlink()
                elif a.src.is_dir():
                    a.src.rmdir()
                results.append(f"deleted {a.src}")
            elif a.kind == "rmdir":
                assert a.src is not None
                if a.src.is_dir():
                    try:
                        a.src.rmdir()
                        results.append(f"rmdir {a.src}")
                    except OSError as exc:
                        # Non-empty: leave it and warn.
                        results.append(f"skipped rmdir {a.src}: {exc}")
            elif a.kind == "write_provenance":
                assert a.dst is not None
                content = build_provenance(plan)
                a.dst.write_text(content)
                results.append(f"wrote {a.dst}")
            elif a.kind == "noop":
                results.append(f"noop: {a.notes}")
            else:
                raise RuntimeError(f"unknown action kind: {a.kind}")
        except Exception as exc:
            results.append(f"FAILED {a.kind} src={a.src} dst={a.dst}: {exc}")
            raise
    return results
