#!/usr/bin/env bash
# sync_travel_bundle.sh — sync the lrgeegfc working data set to/from a travel drive.
#
# The bundle = data/ minus pre-rendered figures and superseded raw archives (~78 GB).
# It carries everything the paper needs offline:
#   - FC caches (imcoh) + LRG caches (imcoh_lrg, lrg, corr)
#   - the matched-strength surrogate nulls  <-- the expensive-to-recompute ones
#   - raw sEEG timeseries (data/raw, live patients only)
#   - preprint / reports / audit tables
# It DROPS: data/outputs/figures/ (regenerate from caches) and
#           data/raw/stereoeeg_patients/archive/ (superseded pre-revendor data).
#
# Usage:
#   sync_travel_bundle.sh push   <MOUNT> [REPO]   # desktop  -> drive
#   sync_travel_bundle.sh pull   <MOUNT> [REPO]   # drive    -> laptop repo
#   sync_travel_bundle.sh verify <MOUNT> [REPO]   # is every bundle file present locally?
# Add -n anywhere for a dry run (rsync --dry-run):
#   sync_travel_bundle.sh push -n /media/giulio/TRAVEL
#
# <MOUNT> = the mounted drive root, e.g. /media/giulio/TRAVEL
# [REPO]  = repo root; defaults to the desktop path, override on the laptop.

set -euo pipefail

REPO_DEFAULT="/home/giulio/Documents/research/neural_networks/lrgeegfc"
BUNDLE_NAME="lrgeegfc_bundle"

EXCLUDES=(
  --exclude 'outputs/figures/'
  --exclude 'raw/stereoeeg_patients/archive/'
  --exclude '.trash_*/'
  --exclude '__pycache__/'
  --exclude '*.pyc'
  --exclude '.DS_Store'
  --exclude '*.tmp'
  --exclude '*~'
)

die() { echo "ERROR: $*" >&2; exit 1; }

# ---- parse args ----------------------------------------------------------
MODE="${1:-}"; shift || true
[[ -n "$MODE" ]] || die "missing mode (push|pull|verify)"

DRY=""
MOUNT=""
REPO="$REPO_DEFAULT"
_pos=0
for arg in "$@"; do
  case "$arg" in
    -n|--dry-run) DRY="--dry-run" ;;
    *)
      if [[ $_pos -eq 0 ]]; then MOUNT="$arg"; else REPO="$arg"; fi
      _pos=$((_pos+1)) ;;
  esac
done

[[ -n "$MOUNT" ]] || die "missing <MOUNT> (e.g. /media/giulio/TRAVEL)"
command -v rsync >/dev/null || die "rsync not installed"
[[ -d "$MOUNT" ]] || die "mount point '$MOUNT' does not exist — is the drive plugged in and mounted?"

DEST="$MOUNT/$BUNDLE_NAME"
SRC_DATA="$REPO/data"

# exFAT-safe: the WD travel drive can't store unix perms/owner/group, so don't
# try to preserve them (would error on every file and, under `set -e`, abort the
# run). --modify-window=1 absorbs exFAT's 2s timestamp granularity on re-sync so
# unchanged files aren't needlessly re-copied.
RSYNC_FLAGS=(-rltD --no-perms --no-owner --no-group --modify-window=1 \
             --info=progress2 --human-readable --partial $DRY)

# ---- helpers -------------------------------------------------------------
count_and_bytes() {  # $1 = dir ; prints "<files> <bytes>" honoring EXCLUDES
  local dir="$1"
  local n b
  n=$(find "$dir" -type f \
        -not -path '*/outputs/figures/*' \
        -not -path '*/raw/stereoeeg_patients/archive/*' \
        -not -path '*/.trash_*/*' \
        -not -path '*/__pycache__/*' 2>/dev/null | wc -l)
  b=$(du -sb --exclude=outputs/figures \
        --exclude=raw/stereoeeg_patients/archive \
        --exclude='.trash_*' "$dir" 2>/dev/null | cut -f1)
  echo "$n $b"
}

# ---- modes ---------------------------------------------------------------
case "$MODE" in
  push)
    [[ -d "$SRC_DATA" ]] || die "source '$SRC_DATA' not found (wrong REPO?)"
    echo ">> PUSH  $SRC_DATA/  ->  $DEST/data/   ${DRY:+(dry run)}"
    mkdir -p "$DEST/data"
    rsync "${RSYNC_FLAGS[@]}" "${EXCLUDES[@]}" "$SRC_DATA/" "$DEST/data/"
    if [[ -z "$DRY" ]]; then
      read -r nfiles nbytes < <(count_and_bytes "$SRC_DATA")
      {
        echo "lrgeegfc travel bundle"
        echo "created:   $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "from host: $(hostname)"
        echo "source:    $SRC_DATA"
        echo "files:     $nfiles"
        echo "bytes:     $nbytes  ($(numfmt --to=iec "$nbytes" 2>/dev/null || echo "$nbytes"))"
        echo "excluded:  outputs/figures/, raw/stereoeeg_patients/archive/, caches junk"
      } > "$DEST/MANIFEST.txt"
      cp -f "$0" "$DEST/sync_travel_bundle.sh" 2>/dev/null || true
      cat > "$DEST/RESTORE.md" <<'EOF'
# Restore on the laptop

1. Clone the repo and install it in the `lapbrain` conda env:
       git clone <repo-url> lrgeegfc && cd lrgeegfc
       conda activate lapbrain      # or: pip install -e .
2. Plug in this drive; note its mount point (e.g. /media/<you>/TRAVEL).
3. Pull the data into the repo (data/ must sit at <repo>/data/):
       bash /media/<you>/TRAVEL/lrgeegfc_bundle/sync_travel_bundle.sh \
            pull /media/<you>/TRAVEL /path/to/lrgeegfc
4. Verify:
       bash .../sync_travel_bundle.sh verify /media/<you>/TRAVEL /path/to/lrgeegfc

Alternative to step 3 (put data anywhere and point the code at it):
       export LRGEEGFC_DATA_ROOT=/some/other/place/data

Regenerate figures from caches (they were intentionally not bundled):
   the FC/LRG/surrogate caches are all present, so figure scripts and
   `lrg-eegfc plot ...` run without recompute.

Sync work BACK after the trip: run `push` from the laptop, then `pull`
on the desktop (swap REPO for each machine's repo path).
EOF
      echo ">> wrote MANIFEST.txt, RESTORE.md, sync_travel_bundle.sh to $DEST/"
      echo ">> $nfiles files, $(numfmt --to=iec "$nbytes" 2>/dev/null || echo "$nbytes")"
    fi
    ;;

  pull)
    [[ -d "$DEST/data" ]] || die "no bundle at '$DEST/data' — wrong mount?"
    echo ">> PULL  $DEST/data/  ->  $SRC_DATA/   ${DRY:+(dry run)}"
    mkdir -p "$SRC_DATA"
    rsync "${RSYNC_FLAGS[@]}" "${EXCLUDES[@]}" "$DEST/data/" "$SRC_DATA/"
    echo ">> done. data root: $SRC_DATA"
    ;;

  verify)
    [[ -d "$DEST/data" ]] || die "no bundle at '$DEST/data'"
    [[ -d "$SRC_DATA" ]]  || die "no local data at '$SRC_DATA'"
    # Merge-safe check. A plain file/byte total over both trees is wrong here:
    # the local data/ legitimately holds files that were never in the bundle
    # (old-layout caches from a prior era, freshly rendered figures), so the
    # totals always DIFFER even when the copy is perfect. Instead ask the only
    # question that matters after a pull — "is every bundle file present and
    # identical locally?" — by running a dry-run rsync in the restore direction
    # and counting what would still transfer. Zero pending == the bundle fully
    # landed. Extra local files are ignored (no --delete anywhere in this tool).
    read -r bf bb < <(count_and_bytes "$DEST/data")
    printf "bundle (%s): %s files, %s bytes\n" "$DEST/data" "$bf" "$bb"
    local_pending() {  # dry-run bundle -> local; itemize; count real transfers
      # grep -c exits 1 on zero matches (the success case here), which would
      # trip `set -e`/pipefail — force a clean exit and let the count stand.
      rsync -rltD --no-perms --no-owner --no-group --modify-window=1 \
            --itemize-changes --dry-run "${EXCLUDES[@]}" \
            "$DEST/data/" "$SRC_DATA/" 2>/dev/null | grep -Ec '^(>f|cd)' || true
    }
    pending=$(local_pending)
    if [[ "$pending" -eq 0 ]]; then
      echo ">> MATCH — every bundle file is present and identical under $SRC_DATA"
      echo "   (extra local files outside the bundle are expected and ignored)"
    else
      echo ">> DIFFER — $pending bundle file(s) missing or out of date locally."
      echo "   Re-run pull, or list the deltas with:"
      echo "     rsync -rltD --no-perms --modify-window=1 -ni ${EXCLUDES[*]} \\"
      echo "       '$DEST/data/' '$SRC_DATA/'"
    fi
    ;;

  *) die "unknown mode '$MODE' (push|pull|verify)" ;;
esac
