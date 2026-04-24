#!/usr/bin/env bash
# migrate_data.sh — Create the new data/ directory layout.
#
# Phase 1 (default): creates symlinks from new paths to old locations.
# Phase 2 (--finalize): replaces symlinks with real directory moves.
#
# Idempotent: safe to run multiple times.

set -euo pipefail

DATA="$(cd "$(dirname "$0")/.." && pwd)/data"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

link_if_missing() {
    local target="$1" link_name="$2"
    # Resolve the target relative to the link's parent directory
    local link_dir
    link_dir="$(dirname "$link_name")"
    local abs_target="$link_dir/$target"
    if [ -L "$link_name" ]; then
        echo "  [skip] $link_name (symlink exists)"
    elif [ -d "$link_name" ]; then
        echo "  [skip] $link_name (real dir exists)"
    elif [ -d "$abs_target" ]; then
        ln -s "$target" "$link_name"
        echo "  [link] $link_name -> $target"
    else
        echo "  [miss] $abs_target does not exist, skipping"
    fi
}

finalize_link() {
    local link_name="$1"
    if [ -L "$link_name" ]; then
        local target
        target="$(readlink -f "$link_name")"
        rm "$link_name"
        mv "$target" "$link_name"
        echo "  [move] $target -> $link_name"
    elif [ -d "$link_name" ]; then
        echo "  [skip] $link_name (already real dir)"
    else
        echo "  [miss] $link_name does not exist"
    fi
}

# ---------------------------------------------------------------------------
# Phase 1: Symlinks
# ---------------------------------------------------------------------------

create_symlinks() {
    echo "=== Creating new data/ layout with symlinks ==="

    # Section 1: Raw data
    mkdir -p "$DATA/raw"
    link_if_missing "../stereoeeg_patients" "$DATA/raw/stereoeeg_patients"

    # Section 2: Cache
    mkdir -p "$DATA/cache"
    link_if_missing "../corr_cache"              "$DATA/cache/corr"
    link_if_missing "../msc_cache"               "$DATA/cache/msc"
    link_if_missing "../msc_cache_dev"            "$DATA/cache/msc_dev"
    link_if_missing "../lrg_cache"               "$DATA/cache/lrg"
    link_if_missing "../lrg_cache_crema"          "$DATA/cache/lrg_crema"
    link_if_missing "../cleaned_corr_cache"       "$DATA/cache/cleaned_corr"
    link_if_missing "../fc_fig_cache"             "$DATA/cache/fc_fig"
    link_if_missing "../metric_concordance"       "$DATA/cache/metric_concordance"
    link_if_missing "../surrogate_validation"     "$DATA/cache/surrogate_validation"

    # ImCoh (imaginary coherence — volume-conduction immune FC)
    link_if_missing "../imcoh_cache"             "$DATA/cache/imcoh"
    link_if_missing "../imcoh_lrg_cache"         "$DATA/cache/imcoh_lrg"

    # Experimental (failed debiasing, kept for reference)
    link_if_missing "../bipolar"                 "$DATA/cache/bipolar"
    link_if_missing "../rescaled"                "$DATA/cache/rescaled"

    # Section 3: Reports (by research topic)
    mkdir -p "$DATA/reports"
    link_if_missing "../wp0_metric_exploration"   "$DATA/reports/metric_exploration"
    link_if_missing "../wp_scalar_vi"             "$DATA/reports/scalar_vi"
    link_if_missing "../wp_spatial"               "$DATA/reports/spatial"

    # Merge wp1 dirs into clinical_application
    if [ ! -d "$DATA/reports/clinical_application" ] && [ ! -L "$DATA/reports/clinical_application" ]; then
        mkdir -p "$DATA/reports/clinical_application"
        # Copy contents from wp1 dirs (they're small — figures + text)
        for d in "$DATA/wp1_report_figures" "$DATA/wp1_fixes" "$DATA/wp1_viz_investigation"; do
            if [ -d "$d" ]; then
                cp -rn "$d"/* "$DATA/reports/clinical_application/" 2>/dev/null || true
                echo "  [merge] $d -> reports/clinical_application/"
            fi
        done
    else
        echo "  [skip] reports/clinical_application (exists)"
    fi

    # Merge coclassification dirs
    if [ ! -d "$DATA/reports/coclassification" ] && [ ! -L "$DATA/reports/coclassification" ]; then
        mkdir -p "$DATA/reports/coclassification"
        for d in "$DATA/wp_coclassification_investigation" "$DATA/wp_coclassification_diagnosis" "$DATA/wp_coclassification_v2"; do
            if [ -d "$d" ]; then
                cp -rn "$d"/* "$DATA/reports/coclassification/" 2>/dev/null || true
                echo "  [merge] $d -> reports/coclassification/"
            fi
        done
    else
        echo "  [skip] reports/coclassification (exists)"
    fi

    # ImCoh investigation reports
    link_if_missing "../imcoh_figures"            "$DATA/reports/imcoh"
    link_if_missing "../imcoh_unanimity"          "$DATA/reports/imcoh_unanimity"
    link_if_missing "../imcoh_vi"                 "$DATA/reports/imcoh_vi"

    # Section 4: Outputs
    mkdir -p "$DATA/outputs"
    link_if_missing "../figures" "$DATA/outputs/figures"
    link_if_missing "../tables"  "$DATA/outputs/tables"

    echo ""
    echo "=== Done. New layout created with symlinks. ==="
    echo "Old directories are untouched. Run with --finalize to replace symlinks with moves."
}

# ---------------------------------------------------------------------------
# Phase 2: Finalize (replace symlinks with real moves)
# ---------------------------------------------------------------------------

finalize() {
    echo "=== Finalizing: replacing symlinks with real directory moves ==="

    finalize_link "$DATA/raw/stereoeeg_patients"

    for sub in corr msc msc_dev lrg lrg_crema cleaned_corr fc_fig metric_concordance surrogate_validation imcoh imcoh_lrg bipolar rescaled; do
        finalize_link "$DATA/cache/$sub"
    done

    finalize_link "$DATA/reports/metric_exploration"
    finalize_link "$DATA/reports/scalar_vi"
    finalize_link "$DATA/reports/spatial"
    finalize_link "$DATA/reports/imcoh"
    finalize_link "$DATA/reports/imcoh_unanimity"
    finalize_link "$DATA/reports/imcoh_vi"

    finalize_link "$DATA/outputs/figures"
    finalize_link "$DATA/outputs/tables"

    # Clean up now-empty old directories
    echo ""
    echo "Checking for empty old directories to remove..."
    for old in stereoeeg_patients corr_cache msc_cache msc_cache_dev lrg_cache lrg_cache_crema \
               cleaned_corr_cache fc_fig_cache metric_concordance surrogate_validation \
               imcoh_cache imcoh_lrg_cache bipolar rescaled \
               wp0_metric_exploration wp_scalar_vi wp_spatial figures tables \
               imcoh_figures imcoh_unanimity imcoh_vi \
               wp1_report_figures wp1_fixes wp1_viz_investigation \
               wp_coclassification_investigation wp_coclassification_diagnosis wp_coclassification_v2; do
        if [ -d "$DATA/$old" ] && [ -z "$(ls -A "$DATA/$old" 2>/dev/null)" ]; then
            rmdir "$DATA/$old"
            echo "  [rmdir] $old (was empty)"
        fi
    done

    echo ""
    echo "=== Finalization complete. ==="
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

case "${1:-}" in
    --finalize)
        finalize
        ;;
    *)
        create_symlinks
        ;;
esac
