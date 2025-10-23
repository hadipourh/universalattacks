#!/usr/bin/env bash

# Orthros Correlation Analysis Visualization Script
# Automatically runs Python visualization after CUDA analysis completes

set -eo pipefail

SCRIPT_NAME=$(basename "$0")
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
DEFAULT_PATTERNS=(
    "present_*r_correlations_*.csv"
    "orthros_*r_*_dlct.csv"
    "orthros_*r_keys_*.csv"
)
DEFAULT_THEME="academic"
PYTHON_BIN=${PYTHON_BIN:-python3}
PIP_BIN=${PIP_BIN:-pip3}
VENV_PATH="${VIRTUAL_ENV:-}"
REQUIREMENTS_FILE="$SCRIPT_DIR/requirements.txt"
VISUALIZER_SCRIPT="$SCRIPT_DIR/src/visualizer.py"

show_header() {
    cat <<'EOF'
Orthros Correlation Visualization System
=======================================
EOF
}

print_usage() {
    cat <<EOF
Usage: $SCRIPT_NAME [options]

Options:
  -i, --input PATTERN       Glob pattern or CSV file to visualize (can repeat)
  -o, --output-dir DIR      Directory to write generated figures (default: alongside CSV)
  -t, --theme NAME          Visualization theme to pass to visualizer (default: $DEFAULT_THEME)
  -r, --requirements PATH   Custom requirements.txt to install dependencies
  --python PATH             Python interpreter to use (default: $PYTHON_BIN)
  --pip PATH                Pip executable to use (default: $PIP_BIN)
  --dry-run                 Show what would run without executing visualizer
  --no-install              Skip automatic dependency installation
  --list                    List matched input files and exit
  -h, --help                Show this message and exit

Environment variables:
  VIRTUAL_ENV               If set, dependencies install into the active virtualenv
  PYTHON_BIN, PIP_BIN       Override interpreter/pip discovery without flags
EOF
}

log() {
    printf '%s %s\n' "[$(date '+%H:%M:%S')]" "$*"
}

warn() {
    printf '%s %s\n' "[$(date '+%H:%M:%S')]" "WARNING: $*" >&2
}

err() {
    printf '%s %s\n' "[$(date '+%H:%M:%S')]" "ERROR: $*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || err "Required command '$1' not found in PATH"
}

install_dependencies() {
    local req_file=${1:-$REQUIREMENTS_FILE}
    if [ ! -f "$req_file" ]; then
        warn "requirements file '$req_file' not found; skipping installation"
        return
    fi
    log "Installing Python dependencies from $req_file"
    "$PIP_BIN" install -r "$req_file"
}

resolve_inputs() {
    local patterns=()
    local files=()
    patterns+=("$@")
    if [ ${#patterns[@]} -eq 0 ]; then
        patterns=("${DEFAULT_PATTERNS[@]}")
    fi
    for pattern in "${patterns[@]}"; do
        for file in $pattern; do
            [ -f "$file" ] && files+=("$file")
        done
    done
    printf '%s\n' "${files[@]}"
}

ensure_visualizer() {
    [ -f "$VISUALIZER_SCRIPT" ] || err "Visualizer script not found at $VISUALIZER_SCRIPT"
}

generate_output_path() {
    local input_file=$1
    local output_dir=${2:-}
    local stem
    stem=$(basename "$input_file" .csv)
    if [ -n "$output_dir" ]; then
        mkdir -p "$output_dir"
        printf '%s/%s_analysis.png\n' "$output_dir" "$stem"
    else
        printf '%s_analysis.png\n' "$stem"
    fi
}

run_visualizer() {
    local csv=$1
    local output=$2
    local theme=$3
    local dry_run=$4
    if [ "$dry_run" = "true" ]; then
        log "[dry-run] $PYTHON_BIN $VISUALIZER_SCRIPT \"$csv\" --output \"$output\" --theme $theme --export-data"
        return 0
    fi
    "$PYTHON_BIN" "$VISUALIZER_SCRIPT" "$csv" --output "$output" --theme "$theme" --export-data
}

main() {
    local inputs=()
    local output_dir=""
    local theme="$DEFAULT_THEME"
    local custom_requirements=""
    local install_deps="true"
    local dry_run="false"
    local list_only="false"

    while [ $# -gt 0 ]; do
        case "$1" in
            -i|--input)
                shift
                [ $# -gt 0 ] || err "--input requires a pattern argument"
                inputs+=("$1")
                ;;
            -o|--output-dir)
                shift
                [ $# -gt 0 ] || err "--output-dir requires a directory argument"
                output_dir="$1"
                ;;
            -t|--theme)
                shift
                [ $# -gt 0 ] || err "--theme requires a name argument"
                theme="$1"
                ;;
            -r|--requirements)
                shift
                [ $# -gt 0 ] || err "--requirements requires a path argument"
                custom_requirements="$1"
                ;;
            --python)
                shift
                [ $# -gt 0 ] || err "--python requires a path argument"
                PYTHON_BIN="$1"
                ;;
            --pip)
                shift
                [ $# -gt 0 ] || err "--pip requires a path argument"
                PIP_BIN="$1"
                ;;
            --dry-run)
                dry_run="true"
                ;;
            --no-install)
                install_deps="false"
                ;;
            --list)
                list_only="true"
                ;;
            -h|--help)
                print_usage
                exit 0
                ;;
            --)
                shift
                break
                ;;
            *)
                err "Unknown argument: $1"
                ;;
        esac
        shift
    done

    show_header
    require_command "$PYTHON_BIN"

    if [ "$install_deps" = "true" ]; then
        require_command "$PIP_BIN"
        if [ -n "$custom_requirements" ]; then
            install_dependencies "$custom_requirements"
        else
            install_dependencies "$REQUIREMENTS_FILE"
        fi
    else
        log "Skipping dependency installation (--no-install)"
    fi

    ensure_visualizer

    # Bash/zsh compatible array handling
    matched_files=()
    while IFS= read -r line; do
        matched_files+=("$line")
    done < <(resolve_inputs "${inputs[@]}")

    if [ ${#matched_files[@]} -eq 0 ]; then
        warn "No correlation data files found. Searched patterns:"
        for pattern in "${inputs[@]:-${DEFAULT_PATTERNS[@]}}"; do
            echo "  - $pattern"
        done
        err "Run the CUDA/OpenMP analysis first to generate data files"
    fi

    log "Found ${#matched_files[@]} correlation data file(s)"

    if [ "$list_only" = "true" ]; then
        printf '%s\n' "${matched_files[@]}"
        exit 0
    fi

    local success_count=0
    local failure_count=0
    local outputs_generated=()

    for csv_file in "${matched_files[@]}"; do
        log "Processing: $csv_file"

        output_path=$(generate_output_path "$csv_file" "$output_dir")
        if run_visualizer "$csv_file" "$output_path" "$theme" "$dry_run"; then
            log "✓ Visualization created: $output_path"
            outputs_generated+=("$output_path")
            success_count=$((success_count + 1))
        else
            warn "✗ Failed to create visualization for $csv_file"
            failure_count=$((failure_count + 1))
        fi
    done

    echo ""
    log "Visualization processing completed!"
    log "Successful: $success_count | Failed: $failure_count"

    if [ ${#outputs_generated[@]} -gt 0 ]; then
        echo ""
        log "Generated files:"
        printf '  %s\n' "${outputs_generated[@]}"
    fi

    if [ "$dry_run" = "true" ]; then
        echo ""
        log "Dry-run mode: no files were created"
    else
        echo ""
        log "To view results:"
        echo "  - Open .pdf files for publication-quality figures"
        echo "  - Open .png files for presentations"
        echo "  - Read _report.txt files for detailed statistical analysis"
    fi
}

main "$@"