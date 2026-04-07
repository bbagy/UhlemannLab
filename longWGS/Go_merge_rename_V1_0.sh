#!/usr/bin/env bash
# Go_merge_rename_V1_0.sh
# Step 1: Merge per-barcode fastq.gz files into single files
# Step 2: Rename merged outputs from old names to new names using a map file

set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
MAP_FILE=""
DRYRUN=false
MERGE_UNTIL_SIZE=""
MERGE_UNTIL_BYTES=0

usage() {
    cat <<'EOF'
Usage:
  bash Go_merge_rename_V1_0.sh -i INPUT_DIR -o OUTPUT_DIR -m MAP_FILE [-n] [--merge-until-size 200M]
  bash Go_merge_rename_V1_0.sh -h

Description:
  Merge per-barcode FASTQ files into one FASTQ per input directory, then rename
  the merged outputs from old names to new names.
  If --merge-until-size is given, files are merged whole-file at a time until
  the cumulative compressed size reaches or exceeds the requested limit.
  Files are never cut mid-gzip file.

Options:
  -i INPUT_DIR   directory containing per-barcode subdirectories (e.g. fastq_pass/)
  -o OUTPUT_DIR  output directory for merged fastq.gz files
  -m MAP_FILE    rename table
  --merge-until-size SIZE
                 merge whole .fastq.gz files until SIZE is reached
                 examples: 200M, 500M, 1G
  -n             dry-run only
  -h             show this help

Map file formats:
  1. Recommended header format:
       old    new
       barcode01 KP0011
       barcode02 KP0063

  2. Barcode-compatible header format:
       barcode  new
       barcode01 KP0011
       barcode02 KP0063

  3. Headerless two-column format:
       barcode01 KP0011
       barcode02 KP0063
     This is interpreted by default as:
       old_name new_name

  4. Legacy sample/barcode format is still supported:
       KP0011 barcode01
       KP0063 barcode02

Recognized header names:
  old-name column:
    old, old_name, current, current_name, source, from, barcode
  new-name column:
    new, new_name, sample, sample_name, target, to

Examples:
  bash Go_merge_rename_V1_0.sh -i fastq_pass -o merged_fastqs -m rename_map.tsv -n
  bash Go_merge_rename_V1_0.sh -i fastq_pass -o merged_fastqs -m rename_map.tsv
  bash Go_merge_rename_V1_0.sh -i fastq_pass -o merged_fastqs -m rename_map.tsv --merge-until-size 500M
EOF
    exit "${1:-1}"
}

parse_size_to_bytes() {
    local raw="${1:-}"
    local num suffix
    if [[ ! "$raw" =~ ^([0-9]+)([KMG]?)$ ]]; then
        echo "ERROR: invalid size '$raw'. Use values like 200M, 500M, 1G, 800K." >&2
        exit 1
    fi
    num="${BASH_REMATCH[1]}"
    suffix="${BASH_REMATCH[2]}"
    case "$suffix" in
        K) echo $((num * 1024)) ;;
        M) echo $((num * 1024 * 1024)) ;;
        G) echo $((num * 1024 * 1024 * 1024)) ;;
        "") echo "$num" ;;
        *)
            echo "ERROR: unsupported size suffix in '$raw'." >&2
            exit 1
            ;;
    esac
}

file_size_bytes() {
    local f="$1"
    if stat -c %s "$f" >/dev/null 2>&1; then
        stat -c %s "$f"
    else
        stat -f %z "$f"
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i)
            INPUT_DIR="${2:-}"
            shift 2
            ;;
        -o)
            OUTPUT_DIR="${2:-}"
            shift 2
            ;;
        -m)
            MAP_FILE="${2:-}"
            shift 2
            ;;
        -n)
            DRYRUN=true
            shift
            ;;
        --merge-until-size)
            [[ $# -ge 2 ]] || { echo "ERROR: --merge-until-size requires a value." >&2; exit 1; }
            MERGE_UNTIL_SIZE="$2"
            MERGE_UNTIL_BYTES="$(parse_size_to_bytes "$2")"
            shift 2
            ;;
        -h|--help)
            usage 0
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            usage 1
            ;;
    esac
done

[ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$MAP_FILE" ] && usage

[ -d "$INPUT_DIR" ] || { echo "ERROR: INPUT_DIR not found: $INPUT_DIR"; exit 1; }
[ -f "$MAP_FILE"  ] || { echo "ERROR: MAP_FILE not found: $MAP_FILE";   exit 1; }

if $DRYRUN; then
    echo "=== DRY-RUN MODE (pass without -n to apply) ==="
else
    echo "=== EXECUTE MODE ==="
fi
if [[ -n "$MERGE_UNTIL_SIZE" ]]; then
    echo "=== MERGE LIMIT: $MERGE_UNTIL_SIZE (${MERGE_UNTIL_BYTES} bytes) ==="
fi
echo ""

# ── load map: MAP[old_name]=new_name ─────────────────────────────────────────
declare -A MAP

normalize_col() {
    echo "$1" | tr '[:upper:]' '[:lower:]' | tr -d ' _-'
}

is_barcode_name() {
    [[ "$1" =~ ^barcode[0-9]+$ ]]
}

header_old_idx=-1
header_new_idx=-1
line_no=0

while IFS=$'\t ' read -r c1 c2 c3 rest; do
    ((line_no += 1))
    [[ -z "${c1:-}" && -z "${c2:-}" ]] && continue
    [[ "${c1:-}" == \#* ]] && continue

    if [[ $line_no -eq 1 ]]; then
        cols=("${c1:-}" "${c2:-}" "${c3:-}")
        for idx in 0 1 2; do
            col="${cols[$idx]}"
            norm="$(normalize_col "$col")"
            case "$norm" in
                old|oldname|current|currentname|source|from|barcode)
                    header_old_idx=$idx
                    ;;
                new|newname|sample|samplename|target|to)
                    header_new_idx=$idx
                    ;;
            esac
        done

        if [[ $header_old_idx -ge 0 && $header_new_idx -ge 0 ]]; then
            continue
        fi
    fi

    vals=("${c1:-}" "${c2:-}" "${c3:-}")
    old_name=""
    new_name=""

    if [[ $header_old_idx -ge 0 && $header_new_idx -ge 0 ]]; then
        old_name="${vals[$header_old_idx]}"
        new_name="${vals[$header_new_idx]}"
    else
        old_name="${c1:-}"
        new_name="${c2:-}"

        if [[ -n "${c1:-}" && -n "${c2:-}" ]] && is_barcode_name "$c2" && ! is_barcode_name "$c1"; then
            old_name="$c2"
            new_name="$c1"
        fi
    fi

    old_name="${old_name//$'\r'/}"
    new_name="${new_name//$'\r'/}"

    [[ -z "$old_name" || -z "$new_name" ]] && continue
    MAP["$old_name"]="$new_name"
done < "$MAP_FILE"
echo "Loaded ${#MAP[@]} rename mappings (old_name → new_name)"
echo ""

# ── Step 1: Merge + rename immediately ───────────────────────────────────────
echo "── Step 1: Merge + rename ───────────────────────────"
$DRYRUN || mkdir -p "$OUTPUT_DIR"

for barcode in $(ls "$INPUT_DIR" | grep barcode | sort); do
    src_dir="$INPUT_DIR/$barcode"
    [ -d "$src_dir" ] || continue

    mapfile -t fastq_files < <(find "$src_dir" -maxdepth 1 -type f -name "*.fastq.gz" | sort)
    n_files="${#fastq_files[@]}"
    if [ "$n_files" -eq 0 ]; then
        echo "  [SKIP] $barcode : no fastq.gz files"
        continue
    fi

    target_name="${MAP[$barcode]:-$barcode}"
    out="$OUTPUT_DIR/${target_name}.fastq.gz"
    if [ -f "$out" ]; then
        echo "  [SKIP-exists] $out"
        continue
    fi

    selected_files=()
    selected_count=0
    cumulative_bytes=0

    for f in "${fastq_files[@]}"; do
        selected_files+=("$f")
        selected_count=$((selected_count + 1))
        cumulative_bytes=$((cumulative_bytes + $(file_size_bytes "$f")))
        if [[ "$MERGE_UNTIL_BYTES" -gt 0 && "$cumulative_bytes" -ge "$MERGE_UNTIL_BYTES" ]]; then
            break
        fi
    done

    if [[ "$selected_count" -eq "$n_files" ]]; then
        echo "  MERGE: $barcode -> ${target_name}.fastq.gz (${selected_count}/${n_files} files, max possible; bytes=$cumulative_bytes)"
    else
        echo "  MERGE: $barcode -> ${target_name}.fastq.gz (${selected_count}/${n_files} files, bytes=$cumulative_bytes, limit=$MERGE_UNTIL_BYTES)"
    fi

    if ! $DRYRUN; then
        cat "${selected_files[@]}" > "$out"
    fi
done
echo ""

if $DRYRUN; then
    echo "Dry-run complete. Re-run without -n to apply."
else
    echo "Done."
    echo ""
    echo "Merged files in final names: $OUTPUT_DIR"
    ls "$OUTPUT_DIR"/*.fastq.gz 2>/dev/null | sort | xargs -I{} basename {}
fi
