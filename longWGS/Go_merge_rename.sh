#!/usr/bin/env bash
# Go_merge_rename.sh
# Step 1: Merge per-barcode fastq.gz files into single files
# Step 2: Rename merged outputs from old names to new names using a map file

set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
MAP_FILE=""
DRYRUN=false

usage() {
    cat <<'EOF'
Usage:
  bash Go_merge_rename.sh -i INPUT_DIR -o OUTPUT_DIR -m MAP_FILE [-n]
  bash Go_merge_rename.sh -h

Description:
  Merge per-barcode FASTQ files into one FASTQ per input directory, then rename
  the merged outputs from old names to new names.

Options:
  -i INPUT_DIR   directory containing per-barcode subdirectories (e.g. fastq_pass/)
  -o OUTPUT_DIR  output directory for merged fastq.gz files
  -m MAP_FILE    rename table
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
  bash Go_merge_rename.sh -i fastq_pass -o merged_fastqs -m rename_map.tsv -n
  bash Go_merge_rename.sh -i fastq_pass -o merged_fastqs -m rename_map.tsv
EOF
    exit "${1:-1}"
}

while getopts ":i:o:m:nh" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        m) MAP_FILE="$OPTARG" ;;
        n) DRYRUN=true ;;
        h) usage 0 ;;
        *) usage ;;
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

    [[ -z "$old_name" || -z "$new_name" ]] && continue
    MAP["$old_name"]="$new_name"
done < "$MAP_FILE"
echo "Loaded ${#MAP[@]} rename mappings (old_name → new_name)"
echo ""

# ── Step 1: Merge ─────────────────────────────────────────────────────────────
echo "── Step 1: Merge ────────────────────────────────────"
$DRYRUN || mkdir -p "$OUTPUT_DIR"

for barcode in $(ls "$INPUT_DIR" | grep barcode | sort); do
    src_dir="$INPUT_DIR/$barcode"
    [ -d "$src_dir" ] || continue

    n_files=$(ls "$src_dir"/*.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
    if [ "$n_files" -eq 0 ]; then
        echo "  [SKIP] $barcode : no fastq.gz files"
        continue
    fi

    out="$OUTPUT_DIR/${barcode}.fastq.gz"
    if [ -f "$out" ]; then
        echo "  [SKIP-exists] $out"
        continue
    fi

    echo "  MERGE: $src_dir/*.fastq.gz (${n_files} files) → $out"
    $DRYRUN || cat "$src_dir"/*.fastq.gz > "$out"
done
echo ""

# ── Step 2: Rename ────────────────────────────────────────────────────────────
echo "── Step 2: Rename ───────────────────────────────────"
for old_name in "${!MAP[@]}"; do
    new_name="${MAP[$old_name]}"
    src="$OUTPUT_DIR/${old_name}.fastq.gz"
    dst="$OUTPUT_DIR/${new_name}.fastq.gz"

    if [ ! -f "$src" ]; then
        echo "  [SKIP-notfound] $src"
        continue
    fi
    if [ -f "$dst" ]; then
        echo "  [SKIP-exists]   $dst"
        continue
    fi

    echo "  RENAME: ${old_name}.fastq.gz  →  ${new_name}.fastq.gz"
    $DRYRUN || mv "$src" "$dst"
done
echo ""

if $DRYRUN; then
    echo "Dry-run complete. Re-run without -n to apply."
else
    echo "Done."
    echo ""
    echo "Merged & renamed files in: $OUTPUT_DIR"
    ls "$OUTPUT_DIR"/*.fastq.gz 2>/dev/null | sort | xargs -I{} basename {}
fi
