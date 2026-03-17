#!/usr/bin/env bash
# Go_merge_rename.sh
# Step 1: Merge per-barcode fastq.gz files into single files
# Step 2: Rename barcode names → sample names using a map file
#
# Usage:
#   bash Go_merge_rename.sh -i INPUT_DIR -o OUTPUT_DIR -m MAP_FILE [-n]
#
#   -i INPUT_DIR   : directory containing barcode subdirs (e.g. fastq_pass/)
#   -o OUTPUT_DIR  : output directory for merged fastq.gz files
#   -m MAP_FILE    : two-column TSV (sample  barcode)
#   -n             : dry-run (show what would be done, no actual changes)
#
# Map file format (tab or space separated):
#   KP0011  barcode01
#   KP0063  barcode02
#   ...
#
# Example (dry-run):
#   bash Go_merge_rename.sh -i fastq_pass -o merged_fastqs -m samples_barcodes.txt -n
# Example (run):
#   bash Go_merge_rename.sh -i fastq_pass -o merged_fastqs -m samples_barcodes.txt

set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
MAP_FILE=""
DRYRUN=false

usage() {
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR -m MAP_FILE [-n]"
    exit 1
}

while getopts "i:o:m:n" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        m) MAP_FILE="$OPTARG" ;;
        n) DRYRUN=true ;;
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

# ── load map: MAP[barcode]=sample ─────────────────────────────────────────────
declare -A MAP
while IFS=$'\t ' read -r sample barcode rest; do
    [[ -z "$sample" || -z "$barcode" ]] && continue
    [[ "$sample" == \#* ]] && continue
    MAP["$barcode"]="$sample"
done < "$MAP_FILE"
echo "Loaded ${#MAP[@]} barcode→sample mappings"
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
for barcode in "${!MAP[@]}"; do
    sample="${MAP[$barcode]}"
    src="$OUTPUT_DIR/${barcode}.fastq.gz"
    dst="$OUTPUT_DIR/${sample}.fastq.gz"

    if [ ! -f "$src" ]; then
        echo "  [SKIP-notfound] $src"
        continue
    fi
    if [ -f "$dst" ]; then
        echo "  [SKIP-exists]   $dst"
        continue
    fi

    echo "  RENAME: ${barcode}.fastq.gz  →  ${sample}.fastq.gz"
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
