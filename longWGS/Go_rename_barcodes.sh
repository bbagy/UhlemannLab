#!/usr/bin/env bash
# Go_rename_barcodes.sh
# Rename longWGS outputs from old names to new names in a result directory.
#
# Usage:
#   bash Go_rename_barcodes.sh <map_file> <longwgs_outdir> [--execute]
#
#   map_file     : rename table with either:
#                    1) header columns (recommended), e.g. old_name  new_name
#                    2) two columns without header, interpreted by default as:
#                         old_name  new_name
#                  Legacy sample/barcode maps are still auto-detected.
#   longwgs_outdir : path to longWGS output dir (contains 1_QC, 2_quast, ...)
#   --execute    : actually rename (default: dry-run only)
#
# Example (dry-run):
#   bash Go_rename_barcodes.sh rename_map.tsv 1_longWGS_out
# Example (real rename):
#   bash Go_rename_barcodes.sh rename_map.tsv 1_longWGS_out --execute

set -euo pipefail

show_help() {
    cat <<'EOF'
Usage:
  bash Go_rename_barcodes.sh <map_file> <longwgs_outdir> [--execute]

Description:
  Rename longWGS output files/directories from old names to new names.
  Default mode is dry-run. Add --execute to actually rename files.

Map file formats:
  1. Recommended header format:
       old    new
       KP0011 KP0011_new
       KP0063 KP0063_new

  2. Barcode-compatible header format:
       barcode  new
       barcode01 KP0011
       barcode02 KP0063

  3. Headerless two-column format:
       KP0011 KP0011_new
       KP0063 KP0063_new
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
  bash Go_rename_barcodes.sh rename_map.tsv 1_longWGS_out
  bash Go_rename_barcodes.sh rename_map.tsv 1_longWGS_out --execute
EOF
}

# ── args ─────────────────────────────────────────────────────────────────────
MAP_FILE="${1:-}"
OUTDIR="${2:-}"
EXECUTE=false
[[ "${3:-}" == "--execute" ]] && EXECUTE=true

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
    exit 0
fi

if [[ -z "$MAP_FILE" || -z "$OUTDIR" ]]; then
    show_help
    exit 1
fi
if [[ ! -f "$MAP_FILE" ]]; then
    echo "ERROR: map file not found: $MAP_FILE"; exit 1
fi
if [[ ! -d "$OUTDIR" ]]; then
    echo "ERROR: output dir not found: $OUTDIR"; exit 1
fi

if $EXECUTE; then
    echo "=== EXECUTE MODE ==="
else
    echo "=== DRY-RUN MODE (pass --execute to apply) ==="
fi
echo ""

# ── rename helper ─────────────────────────────────────────────────────────────
do_rename() {
    local src="$1" dst="$2"
    if [[ ! -e "$src" ]]; then
        echo "  [SKIP-notfound] $src"
        return
    fi
    if [[ -e "$dst" ]]; then
        echo "  [SKIP-exists]   $dst already exists"
        return
    fi
    echo "  RENAME: $(basename "$src")  →  $(basename "$dst")"
    $EXECUTE && mv "$src" "$dst"
}

# ── read mapping ──────────────────────────────────────────────────────────────
declare -A MAP   # MAP[old_name]=new_name

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
        # Default no-header format: old_name new_name
        old_name="${c1:-}"
        new_name="${c2:-}"

        # Backward-compatible auto-detection for legacy sample/barcode tables:
        # if the second column looks like barcodeXX and the first does not,
        # interpret as new_name old_name (sample barcode).
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

# ─────────────────────────────────────────────────────────────────────────────
# Per-subdirectory rename logic
# ─────────────────────────────────────────────────────────────────────────────

for old_name in "${!MAP[@]}"; do
    new_name="${MAP[$old_name]}"
    echo "── $old_name  →  $new_name ──────────────────────────"

    # 1_QC/  :  old_name.clean.fastq.gz  (or .fastq.gz / .fq.gz)
    QC_DIR="$OUTDIR/1_QC"
    if [[ -d "$QC_DIR" ]]; then
        for ext in ".clean.fastq.gz" ".fastq.gz" ".fq.gz" ".clean.fq.gz"; do
            src="$QC_DIR/${old_name}${ext}"
            dst="$QC_DIR/${new_name}${ext}"
            [[ -e "$src" ]] && do_rename "$src" "$dst"
        done
    fi

    # 2_quast/old_name/
    D="$OUTDIR/2_quast"
    [[ -d "$D" ]] && do_rename "$D/$old_name" "$D/$new_name"

    # 3_autocycler/old_name/
    D="$OUTDIR/3_autocycler"
    [[ -d "$D" ]] && do_rename "$D/$old_name" "$D/$new_name"

    # 4_medaka/  :  old_name_final_assembly.fasta
    D="$OUTDIR/4_medaka"
    if [[ -d "$D" ]]; then
        for ext in "_final_assembly.fasta" "_final_assembly.fa"; do
            src="$D/${old_name}${ext}"
            dst="$D/${new_name}${ext}"
            [[ -e "$src" ]] && do_rename "$src" "$dst"
        done
    fi

    # 5_checkm2/old_name/
    D="$OUTDIR/5_checkm2"
    [[ -d "$D" ]] && do_rename "$D/$old_name" "$D/$new_name"

    # 6_coverage/old_name/  (contains old_name.sorted.bam etc.)
    D="$OUTDIR/6_coverage"
    if [[ -d "$D/$old_name" ]]; then
        # rename files inside first
        for f in "$D/$old_name/${old_name}".*; do
            [[ -e "$f" ]] || continue
            fname="$(basename "$f")"
            newfname="${new_name}${fname#${old_name}}"
            do_rename "$f" "$D/$old_name/$newfname"
        done
        # rename the directory itself
        do_rename "$D/$old_name" "$D/$new_name"
    fi

    # 7_bakta/old_name/
    D="$OUTDIR/7_bakta"
    if [[ -d "$D/$old_name" ]]; then
        # rename files inside bakta dir (prefix-named files)
        for f in "$D/$old_name/${old_name}".*; do
            [[ -e "$f" ]] || continue
            fname="$(basename "$f")"
            newfname="${new_name}${fname#${old_name}}"
            do_rename "$f" "$D/$old_name/$newfname"
        done
        do_rename "$D/$old_name" "$D/$new_name"
    fi

    # 8_Bandage_image/  :  old_name_bandage.png, old_name.png, old_name*.* etc.
    D="$OUTDIR/8_Bandage_image"
    if [[ -d "$D" ]]; then
        for f in "$D/${old_name}"*; do
            [[ -e "$f" ]] || continue
            fname="$(basename "$f")"
            suffix="${fname#${old_name}}"
            do_rename "$f" "$D/${new_name}${suffix}"
        done
    fi

    echo ""
done

# ─────────────────────────────────────────────────────────────────────────────
# checkm2_coverage_summary.xlsx  — warn user (content not auto-patched)
# ─────────────────────────────────────────────────────────────────────────────
XLSX="$OUTDIR/checkm2_coverage_summary.xlsx"
if [[ -f "$XLSX" ]]; then
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "NOTE: $XLSX exists."
    echo "  File/dir names will be renamed, but the XLSX"
    echo "  still contains barcode names internally."
    echo "  Re-run the pipeline summary rule, or edit manually."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
fi

echo ""
if $EXECUTE; then
    echo "Done. All renames applied."
else
    echo "Dry-run complete. Re-run with --execute to apply."
fi
