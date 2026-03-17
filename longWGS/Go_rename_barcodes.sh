#!/usr/bin/env bash
# rename_barcodes.sh
# Rename barcode-named outputs → sample names in longWGS result directory
#
# Usage:
#   bash rename_barcodes.sh <map_file> <longwgs_outdir> [--execute]
#
#   map_file     : two-column TSV (sample  barcode)  e.g. samples_barcodes_01-19.txt
#   longwgs_outdir : path to longWGS output dir (contains 1_QC, 2_quast, ...)
#   --execute    : actually rename (default: dry-run only)
#
# Example (dry-run):
#   bash rename_barcodes.sh samples_barcodes_01-19.txt 1_longWGS_out
# Example (real rename):
#   bash rename_barcodes.sh samples_barcodes_01-19.txt 1_longWGS_out --execute

set -euo pipefail

# ── args ─────────────────────────────────────────────────────────────────────
MAP_FILE="${1:-}"
OUTDIR="${2:-}"
EXECUTE=false
[[ "${3:-}" == "--execute" ]] && EXECUTE=true

if [[ -z "$MAP_FILE" || -z "$OUTDIR" ]]; then
    echo "Usage: bash $0 <map_file> <longwgs_outdir> [--execute]"
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
declare -A MAP   # MAP[barcode]=sample
while IFS=$'\t ' read -r sample barcode rest; do
    [[ -z "$sample" || -z "$barcode" ]] && continue
    [[ "$sample" == \#* ]] && continue
    MAP["$barcode"]="$sample"
done < "$MAP_FILE"

echo "Loaded ${#MAP[@]} barcode→sample mappings"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Per-subdirectory rename logic
# ─────────────────────────────────────────────────────────────────────────────

for barcode in "${!MAP[@]}"; do
    sample="${MAP[$barcode]}"
    echo "── $barcode  →  $sample ──────────────────────────"

    # 1_QC/  :  barcode01.clean.fastq.gz  (or .fastq.gz / .fq.gz)
    QC_DIR="$OUTDIR/1_QC"
    if [[ -d "$QC_DIR" ]]; then
        for ext in ".clean.fastq.gz" ".fastq.gz" ".fq.gz" ".clean.fq.gz"; do
            src="$QC_DIR/${barcode}${ext}"
            dst="$QC_DIR/${sample}${ext}"
            [[ -e "$src" ]] && do_rename "$src" "$dst"
        done
    fi

    # 2_quast/barcode01/
    D="$OUTDIR/2_quast"
    [[ -d "$D" ]] && do_rename "$D/$barcode" "$D/$sample"

    # 3_autocycler/barcode01/
    D="$OUTDIR/3_autocycler"
    [[ -d "$D" ]] && do_rename "$D/$barcode" "$D/$sample"

    # 4_medaka/  :  barcode01_final_assembly.fasta
    D="$OUTDIR/4_medaka"
    if [[ -d "$D" ]]; then
        for ext in "_final_assembly.fasta" "_final_assembly.fa"; do
            src="$D/${barcode}${ext}"
            dst="$D/${sample}${ext}"
            [[ -e "$src" ]] && do_rename "$src" "$dst"
        done
    fi

    # 5_checkm2/barcode01/
    D="$OUTDIR/5_checkm2"
    [[ -d "$D" ]] && do_rename "$D/$barcode" "$D/$sample"

    # 6_coverage/barcode01/  (contains barcode01.sorted.bam etc.)
    D="$OUTDIR/6_coverage"
    if [[ -d "$D/$barcode" ]]; then
        # rename files inside first
        for f in "$D/$barcode/${barcode}".*; do
            [[ -e "$f" ]] || continue
            fname="$(basename "$f")"
            newfname="${sample}${fname#${barcode}}"
            do_rename "$f" "$D/$barcode/$newfname"
        done
        # rename the directory itself
        do_rename "$D/$barcode" "$D/$sample"
    fi

    # 7_bakta/barcode01/
    D="$OUTDIR/7_bakta"
    if [[ -d "$D/$barcode" ]]; then
        # rename files inside bakta dir (prefix-named files)
        for f in "$D/$barcode/${barcode}".*; do
            [[ -e "$f" ]] || continue
            fname="$(basename "$f")"
            newfname="${sample}${fname#${barcode}}"
            do_rename "$f" "$D/$barcode/$newfname"
        done
        do_rename "$D/$barcode" "$D/$sample"
    fi

    # 8_Bandage_image/  :  barcode01_bandage.png, barcode01.png, barcode01*.* etc.
    D="$OUTDIR/8_Bandage_image"
    if [[ -d "$D" ]]; then
        for f in "$D/${barcode}"*; do
            [[ -e "$f" ]] || continue
            fname="$(basename "$f")"
            suffix="${fname#${barcode}}"
            do_rename "$f" "$D/${sample}${suffix}"
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
