#!/usr/bin/env bash
set -euo pipefail
usage(){ echo "Usage: $0 -i INPUT -o OUTPUT -d DB -s SNAKEDIR [-p 0|1]"; exit 1; }

INPUT=""; OUTPUT=""; DB=""; SNAKEDIR=""; PORECHOP=0
while getopts "i:o:d:s:p:" opt; do
  case $opt in
    i) INPUT="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    d) DB="$OPTARG" ;;
    s) SNAKEDIR="$OPTARG" ;;
    p) PORECHOP="$OPTARG" ;;
    *) usage ;;
  esac
done
[ -z "$INPUT" ] && usage; [ -z "$OUTPUT" ] && usage; [ -z "$DB" ] && usage; [ -z "$SNAKEDIR" ] && usage

WORKDIR="$(pwd)"

CMD="snakemake --snakefile /pipeline/Go_longWGS.smk \
  --config Fastq_DIRS=\"$INPUT\" output=\"$OUTPUT\" database=\"/db\" threads=4 do_porechop=$PORECHOP \
  --cores 8 --rerun-incomplete"

run(){ docker run --rm -it \
  -u "$(id -u):$(id -g)" \
  -v "$WORKDIR":/work \
  -v "$DB":/db \
  -v "$SNAKEDIR":/pipeline \
  -w /work \
  longwgs bash -lc "$1"; }

set +e
run "$CMD" 2>&1 | tee longwgs.log
rc=${PIPESTATUS[0]}
set -e

if [ $rc -ne 0 ] && grep -qiE "lock|unlock|LockException|cannot be locked" longwgs.log; then
  echo "[Go_longWGS] Detected lock issue -> running --unlock then retry..."
  run "$CMD --unlock"
  
  set +e
  run "$CMD" 2>&1 | tee -a longwgs.log
  rc=${PIPESTATUS[0]}
  set -e
fi

# -----------------------------
# Bandage post-processing (only if success)
# -----------------------------
if [ $rc -eq 0 ]; then
  echo "[Go_longWGS] Pipeline finished -> running Bandage images..."

  BANDAGE_DIR="$OUTPUT/8_Bandage_image"
  mkdir -p "$BANDAGE_DIR"

  for gfa in "$OUTPUT"/3_autocycler/*/autocycler_out/consensus_assembly.gfa; do
    [ -f "$gfa" ] || continue

    sample=$(basename "$(dirname "$(dirname "$gfa")")")
    png="$BANDAGE_DIR/${sample}.consensus_assembly.png"

    echo "  Bandage: $gfa -> $png"
    run "Bandage image \"$gfa\" \"$png\""
  done
fi

exit $rc
