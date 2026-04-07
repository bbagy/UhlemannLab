#!/usr/bin/env bash
set -euo pipefail

usage(){
  echo "Usage: $0 -i MAGS_DIR -o OUTPUT_DIR -g GTDBTK_DATA_DIR -e EGGNOG_DATA_DIR -k KOFAM_SCAN_DIR [-s MAGS_PIPELINE_DIR] [-c CORES] [-j JOBS] [-m IMAGE] [-n] [-K]"
  exit 1
}

abs_path(){
  local p="$1"
  if [ -d "$p" ]; then
    (cd "$p" && pwd -P)
  elif [ -e "$p" ]; then
    (cd "$(dirname "$p")" && printf "%s/%s\n" "$(pwd -P)" "$(basename "$p")")
  else
    return 1
  fi
}

MAGS_DIR=""
OUTPUT_DIR=""
GTDBTK_DATA_DIR=""
EGGNOG_DATA_DIR=""
KOFAM_SCAN_DIR=""
SNAKEDIR=""
CORES=8
JOBS=4
IMAGE="mags-annotation:1.0"
DRYRUN=0
KEEP_GOING=0

while getopts "i:o:g:e:k:s:c:j:m:nK" opt; do
  case "$opt" in
    i) MAGS_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    g) GTDBTK_DATA_DIR="$OPTARG" ;;
    e) EGGNOG_DATA_DIR="$OPTARG" ;;
    k) KOFAM_SCAN_DIR="$OPTARG" ;;
    s) SNAKEDIR="$OPTARG" ;;
    c) CORES="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    m) IMAGE="$OPTARG" ;;
    n) DRYRUN=1 ;;
    K) KEEP_GOING=1 ;;
    *) usage ;;
  esac
done

[ -z "$MAGS_DIR" ] && usage
[ -z "$OUTPUT_DIR" ] && usage
[ -z "$GTDBTK_DATA_DIR" ] && usage
[ -z "$EGGNOG_DATA_DIR" ] && usage
[ -z "$KOFAM_SCAN_DIR" ] && usage

MAGS_DIR_ABS="$(abs_path "$MAGS_DIR")" || { echo "[Go_MAGs_Annotation] MAGS_DIR not found: $MAGS_DIR"; exit 1; }
GTDBTK_DATA_DIR_ABS="$(abs_path "$GTDBTK_DATA_DIR")" || { echo "[Go_MAGs_Annotation] GTDBTK_DATA_DIR not found: $GTDBTK_DATA_DIR"; exit 1; }
EGGNOG_DATA_DIR_ABS="$(abs_path "$EGGNOG_DATA_DIR")" || { echo "[Go_MAGs_Annotation] EGGNOG_DATA_DIR not found: $EGGNOG_DATA_DIR"; exit 1; }
KOFAM_SCAN_DIR_ABS="$(abs_path "$KOFAM_SCAN_DIR")" || { echo "[Go_MAGs_Annotation] KOFAM_SCAN_DIR not found: $KOFAM_SCAN_DIR"; exit 1; }

if ! compgen -G "$MAGS_DIR_ABS/*.fa" >/dev/null && ! compgen -G "$MAGS_DIR_ABS/*.fasta" >/dev/null; then
  echo "[Go_MAGs_Annotation][FATAL] No *.fa or *.fasta MAG files found in: $MAGS_DIR_ABS"
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
PIPELINE_DIR="$SCRIPT_DIR"
if [ -n "$SNAKEDIR" ]; then
  [ -d "$SNAKEDIR" ] || { echo "[Go_MAGs_Annotation] SNAKEDIR not found: $SNAKEDIR"; exit 1; }
  PIPELINE_DIR="$(cd "$SNAKEDIR" && pwd -P)"
fi

SNAKEFILE_NAME="Go_MAGs_Annotation.smk"
if [ ! -f "$PIPELINE_DIR/$SNAKEFILE_NAME" ] && [ -f "$PIPELINE_DIR/workflow/Go_MAGs_Annotation_V1.smk" ]; then
  SNAKEFILE_NAME="workflow/Go_MAGs_Annotation_V1.smk"
fi
if [ ! -f "$PIPELINE_DIR/$SNAKEFILE_NAME" ]; then
  echo "[Go_MAGs_Annotation][FATAL] Snakefile not found: $PIPELINE_DIR/$SNAKEFILE_NAME"
  exit 1
fi

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
  echo "[Go_MAGs_Annotation][FATAL] Docker image not found locally: $IMAGE"
  echo "[Go_MAGs_Annotation] Build example:"
  echo "  cd \"$SCRIPT_DIR/docker/MAGs\" && docker build -f Dockerfile.annotation -t mags-annotation:1.0 ."
  exit 1
fi

WORKDIR="$(pwd)"

run(){
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$WORKDIR":/work \
    -v "$PIPELINE_DIR":/pipeline:ro \
    -v "$MAGS_DIR_ABS":/mags:ro \
    -v "$GTDBTK_DATA_DIR_ABS":/db/gtdbtk:ro \
    -v "$EGGNOG_DATA_DIR_ABS":/db/eggnog:ro \
    -v "$KOFAM_SCAN_DIR_ABS":/db/kofam_scan:ro \
    -w /work \
    "$IMAGE" \
    "$@"
}

BASE_ARGS=(
  snakemake
  --snakefile "/pipeline/$SNAKEFILE_NAME"
  --config
  mags_dir=/mags
  output_dir="$OUTPUT_DIR"
  gtdbtk_data_dir=/db/gtdbtk
  gtdbtk_mash_db=/db/gtdbtk/mash_db/genomic_mash_db.msh
  eggnog_data_dir=/db/eggnog
  kofam_scan_dir=/db/kofam_scan
  kofam_exec=exec_annotation
  gtdbtk_threads="$CORES"
  prokka_threads="$CORES"
  eggnog_threads="$CORES"
  kofamscan_threads="$CORES"
  --cores "$CORES"
  --jobs "$JOBS"
  --latency-wait 60
  --rerun-incomplete
)

if [ "$DRYRUN" -eq 1 ]; then
  BASE_ARGS+=(--dry-run)
fi
if [ "$KEEP_GOING" -eq 1 ]; then
  BASE_ARGS+=(--keep-going)
fi

set +e
run "${BASE_ARGS[@]}" 2>&1 | tee mags_annotation.log
rc=${PIPESTATUS[0]}
set -e

if [ "$rc" -ne 0 ] && grep -qiE "lock|unlock|LockException|cannot be locked" mags_annotation.log; then
  echo "[Go_MAGs_Annotation] Detected lock issue -> running --unlock then retry..."
  run "${BASE_ARGS[@]}" --unlock
  set +e
  run "${BASE_ARGS[@]}" 2>&1 | tee -a mags_annotation.log
  rc=${PIPESTATUS[0]}
  set -e
fi

exit "$rc"
