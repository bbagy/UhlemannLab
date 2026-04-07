#!/usr/bin/env bash
set -euo pipefail

usage(){
  echo "Usage: $0 -i FASTQ_DIR -o OUTPUT_DIR -d KRAKEN2_DB [-s SNAKEDIR] [-c CORES] [-j JOBS] [-m IMAGE] [-n] [-K]"
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

FASTQ_DIR=""
OUTPUT_DIR=""
DB=""
SNAKEDIR=""
CORES=8
JOBS=4
IMAGE="kbracken:1.0"
DRYRUN=0
KEEP_GOING=0

while getopts "i:o:d:s:c:j:m:nK" opt; do
  case "$opt" in
    i) FASTQ_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    d) DB="$OPTARG" ;;
    s) SNAKEDIR="$OPTARG" ;;
    c) CORES="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    m) IMAGE="$OPTARG" ;;
    n) DRYRUN=1 ;;
    K) KEEP_GOING=1 ;;
    *) usage ;;
  esac
done

[ -z "$FASTQ_DIR" ] && usage
[ -z "$OUTPUT_DIR" ] && usage
[ -z "$DB" ] && usage

FASTQ_DIR_ABS="$(abs_path "$FASTQ_DIR")" || { echo "[KBracken] FASTQ_DIR not found: $FASTQ_DIR"; exit 1; }
DB_ABS="$(abs_path "$DB")" || { echo "[KBracken] DB not found: $DB"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
PIPELINE_DIR="$SCRIPT_DIR"
if [ -n "$SNAKEDIR" ]; then
  [ -d "$SNAKEDIR" ] || { echo "[KBracken] SNAKEDIR not found: $SNAKEDIR"; exit 1; }
  PIPELINE_DIR="$(cd "$SNAKEDIR" && pwd -P)"
fi

SNAKEFILE_NAME="Go_KBracken.smk"
if [ ! -f "$PIPELINE_DIR/$SNAKEFILE_NAME" ] && [ -f "$PIPELINE_DIR/Go_KBracken_V1.smk" ]; then
  SNAKEFILE_NAME="Go_KBracken_V1.smk"
fi
if [ ! -f "$PIPELINE_DIR/$SNAKEFILE_NAME" ]; then
  echo "[KBracken][FATAL] Snakefile not found: $PIPELINE_DIR/$SNAKEFILE_NAME"
  exit 1
fi

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
  echo "[KBracken][FATAL] Docker image not found locally: $IMAGE"
  echo "[KBracken] Build example:"
  echo "  cd \"$SCRIPT_DIR/docker/KBracken\" && docker build --network=host -t kbracken:1.0 ."
  exit 1
fi

WORKDIR="$(pwd)"

run(){
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$WORKDIR":/work \
    -v "$PIPELINE_DIR":/pipeline:ro \
    -v "$FASTQ_DIR_ABS":/fastq:ro \
    -v "$DB_ABS":/db/kraken2:ro \
    -w /work \
    "$IMAGE" \
    "$@"
}

BASE_ARGS=(
  snakemake
  --snakefile "/pipeline/$SNAKEFILE_NAME"
  --config
  fastq_dir=/fastq
  output_dir="$OUTPUT_DIR"
  db=/db/kraken2
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
run "${BASE_ARGS[@]}" 2>&1 | tee kbracken.log
rc=${PIPESTATUS[0]}
set -e

if [ "$rc" -ne 0 ] && grep -qiE "lock|unlock|LockException|cannot be locked" kbracken.log; then
  echo "[KBracken] Detected lock issue -> running --unlock then retry..."
  run "${BASE_ARGS[@]}" --unlock
  set +e
  run "${BASE_ARGS[@]}" 2>&1 | tee -a kbracken.log
  rc=${PIPESTATUS[0]}
  set -e
fi

exit "$rc"
