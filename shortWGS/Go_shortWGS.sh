#!/usr/bin/env bash
set -euo pipefail

usage(){
  echo "Usage: $0 -i FASTQ_DIR -o OUTPUT_DIR -d WGS_DB_DIR -k KRAKEN_DB_DIR -r GOWGS_DIR [-c CORES] [-m IMAGE] [-n]"
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
WGS_DB_DIR=""
KRAKEN_DB_DIR=""
GOWGS_DIR=""
CORES=8
IMAGE="shortwgs:1.0"
DRYRUN=0

while getopts "i:o:d:k:r:c:m:n" opt; do
  case "$opt" in
    i) FASTQ_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    d) WGS_DB_DIR="$OPTARG" ;;
    k) KRAKEN_DB_DIR="$OPTARG" ;;
    r) GOWGS_DIR="$OPTARG" ;;
    c) CORES="$OPTARG" ;;
    m) IMAGE="$OPTARG" ;;
    n) DRYRUN=1 ;;
    *) usage ;;
  esac
done

[ -z "$FASTQ_DIR" ] && usage
[ -z "$OUTPUT_DIR" ] && usage
[ -z "$WGS_DB_DIR" ] && usage
[ -z "$KRAKEN_DB_DIR" ] && usage
[ -z "$GOWGS_DIR" ] && usage

FASTQ_DIR_ABS="$(abs_path "$FASTQ_DIR")" || { echo "[Go_shortWGS] FASTQ_DIR not found: $FASTQ_DIR"; exit 1; }
WGS_DB_DIR_ABS="$(abs_path "$WGS_DB_DIR")" || { echo "[Go_shortWGS] WGS_DB_DIR not found: $WGS_DB_DIR"; exit 1; }
KRAKEN_DB_DIR_ABS="$(abs_path "$KRAKEN_DB_DIR")" || { echo "[Go_shortWGS] KRAKEN_DB_DIR not found: $KRAKEN_DB_DIR"; exit 1; }
GOWGS_DIR_ABS="$(abs_path "$GOWGS_DIR")" || { echo "[Go_shortWGS] GOWGS_DIR not found: $GOWGS_DIR"; exit 1; }

WORKDIR="$(pwd)"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"

run(){
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$WORKDIR":/work \
    -v "$SCRIPT_DIR":/pipeline:ro \
    -v "$FASTQ_DIR_ABS":/fastq:ro \
    -v "$WGS_DB_DIR_ABS":/db/wgs_db:ro \
    -v "$KRAKEN_DB_DIR_ABS":/db/kraken2:ro \
    -v "$GOWGS_DIR_ABS":/home/uhlemann/heekuk_path/GoWGS:ro \
    -w /work \
    "$IMAGE" \
    "$@"
}

BASE_ARGS=(
  snakemake
  --snakefile /pipeline/Go_shortWGS_V9_docker.smk
  --cores "$CORES"
  --rerun-incomplete
  --config
  fastq_dir=/fastq
  output_dir="$OUTPUT_DIR"
  wgs_db=/db/wgs_db
  kraken_db=/db/kraken2
)

if [ "$DRYRUN" -eq 1 ]; then
  BASE_ARGS+=(--dry-run)
fi

set +e
run "${BASE_ARGS[@]}" 2>&1 | tee shortwgs.log
rc=${PIPESTATUS[0]}
set -e

if [ "$rc" -ne 0 ] && grep -qiE "lock|unlock|LockException|cannot be locked" shortwgs.log; then
  echo "[Go_shortWGS] Detected lock issue -> running --unlock then retry..."
  run "${BASE_ARGS[@]}" --unlock

  set +e
  run "${BASE_ARGS[@]}" 2>&1 | tee -a shortwgs.log
  rc=${PIPESTATUS[0]}
  set -e
fi

exit "$rc"
