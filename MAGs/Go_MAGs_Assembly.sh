#!/usr/bin/env bash
set -euo pipefail

usage(){
  echo "Usage: $0 -i QC_HOST_FILTERED_FASTQ_DIR -o OUTPUT_DIR [-d CHECKM_DATA_DIR] [-s MAGS_DIR] [-c CORES] [-j JOBS] [-t MEGAHIT_THREADS] [-M MEGAHIT_MEMORY] [-b BINNING_TOOLS] [-m IMAGE] [-n] [-K]"
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
CHECKM_DATA_DIR=""
MAGS_DIR=""
CORES=8
JOBS=4
MEGAHIT_THREADS=""
MEGAHIT_MEMORY=128000
BINNING_TOOLS="concoct,metabat2,maxbin2"
IMAGE="mags-assembly:1.0"
DRYRUN=0
KEEP_GOING=0

while getopts "i:o:d:s:c:j:t:M:b:m:nK" opt; do
  case "$opt" in
    i) FASTQ_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    d) CHECKM_DATA_DIR="$OPTARG" ;;
    s) MAGS_DIR="$OPTARG" ;;
    c) CORES="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    t) MEGAHIT_THREADS="$OPTARG" ;;
    M) MEGAHIT_MEMORY="$OPTARG" ;;
    b) BINNING_TOOLS="$OPTARG" ;;
    m) IMAGE="$OPTARG" ;;
    n) DRYRUN=1 ;;
    K) KEEP_GOING=1 ;;
    *) usage ;;
  esac
done

[ -z "$FASTQ_DIR" ] && usage
[ -z "$OUTPUT_DIR" ] && usage

FASTQ_DIR_ABS="$(abs_path "$FASTQ_DIR")" || { echo "[Go_MAGs_Assembly] FASTQ_DIR not found: $FASTQ_DIR"; exit 1; }
CHECKM_DATA_DIR_ABS=""
if [ -n "$CHECKM_DATA_DIR" ]; then
  CHECKM_DATA_DIR_ABS="$(abs_path "$CHECKM_DATA_DIR")" || { echo "[Go_MAGs_Assembly] CHECKM_DATA_DIR not found: $CHECKM_DATA_DIR"; exit 1; }
fi
if ! compgen -G "$FASTQ_DIR_ABS/*_R1_nohuman.fastq.gz" >/dev/null; then
  echo "[Go_MAGs_Assembly][FATAL] No *_R1_nohuman.fastq.gz files found in: $FASTQ_DIR_ABS"
  echo "[Go_MAGs_Assembly] Use the QC output host_filtered_fastq directory as -i."
  exit 1
fi
for r1 in "$FASTQ_DIR_ABS"/*_R1_nohuman.fastq.gz; do
  r2="${r1/_R1_nohuman.fastq.gz/_R2_nohuman.fastq.gz}"
  if [ ! -f "$r2" ]; then
    echo "[Go_MAGs_Assembly][FATAL] R2 pair missing for: $r1"
    echo "[Go_MAGs_Assembly] Expected: $r2"
    exit 1
  fi
done

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
  echo "[Go_MAGs_Assembly][FATAL] Docker image not found locally: $IMAGE"
  echo "[Go_MAGs_Assembly] Build example:"
  echo "  cd \"$(cd "$(dirname "$0")" && pwd -P)/docker/MAGs\" && docker build -f Dockerfile.assembly -t mags-assembly:1.0 ."
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
PIPELINE_DIR="$SCRIPT_DIR"
if [ -n "$MAGS_DIR" ]; then
  [ -d "$MAGS_DIR" ] || { echo "[Go_MAGs_Assembly] MAGS_DIR not found: $MAGS_DIR"; exit 1; }
  PIPELINE_DIR="$(cd "$MAGS_DIR" && pwd -P)"
fi

SNAKEFILE_NAME="Go_MAGs_Assembly.smk"
if [ ! -f "$PIPELINE_DIR/$SNAKEFILE_NAME" ] && [ -f "$PIPELINE_DIR/workflow/Go_MAGs_Assembly_V1.smk" ]; then
  SNAKEFILE_NAME="workflow/Go_MAGs_Assembly_V1.smk"
fi
if [ ! -f "$PIPELINE_DIR/$SNAKEFILE_NAME" ]; then
  echo "[Go_MAGs_Assembly][FATAL] Snakefile not found: $PIPELINE_DIR/$SNAKEFILE_NAME"
  exit 1
fi

if [ -z "$MEGAHIT_THREADS" ]; then
  MEGAHIT_THREADS="$CORES"
fi

WORKDIR="$(pwd)"

run(){
  local -a docker_args=(
    docker run --rm
    -u "$(id -u):$(id -g)"
    -v "$WORKDIR":/work
    -v "$PIPELINE_DIR":/pipeline:ro
    -v "$FASTQ_DIR_ABS":/fastq:ro
    -w /work
  )
  if [ -n "$CHECKM_DATA_DIR_ABS" ]; then
    docker_args+=(-v "$CHECKM_DATA_DIR_ABS":/db/checkm_data:ro)
  fi
  docker_args+=("$IMAGE")
  "${docker_args[@]}" "$@"
}

BASE_ARGS=(
  snakemake
  --snakefile "/pipeline/$SNAKEFILE_NAME"
  --config
  fastq_dir=/fastq
  output_dir="$OUTPUT_DIR"
  megahit_threads="$MEGAHIT_THREADS"
  megahit_memory="$MEGAHIT_MEMORY"
  binning_tools="$BINNING_TOOLS"
  --cores "$CORES"
  --jobs "$JOBS"
  --latency-wait 60
  --rerun-incomplete
)
if [ -n "$CHECKM_DATA_DIR_ABS" ]; then
  BASE_ARGS+=(checkm_data_dir=/db/checkm_data)
fi

if [ "$DRYRUN" -eq 1 ]; then
  BASE_ARGS+=(--dry-run)
fi
if [ "$KEEP_GOING" -eq 1 ]; then
  BASE_ARGS+=(--keep-going)
fi

set +e
run "${BASE_ARGS[@]}" 2>&1 | tee mags_assembly.log
rc=${PIPESTATUS[0]}
set -e

if [ "$rc" -ne 0 ] && grep -qiE "lock|unlock|LockException|cannot be locked" mags_assembly.log; then
  echo "[Go_MAGs_Assembly] Detected lock issue -> running --unlock then retry..."
  run "${BASE_ARGS[@]}" --unlock
  set +e
  run "${BASE_ARGS[@]}" 2>&1 | tee -a mags_assembly.log
  rc=${PIPESTATUS[0]}
  set -e
fi

exit "$rc"
