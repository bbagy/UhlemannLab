#!/usr/bin/env bash
set -euo pipefail

usage(){
  echo "Usage: $0 -i FASTQ_DIR -o OUTPUT -d REF_FASTA [-g GFF] [-c CORES] [-t THREADS] [-p 1|2] [-m IMAGE] [-n]"
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

FASTQ_DIR=""; OUTPUT=""; REF_FASTA=""; GFF=""
CORES=8
THREADS=4
PAIRED=2
IMAGE="pf-snake:1.0"
DRYRUN=0

while getopts "i:o:d:g:c:t:p:m:n" opt; do
  case $opt in
    i) FASTQ_DIR="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    d) REF_FASTA="$OPTARG" ;;
    g) GFF="$OPTARG" ;;
    c) CORES="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    p) PAIRED="$OPTARG" ;;
    m) IMAGE="$OPTARG" ;;
    n) DRYRUN=1 ;;
    *) usage ;;
  esac
done

[ -z "$FASTQ_DIR" ] && usage
[ -z "$OUTPUT" ] && usage
[ -z "$REF_FASTA" ] && usage

FASTQ_DIR_ABS="$(abs_path "$FASTQ_DIR")" || { echo "[Go_PFsnake] FASTQ_DIR not found: $FASTQ_DIR"; exit 1; }
REF_FASTA_ABS="$(abs_path "$REF_FASTA")" || { echo "[Go_PFsnake] REF_FASTA not found: $REF_FASTA"; exit 1; }

if [ "$PAIRED" != "1" ] && [ "$PAIRED" != "2" ]; then
  echo "[Go_PFsnake] -p must be 1 (SE) or 2 (PE)"
  exit 1
fi

WORKDIR="$(pwd)"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"

run(){
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$WORKDIR":/work \
    -v "$SCRIPT_DIR":/pipeline:ro \
    -v "$FASTQ_DIR_ABS":/fastq:ro \
    -v "$REF_FASTA_ABS":/ref/ref.fasta:ro \
    -w /work \
    "$IMAGE" \
    "$@"
}

BASE_ARGS=(
  snakemake
  --snakefile /pipeline/Go_Pfsnake_V9.2_docker.smk
  --cores "$CORES"
  --rerun-incomplete
  --config
  Fastq_DIRS=/fastq
  output="$OUTPUT"
  hostDB=/ref/ref.fasta
  threads="$THREADS"
  paired="$PAIRED"
)

if [ "$DRYRUN" -eq 1 ]; then
  BASE_ARGS+=(--dry-run)
fi

if [ -n "$GFF" ]; then
  GFF_ABS="$(abs_path "$GFF")" || { echo "[Go_PFsnake] GFF not found: $GFF"; exit 1; }
  BASE_ARGS+=(gff=/work/.pfsnake.gff)
  run_with_gff(){
    docker run --rm \
      -u "$(id -u):$(id -g)" \
      -v "$WORKDIR":/work \
      -v "$SCRIPT_DIR":/pipeline:ro \
      -v "$FASTQ_DIR_ABS":/fastq:ro \
      -v "$REF_FASTA_ABS":/ref/ref.fasta:ro \
      -v "$GFF_ABS":/work/.pfsnake.gff:ro \
      -w /work \
      "$IMAGE" \
      "$@"
  }
else
  run_with_gff(){
    run "$@"
  }
fi

set +e
run_with_gff "${BASE_ARGS[@]}" 2>&1 | tee pfsnake.log
rc=${PIPESTATUS[0]}
set -e

if [ $rc -ne 0 ] && grep -qiE "lock|unlock|LockException|cannot be locked" pfsnake.log; then
  echo "[Go_PFsnake] Detected lock issue -> running --unlock then retry..."
  run_with_gff "${BASE_ARGS[@]}" --unlock

  set +e
  run_with_gff "${BASE_ARGS[@]}" 2>&1 | tee -a pfsnake.log
  rc=${PIPESTATUS[0]}
  set -e
fi

exit $rc
