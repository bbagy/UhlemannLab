#!/usr/bin/env bash
set -euo pipefail

usage(){
  echo "Usage: $0 -i FASTQ_DIR -o OUTPUT_DIR -n NUCLEOTIDE_DB -p PROTEIN_DB [-b METAPHLAN_DB] [-I METAPHLAN_INDEX] [-s SNAKEDIR] [-c CORES] [-j JOBS] [-t THREADS] [-m IMAGE] [-x] [-K] [--run-musicc] [--skip-gene-norm] [--skip-path-split] [--skip-pathcoverage]"
  echo "  Recommended MetaPhlAn input: -b /path/to/metaphlan_db_dir -I mpa_vJan25_CHOCOPhlAnSGB_202503"
  echo "  Backward-compatible shortcut: -b /path/to/mpa_vJan25_CHOCOPhlAnSGB_202503.pkl"
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
NUCLEOTIDE_DB=""
PROTEIN_DB=""
METAPHLAN_DB=""
METAPHLAN_INDEX=""
SNAKEDIR=""
CORES=8
JOBS=4
THREADS=4
IMAGE="humann:1.0"
DRYRUN=0
KEEP_GOING=0
RUN_GENE_NORM=1
RUN_PATH_SPLIT=1
RUN_PATHCOVERAGE=1
RUN_MUSICC=0

while [ "$#" -gt 0 ]; do
  case "$1" in
    -i)
      FASTQ_DIR="${2:-}"
      shift 2
      ;;
    -o)
      OUTPUT_DIR="${2:-}"
      shift 2
      ;;
    -n)
      NUCLEOTIDE_DB="${2:-}"
      shift 2
      ;;
    -p)
      PROTEIN_DB="${2:-}"
      shift 2
      ;;
    -b)
      METAPHLAN_DB="${2:-}"
      shift 2
      ;;
    -I|--metaphlan-index)
      METAPHLAN_INDEX="${2:-}"
      shift 2
      ;;
    -s)
      SNAKEDIR="${2:-}"
      shift 2
      ;;
    -c)
      CORES="${2:-}"
      shift 2
      ;;
    -j)
      JOBS="${2:-}"
      shift 2
      ;;
    -t)
      THREADS="${2:-}"
      shift 2
      ;;
    -m)
      IMAGE="${2:-}"
      shift 2
      ;;
    -x)
      DRYRUN=1
      shift
      ;;
    -K)
      KEEP_GOING=1
      shift
      ;;
    --run-musicc)
      RUN_MUSICC=1
      shift
      ;;
    --skip-gene-norm)
      RUN_GENE_NORM=0
      shift
      ;;
    --skip-path-split)
      RUN_PATH_SPLIT=0
      shift
      ;;
    --skip-pathcoverage)
      RUN_PATHCOVERAGE=0
      shift
      ;;
    --help|-h)
      usage
      ;;
    *)
      usage
      ;;
  esac
done

[ -z "$FASTQ_DIR" ] && usage
[ -z "$OUTPUT_DIR" ] && usage
[ -z "$NUCLEOTIDE_DB" ] && usage
[ -z "$PROTEIN_DB" ] && usage

FASTQ_DIR_ABS="$(abs_path "$FASTQ_DIR")" || { echo "[Humann] FASTQ_DIR not found: $FASTQ_DIR"; exit 1; }
NUCLEOTIDE_DB_ABS="$(abs_path "$NUCLEOTIDE_DB")" || { echo "[Humann] NUCLEOTIDE_DB not found: $NUCLEOTIDE_DB"; exit 1; }
PROTEIN_DB_ABS="$(abs_path "$PROTEIN_DB")" || { echo "[Humann] PROTEIN_DB not found: $PROTEIN_DB"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
PIPELINE_DIR="$SCRIPT_DIR"
if [ -n "$SNAKEDIR" ]; then
  [ -d "$SNAKEDIR" ] || { echo "[Humann] SNAKEDIR not found: $SNAKEDIR"; exit 1; }
  PIPELINE_DIR="$(cd "$SNAKEDIR" && pwd -P)"
fi

SNAKEFILE_NAME="Go_Humann.smk"
if [ ! -f "$PIPELINE_DIR/$SNAKEFILE_NAME" ] && [ -f "$PIPELINE_DIR/Go_Humann_V1.smk" ]; then
  SNAKEFILE_NAME="Go_Humann_V1.smk"
fi
if [ ! -f "$PIPELINE_DIR/$SNAKEFILE_NAME" ]; then
  echo "[Humann][FATAL] Snakefile not found: $PIPELINE_DIR/$SNAKEFILE_NAME"
  exit 1
fi

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
  echo "[Humann][FATAL] Docker image not found locally: $IMAGE"
  echo "[Humann] Build example:"
  echo "  cd \"$SCRIPT_DIR\" && docker build --network=host -t humann:1.0 ."
  exit 1
fi

WORKDIR="$(pwd)"
CONTAINER_HOME="/work/.codex_home_humann"

DOCKER_ARGS=(
  docker run --rm
  -u "$(id -u):$(id -g)"
  -v "$WORKDIR":/work
  -v "$PIPELINE_DIR":/pipeline:ro
  -v "$FASTQ_DIR_ABS":/fastq:ro
  -v "$NUCLEOTIDE_DB_ABS":/db/nucleotide:ro
  -v "$PROTEIN_DB_ABS":/db/protein:ro
  -e HOME="$CONTAINER_HOME"
  -e XDG_CACHE_HOME="$CONTAINER_HOME/.cache"
  -w /work
)

if [ -n "$METAPHLAN_DB" ]; then
  METAPHLAN_DB_ABS="$(abs_path "$METAPHLAN_DB")" || { echo "[Humann] METAPHLAN_DB not found: $METAPHLAN_DB"; exit 1; }
  if [ -f "$METAPHLAN_DB_ABS" ]; then
    case "$METAPHLAN_DB_ABS" in
      *.pkl)
        if [ -z "$METAPHLAN_INDEX" ]; then
          METAPHLAN_INDEX="$(basename "$METAPHLAN_DB_ABS" .pkl)"
        fi
        METAPHLAN_DB_ABS="$(dirname "$METAPHLAN_DB_ABS")"
        ;;
      *)
        echo "[Humann] METAPHLAN_DB file must be a MetaPhlAn .pkl index file: $METAPHLAN_DB"
        exit 1
        ;;
    esac
  fi
  DOCKER_ARGS+=(-v "$METAPHLAN_DB_ABS":/db/metaphlan:ro)
  DOCKER_ARGS+=(-e METAPHLAN_DB_DIR=/db/metaphlan)
fi

run(){
  "${DOCKER_ARGS[@]}" "$IMAGE" "$@"
}

mkdir -p "$WORKDIR/.codex_home_humann/.cache"

BASE_ARGS=(
  snakemake
  --snakefile "/pipeline/$SNAKEFILE_NAME"
  --config
  fastq_dir=/fastq
  output_dir="$OUTPUT_DIR"
  nucleotide_db=/db/nucleotide
  protein_db=/db/protein
  metaphlan_db="${METAPHLAN_DB:+/db/metaphlan}"
  metaphlan_index="$METAPHLAN_INDEX"
  humann_threads="$THREADS"
  run_musicc="$RUN_MUSICC"
  run_gene_norm="$RUN_GENE_NORM"
  run_path_split="$RUN_PATH_SPLIT"
  run_pathcoverage_merge="$RUN_PATHCOVERAGE"
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
run "${BASE_ARGS[@]}" 2>&1 | tee humann.log
rc=${PIPESTATUS[0]}
set -e

if [ "$rc" -ne 0 ] && grep -qiE "lock|unlock|LockException|cannot be locked" humann.log; then
  echo "[Humann] Detected lock issue -> running --unlock then retry..."
  run "${BASE_ARGS[@]}" --unlock
  set +e
  run "${BASE_ARGS[@]}" 2>&1 | tee -a humann.log
  rc=${PIPESTATUS[0]}
  set -e
fi

exit "$rc"
