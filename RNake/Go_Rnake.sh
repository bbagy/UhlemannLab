#!/usr/bin/env bash
set -euo pipefail

usage(){
  echo "Usage: $0 -i READ_DIR -o PROJECT -g GENOME -a GFF [-c CORES] [-m IMAGE] [-n]"
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

READ_DIR=""; PROJECT=""; GENOME=""; GFF=""
CORES=8
IMAGE="rnake:latest"
DRYRUN=0

while getopts "i:o:g:a:c:m:n" opt; do
  case $opt in
    i) READ_DIR="$OPTARG" ;;
    o) PROJECT="$OPTARG" ;;
    g) GENOME="$OPTARG" ;;
    a) GFF="$OPTARG" ;;
    c) CORES="$OPTARG" ;;
    m) IMAGE="$OPTARG" ;;
    n) DRYRUN=1 ;;
    *) usage ;;
  esac
done

[ -z "$READ_DIR" ] && usage
[ -z "$PROJECT" ] && usage
[ -z "$GENOME" ] && usage
[ -z "$GFF" ] && usage

READ_DIR_ABS="$(abs_path "$READ_DIR")" || { echo "[Go_Rnake] READ_DIR not found: $READ_DIR"; exit 1; }
GENOME_ABS="$(abs_path "$GENOME")" || { echo "[Go_Rnake] GENOME not found: $GENOME"; exit 1; }
GFF_ABS="$(abs_path "$GFF")" || { echo "[Go_Rnake] GFF not found: $GFF"; exit 1; }

WORKDIR="$(pwd)"

run(){
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$WORKDIR":/work \
    -v "$READ_DIR_ABS":/reads:ro \
    -v "$GENOME_ABS":/refs/genome.fna:ro \
    -v "$GFF_ABS":/refs/annotation.gff:ro \
    -w /work \
    "$IMAGE" \
    "$@"
}

BASE_ARGS=(
  --cores "$CORES"
  --rerun-incomplete
  --config
  project="$PROJECT"
  read_dir=/reads
  genome=/refs/genome.fna
  gff=/refs/annotation.gff
)

if [ "$DRYRUN" -eq 1 ]; then
  BASE_ARGS+=(--dry-run)
fi

set +e
run "${BASE_ARGS[@]}" 2>&1 | tee rnake.log
rc=${PIPESTATUS[0]}
set -e

if [ $rc -ne 0 ] && grep -qiE "lock|unlock|LockException|cannot be locked" rnake.log; then
  echo "[Go_Rnake] Detected lock issue -> running --unlock then retry..."
  run "${BASE_ARGS[@]}" --unlock

  set +e
  run "${BASE_ARGS[@]}" 2>&1 | tee -a rnake.log
  rc=${PIPESTATUS[0]}
  set -e
fi

exit $rc
