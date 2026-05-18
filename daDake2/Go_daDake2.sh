#!/usr/bin/env bash
# Go_daDake2.sh — DADA2 amplicon pipeline wrapper
# Types: standard_V3V4 | zymo_V1V2 | illumina_ITS
set -euo pipefail

# Resolve symlinks so SCRIPT_DIR points to the actual .sh location
_SELF="$(readlink -f "${BASH_SOURCE[0]}" 2>/dev/null \
         || realpath  "${BASH_SOURCE[0]}" 2>/dev/null \
         || echo      "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(cd "$(dirname "$_SELF")" && pwd -P)"

abs_path(){
  local p="$1"
  if   [ -d "$p" ]; then (cd "$p" && pwd -P)
  elif [ -e "$p" ]; then (cd "$(dirname "$p")" && printf "%s/%s\n" "$(pwd -P)" "$(basename "$p")")
  else echo "$p"
  fi
}

usage(){
  cat <<EOF
Usage: $(basename "$0") -t TYPE -i DIRS -d DB [options]

Required:
  -t TYPE     standard_V3V4 | zymo_V1V2 | illumina_ITS
  -i DIRS     comma-separated project directories (e.g. "ProjA,ProjB")
  -d DB       taxonomy DB path (SILVA or UNITE .fa.gz)

Optional:
  -s SNAKEPATH directory containing a Snakefile, or a specific .smk file
  -c CORES    CPU cores (default: 4)
  -m BYTES    min FASTQ size to include (default: 10000)
  -D          delete filtered FASTQ dirs after successful run
  -n          dry-run
  -K          keep-going on sample failure
  -h          show this help
EOF
  exit 1
}

TYPE=""
INPUT_DIRS=""
DB=""
SNAKEDIR=""
CORES=4
MIN_BYTES=10000
DRYRUN=0
KEEP_GOING=0
DELETE_FILT=0

while getopts "t:i:d:s:c:m:DnKh" opt; do
  case "$opt" in
    t) TYPE="$OPTARG" ;;
    i) INPUT_DIRS="$OPTARG" ;;
    d) DB="$OPTARG" ;;
    s) SNAKEDIR="$OPTARG" ;;
    c) CORES="$OPTARG" ;;
    m) MIN_BYTES="$OPTARG" ;;
    D) DELETE_FILT=1 ;;
    n) DRYRUN=1 ;;
    K) KEEP_GOING=1 ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "$TYPE" ]]       && { echo "[ERROR] -t TYPE is required"; usage; }
[[ -z "$INPUT_DIRS" ]] && { echo "[ERROR] -i DIRS is required"; usage; }
[[ -z "$DB" ]]         && { echo "[ERROR] -d DB is required";   usage; }

case "$TYPE" in
  standard_V3V4|zymo_V1V2|illumina_ITS) ;;
  *) echo "[ERROR] -t must be: standard_V3V4 | zymo_V1V2 | illumina_ITS"; exit 1 ;;
esac

# smk 위치 결정: -s 옵션 > 스크립트와 같은 디렉토리
if [[ -n "$SNAKEDIR" ]]; then
  SNAKEPATH="$(abs_path "$SNAKEDIR")"
  if [[ -d "$SNAKEPATH" ]]; then
    PIPELINE_DIR="$SNAKEPATH"
    SMK="${PIPELINE_DIR}/Go_daDake2.smk"
  elif [[ -f "$SNAKEPATH" ]]; then
    PIPELINE_DIR="$(cd "$(dirname "$SNAKEPATH")" && pwd -P)"
    SMK="$SNAKEPATH"
  else
    echo "[ERROR] -s path not found: $SNAKEDIR"
    exit 1
  fi
else
  PIPELINE_DIR="$SCRIPT_DIR"
  SMK="${PIPELINE_DIR}/Go_daDake2.smk"
fi

SCRIPTS_DIR="${PIPELINE_DIR}/scripts"

[[ -f "$SMK" ]] || {
  echo "[ERROR] Snakefile not found: $SMK"
  echo "        Use -s with a directory or a specific .smk file"
  exit 1
}

DB_ABS="$(abs_path "$DB")"
[[ -f "$DB_ABS" ]] || { echo "[ERROR] DB not found: $DB_ABS"; exit 1; }

RUN_DATE=$(date +%y%m%d)
CFG=$(mktemp /tmp/daDake2_cfg.XXXXXX.yaml)
trap "rm -f '$CFG'" EXIT

cat > "$CFG" <<YAML
type:            "${TYPE}"
project_dirs:    "${INPUT_DIRS}"
db:              "${DB_ABS}"
sdb:             ""
date:            "${RUN_DATE}"
min_fastq_bytes: ${MIN_BYTES}
scripts_dir:     "${SCRIPTS_DIR}"
YAML

echo "[Go_daDake2] type=${TYPE}  date=${RUN_DATE}  cores=${CORES}"
echo "[Go_daDake2] projects: ${INPUT_DIRS}"
echo "[Go_daDake2] db: ${DB_ABS}"
echo "[Go_daDake2] smk: ${SMK}"

SNAKE_ARGS=(
  snakemake
  --snakefile "$SMK"
  --configfile "$CFG"
  --cores "$CORES"
  --rerun-incomplete
  --rerun-triggers mtime
  --printshellcmds
)
[[ "$DRYRUN"    -eq 1 ]] && SNAKE_ARGS+=(--dry-run)
[[ "$KEEP_GOING" -eq 1 ]] && SNAKE_ARGS+=(--keep-going)

on_interrupt(){
  echo "[Go_daDake2] Interrupt received — stopping."
  pkill -TERM -P $$ >/dev/null 2>&1 || true
  exit 130
}
trap on_interrupt INT TERM

set +e
"${SNAKE_ARGS[@]}" 2>&1 | tee dada2_pipeline.log
RC=${PIPESTATUS[0]}
set -e

# auto-unlock on lock error
if [[ "$RC" -ne 0 ]] && grep -qiE "lock|LockException" dada2_pipeline.log; then
  echo "[Go_daDake2] Lock detected — unlocking and retrying..."
  snakemake --snakefile "$SMK" --configfile "$CFG" --unlock
  set +e
  "${SNAKE_ARGS[@]}" 2>&1 | tee -a dada2_pipeline.log
  RC=${PIPESTATUS[0]}
  set -e
fi

# cleanup filtered FASTQs on success
if [[ "$RC" -eq 0 && "$DELETE_FILT" -eq 1 ]]; then
  IFS=',' read -ra _PROJS <<< "$INPUT_DIRS"
  for _proj in "${_PROJS[@]}"; do
    _proj="${_proj// /}"
    for _d in "3_DADA2_filtered" "4_DADA2_filtered" "3_path.cut"; do
      _target="${_proj}_dada2/${_d}"
      if [[ -d "$_target" ]]; then
        rm -rf "$_target"
        echo "[Go_daDake2] Deleted: $_target"
      fi
    done
  done
fi

exit "$RC"
