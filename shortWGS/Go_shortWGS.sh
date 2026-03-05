#!/usr/bin/env bash
set -euo pipefail

usage(){
  echo "Usage: $0 -i FASTQ_DIR -o OUTPUT_DIR -d WGS_DB_DIR -k KRAKEN_DB_DIR -r GOWGS_DIR [-s SNAKEDIR] [-c CORES] [-m IMAGE] [-n] [-K] [-P 0|1]"
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

prefilter_core_init(){
  local input_dir="$1"
  local bad_dir="$2"
  local bad_log="$3"
  local marker="$4"
  if [ ! -d "$input_dir" ]; then
    echo "[Go_shortWGS] FASTQ_DIR not found: $input_dir"
    exit 1
  fi
  mkdir -p "$bad_dir"
  [ -f "$bad_log" ] || echo -e "file_name\tfailure_reason" > "$bad_log"
  if [ -f "$marker" ] && ! find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) -newer "$marker" | grep -q .; then
    local n
    n=$(find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | wc -l | tr -d ' ')
    echo "[Go_shortWGS][Prefilter] PASS (cached): no newer FASTQ detected (n=$n)"
    return 0
  fi
  return 1
}

prefilter_move_bad(){
  local f="$1"
  local reason="$2"
  local bad_dir="$3"
  local bad_log="$4"
  local dst="$bad_dir/$(basename "$f")"
  if [ -e "$dst" ]; then
    dst="$bad_dir/$(basename "$f").$(date +%Y%m%d_%H%M%S)"
  fi
  mv "$f" "$dst"
  echo -e "$(basename "$f")\t$reason" >> "$bad_log"
  echo "[Go_shortWGS][Prefilter] FAIL $(basename "$f") -> $reason"
}

prefilter_illumina_pairs(){
  local input_dir="$1"
  local bad_dir
  local bad_log
  local marker
  local -a files=()
  local f b sample
  local bad_count=0
  local valid_count=0
  local size_warn=0
  local tiny_threshold=$((1 * 1024 * 1024))

  bad_dir="$(dirname "$input_dir")/0_bad_fastqs"
  bad_log="$bad_dir/moved_bad_fastqs.tsv"
  marker="$bad_dir/DONE.txt"

  if prefilter_core_init "$input_dir" "$bad_dir" "$bad_log" "$marker"; then
    return
  fi

  mapfile -t files < <(find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)
  if [ "${#files[@]}" -eq 0 ]; then
    echo "[Go_shortWGS][Prefilter][FATAL] No *.fastq.gz/*.fq.gz in $input_dir"
    exit 1
  fi

  echo "[Go_shortWGS][Prefilter] START n=${#files[@]} pairs_required=1"

  for f in "${files[@]}"; do
    if ! gzip -t "$f" >/dev/null 2>&1; then
      prefilter_move_bad "$f" "gzip_integrity_failed" "$bad_dir" "$bad_log"
      bad_count=$((bad_count + 1))
    fi
  done

  declare -A r1_map=()
  declare -A r2_map=()
  declare -A unknown_map=()
  mapfile -t files < <(find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)
  for f in "${files[@]}"; do
    b="$(basename "$f")"
    if [[ "$b" =~ _R1_001\.(fastq|fq)\.gz$ ]]; then
      sample="${b%_R1_001.fastq.gz}"; sample="${sample%_R1_001.fq.gz}"; r1_map["$sample"]="$f"
    elif [[ "$b" =~ _R2_001\.(fastq|fq)\.gz$ ]]; then
      sample="${b%_R2_001.fastq.gz}"; sample="${sample%_R2_001.fq.gz}"; r2_map["$sample"]="$f"
    elif [[ "$b" =~ _R1\.(fastq|fq)\.gz$ ]]; then
      sample="${b%_R1.fastq.gz}"; sample="${sample%_R1.fq.gz}"; r1_map["$sample"]="$f"
    elif [[ "$b" =~ _R2\.(fastq|fq)\.gz$ ]]; then
      sample="${b%_R2.fastq.gz}"; sample="${sample%_R2.fq.gz}"; r2_map["$sample"]="$f"
    elif [[ "$b" =~ \.R1\.(fastq|fq)\.gz$ ]]; then
      sample="${b%.R1.fastq.gz}"; sample="${sample%.R1.fq.gz}"; r1_map["$sample"]="$f"
    elif [[ "$b" =~ \.R2\.(fastq|fq)\.gz$ ]]; then
      sample="${b%.R2.fastq.gz}"; sample="${sample%.R2.fq.gz}"; r2_map["$sample"]="$f"
    else
      sample="${b%.fastq.gz}"; sample="${sample%.fq.gz}"; unknown_map["$sample"]="$f"
    fi
  done

  for sample in "${!unknown_map[@]}"; do
    [ -e "${unknown_map[$sample]}" ] || continue
    prefilter_move_bad "${unknown_map[$sample]}" "paired_name_not_recognized" "$bad_dir" "$bad_log"
    bad_count=$((bad_count + 1))
  done
  for sample in "${!r1_map[@]}"; do
    if [ -z "${r2_map[$sample]:-}" ] && [ -e "${r1_map[$sample]}" ]; then
      prefilter_move_bad "${r1_map[$sample]}" "missing_pair_R2" "$bad_dir" "$bad_log"
      bad_count=$((bad_count + 1))
    fi
  done
  for sample in "${!r2_map[@]}"; do
    if [ -z "${r1_map[$sample]:-}" ] && [ -e "${r2_map[$sample]}" ]; then
      prefilter_move_bad "${r2_map[$sample]}" "missing_pair_R1" "$bad_dir" "$bad_log"
      bad_count=$((bad_count + 1))
    fi
  done

  mapfile -t files < <(find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)
  for f in "${files[@]}"; do
    valid_count=$((valid_count + 1))
    if [ "$(stat -c %s "$f" 2>/dev/null || stat -f %z "$f")" -lt "$tiny_threshold" ]; then
      size_warn=$((size_warn + 1))
    fi
  done
  if [ "$valid_count" -eq 0 ]; then
    echo "[Go_shortWGS][Prefilter][FATAL] no valid FASTQ left after filtering."
    exit 1
  fi

  {
    echo "prefilter_done_at=$(date '+%Y-%m-%d %H:%M:%S')"
    echo "input_dir=$input_dir"
    echo "valid_fastq=$valid_count"
    echo "moved_bad_fastq=$bad_count"
    echo "small_file_warn=$size_warn"
    echo "pairs_required=1"
  } > "$marker"
  echo "[Go_shortWGS][Prefilter] DONE valid=$valid_count bad=$bad_count warn_small=$size_warn"
}

FASTQ_DIR=""
OUTPUT_DIR=""
WGS_DB_DIR=""
KRAKEN_DB_DIR=""
GOWGS_DIR=""
SNAKEDIR=""
CORES=8
IMAGE="shortwgs:1.0"
DRYRUN=0
KEEP_GOING=0
SHOW_PROGRESS=1
PROGRESS_INTERVAL=60
PROGRESS_PID=""

while getopts "i:o:d:k:r:s:c:m:nKP:" opt; do
  case "$opt" in
    i) FASTQ_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    d) WGS_DB_DIR="$OPTARG" ;;
    k) KRAKEN_DB_DIR="$OPTARG" ;;
    r) GOWGS_DIR="$OPTARG" ;;
    s) SNAKEDIR="$OPTARG" ;;
    c) CORES="$OPTARG" ;;
    m) IMAGE="$OPTARG" ;;
    n) DRYRUN=1 ;;
    K) KEEP_GOING=1 ;;
    P) SHOW_PROGRESS="$OPTARG" ;;
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

prefilter_illumina_pairs "$FASTQ_DIR_ABS"

total_samples(){
  local d="$1"
  local n
  n=$(find "$d" -maxdepth 1 -type f \( -name "*_R1_001.fastq.gz" -o -name "*_R1.fastq.gz" -o -name "*.R1.fastq.gz" \) | wc -l | tr -d ' ')
  if [ "$n" -eq 0 ]; then
    n=$(find "$d" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | wc -l | tr -d ' ')
  fi
  echo "$n"
}

fmt_stage(){
  local n="$1"; local total="$2"
  if [ "$total" -gt 0 ] && [ "$n" -ge "$total" ]; then
    echo "complete"
  else
    echo "${n}/${total}"
  fi
}

progress_snapshot(){
  local stage="$1"
  local now total bad fastp kraken mlst tetyper
  local args_done plas_done summary_done
  now="$(date '+%H:%M:%S')"
  total="$(total_samples "$FASTQ_DIR_ABS")"
  bad="$(awk 'NR>1{n++} END{print n+0}' "$(dirname "$FASTQ_DIR_ABS")/0_bad_fastqs/moved_bad_fastqs.tsv" 2>/dev/null || echo 0)"
  fastp=$(find "$OUTPUT_DIR/1_fastp_out" -maxdepth 1 -type f -name "*.done.txt" 2>/dev/null | wc -l | tr -d ' ')
  kraken=$(find "$OUTPUT_DIR/2a_kraken2" -maxdepth 1 -type f -name "*.top_species.txt" 2>/dev/null | wc -l | tr -d ' ')
  mlst=$(find "$OUTPUT_DIR/2_ST_srst2_out" -maxdepth 1 -type f -name "*.mlst.done.txt" 2>/dev/null | wc -l | tr -d ' ')
  tetyper=$(find "$OUTPUT_DIR/5_TETyper" -maxdepth 1 -type f -name "*_summary.txt" 2>/dev/null | wc -l | tr -d ' ')
  args_done="0"; [ -f "$OUTPUT_DIR/3_ARGs_srst2_out/args_srst2.done.txt" ] && args_done="1"
  plas_done="0"; [ -f "$OUTPUT_DIR/4_Plasmid_srst2_out/plasmid_srst2.done.txt" ] && plas_done="1"
  summary_done="0"; [ -f "$OUTPUT_DIR/summary_report.html" ] && summary_done="1"
  printf "[Go_shortWGS][Progress][%s][%s] total=%s bad=%s fastp=%s kraken=%s mlst=%s tetyper=%s args=%s plas=%s summary=%s\n" \
    "$stage" "$now" "$total" "$bad" \
    "$(fmt_stage "$fastp" "$total")" "$(fmt_stage "$kraken" "$total")" "$(fmt_stage "$mlst" "$total")" "$(fmt_stage "$tetyper" "$total")" \
    "$args_done" "$plas_done" "$summary_done"
}

start_progress_monitor(){
  local stage="$1"
  local current=""
  local current_key=""
  local last=""
  local last_key=""
  if [ "$SHOW_PROGRESS" != "1" ] || [ "$DRYRUN" -eq 1 ]; then
    return
  fi
  current="$(progress_snapshot "$stage")"
  current_key="$(printf '%s\n' "$current" | sed -E 's/\[[0-9]{2}:[0-9]{2}:[0-9]{2}\]//')"
  echo "$current"
  last="$current"
  last_key="$current_key"
  (
    while true; do
      sleep "$PROGRESS_INTERVAL"
      current="$(progress_snapshot "$stage")"
      current_key="$(printf '%s\n' "$current" | sed -E 's/\[[0-9]{2}:[0-9]{2}:[0-9]{2}\]//')"
      if [ "$current_key" != "$last_key" ]; then
        echo "$current"
        last="$current"
        last_key="$current_key"
      fi
    done
  ) &
  PROGRESS_PID=$!
}

stop_progress_monitor(){
  if [ -n "${PROGRESS_PID:-}" ]; then
    kill "$PROGRESS_PID" >/dev/null 2>&1 || true
    wait "$PROGRESS_PID" 2>/dev/null || true
    PROGRESS_PID=""
  fi
}

WORKDIR="$(pwd)"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
PIPELINE_DIR="$SCRIPT_DIR"
if [ -n "$SNAKEDIR" ]; then
  [ -d "$SNAKEDIR" ] || { echo "[Go_shortWGS] SNAKEDIR not found: $SNAKEDIR"; exit 1; }
  PIPELINE_DIR="$(cd "$SNAKEDIR" && pwd -P)"
fi
SNAKEFILE_NAME="Go_shortWGS_V9_docker.smk"
if [ ! -f "$PIPELINE_DIR/$SNAKEFILE_NAME" ]; then
  echo "[Go_shortWGS][FATAL] Snakefile not found: $PIPELINE_DIR/$SNAKEFILE_NAME"
  exit 1
fi
echo "[Go_shortWGS] Using Snakefile: $PIPELINE_DIR/$SNAKEFILE_NAME"

run(){
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$WORKDIR":/work \
    -v "$PIPELINE_DIR":/pipeline:ro \
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
  --snakefile "/pipeline/$SNAKEFILE_NAME"
  --cores "$CORES"
  --rerun-incomplete
  --rerun-triggers mtime
  --config
  fastq_dir=/fastq
  output_dir="$OUTPUT_DIR"
  wgs_db=/db/wgs_db
  kraken_db=/db/kraken2
)

if [ "$DRYRUN" -eq 1 ]; then
  BASE_ARGS+=(--dry-run)
fi
if [ "$KEEP_GOING" -eq 1 ]; then
  BASE_ARGS+=(--keep-going)
fi

on_interrupt(){
  echo "[Go_shortWGS] Interrupt received. Stopping running jobs..."
  stop_progress_monitor || true
  pkill -TERM -P $$ >/dev/null 2>&1 || true
  sleep 1
  pkill -KILL -P $$ >/dev/null 2>&1 || true
  exit 130
}
trap on_interrupt INT TERM

set +e
start_progress_monitor "main"
run "${BASE_ARGS[@]}" 2>&1 | tee shortwgs.log
rc=${PIPESTATUS[0]}
stop_progress_monitor
set -e

if [ "$rc" -ne 0 ] && grep -qiE "lock|unlock|LockException|cannot be locked" shortwgs.log; then
  echo "[Go_shortWGS] Detected lock issue -> running --unlock then retry..."
  run "${BASE_ARGS[@]}" --unlock

  set +e
  start_progress_monitor "retry"
  run "${BASE_ARGS[@]}" 2>&1 | tee -a shortwgs.log
  rc=${PIPESTATUS[0]}
  stop_progress_monitor
  set -e
fi

exit "$rc"
