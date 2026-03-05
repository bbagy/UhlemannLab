#!/usr/bin/env bash
set -euo pipefail

usage(){
  echo "Usage: $0 -i READ_DIR -o PROJECT -g GENOME -a GFF [-s SNAKEDIR] [-c CORES] [-m IMAGE] [-n] [-K] [-P 0|1]"
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

prefilter_rnaseq_fastqs(){
  local input_dir="$1"
  local bad_dir
  local bad_log
  local marker
  local -a files=()
  local f dst
  local bad_count=0
  local valid_count=0

  if [ ! -d "$input_dir" ]; then
    echo "[Go_Rnake] READ_DIR not found: $input_dir"
    exit 1
  fi

  bad_dir="$(dirname "$input_dir")/0_bad_fastqs"
  bad_log="$bad_dir/moved_bad_fastqs.tsv"
  marker="$bad_dir/DONE.txt"
  mkdir -p "$bad_dir"
  [ -f "$bad_log" ] || echo -e "file_name\tfailure_reason" > "$bad_log"

  if [ -f "$marker" ] && ! find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) -newer "$marker" | grep -q .; then
    local n
    n=$(find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | wc -l | tr -d ' ')
    echo "[Go_Rnake][Prefilter] PASS (cached): no newer FASTQ detected (n=$n)"
    return
  fi

  mapfile -t files < <(find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)
  if [ "${#files[@]}" -eq 0 ]; then
    echo "[Go_Rnake][Prefilter][FATAL] No *.fastq.gz/*.fq.gz in $input_dir"
    exit 1
  fi
  echo "[Go_Rnake][Prefilter] START n=${#files[@]}"

  for f in "${files[@]}"; do
    if ! gzip -t "$f" >/dev/null 2>&1; then
      dst="$bad_dir/$(basename "$f")"
      if [ -e "$dst" ]; then
        dst="$bad_dir/$(basename "$f").$(date +%Y%m%d_%H%M%S)"
      fi
      mv "$f" "$dst"
      echo -e "$(basename "$f")\tgzip_integrity_failed" >> "$bad_log"
      echo "[Go_Rnake][Prefilter] FAIL $(basename "$f") -> gzip_integrity_failed"
      bad_count=$((bad_count + 1))
    fi
  done

  valid_count=$(find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | wc -l | tr -d ' ')
  if [ "$valid_count" -eq 0 ]; then
    echo "[Go_Rnake][Prefilter][FATAL] no valid FASTQ left after filtering."
    exit 1
  fi

  {
    echo "prefilter_done_at=$(date '+%Y-%m-%d %H:%M:%S')"
    echo "input_dir=$input_dir"
    echo "valid_fastq=$valid_count"
    echo "moved_bad_fastq=$bad_count"
  } > "$marker"
  echo "[Go_Rnake][Prefilter] DONE valid=$valid_count bad=$bad_count"
}

READ_DIR=""; PROJECT=""; GENOME=""; GFF=""; SNAKEDIR=""
CORES=8
IMAGE="rnake:latest"
DRYRUN=0
KEEP_GOING=0
SHOW_PROGRESS=1
PROGRESS_INTERVAL=60
PROGRESS_PID=""

while getopts "i:o:g:a:s:c:m:nKP:" opt; do
  case $opt in
    i) READ_DIR="$OPTARG" ;;
    o) PROJECT="$OPTARG" ;;
    g) GENOME="$OPTARG" ;;
    a) GFF="$OPTARG" ;;
    s) SNAKEDIR="$OPTARG" ;;
    c) CORES="$OPTARG" ;;
    m) IMAGE="$OPTARG" ;;
    n) DRYRUN=1 ;;
    K) KEEP_GOING=1 ;;
    P) SHOW_PROGRESS="$OPTARG" ;;
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

prefilter_rnaseq_fastqs "$READ_DIR_ABS"

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
  local outdir now total bad trim map counts merged
  outdir="${PROJECT}_RNAseq_output"
  now="$(date '+%H:%M:%S')"
  total="$(total_samples "$READ_DIR_ABS")"
  bad="$(awk 'NR>1{n++} END{print n+0}' "$(dirname "$READ_DIR_ABS")/0_bad_fastqs/moved_bad_fastqs.tsv" 2>/dev/null || echo 0)"
  trim=$(find "$outdir/1_trim" -maxdepth 1 -type f -name "*.R1.paired.output.fastq.gz" 2>/dev/null | wc -l | tr -d ' ')
  map=$(find "$outdir/3_bowtie2_files" -maxdepth 1 -type f -name "*.sam" 2>/dev/null | wc -l | tr -d ' ')
  counts=$(find "$outdir/5_counts" -maxdepth 1 -type f -name "*.counts" 2>/dev/null | wc -l | tr -d ' ')
  merged="0"; [ -f "$outdir/merged_counts_with_gene_names.csv" ] && merged="1"
  printf "[Go_Rnake][Progress][%s][%s] total=%s bad=%s trim=%s map=%s counts=%s merged=%s\n" \
    "$stage" "$now" "$total" "$bad" \
    "$(fmt_stage "$trim" "$total")" "$(fmt_stage "$map" "$total")" "$(fmt_stage "$counts" "$total")" "$merged"
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
  [ -d "$SNAKEDIR" ] || { echo "[Go_Rnake] SNAKEDIR not found: $SNAKEDIR"; exit 1; }
  PIPELINE_DIR="$(cd "$SNAKEDIR" && pwd -P)"
fi
SNAKEFILE_NAME="Go_bacteriaRNake_paired_V4.smk"
if [ ! -f "$PIPELINE_DIR/$SNAKEFILE_NAME" ]; then
  echo "[Go_Rnake][FATAL] Snakefile not found: $PIPELINE_DIR/$SNAKEFILE_NAME"
  exit 1
fi
echo "[Go_Rnake] Using Snakefile: $PIPELINE_DIR/$SNAKEFILE_NAME"

run(){
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$WORKDIR":/work \
    -v "$PIPELINE_DIR":/pipeline:ro \
    -v "$READ_DIR_ABS":/reads:ro \
    -v "$GENOME_ABS":/refs/genome.fna:ro \
    -v "$GFF_ABS":/refs/annotation.gff:ro \
    -w /work \
    "$IMAGE" \
    "$@"
}

BASE_ARGS=(
  --snakefile "/pipeline/$SNAKEFILE_NAME"
  --cores "$CORES"
  --rerun-incomplete
  --rerun-triggers mtime
  --config
  project="$PROJECT"
  read_dir=/reads
  genome=/refs/genome.fna
  gff=/refs/annotation.gff
)

if [ "$DRYRUN" -eq 1 ]; then
  BASE_ARGS+=(--dry-run)
fi
if [ "$KEEP_GOING" -eq 1 ]; then
  BASE_ARGS+=(--keep-going)
fi

on_interrupt(){
  echo "[Go_Rnake] Interrupt received. Stopping running jobs..."
  stop_progress_monitor || true
  pkill -TERM -P $$ >/dev/null 2>&1 || true
  sleep 1
  pkill -KILL -P $$ >/dev/null 2>&1 || true
  exit 130
}
trap on_interrupt INT TERM

set +e
start_progress_monitor "main"
run "${BASE_ARGS[@]}" 2>&1 | tee rnake.log
rc=${PIPESTATUS[0]}
stop_progress_monitor
set -e

if [ $rc -ne 0 ] && grep -qiE "lock|unlock|LockException|cannot be locked" rnake.log; then
  echo "[Go_Rnake] Detected lock issue -> running --unlock then retry..."
  run "${BASE_ARGS[@]}" --unlock

  set +e
  start_progress_monitor "retry"
  run "${BASE_ARGS[@]}" 2>&1 | tee -a rnake.log
  rc=${PIPESTATUS[0]}
  stop_progress_monitor
  set -e
fi

exit $rc
