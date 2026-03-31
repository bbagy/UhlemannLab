#!/usr/bin/env bash
set -euo pipefail
usage(){ echo "Usage: $0 -i INPUT -o OUTPUT -d DB [-m MAP_FILE] [-s SNAKEDIR] [-p 0|1] [-n] [-K]"; exit 1; }

INPUT=""; OUTPUT=""; DB=""; SNAKEDIR=""; PORECHOP=0
MAP_FILE=""
DRYRUN=0
KEEP_GOING=0
PROGRESS_INTERVAL=60
PROGRESS_PID=""
DEFAULT_THREADS=2
DEFAULT_CORES=4
DEFAULT_AUTOCYCLER_JOBS=1
PLASSEMBLER_DB_PATH="${PLASSEMBLER_DB_PATH:-/db/plassembler_db}"
MEDAKA_DATA_PATH="${MEDAKA_DATA_PATH:-/db/medaka_models}"

is_nonneg_int() {
  [[ "${1:-}" =~ ^[0-9]+$ ]]
}
while getopts "i:o:d:s:p:nK" opt; do
  case $opt in
    i) INPUT="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    d) DB="$OPTARG" ;;
    s) SNAKEDIR="$OPTARG" ;;
    p) PORECHOP="$OPTARG" ;;
    n) DRYRUN=1 ;;
    K) KEEP_GOING=1 ;;
    *) usage ;;
  esac
done
[ -z "$INPUT" ] && usage; [ -z "$OUTPUT" ] && usage; [ -z "$DB" ] && usage

file_size_bytes() {
  local f="$1"
  if stat -c %s "$f" >/dev/null 2>&1; then
    stat -c %s "$f"
  else
    stat -f %z "$f"
  fi
}

human_size() {
  local b="$1"
  awk -v B="$b" 'BEGIN{
    split("B KB MB GB TB", u, " ");
    i=1;
    while (B>=1024 && i<5) { B/=1024; i++ }
    printf "%.2f %s", B, u[i];
  }'
}

preflight_check_fastqs() {
  local input_dir="$1"
  local -a fastqs=()
  local -a current_fastqs=()
  local total_bytes=0
  local bad_count=0
  local valid_count=0
  local tiny_count=0
  local bad_dir
  local moved_path
  local marker
  local bad_log
  local current_count
  local idx
  local msg
  local tiny_threshold=$((5 * 1024 * 1024))
  local f sz

  if [ ! -d "$input_dir" ]; then
    echo "[Go_longWGS] INPUT directory not found: $input_dir"
    exit 1
  fi

  bad_dir="$(dirname "$input_dir")/0_bad_fastqs"
  mkdir -p "$bad_dir"
  marker="$bad_dir/DONE.txt"
  bad_log="$bad_dir/moved_bad_fastqs.tsv"

  if [ ! -f "$bad_log" ]; then
    echo -e "file_name\tfailure_reason" > "$bad_log"
  fi

  mapfile -t current_fastqs < <(find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)
  current_count="${#current_fastqs[@]}"
  if [ "$current_count" -eq 0 ]; then
    echo "[Go_longWGS] No *.fastq.gz/*.fq.gz files found in: $input_dir"
    exit 1
  fi

  # Skip re-check if prefilter already done and no FASTQ was updated after marker.
  if [ -f "$marker" ] && ! find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) -newer "$marker" | grep -q .; then
    echo "[Go_longWGS][Prefilter] PASS (cached): no newer FASTQ detected"
    echo "[Go_longWGS][Prefilter] marker: $marker"
    echo "[Go_longWGS][Prefilter] current FASTQ files: $current_count"
    return
  fi

  fastqs=("${current_fastqs[@]}")

  echo "[Go_longWGS][Prefilter] START: ${#fastqs[@]} FASTQ files"
  idx=0
  for f in "${fastqs[@]}"; do
    idx=$((idx + 1))
    msg="[${idx}/${#fastqs[@]}] $(basename "$f")"
    echo "[Go_longWGS][Prefilter] CHECK $msg"

    if ! gzip -t "$f" >/dev/null 2>&1; then
      moved_path="$bad_dir/$(basename "$f")"
      if [ -e "$moved_path" ]; then
        moved_path="$bad_dir/$(basename "$f").$(date +%Y%m%d_%H%M%S)"
      fi
      mv "$f" "$moved_path"
      echo -e "$(basename "$f")\tgzip_integrity_failed" >> "$bad_log"
      echo "[Go_longWGS][Prefilter] FAIL $msg -> moved to: $moved_path"
      bad_count=$((bad_count + 1))
      echo "[Go_longWGS][Prefilter] PROGRESS: ok=$valid_count bad=$bad_count warn_small=$tiny_count"
      continue
    fi

    valid_count=$((valid_count + 1))
    sz=$(file_size_bytes "$f")
    total_bytes=$((total_bytes + sz))
    if [ "$sz" -lt "$tiny_threshold" ]; then
      echo "[Go_longWGS][Prefilter] WARN $msg size=$(human_size "$sz") (< 5.00 MB)"
      tiny_count=$((tiny_count + 1))
    else
      echo "[Go_longWGS][Prefilter] PASS $msg size=$(human_size "$sz")"
    fi
    echo "[Go_longWGS][Prefilter] PROGRESS: ok=$valid_count bad=$bad_count warn_small=$tiny_count"
  done

  echo "[Go_longWGS][Prefilter] DONE"
  echo "[Go_longWGS][Prefilter] total input size: $(human_size "$total_bytes") (${total_bytes} bytes)"
  echo "[Go_longWGS][Prefilter] valid FASTQ files: $valid_count / ${#fastqs[@]}"
  if [ "$tiny_count" -gt 0 ]; then
    echo "[Go_longWGS][Prefilter] small-file warnings: $tiny_count"
  fi
  if [ "$bad_count" -gt 0 ]; then
    echo "[Go_longWGS][Prefilter] moved bad FASTQ files: $bad_count (dir: $bad_dir)"
  fi
  if [ "$valid_count" -eq 0 ]; then
    echo "[Go_longWGS][Prefilter][FATAL] no valid FASTQ files remain after integrity filtering."
    exit 1
  fi

  {
    echo "prefilter_done_at=$(date '+%Y-%m-%d %H:%M:%S')"
    echo "input_dir=$input_dir"
    echo "valid_fastq=$valid_count"
    echo "moved_bad_fastq=$bad_count"
    echo "small_file_warn=$tiny_count"
  } > "$marker"
  echo "[Go_longWGS][Prefilter] wrote marker: $marker"
}

preflight_check_fastqs "$INPUT"

WORKDIR="$(pwd)"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
PIPELINE_DIR="$SCRIPT_DIR"
if [ -n "$SNAKEDIR" ]; then
  if [ ! -d "$SNAKEDIR" ]; then
    echo "[Go_longWGS] SNAKEDIR not found: $SNAKEDIR"
    exit 1
  fi
  PIPELINE_DIR="$(cd "$SNAKEDIR" && pwd -P)"
fi

dir_count() {
  local d="$1"
  if [ -d "$d" ]; then
    find "$d" -mindepth 1 -maxdepth 1 2>/dev/null | wc -l | tr -d ' '
  else
    echo 0
  fi
}

stage_done_count() {
  local stage="$1"
  case "$stage" in
    qc)
      find "$OUTPUT/1_QC" -maxdepth 1 -type f -name "*.clean.fastq.gz" 2>/dev/null | wc -l | tr -d ' '
      ;;
    quast)
      find "$OUTPUT/2_quast" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' '
      ;;
    autocycler)
      find "$OUTPUT/3_autocycler" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' '
      ;;
    medaka)
      find "$OUTPUT/4_medaka" -maxdepth 1 -type f -name "*_final_assembly.fasta" 2>/dev/null | wc -l | tr -d ' '
      ;;
    checkm2)
      find "$OUTPUT/5_checkm2" -mindepth 2 -maxdepth 2 -type f -name "DONE.txt" 2>/dev/null | wc -l | tr -d ' '
      ;;
    coverage)
      find "$OUTPUT/6_coverage" -mindepth 2 -maxdepth 2 -type f -name "DONE.txt" 2>/dev/null | wc -l | tr -d ' '
      ;;
    kraken2)
      find "$OUTPUT/7a_kraken2" -maxdepth 1 -type f -name "*.done.txt" 2>/dev/null | wc -l | tr -d ' '
      ;;
    bakta)
      find "$OUTPUT/7_bakta" -mindepth 2 -maxdepth 2 -type f -name "DONE.txt" 2>/dev/null | wc -l | tr -d ' '
      ;;
    *)
      echo 0
      ;;
  esac
}

stage_fail_count() {
  local stage="$1"
  local fail_log="$OUTPUT/0_failed_samples.tsv"
  if [ ! -f "$fail_log" ]; then
    echo 0
    return
  fi
  awk -F'\t' -v s="$stage" 'NR>1 && $2==s {seen[$1]=1} END{for (k in seen) n++; print n+0}' "$fail_log"
}

stage_status() {
  local stage="$1"
  local total="$2"
  local done
  local fail
  local status
  done="$(stage_done_count "$stage")"
  fail="$(stage_fail_count "$stage")"

  is_nonneg_int "$total" || total=0
  is_nonneg_int "$done" || done=0
  is_nonneg_int "$fail" || fail=0

  if [ "$total" -gt 0 ] && [ $((done + fail)) -ge "$total" ]; then
    status="complete"
  else
    status="${done}/${total}"
  fi

  if [ "$fail" -gt 0 ]; then
    printf "%s (fail %s)" "$status" "$fail"
  else
    printf "%s" "$status"
  fi
}

bad_fastq_count() {
  local bad_log
  bad_log="$(dirname "$INPUT")/0_bad_fastqs/moved_bad_fastqs.tsv"
  if [ -f "$bad_log" ]; then
    awk 'NR>1{n++} END{print n+0}' "$bad_log"
  else
    echo 0
  fi
}

valid_fastq_count() {
  local marker
  local n
  marker="$(dirname "$INPUT")/0_bad_fastqs/DONE.txt"
  if [ -f "$marker" ]; then
    n="$(awk -F'=' '/^valid_fastq=/{print $2}' "$marker" | tail -n 1)"
    is_nonneg_int "$n" || n=0
    echo "$n"
  else
    n="$(find "$INPUT" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | wc -l | tr -d ' ')"
    is_nonneg_int "$n" || n=0
    echo "$n"
  fi
}

progress_snapshot() {
  local stage="$1"
  local now
  local total
  local bad
  local s_qc s_quast s_autocycler s_medaka s_checkm2 s_coverage s_kraken2 s_bakta
  now="$(date '+%H:%M:%S')"
  total="$(valid_fastq_count)"
  bad="$(bad_fastq_count)"
  s_qc="$(stage_status qc "$total")"
  s_quast="$(stage_status quast "$total")"
  s_autocycler="$(stage_status autocycler "$total")"
  s_medaka="$(stage_status medaka "$total")"
  s_checkm2="$(stage_status checkm2 "$total")"
  s_coverage="$(stage_status coverage "$total")"
  s_kraken2="$(stage_status kraken2 "$total")"
  s_bakta="$(stage_status bakta "$total")"
  printf "[Go_longWGS][Progress][%s][%s] total=%s bad=%s 1_QC=%s 2_quast=%s 3_autocycler=%s 4_medaka=%s 5_checkm2=%s 6_coverage=%s 7a_kraken2=%s 7_bakta=%s\n" \
    "$stage" "$now" "$total" "$bad" "$s_qc" "$s_quast" "$s_autocycler" "$s_medaka" "$s_checkm2" "$s_coverage" "$s_kraken2" "$s_bakta"
}

start_progress_monitor() {
  local stage="$1"
  local current=""
  local current_key=""
  local last=""
  local last_key=""
  if [ "$DRYRUN" -eq 1 ]; then
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

stop_progress_monitor() {
  if [ -n "${PROGRESS_PID:-}" ]; then
    kill "$PROGRESS_PID" >/dev/null 2>&1 || true
    wait "$PROGRESS_PID" 2>/dev/null || true
    PROGRESS_PID=""
  fi
}

SNAKEFILE_NAME="Go_longWGS_V6_1_docker.smk"
if [ ! -f "$PIPELINE_DIR/$SNAKEFILE_NAME" ]; then
  echo "[Go_longWGS][FATAL] Snakefile not found: $PIPELINE_DIR/$SNAKEFILE_NAME"
  exit 1
fi
echo "[Go_longWGS] Using Snakefile: $PIPELINE_DIR/$SNAKEFILE_NAME"
mkdir -p "$DB/medaka_models"

CMD="snakemake --snakefile /pipeline/$SNAKEFILE_NAME \
  --config Fastq_DIRS=\"$INPUT\" output=\"$OUTPUT\" database=\"/db\" \
threads=$DEFAULT_THREADS jobs=$DEFAULT_AUTOCYCLER_JOBS do_porechop=$PORECHOP \
qc_threads=$DEFAULT_THREADS medaka_threads=$DEFAULT_THREADS quast_threads=$DEFAULT_THREADS \
checkm2_threads=$DEFAULT_THREADS cov_threads=$DEFAULT_THREADS bakta_threads=$DEFAULT_THREADS \
kraken_threads=$DEFAULT_THREADS \
  --cores $DEFAULT_CORES --rerun-incomplete --rerun-triggers mtime"

if [ "$DRYRUN" -eq 1 ]; then
  CMD="$CMD --dry-run"
fi
if [ "$KEEP_GOING" -eq 1 ]; then
  CMD="$CMD --keep-going"
fi

run(){ docker run --rm \
  --network host \
  -u "$(id -u):$(id -g)" \
  -v "$WORKDIR":/work \
  -v "$DB":/db \
  -v "$PIPELINE_DIR":/pipeline:ro \
  -w /work \
  -e XDG_CACHE_HOME=/tmp \
  -e PLASSEMBLER_DB="$PLASSEMBLER_DB_PATH" \
  -e MEDAKA_DATA="$MEDAKA_DATA_PATH" \
  longwgs $1; }

on_interrupt() {
  echo "[Go_longWGS] Interrupt received. Stopping running jobs..."
  stop_progress_monitor || true
  pkill -TERM -P $$ >/dev/null 2>&1 || true
  sleep 1
  pkill -KILL -P $$ >/dev/null 2>&1 || true
  exit 130
}
trap on_interrupt INT TERM

set +e
start_progress_monitor "main"
run "$CMD" 2>&1 | tee longwgs.log
rc=${PIPESTATUS[0]}
stop_progress_monitor
set -e

if [ $rc -ne 0 ] && grep -qiE "lock|unlock|LockException|cannot be locked" longwgs.log; then
  echo "[Go_longWGS] Detected lock issue -> running --unlock then retry..."
  run "$CMD --unlock"
  
  set +e
  start_progress_monitor "retry"
  run "$CMD" 2>&1 | tee -a longwgs.log
  rc=${PIPESTATUS[0]}
  stop_progress_monitor
  set -e
fi

# -----------------------------
# Bandage post-processing (run for available GFA files)
# -----------------------------
if [ "$DRYRUN" -eq 0 ]; then
  BANDAGE_DIR="$OUTPUT/8_Bandage_image"

  if [ ! -d "$OUTPUT/3_autocycler" ]; then
    echo "[Go_longWGS][Bandage] No autocycler output directory found. Skipping."
  elif ! command -v Bandage >/dev/null 2>&1; then
    echo "[Go_longWGS][WARN] Host Bandage not found in PATH. Skipping Bandage image generation."
  else
    shopt -s nullglob
    gfas=("$OUTPUT"/3_autocycler/*/autocycler_out/consensus_assembly.gfa)
    shopt -u nullglob

    if [ "${#gfas[@]}" -eq 0 ]; then
      echo "[Go_longWGS][Bandage] No consensus_assembly.gfa files found. Skipping."
    else
      echo "[Go_longWGS] Running Bandage images for available assemblies..."
      mkdir -p "$BANDAGE_DIR"
      for gfa in "${gfas[@]}"; do
        sample=$(basename "$(dirname "$(dirname "$gfa")")")
        gfa_copy="$BANDAGE_DIR/${sample}.consensus_assembly.gfa"
        png="$BANDAGE_DIR/${sample}.consensus_assembly.png"
        if [ ! -s "$gfa_copy" ] || [ "$gfa" -nt "$gfa_copy" ]; then
          cp -f "$gfa" "$gfa_copy"
          echo "  Bandage: copied GFA -> $gfa_copy"
        fi
        if [ -s "$png" ]; then
          echo "  Bandage: skip existing $png"
          continue
        fi
        echo "  Bandage: $gfa -> $png"
        Bandage image "$gfa" "$png" || echo "[Go_longWGS][WARN] Bandage failed for $sample"
      done
    fi
  fi
fi

exit $rc
