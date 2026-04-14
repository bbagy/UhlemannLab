#!/usr/bin/env bash
# =============================================================================
# MEMO (2026-04-11)
# =============================================================================
# 목표: Docker 기반 HuMAnN3 + MetaPhlAn 4 파이프라인 세팅
#
# 현재 상태:
#   - Docker image humann:1.0 빌드 완료
#       - MetaPhlAn 4.1.0 설치
#       - Python wrapper로 metaphlan --version exit 0 보장
#   - Snakemake workflow (Go_Humann_V1.smk) 수정 완료
#       - --db_dir → --bowtie2db 수정
#       - Go_Humann.sh: site-packages 이중 마운트 제거
#
# R1 R2 따로 해도 나중에 합칠수 있음.

# 남은 이슈:
#   - Docker로 하나씩 샘플 돌려보면서 실제 에러 확인 필요
#   - log 파일 (humann3_out/6_logs/<sample>.humann.log) 반드시 확인
#   - 에러 지속 시 humann_temp/<sample>.log 도 확인
#
# 임시 대안:
#   - 급할 때는 예전 conda 방식 (humann3 env + MetaPhlAn3_v3 DB) 사용 가능
#   - DB: /media/uhlemann/core4/DB/humann_db/humann3/MetaPhlAn3_v3
#
# 다음 할 일:
#   1. test/ 샘플 2쌍으로 Docker 테스트 → direct humann stderr 확인
#   2. 에러 해결 후 전체 샘플 실행
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
cd /media/uhlemann/core5/01_MG/20260409_DEAPIM30

inputDIR="DEAPIM30_QC/test"
outDIR="humann3_out_test"
chocophlanDB="/media/uhlemann/core4/DB/humann_db/humann3/chocophlan"
uniref90DB="/media/uhlemann/core4/DB/humann_db/humann3/uniref"
metaphlanDB="/media/uhlemann/core4/DB/humann_db/metaphlan4"
metaphlanIndex="mpa_vJan25_CHOCOPhlAnSGB_202503"
image="humann:1.0"


# ---------------------------------------------------------------------------
# Example: full Go_Humann.sh run
# 아래 형식으로 test 입력을 wrapper로 실행 가능
# ---------------------------------------------------------------------------

cat <<'EOF'
Go_Humann.sh \
   -i DEAPIM30_QC/test \
   -o humann3_out_test \
   -n /media/uhlemann/core4/DB/humann_db/humann3/chocophlan \
   -p /media/uhlemann/core4/DB/humann_db/humann3/uniref \
   -b /media/uhlemann/core4/DB/humann_db/metaphlan4_vJun23 \
   -I mpa_vJun23_CHOCOPhlAnSGB_202307 \
   -c 4 \
   -j 1 \
   -t 2 \
   -K
EOF

# 실제 실행이 필요하면 아래 블록을 사용
#
# Go_Humann.sh \
#   -i "$inputDIR" \
#   -o "$outDIR" \
#   -n "$chocophlanDB" \
#   -p "$uniref90DB" \
#   -b "$metaphlanDB" \
#   -I "$metaphlanIndex" \
#   -c 4 \
#   -j 1 \
#   -t 3 \
#   -K

# ---------------------------------------------------------------------------
# Direct docker run for one sample
# Snakemake/Go_Humann.sh를 거치지 않고 humann만 바로 실행해서 에러 확인
# ---------------------------------------------------------------------------

sample="DPM10002_S202"
mkdir -p humann_direct_debug

zcat \
  "$inputDIR/${sample}_R1_nohuman.fastq.gz" \
  "$inputDIR/${sample}_R2_nohuman.fastq.gz" \
  > "humann_direct_debug/${sample}.fastq"

echo "[check] image=$image"
docker run --rm "$image" which humann
docker run --rm "$image" python -c "import humann; print(humann.__file__)"

docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$(pwd)":/work \
  -v "$chocophlanDB":/db/chocophlan:ro \
  -v "$uniref90DB":/db/uniref:ro \
  -v "$metaphlanDB":/db/metaphlan:ro \
  -w /work \
  "$image" \
  humann \
  --input "humann_direct_debug/${sample}.fastq" \
  --output "humann_direct_debug/${sample}_out" \
  --threads 3 \
  --nucleotide-database /db/chocophlan \
  --protein-database /db/uniref \
  --metaphlan-options "--bowtie2db /db/metaphlan --index $metaphlanIndex" \
  2>&1 | tee "humann_direct_debug/${sample}.direct.log"

tail -n 100 "humann_direct_debug/${sample}.direct.log" || true
