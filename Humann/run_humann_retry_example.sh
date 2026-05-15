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


docker run --rm \
     -u "$(id -u):$(id -g)" \
     -v /media/uhlemann/core4/DB/humann_db/metaphlan4_vJun23:/db \
     humann:1.0 \
     metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202307 --bowtie2db /db
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v /media/uhlemann/core4/DB/humann_db/metaphlan4_vJun23:/db \
    humann:1.0 \
    metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202307 --bowtie2db /db



 docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v /media/uhlemann/core4/DB/humann_db/metaphlan4_vJun23:/db \
    humann:1.0 \
    bash -lc 'export PATH=/opt/conda/envs/humann/bin:$PATH; which bowtie2-build; /opt/conda/envs/humann/bin/bowtie2-build --version; /opt/conda/envs/humann/bin/bowtie2-build --large-index -f /db/mpa_vJun23_CHOCOPhlAnSGB_202307.fna /db/mpa_vJun23_CHOCOPhlAnSGB_202307'


# ---------------------------------------------------------------------------
# Server-side diagnostics for MetaPhlAn DB / bowtie2-build failure
# ---------------------------------------------------------------------------
#
# This workflow expects MetaPhlAn DB as:
#   --bowtie2db /db/metaphlan --index <basename>
# i.e. an already prepared DB directory, not a rebuild from .fna during the
# normal Go_Humann.sh path.
#
# The bowtie2-build failure seen on /db/mpa_vJun23_CHOCOPhlAnSGB_202307.fna
# should be diagnosed separately from the main HUMAnN wrapper.
#
# 1. Check image toolchain versions inside the exact container.
#
cat <<'EOF'
docker run --rm humann:1.0 \
  bash -lc 'export PATH=/opt/conda/envs/humann/bin:$PATH; \
    echo "[humann]"; humann --version; \
    echo "[bowtie2-build]"; bowtie2-build --version | head -n 1; \
    echo "[metaphlan]"; metaphlan --version || true; \
    echo "[python modules]"; python -c "import humann, metaphlan; print(humann.__file__); print(metaphlan.__file__)"'
EOF
#
# 2. Confirm mounted DB path, free space, inode usage, and write test.
#
cat <<'EOF'
docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v /media/uhlemann/core4/DB/humann_db/metaphlan4_vJun23:/db \
  humann:1.0 \
  bash -lc 'set -euo pipefail; \
    echo "[mount]"; ls -ld /db; realpath /db; \
    echo "[space]"; df -h /db; df -i /db; \
    echo "[write test]"; touch /db/.codex_write_test && ls -lh /db/.codex_write_test && rm /db/.codex_write_test; \
    echo "[db files]"; ls -lh /db | sed -n "1,20p"'
EOF
#
# 3. If rebuild is unavoidable, write to a dedicated output prefix and keep
#    stdout/stderr in a log file for the exact failure point.
#
cat <<'EOF'
docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v /media/uhlemann/core4/DB/humann_db/metaphlan4_vJun23:/db \
  humann:1.0 \
  bash -lc 'set -euo pipefail; \
    export PATH=/opt/conda/envs/humann/bin:$PATH; \
    mkdir -p /db/tmp_bt2; \
    export TMPDIR=/db/tmp_bt2; \
    bowtie2-build --large-index -f \
      /db/mpa_vJun23_CHOCOPhlAnSGB_202307.fna \
      /db/mpa_vJun23_CHOCOPhlAnSGB_202307 \
      2>&1 | tee /db/bowtie2-build.mpa_vJun23.log'
EOF
#
# 4. For the main pipeline, prefer a MetaPhlAn DB directory that already
#    contains the matching .pkl and .bt2/.bt2l files for the selected index.


  docker run --rm humann:1.0 \
    bash -lc 'export PATH=/opt/conda/envs/humann/bin:$PATH;bowtie2-build --version | head -n 1'

 grep -n "bowtie2" Dockerfile
  docker run --rm humann:1.0 bash -lc 'export PATH=/opt/conda/envs/humann/bin:$PATH; conda list bowtie2 || micromamba list -n humann bowtie2'


    docker run --rm humann:1.0 \
    bash -lc 'export PATH=/opt/conda/envs/humann/bin:$PATH; \
      echo "[which]"; which bowtie2-build; \
      echo "[ls]"; ls -l "$(which bowtie2-build)"; \
      echo "[conda list]"; micromamba list -n humann | grep -E "^bowtie2|^metaphlan|^humann"; \
      echo "[version]"; bowtie2-build --version | head -n 3'



  docker run --rm humann:1.0 \
    bash -lc 'export PATH=/opt/conda/envs/humann/bin:$PATH; \
      echo "[which]"; which bowtie2-build; \
      echo "[type]"; type -a bowtie2-build; \
      echo "[ls]"; ls -l "$(which bowtie2-build)"; \
      echo "[find]"; find / -name bowtie2-build 2>/dev/null'

  docker run --rm humann:1.0 \
    bash -lc '/opt/conda/envs/humann/bin/bowtie2-build --version | head -n 3'


  docker run --rm humann:1.0 \
    bash -lc 'export PATH=/opt/conda/envs/humann/bin:$PATH;export LD_LIBRARY_PATH=/opt/conda/envs/humann/lib:$LD_LIBRARY_PATH; bowtie2-build --version | head -n 3'


# =============================================================================
# MEMO (2026-04-16) next retry strategy
# =============================================================================
# Now that bowtie2-build runs correctly as 2.5.5 inside humann:1.0 when:
#   export PATH=/opt/conda/envs/humann/bin:$PATH
#   export LD_LIBRARY_PATH=/opt/conda/envs/humann/lib:$LD_LIBRARY_PATH
#
# the next reasonable retry is NOT manual rebuild first.
# Try official MetaPhlAn install/download first, because:
#   - previous manual rebuild attempts were started mainly because bowtie2
#     version/runtime was broken inside the container
#   - that runtime issue is now fixed
#   - main Go_Humann.sh workflow expects a prepared MetaPhlAn DB directory
#     plus index basename anyway
#
# Recommended order:
#   1. Retry `metaphlan --install --index ... --bowtie2db /db`
#   2. If install succeeds, use that DB directly with:
#        -b /media/uhlemann/core4/DB/humann_db/metaphlan4_vJun23
#        -I mpa_vJun23_CHOCOPhlAnSGB_202307
#   3. Only if official install still fails, fall back to manual bowtie2-build
#
# Official MetaPhlAn install retry:
cat <<'EOF'
docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v /media/uhlemann/core4/DB/humann_db/metaphlan4_vJun23:/db \
  humann:1.0 \
  bash -lc 'set -euo pipefail; \
    export PATH=/opt/conda/envs/humann/bin:$PATH; \
    export LD_LIBRARY_PATH=/opt/conda/envs/humann/lib:$LD_LIBRARY_PATH; \
    metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202307 --bowtie2db /db \
      2>&1 | tee /db/metaphlan_install.mpa_vJun23.log'
EOF
#
# Post-install check:
cat <<'EOF'
docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v /media/uhlemann/core4/DB/humann_db/metaphlan4_vJun23:/db \
  humann:1.0 \
  bash -lc 'set -euo pipefail; \
    ls -lh /db/mpa_vJun23_CHOCOPhlAnSGB_202307*'
EOF

# =============================================================================
# MEMO (2026-04-27) 디스크 정리 및 root 용량 조사
# =============================================================================
# 상황:
#   - root (/) 621GB used / 747GB (88%)
#   - Docker image prune 0B 회수 (<none> 이미지들이 humann:1.0과 레이어 공유)
#   - diamond 임시파일 생성 실패 원인 조사 중
#
# 1. 오래된 Docker 이미지 삭제 (OrthoVenn3 + islandpath, ~11GB)
docker rmi \
  lufang0411/orthovenn3-api:latest \
  leeoluo/orthovenn3-front:latest \
  lufang0411/orthovenn3-mysql:latest \
  brinkmanlab/islandpath:1.0.0 \
  quay.io/biocontainers/islandpath:1.0.6--hdfd78af_0

# 2. 멈춘 컨테이너 정리
docker container prune -f

# 3. root 파티션 용량 범인 찾기
df -h /
du -sh /home/* /var/log /var/lib /opt /tmp /root /snap 2>/dev/null | sort -rh | head -15
