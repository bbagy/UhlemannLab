#############################################
# Go_MAGs_Assembly_V1.smk
# 20260406
# Heekuk Park
#
# Pipeline:
#   QC nohuman FASTQ -> MEGAHIT -> contig filtering -> bowtie2 mapping
#   -> coverage -> MetaBAT2/MaxBin2/CONCOCT -> DAS_Tool -> CheckM v1
#
# Input:
#   FASTQ_DIR/
#     <sample>_R1_nohuman.fastq.gz
#     <sample>_R2_nohuman.fastq.gz
#
# OUTDIR/
#   1_MEGAHIT/         ({sample}/final.contigs.fa)
#   2_Contigs/         ({sample}_contigs.fa, {sample}_filtered_contigs.fa)
#   3_Binning/         (MetaBAT2, MaxBin2, CONCOCT outputs)
#   4_Mapping/         (bowtie2 contig indexes)
#   5_Coverage/        ({sample}_depth.txt, {sample}_depth_fixed.txt)
#   6_Mapped_Reads/    ({sample}.bam)
#   7_DAS_tool_out/    (DAS_Tool bin integration)
#   8_checkM_summary/  ({sample}_checkm_summary.csv)
#   9_Final_MAGs/      (flat final MAG FASTA files for annotation input)
#   assembly_summary.csv
#   checkm_summary_combined.csv
#
# Example run:
# fastq_dir="mags_qc_out/host_filtered_fastq"
# output_dir="mags_assembly_out"
#
# snakemake --snakefile /home/uhlemann/heekuk_path/MAGs/workflow/Go_MAGs_Assembly_V1.smk \
#   --config fastq_dir="$fastq_dir" output_dir="$output_dir" \
#   megahit_threads=24 megahit_memory=128000 binning_tools="concoct,metabat2,maxbin2" \
#   --cores 24 --jobs 4 --latency-wait 60 --rerun-incomplete --keep-going
#############################################

import glob
import subprocess
import os
import re
import pandas as pd

# Config 값 가져오기
FASTQ_DIR = config["fastq_dir"]
OUTPUT_DIR = config["output_dir"]
MEGAHIT_THREADS = int(config.get("megahit_threads", 8))
MEGAHIT_MEMORY  = int(config.get("megahit_memory", 128000))
CHECKM_DATA_DIR = str(config.get("checkm_data_dir", "")).strip()
binning_tools_str = config.get("binning_tools", "concoct,metabat2,maxbin2")
BINNING_TOOLS = [x.strip() for x in binning_tools_str.split(",") if x.strip()]


# Output 디렉토리 생성 (순차 번호 부여)
os.makedirs(f"{OUTPUT_DIR}/1_MEGAHIT", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/2_Contigs", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/3_Binning", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/4_Mapping", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/5_Coverage", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/6_Mapped_Reads", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/7_DAS_tool_out", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/8_checkM_summary", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/logs", exist_ok=True)

# 샘플별 FASTQ 파일 검색
fastq_files = sorted(glob.glob(f"{FASTQ_DIR}/*_R1_nohuman.fastq.gz"))
# samples = list(set(["_".join(os.path.basename(f).split("_")[:1]) for f in fastq_files]))

samples = list(set([
    re.sub(r'_R1_nohuman\.fastq\.gz$', '', os.path.basename(f))
    for f in fastq_files
]))


###############################################################################
# rule all
###############################################################################
rule all:
    input:
        expand(f"{OUTPUT_DIR}/2_Contigs/{{sample}}_contigs.fa", sample=samples),  # ✅ 추가
        expand(f"{OUTPUT_DIR}/2_Contigs/{{sample}}_filtered_contigs.fa", sample=samples),
        expand(f"{OUTPUT_DIR}/5_Coverage/{{sample}}_depth.txt", sample=samples),
        expand(f"{OUTPUT_DIR}/5_Coverage/{{sample}}_depth_fixed.txt", sample=samples),
        expand(f"{OUTPUT_DIR}/6_Mapped_Reads/{{sample}}.bam", sample=samples),
        expand(f"{OUTPUT_DIR}/8_checkM_summary/{{sample}}_checkm_summary.csv", sample=samples),
        f"{OUTPUT_DIR}/assembly_summary.csv",
        f"{OUTPUT_DIR}/checkm_summary_combined.csv",
        f"{OUTPUT_DIR}/9_Final_MAGs/DONE.txt"
        


        
        

###############################################################################
# 모니터링 규칙들 (swap, cache etc.)
###############################################################################
rule monitor_swap:
    output:
        # 플래그 파일(메모리 감시 완료 표시)
        flag = f"{OUTPUT_DIR}/logs/swap_monitor.flag",
        # 로그 파일(메모리 감시 결과 기록)
        log = f"{OUTPUT_DIR}/logs/swap_monitor.log"
    shell:
        """
        echo "🚀 Checking Swap & CPU usage at $(date '+%Y-%m-%d %H:%M:%S')" >> {output.log}

        # ✅ 현재 Swap 사용량(%) 확인 (Swap이 사용 중인지 확인)
        SWAP_TOTAL=$(free | awk '/Swap:/ {{print $2}}')
        SWAP_USED=$(free | awk '/Swap:/ {{print $3}}')

        if [ "$SWAP_TOTAL" -gt 0 ]; then
            SWAP_PERCENT=$(echo "scale=2; ($SWAP_USED/$SWAP_TOTAL)*100" | bc)
        else
            SWAP_PERCENT=0
        fi

        # ✅ 현재 사용 가능한 RAM 용량 (단위: KB)
        MEM_AVAILABLE=$(free | awk '/Mem:/ {{print $7}}')

        # ✅ 현재 CPU 사용률 (%) (user + system)
        CPU_LOAD=$(top -bn1 | grep "Cpu(s)" | awk '{{print $2 + $4}}')

        # ✅ Snakemake가 실행 중인 프로세스 개수 확인
        #   (주의: grep \"snakemake\"에 의해 grep 본인도 잡힐 수 있으므로, 실제 환경에 맞게 수정 가능)
        SNAKEMAKE_RUNNING=$(ps -ef | grep "snakemake" | wc -l)

        echo "🔍 Swap: $SWAP_PERCENT%, Free RAM: $MEM_AVAILABLE KB, CPU Load: $CPU_LOAD%, Running Jobs: $SNAKEMAKE_RUNNING" >> {output.log}

        # 🔥 Swap 초기화 실행 조건 (안전한 상태에서만 실행)
        if [ "$(echo "$SWAP_PERCENT >= 90" | bc -l)" -eq 1 ] && \
           [ "$MEM_AVAILABLE" -gt 5000000 ] && \
           [ "$SNAKEMAKE_RUNNING" -lt 10 ] && \
           [ "$(echo "$CPU_LOAD < 90" | bc -l)" -eq 1 ]; then
            echo "⚠️ High Swap usage ($SWAP_PERCENT%) detected. Resetting swap..." >> {output.log}

            echo "⚠️ Swap reset skipped inside Docker at $(date '+%Y-%m-%d %H:%M:%S')" >> {output.log}

            # 🕒 CPU 부하 방지를 위해 10초 대기
            sleep 10
        else
            echo "✅ Swap usage is $SWAP_PERCENT%. No reset needed." >> {output.log}
        fi

        # 🔍 Snakemake가 실행 완료를 감지할 수 있도록 플래그 파일 생성
        touch {output.flag}
        """

rule manage_swap:
    output:
        log = f"{OUTPUT_DIR}/logs/swap_management.log"  # ✅ Swap 초기화 로그 기록 파일
    shell:
        """
        echo "🚀 Managing swap space..." >> {output.log}

        # 🔍 Swap 초기화 실행 시간을 기록하는 파일 경로
        LAST_RESET_FILE="{OUTPUT_DIR}/last_swap_reset.txt"

        # 🔍 이전 Swap 초기화 시간을 확인하여 일정 시간(30분) 이내이면 실행하지 않음
        if [ -f "$LAST_RESET_FILE" ]; then
            LAST_RESET=$(cat "$LAST_RESET_FILE")  # ✅ 마지막 Swap 초기화 시간(초 단위) 가져오기
            NOW=$(date +%s)  # ✅ 현재 시간(초 단위)
            DIFF=$((NOW - LAST_RESET))  # ✅ 마지막 초기화 이후 경과 시간 계산

            # ⏳ 만약 30분(1800초) 이내라면 초기화 스킵
            if [ "$DIFF" -lt 1800 ]; then
                echo "⏳ Swap reset skipped (Last reset was $DIFF seconds ago)" >> {output.log}
                exit 0
            fi
        fi

        # ✅ 현재 RAM 사용률(%) 확인 (수정된 버전)
        RAM_TOTAL=$(free | awk '/Mem:/ {{print $2}}')
        RAM_USED=$(free | awk '/Mem:/ {{print $3}}')
        RAM_PERCENT=$(echo "scale=2; ($RAM_USED/$RAM_TOTAL)*100" | bc)

        # ✅ 현재 I/O 부하 확인 (디스크 입출력 사용량)
        IO_USAGE=$(iostat -d 1 2 | awk 'NR>3 {{sum+=$2}} END {{print sum}}')

        # ⚠️ RAM 사용률이 95% 이상이면 Swap 초기화를 하지 않음 (오히려 RAM 부족으로 크래시 발생 가능)
        if (( $(echo "$RAM_PERCENT >= 95" | bc -l) )); then
            echo "⚠️ High RAM usage ($RAM_PERCENT%). Skipping swap reset." >> {output.log}
            exit 0
        fi

        # ⚠️ 디스크 I/O 부하가 50% 이상이면 Swap 초기화를 하지 않음 (I/O 부하가 높은 상태에서 실행하면 성능 저하 가능)
        if (( $(echo "$IO_USAGE >= 50" | bc -l) )); then
            echo "⚠️ High I/O load ($IO_USAGE%). Skipping swap reset." >> {output.log}
            exit 0
        fi

        # Docker-safe no-op: do not modify host swap from inside the container.
        echo "⚠️ Swap reset skipped inside Docker" >> {output.log}

        # 🔍 초기화 시간을 기록하여 30분 이내 중복 실행 방지
        echo "$(date +%s)" > "$LAST_RESET_FILE"

        echo "✅ Swap has been reset (RAM: $RAM_PERCENT%, I/O: $IO_USAGE%)" >> {output.log}
        """
      
        
rule monitor_cache:
    output:
        flag = f"{OUTPUT_DIR}/logs/cache_monitor.flag"  # ✅ 캐시 정리 완료 여부를 기록하는 플래그 파일
    shell:
        """
        echo "🚀 Checking Cache Usage..."

        # ✅ 현재 캐시 사용량(GB)
        CACHE_USED=$(free -g | awk '/Mem/ {{print $6}}')

        # ✅ 사용 가능한 메모리(GB)
        MEM_AVAILABLE=$(free -g | awk '/Mem/ {{print $7}}')

        # ✅ 현재 CPU 부하 확인
        CPU_LOAD=$(top -bn1 | grep "Cpu(s)" | awk '{{print $2 + $4}}')

        # ✅ 현재 I/O 부하 확인 (디스크 입출력 사용량)
        IO_USAGE=$(iostat -d 1 2 | awk 'NR>3 {{sum+=$2}} END {{print sum}}')

        # ✅ Snakemake 실행 중인 프로세스 개수 확인
        SNAKEMAKE_RUNNING=$(ps -ef | grep "snakemake" | wc -l)

        # 🔍 최근 캐시 정리 시간 확인
        LAST_CLEAR_FILE="{OUTPUT_DIR}/last_cache_clear.txt"

        if [ -f "$LAST_CLEAR_FILE" ]; then
            LAST_CLEAR=$(cat "$LAST_CLEAR_FILE")
            NOW=$(date +%s)
            DIFF=$((NOW - LAST_CLEAR))

            # ⏳ 30분(1800초) 이내라면 캐시 정리 건너뜀
            if [ "$DIFF" -lt 1800 ]; then
                echo "⏳ Cache reset skipped (Last reset was $DIFF seconds ago)"
                exit 0
            fi
        fi

        # 🔥 캐시 정리 실행 조건 (안전한 상태에서만 실행)
        if [ "$CACHE_USED" -ge 50 ] && [ "$MEM_AVAILABLE" -lt 10 ] && \
           [ "$IO_USAGE" -lt 50 ] && [ "$SNAKEMAKE_RUNNING" -lt 10 ] && \
           [ "$(echo "$CPU_LOAD < 90" | bc)" -eq 1 ]; then
            echo "⚠️ High Cache usage ($CACHE_USED GB). Clearing cache..."

            echo "⚠️ Cache clearing skipped inside Docker."

            # 🔍 캐시 정리 실행 시간을 기록하여 중복 실행 방지
            echo "$(date +%s)" > "$LAST_CLEAR_FILE"

            # 🕒 CPU 부하 방지를 위해 10초 대기
            sleep 10
        else
            echo "✅ Cache usage is $CACHE_USED GB. No clearing needed."
        fi

        # 🔍 Snakemake가 실행 완료를 감지할 수 있도록 플래그 파일 생성
        touch {output.flag}
        """


###############################################################################
# 1) MEGAHIT assembly
###############################################################################
rule megahit_assembly:
    input:
        R1 = f"{FASTQ_DIR}/{{sample}}_R1_nohuman.fastq.gz",
        R2 = f"{FASTQ_DIR}/{{sample}}_R2_nohuman.fastq.gz"
    output:
        contigs = f"{OUTPUT_DIR}/1_MEGAHIT/{{sample}}/final.contigs.fa"
    shell:
        """
        echo "Checking if MEGAHIT output directory exists: {OUTPUT_DIR}/1_MEGAHIT/{wildcards.sample}"
        
        if [ -d "{OUTPUT_DIR}/1_MEGAHIT/{wildcards.sample}" ]; then
            if [ -f "{OUTPUT_DIR}/1_MEGAHIT/{wildcards.sample}/options.json" ]; then
                echo "Previous MEGAHIT run detected, resuming with --continue..."
                megahit --continue -o {OUTPUT_DIR}/1_MEGAHIT/{wildcards.sample} --num-cpu-threads {MEGAHIT_THREADS} --memory {MEGAHIT_MEMORY}
            else
                echo "Incomplete MEGAHIT run detected, removing and restarting..."
                rm -rf {OUTPUT_DIR}/1_MEGAHIT/{wildcards.sample}
                megahit -1 {input.R1} -2 {input.R2} \
                        -o {OUTPUT_DIR}/1_MEGAHIT/{wildcards.sample} --num-cpu-threads {MEGAHIT_THREADS} --memory {MEGAHIT_MEMORY}
            fi
        else
            echo "No previous MEGAHIT run detected, starting fresh..."
            megahit -1 {input.R1} -2 {input.R2} \
                    -o {OUTPUT_DIR}/1_MEGAHIT/{wildcards.sample} --num-cpu-threads {MEGAHIT_THREADS} --memory {MEGAHIT_MEMORY}
        fi

        # ✅ MEGAHIT 완료 후 컨티그 개수 확인
        num_contigs=$(grep -c "^>" {output.contigs} || echo 0)

        if [ "$num_contigs" -lt 2 ]; then
            echo "⚠️ Skipping MEGAHIT for {wildcards.sample} (only $num_contigs contig(s) found)"
            
            # ✅ `skip.log` 업데이트 (flock 사용하여 동시 접근 방지)
            (
                flock -x 200
                mkdir -p {OUTPUT_DIR}
                echo "{wildcards.sample} - MEGAHIT skipped (only $num_contigs contigs found)" >> {OUTPUT_DIR}/skip.log
                touch {OUTPUT_DIR}/skip.log
            ) 200> {OUTPUT_DIR}/skip.log.lock
            
            exit 0
        fi

        echo "✅ MEGAHIT completed successfully for {wildcards.sample} with $num_contigs contigs."
        """

###############################################################################
# 2) Copy contigs -> 2_Contigs
###############################################################################
rule copy_contigs:
    input:
        contigs = f"{OUTPUT_DIR}/1_MEGAHIT/{{sample}}/final.contigs.fa"
    output:
        copied_contigs = f"{OUTPUT_DIR}/2_Contigs/{{sample}}_contigs.fa"
    shell:
        """
        echo "🔍 Checking if final.contigs.fa exists for sample: {wildcards.sample}"
        
        # ✅ 파일 존재 여부 확인
        while [ ! -f {input.contigs} ]; do
            echo "⏳ Waiting for {input.contigs} to be generated..."
            sleep 10  # 10초 대기 후 다시 체크
        done

        if [ ! -s {input.contigs} ]; then
            echo "⚠️ WARNING: {input.contigs} is empty. Skipping copy..."
            touch {output.copied_contigs}
            exit 0
        fi

        echo "📌 Copying {input.contigs} to {output.copied_contigs}"
        cp {input.contigs} {output.copied_contigs}

        # ✅ 복사 확인
        if [ ! -s {output.copied_contigs} ]; then
            echo "❌ ERROR: {output.copied_contigs} is empty after copy!"
            exit 1
        fi

        echo "✅ Successfully copied: {output.copied_contigs}"
        """

###############################################################################
# 3) Filter contigs >=1kb
###############################################################################
rule filter_contigs:
    input:
        contigs = f"{OUTPUT_DIR}/2_Contigs/{{sample}}_contigs.fa"
    output:
        filtered_contigs = f"{OUTPUT_DIR}/2_Contigs/{{sample}}_filtered_contigs.fa"
    shell:
        """
        awk 'BEGIN {{seq=""; header=""}}
             /^>/ {{if (length(seq) >= 1000) print header "\\n" seq; header=$0; seq=""}}
             !/^>/ {{seq=seq""$0}}
             END {{if (length(seq) >= 1000) print header "\\n" seq}}' {input.contigs} > {output.filtered_contigs}

        # ✅ filtered_contigs.fa 파일 크기 확인 (1MB 이하인지 체크)
        FILE_SIZE=$(stat -c%s "{output.filtered_contigs}")

        if [ "$FILE_SIZE" -le 1048576 ]; then
            mkdir -p {OUTPUT_DIR}
            echo "⚠️ {wildcards.sample} - Skipping (Filtered contigs < 1MB)" >> {OUTPUT_DIR}/skip.log
            echo "⚠️ {wildcards.sample} - Skipped due to small contigs file size"
            touch {output.filtered_contigs}  # 빈 파일 생성 (Snakemake 오류 방지)
            exit 0
        fi
        """

###############################################################################
# 4) Bowtie2 mapping -> BAM
###############################################################################
rule map_reads:
    input:
        R1 = f"{FASTQ_DIR}/{{sample}}_R1_nohuman.fastq.gz",
        R2 = f"{FASTQ_DIR}/{{sample}}_R2_nohuman.fastq.gz",
        contigs = f"{OUTPUT_DIR}/2_Contigs/{{sample}}_filtered_contigs.fa"
    output:
        bam = f"{OUTPUT_DIR}/6_Mapped_Reads/{{sample}}.bam"
    run:
        import os

        skip_file = f"{OUTPUT_DIR}/skip.log"
        skip_samples = set()

        if os.path.exists(skip_file):
            with open(skip_file) as f:
                skip_samples = {line.strip().split(" - ")[0] for line in f}

        if wildcards.sample in skip_samples:
            print(f"⚠️ Skipping mapping: {wildcards.sample} is listed in skip.log")
            os.makedirs(os.path.dirname(output.bam), exist_ok=True)
            open(output.bam, "w").close()
            return

        shell("""
        echo "🚀 Mapping reads to contigs with Bowtie2 for {wildcards.sample}..."

        # ✅ contigs 파일 존재/비어있는지 확인
        if [ ! -s {input.contigs} ]; then
            echo "⚠️ Skipping mapping: contigs file is missing or empty for {wildcards.sample}"
            echo "{wildcards.sample} - Skipped mapping (empty contigs)" >> {OUTPUT_DIR}/skip.log
            touch {output.bam}
            exit 0
        fi

        mkdir -p {OUTPUT_DIR}/4_Mapping

        if [ ! -f {OUTPUT_DIR}/4_Mapping/{wildcards.sample}.1.bt2 ]; then
            bowtie2-build {input.contigs} {OUTPUT_DIR}/4_Mapping/{wildcards.sample}
        fi

        bowtie2 -x {OUTPUT_DIR}/4_Mapping/{wildcards.sample} -1 {input.R1} -2 {input.R2} | \
        samtools view -Sb - > {OUTPUT_DIR}/6_Mapped_Reads/{wildcards.sample}.unsorted.bam

        echo "🔄 Checking if BAM file is sorted..."
        if samtools view -H {OUTPUT_DIR}/6_Mapped_Reads/{wildcards.sample}.unsorted.bam | grep -q '@HD.*SO:coordinate'; then
            echo "✅ BAM file is already sorted!"
            mv {OUTPUT_DIR}/6_Mapped_Reads/{wildcards.sample}.unsorted.bam {output.bam}
        else
            echo "🔄 Sorting BAM file..."
            samtools sort -o {output.bam} {OUTPUT_DIR}/6_Mapped_Reads/{wildcards.sample}.unsorted.bam
            rm {OUTPUT_DIR}/6_Mapped_Reads/{wildcards.sample}.unsorted.bam
        fi

        # ✅ BAM 파일이 비어있는 경우 추가적으로 skip 처리
        if [ ! -s {output.bam} ]; then
            mkdir -p {OUTPUT_DIR}
            echo "⚠️ {wildcards.sample} - Skipped (Empty BAM file)" >> {OUTPUT_DIR}/skip.log
            exit 0
        fi
        """)


###############################################################################
# 5) jgi_summarize_bam_contig_depths -> coverage
###############################################################################
rule generate_abundance:
    input:
        bam = f"{OUTPUT_DIR}/6_Mapped_Reads/{{sample}}.bam"
    output:
        depth_raw = f"{OUTPUT_DIR}/5_Coverage/{{sample}}_depth.txt",
        depth_fixed = f"{OUTPUT_DIR}/5_Coverage/{{sample}}_depth_fixed.txt"
    log:
        f"{OUTPUT_DIR}/6_Mapped_Reads/{{sample}}_generate_abundance.log"
    shell:
        """
        echo "🚀 Checking BAM file for {wildcards.sample}..." | tee -a {log}

        # 🔍 BAM 파일 확인 (파일이 비어 있으면 스킵)
        if [ ! -s {input.bam} ]; then
            echo "⚠️ BAM file for {wildcards.sample} is empty! Skipping." | tee -a {log}
            mkdir -p {OUTPUT_DIR}
            echo "{wildcards.sample} - Skipped (Empty BAM file)" >> {OUTPUT_DIR}/skip.log
            touch {output.depth_raw}
            touch {output.depth_fixed}
            exit 0
        fi

        # 🔍 BAM 파일이 사용 중인지 확인
        if lsof {input.bam} >/dev/null 2>&1; then
            echo "❌ ERROR: BAM file {input.bam} is in use by another process!" | tee -a {log}
            exit 1
        fi

        # 🔄 기존 임시 파일 삭제 (정렬 실패 방지)
        echo "🔄 Cleaning up temporary BAM files..." | tee -a {log}
        rm -f {OUTPUT_DIR}/6_Mapped_Reads/{wildcards.sample}.bam.tmp.*.bam

        # 🔄 BAM 파일 정렬 여부 확인 및 정렬
        if samtools view -H {input.bam} | grep -q '@HD.*SO:coordinate'; then
            echo "✅ BAM file is already sorted." | tee -a {log}
        else
            echo "🔄 Sorting BAM file..." | tee -a {log}
            samtools sort -T {OUTPUT_DIR}/6_Mapped_Reads/{wildcards.sample}.bam.tmp -o {OUTPUT_DIR}/6_Mapped_Reads/{wildcards.sample}.sorted.bam {input.bam}

            if [ ! -s {OUTPUT_DIR}/6_Mapped_Reads/{wildcards.sample}.sorted.bam ]; then
                echo "❌ ERROR: Sorting failed for {wildcards.sample}!" | tee -a {log}
                exit 1
            fi

            mv {OUTPUT_DIR}/6_Mapped_Reads/{wildcards.sample}.sorted.bam {input.bam}
        fi

        # 🔍 Contig 개수 확인 (2개 미만이면 스킵)
        num_contigs=$(awk 'BEGIN {{ count=0 }} /^>/ {{ count++ }} END {{ print count }}' {OUTPUT_DIR}/2_Contigs/{wildcards.sample}_filtered_contigs.fa)

        if [ "$num_contigs" -lt 2 ]; then
            echo "⚠️ Skipping abundance calculation for {wildcards.sample} (only $num_contigs contig(s) found)" | tee -a {log}
            mkdir -p {OUTPUT_DIR}
            echo "{wildcards.sample} - Skipped (only $num_contigs found)" >> {OUTPUT_DIR}/skip.log
            touch {output.depth_raw}
            touch {output.depth_fixed}
            exit 0
        fi

        echo "🚀 Generating raw abundance file for {wildcards.sample}" | tee -a {log}
        jgi_summarize_bam_contig_depths --outputDepth {output.depth_raw} {input.bam} 2>> {log}

        if [ ! -s {output.depth_raw} ]; then
            echo "❌ ERROR: jgi_summarize_bam_contig_depths failed for {wildcards.sample}!" | tee -a {log}
            exit 1
        fi

        echo "🔍 Transforming coverage file for CONCOCT" | tee -a {log}
        awk -v sample="Sample_{wildcards.sample}" 'NR==1 {{print $1, sample}} NR>1 {{print $1, $4}}' OFS="\\t" {output.depth_raw} > {output.depth_fixed}

        if [ ! -s {output.depth_fixed} ]; then
            echo "❌ ERROR: Transformed coverage file is empty! Skipping {wildcards.sample}." | tee -a {log}
            exit 1
        fi

        echo "✅ Fixed coverage file saved: {output.depth_fixed}" | tee -a {log}

        # ✅ Snakemake가 감지하도록 2초 대기
        sleep 2
        """


###############################################################################
# 6) MetaBAT2 / MaxBin2 / CONCOCT (binning)
###############################################################################
# 6️⃣ MetaBAT2 Binning
rule metabat2:
    input:
        contigs = f"{OUTPUT_DIR}/2_Contigs/{{sample}}_filtered_contigs.fa",
        abundance = f"{OUTPUT_DIR}/5_Coverage/{{sample}}_depth.txt"
    output:
        bins = directory(f"{OUTPUT_DIR}/3_Binning/MetaBAT2/{{sample}}")
    shell:
        """
        echo "🚀 Running MetaBAT2 for {wildcards.sample}"

        # Contig 개수 확인
        num_contigs=$(awk '/^>/ {{count+=1}} END {{print count}}' {input.contigs})

        if [ "$num_contigs" -lt 2 ]; then
            echo "⚠️ Skipping MetaBAT2 for {wildcards.sample} (only $num_contigs contig(s) found)"
            echo "{wildcards.sample} - MetaBAT2 skipped (only $num_contigs contig(s) found)" >> {OUTPUT_DIR}/skip.log
            mkdir -p {output.bins}  # 빈 폴더 생성하여 Snakemake 오류 방지
            touch {output.bins}/bin_dummy.txt  # 🚀 MetaBAT2가 필요로 하는 빈 파일 생성
            exit 0
        fi

        # Read Depth 파일이 존재하는 경우 실행
        if [ -s {input.abundance} ]; then
            echo "📌 Using read depth for MetaBAT2"
            metabat2 -i {input.contigs} -a {input.abundance} -o {output.bins}/bin
        else
            echo "⚠️ No read depth file found for {wildcards.sample}. Running MetaBAT2 without read depth."
            mkdir -p {OUTPUT_DIR}
            echo "{wildcards.sample} - MetaBAT2 running without read depth" >> {OUTPUT_DIR}/skip.log
            metabat2 -i {input.contigs} -o {output.bins}/bin
        fi
        """

# 7️⃣ MaxBin2 Binning
rule maxbin2:
    input:
        contigs = f"{OUTPUT_DIR}/2_Contigs/{{sample}}_filtered_contigs.fa",
        abundance = f"{OUTPUT_DIR}/5_Coverage/{{sample}}_depth.txt"
    output:
        bins = directory(f"{OUTPUT_DIR}/3_Binning/MaxBin2/{{sample}}")
    shell:
        """
        echo "🚀 Running MaxBin2 for {wildcards.sample}"

        num_contigs=$(awk 'BEGIN {{count=0}} /^>/ {{count++}} END {{print count}}' {input.contigs})

        if [ "$num_contigs" -lt 2 ]; then
            echo "⚠️ Skipping MaxBin2 for {wildcards.sample} (only $num_contigs contig(s) found)"
            echo "{wildcards.sample} - MaxBin2 skipped (only $num_contigs contig(s) found)" >> {OUTPUT_DIR}/skip.log
            mkdir -p {output.bins}
            touch {output.bins}/bin_dummy.txt
            exit 0
        fi

        OUT_PREFIX={output.bins}/maxbin_out

        if [ -e "$OUT_PREFIX" ]; then
            echo "⚠️ Existing output prefix may interfere. Removing..."
            rm -rf "$OUT_PREFIX"
        fi

        mkdir -p {output.bins}

        if [ -s {input.abundance} ]; then
            echo "📌 Using read depth for MaxBin2"
            run_MaxBin.pl -contig {input.contigs} -out "$OUT_PREFIX" -abund {input.abundance} -thread {MEGAHIT_THREADS} || {{
                echo "⚠️ MaxBin2 failed (possibly due to lack of marker genes)"
                echo "{wildcards.sample} - MaxBin2 skipped (no marker genes or internal failure)" >> {OUTPUT_DIR}/skip.log
                touch {output.bins}/bin_dummy.txt
                exit 0
            }}
        else
            echo "⚠️ No read depth file found. Skipping MaxBin2."
            echo "{wildcards.sample} - MaxBin2 skipped (no read depth file)" >> {OUTPUT_DIR}/skip.log
            touch {output.bins}/bin_dummy.txt
            exit 0
        fi
        """

# 8️⃣ CONCOCT Binning
rule concoct:
    input:
        contigs = f"{OUTPUT_DIR}/2_Contigs/{{sample}}_filtered_contigs.fa",
        coverage = f"{OUTPUT_DIR}/5_Coverage/{{sample}}_depth_fixed.txt"
    output:
        bins = directory(f"{OUTPUT_DIR}/3_Binning/CONCOCT/{{sample}}/fasta_bins")
    shell:
        """
        echo "🚀 Running CONCOCT for {wildcards.sample}"

        # ✅ 컨티그 개수 확인
        num_contigs=$(awk 'BEGIN {{count=0}} /^>/ {{count++}} END {{print count}}' {input.contigs})
        if [ "$num_contigs" -lt 2 ]; then
            echo "⚠️ Skipping CONCOCT for {wildcards.sample} (only $num_contigs found)"
            mkdir -p {OUTPUT_DIR}
            echo "{wildcards.sample} - CONCOCT skipped (only $num_contigs found)" >> {OUTPUT_DIR}/skip.log
            mkdir -p {output.bins}
            touch {output.bins}/.snakemake_success
            exit 0
        fi

        # ✅ Read Depth 파일 존재 확인
        if [ -s {input.coverage} ]; then
            echo "📌 Using transformed read depth for CONCOCT"
            mkdir -p {OUTPUT_DIR}/3_Binning/CONCOCT/{wildcards.sample}

            concoct --composition_file {input.contigs} \
                    --coverage_file {input.coverage} \
                    -b {OUTPUT_DIR}/3_Binning/CONCOCT/{wildcards.sample}/ \
                    --threads {MEGAHIT_THREADS}

            # ✅ `bin_clustering_gt1000.csv` 존재 확인 (이제 bin_ 접두사 포함)
            if [ -f {OUTPUT_DIR}/3_Binning/CONCOCT/{wildcards.sample}/clustering_gt1000.csv ]; then
                echo "🚀 Extracting FASTA bins for {wildcards.sample}"
                mkdir -p {output.bins}

                extract_fasta_bins.py {input.contigs} \
                                      {OUTPUT_DIR}/3_Binning/CONCOCT/{wildcards.sample}/clustering_gt1000.csv \
                                      --output_path {output.bins}

                # ✅ Snakemake가 폴더를 인식할 수 있도록 빈 파일 추가
                touch {output.bins}/.snakemake_success
            else
                echo "⚠️ No clustering_gt1000.csv found. CONCOCT did not generate bin assignments."
                mkdir -p {output.bins}
                touch {output.bins}/.snakemake_success
                exit 0
            fi
        else
            echo "⚠️ No read depth file found. Skipping CONCOCT for {wildcards.sample}"
            mkdir -p {OUTPUT_DIR}
            echo "{wildcards.sample} - CONCOCT skipped (no read depth file)" >> {OUTPUT_DIR}/skip.log
            mkdir -p {output.bins}
            touch {output.bins}/.snakemake_success
            exit 0
        fi

        # ✅ Snakemake가 `fasta_bins` 디렉토리를 인식하도록 5초 대기
        sleep 5

        if [ -z "$(ls -A {output.bins})" ]; then
            echo "❌ ERROR: No FASTA bins were created for {wildcards.sample}."
            exit 1
        fi

        echo "✅ CONCOCT completed successfully for {wildcards.sample}"
        """


###############################################################################
# 7) QUAST -> assembly_summary.csv
###############################################################################
rule assembly_summary:
    input:
        contigs = expand(f"{OUTPUT_DIR}/2_Contigs/{{sample}}_filtered_contigs.fa", sample=samples)
    output:
        summary = f"{OUTPUT_DIR}/assembly_summary.csv"
    shell:
        """
        echo "🚀 Running QUAST for assembly evaluation..."
        quast -o {OUTPUT_DIR}/1_MEGAHIT/quast_report {input.contigs}

        # ✅ QUAST 실행 후 결과 파일 확인
        QUAST_REPORT={OUTPUT_DIR}/1_MEGAHIT/quast_report/transposed_report.tsv
        if [ ! -f "$QUAST_REPORT" ]; then
            echo "❌ ERROR: QUAST output not found! ($QUAST_REPORT)"
            exit 1
        fi

        # ✅ 이전 assembly_summary.csv 삭제 후 새로 생성
        echo "Sample,N50,Total_contigs" > {output.summary}

        # ✅ N50 및 Total Contigs 추출하여 CSV 생성
        awk 'NR>1 {{print $1","$17","$14}}' "$QUAST_REPORT" >> {output.summary}

        # ✅ Snakemake에서 파일 확인하도록 추가 대기 (최대 10초)
        sleep 10

        # ✅ 파일이 정상적으로 생성되었는지 확인
        if [ ! -s {output.summary} ]; then
            echo "❌ ERROR: assembly_summary.csv is empty! Check QUAST output."
            exit 1
        fi

        echo "✅ Assembly summary successfully generated: {output.summary}"
        """

###############################################################################
# 8) gather_binning_tsv (각 binning 결과 -> contig-bin mapping TSV)
###############################################################################

rule gather_binning_tsv:
    input:
        metabat2_dir = directory(f"{OUTPUT_DIR}/3_Binning/MetaBAT2/{{sample}}"),
        maxbin2_dir  = directory(f"{OUTPUT_DIR}/3_Binning/MaxBin2/{{sample}}"),
        concoct_dir  = directory(f"{OUTPUT_DIR}/3_Binning/CONCOCT/{{sample}}/fasta_bins")
    output:
        metabat2_tsv = f"{OUTPUT_DIR}/7_DAS_tool_out/{{sample}}/metabat2_bins.tsv",
        maxbin2_tsv  = f"{OUTPUT_DIR}/7_DAS_tool_out/{{sample}}/maxbin2_bins.tsv",
        concoct_tsv  = f"{OUTPUT_DIR}/7_DAS_tool_out/{{sample}}/concoct_bins.tsv"
    run:
        import os

        # 🚀 `skip.log`에서 제외할 샘플 가져오기
        SKIP_FILE = f"{OUTPUT_DIR}/skip.log"
        skip_samples = set()

        if os.path.exists(SKIP_FILE):
            with open(SKIP_FILE, "r") as f:
                for line in f:
                    sample_name = line.strip().split(" - ")[0]
                    skip_samples.add(sample_name)

        # ✅ 현재 실행 중인 sample이 `skip.log`에 있는지 확인
        if wildcards.sample in skip_samples:
            print(f"⚠️ Skipping {wildcards.sample} (found in skip.log)")
            shell(f"touch {output.metabat2_tsv} {output.maxbin2_tsv} {output.concoct_tsv}")  # 빈 파일 생성
            return  # ✅ `exit(0)` 대신 `return`을 사용하여 Snakemake가 정상적으로 진행하도록 함

        # ✅ 디렉토리 존재 여부 확인 후, 없으면 즉시 종료
        if not os.path.exists(input.metabat2_dir) and not os.path.exists(input.maxbin2_dir) and not os.path.exists(input.concoct_dir):
            print(f"⚠️ No binning directories found for {wildcards.sample}, skipping.")
            shell(f"touch {output.metabat2_tsv} {output.maxbin2_tsv} {output.concoct_tsv}")  # 빈 파일 생성
            return  # ✅ Snakemake가 정상적으로 진행하도록 `return` 사용

        shell(r"""
        set -e

        mkdir -p {OUTPUT_DIR}/7_DAS_tool_out/{wildcards.sample}

        # 1) MetaBAT2
        rm -f {output.metabat2_tsv}
        touch {output.metabat2_tsv}
        for binfa in {input.metabat2_dir}/*.fa; do
            [ -f "$binfa" ] || continue
            binname=$(basename "$binfa" .fa)
            awk -v BN="$binname" '/^>/{{sub(/^>/,"",$1); print $1"\t"BN}}' "$binfa" >> {output.metabat2_tsv}
        done

        # 2) MaxBin2
        rm -f {output.maxbin2_tsv}
        touch {output.maxbin2_tsv}
        for binfa in {input.maxbin2_dir}/*.fasta; do
            [ -f "$binfa" ] || continue
            binname=$(basename "$binfa" .fasta)
            awk -v BN="$binname" '/^>/{{sub(/^>/,"",$1); print $1"\t"BN}}' "$binfa" >> {output.maxbin2_tsv}
        done

        # 3) CONCOCT
        rm -f {output.concoct_tsv}
        touch {output.concoct_tsv}
        if ls {input.concoct_dir}/*.fa 1> /dev/null 2>&1; then
            for binfa in {input.concoct_dir}/*.fa; do
                [ -f "$binfa" ] || continue
                binname=$(basename "$binfa" .fa)
                awk -v BN="concoct_bin_$binname" '/^>/{{sub(/^>/,"",$1); print $1"\t"BN}}' "$binfa" >> {output.concoct_tsv}
            done
        else
            echo "⚠️ No CONCOCT bins found for {wildcards.sample}, skipping."
            rm -f {output.concoct_tsv}  # 빈 파일 제거
        fi

        echo "✅ gather_binning_tsv done for sample {wildcards.sample}."
        """)





###############################################################################
# 9) DAS_Tool -> bin 통합
###############################################################################
rule das_tool:
    input:
        metabat2_tsv = f"{OUTPUT_DIR}/7_DAS_tool_out/{{sample}}/metabat2_bins.tsv",
        maxbin2_tsv  = f"{OUTPUT_DIR}/7_DAS_tool_out/{{sample}}/maxbin2_bins.tsv",
        concoct_tsv  = f"{OUTPUT_DIR}/7_DAS_tool_out/{{sample}}/concoct_bins.tsv",
        contigs = f"{OUTPUT_DIR}/2_Contigs/{{sample}}_filtered_contigs.fa"
    output:
        final_bins = directory(f"{OUTPUT_DIR}/7_DAS_tool_out/{{sample}}/DAS_Tool_run_DASTool_bins")
    params:
        prefix = f"{OUTPUT_DIR}/7_DAS_tool_out/{{sample}}/DAS_Tool_run"
    log:
        f"{OUTPUT_DIR}/7_DAS_tool_out/{{sample}}/DAS_Tool_run.log"
    shell:
        r"""
        set -e
        echo "🚀 Running DAS_Tool for {wildcards.sample}" > {log}

        # 입력 파일 확인 (Binning 결과가 없으면 DAS_Tool을 실행하지 않음)
        if [ ! -s {input.metabat2_tsv} ] && [ ! -s {input.maxbin2_tsv} ] && [ ! -s {input.concoct_tsv} ]; then
            echo "⚠️ No valid binning results for {wildcards.sample}. Skipping DAS_Tool." >> {log}
            mkdir -p {output.final_bins}
            touch {output.final_bins}/empty
            echo "{wildcards.sample} - DAS_Tool skipped (No valid binning results)" >> {OUTPUT_DIR}/skip.log
            exit 0
        fi

        # DAS_Tool 실행 (실패해도 Snakemake 중단 방지)
        if ! DAS_Tool \
            -i {input.metabat2_tsv},{input.maxbin2_tsv},{input.concoct_tsv} \
            -l metabat2,maxbin2,concoct \
            -c {input.contigs} \
            -o {params.prefix} \
            --write_bins \
            --search_engine diamond \
            --threads {MEGAHIT_THREADS} >> {log} 2>&1; then
            echo "❌ DAS_Tool failed for {wildcards.sample}" >> {log}
            echo "{wildcards.sample} - DAS_Tool skipped (DAS_Tool execution failed)" >> {OUTPUT_DIR}/skip.log
            mkdir -p {output.final_bins}
            touch {output.final_bins}/empty
            exit 0
        fi

        # DAS_Tool 실행 후 결과 폴더 확인
        if [ -d "{params.prefix}_DASTool_bins" ]; then
            echo "✅ DAS_Tool output folder found: {params.prefix}_DASTool_bins" >> {log}
        else
            echo "❌ DAS_Tool output folder not found! Skipping sample." >> {log}
            mkdir -p {OUTPUT_DIR}
            echo "{wildcards.sample} - DAS_Tool skipped (No bins generated)" >> {OUTPUT_DIR}/skip.log
            mkdir -p {output.final_bins}
            touch {output.final_bins}/empty
            exit 0
        fi
        """

###############################################################################
# 10) CheckM -> 최종 bin 품질평가 (CheckM v1)
###############################################################################
rule checkm_final:
    input:
        final_bins = f"{OUTPUT_DIR}/7_DAS_tool_out/{{sample}}/DAS_Tool_run_DASTool_bins"
    output:
        summary_csv = f"{OUTPUT_DIR}/8_checkM_summary/{{sample}}_checkm_summary.csv",
        log = f"{OUTPUT_DIR}/8_checkM_summary/{{sample}}_checkm_summary.log",
        lineage_ms = f"{OUTPUT_DIR}/8_checkM_summary/{{sample}}_temp/lineage.ms"
    run:
        import os

        # 🚀 `skip.log`에서 제외할 샘플 가져오기
        SKIP_FILE = f"{OUTPUT_DIR}/skip.log"
        skip_samples = set()

        if os.path.exists(SKIP_FILE):
            with open(SKIP_FILE, "r") as f:
                for line in f:
                    sample_name = line.strip().split(" - ")[0]
                    skip_samples.add(sample_name)

        # ✅ 현재 실행 중인 샘플이 `skip.log`에 있는지 확인
        if wildcards.sample in skip_samples:
            print(f"⚠️ Skipping {wildcards.sample} (found in skip.log)")

            # 빈 결과 파일 생성 (Snakemake가 정상적으로 진행되도록)
            for out_file in output:
                os.makedirs(os.path.dirname(out_file), exist_ok=True)
                open(out_file, "w").close()  # 빈 파일 생성
            return  # ✅ `exit(0)` 대신 `return` 사용하여 정상 종료

        shell("""
        set -e
        TEMP_DIR="{OUTPUT_DIR}/8_checkM_summary/{wildcards.sample}_temp"
        mkdir -p "$TEMP_DIR"

        echo "🚀 Running CheckM v1 on final DAS_Tool bins for sample {wildcards.sample}" | tee -a {output.log}

        # ✅ 입력 디렉토리 확인
        if [ ! -d "{input.final_bins}" ] || [ -z "$(ls -A {input.final_bins})" ]; then
            echo "⚠️ No DAS_Tool results found for {wildcards.sample}, skipping CheckM." | tee -a {output.log}
            
            # ✅ CheckM 스킵된 샘플을 skip.log에 추가
            echo "{wildcards.sample} - CheckM skipped (No DAS_Tool results)" >> {OUTPUT_DIR}/skip.log
            
            # 빈 결과 파일 생성 (Snakemake가 정상적으로 진행되도록)
            touch {output.summary_csv} {output.lineage_ms}
            exit 0
        fi

        # ✅ CheckM 실행
        if [ -n "{CHECKM_DATA_DIR}" ]; then
            checkm data setRoot "{CHECKM_DATA_DIR}" >> {output.log} 2>&1
        fi
        checkm lineage_wf -x fa "{input.final_bins}" "$TEMP_DIR" -t {MEGAHIT_THREADS} --quiet 2>> {output.log}

        # ✅ lineage.ms 파일 확인
        if [ ! -s "{output.lineage_ms}" ]; then
            echo "❌ ERROR: lineage.ms file not found or empty! Possible RAM issue." | tee -a {output.log}

            # ✅ lineage.ms 없으면 즉시 종료 (RAM 부족 가능성 높음)
            exit 1
        fi

        # ✅ CheckM 품질 평가
        checkm qa "{output.lineage_ms}" "$TEMP_DIR" -o 2 > {output.summary_csv} 2>> {output.log}

        # ✅ 결과 검증
        if [ ! -s "{output.summary_csv}" ]; then
            echo "❌ CheckM v1 summary is empty! Skipping {wildcards.sample}." | tee -a {output.log}

            # ✅ CheckM 결과 없으면 skip.log에 추가
            echo "{wildcards.sample} - CheckM skipped (CheckM summary empty)" >> {OUTPUT_DIR}/skip.log
            
            exit 1
        fi

        echo "✅ CheckM v1 done: {output.summary_csv}" | tee -a {output.log}
        """)



###############################################################################
# 11) CheckM1 통합 summary 생성
###############################################################################
rule checkm_summary:
    """
    CheckM1 결과를 단일 summary 파일로 저장 (로그 제거 후 정리)
    """
    input:
        checkm1 = expand(f"{OUTPUT_DIR}/8_checkM_summary/{{sample}}_checkm_summary.csv", sample=samples)
    output:
        summary_csv = f"{OUTPUT_DIR}/checkm_summary_combined.csv"
    run:
        import pandas as pd

        checkm1_files = input.checkm1
        combined_df = pd.DataFrame()

        for f in checkm1_files:
            try:
                sample_name = f.split("/")[-1].replace("_checkm_summary.csv", "")

                with open(f, "r") as file:
                    lines = file.readlines()

                    # ✅ "Bin Id" 또는 "Bin_ID"가 포함된 첫 번째 줄을 찾기
                    start_line = next((i for i, line in enumerate(lines) if "Bin Id" in line or "Bin_ID" in line), None)

                    if start_line is None:
                        print(f"⚠️ Warning: No 'Bin Id' found in {f}. Skipping this file.")
                        continue

                # ✅ "Bin Id" 이후 데이터만 읽기
                df = pd.read_csv(f, sep="\s+", skiprows=start_line, engine="python")

                # ✅ 컬럼이 예상되는 개수와 다를 경우 로그 출력
                if df.shape[1] != 25:
                    print(f"⚠️ Warning: Unexpected number of columns in {f}. Columns: {df.shape[1]}")

                # ✅ 예상되는 컬럼 목록
                expected_columns = [
                    "Bin_ID", "Marker_lineage", "#_genomes", "#_markers", "#_marker_sets",
                    "CheckM1_Completeness", "CheckM1_Contamination", "CheckM1_Strain_Heterogeneity",
                    "Genome_size_bp", "#_ambiguous_bases", "#_scaffolds", "#_contigs",
                    "N50_scaffolds", "N50_contigs", "Mean_scaffold_length_bp", "Mean_contig_length_bp",
                    "Longest_scaffold_bp", "Longest_contig_bp", "GC", "GC_std_scaffolds_1kbp",
                    "Coding_density", "Translation_table", "#_predicted_genes", "0", "1", "2", "3", "4", "5+", "extra"
                ]

                # ✅ 컬럼 개수가 다를 경우 자동으로 맞춤
                if len(df.columns) > len(expected_columns):
                    df = df.iloc[:, :len(expected_columns)]
                elif len(df.columns) < len(expected_columns):
                    print(f"⚠️ Warning: {f} has missing columns. Skipping this file.")
                    continue

                df.columns = expected_columns

                # ✅ 샘플명 추가
                df.insert(0, "Sample", sample_name)

                # ✅ 데이터 타입 변환 (숫자로 변환 가능한 항목 변환)
                numeric_cols = df.columns[3:]  # Sample, Bin_ID, Marker_lineage 제외
                df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors="coerce")

                # ✅ "Bin_ID"가 "-" 혹은 이상한 값인 경우 제거
                df = df[df["Bin_ID"].notnull()]
                df = df[~df["Bin_ID"].astype(str).str.contains("-{5,}")]  # "-----" 같은 이상한 값 제거

                # ✅ "CheckM1_Completeness" 값이 없는 경우 제거
                df = df[df["CheckM1_Completeness"].notnull()]

                # ✅ "Marker_lineage"가 없는 경우 제거
                df = df[df["Marker_lineage"].notnull()]

                combined_df = pd.concat([combined_df, df], ignore_index=True)

            except Exception as e:
                print(f"❌ Error reading CheckM1 file {f}: {e}")

        # ✅ 최종 CSV 저장
        combined_df.to_csv(output.summary_csv, index=False)
        print(f"✅ Cleaned CheckM1 summary saved: {output.summary_csv}")


###############################################################################
# 12) Collect final MAG FASTA files for annotation input
###############################################################################
rule collect_final_mags:
    input:
        final_bins = expand(f"{OUTPUT_DIR}/7_DAS_tool_out/{{sample}}/DAS_Tool_run_DASTool_bins", sample=samples),
        checkm_summary = f"{OUTPUT_DIR}/checkm_summary_combined.csv"
    output:
        done = f"{OUTPUT_DIR}/9_Final_MAGs/DONE.txt"
    shell:
        r"""
        set -euo pipefail
        outdir="{OUTPUT_DIR}/9_Final_MAGs"
        rm -rf "$outdir"
        mkdir -p "$outdir"

        for d in {input.final_bins}; do
            [ -d "$d" ] || continue
            sample=$(basename "$(dirname "$d")")
            for fa in "$d"/*.fa "$d"/*.fasta; do
                [ -f "$fa" ] || continue
                base=$(basename "$fa")
                cp "$fa" "$outdir/${sample}__${base}"
            done
        done

        n=$(find "$outdir" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) | wc -l | tr -d ' ')
        echo "final_mag_count=$n" > {output.done}
        if [ "$n" -eq 0 ]; then
            echo "⚠️ No final MAG FASTA files were collected into $outdir" >&2
        else
            echo "✅ Collected $n final MAG FASTA files into $outdir"
        fi
        """
