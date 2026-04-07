#############################################
# Go_MAGs_Annotation_V1.smk
# 20260406
# Heekuk Park
#
# Pipeline:
#   MAG fasta -> GTDB-Tk taxonomy/tree -> Prokka -> EggNOG-mapper -> KofamScan
#   -> EggNOG summary tables -> KofamScan summary tables
#
# Input:
#   MAGs_DIR/
#     <mag>.fa or <mag>.fasta
#
# OUTDIR/
#   0_Input_MAGs/      (normalized *.fa copies)
#   1_GTDB_Tk/        (GTDB-Tk classify and de novo tree outputs)
#   2_Prokka/         ({sample}.gff, {sample}.faa)
#   3_EggNOG/         ({sample}_eggnog.emapper.annotations)
#   4_KEGG/           ({sample}_kofamscan.txt)
#   eggnog_summary.presence.csv
#   eggnog_summary.count.csv
#   eggnog_summary.score.csv
#   kofamscan_summary.presence.csv
#   kofamscan_summary.count.csv
#   kofamscan_summary.score.csv
#
# Example run:
# mags_dir="mags_assembly_out/7_DAS_tool_out/all_bins"
# output_dir="mags_annotation_out"
#
# snakemake --snakefile /home/uhlemann/heekuk_path/MAGs/workflow/Go_MAGs_Annotation_V1.smk \
#   --config mags_dir="$mags_dir" output_dir="$output_dir" \
#   gtdbtk_data_dir="/db/gtdbtk" gtdbtk_mash_db="/db/gtdbtk/mash_db/genomic_mash_db.msh" \
#   eggnog_data_dir="/db/eggnog" kofam_scan_dir="/db/kofam_scan" \
#   --cores 24 --jobs 4 --latency-wait 60 --rerun-incomplete --keep-going
#############################################

import os
import glob

# Config 값 가져오기
MAGs_DIR = config["mags_dir"]
OUTPUT_DIR = config["output_dir"]
GTDBTK_DATA_DIR = config["gtdbtk_data_dir"]
GTDBTK_MASH_DB = config.get("gtdbtk_mash_db", os.path.join(GTDBTK_DATA_DIR, "mash_db", "genomic_mash_db.msh"))
EGGNOG_DATA_DIR = config["eggnog_data_dir"]
KOFAM_SCAN_DIR = config["kofam_scan_dir"]
KOFAM_PROFILE_DIR = config.get("kofam_profile_dir", os.path.join(KOFAM_SCAN_DIR, "profiles"))
KOFAM_KO_LIST = config.get("kofam_ko_list", os.path.join(KOFAM_SCAN_DIR, "ko_list"))
KOFAM_EXEC = config.get("kofam_exec", "exec_annotation")

GTDBTK_THREADS = int(config.get("gtdbtk_threads", 4))
PROKKA_THREADS = int(config.get("prokka_threads", 4))
EGGNOG_THREADS = int(config.get("eggnog_threads", 4))
KEGG_THREADS = int(config.get("kofamscan_threads", 4))
NORM_MAGS_DIR = f"{OUTPUT_DIR}/0_Input_MAGs"

# sample 리스트 만들기
mags = glob.glob(os.path.join(MAGs_DIR, "*.fa")) + glob.glob(os.path.join(MAGs_DIR, "*.fasta"))
MAGS_BY_SAMPLE = {
    os.path.splitext(os.path.basename(mag))[0]: mag
    for mag in mags
}
samples = sorted([
    os.path.splitext(os.path.basename(mag))[0] for mag in mags
])

print("✅ Loaded samples:", samples)

# Output 디렉토리 생성 (리넘버링 완료)
os.makedirs(f"{OUTPUT_DIR}/1_GTDB_Tk", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/2_Prokka", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/3_EggNOG", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/4_KEGG", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/logs", exist_ok=True)
os.makedirs(NORM_MAGS_DIR, exist_ok=True)

# 최종 output 정의
rule all:
    input:
        expand(f"{NORM_MAGS_DIR}/{{sample}}.fa", sample=samples),
        f"{OUTPUT_DIR}/1_GTDB_Tk/infer/infer/gtdbtk.bac120.decorated.tree-table",
        f"{OUTPUT_DIR}/1_GTDB_Tk/infer/infer/gtdbtk.bac120.decorated.tree",
        f"{OUTPUT_DIR}/1_GTDB_Tk/classify/gtdbtk.bac120.summary.tsv",
        #f"{OUTPUT_DIR}/1_GTDB_Tk/classify/gtdbtk.bac120.tree.mapping.tsv",
        expand("{output_dir}/2_Prokka/{sample}/{sample}.gff", sample=samples, output_dir=OUTPUT_DIR),
        expand("{output_dir}/3_EggNOG/{sample}_eggnog.emapper.annotations", sample=samples, output_dir=OUTPUT_DIR),
        expand("{output_dir}/4_KEGG/{sample}_kofamscan.txt", sample=samples, output_dir=OUTPUT_DIR),

        f"{OUTPUT_DIR}/eggnog_summary.presence.csv",
        f"{OUTPUT_DIR}/eggnog_summary.count.csv",
        f"{OUTPUT_DIR}/eggnog_summary.score.csv",
        f"{OUTPUT_DIR}/kofamscan_summary.presence.csv",
        f"{OUTPUT_DIR}/kofamscan_summary.count.csv",
        f"{OUTPUT_DIR}/kofamscan_summary.score.csv"
        
        

rule normalize_mags:
    input:
        mag = lambda wildcards: MAGS_BY_SAMPLE[wildcards.sample]
    output:
        mag = f"{NORM_MAGS_DIR}/{{sample}}.fa"
    shell:
        """
        cp {input.mag} {output.mag}
        """


###############################################################################
# 0️⃣ 모니터링 규칙들 (swap, cache etc.)
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

        # ⚠️ 디스크 I/O 부하가 80% 이상이면 Swap 초기화를 하지 않음 (I/O 부하가 높은 상태에서 실행하면 성능 저하 가능)
        if (( $(echo "$IO_USAGE >= 80" | bc -l) )); then
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
# 1️⃣ GTDB-Tk 실행 (계통 분석)
###############################################################################
rule gtdbtk_classify_all:
    input:
        bin_fastas = expand(f"{NORM_MAGS_DIR}/{{sample}}.fa", sample=samples)
    output:
        summary = f"{OUTPUT_DIR}/1_GTDB_Tk/classify/gtdbtk.bac120.summary.tsv",
        #tree_mapping = temp(f"{OUTPUT_DIR}/1_GTDB_Tk/classify/gtdbtk.bac120.tree.mapping.tsv")

    params:
        threads = GTDBTK_THREADS,
        gtdbtk_db = GTDBTK_DATA_DIR,
        mash_db = GTDBTK_MASH_DB
    shell:
        """
        set -e
        mkdir -p {OUTPUT_DIR}/1_GTDB_Tk

        echo "🚀 Running GTDB-Tk batch classification with {params.threads} threads"

        export GTDBTK_DATA_PATH={params.gtdbtk_db}

        gtdbtk classify_wf \
            --genome_dir {NORM_MAGS_DIR} \
            --out_dir {OUTPUT_DIR}/1_GTDB_Tk \
            --extension fa \
            --cpus {params.threads} \
            --mash_db {params.mash_db}

        echo "✅ GTDB-Tk batch classification completed"
        """



rule gtdbtk_de_novo_tree:
    input:
        bin_fastas = expand(f"{NORM_MAGS_DIR}/{{sample}}.fa", sample=samples)
    output:
        tree_table = f"{OUTPUT_DIR}/1_GTDB_Tk/infer/infer/gtdbtk.bac120.decorated.tree-table",
        tree = f"{OUTPUT_DIR}/1_GTDB_Tk/infer/infer/gtdbtk.bac120.decorated.tree"
    params:
        threads = GTDBTK_THREADS,
        gtdbtk_db = GTDBTK_DATA_DIR,
        outgroup = "p__Chloroflexota"
    shell:
        """
        set -e
        mkdir -p {OUTPUT_DIR}/1_GTDB_Tk/infer

        echo "🚀 Running GTDB-Tk de novo tree with {params.threads} threads"

        export GTDBTK_DATA_PATH={params.gtdbtk_db}

        gtdbtk de_novo_wf \
            --genome_dir {NORM_MAGS_DIR} \
            --out_dir {OUTPUT_DIR}/1_GTDB_Tk/infer \
            --extension fa \
            --cpus {params.threads} \
            --bacteria \
            --outgroup_taxon {params.outgroup}

        echo "🌳 Rooting tree..."
        gtdbtk root \
            --input_tree {OUTPUT_DIR}/1_GTDB_Tk/infer/infer/intermediate_results/gtdbtk.bac120.unrooted.tree \
            --output_tree {output.tree} \
            --outgroup_taxon {params.outgroup}

        cp {OUTPUT_DIR}/1_GTDB_Tk/infer/infer/intermediate_results/gtdbtk.bac120.tree.log \
           {output.tree_table}

        echo "✅ GTDB-Tk de novo tree and rooting completed"
        """




# classify_wf no tree
# de_novo_wf yes tree





###############################################################################
# 2️⃣ Prokka 실행
###############################################################################
rule prokka_annotation:
    input:
        bin_fasta = f"{NORM_MAGS_DIR}/{{sample}}.fa"
    output:
        gff = f"{OUTPUT_DIR}/2_Prokka/{{sample}}/{{sample}}.gff",
        faa = f"{OUTPUT_DIR}/2_Prokka/{{sample}}/{{sample}}.faa"
    params:
        kingdom = "Bacteria",
        threads = PROKKA_THREADS
    shell:
        """
        set -e
        mkdir -p $(dirname {output.gff})

        echo "🚀 Running optimized Prokka on {wildcards.sample} with {params.threads} threads"

        prokka --outdir $(dirname {output.gff}) \
               --prefix {wildcards.sample} \
               --kingdom {params.kingdom} \
               --mincontiglen 200 \
               --cpus {params.threads} \
               --force \
               {input.bin_fasta}

        echo "✅ Prokka annotation completed: $(dirname {output.gff})"
        """


###############################################################################
# 3️⃣ EggNOG-mapper 실행
###############################################################################
rule eggnog_mapper:
    input:
        faa = f"{OUTPUT_DIR}/2_Prokka/{{sample}}/{{sample}}.faa"
    output:
        annotations = f"{OUTPUT_DIR}/3_EggNOG/{{sample}}_eggnog.emapper.annotations"
    params:
        threads = EGGNOG_THREADS,
        eggnog_db = EGGNOG_DATA_DIR
    shell:
        """
        set -e
        mkdir -p {OUTPUT_DIR}/3_EggNOG
        mkdir -p {OUTPUT_DIR}/3_EggNOG/tmp_{wildcards.sample}  # 👈 temp_dir 미리 생성

        echo "🚀 Running EggNOG-mapper on {wildcards.sample} with {params.threads} threads"

        emapper.py \
            -i {input.faa} \
            --output {OUTPUT_DIR}/3_EggNOG/{wildcards.sample}_eggnog \
            --cpu {params.threads} \
            --data_dir {params.eggnog_db} \
            --dmnd_db {params.eggnog_db}/eggnog_proteins.dmnd \
            --temp_dir {OUTPUT_DIR}/3_EggNOG/tmp_{wildcards.sample} \
            --usemem \
            --override

        echo "✅ EggNOG-mapper completed: {output.annotations}"
        """

rule eggnog_summary:
    input:
        eggnog_files = lambda wildcards: sorted(glob.glob(f"{OUTPUT_DIR}/3_EggNOG/*emapper.annotations"))
    output:
        presence_out = f"{OUTPUT_DIR}/eggnog_summary.presence.csv",
        count_out = f"{OUTPUT_DIR}/eggnog_summary.count.csv",
        score_out = f"{OUTPUT_DIR}/eggnog_summary.score.csv"
    run:
        import pandas as pd
        import os
        from collections import defaultdict

        presence_dict = defaultdict(dict)
        count_dict = defaultdict(dict)
        score_dict = defaultdict(dict)
        feature_info = {}

        for file in input.eggnog_files:
            sample = os.path.basename(file).replace(".emapper.annotations", "")
            df = pd.read_csv(file, sep="\t", skiprows=4, low_memory=False)

            for _, row in df.iterrows():
                ko = str(row.get("KEGG_ko", ""))
                ec = str(row.get("EC", ""))
                go = str(row.get("GOs", ""))
                pw = str(row.get("KEGG_Pathway", ""))
                if pd.isna(ko) and pd.isna(ec) and pd.isna(go) and pd.isna(pw):
                    continue
                feature_id = f"{ko}|{ec}|{go}|{pw}"

                # Presence
                presence_dict[feature_id][sample] = 1

                # Count
                count_dict[feature_id][sample] = count_dict[feature_id].get(sample, 0) + 1

                # Score: 여기선 Score가 따로 없기 때문에 count와 동일하게 처리 (원하면 수정 가능)
                prev_score = score_dict[feature_id].get(sample, 0)
                score_dict[feature_id][sample] = max(prev_score, 1)

                if feature_id not in feature_info:
                    feature_info[feature_id] = (ko, ec, go, pw)

        samples_all = sorted(set(os.path.basename(f).replace(".emapper.annotations", "") for f in input.eggnog_files))

        def write_table(data_dict, outfile, prefix):
            df = pd.DataFrame.from_dict(data_dict, orient='index').fillna(0)
            df = df[samples_all]  # consistent column order
            df.index.name = "Feature_ID"
            df.insert(0, "KO", [feature_info[i][0] for i in df.index])
            df.insert(1, "EC", [feature_info[i][1] for i in df.index])
            df.insert(2, "GO", [feature_info[i][2] for i in df.index])
            df.insert(3, "KEGG_Pathway", [feature_info[i][3] for i in df.index])
            df.index = [f"{prefix}{i+1:05d}" for i in range(len(df))]
            df.to_csv(outfile, index_label="Feature_ID", encoding="utf-8")

        write_table(presence_dict, output.presence_out, prefix="presence")
        write_table(count_dict, output.count_out, prefix="count")
        write_table(score_dict, output.score_out, prefix="score")

        
# cd /db/eggnog
# download_eggnog_data.py --data_dir .
###############################################################################
#4️⃣ KofamScan 실행
###############################################################################

rule clean_faa_headers:
    input:
        faa = f"{OUTPUT_DIR}/2_Prokka/{{sample}}/{{sample}}.faa"
    output:
        cleaned_faa = f"{OUTPUT_DIR}/2_Prokka/{{sample}}/{{sample}}.clean.faa"
    shell:
        r"""
        awk 'BEGIN {{i=0}} /^>/ {{print ">clean_" ++i; next}} !/^>/ && NF > 0 {{print}}' {input.faa} > {output.cleaned_faa}
        """


rule prepare_ko_list:
    input:
        original = KOFAM_KO_LIST
    output:
        filtered = f"{OUTPUT_DIR}/4_KEGG/ko_list.filtered.tsv"
    shell:
        """
        awk 'NR==1 || $1 ~ /^K[0-9]{{5}}$/' {input.original} > {output.filtered}
        """



rule kofamscan_annotation:
    input:
        faa = f"{OUTPUT_DIR}/2_Prokka/{{sample}}/{{sample}}.clean.faa"
    output:
        kofamscan_results = f"{OUTPUT_DIR}/4_KEGG/{{sample}}_kofamscan.txt"
    params:
        profile = KOFAM_PROFILE_DIR,
        ko_list = KOFAM_KO_LIST,
        kofam_exec = KOFAM_EXEC,
        threads = KEGG_THREADS
    resources:
        kofamscan_slots=1
    shell:
        """
        set -e
        mkdir -p {OUTPUT_DIR}/4_KEGG

        echo "▶ CHECK ko_list: '{params.ko_list}'" >&2
        echo "▶ CHECK input faa: '{input.faa}'" >&2

        if [ ! -s {input.faa} ]; then
            echo "⚠️ {input.faa} is empty. Skipping sample {wildcards.sample}." >&2
            touch {output.kofamscan_results}
            exit 0
        fi

        echo "🚀 Running KofamScan on {wildcards.sample} with {params.threads} threads"

        {params.kofam_exec} \
            -o {output.kofamscan_results} \
            -f detail-tsv \
            --cpu {params.threads} \
            --profile {params.profile} \
            --ko-list {params.ko_list} \
            {input.faa} || {{
                echo "❌ KofamScan failed for {wildcards.sample}" >&2
                rm -f {output.kofamscan_results}
                exit 1
            }}

        echo "✅ KofamScan completed: {output.kofamscan_results}"
        """


rule kofamscan_summary:
    input:
        kofamscan_results=lambda wildcards: sorted(
            glob.glob(f"{OUTPUT_DIR}/4_KEGG/*_kofamscan.txt")
        )
    output:
        presence=f"{OUTPUT_DIR}/kofamscan_summary.presence.csv",
        count_out=f"{OUTPUT_DIR}/kofamscan_summary.count.csv",
        score_out=f"{OUTPUT_DIR}/kofamscan_summary.score.csv"
    run:
        import pandas as pd
        import os
        from collections import defaultdict

        presence_dict = defaultdict(dict)
        count_dict = defaultdict(dict)
        score_dict = defaultdict(dict)
        ko_info_dict = {}

        for file in input.kofamscan_results:
            sample = os.path.basename(file).replace("_kofamscan.txt", "")
            with open(file) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue

                    parts = line.strip().split("\t")
                    if len(parts) < 6:
                        continue

                    gene = parts[0]
                    ko = parts[1]

                    try:
                        thrshld = float(parts[2])
                        score = float(parts[3])
                    except ValueError:
                        continue  # KO 자리나 수치가 깨졌을 경우 skip

                    evalue = parts[4]
                    definition = parts[5].strip('"')

                    # presence
                    presence_dict[ko][sample] = 1

                    # count
                    count_dict[ko][sample] = count_dict[ko].get(sample, 0) + 1

                    # score (max per KO per sample)
                    prev_score = score_dict[ko].get(sample, 0)
                    score_dict[ko][sample] = max(prev_score, score)

                    # KO 정의 저장
                    if ko not in ko_info_dict:
                        ko_info_dict[ko] = definition

        samples_all = sorted(set(os.path.basename(f).replace("_kofamscan.txt", "") for f in input.kofamscan_results))

            
        def write_table(data_dict, outfile, prefix):
            df = pd.DataFrame.from_dict(data_dict, orient='index').fillna(0)

            # 모든 샘플을 column으로 강제 포함 (없으면 0으로 추가)
            for sample in samples_all:
                if sample not in df.columns:
                    df[sample] = 0

            df = df[samples_all]  # column 순서 맞추기
            df.insert(0, "KO", df.index)
            df.insert(1, "Definition", df["KO"].map(ko_info_dict))
            df.index = [f"{prefix}{i+1:05d}" for i in range(len(df))]
            df.to_csv(outfile, sep=',', index_label="Feature", encoding="utf-8")

        write_table(presence_dict, output.presence, prefix="presence")
        write_table(count_dict, output.count_out, prefix="count")
        write_table(score_dict, output.score_out, prefix="score")




        
# 다운로드 (프로필 + KO ID 리스트)
# cd /db/kofam_scan
# git clone https://github.com/takaram/kofam_scan.git
# wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
# wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz

# 압축 해제
# tar -xzf profiles.tar.gz
# gunzip ko_list.gz

###############################################################################
# 5️⃣
###############################################################################
