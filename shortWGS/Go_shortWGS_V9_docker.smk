###############################################
# GoBacterialWGS.smk  (WGS_DB2 구조 전제)
# Flow:
# 1) fastp QC
# 2) genotyping (srst2)
#    a) args_srst2 (배치)
#    b) plas_srst2 (배치)
#    c) kraken2 → (top species) → per-sample mlst_srst2
# 3) TETyper (per-sample → summary)
# 4) MultiQC
# 5) Summary (Rmd: 20250930_Summary_WGS_tem_v2.Rmd)
# 6) summary_master.csv (통합 요약 CSV)
###############################################

import os, re, glob, json, subprocess, shlex, datetime, csv
from pathlib import Path

# ---------- Config ----------
FASTQ_DIR  = config["fastq_dir"]
OUTPUT_DIR = config["output_dir"]
DB_DIR     = config["wgs_db"]             # e.g., /media/uhlemann/core4/DB/WGS_DB2
KRAKEN_DB  = config["kraken_db"]          # e.g., /media/.../kraken2DB/k2_pluspfp_...

# DB: WGS_DB2 내부 고정 경로
CARD_DIR        = f"{DB_DIR}/CARD"
PLASMID_DIR     = f"{DB_DIR}/PlasmidFinder"
MLST_DB_ROOT    = f"{DB_DIR}/MLST"
TETYPE_REF_DIR  = f"{DB_DIR}/tetyper"

# TETyper 내부
STRUCT_PROF = f"{TETYPE_REF_DIR}/struct_profiles.txt"
SNP_PROF    = f"{TETYPE_REF_DIR}/snp_profiles.txt"
SHOW_REGION = config.get("show_region", "7202-8083")

# Global skip log
SKIP_LOG = f"{OUTPUT_DIR}/skip.log"

# ---------- Helpers ----------
def _first_fasta_in(dirpath):
    fas = sorted(glob.glob(os.path.join(dirpath, "*.fasta")))
    if not fas:
        raise FileNotFoundError(f"No .fasta in: {dirpath}")
    return fas[0]

ARGS_DB     = _first_fasta_in(CARD_DIR)         # e.g., CARD_v3.*.fasta
PLASMID_DB  = _first_fasta_in(PLASMID_DIR)      # e.g., PlasmidFinder*.fasta

def _find_mlst_definitions(species_dir):
    """
    'profiles_csv' (dir), 'profiles_csv' (file), 'profiles.csv' (file)
    순으로 우선 탐색. 없으면 *.csv(이름에 profile 포함 선호) fallback.
    """
    c1 = os.path.join(species_dir, "profiles_csv")
    if os.path.isdir(c1):  return c1
    if os.path.isfile(c1): return c1
    c2 = os.path.join(species_dir, "profiles.csv")
    if os.path.isfile(c2): return c2
    csvs = sorted(glob.glob(os.path.join(species_dir, "*.csv")))
    if csvs:
        preferred = [p for p in csvs if "profile" in os.path.basename(p).lower()]
        return preferred[-1] if preferred else csvs[-1]
    return None

# Kraken2 종명 표준화 → MLST 코드 매핑
def parse_kraken_species(raw: str) -> str:
    if not raw:
        return None
    s = raw.strip()
    s = re.sub(r'^[ksfpgco]__', '', s, flags=re.IGNORECASE)   # k__/p__/c__/o__/f__/g__/s__
    s = re.split(r'[\(\[\{,;/]| strain | subgroup | subsp\. | subspecies | biovar | serovar ', s, 1, flags=re.IGNORECASE)[0]
    s = s.replace('_', ' ')
    s = re.sub(r'\s+', ' ', s).strip()
    abbrev = re.match(r'^([A-Z])\.\s*([a-z]+)$', s)
    if abbrev:
        initial = abbrev.group(1).upper()
        species = abbrev.group(2).lower()
        initial_to_genus = {
            'A': 'Acinetobacter',
            'C': 'Citrobacter',
            'E': 'Enterococcus',
            'K': 'Klebsiella',
            'P': 'Pseudomonas',
            'S': 'Staphylococcus',
        }
        genus = initial_to_genus.get(initial)
        if genus:
            s = f'{genus} {species}'
    s = ' '.join(w.capitalize() for w in s.split())
    if len(s.split()) < 2:
        return None
    genus, species = s.split()[:2]
    return f"{genus} {species}"

def norm_species_to_code(std_name: str) -> str:
    if not std_name:
        return None
    key = std_name.strip().lower()
    base = {
        "acinetobacter baumannii": "a_baumannii",
        "citrobacter amalonaticus": "c_amalonaticus",
        "citrobacter freundii": "c_fruendii",
        "enterobacter cloacae": "e_cloacae",
        "escherichia coli": "e_coli",
        "enterococcus faecalis": "e_faecalis",
        "enterococcus faecium": "e_faecium",
        "klebsiella aerogenes": "k_aerogenes",
        "klebsiella oxytoca": "k_oxytoca",
        "klebsiella pneumoniae": "k_pneumoniae",
        "pseudomonas aeruginosa": "p_aeruginosa",
        "staphylococcus aureus": "s_aureus",
        "staphylococcus hominis": "s_hominis",
    }
    synonyms = {
        "enterococcus faecum": "enterococcus faecium",
        "enterococcus faeceum": "enterococcus faecium",
        "enterococcus faecium ": "enterococcus faecium",
        "enterococcus faecallis": "enterococcus faecalis",
        "enterobacter aerogenes": "klebsiella aerogenes",
        
        "enterobacter cloacae complex": "enterobacter cloacae",
        "enterobacter hormaechei": "enterobacter cloacae",
        "enterobacter asburiae": "enterobacter cloacae",
        "enterobacter kobei": "enterobacter cloacae",
        "enterobacter ludwigii": "enterobacter cloacae",
        "enterobacter mori": "enterobacter cloacae",
        "enterobacter xiangfangensis": "enterobacter cloacae",
    }
    if key in synonyms:
        key = synonyms[key]
    return base.get(key)

# ---------- discover samples ----------
def _find_samples(fastq_dir):
    # 다양한 naming pattern 지원
    patterns = [
        "*_R1_001.fastq.gz",
        "*_R1.fastq.gz",
        "*.R1.fastq.gz"
    ]
    r1s = []
    for pat in patterns:
        r1s += sorted(glob.glob(os.path.join(fastq_dir, pat)))

    samples = set()
    for r1 in r1s:
        b = os.path.basename(r1)
        # _R1, _R1_001, .R1 모두 인식
        m = re.match(r"(.+?)(?:[_\.]R1(?:_001)?)\.fastq\.gz$", b)
        if m:
            samples.add(m.group(1))

    return sorted(samples)

SAMPLES = _find_samples(FASTQ_DIR)
if not SAMPLES:
    raise ValueError(
        f"[FATAL] No samples found under {FASTQ_DIR}. "
        f"Expecting files like <SAMPLE>_R1.fastq.gz / <SAMPLE>.R1.fastq.gz / <SAMPLE>_R1_001.fastq.gz"
    )

# ---------- Targets ----------
rule all:
    input:
        # 1) fastp → multiqc
        fastp_done = expand(f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.done.txt", sample=SAMPLES),
        multiqc    = f"{OUTPUT_DIR}/multiqc_report.html",

        # 2) ARG/PLASMID (배치 실행)
        args_done  = f"{OUTPUT_DIR}/3_ARGs_srst2_out/args_srst2.done.txt",
        plas_done  = f"{OUTPUT_DIR}/4_Plasmid_srst2_out/plasmid_srst2.done.txt",

        # 3) kraken2(샘플별) → mlst(샘플별) → 집계 done
        kraken_top = expand(f"{OUTPUT_DIR}/2a_kraken2/{{sample}}.top_species.txt", sample=SAMPLES),
        mlst_all   = f"{OUTPUT_DIR}/2_ST_srst2_out/mlst_srst2.done.txt",
        mlst_table = f"{OUTPUT_DIR}/2_ST_srst2_out/mlst_master.csv",

        # 4) TETyper
        tetyper_all = f"{OUTPUT_DIR}/5_TETyper/tetyper.done.txt",

        # 5) Master CSV
        master_csv = f"{OUTPUT_DIR}/summary_master.csv",
        
        # 6) Summary (Rmd)
        summary = f"{OUTPUT_DIR}/summary_report.html",

        # 로그
        skiplog = SKIP_LOG

# ---------- Utils ----------
rule init_skip_log:
    output: SKIP_LOG
    shell:  "mkdir -p {OUTPUT_DIR} && : > {output}"

# ---------- 1) fastp QC ----------
rule fastp_qc:
    input:
        r1 = lambda wc: next((p for p in [
            f"{FASTQ_DIR}/{wc.sample}_R1_001.fastq.gz",
            f"{FASTQ_DIR}/{wc.sample}_R1.fastq.gz",
            f"{FASTQ_DIR}/{wc.sample}.R1.fastq.gz"
        ] if os.path.exists(p)), ""),
        r2 = lambda wc: next((p for p in [
            f"{FASTQ_DIR}/{wc.sample}_R2_001.fastq.gz",
            f"{FASTQ_DIR}/{wc.sample}_R2.fastq.gz",
            f"{FASTQ_DIR}/{wc.sample}.R2.fastq.gz"
        ] if os.path.exists(p)), ""),
    output:
        r1_out = f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.clean.R1.fastq.gz",
        r2_out = f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.clean.R2.fastq.gz",
        json   = f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.json",
        html   = f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.html",
        done   = f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.done.txt"
    params:
        skip_log = SKIP_LOG
    shell:
        r'''
        set -euo pipefail
        mkdir -p {OUTPUT_DIR}/1_fastp_out

        if [ ! -s "{input.r1}" ] || [ ! -s "{input.r2}" ]; then
            ( flock -x 200; echo "{wildcards.sample}\tFASTQ_NOT_FOUND" >> {params.skip_log}; ) 200> {params.skip_log}.lock
            : > {output.r1_out}; : > {output.r2_out}; : > {output.json}; : > {output.html}; : > {output.done}
            exit 0
        fi

        R1_SIZE=$(stat -c %s {input.r1}); R2_SIZE=$(stat -c %s {input.r2})
        if [ "$R1_SIZE" -lt 10000 ] || [ "$R2_SIZE" -lt 10000 ]; then
            ( flock -x 200; echo "{wildcards.sample}\tTOO_SMALL" >> {params.skip_log}; ) 200> {params.skip_log}.lock
            : > {output.r1_out}; : > {output.r2_out}; : > {output.json}; : > {output.html}; : > {output.done}
            exit 0
        fi

        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1_out} -O {output.r2_out} \
              -h {output.html} -j {output.json}
        touch {output.done}
        '''

# ---------- 4) MultiQC (fastp 요약) ----------
rule multiqc:
    input:
        fastp_json = expand(f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.json", sample=SAMPLES),
        skip       = SKIP_LOG
    output:
        html = f"{OUTPUT_DIR}/multiqc_report.html"
    shell:
        r'''
        mkdir -p {OUTPUT_DIR}
        rm -f {output.html}
        rm -rf {OUTPUT_DIR}/multiqc_report_data {OUTPUT_DIR}/multiqc_data
        multiqc {OUTPUT_DIR}/1_fastp_out -o {OUTPUT_DIR} -n multiqc_report.html
        '''

# ---------- 2) genotyping (srst2) ----------
# 2a) ARGs (배치)
rule args_srst2:
    input:
        r1 = expand(f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.clean.R1.fastq.gz", sample=SAMPLES),
        r2 = expand(f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.clean.R2.fastq.gz", sample=SAMPLES),
        gene_db = ARGS_DB
    output:
        done = f"{OUTPUT_DIR}/3_ARGs_srst2_out/args_srst2.done.txt"
    params:
        outprefix      = f"{OUTPUT_DIR}/3_ARGs_srst2_out/__genes__CARD",
        # 결과 파일: 접미사 유무 모두 허용
        genes_exact    = f"{OUTPUT_DIR}/3_ARGs_srst2_out/__genes__CARD__genes__CARD__results.txt",
        genes_wild     = f"{OUTPUT_DIR}/3_ARGs_srst2_out/__genes__CARD__genes__CARD_*__results.txt",
        full_exact     = f"{OUTPUT_DIR}/3_ARGs_srst2_out/__genes__CARD__fullgenes__CARD__results.txt",
        full_wild      = f"{OUTPUT_DIR}/3_ARGs_srst2_out/__genes__CARD__fullgenes__CARD_*__results.txt"
    shell:
        r'''
        set -euo pipefail
        mkdir -p "{OUTPUT_DIR}/3_ARGs_srst2_out"

        srst2 --threads 6 --input_pe {input.r1} {input.r2} \
            --forward .clean.R1 --reverse .clean.R2 \
            --output "{params.outprefix}" --log --gene_db "{input.gene_db}" || true

        if compgen -G "{params.genes_exact}" > /dev/null || \
           compgen -G "{params.genes_wild}"  > /dev/null || \
           compgen -G "{params.full_exact}"  > /dev/null || \
           compgen -G "{params.full_wild}"   > /dev/null; then
            touch "{output.done}"
        else
            echo "[WARNING] ARGs results missing/empty" >&2
            ls -l "{OUTPUT_DIR}/3_ARGs_srst2_out" >&2 || true
            LOGF="{params.outprefix}.log"
            if [ -s "$LOGF" ]; then
              ( echo "== tail -n 200 $LOGF ==" >&2; tail -n 200 "$LOGF" >&2 ) || true
            fi
            exit 1
        fi
        '''

# 2b) Plasmid (배치)
rule plas_srst2:
    input:
        r1 = expand(f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.clean.R1.fastq.gz", sample=SAMPLES),
        r2 = expand(f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.clean.R2.fastq.gz", sample=SAMPLES),
        gene_db = PLASMID_DB
    output:
        done = f"{OUTPUT_DIR}/4_Plasmid_srst2_out/plasmid_srst2.done.txt"
    params:
        outprefix      = f"{OUTPUT_DIR}/4_Plasmid_srst2_out/__genes__PlasmidFinder",
        genes_exact    = f"{OUTPUT_DIR}/4_Plasmid_srst2_out/__genes__PlasmidFinder__genes__PlasmidFinder__results.txt",
        genes_wild     = f"{OUTPUT_DIR}/4_Plasmid_srst2_out/__genes__PlasmidFinder__genes__PlasmidFinder_*__results.txt",
        full_exact     = f"{OUTPUT_DIR}/4_Plasmid_srst2_out/__genes__PlasmidFinder__fullgenes__PlasmidFinder__results.txt",
        full_wild      = f"{OUTPUT_DIR}/4_Plasmid_srst2_out/__genes__PlasmidFinder__fullgenes__PlasmidFinder_*__results.txt"
    shell:
        r'''
        set -euo pipefail
        mkdir -p "{OUTPUT_DIR}/4_Plasmid_srst2_out"

        srst2 --threads 6 --input_pe {input.r1} {input.r2} \
            --forward .clean.R1 --reverse .clean.R2 \
            --output "{params.outprefix}" --log --gene_db "{input.gene_db}" || true

        if compgen -G "{params.genes_exact}" > /dev/null || \
           compgen -G "{params.genes_wild}"  > /dev/null || \
           compgen -G "{params.full_exact}"  > /dev/null || \
           compgen -G "{params.full_wild}"   > /dev/null; then
            touch "{output.done}"
        else
            echo "[WARNING] Plasmid results missing/empty" >&2
            ls -l "{OUTPUT_DIR}/4_Plasmid_srst2_out" >&2 || true
            LOGF="{params.outprefix}.log"
            if [ -s "$LOGF" ]; then
              ( echo "== tail -n 200 $LOGF ==" >&2; tail -n 200 "$LOGF" >&2 ) || true
            fi
            exit 1
        fi
        '''

# 2c) Kraken2 → top species (샘플별)
rule kraken2_classify:
    input:
        r1 = f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.clean.R1.fastq.gz",
        r2 = f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.clean.R2.fastq.gz"
    output:
        kreport = f"{OUTPUT_DIR}/2a_kraken2/{{sample}}.kreport",
        out     = f"{OUTPUT_DIR}/2a_kraken2/{{sample}}.kraken",
        done    = f"{OUTPUT_DIR}/2a_kraken2/{{sample}}.done.txt"
    params:
        db = KRAKEN_DB
    shell:
        r'''
        mkdir -p {OUTPUT_DIR}/2a_kraken2
        if [ ! -s {input.r1} ] || [ ! -s {input.r2} ]; then
            : > {output.kreport}; : > {output.out}; touch {output.done}; exit 0
        fi
        kraken2 --db {params.db} --paired --gzip-compressed \
                --report {output.kreport} --output {output.out} \
                {input.r1} {input.r2}
        touch {output.done}
        '''

rule kraken2_pick_species:
    input:
        kreport = f"{OUTPUT_DIR}/2a_kraken2/{{sample}}.kreport",
        done    = f"{OUTPUT_DIR}/2a_kraken2/{{sample}}.done.txt"
    output:
        top = f"{OUTPUT_DIR}/2a_kraken2/{{sample}}.top_species.txt"
    run:
        top_name = "unknown"; top_pct = -1.0
        with open(input.kreport) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 6: continue
                try: pct = float(p[0].strip())
                except: continue
                if p[3].strip() == 'S' and pct > top_pct:
                    top_pct = pct; top_name = p[5].strip()
        with open(output.top, "w") as w: w.write(top_name + "\n")

# 2d) per-sample MLST (kraken 결과 기반)
rule mlst_srst2_per_sample:
    input:
        r1  = f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.clean.R1.fastq.gz",
        r2  = f"{OUTPUT_DIR}/1_fastp_out/{{sample}}.clean.R2.fastq.gz",
        top = f"{OUTPUT_DIR}/2a_kraken2/{{sample}}.top_species.txt"
    output:
        result = f"{OUTPUT_DIR}/2_ST_srst2_out/{{sample}}.mlst.txt",
        logf   = f"{OUTPUT_DIR}/2_ST_srst2_out/{{sample}}.mlst.log",
        done   = f"{OUTPUT_DIR}/2_ST_srst2_out/{{sample}}.mlst.done.txt"
    run:
        import os, glob, shlex, subprocess, datetime, pathlib

        os.makedirs(f"{OUTPUT_DIR}/2_ST_srst2_out", exist_ok=True)

        with open(output.logf, "w") as LOG:
            def w(msg): print(msg, file=LOG, flush=True)

            w(f"[MLST] {wildcards.sample} @ {datetime.datetime.now()}")
            w(f"FASTQs: {input.r1}, {input.r2}")

            # 빈 FASTQ 가드
            if (not os.path.exists(input.r1) or os.path.getsize(input.r1) == 0 or
                not os.path.exists(input.r2) or os.path.getsize(input.r2) == 0):
                w("[SKIP] empty FASTQ")
                pathlib.Path(output.result).touch()
                pathlib.Path(output.done).touch()
                return

            # Kraken2 top species -> 표준화 -> MLST 코드
            with open(input.top) as f:
                raw_species = f.read().strip()
            std_species = parse_kraken_species(raw_species)
            mlst_code   = norm_species_to_code(std_species)
            w(f"Top species (kraken2): '{raw_species}'")
            w(f"Normalized species:    '{std_species}'")
            w(f"MLST code:             '{mlst_code}'")

            if not mlst_code:
                w(f"[SKIP] unmapped species: '{raw_species}'")
                pathlib.Path(output.result).touch()
                pathlib.Path(output.done).touch()
                return

            # MLST DB 찾기
            species_dir = os.path.join(MLST_DB_ROOT, mlst_code)
            w(f"MLST_DB_ROOT: {MLST_DB_ROOT}")
            w(f"species dir: {species_dir}")

            fasta_candidates = (
                glob.glob(os.path.join(species_dir, "*.fasta")) +
                glob.glob(os.path.join(species_dir, "*.fa")) +
                glob.glob(os.path.join(species_dir, "*.fna"))
            )
            w(f"found fasta candidates: {fasta_candidates}")
            if not fasta_candidates:
                w(f"[SKIP] no fasta in '{species_dir}'")
                pathlib.Path(output.result).touch()
                pathlib.Path(output.done).touch()
                return
            db_fa = sorted(fasta_candidates)[0]

            # profiles_csv / profiles.csv / 임의 *.csv 탐색
            prof_path = _find_mlst_definitions(species_dir)
            w(f"profiles_csv path: {prof_path}")
            if not prof_path or not os.path.exists(prof_path):
                w(f"[SKIP] MLST definitions missing for '{mlst_code}' under '{species_dir}'")
                pathlib.Path(output.result).touch()
                pathlib.Path(output.done).touch()
                return

            outprefix = os.path.join(
                OUTPUT_DIR, "2_ST_srst2_out",
                f"{wildcards.sample}__mlst__{mlst_code}"
            )

            cmd = (
                "srst2 "
                f"--input_pe {input.r1} {input.r2} "
                f"--mlst_db {db_fa} "
                f"--mlst_definitions {prof_path} "
                "--mlst_delimiter '_' "
                f"--output {outprefix} "
                "--forward .clean.R1 --reverse .clean.R2"
            )
            w(f"CMD: {cmd}")
            rc = subprocess.call(shlex.split(cmd), stdout=LOG, stderr=LOG)

            # 결과 파일(가변 네이밍) 가장 최신 것 선택
            res = sorted(glob.glob(outprefix + "*__results.txt"), key=os.path.getmtime, reverse=True)
            if rc == 0 and res and os.path.getsize(res[0]) > 0:
                os.replace(res[0], output.result)
                w(f"[OK] wrote {output.result} <- {res[0]}")
            else:
                w(f"[WARN] srst2 failed or empty (rc={rc}); write empty result")
                pathlib.Path(output.result).touch()

            pathlib.Path(output.done).touch()

rule mlst_srst2:
    input: expand(f"{OUTPUT_DIR}/2_ST_srst2_out/{{sample}}.mlst.done.txt", sample=SAMPLES)
    output: f"{OUTPUT_DIR}/2_ST_srst2_out/mlst_srst2.done.txt"
    shell:  "touch {output}"


rule make_mlst_table:
    input:
        mlst_done = f"{OUTPUT_DIR}/2_ST_srst2_out/mlst_srst2.done.txt"
    output:
        table = f"{OUTPUT_DIR}/2_ST_srst2_out/mlst_master.csv"
    run:
        import os, glob, csv, re

        os.makedirs(os.path.dirname(output.table), exist_ok=True)

        META = {
            "sample","st","alleles","profile","scheme","species",
            "mismatches","uncertainty","depth","maxmaf"
        }

        # ── MLST 코드(s_aureus 등) 추출
        def infer_mlst_code(sample):
            base = os.path.join(OUTPUT_DIR, "2_ST_srst2_out")
            patterns = [
                os.path.join(base, f"{sample}__mlst__*__*.sorted.bam"),
                os.path.join(base, f"{sample}__mlst__*__*.bam"),
                os.path.join(base, f"{sample}__mlst__*__*__results.txt"),
                os.path.join(base, f"{sample}__mlst__*"),
            ]
            hits = []
            for pat in patterns:
                hits.extend(glob.glob(pat))
            if not hits:
                return "NA"
            hit = max(hits, key=lambda p: os.path.getmtime(p))
            m = re.search(r"__mlst__(.+?)__", os.path.basename(hit))
            if not m:
                m = re.search(r"__mlst__(.+?)__", hit)
            return m.group(1) if m else "NA"

        def read_tsv_rows(path):
            rows = []
            with open(path) as fh:
                lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
            if not lines:
                return [], []
            header = lines[0].split("\t")
            for ln in lines[1:]:
                row = ln.split("\t")
                rec = {header[i]: (row[i] if i < len(row) else "") for i in range(len(header))}
                rows.append(rec)
            return header, rows

        def val(rec, name):
            if name in rec:
                return rec[name]
            name_low = name.strip().lower()
            for k in rec.keys():
                if k.strip().lower() == name_low:
                    return rec[k]
            return ""

        # ── 최종 헤더 (MLST_Scheme 제거)
        fieldnames = [
            "Sample",
            "MLST_DB_Code",
            "MLST_ST",
            "Allele1","Allele2","Allele3","Allele4","Allele5","Allele6","Allele7",
            "mismatches","uncertainty","depth","maxMAF"
        ]

        with open(output.table, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=fieldnames)
            w.writeheader()

            mlst_paths = sorted(glob.glob(f"{OUTPUT_DIR}/2_ST_srst2_out/*.mlst.txt"))
            for p in mlst_paths:
                header, rows = read_tsv_rows(p)
                if not header or not rows:
                    continue

                header_low = [h.strip().lower() for h in header]
                loci_in_file = [header[i] for i, h in enumerate(header_low) if h not in META]

                sample_base = os.path.basename(p).replace(".mlst.txt", "")
                mlst_code = infer_mlst_code(sample_base)

                for idx, rec in enumerate(rows, start=1):
                    sm = val(rec, "Sample").strip() or sample_base
                    sample_out = f"{sm}_{idx}" if len(rows) > 1 else sm

                    st      = val(rec, "ST") or "NA"
                    mism    = val(rec, "mismatches") or "NA"
                    uncrt   = val(rec, "uncertainty") or "NA"
                    depth   = val(rec, "depth") or "NA"
                    maxmaf  = val(rec, "maxMAF") or "NA"

                    # Allele 토큰 구성
                    tokens = []
                    for locus in loci_in_file:
                        v = rec.get(locus, "") or val(rec, locus)
                        v = v.strip()
                        if v:
                            tokens.append(f"{locus}({v})")

                    # 7칸 고정
                    if len(tokens) < 7:
                        tokens += ["NA"] * (7 - len(tokens))
                    elif len(tokens) > 7:
                        tokens = tokens[:7]

                    row_out = {
                        "Sample": sample_out,
                        "MLST_DB_Code": mlst_code,
                        "MLST_ST": st,
                        "Allele1": tokens[0],
                        "Allele2": tokens[1],
                        "Allele3": tokens[2],
                        "Allele4": tokens[3],
                        "Allele5": tokens[4],
                        "Allele6": tokens[5],
                        "Allele7": tokens[6],
                        "mismatches": mism,
                        "uncertainty": uncrt,
                        "depth": depth,
                        "maxMAF": maxmaf,
                    }
                    w.writerow(row_out)

# ---------- 3) TETyper ----------
TETY_ENV = config.get("tetyper_env", "tetyper")
TETY_CMD = config.get("tetyper_cmd", f"conda run -n {TETY_ENV} python /home/uhlemann/.local/bin/TETyper.py")

# ---------- 3) TETyper ----------
TETY_ENV = config.get("tetyper_env", "tetyper")
TETY_CMD = config.get("tetyper_cmd", f"conda run -n {TETY_ENV} python /home/uhlemann/.local/bin/TETyper.py")

rule tetyper:
    input:
        r1     = lambda wc: f"{OUTPUT_DIR}/1_fastp_out/{wc.sample}.clean.R1.fastq.gz",
        r2     = lambda wc: f"{OUTPUT_DIR}/1_fastp_out/{wc.sample}.clean.R2.fastq.gz",
        ref    = f"{DB_DIR}/tetyper/Tn4401b-1.fasta",
        struct = f"{DB_DIR}/tetyper/struct_profiles.txt",
        snp    = f"{DB_DIR}/tetyper/snp_profiles.txt",
        args_done = f"{OUTPUT_DIR}/3_ARGs_srst2_out/args_srst2.done.txt"
    output:
        summary = f"{OUTPUT_DIR}/5_TETyper/{{sample}}_summary.txt"
    params:
        prefix      = lambda wc: f"{OUTPUT_DIR}/5_TETyper/{wc.sample}",
        show_region = SHOW_REGION,
        skip_log    = SKIP_LOG,
        tety_cmd    = TETY_CMD,
        tety_env    = TETY_ENV,
        tmpdir      = lambda wc: f"{OUTPUT_DIR}/5_TETyper/tmp/{wc.sample}",
        threads     = 4,
        mem_gb      = 16,
        args_glob1  = lambda wc: f"{OUTPUT_DIR}/3_ARGs_srst2_out/__genes__CARD__genes__CARD_*__results.txt",
        args_glob2  = lambda wc: f"{OUTPUT_DIR}/3_ARGs_srst2_out/__genes__CARD__fullgenes__CARD_*__results.txt"
    log:
        tlog = f"{OUTPUT_DIR}/5_TETyper/{{sample}}.tetyper.log"
    run:
        import os, re, csv, glob, pathlib

        os.makedirs(f"{OUTPUT_DIR}/5_TETyper", exist_ok=True)
        os.makedirs(params.tmpdir, exist_ok=True)

        def log_msg(msg: str):
            with open(log.tlog, "a") as L:
                L.write(msg.rstrip() + "\n")

        def norm_sample(s: str) -> str:
            if not s: return ""
            s = s.strip()
            s = re.sub(r"\.fastq(\.gz)?$", "", s, flags=re.I)
            s = re.sub(r"\.clean$", "", s, flags=re.I)
            s = re.sub(r"[_\.]R[12](?:_001)?$", "", s, flags=re.I)
            return s

        def same_sample(a: str, b: str) -> bool:
            if not a or not b: return False
            return (a == b) or a.startswith(b) or b.startswith(a)

        def looks_like_kpc(text: str) -> bool:
            if not text: return False
            t = re.sub(r"[^A-Za-z0-9]+", "", text).lower()
            return "kpc" in t

        def sample_has_kpc(sample):
            arg_files = sorted(glob.glob(params.args_glob1)) + sorted(glob.glob(params.args_glob2))
            if not arg_files:
                log_msg("[KPC] No ARG result files found")
                return False
            sm_ref = norm_sample(sample)
            col_candidates = ["gene","allele","amr_gene","best_hit_aro","aro","name","hit","arg"]
            for p in arg_files:
                try:
                    with open(p, newline="") as fh:
                        reader = csv.reader(fh, delimiter="\t")
                        header = next(reader, [])
                        hlow   = [h.strip().lower() for h in header]
                        i_sample = hlow.index("sample") if "sample" in hlow else None
                        idxs = [hlow.index(c) for c in col_candidates if c in hlow]
                        for row in reader:
                            if not row: continue
                            s_raw = row[i_sample] if (i_sample is not None and i_sample < len(row)) else ""
                            s     = norm_sample(s_raw) or sm_ref
                            if not same_sample(s, sm_ref): continue
                            hit = any(i < len(row) and looks_like_kpc(row[i]) for i in idxs)
                            if not hit and looks_like_kpc(" ".join(row)):
                                hit = True
                            if hit:
                                log_msg(f"[KPC] HIT: sample={sample} file={os.path.basename(p)}")
                                return True
                except Exception as e:
                    log_msg(f"[KPC] Read fail: {p} ({e})")
                    continue
            log_msg(f"[KPC] NO-HIT for sample={sample} in {len(arg_files)} files")
            return False

        # KPC 없으면: 빈 요약 파일로 종료
        if not sample_has_kpc(wildcards.sample):
            with open(params.skip_log, "a") as sk:
                sk.write(f"{wildcards.sample}\tTETYPER_SKIP_NO_KPC\n")
            pathlib.Path(output.summary).touch()
            return

        # FASTQ 가드
        if (not os.path.exists(input.r1) or os.path.getsize(input.r1) == 0 or
            not os.path.exists(input.r2) or os.path.getsize(input.r2) == 0):
            with open(params.skip_log, "a") as sk:
                sk.write(f"{wildcards.sample}\tTETYPER_EMPTY_FASTQ\n")
            pathlib.Path(output.summary).touch()
            return

        shell(r'''
        set -euo pipefail
        mkdir -p "{OUTPUT_DIR}/5_TETyper" "{params.tmpdir}"

        SPADES_BIN_DIR=$(conda run -n "{params.tety_env}" bash -lc 'dirname "$(which spades.py)"' || true)
        SAMTOOLS_BIN_DIR=$(conda run -n "{params.tety_env}" bash -lc 'dirname "$(which samtools)"' || true)
        BCFTOOLS_BIN_DIR=$(conda run -n "{params.tety_env}" bash -lc 'dirname "$(which bcftools)"' || true)
        BLAST_BIN_DIR=$(conda run -n "{params.tety_env}" bash -lc 'dirname "$(which blastn)"' || true)
        PY_BIN=$(conda run -n "{params.tety_env}" bash -lc 'which python' || true)

        export PATH="$SPADES_BIN_DIR:$SAMTOOLS_BIN_DIR:$BCFTOOLS_BIN_DIR:$BLAST_BIN_DIR:$PATH"
        export TMPDIR="{params.tmpdir}"
        export OMP_NUM_THREADS="{params.threads}"
        export SPADES_MAX_MEMORY="{params.mem_gb}"

        (
          echo "=== TETyper env check ==="
          echo "python: $PY_BIN"
          echo "spades: $(which spades.py || true)"
          echo "samtools: $(which samtools || true)"
          echo "bcftools: $(which bcftools || true)  ($({{ which bcftools >/dev/null 2>&1 && bcftools --version | head -n1; }} || echo 'N/A'))"
          echo "blastn:  $(which blastn || true)"
        ) >> "{log.tlog}" 2>&1 || true

        {params.tety_cmd} \
          --ref "{input.ref}" \
          --fq1 "{input.r1}" \
          --fq2 "{input.r2}" \
          --outprefix "{params.prefix}" \
          --flank_len 5 \
          --struct_profiles "{input.struct}" \
          --snp_profiles "{input.snp}" \
          --show_region "{params.show_region}" >> "{log.tlog}" 2>&1 || true

        # 주의: TETyper가 바로 "{params.prefix}_summary.txt"를 만듭니다.
        # output.summary 와 동일 경로이므로 별도 ln/cp 필요 없음.
        # 다만 혹시 실패했을 경우를 대비해 빈 파일 보장:
        if [ ! -e "{output.summary}" ]; then
          : > "{output.summary}"
        fi
        ''')

rule make_tetyper_summary:
    input:
        txts = expand(f"{OUTPUT_DIR}/5_TETyper/{{sample}}_summary.txt", sample=SAMPLES)
    output:
        summary_json = f"{OUTPUT_DIR}/5_TETyper/tetyper_summary.json"
    run:
        import csv, json, glob, os
        from pathlib import Path

        tdir = Path(f"{OUTPUT_DIR}/5_TETyper")
        rows = []

        for fn in sorted(tdir.glob("*_summary.txt")):
            sample = fn.name.replace("_summary.txt", "")
            # 기본값
            rec = {
                "Sample": sample,
                "Deletions": "N/A",
                "Structural Variant": "N/A",
                "SNPs (Hom)": "N/A",
                "SNPs (Het)": "N/A",
                "Het SNP counts": "N/A",
                "SNP Variant": "N/A",
                "Combined Variant": "N/A",
                "Left flanks": "N/A",
                "Right flanks": "N/A",
                "Left flank counts": "N/A",
                "Right flank counts": "N/A",
                "Region presence": "N/A"  # 예: 7202-8083_presence
            }

            try:
                with open(fn, newline="") as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    for row in reader:
                        # 컬럼 이름들 표준화
                        def get(colnames, default="N/A"):
                            for c in colnames:
                                if c in row and str(row[c]).strip() != "":
                                    return str(row[c]).strip()
                            return default

                        rec["Deletions"]           = get(["Deletions"])
                        rec["Structural Variant"]  = get(["Structural_variant", "Structural Variant"])
                        rec["SNPs (Hom)"]          = get(["SNPs_homozygous", "SNPs (Hom)"])
                        rec["SNPs (Het)"]          = get(["SNPs_heterozygous", "SNPs (Het)"])
                        rec["Het SNP counts"]      = get(["Heterozygous_SNP_counts", "Het SNP counts"])
                        rec["SNP Variant"]         = get(["SNP_variant", "SNP Variant"])
                        rec["Combined Variant"]    = get(["Combined_variant", "Combined Variant"])
                        rec["Left flanks"]         = get(["Left_flanks", "Left flanks"])
                        rec["Right flanks"]        = get(["Right_flanks", "Right flanks"])
                        rec["Left flank counts"]   = get(["Left_flank_counts", "Left flank counts"])
                        rec["Right flank counts"]  = get(["Right_flank_counts", "Right flank counts"])
                        # region presence 컬럼은 이름이 가변적(예: "7202-8083_presence")
                        pres_col = next((c for c in row.keys() if c.endswith("_presence")), None)
                        if pres_col:
                            rec["Region presence"] = str(row[pres_col]).strip()
                        break  # 한 줄만 있으면 충분
            except Exception:
                pass

            # Rmd 호환을 위해 예전 키도 함께 제공(원하면 생략 가능)
            rec["Struct Profile"] = rec["Structural Variant"]
            rec["SNP Profile"]    = rec["SNP Variant"]

            rows.append(rec)

        Path(output.summary_json).write_text(json.dumps(rows, indent=2))

rule tetyper_done:
    input:
        summaries = expand(f"{OUTPUT_DIR}/5_TETyper/{{sample}}_summary.txt", sample=SAMPLES),
        summary_json = f"{OUTPUT_DIR}/5_TETyper/tetyper_summary.json"
    output:
        f"{OUTPUT_DIR}/5_TETyper/tetyper.done.txt"
    shell:
        "touch {output}"
        
        

# ---------- 4) MultiQC ----------
rule multiqc_done_gate:
    input:
        html = f"{OUTPUT_DIR}/multiqc_report.html"
    output:
        f"{OUTPUT_DIR}/multiqc.done.txt"
    shell:
        "touch {output}"


# ---------- 6) 통합 요약 CSV (샘플별 1행) ----------
rule make_master_table:
    input:
        kraken_top = expand(f"{OUTPUT_DIR}/2a_kraken2/{{sample}}.top_species.txt", sample=SAMPLES),
        mlst_done  = f"{OUTPUT_DIR}/2_ST_srst2_out/mlst_srst2.done.txt",
        args_done  = f"{OUTPUT_DIR}/3_ARGs_srst2_out/args_srst2.done.txt",
        plas_done  = f"{OUTPUT_DIR}/4_Plasmid_srst2_out/plasmid_srst2.done.txt",
        tety_json  = f"{OUTPUT_DIR}/5_TETyper/tetyper_summary.json"
    output:
        table = f"{OUTPUT_DIR}/summary_master.csv"
    run:
        import os, glob, json, csv, re
        from collections import defaultdict, Counter

        outdir = os.path.dirname(output.table)
        os.makedirs(outdir, exist_ok=True)

        # 1) Kraken2 top species
        top_species = {}
        for p in input.kraken_top:
            sample = os.path.basename(p).replace(".top_species.txt", "")
            try:
                with open(p) as f:
                    top_species[sample] = f.read().strip() or "unknown"
            except Exception:
                top_species[sample] = "unknown"

        # 2) MLST (유연 파싱)
        mlst_st, mlst_alle = {}, {}
        for sample in SAMPLES:
            r = f"{OUTPUT_DIR}/2_ST_srst2_out/{sample}.mlst.txt"
            st, alle = "NA", "NA"
            if os.path.exists(r) and os.path.getsize(r) > 0:
                try:
                    with open(r) as fh:
                        header = fh.readline().rstrip("\n").split("\t")
                        row    = fh.readline().rstrip("\n").split("\t") if header else []
                    if header and row:
                        hlow = [h.strip().lower() for h in header]
                        hmap = {hlow[i]: i for i in range(len(hlow))}

                        # ST 값
                        if "st" in hmap and hmap["st"] < len(row) and row[hmap["st"]].strip():
                            st = row[hmap["st"]].strip()

                        # 7개 allele만 추출
                        allele_parts = []
                        for col, val in zip(header, row):
                            col_l = col.strip().lower()
                            # allele 이름에 괄호가 있는 "통계 컬럼" 제외하고
                            if "(" in col_l or col_l in ["sample","scheme","st","alleles","flags","st (mlst_check)",
                                                         "coverage","depth","diffs","uncertainty","divergence","length",
                                                         "maxmaf","clusterid","seqid","annotation","mismatches","uncertainty","maxmaf"]:
                                continue
                            # allele 값이 "gene(숫자)" 형태면 그대로 추가
                            if re.match(r"^[A-Za-z0-9_]+$", col) and val.strip().isdigit():
                                allele_parts.append(f"{col}({val.strip()})")
                        if allele_parts:
                            # 상위 7개까지만 (MLST는 항상 7개 locus 기준)
                            alle = ",".join(allele_parts[:7])
                except Exception:
                    pass
            mlst_st[sample], mlst_alle[sample] = st, alle

        # 3) ARGs (CARD) - 배치 통합 결과만
        arg_counts = Counter()
        arg_genes  = defaultdict(set)
        arg_files = sorted(
              glob.glob(f"{OUTPUT_DIR}/3_ARGs_srst2_out/__genes__CARD__genes__CARD_*__results.txt")
            + glob.glob(f"{OUTPUT_DIR}/3_ARGs_srst2_out/__genes__CARD__fullgenes__CARD_*__results.txt")
        )

        for p in arg_files:
            try:
                with open(p) as fh:
                    reader = csv.reader(fh, delimiter="\t")
                    header = next(reader, [])
                    if not header:
                        continue
                    colmap = {h.lower(): i for i, h in enumerate(header)}
                    i_sample = colmap.get("sample", 0)
                    i_gene   = colmap.get("gene", colmap.get("allele"))
                    for row in reader:
                        if not row:
                            continue
                        s = row[i_sample] if i_sample is not None and i_sample < len(row) else None
                        g = row[i_gene]   if i_gene   is not None and i_gene   < len(row) else None
                        if not s or not g:
                            continue
                        core = re.split(r"[|;,\s]", g)[0]
                        arg_genes[s].add(core)
            except Exception:
                pass

        for s, gs in arg_genes.items():
            arg_counts[s] = len(gs)

        # 4) PlasmidFinder (ARGs와 동일 로직)
        plas_counts = Counter()
        plas_reps   = defaultdict(set)

        plas_files = sorted(
              glob.glob(f"{OUTPUT_DIR}/4_Plasmid_srst2_out/__genes__PlasmidFinder__genes__PlasmidFinder*__results.txt")
            + glob.glob(f"{OUTPUT_DIR}/4_Plasmid_srst2_out/__genes__PlasmidFinder__fullgenes__PlasmidFinder*__results.txt")
        )

        for p in plas_files:
            try:
                with open(p) as fh:
                    reader = csv.reader(fh, delimiter="\t")
                    header = next(reader, [])
                    if not header:
                        continue
                    colmap = {h.lower(): i for i, h in enumerate(header)}
                    i_sample = colmap.get("sample", 0)
                    i_gene   = colmap.get("gene", colmap.get("allele"))  # gene 또는 allele

                    for row in reader:
                        if not row:
                            continue
                        s = row[i_sample] if i_sample is not None and i_sample < len(row) else None
                        g = row[i_gene]   if i_gene   is not None and i_gene   < len(row) else None
                        if not s or not g:
                            continue
                        core = re.split(r"[|;,\s]", g)[0]
                        plas_reps[s].add(core)
            except Exception:
                pass

        for s, rs in plas_reps.items():
            plas_counts[s] = len(rs)

        # 5) TETyper 요약(JSON)
        tety_struct = defaultdict(lambda: "NA")
        tety_snp    = defaultdict(lambda: "NA")
        try:
            with open(input.tety_json) as jf:
                items = json.load(jf)
            for it in items:
                s = it.get("Sample")
                if not s:
                    continue
                tety_struct[s] = it.get("Struct Profile", "NA") or "NA"
                tety_snp[s]    = it.get("SNP Profile", "NA") or "NA"
        except Exception:
            pass

        # 6) CSV 작성
        def _top5_str(xset):
            if not xset:
                return "NA"
            xs = sorted(xset)[:5]
            return ",".join(xs)

        with open(output.table, "w", newline="") as out:
            w = csv.writer(out)
            w.writerow([
                "Sample",
                "Kraken2_Top_Species",
                "MLST_ST",
                "MLST_Alleles",
                "ARG_Gene_Count",
                "ARG_Genes_Top5",
                "Plasmid_Rep_Count",
                "Plasmid_Reps_Top5",
                "TETyper_Struct_Profile",
                "TETyper_SNP_Profile",
            ])
            for s in SAMPLES:
                w.writerow([
                    s,
                    top_species.get(s, "unknown"),
                    mlst_st.get(s, "NA"),
                    mlst_alle.get(s, "NA"),
                    arg_counts.get(s, 0),
                    _top5_str(arg_genes.get(s, set())),
                    plas_counts.get(s, 0),
                    _top5_str(plas_reps.get(s, set())),
                    tety_struct.get(s, "NA"),
                    tety_snp.get(s, "NA"),
                ])
                
                
# ---------- 5) Summary (Rmd 호출) ----------
rule summary_report:
    input:
        mlst_table  = f"{OUTPUT_DIR}/2_ST_srst2_out/mlst_master.csv",
        mlst_done   = f"{OUTPUT_DIR}/2_ST_srst2_out/mlst_srst2.done.txt",
        args_done   = f"{OUTPUT_DIR}/3_ARGs_srst2_out/args_srst2.done.txt",
        plas_done   = f"{OUTPUT_DIR}/4_Plasmid_srst2_out/plasmid_srst2.done.txt",
        tetyper_json= f"{OUTPUT_DIR}/5_TETyper/tetyper_summary.json",
        multiqc_gate= f"{OUTPUT_DIR}/multiqc.done.txt",
        master_csv  = f"{OUTPUT_DIR}/summary_master.csv",
    output:
        html = f"{OUTPUT_DIR}/summary_report.html"
    run:
        import os, glob
        summary_html   = os.path.abspath(output.html)
        mlst_master    = os.path.abspath(input.mlst_table)
        tetyper_json   = os.path.abspath(input.tetyper_json)
        master_csv     = os.path.abspath(input.master_csv)


        args_full_list  = sorted(glob.glob(f"{OUTPUT_DIR}/3_ARGs_srst2_out/__genes__CARD__fullgenes__CARD*__results.txt"))
        args_short_list = sorted(glob.glob(f"{OUTPUT_DIR}/3_ARGs_srst2_out/__genes__CARD__genes__CARD*__results.txt"))
        plas_full_list  = sorted(glob.glob(f"{OUTPUT_DIR}/4_Plasmid_srst2_out/__genes__PlasmidFinder__fullgenes__PlasmidFinder*__results.txt"))
        plas_short_list = sorted(glob.glob(f"{OUTPUT_DIR}/4_Plasmid_srst2_out/__genes__PlasmidFinder__genes__PlasmidFinder*__results.txt"))

  
        def to_c(paths):
            return "character(0)" if not paths else "c(" + ",".join(f"'{os.path.abspath(p)}'" for p in paths) + ")"
        args_full_arg  = to_c(args_full_list)
        args_short_arg = to_c(args_short_list)
        plas_full_arg  = to_c(plas_full_list)
        plas_short_arg  = to_c(plas_short_list)

        shell(f"""
          mkdir -p {os.path.dirname(summary_html)}
          Rscript -e "rmarkdown::render('/home/uhlemann/heekuk_path/GoWGS/scripts/20251007_Summary_WGS_tem_v3.Rmd',
            output_file='{summary_html}',
            params=list(
              mlst_master_csv='{mlst_master}',
              args_full_files={args_full_arg},
              args_short_files={args_short_arg},
              plas_full_files={plas_full_arg},
              plas_short_files={plas_short_arg},
              tetyper_json='{tetyper_json}',
              master_csv='{master_csv}'
            ))"
        """)
