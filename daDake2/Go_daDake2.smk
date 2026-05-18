###############################################
# Go_daDake2.smk — DADA2 Amplicon Pipeline
# Types: standard_V3V4 | zymo_V1V2 | illumina_ITS
#
# Flow:
# [ITS only] 1a) cutadapt primer trimming → 3_path.cut/
#            1)  filterAndTrim + quality plots → 3_DADA2_filtered/ or 4_DADA2_filtered/
#            2)  learnErrors → errF.rds, errR.rds
#            3)  dada + merge + chimera removal → seqtab.nochim.rds, seqs.fna
#            4)  assignTaxonomy → tax.rds
#            5)  phyloseq + export CSVs + mapping
#            6)  qiime2 tree
# failed.csv written at startup for tiny/missing FASTQ samples
###############################################

import os, re, glob, csv, datetime
from pathlib import Path

# ── Config ──────────────────────────────────────────────────────────────────
TYPE         = config["type"]
PROJECT_DIRS = [p.strip() for p in config["project_dirs"].split(",") if p.strip()]
DB           = config["db"]
SDB          = config.get("sdb", "")
DATE         = config.get("date", datetime.datetime.now().strftime("%y%m%d"))
MIN_SIZE     = int(config.get("min_fastq_bytes", 10000))
SCRIPTS      = config.get("scripts_dir", os.path.join(os.path.dirname(workflow.snakefile), "scripts"))

# Type-specific parameters
_PARAMS = {
    "standard_V3V4": dict(trimleft_f=20, trimleft_r=21, trunclen_f=250, trunclen_r=220,
                          primer_f="", primer_r="", filt_subdir="3_DADA2_filtered",
                          filt_sfx_r1="_R1_filt.fastq.gz", filt_sfx_r2="_R2_filt.fastq.gz",
                          remove_host=True),
    "zymo_V1V2":     dict(trimleft_f=20, trimleft_r=17, trunclen_f=230, trunclen_r=170,
                          primer_f="", primer_r="", filt_subdir="3_DADA2_filtered",
                          filt_sfx_r1="_R1_filt.fastq.gz", filt_sfx_r2="_R2_filt.fastq.gz",
                          remove_host=True),
    "illumina_ITS":  dict(trimleft_f=0,  trimleft_r=0,  trunclen_f=240, trunclen_r=200,
                          primer_f="GCATCGATGAAGAACGCAG", primer_r="TCCTCCGCTTATTGATATGC",
                          filt_subdir="4_DADA2_filtered",
                          filt_sfx_r1="_1_filt.fastq.gz", filt_sfx_r2="_2_filt.fastq.gz",
                          remove_host=False),
}
P = _PARAMS[TYPE]

def _rc(seq):
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

if TYPE == "illumina_ITS":
    PRIMER_F_RC = _rc(P["primer_f"])
    PRIMER_R_RC = _rc(P["primer_r"])

# ── Preflight: scan each project, write failed.csv ───────────────────────────
def _find_r1s(proj):
    hits = []
    for pat in [f"{proj}/*_L001_R1_001.fastq.gz", f"{proj}/*_R1_001.fastq.gz",
                f"{proj}/*_R1.fastq.gz"]:
        hits += glob.glob(pat)
    return sorted(set(hits))

def _to_r2(r1):
    b = os.path.basename(r1)
    b2 = re.sub(r"(_L001)?_R1(_001)?\.fastq\.gz$",
                lambda m: m.group(0).replace("_R1", "_R2"), b)
    return os.path.join(os.path.dirname(r1), b2)

def _sample_name(r1):
    b = os.path.basename(r1)
    return re.sub(r"(_L001)?_R1(_001)?\.fastq\.gz$", "", b)

PROJECT_SAMPLES = {}   # proj -> [passing sample names]

for _proj in PROJECT_DIRS:
    _dada2_dir = f"{_proj}_dada2"
    os.makedirs(f"{_dada2_dir}/1_out", exist_ok=True)
    os.makedirs(f"{_dada2_dir}/2_rds", exist_ok=True)

    _r1s = _find_r1s(_proj)
    _failed = []
    _passing = []

    for _r1 in _r1s:
        _sn = _sample_name(_r1)
        _r2 = _to_r2(_r1)
        _s1 = os.path.getsize(_r1) if os.path.exists(_r1) else 0
        _s2 = os.path.getsize(_r2) if os.path.exists(_r2) else 0
        _reason = ""
        if not os.path.exists(_r2):         _reason = "R2_missing"
        elif _s1 < MIN_SIZE:                _reason = f"R1_too_small({_s1}B)"
        elif _s2 < MIN_SIZE:                _reason = f"R2_too_small({_s2}B)"
        if _reason:
            _failed.append(dict(sample=_sn, r1=_r1, r2=_r2,
                                r1_bytes=_s1, r2_bytes=_s2, reason=_reason))
        else:
            _passing.append(_sn)

    with open(f"{_dada2_dir}/failed.csv", "w", newline="") as _fh:
        _w = csv.DictWriter(_fh, fieldnames=["sample","r1","r2","r1_bytes","r2_bytes","reason"])
        _w.writeheader()
        _w.writerows(_failed)

    if _failed:
        print(f"[daDake2][{_proj}] {len(_failed)} sample(s) skipped → {_dada2_dir}/failed.csv")
    if not _passing:
        raise ValueError(f"[FATAL] No valid samples in {_proj}. See {_dada2_dir}/failed.csv")

    PROJECT_SAMPLES[_proj] = _passing

# ── Helper: build rule output paths ─────────────────────────────────────────
def _out(proj, *parts):
    return os.path.join(f"{proj}_dada2", *parts)

# ── Final targets ────────────────────────────────────────────────────────────
rule all:
    input:
        expand(f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.asvTable.csv",  proj=PROJECT_DIRS),
        expand(f"{{proj}}_dada2/2_rds/ps.{{proj}}.{DATE}.rds",        proj=PROJECT_DIRS),
        expand(f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.seqs.fna_tree/exported-tree/tree.nwk",
               proj=PROJECT_DIRS),

# ══════════════════════════════════════════════════════════════════════════════
# ITS ONLY: cutadapt primer trimming
# ══════════════════════════════════════════════════════════════════════════════
if TYPE == "illumina_ITS":
    rule cutadapt_primers:
        input:  fastq_dir = lambda wc: wc.proj
        output: done = "{proj}_dada2/3_path.cut/.done"
        log:    "{proj}_dada2/logs/cutadapt.log"
        params:
            cut_dir   = "{proj}_dada2/3_path.cut",
            primer_f  = P["primer_f"],
            primer_r  = P["primer_r"],
            primer_frc = PRIMER_F_RC,
            primer_rrc = PRIMER_R_RC,
        shell:
            r"""
            set -euo pipefail
            mkdir -p "{params.cut_dir}" "$(dirname {log})"

            shopt -s nullglob
            R1S=( {input.fastq_dir}/*_L001_R1_001.fastq.gz \
                  {input.fastq_dir}/*_R1_001.fastq.gz )
            shopt -u nullglob

            for R1 in "${{R1S[@]}}"; do
                R2="${{R1/_R1_/_R2_}}"
                R2="${{R2/_L001_R1_001.fastq.gz/_L001_R2_001.fastq.gz}}"
                R2="${{R2/_R1_001.fastq.gz/_R2_001.fastq.gz}}"
                [ -f "$R2" ] || {{ echo "[cutadapt] SKIP $R1 — R2 missing" >> {log}; continue; }}
                OUT_R1="{params.cut_dir}/$(basename "$R1")"
                OUT_R2="{params.cut_dir}/$(basename "$R2")"
                cutadapt \
                    -g {params.primer_f}  -a {params.primer_rrc} \
                    -G {params.primer_r}  -A {params.primer_frc} \
                    -n 2 --discard-untrimmed \
                    -o "$OUT_R1" -p "$OUT_R2" \
                    "$R1" "$R2" >> {log} 2>&1
            done
            touch {output.done}
            """

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: filterAndTrim + quality plots
# ══════════════════════════════════════════════════════════════════════════════
if TYPE == "illumina_ITS":
    rule filter_trim:
        input:
            cut_done = "{proj}_dada2/3_path.cut/.done",
        output:
            done     = f"{{proj}}_dada2/{P['filt_subdir']}/.done",
            qc_raw   = f"{{proj}}_dada2/{{proj}}.{DATE}.qualityProfiles.pdf",
            qc_filt  = f"{{proj}}_dada2/{{proj}}.{DATE}.qualityProfiles.filt.pdf",
            filt_rds = f"{{proj}}_dada2/2_rds/filter_stats.rds",
        log: "{proj}_dada2/logs/filter_trim.log"
        params:
            dada2_dir  = "{proj}_dada2",
            fastq_dir  = "{proj}_dada2/3_path.cut",
            project    = "{proj}",
            filt_sub   = P["filt_subdir"],
            trimleft_f = P["trimleft_f"],
            trimleft_r = P["trimleft_r"],
            trunclen_f = P["trunclen_f"],
            trunclen_r = P["trunclen_r"],
            type_      = TYPE,
            date       = DATE,
            scripts    = SCRIPTS,
            failed_csv = "{proj}_dada2/failed.csv",
        shell:
            r"""
            set -euo pipefail
            mkdir -p "{params.dada2_dir}/logs"
            Rscript {params.scripts}/01_filter_trim.R \
                --dada2_dir  "{params.dada2_dir}" \
                --fastq_dir  "{params.fastq_dir}" \
                --project    "{params.project}" \
                --filt_sub   "{params.filt_sub}" \
                --trimleft_f {params.trimleft_f} \
                --trimleft_r {params.trimleft_r} \
                --trunclen_f {params.trunclen_f} \
                --trunclen_r {params.trunclen_r} \
                --type       "{params.type_}" \
                --date       "{params.date}" \
                --failed_csv "{params.failed_csv}" \
                > {log} 2>&1
            """
else:
    rule filter_trim:
        input:
            fastq_dir = lambda wc: wc.proj,
        output:
            done     = f"{{proj}}_dada2/{P['filt_subdir']}/.done",
            qc_raw   = f"{{proj}}_dada2/{{proj}}.{DATE}.qualityProfiles.pdf",
            qc_filt  = f"{{proj}}_dada2/{{proj}}.{DATE}.qualityProfiles.filt.pdf",
            filt_rds = f"{{proj}}_dada2/2_rds/filter_stats.rds",
        log: "{proj}_dada2/logs/filter_trim.log"
        params:
            dada2_dir  = "{proj}_dada2",
            fastq_dir  = "{proj}",
            project    = "{proj}",
            filt_sub   = P["filt_subdir"],
            trimleft_f = P["trimleft_f"],
            trimleft_r = P["trimleft_r"],
            trunclen_f = P["trunclen_f"],
            trunclen_r = P["trunclen_r"],
            type_      = TYPE,
            date       = DATE,
            scripts    = SCRIPTS,
            failed_csv = "{proj}_dada2/failed.csv",
        shell:
            r"""
            set -euo pipefail
            mkdir -p "{params.dada2_dir}/logs"
            Rscript {params.scripts}/01_filter_trim.R \
                --dada2_dir  "{params.dada2_dir}" \
                --fastq_dir  "{params.fastq_dir}" \
                --project    "{params.project}" \
                --filt_sub   "{params.filt_sub}" \
                --trimleft_f {params.trimleft_f} \
                --trimleft_r {params.trimleft_r} \
                --trunclen_f {params.trunclen_f} \
                --trunclen_r {params.trunclen_r} \
                --type       "{params.type_}" \
                --date       "{params.date}" \
                --failed_csv "{params.failed_csv}" \
                > {log} 2>&1
            """

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: learnErrors
# ══════════════════════════════════════════════════════════════════════════════
rule learn_errors:
    input:
        filt_done = f"{{proj}}_dada2/{P['filt_subdir']}/.done",
    output:
        errF    = "{proj}_dada2/2_rds/errF.rds",
        errR    = "{proj}_dada2/2_rds/errR.rds",
        pdf_F   = f"{{proj}}_dada2/{{proj}}.{DATE}.splotErrors.errF1.pdf",
        pdf_R   = f"{{proj}}_dada2/{{proj}}.{DATE}.splotErrors.errF2.pdf",
    log: "{proj}_dada2/logs/learn_errors.log"
    params:
        dada2_dir = "{proj}_dada2",
        project   = "{proj}",
        filt_sub  = P["filt_subdir"],
        type_     = TYPE,
        date      = DATE,
        scripts   = SCRIPTS,
        failed_csv = "{proj}_dada2/failed.csv",
    shell:
        r"""
        set -euo pipefail
        Rscript {params.scripts}/02_learn_errors.R \
            --dada2_dir  "{params.dada2_dir}" \
            --project    "{params.project}" \
            --filt_sub   "{params.filt_sub}" \
            --type       "{params.type_}" \
            --date       "{params.date}" \
            --failed_csv "{params.failed_csv}" \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: dada + merge + chimera removal → seqtab.nochim + seqs.fna
# ══════════════════════════════════════════════════════════════════════════════
rule denoise_merge:
    input:
        errF      = "{proj}_dada2/2_rds/errF.rds",
        errR      = "{proj}_dada2/2_rds/errR.rds",
        filt_done = f"{{proj}}_dada2/{P['filt_subdir']}/.done",
    output:
        seqtab  = f"{{proj}}_dada2/2_rds/seqtab.nochim.{{proj}}.{DATE}.rds",
        seqsfna = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.seqs.fna",
        dada_rds = f"{{proj}}_dada2/2_rds/dada_stats.rds",
    log: "{proj}_dada2/logs/denoise_merge.log"
    params:
        dada2_dir = "{proj}_dada2",
        project   = "{proj}",
        filt_sub  = P["filt_subdir"],
        type_     = TYPE,
        date      = DATE,
        scripts   = SCRIPTS,
        failed_csv = "{proj}_dada2/failed.csv",
    shell:
        r"""
        set -euo pipefail
        Rscript {params.scripts}/03_denoise_merge.R \
            --dada2_dir  "{params.dada2_dir}" \
            --project    "{params.project}" \
            --filt_sub   "{params.filt_sub}" \
            --type       "{params.type_}" \
            --date       "{params.date}" \
            --failed_csv "{params.failed_csv}" \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: assignTaxonomy (+ addSpecies for 16S)
# ══════════════════════════════════════════════════════════════════════════════
rule assign_taxonomy:
    input:
        seqtab = f"{{proj}}_dada2/2_rds/seqtab.nochim.{{proj}}.{DATE}.rds",
    output:
        tax = f"{{proj}}_dada2/2_rds/tax.{{proj}}.{DATE}.rds",
    log: "{proj}_dada2/logs/assign_taxonomy.log"
    params:
        dada2_dir   = "{proj}_dada2",
        project     = "{proj}",
        date        = DATE,
        db          = DB,
        sdb         = SDB,
        remove_host = str(P["remove_host"]).upper(),
        type_       = TYPE,
        scripts     = SCRIPTS,
    shell:
        r"""
        set -euo pipefail
        Rscript {params.scripts}/04_taxonomy.R \
            --dada2_dir   "{params.dada2_dir}" \
            --project     "{params.project}" \
            --date        "{params.date}" \
            --db          "{params.db}" \
            --sdb         "{params.sdb}" \
            --remove_host {params.remove_host} \
            --type        "{params.type_}" \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: phyloseq object + export all CSVs + mapping files
# ══════════════════════════════════════════════════════════════════════════════
rule export_tables:
    input:
        seqtab   = f"{{proj}}_dada2/2_rds/seqtab.nochim.{{proj}}.{DATE}.rds",
        tax      = f"{{proj}}_dada2/2_rds/tax.{{proj}}.{DATE}.rds",
        dada_rds = f"{{proj}}_dada2/2_rds/dada_stats.rds",
        filt_rds = f"{{proj}}_dada2/2_rds/filter_stats.rds",
    output:
        ps        = f"{{proj}}_dada2/2_rds/ps.{{proj}}.{DATE}.rds",
        track_csv = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.track.csv",
        asv_csv   = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.asv.csv",
        tax_csv   = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.tax.csv",
        asv_table = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.asvTable.csv",
    log: "{proj}_dada2/logs/export_tables.log"
    params:
        dada2_dir = "{proj}_dada2",
        project   = "{proj}",
        date      = DATE,
        type_     = TYPE,
        scripts   = SCRIPTS,
    shell:
        r"""
        set -euo pipefail
        Rscript {params.scripts}/05_export.R \
            --dada2_dir "{params.dada2_dir}" \
            --project   "{params.project}" \
            --date      "{params.date}" \
            --type      "{params.type_}" \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: qiime2 phylogenetic tree
# ══════════════════════════════════════════════════════════════════════════════
rule qiime_tree:
    input:
        fna = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.seqs.fna",
    output:
        nwk = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.seqs.fna_tree/exported-tree/tree.nwk",
    log: "{proj}_dada2/logs/qiime_tree.log"
    params:
        fna_name = f"{{proj}}.{DATE}.seqs.fna",
        out_dir  = f"{{proj}}_dada2/1_out",
    shell:
        r"""
        set -euo pipefail
        cd "{params.out_dir}"

        FNA="{params.fna_name}"
        TREE_DIR="${{FNA}}_tree"
        mkdir -p "$TREE_DIR"

        conda run -n qiime2 qiime tools import \
            --input-path  "$FNA" \
            --output-path "$TREE_DIR/${{FNA%.fna}}.qza" \
            --type        'FeatureData[Sequence]' >> {log} 2>&1

        conda run -n qiime2 qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences      "$TREE_DIR/${{FNA%.fna}}.qza" \
            --o-alignment      "$TREE_DIR/${{FNA%.fna}}_aligned-rep-seqs.qza" \
            --o-masked-alignment "$TREE_DIR/${{FNA%.fna}}_masked-aligned-rep-seqs.qza" \
            --o-tree           "$TREE_DIR/${{FNA%.fna}}_unrooted-tree.qza" \
            --o-rooted-tree    "$TREE_DIR/${{FNA%.fna}}_rooted-tree.qza" >> {log} 2>&1

        conda run -n qiime2 qiime tools export \
            --input-path  "$TREE_DIR/${{FNA%.fna}}_rooted-tree.qza" \
            --output-path "$TREE_DIR/exported-tree" >> {log} 2>&1
        """
