###############################################
# Go_daDake2.smk — DADA2 Amplicon Pipeline
# Types: standard_V3V4_600 | standard_V3V4_250 | standard_V4 | public_V4_notrim | standard_V1V2 | zymo_V1V2 | illumina_ITS
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
# primer_len_f/r = true biological primer length (no early-cycle buffer);
# amplicon_len   = expected insert length WITHOUT primers. V3V4=403 is measured
#                   from real merged-sequence length distributions on this lab's
#                   own V3V4 runs (20260820_16S_Deb denoise_merge.log: modal
#                   401-405bp) -- the literature figure (~426bp, Klindworth 2013)
#                   overshot by ~24bp, which made FIGARO read needlessly deep into
#                   low-quality tail and tanked mergePairs() success on some
#                   samples (maxMismatch=0 has zero tolerance). V4/V1V2 values are
#                   still literature estimates (Caporaso/EMP, rough V1V2 approx) --
#                   override with -a if a kit spec or your own data differs.
# Both are only used by the FIGARO auto-trunclen path (-A); the default
# trimleft_f/r + trunclen_f/r below are untouched and remain the normal path.
_PARAMS = {
    "standard_V3V4": dict(trimleft_f=20, trimleft_r=21, trunclen_f=240, trunclen_r=240,
                          primer_f="", primer_r="", filt_subdir="3_DADA2_filtered",
                          filt_sfx_r1="_R1_filt.fastq.gz", filt_sfx_r2="_R2_filt.fastq.gz",
                          remove_host=False,
                          primer_len_f=17, primer_len_r=21, amplicon_len=403),
    "standard_V3V4_600": dict(trimleft_f=20, trimleft_r=21, trunclen_f=240, trunclen_r=240,
                               primer_f="", primer_r="", filt_subdir="3_DADA2_filtered",
                               filt_sfx_r1="_R1_filt.fastq.gz", filt_sfx_r2="_R2_filt.fastq.gz",
                               remove_host=False,
                               primer_len_f=17, primer_len_r=21, amplicon_len=403),
    "standard_V3V4_250": dict(trimleft_f=20, trimleft_r=21, trunclen_f=230, trunclen_r=229,
                               primer_f="", primer_r="", filt_subdir="3_DADA2_filtered",
                               filt_sfx_r1="_R1_filt.fastq.gz", filt_sfx_r2="_R2_filt.fastq.gz",
                               remove_host=False,
                               primer_len_f=17, primer_len_r=21, amplicon_len=403),
    "standard_V4":     dict(trimleft_f=19, trimleft_r=20, trunclen_f=230, trunclen_r=220,
                             primer_f="", primer_r="", filt_subdir="3_DADA2_filtered",
                             filt_sfx_r1="_R1_filt.fastq.gz", filt_sfx_r2="_R2_filt.fastq.gz",
                             remove_host=False,
                             primer_len_f=19, primer_len_r=20, amplicon_len=253),
    # V4 (515F/806R), primers already stripped at deposit (common for public ENA/SRA
    # 16S runs) -- trimleft=0 unlike standard_V4, which assumes primers are still in
    # the read. trunclen chosen from the actual quality profile (mean Q>=30 to ~240bp
    # fwd, ~200bp rev before the tail drop), not copied from another type.
    "public_V4_notrim": dict(trimleft_f=0, trimleft_r=0, trunclen_f=240, trunclen_r=200,
                              primer_f="", primer_r="", filt_subdir="3_DADA2_filtered",
                              filt_sfx_r1="_R1_filt.fastq.gz", filt_sfx_r2="_R2_filt.fastq.gz",
                              remove_host=False,
                              primer_len_f=0, primer_len_r=0, amplicon_len=253),
    # V1-V2 (27F/338R, 20bp/19bp), not the Zymo kit primer set -- trimleft_f=23
    # adds the same +3bp early-cycle buffer used for standard_V3V4 (17bp primer ->
    # trimleft 20); trimleft_r=19 matches 338R exactly (no buffer needed there,
    # mirroring standard_V3V4's R). trunclen reused from zymo_V1V2 (same V1-V2
    # amplicon length).
    "standard_V1V2": dict(trimleft_f=23, trimleft_r=19, trunclen_f=230, trunclen_r=170,
                          primer_f="", primer_r="", filt_subdir="3_DADA2_filtered",
                          filt_sfx_r1="_R1_filt.fastq.gz", filt_sfx_r2="_R2_filt.fastq.gz",
                          remove_host=False,
                          primer_len_f=20, primer_len_r=19, amplicon_len=300),
    "zymo_V1V2":     dict(trimleft_f=20, trimleft_r=17, trunclen_f=230, trunclen_r=170,
                          primer_f="", primer_r="", filt_subdir="3_DADA2_filtered",
                          filt_sfx_r1="_R1_filt.fastq.gz", filt_sfx_r2="_R2_filt.fastq.gz",
                          remove_host=False,
                          primer_len_f=20, primer_len_r=17, amplicon_len=300),
    "illumina_ITS":  dict(trimleft_f=0,  trimleft_r=0,  trunclen_f=240, trunclen_r=200,
                          primer_f="GCATCGATGAAGAACGCAG", primer_r="TCCTCCGCTTATTGATATGC",
                          filt_subdir="4_DADA2_filtered",
                          filt_sfx_r1="_1_filt.fastq.gz", filt_sfx_r2="_2_filt.fastq.gz",
                          remove_host=False),
}
P = _PARAMS[TYPE]

# ── FIGARO auto-trunclen (opt-in, -A) ─────────────────────────────────────────
# Computes per-project truncLen from the actual quality profile instead of the
# fixed defaults above. trimLeft always uses the true primer length (no manual
# buffer) since FIGARO's optimization already accounts for real quality decay.
AUTO_TRUNCLEN = str(config.get("auto_trunclen", "false")).lower() == "true"
if AUTO_TRUNCLEN:
    if TYPE == "illumina_ITS":
        raise ValueError("[FATAL] -A (FIGARO auto-trunclen) is not supported for illumina_ITS.")
    AMPLICON_LEN  = int(config.get("amplicon_len", 0)) or P["amplicon_len"]
    # FIGARO's own default (-m 20) picks truncLen combos that hug the minimum
    # overlap exactly -- fine in theory, but mergePairs()'s real minOverlap=12
    # leaves zero margin for indels/quality noise, so real runs can merge
    # almost nothing even though filter retention looks great. Swept on
    # 20260820_16S_Deb (21 samples, standard_V3V4, amplicon_len=403):
    #   20, 30 -> 8/21 samples (DM-17/18/19/20/21/22/23/Pos-6, same set every
    #             time) stuck at merged=0
    #   50     -> all 21 merge, but total final (nonchim) reads = 43,197,
    #             only 82% of the old fixed standard_V3V4 defaults (240/240,
    #             ~36bp real overlap) on the same data (52,541)
    #   40     -> all 21 merge AND total nonchim = 68,037 -- 130% of the old
    #             fixed defaults. Best of the sweep; keep as default.
    MIN_OVERLAP   = int(config.get("min_overlap", 0)) or 40
    FIGARO_ENV    = "figaro_env"
    FIGARO_REPO   = os.path.join(os.path.dirname(workflow.snakefile), ".figaro_src")
    FIGARO_SCRIPT = os.path.join(FIGARO_REPO, "figaro", "figaro.py")
    FIGARO_RUNNER = os.path.join(SCRIPTS, "figaro_compat.py")

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

# ══════════════════════════════════════════════════════════════════════════════
# FIGARO auto-trunclen: self-installing setup + per-project optimization
# ══════════════════════════════════════════════════════════════════════════════
if AUTO_TRUNCLEN:
    def _passing_fastqs(wc):
        """Return complete FASTQ pairs that passed the shared preflight scan."""
        passing = set(PROJECT_SAMPLES[wc.proj])
        files = []
        for r1 in _find_r1s(wc.proj):
            if _sample_name(r1) in passing:
                files.extend([r1, _to_r2(r1)])
        return files

    rule figaro_setup:
        output:
            marker = os.path.join(FIGARO_REPO, ".installed"),
        log: "logs/figaro_setup.log"
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {log})"
            if [ ! -f "{FIGARO_SCRIPT}" ]; then
                git clone --depth 1 https://github.com/Zymo-Research/figaro.git "{FIGARO_REPO}" >> {log} 2>&1
            fi
            if ! conda env list | awk '{{print $1}}' | grep -qx "{FIGARO_ENV}"; then
                conda create -y -c conda-forge -n {FIGARO_ENV} python=3.10 numpy scipy matplotlib >> {log} 2>&1
            fi
            touch {output.marker}
            """

    rule figaro_stage_fastqs:
        input:
            fastqs = _passing_fastqs,
        output:
            fastq_dir = directory("{proj}_dada2/figaro_input_fixedlen"),
        script:
            os.path.join(SCRIPTS, "prepare_figaro_fastqs.py")

    checkpoint figaro_optimize:
        input:
            fastq_dir = rules.figaro_stage_fastqs.output.fastq_dir,
            setup     = os.path.join(FIGARO_REPO, ".installed"),
        output:
            json = "{proj}_dada2/figaro/trimParameters.json",
        log: "{proj}_dada2/logs/figaro.log"
        params:
            outdir       = "{proj}_dada2/figaro",
            amplicon_len = AMPLICON_LEN,
            primer_f_len = P["primer_len_f"],
            primer_r_len = P["primer_len_r"],
            min_overlap  = MIN_OVERLAP,
        shell:
            r"""
            set -euo pipefail
            mkdir -p "{params.outdir}" "$(dirname {log})"
            conda run -n {FIGARO_ENV} python3 "{FIGARO_RUNNER}" "{FIGARO_SCRIPT}" \
                -i "{input.fastq_dir}" \
                -o "{params.outdir}" \
                -a {params.amplicon_len} \
                -f {params.primer_f_len} \
                -r {params.primer_r_len} \
                -m {params.min_overlap} \
                -F illumina \
                > {log} 2>&1
            """

    def _figaro_trunclen(wc):
        """Parse FIGARO's best-scoring candidate into DADA2 truncLen values.
        trimPosition is a raw-read cycle number (primer included, see
        Zymo-Research/figaro trimParameterPrediction.py) and DADA2's own
        truncLen is ALSO measured from the raw read, before trimLeft is
        applied -- dada2::filterAndTrim: "if both truncLen and trimLeft are
        provided, filtered reads will have length truncLen-trimLeft". So
        trimPosition maps directly to truncLen with NO subtraction; DADA2
        does the primer-length subtraction itself via trimLeft. (Previously
        this subtracted primer_len_f/r here too -- a double subtraction that
        silently shrank the real overlap by primer_len_f+primer_len_r and
        collapsed mergePairs() success on several samples.)"""
        import json
        json_path = checkpoints.figaro_optimize.get(proj=wc.proj).output.json
        with open(json_path) as fh:
            best = json.load(fh)[0]
        fwd_pos, rev_pos = best["trimPosition"]
        return dict(trunclen_f=fwd_pos, trunclen_r=rev_pos)

# ── Final targets ────────────────────────────────────────────────────────────
rule all:
    default_target: True
    input:
        expand(f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.asvTable.csv",  proj=PROJECT_DIRS),
        expand(f"{{proj}}_dada2/2_rds/ps.{{proj}}.{DATE}.rds",        proj=PROJECT_DIRS),
        expand(f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.seqs.fna_tree/exported-tree/tree.nwk",
               proj=PROJECT_DIRS),
        expand("zip_projects/{proj}_dada2.zip", proj=PROJECT_DIRS),

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
                  {input.fastq_dir}/*_R1_001.fastq.gz \
                  {input.fastq_dir}/*_R1.fastq.gz )
            shopt -u nullglob

            for R1 in "${{R1S[@]}}"; do
                R2="${{R1/_R1_/_R2_}}"
                R2="${{R2/_L001_R1_001.fastq.gz/_L001_R2_001.fastq.gz}}"
                R2="${{R2/_R1_001.fastq.gz/_R2_001.fastq.gz}}"
                R2="${{R2/_R1.fastq.gz/_R2.fastq.gz}}"
                [ -f "$R2" ] || {{ echo "[cutadapt] SKIP $R1 — R2 missing" >> {log}; continue; }}
                OUT_R1="{params.cut_dir}/$(basename "$R1")"
                OUT_R2="{params.cut_dir}/$(basename "$R2")"
                cutadapt \
                    -g {params.primer_f}  -a {params.primer_rrc} \
                    -G {params.primer_r}  -A {params.primer_frc} \
                    -n 2 \
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
            fastq_dir   = lambda wc: wc.proj,
            figaro_json = (lambda wc: checkpoints.figaro_optimize.get(proj=wc.proj).output.json) \
                          if AUTO_TRUNCLEN else [],
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
            trimleft_f = P["primer_len_f"] if AUTO_TRUNCLEN else P["trimleft_f"],
            trimleft_r = P["primer_len_r"] if AUTO_TRUNCLEN else P["trimleft_r"],
            trunclen_f = (lambda wc: _figaro_trunclen(wc)["trunclen_f"]) if AUTO_TRUNCLEN else P["trunclen_f"],
            trunclen_r = (lambda wc: _figaro_trunclen(wc)["trunclen_r"]) if AUTO_TRUNCLEN else P["trunclen_r"],
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
        fna       = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.seqs.fna",
        asv_table = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.asvTable.csv",
    output:
        nwk = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.seqs.fna_tree/exported-tree/tree.nwk",
    log: "{proj}_dada2/logs/qiime_tree.log"
    params:
        fna      = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.seqs.fna",
        tree_dir = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.seqs.fna_tree",
        stem     = f"{{proj}}.{DATE}.seqs",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})" "{params.tree_dir}"

        conda run -n qiime2 qiime tools import \
            --input-path  "{params.fna}" \
            --output-path "{params.tree_dir}/{params.stem}.qza" \
            --type        'FeatureData[Sequence]' >> {log} 2>&1

        conda run -n qiime2 qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences        "{params.tree_dir}/{params.stem}.qza" \
            --o-alignment        "{params.tree_dir}/{params.stem}_aligned-rep-seqs.qza" \
            --o-masked-alignment "{params.tree_dir}/{params.stem}_masked-aligned-rep-seqs.qza" \
            --o-tree             "{params.tree_dir}/{params.stem}_unrooted-tree.qza" \
            --o-rooted-tree      "{params.tree_dir}/{params.stem}_rooted-tree.qza" >> {log} 2>&1

        conda run -n qiime2 qiime tools export \
            --input-path  "{params.tree_dir}/{params.stem}_rooted-tree.qza" \
            --output-path "{params.tree_dir}/exported-tree" >> {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# Step 7: zip results
# ══════════════════════════════════════════════════════════════════════════════
rule create_zip:
    input:
        nwk = f"{{proj}}_dada2/1_out/{{proj}}.{DATE}.seqs.fna_tree/exported-tree/tree.nwk",
        ps  = f"{{proj}}_dada2/2_rds/ps.{{proj}}.{DATE}.rds",
    output:
        zip = "zip_projects/{proj}_dada2.zip",
    shell:
        r"""
        set -euo pipefail
        mkdir -p zip_projects
        zip -s 3700m -r {output.zip} {wildcards.proj}_dada2 >/dev/null
        echo "[create_zip] {wildcards.proj}_dada2.zip done"
        """
