# Humann

Snakemake + Docker rewrite of the HUMAnN3 execution block used in Uhlemann Lab shotgun workflows.

## Files

- `Go_Humann_V1.smk`: main Snakemake workflow
- `Dockerfile`: dependencies for Snakemake and HUMAnN3
- `Go_Humann.sh`: Docker wrapper
- `scripts/humann_masterlog.py`: merges per-sample HUMAnN logs into `humann3_log.txt`

## Workstation Layout

After `Go_toWorkstation.sh Humann`, files are placed as:

```text
/home/uhlemann*/heekuk_path/
  Go_Humann.sh
  Go_Humann.smk
  scripts/
    humann_masterlog.py
  docker/Humann/
    Dockerfile
```

## Build

```bash
cd /home/uhlemann/heekuk_path/docker/Humann
docker build --network=host -t humann:1.0 .
```

## Direct Snakemake Run

```bash
fastq_dir="host_filtered_fastq"
output_dir="humann_run"
chocophlan="/media/uhlemann/core4/DB/humann_db/humann3/chocophlan"
uniref="/media/uhlemann/core4/DB/humann_db/humann3/uniref"
metaphlan_db="/media/uhlemann/core4/DB/humann_db/metaphlan4"
metaphlan_index="mpa_vJan25_CHOCOPhlAnSGB_202503"

snakemake --snakefile /home/uhlemann/heekuk_path/Go_Humann.smk \
  --config \
  fastq_dir="$fastq_dir" \
  output_dir="$output_dir" \
  nucleotide_db="$chocophlan" \
  protein_db="$uniref" \
  metaphlan_db="$metaphlan_db" \
  metaphlan_index="$metaphlan_index" \
  run_musicc=true \
  humann_threads=4 \
  --cores 8 --jobs 4 \
  --latency-wait 60 --rerun-incomplete
```

## Docker Wrapper Run

```bash
./Go_Humann.sh \
  -i /path/to/host_filtered_fastq \
  -o humann_run \
  -n /media/uhlemann/core4/DB/humann_db/humann3/chocophlan \
  -p /media/uhlemann/core4/DB/humann_db/humann3/uniref \
  -b /media/uhlemann/core4/DB/humann_db/metaphlan4 \
  -I mpa_vJan25_CHOCOPhlAnSGB_202503 \
  -c 8 \
  -j 4 \
  -t 4 \
  --run-musicc \
  -K
```

## Output Layout

```text
output/
  1_humann3_out/
    <sample>_genefamilies.tsv
    <sample>_pathabundance.tsv
    <sample>_pathcoverage.tsv
    <sample>.humann.done
  2_humann3_final_out/
    merged_genefamilies.txt
    merged_pathabundance.txt
    merged_pathcoverage.txt
  3_MetaPhlAn_bug_list/
    <sample>_metaphlan_bugs_list.tsv
  4_kegg-orthology/
    merged_genefamilies_cpm.txt
    merged_genefamilies_uniref90_rxn_cpm.txt
    merged_genefamilies_uniref90_rxn_musicc.txt
    merged_genefamilies_uniref90_rxn_kegg-orthology_cpm.txt
    stratified_out/
      merged_genefamilies_uniref90_rxn_kegg-orthology_cpm_unstratified.txt
      merged_genefamilies_uniref90_rxn_kegg-orthology_cpm_unstratified_filtered.txt
  5_pathabundance_stratified_out/
    merged_pathabundance_unstratified.txt
  6_logs/
    *.log
  7_intermediate/
    <sample>.merged.fastq.gz
    <sample>.input.done
  humann3_log.txt
```

## Useful Configs

```bash
--config \
  fastq_dir=host_filtered_fastq \
  output_dir=humann_run \
  nucleotide_db=/path/to/chocophlan \
  protein_db=/path/to/uniref \
  metaphlan_db=/path/to/metaphlan_db \
  metaphlan_index=mpa_vJan25_CHOCOPhlAnSGB_202503 \
  run_musicc=false \
  humann_threads=4 \
  run_gene_norm=true \
  run_path_split=true \
  run_pathcoverage_merge=true
```

## Notes

- Recommended MetaPhlAn input is `metaphlan_db=<directory>` plus `metaphlan_index=<basename without .pkl>`.
- `Go_Humann.sh` still accepts a `.pkl` path via `-b` and auto-converts it to directory + index for backward compatibility.
- Input discovery is intentionally simple: every `*.fastq`, `*.fq`, `*.fastq.gz`, `*.fq.gz` file is treated as one HUMAnN input unit.
- Single FASTQ inputs are processed as one HUMAnN input unit.
- If inputs are split as `*_R1.fastq.gz` and `*_R2.fastq.gz`, including `*_R1_nohuman.fastq.gz` and `*_R2_nohuman.fastq.gz`, they are merged per sample into `7_intermediate/<sample>.merged.fastq.gz` and then run as one HUMAnN job.
- Per-sample `.humann.done` markers are added so interrupted runs resume cleanly.
- `MetaPhlAn_bug_list` files are extracted from `{sample}_humann_temp/` before that temp directory is removed.
- Internal output directories are numbered in run order: `1_` through `7_`.
- `--run-musicc` enables MUSiCC correction on the regrouped KO table before `humann_rename_table`.
- `--skip-gene-norm` is useful if only raw HUMAnN tables are needed.
- `--skip-path-split` is useful if pathway split tables are not needed for the current project.
- By default the workflow uses standard HUMAnN normalization only. MUSiCC is optional.
