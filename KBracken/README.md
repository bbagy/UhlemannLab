# KBracken

Snakemake + Docker rewrite of `Gobracken2_V4.pl` for Kraken2 + Bracken profiling.

## Files

- `Go_KBracken_V1.smk`: main Snakemake workflow
- `Dockerfile`: dependencies for Snakemake, Kraken2, and Bracken
- `Go_KBracken.sh`: Docker wrapper
- `scripts/merge_mpa_tables.py`: fixed local replacement for `merge_metaphlan_tables.py`
- `scripts/bracken_to_mpa.py`: Bracken table to MPA-style table
- `scripts/fill_bracken_taxonomy.py`: adds Kraken2 taxonomy ranks to merged Bracken table
- `scripts/kraken_masterlog.py`: creates `kraken2_log.txt`

## Workstation Layout

After `Go_toWorkstation.sh KBracken`, files are placed as:

```text
/home/uhlemann*/heekuk_path/
  Go_KBracken.sh
  Go_KBracken.smk
  scripts/
    bracken_to_mpa.py
    fill_bracken_taxonomy.py
    kraken_masterlog.py
    merge_mpa_tables.py
  docker/KBracken/
    Dockerfile
```

## Build

```bash
cd /home/uhlemann/heekuk_path/docker/KBracken
docker build --network=host -t kbracken:1.0 .
```

## Direct Snakemake Run

```bash
fastq_dir="input_fastqs"
output_dir="output"
DB="/media/uhlemann/core4/DB/kraken2DB/k2_pluspfp_16gb_20241228"

snakemake --snakefile /home/uhlemann/heekuk_path/Go_KBracken.smk \
  --config fastq_dir="$fastq_dir" output_dir="$output_dir" db="$DB" \
  --cores 8 --jobs 4 \
  --latency-wait 60 --rerun-incomplete
```

## Docker Wrapper Run

```bash
./Go_KBracken.sh \
  -i /path/to/input_fastqs \
  -o output \
  -d /media/uhlemann/core4/DB/kraken2DB/k2_pluspfp_16gb_20241228 \
  -c 8 \
  -j 4 \
  -K
```

## Input

The workflow auto-detects paired vs single-end FASTQs. It supports:

- `<sample>_R1_001.fastq.gz` / `<sample>_R2_001.fastq.gz`
- `<sample>_R1.fastq.gz` / `<sample>_R2.fastq.gz`
- `<sample>.R1.fastq.gz` / `<sample>.R2.fastq.gz`
- non-gzipped `.fastq` / `.fq`

Sample names are inferred from the full string before an explicit read token:

- `ASB_1_R1_fastq.gz` / `ASB_1_R2_fastq.gz` becomes sample `ASB_1`.
- `ASB_2_R1_fastq.gz` / `ASB_2_R2_fastq.gz` becomes sample `ASB_2`.
- `S1_1.fastq.gz` and `S1_2.fastq.gz` remain separate single-end samples because `_1/_2` alone is not treated as paired-end evidence.
- Paired-end inference requires explicit `R1/R2` or `forward/reverse` tokens.

## Output Layout

```text
output/
  1_out/
    <sample>_out.txt
    <sample>_out.log
  2_report/
    <sample>_report.txt
  3_classified/
    <sample>_classified*.fastq.gz
  4_unclassified/
    <sample>_unclassified*.fastq.gz
  5_mpa_report/
    <sample>_mpa.txt
  6_bracken_out/
    <sample>_bracken.txt
  7_bracken_mpa/
    <sample>_bracken_mpa.txt
  kraken2_mpa.txt
  bracken_mpa.txt
  bracken_mpa_filled.txt
  kraken2_log.txt
```

## Useful Configs

```bash
--config \
  fastq_dir=input_fastqs \
  output_dir=output \
  db=/path/to/kraken2_db \
  kraken_threads=4 \
  bracken_read_len=100 \
  bracken_level=S \
  bracken_threshold=10
```

## Notes

- This version intentionally keeps the `Gobracken2_V4.pl` output naming.
- This version does not depend on MetaPhlAn's `merge_metaphlan_tables.py`; it uses the bundled `scripts/merge_mpa_tables.py` for stable merge behavior.
- Per-sample `.kraken2.done` markers are added inside `1_out/` so Snakemake can track side-effect FASTQ files safely.
- `db="$DB"` is required in the direct Snakemake command. `db="DB"` would pass the literal string `DB`.
