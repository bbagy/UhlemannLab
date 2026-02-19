# PFsnake — *Plasmodium falciparum* WGS pipeline (Illumina)

This pipeline performs standardized Illumina WGS analysis for *P. falciparum*.

---

## What this pipeline does

- QC
- Mapping
- Variant calling
- CNV
- Drug resistance summary
- Final reports (tables / Excel)

---

## Requirements

- Docker
- Input FASTQs (paired-end)

---

## Build docker image

```bash
docker build -t pf-snake:1.0 .
```

---

## Run

Move to your project directory (where FASTQs exist), then run:

```bash
/path/to/PFsnake/run.sh \
  --fastq-dir fastqs \
  --out out_pf \
  --cores 8
```

---

## Inputs

- FASTQ naming convention:
  - `*_R1_*.fastq.gz`
  - `*_R2_*.fastq.gz`

---

## Output structure

Output will be written to:

- `out_pf/`

Main subfolders:
- `1_QC/`
- `2_mapping/`
- `3_variants/`
- `6_cnv/`
- `9_drug_resistance/`

---

## Troubleshooting

### Docker permission denied
If you see:

`permission denied while trying to connect to the docker API`

You likely need to run docker as a user in the docker group.

---

## Maintainer
Heekuk Park
