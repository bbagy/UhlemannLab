#!/usr/bin/env python3
"""Create fixed-read-length paired FASTQs for FIGARO quality modeling."""

import collections
import gzip
import itertools
import os


def records(path):
    with gzip.open(path, "rt") as handle:
        while True:
            record = tuple(itertools.islice(handle, 4))
            if not record:
                return
            if len(record) != 4:
                raise ValueError(f"Incomplete FASTQ record in {path}")
            yield record


def modal_length(paths, sample_reads=10000):
    counts = collections.Counter()
    for path in paths:
        for record in itertools.islice(records(path), sample_reads):
            counts[len(record[1].rstrip("\r\n"))] += 1
    if not counts:
        raise ValueError("No FASTQ reads found while determining modal read length")
    # Illumina adapter trimming truncates primer-dimer/no-insert reads to a fixed
    # short length (commonly 35bp) regardless of sample -- in low-biomass runs
    # this spike can outnumber real amplicon reads, which vary in length and so
    # split their vote across several lengths near the true (near-full-cycle)
    # length. Restrict the mode to lengths near the observed max to avoid
    # picking the adapter-trim artifact.
    max_len = max(counts)
    candidates = {length: n for length, n in counts.items() if length >= 0.9 * max_len}
    return max(candidates, key=candidates.get)


inputs = list(snakemake.input.fastqs)
if len(inputs) % 2:
    raise ValueError("FIGARO staging requires complete R1/R2 FASTQ pairs")

r1s = inputs[0::2]
r2s = inputs[1::2]
target_r1 = modal_length(r1s)
target_r2 = modal_length(r2s)
os.makedirs(snakemake.output.fastq_dir, exist_ok=True)

stats = []
for r1, r2 in zip(r1s, r2s):
    out1 = os.path.join(snakemake.output.fastq_dir, os.path.basename(r1))
    out2 = os.path.join(snakemake.output.fastq_dir, os.path.basename(r2))
    kept = dropped = 0
    with gzip.open(out1, "wt", compresslevel=1) as dst1, gzip.open(
        out2, "wt", compresslevel=1
    ) as dst2:
        for rec1, rec2 in itertools.zip_longest(records(r1), records(r2)):
            if rec1 is None or rec2 is None:
                raise ValueError(f"Different record counts in paired FASTQs: {r1}, {r2}")
            len1 = len(rec1[1].rstrip("\r\n"))
            len2 = len(rec2[1].rstrip("\r\n"))
            if len1 == target_r1 and len2 == target_r2:
                dst1.writelines(rec1)
                dst2.writelines(rec2)
                kept += 1
            else:
                dropped += 1
    if kept == 0:
        raise ValueError(f"No fixed-length read pairs retained for {r1}, {r2}")
    stats.append((os.path.basename(r1), kept, dropped))

with open(os.path.join(snakemake.output.fastq_dir, "preparation.tsv"), "w") as handle:
    handle.write(f"# modal_R1={target_r1}\tmodal_R2={target_r2}\n")
    handle.write("sample_R1\tkept_pairs\tdropped_pairs\n")
    for name, kept, dropped in stats:
        handle.write(f"{name}\t{kept}\t{dropped}\n")
