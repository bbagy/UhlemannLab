#!/usr/bin/env python3

import argparse
import gzip
import re
import shutil
import subprocess
from pathlib import Path
from collections import Counter

import pandas as pd


###############################################################################
# Helpers
###############################################################################

def run_cmd(cmd, log_path=None, cwd=None):
    cmd_str = " ".join([str(x) for x in cmd])
    print(f"[CMD] {cmd_str}")

    if log_path is not None:
        log_path = Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as log:
            p = subprocess.run(
                cmd,
                cwd=cwd,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True
            )
    else:
        p = subprocess.run(cmd, cwd=cwd)

    if p.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd_str}")


def list_fastqs(fastq_dir: Path):
    fastqs = sorted(list(fastq_dir.glob("*.fastq.gz")))
    if len(fastqs) == 0:
        raise FileNotFoundError(f"No *.fastq.gz found in: {fastq_dir}")
    return fastqs


def sample_name_from_fastq(fq: Path):
    # e.g., KpGD42_barcode01.fastq.gz -> KpGD42_barcode01
    name = fq.name
    name = re.sub(r"\.fastq\.gz$", "", name)
    return name


def write_fasta(path: Path, header: str, seq: str):
    with open(path, "w") as f:
        f.write(f">{header}\n")
        f.write(seq.strip() + "\n")


def revcomp(seq: str):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def trim_at_polyC(seq: str, min_c: int = 12):
    """
    Find first poly-C run (>= min_c) and trim sequence BEFORE it.
    If not found, return original.
    """
    m = re.search(rf"C{{{min_c},}}", seq)
    if m:
        return seq[:m.start()]
    return seq


def count_fastq_reads(fastq_gz: Path) -> int:
    """
    Count reads in gzipped FASTQ.
    FASTQ has 4 lines per read.
    """
    n_lines = 0
    with gzip.open(fastq_gz, "rt") as f:
        for _ in f:
            n_lines += 1
    return n_lines // 4


def read_kv_stats(path: Path):
    d = {}
    with open(path) as f:
        for line in f:
            k, v = line.strip().split("\t")
            d[k] = int(v)
    return d


def fasta_length_stats(fa: Path):
    """
    Return basic length stats for a FASTA file.
    Works for multi-line FASTA as well.
    """
    if not fa.exists() or fa.stat().st_size == 0:
        return {
            "n_seqs": 0,
            "min_len": None,
            "mean_len": None,
            "median_len": None,
            "max_len": None
        }

    lengths = []
    seq = []

    with open(fa) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    lengths.append(len("".join(seq)))
                    seq = []
            else:
                seq.append(line)

        if seq:
            lengths.append(len("".join(seq)))

    if len(lengths) == 0:
        return {
            "n_seqs": 0,
            "min_len": None,
            "mean_len": None,
            "median_len": None,
            "max_len": None
        }

    lengths_sorted = sorted(lengths)
    n = len(lengths_sorted)
    mean_len = sum(lengths_sorted) / n

    if n % 2 == 1:
        median_len = lengths_sorted[n // 2]
    else:
        median_len = (lengths_sorted[n // 2 - 1] + lengths_sorted[n // 2]) / 2

    return {
        "n_seqs": n,
        "min_len": lengths_sorted[0],
        "mean_len": round(mean_len, 2),
        "median_len": median_len,
        "max_len": lengths_sorted[-1]
    }


###############################################################################
# Step 1: QC / trimming
###############################################################################

def step_qc_and_trim(
    in_fastq: Path,
    out_fastq: Path,
    qmin: int,
    min_len: int,
    overhangs_5: list,
    overhangs_3: list,
    discard_untrimmed: bool,
    qc_stat_path: Path,
    log_dir: Path,
    force: bool
):
    if out_fastq.exists() and not force:
        return

    tmp1 = out_fastq.with_suffix(".tmp1.fastq.gz")
    tmp2 = out_fastq.with_suffix(".tmp2.fastq.gz")

    # 0) raw read count
    raw_reads = count_fastq_reads(in_fastq)

    # 1) NanoFilt (optional)
    if qmin is None and min_len is None:
        shutil.copyfile(in_fastq, tmp1)
    else:
        cmd = [
            "bash", "-c",
            f"gunzip -c {in_fastq} | NanoFilt "
            f"{'' if qmin is None else f'-q {qmin}'} "
            f"{'' if min_len is None else f'-l {min_len}'} "
            f"| gzip -c > {tmp1}"
        ]
        run_cmd(cmd, log_path=log_dir / f"{in_fastq.stem}.nanofilt.log")

    nanofilt_reads = count_fastq_reads(tmp1)

    # 2) cutadapt trimming (optional)
    if (len(overhangs_5) == 0) and (len(overhangs_3) == 0):
        tmp1.rename(out_fastq)
        clean_reads = count_fastq_reads(out_fastq)

        qc_stat_path.parent.mkdir(parents=True, exist_ok=True)
        with open(qc_stat_path, "w") as f:
            f.write(f"raw_reads\t{raw_reads}\n")
            f.write(f"nanofilt_reads\t{nanofilt_reads}\n")
            f.write(f"clean_reads\t{clean_reads}\n")
        return

    cutadapt_cmd = ["cutadapt", "-j", "1"]

    # ONT에서는 기본적으로 discard를 끄는 게 맞음
    if discard_untrimmed:
        cutadapt_cmd += ["--discard-untrimmed"]

    for s in overhangs_5:
        cutadapt_cmd += ["-g", s]
    for s in overhangs_3:
        cutadapt_cmd += ["-a", s]

    cutadapt_cmd += ["-o", str(tmp2), str(tmp1)]
    run_cmd(cutadapt_cmd, log_path=log_dir / f"{in_fastq.stem}.cutadapt.log")

    tmp1.unlink(missing_ok=True)
    tmp2.rename(out_fastq)

    clean_reads = count_fastq_reads(out_fastq)

    qc_stat_path.parent.mkdir(parents=True, exist_ok=True)
    with open(qc_stat_path, "w") as f:
        f.write(f"raw_reads\t{raw_reads}\n")
        f.write(f"nanofilt_reads\t{nanofilt_reads}\n")
        f.write(f"clean_reads\t{clean_reads}\n")


###############################################################################
# Step 2: minimap2 mapping to anchor reference
###############################################################################

def step_minimap2_anchor_map(
    clean_fastq: Path,
    anchor_fa: Path,
    out_paf: Path,
    threads: int,
    log_dir: Path,
    force: bool
):
    if out_paf.exists() and not force:
        return

    out_paf.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "bash", "-c",
        f"minimap2 -t {threads} -x map-ont {anchor_fa} {clean_fastq} > {out_paf}"
    ]
    run_cmd(cmd, log_path=log_dir / f"{clean_fastq.stem}.minimap2.log")


###############################################################################
# Step 3: Extract flank (and anchorHit+flank) using PAF coordinates
###############################################################################

def parse_fastq_stream(fastq_gz: Path):
    with gzip.open(fastq_gz, "rt") as f:
        while True:
            h = f.readline().strip()
            if not h:
                break
            seq = f.readline().strip()
            _plus = f.readline().strip()
            qual = f.readline().strip()
            name = h.split()[0].replace("@", "")
            yield name, seq, qual


def load_anchor_hits_from_paf(paf: Path):
    hits = {}
    with open(paf) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            qname = parts[0]
            qstart = int(parts[2])
            qend = int(parts[3])
            strand = parts[4]
            aln_len = qend - qstart
            if qname not in hits or aln_len > hits[qname]["aln_len"]:
                hits[qname] = {
                    "strand": strand,
                    "qstart": qstart,
                    "qend": qend,
                    "aln_len": aln_len
                }
    return hits


def step_extract_flank(
    clean_fastq: Path,
    paf: Path,
    out_flank_fa: Path,
    out_anchorhit_plus_flank_fa: Path,
    flank_len: int,
    log_dir: Path,
    force: bool,
    min_polyC: int = 12
):
    """
    Produces two FASTA files:

    1) out_flank_fa
       - sequence AFTER the anchor hit end (qend)
       - anchor NOT included
       - trimmed at polyC (>= min_polyC)

    2) out_anchorhit_plus_flank_fa
       - actual read subsequence from qstart to (qend + flank_len)
       - includes anchor-hit region as it appears in the read
       - NOT trimmed (keeps polyC + primer tail if present)
       - useful for checking ONT error/SNPs inside anchor region
    """
    if out_flank_fa.exists() and out_anchorhit_plus_flank_fa.exists() and not force:
        return

    hits = load_anchor_hits_from_paf(paf)

    n_total = 0
    n_hit = 0
    n_flank = 0
    n_anchor_plus = 0
    n_flank_polyC_trimmed = 0

    out_flank_fa.parent.mkdir(parents=True, exist_ok=True)

    with open(out_flank_fa, "w") as out_flank, open(out_anchorhit_plus_flank_fa, "w") as out_anchor_plus:
        for name, seq, _qual in parse_fastq_stream(clean_fastq):
            n_total += 1
            if name not in hits:
                continue

            n_hit += 1
            h = hits[name]
            L = len(seq)

            qstart = h["qstart"]
            qend = h["qend"]
            strand = h["strand"]

            # orient to forward
            if strand == "-":
                seq = revcomp(seq)
                qstart, qend = (L - qend), (L - qstart)

            anchor_start = qstart
            anchor_end = qend

            # 1) flank-only (after anchor end)
            flank_raw = seq[anchor_end: anchor_end + flank_len].strip()
            flank_trimmed = trim_at_polyC(flank_raw, min_c=min_polyC)

            if len(flank_trimmed) >= 20:
                n_flank += 1
                if flank_trimmed != flank_raw:
                    n_flank_polyC_trimmed += 1
                out_flank.write(f">{name}\n{flank_trimmed}\n")

            # 2) anchorHit + flank (full junction region)
            anchor_plus = seq[anchor_start: anchor_end + flank_len].strip()
            if len(anchor_plus) >= 40:
                n_anchor_plus += 1
                out_anchor_plus.write(f">{name}\n{anchor_plus}\n")

    stat_path = log_dir / f"{clean_fastq.stem}.flank_stats.txt"
    stat_path.parent.mkdir(parents=True, exist_ok=True)
    with open(stat_path, "w") as f:
        f.write(f"total_reads\t{n_total}\n")
        f.write(f"anchor_hits\t{n_hit}\n")
        f.write(f"flank_written\t{n_flank}\n")
        f.write(f"anchorHit_plus_flank_written\t{n_anchor_plus}\n")
        f.write(f"flank_polyC_trimmed\t{n_flank_polyC_trimmed}\n")


###############################################################################
# Step 4: vsearch clustering (on flank-only)
###############################################################################

def step_vsearch_cluster(
    flank_fa: Path,
    out_centroids: Path,
    out_uc: Path,
    out_sizes: Path,
    identity: float,
    log_dir: Path,
    force: bool
):
    if out_centroids.exists() and not force:
        return

    out_centroids.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "vsearch",
        "--cluster_fast", str(flank_fa),
        "--id", str(identity),
        "--centroids", str(out_centroids),
        "--uc", str(out_uc),
        "--sizeout"
    ]
    run_cmd(cmd, log_path=log_dir / f"{flank_fa.stem}.vsearch.log")

    cluster_counts = Counter()
    with open(out_uc) as f:
        for line in f:
            if not (line.startswith("H") or line.startswith("S")):
                continue
            parts = line.strip().split("\t")
            cluster_id = parts[1]
            cluster_counts[cluster_id] += 1

    df = pd.DataFrame({
        "cluster_id": list(cluster_counts.keys()),
        "n_reads": list(cluster_counts.values())
    }).sort_values("n_reads", ascending=False)

    df.to_csv(out_sizes, sep="\t", index=False)


###############################################################################
# Summary
###############################################################################

def read_stats_from_flank_stats(stat_file: Path):
    d = {}
    with open(stat_file) as f:
        for line in f:
            k, v = line.strip().split("\t")
            d[k] = int(v)
    return d


def build_summary(out_dir: Path, readstats_rows, topcluster_rows):
    df1 = pd.DataFrame(readstats_rows)
    df2 = pd.DataFrame(topcluster_rows)

    xlsx = out_dir / "summary.xlsx"
    try:
        import openpyxl  # noqa: F401
        with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
            df1.to_excel(writer, sheet_name="ReadStats", index=False)
            df2.to_excel(writer, sheet_name="TopClusters", index=False)
        print(f"Summary: {xlsx}")
        return
    except ModuleNotFoundError:
        tsv1 = out_dir / "summary.ReadStats.tsv"
        tsv2 = out_dir / "summary.TopClusters.tsv"
        df1.to_csv(tsv1, sep="\t", index=False)
        df2.to_csv(tsv2, sep="\t", index=False)
        print("[WARN] openpyxl not installed. Wrote TSV summaries instead:")
        print(f"  - {tsv1}")
        print(f"  - {tsv2}")


###############################################################################
# Main
###############################################################################

def main():
    ap = argparse.ArgumentParser(
        description="Search TNseq/amplicon flank sequence using ONT reads and known anchor reference."
    )
    ap.add_argument("-i", "--input", required=True, help="Input fastq directory")
    ap.add_argument("-o", "--output", required=True, help="Output directory")

    ap.add_argument("--anchor", required=True, help="Anchor reference sequence for minimap2 (recommend >=50bp)")
    ap.add_argument("--anchor_name", default="anchor", help="Anchor name in FASTA")

    ap.add_argument("--flank_len", type=int, default=300, help="Length downstream of anchor end to extract")
    ap.add_argument("--threads", type=int, default=4)

    ap.add_argument("--qmin", type=int, default=None, help="NanoFilt Q-score min (e.g. 9)")
    ap.add_argument("--min_len", type=int, default=None, help="NanoFilt min length (e.g. 200)")

    ap.add_argument("--overhang5", action="append", default=[], help="5' overhang/adapter to trim (repeatable)")
    ap.add_argument("--overhang3", action="append", default=[], help="3' overhang/adapter to trim (repeatable)")

    ap.add_argument("--discard_untrimmed", action="store_true",
                    help="If set, cutadapt will discard reads that do not contain the adapter/overhang. (NOT recommended for ONT)")

    ap.add_argument("--cluster_id", type=float, default=0.95, help="vsearch clustering identity (0-1)")
    ap.add_argument("--top_n", type=int, default=3, help="How many top clusters to report")
    ap.add_argument("--force", action="store_true", help="Overwrite outputs")

    args = ap.parse_args()

    in_dir = Path(args.input)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    anchor_fa = out_dir / "anchor.fa"
    write_fasta(anchor_fa, args.anchor_name, args.anchor)

    qc_dir = out_dir / "1_QC"
    map_dir = out_dir / "2_minimap2"
    flank_dir = out_dir / "3_flank"
    cluster_dir = out_dir / "4_cluster"
    log_dir = out_dir / "logs"

    qc_dir.mkdir(exist_ok=True)
    map_dir.mkdir(exist_ok=True)
    flank_dir.mkdir(exist_ok=True)
    cluster_dir.mkdir(exist_ok=True)
    log_dir.mkdir(exist_ok=True)

    readstats_rows = []
    topcluster_rows = []

    fastqs = list_fastqs(in_dir)

    for fq in fastqs:
        sample = sample_name_from_fastq(fq)
        print("\n==============================")
        print(f"[SAMPLE] {sample}")
        print("==============================")

        qc_stat_file = log_dir / f"{sample}.qc_stats.txt"

        # 1) QC / trimming
        clean_fq = qc_dir / f"{sample}.clean.fastq.gz"
        step_qc_and_trim(
            in_fastq=fq,
            out_fastq=clean_fq,
            qmin=args.qmin,
            min_len=args.min_len,
            overhangs_5=args.overhang5,
            overhangs_3=args.overhang3,
            discard_untrimmed=args.discard_untrimmed,
            qc_stat_path=qc_stat_file,
            log_dir=log_dir,
            force=args.force
        )

        qc_stats = read_kv_stats(qc_stat_file)

        # 2) minimap2 mapping
        sample_map_dir = map_dir / sample
        sample_map_dir.mkdir(parents=True, exist_ok=True)
        paf = sample_map_dir / f"{sample}.anchor.paf"

        step_minimap2_anchor_map(
            clean_fastq=clean_fq,
            anchor_fa=anchor_fa,
            out_paf=paf,
            threads=args.threads,
            log_dir=log_dir,
            force=args.force
        )

        # 3) flank extraction
        flank_fa = flank_dir / f"{sample}.flank.fa"
        anchor_plus_fa = flank_dir / f"{sample}.anchorHit_plus_flank.fa"

        step_extract_flank(
            clean_fastq=clean_fq,
            paf=paf,
            out_flank_fa=flank_fa,
            out_anchorhit_plus_flank_fa=anchor_plus_fa,
            flank_len=args.flank_len,
            log_dir=log_dir,
            force=args.force
        )

        # flank stats
        stat_file = log_dir / f"{clean_fq.stem}.flank_stats.txt"
        stats = read_stats_from_flank_stats(stat_file)

        raw_reads = qc_stats.get("raw_reads", None)
        clean_reads = qc_stats.get("clean_reads", None)
        nanofilt_reads = qc_stats.get("nanofilt_reads", None)

        qc_keep_rate = None
        if raw_reads is not None and clean_reads is not None and raw_reads > 0:
            qc_keep_rate = round(clean_reads / raw_reads * 100, 2)

        # NEW: fasta length stats
        flank_len_stats = fasta_length_stats(flank_fa)
        anchor_plus_len_stats = fasta_length_stats(anchor_plus_fa)

        readstats_rows.append({
            "sample": sample,
            "raw_fastq": str(fq),
            "clean_fastq": str(clean_fq),

            "raw_reads": raw_reads,
            "nanofilt_reads": nanofilt_reads,
            "clean_reads": clean_reads,
            "qc_keep_rate_%": qc_keep_rate,

            "total_reads_used_for_mapping": stats.get("total_reads", None),
            "anchor_hits": stats.get("anchor_hits", None),
            "flank_written": stats.get("flank_written", None),
            "anchorHit_plus_flank_written": stats.get("anchorHit_plus_flank_written", None),
            "anchor_hit_rate_%": round(stats.get("anchor_hits", 0) / max(stats.get("total_reads", 1), 1) * 100, 2),

            # NEW: flank.fa length stats
            "flank_n_seqs": flank_len_stats["n_seqs"],
            "flank_min_len": flank_len_stats["min_len"],
            "flank_mean_len": flank_len_stats["mean_len"],
            "flank_median_len": flank_len_stats["median_len"],
            "flank_max_len": flank_len_stats["max_len"],

            # NEW: anchorHit_plus_flank.fa length stats
            "anchorPlus_n_seqs": anchor_plus_len_stats["n_seqs"],
            "anchorPlus_min_len": anchor_plus_len_stats["min_len"],
            "anchorPlus_mean_len": anchor_plus_len_stats["mean_len"],
            "anchorPlus_median_len": anchor_plus_len_stats["median_len"],
            "anchorPlus_max_len": anchor_plus_len_stats["max_len"],
        })

        # 4) clustering (on flank-only)
        sample_cluster_dir = cluster_dir / sample
        sample_cluster_dir.mkdir(parents=True, exist_ok=True)

        centroids = sample_cluster_dir / f"{sample}.centroids.fa"
        uc = sample_cluster_dir / f"{sample}.clusters.uc"
        sizes = sample_cluster_dir / f"{sample}.cluster_sizes.tsv"

        step_vsearch_cluster(
            flank_fa=flank_fa,
            out_centroids=centroids,
            out_uc=uc,
            out_sizes=sizes,
            identity=args.cluster_id,
            log_dir=log_dir,
            force=args.force
        )

        if sizes.exists() and sizes.stat().st_size > 0:
            df_sizes = pd.read_csv(sizes, sep="\t")
            if df_sizes.shape[0] > 0 and df_sizes["n_reads"].sum() > 0:
                df_sizes["fraction"] = df_sizes["n_reads"] / df_sizes["n_reads"].sum()
                for _, row in df_sizes.head(args.top_n).iterrows():
                    topcluster_rows.append({
                        "sample": sample,
                        "cluster_id": row["cluster_id"],
                        "n_reads": int(row["n_reads"]),
                        "fraction_%": round(float(row["fraction"]) * 100, 2),
                        "centroids_fasta": str(centroids)
                    })

    build_summary(out_dir, readstats_rows, topcluster_rows)

    print("\nDONE")


if __name__ == "__main__":
    main()