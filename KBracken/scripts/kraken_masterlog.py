#!/usr/bin/env python3
import argparse
import os
import re


CLASSIFIED_RE = re.compile(r"(\d+) sequences classified \((\S+)%\)")
UNCLASSIFIED_RE = re.compile(r"(\d+) sequences unclassified \((\S+)%\)")


def parse_log(path):
    with open(path) as fh:
        lines = [line.strip() for line in fh if line.strip()]
    classified = classified_percent = unclassified = unclassified_percent = None
    for line in lines:
        m = CLASSIFIED_RE.search(line)
        if m:
            classified, classified_percent = int(m.group(1)), m.group(2)
        m = UNCLASSIFIED_RE.search(line)
        if m:
            unclassified, unclassified_percent = int(m.group(1)), m.group(2)
    if classified is None or unclassified is None:
        raise ValueError(f"Could not parse Kraken2 classified/unclassified counts from {path}")
    return classified, classified_percent, unclassified, unclassified_percent


def main():
    parser = argparse.ArgumentParser(description="Create Gobracken2-style master Kraken2 log.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--keep-logs", type=int, default=1)
    parser.add_argument("logs", nargs="+")
    args = parser.parse_args()

    with open(args.output, "w") as out_fh:
        out_fh.write("sample total_reads     classified      unclassified    classified_percent      unclassified_percent\n")
        for log_path in args.logs:
            classified, classified_percent, unclassified, unclassified_percent = parse_log(log_path)
            total = classified + unclassified
            sample = os.path.basename(log_path).removesuffix("_out.log")
            out_fh.write(f"{sample}  {total}  {classified}     {unclassified}   {classified_percent}     {unclassified_percent}\n")
            if not args.keep_logs:
                os.remove(log_path)


if __name__ == "__main__":
    main()
