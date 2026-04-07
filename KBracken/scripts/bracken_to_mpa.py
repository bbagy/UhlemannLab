#!/usr/bin/env python3
import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description="Convert Bracken species output to the MPA-like format used by Gobracken2_V4.pl.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.input, newline="") as in_fh, open(args.output, "w", newline="") as out_fh:
        reader = csv.DictReader(in_fh, delimiter="\t")
        for row in reader:
            name = row.get("name", "").replace(" ", "_")
            new_est_reads = row.get("new_est_reads", "")
            if not name:
                continue
            out_fh.write(f"s__{name}\t{new_est_reads}\n")


if __name__ == "__main__":
    main()
