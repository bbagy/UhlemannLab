#!/usr/bin/env python3
import argparse
import os
from collections import OrderedDict


def sample_name(path):
    name = os.path.basename(path)
    for suffix in [
        "_bracken_mpa.txt",
        "_mpa.txt",
        ".txt",
        ".tsv",
    ]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


def read_mpa(path):
    values = OrderedDict()
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("ID\t"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            values[fields[0]] = fields[1]
    return values


def main():
    parser = argparse.ArgumentParser(description="Merge MPA-style two-column tables with stable KBracken behavior.")
    parser.add_argument("--output", required=True)
    parser.add_argument("tables", nargs="+")
    args = parser.parse_args()

    sample_to_values = OrderedDict()
    taxon_order = OrderedDict()

    for path in args.tables:
        sample = sample_name(path)
        values = read_mpa(path)
        sample_to_values[sample] = values
        for taxon in values:
            taxon_order.setdefault(taxon, None)

    with open(args.output, "w") as out_fh:
        out_fh.write("ID\t" + "\t".join(sample_to_values.keys()) + "\n")
        for taxon in taxon_order:
            row = [taxon]
            for values in sample_to_values.values():
                row.append(values.get(taxon, "0"))
            out_fh.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
