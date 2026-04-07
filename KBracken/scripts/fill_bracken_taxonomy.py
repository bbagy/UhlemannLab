#!/usr/bin/env python3
import argparse
import csv


def load_kraken_taxonomy(path):
    taxonomy = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for fields in reader:
            if not fields or fields[0] == "ID":
                continue
            taxon = fields[0]
            if not all(rank in taxon for rank in ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]):
                continue
            ranks = taxon.split("|")
            if len(ranks) < 7:
                continue
            cleaned = {
                "phylum": ranks[1].removeprefix("p__"),
                "class": ranks[2].removeprefix("c__"),
                "order": ranks[3].removeprefix("o__"),
                "family": ranks[4].removeprefix("f__"),
                "genus": ranks[5].removeprefix("g__"),
                "species": ranks[6].removeprefix("s__").replace("_", " "),
            }
            taxonomy[cleaned["species"]] = cleaned
    return taxonomy


def main():
    parser = argparse.ArgumentParser(description="Append Kraken2 taxonomy ranks to merged Bracken MPA table.")
    parser.add_argument("--kraken-mpa", required=True)
    parser.add_argument("--bracken-mpa", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    kraken_taxonomy = load_kraken_taxonomy(args.kraken_mpa)

    with open(args.bracken_mpa, newline="") as in_fh, open(args.output, "w", newline="") as out_fh:
        for raw in in_fh:
            line = raw.rstrip("\n")
            if line.startswith("#") or line.startswith("ID"):
                out_fh.write(f"{line}\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
                continue
            if not line:
                out_fh.write("\n")
                continue
            fields = line.split("\t")
            species = fields[0].removeprefix("s__").replace("_", " ")
            tax = kraken_taxonomy.get(species)
            if tax:
                out_fh.write(
                    f"{line}\t{tax['phylum']}\t{tax['class']}\t{tax['order']}\t"
                    f"{tax['family']}\t{tax['genus']}\t{tax['species']}\n"
                )
            else:
                out_fh.write(f"{line}\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n")


if __name__ == "__main__":
    main()
