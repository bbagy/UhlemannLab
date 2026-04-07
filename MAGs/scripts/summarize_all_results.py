#!/usr/bin/env python3

# summarize_all_results.py
import pandas as pd
import os
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--gtdb_dir', required=True)
parser.add_argument('--prokka_dir', required=True)
parser.add_argument('--eggnog_dir', required=True)
parser.add_argument('--kofam_dir', required=True)
parser.add_argument('--output_excel', required=True)
args = parser.parse_args()

# 1️⃣ GTDB-Tk summary
summary_files = glob.glob(os.path.join(args.gtdb_dir, "*/gtdbtk.bac120.summary.tsv"))
gtdb_list = [pd.read_csv(f, sep="\t") for f in summary_files if os.path.getsize(f) > 0]
gtdb_df = pd.concat(gtdb_list, ignore_index=True) if gtdb_list else pd.DataFrame()

# 2️⃣ Prokka stats (detailed)
prokka_summaries = []
for d in glob.glob(os.path.join(args.prokka_dir, "*/")):
    sample = os.path.basename(os.path.normpath(d))
    stats_file = os.path.join(d, f"{sample}.txt")
    if os.path.exists(stats_file):
        entry = {"sample": sample}
        with open(stats_file) as f:
            for line in f:
                if ":" in line:
                    key, value = line.strip().split(":", 1)
                    entry[key.strip()] = value.strip()
        prokka_summaries.append(entry)
prokka_df = pd.DataFrame(prokka_summaries)

# 3️⃣ EggNOG hits
eggnog_files = glob.glob(os.path.join(args.eggnog_dir, "*_eggnog.emapper.annotations"))
eggnog_df = pd.DataFrame()
if eggnog_files:
    dfs = []
    for f in eggnog_files:
        try:
            df = pd.read_csv(f, sep='\t', comment='#')
            sample = os.path.basename(f).split("_eggnog.")[0]
            df.insert(0, "sample", sample)
            dfs.append(df)
        except Exception as e:
            print(f"⚠️ Error reading {f}: {e}")
    eggnog_df = pd.concat(dfs, ignore_index=True)

# 4️⃣ KofamScan results
kofam_files = glob.glob(os.path.join(args.kofam_dir, "*_kofamscan.txt"))
kofam_dfs = []
for f in kofam_files:
    try:
        df = pd.read_csv(f, sep='\t', comment='#', header=None)
        df.columns = ["query", "KO", "score", "threshold", "Evalue", "KO_name"]
        sample = os.path.basename(f).split("_kofamscan.")[0]
        df.insert(0, "sample", sample)
        kofam_dfs.append(df)
    except Exception as e:
        print(f"⚠️ Error reading {f}: {e}")
kofam_df = pd.concat(kofam_dfs, ignore_index=True) if kofam_dfs else pd.DataFrame()

# Write to Excel
with pd.ExcelWriter(args.output_excel) as writer:
    gtdb_df.to_excel(writer, sheet_name="GTDB_Tk", index=False)
    prokka_df.to_excel(writer, sheet_name="Prokka_Stats", index=False)
    eggnog_df.to_excel(writer, sheet_name="EggNOG", index=False)
    kofam_df.to_excel(writer, sheet_name="KofamScan", index=False)

print(f"✅ Summary Excel saved to: {args.output_excel}")
