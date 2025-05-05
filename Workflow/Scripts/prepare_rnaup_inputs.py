#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import pandas as pd
from Bio import SeqIO
import re

# --- Input arguments ---
result_file = sys.argv[1]
mirna_file = sys.argv[2]
fasta_file = sys.argv[3]
output_folder = sys.argv[4]

# --- Load data ---
df = pd.read_csv(result_file, sep="\t")
required_cols = {"Start", "End", "Contig", "miRNA"}
if not required_cols.issubset(df.columns):
    raise ValueError(f"Missing columns: {required_cols - set(df.columns)}")

df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
df["End"] = pd.to_numeric(df["End"], errors="coerce")
df = df.dropna(subset=["Start", "End"])

# --- Load sequences ---
mirnas = SeqIO.to_dict(SeqIO.parse(mirna_file, "fasta"))
targets = list(SeqIO.parse(fasta_file, "fasta"))

# --- Prepare output ---
print("Generating RNAup inputs...")
os.makedirs(output_folder, exist_ok=True)
metadata_list = []

for _, row in df.iterrows():
    contig = row["Contig"]
    start = int(row["Start"])
    end = int(row["End"])
    mirna = row["miRNA"]

    try:
        mi_seq = str(mirnas[mirna].seq)
    except KeyError:
        print(f"[!] miRNA not found: {mirna}")
        continue

    found = False
    for record in targets:
        if contig in record.id:
            try:
                range_str = record.id.split(":")[1]
                win_start, win_end = map(int, range_str.split("-"))
            except Exception:
                continue

            if win_start <= start and end <= win_end:
                target_seq = str(record.seq)
                found = True
                break

    if not found:
        print(f"[!] No matching window found for: {contig}:{start}-{end}")
        continue

    safe_contig = re.sub(r"[^a-zA-Z0-9_]", "_", contig)
    seq_id = f"{mirna}_{safe_contig}_{start}_{end}"
    output_path = os.path.join(output_folder, f"{seq_id}.fa")
    with open(output_path, "w") as f:
        f.write(f">{seq_id}\n{mi_seq}&{target_seq}\n")

    metadata_list.append({
        "file": f"{seq_id}.fa",
        "miRNA": mirna,
        "Contig": contig,
        "Start": start,
        "End": end
    })

# --- Save metadata ---
metadata_df = pd.DataFrame(metadata_list)
metadata_df.to_csv(os.path.join(output_folder, "input_metadata.tsv"), sep="\t", index=False)

# --- Completion marker ---
open(os.path.join(output_folder, ".inputs_prepared"), "w").close()
print("Done.")