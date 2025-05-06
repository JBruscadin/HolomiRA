#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd
import re
from pathlib import Path

# --- Inputs ---
rna_dir = sys.argv[1]
finalresults = sys.argv[2]
output_dir = sys.argv[3]
dg_cutoff = float(sys.argv[4])

# --- Create output directory if needed ---
os.makedirs(output_dir, exist_ok=True)

print("Reading finalresults.txt...")
final_df = pd.read_csv(finalresults, sep="\t")

print("Scanning RNAup output files...")
rna_rows = []
for root, _, files in os.walk(rna_dir):
    for f in files:
        if not f.endswith("_rnaup.txt"):
            continue
        file_path = os.path.join(root, f)
        try:
            with open(file_path, encoding="utf-8") as file:
                lines = [line.strip() for line in file if line.strip()]
        except UnicodeDecodeError:
            with open(file_path, encoding="latin1", errors="ignore") as file:
                lines = [line.strip() for line in file if line.strip()]

        if len(lines) < 2:
            continue

        header = lines[0].lstrip(">")
        align_line = lines[1]

        energy_match = re.search(r"\(([-\d\.]+) = ([-\d\.]+) \+ ([-\d\.]+) \+ ([-\d\.]+)\)", align_line)
        coord_match = re.search(r"(\d+),(\d+)\s+:\s+(\d+,\d+)", align_line)

        if not (energy_match and coord_match):
            continue

        dg_total = float(energy_match.group(1))
        dg_binding = float(energy_match.group(2))
        dg_open_target = float(energy_match.group(3))
        dg_open_mirna = float(energy_match.group(4))
        pos1 = int(coord_match.group(1))
        pos2 = int(coord_match.group(2))
        mirna_pairing = coord_match.group(3)

        # Parse filename to extract miRNA, contig, and coordinates
        match = re.match(r"(.+?)_(gnl_X_.+?)_(\d+)_(\d+)_rnaup\.txt", f)
        if not match:
            continue

        mirna = match.group(1)
        contig = match.group(2).replace("gnl_X_", "gnl|X|")
        start = int(match.group(3))
        end = int(match.group(4))

        rna_rows.append({
            "Contig": contig,
            "Start": start,
            "End": end,
            "miRNA": mirna,
            "pos1": pos1,
            "pos2": pos2,
            "mirNA_pairing": mirna_pairing,
            "dG_total": dg_total,
            "dG_binding": dg_binding,
            "dG_opening_target": dg_open_target,
            "dG_opening_miRNA": dg_open_mirna
        })

rna_df = pd.DataFrame(rna_rows)

if rna_df.empty:
    print("[!] RNAup result DataFrame is empty. No valid entries parsed.")
    summary_path = os.path.join(rna_dir, "RNAup_summary_results.tsv")
    rna_df.to_csv(summary_path, sep="\t", index=False)
    final_df["dG_total"] = pd.NA
    final_df["dG_binding"] = pd.NA
    final_df["dG_opening_target"] = pd.NA
    final_df["dG_opening_miRNA"] = pd.NA
    final_df.to_csv(os.path.join(output_dir, "HolomiRA_discarded.tsv"), sep="\t", index=False)
    open(os.path.join(output_dir, "HolomiRA_results.tsv"), "w").close()
    sys.exit(0)

summary_path = os.path.join(rna_dir, "RNAup_summary_results.tsv")
rna_df.to_csv(summary_path, sep="\t", index=False)
print(f"Summary saved: {summary_path}")

# --- Filter: Keep rows where the binding is inside the target region ---
rna_df_filtered = rna_df[(rna_df["pos1"] > 150) & (rna_df["pos2"] < 186)]

# --- Ensure consistent types for merging ---
for col in ["Start", "End"]:
    final_df[col] = pd.to_numeric(final_df[col], errors="coerce").astype("Int64")
    rna_df_filtered[col] = pd.to_numeric(rna_df_filtered[col], errors="coerce").astype("Int64")

print("Merging finalresults with RNAup information...")
merged = pd.merge(final_df, rna_df_filtered, how="left", on=["miRNA", "Contig", "Start", "End"])

# --- Filter based on dG_total threshold ---
valid = merged[merged["dG_total"] <= dg_cutoff]
discarded = merged[~(merged["dG_total"] <= dg_cutoff)]

# --- Save final tables ---
valid.to_csv(os.path.join(output_dir, "HolomiRA_results.tsv"), sep="\t", index=False)
discarded.to_csv(os.path.join(output_dir, "HolomiRA_discarded.tsv"), sep="\t", index=False)

print(f"\nMerge complete!\n   Valid hits:     {len(valid)}\n   Discarded hits: {len(discarded)}")
