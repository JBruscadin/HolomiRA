#!/usr/bin/env python3
import sys
import os
import pandas as pd

# --- Inputs ---
finalresults = sys.argv[1]        # finalresults.txt
fasta_dir = sys.argv[2]           # path to .fna files
output_prefix = sys.argv[3]       # prefix for output (e.g., OUT_DIR/structure/sig_hits)
window = 150                      # nt upstream and downstream

# --- Load table ---
df = pd.read_csv(finalresults, sep="\t", comment="#").dropna(how="all")

# --- Validate required columns ---
required_cols = ["Contig", "Start", "End", "miRNA", "MAG"]
if not all(col in df.columns for col in required_cols):
    print("? ERROR: Missing required columns in finalresults.txt")
    sys.exit(1)

# --- Sanitize positions ---
df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
df["End"] = pd.to_numeric(df["End"], errors="coerce")
df = df.dropna(subset=["Start", "End"])

# --- Remove duplicated GFF target IDs ---
df["Target_ID"] = df.apply(lambda row: f"{row['miRNA']}_{row['MAG']}", axis=1)
df = df.drop_duplicates(subset=["Contig", "Target_ID", "Start", "End"])

# --- Create GFF entries (no strand info) ---
gff_entries = []
for _, row in df.iterrows():
    contig = row["Contig"]
    start = int(row["Start"])
    end = int(row["End"])
    name = row["Target_ID"]

    center = (start + end) // 2
    gff_start = max(1, center - window)
    gff_end = center + window

    gff_line = [
        contig,
        "RNAhybrid",
        "target_site",
        str(gff_start),
        str(gff_end),
        ".",
        ".",
        ".",
        f"ID={name}"
    ]
    gff_entries.append("\t".join(gff_line))

# --- Save GFF ---
gff_out = f"{output_prefix}.gff"
os.makedirs(os.path.dirname(gff_out), exist_ok=True)
with open(gff_out, "w") as f:
    for line in gff_entries:
        f.write(line + "\n")
#print(f"? GFF written: {gff_out}")

# --- Identify .fna files by contig ---
unique_contigs = set(df["Contig"])
contig_to_fasta = {}

for root, _, files in os.walk(fasta_dir):
    for file in files:
        if file.endswith(".fna"):
            full_path = os.path.join(root, file)
            with open(full_path) as f:
                for line in f:
                    if line.startswith(">"):
                        header = line[1:].strip().split()[0]
                        if header not in contig_to_fasta:
                            contig_to_fasta[header] = full_path

# --- Merge FASTA files that match contigs ---
matched_fastas = set()
for contig in unique_contigs:
    if contig in contig_to_fasta:
        matched_fastas.add(contig_to_fasta[contig])
    else:
        print(f" Contig not found in FASTA files: {contig}")

if matched_fastas:
    merged_fasta = f"{output_prefix}_merged.fna"
    seen = set()
    with open(merged_fasta, "w") as out:
        for fasta in matched_fastas:
            with open(fasta) as f:
                for record in f.read().split(">")[1:]:
                    header = record.splitlines()[0].strip().split()[0]
                    if header not in seen:
                        seen.add(header)
                        out.write(">" + record)
    print(f"Merged FASTA: {merged_fasta}")

# --- Run bedtools (without -s to ignore strand) ---
    fasta_out = f"{output_prefix}.fasta"
    cmd = f"bedtools getfasta -fi {merged_fasta} -bed {gff_out} -fo {fasta_out}"
    os.system(cmd)
    print(f"FASTA written: {fasta_out}")
else:
    print("No matching contigs found in .fna files.")
