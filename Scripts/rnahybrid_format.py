#!/usr/bin/env python
import pandas as pd
import sys
import os

# Retrieve command line arguments: MAG_ID, output directory, and input file
MAG_ID = sys.argv[1]
out_dir = sys.argv[2]
input_file = sys.argv[3]

# Define file paths for the necessary files
cds_file = f"{out_dir}/annotation/{MAG_ID}/{MAG_ID}_cds_fiveprime.gff"
tsv_file = f"{out_dir}/annotation/{MAG_ID}/{MAG_ID}.tsv"
id_info_file = f"{out_dir}/annotation/{MAG_ID}/{MAG_ID}_IDs.txt"
output_file_path = f"{out_dir}/rnahybrid/{MAG_ID}_bsites.tsv"

# Header for the output file (includes start_gene and end_gene)
header = "sample\tseq\tstart\tend\tmir\tID\tmfe\tp\tgene\tcds_start\tcds_end\tstart_gene\tend_gene\n"

# If the input file is empty, create a file with only the header and exit
if os.stat(input_file).st_size == 0:
    with open(output_file_path, "w") as output_file:
        output_file.write(header)
    sys.exit()

# Read the input file containing miRNA binding information
a = pd.read_csv(input_file, sep=':', names=['seq', 'position', 'nc', 'mir', 'ncmir', 'mfe', 'p', '1', '2', '3', '4', '5'])
a[['startend', 'strand']] = a['position'].str.split("(", expand=True)
a[['start', 'end']] = a['startend'].str.split("-", expand=True)
a['start'] = a['start'].astype(int)
a['end'] = a['end'].astype(int)

# Read the *_cds_fiveprime.gff file (contains gene coordinates)
id_df = pd.read_csv(cds_file, sep="\t", header=None, names=["seq", "cds_start", "cds_end", "ID", "strand"])
id_df["ID"] = id_df["ID"].str.replace("ID=", "", regex=True).str.strip()
id_df["cds_start"] = id_df["cds_start"].astype(int)
id_df["cds_end"] = id_df["cds_end"].astype(int)

# Read the gene annotation TSV file and filter for gene entries
tsv = pd.read_csv(tsv_file, sep="\t")
tsv = tsv[tsv['ftype'] == 'gene'][["locus_tag", "gene"]].drop_duplicates()
tsv["locus_tag"] = tsv["locus_tag"].str.strip()
tsv["gene"] = tsv["gene"].fillna("")

# Create a list to store merged data from miRNA bindings and gene coordinates
merged_data = []
for _, row in a.iterrows():
    matches = id_df[
        (id_df["seq"] == row["seq"]) &
        (id_df["cds_start"] <= row["end"]) &
        (id_df["cds_end"] >= row["start"])
    ]
    for _, match in matches.iterrows():
        merged_data.append({
            "sample": MAG_ID,
            "seq": row["seq"],
            "start": row["start"],
            "end": row["end"],
            "mir": row["mir"],
            "ID": match["ID"],
            "mfe": row["mfe"],
            "p": row["p"],
            "cds_start": match["cds_start"],
            "cds_end": match["cds_end"]
        })

# If no matches are found, create an empty file with only the header and exit
if not merged_data:
    with open(output_file_path, "w") as output_file:
        output_file.write(header)
    sys.exit()

# Create DataFrame from merged data and merge with gene annotations
a_id = pd.DataFrame(merged_data)
final_data = a_id.merge(tsv, how="left", left_on="ID", right_on="locus_tag")
final_data = final_data[["sample", "seq", "start", "end", "mir", "ID", "mfe", "p", "gene", "cds_start", "cds_end"]]

# Read the {MAG_ID}_ID.txt file to extract gene start and end positions
id_info = pd.read_csv(id_info_file, sep="\t", header=None, names=["seq", "start_gene", "end_gene", "ID_raw", "strand"])
id_info["ID"] = id_info["ID_raw"].str.replace("ID=", "", regex=True).str.strip()
id_info = id_info[["ID", "start_gene", "end_gene"]]

# Merge the gene start and end info into final_data based on matching ID
final_data = final_data.merge(id_info, on="ID", how="left")

# Save only the final merged data to the output file
final_data.to_csv(output_file_path, sep="\t", index=False)
