# -*- coding: utf-8 -*-
import os
import glob
import pandas as pd
import sys

# --- Get the output directory from the command line arguments ---
out_dir = sys.argv[1]

# --- Input file: HolomiRA results ---
input_file = os.path.join(out_dir, "final_results/HolomiRA_results.tsv")

# --- Ensure the input file exists ---
if not os.path.exists(input_file):
    print(f"Error: {input_file} does not exist.")
    sys.exit(1)

# --- Read the HolomiRA results to extract MAGs and miRNA names ---
holomira_results = pd.read_csv(input_file, sep="\t")
MAG_names = set(holomira_results["MAG"].dropna().unique())  # Ensure no NaN values
miRNA_names = set(holomira_results["miRNA"].dropna().unique())

# --- Identify all environments dynamically from the folder names ---
environment_dirs = glob.glob(os.path.join(out_dir, "function/miRNA_*"))
environments = [os.path.basename(d).replace("miRNA_", "") for d in environment_dirs]

# --- Initialize lists to collect all concatenated files from all environments ---
all_miRNA_files = []
all_MAG_files = []

for environment in environments:
    # Define paths for the current environment
    by_miRNA_dir = os.path.join(out_dir, f"function/miRNA_{environment}")
    by_MAG_dir = os.path.join(out_dir, f"function/MAGs_{environment}")

    # Ensure the directories exist
    for dir_to_create in [by_miRNA_dir, by_MAG_dir]:
        if not os.path.exists(dir_to_create):
            os.makedirs(dir_to_create)

    # Merge the affected CDS sequences for each miRNA_name in the current environment
    for miRNA_name in miRNA_names:
        miRNA_files = glob.glob(os.path.join(out_dir, f"annotation/*/temp_*_affected_cds_miRNA_{miRNA_name}_{environment}.fasta"))

        if miRNA_files:
            output_file_path = os.path.join(by_miRNA_dir, f"concatenated_by_{miRNA_name}_{environment}.fasta")
            with open(output_file_path, "w") as outfile:
                for fname in miRNA_files:
                    with open(fname) as infile:
                        outfile.write(infile.read())
            all_miRNA_files.append(output_file_path)  # Collect the miRNA files
        else:
            print(f"No files found for miRNA: {miRNA_name} in environment: {environment}")

    # Check for the presence of contig files for each MAG before merging sequences
    contig_files = glob.glob(os.path.join(by_MAG_dir, "contigs_for_*.txt"))
    mag_contig_present = set(
        [os.path.basename(f).replace("contigs_for_", "").replace(".txt", "") for f in contig_files]
    )

    # Merge the affected CDS sequences for each MAG in the current environment
    for MAG_name in MAG_names:
        MAG_files = glob.glob(os.path.join(out_dir, f"annotation/*/temp_{MAG_name}_affected_cds_miRNA_*_{environment}.fasta"))

        if MAG_name in mag_contig_present:
            if MAG_files:
                output_file_path = os.path.join(by_MAG_dir, f"concatenated_by_{MAG_name}_{environment}.fasta")
                with open(output_file_path, "w") as outfile:
                    for fname in MAG_files:
                        with open(fname) as infile:
                            outfile.write(infile.read())
                all_MAG_files.append(output_file_path)  # Collect the MAG files
            else:
                print(f"No files found for MAG: {MAG_name} in environment: {environment}")
        else:
            print(f"Skipping MAG: {MAG_name} as no contigs are present in {by_MAG_dir} for environment: {environment}")

# --- Create the general control files for miRNA and MAG concatenated results (no environment separation) ---
miRNA_control_file = os.path.join(out_dir, "function/temp_concat_end_miRNA")
MAG_control_file = os.path.join(out_dir, "function/temp_concat_end_MAG")

# ---  Write the list of all concatenated miRNA fasta files to the general control file ---
with open(miRNA_control_file, "w") as control_file:
    if all_miRNA_files:
        for file_path in all_miRNA_files:
            control_file.write(file_path + "\n")
    else:
        print("No concatenated miRNA files found to write in general control file.")

# ---  Write the list of all concatenated MAG fasta files to the general control file ---
with open(MAG_control_file, "w") as control_file:
    if all_MAG_files:
        for file_path in all_MAG_files:
            control_file.write(file_path + "\n")
    else:
        print("No concatenated MAG files found to write in general control file.")