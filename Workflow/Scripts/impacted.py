#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd
import subprocess
import warnings

# ---  Get input file and output directory from command-line arguments ---
if len(sys.argv) < 3:
    sys.exit(1)

input_file = sys.argv[1]
out_dir = sys.argv[2]

# ---  Create output directories if they do not exist ---
output_dir_function = os.path.join(out_dir, "function")
os.makedirs(output_dir_function, exist_ok=True)

# ---  Dictionaries to store contigs per miRNA and MAGs per environment ---
contig_sets_by_mirna = {}
mag_sets_by_environment = {}

print(f"Reading input file: {input_file}")

# ---  Check if the input file exists ---
if not os.path.exists(input_file):
    print(f"ERROR: Input file {input_file} does not exist.")
    sys.exit(1)

# --- Read input file ---
with open(input_file, 'r', encoding='utf-8', errors='replace') as file:
    header = file.readline().strip().split('\t')  # Read and identify header
    print(f"Detected header ({len(header)} columns): {header}")

    for line_number, line in enumerate(file, start=2):
        columns = line.strip().split('\t')

        # Ensure the line has exactly 16 columns
        if len(columns) != 22:
            print(f"ERROR: Line {line_number} has {len(columns)} columns instead of 22. Skipping: {columns}")
            continue

        try:
            mag_name     = columns[0]  # MAG name
            contig_value = columns[1]  # Contig name
            start_gene   = int(columns[11])  # Gene start coordinate
            end_gene     = int(columns[12])  # Gene end coordinate
            mirna_name   = columns[4]  # miRNA name
            environment  = columns[14]  # Environment
        except ValueError as e:
            print(f"ERROR: Line {line_number} contains invalid values! {e} -> {columns}")
            continue

        # Store data for miRNAs
        contig_sets_by_mirna.setdefault(mirna_name, {}).setdefault(environment, []).append(
            (mag_name, contig_value, start_gene, end_gene)
        )

        # Store data for MAGs
        mag_sets_by_environment.setdefault(environment, {}).setdefault(mag_name, []).append(
            (contig_value, start_gene, end_gene)
        )

# --- EXTRACTING FASTA SEQUENCES FOR EACH MAG ---

all_affected_cds_files = {}

MAG_folders = glob.glob(os.path.join(out_dir, "annotation/*"))

for MAG_folder in MAG_folders:
    MAG_name = os.path.basename(MAG_folder)

    gff_file = os.path.join(MAG_folder, f"{MAG_name}_cds.gff")
    if not os.path.exists(gff_file):
        print(f"WARNING: GFF file {gff_file} not found. Skipping {MAG_name}.")
        continue

    df = pd.read_csv(
        gff_file, sep='\t', comment='#', header=None,
        names=["contig", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    )
    df["contig"] = df["contig"].str.strip()

    for environment, mag_dict in mag_sets_by_environment.items():
        if MAG_name not in mag_dict:
            continue

        output_dir_env = os.path.join(out_dir, f"function/MAGs_{environment}")
        os.makedirs(output_dir_env, exist_ok=True)

        bed_file_path = os.path.join(output_dir_env, f"temp_{MAG_name}_{environment}.bed")
        contig_for_path = os.path.join(output_dir_env, f"contigs_for_{MAG_name}.txt")

        # Create BED file and contigs_for_* only with data from the respective MAG
        with open(bed_file_path, 'w') as bed_file, open(contig_for_path, 'w') as contig_file:
            for contig, start_gene, end_gene in mag_dict[MAG_name]:
                bed_file.write(f"{contig}\t{start_gene-1}\t{end_gene}\n")  # BED format is 0-based
                contig_file.write(f"{contig} {start_gene} {end_gene}\n")  # Include coordinates in contig_for_*

        input_fasta = os.path.join(MAG_folder, f"{MAG_name}.fna")
        output_fa = os.path.join(output_dir_env, f"temp_{MAG_name}_affected_cds.fasta")

        warnings.filterwarnings("ignore")
        cmd = f"bedtools getfasta -fi {input_fasta} -bed {bed_file_path} -fo {output_fa}"
        subprocess.call(cmd, shell=True)

        all_affected_cds_files.setdefault(environment, {}).setdefault(MAG_name, []).append(output_fa)

# --- Merge FASTA files per MAG — avoiding duplicated headers ---
for environment, mag_dict in all_affected_cds_files.items():
    for mag_name, fasta_files in mag_dict.items():
        mag_output_fa = os.path.join(out_dir, f"function/MAGs_{environment}/merged_{mag_name}_{environment}.fasta")
        seen_headers = set()

        with open(mag_output_fa, 'w') as merged_file:
            for fasta_file in fasta_files:
                if not os.path.exists(fasta_file):
                    continue

                with open(fasta_file, 'r') as f:
                    sequence_data = f.read().split('>')[1:]  # skip empty first part

                    for seq in sequence_data:
                        header, sequence = seq.split("\n", 1)
                        header = header.strip()
                        sequence = sequence.replace("\n", "")

                        if header not in seen_headers:
                            merged_file.write(f">{header}\n{sequence}\n")
                            seen_headers.add(header)
                        else:
                            print(f"Duplicate sequence skipped in {mag_name}: {header}")

# ---  Function to match contig:start-end with tolerance of ±1 position ---
def matches_header(header, contig, start, end):
    try:
        contig_h, coords = header.split(":")
        start_h, end_h = map(int, coords.split("-"))
        return (
            contig == contig_h and
            abs(start_h - start) <= 1 and
            abs(end_h - end) <= 1
        )
    except Exception:
        return False

# ---  Filter and merge FASTA files per miRNA ---
for mirna_name, env_dict in contig_sets_by_mirna.items():
    for environment, mirna_contigs in env_dict.items():
        filtered_sequences = []
        already_added_headers = set()
        mirna_contigs_file_path = os.path.join(out_dir, f"function/miRNA_{environment}/contigs_for_{mirna_name}.txt")

        os.makedirs(os.path.dirname(mirna_contigs_file_path), exist_ok=True)

        print(f"Processing miRNA: {mirna_name} in {environment}")

        # Save expected contig+coords
        contig_coords_list = []
        with open(mirna_contigs_file_path, 'w') as contig_file:
            for mag_name, contig, start_gene, end_gene in mirna_contigs:
                contig_file.write(f"{contig} {start_gene} {end_gene}\n")
                contig_coords_list.append((contig, start_gene, end_gene))

        # Scan MAG FASTAs
        for mag_name, _, _, _ in mirna_contigs:
            mag_fasta_path = os.path.join(out_dir, f"function/MAGs_{environment}/merged_{mag_name}_{environment}.fasta")

            if not os.path.exists(mag_fasta_path):
                continue

            #print(f"Scanning FASTA: {mag_fasta_path}")

            with open(mag_fasta_path, 'r') as fasta_file:
                sequence_data = fasta_file.read().split('>')[1:]  # Skip first empty item

                for seq in sequence_data:
                    header, sequence = seq.split("\n", 1)
                    header = header.strip()
                    sequence = sequence.replace("\n", "")

                    for contig, start_gene, end_gene in contig_coords_list:
                        if matches_header(header, contig, start_gene, end_gene):
                            if header not in already_added_headers:
                                filtered_sequences.append(f">{header}\n{sequence}\n")
                                already_added_headers.add(header)
                                #print(f"Added unique match: {header}")
                            else:
                                #print(f"Skipping duplicate: {header}")
                            break  # Once matched, stop checking this sequence

        # Save FASTA if matches found
        if filtered_sequences:
            mirna_output_fa = os.path.join(out_dir, f"function/miRNA_{environment}/merged_miRNA_{mirna_name}_{environment}.fasta")
            with open(mirna_output_fa, 'w') as merged_file:
                merged_file.writelines(filtered_sequences)

            #print(f"Saved {len(filtered_sequences)} unique sequences to {mirna_output_fa}")
        else:
            #print(f"WARNING: No matching sequences found for {mirna_name} in {environment}")


# ---  Cleanup temporary files ---
for temp_file in glob.glob(os.path.join(out_dir, "function/*/temp_*.fasta")) + \
                  glob.glob(os.path.join(out_dir, "function/MAGs_*/temp_*.bed")) + \
                  glob.glob(os.path.join(out_dir, "annotation/*/temp_*.gff")):
    os.remove(temp_file)

# ---  Control - temp_merged_affected_cds.fasta "Done!" ---
final_status_file = os.path.join(out_dir, "function/temp_merged_affected_cds.fasta")
with open(final_status_file, 'w') as f:
    f.write("Done!\n")