#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd
import subprocess
import warnings

# --- Get input file and output directory from command-line arguments ---
if len(sys.argv) < 3:
    print("Usage: impacted.py <input_file> <output_dir>")
    sys.exit(1)

input_file = sys.argv[1]
out_dir = sys.argv[2]

# --- Create output directory if it doesn't exist ---
output_dir_function = os.path.join(out_dir, "function")
os.makedirs(output_dir_function, exist_ok=True)

# --- Initialize dictionaries to store contigs per miRNA and MAGs per environment ---
contig_sets_by_mirna = {}
mag_sets_by_environment = {}

print(f"Reading input file: {input_file}")

# --- Check file existence ---
if not os.path.exists(input_file):
    print(f"ERROR: Input file {input_file} does not exist.")
    sys.exit(1)

# --- Read the input TSV file using pandas with python engine ---
try:
    df_input = pd.read_csv(input_file, sep='\t', engine='python')
    if df_input.shape[1] != 22:
        print(f"ERROR: Expected 22 columns, but found {df_input.shape[1]}.")
        print(f"Columns detected: {list(df_input.columns)}")
        sys.exit(1)
    print(f"Detected header ({df_input.shape[1]} columns): {list(df_input.columns)}")
except Exception as e:
    print(f"ERROR while reading input file: {e}")
    sys.exit(1)

# --- Parse each row ---
for line_number, row in df_input.iterrows():
    try:
        mag_name = row['MAG']
        contig_value = row['Contig']
        start_gene = int(row['start_gene']) if 'start_gene' in row else int(row[11])
        end_gene = int(row['end_gene']) if 'end_gene' in row else int(row[12])
        mirna_name = row['miRNA']
        environment = row['Environment']
    except Exception as e:
        print(f"ERROR: Line {line_number + 2} contains invalid data: {e}")
        continue

    contig_sets_by_mirna.setdefault(mirna_name, {}).setdefault(environment, []).append(
        (mag_name, contig_value, start_gene, end_gene)
    )
    mag_sets_by_environment.setdefault(environment, {}).setdefault(mag_name, []).append(
        (contig_value, start_gene, end_gene)
    )

# --- Extract affected CDS regions with bedtools ---
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

        with open(bed_file_path, 'w') as bed_file, open(contig_for_path, 'w') as contig_file:
            for contig, start_gene, end_gene in mag_dict[MAG_name]:
                bed_file.write(f"{contig}\t{start_gene - 1}\t{end_gene}\n")
                contig_file.write(f"{contig} {start_gene} {end_gene}\n")

        input_fasta = os.path.join(MAG_folder, f"{MAG_name}.fna")
        output_fa = os.path.join(output_dir_env, f"temp_{MAG_name}_affected_cds.fasta")

        warnings.filterwarnings("ignore")
        cmd = f"bedtools getfasta -fi {input_fasta} -bed {bed_file_path} -fo {output_fa}"
        subprocess.call(cmd, shell=True)

        all_affected_cds_files.setdefault(environment, {}).setdefault(MAG_name, []).append(output_fa)

# --- Merge FASTA files per MAG avoiding duplicates ---
for environment, mag_dict in all_affected_cds_files.items():
    for mag_name, fasta_files in mag_dict.items():
        mag_output_fa = os.path.join(out_dir, f"function/MAGs_{environment}/merged_{mag_name}_{environment}.fasta")
        seen_headers = set()

        with open(mag_output_fa, 'w') as merged_file:
            for fasta_file in fasta_files:
                if not os.path.exists(fasta_file):
                    continue
                with open(fasta_file, 'r') as f:
                    sequence_data = f.read().split('>')[1:]
                    for seq in sequence_data:
                        header, sequence = seq.split("\n", 1)
                        header = header.strip()
                        sequence = sequence.replace("\n", "")
                        if header not in seen_headers:
                            merged_file.write(f">{header}\n{sequence}\n")
                            seen_headers.add(header)
                        else:
                            print(f"Duplicate sequence skipped in {mag_name}: {header}")

# --- Header matching function with ±1 tolerance ---
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

# --- Filter and merge per miRNA ---
for mirna_name, env_dict in contig_sets_by_mirna.items():
    for environment, mirna_contigs in env_dict.items():
        filtered_sequences = []
        already_added_headers = set()
        mirna_dir = os.path.join(out_dir, f"function/miRNA_{environment}")
        os.makedirs(mirna_dir, exist_ok=True)
        contig_coords_list = []

        contig_path = os.path.join(mirna_dir, f"contigs_for_{mirna_name}.txt")
        with open(contig_path, 'w') as contig_file:
            for mag_name, contig, start_gene, end_gene in mirna_contigs:
                contig_file.write(f"{contig} {start_gene} {end_gene}\n")
                contig_coords_list.append((contig, start_gene, end_gene))

        for mag_name, _, _, _ in mirna_contigs:
            mag_fasta_path = os.path.join(out_dir, f"function/MAGs_{environment}/merged_{mag_name}_{environment}.fasta")
            if not os.path.exists(mag_fasta_path):
                continue

            with open(mag_fasta_path, 'r') as fasta_file:
                sequence_data = fasta_file.read().split('>')[1:]
                for seq in sequence_data:
                    header, sequence = seq.split("\n", 1)
                    header = header.strip()
                    sequence = sequence.replace("\n", "")
                    for contig, start_gene, end_gene in contig_coords_list:
                        if matches_header(header, contig, start_gene, end_gene):
                            if header not in already_added_headers:
                                filtered_sequences.append(f">{header}\n{sequence}\n")
                                already_added_headers.add(header)
                                print(f"Added unique match: {header}")
                            else:
                                print(f"Skipping duplicate: {header}")
                            break

        if filtered_sequences:
            output_fasta = os.path.join(mirna_dir, f"merged_miRNA_{mirna_name}_{environment}.fasta")
            with open(output_fasta, 'w') as merged_file:
                merged_file.writelines(filtered_sequences)
            print(f"Saved {len(filtered_sequences)} unique sequences to {output_fasta}")
        else:
            print(f"WARNING: No matching sequences found for {mirna_name} in {environment}")

# --- Cleanup temporary files ---
for temp_file in glob.glob(os.path.join(out_dir, "function/*/temp_*.fasta")) + \
                  glob.glob(os.path.join(out_dir, "function/MAGs_*/temp_*.bed")) + \
                  glob.glob(os.path.join(out_dir, "annotation/*/temp_*.gff")):
    os.remove(temp_file)

# --- Final status file ---
final_status_file = os.path.join(out_dir, "function/temp_merged_affected_cds.fasta")
with open(final_status_file, 'w') as f:
    f.write("Done!\n")
