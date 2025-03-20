#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd
import subprocess
import warnings

# Get input file and output directory from command-line arguments
if len(sys.argv) < 3:
    print("Usage: python impacted.py <input_file> <output_directory>")
    sys.exit(1)

input_file = sys.argv[1]
out_dir = sys.argv[2]

# Create output directories if they do not exist
output_dir_function = os.path.join(out_dir, "function")
os.makedirs(output_dir_function, exist_ok=True)

# Dictionaries to store contigs per miRNA and MAGs per environment
contig_sets_by_mirna = {}
mag_sets_by_environment = {}

print(f"?? Reading input file: {input_file}")

# Check if the input file exists
if not os.path.exists(input_file):
    print(f"?? ERROR: Input file {input_file} does not exist.")
    sys.exit(1)

# Read input file with UTF-8 encoding
with open(input_file, 'r', encoding='utf-8', errors='replace') as file:
    header = file.readline().strip().split('\t')  # Read and identify header
    print(f"?? Detected header ({len(header)} columns): {header}")

    for line_number, line in enumerate(file, start=2):
        columns = line.strip().split('\t')

        # Ensure the line has exactly 15 columns
        if len(columns) != 15:
            print(f"?? ERROR: Line {line_number} has {len(columns)} columns instead of 15. Skipping: {columns}")
            continue

        try:
            mag_name     = columns[0]  # MAG name
            contig_value = columns[1]  # Contig name
            start_gene   = int(columns[11])  # Gene start coordinate
            end_gene     = int(columns[12])  # Gene end coordinate
            mirna_name   = columns[4]  # miRNA name
            environment  = columns[14]  # Environment
        except ValueError as e:
            print(f"?? ERROR: Line {line_number} contains invalid values! {e} -> {columns}")
            continue

        # Store data for miRNAs
        contig_sets_by_mirna.setdefault(mirna_name, {}).setdefault(environment, []).append(
            (mag_name, contig_value, start_gene, end_gene)
        )

        # Store data for MAGs
        mag_sets_by_environment.setdefault(environment, {}).setdefault(mag_name, []).append(
            (contig_value, start_gene, end_gene)
        )

######################
# ?? EXTRACTING FASTA SEQUENCES FOR EACH MAG
######################
all_affected_cds_files = {}

MAG_folders = glob.glob(os.path.join(out_dir, "annotation/*"))

for MAG_folder in MAG_folders:
    MAG_name = os.path.basename(MAG_folder)

    gff_file = os.path.join(MAG_folder, f"{MAG_name}_cds.gff")
    if not os.path.exists(gff_file):
        print(f"?? WARNING: GFF file {gff_file} not found. Skipping {MAG_name}.")
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

        # Criar BED file e contigs_for_* apenas com dados do próprio MAG
        with open(bed_file_path, 'w') as bed_file, open(contig_for_path, 'w') as contig_file:
            for contig, start_gene, end_gene in mag_dict[MAG_name]:
                bed_file.write(f"{contig}\t{start_gene-1}\t{end_gene}\n")  # BED format is 0-based
                contig_file.write(f"{contig} {start_gene} {end_gene}\n")  # Incluir coordenadas no contig_for_*

        input_fasta = os.path.join(MAG_folder, f"{MAG_name}.fna")
        output_fa = os.path.join(output_dir_env, f"temp_{MAG_name}_affected_cds.fasta")

        warnings.filterwarnings("ignore")
        cmd = f"bedtools getfasta -fi {input_fasta} -bed {bed_file_path} -fo {output_fa}"
        subprocess.call(cmd, shell=True)

        all_affected_cds_files.setdefault(environment, {}).setdefault(MAG_name, []).append(output_fa)

# **Merge FASTA files per MAG**
for environment, mag_dict in all_affected_cds_files.items():
    for mag_name, fasta_files in mag_dict.items():
        mag_output_fa = os.path.join(out_dir, f"function/MAGs_{environment}/merged_{mag_name}_{environment}.fasta")
        with open(mag_output_fa, 'w') as merged_file:
            for fasta_file in fasta_files:
                if os.path.exists(fasta_file):
                    with open(fasta_file, 'r') as f:
                        merged_file.write(f.read())

# **Merge FASTA files per miRNA**
for mirna_name, env_dict in contig_sets_by_mirna.items():
    for environment, mirna_contigs in env_dict.items():
        mirna_fasta_files = []

        # Get affected MAGs for this miRNA
        for mag_name, contig, start_gene, end_gene in mirna_contigs:
            mag_fasta_path = os.path.join(out_dir, f"function/MAGs_{environment}/merged_{mag_name}_{environment}.fasta")
            if os.path.exists(mag_fasta_path):
                mirna_fasta_files.append(mag_fasta_path)

        # Merge affected MAGs into one miRNA FASTA
        if mirna_fasta_files:
            mirna_output_fa = os.path.join(out_dir, f"function/miRNA_{environment}/merged_miRNA_{mirna_name}_{environment}.fasta")
            os.makedirs(os.path.dirname(mirna_output_fa), exist_ok=True)
            with open(mirna_output_fa, 'w') as merged_file:
                for fasta_file in mirna_fasta_files:
                    with open(fasta_file, 'r') as f:
                        merged_file.write(f.read())

# Cleanup temporary files
for temp_file in glob.glob(os.path.join(out_dir, "function/*/temp_*.fasta")) + \
                  glob.glob(os.path.join(out_dir, "function/MAGs_*/temp_*.bed")) + \
                  glob.glob(os.path.join(out_dir, "annotation/*/temp_*.gff")):
    os.remove(temp_file)

# control -  temp_merged_affected_cds.fasta "Done!"
final_status_file = os.path.join(out_dir, "function/temp_merged_affected_cds.fasta")
with open(final_status_file, 'w') as f:
    f.write("Done!\n")
            
