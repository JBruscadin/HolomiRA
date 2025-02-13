# -*- coding: utf-8 -*-

# Import the necessary modules
import subprocess
import pandas as pd
import glob
import os
import sys
import warnings

# Get the input and output file paths from the command line arguments
input_file = sys.argv[1]
out_dir = sys.argv[2]

# Check if the output folder exists, and if not, create it
output_dir_mirna_base = os.path.join(out_dir, "function")
if not os.path.exists(output_dir_mirna_base):
    os.makedirs(output_dir_mirna_base)

# Initialize dictionaries to store contigs per miRNA and MAGs by Environment
contig_sets_by_mirna = {}
mag_sets_by_environment = {}

# Open the input file and read its lines
with open(input_file, 'r') as file:
    # Use enumerate to keep track of the line number
    for line_number, line in enumerate(file):
        # Skip the header in the first line
        if line_number == 0:
            continue
        # Split the line into columns using tab as the separator
        columns = line.strip().split('\t')
        # Extract relevant columns
        contig_value = columns[1]
        mirna_name = columns[4]
        mag_name = columns[0]
        environment = columns[10]

        # Handle miRNA related files and directories
        if mirna_name not in contig_sets_by_mirna:
            contig_sets_by_mirna[mirna_name] = {}
        
        # Add environment-specific contigs for miRNA
        if environment not in contig_sets_by_mirna[mirna_name]:
            contig_sets_by_mirna[mirna_name][environment] = set()
        contig_sets_by_mirna[mirna_name][environment].add(contig_value)

        # Handle MAG-related files and directories
        if environment not in mag_sets_by_environment:
            mag_sets_by_environment[environment] = {}
        if mag_name not in mag_sets_by_environment[environment]:
            mag_sets_by_environment[environment][mag_name] = set()
        mag_sets_by_environment[environment][mag_name].add(contig_value)

# Create miRNA and MAG directories and save files
for mirna_name, env_dict in contig_sets_by_mirna.items():
    for environment, contig_set in env_dict.items():
        output_dir_mirna = os.path.join(out_dir, f"function/miRNA_{environment}")
        if not os.path.exists(output_dir_mirna):
            os.makedirs(output_dir_mirna)
        filename = os.path.join(output_dir_mirna, f"contigs_for_{mirna_name}.txt")
        with open(filename, 'w') as output_file:
            output_file.write('\n'.join(contig_set))

for environment, mag_dict in mag_sets_by_environment.items():
    output_dir_mag = os.path.join(out_dir, f"function/MAGs_{environment}")
    if not os.path.exists(output_dir_mag):
        os.makedirs(output_dir_mag)
    
    for mag_name, contig_set in mag_dict.items():
        filename = os.path.join(output_dir_mag, f"contigs_for_{mag_name}.txt")
        with open(filename, 'w') as output_file:
            output_file.write('\n'.join(contig_set))

# Initialize a list to store the paths of all affected CDS files
all_affected_cds_files = []

# List all animal-specific subdirectories
MAG_folders = glob.glob(os.path.join(out_dir, "annotation/*"))

# Iterate over the MAG-specific subdirectories
for MAG_folder in MAG_folders:
    # Extract the MAG name from the folder path
    MAG_name = os.path.basename(MAG_folder)

    # Read the GFF file into a Pandas DataFrame
    gff_file = os.path.join(MAG_folder, f"{MAG_name}_cds.gff")
    df = pd.read_csv(
        gff_file, sep='\t', comment='#', header=None,
        names=["contig", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    )
    df["contig"] = df["contig"].str.strip()  # Normalize contig names

    # Process each environment
    for environment, mag_dict in mag_sets_by_environment.items():
        if MAG_name not in mag_dict:
            continue

        # List contig list files (contigs_for_*.txt)
        contig_list_files = glob.glob(os.path.join(out_dir, f"function/*_{environment}/contigs_for_*.txt"))

        # Process each contig list file
        for contig_list_file in contig_list_files:
            # Extract the wildcard part (*) from the filename
            wildcard_part = os.path.basename(contig_list_file).replace("contigs_for_", "").replace(".txt", "")

            # Read the list of contigs from the file
            with open(contig_list_file, 'r') as f:
                contigs_to_keep = [line.strip() for line in f]

            # Filter the DataFrame based on the list of contigs
            filtered_df = df[df["contig"].isin(contigs_to_keep)]

            if not filtered_df.empty:  # Only proceed if there are matches
                # Save the filtered GFF entries to a new GFF file
                filtered_gff_file = os.path.join(
                    MAG_folder,
                    f"temp_{MAG_name}_filtered_miRNA_{wildcard_part}_{environment}.gff"
                )
                filtered_df.to_csv(filtered_gff_file, sep='\t', index=False, header=False)

                # Create the command to extract affected CDS sequences
                input_fasta = os.path.join(MAG_folder, f"{MAG_name}.fna")
                output_fa = os.path.join(
                    MAG_folder,
                    f"temp_{MAG_name}_affected_cds_miRNA_{wildcard_part}_{environment}.fasta"
                )

                warnings.filterwarnings("ignore")
                cmd = f"bedtools getfasta -fi {input_fasta} -bed {filtered_gff_file} -fo {output_fa}"
                subprocess.call(cmd, shell=True)

                # Append the path of the generated affected CDS file to the list
                all_affected_cds_files.append(output_fa)
            else:
                print(f"No matching contigs found in {contig_list_file} for MAG {MAG_name}")

# After processing all files, merge the affected CDS files
output_dir = os.path.join(out_dir, "function/")
merged_affected_cds_file = os.path.join(output_dir, "temp_merged_affected_cds.fasta")
with open(merged_affected_cds_file, 'w') as output_file:
    for affected_cds_file in all_affected_cds_files:
        with open(affected_cds_file, 'r') as input_file:
            output_file.write(input_file.read())
