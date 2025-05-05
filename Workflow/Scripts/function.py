# -*- coding: utf-8 -*-

# Import necessary modules
import subprocess
import pandas as pd
import glob
import os
import sys

# Get input and output file paths from command line arguments
input_file = sys.argv[1]
out_dir = sys.argv[2]

# Check if the output folder exists, and if not, create it
output_dir = os.path.join(out_dir, "function/miRNA")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Initialize a dictionary to store sets of unique contigs per miRNA
contig_sets_by_mirna = {}

# Open the input file and read its lines
with open(input_file, 'r') as file:
    # Use enumerate to keep track of the line number
    for line_number, line in enumerate(file):
        # Skip the header in the first line
        if line_number == 0:
            continue
        # Split the line into columns using tab as the separator
        columns = line.strip().split('\t')
        # Get the value from the "Contig" column (column 1)
        contig_value = columns[1]
        # Get the miRNA value (column 4)
        mirna_name = columns[4]
        # Check if the miRNA is already in the dictionary, if not, create an empty set
        if mirna_name not in contig_sets_by_mirna:
            contig_sets_by_mirna[mirna_name] = set()
        # Add the contig to the corresponding miRNA's set
        contig_sets_by_mirna[mirna_name].add(contig_value)

# Save the sets of unique contigs in separate files for each miRNA
for mirna_name, contig_set in contig_sets_by_mirna.items():
    filename = os.path.join(output_dir, f"contigs_for_{mirna_name}.txt")
    with open(filename, 'w') as output_file:
        output_file.write('\n'.join(contig_set))

# List all animal-specific subdirectories
MAG_folders = glob.glob(os.path.join(out_dir, "annotation/*"))

# Iterate over the MAG-specific subdirectories
for MAG_folder in MAG_folders:
    # Extract the MAG name from the folder path
    MAG_name = os.path.basename(MAG_folder)

    # Read the GFF file into a Pandas DataFrame
    gff_file = os.path.join(MAG_folder, f"{MAG_name}_cds.gff")
    df = pd.read_csv(gff_file, sep='\t', comment='#', header=None, names=["contig", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])

    # SList contig list files (contigs_for_*.txt)
    contig_list_files = glob.glob(os.path.join(output_dir, "contigs_for_*.txt"))

    # Process each contig list file
    for contig_list_file in contig_list_files:
        # Extract the wildcard part (*) from the filename
        wildcard_part = os.path.basename(contig_list_file).replace("contigs_for_", "").replace(".txt", "")

        # Read the list of contigs from the file
        with open(contig_list_file, 'r') as f:
            contigs_to keep = [line.strip() for line in f]

        # Filter the DataFrame based on the list of contigs
        filtered_df = df[df["contig"].isin(contigs_to_keep)]

        # Save the filtered GFF entries to a new GFF file
        filtered_gff_file = os.path.join(MAG_folder, f"temp_{MAG_name}_filtered_miRNA_{wildcard_part}.gff")
        filtered_df.to_csv(filtered_gff_file, sep='\t', index=False, header=False)

        # Create the command to extract affected CDS sequences
        input_fasta = os.path.join(MAG_folder, f"{MAG_name}.fna")
        input_bed = filtered_gff_file
        output_fa = os.path.join(MAG_folder, f"temp_{MAG_name}_affected_cds_miRNA_{wildcard_part}.fasta")

        cmd = f"bedtools getfasta -fi {input_fasta} -bed {input_bed} -fo {output_fa}"
        subprocess.call(cmd, shell=True)
