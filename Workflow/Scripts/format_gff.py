# -*- coding: utf-8 -*-
import os
import sys

# Check if the correct number of arguments was provided
if len(sys.argv) != 2:
    print("Usage: python format_gff.py OUT_DIR")
    sys.exit(1)

# Get the value of OUT_DIR from the command line arguments
OUT_DIR = sys.argv[1]

# Set the root directory where the folders with *_cds.gff files are located
root_dir = os.path.join(OUT_DIR, "annotation")

# Iterate over the folders within the root directory
for dir_path, dir_names, file_names in os.walk(root_dir):
    for dir_name in dir_names:
        dir_full_path = os.path.join(dir_path, dir_name)

        # Execute the command for each *_cds.gff file within the directory
        for file_name in os.listdir(dir_full_path):
            if file_name.endswith("_cds.gff"):
                gff_file_path = os.path.join(dir_full_path, file_name)
                output_file_path = os.path.join(dir_full_path, file_name.replace("_cds.gff", "_IDs.txt"))

                with open(gff_file_path, "r") as gff_file, open(output_file_path, "w") as output_file:
                    for line in gff_file:
                        if not line.startswith("#"):
                            columns = line.strip().split("\t")
                            if len(columns) >= 9:
                                column_values = columns[0], columns[3], columns[4], columns[6], columns[8]
                                value = column_values[4].split(";", 1)[0]
                                output_line = "\t".join(column_values[:4] + [value]) + "\n"
                                output_file.write(output_line)
