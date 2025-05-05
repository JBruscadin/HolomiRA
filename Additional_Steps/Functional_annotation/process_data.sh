#!/bin/bash

if [ "$1" == "--help" ]; then
    echo "Usage: ./process_data.sh <subsystem_level_X.xls> <counts|relab> <sample_type> <dataset_name>"
    echo ""
    echo "Arguments:"
    echo "  <subsystem_level_X.xls>  : The input file containing subsystem information (output from Superfocus)."
    echo "  <counts|relab>           : The type of data to be processed (counts for raw counts, relab for relative abundance)."
    echo "  <sample_type>            : The type of sample (e.g., 'Mags', 'Genes', 'miRNA)'."
    echo "  <dataset_name>           : The name of the organism to be analyzed (e.g., 'Human', 'Ruminants')"
    echo "                           : or the conditions of the same organism (e.g., 'decrease_RME', 'increase_RME')."
    exit 0
fi



# Check if the correct number of arguments was passed
if [ "$#" -ne 4 ]; then
    echo "Usage: ./process_data.sh <subsystem_level_X.xls> <counts|relab> <sample_type> <dataset_name>"
    exit 1
fi

# Capture arguments
input_file=$1
method=$2
sample_type=$3
dataset_name=$4

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Input file not found: $input_file"
    exit 1
fi

# Extract the level from the file name (level_1, level_2, or level_3)
level=$(echo "$input_file" | grep -oP 'level_\d')
if [ -z "$level" ]; then
    echo "Invalid file name. The file name must contain 'level_1', 'level_2', or 'level_3'."
    exit 1
fi

# Create a directory for temporary files if it doesn't exist
mkdir -p Temp/${level}_${sample_type}_${dataset_name}

# Process column with pathway names
tail -n +5 "$input_file" | awk -F'\t' '{print $1}' | sed 's/ /_/g' > Temp/${level}_${sample_type}_${dataset_name}/output_file_${level}_${sample_type}_${dataset_name}.txt
awk '{for (i=1; i<=NF; i++) a[i, NR] = $i} NF > p { p = NF } END { for (i=1; i<=p; i++) { for (j=1; j<=NR; j++) printf "%s ", a[i, j]; print "" }}' Temp/${level}_${sample_type}_${dataset_name}/output_file_${level}_${sample_type}_${dataset_name}.txt > Temp/${level}_${sample_type}_${dataset_name}/output_file_${level}_${sample_type}_${dataset_name}_transposed.txt

# Generate file with column names
tail -n +5 "$input_file" | awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i ~ /%/) printf "%s\n", $i}' > Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_columns_names.txt

# Transpose the original data
tail -n +5 "$input_file" | awk -F '\t' '{for (i=1; i<=NF; i++) a[i, NR] = $i} NF > p { p = NF } END { for (i=1; i<=p; i++) { for (j=1; j<=NR; j++) printf "%s ", a[i, j]; print "" }}' > Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_transposed.txt

# Processing logic based on the chosen method
if [ "$method" == "relab" ]; then
    # Relative Abundance method
    cp Temp/${level}_${sample_type}_${dataset_name}/output_file_${level}_${sample_type}_${dataset_name}_transposed.txt Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_filter_relab.txt
    grep -w -f Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_columns_names.txt Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_transposed.txt >> Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_filter_relab.txt
    sed 's/% //g' Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_filter_relab.txt > Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_filter_relab_cleaned.txt
    awk -F ' ' '{for (i=1; i<=NF; i++) a[i, NR] = $i} NF > p { p = NF } END { for (i=1; i<=p; i++) { for (j=1; j<=NR; j++) printf "%s ", a[i, j]; print "" }}' Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_filter_relab_cleaned.txt > Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_final_relab.txt
    sed "s/'//g; s/&#//g" Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_final_relab.txt  > ${level}_${sample_type}_${dataset_name}_relab.txt
    echo "Relative Abundance processing completed for ${level} with sample type ${sample_type} of ${dataset_name}."

elif [ "$method" == "counts" ]; then
    # Counts method
    cp Temp/${level}_${sample_type}_${dataset_name}/output_file_${level}_${sample_type}_${dataset_name}_transposed.txt Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_filter_counts.txt
    grep -w -v -f Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_columns_names.txt Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_transposed.txt >> Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_filter_counts.txt
    sed -i '2d' Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_filter_counts.txt
    awk -F ' ' '{for (i=1; i<=NF; i++) a[i, NR] = $i} NF > p { p = NF } END { for (i=1; i<=p; i++) { for (j=1; j<=NR; j++) printf "%s ", a[i, j]; print "" }}' Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_filter_counts.txt > Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_final_counts.txt
    sed "s/'//g; s/&#//g" Temp/${level}_${sample_type}_${dataset_name}/${level}_${sample_type}_${dataset_name}_final_counts.txt  > ${level}_${sample_type}_${dataset_name}_counts.txt
    echo "Counts processing completed for ${level} with sample type ${sample_type} of ${dataset_name}."

else
    echo "Invalid option. Choose 'counts' or 'relab'."
    exit 1
fi

# Clean up intermediate files (optional)
# Uncomment the following line if you want to remove the Temp directory after the script runs:
# rm -r Temp
