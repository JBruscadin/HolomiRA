from typing import Dict
import os

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import seaborn as sns
import shutil

## HolomiRA: Venn diagram

comparation = input(
    "Enter 1 to compare groups between samples and 2 to compare species: "
)

def ensure_valid_directory(final_results_dir: str) -> str:
    """
    Ensure a valid directory path is provided, allowing the user to overwrite if it exists 
    or specify a new directory.

    Args:
        final_results_dir (str): The full path to the directory to be created or overwritten.

    Returns:
        str: The full path to the valid directory.
    """
    while True:
        # Check if the directory already exists
        if os.path.exists(final_results_dir):
            # Warn the user and ask if they want to overwrite
            user_input = input(f"Warning: The directory '{final_results_dir}' already exists. Do you want to overwrite it? (yes/no): ").strip().lower()

            if user_input in ['yes', 'y']:
                # Remove the existing directory and create a new one
                shutil.rmtree(final_results_dir)
                os.makedirs(final_results_dir)
                print(f"Directory '{final_results_dir}' has been overwritten.")
                break  # Exit the loop as the path is valid and processed
            else:
                # Ask the user to enter a new path
                final_results_dir = input("Please enter a new directory path: ").strip()
        else:
            # Create the directory if it doesn't exist
            os.makedirs(final_results_dir)
            print(f"Directory '{final_results_dir}' has been created.")
            break  # Exit the loop as the directory is successfully created

    return final_results_dir

def create_combined_venn_diagram(df):

    def extract_up_to_level(taxonomy_string, level_prefix):
            """
            Extracts the part of the taxonomy string up to and including the specified level prefix.
            If the prefix is not found, returns None.
            """
            parts = taxonomy_string.split(';')
            result = []
            for part in parts:
                result.append(part)
                if part.startswith(level_prefix):
                    parts = part
                    return parts
            # if level_prefix in result[-1]:  # Ensure the level is in the taxonomy string
            #     return ';'.join(result)  # Join the extracted parts back into a string
            # else:
            #     return None
    # List of valid taxonomic levels
    valid_taxonomic_levels = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    # Validate the user input
    while True:
        # Ask user to input the desired taxonomy level (e.g., d__, p__, c__)
        desired_level = input("Enter the taxonomy level (e.g., 'd__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__'): ")
        # Validate the user input
        if desired_level in valid_taxonomic_levels:
            break  # Exit the loop if the input is valid
        else:
            print(
                f"Invalid taxonomy level: {desired_level}. "
                f"Must be one of {valid_taxonomic_levels}. Please try again."
            )
    df['Taxonomy'] = df['Taxonomy'].apply(lambda x: extract_up_to_level(x, desired_level))

    # Mapping of taxonomy prefixes to their full names
    taxonomy_level_names = {
        'd__': 'Domain',
        'p__': 'Phylum',
        'c__': 'Class',
        'o__': 'Order',
        'f__': 'Family',
        'g__': 'Genus',
        's__': 'Species'
    }

    # Get the full name of the taxonomy level based on the user input
    taxonomy_name = taxonomy_level_names.get(desired_level, "Unknown Level")

    # Get the unique values in the 'Environment' column
    unique_environments = df["Environment"].unique()
    # Create sets for unique miRNA, Gene, and Taxonomy for each group
    mirna_sets = [
            set(df[df["Environment"] == env]["miRNA"]) for env in unique_environments
    ]
    gene_sets = [
        set(df[df["Environment"] == env]["Gene"]) for env in unique_environments
    ]
    # Discard rows where the taxonomy level is missing (i.e., where extract_up_to_level returned None)
    # df = df.dropna(subset=['Taxonomy'])
    taxonomy_sets = [
        set(df[df["Environment"] == env]["Taxonomy"]) for env in unique_environments
    ]
    # Open a file to save the summary of unique and shared elements
    with open(f"{final_results_dir}/venn_summary.txt", "w") as f:
        # Create a single figure with 3 subplots
        _, axs = plt.subplots(
            1, 3, figsize=(18, 6)
        )  # Set the figure size for the combined plot
        # Create Venn diagrams based on the number of Groups
        num_Groups = len(unique_environments)
        if num_Groups == 2:
            # Venn for miRNA
            venn2(
                [mirna_sets[0], mirna_sets[1]],
                set_labels=unique_environments,
                set_colors=[palette[0], palette[1]],
                ax=axs[0],
            )
            axs[0].set_title("Unique miRNA", fontsize=14)
            # Write miRNA unique/shared sets to file
            f.write("### miRNA\n")
            f.write(
                f"Unique to {unique_environments[0]}: {mirna_sets[0] - mirna_sets[1]}\n"
            )
            f.write(
                f"Unique to {unique_environments[1]}: {mirna_sets[1] - mirna_sets[0]}\n"
            )
            f.write(f"Shared: {mirna_sets[0] & mirna_sets[1]}\n\n")
            # Venn for Gene
            venn2(
                [gene_sets[0], gene_sets[1]],
                set_labels=unique_environments,
                set_colors=[palette[0], palette[1]],
                ax=axs[1],
            )
            axs[1].set_title("Unique Gene", fontsize=14)
            # Write Gene unique/shared sets to file
            f.write("### Gene\n")
            f.write(
                f"Unique to {unique_environments[0]}: {gene_sets[0] - gene_sets[1]}\n"
            )
            f.write(
                f"Unique to {unique_environments[1]}: {gene_sets[1] - gene_sets[0]}\n"
            )
            f.write(f"Shared: {gene_sets[0] & gene_sets[1]}\n\n")
            # Venn for Taxonomy
            venn2(
                [taxonomy_sets[0], taxonomy_sets[1]],
                set_labels=unique_environments,
                set_colors=[palette[0], palette[1]],
                ax=axs[2],
            )
            axs[2].set_title(f"Unique Taxonomy ({taxonomy_name})", fontsize=14)
            # Write Taxonomy unique/shared sets to file
            f.write("### Taxonomy\n")
            f.write(
                f"Unique to {unique_environments[0]}: {taxonomy_sets[0] - taxonomy_sets[1]}\n"
            )
            f.write(
                f"Unique to {unique_environments[1]}: {taxonomy_sets[1] - taxonomy_sets[0]}\n"
            )
            f.write(f"Shared: {taxonomy_sets[0] & taxonomy_sets[1]}\n\n")
        elif num_Groups == 3:
            # Venn for miRNA
            venn3(
                [mirna_sets[0], mirna_sets[1], mirna_sets[2]],
                set_labels=unique_environments,
                set_colors=[palette[0], palette[1], palette[2]],
                ax=axs[0],
            )
            axs[0].set_title("Unique miRNA", fontsize=14)
            # Write miRNA unique/shared sets to file
            f.write("### miRNA\n")
            f.write(
                f"Unique to {unique_environments[0]}: {mirna_sets[0] - mirna_sets[1] - mirna_sets[2]}\n"
            )
            f.write(
                f"Unique to {unique_environments[1]}: {mirna_sets[1] - mirna_sets[0] - mirna_sets[2]}\n"
            )
            f.write(
                f"Unique to {unique_environments[2]}: {mirna_sets[2] - mirna_sets[0] - mirna_sets[1]}\n"
            )
            f.write(f"Shared: {mirna_sets[0] & mirna_sets[1] & mirna_sets[2]}\n\n")
            # Venn for Gene
            venn3(
                [gene_sets[0], gene_sets[1], gene_sets[2]],
                set_labels=unique_environments,
                set_colors=[palette[0], palette[1], palette[2]],
                ax=axs[1],
            )
            axs[1].set_title("Unique Gene", fontsize=14)
            # Write Gene unique/shared sets to file
            f.write("### Gene\n")
            f.write(
                f"Unique to {unique_environments[0]}: {gene_sets[0] - gene_sets[1] - gene_sets[2]}\n"
            )
            f.write(
                f"Unique to {unique_environments[1]}: {gene_sets[1] - gene_sets[0] - gene_sets[2]}\n"
            )
            f.write(
                f"Unique to {unique_environments[2]}: {gene_sets[2] - gene_sets[0] - gene_sets[1]}\n"
            )
            f.write(f"Shared: {gene_sets[0] & gene_sets[1] & gene_sets[2]}\n\n")
            # Venn for Taxonomy
            venn3(
                [taxonomy_sets[0], taxonomy_sets[1], taxonomy_sets[2]],
                set_labels=unique_environments,
                set_colors=[palette[0], palette[1], palette[2]],
                ax=axs[2],
            )
            axs[2].set_title(f"Unique Taxonomy ({taxonomy_name})", fontsize=14)
            # Write Taxonomy unique/shared sets to file
            f.write("### Taxonomy\n")
            f.write(
                f"Unique to {unique_environments[0]}: {taxonomy_sets[0] - taxonomy_sets[1] - taxonomy_sets[2]}\n"
            )
            f.write(
                f"Unique to {unique_environments[1]}: {taxonomy_sets[1] - taxonomy_sets[0] - taxonomy_sets[2]}\n"
            )
            f.write(
                f"Unique to {unique_environments[2]}: {taxonomy_sets[2] - taxonomy_sets[0] - taxonomy_sets[1]}\n"
            )
            f.write(
                f"Shared: {taxonomy_sets[0] & taxonomy_sets[1] & taxonomy_sets[2]}\n\n"
            )
    # Adjust the layout to ensure titles are aligned
    plt.subplots_adjust(top=0.85, wspace=0.3)
    # Set a super title for the whole figure (optional)
    plt.suptitle("Venn Diagrams for miRNA, Gene, and Taxonomy", fontsize=16)
    # Save the combined figure with the three Venn diagrams
    plt.savefig(f"{final_results_dir}/Venn_diagram_combined.png")
    plt.close()  # Close the plot to avoid display issues


if comparation == "1":
    input_file = input("Please enter the path of the directory with the file: ")
    out_dir = input(
        "Please enter the path to the output directory (leave blank to create one): "
    )

    # Read the all_results.tsv file into a pandas DataFrame
    df = pd.read_csv(input_file, sep="\t")

    # Create the output directory if it doesn't exist
    final_results_dir = os.path.join(out_dir, "results_venn_diagram")
    if not os.path.exists(final_results_dir):
        os.makedirs(final_results_dir)

    final_results_dir = ensure_valid_directory(final_results_dir)

    # Extract the "deep" color palette from Seaborn
    palette = sns.color_palette("deep")

    # Get the unique values in the 'Environment' column
    unique_environments = df["Environment"].unique()

    # Check the number of unique environments
    num_environments = len(unique_environments)

    # Check the number of environments and call the function accordingly
    if num_environments == 1:
        print("Error: Only one group found.")
    elif num_environments > 4:
        print("Error: More than four groups found.")
    else:
        create_combined_venn_diagram(
            df
        )
        print(f"Combined Venn diagrams generated and saved successfully in {out_dir}")
        print(
            f"Venn summary saved successfully in {final_results_dir}/venn_summary.txt."
        )

elif comparation == "2":

    input_file1 = input("Please enter the path of the directory of the first sample: ")
    input_file2 = input("Please enter the path of the directory of the second sample: ")

    out_dir = input(
        "Please enter the path to the output directory (leave blank to create one): "
    )

    df1 = pd.read_csv(input_file1, sep="\t")
    df1['miRNA'] = df1['miRNA'].str.split('-', n=1).str[1]

    df2 = pd.read_csv(input_file2, sep="\t")
    df2['miRNA'] = df2['miRNA'].str.split('-', n=1).str[1]

    # Concatenate the two DataFrames by rows
    df = pd.concat([df1, df2], ignore_index=True)

    # Get the unique values in the 'Environment' column
    unique_groups = df["Environment"].unique()

    # Create the output directory if it doesn't exist
    final_results_dir = os.path.join(out_dir, "results_venn_diagram")
    final_results_dir = ensure_valid_directory(final_results_dir)

    # Extract the "deep" color palette from Seaborn
    palette = sns.color_palette("deep")

    # Check the number of Groups and call the function accordingly
    if len(unique_groups) == 1:
        print("Error: Only one environment found.")
    elif len(unique_groups) > 3:
        print("Error: More than three environments found.")
    else:
        create_combined_venn_diagram(
            df
        )
        print(f"Combined Venn diagrams generated and saved successfully in {out_dir}")
        print(
            f"Venn summary saved successfully in {final_results_dir}/venn_summary.txt."
        )

else:
    print("Error: variable is not 1 or 2")
