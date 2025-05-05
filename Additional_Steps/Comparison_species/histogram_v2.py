import os
import shutil

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt


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
            user_input = input(
                f"Warning: The directory '{final_results_dir}' already exists. "
                "Do you want to overwrite it? (yes/no): ").strip().lower()

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
##HolomiRA: Histograms

comparation = input(
    "Enter 1 to compare groups between samples and 2 to compare species: "
)

if comparation == "1":
    input_file = input("Please enter the path of the directory with the file: ")
    out_dir = input(
        "Please enter the path to the output directory (leave blank to create one): "
    )

    if not out_dir:
        out_dir = os.path.join(os.getcwd(), "output")

    # Create the output directory if it doesn't exist
    final_results_dir = os.path.join(out_dir, "results_histogram")
    final_results_dir = ensure_valid_directory(final_results_dir)

    df = pd.read_csv(input_file, sep='\t')

     # Calculate the total number of unique miRNA, gene, MAG and Taxonomy for each Environment
    counts = df.groupby(['Environment'])[['miRNA', 'Gene', 'MAG', 'Taxonomy']].nunique()

    # Reset index to use 'Environment' as column
    counts.reset_index(inplace=True)

    # Unpivot the data to use 'Variable' for the types (miRNA, Gene, MAG, Taxonomy)
    counts = pd.melt(counts, id_vars='Environment', var_name='Variable', value_name='Counts')

    # Plot the barplot
    plt.figure(figsize=(9, 6))
    sns.set(style="white")  # No grid lines

    # Create the plot using Seaborn
    plt.figure(figsize=(12, 6))
    sns.barplot(data=counts, x='Variable', y='Counts', hue='Environment', palette='deep')
    plt.title('Unique Number of miRNAs, Genes, MAGs, and Taxonomies', fontsize='16')
    plt.xlabel('Category', fontsize=14)
    plt.ylabel('Unique Count', fontsize=14)

    # Add count numbers on top of each bar
    for p in plt.gca().patches:
        height = int(p.get_height())
        if height > 0:
                plt.gca().annotate(f'{height}', (p.get_x() + p.get_width() / 2., height / 2.),
                        ha='center', va='center', fontsize=12, color='black', xytext=(0, 5),
                        textcoords='offset points')


    # Save the plot
    plt.savefig(f'{final_results_dir}/histogram.png')
    print("Histograms generated and saved successfully.")

elif comparation == "2":

    input_file1 = input("Please enter the path of the directory of the first sample: ")
    input_file2 = input("Please enter the path of the directory of the second sample: ")

    out_dir = input(
        "Please enter the path to the output directory (leave blank to create one): "
    )

    if not out_dir:
        out_dir = os.path.join(os.getcwd(), "output")

    # Create the output directory if it doesn't exist
    final_results_dir = os.path.join(out_dir, "results_histogram")
    final_results_dir = ensure_valid_directory(final_results_dir)

    df1 = pd.read_csv(input_file1, sep="\t")

    df2 = pd.read_csv(input_file2, sep="\t")

    # Concatenate the two DataFrames by rows
    df = pd.concat([df1, df2], ignore_index=True)

     # Calculate the total number of unique miRNA, gene, MAG and Taxonomy for each Environment
    counts = df.groupby(['Environment'])[['miRNA', 'Gene', 'MAG', 'Taxonomy']].nunique()

    # Reset index to use 'Environment' as column
    counts.reset_index(inplace=True)

    # Unpivot the data to use 'Variable' for the types (miRNA, Gene, MAG, Taxonomy)
    counts = pd.melt(counts, id_vars='Environment', var_name='Variable', value_name='Counts')

    # fig, axes = plt.subplots(2, 2, tight_layout=True) 

    # # Flatten the axes array to simplify indexing
    # axes = axes.flatten()

    # # Loop through each category to create individual plots
    # for i, category in enumerate(counts['Variable'].unique()):
    #     # Filter data for the current category
    #     category_data = counts[counts['Variable'] == category]

    #     # Create the barplot on the corresponding subplot
    #     sns.barplot(data=category_data, x='Variable', y='Counts', hue='Environment', ax=axes[i])

    #     # Customize the plot
    #     axes[i].set_title(f'Counts for {category}', fontsize=14)
    #     axes[i].set_xlabel('')
    #     axes[i].set_ylabel('Counts', fontsize=12)

    #     # Remove legend from the subplot
    #     axes[i].legend_.remove()

    #     # Add count numbers in the middle of each bar
    #     for p in axes[i].patches:
    #         height = int(p.get_height())
    #         if height > 0:
    #             axes[i].annotate(f'{height}', (p.get_x() + p.get_width() / 2., height / 2.),
    #                     ha='center', va='center', fontsize=12, color='black', xytext=(0, 5),
    #                     textcoords='offset points')
                
    # # Add a single legend on the right
    # handles, labels = axes[0].get_legend_handles_labels()  # Get legend info from the first subplot
    # fig.legend(handles, labels, loc='upper left', title='Environment', fontsize=12, title_fontsize=14, bbox_to_anchor=(1.04, 1))

    # Plot the barplot
    plt.figure(figsize=(9, 6))
    sns.set(style="white")  # No grid lines

    # Create the plot using Seaborn
    plt.figure(figsize=(12, 6))
    sns.barplot(data=counts, x='Variable', y='Counts', hue='Environment', palette='deep')
    plt.title('Unique Number of miRNAs, Genes, MAGs, and Taxonomies', fontsize='16')
    plt.xlabel('Category', fontsize=14)
    plt.ylabel('Unique Count', fontsize=14)

    # Add count numbers on top of each bar
    for p in plt.gca().patches:
        height = p.get_height()
        plt.gca().annotate(f'{int(p.get_height())}', (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='center', fontsize=12, color='black', xytext=(0, 5),
                    textcoords='offset points')

    # Save the combined plot
    plt.savefig(f'{out_dir}/Grouped_Histograms.png', bbox_inches="tight")
    print("Grouped histograms generated and saved successfully.")