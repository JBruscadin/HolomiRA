import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import sys
import os

## HolomiRA: Venn diagram
input_file = sys.argv[1]
out_dir = sys.argv[2]

# --- Read the all_results.tsv file into a pandas DataFrame ---
df = pd.read_csv(input_file, sep='\t')

# --- Get the unique values in the 'Environment' column ---
unique_environments = df['Environment'].unique()

# --- Create the output directory if it doesn't exist ---
final_results_dir = os.path.join(out_dir, 'final_results')
if not os.path.exists(final_results_dir):
    os.makedirs(final_results_dir)

# --- Extract the "deep" color palette from Seaborn ---
palette = sns.color_palette("deep")

# --- Function to create combined Venn diagrams with aligned titles ---
def create_combined_venn_diagram(df, unique_environments):
    # Create sets for unique miRNA, Gene, and Taxonomy for each environment
    mirna_sets = [set(df[df['Environment'] == env]['miRNA']) for env in unique_environments]
    gene_sets = [set(df[df['Environment'] == env]['Gene']) for env in unique_environments]
    taxonomy_sets = [set(df[df['Environment'] == env]['Taxonomy']) for env in unique_environments]

    # Open a file to save the summary of unique and shared elements
    with open(f'{final_results_dir}/venn_summary.txt', 'w') as f:
        # Create a single figure with 3 subplots
        fig, axs = plt.subplots(1, 3, figsize=(18, 6))  # Set the figure size for the combined plot
        
        # Create Venn diagrams based on the number of environments
        num_environments = len(unique_environments)
        if num_environments == 2:
            # Venn for miRNA
            venn = venn2([mirna_sets[0], mirna_sets[1]], set_labels=unique_environments, set_colors=[palette[0], palette[1]], ax=axs[0])
            axs[0].set_title('Unique miRNA', fontsize=14)

            # Write miRNA unique/shared sets to file
            f.write('### miRNA\n')
            f.write(f"Unique to {unique_environments[0]}: {mirna_sets[0] - mirna_sets[1]}\n")
            f.write(f"Unique to {unique_environments[1]}: {mirna_sets[1] - mirna_sets[0]}\n")
            f.write(f"Shared: {mirna_sets[0] & mirna_sets[1]}\n\n")

            # Venn for Gene
            venn = venn2([gene_sets[0], gene_sets[1]], set_labels=unique_environments, set_colors=[palette[0], palette[1]], ax=axs[1])
            axs[1].set_title('Unique Gene', fontsize=14)

            # Write Gene unique/shared sets to file
            f.write('### Gene\n')
            f.write(f"Unique to {unique_environments[0]}: {gene_sets[0] - gene_sets[1]}\n")
            f.write(f"Unique to {unique_environments[1]}: {gene_sets[1] - gene_sets[0]}\n")
            f.write(f"Shared: {gene_sets[0] & gene_sets[1]}\n\n")

            # Venn for Taxonomy
            venn = venn2([taxonomy_sets[0], taxonomy_sets[1]], set_labels=unique_environments, set_colors=[palette[0], palette[1]], ax=axs[2])
            axs[2].set_title('Unique Taxonomy', fontsize=14)

            # Write Taxonomy unique/shared sets to file
            f.write('### Taxonomy\n')
            f.write(f"Unique to {unique_environments[0]}: {taxonomy_sets[0] - taxonomy_sets[1]}\n")
            f.write(f"Unique to {unique_environments[1]}: {taxonomy_sets[1] - taxonomy_sets[0]}\n")
            f.write(f"Shared: {taxonomy_sets[0] & taxonomy_sets[1]}\n\n")

        elif num_environments == 3:
            # Venn for miRNA
            venn = venn3([mirna_sets[0], mirna_sets[1], mirna_sets[2]], set_labels=unique_environments, set_colors=[palette[0], palette[1], palette[2]], ax=axs[0])
            axs[0].set_title('Unique miRNA', fontsize=14)

            # Write miRNA unique/shared sets to file
            f.write('### miRNA\n')
            f.write(f"Unique to {unique_environments[0]}: {mirna_sets[0] - mirna_sets[1] - mirna_sets[2]}\n")
            f.write(f"Unique to {unique_environments[1]}: {mirna_sets[1] - mirna_sets[0] - mirna_sets[2]}\n")
            f.write(f"Unique to {unique_environments[2]}: {mirna_sets[2] - mirna_sets[0] - mirna_sets[1]}\n")
            f.write(f"Shared: {mirna_sets[0] & mirna_sets[1] & mirna_sets[2]}\n\n")

            # Venn for Gene
            venn = venn3([gene_sets[0], gene_sets[1], gene_sets[2]], set_labels=unique_environments, set_colors=[palette[0], palette[1], palette[2]], ax=axs[1])
            axs[1].set_title('Unique Gene', fontsize=14)

            # Write Gene unique/shared sets to file
            f.write('### Gene\n')
            f.write(f"Unique to {unique_environments[0]}: {gene_sets[0] - gene_sets[1] - gene_sets[2]}\n")
            f.write(f"Unique to {unique_environments[1]}: {gene_sets[1] - gene_sets[0] - gene_sets[2]}\n")
            f.write(f"Unique to {unique_environments[2]}: {gene_sets[2] - gene_sets[0] - gene_sets[1]}\n")
            f.write(f"Shared: {gene_sets[0] & gene_sets[1] & gene_sets[2]}\n\n")

            # Venn for Taxonomy
            venn = venn3([taxonomy_sets[0], taxonomy_sets[1], taxonomy_sets[2]], set_labels=unique_environments, set_colors=[palette[0], palette[1], palette[2]], ax=axs[2])
            axs[2].set_title('Unique Taxonomy', fontsize=14)

            # Write Taxonomy unique/shared sets to file
            f.write('### Taxonomy\n')
            f.write(f"Unique to {unique_environments[0]}: {taxonomy_sets[0] - taxonomy_sets[1] - taxonomy_sets[2]}\n")
            f.write(f"Unique to {unique_environments[1]}: {taxonomy_sets[1] - taxonomy_sets[0] - taxonomy_sets[2]}\n")
            f.write(f"Unique to {unique_environments[2]}: {taxonomy_sets[2] - taxonomy_sets[0] - taxonomy_sets[1]}\n")
            f.write(f"Shared: {taxonomy_sets[0] & taxonomy_sets[1] & taxonomy_sets[2]}\n\n")

    # Adjust the layout to ensure titles are aligned
    plt.subplots_adjust(top=0.85, wspace=0.3)

    # Set a super title for the whole figure (optional)
    plt.suptitle('Venn Diagrams for miRNA, Gene, and Taxonomy', fontsize=16)

    # Save the combined figure with the three Venn diagrams
    plt.savefig(f'{out_dir}/plots/Venn_diagram_combined.png')
    plt.close()  # Close the plot to avoid display issues

# --- Check the number of unique environments
num_environments = len(unique_environments)

if num_environments == 1:
    print("Error: Only one environment found.")
elif num_environments > 3:
    print("Error: More than three environments found.")
else:
    create_combined_venn_diagram(df, unique_environments)
    print(f"Combined Venn diagrams generated and saved successfully in {out_dir}/plots/")
    print(f"Venn summary saved successfully in {final_results_dir}/venn_summary.txt.")