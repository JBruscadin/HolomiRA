import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import sys
import os

# --- Input arguments ---
input_file = sys.argv[1]
out_dir = sys.argv[2]

# --- Output directories ---
final_results_dir = os.path.join(out_dir, 'final_results')
plots_dir = os.path.join(out_dir, 'plots')
os.makedirs(final_results_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

# --- Read the input file ---
df = pd.read_csv(input_file, sep='\t')

# --- Get unique environments ---
unique_environments = df['Environment'].unique()
num_environments = len(unique_environments)

# --- Color palette ---
palette = sns.color_palette("deep")


# --- Function to generate a blank figure with an informative message ---
def generate_blank_figure(message, output_path):
    plt.figure(figsize=(8, 6))
    plt.text(0.5, 0.5, message, fontsize=14, ha='center', va='center')
    plt.axis('off')
    plt.title("Venn Diagram Not Generated", fontsize=16)
    plt.savefig(output_path)
    plt.close()
    print(f"Blank figure saved with message: '{message}'")


# --- Function to create combined Venn diagrams ---
def create_combined_venn_diagram(df, unique_environments):
    mirna_sets = [set(df[df['Environment'] == env]['miRNA']) for env in unique_environments]
    gene_sets = [set(df[df['Environment'] == env]['Gene']) for env in unique_environments]
    taxonomy_sets = [set(df[df['Environment'] == env]['Taxonomy']) for env in unique_environments]

    with open(f'{final_results_dir}/venn_summary.txt', 'w') as f:
        fig, axs = plt.subplots(1, 3, figsize=(18, 6))

        if num_environments == 2:
            venn2([mirna_sets[0], mirna_sets[1]], set_labels=unique_environments, set_colors=[palette[0], palette[1]], ax=axs[0])
            axs[0].set_title('Unique miRNA', fontsize=14)
            f.write('### miRNA\n')
            f.write(f"Unique to {unique_environments[0]}: {mirna_sets[0] - mirna_sets[1]}\n")
            f.write(f"Unique to {unique_environments[1]}: {mirna_sets[1] - mirna_sets[0]}\n")
            f.write(f"Shared: {mirna_sets[0] & mirna_sets[1]}\n\n")

            venn2([gene_sets[0], gene_sets[1]], set_labels=unique_environments, set_colors=[palette[0], palette[1]], ax=axs[1])
            axs[1].set_title('Unique Gene', fontsize=14)
            f.write('### Gene\n')
            f.write(f"Unique to {unique_environments[0]}: {gene_sets[0] - gene_sets[1]}\n")
            f.write(f"Unique to {unique_environments[1]}: {gene_sets[1] - gene_sets[0]}\n")
            f.write(f"Shared: {gene_sets[0] & gene_sets[1]}\n\n")

            venn2([taxonomy_sets[0], taxonomy_sets[1]], set_labels=unique_environments, set_colors=[palette[0], palette[1]], ax=axs[2])
            axs[2].set_title('Unique Taxonomy', fontsize=14)
            f.write('### Taxonomy\n')
            f.write(f"Unique to {unique_environments[0]}: {taxonomy_sets[0] - taxonomy_sets[1]}\n")
            f.write(f"Unique to {unique_environments[1]}: {taxonomy_sets[1] - taxonomy_sets[0]}\n")
            f.write(f"Shared: {taxonomy_sets[0] & taxonomy_sets[1]}\n\n")

        elif num_environments == 3:
            venn3([mirna_sets[0], mirna_sets[1], mirna_sets[2]], set_labels=unique_environments, set_colors=[palette[0], palette[1], palette[2]], ax=axs[0])
            axs[0].set_title('Unique miRNA', fontsize=14)
            f.write('### miRNA\n')
            f.write(f"Unique to {unique_environments[0]}: {mirna_sets[0] - mirna_sets[1] - mirna_sets[2]}\n")
            f.write(f"Unique to {unique_environments[1]}: {mirna_sets[1] - mirna_sets[0] - mirna_sets[2]}\n")
            f.write(f"Unique to {unique_environments[2]}: {mirna_sets[2] - mirna_sets[0] - mirna_sets[1]}\n")
            f.write(f"Shared: {mirna_sets[0] & mirna_sets[1] & mirna_sets[2]}\n\n")

            venn3([gene_sets[0], gene_sets[1], gene_sets[2]], set_labels=unique_environments, set_colors=[palette[0], palette[1], palette[2]], ax=axs[1])
            axs[1].set_title('Unique Gene', fontsize=14)
            f.write('### Gene\n')
            f.write(f"Unique to {unique_environments[0]}: {gene_sets[0] - gene_sets[1] - gene_sets[2]}\n")
            f.write(f"Unique to {unique_environments[1]}: {gene_sets[1] - gene_sets[0] - gene_sets[2]}\n")
            f.write(f"Unique to {unique_environments[2]}: {gene_sets[2] - gene_sets[0] - gene_sets[1]}\n")
            f.write(f"Shared: {gene_sets[0] & gene_sets[1] & gene_sets[2]}\n\n")

            venn3([taxonomy_sets[0], taxonomy_sets[1], taxonomy_sets[2]], set_labels=unique_environments, set_colors=[palette[0], palette[1], palette[2]], ax=axs[2])
            axs[2].set_title('Unique Taxonomy', fontsize=14)
            f.write('### Taxonomy\n')
            f.write(f"Unique to {unique_environments[0]}: {taxonomy_sets[0] - taxonomy_sets[1] - taxonomy_sets[2]}\n")
            f.write(f"Unique to {unique_environments[1]}: {taxonomy_sets[1] - taxonomy_sets[0] - taxonomy_sets[2]}\n")
            f.write(f"Unique to {unique_environments[2]}: {taxonomy_sets[2] - taxonomy_sets[0] - taxonomy_sets[1]}\n")
            f.write(f"Shared: {taxonomy_sets[0] & taxonomy_sets[1] & taxonomy_sets[2]}\n\n")

        plt.subplots_adjust(top=0.85, wspace=0.3)
        plt.suptitle('Venn Diagrams for miRNA, Gene, and Taxonomy', fontsize=16)
        plt.savefig(f'{plots_dir}/Venn_diagram_combined.png')
        plt.close()


# --- Execute based on number of environments ---
if num_environments == 1:
    msg = "Only one environment detected.\nAt least 2 are required to generate a Venn diagram."
    generate_blank_figure(msg, f'{plots_dir}/Venn_diagram_combined.png')

elif num_environments > 3:
    msg = "More than three environments detected.\nVenn diagram generation supports up to 3."
    generate_blank_figure(msg, f'{plots_dir}/Venn_diagram_combined.png')

else:
    create_combined_venn_diagram(df, unique_environments)
    print(f"Venn diagrams and summary saved in {plots_dir}/ and {final_results_dir}/")
