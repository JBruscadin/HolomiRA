import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import warnings
import sys

##HolomiRA:  Venn diagram
input_file=sys.argv[1]
out_dir=sys.argv[2]

# Read the all_results.tsv file into a pandas DataFrame
df = pd.read_csv(input_file, sep='\t')

# Get the unique values in the 'Environment' column
unique_environments = df['Environment'].unique()

def create_venn_diagrams(df, unique_environments):
    # Create sets for unique miRNA, Gene, and MAG for each environment
    mirna_sets = [set(df[df['Environment'] == env]['miRNA']) for env in unique_environments]
    gene_sets = [set(df[df['Environment'] == env]['Gene']) for env in unique_environments]
    taxonomy_sets = [set(df[df['Environment'] == env]['Taxonomy']) for env in unique_environments]

    # Create Venn diagrams based on the number of environments
    num_environments = len(unique_environments)
    if num_environments == 2:
        plt.figure(figsize=(10, 6))
        venn2([mirna_sets[0], mirna_sets[1]], set_labels=unique_environments)
        plt.title('Unique miRNA')
        plt.savefig(f'{out_dir}/plots/Venn_diagram_miRNA.png')

        plt.figure(figsize=(10, 6))
        venn2([gene_sets[0], gene_sets[1]], set_labels=unique_environments)
        plt.title('Unique Gene')
        plt.savefig(f'{out_dir}/plots/Venn_diagram_gene.png')

        plt.figure(figsize=(10, 6))
        venn2([taxonomy_sets[0], taxonomy_sets[1]], set_labels=unique_environments)
        plt.title('Unique Taxonomy')
        plt.savefig(f'{out_dir}/plots/Venn_diagram_taxonomy.png')

    elif num_environments == 3:
        plt.figure(figsize=(10, 6))
        venn3([mirna_sets[0], mirna_sets[1], mirna_sets[2]], set_labels=unique_environments)
        plt.title('Unique miRNA')
        plt.savefig(f'{out_dir}/plots/Venn_diagram_miRNA.png')

        plt.figure(figsize=(10, 6))
        venn3([gene_sets[0], gene_sets[1], gene_sets[2]], set_labels=unique_environments)
        plt.title('Unique Gene')
        plt.savefig(f'{out_dir}/plots/Venn_diagram_gene.png')

        plt.figure(figsize=(10, 6))
        venn3([taxonomy_sets[0], taxonomy_sets[1], taxonomy_sets[2]], set_labels=unique_environments)
        plt.title('Unique Taxonomy')
        plt.savefig(f'{out_dir}/plots/Venn_diagram_taxonomy.png')

# Check the number of unique environments
num_environments = len(unique_environments)

# Check the number of environments and call the function accordingly
if num_environments == 1:
    print("Error: Only one environment found.")
elif num_environments > 3:
    print("Error: More than three environments found.")
else:
    create_venn_diagrams(df, unique_environments)
    print("Venn diagrams generated and saved successfully.")

