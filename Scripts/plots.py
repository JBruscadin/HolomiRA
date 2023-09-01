#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
import os
from os import system
import seaborn as sns
import matplotlib.pyplot as plt
#from matplotlib_venn import venn2, venn3
import warnings

out_dir=sys.argv[1]
list0=sys.argv[2]
file=sys.argv[3]

list=pd.read_csv(list0, sep="\t")

#### plots ####
##HolomiRA: summaries plots and tables

## Summaries tables:
# Two output files are expected for each environment provided by the user.
# One summarizy by miRNA and another by Taxonomy

# Read the all_results.tsv file into a pandas DataFrame
df = pd.read_csv(file, sep='\t')
df=df.drop_duplicates()

unique_environments = df['Environment'].unique()   #user defined
#unique_environments=environment

# Loop over each unique environment and perform the grouping and calculations
for env in unique_environments:
    # Create a subset DataFrame for the current environment
    subset_df = df[df['Environment'] == env]

    # Group by 'MAG' and calculate the required information for each SamPvaluele_id
    grouped_data_taxonomy = subset_df.groupby('MAG').agg({
        'miRNA': ['nunique', lambda x: ', '.join(x.dropna().unique())],
        'Gene': ['nunique', lambda x: ', '.join(x.dropna().unique())]
    }).reset_index()

    # Flatten the multi-level column index
    grouped_data_taxonomy.columns = ['MAG', 'num_unique_miRNAs', 'unique_miRNAs', 'num_unique_genes', 'unique_genes']

    # Reorder the columns as needed
    grouped_data_taxonomy = grouped_data_taxonomy[['MAG', 'num_unique_miRNAs', 'num_unique_genes', 'unique_miRNAs', 'unique_genes']]

    # Save the result to a new file named "result_table_summary_taxonomy_<environment>.tsv"
    grouped_data_taxonomy.to_csv(f'MAG_result_table_summary_taxonomy_{env}.tsv', sep='\t', index=False)

    # Group by 'miRNA' and calculate the required information for each miRNA
    grouped_data_miRNA = subset_df.groupby('miRNA').agg({
        'MAG': ['nunique', lambda x: ', '.join(x.dropna().unique())],
        'Gene': ['nunique', lambda x: ', '.join(x.dropna().unique())]
    }).reset_index()

    # Flatten the multi-level column index
    grouped_data_miRNA.columns = ['miRNA', 'num_unique_taxa', 'unique_taxa', 'num_unique_genes', 'unique_genes']

    # Reorder the columns as needed
    grouped_data_miRNA = grouped_data_miRNA[['miRNA', 'num_unique_taxa', 'num_unique_genes', 'unique_taxa', 'unique_genes']]

    # Save the result to a new file named "result_table_summary_miRNA_<environment>.tsv"
    grouped_data_miRNA.to_csv(f'MAG_result_table_summary_miRNA_{env}.tsv', sep='\t', index=False)

print("Summary tables generated and saved successfully.")

## PLOTs

####plot 1
# Calculate the total number of unique miRNA, gene, and MAG for each Environment
counts = df.groupby(['Environment'])[['miRNA', 'Gene', 'MAG']].nunique()

# Reset index to use 'Environment' as column
counts.reset_index(inplace=True)

# Unpivot the data to use 'Variable' for the types (miRNA, Gene, MAG)
counts = pd.melt(counts, id_vars='Environment', var_name='Variable', value_name='Counts')

# Create the plot using Seaborn
plt.figure(figsize=(12, 6))
sns.barplot(data=counts, x='Variable', y='Counts', hue='Environment')
plt.title('')
plt.xlabel('Variable')
plt.ylabel('Counts')

# Add count numbers on top of each bar
for p in plt.gca().patches:
    height = p.get_height()
    plt.gca().annotate(f'{int(height)}', (p.get_x() + p.get_width() / 2., height),
                       ha='center', va='bottom')

# Save the plot
plt.savefig('MAG_Histograms.png')
print("Histograms generated and saved successfully.")

###### veen diagram

###### veen diagram
def create_venn_diagrams(df, unique_environments):
    # Create sets for unique miRNA, Gene, and MAG for each environment
    mirna_sets = [set(df[df['Environment'] == env]['miRNA']) for env in unique_environments]
    gene_sets = [set(df[df['Environment'] == env]['Gene']) for env in unique_environments]
    taxonomy_sets = [set(df[df['Environment'] == env]['MAG']) for env in unique_environments]

    # Create Venn diagrams based on the number of environments
    num_environments = len(unique_environments)
    if num_environments == 2:
        plt.figure(figsize=(10, 6))
        venn2([mirna_sets[0], mirna_sets[1]], set_labels=unique_environments)
        plt.title('Unique miRNA')
        plt.savefig('MAG_venn_diagram_miRNA.png')

        plt.figure(figsize=(10, 6))
        venn2([gene_sets[0], gene_sets[1]], set_labels=unique_environments)
        plt.title('Unique Gene')
        plt.savefig('MAG_venn_diagram_gene.png')

        plt.figure(figsize=(10, 6))
        venn2([taxonomy_sets[0], taxonomy_sets[1]], set_labels=unique_environments)
        plt.title('Unique MAG')
        plt.savefig('venn_diagram_taxonomy.png')

    elif num_environments == 3:
        plt.figure(figsize=(10, 6))
        venn3([mirna_sets[0], mirna_sets[1], mirna_sets[2]], set_labels=unique_environments)
        plt.title('Unique miRNA')
        plt.savefig('MAG_venn_diagram_miRNA.png')

        plt.figure(figsize=(10, 6))
        venn3([gene_sets[0], gene_sets[1], gene_sets[2]], set_labels=unique_environments)
        plt.title('Unique Gene')
        plt.savefig('MAG_venn_diagram_gene.png')

        plt.figure(figsize=(10, 6))
        venn3([taxonomy_sets[0], taxonomy_sets[1], taxonomy_sets[2]], set_labels=unique_environments)
        plt.title('Unique MAG')
        plt.savefig('MAG_venn_diagram_taxonomy.png')

# Check the number of unique environments
num_environments = len(unique_environments)

# Check the number of environments and call the function accordingly
if num_environments == 1:
    print("Error: Only one environment found.")
elif num_environments > 3:
    print("Error: More than three environments found.")
else:
    create_venn_diagrams(df, unique_environments)
    print("Veen diagrams generated and saved successfully.")

### plot top 10

def generate_plots(env, file_prefix, output_prefix):
    try:
        # Read the tables into pandas DataFrames
        table1 = pd.read_csv(f"{file_prefix}_{env}.tsv", sep="\t")
        table2 = pd.read_csv(f"{output_prefix}_{env}.tsv", sep="\t")

        # Generate plots for table 1
        fig, axes = plt.subplots(2, 1, figsize=(8, 10))

       # Plot top 10 mirnas based on num_unique_taxa
        table1_sorted_taxa = table1.sort_values("num_unique_taxa")
        top10_taxa = pd.concat([table1_sorted_taxa.head(10), table1_sorted_taxa.tail(10)])
        axes[0].barh(top10_taxa["miRNA"], top10_taxa["num_unique_taxa"])
        axes[0].set_xlabel("Number of Unique Taxa")
        axes[0].set_ylabel("miRNA")
        axes[0].set_title("Top 10 miRNAs by Unique Taxa")

        # Plot top 10 mirnas based on num_unique_genes
        table1_sorted_genes = table1.sort_values("num_unique_genes")
        top10_genes = pd.concat([table1_sorted_genes.head(10), table1_sorted_genes.tail(10)])
        axes[1].barh(top10_genes["miRNA"], top10_genes["num_unique_genes"])
        axes[1].set_xlabel("Number of Unique Genes")
        axes[1].set_ylabel("miRNA")
        axes[1].set_title("Top 10 miRNAs by Unique Genes")

        # Adjust the spacing between subplots
        plt.tight_layout()


        # Adjust the spacing between subplots
        plt.tight_layout()

        # Save the plot as png with 300 PPI
        plt.savefig(f"{env}_Top_10_miRNAs_MAG.png", dpi=300)
        plt.close()

        # Generate plots for table 2
        fig, axes = plt.subplots(2, 1, figsize=(8, 10))

       # Plot top 10 taxonomy based on num_unique_mir
      # Plot top 10 taxonomy based on num_unique_mir
        table2_sorted_mir = table2.sort_values("num_unique_miRNAs")
        top10_mir = pd.concat([table2_sorted_mir.head(10), table2_sorted_mir.tail(10)])
        axes[0].barh(top10_mir["MAG"], top10_mir["num_unique_miRNAs"])
        axes[0].set_xlabel("Number of Unique miRNAs")
        axes[0].set_ylabel("MAG")
        axes[0].set_title("Top 10 MAG by Unique miRNAs")

        # Plot top 10 taxonomy based on num_unique_genes
        table2_sorted_genes = table2.sort_values("num_unique_genes")
        top10_genes = pd.concat([table2_sorted_genes.head(10), table2_sorted_genes.tail(10)])
        axes[1].barh(top10_genes["MAG"], top10_genes["num_unique_genes"])
        axes[1].set_xlabel("Number of Unique Genes")
        axes[1].set_ylabel("MAG")
        axes[1].set_title("Top 10 MAG by Unique Genes")

        # Adjust the spacing between subplots
        plt.tight_layout()

        # Save the plot as JPG with 300 PPI
        plt.savefig(f"{env}_Top_10_taxonomy_MAG.png", dpi=300)
        plt.close()

        print(f"Plots Top 10 for {env} generated and saved successfully.")
    except Exception as e:
        # Log the error
        print(f"Error: {e}")

# Loop over each unique environment and perform the grouping and calculations
for env in unique_environments:
    # Call the function to generate plots for the current environment
    generate_plots(env, 'MAG_result_table_summary_miRNA', 'MAG_result_table_summary_taxonomy')

