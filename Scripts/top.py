import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import warnings
import sys

##HolomiRA: Top 10

# Read the all_results.tsv file into a pandas DataFrame

input_file=sys.argv[1]
out_dir=sys.argv[2]

df = pd.read_csv(input_file, sep='\t')

# Get the unique values in the 'Environment' column
unique_environments = df['Environment'].unique()

def generate_plots(env, file_prefix, output_prefix):
    try:
        # Read the tables into pandas DataFrames
        table1 = pd.read_csv(f"{file_prefix}_{env}.tsv", sep="\t")
        table2 = pd.read_csv(f"{output_prefix}_{env}.tsv", sep="\t")
        
        # Generate plots for table 1
        fig, axes = plt.subplots(2, 1, figsize=(8, 10))

       # Plot top 10 mirnas based on num_unique_MAG
        table1_sorted_taxa = table1.sort_values("num_unique_MAG")
        top10_taxa = pd.concat([table1_sorted_taxa.head(10), table1_sorted_taxa.tail(10)])
        axes[0].barh(top10_taxa["miRNA"], top10_taxa["num_unique_MAG"])
        axes[0].set_xlabel("Number of MAG")
        axes[0].set_ylabel("miRNA")
        axes[0].set_title("Top 10 miRNAs by MAG")

        # Plot top 10 mirnas based on num_unique_genes
        table1_sorted_genes = table1.sort_values("num_unique_genes")
        top10_genes = pd.concat([table1_sorted_genes.head(10), table1_sorted_genes.tail(10)])
        axes[1].barh(top10_genes["miRNA"], top10_genes["num_unique_genes"])
        axes[1].set_xlabel("Number of Unique Genes")
        axes[1].set_ylabel("miRNA")
        axes[1].set_title("Top 10 miRNAs by Unique Genes")

        # Adjust the spacing between subplots
        plt.tight_layout(pad=1.0) 

        # Save the plot as png with 300 PPI
        plt.savefig(f"{out_dir}/plots/{env}_Top_10_miRNAs.png", dpi=300)
        plt.close()

        # Generate plots for table 2
        fig, axes = plt.subplots(2, 1, figsize=(8, 10))

       # Plot top 10 MAG based on num_unique_mir
        table2_sorted_mir = table2.sort_values("num_unique_miRNAs")
        top10_mir = pd.concat([table2_sorted_mir.head(10), table2_sorted_mir.tail(10)])
        axes[0].barh(top10_mir["MAG"], top10_mir["num_unique_miRNAs"])
        axes[0].set_xlabel("Number of Unique miRNAs")
        axes[0].set_ylabel("MAG")
        axes[0].set_title("Top 10 MAG by Unique miRNAs")

        # Plot top 10 MAG based on num_unique_genes
        table2_sorted_genes = table2.sort_values("num_unique_genes")
        top10_genes = pd.concat([table2_sorted_genes.head(10), table2_sorted_genes.tail(10)])
        axes[1].barh(top10_genes["MAG"], top10_genes["num_unique_genes"])
        axes[1].set_xlabel("Number of Unique Genes")
        axes[1].set_ylabel("MAG")
        axes[1].set_title("Top 10 MAG by Unique Genes")

        # Adjust the spacing between subplots
        plt.tight_layout(pad=1.0) 

        # Save the plot as JPG with 300 PPI
        plt.savefig(f"{out_dir}/plots/{env}_Top_10_MAG.png", dpi=300)
        plt.close()

        print(f"Plots Top 10 for {env} generated and saved successfully.")
    except Exception as e:
        # Log the error
        print(f"Error: {e}")

# Loop over each unique environment and perform the grouping and calculations
for env in unique_environments:
    # Call the function to generate plots for the current environment
    generate_plots(env, f'{out_dir}/plots/MAG_result_table_summary_miRNA', f'{out_dir}/plots/MAG_result_table_summary_taxonomy')



