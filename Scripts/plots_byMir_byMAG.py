import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# Input file and output directory from command line arguments
input_file = sys.argv[1]
out_dir = sys.argv[2]

# Load the main data
df = pd.read_csv(input_file, sep='\t')

# Get unique environments
unique_environments = df['Environment'].unique()

def generate_plots(env, file_prefix, output_prefix):
    try:
        # Read the tables into pandas DataFrames
        table1 = pd.read_csv(f"{file_prefix}_{env}.tsv", sep="\t")
        table2 = pd.read_csv(f"{output_prefix}_{env}.tsv", sep="\t")
        
        # Prepare data for plotting, sorting values, and getting top 20
        df_taxa_top20 = table1[['miRNA', 'num_unique_MAG']].sort_values(by='num_unique_MAG', ascending=False).iloc[:20,]
        df_genes_top20 = table1[['miRNA', 'num_unique_genes']].sort_values(by='num_unique_genes', ascending=False).iloc[:20,]
        df_mir_top20 = table2[['MAG', 'num_unique_miRNAs']].sort_values(by='num_unique_miRNAs', ascending=False).iloc[:20,]
        df_mag_genes_top20 = table2[['MAG', 'num_unique_genes']].sort_values(by='num_unique_genes', ascending=False).iloc[:20,]

        # Create a 2x2 subplot layout
        fig, axs = plt.subplots(2, 2, figsize=(16, 12))  # 2x2 grid with larger size
        sns.set(style="white")  # No grid lines

        # Plot 1: Top 20 miRNAs by number of MAGs
        sns.barplot(x='num_unique_MAG', y='miRNA', data=df_taxa_top20, color='#55A868', ax=axs[0, 0])
        axs[0, 0].set_title('Top 20 miRNAs by number of MAGs', fontsize=14)
        axs[0, 0].set_xlabel('MAG Count', fontsize=12)
        axs[0, 0].set_ylabel('miRNA', fontsize=12)
        for spine in axs[0, 0].spines.values():
            spine.set_linewidth(2)

        # Plot 2: Top 20 miRNAs by number of target genes
        sns.barplot(x='num_unique_genes', y='miRNA', data=df_genes_top20, color='#DD8452', ax=axs[0, 1])
        axs[0, 1].set_title('Top 20 miRNAs by number of target genes', fontsize=14)
        axs[0, 1].set_xlabel('Gene Count', fontsize=12)
        axs[0, 1].set_ylabel('miRNA', fontsize=12)
        for spine in axs[0, 1].spines.values():
            spine.set_linewidth(2)

        # Plot 3: Top 20 MAGs by number of miRNAs
        sns.barplot(x='num_unique_miRNAs', y='MAG', data=df_mir_top20, color='#4C72B0', ax=axs[1, 0])
        axs[1, 0].set_title('Top 20 MAGs by number of miRNAs', fontsize=14)
        axs[1, 0].set_xlabel('miRNA Count', fontsize=12)
        axs[1, 0].set_ylabel('MAG', fontsize=12)
        for spine in axs[1, 0].spines.values():
            spine.set_linewidth(2)

        # Plot 4: Top 20 MAGs by number of target genes
        sns.barplot(x='num_unique_genes', y='MAG', data=df_mag_genes_top20, color='#DD8452', ax=axs[1, 1])
        axs[1, 1].set_title('Top 20 MAGs by number of target genes', fontsize=14)
        axs[1, 1].set_xlabel('Gene Count', fontsize=12)
        axs[1, 1].set_ylabel('MAG', fontsize=12)
        for spine in axs[1, 1].spines.values():
            spine.set_linewidth(2)

        # Adjust layout for more space between subplots
        plt.subplots_adjust(hspace=0.8, wspace=0.6)

        # Further adjust layout to avoid overlap
        plt.tight_layout()

        # Save the figure as PNG with 300 DPI
        plt.savefig(f"{out_dir}/plots/{env}_Top_20_miRNAs_and_MAGs.png", dpi=300)
        plt.close()

        print(f"Top 20 plots for {env} generated and saved successfully.")

    except Exception as e:
        print(f"Error while generating plots for {env}: {e}")

# Loop over each unique environment and generate the plots
for env in unique_environments:
    generate_plots(env, f'{out_dir}/final_results/MAG_result_table_summary_miRNA', f'{out_dir}/final_results/MAG_result_table_summary_taxonomy')
