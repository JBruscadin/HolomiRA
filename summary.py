import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import warnings
import sys

##HolomiRA: summaries tables
input_file=sys.argv[1]
out_dir=sys.argv[2]

###Summary Tables:
#Two output files are expected for each environment provided by the user.
#One summarizes data by miRNA, and another summarizes data by MAG.
    
# Read the final file into a pandas DataFrame
df = pd.read_csv(input_file, sep='\t')

# Get the unique values in the 'Environment' column
unique_environments = df['Environment'].unique()

# Loop over each unique environment and perform the grouping and calculations
for env in unique_environments:
    # Create a subset DataFrame for the current environment
    subset_df = df[df['Environment'] == env]

    # Group by 'MAG' and calculate the required information for each MAG
    grouped_data_taxonomy = subset_df.groupby('MAG').agg({
        'Taxonomy': ['nunique', lambda x: ', '.join(x.dropna().unique())],
        'miRNA': ['nunique', lambda x: ', '.join(x.dropna().unique())],
        'Gene': ['nunique', lambda x: ', '.join(x.dropna().unique())]
    }).reset_index()

    # Flatten the multi-level column index
    grouped_data_taxonomy.columns = ['MAG', 'num_Taxonomy', 'Taxonomy', 'num_unique_miRNAs', 'unique_miRNAs', 'num_unique_genes', 'unique_genes']

    # Reorder the columns as needed
    grouped_data_taxonomy = grouped_data_taxonomy[['MAG', 'num_Taxonomy', 'Taxonomy', 'num_unique_miRNAs', 'num_unique_genes', 'unique_miRNAs', 'unique_genes']]

    # Save the result to a new file named "MAG_result_table_summary_taxonomy_<environment>.tsv"
    grouped_data_taxonomy.to_csv(f'{out_dir}/plots/MAG_result_table_summary_taxonomy_{env}.tsv', sep='\t', index=False)

    # Group by 'miRNA' and calculate the required information for each miRNA
    grouped_data_miRNA = subset_df.groupby('miRNA').agg({
        'Taxonomy': ['nunique', lambda x: ', '.join(x.dropna().unique())],
        'MAG': ['nunique', lambda x: ', '.join(x.dropna().unique())],
        'Gene': ['nunique', lambda x: ', '.join(x.dropna().unique())]
    }).reset_index()

    # Flatten the multi-level column index
    grouped_data_miRNA.columns = ['miRNA', 'num_unique_Taxa', 'unique_Taxa', 'num_unique_MAG', 'unique_MAG', 'num_unique_genes', 'unique_genes']

    # Reorder the columns as needed
    grouped_data_miRNA = grouped_data_miRNA[['miRNA', 'num_unique_Taxa','num_unique_MAG', 'num_unique_genes', 'unique_Taxa', 'unique_MAG', 'unique_genes']]

    # Save the result to a new file named "MAG_result_table_summary_miRNA_<environment>.tsv"
    grouped_data_miRNA.to_csv(f'{out_dir}/plots/MAG_result_table_summary_miRNA_{env}.tsv', sep='\t', index=False)

print("Summary tables generated and saved successfully.")



