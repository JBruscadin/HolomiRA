import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import warnings
import sys
##HolomiRA: Histograms

# Read the all_results.tsv file into a pandas DataFrame
#df = pd.read_csv('HolomiRA_results.tsv', sep='\t')

input_file=sys.argv[1]
out_dir=sys.argv[2]

df = pd.read_csv(input_file, sep='\t')

# Calculate the total number of unique miRNA, gene, MAG and Taxonomy for each Environment
counts = df.groupby(['Environment'])[['miRNA', 'Gene', 'MAG', 'Taxonomy']].nunique()

# Reset index to use 'Environment' as column
counts.reset_index(inplace=True)

# Unpivot the data to use 'Variable' for the types (miRNA, Gene, MAG, Taxonomy)
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
plt.savefig(f'{out_dir}/plots/MAG_Histograms.png')
print("Histograms generated and saved successfully.")
