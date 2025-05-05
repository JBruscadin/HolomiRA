#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# ---  Input arguments ---
input_file = sys.argv[1] 
out_dir = sys.argv[2]

# --- Read input data ---
df = pd.read_csv(input_file, sep='\t')

# --- Calculate unique counts per Environment ---
counts = df.groupby('Environment')[['miRNA', 'Gene', 'MAG', 'Taxonomy']].nunique().reset_index()

# --- Reshape for plotting ---
counts = pd.melt(counts, id_vars='Environment', var_name='Variable', value_name='Counts')

# --- Plot ---
plt.figure(figsize=(12, 6))
sns.set(style="white")

barplot = sns.barplot(data=counts, x='Variable', y='Counts', hue='Environment', palette='deep')
plt.title('Unique Number of miRNAs, Genes, MAGs, and Taxonomies', fontsize=16)
plt.xlabel('Category', fontsize=14)
plt.ylabel('Unique Count', fontsize=14)

# --- Add count annotations above bars (only if > 0) ---
for p in barplot.patches:
    height = p.get_height()
    if height > 0:
        barplot.annotate(
            f'{int(height)}',
            (p.get_x() + p.get_width() / 2., height),
            ha='center', va='bottom',
            fontsize=12, color='black',
            xytext=(0, 5), textcoords='offset points'
        )

plt.tight_layout()
plt.savefig(f"{out_dir}/plots/MAG_Histograms.png")
print("Histogram plot saved successfully.")