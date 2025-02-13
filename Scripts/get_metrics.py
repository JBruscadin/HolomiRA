#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
import os
from os import system

out_dir=sys.argv[1]
id=sys.argv[2]
file=sys.argv[3]
list0=sys.argv[4]

list=pd.read_csv(list0, sep="\t")

for index,row in list.iterrows():
    MAG_ID=row["SampleID"]
    #print(index)

    #print(f"Creating {MAG_ID} final results file")
    a=pd.read_csv(f'{out_dir}/rnahybrid/{MAG_ID}_bsites.tsv', sep="\t").drop_duplicates()
    taxon=pd.read_csv(id, sep="\t", names=["sample", "Taxonomy", "Environment"])
    a_taxon=pd.merge(a, taxon, how="inner", on="sample").drop_duplicates()
    names={ 'sample' : 'MAG',
            'seq' : 'Contig',
            'start': 'Start',
            'end': 'End',
            'mir': 'miRNA',
            'ID' : 'Locus_tag',
            'mfe': 'MFE',
            'p': 'Pvalue',
            'gene': 'Gene',
            }

    a_taxon.rename(columns=names,inplace=True)
    a_taxon.to_csv(f"{out_dir}/rnahybrid/{MAG_ID}_finalresults.tsv", sep="\t", index=None)
