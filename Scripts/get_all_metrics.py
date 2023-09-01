#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
import os
from os import system

out_dir=sys.argv[1]
list0=sys.argv[2]
file=sys.argv[3]

list=pd.read_csv(list0, sep="\t")


final = pd.DataFrame([])
print('------- Start obtaining metrics for all results ----------------')

for index,row in list.iterrows():
    MAG_ID=row["SampleID"]
    path=f'{out_dir}/rnahybrid/{MAG_ID}_finalresults.tsv'
    indiv=pd.read_csv(path, sep="\t")
    final=pd.concat([indiv, final])
final.drop_duplicates().to_csv(f'{out_dir}/final_results/HolomiRA_results.tsv', sep="\t", index=None)


