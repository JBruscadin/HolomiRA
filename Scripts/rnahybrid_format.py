#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys

MAG_ID = sys.argv[1]
out_dir=sys.argv[2]
input_file=sys.argv[3]

a=pd.read_csv(input_file, sep=':', names=['seq','position','nc','mir','ncmir','mfe','p','1','2','3','4','5'])
#a[['startend','strand']]=a.position.str.split("\(", expand=True)
a[['startend', 'strand']] = a['position'].str.split("\(", expand=True)
#a[['start','end']]=a.startend.str.split("-", expand=True)
a[['start', 'end']] = a['startend'].str.split("-", expand=True)
id_df=pd.read_csv(f"{out_dir}/annotation/{MAG_ID}/{MAG_ID}_IDs.txt", sep="\t", header=None, names=["seq","1","2", "ID", "3"])
##a_id=pd.merge(a, id_df[["seq","ID"]], how="inner", on="seq")
a_id = pd.merge(a, id_df[["seq", "ID"]], how="inner", on="seq")

a_id.loc[:, 'ID'] = a_id['ID'].str.replace(r"ID=", "")
##b=a_id[["seq","start","end","mir","mfe","p", "ID"]]
#b.loc[:, 'ID'] = b['ID'].str.replace(r"ID=", "")
#b['ID']=b['ID'].str.replace(r"ID=","")
tsv=pd.read_csv(f"{out_dir}/annotation/{MAG_ID}/{MAG_ID}.tsv", sep="\t")
tsv = tsv[tsv['ftype'] == 'gene']
tsv['gene'] = tsv['gene'].str.split('_').str[0] #removing duplicates
tsv=tsv.drop_duplicates()

#tsv_b=pd.merge(b, tsv, how="inner", left_on="ID", right_on="locus_tag")
tsv_b=pd.merge(a_id, tsv, how="inner", left_on="ID", right_on="locus_tag")
tsv_b["sample"]=MAG_ID
bsites=tsv_b[["sample", "seq", "start", "end", "mir", "ID", "mfe", "p", "gene"]]
bsites.to_csv(f"{out_dir}/rnahybrid/{MAG_ID}_bsites.tsv", sep="\t", index=False)

