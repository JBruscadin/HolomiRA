#!/usr/bin/env python
import pybedtools
import sys
from pybedtools import BedTool
from pybedtools import featurefuncs
import os
import pysam

#include bash variables
sample = sys.argv[1]
fasta_dir = sys.argv[2]
upstream = int(sys.argv[3])
downstream = int(sys.argv[4])
out_dir=sys.argv[5]
gff=sys.argv[6]


a=BedTool(gff)
fasta=f"{fasta_dir}/{sample}.fa"


b=a.each(pybedtools.featurefuncs.five_prime, upstream, downstream, add_to_name=None, genome=None).saveas(f"{out_dir}/annotation/{sample}/{sample}_cds_fiveprime.gff")

input_gff = f"{out_dir}/annotation/{sample}/{sample}_cds_fiveprime.gff"
fasta = f"{out_dir}/annotation/{sample}/{sample}.fna"
output_fasta = f"{out_dir}/target_fasta/{sample}_filtered.fa"
warnings_file = "warnings.log"

command = f"bedtools getfasta -s -fo {output_fasta} -fi {fasta} -bed {input_gff} 2> {warnings_file}"
os.system(command)

input_gff2=gff
output_fasta = f"{out_dir}/target_fasta/{sample}_CDS.fa"
warnings_file = "warnings2.log"

command2 = f"bedtools getfasta -s -fo {output_fasta} -fi {fasta} -bed {input_gff2} 2> {warnings_file}"
os.system(command2)
