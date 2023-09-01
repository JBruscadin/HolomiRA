import pandas as pd
import os

configfile: "config.yaml"
FASTA_DIR=config["fasta_dir"]
OUT_DIR=config["out_dir"]
UPS=config["upstream"]
DWNS=config["downstream"]
REF_MIR=config["ref_mir"]
SEED=config["seed"]
ENERGY=config["energy"]
PVALUE=config["pvalue"]
ID=config["id"]
ENV=config["environment"]
sample_tab=pd.read_csv(config["sample_tab"], header=0, sep = "\t")
sample=sample_tab["SampleID"].drop_duplicates().to_list()


rule all:
	input:  expand("{out_dir}/rnahybrid/{sample}_finalresults.tsv", out_dir=OUT_DIR, sample=sample),
		OUT_DIR+"/final_results/HolomiRA_results.tsv",
		expand("{out_dir}/plots/MAG_result_table_summary_miRNA_{env}.tsv", out_dir=OUT_DIR, env=ENV),
		OUT_DIR+"/plots/MAG_Histograms.png",
		OUT_DIR+"/plots/Venn_diagram_taxonomy.png",
		expand("{out_dir}/plots/{env}_Top_10_MAG.png", out_dir=OUT_DIR, env=ENV),
		
rule annotate_prokka:
	input: FASTA_DIR+"{sample}.fa"
	output:
		OUT_DIR+"/annotation/{sample}/{sample}.gff"

	params: out_dir=OUT_DIR
	conda: "Envs/prokka.yml"
	shell: 
		"""prokka --quiet --outdir {params.out_dir}/annotation/{wildcards.sample} --prefix {wildcards.sample} --addgenes --cpus 12 {input} --force """


rule cds:
	input: 
		OUT_DIR+"/annotation/{sample}/{sample}.gff"
	output: 
		OUT_DIR+"/annotation/{sample}/{sample}_cds.gff"
	shell: """ awk '$3 == "CDS" {{print $0}}' {input} > {output} """
	

rule id:
	input:	OUT_DIR+"/annotation/{sample}/{sample}_cds.gff"
	output:	OUT_DIR+"/annotation/{sample}/{sample}_IDs.txt"
	shell:	""" grep -v '^#' {input} | awk '{{print $1, $4, $5,$7,$9}}' | awk '{{split($5,a,/;/); print $1,$2,$3,a[1],$4}}' OFS="\t"> {output} """



rule five_prime:
	input: OUT_DIR+"/annotation/{sample}/{sample}_cds.gff"
	params:
		fasta=FASTA_DIR,
		upstream=UPS,
		downstream=DWNS,
		out_dir=OUT_DIR
	output: OUT_DIR+"/target_fasta/{sample}_filtered.fa",
		OUT_DIR+"/target_fasta/{sample}_CDS.fa"
	
	conda: "Envs/fiveprime.yml"
	shell: """ python get_fiveprime_fa.py {wildcards.sample} {params.fasta} {params.upstream} {params.downstream} {params.out_dir} {input}"""


rule find_targets:
	input: 
		fasta=OUT_DIR+"/target_fasta/{sample}_filtered.fa",
		ref_mir=REF_MIR
	output: 
		OUT_DIR+"/rnahybrid/{sample}putative_targets.tsv"
	params: 
		seed='-f '+SEED if SEED!='NA' else [],
		e="-e "+str(ENERGY),
		p="-p "+str(PVALUE)
	conda: "Envs/rnahybrid.yml"
	shell: """ RNAhybrid -s 3utr_human -c -t {input.fasta} -q {input.ref_mir} {params.seed} {params.e} {params.p} > {output} """


rule format_rnahybrid:
	input: OUT_DIR+"/rnahybrid/{sample}putative_targets.tsv"
	conda: "Envs/formatOutputs.yml"
	params: out_dir=OUT_DIR
	output: OUT_DIR+"/rnahybrid/{sample}_bsites.tsv"
	shell: """ python rnahybrid_format.py {wildcards.sample} {params.out_dir} {input}"""


rule get_indiv_metrics:
	input: file=OUT_DIR+"/rnahybrid/{sample}_bsites.tsv",
		list=config["sample_tab"]
	output: OUT_DIR+"/rnahybrid/{sample}_finalresults.tsv"
	conda: "Envs/formatOutputs.yml"
	params: out_dir=OUT_DIR,
		id=ID
	shell: """   python get_metrics.py {params.out_dir} {params.id} {input.file} {input.list} """

rule get_global_metrics:
	input: files=expand("{out_dir}/rnahybrid/{sample}_finalresults.tsv", out_dir=OUT_DIR, sample=sample),
		list=config["sample_tab"],
	output: OUT_DIR+"/final_results/HolomiRA_results.tsv"
        params: out_dir=OUT_DIR
	conda: "Envs/formatOutputs.yml"
	shell: """   python get_all_metrics.py {params.out_dir} {input.list} {input.files} """


rule summary:
	input: OUT_DIR+"/final_results/HolomiRA_results.tsv"
	output:expand("{out_dir}/plots/MAG_result_table_summary_miRNA_{env}.tsv", out_dir=OUT_DIR, env=ENV)
	conda: "Envs/plots.yml"
	params: out_dir=OUT_DIR
	shell: """ python summary.py {params.out_dir}/final_results/HolomiRA_results.tsv {params.out_dir} """

rule plt_histogram:
        input: OUT_DIR+"/final_results/HolomiRA_results.tsv"
        output: OUT_DIR+"/plots/MAG_Histograms.png"
        conda: "Envs/plots.yml"
        params: out_dir=OUT_DIR
        shell: """ python histogram.py {params.out_dir}/final_results/HolomiRA_results.tsv {params.out_dir} """

rule plt_venn:
        input: OUT_DIR+"/final_results/HolomiRA_results.tsv"
        output: OUT_DIR+"/plots/Venn_diagram_taxonomy.png"
        conda: "Envs/plots.yml"
        params: out_dir=OUT_DIR
        shell: """ python venndiagram.py {params.out_dir}/final_results/HolomiRA_results.tsv {params.out_dir} """

rule plt_top:
        input: OUT_DIR+"/final_results/HolomiRA_results.tsv"
        output: expand("{out_dir}/plots/{env}_Top_10_MAG.png", out_dir=OUT_DIR, env=ENV)
        conda: "Envs/plots.yml"
        params: out_dir=OUT_DIR
        shell: """ python top.py {params.out_dir}/final_results/HolomiRA_results.tsv {params.out_dir} """

