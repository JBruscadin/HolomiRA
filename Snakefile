import pandas as pd
import os
import time

configfile: "/home/externo/tainacardoso/HolomiRA/config.yaml"
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
    input:
        expand(OUT_DIR + "/rnahybrid/{sample}_putative_targets.tsv", sample=sample),
        expand(OUT_DIR + "/rnahybrid/{sample}_bsites.tsv", sample=sample),
        expand("{out_dir}/rnahybrid/{sample}_finalresults.tsv", out_dir=OUT_DIR, sample=sample),
        OUT_DIR + "/final_results/HolomiRA_results.tsv",
        expand("{out_dir}/final_results/MAG_result_table_summary_miRNA_{env}.tsv", out_dir=OUT_DIR, env=ENV),
        OUT_DIR + "/plots/MAG_Histograms.png",
        OUT_DIR + "/plots/Venn_diagram_combined.png",
        expand("{out_dir}/plots/{env}_Top_20_miRNAs_and_MAGs.png", out_dir=OUT_DIR, env=ENV),
        config["out_dir"] + "/function/temp_merged_affected_cds.fasta",
        config["out_dir"] + "/function/temp_concat_end_MAG",
        config["out_dir"] + "/function/temp_concat_end_miRNA",
        config["out_dir"] + "/function/output_all_levels_and_function_done_mags.txt",
        config["out_dir"] + "/function/output_all_levels_and_function_done_mirna.txt"

        
rule annotate_prokka:
	input: FASTA_DIR+"{sample}.fa"
	output:
		OUT_DIR+"/annotation/{sample}/{sample}.gff"

	params: out_dir=OUT_DIR
	conda: "Envs/prokka.yml"
	shell: 
		"""prokka --quiet --outdir {params.out_dir}/annotation/{wildcards.sample} --prefix {wildcards.sample} --addgenes --centre X --compliant --cpus 12 {input} --force """

rule cds:
    input:
        expand("{out_dir}/annotation/{sample}/{sample}.gff", out_dir=OUT_DIR, sample=sample)
    output:
        OUT_DIR+"/annotation/{sample}/{sample}_cds.gff"
    shell:
        """ awk '$3 == "CDS" {{print $0}}' {input} > {output} """

rule id:
	input:	OUT_DIR+"/annotation/{sample}/{sample}_cds.gff"
	output:	OUT_DIR+"/annotation/{sample}/{sample}_IDs.txt"
	shell:	""" grep -v '^#' {input} | awk '{{print $1, $4, $5,$7,$9}}' | awk '{{split($5,a,/;/); print $1,$2,$3,a[1],$4}}' OFS="\t"> {output} """

rule five_prime:
    input: OUT_DIR + "/annotation/{sample}/{sample}_IDs.txt"
    params: fasta=FASTA_DIR, upstream=UPS, downstream=DWNS, out_dir=OUT_DIR
    output: OUT_DIR + "/target_fasta/{sample}_filtered.fa", OUT_DIR + "/target_fasta/{sample}_CDS.fa"
    conda: "Envs/fiveprime.yml"
    shell: "python Scripts/get_fiveprime.py {wildcards.sample} {params.fasta} {params.upstream} {params.downstream} {params.out_dir} {input}"

rule find_targets:
	input: 
		fasta=OUT_DIR+"/target_fasta/{sample}_filtered.fa",
		ref_mir=REF_MIR
	output: 
		OUT_DIR+"/rnahybrid/{sample}_putative_targets.tsv"
	params: 
		seed='-f '+SEED if SEED!='NA' else [],
		e="-e "+str(ENERGY),
		p="-p "+str(PVALUE)
	conda: "Envs/rnahybrid.yml"
	shell: """ RNAhybrid -s 3utr_human -c -t {input.fasta} -q {input.ref_mir} {params.seed} {params.e} {params.p} > {output} """

checkpoint format_rnahybrid:
    input:
        "{out_dir}/rnahybrid/{sample}_putative_targets.tsv"
    output:
        "{out_dir}/rnahybrid/{sample}_bsites.tsv"
    conda: "Envs/formatOutputs.yml"
    params: out_dir=OUT_DIR
    shell:
        """
        python Scripts/rnahybrid_format.py {wildcards.sample} {params.out_dir} {input} > {output}
        """
rule aggregate_bsites:
    input:
        expand("{out_dir}/rnahybrid/{sample}_bsites.tsv", out_dir=OUT_DIR, sample=sample)
    output:
        OUT_DIR + "/rnahybrid/all_bsites.txt"
    run:
        # Concatene todos os arquivos {sample}_bsites.tsv em um único arquivo
        with open(output[0], 'w') as output_file:
            for input_file in input:
                with open(input_file, 'r') as input_file:
                    output_file.write(input_file.read())
                    
rule get_indiv_metrics:
	input: OUT_DIR + "/rnahybrid/all_bsites.txt",
		list=config["sample_tab"]
	output: OUT_DIR+"/rnahybrid/{sample}_finalresults.tsv"
	conda: "Envs/formatOutputs.yml"
	params: out_dir=OUT_DIR,
		id=ID
	shell: """   python Scripts/get_metrics.py {params.out_dir} {params.id} {input[0]} {input.list} """



rule aggregate_finalresults:
    input:
        expand("{out_dir}/rnahybrid/{sample}_finalresults.tsv", out_dir=OUT_DIR, sample=sample)
    output:
        OUT_DIR + "/rnahybrid/finalresults.txt"
    run:
        if not os.path.exists(output[0]):
            with open(output[0], 'w') as output_file:
                for input_file in input:
                    with open(input_file, 'r') as input_file:
                        output_file.write(input_file.read())

rule get_global_metrics:
	input: files= OUT_DIR + "/rnahybrid/finalresults.txt",
		list=config["sample_tab"],
	output: OUT_DIR+"/final_results/HolomiRA_results.tsv"
        params: out_dir=OUT_DIR
	conda: "Envs/formatOutputs.yml"
	shell: """   python Scripts/get_all_metrics.py {params.out_dir} {input.list} {input.files} 
	touch {output} 
	"""

rule summary:
	input: OUT_DIR+"/final_results/HolomiRA_results.tsv"
	output:expand("{out_dir}/final_results/MAG_result_table_summary_miRNA_{env}.tsv", out_dir=OUT_DIR, env=ENV)
	conda: "Envs/plots.yml"
	params: out_dir=OUT_DIR
	shell: """ python Scripts/summary.py {params.out_dir}/final_results/HolomiRA_results.tsv {params.out_dir} """

rule plt_histogram:
        input: OUT_DIR+"/final_results/HolomiRA_results.tsv"
        output: OUT_DIR+"/plots/MAG_Histograms.png"
        conda: "Envs/plots.yml"
        params: out_dir=OUT_DIR
        shell: """ python Scripts/histogram.py {params.out_dir}/final_results/HolomiRA_results.tsv {params.out_dir} """

rule plt_venn:
        input: OUT_DIR+"/final_results/HolomiRA_results.tsv"
        output: OUT_DIR+"/plots/Venn_diagram_combined.png"
        conda: "Envs/plots.yml"
        params: out_dir=OUT_DIR
        shell: """ python Scripts/venndiagram.py {params.out_dir}/final_results/HolomiRA_results.tsv {params.out_dir} """

rule plt_top:
        input: OUT_DIR+"/final_results/HolomiRA_results.tsv"
        output: OUT_DIR+"/plots/{env}_Top_20_miRNAs_and_MAGs.png"
        conda: "Envs/plots.yml"
        params: out_dir=OUT_DIR
        shell: """ python Scripts/plots_byMir_byMAG.py {params.out_dir}/final_results/HolomiRA_results.tsv {params.out_dir} """

rule impacted:
    input:
        config["out_dir"] + "/final_results/HolomiRA_results.tsv",
    output:
        config["out_dir"] + "/function/temp_merged_affected_cds.fasta"
    params:
        out_dir=config["out_dir"]
    shell:
        """python Scripts/impacted.py {input} {params.out_dir}"""

rule prep_superfocus:
    input: config["out_dir"] + "/function/temp_merged_affected_cds.fasta"
    output: miRNA=config["out_dir"] + "/function/temp_concat_end_miRNA", MAG=config["out_dir"] + "/function/temp_concat_end_MAG"
    params: out_dir=config["out_dir"]
    shell: "python Scripts/prep_superfocus.py {params.out_dir}"

rule run_superfocus_MAG:
    input:
        config["out_dir"] + "/function/temp_concat_end_MAG"  
    output:
        config["out_dir"] + "/function/output_all_levels_and_function_done_mags.txt"
    params:
        out_dir = config["out_dir"]
    shell:
        """
       
        for dir in {params.out_dir}/function/MAGs_*; do
            # Verifica se é realmente um diretório
            if [ -d "$dir" ]; then
                echo "Processando $dir"
                superfocus -q $dir -dir $dir -a diamond
            fi
        done
        
        touch {output}
        """

rule run_superfocus_miRNA:
    input:
        config["out_dir"] + "/function/temp_concat_end_miRNA"  
    output:
        config["out_dir"] + "/function/output_all_levels_and_function_done_mirna.txt"
    params:
        out_dir = config["out_dir"]
    shell:
        """
        
        for dir in {params.out_dir}/function/miRNA_*; do
            
            if [ -d "$dir" ]; then
                echo "Processing $dir"
                superfocus -q $dir -dir $dir -a diamond
            fi
        done
        
        touch {output}
        """
