# HolomiRA - Holobiome miRNA Affinity Predictor

## Snakemake Workflow

HolomiRA is a comprehensive tool for miRNA binding site prediction in metagenomic assemblies. This repository contains the Snakemake workflow implementation for HolomiRA. Follow this guide to set up and run HolomiRA using Snakemake.

## Snakemake Installation
Snakemake recommends using mamba as the environment manager. If you don't have mamba installed, you can do so with the following command:
```shell
$ conda install -n base -c conda-forge mamba
```
Then, create the Snakemake environment using mamba:
```shell
$ conda activate base
$ mamba create -c conda-forge -c bioconda -n HolomiRA snakemake
```
To activate the Snakemake environment, run:

```shell
$  conda activate HolomiRA
```
For more details about Snakemake, refer to the [official documentation](https://snakemake.readthedocs.io/en/stable/index.html).

## HolomiRA Downloading and Setup 
**1-** Clone the HoloMirA repository and set up the environment:

```bash
git clone https://github.com/JBruscadin/HolomiRA.git
```

**2-** Customize HoloMirA by editing the configuration file - [config.yaml](https://github.com/JBruscadin/HolomiRA/blob/main/config.yaml) - with your preferred text editor. 

* **fasta_dir:** Path to the directory containing your fasta file.
* **out_dir:**  Path to the desired output directory.
* **sample_tab:** File listing all MAGs to be analyzed. 
* **upstream:** Upstream distance in base pairs before the CDS start position for miRNA binding site search.
* **downstream:** Downstream distance in base pairs after the CDS start position for miRNA binding site search.
* **ref_mir:** Host reference miRNA sequences.
* **seed** - RNAHybrid parameter. Comma-separated start nucleotide of the seed region and its length (e.g., 2,8 indicates a seed starting at nucleotide 2 with a length of 8).
* **energy** - RNAHybrid parameter. Maximum allowed minimum free energy for miRNA binding.
* **pvalue** -  RNAHybrid parameter. P-value threshold for target-miRNA interaction. 
* **id** - Tab-separated file containing MAG IDs, taxonomy, and environment/tissue info.
* **environment** - Python list of studied environments/tissues (e.g., ["feces", "rumen"]).
  
Learn more about RNAHybrid parameters in the tool’s  [manual](https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid?id=rnahybrid_manual_manual).

To run the enrichment analyses performed by Super-Focus, the prebuilt databases (clusters) can be downloaded in the developer's [Github page](/https://github.com/metageni/SUPER-FOCUS). 

## Running HolomiRA

Execute HoloMirA using Snakemake:
```bash
$ cd HolomiRA
$ snakemake -s Workflow/Snakefile --use-conda --cores N 
```
If your server supports a cluster system for parallel job execution, you can use these example commands (customize according to your resources):

Examples for different clusters:
```bash
$snakemake -s Snakefile --cluster 'sbatch -t 60 --mem=2g -c 1' -j 10
$snakemake -s Snakefile --cluster 'qsub -cwd -N HoloMira' -j 10
```
For more information about cluster execution in Snakemake, refer to the [documentation]( https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

## Processing steps

**Step 1** - Predict coding sequences from microbial genomes, implemented by Prodigal

**Step 2** - Select candidate target regions according to upstream and downstream distances described in the config file

**Step 3** - Search for host miRNAs binding sites in the candidate regions, implemented by RNAHybrid

**Step 4** - Get summarized metrics and plots from the results

**Step 5** - Functional Enrichment analysis, performed with Super-focus

**Step 6** - Comparative analysis between sample groups defined in the config file

**Additional Step** - Comparative analysis between different results sets (_i.e.,_ different host genomes)

## Output files

Five folders will be created in your output directory, and the output files in each folder are listed below:

**1 Annotation** 
* It contains directories named with sample prefixes, storing various files from Prokka-implemented annotation step. 
* In addition, HolomiRA creates new files:
  * Prefix_cds.gff contains contigs located in coding regions;
  * Prefix_cds_fiveprime.gff have the regions around the CDS start created based on the user-defined genomic ranges (input for miRNA target prediction);
  * Prefix_ID.txt:
    
**2 Target Fasta** 
* It contains sequence files in fasta format for all samples obtained from filtering the MAGs sequences for the regions obtained in the Annotation Step:
    * Prefix_CDS.fa is produced from filtering MAGs for the regions from Prefix_cds.gff;
    * Prefix_filtered.fa is produced from filtering MAGs for the regions from Prefix_cds_fiveprime.gff.

**3 Rnahybrid** 
* Prefixputative_targets.tsv
* Prefix_bsites.tsv
* Prefix_finalresults.tsv

**4 Final results** 
* HolomiRA_results.tsv
* MAG_result_table_summary_miRNA_{phenotype}.tsv
* MAG_result_table_summary_taxonomy_{phenotype}.tsv
* venn_summary.txt
  
**5 Plots**
* {phenotype}_Top_20_miRNAs_and_MAGs.png
* MAG_Histograms.png
* Venn_diagram_combined.png

## Acknowledgments
This research has been facilitated by the financial support of Brazillian institutions, namely FAPESP (Foundation for Research Support of the State of São Paulo), CAPES (Coordination for the Improvement of Higher Education Personnel), and CNPq (National Council for Scientific and Technological Development). The indispensable infrastructure and resources provided by UFSCar (Federal University of São Carlos) and EMBRAPA (Brazilian Agricultural Research Corporation), with particular emphasis on the Embrapa Pecuária Sudeste and Embrapa Informática Agropecuária units, played an integral role in enabling the execution of this work.


