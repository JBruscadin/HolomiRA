# HolomiRA - Holobiome miRNA Affinity Predictor

## Snakemake Workflow

HoloMirA is a comprehensive tool for miRNA binding site prediction in metagenomic assemblies. This repository contains the Snakemake workflow implementation for HoloMirA. Follow this guide to set up and run HoloMirA using Snakemake.

## Snakemake Installation
Snakemake recommends using mamba as the environment manager. If you don't have mamba installed, you can do so with the following command:
```shell
$ conda install -n base -c conda-forge mamba
```
Then, create the Snakemake environment using mamba:
```shell
$ conda activate base
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
To activate the Snakemake environment, run:

```shell
$  conda activate snakemake
```
For more details about Snakemake, refer to the [official documentation](https://snakemake.readthedocs.io/en/stable/index.html).

## HolomiRA Downloading and Setup 
**1-** Clone the HoloMirA repository and set up the environment:

```bash
git clone https://github.com/JBruscadin/HolomiRA.git
```

**2-** Customize HoloMirA by editing the configuration file - [config.yaml] (https://github.com/JBruscadin/HolomiRA/blob/main/config.yaml) - with your preferred text editor. 

* **fasta_dir:** Path to the directory containing your fasta file.
* **out_dir:**  Path to the desired output directory.
* **sample_tab:** File listing all MAGs to be analyzed. 
* **upstream:** Upstream distance in base pairs before the CDS start position for miRNA binding site search.
* **downstream:** Downstream distance in base pairs after the CDS start position for miRNA binding site search.
* **ref_mir:** Host reference miRNA sequences.
* **seed** - RNAHybrid parameter. Start nucleotide of the seed region and its length (e.g., 2,8 indicates a seed starting at nucleotide 2 with a length of 8).
* **energy** - RNAHybrid parameter. Maximum allowed minimum free energy for miRNA binding.
* **pvalue** -  RNAHybrid parameter. P-value threshold for target-miRNA interaction. 
* **id** - Tab-separated file containing MAG IDs, taxonomy, and environment/tissue info. See [example](path_to_example_file) 
* **environment** - Python list of studied environments/tissues (e.g., ["fezes", "rumen"]).
  
Learn more about RNAHybrid parameters in the tool’s  [manual](https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid?id=rnahybrid_manual_manual).

## Running HolomiRA

Execute HoloMirA using Snakemake:
```bash
$snakemake -s Snakefile --use-conda --cores N 
```
If your server supports a cluster system for parallel job execution, you can use these example commands (customize according to your resources):

Examples for different clusters:
```bash
$snakemake -s Snakefile --cluster 'sbatch -t 60 --mem=2g -c 1' -j 10
$snakemake -s Snakefile --cluster 'qsub -cwd -N HoloMira' -j 10
```
For more information about cluster execution in Snakemake, refer to the [documentation]( https://snakemake.readthedocs.io/en/stable/executing/cluster.html).




