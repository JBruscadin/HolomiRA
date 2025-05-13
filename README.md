# HolomiRA - Holobiome miRNA Affinity Predictor

HolomiRA is a Snakemake-based pipeline for predicting host miRNA binding sites in microbial genomes. It automates all steps, from genome annotation to interaction prediction and functional enrichment, offering a user-friendly and fully reproducible workflow.

## HolomiRA Downloading and Setup 
1. Install Snakemake and Mamba

We recommend using mamba for efficient environment management:

```shell
conda install -n base -c conda-forge mamba
```

2. Clone the repository and set up the environment

```shell
git clone https://github.com/JBruscadin/HolomiRA.git
cd HolomiRA
mamba env create -f Workflow/Envs/HolomiRA_versions.yml
conda activate HolomiRA
```

For more details about Snakemake, refer to the [official documentation](https://snakemake.readthedocs.io/en/stable/index.html).

## Configuration

Edit the file: Config/config.yaml

Required parameters:

* **fasta_dir:** Path to the directory containing genomes'fasta files.
* **ref_mir:** Path for the host miRNA sequences'fasta.
* **sample_tab:** Sample list (one sample ID per line) 
* **id:** Metadata table with genome ID, taxonomy, etc.


Optional parameters (with defaults):
* **out_dir:**  Output directory (default: Results/)
* **upstream:** nt upstream of CDS start for binding site search (default: 15)
* **downstream:** nt downstream of CDS start (default: 20)
* **seed:** Comma-separated seed start and length (e.g., 2,8;1,7); (default: NA)
* **energy:** RNAHybrid energy cutoff (default: -20)
* **pvalue:** RNAHybrid p-value threshold (default: 0.01)
* **DGopen_cutoff:** RNAup ΔG total cutoff for accessibility (default: -10)


**HolomiRA uses DIAMOND + MMseqs2 clusters. Download and format**

```bash
mkdir Superfocus
cd Superfocus
wget https://figshare.com/ndownloader/files/54459941?private_link=fb65bc6be0fb68ebbaf2
unzip '54459941?private_link=fb65bc6be0fb68ebbaf2' -d 90_clusters/
superfocus_downloadDB -i ./Superfocus -a diamond -c 90
```

For other versions and databases - [SUPER-FOCUS](https://github.com/metageni/SUPER-FOCUS).

## Running HolomiRA

Execute HolomiRA:
```bash
cd HolomiRA
snakemake -s Workflow/Snakefile --use-conda --cores N 
```
If your server supports a cluster system for parallel job execution, you can use these example commands (customize according to your resources):

Examples for different clusters:
```bash
snakemake -s Snakefile --cluster 'sbatch -t 60 --mem=2g -c 1' -j 10
snakemake -s Snakefile --cluster 'qsub -cwd -N HoloMira' -j 10
```
For more information about cluster execution in Snakemake, refer to the [documentation](https://snakemake.readthedocs.io/en/v7.19.1/executing/cluster.html).

* **Attention**: If error is given during SUPER-FOCUS analysis, please delete *m8 and its corresponding *.fasta and run snakemake again.
  
## Processing steps

**Step 1**: Predict CDS using Prokka
**Step 2**: Extract target candidate regions
**Step 3**: Predict miRNA binding (RNAHybrid + RNAup)
**Step 4**: Summarize and visualize interactions
**Step 5**: Perform functional enrichment (SUPER-FOCUS)


**Additional Steps** 

-- Additional Step 1 - Functional Enrichment Analysis (Extra Visualizations)

Usage:


./process_data.sh <subsystem_level_X.xls> <counts|relab> <sample_type>
<dataset_name>

Rscript analyse_by_dataset.R <input_file> <dataset_name> <sample_type>
<abundance_threshold> <prevalence_threshold> <top> <color_data>

Rscript comparative_statistical_analysis.R <input_data1> <dataset_name1> <input_data2>
<dataset_name2> <sample_type> <abundance_threshold> <prevalence_threshold> <padj>
<log_threshold> <color_data1> <color_data2> <heatmap_color>

Rscript exclusive_functions.R <input_file> <dataset_name> <sample_type> <top>
<color_graph>

** This step uses default values for filtering: an abundance threshold of 0.01, a prevalence threshold of 0.2, and adjusted p-value ≤ 0.05 and a |log2 fold change| of at least 1.

Practical exemplo:

```bash
cd HolomiRA
mkdir ./Results/Functional_analysis
cd ./Results/Functional_analysis

#process_data
bash ../../Additional_Steps/Functional_annotation/process_data.sh ../../Results/function/microbial genomes_Feces/output_subsystem_level_3.xls relab MAG Feces
bash ../../Additional_Steps/Functional_annotation/process_data.sh ../../Results/function/microbial genomes_Rumen/output_subsystem_level_3.xls relab MAG Rumen

#analyse_by_dataset
Rscript ../../Additional_Steps/Functional_annotation/analyse_by_dataset.R --input_file level_3_MAG_Feces_relab.txt --dataset_name Feces --sample_type MAG --abundance_threshold 0.01 --prevalence_threshold 0.2 --top 10 --color_data blue  
Rscript ../../Additional_Steps/Functional_annotation/analyse_by_dataset.R level_3_MAG_Rumen_relab.txt --dataset_name Rumen --sample_type MAG --abundance_threshold 0.01 --prevalence_threshold 0.2 --top 10 --color_data red

#comparative_statistical_analysis
Rscript ../../Additional_Steps/Functional_annotation/comparative_statistical_analysis.R --input_data1 level_3_MAG_Feces_relab.txt --input_data2 level_3_MAG_Rumen_relab.txt --dataset_name1 Feces --dataset_name2 Rumen --sample_type MAG --abundance_threshold 0.01 --prevalence_threshold 0.2 --padj 0.05 --log_threshold 1 --color_data1 "#FF5733" --color_data2 "#33FF57" --heatmap_color "purple"

#exclusive_functions
Rscript ../../Additional_Steps/Functional_annotation/exclusive_functions.R --input_file exclusive_functions_MAG_level_3_decrease_RME.txt --dataset_name decrease_RME --sample_type MAG --top 10 --color_graph green
Rscript ../../Additional_Steps/Functional_annotation/exclusive_functions.R --input_file exclusive_functions_MAG_level_3_increase_RME.txt --dataset_name increase_RME  --sample_type MAG --top 10 --color_graph red

```
For complete description of parameters, please see [Functional annotation documentation](Additional_Steps/Functional_annotation/Documentation).

-- Additional Step 2 - Comparative Analysis Between Results
```bash
cd HolomiRA
python Additional_Steps/Comparison_species/histogram.py
python Additional_Steps/Comparison_species/venndiagram.py
```

## Output files

Each subfolder in Results/ corresponds to a specific step. Example contents:

* **Annotation/**: Prokka results, GFF, CDS coordinates
* **Target_Fasta/**: CDS and filtered sequences
* **Rnahybrid/**: RNAHybrid output (putative target sites)
* **Structure/**: Pre-RNAup formatting files
* **RNAup/**: Accessibility results
* **Final_results/**: HolomiRA results + summary tables
* **Plots/**: miRNA-target genome visuals
* **Function/**: SuperFocus output by phenotype

When running Additional Step 1, these files are added:

* **Functional_analysis/**: Differential results and heatmaps

When running Additional Step 2, these files are added: miRNA-target genome visuals

## Acknowledgments

This research has been facilitated by the financial support of Brazillian institutions, namely FAPESP (Foundation for Research Support of the State of São Paulo), CAPES (Coordination for the Improvement of Higher Education Personnel), and CNPq (National Council for Scientific and Technological Development). The indispensable infrastructure and resources provided by UFSCar (Federal University of São Carlos) and EMBRAPA (Brazilian Agricultural Research Corporation), with particular emphasis on the Embrapa Pecuária Sudeste and Embrapa Informática Agropecuária units, played an integral role in enabling the execution of this work.

## Citation

If you use HolomiRA in your research, please cite:

Bruscadin, J. et al. HolomiRA: a reproducible pipeline for miRNA binding site prediction in microbial genomes. (2024) GitHub: https://github.com/JBruscadin/HolomiRA




