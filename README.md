# HolomiRA - Holobiome miRNA Affinity Predictor

HolomiRA is a Snakemake-based pipeline for predicting host miRNA binding sites in microbial genomes. It automates all steps, from genome annotation to interaction prediction and functional enrichment, offering a user-friendly and fully reproducible workflow.

## Installation and Setup
1. Clone the repository

```shell
git clone https://github.com/JBruscadin/HolomiRA.git
cd HolomiRA
```

2. Create and activate the HolomiRA environment

HolomiRA uses Conda environments for dependency management.
If necessary, install Conda (≥ 24.7.0) (and optionally Mamba) before creating and activating the HolomiRA environment.

```shell
conda env create -f Workflow/Envs/HolomiRA_versions.yml
conda activate HolomiRA
```

If Mamba is available, we recommend using it for faster and more reliable environment resolution.

```shell
mamba env create -f Workflow/Envs/HolomiRA_versions.yml
conda activate HolomiRA
```

3. Check software versions

Check that all required tools are correctly installed and meet the minimum recommended versions:

* Snakemake	≥ 7.32.3
* Conda	≥ 24.1.0
* Mamba	≥ 1.5.7 (if used)
* Prokka	1.14.6
* Prodigal	2.6.3
* RNAHybrid	2.1.2
* RNAup / ViennaRNA	2.5.1
* SUPER-FOCUS	1.4.1
* DIAMOND	2.1.8
* unzip	Any recent version (tested with unzip 6.0)

If any tool does not meet the recommended version, please update your Conda environment or install the correct version manually.

For more details about Snakemake, refer to the [official documentation](https://snakemake.readthedocs.io/en/stable/index.html).

## Configuration

Edit the file: 

```shell
Config/config.yaml
```

Required parameters:

* **fasta_dir:** Path to the directory containing genome FASTA files.
* **ref_mir:** Path to the host miRNA sequences FASTA file.
* **sample_tab:** Sample list (one sample ID per line). 
* **id:** Metadata table with genome ID, taxonomy, etc.


Optional parameters (with defaults):
* **out_dir:**  Output directory (default: Results/)
* **upstream:** nt upstream of CDS start for binding site search (default: 15)
* **downstream:** nt downstream of CDS start (default: 20)
* **seed:** Comma-separated seed start and length (e.g., 2,8;1,7); (default: NA)
* **energy:** RNAHybrid energy cutoff (default: -20)
* **pvalue:** RNAHybrid p-value threshold (default: 0.01)
* **DGopen_cutoff:** RNAup ΔG total cutoff for accessibility (default: -10)


**SuperFocus Database Preparation**

HolomiRA uses DIAMOND + MMseqs2 clusters from SUPER-FOCUS for functional annotation.

We recommend downloading the 90_clusters.zip file manually via your web browser from the following link:

```bash
https://doi.org/10.6084/m9.figshare.29040107.v1
```

After downloading, place the file inside the Superfocus directory and unzip it:

```bash
mkdir -p Superfocus
mv /path/to/90_clusters.zip Superfocus/
cd Superfocus
unzip 90_clusters.zip -d 90_clusters
```

3. Download and format the SUPER-FOCUS database:
```bash
superfocus_downloadDB   --input 90_clusters/   -d $CONDA_PREFIX/lib/python3.9/site-packages/superfocus_app/   --aligner diamond   --clusters 90

```

For other versions and databases - [SUPER-FOCUS](https://github.com/metageni/SUPER-FOCUS).

## Running HolomiRA

From the root directory of the HolomiRA repository (where the Workflow/Snakefile is located), run:

```bash
snakemake -s Workflow/Snakefile --conda-frontend conda --cores N 
```
Where *N* is the number of CPU cores you want to use.

If your server supports a cluster system for parallel job execution, you can use these example commands (customize according to your resources):

Examples for different clusters:
```bash
snakemake -s Snakefile --cluster 'sbatch -t 60 --mem=2g -c 1' -j 10
snakemake -s Snakefile --cluster 'qsub -cwd -N HoloMira' -j 10
```
For more information about cluster execution in Snakemake, refer to the [documentation](https://snakemake.readthedocs.io/en/v7.19.1/executing/cluster.html).

* **Note 1**: RNAup can be memory-intensive when analyzing long sequences. If you encounter segmentation faults or buffer overflow errors (core dumped), try running the analysis on a machine with more available RAM.
 **Note 2**: If you encounter errors during SUPER-FOCUS steps, please delete the affected .m8 and .fasta files and re-run Snakemake.

## Workflow Steps

* **Step 1**: Predict CDS using Prokka
* **Step 2**: Extract target candidate regions
* **Step 3**: Predict miRNA binding (RNAHybrid + RNAup)
* **Step 4**: Summarize and visualize interactions
* **Step 5**: Perform functional enrichment (SUPER-FOCUS)


**Additional Steps** 

-- Additional Step 1 - Functional Enrichment Analysis (Extra Visualizations)

For complete description of parameters, please see [Functional annotation documentation](Additional_Steps/Functional_annotation/Documentation).

** This step uses default values for filtering: an abundance threshold of 0.01, a prevalence threshold of 0.2, and adjusted p-value ≤ 0.05 and a |log2 fold change| of at least 1.

Example workflow:

```bash
cd HolomiRA
mkdir ./Results/Functional_analysis
cd ./Results/Functional_analysis

#process_data
bash ../../Additional_Steps/Functional_annotation/process_data.sh ../../Results/function/MAGs_Feces/output_subsystem_level_3.xls relab MAG Feces
bash ../../Additional_Steps/Functional_annotation/process_data.sh ../../Results/function/MAGs_Rumen/output_subsystem_level_3.xls relab MAG Rumen

#analyse_by_dataset
Rscript ../../Additional_Steps/Functional_annotation/analyse_by_dataset.R level_3_MAG_Feces_relab.txt Feces MAG 0.01 0.2 10 blue  
Rscript ../../Additional_Steps/Functional_annotation/analyse_by_dataset.R level_3_MAG_Rumen_relab.txt Rumen MAG 0.01 0.2 10 red

#comparative_statistical_analysis
Rscript ../../Additional_Steps/Functional_annotation/comparative_statistical_analysis.R level_3_MAG_Feces_relab.txt Feces level_3_MAG_Rumen_relab.txt Rumen MAG 0.01 0.2 0.05 1 "#FF5733" "#33FF57" "purple"

#exclusive_functions
Rscript ../../Additional_Steps/Functional_annotation/exclusive_functions.R exclusive_functions_MAG_level_3_Feces.txt Feces MAG 10 green
Rscript ../../Additional_Steps/Functional_annotation/exclusive_functions.R exclusive_functions_MAG_level_3_Rumen.txt Rumen MAG 10 red

```

-- Additional Step 2 - Comparative Analysis Between Results
For cross-sample or species comparisons

```bash
cd ../../HolomiRA
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




