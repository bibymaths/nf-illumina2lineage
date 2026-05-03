This guide will help you get the pipeline up and running as quickly as possible. To use your own data, you must override
the input parameters or modify the downloadData process, as the current script hardcodes a download from OSF.

## Step 1: Clone the Repository

```bash
git clone https://github.com/bibymaths/nf-illumina2lineage.git
cd nf-illumina2lineage
```

## Step 2: Install Mamba (if not already installed)

```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh"
bash Mambaforge-Linux-x86_64.sh
```

## Step 3: Create and Activate the Environment

```bash
conda env create -f environment.yml
conda activate illumina2lineage
python -m pip install -r requirements.txt
```

## Set up Java and Nexflow if not already installed

```bash 
curl -s https://get.sdkman.io | bash 
```  

In a new terminal:

```bash 
sdk install java 17.0.10-tem  
java -version
```  

```bash 
curl -s https://get.nextflow.io | bash  
chmod +x nextflow 
mkdir -p $HOME/.local/bin/
mv nextflow $HOME/.local/bin/ 
nextflow info
``` 

Before you launch Nextflow, unset the Conda‐JDK variables and point to SDKMAN’s inside the conda environment.

```bash 
unset JAVA_CMD JAVA_HOME
export JAVA_HOME="$HOME/.sdkman/candidates/java/current"
export JAVA_CMD="$JAVA_HOME/bin/java" 
```

Nextflow requires Bash 3.2 (or later) and Java 17 (or later, up to 24) to be installed.
More information on [Nextflow installation](https://www.nextflow.io/docs/latest/install.html).

## Step 4: Run the Pipeline with Nextflow

```bash
nextflow run main.nf
```  

## Input Data Requirements

- Illumina paired-end FASTQ files
- SARS-CoV-2 reference genome (automatically downloaded by the pipeline)

## Output Overview

- Cleaned FASTQ files and QC reports
- Sorted, indexed BAM files
- VCF variant files
- Consensus sequences (FASTA)
- Pangolin lineage assignments
- Phylogenetic tree (Newick format)

For advanced configuration, see the [parameters](parameters.md) section.

> 💡 Tip: Use `multiqc` to summarize all QC results in one place.

You're now ready to use the pipeline!
