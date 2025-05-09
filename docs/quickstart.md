# Quickstart Guide

This guide will help you get the pipeline up and running as quickly as possible.

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
conda update -y conda
mamba env create -p ./envs/projectSARS --file environment.yaml
mamba activate ./envs/projectSARS
```

## Step 4: Run the Pipeline with Nextflow
```bash
nextflow run main.nf
```

To use a Docker profile:
```bash
nextflow run main.nf -profile docker
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
