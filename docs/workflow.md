# Pipeline Workflow Diagram and Overview

This section outlines the logical structure of the pipeline and provides a visual overview of each step.

## Workflow Summary
The pipeline consists of the following key stages:

1. **Environment Setup**: Initialize and configure the bioinformatics environment.
2. **Data Preparation**: Download SARS-CoV-2 reference genome and raw Illumina sequencing data.
3. **Quality Control**: Assess and clean raw sequencing reads using `fastqc`, `fastp`, and `multiqc`.
4. **Mapping**: Align cleaned reads to the SARS-CoV-2 reference genome using `minimap2` and process alignments with `samtools`.
5. **Primer Clipping**: Remove primer sequences from alignments using `bamclipper`.
6. **Variant Calling**: Call genetic variants with `freebayes`.
7. **Filtering & Masking**: Post-process variant calls using R scripts and `vcfR`.
8. **Consensus Generation**: Generate consensus sequences using `bcftools`.
9. **Lineage Annotation**: Annotate sequences using `pangolin`.
10. **Phylogenetic Analysis**: Perform multiple sequence alignment with `mafft` and infer phylogeny using `iqtree`.

## Workflow Diagram

![Pipeline Diagram](files/methodSARS.svg)

> **Note**: Each step is implemented as a separate Nextflow `process` and connected logically in `main.nf`.

For process-specific logic and input/output, see the [Process Details](processes.md) section.
