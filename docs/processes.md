# Process Details

This section describes each Nextflow `process` in the pipeline, including their inputs, outputs, and core functionality.

---

## `downloadData`
**Purpose**: Download Illumina sequencing data archive and extract contents.
- **Input**: None
- **Output**: Raw `.fastq.gz` files

---

## `referenceGenome`
**Purpose**: Download SARS-CoV-2 reference genome (NC_045512.2) using NCBI Entrez Direct.
- **Input**: None
- **Output**: `reference.fasta`

---

## `qc`
**Purpose**: Quality control and cleaning of raw sequencing reads.
- **Input**: Paired-end FASTQ files
- **Output**:
  - Cleaned FASTQ files (`pair*.R1/2.clean.fastq.gz`)
  - FastQC reports (`.html`, `.zip`)
  - Fastp reports (`.json`, `.html`)
  - MultiQC summary

---

## `mapping`
**Purpose**: Align reads to reference genome and produce sorted, indexed BAM files.
- **Input**: Cleaned FASTQ files, `reference.fasta`
- **Output**: Sorted and indexed BAM files (`.sorted.bam`, `.bai`)

---

## `primerClipping`
**Purpose**: Clip primer sequences using a CleanPlex BEDPE file and `bamclipper`.
- **Input**: Sorted BAM files
- **Output**: Primer-clipped BAM files (`.primerclipped.bam`)

---

## `variantCalling`
**Purpose**: Call variants from aligned BAM files using `freebayes`.
- **Input**: Primer-clipped BAM files, `reference.fasta`
- **Output**: Raw VCF files (`freebayes-illumina*.vcf`)

---

## `consensusGeneration`
**Purpose**: Generate consensus sequences using `bcftools` from VCFs.
- **Input**: VCF files, `reference.fasta`
- **Output**: FASTA files for consensus sequences (`consensus-*.fasta`, `consensus-seqs.fasta`)

---

## `pangolinLineage`
**Purpose**: Assign SARS-CoV-2 lineages using `pangolin`.
- **Input**: Combined consensus FASTA file
- **Output**: `lineage_report.csv`, optional TSV and summary files

---

## `consensusQC`
**Purpose**: Assess consensus quality using `president`.
- **Input**: Consensus FASTA, `reference.fasta`
- **Output**: QC summary in `output/` folder

---

## `phylogeny`
**Purpose**: Perform multiple sequence alignment and build phylogenetic tree.
- **Input**: Consensus FASTA file
- **Output**: `alignment.fasta`, IQ-TREE outputs (e.g., `.treefile`, `.log`)


> For a visual overview, see the [Workflow](workflow.md) section.