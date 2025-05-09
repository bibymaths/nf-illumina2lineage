# nf-illumina2lineage

A reproducible and modular **Nextflow pipeline** for **SARS-CoV-2 genome assembly and lineage analysis** from Illumina paired-end sequencing data.

## Overview

This pipeline automates:
- Read quality control
- Reference-based mapping
- Primer clipping
- Variant calling
- Consensus generation
- Lineage assignment
- Phylogenetic analysis

It is based on best-practice tools and developed as part of the *SARS-2 Bioinformatics & Data Science* course by Freie Universität Berlin and the Robert Koch Institute.

## Quickstart

```bash
git clone https://github.com/bibymaths/nf-illumina2lineage.git
cd nf-illumina2lineage

# Setup environment
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh"
bash Mambaforge-Linux-x86_64.sh
mamba env create -p ./envs/projectSARS --file environment.yaml
mamba activate ./envs/projectSARS

# Run the pipeline
nextflow run main.nf -profile docker
````

> 💡 See [docs/quickstart.md](docs/quickstart.md) for full details.

## Inputs

* Illumina paired-end `.fastq.gz` files
* SARS-CoV-2 reference genome (downloaded automatically)

## Outputs

* QC reports: FastQC, Fastp, MultiQC
* BAM & VCF files
* Consensus FASTA sequences
* Pangolin lineage annotations
* Phylogenetic tree (.treefile)

For a full output structure, see [docs/outputs.md](docs/outputs.md).

## Dependencies

Managed via `mamba` or `Docker`:

* QC: `fastqc`, `fastp`, `multiqc`
* Mapping: `minimap2`, `samtools`, `bamclipper`
* Variant Calling: `freebayes`, `vcftools`, `bcftools`
* Consensus: `vcfR`, `bcftools`, `president`
* Lineage & MSA: `pangolin`, `mafft`, `iqtree`

## Documentation

Complete documentation is available under the `docs/` folder and rendered via [MkDocs](https://www.mkdocs.org/). Includes:

* [Pipeline overview](docs/workflow.md)
* [Process details](docs/processes.md)
* [Parameters](docs/parameters.md)
* [Container usage](docs/containers.md)
* [Lineage QC](docs/lineage_qc.md)

## License

This project is licensed under the **BSD 3-Clause License**. See [LICENSE](LICENSE).

## Author

**Abhinav Mishra**
[mishraabhinav36@gmail.com](mailto:mishraabhinav36@gmail.com)

## Acknowledgments

Developed during the [SARS-2 Bioinformatics & Data Science](https://github.com/rki-mf1/2023-SC2-Data-Science) course at FU Berlin & RKI, under guidance of **Max von Kleist** and **Martin Hölzer**.

---
