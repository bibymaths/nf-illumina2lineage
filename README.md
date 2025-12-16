# nf-illumina2lineage

![GitHub last commit](https://img.shields.io/github/last-commit/bibymaths/nf-illumina2lineage)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.6-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15376065.svg)](https://doi.org/10.5281/zenodo.15376065)

A Nextflow pipeline for SARS-CoV-2 genome assembly and lineage assignment from Illumina paired-end sequencing data.

---

## Quickstart

```bash
git clone https://github.com/bibymaths/nf-illumina2lineage.git
cd nf-illumina2lineage
nextflow run main.nf
````

---

## Documentation

All documentation lives in `docs/` and is rendered [here](https://bibymaths.github.io/nf-illumina2lineage/).

* Getting started: `docs/quickstart.md`
* Workflow overview: `docs/workflow.md`
* Processes: `docs/processes.md`
* Parameters: `docs/parameters.md`
* Outputs: `docs/outputs.md`
* Lineage QC: `docs/lineage_qc.md`

---

## Inputs

* Illumina paired-end FASTQ.gz
* SARS-CoV-2 reference genome (auto-downloaded)

---

## Outputs

* QC reports (FastQC, Fastp, MultiQC)
* BAM and VCF files
* Consensus FASTA
* Pangolin lineage assignments
* Phylogenetic tree

See `docs/outputs.md` for details.

---

## License

BSD 3-Clause License. See `LICENSE`.

---

## Author

Abhinav Mishra
Email: [mishraabhinav@gmail.com](mailto:mishraabhinav@gmail.com)

---

## Acknowledgments

Developed during the SARS-2 Bioinformatics and Data Science course
at Freie Universität Berlin and the Robert Koch Institute.

---