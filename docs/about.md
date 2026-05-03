This pipeline was developed as part of the **SARS-2 Bioinformatics & Data Science** course offered jointly by **Freie
Universität Berlin** and the **Robert Koch Institute**.

### Objectives

The aim was to create a reproducible, modular, and fully-documented pipeline for assembling and analyzing SARS-CoV-2
genomes from Illumina paired-end sequencing data using best-practice bioinformatics tools.

### Author

**Abhinav Mishra**  
Email: [mishraabhinav@gmail.com](mailto:mishraabhinav@gmail.com)

### Instructors

- **Max von Kleist**, Robert Koch Institute
- **Martin Hölzer**, Freie Universität Berlin

### Course Info

You can view on [GitHub](https://github.com/rki-mf1/2023-SC2-Data-Science)

### Device and Environment Tested On

#### Report details
- **Date generated:**                              2025-12-16 20:43:27

#### Hardware Information:
- **Hardware Model:**                              Dell Inc. OptiPlex 3050
- **Memory:**                                      32.0 GiB
- **Processor:**                                   Intel® Core™ i7-6700 × 8
- **Graphics:**                                    Intel® HD Graphics 530 (SKL GT2)
- **Disk Capacity:**                               1.0 TB

#### Software Information:
- **Firmware Version:**                            1.31.0
- **OS Name:**                                     Fedora Linux 43 (Workstation Edition)
- **OS Build:**                                    (null)
- **OS Type:**                                     64-bit
- **GNOME Version:**                               49
- **Windowing System:**                            Wayland
- **Kernel Version:**                              Linux 6.17.11-300.fc43.x86_64

### Execution Results

Below are the two execution result tables in Markdown format.

---

##### 1. Resource Usage by Process

| Process Name | Count | Avg Duration | Max Duration | Avg Memory (RSS) | Max Memory (RSS) | Avg CPU |
| --- | --- | --- | --- | --- | --- | --- |
| **downloadData** | 1 | 26.4s | 26.4s | 27.0 MB | 27.0 MB | 30.1% |
| **referenceGenome** | 1 | 889ms | 889ms | 12.0 MB | 12.0 MB | 8.7% |
| **qc** | 2 | 2m 45s | 2m 45s | 2.45 GB | 2.8 GB | 208.7% |
| **mapping** | 2 | 19.4s | 19.5s | 722.5 MB | 727.0 MB | 312.1% |
| **primerClipping** | 2 | 1.8s | 2.0s | 50.1 MB | 50.9 MB | 100.9% |
| **variantCalling** | 2 | 1.3s | 1.3s | 10.4 MB | 10.5 MB | 99.1% |
| **consensusGeneration** | 2 | 370ms | 385ms | 3.0 MB | 3.0 MB | 3632.5% |
| **mergeConsensus** | 1 | 225ms | 225ms | 0 MB | 0 MB | 0% |
| **phylogeny** | 1 | 1.1s | 1.1s | 6.8 MB | 6.8 MB | 0% |

---

##### 2. Detailed Task Execution Table

| Task ID | Process | Status | Realtime | %CPU | Peak Memory |
| --- | --- | --- | --- | --- | --- |
| 1 | referenceGenome | COMPLETED | 784ms | 8.7% | 12 MB |
| 2 | downloadData | COMPLETED | 26.3s | 30.1% | 27 MB |
| 3 | qc (200408_20-04246...) | COMPLETED | 1m 4s | 162.3% | 2.1 GB |
| 5 | qc (200423_20-04411...) | COMPLETED | 2m 44s | 255.0% | 2.8 GB |
| 6 | mapping (200408_20-04246...) | COMPLETED | 19.4s | 315.6% | 718 MB |
| 8 | mapping (200423_20-04411...) | COMPLETED | 19.5s | 308.5% | 727 MB |
| 9 | primerClipping (200408...) | COMPLETED | 2.0s | 102.6% | 50.9 MB |
| 10 | primerClipping (200423...) | COMPLETED | 1.6s | 99.1% | 49.3 MB |
| 11 | variantCalling (200408...) | COMPLETED | 1.3s | 98.8% | 10.4 MB |
| 12 | variantCalling (200423...) | COMPLETED | 1.2s | 99.3% | 10.5 MB |
| 13 | consensusGeneration | COMPLETED | 385ms | 4538.1% | 3 MB |
| 14 | mergeConsensus | COMPLETED | 225ms | 0% | 0 MB |
| 15 | phylogeny | COMPLETED | 1.1s | 0% | 6.8 MB |


### Acknowledgments

Special thanks to the course instructors and the contributing open-source bioinformatics community whose tools made this
project possible.

### License

This work is licensed under the BSD 3-Clause License.