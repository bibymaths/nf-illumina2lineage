This section lists the configurable parameters used in the pipeline. These can be adjusted using a custom configuration
file or passed directly via the command line.

---

## Core Parameters

| Parameter   | Default           | Description                                      |
|-------------|-------------------|--------------------------------------------------|
| `reads`     | `*.fastq.gz`      | Glob pattern to match input FASTQ files          |
| `reference` | `reference.fasta` | Reference genome in FASTA format                 |
| `threads`   | `8`               | Number of threads to use for multithreaded tools |

---

## Default Configuration

```bash
nextflow run main.nf -profile pixi
``` 

or 

```bash
pixi run run-pipeline
``` 

> 💡 All parameters can also be set directly via `-params` flag or in `nextflow.config`.
