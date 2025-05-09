# Output Files and Directory Structure

The pipeline generates multiple outputs at each stage of analysis. Below is a summary of key files and their locations.

---

## Directory Structure
```
results/
├── qc/
│   ├── *.fastqc.html           # FastQC reports
│   ├── *.fastp.html/json       # Fastp reports
│   └── multiqc_report.html     # Summary report
├── alignments/
│   ├── *.sam                   # Raw alignments
│   ├── *.sorted.bam            # Sorted BAM files
│   └── *.bai                   # BAM index files
├── clipped/
│   └── *.primerclipped.bam     # Primer-removed alignments
├── variants/
│   └── *.vcf                   # Raw variant calls
├── consensus/
│   └── *.fasta                 # Consensus genome sequences
├── phylogeny/
│   ├── alignment.fasta         # MSA file
│   └── *.treefile              # Phylogenetic tree
├── lineage/
│   └── lineage_report.csv      # Pangolin lineage calls
└── qc_consensus/
    └── output/                 # PRESIDENT quality reports
```

---

## Output Types

| Stage              | File Type(s)        | Description |
|-------------------|---------------------|-------------|
| Quality Control    | `.html`, `.json`   | QC and trimming reports |
| Mapping            | `.bam`, `.bai`     | Aligned reads, indexed |
| Primer Clipping    | `.primerclipped.bam`| Primer-trimmed alignments |
| Variant Calling    | `.vcf`             | Raw variant calls |
| Consensus          | `.fasta`           | Consensus genome sequences |
| Lineage Annotation | `.csv`             | Pangolin lineage table |
| Phylogenetic Tree  | `.fasta`, `.treefile` | MSA + tree structure |
| QC (Consensus)     | `output/` files    | Nucleotide identity, ambiguities |

> 📁 All files are written into structured subdirectories under `results/` by default.

You can modify the output paths via `nextflow.config` or process-specific directives.
