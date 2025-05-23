# Lineage Annotation and Consensus Quality Control

This section explains how the pipeline handles lineage annotation and consensus sequence quality assessment.

---

## Lineage Annotation with Pangolin

`pangolin` is used to assign SARS-CoV-2 lineages based on the consensus sequences generated by the pipeline.

### Input
- `consensus-seqs.fasta` (combined consensus sequences for all samples)

### Output
- `lineage_report.csv` — Pangolin lineage calls with sample metadata
- Additional summary files: `summary.csv`, `pangolin.json` (optional)

### Usage
Pangolin is executed via the following command in the pipeline:
```bash
pangolin -t 8 consensus-seqs.fasta
```

---

## Consensus QC with PRESIDENT

[PRESIDENT](https://github.com/rki-mf1/president) evaluates consensus sequences by comparing them to the reference genome.

### Checks Performed
- Pairwise nucleotide identity
- Ambiguous base counts
- Masked regions

### Input
- `reference.fasta`
- `consensus-seqs.fasta`

### Output
- Summary reports written to the `output/` directory
- `.html` and tabular output formats depending on PRESIDENT version

### Usage
Executed as part of the pipeline:
```bash
president -r reference.fasta -q consensus-seqs.fasta -t 8 -a -p output/ -f consensus_
```

> 📈 These results help assess sequence reliability before downstream analysis like phylogenetics.