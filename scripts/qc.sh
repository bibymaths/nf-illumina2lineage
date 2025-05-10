#!/usr/bin/env bash
set -euo pipefail

# Usage: qc.sh <input_dir> <output_dir> <threads>
# Example: ./qc.sh data datafiles 8

indir=${1:-data}
outdir=${2:-datafiles}
threads=${3:-8}

mkdir -p "$outdir"

# Find all R1 files, sorted
mapfile -t r1_files < <(ls "${indir}"/*_R1_001.fastq.gz 2>/dev/null | sort)
if [ "${#r1_files[@]}" -eq 0 ]; then
  echo "No FASTQ R1 files found in ${indir}" >&2
  exit 1
fi

i=1
for r1 in "${r1_files[@]}"; do
  r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
  sample="pair${i}"

  echo "[$(date +'%T')] QC on $sample:"
  echo "  R1: $r1"
  echo "  R2: $r2"

  # 1) Initial FastQC
  fastqc -t "$threads" "$r1" "$r2" -o "$outdir"

  # 2) fastp trimming & QC
  fastp \
    --thread "$threads" \
    --detect_adapter_for_pe \
    --overrepresentation_analysis \
    --correction \
    --qualified_quality_phred 20 \
    --cut_right \
    --html "${outdir}/${sample}.fastp.html" \
    --json "${outdir}/${sample}.fastp.json" \
    -i "$r1" -I "$r2" \
    -o "${outdir}/${sample}.R1.clean.fastq.gz" \
    -O "${outdir}/${sample}.R2.clean.fastq.gz"

  # 3) FastQC on cleaned reads
  fastqc -t "$threads" \
    "${outdir}/${sample}.R1.clean.fastq.gz" \
    "${outdir}/${sample}.R2.clean.fastq.gz" \
    -o "$outdir"

  i=$((i+1))
done

# 4) MultiQC over all reports in outdir
multiqc "$outdir" -o "$outdir"

echo "QC complete: processed ${i-1} samples."
