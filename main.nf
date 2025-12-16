#!/usr/bin/env nextflow

// Description: SARS-CoV-2 Amplicon Sequencing Analysis Pipeline
// Author: Abhinav Mishra
// License: BSD-3-Clause
// Source: https://github.com/bibymaths/nf-illumina2lineage

// Enable DSL2 syntax for modular workflow definitions
nextflow.enable.dsl=2

// ===================================================================================
// WORKFLOW DEFINITION
// ===================================================================================
workflow {

    // 1. Data Acquisition
    // Calls the downloadData process to fetch raw sequencing files
    raw_list_ch = downloadData()

    // Flatten the list of files emitted by downloadData into individual items
    raw_reads_ch = raw_list_ch
        .flatMap { it }

    // Fetch the reference genome (NC_045512.2)
    ref_ch = referenceGenome()

    // 2. Channel Manipulation: Prepare Read Pairs
    // This logic converts a stream of individual files into paired tuples: [sample_id, R1, R2]
    reads_ch = raw_reads_ch
      .map { file ->
        // Extract Sample ID by removing the suffix (e.g., SampleA_R1_001.fastq.gz -> SampleA)
        def id = file.name.replaceFirst(/_R[12]_001\.fastq\.gz$/, '')
        tuple(id, file)
      }
      // Group files by Sample ID (result: [id, [file1, file2]])
      .groupTuple()
      // explicit sorting/assignment to ensure R1 and R2 are correctly identified
      .map { sid, files ->
        def r1 = files.find{ it.name.contains('_R1_') }
        def r2 = files.find{ it.name.contains('_R2_') }
        tuple(sid, r1, r2)
      }

    // 3. Quality Control
    // Run FastQC, Fastp (trimming), and MultiQC
    // Returns cleaned reads (qc_reads) and reports
    def ( qc_reads, qc_html_ch, qc_json_ch ) = qc(reads_ch)

    // 4. Alignment
    // Map cleaned reads to reference using Minimap2
    mapping_ch = mapping(qc_reads, ref_ch)

    // 5. Primer Clipping
    // Remove amplicon primer sequences using BamClipper
    primer_ch = primerClipping(mapping_ch)

    // 6. Variant Calling
    // Call variants using Freebayes on the clipped BAM files
    vcf_ch = variantCalling(primer_ch, ref_ch)

    // 7. Consensus Generation
    // Extract VCF files from the channel for collection if needed
    vcf_list = vcf_ch.map { it[1] }.collect()

    // Generate consensus FASTA sequence for each sample individually
    individual_consensus_ch = consensusGeneration(vcf_ch, ref_ch)

    // Merge all individual consensus sequences into a single multi-FASTA file
    merged_consensus_ch = mergeConsensus(individual_consensus_ch.collect())

    // 8. Lineage Assignment
    // Run Pangolin on the merged consensus file
    pangolinLineage(merged_consensus_ch)

    // 9. Downstream Analysis
    // Run QC on consensus sequences (President) and build Phylogenetic tree (IQ-TREE)
    consensusQC(merged_consensus_ch, ref_ch)
    phylogeny(merged_consensus_ch)
}

// ===================================================================================
// PROCESS DEFINITIONS
// ===================================================================================

/*
 * Process: downloadData
 * Purpose: Downloads the raw Illumina sequencing data archive.
 * Tools: aria2c (download), tar/pigz (extraction).
 */
process downloadData {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    output:
      path params.reads, emit: raw_list_ch

    script:
    """
    set -euo pipefail

    # Download archive using aria2c (multi-connection download)
    aria2c -x16 -s16 -d . -o illumina-amplicon-capture-wgs.tar.gz \\
      https://osf.io/qu3bh/download

    # Extract the tar.gz archive
    tar -xzf illumina-amplicon-capture-wgs.tar.gz -C .

    # Organize files: Move all FASTQ files from nested folders to the current directory
    find . -mindepth 2 -type f -name "*.fastq.gz" -exec mv {} . \\;

    # Cleanup archive to save space
    rm -rf illumina-amplicon-capture-wgs.tar.gz __MACOSX
    """
}

/*
 * Process: referenceGenome
 * Purpose: Downloads the SARS-CoV-2 reference genome (NC_045512.2).
 * Tools: wget (using NCBI E-utilities).
 */
process referenceGenome {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    output:
      path "reference.fasta"

    script:
    """
    set -euo pipefail

    echo "Downloading SARS-CoV-2 Reference..."

    # Fetch FASTA directly from NCBI Nuccore database
    wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text" -O reference.fasta

    # Validation: Ensure the file is not empty
    if [ ! -s reference.fasta ]; then
        echo "ERROR: reference.fasta is empty. Download failed."
        exit 1
    fi

    # Log the first few lines for verification
    head -n 2 reference.fasta
    """
}

/*
 * Process: qc
 * Purpose: Quality control, adapter trimming, and reporting.
 * Tools: FastQC (pre/post), Fastp (trimming), MultiQC (aggregation).
 */
process qc {
    tag "$sample_id"
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    input:
      tuple val(sample_id), path(r1), path(r2)

    output:
      tuple val(sample_id), path("${sample_id}.R1.clean.fastq.gz"), path("${sample_id}.R2.clean.fastq.gz"), emit: qc_reads
      path "${sample_id}.fastp.html", emit: qc_html
      path "${sample_id}.fastp.json", emit: qc_json

    shell:
    """
    set -euo pipefail

    # 1. Pre-trimming QC
    fastqc -t ${task.cpus} "$r1" "$r2"

    # 2. Trimming and Cleaning with Fastp
    # --detect_adapter_for_pe: Auto-detect adapters
    # --cut_right: Sliding window quality filtering
    # --qualified_quality_phred 20: Minimum quality score
    fastp --detect_adapter_for_pe \\
          --overrepresentation_analysis \\
          --correction \\
          --qualified_quality_phred 20 \\
          --cut_right \\
          --thread ${task.cpus} \\
          --html  ${sample_id}.fastp.html \\
          --json  ${sample_id}.fastp.json \\
          -i "$r1" -I "$r2" \\
          -o ${sample_id}.R1.clean.fastq.gz \\
          -O ${sample_id}.R2.clean.fastq.gz

    # 3. Post-trimming QC
    fastqc -t ${task.cpus} \\
           ${sample_id}.R1.clean.fastq.gz \\
           ${sample_id}.R2.clean.fastq.gz

    # 4. Aggregate reports into MultiQC
    multiqc . -o .
    """
}

/*
 * Process: mapping
 * Purpose: Align cleaned reads to the reference genome.
 * Tools: Minimap2 (alignment), Samtools (sort/index).
 */
process mapping {
    tag "$sample_id"
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    input:
      tuple val(sample_id), path(r1), path(r2)
      path(ref)

    output:
      tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam_files

    script:
    """
    # Align reads using Minimap2 in short-read mode (-x sr)
    minimap2 -x sr -t ${task.cpus} -a -o ${sample_id}.sam ${ref} ${r1} ${r2}

    # Convert SAM to sorted BAM
    samtools view -bS ${sample_id}.sam | samtools sort -o ${sample_id}.sorted.bam

    # Index the sorted BAM for downstream tools
    samtools index ${sample_id}.sorted.bam
    """
}

/*
 * Process: primerClipping
 * Purpose: Remove PCR primer sequences from BAM files to prevent false positives.
 * Tools: bamclipper.sh, python (string formatting), awk (filtering).
 * Note: Uses CleanPlex BEDPE definitions.
 */
process primerClipping {
  tag "$sample_id"
  publishDir "${params.intermediate}/${task.process}", mode: 'copy'

  input:
    tuple val(sample_id), path(bam_file)

  output:
    tuple val(sample_id), path("${sample_id}.primerclipped.bam"), emit: primer_bam

  script:
  """
  set -euo pipefail

  # 1. Prepare Primer Bedfile
  # Download CleanPlex amplicon definitions
  wget -qO cleanplex.amplicons.bedpe -L https://osf.io/4nztj/download

  # Fix Windows line endings (\r\n) using Python
  python3 - <<'PY'
from pathlib import Path
p = Path("cleanplex.amplicons.bedpe")
b = p.read_bytes()
b = b.replace(b"\\r\\n", b"\\n").replace(b"\\r", b"\\n")
p.write_bytes(b)
PY

  # Filter BEDPE to keep only SARS-CoV-2 entries and normalize chromosome names
  awk -v OFS='\\t' '
    NF>=6 {
      gsub("NM_003194","NC_045512.2", \$1);
      gsub("NM_003194","NC_045512.2", \$4);
      if (\$1=="NC_045512.2" && \$4=="NC_045512.2") print \$1,\$2,\$3,\$4,\$5,\$6
    }
  ' cleanplex.amplicons.bedpe > SARSCoV2.amplicons.bedpe

  # 2. Prepare BAM for Clipper
  # Localize and re-index input BAM to ensure stable paths for the script
  cp "${bam_file}" input.bam
  samtools index input.bam

  # 3. Run BamClipper
  # Soft-clips primers based on the BEDPE file
  bamclipper.sh \
      -b input.bam \
      -p SARSCoV2.amplicons.bedpe \
      -n ${task.cpus} \
      -u 10 -d 10

  # Rename output to match pipeline convention
  mv input.primerclipped.bam "${sample_id}.primerclipped.bam"

  # Index the clipped BAM
  samtools index "${sample_id}.primerclipped.bam"
  """
}

/*
 * Process: variantCalling
 * Purpose: Identify genetic variants (SNPs/Indels) from the alignments.
 * Tool: Freebayes.
 */
process variantCalling {
  input:
    tuple val(sample_id), path(bam)
    path ref
  output:
    tuple val(sample_id), path("${sample_id}.vcf")
  script:
  """
  # -f: Reference fasta
  freebayes -f ${ref} "${bam}" > ${sample_id}.vcf
  """
}

/*
 * Process: consensusGeneration
 * Purpose: Generate a consensus FASTA sequence from the VCF and Reference.
 * Tools: bcftools (norm, index, consensus), sed.
 */
process consensusGeneration {
    tag "$sample_id"
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    input:
      tuple val(sample_id), path(vcf)
      path ref

    output:
      path "${sample_id}.consensus.fasta", emit: consensus_fasta

    script:
    """
    set -euo pipefail

    # [CRITICAL FIX]: Normalize VCF
    # Freebayes can produce overlapping variants that break 'bcftools consensus'.
    # 'bcftools norm' merges these multiallelic sites.
    bcftools norm -f "${ref}" -m -both "${vcf}" -o normalized.vcf

    # Compress and index the normalized VCF for bcftools consensus
    bcftools view normalized.vcf -Oz -o masked-strict.vcf.gz
    bcftools index masked-strict.vcf.gz

    # Apply variants to the reference genome to create consensus
    bcftools consensus -f "${ref}" masked-strict.vcf.gz -o "${sample_id}.consensus.fasta"

    # Rename the header in the FASTA file from Ref ID to Sample ID
    sed -i "s/NC_045512.2/${sample_id}/g" "${sample_id}.consensus.fasta"
    """
}

/*
 * Process: mergeConsensus
 * Purpose: Concatenate all individual consensus sequences into one file.
 * Tool: cat.
 */
process mergeConsensus {
  publishDir "${params.intermediate}/${task.process}", mode: 'copy'

  input:
    path fasta_files

  output:
    path "consensus-seqs.fasta"

  script:
  """
  set -euo pipefail

  # List files for debugging logs, allow failure if empty list (|| true)
  ls -lh *.fasta || true

  # Concatenate all input FASTAs
  cat *.fasta > consensus-seqs.fasta

  # Debugging: check size and headers
  echo "[DEBUG] merged file size:"
  ls -lh consensus-seqs.fasta

  echo "[DEBUG] headers:"
  grep '^>' consensus-seqs.fasta || true
  """
}

/*
 * Process: pangolinLineage
 * Purpose: Assign SARS-CoV-2 lineages (e.g., B.1.1.7) to the sequences.
 * Tool: Pangolin.
 * Dependencies: Requires 'bioconda::pangolin' and 'bioconda::pangolin-data'.
 */
process pangolinLineage {
    conda 'bioconda::pangolin bioconda::pangolin-data'

    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    input:
      path consensus_files

    output:
      path "lineage_report.csv"

    script:
    """
    echo "Running Pangolin on ${consensus_files}..."

    # Run Pangolin with specified threads
    # Output is a CSV file containing lineage assignments
    pangolin ${consensus_files} --threads ${task.cpus} --outfile lineage_report.csv
    """
}

/*
 * Process: consensusQC
 * Purpose: Quality Control of consensus sequences (check for Ns, etc.).
 * Tool: President.
 */
process consensusQC {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    input:
      path consensus_files
      path ref

    output:
      path "consensus_*"

    script:
    """
    mkdir -p results/output

    # Run President QC
    # -a: Calculate ALignments
    # -p: Output path
    # -f: Output file prefix
    president -r "${ref}" -q "${consensus_files}" -t ${task.cpus} \\
      -a -p . -f consensus_
    """
}

/*
 * Process: phylogeny
 * Purpose: Create a phylogenetic tree from the consensus sequences.
 * Tools: MAFFT (alignment), IQ-TREE (tree building).
 */
process phylogeny {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    input:
      path consensus_files

    output:
      path "phylo.treefile", optional: true

    script:
    """
    mkdir -p results

    # 1. Merge (redundant if input is already merged, but ensures single file)
    cat ${consensus_files} > all_sequences.fasta

    # 2. Safety Check: Skip if file is empty or too small (<1KB)
    if [ ! -s all_sequences.fasta ] || [ \$(stat -c%s all_sequences.fasta) -lt 1000 ];
    then
        echo "ERROR: Combined consensus file is empty. Skipping tree generation."
        exit 0
    fi

    # 3. Multiple Sequence Alignment (MSA) using MAFFT
    mafft --thread ${task.cpus} all_sequences.fasta > alignment.fasta

    # 4. Build Tree using IQ-TREE (Model: GTR)
    if [ -s alignment.fasta ];
    then
        iqtree -s alignment.fasta -nt ${task.cpus} -pre phylo -m GTR
    else
        echo "ERROR: Alignment failed or is empty."
    fi
    """
}