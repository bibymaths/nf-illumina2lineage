#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl=2


workflow {

    raw_list_ch = downloadData()

    raw_reads_ch = raw_list_ch
        .flatMap { it }

    ref_ch = referenceGenome()

    reads_ch = raw_reads_ch
      .map { file ->
        def id = file.name.replaceFirst(/_R[12]_001\.fastq\.gz$/, '')
        tuple(id, file)
      }
      .groupTuple()
      .map { sid, files ->
        def r1 = files.find{ it.name.contains('_R1_') }
        def r2 = files.find{ it.name.contains('_R2_') }
        tuple(sid, r1, r2)
      }

    def ( qc_reads, qc_html_ch, qc_json_ch ) = qc(reads_ch)

    mapping_ch = mapping(qc_reads, ref_ch)
    primer_ch    = primerClipping(mapping_ch)
    vcf_ch = variantCalling(primer_ch, ref_ch)

    // 1. Prepare inputs for QC: extract just the files and collect them into a list
    vcf_list = vcf_ch.map { it[1] }.collect()

    // Generate per-sample consensus
    individual_consensus_ch = consensusGeneration(vcf_ch, ref_ch)

    // Merge
    merged_consensus_ch = mergeConsensus(individual_consensus_ch.collect())

    // 3. Run Pangolin on the merged file
    pangolinLineage(merged_consensus_ch)

    // 4. Update consensusQC and phylogeny to use the merged file as well
    consensusQC(merged_consensus_ch, ref_ch)
    phylogeny(merged_consensus_ch)
}

process downloadData {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    output:
      path params.reads, emit: raw_list_ch

    script:
    """
    set -euo pipefail

    # Download into current directory
    aria2c -x16 -s16 -d . -o illumina-amplicon-capture-wgs.tar.gz \\
      https://osf.io/qu3bh/download

    # Unpack into cwd
    # pigz -dc illumina-amplicon-capture-wgs.tar.gz | tar -x -C .

    tar -xzf illumina-amplicon-capture-wgs.tar.gz -C .

    # Move any fastq.gz out of nested dirs into cwd
    find . -mindepth 2 -type f -name "*.fastq.gz" -exec mv {} . \\;

    # Clean up
    rm -rf illumina-amplicon-capture-wgs.tar.gz __MACOSX
    """
}

process referenceGenome {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    output:
      path "reference.fasta"

    script:
    """
    set -euo pipefail

    echo "Downloading SARS-CoV-2 Reference..."

    # Use a direct, stable link from NCBI (E-utilities)
    wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text" -O reference.fasta

    # SAFETY CHECK: Fail if the file is empty
    if [ ! -s reference.fasta ]; then
        echo "ERROR: reference.fasta is empty. Download failed."
        exit 1
    fi

    # Check the first few lines to ensure it looks like a FASTA
    head -n 2 reference.fasta
    """
}


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

    # raw QC
    fastqc -t ${task.cpus} "$r1" "$r2"

    # trimming & per-sample QC
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

    # QC on cleaned
    fastqc -t ${task.cpus} \\
           ${sample_id}.R1.clean.fastq.gz \\
           ${sample_id}.R2.clean.fastq.gz

    # aggregate
    multiqc . -o .
    """
}


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
    minimap2 -x sr -t ${task.cpus} -a -o ${sample_id}.sam ${ref} ${r1} ${r2}
    samtools view -bS ${sample_id}.sam | samtools sort -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}

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

  # Download BEDPE
  wget -qO cleanplex.amplicons.bedpe -L https://osf.io/4nztj/download

  # Python script to fix line endings
  python3 - <<'PY'
from pathlib import Path
p = Path("cleanplex.amplicons.bedpe")
b = p.read_bytes()
b = b.replace(b"\\r\\n", b"\\n").replace(b"\\r", b"\\n")
p.write_bytes(b)
PY

  # AWK processing to filter for SARS-CoV-2
  awk -v OFS='\\t' '
    NF>=6 {
      gsub("NM_003194","NC_045512.2", \$1);
      gsub("NM_003194","NC_045512.2", \$4);
      if (\$1=="NC_045512.2" && \$4=="NC_045512.2") print \$1,\$2,\$3,\$4,\$5,\$6
    }
  ' cleanplex.amplicons.bedpe > SARSCoV2.amplicons.bedpe

  # FIX: Localize and Index the input BAM
  # bamclipper requires the index, and symlinks from other dirs can sometimes cause path issues
  cp "${bam_file}" input.bam
  samtools index input.bam

  # Run bamclipper on the local indexed file
  # This produces 'input.primerclipped.bam'
  bamclipper.sh \
      -b input.bam \
      -p SARSCoV2.amplicons.bedpe \
      -n ${task.cpus} \
      -u 10 -d 10

  # Rename the output to the expected sample ID format
  mv input.primerclipped.bam "${sample_id}.primerclipped.bam"

  # Index the result (good practice for downstream steps)
  samtools index "${sample_id}.primerclipped.bam"
  """
}

process variantCalling {
  input:
    tuple val(sample_id), path(bam)
    path ref
  output:
    tuple val(sample_id), path("${sample_id}.vcf")
  script:
  """
  freebayes -f ${ref} "${bam}" > ${sample_id}.vcf
  """
}

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

    # [CRITICAL FIX]
    # Freebayes produces overlapping variants that crash 'bcftools consensus'.
    # We MUST normalize first to fix these overlaps.
    bcftools norm -f "${ref}" -m -both "${vcf}" -o normalized.vcf

    # Now run your simple steps using the FIXED vcf
    bcftools view normalized.vcf -Oz -o masked-strict.vcf.gz
    bcftools index masked-strict.vcf.gz

    bcftools consensus -f "${ref}" masked-strict.vcf.gz -o "${sample_id}.consensus.fasta"

    sed -i "s/NC_045512.2/${sample_id}/g" "${sample_id}.consensus.fasta"
    """
}

process mergeConsensus {
  publishDir "${params.intermediate}/${task.process}", mode: 'copy'

  input:
    path fasta_files

  output:
    path "consensus-seqs.fasta"

  script:
  """
  set -euo pipefail
  ls -lh *.fasta || true
  cat *.fasta > consensus-seqs.fasta

  echo "[DEBUG] merged file size:"
  ls -lh consensus-seqs.fasta

  echo "[DEBUG] headers:"
  grep '^>' consensus-seqs.fasta || true
  """
}

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

    # Run Pangolin
    # -t: threads
    # --outfile: specifies the output CSV name
    pangolin ${consensus_files} --threads ${task.cpus} --outfile lineage_report.csv
    """
}

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
    president -r "${ref}" -q "${consensus_files}" -t ${task.cpus} \\
      -a -p . -f consensus_
    """
}

process phylogeny {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    input:
      path consensus_files

    output:
      path "phylo.treefile", optional: true

    script:
    """
    mkdir -p results

    # 1. Merge all consensus files into one
    cat ${consensus_files} > all_sequences.fasta

    # 2. Check if file is empty or too small (less than 1KB)
    if [ ! -s all_sequences.fasta ] || [ \$(stat -c%s all_sequences.fasta) -lt 1000 ]; then
        echo "ERROR: Combined consensus file is empty. Skipping tree generation."
        exit 0
    fi

    # 3. Align
    mafft --thread ${task.cpus} all_sequences.fasta > alignment.fasta

    # 4. Build Tree (Check if alignment succeeded)
    if [ -s alignment.fasta ]; then
        iqtree -s alignment.fasta -nt ${task.cpus} -pre phylo -m GTR
    else
        echo "ERROR: Alignment failed or is empty."
    fi
    """
}