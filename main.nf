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

    // 2. Define where your R script is located
    // Make sure 'scripts/vcf_qc_masking.R' exists in your project folder!
    r_script = file("${projectDir}/scripts/vcf_qc_masking.R")

    // 3. Run the QC
    qc_masked = vcfMaskingQC(vcf_list, ref_ch, r_script)

    // 1. Generate consensus for each sample in parallel
    individual_consensus_ch = consensusGeneration(vcf_ch, ref_ch)

    // 2. Merge all individual FASTAs into one file for Pangolin
    merged_consensus_ch = individual_consensus_ch.collectFile(name: 'consensus-seqs.fasta')

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
      path params.reference

    script:
    """
    set -euo pipefail
    mkdir -p \$(dirname ${params.reference})

    wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O - | bash
    export PATH=\$HOME/edirect:\$PATH

    esearch -db nucleotide -query "NC_045512.2" | efetch -format fasta \\
      > ${params.reference}
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
  freebayes -f ${ref} --min-alternate-count 10 \
    --min-alternate-fraction 0.1 --min-coverage 20 \
    --pooled-continuous --haplotype-length -1 "${bam}" \
    > ${sample_id}.vcf
  """
}

process vcfMaskingQC {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    input:
      path vcf_files
      path "reference.fasta"
      path r_script

    output:
      path "pair*/masked-strict.vcf", emit: masked_vcfs
      path "pair*/*.pdf",             emit: qc_plots

    script:
    """
    mkdir -p pair1 pair2 pair3 scripts

    cp "${r_script}" scripts/vcf_qc_masking.R

    # FIX: Safely clean the header using a temporary file
    # This keeps only the ID (first column) and removes the description
    cut -d ' ' -f 1 reference.fasta > reference.tmp
    mv reference.tmp reference.fasta

    # Dynamically assign the input VCFs
    i=1
    for vcf in \$(ls *.vcf | sort); do
        if [ "\$i" -gt 3 ]; then break; fi
        cp "\$vcf" "pair\$i/freebayes-illumina.vcf"
        i=\$((i+1))
    done

    Rscript scripts/vcf_qc_masking.R
    """
}

process consensusGeneration {
    tag "$sample_id"
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    input:
      tuple val(sample_id), path(vcf)
      path ref

    output:
      path "${sample_id}.consensus.fasta"

    script:
    """
    set -u
    # set -e is removed to allow conditional handling, but we must handle pipe errors manually

    # FIX: Add '|| true' so if grep finds 0 lines (exit code 1), the script continues
    VAR_COUNT=\$(grep -v "^#" "${vcf}" | wc -l || true)

    if [ "\$VAR_COUNT" -eq 0 ]; then
        echo "No variants found for ${sample_id}. Using reference as consensus."
        sed "s/^>NC_045512.2/>${sample_id}/" "${ref}" > "${sample_id}.consensus.fasta"
    else
        # Normalize
        bcftools norm -f "${ref}" -m -both -o normalized.vcf "${vcf}"
        bgzip -c normalized.vcf > normalized.vcf.gz
        bcftools index normalized.vcf.gz

        # Try generating consensus
        if bcftools consensus -f "${ref}" normalized.vcf.gz -o temp.fasta; then
            sed "s/^>NC_045512.2/>${sample_id}/" temp.fasta > "${sample_id}.consensus.fasta"
        else
            echo "Standard consensus failed. Falling back to SNP-only mode."
            bcftools view -v snps normalized.vcf.gz -Oz -o snps_only.vcf.gz
            bcftools index snps_only.vcf.gz

            if bcftools consensus -f "${ref}" snps_only.vcf.gz -o temp_snps.fasta; then
                 sed "s/^>NC_045512.2/>${sample_id}_SNPsOnly/" temp_snps.fasta > "${sample_id}.consensus.fasta"
            else
                 echo "All methods failed. Using reference."
                 sed "s/^>NC_045512.2/>${sample_id}_FAILED/" "${ref}" > "${sample_id}.consensus.fasta"
            fi
        fi
    fi
    """
}

process pangolinLineage {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'
    input:
      path consensus_files

    output:
      path "results/lineage_report.csv"

    script:
    """
    mkdir -p results
    pangolin ${consensus_files} --threads ${task.cpus} -o results
    """
}

process consensusQC {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'
    input:
      path consensus_files
      path ref

    output:
      path "results/output/"

    script:
    """
    mkdir -p results/output
    # FIX: Removed backslash before consensus_files
    president -r "${ref}" -q "${consensus_files}" -t ${task.cpus} \\
      -a -p results/output/ -f consensus_
    """
}

process phylogeny {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'
    input:
      path consensus_files

    output:
      path "results/phylo.treefile"

    script:
    """
    mkdir -p results

    # Align sequences
    mafft "${consensus_files}" > alignment.fasta
    # Build phylogenetic tree
    iqtree -s alignment.fasta -nt ${task.cpus} -pre results/phylo
    """
}
