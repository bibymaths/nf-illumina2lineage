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
    vcf_ch       = variantCalling(primer_ch, ref_ch)
    qc_masked    = vcfMaskingQC(vcf_ch, ref_ch)
    consensus_ch = consensusGeneration(vcf_ch, ref_ch)
    pangolinLineage(consensus_ch)
    consensusQC(consensus_ch, ref_ch)
    phylogeny(consensus_ch)
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
    pigz -dc illumina-amplicon-capture-wgs.tar.gz | tar -x -C .

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
      tuple val(sample_id), path(bam)

    output:
      path "${sample_id}.primerclipped.bam"

    script:
    """
    wget -qO cleanplex.amplicons.bedpe https://osf.io/4nztj/download
    sed 's/NM_003194/NC_045512.2/g' cleanplex.amplicons.bedpe > SARSCoV2.amplicons.bedpe

    bamclipper.sh -b \$bam -p SARSCoV2.amplicons.bedpe -n ${task.cpus} -u 10 -d 10 \\
      -o ${sample_id}.primerclipped.bam
    """
}


process variantCalling {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'
    input:
      path primer_files
      path ref

    output:
      path "datafiles/freebayes-illumina*.vcf"

    script:
    """
    mkdir -p datafiles
    i=1
    for bam in datafiles/*.primerclipped.bam; do
      freebayes -f $ref --min-alternate-count 10 \\
        --min-alternate-fraction 0.1 --min-coverage 20 \\
        --pooled-continuous --haplotype-length -1 \$bam \\
        > datafiles/freebayes-illumina\$i.vcf
      i=\$((i+1))
    done
    """
}


process vcfMaskingQC {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'
    input:
      path vcf_files
      path ref

    output:
      path "datafiles/pair*/masked-strict.vcf"
      path "datafiles/pair*/qc_*.pdf"

    script:
    """
    mkdir -p pair1 pair2 pair3 datafiles
    cp datafiles/freebayes-illumina1.vcf pair1/
    cp datafiles/freebayes-illumina2.vcf pair2/
    cp datafiles/freebayes-illumina3.vcf pair3/
    cp \$ref reference.fasta
    Rscript scripts/vcf_qc_masking.R
    mv pair*/masked-strict.vcf datafiles/pair*/
    mv pair*/*.pdf          datafiles/pair*/
    """
}


process consensusGeneration {

    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    input:
      path vcf_files
      path ref

    output:
      path "results/consensus-seqs.fasta"

    script:
    """
    mkdir -p results
    i=1
    for vcf in datafiles/freebayes-illumina*.vcf; do
      bgzip -c \$vcf > masked-strict\$i.vcf.gz
      bcftools index masked-strict\$i.vcf.gz
      bcftools consensus -f $ref masked-strict\$i.vcf.gz \\
        -o consensus-illumina-qc-strict\$i.fasta
      sed -i "s/NC_045512.2/Consensus-Illumina-PE\$i/g" consensus-illumina-qc-strict\$i.fasta
      i=\$((i+1))
    done
    cat consensus-illumina-qc-strict*.fasta > results/consensus-seqs.fasta
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
    pangolin --update
    pangolin --update-data
    pangolin -t ${task.cpus} \$consensus_files -o results
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
    president -r $ref -q \$consensus_files -t ${task.cpus} \\
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
    mafft --thread ${task.cpus} \$consensus_files > alignment.fasta
    iqtree -nt ${task.cpus} -s alignment.fasta --prefix results/phylo
    """
}
