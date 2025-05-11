#!/usr/bin/env nextflow
nextflow.enable.dsl=2


workflow {

    //
    // 1) Download creates the FASTQs at runtime
    //    It emits _one_ List<File> containing all the .fastq.gz files.
    //
    raw_list_ch = downloadData()

    //
    // 2) Flatten that List<File> into individual File messages
    //
    raw_reads_ch = raw_list_ch
        .flatMap { it }    // now each message is a single File

    //
    // 3) Fetch the reference (independent)
    //
    ref_ch = referenceGenome()

    //
    // 4) Take each File, strip off its _R1_ or _R2_ suffix to get the sample ID,
    //    then group into two-element lists `[R1_file, R2_file]`.
    //
    raw_pairs = raw_reads_ch
        .map { file ->
            // e.g. "sample_X_L001_R1_001.fastq.gz" → "sample_X_L001"
            def id = file.name.replaceFirst(/_R[12]_001\.fastq\.gz$/, '')
            tuple( id, file )
        }
        .groupTuple()   // emits exactly ( id, [r1_file, r2_file] )
   //
   // 4) Turn (id, [r1,r2]) → ( id, r1, r2 )
   //
    reads_ch = raw_pairs
      .map { sample_id, files ->
           // locate R1 versus R2
           def r1 = files.find{ it.name.contains('_R1_') }
           def r2 = files.find{ it.name.contains('_R2_') }
           tuple( sample_id, r1, r2 )
       }
    //
    // 5) QC now sees a real (id, r1, r2) tuple and will run.
    //
    qc_out = qc(reads_ch)

    // 1. Map R1 to (sample_id, R1 file)
    r1_ch = qc_out.cleaned_r1.map { file ->
        def sid = file.getBaseName().replaceFirst(/\\.R1\\.clean$/, '')
        tuple(sid, file)
    }

    // 2. Map R2 to (sample_id, R2 file)
    r2_ch = qc_out.cleaned_r2.map { file ->
        def sid = file.getBaseName().replaceFirst(/\\.R2\\.clean$/, '')
        tuple(sid, file)
    }

    // ✅ 3. Use `.join()` (not `.combine()`) inside workflow
    paired_ch = r1_ch.join(r2_ch)
        .map { a, b ->
            assert a[0] == b[0]
            tuple(a[0], a[1], b[1])
        }

    // 4. Add reference to each tuple
    mapping_input_ch = paired_ch.map { sid, r1, r2 -> tuple(sid, r1, r2, ref_ch) }

    // 5. Run mapping
    mapping_ch = mapping(mapping_input_ch)


    primer_ch    = primerClipping(mapping_ch.bam_files)
    vcf_ch       = variantCalling(primer_ch, ref_ch)
    qc_masked    = vcfMaskingQC(vcf_ch, ref_ch)
    consensus_ch = consensusGeneration(vcf_ch, ref_ch)
    pangolinLineage(consensus_ch)
    consensusQC(consensus_ch, ref_ch)
    phylogeny(consensus_ch)
}

/***************************************************************************
 * Download raw FASTQs into data/*.fastq.gz
 ***************************************************************************/
process downloadData {
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    output:
      path "data/*.fastq.gz", emit: raw_list_ch

    script:
    """
    set -euo pipefail
    mkdir -p data

    # parallel fetch + extract
    aria2c -x16 -s16 -d data -o illumina-amplicon-capture-wgs.tar.gz \\
      https://osf.io/qu3bh/download
    pigz -dc data/illumina-amplicon-capture-wgs.tar.gz | tar -x -C data

    # flatten any nested dirs
    find data/ -name "*.fastq.gz" -exec mv {} data/ \\;
    rm -rf data/illumina-amplicon-capture-wgs* data/__MACOSX || true
    """
}



/***************************************************************************
 * Fetch reference genome
 ***************************************************************************/
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



/***************************************************************************
 * QC + trimming
 ***************************************************************************/
process qc {
    tag "$sample_id"

    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    input:
      tuple val(sample_id), path(r1), path(r2)

    output:
      path "datafiles/${sample_id}.R1.clean.fastq.gz", emit: cleaned_r1
      path "datafiles/${sample_id}.R2.clean.fastq.gz", emit: cleaned_r2
      path "datafiles/${sample_id}.fastp.html",        emit: qc_html
      path "datafiles/${sample_id}.fastp.json",        emit: qc_json

    shell:
    """
    set -euo pipefail
    mkdir -p datafiles

    # raw QC
    fastqc -t ${task.cpus} "$r1" "$r2" -o datafiles

    # trimming & per-sample QC
    fastp --detect_adapter_for_pe \
          --overrepresentation_analysis \
          --correction \
          --qualified_quality_phred 20 \
          --cut_right \
          --thread ${task.cpus} \
          --html  datafiles/${sample_id}.fastp.html \
          --json  datafiles/${sample_id}.fastp.json \
          -i "$r1" -I "$r2" \
          -o datafiles/${sample_id}.R1.clean.fastq.gz \
          -O datafiles/${sample_id}.R2.clean.fastq.gz

    # QC on cleaned
    fastqc -t ${task.cpus} \\
           datafiles/${sample_id}.R1.clean.fastq.gz \\
           datafiles/${sample_id}.R2.clean.fastq.gz \\
           -o datafiles

    # aggregate
    multiqc datafiles -o datafiles
    """
}



/***************************************************************************
 * Mapping
 ***************************************************************************/
process mapping {
    tag "$sample_id"
    publishDir "${params.intermediate}/${task.process}", mode: 'copy'

    input:
      tuple val(sample_id), path(r1), path(r2), path(ref)

    output:
      tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam_files

    script:
    """
    minimap2 -x sr -t ${task.cpus} -a -o ${sample_id}.sam \$ref \$r1 \$r2
    samtools view -bS ${sample_id}.sam | samtools sort -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}

/***************************************************************************
 * Primer clipping
 ***************************************************************************/
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

/***************************************************************************
 * Variant calling
 ***************************************************************************/
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



/***************************************************************************
 * Masking QC
 ***************************************************************************/
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



/***************************************************************************
 * Consensus generation
 ***************************************************************************/
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



/***************************************************************************
 * Lineage assignment
 ***************************************************************************/
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



/***************************************************************************
 * Consensus QC
 ***************************************************************************/
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



/***************************************************************************
 * Phylogenetic analysis
 ***************************************************************************/
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
