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

    //
    // 6) The rest of your pipeline:
    //
    mapping_ch   = mapping(qc_out.cleaned_r1, ref_ch)
    primer_ch    = primerClipping(mapping_ch)
    vcf_ch       = variantCalling(primer_ch, ref_ch)
    qc_masked    = vcfMaskingQC(vcf_ch, ref_ch)
    consensus_ch = consensusGeneration(vcf_ch, ref_ch)
    pangolinLineage(consensus_ch)
    consensusQC(consensus_ch, ref_ch)
    phylogeny(consensus_ch)
}


// … your existing processes (downloadData, referenceGenome, qc, mapping, …) follow unchanged.



/***************************************************************************
 * Download raw FASTQs into data/*.fastq.gz
 ***************************************************************************/
process downloadData {
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
    input:
      path cleaned_r1
      path ref

    output:
      path "datafiles/minimap2-illumina*.sorted.bam"

    script:
    """
    mkdir -p datafiles
    i=1
    for R1 in datafiles/*R1.clean.fastq.gz; do
      R2=\${R1/R1.clean.fastq.gz/R2.clean.fastq.gz}
      minimap2 -x sr -t ${task.cpus} -a -o minimap2-illumina\$i.sam $ref \$R1 \$R2
      samtools view -bS minimap2-illumina\$i.sam | samtools sort -o datafiles/minimap2-illumina\$i.sorted.bam
      samtools index    datafiles/minimap2-illumina\$i.sorted.bam
      i=\$((i+1))
    done
    """
}



/***************************************************************************
 * Primer clipping
 ***************************************************************************/
process primerClipping {
    input:
      path bam_files

    output:
      path "datafiles/*.primerclipped.bam"

    script:
    """
    mkdir -p datafiles
    wget -qO cleanplex.amplicons.bedpe https://osf.io/4nztj/download
    sed 's/NM_003194/NC_045512.2/g' cleanplex.amplicons.bedpe > SARSCoV2.amplicons.bedpe
    for bam in datafiles/*.sorted.bam; do
      bamclipper.sh -b \$bam -p SARSCoV2.amplicons.bedpe -n ${task.cpus}
    done
    """
}



/***************************************************************************
 * Variant calling
 ***************************************************************************/
process variantCalling {
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
