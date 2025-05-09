// include { banner } from './scripts/banner.nf'

params.reads     = "*.fastq.gz"
params.reference = "data/reference.fasta"

workflow {
    // 1) Download data and reference
    downloadData()
    ref_ch = referenceGenome()

    // 2) Pair FASTQs *after* they exist
    reads_ch = Channel
        .fromFilePairs( 'data/*_R{1,2}_001.fastq.gz', flat: true )
        .map { sample, files -> tuple(sample, files) }

    // 3) QC — now reads_ch is non‐empty
    qc_out = qc(reads_ch)

    // 4) Downstream steps
    mapping_ch  = mapping(qc_out.cleaned_r1, ref_ch)
    primer_ch   = primerClipping(mapping_ch)
    vcf_ch      = variantCalling(primer_ch, ref_ch)
    qc_masked   = vcfMaskingQC(vcf_ch, ref_ch)
    consensus_ch = consensusGeneration(vcf_ch, ref_ch)
    pangolinLineage(consensus_ch)
    consensusQC(consensus_ch, ref_ch)
    phylogeny(consensus_ch)
}


// Download raw data
process downloadData {
    output:
    path "data/*.fastq.gz"

    // Require aria2c and pigz to be installed in your env
    script:
    """
    #!/usr/bin/env bash
    set -e

    mkdir -p data

    # 1) Download with aria2c (16 connections)
    aria2c \\
      -x16 \\
      -s16 \\
      -d data \\
      -o illumina-amplicon-capture-wgs.tar.gz \\
      https://osf.io/qu3bh/download

    # 2) Extract in parallel with pigz
    #    Use pigz to decompress .gz input and tar to extract
    pigz -dc data/illumina-amplicon-capture-wgs.tar.gz | tar -x -C data

    # 3) Move FASTQs to top-level data/
    find data/ -name "*.fastq.gz" -exec mv {} data/ \\;

    # 4) Clean up
    rm -rf data/illumina-amplicon-capture-wgs* data/__MACOSX || true
    """
}


// Fetch reference genome
process referenceGenome {

    shell:
    'bash'

    output:
    path "data/reference.fasta"

    script:
    '''
    mkdir -p data

    wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O - | bash
    export PATH=$HOME/edirect:$PATH

    esearch -db nucleotide -query "NC_045512.2" | efetch -format fasta > data/reference.fasta
    '''
}


// Quality control and trimming
process qc {
    tag "$sample_id"
    input:
    tuple val(sample_id), path(reads)   // reads is a list of two files

    output:
    path "datafiles/${sample_id}.R1.clean.fastq.gz", emit: cleaned_r1
    path "datafiles/${sample_id}.R2.clean.fastq.gz", emit: cleaned_r2
    path "datafiles/${sample_id}.fastp.html",      emit: qc_html

    script:
    """
    mkdir -p datafiles

    # Assign files
    R1=${reads[0]}
    R2=${reads[1]}

    # QC
    fastqc -t 8 \$R1 \$R2

    fastp --detect_adapter_for_pe \
          --correction \
          --qualified_quality_phred 20 \
          --cut_right \
          --thread 8 \
          --html  datafiles/${sample_id}.fastp.html \
          --json  datafiles/${sample_id}.fastp.json \
          -i \$R1 -I \$R2 \
          -o  datafiles/${sample_id}.R1.clean.fastq.gz \
          -O  datafiles/${sample_id}.R2.clean.fastq.gz

    fastqc -t 8 datafiles/${sample_id}.R1.clean.fastq.gz datafiles/${sample_id}.R2.clean.fastq.gz

    multiqc -o datafiles .

    """
}


// Align reads to reference
process mapping {
    input:
    path cleaned_r1
    path ref

    output:
    path "datafiles/minimap2-illumina*.sorted.bam"

    script:
    '''
    mkdir -p datafiles
    i=1
    for R1 in $(ls datafiles/pair*.R1.clean.fastq.gz); do
      R2="${R1/R1.clean.fastq.gz/R2.clean.fastq.gz}"
      minimap2 -x sr -t 8 -a -o minimap2-illumina$i.sam $ref $R1 $R2
      samtools view -bS minimap2-illumina$i.sam | samtools sort -o datafiles/minimap2-illumina$i.sorted.bam
      samtools index datafiles/minimap2-illumina$i.sorted.bam
      i=$((i+1))
    done
    '''
}

// Clip primers
process primerClipping {
    input:
    path mapping_ch

    output:
    path "datafiles/*.primerclipped.bam"

    script:
    """
    mkdir -p datafiles
    wget -O cleanplex.amplicons.bedpe https://osf.io/4nztj/download
    sed 's/NM_003194/NC_045512.2/g' cleanplex.amplicons.bedpe > SARSCoV2.amplicons.bedpe
    for bam in datafiles/*.sorted.bam; do
      bamclipper.sh -b $bam -p SARSCoV2.amplicons.bedpe -n 8
    done
    """
}

// Call variants
process variantCalling {
    input:
    path primer_ch
    path ref

    output:
    path "datafiles/freebayes-illumina*.vcf"

    script:
    '''
    mkdir -p datafiles
    i=1
    for bam in datafiles/*.primerclipped.bam; do
      freebayes -f $ref --min-alternate-count 10 --min-alternate-fraction 0.1 \
        --min-coverage 20 --pooled-continuous --haplotype-length -1 $bam > datafiles/freebayes-illumina$i.vcf
      i=$((i+1))
    done
    '''
}

// QC & mask VCF
process vcfMaskingQC {
    input:
    path vcf_ch
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
    cp $ref reference.fasta
    Rscript scripts/vcf_qc_masking.R
    mv pair*/masked-strict.vcf datafiles/pair*/
    mv pair*/*.pdf datafiles/pair*/
    """
}

// Build consensus
process consensusGeneration {
    input:
    path vcf_ch
    path ref

    output:
    path "results/consensus-seqs.fasta"

    script:
    '''
    mkdir -p results
    i=1
    for vcf in datafiles/freebayes-illumina*.vcf; do
      bgzip -c $vcf > masked-strict$i.vcf.gz
      bcftools index masked-strict$i.vcf.gz
      bcftools consensus -f $ref masked-strict$i.vcf.gz -o consensus-illumina-qc-strict$i.fasta
      sed -i "s/NC_045512.2/Consensus-Illumina-PE$i/g" consensus-illumina-qc-strict$i.fasta
      i=$((i+1))
    done
    cat consensus-illumina-qc-strict*.fasta > results/consensus-seqs.fasta
    '''
}

// Assign lineage
process pangolinLineage {
    input:
    path consensus_ch

    output:
    path "results/lineage_report.csv"

    script:
    """
    mkdir -p results
    pangolin --update
    pangolin --update-data
    pangolin -t 8 $consensus_ch -o results
    """
}

// QC consensus
process consensusQC {
    input:
    path consensus_ch
    path ref

    output:
    path "results/output/"

    script:
    """
    mkdir -p results/output
    president -r $ref -q $consensus_ch -t 8 -a -p results/output/ -f consensus_
    """
}

// Phylogenetic analysis
process phylogeny {
    input:
    path consensus_ch

    output:
    path "results/phylo.treefile"

    script:
    """
    mkdir -p results
    mafft --thread 8 $consensus_ch > alignment.fasta
    iqtree -nt 8 -s alignment.fasta --prefix results/phylo
    """
}
