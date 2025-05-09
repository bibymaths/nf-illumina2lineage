include {} from './scripts/banner.nf'

params.reads = "*.fastq.gz"
params.reference = "data/reference.fasta"

workflow {
    Channel
        .from(downloadData())
        .set { reads_ch }

    Channel
        .from(referenceGenome())
        .set { ref_ch }

    qc(reads_ch)
    mapping(qc.cleaned_r1, ref_ch)
    primerClipping(mapping.out)
    variantCalling(primerClipping.out, ref_ch)
    vcfMaskingQC(variantCalling.out, ref_ch)
    consensusGeneration(variantCalling.out, ref_ch)
    pangolinLineage(consensusGeneration.out)
    consensusQC(consensusGeneration.out, ref_ch)
    phylogeny(consensusGeneration.out)
}

process downloadData {
    output:
    path "data/illumina*.fastq.gz"

    script:
    """
    mkdir -p data
    wget -O data/illumina-amplicon-capture-wgs.tar.gz https://osf.io/qu3bh/download
    tar -xf data/illumina-amplicon-capture-wgs.tar.gz -C data
    mv data/illumina-amplicon-capture-wgs/* data/
    rm -rf data/illumina-amplicon-capture-wgs*
    """
}

process referenceGenome {
    output:
    path "data/reference.fasta"

    script:
    """
    mkdir -p data
    sh -c "\$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
    export PATH=\$HOME/edirect:\$PATH
    esearch -db nucleotide -query "NC_045512.2" | efetch -format fasta > data/reference.fasta
    """
}

process qc {
    input:
    path reads

    output:
    path "datafiles/pair*.R1.clean.fastq.gz", emit: cleaned_r1
    path "datafiles/pair*.R2.clean.fastq.gz", emit: cleaned_r2
    path "datafiles/*.html", emit: qc_html

    script:
    """
    mkdir -p datafiles
    i=1
    for ((j=1;j<=6;j+=2)); do
      R1=\$(ls *R1_001.fastq.gz | sed -n "\${j}p")
      R2=\$(ls *R2_001.fastq.gz | sed -n "\$((j+1))p")
      fastqc -t 8 \$R1 \$R2
      fastp --detect_adapter_for_pe --correction --qualified_quality_phred 20 --cut_right --thread 8 \
        --html datafiles/pair\$i.fastp.html --json datafiles/pair\$i.fastp.json \
        -i \$R1 -I \$R2 \
        -o datafiles/pair\$i.R1.clean.fastq.gz -O datafiles/pair\$i.R2.clean.fastq.gz
      fastqc -t 8 datafiles/pair\$i.R{1,2}.clean.fastq.gz
      i=\$((i+1))
    done
    multiqc .
    """
}

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
    for R1 in \$(ls datafiles/pair*.R1.clean.fastq.gz); do
      R2=\"\${R1/R1/R2}\"
      minimap2 -x sr -t 8 -a -o minimap2-illumina\$i.sam \$ref \$R1 \$R2
      samtools view -bS minimap2-illumina\$i.sam | samtools sort -o datafiles/minimap2-illumina\$i.sorted.bam
      samtools index datafiles/minimap2-illumina\$i.sorted.bam
      i=\$((i+1))
    done
    """
}

process primerClipping {
    input:
    path bams

    output:
    path "datafiles/*.primerclipped.bam"

    script:
    """
    mkdir -p datafiles
    wget -O cleanplex.amplicons.bedpe https://osf.io/4nztj/download
    sed 's/NM_003194/NC_045512.2/g' cleanplex.amplicons.bedpe > SARSCoV2.amplicons.bedpe
    for bam in datafiles/*.sorted.bam; do
      bamclipper.sh -b \$bam -p SARSCoV2.amplicons.bedpe -n 8
    done
    """
}

process variantCalling {
    input:
    path clipped
    path ref

    output:
    path "datafiles/freebayes-illumina*.vcf"

    script:
    """
    mkdir -p datafiles
    i=1
    for bam in datafiles/*.primerclipped.bam; do
      freebayes -f \$ref --min-alternate-count 10 --min-alternate-fraction 0.1 \
        --min-coverage 20 --pooled-continuous --haplotype-length -1 \$bam > datafiles/freebayes-illumina\$i.vcf
      i=\$((i+1))
    done
    """
}

process vcfMaskingQC {
    input:
    path vcfs
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
    mv pair*/*.pdf datafiles/pair*/
    """
}

process consensusGeneration {
    input:
    path vcfs
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
      bcftools consensus -f \$ref masked-strict\$i.vcf.gz -o consensus-illumina-qc-strict\$i.fasta
      sed -i "s/NC_045512.2/Consensus-Illumina-PE\$i/g" consensus-illumina-qc-strict\$i.fasta
      i=\$((i+1))
    done
    cat consensus-illumina-qc-strict*.fasta > results/consensus-seqs.fasta
    """
}

process pangolinLineage {
    input:
    path consensus

    output:
    path "results/lineage_report.csv"

    script:
    """
    mkdir -p results
    pangolin --update
    pangolin --update-data
    pangolin -t 8 \$consensus -o results
    """
}

process consensusQC {
    input:
    path consensus
    path ref

    output:
    path "results/output/"

    script:
    """
    mkdir -p results/output
    president -r \$ref -q \$consensus -t 8 -a -p results/output/ -f consensus_
    """
}

process phylogeny {
    input:
    path consensus

    output:
    path "results/phylo.treefile"

    script:
    """
    mkdir -p results
    mafft --thread 8 \$consensus > alignment.fasta
    iqtree -nt 8 -s alignment.fasta --prefix results/phylo
    """
}
