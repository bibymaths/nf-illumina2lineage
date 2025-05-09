// main.nf

params.reads = "*.fastq.gz"
params.reference = "reference.fasta"

workflow {
    downloadData()
    referenceGenome()
    qc()
    mapping()
    primerClipping()
    variantCalling()
    consensusGeneration()
    pangolinLineage()
    consensusQC()
    phylogeny()
}

process downloadData {
    output:
    path "illumina*.fastq.gz"

    script:
    """
    wget -O illumina-amplicon-capture-wgs.tar.gz https://osf.io/qu3bh/download
    tar -xf illumina-amplicon-capture-wgs.tar.gz
    mv illumina-amplicon-capture-wgs/* ./
    rm -rf illumina-amplicon-capture-wgs*
    """
}

process referenceGenome {
    output:
    path "reference.fasta"

    script:
    """
    sh -c \"\$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)\"
    export PATH=\$HOME/edirect:\$PATH
    esearch -db nucleotide -query \"NC_045512.2\" | efetch -format fasta > reference.fasta
    """
}

process qc {
    input:
    path reads from downloadData.out

    output:
    path "pair*.R1.clean.fastq.gz"
    path "pair*.R2.clean.fastq.gz"
    path "*.html", emit: qc_html

    script:
    """
    i=1
    for ((j=1;j<=6;j+=2)); do
      R1=\$(ls *R1_001.fastq.gz | sed -n \"\${j}p\")
      R2=\$(ls *R2_001.fastq.gz | sed -n \"\$((j+1))p\")
      fastqc -t 8 \$R1 \$R2
      fastp --detect_adapter_for_pe --correction --qualified_quality_phred 20 --cut_right --thread 8 \\
        --html pair\$i.fastp.html --json pair\$i.fastp.json \\
        -i \$R1 -I \$R2 \\
        -o pair\$i.R1.clean.fastq.gz -O pair\$i.R2.clean.fastq.gz
      fastqc -t 8 pair\$i.R{1,2}.clean.fastq.gz
      i=\$((i+1))
    done
    multiqc .
    """
}

process mapping {
    input:
    path cleaned from qc.out
    path ref from referenceGenome.out

    output:
    path "minimap2-illumina*.sorted.bam"

    script:
    """
    i=1
    for R1 in pair*.R1.clean.fastq.gz; do
      R2=\"\${R1/R1/R2}\"
      minimap2 -x sr -t 8 -a -o minimap2-illumina\$i.sam \$ref \$R1 \$R2
      samtools view -bS minimap2-illumina\$i.sam | samtools sort -o minimap2-illumina\$i.sorted.bam
      samtools index minimap2-illumina\$i.sorted.bam
      i=\$((i+1))
    done
    """
}

process primerClipping {
    input:
    path bams from mapping.out

    output:
    path "*.primerclipped.bam"

    script:
    """
    wget -O cleanplex.amplicons.bedpe https://osf.io/4nztj/download
    sed 's/NM_003194/NC_045512.2/g' cleanplex.amplicons.bedpe > SARSCoV2.amplicons.bedpe
    for bam in *.sorted.bam; do
      bamclipper.sh -b \$bam -p SARSCoV2.amplicons.bedpe -n 8
    done
    """
}

process variantCalling {
    input:
    path clipped from primerClipping.out
    path ref from referenceGenome.out

    output:
    path "freebayes-illumina*.vcf"

    script:
    """
    i=1
    for bam in *.primerclipped.bam; do
      freebayes -f \$ref --min-alternate-count 10 --min-alternate-fraction 0.1 \\
        --min-coverage 20 --pooled-continuous --haplotype-length -1 \$bam > freebayes-illumina\$i.vcf
      i=\$((i+1))
    done
    """
}

process consensusGeneration {
    input:
    path vcfs from variantCalling.out
    path ref from referenceGenome.out

    output:
    path "consensus-seqs.fasta"

    script:
    """
    i=1
    for vcf in freebayes-illumina*.vcf; do
      bgzip -c \$vcf > masked-strict\$i.vcf.gz
      bcftools index masked-strict\$i.vcf.gz
      bcftools consensus -f \$ref masked-strict\$i.vcf.gz -o consensus-illumina-qc-strict\$i.fasta
      sed -i \"s/NC_045512.2/Consensus-Illumina-PE\$i/g\" consensus-illumina-qc-strict\$i.fasta
      i=\$((i+1))
    done
    cat consensus-illumina-qc-strict*.fasta > consensus-seqs.fasta
    """
}

process pangolinLineage {
    input:
    path consensus from consensusGeneration.out

    script:
    """
    pangolin --update
    pangolin --update-data
    pangolin -t 8 \$consensus
    """
}

process consensusQC {
    input:
    path consensus from consensusGeneration.out
    path ref from referenceGenome.out

    output:
    path "output/"

    script:
    """
    president -r \$ref -q \$consensus -t 8 -a -p output/ -f consensus_
    """
}

process phylogeny {
    input:
    path consensus from consensusGeneration.out

    script:
    """
    mafft --thread 8 \$consensus > alignment.fasta
    iqtree -nt 8 -s alignment.fasta --prefix phylo
    """
}
