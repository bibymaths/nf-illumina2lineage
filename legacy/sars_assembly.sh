#!/bin/bash

# SARS-CoV-2 genome assembly script

# Step 1: Preliminary Setup

# Install environment manager (mamba=1.4.2 and conda=23.3.1)
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh"
bash Mambaforge-Linux-x86_64.sh
conda update -y conda

# Add necessary channels for package installations
conda config --add channels default
conda config --add channels bioconda
conda config --add channels conda-forge

# Create and activate an environment with required bioinformatics tools
mamba create -y -p ~/envs/projectSARS multiqc fastqc fastp minimap2 samtools bcftools igv pangolin president bamclipper freebayes vcftools vcflib mafft bedtools iqtree jalview gnuplot
mamba activate ~/envs/projectSARS

# Download data and prepare project directory structure
mkdir -p ~/sars-project
cd ~/sars-project

wget --no-check-certificate https://osf.io/qu3bh/download -O illumina-amplicon-capture-wgs.tar.gz
tar -xf illumina-amplicon-capture-wgs.tar.gz
mv illumina-amplicon-capture-wgs/* ./
rm -rf illumina-amplicon-capture-wgs illumina-amplicon-capture-wgs.tar.gz

# Download the reference genome using NCBI Entrez Direct utilities
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
export PATH=${HOME}/edirect:${PATH}
esearch -db nucleotide -query "NC_045512.2" | efetch -format fasta > reference.fasta

# Step 2: Quality Control

# Define sample names for paired-end (PE) sequencing data
ILLUMINA_SAMPLE1='200408_20-04246_A_S1_L000_R1_001.fastq.gz'
ILLUMINA_SAMPLE2='200408_20-04246_A_S1_L000_R2_001.fastq.gz'
ILLUMINA_SAMPLE3='200422_20-04444_A_S1_L000_R1_001.fastq.gz'
ILLUMINA_SAMPLE4='200422_20-04444_A_S1_L000_R2_001.fastq.gz'
ILLUMINA_SAMPLE5='200423_20-04411_A_S1_L000_R1_001.fastq.gz'
ILLUMINA_SAMPLE6='200423_20-04411_A_S1_L000_R2_001.fastq.gz'

# Perform QC analysis and preprocessing for each PE group using fastqc and fastp
for i in {1..3}; do
  fastqc -t 8 ${!ILLUMINA_SAMPLE$((2*i-1))} ${!ILLUMINA_SAMPLE$((2*i))}
  fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --qualified_quality_phred 20 --cut_right --thread 8 \
    --html pair${i}.fastp.html --json pair${i}.fastp.json \
    -i ${!ILLUMINA_SAMPLE$((2*i-1))} -I ${!ILLUMINA_SAMPLE$((2*i))} \
    -o pair${i}.R1.clean.fastq.gz -O pair${i}.R2.clean.fastq.gz
  fastqc -t 8 pair${i}.R{1,2}.clean.fastq.gz
  multiqc .
done

# Step 3: Mapping & Visualization

# Align Illumina PE reads with the reference genome to create SAM files
for i in {1..3}; do
  minimap2 -x sr -t 8 -a -o minimap2-illumina${i}.sam reference.fasta pair${i}.R1.clean.fastq.gz pair${i}.R2.clean.fastq.gz

  # Convert SAM to BAM, sort, index, and prepare for visualization
  samtools view -bS minimap2-illumina${i}.sam > minimap2-illumina${i}.bam
  samtools sort minimap2-illumina${i}.bam > minimap2-illumina${i}.sorted.bam
  samtools index minimap2-illumina${i}.sorted.bam

  # Open IGV viewer to visualize original and sorted BAM files
  igv &
done

# Step 4: Primer Clipping

# Download and preprocess CleanPlex amplicons BEDPE file
wget --no-check-certificate https://osf.io/4nztj/download -O cleanplex.amplicons.bedpe
sed 's/NM_003194/NC_045512.2/g' cleanplex.amplicons.bedpe > SARSCoV2.amplicons.bedpe

# Clip primer sequences for each sorted BAM file
for i in {1..3}; do
  bamclipper.sh -b minimap2-illumina${i}.sorted.bam -p SARSCoV2.amplicons.bedpe -n 8
done

# Step 5: Variant Calling

# Call variants for each BAM file with FreeBayes
for i in {1..3}; do
  freebayes -f reference.fasta --min-alternate-count 10 --min-alternate-fraction 0.1 \
    --min-coverage 20 --pooled-continuous --haplotype-length -1 minimap2-illumina${i}.sorted.primerclipped.bam \
    > freebayes-illumina${i}.vcf

done

# Step 6: Filtering & Masking

# Perform quality control, visualization, and write filtered VCF files using R
Rscript -e "install.packages('vcfR'); install.packages('ape'); library(ape); library(vcfR);
vcf1 <- read.vcfR('/Volumes/EDEN/sars-project/pair1/freebayes-illumina.vcf', verbose = FALSE);
vcf2 <- read.vcfR('/Volumes/EDEN/sars-project/pair2/freebayes-illumina.vcf', verbose = FALSE);
vcf3 <- read.vcfR('/Volumes/EDEN/sars-project/pair3/freebayes-illumina.vcf', verbose = FALSE);
dna <- ape::read.dna('/Volumes/EDEN/sars-project/reference.fasta', format = 'fasta');
chrom1 <- create.chromR(name='Illumina-PE1 | 200408,A,20-04246,CleanPlex SARS-CoV-2,IQ', vcf=vcf1, seq=dna);
chrom2 <- create.chromR(name='Illumina-PE2 | 200422,A,20-04444,Nextera Flex,IQ', vcf=vcf2, seq=dna);
chrom3 <- create.chromR(name='Illumina-PE3 | 200423,A,20-04411,Nextera_XT,NX', vcf=vcf3, seq=dna);
chrom1 <- masker(chrom1);
chrom2 <- masker(chrom2);
chrom3 <- masker(chrom3, min_QUAL = 1);
write.vcf(chrom1, file = '~/sars-project/pair1/masked-strict.vcf', mask = TRUE);
write.vcf(chrom2, file = '~/sars-project/pair2/masked-strict.vcf', mask = TRUE);
write.vcf(chrom3, file = '~/sars-project/pair3/masked-strict.vcf', mask = TRUE);"

# Step 7: Consensus Generation

# Generate consensus sequences for each dataset and combine into a single FASTA file
for i in {1..3}; do
  bcftools view freebayes-illumina${i}.vcf -Oz -o masked-strict${i}.vcf.gz
  bcftools index masked-strict${i}.vcf.gz
  bcftools consensus -f reference.fasta masked-strict${i}.vcf.gz -o consensus-illumina-qc-strict${i}.fasta
  sed -i "s/NC_045512.2/Consensus-Illumina-PE${i}/g" consensus-illumina-qc-strict${i}.fasta
done

# Combine consensus sequences into a single FASTA file
cat consensus-illumina-qc-strict{1,2,3}.fasta > consensus-seqs.fasta

# Step 8: Lineage Annotation

# Assign SARS-CoV-2 lineages using Pangolin
pangolin --update
pangolin --update-data
pangolin -t 8 consensus-seqs.fasta

# Step 9: Consensus QC

# Evaluate consensus sequences using PRESIDENT for pairwise nucleotide identity and ambiguous base reporting
president -r reference.fasta -q consensus-seqs.fasta -t 8 -a -p output/ -f consensus_

# Step 10: Phylogeny & MSA

# Perform multiple sequence alignment and phylogenetic analysis
mafft --thread 8 consensus-seqs.fasta > alignment.fasta
jalview -open alignment.fasta
iqtree -nt 8 -s alignment.fasta --prefix phylo

# Notify user
echo "Pipeline execution completed."
