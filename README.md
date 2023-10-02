# SARS-CoV-2 genome assembly from Illumina reads 
Course: SARS-2 Bioinformatics & Data Science <br>
Intructors: Max von Kleist, Martin Hölzer <br>
Institution: Freie Universität Berlin, Robert-Koch Institute <br> 

![alt_text](methodSARS.svg) 

### Data & File description  
![alt_text](filedesc.png) 

### Step 1 - Preliminary  

Installing environment manager: _mamba=1.4.2_ and _conda=23.3.1 _ 
 
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh" 
bash Mambaforge-Linux-x86_64.sh 
conda update conda  
```
 
Adding channels
 
```
conda config --add channels default
conda config --add channels bioconda
conda config --add channels conda-forge  
```

Creating environment with tools

```
mamba create -y -p envs/projectSARS  multiqc, fastqc, fastp, minimap2, samtools, bcftools, igv, pangolin, president, bamclipper, freebayes, vcftools, vcflib, mafft, bedtools, iqtree, jalview gnuplot
mamba activate /home/abhinavmishra/envs/projectSARS  
```

Downloading data, and creating three folder for PE data

```
mkdir sars-project
cd sars-project

wget --no-check-certificate https://osf.io/qu3bh/download -O illumina-amplicon-capture-wgs.tar.gz  
tar -xf illumina-amplicon-capture-wgs.tar.gz 

mv /home/abhinavmishra/sars-project/illumina-amplicon-capture-wgs/* /home/abhinavmishra/sars-project/ 
sudo rm -rf illumina-amplicon-capture-wgs 
rm illumina-amplicon-capture-wgs.tar.gz   
```

Download reference genome 

```
## NCBI Entrez Direct UNIX E-utilities 
## download reference genome, and index
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
export PATH=${HOME}/edirect:${PATH}
esearch -db nucleotide -query "NC_045512.2" | efetch -format fasta > reference.fasta 
```
 
### Step 2 - Quality Control

```
## Declaring some environment variables
ILLUMINA_SAMPLE1='200408_20-04246_A_S1_L000_R1_001.fastq.gz'
ILLUMINA_SAMPLE2='200408_20-04246_A_S1_L000_R2_001.fastq.gz' 
ILLUMINA_SAMPLE3='200422_20-04444_A_S1_L000_R1_001.fastq.gz'
ILLUMINA_SAMPLE4='200422_20-04444_A_S1_L000_R2_001.fastq.gz' 
ILLUMINA_SAMPLE5='200423_20-04411_A_S1_L000_R1_001.fastq.gz'
ILLUMINA_SAMPLE6='200423_20-04411_A_S1_L000_R2_001.fastq.gz'   
```

For each PE group/pair, we will do QC analysis and preprocessing using fastqc+fastp 

```
# --cut_right  = 5’->3’ : meanQ <20 (drop bases-stop)
# --correction = overlap correction - base quality, > Q20

## PE 1
fastqc -t 8 $ILLUMINA_SAMPLE1 $ILLUMINA_SAMPLE2 
fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --qualified_quality_phred 20 --cut_right --thread 8 --html pair1.fastp.html --json pair1.fastp.json -i $ILLUMINA_SAMPLE1 -I $ILLUMINA_SAMPLE2 -o pair1.R1.clean.fastq.gz -O pair1.R2.clean.fastq.gz    
fastqc -t 8 pair1.R{1,2}.clean.fastq.gz  
 
## PE 2 
fastqc -t 8 $ILLUMINA_SAMPLE3 $ILLUMINA_SAMPLE4  
fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --qualified_quality_phred 20 --cut_right --thread 8 --html pair2.fastp.html --json pair2.fastp.json -i $ILLUMINA_SAMPLE3 -I $ILLUMINA_SAMPLE4 -o pair2.R1.clean.fastq.gz -O pair2.R2.clean.fastq.gz   
fastqc -t 8 pair2.R{1,2}.clean.fastq.gz    

## PE 3 
fastqc -t 8 $ILLUMINA_SAMPLE5 $ILLUMINA_SAMPLE6  
fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --qualified_quality_phred 20 --cut_right --thread 8 --html pair3.fastp.html --json pair3.fastp.json -i $ILLUMINA_SAMPLE5 -I $ILLUMINA_SAMPLE6 -o pair3.R1.clean.fastq.gz -O pair3.R2.clean.fastq.gz
fastqc -t 8 pair3.R{1,2}.clean.fastq.gz  
   
```

Comparing _fastqc_ and _fastp_ reports 

```
# recrusive search on the html reports in the directory
# modify filenames according to readability to user in the report
multiqc .
```

### Step 3 - Mapping & Visualisation
 
Aligning Illumina PE-reads with reference that runs a sequence alignment program that aligns DNA, creating _.sam _ files  

```
# alignment score >= 80
# 20,000 >= bandwidth >= 500
minimap2 -x sr -t 8 -a -o minimap2-illumina.sam reference.fasta pair1.R1.clean.fastq.gz pair1.R2.clean.fastq.gz  
minimap2 -x sr -t 8 -a -o minimap2-illumina.sam reference.fasta pair2.R1.clean.fastq.gz pair2.R2.clean.fastq.gz 
minimap2 -x sr -t 8 -a -o minimap2-illumina.sam reference.fasta pair3.R1.clean.fastq.gz pair3.R2.clean.fastq.gz  
```

For each pair, convert _.sam_ to _.bam_ and processing mapping reads 

```
for SAM in minimap2-illumina.sam; do
    BN=$(basename $SAM.sam); 
    samtools view -bS $SAM > $BN.bam; 
    samtools sort $BN.bam > $BN.sorted.bam; 
    samtools index $BN.sorted.bam;  
    done  
```
 
The loop above indexes and sorts, then visualize the mapped reads in _igv_ viewer for both: original and sorted _.bam_ files 

```
igv&
```

  
### Step 4 - Primer Clipping   

We checked and related to cleanplex amplicons _.bedpe_ file 

```
wget --no-check-certificate https://osf.io/4nztj/download -O bed_files/cleanplex.amplicons.bedpe  
## Checking the difference of locations and modifying after
bedtools pairtobed -a bed_files/cleanplex.amplicons.bedpe -b bed_files/SARSCoV2.amplicon.bed -type neither 
bedtools pairtobed -a bed_files/cleanplex.amplicons.bedpe -b bed_files/SARSCoV2.ampInsert.bed -type neither 
bedtools pairtobed -a bed_files/cleanplex.amplicons.bedpe -b bed_files/SARSCoV2.ampInsert.shift1.bed -type neither 
```

After matching the _.FASTA_ header and _ID_s 

```
head reference.fasta 
head bed_files/cleanplex.amplicons.bedpe 
sed 's/NM_003194/NC_045512.2/g' bed_files/cleanplex.amplicons.bedpe > bed_files/SARSCoV2.amplicons.bedpe
```

Finally, clipping the primer sequences

```
# For every PE folder 
bamclipper.sh -b minimap2-illumina.sorted.bam -p SARSCoV2.amplicons.bedpe -n 8 
```

### Step 5 - Variant Calling 

_--min-alternate-count_ (N reads in a sample support an allele)
_--min-alternate-fraction_ (N fraction of observations in supporting an allele)
_--min-coverage_, _--pooled-continuous_ (number of samples in the pool)
_--haplotype-length_ (short reads issue - avoid base quality recalibration )

```
# For every PE folder
freebayes -f reference.fasta --min-alternate-count 10 --min-alternate-fraction 0.1 --min-coverage 20 --pooled-continuous --haplotype-length -1 minimap2-illumina.sorted.primerclipped.bam > freebayes-illumina.vcf 
```

PE3 didn't had any variants so a quorum effort was done that gave a lot of unknown, unphased samples using _--report-monomorphic_ (loci which appear to bemonomorphic, and alleles, even those not present in called genotypes)

```
# For PE3
freebayes -f reference.fasta -b minimap2-illumina.sorted.primerclipped.bam --report-monomorphic --pooled-continuous > freebayes-illumina.vcf
```

### Step 6 - Filtering & Masking

A small _R_ Script using _vcfR_ package was used for quality control, visualization, and writing out the new _.vcf_ file. It does annotation for each PE data for analysis and composite plots, with masking low-coverage region based on quality and depth. 

```r
install.packages("vcfR")
install.packages("ape")
library(ape)
library(vcfR)

## Author: Abhinav Mishra, FU Berlin 
## A script to visulize, quality control, and writing out the VCF files
## After the masking, what happens 
## *../pair{1,2,3}/freebayes-illumina.vcf => *../pair{1,2,3}/masked-strict.vcf

  
## Load VCF data 

vcf1 <- read.vcfR("/Volumes/EDEN/sars-project/pair1/freebayes-illumina.vcf", verbose = FALSE)
vcf2 <- read.vcfR("/Volumes/EDEN/sars-project/pair2/freebayes-illumina.vcf", verbose = FALSE)
vcf3 <- read.vcfR("/Volumes/EDEN/sars-project/pair3/freebayes-illumina.vcf", verbose = FALSE)
 
## Load reference genome 

dna <- ape::read.dna("/Volumes/EDEN/sars-project/reference.fasta", format = "fasta") 
 
## Create chromR object for analysis, with proper annotation  
## for the datsets, PE : Paired end
chrom1 <- create.chromR(name='Illumina-PE1 | 200408,A,20-04246,CleanPlex SARS-CoV-2,IQ', vcf=vcf1, seq=dna) 
chrom2 <- create.chromR(name='Illumina-PE2 | 200422,A,20-04444,Nextera Flex,IQ', vcf=vcf2, seq=dna) 
chrom3 <- create.chromR(name='Illumina-PE3 | 200423,A,20-04411,Nextera_XT,NX', vcf=vcf3, seq=dna) 
 
## Quantitative plots
plot(chrom1, main = "Illumina-PE1", sub = "Before Masking")  
plot(chrom2, main = "Illumina-PE2", sub = "Before Masking") 
plot(chrom3, main = "Illumina-PE3", sub = "Before Masking") 

## Masking out data that we do not have high confidence in 
## based on quality and depth
 
chrom1 <-masker(chrom1)  
chrom2 <-masker(chrom2) 
chrom3 <-masker(chrom3, min_QUAL = 1) 
 
## Check again 

plot(chrom1, main = "Illumina-PE1", sub = "After Masking") 
plot(chrom2, main = "Illumina-PE2", sub = "After Masking") 
plot(chrom3, main = "Illumina-PE3", sub = "After Masking") 
 
## Adding variant counts per window 

chrom1 <- proc.chromR(chrom1, verbose = FALSE) 
chrom2 <- proc.chromR(chrom2, verbose = FALSE) 
chrom3 <- proc.chromR(chrom3, verbose = FALSE) 
 
## Check one more time 

plot(chrom1, main = "Illumina-PE1", sub = "Adding Variant counts")  
plot(chrom2, main = "Illumina-PE2", sub = "Adding Variant counts")
plot(chrom3, main = "Illumina-PE3", sub = "Adding Variant counts")
 
## Composite plots-summary  
## PE3 didn't had any variants   

chromoqc(chrom1)  
chromoqc(chrom2) 
chromoqc(chrom3) 
 
## Writing out new VCF files with QC checks,  
## filtering for low-coverage, masking them    

write.vcf(chrom1, file = "~/sars-project/pair1/masked-strict.vcf", mask = TRUE) 
write.vcf(chrom2, file = "~/sars-project/pair2/masked-strict.vcf", mask = TRUE) 
write.vcf(chrom2, file = "~/sars-project/pair3/masked-strict.vcf", mask = TRUE)
```

For generating _.vcf_ file statistics for interpretation and reassessment 

```
## mean depth per individual
vcftools --gzvcf masked-strict.vcf.gz --depth --out filterVCF/stats   
## mean depth per site
vcftools --gzvcf masked-strict.vcf.gz --site-mean-depth --out filterVCF/stats   
## proportion of missing data per site 
vcftools --gzvcf masked-strict.vcf.gz --missing-site --out filterVCF/stats   
## site quality
vcftools --gzvcf masked-strict.vcf.gz --site-quality --out filterVCF/stats 

vcf-compare pair1/masked.vcf.gz pair2/masked.vcf.gz pair3/masked.vcf.gz 

##  Venn diagrams -  numbers of positions contained in one but not the other files; two but not the other files, etc. 
vcf-compare pair1/masked.vcf.gz pair2/masked.vcf.gz pair3/masked.vcf.gz | grep ^VN | cut -f 2- > filterVCF/fullstattree.txt 
## returns a tree with basic stats e.g. SNPs, INDELS etc.
vcf-stats pair1/masked.vcf.gz >> filterVCF/fullstattree.txt 
## Append TSTV ratio to stats file - all 
cat pair1/masked.vcf | vcf-tstv >> pair1/filterVCF/fullstatstree.txt
```

For generating some plots using a perl [script](https://github.com/ENCODE-DCC/chip-seq-pipeline/blob/master/dnanexus/bam2tagalign/resources/usr/local/bin/samtools-1.0/bin/plot-bamstats) from [Petr Danecek](https://www.sanger.ac.uk/person/danecek-petr/) 

```
#For every PE
samtools stats minimap2-illumina.sorted.primerclipped.bam > bamstats.bam.bc
perl bamstats.pl -p gnuplots/ bamstats.bam.bc
```

### Step 7 - Consensus  

After making the variant calls, normalising indels and filtering low-coverage

```
#PE1
bcftools view masked-strict.vcf -Oz -o masked-strict.vcf.gz 
bcftools index masked-strict.vcf.gz 
bcftools consensus -f reference.fasta masked-strict.vcf.gz -o consensus-illumina-qc-strict.fasta 
sed -i 's/NC_045512.2/Consensus-Illumina-PE1 | 200408,A,20-04246,CleanPlex SARS-CoV-2,IQ/g' consensus-illumina-qc-strict.fasta  
 
#PE2
bcftools view masked-strict.vcf -Oz -o masked-strict.vcf.gz 
bcftools index masked-strict.vcf.gz 
bcftools consensus -f reference.fasta masked-strict.vcf.gz -o consensus-illumina-qc-strict.fasta 
sed -i 's/NC_045512.2/Consensus-Illumina-PE2 | 200422,A,20-04444,Nextera Flex,IQ/g' consensus-illumina-qc-strict.fasta  
 
#PE3
bcftools view masked-strict.vcf -Oz -o masked-strict.vcf.gz 
bcftools index masked-strict.vcf.gz 
bcftools consensus -f reference.fasta masked-strict.vcf.gz -o consensus-illumina-qc-strict.fasta 
sed -i 's/NC_045512.2/Consensus-Illumina-PE3 | 200423,A,20-04411,Nextera_XT,NX/g' consensus-illumina-qc-strict.fasta 
```

Put all three consensus from PE datsets into one _.fasta_ file 

### Step 8 - Lineage Annotation 

For the dynamic nomenclature of SARS-CoV-2 lineages a.k.a. pangloin nomenclature, and assigning most likely lineage (Pango lineage) to query sequences

```
pangolin --update 
pangolin --update-data 
pangolin -t 8 consensus-seqs.fasta 
```

### Step 9 - Consensus QC 

Based on calculating pairwise nucleotide identity and reporting ambiguous _N_'s

```
# use -a for storing the alignment
president -r reference.fasta -q consensus-seqs.fasta -t 8 -a -p output/ -f consensus_ 
```

### Step 10 - Phylogeny & MSA

For checking the positions of the mismatches, and distance matrix between each and every query and reference sequences

```
mafft --thread 8 consensus-seqs.fasta > alignment.fasta
jalview -open alignment.fasta 
iqtree -nt 8 -s alignment.fasta --prefix phylo 
```

The .treefile is visuzalized in [IROKI](https://www.iroki.net/viewer).
