# SARS-CoV-2 genome assembly from Illumina reads 
Course: SARS-2 Bioinformatics & Data Science <br>
Intructors: Max von Kleist, Martin Hölzer <br>
Institution: Freie Universität Berlin, Robert-Koch Institute <br> 

![alt_text](methodSARS.svg) 

### Data & File description  
![alt_text](filedesc.png) 

### Step 1 - Preliminary  

Installing environment manager: _mamba=1.4.2_ and _conda=23.3.1 _ 
 
```shell
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh" 
shell Mambaforge-Linux-x86_64.sh 
conda update conda  
```
 
Adding channels
 
```shell
conda config --add channels default
conda config --add channels bioconda
conda config --add channels conda-forge  
```

Creating environment with tools

```shell
mamba create -y -p envs/projectSARS  multiqc, fastqc, fastp, minimap2, samtools, bcftools, igv, pangolin, president, bamclipper, freebayes, vcftools, vcflib, mafft, bedtools, iqtree, jalview gnuplot
mamba activate /home/abhinavmishra/envs/projectSARS  
```

Downloading data, and creating three folder for PE data

```shell
mkdir sars-project
cd sars-project

wget --no-check-certificate https://osf.io/qu3bh/download -O illumina-amplicon-capture-wgs.tar.gz  
tar -xf illumina-amplicon-capture-wgs.tar.gz 

mv /home/abhinavmishra/sars-project/illumina-amplicon-capture-wgs/* /home/abhinavmishra/sars-project/ 
sudo rm -rf illumina-amplicon-capture-wgs 
rm illumina-amplicon-capture-wgs.tar.gz   
```

Download reference genome 

```shell
## NCBI Entrez Direct UNIX E-utilities 
## download reference genome, and index
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
export PATH=${HOME}/edirect:${PATH}
esearch -db nucleotide -query "NC_045512.2" | efetch -format fasta > reference.fasta 
```
 
### Step 2 - Quality Control

```shell
## Declaring some environment variables
ILLUMINA_SAMPLE1='200408_20-04246_A_S1_L000_R1_001.fastq.gz'
ILLUMINA_SAMPLE2='200408_20-04246_A_S1_L000_R2_001.fastq.gz' 
ILLUMINA_SAMPLE3='200422_20-04444_A_S1_L000_R1_001.fastq.gz'
ILLUMINA_SAMPLE4='200422_20-04444_A_S1_L000_R2_001.fastq.gz' 
ILLUMINA_SAMPLE5='200423_20-04411_A_S1_L000_R1_001.fastq.gz'
ILLUMINA_SAMPLE6='200423_20-04411_A_S1_L000_R2_001.fastq.gz'   
```

For each PE group/pair, we will do QC analysis and preprocessing using fastqc+fastp 

```shell
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

```shell
# recrusive search on the html reports in the directory
# modify filenames according to readability to user in the report
multiqc .
```

### Step 3 - Mapping & Visualisation
 
Aligning Illumina PE-reads with reference that runs a sequence alignment program that aligns DNA, creating _.sam _ files  

```shell
# alignment score >= 80
# 20,000 >= bandwidth >= 500
minimap2 -x sr -t 8 -a -o minimap2-illumina.sam reference.fasta pair1.R1.clean.fastq.gz pair1.R2.clean.fastq.gz  
minimap2 -x sr -t 8 -a -o minimap2-illumina.sam reference.fasta pair2.R1.clean.fastq.gz pair2.R2.clean.fastq.gz 
minimap2 -x sr -t 8 -a -o minimap2-illumina.sam reference.fasta pair3.R1.clean.fastq.gz pair3.R2.clean.fastq.gz  
```

For each pair, convert _.sam_ to _.bam_ and processing mapping reads 

```shell
for SAM in minimap2-illumina.sam; do
    BN=$(basename $SAM.sam); 
    samtools view -bS $SAM > $BN.bam; 
    samtools sort $BN.bam > $BN.sorted.bam; 
    samtools index $BN.sorted.bam;  
    done  
```
 
The loop above indexes and sorts, then visualize the mapped reads in _igv_ viewer for both: original and sorted _.bam_ files 

```shell
igv&
```

  
### Step 4 - Primer Clipping   

We checked and related to cleanplex amplicons _.bedpe_ file 

```shell
wget --no-check-certificate https://osf.io/4nztj/download -O bed_files/cleanplex.amplicons.bedpe  
## Checking the difference of locations and modifying after
bedtools pairtobed -a bed_files/cleanplex.amplicons.bedpe -b bed_files/SARSCoV2.amplicon.bed -type neither 
bedtools pairtobed -a bed_files/cleanplex.amplicons.bedpe -b bed_files/SARSCoV2.ampInsert.bed -type neither 
bedtools pairtobed -a bed_files/cleanplex.amplicons.bedpe -b bed_files/SARSCoV2.ampInsert.shift1.bed -type neither 
```

After matching the _.FASTA_ header and _ID_s 

```shell
head reference.fasta 
head bed_files/cleanplex.amplicons.bedpe 
sed 's/NM_003194/NC_045512.2/g' bed_files/cleanplex.amplicons.bedpe > bed_files/SARSCoV2.amplicons.bedpe
```

Finally, clipping the primer sequences

```shell
# For every PE folder 
bamclipper.sh -b minimap2-illumina.sorted.bam -p SARSCoV2.amplicons.bedpe -n 8 
```

### Step 5 - Variant Calling 

_--min-alternate-count_ (N reads in a sample support an allele)
_--min-alternate-fraction_ (N fraction of observations in supporting an allele)
_--min-coverage_, _--pooled-continuous_ (number of samples in the pool)
_--haplotype-length_ (short reads issue - avoid base quality recalibration )

```shell
# For every PE folder
freebayes -f reference.fasta --min-alternate-count 10 --min-alternate-fraction 0.1 --min-coverage 20 --pooled-continuous --haplotype-length -1 minimap2-illumina.sorted.primerclipped.bam > freebayes-illumina.vcf 
```

PE3 didn't had any variants so a quorum effort was done that gave a lot of unknown, unphased samples using _--report-monomorphic_ (loci which appear to bemonomorphic, and alleles, even those not present in called genotypes)

```shell
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

```shell
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

```shell
#For every PE
samtools stats minimap2-illumina.sorted.primerclipped.bam > bamstats.bam.bc
perl bamstats.pl -p gnuplots/ bamstats.bam.bc
```

### Step 7 - Consensus  

After making the variant calls, normalising indels and filtering low-coverage

```shell
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

```shell
pangolin --update 
pangolin --update-data 
pangolin -t 8 consensus-seqs.fasta 
```

### Step 9 - Consensus QC 

Based on calculating pairwise nucleotide identity and reporting ambiguous _N_'s

```shell
# use -a for storing the alignment
president -r reference.fasta -q consensus-seqs.fasta -t 8 -a -p output/ -f consensus_ 
```

### Step 10 - Phylogeny & MSA

For checking the positions of the mismatches, and distance matrix between each and every query and reference sequences

```shell
mafft --thread 8 consensus-seqs.fasta > alignment.fasta
jalview -open alignment.fasta 
iqtree -nt 8 -s alignment.fasta --prefix phylo 
```

The .treefile is visuzalized in [IROKI](https://www.iroki.net/viewer). 

## References 

1. bamstats.pl script, 2012-2014 Genome Research Ltd (Author: Petr Danecek <pd3@sanger.ac.uk>)
2. Anon, 2020. Anaconda Software Distribution, Anaconda Inc. Available at: https://docs.anaconda.com/.
3. Philip Ewels, Måns Magnusson, Sverker Lundin, Max Käller, MultiQC: summarize analysis results for
multiple tools and samples in a single report, Bioinformatics, Volume 32, Issue 19, October 2016, Pages
3047–3048, https://doi.org/10.1093/bioinformatics/btw354
4. Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online].
5. Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu, fastp: an ultra-fast all-in-one FASTQ preprocessor,
Bioinformatics, Volume 34, Issue 17, September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
6. Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34,
Issue 18, September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191
7. Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O
Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li,
Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021,
giab008, https://doi.org/10.1093/gigascience/giab008
8. Robinson, J. T., Thorvaldsdóttir, H., Winckler, W., Guttman, M., Lander, E. S., Getz, G., &
Mesirov, J. P. (2011). Integrative genomics viewer. Nature biotechnology, 29(1), 24–26. https://doi.org/10.1038/nbt.1754
9. O'Toole Á, Hill V, Pybus OG et al. Tracking the international spread of SARS-CoV-2 lineages B.1.1.7 and
B.1.351/501Y-V2 [version 1; peer review: 3 approved]. Wellcome Open Res 2021, 6:121 (https://doi.org/10.12688/wellcomeopenres.16661.1)
10. Áine O’Toole, Emily Scher, Anthony Underwood, Ben Jackson, Verity Hill, John T McCrone, Rachel
Colquhoun, Chris Ruis, Khalil Abu-Dahab, Ben Taylor, Corin Yeats, Louis du Plessis, Daniel Maloney, Nathan
Medd, Stephen W Attwood, David M Aanensen, Edward C Holmes, Oliver G Pybus, Andrew Rambaut,
Assignment of epidemiological lineages in an emerging pandemic using the pangolin tool, Virus Evolution,
Volume 7, Issue 2, December 2021, veab064, https://doi.org/10.1093/ve/veab064
11. Rambaut, A., Holmes, E.C., O’Toole, Á. et al. A dynamic nomenclature proposal for SARS-CoV-2 lineages to
assist genomic epidemiology. Nat Microbiol 5, 1403–1407 (2020). https://doi.org/10.1038/s41564-020-0770-5
12. Tool for QC with consensus sequences https://github.com/rki-mf1/president
13. Au, C., Ho, D., Kwong, A. et al. BAMClipper: removing primers from alignments to
minimize false-negative mutations in amplicon next-generation sequencing. Sci Rep 7, 1567
(2017). https://doi.org/10.1038/s41598-017-01703-6
14. Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing.
arXiv preprint arXiv:1207.3907 [q-bio.GN] 2012
15. Vcflib and tools for processing the VCF variant call format. Erik Garrison, Zev N.
Kronenberg, Eric T. Dawson, Brent S. Pedersen, Pjotr Prins. bioRxiv 2021.05.21.445151; doi:
https://doi.org/10.1101/2021.05.21.445151
16. Petr Danecek, Adam Auton, Goncalo Abecasis, Cornelis A. Albers, Eric Banks, Mark A. DePristo, Robert E. Handsaker,
Gerton Lunter, Gabor T. Marth, Stephen T. Sherry, Gilean McVean, Richard Durbin, 1000 Genomes Project Analysis
Group, The variant call format and VCFtools, Bioinformatics, Volume 27, Issue 15, August 2011, Pages 2156–2158, https://doi.org/10.1093/bioinformatics/btr330
17. Waterhouse, A.M., Procter, J.B., Martin, D.M.A, Clamp, M., Barton, G.J (2009), "Jalview version 2: A Multiple Sequence
Alignment and Analysis Workbench,” Bioinformatics 25 (9) 1189-1191 doi: 10.1093/bioinformatics/btp033
18. Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2015). IQ-TREE: a fast and effective stochastic algorithm
for estimating maximum-likelihood phylogenies. Molecular biology and evolution, 32(1), 268–274. https://doi.org/10.1093/molbev/msu300
19. Bui Quang Minh, Heiko A Schmidt, Olga Chernomor, Dominik Schrempf, Michael D Woodhams, Arndt von Haeseler,
Robert Lanfear, IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era, Molecular
Biology and Evolution, Volume 37, Issue 5, May 2020, Pages 1530–1534, https://doi.org/10.1093/molbev/msaa015
20. Aaron R. Quinlan, Ira M. Hall, BEDTools: a flexible suite of utilities for comparing genomic features,
Bioinformatics, Volume 26, Issue 6, March 2010, Pages 841–842, https://doi.org/10.1093/bioinformatics/btq033
21. Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in
performance and usability. Molecular biology and evolution, 30(4), 772–780. https://doi.org/10.1093/molbev/mst010
22. R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical
Computing, Vienna, Austria. URL https://www.R-project.org/.
23. Knaus, B.J. and Grünwald, N.J. (2017), vcfr: a package to manipulate and visualize variant call format data in R.
Mol Ecol Resour, 17: 44-53. https://doi.org/10.1111/1755-0998.12549
24. National Center for Biotechnology Information (NCBI)[Internet]. Bethesda (MD): National Library of
Medicine (US), National Center for Biotechnology Information; [1988] – [cited 2023 Sep 29]. Available from:
https://www.ncbi.nlm.nih.gov/

### Device Info 

The analysis and results was done and generated on 

**OS**          Fedora Linux 38 <br>
**Kernel**      Linux 6.4.15-200.fc38.x86_64 <br>
**Processor**   Intel i5-8250U (8 slots), with CUDA support <br>
**Graphics**    UHD 620 (KBL GT2) <br>
**Memory**      8 GB 
