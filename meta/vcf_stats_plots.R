install.packages("vcfR")
install.packages("ape")
library(ape)
library(vcfR)

## Author: Abhinav Mishra, FU Berlin 
## A script to visulize, quality control, and writing out the VCF files  
  
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
