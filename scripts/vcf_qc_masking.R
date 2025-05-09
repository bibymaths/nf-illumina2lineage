# Load required packages
if (!requireNamespace("vcfR", quietly = TRUE)) install.packages("vcfR")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")

library(ape)
library(vcfR)

## Author: Abhinav Mishra, FU Berlin
## Visualize, QC, and mask VCF files for SARS-CoV-2

# Load data
vcf1 <- read.vcfR("pair1/freebayes-illumina.vcf", verbose = FALSE)
vcf2 <- read.vcfR("pair2/freebayes-illumina.vcf", verbose = FALSE)
vcf3 <- read.vcfR("pair3/freebayes-illumina.vcf", verbose = FALSE)

dna <- ape::read.dna("reference.fasta", format = "fasta")

# Create chromR objects
chrom1 <- create.chromR(name='Illumina-PE1', vcf=vcf1, seq=dna)
chrom2 <- create.chromR(name='Illumina-PE2', vcf=vcf2, seq=dna)
chrom3 <- create.chromR(name='Illumina-PE3', vcf=vcf3, seq=dna)

# Plot: Before Masking
pdf("pair1/qc_before_masking.pdf")
plot(chrom1, main = "Illumina-PE1", sub = "Before Masking")
dev.off()

pdf("pair2/qc_before_masking.pdf")
plot(chrom2, main = "Illumina-PE2", sub = "Before Masking")
dev.off()

pdf("pair3/qc_before_masking.pdf")
plot(chrom3, main = "Illumina-PE3", sub = "Before Masking")
dev.off()

# Mask low-quality regions
chrom1 <- masker(chrom1)
chrom2 <- masker(chrom2)
chrom3 <- masker(chrom3, min_QUAL = 1)

# Plot: After Masking
pdf("pair1/qc_after_masking.pdf")
plot(chrom1, main = "Illumina-PE1", sub = "After Masking")
dev.off()

pdf("pair2/qc_after_masking.pdf")
plot(chrom2, main = "Illumina-PE2", sub = "After Masking")
dev.off()

pdf("pair3/qc_after_masking.pdf")
plot(chrom3, main = "Illumina-PE3", sub = "After Masking")
dev.off()

# Process for variant window metrics
chrom1 <- proc.chromR(chrom1, verbose = FALSE)
chrom2 <- proc.chromR(chrom2, verbose = FALSE)
chrom3 <- proc.chromR(chrom3, verbose = FALSE)

# Plot: Variant counts
pdf("pair1/qc_variant_counts.pdf")
plot(chrom1, main = "Illumina-PE1", sub = "Variant Counts")
dev.off()

pdf("pair2/qc_variant_counts.pdf")
plot(chrom2, main = "Illumina-PE2", sub = "Variant Counts")
dev.off()

pdf("pair3/qc_variant_counts.pdf")
plot(chrom3, main = "Illumina-PE3", sub = "Variant Counts")
dev.off()

# Composite plots
pdf("pair1/chromoqc_summary.pdf")
chromoqc(chrom1)
dev.off()

pdf("pair2/chromoqc_summary.pdf")
chromoqc(chrom2)
dev.off()

pdf("pair3/chromoqc_summary.pdf")
chromoqc(chrom3)
dev.off()

# Write masked VCFs
write.vcf(chrom1, file = "pair1/masked-strict.vcf", mask = TRUE)
write.vcf(chrom2, file = "pair2/masked-strict.vcf", mask = TRUE)
write.vcf(chrom3, file = "pair3/masked-strict.vcf", mask = TRUE)