# SARS-CoV-2 genome assembly from Illumina reads 
Course: [SARS-2 Bioinformatics & Data Science](https://github.com/rki-mf1/2023-SC2-Data-Science) <br>
Intructors: Max von Kleist, Martin Hölzer <br>
Institution: Freie Universität Berlin, Robert-Koch Institute <br> 

![alt_text](methodSARS.svg) 

### Data & File description  
![alt_text](filedesc.png) 

### Step 1 - Preliminary  

Installing environment manager: _mamba=1.4.2_ and _conda=23.3.1 _ 

Adding channels
 
Creating environment with tools

Downloading data, and creating three folder for PE data

Download reference genome 
 
### Step 2 - Quality Control

For each PE group/pair, we will do QC analysis and preprocessing using fastqc+fastp 

Comparing _fastqc_ and _fastp_ reports 

### Step 3 - Mapping & Visualisation
 
Aligning Illumina PE-reads with reference that runs a sequence alignment program that aligns DNA, creating _.sam _ files  

For each pair, convert _.sam_ to _.bam_ and processing mapping reads 

The loop above indexes and sorts, then visualize the mapped reads in _igv_ viewer for both: original and sorted _.bam_ files 
  
### Step 4 - Primer Clipping   

We checked and related to cleanplex amplicons _.bedpe_ file 

After matching the _.FASTA_ header and _ID_s 

Finally, clipping the primer sequences

### Step 5 - Variant Calling 

PE3 didn't had any variants so a quorum effort was done that gave a lot of unknown, unphased samples using _--report-monomorphic_ (loci which appear to bemonomorphic, and alleles, even those not present in called genotypes)

### Step 6 - Filtering & Masking

A small _R_ Script using _vcfR_ package was used for quality control, visualization, and writing out the new _.vcf_ file. It does annotation for each PE data for analysis and composite plots, with masking low-coverage region based on quality and depth. 

For generating _.vcf_ file statistics for interpretation and reassessment 

For generating some plots using a perl [script](https://github.com/ENCODE-DCC/chip-seq-pipeline/blob/master/dnanexus/bam2tagalign/resources/usr/local/bin/samtools-1.0/bin/plot-bamstats) from [Petr Danecek](https://www.sanger.ac.uk/person/danecek-petr/) 

### Step 7 - Consensus  

After making the variant calls, normalising indels and filtering low-coverage

Put all three consensus from PE datsets into one _.fasta_ file 

### Step 8 - Lineage Annotation 

For the dynamic nomenclature of SARS-CoV-2 lineages a.k.a. pangloin nomenclature, and assigning most likely lineage (Pango lineage) to query sequences

### Step 9 - Consensus QC 

Based on calculating pairwise nucleotide identity and reporting ambiguous _N_'s

### Step 10 - Phylogeny & MSA

For checking the positions of the mismatches, and distance matrix between each and every query and reference sequences

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
