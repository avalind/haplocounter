Haplocounter
---

Haplocounter is a tool for studying XCI using 10x scRNAseq data.
The basic idea is given a vcf of either

- paternal SNPs in the non-PAR region that are heterozygous the sample
- maternal (sample with 10x data) & transmitted het SNPS in the nonPAR region

as well as a genome aligned .bam file from the cellranger pipeline, haplocounter
tallies for each cellular barcode (bam CB-tag) the number of times a variant nucleotide
is expressed for each site in the vcf file, and then uses this data to classify
each CB as to which X-chromosome that was inactivated.

When preparing a .vcf-file   
