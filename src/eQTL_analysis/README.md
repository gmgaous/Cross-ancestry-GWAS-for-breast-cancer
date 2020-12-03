The pipeline for the eQTL analysis using data from TCGA breast cancer
D. Huo

1) Extract genotype data from TCGA for the selected loci using gtool

	TCGA_genotype_extraction.pbs

2) Select the genes 2MB around the 6 GWAS loci using Stata 

	 Gene_select_4eQTL.do

3) Conduct eQTL analysis of nearby genes and SNPs in the 6 GWAS loci for Blacks and Whites seperately using R:
	qsub_eQTL2.pbs
	eQTL_analysis2.R


