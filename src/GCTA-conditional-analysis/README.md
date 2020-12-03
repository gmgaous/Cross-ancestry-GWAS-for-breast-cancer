
The pipeline for the conditional analysis to identify additional signals of loci

1) Download and prepare the reference panel genotype data in Linux shell

	data_prepare1.pbs

2) Update the SNP IDs of reference panel data to be in the format chr:position:A1:A2 so to match with study genotype data

Stata do.file for generating new ID from plink's .bim file:
	Generation_updated_SNPID_conditional_analysis.do

Then use plink to update the SNP IDs
	plink –bfile rs17024629.genotypes.chr1 --update-map rs1702_change_rsid.txt --update-name –make-bed -- rs1702_newID  --noweb

3) Split the reference panel data into AFR and EUR subsets

	plink --bfile rs1702_newID --keep KGP_AFR.txt  --make-bed –out rs1702_newID_AFR_subset --noweb
	plink --bfile rs1702_newID --keep KGP_EUR.txt  --make-bed –out rs1702_newID_EUR_subset --noweb

4) Run the GCTA conditional analysis using the command below, conditioning on the lead index SNP

	gcta --bfile rs1702_newID_AFR_subset --cojo-file rs1702_AFR_500kb_betachanged_chk.txt --cojo-cond rs1702_only.txt --out rs1702_gcta_swapped_AFR_betachanged_chk
	gcta --bfile rs1702_newID_EUR_subset --cojo-file rs1702_EUR_500kb_betachanged_chk.txt --cojo-cond rs1702_only.txt --out rs1702_gcta_swapped_EUR_betachanged_chk

5) Run the GCTA conditional analysis using the command below, conditioning on two top SNPs (lead index SNP and the next most highly significant SNP)

	gcta --bfile rs1702_newID_AFR_subset --cojo-file rs1702_AFR_500kb_betachanged.txt --cojo-cond rs1702_and_top_SNP_only.txt --out rs1702plus_top_cond_AFR_results
	gcta --bfile rs1702_newID_EUR_subset --cojo-file rs1702_EUR_500kb_betachanged.txt --cojo-cond rs1702_and_top_SNP_only.txt --out rs1702plus_top_cond_EUR_results

repeat steps 5 until no SNPs have conditional p value <10e-4

then repeat steps 1-

6) Generation of pooled ORs, p values and CIs from AFR and EUR estimates using meta-analysis method

	Meta-analyss_5AAs_BCAC_white_cond_regression.do 

