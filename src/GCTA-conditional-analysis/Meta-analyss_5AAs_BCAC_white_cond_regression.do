
*Modification of the cma.cojo results output of GCTA and generation of datatset
gen beta_numeric=real( bc)
gen SE_numeric=real( bc_se)
gen pvalue_numeric=real( pc )
save "U:\Samples 5AAs and BCAC\Conditional analysis final GCTA results\rs1702_AFR_cond_results.dta"



**generation of parameters from conditional regression output for merge
use "U:\Samples 5AAs and BCAC\Conditional analysis final GCTA results\rs1702_EUR_cond_results.dta"
drop chr bp refa b se p n freq_geno bc bc_se pc
gen sample_size=228951
rename (freq beta_numeric SE_numeric pvalue_numeric) (freq_EUR beta_numeric_EUR SE_numeric_EUR pvalue_numeric_EUR)
save "U:\Samples 5AAs and BCAC\Conditional analysis final GCTA results\rs1702_EUR_cond_results_for merge.dta"



****Overall
gen beta_numeric=real( bc)
gen SE_numeric=real( bc_se)
gen pvalue_numeric=real( pc )
drop chr bp refa b se p n freq_geno bc bc_se pc
gen sample_size=228951
rename (freq beta_numeric SE_numeric pvalue_numeric) (freq_EUR beta_numeric_EUR SE_numeric_EUR pvalue_numeric_EUR)


gen beta_numeric=real( bc)
gen SE_numeric=real( bc_se)
gen pvalue_numeric=real( pc )
drop chr bp refa b se p n freq_geno bc bc_se pc
gen sample_size=19433
rename (freq beta_numeric SE_numeric pvalue_numeric) (freq_AFR beta_numeric_AFR SE_numeric_AFR pvalue_numeric_AFR)


****ERnegative
gen beta_numeric=real( bc)
gen SE_numeric=real( bc_se)
gen pvalue_numeric=real( pc )
drop chr bp refa b se p n freq_geno bc bc_se pc
gen sample_size=127442
rename (freq beta_numeric SE_numeric pvalue_numeric) (freq_EUR beta_numeric_EUR SE_numeric_EUR pvalue_numeric_EUR)


gen beta_numeric=real( bc)
gen SE_numeric=real( bc_se)
gen pvalue_numeric=real( pc )
drop chr bp refa b se p n freq_geno bc bc_se pc
gen sample_size=12857
rename (freq beta_numeric SE_numeric pvalue_numeric) (freq_AFR beta_numeric_AFR SE_numeric_AFR pvalue_numeric_AFR)



use "U:\Samples 5AAs and BCAC\Conditional analysis final GCTA results\rs1702_EUR_cond_results.dta"
drop chr bp refa b se p n freq_geno bc bc_se pc
gen sample_size=19433
rename (freq beta_numeric SE_numeric pvalue_numeric) (freq_AFR beta_numeric_AFR SE_numeric_AFR pvalue_numeric_AFR)


**merge with European data
merge 1:1 snp using "U:\Samples 5AAs and BCAC\Conditional analysis final GCTA results\GCTA_AFR_EUR_results_for_pooling_metaanalyss\rs1702_AFR_cond_results_for merge.dta"
save "U:\Samples 5AAs and BCAC\Conditional analysis final GCTA results\GCTA_AFR_EUR_results_for_pooling_metaanalyss\rs1702_AFR_EUR_cond_results_merged.dta"



*****Calculation of pooled beta coefficients and p values 5AAs and BCAC White results
gen Z_AFR=-invnormal( pvalue_numeric_AFR/2)
gen Z_EUR=-invnormal( pvalue_numeric_EUR/2)
gen weight_AFR=19433^0.5
gen weight_EUR=228951^0.5
gen Z_pooled=(weight_AFR*Z_AFR + weight_EUR* Z_EUR)/((weight_AFR^2 + weight_EUR^2)^0.5)
gen pvalue_pooled=2*(1-normal(abs( Z_pooled)))
*Calculatio. of pvalue_pooled and CI
gen SE_pooled=(1/(1/SE_numeric_EUR^2 + 1/SE_numeric_AFR^2))^0.5
gen weight_for_beta_AFR=1/(SE_numeric_AFR^2)
gen weight_for_beta_EUR=1/(SE_numeric_EUR^2)
gen beta_pooled=(weight_for_beta_AFR*beta_numeric_AFR + weight_for_beta_EUR* beta_numeric_EUR)/( weight_for_beta_AFR+ weight_for_beta_EUR)
gen OR_pooled=exp( beta_pooled)
gen upper_CIOR_pooled=exp( beta_pooled+1.96* SE_pooled)
gen lower_CIOR_pooled=exp( beta_pooled-1.96* SE_pooled)

*Generation of ORs and CIs for AAs and Whites
gen OR_AFR=exp( beta_numeric_AFR)
gen OR_EUR=exp( beta_numeric_EUR)
gen upper_CIOR_AFR=exp( beta_numeric_AFR+ 1.96*SE_numeric_AFR)
gen lower_CIOR_AFR=exp( beta_numeric_AFR- 1.96*SE_numeric_AFR)
gen upper_CIOR_EUR=exp( beta_numeric_EUR+ 1.96*SE_numeric_EUR)
gen lower_CIOR_EUR=exp( beta_numeric_EUR- 1.96*SE_numeric_EUR)


list snp freq_EUR freq_AFR OR_AFR lower_CIOR_AFR upper_CIOR_AFR OR_EUR lower_CIOR_EUR upper_CIOR_EUR if pvalue_pooled <0.0001
list snp pvalue_numeric_EUR pvalue_numeric_AFR OR_pooled lower_CIOR_pooled upper_CIOR_pooled pvalue_pooled if pvalue_pooled <0.0001

   
 


































***Step 5 conditional anaylsis pipeline - Harmization of the IDs of IKG (USING THE BIM file and the summary statistics file to be inputted into GCTA)

**European ancestry
import delimited "U:\Samples 5AAs and BCAC\rs6038_newID_EUR_subset.bim", clear
save "U:\Samples 5AAs and BCAC/rs6038_harmizatn_1kgID_others_Jan\rs6038_newID_EUR_bimSTATA.dta"
duplicates drop v4, force
rename v4 position
save "U:\Samples 5AAs and BCAC\rs6038_harmizatn_1kgID_others_Jan\rs6038_newID_EUR_bimSTATA.dta", replace
use "U:\Samples 5AAs and BCAC\rs6038_overall_GCTAcojo_EUR.dta", clear    /* this is the dataset for the 8 column file for GCTA*/
split snp_id, p(:)
drop snp_id1 snp_id3 snp_id4
destring snp_id2, replace
rename snp_id2 position
duplicates drop position, force
save U:\Samples 5AAs and BCAC\rs6038_harmizatn_1kgID_others_Jan\rs6038_overall_GCTAcojo_EUR_STATA.dta", replace
merge 1:1 position using "U:\Samples 5AAs and BCAC\rs6038_harmizatn_1kgID_others_Jan\rs6038_newID_EUR_bimSTATA.dta"
save "U:\Samples 5AAs and BCAC/rs6038_harmizatn_1kgID_others_Jan\rs6038_EUR_bim_500_merged.dta"
keep if _merge==3
list a1 a0 v5 v6 in 1/10
gen match=1 if (a1 ==v5) & (a0 ==v6)
tab match
replace a1=v5 if match==.
replace a0=v6 if match==.
replace b=-b if match==.    /* In  addition to swapping the effect and ref allele, beta sign also has to be changed*/
keep a1 a0 freq b se p sample_size v2
order v2 
save "U:\Samples 5AAs and BCAC\rs6038_harmizatn_1kgID_others_Jan\rs6038_EUR_500kb_betachanged.dta"





**African ancestry
import delimited "U:\Samples 5AAs and BCAC/rs6038_newID_AFR_subset.bim", clear
save "U:\Samples 5AAs and BCAC/rs6038_harmizatn_1kgID_others_Jan\rs6038_newID_AFR_bimSTATA.dta"
duplicates drop v4, force
rename v4 position
save "U:\Samples 5AAs and BCAC\rs6038_harmizatn_1kgID_others_Jan\rs6038_newID_AFR_bimSTATA.dta", replace
use "U:\Samples 5AAs and BCAC\rs6038_overall_GCTAcojo_AFR.dta", clear    /* this is the dataset for the 8 column file for GCTA*/
split markername, p(:)
drop markername1 markername3 markername4 markername5 markername6
destring markername2, replace
rename markername2 position
duplicates drop position, force
save U:\Samples 5AAs and BCAC/rs6038_harmizatn_1kgID_others_Jan\rs6038_overall_GCTAcojo_AFR.dta", replace
merge 1:1 position using "U:\Samples 5AAs and BCAC\rs6038_harmizatn_1kgID_others_Jan\rs6038_newID_AFR_bimSTATA.dta"
save "U:\Samples 5AAs and BCAC/rs6038_harmizatn_1kgID_others_Jan\rs6038_AFR_bim_500_merged.dta"
keep if _merge==3
gen allele_match=1 if (eff_allele == v5) & (ref_allele== v6)
tab allele_match, m
replace eff_allele=v5 if allele_match ==.
replace ref_allele =v6 if allele_match ==.
replace effect =-(effect) if allele_match==.
keep v2 eff_allele ref_allele freq1 effect stderr pvalue sample_size
order v2
save "U:\Samples 5AAs and BCAC\rs6038_harmizatn_1kgID_others_Jan\rs6038_AFR_500_betachanged.dta"

