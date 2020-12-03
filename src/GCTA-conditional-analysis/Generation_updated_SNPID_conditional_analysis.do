***Conditional analysis Six novel loci


*****rs1702

import delimited "U:\Samples 5AAs and BCAC\rs17024629.genotypes.chr1.bim", clear

tostring v4, gen (postion_string)

gen markerID="1" + ":" + postion_string + ":" + v6 + ":" + v5

save "U:\Samples 5AAs and BCAC\rs1702_bim_converted.dta"

keep v2 markerID

duplicates drop v2, force

save "U:\Samples 5AAs and BCAC\rs1702_update_SNP.dta", replace

export delimited using "U:\Samples 5AAs and BCAC\rs1702_change_rsid.txt", delimiter(tab) replace





*****rs6793

import delimited "U:\Samples 5AAs and BCAC\rs67931591.genotypes.chr1.bim", clear

tostring v4, gen (postion_string)

gen markerID="1" + ":" + postion_string + ":" + v6 + ":" + v5

save "U:\Samples 5AAs and BCAC\rs6793_bim_converted.dta"

keep v2 markerID

duplicates drop v2, force

save "U:\Samples 5AAs and BCAC\rs6793_update_SNP.dta", replace

export delimited using "U:\Samples 5AAs and BCAC\rs6793_change_rsid.txt", delimiter(tab) replace



*****rs2522

import delimited "U:\Samples 5AAs and BCAC\rs2522057.genotypes.chr5.bim", clear

tostring v4, gen (postion_string)

gen markerID="5" + ":" + postion_string + ":" + v6 + ":" + v5

save "U:\Samples 5AAs and BCAC\rs2522_bim_converted.dta"

keep v2 markerID

duplicates drop v2, force

save "U:\Samples 5AAs and BCAC\rs2522_update_SNP.dta", replace

export delimited using "U:\Samples 5AAs and BCAC\rs2522_change_rsid.txt", delimiter(tab) replace




*****rs1637

import delimited "U:\Samples 5AAs and BCAC\rs1637365.genotypes.chr7.bim", clear

tostring v4, gen (postion_string)

gen markerID="7" + ":" + postion_string + ":" + v6 + ":" + v5

save "U:\Samples 5AAs and BCAC\rs1637_bim_converted.dta"

keep v2 markerID

duplicates drop v2, force

save "U:\Samples 5AAs and BCAC\rs1637_update_SNP.dta", replace

export delimited using "U:\Samples 5AAs and BCAC\rs1637_change_rsid.txt", delimiter(tab) replace



******rs1869

import delimited "U:\Samples 5AAs and BCAC\rs1869959.genotypes.chr15.bim", clear

tostring v4, gen (postion_string)

gen markerID="15" + ":" + postion_string + ":" + v6 + ":" + v5

save "U:\Samples 5AAs and BCAC\rs1869_bim_converted.dta"

keep v2 markerID

duplicates drop v2, force

save "U:\Samples 5AAs and BCAC\rs1869_update_SNP.dta", replace

export delimited using "U:\Samples 5AAs and BCAC\rs1869_change_rsid.txt", delimiter(tab) replace




*****rs1813

import delimited "U:\Samples 5AAs and BCAC\rs181337095.genotypes.chr15.bim", clear

tostring v4, gen (postion_string)

gen markerID="15" + ":" + postion_string + ":" + v6 + ":" + v5

save "U:\Samples 5AAs and BCAC\rs1813_bim_converted.dta"

keep v2 markerID

duplicates drop v2, force

save "U:\Samples 5AAs and BCAC\rs1813_update_SNP.dta", replace

export delimited using "U:\Samples 5AAs and BCAC\rs1813_change_rsid.txt", delimiter(tab) replace




****Generation of datafile to impute for the gcta-cojo 8 column
import delimited "U:\AABC Onco splits used for 23 Xsomes\ER NEGATIVE\meta_analysis_5AAs_ERneg_filtered_chr1.txt"
split markername, p(:)
drop markername1 markername3 markername4 markername5 markername6
destring markername2, replace
keep if markername2<215830292 & markername2>214830292
list markername allele1 allele2 if markername2==215330292
gen eff_allele=upper( allele1)
gen ref_allele=upper(allele2)
gen sample_size=12857
keep markername eff_allele ref_allele freq1 effect stderr pvalue sample_size
order markername eff_allele ref_allele freq1 effect stderr pvalue sample_size
save "U:\Samples 5AAs and BCAC\rs6793_ERneg_GCTAcojo_AFR.dta"
export delimited using "U:\Samples 5AAs and BCAC\rs6793_500kb_AFR_GCTA.txt", delimiter(tab) replace







