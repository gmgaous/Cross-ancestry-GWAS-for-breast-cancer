*** eqtl analysis
*** select genes 2MB around the 6 GWAS loci
*** D Huo, B Adedokun 


** Extraction of gene positions form Ensembl Human file (Homo_sapiens.GRCh37.87.chr.gtf)
** import delimited "Y:\ROOT\Tunde\eQTL\Homo_sapiens.GRCh37.87.chr.gtf", delimiter(tab) clear
import delimited "Y:\ROOT\Tunde\eQTL\Homo_sapiens.GRCh37.87.chr.gtf", rowrange(6) clear 
keep if v3=="gene"

gen loci ="1p13" if (v4>108179756 & v4<112179756) & v1=="1" 
replace loci ="1q41" if (v4>213330292 & v4<217330292) & v1=="1" 
replace loci ="5q31" if (v4>129801947 & v4<133801947) & v1=="5" 
replace loci ="7q11" if (v4>72359358 & v4<76359358) & v1=="7" 
replace loci ="15q24" if (v4>73147332 & v4<77147332) & v1=="15" 
replace loci ="15q26" if (v4>98907094 & v4<102907094) & v1=="15"  
tab loci 
keep if loci!=""

split v9, p(;)

split v91, p() 
split v92, p() 
split v93, p() 
split v94, p() 
split v95, p() 

gen chr="chr"+v1 

*Generation of ENS ID (Combination of ENS number and gene version number of the for ENSG00000123.7 where the number after the dot is the gene version number)
gen version = subinstr(v922, char(34), "", .)
gen gen_ensembl= subinstr(v912, char(34), "", .)
gen gene_ID= gen_ensembl + "." + version

gen gene_name = subinstr(v932, char(34), "", .)
gen gene_biotype  = subinstr(v952, char(34), "", .)

ren v4 gene_start 
ren v5 gene_end 
keep chr loci gene_ID gene_name gene_biotype gene_start gene_end

duplicates tag gene_name , gen(dup_genename)
list if dup_genename >0
drop if dup_genename>0

drop dup_genename

export delimited using "Y:\ROOT\Tunde\eQTL\Gene_position_6Loci.csv", delimiter(tab) replace

