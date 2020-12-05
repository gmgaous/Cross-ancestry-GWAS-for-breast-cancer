Summary statistics were from meta-analysis for overall, ER-positive, and ER-negative. Only top 10,000 SNPs included in the csv files.
The meta-analysis combined GWAS results for one BCAC study on women of European ancestry (BCAC-white) and five studies on women of African Ancestry (BCAC-black, AABC, Ghana, Amber, ROOT) 
Based on meta-analysis results, we selected the top SNPs with MinFreq >0.005 and MaxFreq <0.995 (see below for definition of MinFreq and MaxFreq).
The summary statistics csv files include the following columns:
# MarkerName    - this is the marker name
# Allele1   - the first allele for this marker in the first file (listed in the metal command file) where it occurs
# Allele2   - the second allele for this marker in the first file where it occurs
# Freq1       - weighted average of frequency for allele 1 across all studies
# FreqSE      - corresponding standard error for allele frequency estimate
# MinFreq     - minimum frequency for allele 1 across all studies
# MaxFreq     - maximum frequency for allele 1 across all studies
# Effect    - overall estimated effect size for allele1
# StdErr    - overall standard error for effect size estimate
# P-value   - meta-analysis p-value
# Direction - summary of effect direction for each study, with one '+' or '-' per study
# HetISq    - I^2 statistic which measures heterogeneity on scale of 0-100%
# HetChiSq  - chi-squared statistic in simple test of heterogeneity
# df        - degrees of freedom for heterogeneity statistic
# HetPVal   - P-value for heterogeneity statistic
