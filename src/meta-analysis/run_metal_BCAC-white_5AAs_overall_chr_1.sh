#!/bin/bash
module load gcc/6.2.0
module load zlib/1.2.8
module load metal/2011.03.25

rm METAANALYSIS1.TBL
rm METAANALYSIS1.TBL.info
rm meta_analysis_BCAC-white_5AAs_ALL_filtered_chr1.txt
rm meta_analysis_BCAC-white_5AAs_ALL_filtered_chr1_info.txt

metal < metal_BCAC-white_5AAs_overall_chr1_code.txt  #metal_AABC_BCAC_ALL_filtered_chr1_code.txt #metal_AABC_BCAC_chr1.txt

mv METAANALYSIS1.TBL meta_analysis_BCAC-white_5AAs_ALL_filtered_chr1.txt
mv METAANALYSIS1.TBL.info meta_analysis_BCAC-white_5AAs_ALL_filtered_chr1_info.txt
echo "finish chr1"
