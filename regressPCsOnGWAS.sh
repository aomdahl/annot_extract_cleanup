#!/bin/bash
#Quick script for running many regressions of GWAS vs annotations at once.
#Inputs are- 
# 1)list of GWAS paths to run against,
# 2) path to annotations to regress against, requires 
# 3) Extension in output file to save as (usually cc for coding common or ncc for noncoding common)
# 4) Number of top correlating PCs to keep 
ml gcc/5.5.0
ml R
set -e
ANNOT=$2
EXT=$3
QUERY=$1
TOPN=$4
for i in `cat ${QUERY}`; do
    t=`(basename $i| cut -f 1 -d "_")`
    echo $t
    Rscript /work-zfs/abattle4/ashton/prs_dev/prs_tools/annotation_analysis/top_annots_script.R  --sum_stats $i --annots ${ANNOT} --topn ${TOPN} --output ${t}.${EXT}
done
~                                                                                                                                                        
~               
