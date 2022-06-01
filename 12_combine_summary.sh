#---
# title: Combine all summary outputs files
#---

#!/bin/bash
######################
# Variables
######################
dirin=${1}
dirout=${2}
type=${3}
DEPTH=${4}

head -1 ${dirout}/Mix_${type}/summary_error_maf_proba_${DEPTH}_1.txt > ${dirout}/Mix_${type}/summary_error_maf_proba_${DEPTH}.txt
tail -n +2 -q ${dirout}/Mix_${type}/summary_error_maf_proba_${DEPTH}_*.txt >> ${dirout}/Mix_${type}/summary_error_maf_proba_${DEPTH}.txt
rm  ${dirout}/Mix_${type}/summary_error_maf_proba_${DEPTH}_*.txt
head -1 ${dirout}/Mix_${type}/summary_error_maf_${DEPTH}_1.txt > ${dirout}/Mix_${type}/summary_error_maf_${DEPTH}.txt
tail -n +2 -q ${dirout}/Mix_${type}/summary_error_maf_${DEPTH}_*.txt >> ${dirout}/Mix_${type}/summary_error_maf_${DEPTH}.txt
rm  ${dirout}/Mix_${type}/summary_error_maf_${DEPTH}_*.txt
head -1 ${dirout}/Mix_${type}/summary_error_proba_${DEPTH}_1.txt > ${dirout}/Mix_${type}/summary_error_proba_${DEPTH}.txt
tail -n +2 -q ${dirout}/Mix_${type}/summary_error_proba_${DEPTH}_*.txt >> ${dirout}/Mix_${type}/summary_error_proba_${DEPTH}.txt
rm  ${dirout}/Mix_${type}/summary_error_proba_${DEPTH}_*.txt
head -1 ${dirout}/Mix_${type}/summary_error_${DEPTH}_1.txt > ${dirout}/Mix_${type}/summary_error_${DEPTH}.txt
tail -n +2 -q ${dirout}/Mix_${type}/summary_error_${DEPTH}_*.txt >> ${dirout}/Mix_${type}/summary_error_${DEPTH}.txt
rm  ${dirout}/Mix_${type}/summary_error_${DEPTH}_*.txt
head -1 ${dirout}/Mix_${type}/summary_maf_proba_${DEPTH}_1.txt > ${dirout}/Mix_${type}/summary_maf_proba_${DEPTH}.txt
tail -n +2 -q ${dirout}/Mix_${type}/summary_maf_proba_${DEPTH}_*.txt >> ${dirout}/Mix_${type}/summary_maf_proba_${DEPTH}.txt
rm  ${dirout}/Mix_${type}/summary_maf_proba_${DEPTH}_*.txt
head -1 ${dirout}/Mix_${type}/summary_calib_${DEPTH}_1.txt > ${dirout}/Mix_${type}/summary_calib_${DEPTH}.txt
tail -n +2 -q ${dirout}/Mix_${type}/summary_calib_${DEPTH}_*.txt >> ${dirout}/Mix_${type}/summary_calib_${DEPTH}.txt
rm  ${dirout}/Mix_${type}/summary_calib_${DEPTH}_*.txt
head -1 ${dirout}/Mix_${type}/summary_proba_${DEPTH}_1.txt > ${dirout}/Mix_${type}/summary_proba_${DEPTH}.txt
tail -n +2 -q ${dirout}/Mix_${type}/summary_proba_${DEPTH}_*.txt >> ${dirout}/Mix_${type}/summary_proba_${DEPTH}.txt
rm  ${dirout}/Mix_${type}/summary_proba_${DEPTH}_*.txt
