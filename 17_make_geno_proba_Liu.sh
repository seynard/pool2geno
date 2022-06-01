#---
# title: Queen genotype reconstruction from offspring data in publicaly available data from Liu et al. (2015) allowing for boostrapping 
#---

#! /bin/bash
######################
# Variables
######################
script=${1}
dirin=${2}
dirout=${3}
prefix=${4}
col=${5}
nind_test=${6}
nboot=${7}
seq_error=${8}
data=${9}
type=${10}

IFS=', ' read -r -a col <<< "$col"
IFS=', ' read -r -a nind_test <<< "$nind_test"
for c in ${col[@]}
do
for n in ${nind_test[@]}
do
for b in $(seq 1 1 ${nboot})
do
sbatch -W --wrap="${script}/17_1_make_geno_Liu.sh ${dirin} ${dirout} ${prefix} ${c} ${n} ${b} ${type} ${data}"
sbatch -o ${dirout}/log/male_2_queen_Liu.out -e ${dirout}/log/male_2_queen_Liu.err --wrap="python ${script}/0_male_proba_geno.py ${dirout}/${prefix} ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt ${seq_error}"
done
done
done
