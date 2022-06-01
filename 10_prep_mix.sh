#---
# title: Make homogeneous populations out of all our simulations
#---

#!/bin/bash
######################
# Variables
######################
dirin=${1}
dirout=${2}
script=${3}
type=${4}
DEPTH=${5}
sce_list=${6}
n_col_snpid=${7}

mkdir -p ${dirin}/Mix_${type}
mkdir -p ${dirout}/Mix_${type}
if [ ${i} = 30 ]
then
Rscript ${script}/10_1_choice_simul_mix.r ${dirin} ${dirout} ${type} ${sce_list} ${dirin}/Mix_${type} ${DEPTH}
else
dirinx=${dirin}/depth${DEPTH}
diroutx=${dirout}/depth${DEPTH}
Rscript ${script}/10_1_choice_simul_mix.r ${dirinx} ${diroutx} ${type} ${sce_list} ${dirin}/Mix_${type} ${DEPTH}
fi
