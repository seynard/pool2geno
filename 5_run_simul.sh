#---
# title: Generate queen genotype and run model AM
#---

#!/bin/bash
######################
# Variables
######################
dirin=${1}
dirout=${2}
script=${3}
nb_marker=${4}
simul_number=${5}
freq=${6}
n_col=${7}
n_col_prop=${8}
depth=${9}
distrib_q=${10}
distrib_d=${11}

mkdir -p ${dirin}/simul${simul_number}
mkdir -p ${dirout}/simul${simul_number}
cp ${dirin}/${freq} ${dirin}/simul${simul_number}
dirin_simul=''${dirin}'/simul'${simul_number}'/'
dirout_simul=''${dirout}'/simul'${simul_number}'/'
mkdir -p ${dirout_simul}/log
ncpu=${12}
n_col_snpid=${13}
n_pop=${14}
pop_id=${15}
nbmales=${16}
module load bioinfo/plink-v1.90b5.3

Rscript ${script}/5_1_simul.r ${dirin_simul} ${n_col_prop} ${nb_marker} ${depth} ${distrib_q} ${distrib_d} ${freq} ${n_pop} ${pop_id} ${nbmales}> ${dirout_simul}/log/simul.Rout 2>&1
echo 'simul.r done'

############################################################################################
# Run model AM
############################################################################################
for i in $(seq 1 1 ${n_col})
do
	header="CHROM POS D X"
	sed -i "1s/.*/${header}/" ${dirin_simul}/sim_depth_count${i}.txt
	.local/bin/qg_pool --Fmatrix ${dirin}/sim_freq.txt ${dirin_simul}/sim_depth_count${i}.txt -o ${dirout_simul}
done
