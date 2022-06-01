#---
# title: Create input files for simulations on real data (haploid individuals to create parents of the pool)
#---

#!/bin/bash
######################
# Variables
######################
script=${1}
dirin=${2}
dirout=${3}
nbmales=${4}
depth=${5}
nboffspring=${6}
nbsimul=${7}
prefix=${8}

mkdir -p ${dirin}/${prefix}${nbsimul}
mkdir -p ${dirin}/${prefix}${nbsimul}/log
mkdir -p ${dirout}/${prefix}${nbsimul}
mkdir -p ${dirout}/${prefix}${nbsimul}/log
sample=($(head -n1 ${dirin}/geno_males.txt))
awk -v i=${nbsimul} -F' ' '{if($1==i) print}' ${dirin}/${prefix}.txt > ${dirin}/${prefix}${nbsimul}/tmp.txt 
n_col=$(wc -l < ${dirin}/${prefix}${nbsimul}/tmp.txt)
for x in $(seq 1 1 ${n_col})
do
	id=($(sed "${x}q;d" ${dirin}/${prefix}${nbsimul}/tmp.txt))
	id=("${id[@]:1}")
	sample2='$1" "$2" "$3" "$4'
	for a in ${id[@]}
	do
		n=$(echo ${sample[@]} | tr -s " " "\n" | grep -w -n ${a} | cut -d":" -f 1)
		sample2+='" "$'${n}
	done
	awk -f <(echo "{print ${sample2}}") ${dirin}/geno_males.txt > ${dirin}/${prefix}${nbsimul}/geno_males_simul_${x}.txt
	Rscript ${script}/7_1_simul_data.r ${dirin}/${prefix}${nbsimul} ${x} ${depth} ${nboffspring} ${nbmales}
done
ls -v ${dirin}/${prefix}${nbsimul}/queen_simul_data*.txt | xargs paste -d' ' > ${dirin}/${prefix}${nbsimul}/Geno_queen_simul_data.txt
ls -v ${dirin}/${prefix}${nbsimul}/freq_pop_simul_data*.txt | xargs paste -d' ' > ${dirin}/${prefix}${nbsimul}/Freq_pop_simul_data.txt
ls -v ${dirin}/${prefix}${nbsimul}/freq_drone_simul_data*.txt | xargs paste -d' ' > ${dirin}/${prefix}${nbsimul}/Freq_drone_simul_data.txt
ls -v ${dirin}/${prefix}${nbsimul}/count_ref_simul_data*.txt | xargs paste -d' ' > ${dirin}/${prefix}${nbsimul}/Count_ref_simul_data.txt
ls -v ${dirin}/${prefix}${nbsimul}/depth_simul_data*.txt | xargs paste -d' ' > ${dirin}/${prefix}${nbsimul}/Depth_simul_data.txt
rm ${dirin}/${prefix}${nbsimul}/queen_simul_data*.txt ${dirin}/${prefix}${nbsimul}/geno_males_simul*.txt ${dirin}/${prefix}${nbsimul}/freq_pop_simul_data*.txt ${dirin}/${prefix}${nbsimul}/freq_drone_simul_data*.txt ${dirin}/${prefix}${nbsimul}/count_ref_simul_data*.txt ${dirin}/${prefix}${nbsimul}/depth_simul_data*.txt
