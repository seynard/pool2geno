#---
# title: Prepare Bam files for real data of pool sequence experiment
#---

#!/bin/bash

######################
# Variables
######################
dir_in=$1
dir_out=$2
script=$3
BamList_name=$4

bam_list=($(find ${dir_in} -name *.bam))

if [ -e ${dir_out}/${BamList_name}.bz2 ]
then
    bzip2 -d ${dir_out}/${BamList_name}.bz2
	for i in ${bam_list[@]} 
	do
	if [[ ${i} != *'ToMerge'* ]]
	then
	n=($(echo ${i} | rev | cut -d '/' -f 1 | rev | cut -d'_' -f 1))
	l=($(bzgrep -R ${i} ${dir_out}/BamList))
	if [ ${#l[@]} -eq 0 ]
	then
	echo "${i} ${n}" >> ${dir_out}/${BamList_name}
	fi
	fi
	done
elif [ -e ${dir_out}/${BamList_name} ]
then
	for i in ${bam_list[@]} 
	do
	if [[ ${i} != *'ToMerge'* ]]
	then
	n=($(echo ${i} | rev | cut -d '/' -f 1 | rev | cut -d'_' -f 1))
	l=($(grep -R ${i} ${dir_sonia}/BamList))
	if [ ${#l[@]} -eq 0 ]
	then
	echo "${i} ${n}" >> ${dir_out}/${BamList_name}
	fi
	fi
	done	
else
	for i in ${bam_list[@]} 
	do
	if [[ ${i} != *'ToMerge'* ]]
	then
	n=($(echo ${i} | rev | cut -d '/' -f 1 | rev | cut -d'_' -f 1))
	printf "%s\n" "${i} ${n}" >> ${dir_out}/${BamList_name}
	fi
	done
fi
Rscript ${script}/13_1_1_recode_dup.r ${dir_out} ${BamList_name}
