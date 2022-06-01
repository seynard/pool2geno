#---
# title: Estimate average sequencing depth per individual 
#---

#! /bin/bash
######################
# Variables
######################
depth_in=${1}
depth_out=${2}

ncol=$(awk '{print NF;exit}' ${depth_in})
unset total
for i in $(seq 5 1 ${ncol})
do
echo ${i}
awk -v i=${i} '{ total += $i } END { print total/NR }' ${depth_in} >> ${depth_out}
done


