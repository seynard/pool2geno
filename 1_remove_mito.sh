#---
# title: Remove markers on mitochondria
#---

#!/bin/bash
######################
# Variables
######################
dir_in=${1}
vcf_file=${2}
mito=${3}

awk '$5 !~ /,/' < ${vcf_file} | grep -v ${mito} > ${dir_in}/vcf.vcf
