#---
# title: Sample groups of haploid individuals matching our simulation parameters
#---

#!/bin/bash
######################
# Variables
######################
script=${1}
dirout=${2}
dirin=${3}
vcf_file=${4} 
nbmales=${5}
npop=${6}

sample=$(awk '{print $1}' ${dirin}/panel_seqapipop/list_seqapipop.txt)
sample=`echo ${sample} | sed 's/ /,/g'`
bcftools view -s ${sample} ${vcf_file} > ${dirin}/vcf_sample.vcf

sample=$(bcftools query -l ${dirin}/vcf_sample.vcf)
header='CHROM POS REF ALT '${sample}
bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' ${dirin}/vcf_sample.vcf > ${dirin}/geno_males.txt
(echo ${header} && cat ${dirin}/geno_males.txt) > ${dirin}/geno_males.tmp && mv ${dirin}/geno_males.tmp ${dirin}/geno_males.txt
sed -i -e 's:0/0:2:g' -e 's:0|0:2:g' -e 's:0/1:1:g' -e 's:0|1:1:g' -e 's:1/1:0:g' -e 's:1|1:0:g' -e 's:./.:-9:g' -e 's:\.:-9:g' ${dirin}/geno_males.txt
sed -i -e 's:-91:\.1:g' ${dirin}/geno_males.txt
Rscript ${script}/6_1_choice_simul_data.r ${dirin} seqapipop_sample ${nbmales} ${npop}
