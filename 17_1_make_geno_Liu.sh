#---
# title: Bootstrap sampling from publicaly available data from Liu et al. (2015)
#---

#! /bin/bash
######################
# Variables
######################
dirin=${1}
dirout=${2}
prefix=${3}
c=${4}
n=${5}
b=${6}
type=${7}
data=${8}

if [ ${type} = 'cm' ]
then
if [ ${data} != '50k' ]
then 
bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' ${dirin}/${prefix}/vcf_${prefix}_cm_${c}_${n}_${b}.vcf.gz > ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt
else
bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' ${dirin}/${prefix}/vcf_${prefix}_cm_${c}_${n}_${b}_${data}.vcf.gz > ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt
fi
else
if [ ${data} != '50k' ]
then 
bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' ${dirin}/${prefix}/vcf_${prefix}_${c}_${n}_${b}.vcf > ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt
else
bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' ${dirin}/${prefix}/vcf_${prefix}_${c}_${n}_${b}_${data}.vcf > ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt
fi
fi
ncol=$(awk '{print NF;exit}' ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt)
ncol=$((ncol-4))
header="CHROM POS REF ALT " 
for i in $(seq 1 1 ${ncol});do header+=${c}"_"${n}"_"${b}"__ ";done
(echo ${header} && cat ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt) > ${dirin}/tmp_vcf.txt && mv ${dirin}/tmp_vcf.txt ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt
sed -i -e 's,0/0,0,g' -e 's,0|0,0,g' -e 's,1/1,2,g' -e 's,1|1,2,g' -e 's,0/1,1,g' -e 's,0|1,1,g' -e 's,1/0,1,g' -e 's,1|0,1,g' -e 's,./.,-9,g' -e 's,.|.,-9,g' ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt
if [ ${type} = 'cm' ]
then
convers=${dirin}/${prefix}/conversion_${data}.bim
awk 'NR==FNR{a[$1,$2]=$0;next} ($1,$8) in a{print $1" "$4, a[$1,$8]}' ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt ${convers} > ${dirin}/${prefix}/tmp_file.txt
cut -d' ' -f3,4 --complement ${dirin}/${prefix}/tmp_file.txt > ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt
header="CHROM POS REF ALT " 
for i in $(seq 1 1 ${ncol});do header+=${c}"_"${n}"_"${b}"__ ";done
(echo ${header} && cat ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt) > ${dirin}/tmp_vcf.txt && mv ${dirin}/tmp_vcf.txt ${dirin}/${prefix}/geno_${c}_${n}_${b}.txt
fi






