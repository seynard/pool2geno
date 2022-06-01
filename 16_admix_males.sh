#---
# title: Run ADMIXTURE on queen offspring and queen genotype reconstructed from these offspring 
#---

#! /bin/bash
######################
# Variables
######################
script=${1}
dirin=${2}
dirout=${3}
freq=${4}
n_col_snpid=${5}
ncpu=${6}
n_pop=${7}
pop_id=${8}
prefix=${9}
snp_list=${10}
impute=${11}

dirin_data=''${dirin}'/'${prefix}'/'
dirout_data=''${dirout}'/'${prefix}'/'

############################################################################################
# Reconstruct queen genotype using offspring information
############################################################################################
n=$(ls -l ${dirout_data}/*QueenGenoFromMale.txt | wc -l)
for i in $(seq 1 1 ${n})
do
sbatch -W --mem=10G --wrap="Rscript ${script}/16_1_GenoFromMales.r ${dirin_data} ${dirout_data} ${prefix} ${i}"
done
cp ${dirout_data}/GenoFromMale_1.txt ${dirout_data}/tmp_geno.txt
for i in $(seq 2 1 ${n})
do
cut -d' ' -f3 ${dirout_data}/GenoFromMale_${i}.txt > ${dirout_data}/tmp_geno2.txt
paste -d' ' ${dirout_data}/tmp_geno.txt ${dirout_data}/tmp_geno2.txt > ${dirout_data}/tmp.txt && mv ${dirout_data}/tmp.txt ${dirout_data}/tmp_geno.txt 
done
mv ${dirout_data}/tmp_geno.txt ${dirout_data}/GenoFromMale.txt 
rm ${dirout_data}/GenoFromMale_*.txt

############################################################################################
# ADMIXTURE for the queen genotypes reconstructed from offspring (50 000 reference SNPs)
############################################################################################
Ncol=$(awk '{print NF;exit}' ${dirout_data}/GenoFromMale.txt)
n_col=$((${Ncol}-2))
cp ${dirin}/Unif_k_pop.fam ${dirin_data}/${prefix}_queen.tfam
for i in $(seq 1 1 ${n_col})
do
line="- col${i} 0 0 0 -9"
echo ${line} >> ${dirin_data}/${prefix}_queen.tfam
done
awk '{print $1}' ${dirin_data}/${prefix}_queen.tfam > ${dirin_data}/${prefix}_queen.pop
cp ${dirin_data}/${prefix}_queen.tfam ${dirout_data}/${prefix}_queen.tfam 
cp ${dirin_data}/${prefix}_queen.pop ${dirout_data}/${prefix}_queen.pop 
grep -Fwf ${dirin_data}/snp_list.txt ${dirout_data}/GenoFromMale.txt > ${dirout_data}/GenoFromMale_sub.txt 
cp ${dirin_data}/geno_ref_sub.txt ${dirout_data}/geno_queen.txt
tail -n +2 ${dirout_data}/GenoFromMale_sub.txt > ${dirout_data}/tmp_geno.txt
sort -t' ' -V -k1,1 -k2,2 ${dirout_data}/tmp_geno.txt > ${dirout_data}/tmp_geno2.txt 
cut -d' ' -f1,2 --complement ${dirout_data}/tmp_geno2.txt > ${dirout_data}/tmp_geno.txt
paste -d' ' ${dirout_data}/geno_queen.txt ${dirout_data}/tmp_geno.txt > ${dirout_data}/tmp.txt && mv ${dirout_data}/tmp.txt ${dirout_data}/geno_queen.txt 
sed -i -e 's|0|0 0|g' -e 's|1|0 1|g' -e 's|2|1 1|g' -e 's|-9|-9 -9|g' -e 's|,| |g' ${dirout_data}/geno_queen.txt 
Rscript ${script}/0_make_ped.r ${dirout_data} ${dirin_data} geno_queen.txt allele_id_sub.txt
paste -d' ' ${dirin_data}/${prefix}.map ${dirout_data}/geno_queen.txt > ${dirout_data}/tmp.txt && mv ${dirout_data}/tmp.txt ${dirout_data}/geno_queen.txt
mv ${dirout_data}/geno_queen.txt ${dirout_data}/${prefix}_queen.tped
cp ${dirin_data}/${prefix}_queen.tfam ${dirout_data}/${prefix}_queen.tfam
cp ${dirin_data}/${prefix}_queen.pop ${dirout_data}/${prefix}_queen.pop
plink --tfile ${dirout_data}/${prefix}_queen --make-bed --out ${dirout_data}/${prefix}_queen
cd ${dirout_data}
admixture ${dirout_data}/${prefix}_queen.bed -B ${n_pop} --supervised
cd

############################################################################################
# ADMIXTURE for male offspring of the queen (50 000 reference SNPs)
############################################################################################
cp ${dirin}/geno_males_${prefix}.txt ${dirin_data}/geno_males_${prefix}.txt)
Ncol=$(awk '{print NF;exit}' ${dirin_data}/geno_males_${prefix}.txt)
n_col=$((${Ncol}-${n_col_snpid}))
cp ${dirin}/Unif_k_pop.fam ${dirin_data}/${prefix}_males.tfam
for i in $(seq 1 1 ${n_col})
do
	line="- col${i} 0 0 0 -9"
	echo ${line} >> ${dirin_data}/${prefix}_males.tfam
done
awk '{print $1}' ${dirin_data}/${prefix}_males.tfam > ${dirin_data}/${prefix}_males.pop
cp ${dirin_data}/${prefix}_males.tfam ${dirout_data}/${prefix}_males.tfam 
cp ${dirin_data}/${prefix}_males.pop ${dirout_data}/${prefix}_males.pop 
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_data}/geno_males_${prefix}.txt
grep -Fwf ${dirin_data}/snp_list.txt ${dirin_data}/geno_males_${prefix}.txt > ${dirin_data}/geno_males_${prefix}_sub.txt
tail -n +2 ${dirin_data}/geno_males_${prefix}_sub.txt > ${dirout_data}/tmp_geno.txt
sort -V -k1,1 -k2,2 ${dirout_data}/tmp_geno.txt > ${dirout_data}/tmp_geno2.txt 
cut -d' ' -f1,2,3,4 --complement ${dirout_data}/tmp_geno2.txt > ${dirout_data}/tmp_geno.txt
paste -d' ' ${dirin_data}/geno_ref_sub.txt ${dirout_data}/tmp_geno.txt > ${dirout_data}/tmp.txt && mv ${dirout_data}/tmp.txt ${dirout_data}/${prefix}_males.txt 
sed -i -e 's:-9:-9 -9:g' -e 's:0:0 0:g' -e 's:1:0 1:g' -e 's:2:1 1:g' ${dirout_data}/${prefix}_males.txt
Rscript ${script}/0_make_ped.r ${dirout_data} ${dirin_data} ${prefix}_males.txt allele_id_sub.txt
paste -d' ' ${dirin_data}/${prefix}.map ${dirout_data}/${prefix}_males.txt > ${dirout_data}/tmp.txt && mv ${dirout_data}/tmp.txt ${dirout_data}/${prefix}_males.txt
mv ${dirout_data}/${prefix}_males.txt ${dirout_data}/${prefix}_males.tped
plink --tfile ${dirout_data}/${prefix}_males --recode --make-bed --out ${dirout_data}/${prefix}_males
cd ${dirout_data}
admixture ${dirout_data}/${prefix}_males.bed -B ${n_pop} --supervised 
cd
rm ${dirout_data}/tmp*.txt

