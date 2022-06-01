#---
# title: Use the diversity panel from Wragg et al 2021 (https://doi.org/10.1101/2021.09.20.460798) to define reference populations and estimate parameters for these populations (genotype, allele frequencies)
#---

#! /bin/bash
######################
# Variables
######################
script=${1}
dirin=${2}
dir_save=${3}
dirout=${4}
fasta=${5}
vcf_file=${6}
n_pop=${7}
ncpu=${8}
pop_id=${9}
unif_threshold=${10}
mkdir -p ${dirin}/temp
mkdir -p ${dirin}/log
mkdir -p ${dirin}/panel_seqapipop

############################################################################################
# Definition of genetic composition for the individuals in the diversity panel, under the assumption of 3 reference populations
############################################################################################
sbatch -W --wrap="module load bioinfo/plink-v1.90b5.3; plink --vcf ${dirin}/vcf.vcf --keep-allele-order --make-bed --out ${dirin}/seqapipop"
awk '{print $1" "$1}' ${dirin}/list_seqapipop.txt > ${dirin}/sample_seqapipop.txt 
sbatch -W --wrap="module load bioinfo/plink-v1.90b5.3; plink --bfile ${dirin}/seqapipop --keep ${dirin}/sample_seqapipop.txt --keep-allele-order --make-bed --out ${dirin}/seqapipop_sample"
K=${n_pop}
cd ${dirin}	
sbatch -W --mem=20G -J admix_seqapipop_sample -o ${dirin}/log/admix_seqapipop_sample_${K}.out -e ${dirin}/log/admix_seqapipop_sample_${K}.err --wrap="admixture --cv -j${ncpu} seqapipop_sample.bed ${K}" # run ADMXITURE for k=3 populations, this allows for the estimation of reference allele frequencies in each of the 3 reference population defined
cd
sbatch -J plot_admixture -o ${dirin}/log/plot_admixture_seqapipop.out -e ${dirin}/log/plot_admixture_seqapipop.err --wrap="Rscript ${script}/2_1_admix_plot.r ${dirin} seqapipop_sample ${pop} apriori_pop.txt" # plot results from ADMIXTURE knowing pop of origin a priori

############################################################################################
# Definition of reference individuals, under the assumption of 3 reference populations, for a defined threshold
############################################################################################
sbatch -W --mem=20G -J prep_admixture_${n_pop} -o ${dirin}/log/prep_admixture_${n_pop}.out -e ${dirin}/log/prep_admixture_${n_pop}.err --wrap="Rscript ${script}/2_2_prep_admix_k.r ${dirin} seqapipop_sample ${pop_id} ${unif_threshold}" # choose reference individuals for a defined genetic composition threshold
sample=($(awk '{print $1}' ${dirin}/Unif_k.fam))
sample2=''
for i in ${sample[@]}
do
sample2+=','${i}
done
sample2="${sample2:1}"
sbatch -W --wrap="module load bioinfo/bcftools-1.6; bcftools view -s ${sample2} ${vcf_file} > ${dirin}/vcf_subsp.vcf" # extract only reference individuals from vcf
IFS=', ' read -r -a popID <<< "${pop_id}"
for i in ${popID[@]}
do
samplei=$(awk -v var=${i} -F' ' '{if($1==var)print $2}' ${dirin}/Unif_k_pop.fam )
samplei=$(echo $samplei | sed 'y/ /,/')
sbatch -W --wrap="module load bioinfo/bcftools-1.6; bcftools view -s ${samplei} ${vcf_file} > ${dirin}/vcf_${i}.vcf" # make a vcf for each of the reference population (in our case 3)
cp ${dirin}/vcf_${i}.vcf ${dirin}/vcf_${i}_recode.vcf
bcftools annotate --set-id +'%CHROM:%POS' ${dirin}/vcf_${i}_recode.vcf > ${dirin}/tmp.vcf; mv ${dirin}/tmp.vcf ${dirin}/vcf_${i}_recode.vcf
bgzip -c ${dirin}/vcf_${i}.vcf > ${dirin}/vcf_${i}.vcf.gz
bgzip -c ${dirin}/vcf_${i}_recode.vcf > ${dirin}/vcf_${i}_recode.vcf.gz
mv ${dirin}/vcf_${i}* ${dirin}/panel_seqapipop
done
awk '{print $1}' ${dirin}/Unif_k_pop.fam > ${dirin}/vcf_subsp.pop
cp ${dirin}/seqapipop_sample.${n_pop}.P ${dirin}/freq_admix${n_pop}.txt
tail -n +2 ${dirin}/snp_pos.txt > ${dirin}/snp_tmp.txt
paste --delimiter=' ' ${dirin}/snp_tmp.txt ${dirin}/freq_admix${n_pop}.txt > ${dirin}/temp.txt && mv ${dirin}/temp.txt ${dirin}/freq_admix${n_pop}.txt
rm ${dirin}/snp_tmp.txt
col_id=$(cat ${dirin}/ColNames.txt)
header='CHROM POS '${col_id}
header=$(echo ${header} | sed 's/\n/ /g')
awk -v v="$header" 'NR==1{print v}1' ${dirin}/freq_admix${n_pop}.txt > ${dirin}/temp.txt && mv ${dirin}/temp.txt ${dirin}/freq_admix${n_pop}.txt
sbatch -J freq_plot -o ${dirin}/log/freq_plot.out -e ${dirin}/log/freq_plot.err --wrap="Rscript ${script}/2_3_freq_plot.r ${dirin} freq_admix${n_pop}.txt" # estimate reference allele frequency for each population independently
sbatch -W --wrap="bgzip -c ${dirin}/vcf_subsp.vcf > ${dirin}/vcf_subsp.vcf.gz"
mv ${dirin}/{*seqapipop*,vcf_subsp*,vcf_sample*,plot_freq*}  ${dirin}/panel_seqapipop
