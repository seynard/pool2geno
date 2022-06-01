#---
# title: Combine all summary outputs files
#---

#!/bin/bash
######################
# Variables
######################
script=${1}
dir_out=${2}
dir_save=${3}
dir_in=${4}
dir_popoolation=${5}
fasta=${6}
vcf_file=${7}
vcf_males=${8}
pool_size_def=${9}
nbjobs=${10}

mkdir -p ${dir_in}/temp
mkdir -p ${dir_in}/temp/pileup_mosar
mkdir -p ${dir_in}/temp/count_mosar
mkdir -p ${dir_in}/log

############################################################################################
# From BAM file to allele count, depth and frequencies
############################################################################################
# Colony by colony create depth and count files, rerun when new sequences arrive
if [ -e ${dir_in}/temp/depth_count/depth_mosar.txt.zip ]
then
	sbatch -W --wrap="module load bioinfo/bcftools-1.6; bcftools query -f '%CHROM %POS \n' ${vcf_file} > ${dir_in}/snp_pos.txt; bcftools query -f '%CHROM %POS %REF %ALT \n' ${vcf_file} > ${dir_in}/allele_id.txt"
	unzip ${dir_in}/temp/depth_count/depth_mosar.txt.zip -d ${dir_in}
	unzip ${dir_in}/temp/depth_count/count_ref_mosar.txt.zip -d ${dir_in}
	unzip ${dir_in}/temp/depth_count/count_alt_mosar.txt.zip -d ${dir_in}
else
	sbatch -W --wrap="module load bioinfo/bcftools-1.6; bcftools query -f '%CHROM %POS \n' ${vcf_file} > ${dir_in}/snp_pos.txt; bcftools query -f '%CHROM %POS %REF %ALT \n' ${vcf_file} > ${dir_in}/allele_id.txt; cp ${dir_in}/allele_id.txt ${dir_in}/depth_mosar.txt"
	header="CHROM POS REF ALT "
	awk -v x="${header}" 'NR==1{print x} 1' ${dir_in}/depth_mosar.txt > ${dir_in}/tmp && mv ${dir_in}/tmp ${dir_in}/depth_mosar.txt; cp ${dir_in}/depth_mosar.txt ${dir_in}/count_ref_mosar.txt; cp ${dir_in}/depth_mosar.txt ${dir_in}/count_alt_mosar.txt
fi
sbatch -W -J init -o ${dir_in}/log/init.out -e ${dir_in}/log/init.err --wrap="sh ${script}/13_1_file_init.sh ${dir_save}/BAM_MOSAR_HAV31/ ${dir_in} ${script} BamList_mosar"
sbatch -W --wrap="grep -v 'done' ${dir_in}/BamList_mosar > ${dir_in}/Bam_not_done_mosar"
bam_list2=($(awk -F' ' '{print $1}' ${dir_in}/Bam_not_done_mosar))
rm ${script}/pileup_array2_mosar*
x=1
for i in ${bam_list2[@]} 
do
n=$(awk -F ' ' 'NR==num{print $2}' num=${x} ${dir_in}/Bam_not_done_mosar)
x=$[$x + 1]
echo "sh "${script}"/13_2_make_pileup.sh "${dir_in}" "${dir_popoolation}" "${fasta}" "${vcf_file}" "${script}" "${i}" "${n}" BamList_mosar;" >> ${script}/pileup_array2_mosar.sh
done
split -l ${nbjobs} ${script}/pileup_array2_mosar.sh ${script}/pileup_array2_mosar_split_
pileup_list=(${script}/pileup_array2_mosar_split_*)
for i in ${pileup_list[@]}
do
echo -e "#! /bin/bash\n$(cat ${i})" > ${i}
sbatch ${i}
done
n_job_run=$(squeue -u seynard | grep 'pileup' | wc -l)
while [ ${n_job_run} != 0 ]
do
sleep 10m 
echo ${n_job_run}
n_job_run=$(squeue -u seynard | grep 'pileup' | wc -l)
done
mv ${dir_in}/BS*.zip ${dir_in}/temp/pileup_mosar

count_list=(${dir_in}/BS*.count)
last=$(echo ${count_list[@]:(-1)}|rev|cut -d'/' -f 1|rev | cut -d'.' -f 1)
for i in ${count_list[@]} 
do
n=($(echo ${i}|rev| cut -d '/' -f 1|rev | cut -d'.' -f 1))
if grep ${i} ${dir_in}/depth_mosar.txt
then
echo 'already done'
else
echo 'to do'
if [ $n != $last ] 
then
sbatch -J synctocount -o ${dir_in}/log/synctocount_mosar_${n}.out -e ${dir_in}/log/synctocount_mosar_${n}.err --wrap="python ${script}/13_3_synctocount.py ${i} ${dir_in}" # tailor made script to correct for sequencing error and count alleles 
else
sbatch -W -J synctocount -o ${dir_in}/log/synctocount_mosar_${n}.out -e ${dir_in}/log/synctocount_mosar_${n}.err --wrap="python ${script}/13_3_synctocount.py ${i} ${dir_in}"
fi
fi
done
for i in ${count_list[@]} 
do
n=($(echo ${i}|rev| cut -d '/' -f 1|rev | cut -d'.' -f 1))
sbatch -W --wrap="zip -j ${dir_in}/temp/count_mosar/${n}.count.zip ${dir_in}/${n}.count"
done
rm ${dir_in}/BS*.count

# Merge files across all colonies
ls -1 ${dir_in}/tmp_depth_BS* | split -l 100 -d - lists_mosar
for list in lists_mosar*; do paste $(cat $list) > ${dir_in}/merge${list##lists}; done
sbatch -W --wrap="paste ${dir_in}/merge_mosar* > ${dir_in}/depth_mosar.tmp; paste -d ' ' ${dir_in}/depth_mosar.txt ${dir_in}/depth_mosar.tmp > ${dir_in}/depth2_mosar.txt; mv ${dir_in}/depth2_mosar.txt ${dir_in}/depth_mosar.txt; sed -i 's/  / /g' ${dir_in}/depth_mosar.txt; sed -i 's/\t/ /g' ${dir_in}/depth_mosar.txt"
rm ${dir_in}/merge* lists* ${dir_in}/depth_mosar.tmp
ls -1 ${dir_in}/tmp_count_ref_BS* | split -l 100 -d - lists_mosar
for list in lists_mosar*; do paste $(cat $list) > ${dir_in}/merge${list##lists}; done
sbatch -W --wrap="paste ${dir_in}/merge* > ${dir_in}/count_ref_mosar.tmp; paste -d ' ' ${dir_in}/count_ref_mosar.txt ${dir_in}/count_ref_mosar.tmp > ${dir_in}/count_ref2_mosar.txt; mv ${dir_in}/count_ref2_mosar.txt ${dir_in}/count_ref_mosar.txt; sed -i 's/  / /g' ${dir_in}/count_ref_mosar.txt; sed -i 's/\t/ /g' ${dir_in}/count_ref_mosar.txt"
rm ${dir_in}/merge* lists* ${dir_in}/count_ref_mosar.tmp
ls -1 ${dir_in}/tmp_count_alt_BS* | split -l 100 -d - lists_mosar
for list in lists_mosar*; do paste $(cat $list) > ${dir_in}/merge${list##lists}; done
sbatch -W --wrap="paste ${dir_in}/merge* > ${dir_in}/count_alt_mosar.tmp; paste -d ' ' ${dir_in}/count_alt_mosar.txt ${dir_in}/count_alt_mosar.tmp > ${dir_in}/count_alt2_mosar.txt; mv ${dir_in}/count_alt2_mosar.txt ${dir_in}/count_alt_mosar.txt; sed -i 's/  / /g' ${dir_in}/count_alt_mosar.txt; sed -i 's/\t/ /g' ${dir_in}/count_alt_mosar.txt"
rm ${dir_in}/merge* lists* ${dir_in}/count_alt_mosar.tmp
rm ${dir_in}/tmp_depth_BS*.txt ${dir_in}/tmp_count_ref_BS*.txt ${dir_in}/tmp_count_alt_BS*.txt

# Extract genotypes of individual males offspring of the queen
bcftools query -l ${vcf_file} > ${dir_in}/sample_males.txt
grep 'BS' ${dir_in}/sample_males.txt > ${dir_in}/tmp_BS && mv ${dir_in}/tmp_BS ${dir_in}/sample_males.txt
sbatch -W --wrap="module load bioinfo/bcftools-1.6; bcftools view -S ${dir_in}/sample_males.txt ${vcf_file} > ${vcf_males}"
sbatch -W --wrap="module load bioinfo/bcftools-1.6; bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' ${vcf_males} > ${dir_in}/tmp.vcf && mv ${dir_in}/tmp.vcf ${dir_in}/geno_males_mosar.txt"
grep -Fwf ${dir_in}/snp_pos.txt ${dir_in}/geno_males_mosar.txt > ${dir_in}/tmp.txt && mv ${dir_in}/tmp.txt ${dir_in}/geno_males_mosar.txt 
rm ${vcf_males}
sample=`cat ${dir_in}/sample_males.txt`
sample=$(tr '\n' ' ' <<<${sample})
header='CHROM POS REF ALT '${sample}
awk -v v="${header}" 'NR==1{print v}1' ${dir_in}/geno_males_mosar.txt > ${dir_in}/tmp.txt && mv ${dir_in}/tmp.txt ${dir_in}/geno_males_mosar.txt
sbatch -W --wrap="sed -i -e 's:0/0:2:g' -e 's:0|0:2:g' -e 's:0/1:-9:g' -e 's:0|1:-9:g' -e 's:1/1:0:g' -e 's:1|1:0:g' -e 's:2/2:0:g' -e 's:2|2:0:g' -e 's:./.:-9:g' -e 's:\.:-9:g' ${dir_in}/geno_males_mosar.txt"
sbatch -W --wrap="sed -i -e 's:-91:\.1:g' ${dir_in}/geno_males_mosar.txt"
