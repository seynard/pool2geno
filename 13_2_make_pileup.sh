#---
# title: Script to produce pileup (prerequisite to perform allele frequency estimation from poolseq with Popoolation2)
#---

#! /bin/bash
######################
# Variables
######################
dir_out=${1}
dir_popoolation=${2}
fasta=${3}
vcf_position=${4}
script=${5}
bam_file=${6}
bam_list=${8}

n=${7}
output_pileup=${n}.pileup
sbatch -W -J pileup_${n} -o ${dir_out}/log/pileup_${n}.out -e ${dir_out}/log/pileup_${n}.err --wrap="module load bioinfo/samtools-1.8; samtools mpileup -I -l ${dir_out}/snp_pos.txt -f ${fasta} -C 50 -q 20 -Q 20 ${bam_file} -o ${dir_out}/${output_pileup}" # run pileup program from Popoolation2
output_sync=${n}.sync
sbatch --mem=40G -W -J sync_${n} -o ${dir_out}/log/sync_${n}.out -e ${dir_out}/log/sync_${n}.err --wrap="java -ea -Xmx10g -jar ${dir_popoolation}/mpileup2sync.jar --fastq-type sanger --min-qual 20 --input ${dir_out}/${output_pileup} --output ${dir_out}/${output_sync}" # run sync, pileup output interpretation program
sed -i "s/a/A/g;s/t/T/g;s/c/C/g;s/g/G/g" ${dir_out}/${n}.sync
echo 'DONE pileup and sync'
if [ -e ${dir_out}/${n}.count ]	
then
echo 'colony ${n} already in' 
else 
python ${script}/13_2_1_combine_col.py ${dir_out}/allele_id.txt ${dir_out}/${n} ${dir_out}/${n}.count
awk '$2=="'"${n}"'"{$(NF+2)=" done"} 1' ${dir_out}/${bam_list} > ${dir_out}/tmp && mv ${dir_out}/tmp ${dir_out}/${bam_list}
fi
rm ${dir_out}/${n}.sync
zip -j ${dir_out}/${n}.pileup.zip ${dir_out}/${n}.pileup
rm ${dir_out}/${n}.pileup


