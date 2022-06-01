#---
# title: Run models on real data 
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

mkdir -p ${dirin}/${prefix}
mkdir -p ${dirout}/${prefix}
dirin_data=''${dirin}'/'${prefix}'/'
dirout_data=''${dirout}'/'${prefix}'/'
mkdir -p ${dirout_data}/log
cp ${dirin}/*${prefix}* ${dirin_data}
cp ${dirin}/${freq} ${dirin_data}
cp ${snp_list} ${dirin_data}/snp_list.txt
header="CHROM,POS,${pop_id}"
awk -f ${script}/0_awk_cut.awk -v c=${header} ${dirin_data}/${freq} > ${dirin_data}/tmp && mv ${dirin_data}/tmp ${dirin_data}/${freq}

awk '{print$1" "$2}' ${dirin_data}/${freq} > ${dirin_data}/snp.txt
n_row=$(wc -l < ${dirin_data}/depth_${prefix}.txt)
nb_marker=$((${n_row}-1))
chr=($(cut -f1 -d ' ' ${dirin_data}/${freq} | sort | uniq ))
chr=("${chr[@]:1}")

sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_data}/depth_${prefix}.txt
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_data}/count_alt_${prefix}.txt
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_data}/count_ref_${prefix}.txt
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_data}/snp.txt
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_data}/${freq}
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_data}/snp_list.txt
chr=($(seq 1 1 16))

############################################################################################
# Run model AM (base don 50 000 reference SNPs)
############################################################################################
Ncol=$(awk '{print NF;exit}' ${dirin_data}/depth_${prefix}.txt)
n_col=$((${Ncol}-${n_col_snpid}))
a=${n_col_snpid}
n=1
for i in $(seq $((${a}+${n})) 1 ${Ncol})
do
	header="CHROM POS D X"
	awk '{print $1" "$2" "$"'"${i}"'"}' ${dirin_data}/depth_${prefix}.txt > ${dirin_data}/tmp_depth
	awk '{print $'${i}'}' ${dirin_data}/count_ref_${prefix}.txt > ${dirin_data}/tmp_count	
	paste -d' ' ${dirin_data}/tmp_depth ${dirin_data}/tmp_count > ${dirin_data}/tmp.txt && mv ${dirin_data}/tmp.txt ${dirin_data}/sim_depth_count${n}.txt 	
	sed -i "1s/.*/${header}/" ${dirin_data}/sim_depth_count${n}.txt
	sbatch -J ${i}_qg_pool${n} --mem=50G --wrap=".local/bin/qg_pool --Fmatrix ${dirin_data}/${freq} ${dirin_data}/sim_depth_count${n}.txt --extract ${dirin_data}/snp_list.txt -o ${dirout_data}"
	let "n+=1"
done
rm ${dirin_data}/tmp_count ${dirin_data}/tmp_depth

############################################################################################
# Run model HP
############################################################################################
sbatch -J ${prefix}_hom --mem=200G -c ${ncpu} -W --wrap="python ${script}/beethoven/model_genoqueen_hom.py ${dirin_data} depth_${prefix}.txt count_ref_${prefix}.txt ${n_col_snpid} ${ncpu} sim_model1 100000 >> ${dirout_data}/log/mod1.out"
mv ${dirin_data}/sim_model1* ${dirout_data}/

