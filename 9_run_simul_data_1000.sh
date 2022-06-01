#---
# title: Run simulation on real data for subset of 1000 markers (haploid individuals to create parents of the pool)
#---

#!/bin/bash
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
prefix0=${10}
impute=${11}

mkdir -p ${dirin}/${prefix}
mkdir -p ${dirout}/${prefix}
cp ${dirin}/${freq} ${dirin}/${prefix}
dirin_simul=''${dirin}'/'${prefix}'/'
dirout_simul=''${dirout}'/'${prefix}'/'
header="CHROM,POS,${pop_id}"
awk -f ${script}/0_awk_cut.awk -v c=${header} ${dirin_simul}/${freq} > ${dirin_simul}/tmp && mv ${dirin_simul}/tmp ${dirin_simul}/${freq}
mkdir -p ${dirout_simul}/log

cp ${dirin}/${prefix0}/Depth_simul_data.txt ${dirin_simul}/
cp ${dirin}/${prefix0}/Count_ref_simul_data.txt ${dirin_simul}/
cp ${dirin}/${prefix0}/Geno_queen_simul_data.txt ${dirin_simul}/
cp ${dirin}/${prefix0}/Freq_drone_simul_data.txt ${dirin_simul}/
cp ${dirin}/${prefix0}/Freq_pop_simul_data.txt ${dirin_simul}/

chr=($(cut -f1 -d ' ' ${dirin}/${freq} | sort | uniq ))
chr=("${chr[@]:1}")
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_simul}/${freq}
awk '{print$1" "$2}' ${dirin_simul}/${freq} > ${dirin_simul}/snp.txt

h=$(head -n1 ${dirin_simul}/Depth_simul_data.txt)
if [[ ! ${h} == *"CHROM"* ]]
then
echo -e "CHROM POS REF ALT\n$(cat ${dirin}/allele_id.txt)" > ${dirin_simul}/tmp
paste -d' ' ${dirin_simul}/tmp ${dirin_simul}/Depth_simul_data.txt > ${dirin_simul}/tmp2 && mv ${dirin_simul}/tmp2 ${dirin_simul}/Depth_simul_data.txt
paste -d' ' ${dirin_simul}/tmp ${dirin_simul}/Count_ref_simul_data.txt > ${dirin_simul}/tmp2 && mv ${dirin_simul}/tmp2 ${dirin_simul}/Count_ref_simul_data.txt
paste -d' ' ${dirin_simul}/tmp ${dirin_simul}/Geno_queen_simul_data.txt > ${dirin_simul}/tmp2 && mv ${dirin_simul}/tmp2 ${dirin_simul}/Geno_queen_simul_data.txt
paste -d' ' ${dirin_simul}/tmp ${dirin_simul}/Freq_drone_simul_data.txt > ${dirin_simul}/tmp2 && mv ${dirin_simul}/tmp2 ${dirin_simul}/Freq_drone_simul_data.txt
paste -d' ' ${dirin_simul}/tmp ${dirin_simul}/Freq_pop_simul_data.txt > ${dirin_simul}/tmp2 && mv ${dirin_simul}/tmp2 ${dirin_simul}/Freq_pop_simul_data.txt
rm ${dirin_simul}/tmp
fi

############################################################################################
# Subset input files to only keep 1000 markers
############################################################################################
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_simul}/Depth_simul_data.txt
grep -Fwf ${dirin_simul}/snp.txt ${dirin_simul}/Depth_simul_data.txt > ${dirin_simul}/tmp && mv ${dirin_simul}/tmp ${dirin_simul}/Depth_simul_data.txt
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_simul}/Count_ref_simul_data.txt
grep -Fwf ${dirin_simul}/snp.txt ${dirin_simul}/Count_ref_simul_data.txt > ${dirin_simul}/tmp && mv ${dirin_simul}/tmp ${dirin_simul}/Count_ref_simul_data.txt
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_simul}/Geno_queen_simul_data.txt
grep -Fwf ${dirin_simul}/snp.txt ${dirin_simul}/Geno_queen_simul_data.txt > ${dirin_simul}/tmp && mv ${dirin_simul}/tmp ${dirin_simul}/Geno_queen_simul_data.txt
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_simul}/Freq_drone_simul_data.txt
grep -Fwf ${dirin_simul}/snp.txt ${dirin_simul}/Freq_drone_simul_data.txt > ${dirin_simul}/tmp && mv ${dirin_simul}/tmp ${dirin_simul}/Freq_drone_simul_data.txt
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_simul}/Freq_pop_simul_data.txt
grep -Fwf ${dirin_simul}/snp.txt ${dirin_simul}/Freq_pop_simul_data.txt > ${dirin_simul}/tmp && mv ${dirin_simul}/tmp ${dirin_simul}/Freq_pop_simul_data.txt
python ${script}/0_recode.py ${dirin_simul} Depth_simul_data.txt Count_ref_simul_data.txt Depth_simul_data_recode.txt Count_ref_simul_data_recode.txt
n_row=$(wc -l < ${dirin_simul}/Depth_simul_data.txt)
nb_marker=$((${n_row}-1))

############################################################################################
# Run model AM
############################################################################################
Ncol=$(awk '{print NF;exit}' ${dirin_simul}/Geno_queen_simul_data.txt)
n_col=$((${Ncol}-${n_col_snpid}))
a=${n_col_snpid}
n=1
for i in $(seq $((${a}+${n})) 1 ${Ncol})
do
	header="CHROM POS D X"
	awk '{print $1" "$2" "$"'"${i}"'"}' ${dirin_simul}/Depth_simul_data_recode.txt > ${dirin_simul}/tmp_depth
	awk '{print $'${i}'}' ${dirin_simul}/Count_ref_simul_data_recode.txt > ${dirin_simul}/tmp_count	
	paste -d' ' ${dirin_simul}/tmp_depth ${dirin_simul}/tmp_count > ${dirin_simul}/tmp.txt && mv ${dirin_simul}/tmp.txt ${dirin_simul}/sim_depth_count${n}.txt 	
	sed -i "1s/.*/${header}/" ${dirin_simul}/sim_depth_count${n}.txt
	.local/bin/qg_pool --Fmatrix ${dirin_simul}/${freq} ${dirin_simul}/sim_depth_count${n}.txt -o ${dirout_simul}
	let "n+=1"
done
echo 'model2 done'
rm ${dirin_simul}/tmp_count ${dirin_simul}/tmp_depth

############################################################################################
# Prep admixture (1000 markers)
############################################################################################
cp ${dirin_simul}/Geno_queen_simul_data.txt ${dirin_simul}/geno_queen.txt
cp ${dirin}/allele_id.txt ${dirin_simul}/allele_id.txt
cp ${dirin}/geno_ref.txt ${dirin_simul}/geno_ref.txt
chr=($(cut -f1 -d  ' ' ${dirin}/allele_id.txt | sort | uniq ))
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_simul}/allele_id.txt
sed -i -e "s/${chr[0]}/1/g" -e "s/${chr[1]}/2/g" -e "s/${chr[2]}/3/g" -e "s/${chr[3]}/4/g" -e "s/${chr[4]}/5/g" -e "s/${chr[5]}/6/g" -e "s/${chr[6]}/7/g" -e "s/${chr[7]}/8/g" -e "s/${chr[8]}/9/g" -e "s/${chr[9]}/10/g" -e "s/${chr[10]}/11/g" -e "s/${chr[11]}/12/g" -e "s/${chr[12]}/13/g" -e "s/${chr[13]}/14/g" -e "s/${chr[14]}/15/g" -e "s/${chr[15]}/16/g" ${dirin_simul}/geno_ref.txt
grep -Fwf ${dirin_simul}/snp.txt ${dirin_simul}/allele_id.txt > ${dirin_simul}/tmp.txt && mv ${dirin_simul}/tmp.txt ${dirin_simul}/allele_id.txt
sed 's: :,:g' ${dirin_simul}/snp.txt > ${dirin_simul}/snp2.txt
grep -Fwf ${dirin_simul}/snp2.txt ${dirin_simul}/geno_ref.txt > ${dirin_simul}/tmp.txt && mv ${dirin_simul}/tmp.txt ${dirin_simul}/geno_ref.txt
sort -V -k1,1 -k2,2 ${dirin_simul}/geno_ref.txt > ${dirin_simul}/tmp_geno.txt && mv ${dirin_simul}/tmp_geno.txt ${dirin_simul}/geno_ref.txt
sort -V -k1,1 -k2,2 ${dirin_simul}/allele_id.txt > ${dirin_simul}/tmp_geno.txt && mv ${dirin_simul}/tmp_geno.txt ${dirin_simul}/allele_id.txt
awk -F',' '{print $1" "$1":"$2" "0" "$2}' ${dirin_simul}/geno_ref.txt > ${dirin_simul}/sim_admix.map
cut -d',' -f1,2,3,4 --complement ${dirin_simul}/geno_ref.txt > ${dirin_simul}/tmp.txt && mv ${dirin_simul}/tmp.txt ${dirin_simul}/geno_ref.txt
Ncol=$(awk '{print NF;exit}' ${dirin_simul}/Geno_queen_simul_data.txt)
n_col=$((${Ncol}-${n_col_snpid}))
cp ${dirin}/Unif_k_pop.fam ${dirin_simul}/sim.tfam
for i in $(seq 1 1 ${n_col})
do
	line="- col${i} 0 0 0 -9"
	echo ${line} >> ${dirin_simul}/sim.tfam
done
awk '{print $1}' ${dirin_simul}/sim.tfam > ${dirin_simul}/sim.pop
cp ${dirin_simul}/sim.tfam ${dirout_simul}/sim.tfam 
cp ${dirin_simul}/sim.pop ${dirout_simul}/sim.pop 

############################################################################################
# Run ADMIXTURE on simulated real data (1000 markers)
############################################################################################
cp ${dirin_simul}/geno_ref.txt ${dirout_simul}/geno_admix.txt
tail -n +2 ${dirin_simul}/geno_queen.txt > ${dirout_simul}/tmp_geno.txt
sort -V -k1,1 -k2,2 ${dirout_simul}/tmp_geno.txt > ${dirout_simul}/tmp_geno2.txt 
cut -d' ' -f1,2,3,4 --complement ${dirout_simul}/tmp_geno2.txt > ${dirout_simul}/tmp_geno.txt
paste -d' ' ${dirout_simul}/geno_admix.txt ${dirout_simul}/tmp_geno.txt > ${dirout_simul}/tmp.txt && mv ${dirout_simul}/tmp.txt ${dirout_simul}/geno_admix.txt 
sed -i -e 's:NA:-9:g' -e 's:0:0:g' -e 's:1:2:g' -e 's:0.5:1:g' -e 's:,: :g' ${dirout_simul}/geno_admix.txt
sed -i -e 's:-9:-9 -9:g' -e 's:0:0 0:g' -e 's:1:0 1:g' -e 's:2:1 1:g' ${dirout_simul}/geno_admix.txt
Rscript ${script}/0_make_ped.r ${dirout_simul} ${dirin_simul} geno_admix.txt allele_id.txt
paste -d' ' ${dirin_simul}/sim_admix.map ${dirout_simul}/geno_admix.txt > ${dirout_simul}/tmp.txt && mv ${dirout_simul}/tmp.txt ${dirout_simul}/geno_admix.txt
mv ${dirout_simul}/geno_admix.txt ${dirout_simul}/sim.tped
plink --tfile ${dirout_simul}/sim --recode --make-bed --out ${dirout_simul}/sim
cd ${dirout_simul}
admixture ${dirout_simul}/sim.bed -B ${n_pop} --supervised 
cd
