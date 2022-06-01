#---
# title: Subset the file containing allele frequencies for the reference populations based on MAF and variance to select a subset of markers for which allele frequencies in the reference population are interesting for our study
#---

######################
# libraries
######################
installpackages<-function(package_name){if(package_name %in% rownames(installed.packages()) == FALSE) {install.packages(package_name)}}
package_list<-c('data.table','dplyr','tidyr','matrixStats')
for(i in 1:length(package_list)){
	installpackages(package_list[i])
	library(package_list[i],character.only=T)}
######################
# Variables
######################
args<-commandArgs(TRUE)
dirin<-args[1]
pop<-args[2]
pop=unlist(strsplit(pop,','))
n_pop=length(pop)

# load the file containing allele frequencies for the reference population for all the markers
freq<-fread(paste0(dirin,'/freq_admix',n_pop,'.txt'),data.table=F)
# estimate MAF for each marker
for(i in 1:length(pop)){freq[,paste0('MAF_',pop[i])]<-NA
freq[freq[,pop[i]]<0.5,paste0('MAF_',pop[i])]<-freq[freq[,pop[i]]<0.5,pop[i]]
freq[freq[,pop[i]]>=0.5,paste0('MAF_',pop[i])]<-(1-freq[freq[,pop[i]]>=0.5,pop[i]])}
# estimate variance for each marker
freq$var<-rowVars(as.matrix(freq[,pop]))
x<-grep('MAF_',colnames(freq))
freq$meanmaf<-rowMeans(as.matrix(freq[,x]))
freq$lowmaf<-apply(freq[,x],1,FUN=min)

# selection 1: on variance
freq.v<-subset(freq,freq$var>0.1)
fv<-freq.v[,c('CHROM','POS',pop)]
write.table(fv,row.names=F,quote=F,file=paste0(dirin,'/freq_admix_v.txt'))
# selection 2: on average MAF
freq.m<-subset(freq,freq$meanmaf>0.1)
fm<-freq.m[,c('CHROM','POS',pop)]
write.table(fm,row.names=F,quote=F,file=paste0(dirin,'/freq_admix_m.txt'))
# selection 3: on lowest MAF
freq.l<-subset(freq,freq$lowmaf>0.1)
fl<-freq.l[,c('CHROM','POS',pop)]
write.table(fl,row.names=F,quote=F,file=paste0(dirin,'/freq_admix_l.txt'))
# selection 3: on variance and average MAF USED IN OUR ANALYSIS
freq.mv<-subset(freq,freq$var>0.1 & freq$meanmaf>0.1)
fmv<-freq.mv[,c('CHROM','POS',pop)]
write.table(fmv,row.names=F,quote=F,file=paste0(dirin,'/freq_admix_mv.txt'))
# selection 3: on variance and lowest MAF
freq.lv<-subset(freq,freq$var>0.1 & freq$lowmaf>0.1)
flv<-freq.lv[,c('CHROM','POS',pop)]
write.table(flv,row.names=F,quote=F,file=paste0(dirin,'/freq_admix_lv.txt'))
