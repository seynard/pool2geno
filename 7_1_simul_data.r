#---
# title: Select the subset of markers used for our simulation
#---

######################
# libraries
######################
installpackages<-function(package_name){if(package_name %in% rownames(installed.packages()) == FALSE) {install.packages(package_name)}}
package_list<-c('data.table')
for(i in 1:length(package_list)){
	installpackages(package_list[i])
	library(package_list[i],character.only=T)}
######################
# Variables
######################
args<-commandArgs(TRUE)
dir<-args[1]
simul<-as.numeric(args[2])
depth<-as.numeric(args[3])
nboffspring<-as.numeric(args[4])
nbmales<-as.numeric(args[5])

geno_file<-fread(paste0(dir,'/geno_males_simul_',simul,'.txt'),header=T,data.table=F)
n<-which.min(colnames(geno_file)%in%c('CHROM','POS','REF','ALT'))
# Queen Q genotype
Queen<-geno_file[,c(n,(n+1))]
Queen$geno<-paste0(Queen[,1],'|',Queen[,2])
Queen$geno[is.na(Queen[,1]) | is.na(Queen[,2])]<-NA
Q<-data.frame(Queen$geno)
colnames(Q)<-'geno'
Q$f<-NA
Q$f[Q$geno=='2|2']<-1
Q$f[Q$geno=='0|2']<-0.5
Q$f[Q$geno=='2|0']<-0.5
Q$f[Q$geno=='0|0']<-0
Q$geno<-NULL
colnames(Q)<-paste0('queen',simul)
write.table(Q,paste0(dir,'/queen_simul_data',simul,'.txt'),col.names=T,row.names=F,quote=F)

# Frequencies in sperm
f_drone<-geno_file[,c((n+2):ncol(geno_file))]
f_drone$freq<-(rowSums(f_drone=='2'))/nbmales
Freq_drone<-data.frame(f_drone$freq)
colnames(Freq_drone)<-paste0('pop',simul)
write.table(Freq_drone,paste0(dir,'/freq_drone_simul_data',simul,'.txt'),col.names=T,row.names=F,quote=F)

# Frequencies in pool descending from queen x sperm
if(nboffspring==0){
	freq<-(Q[,1]+Freq_drone)/2
}else{
	ind<-list()
	for(i in 1:nboffspring){
		d<-sample(nbmales,1)
		d<-geno_file[,d]
		q<-sample(c(1,2),1)
		q<-Queen[,q]
		ind[[i]]<-paste0(d,'|',q)
		ind[[i]][ind[[i]]=='0|NA']<-NA
		ind[[i]][ind[[i]]=='1|NA']<-NA
		ind[[i]][ind[[i]]=='NA|0']<-NA
		ind[[i]][ind[[i]]=='NA|1']<-NA
		ind[[i]][ind[[i]]=='NA|NA']<-NA
		ind[[i]][ind[[i]]=='1|0']<-'0|1'
		}
	Geno<-do.call(cbind,ind)
	freq<-(rowSums(Geno=='0|0')+0.5*rowSums(Geno=='0|1'))/nboffspring
}

# Allele count
C<-ceiling(freq*depth)
C<-data.frame(C)
colnames(C)<-paste0('pop',simul)
write.table(C,paste0(dir,'/count_ref_simul_data',simul,'.txt'),col.names=T,row.names=F,quote=F)

# Sequencing depth
D<-data.frame(rep(depth,nrow(Q)))
colnames(D)<-paste0('pop',simul)
write.table(D,paste0(dir,'/depth_simul_data',simul,'.txt'),col.names=T,row.names=F,quote=F)

# Frequencies in pool
Freq<-data.frame(freq)
colnames(Freq)<-paste0('pop',simul)
write.table(Freq,paste0(dir,'/freq_pop_simul_data',simul,'.txt'),col.names=T,row.names=F,quote=F)

