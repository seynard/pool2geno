#---
# title: Select the subset of markers used for our simulation
#---

######################
# libraries
######################
installpackages<-function(package_name){if(package_name %in% rownames(installed.packages()) == FALSE) {install.packages(package_name)}}
package_list<-c('data.table','tidyverse')
for(i in 1:length(package_list)){
	installpackages(package_list[i])
	library(package_list[i],character.only=T)}
######################
# Variables
######################
args<-commandArgs(TRUE)
dirin<-args[1]
dirout<-args[2]
type<-args[3]
sce_list<-args[4] 
sce_list<-unlist(strsplit(sce_list,','))
dirin2<-paste0('data/Mix_',type)

# Read all Q vectors estimated by the AM model
Q<-list()
for(i in 1:length(sce_list)){
	if(grepl('1000',type)){lq<-list.files(path=paste0(dirout,'/simul_data',sce_list[i],'_1000'),pattern='_st_het.Q')}else{lq<-list.files(path=paste0(dirout,'/',type,sce_list[i]),pattern='_st_het.Q')}
	q<-list()
	for(j in 1:length(lq)){
		if(grepl('1000',type)){x<-fread(paste0(dirout,'/simul_data',sce_list[i],'_1000/',lq[j]))}else{x<-fread(paste0(dirout,'/',type,sce_list[i],'/',lq[j]))}
		n<-gsub('sim_depth_count','',lq[j])
		n<-gsub('_st_het.Q','',n)
		x$id<-paste0('simul',sce_list[i],'_col_',n)
		x<-spread(x,V1,V2)
		q[[j]]<-x
	}
Q[[i]]<-do.call(rbind,q)
}	
Q<-do.call(rbind,Q)

# Make groups based on Q vectors 
Q$k<-NA
Q$k[is.na(Q$k)]<-'LMC'
Q$k[Q$Mellifera+Q$Caucasica>=0.8 & Q$Mellifera+Q$Ligustica_Carnica<0.85]<-'MC'
Q$k[Q$Ligustica_Carnica+Q$Mellifera>=0.8 & Q$Mellifera+Q$Caucasica<0.85]<-'LM'
Q$k[Q$Ligustica_Carnica+Q$Caucasica>=0.8 & Q$Caucasica+Q$Mellifera<0.85]<-'LC'
Q$k[Q$Mellifera>=0.75]<-'M'
Q$k[Q$Caucasica>=0.75]<-'C'
Q$k[Q$Ligustica_Carnica>=0.75]<-'L'

# List colonies of each scenarios grouped together and make input files (depth and count files) for further analysis
for(i in 1:length(unique(Q$k))){
	q_k<-subset(Q$id,Q$k==unique(Q$k)[i])
	d<-list()
	c<-list()
	for(j in 1:length(q_k)){
		a<-unlist(strsplit(q_k[j],'_'))
		sim_id<-gsub('simul','',a[1])
		col_id<-a[3]
		if(grepl('1000',type)){x<-fread(paste0(dirin,'/simul_data',sim_id,'_1000/sim_depth_count',col_id,'.txt'))}else{x<-fread(paste0(dirin,'/',type,sim_id,'/sim_depth_count',col_id,'.txt'))}
		d[[j]]<-x[,c('CHROM','POS','D')]
		colnames(d[[j]])[3]<-q_k[j]
		c[[j]]<-x[,c('CHROM','POS','X')]
		colnames(c[[j]])[3]<-q_k[j]
		}
	depth<-Reduce(function(x,y) merge(x,y,by=c('CHROM','POS'),all=T),d)
	count<-Reduce(function(x,y) merge(x,y,by=c('CHROM','POS'),all=T),c)
	write.table(paste0(dirin2,'/depth',d,'_',i,'.txt'),col.names=T,row.names=F,quote=F)
	write.table(paste0(dirin2,'/count',c,'_',i,'.txt'),col.names=T,row.names=F,quote=F)
}



