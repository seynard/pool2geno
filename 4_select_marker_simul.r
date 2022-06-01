#---
# title: Select the subset of markers used for our simulation
#---

######################
# libraries
######################
installpackages<-function(package_name){if(package_name %in% rownames(installed.packages()) == FALSE) {install.packages(package_name)}}
package_list<-c('data.table','dplyr','tidyr','gtools')
for(i in 1:length(package_list)){
	installpackages(package_list[i])
	library(package_list[i],character.only=T)}
######################
# Variables
######################
args<-commandArgs(TRUE)
dir<-args[1]
nb_marker<-as.numeric(args[2])
n_pop=as.numeric(args[4])
pop_id=args[5]
pop_ID<-unlist(strsplit(pop_id,','))

freq<-as_tibble(fread(paste0(dir,'/',args[3]),data.table=F))
freq<-freq%>% sample_n(nb_marker,replace=F)
list_marker<-paste0(freq$CHROM,'_',freq$POS)
freq_marker<-freq[order(freq$CHROM,freq$POS),]
n<-2+n_pop
freq_marker<-freq_marker[,1:n]
col_order<-c('CHROM','POS',pop_ID)
freq_marker<-freq_marker[,col_order]
write.table(freq_marker,paste0(dir,'/sim_freq.txt'),col.names=T,row.names=F,quote=F)
