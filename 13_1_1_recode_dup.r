#---
# title: Recode duplicated Bam file when running pool sequencing analysis
#---

######################
# Variables
######################
args<-commandArgs(TRUE)
dir_out<-args[1]
BamList_name<-args[2]

dat<-read.table(paste0(dir_out,'/',BamList_name),sep=' ',header=F,fill=T)
n<-which(duplicated(dat[,2]))
dat[,2]<-as.character(dat[,2])
if(length(n)>0){
for (i in 1:length(n)){
	dat[n[i],2]<-paste0(dat[n[i],2],'_bis')
	}
}
write.table(dat,paste0(dir_out,'/',BamList_name),sep=' ',col.names=F,row.names=F,quote=F)
