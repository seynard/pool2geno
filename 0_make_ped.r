#---
# title: Convert file to .ped (PLINK)
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
dirout<-args[1]
dirin<-args[2]
dat<-args[3]
id<-args[4]

data<-fread(paste0(dirout,'/',dat),data.table=F,header=F)
id<-fread(paste0(dirin,'/',id),data.table=F,header=F)
if(ncol(data)==2){data<-paste(data[,1],data[,2],sep=' ')
}else {data<-apply(data[,1:ncol(data)],1,paste,collapse=" ")}
Dat<-list()
for(i in 1:length(data)){
d<-unlist(strsplit(data[i],' '))
ra<-id[i,]
d[d=='0']<-ra[,4]
d[d=='1']<-ra[,3]
d[d=='-9']<-'0'
Dat[[i]]<-paste0(d,collapse=' ')
}
D<-do.call(rbind,Dat)
write.table(D,paste0(dirout,'/',dat),col.names=F,row.names=F,quote=F)


