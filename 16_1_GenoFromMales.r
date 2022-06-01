#---
# title: Combine information from reconstructed queen genotype using offspring information
#---

######################
# libraries
######################
installpackages<-function(package_name){if(package_name %in% rownames(installed.packages()) == FALSE) {install.packages(package_name)}}
package_list<-c('data.table','tidyverse','reshape2','gridExtra')
for(i in 1:length(package_list)){
	installpackages(package_list[i])
	library(package_list[i],character.only=T)}
######################
# Variables
######################
args<-commandArgs(TRUE)
dirin<-args[1]
dirout<-args[2]
prefix<-args[3]
i<-as.numeric(args[4])

list_geno<-list.files(path=dirout,pattern='QueenGenoFromMale.txt')

if(prefix=='mosar'){col<-paste0(gsub('-QueenGenoFromMale.txt','',list_geno[i]),'-F')}else{col<-gsub('-QueenGenoFromMale.txt','',list_geno[i])}
print(col)
males2queen<-fread(paste0(dirout,'/',list_geno[i]),header=F,sep='\t',data.table=F)
a<-strsplit(males2queen[,1],"\\[")
males2queen<-do.call(rbind,a)
a<-colsplit(males2queen[,1]," ",c('CHROM','POS'))
b<-colsplit(males2queen[,2],"\\]",c('v1','v2'))
b$v1<-gsub(", ",",",b$v1)
b$v2<-trimws(b$v2,which="both",whitespace="[ \t\r\n]")
c<-colsplit(b$v2," ",c('P_AA','P_AR','P_RR'))
males2queen<-data.frame('CHROM'=a$CHROM,'POS'=a$POS,'geno_males'=b$v1,'P_AA'=c$P_AA,'P_AR'=c$P_AR,'P_RR'=c$P_RR)
males2queen$CHROM<-as.character(males2queen$CHROM)
males2queen$CHROM[males2queen$CHROM=='NC_037638.1']<-1
males2queen$CHROM[males2queen$CHROM=='NC_037639.1']<-2
males2queen$CHROM[males2queen$CHROM=='NC_037640.1']<-3
males2queen$CHROM[males2queen$CHROM=='NC_037641.1']<-4
males2queen$CHROM[males2queen$CHROM=='NC_037642.1']<-5
males2queen$CHROM[males2queen$CHROM=='NC_037643.1']<-6
males2queen$CHROM[males2queen$CHROM=='NC_037644.1']<-7
males2queen$CHROM[males2queen$CHROM=='NC_037645.1']<-8
males2queen$CHROM[males2queen$CHROM=='NC_037646.1']<-9
males2queen$CHROM[males2queen$CHROM=='NC_037647.1']<-10
males2queen$CHROM[males2queen$CHROM=='NC_037648.1']<-11
males2queen$CHROM[males2queen$CHROM=='NC_037649.1']<-12
males2queen$CHROM[males2queen$CHROM=='NC_037650.1']<-13
males2queen$CHROM[males2queen$CHROM=='NC_037651.1']<-14
males2queen$CHROM[males2queen$CHROM=='NC_037652.1']<-15
males2queen$CHROM[males2queen$CHROM=='NC_037653.1']<-16
males2queen$CHROM<-as.numeric(males2queen$CHROM)
males2queen$maximum_column<-colnames(males2queen)[4:6][apply(males2queen[4:6],1,which.max)]
males2queen$maximum_column[males2queen$maximum_column=='P_AA']<-0
males2queen$maximum_column[males2queen$maximum_column=='P_AR']<-1
males2queen$maximum_column[males2queen$maximum_column=='P_RR']<-2
DF<-males2queen[,c('CHROM','POS','maximum_column')]
colnames(DF)[3]<-col
write.table(DF,paste0(dirout,'/GenoFromMale_',i,'.txt'),col.names=T,row.names=F,sep=' ',quote=F)




