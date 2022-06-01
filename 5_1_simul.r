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
n_col<-args[2]
nb_marker<-as.numeric(args[3])
depth<-as.numeric(args[4])
distrib_q<-args[5]
distrib_d<-args[6]
n_pop=as.numeric(args[8])
pop_id=args[9]
pop_ID<-unlist(strsplit(pop_id,','))
n_males<-as.numeric(args[10])

# Describe parameters of the simulation
scenario<-data.frame('nb_col_sub_sp'=unlist(strsplit(n_col,'~')),'distrib_male_sub_sp'=unlist(strsplit(distrib_d,'~')),'distrib_queen_sub_sp'=unlist(strsplit(distrib_q,'~')))
scenario$nb_col_sub_sp<-as.numeric(as.character(scenario$nb_col_sub_sp))
scenario$distrib_male_sub_sp<-as.character(scenario$distrib_male_sub_sp)
scenario$distrib_queen_sub_sp<-as.character(scenario$distrib_queen_sub_sp)
freq_marker<-fread(paste0(dir,'sim_freq.txt'),header=T,data.table=F)

NCOL<-1
x_depth<-list()
x_count<-list()
for(n in 1:nrow(scenario)){
	X_depth<-list()
	X_count<-list()
	n_col<-scenario$nb_col_sub_sp[n]
	if(n_col>0){
	for(i in 1:n_col){
		# Queen Q matrix and genotype
		Q<-rdirichlet(1,as.numeric(unlist(strsplit(scenario$distrib_queen_sub_sp[n],'_'))))
		print(paste0('Q matrix queen (',pop_id,') :',Q[1],' ',Q[2],' ',Q[3]))
		alpha_queen_lig<-as.numeric(Q[1])
		alpha_queen_mel<-as.numeric(Q[2])
		alpha_queen_cau<-as.numeric(Q[3])
		alpha_queen<-c(alpha_queen_lig,alpha_queen_mel,alpha_queen_cau)
		X1<-apply(rmultinom(nb_marker,1,alpha_queen),2,function(x) {which(x==1)})
		X2<-apply(rmultinom(nb_marker,1,alpha_queen),2,function(x) {which(x==1)})
		ZZ<-data.frame(freq_marker$CHROM,freq_marker$POS,X1,X2)
		F_queen1<-unlist(sapply(1:nb_marker,function(x){freq_marker[x,(X1+2)[x]]}))
		F_queen2<-unlist(sapply(1:nb_marker,function(x){freq_marker[x,(X2+2)[x]]}))
		F_queen<-data.frame('CHROM'=freq_marker$CHROM,'POS'=freq_marker$POS,'freq_queen'=rowMeans(cbind(F_queen1,F_queen2)))
		g1<-(rbinom(nb_marker,1,F_queen1))
		g2<-(rbinom(nb_marker,1,F_queen2))
		G<-(g1+g2)/2
		G<-data.frame(freq_marker$CHROM,freq_marker$POS,G)
		colnames(G)<-c('CHROM','POS','geno_queen')

		# Drones Q matrix and frequencies
		if(distrib_d=='0_0_0'){D<-Q
		}else{D<-rdirichlet(1,as.numeric(unlist(strsplit(scenario$distrib_male_sub_sp[n],'_'))))}
		print(paste0('Q matrix drones (',pop_id,') :',D[1],' ',D[2],' ',D[3]))
		alpha_drone_lig<-as.numeric(D[1])
		alpha_drone_mel<-as.numeric(D[2])
		alpha_drone_cau<-as.numeric(D[3])
		alpha_drone<-c(alpha_drone_lig,alpha_drone_mel,alpha_drone_cau)
		F_drone<-(as.matrix(freq_marker[,3:ncol(freq_marker)])%*%alpha_drone)/sum(alpha_drone)
		G_drone<-rbinom(nb_marker,n_males,F_drone)/n_males
		F_drone<-data.frame(freq_marker$CHROM,freq_marker$POS,G_drone)
		colnames(F_drone)<-c('CHROM','POS','freq_drone')
	
		# Allele count and depth
		x<-rbinom(nb_marker,depth,(F_drone$freq_drone+G$geno_queen)/2)
		x<-data.frame(freq_marker$CHROM,freq_marker$POS,depth,x)
		colnames(x)<-c('CHROM','POS','depth','count')
		X_depth[[NCOL]]<-x$depth
		X_count[[NCOL]]<-x$count

		write.table(F_drone,paste0(dir,'sim_freq_drones',NCOL,'.txt'),col.names=T,row.names=F,quote=F)
		write.table(F_queen,paste0(dir,'sim_freq_queen',NCOL,'.txt'),col.names=T,row.names=F,quote=F)
		write.table(G,paste0(dir,'sim_queen_geno',NCOL,'.txt'),col.names=T,row.names=F,quote=F)
		write.table(x,paste0(dir,'sim_depth_count',NCOL,'.txt'),col.names=T,row.names=F,quote=F)
		write.table(ZZ,paste0(dir,'sim_zz',NCOL,'.txt'),col.names=T,row.names=F,quote=F)
		NCOL<-NCOL+1
	}}else{X_depth[[i]]<-NA; X_count[[i]]<-NA}
x_depth[[n]]<-do.call(cbind,X_depth)
x_count[[n]]<-do.call(cbind,X_count)
x_depth[[n]]<-data.frame(as.character(freq_marker$CHROM),freq_marker$POS,x_depth[[n]])
colnames(x_depth[[n]])<-c('CHROM','POS',rep('depth',n_col))
x_count[[n]]<-data.frame(as.character(x$CHROM),x$POS,x_count[[n]])
colnames(x_count[[n]])<-c('CHROM','POS',rep('count',n_col))
}
X_d<-Reduce(function(x,y) merge(x,y,by=c('CHROM','POS')),x_depth)
X_c<-Reduce(function(x,y) merge(x,y,by=c('CHROM','POS')),x_count)
X_d<-X_d[colSums(!is.na(X_d))>0]
X_c<-X_c[colSums(!is.na(X_c))>0]
print(scenario)
write.table(X_d,paste0(dir,'sim_depth.txt'),col.names=T,row.names=F,quote=F)
write.table(X_c,paste0(dir,'sim_count.txt'),col.names=T,row.names=F,quote=F)







