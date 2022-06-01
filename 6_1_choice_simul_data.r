#---
# title: Select the subset of markers used for our simulation
#---

######################
# libraries
######################
installpackages<-function(package_name){if(package_name %in% rownames(installed.packages()) == FALSE) {install.packages(package_name)}}
package_list<-c('data.table','ggplot2','factoextra','NbClust')
for(i in 1:length(package_list)){
	installpackages(package_list[i])
	library(package_list[i],character.only=T)}
######################
# Variables
######################
args<-commandArgs(TRUE)
dirin<-args[1]
type<-args[2]
n_males<-as.numeric(args[3])
npop<-as.numeric(args[4])


choice_scenario<-fread(paste0(dirin,'/choice_scenario_simul_data.txt'),head=F,data.table=T)
n_scenario<-nrow(choice_scenario)

# Load reference data set (Wragg et al. 2022)
pca<-fread(paste0(dirin,'/panel_seqapipop/',type,'.',npop,'.Q'),header=F,data.table=F)
fam<-fread(paste0(dirin,'/panel_seqapipop/',type,'.fam'),header=F,data.table=F)
pc<-cbind(fam[,1],pca)
colnames(pc)<-c('ind','Q1','Q2','Q3')
pc$group<-NA
pc$group[grep('AOC',pc$ind)]<-'AOC'
pc$group[grep('BER',pc$ind)]<-'BER'
pc$group[grep('BR',pc$ind)]<-'BR'
pc$group[grep('Buck',pc$ind)]<-'Buck'
pc$group[grep('CAR',pc$ind)]<-'CAR'
pc$group[grep('CAU',pc$ind)]<-'CAU'
pc$group[grep('CHI',pc$ind)]<-'CHI'
pc$group[grep('CIReine',pc$ind)]<-'CIReine'
pc$group[grep('DAN',pc$ind)]<-'DAN'
pc$group[grep('ESP',pc$ind)]<-'ESP'
pc$group[grep('FL',pc$ind)]<-'FL'
pc$group[grep('ITA',pc$ind)]<-'ITA'
pc$group[grep('ITSAP',pc$ind)]<-'ITSAP'
pc$group[grep('JFM',pc$ind)]<-'JFM'
pc$group[grep('JSC',pc$ind)]<-'JSC'
pc$group[grep('KF',pc$ind)]<-'KF'
pc$group[grep('NM',pc$ind)]<-'NM'
pc$group[grep('OUE',pc$ind)]<-'OUE'
pc$group[grep('PacBio',pc$ind)]<-'OUE'
pc$group[grep('POL',pc$ind)]<-'POL'
pc$group[grep('POR',pc$ind)]<-'POR'
pc$group[grep('SLO',pc$ind)]<-'SLO'
pc$group[grep('SOL',pc$ind)]<-'SOL'
pc$group[grep('Sar',pc$ind)]<-'Sar'
pc$group[grep('Sav',pc$ind)]<-'Sav'
pc$group[grep('TES',pc$ind)]<-'TES'
pc$group[grep('Ticino',pc$ind)]<-'Ticino'
pc$group[grep('UK',pc$ind)]<-'UK'
pc$group[grep('XC',pc$ind)]<-'XC'
pc$group[grep('YC',pc$ind)]<-'YC'

# Make groups of individuals based on their Q vectors 
pc$k<-NA
pc$k[is.na(pc$k)]<-'123'
pc$k[pc$Q1+pc$Q2>=0.8 & pc$Q1+pc$Q3<0.85]<-'12'
pc$k[pc$Q1+pc$Q3>=0.8 & pc$Q1+pc$Q2<0.85]<-'13'
pc$k[pc$Q2+pc$Q3>=0.8 & pc$Q2+pc$Q1<0.85]<-'23'
pc$k[pc$Q1>=0.75]<-'1'
pc$k[pc$Q2>=0.75]<-'2'
pc$k[pc$Q3>=0.75]<-'3'
t<-data.frame(table(pc$group,pc$k))
t<-subset(t,t[,3]>0)
colnames(t)<-c('group','k','freq')
pc$k_id<-NA
for(i in 1:length(unique(pc$k))){
	tk<-t[t$k==unique(pc$k)[i],]
	if(tk$group[which.max(tk$freq)]=='OUE'){pc$k_id[pc$k==unique(pc$k)[i]]<-'M'
	}else if(tk$group[which.max(tk$freq)]=='CAU'){pc$k_id[pc$k==unique(pc$k)[i]]<-'C'
	}else if('ITA' %in% tk$group){pc$k_id[pc$k==unique(pc$k)[i]]<-'L'
	}else if(tk$group[which.max(tk$freq)]=='AOC'){pc$k_id[pc$k==unique(pc$k)[i]]<-'LM'
	}else if('Sav'%in%tk$group & 'JSC'%in%tk$group){pc$k_id[pc$k==unique(pc$k)[i]]<-'MC'
	}else if('JFM'%in%tk$group & 'ITSAP'%in%tk$group){pc$k_id[pc$k==unique(pc$k)[i]]<-'LC'
	}else {pc$k_id[pc$k==unique(pc$k)[i]]<-'LMC'}}	
pc$k_id<-factor(pc$k_id,levels=c('L','M','C','LC','LM','MC','LMC'))

# Choose individuals to perform simulation with under the acenario of interest
pc$ind<-as.character(pc$ind)
list_ind<-matrix(ncol=(n_males+2+1))
colnames(list_ind)<-c('scenario','Q1','Q2',paste0('M',seq(1:n_males)))
for(i in 1:n_scenario){
	x<-as.character(unlist(choice_scenario[i,1]))
	q<-as.character(unlist(choice_scenario[i,2]))
	d<-as.character(unlist(choice_scenario[i,3]))
	n<-as.character(unlist(choice_scenario[i,4]))
	if(grepl('~',n)){
		n<-as.character(unlist(strsplit(n,'~')))
		q<-as.character(unlist(strsplit(q,'~')))
		d<-as.character(unlist(strsplit(d,'~')))
		D<-data.frame(q,d,n)
		for(k in 1:nrow(D)){
			n<-as.numeric(as.character(D[k,'n']))
			q<-as.character(D[k,'q'])
			d<-as.character(D[k,'d'])
			for(j in 1:n){
				if(q==d){
					if(q==0 & d==0){L[[j]]<-sample(pc$ind,(n_males+2),replace=F)
					}else if (q=='x' & d=='x'){
						g<-sample(unique(pc$k_id),1,replace=F)
						while(length(pc$ind[pc$k_id==g])<n_males+2){g<-sample(unique(pc$k_id),1,replace=F)}
						lq<-sample(pc$ind[pc$k_id==g],2,replace=F)
						ld<-sample(pc$ind[pc$k_id==g & !(pc$ind%in%lq)],n_males,replace=F)				
						list_ind<-rbind(list_ind,c(x,lq,ld))
					}else if(nchar(q)==1 & nchar(d)==1){
						p<-unlist(strsplit(q,''))
						list_ind<-rbind(list_ind,c(x,sample(pc$ind[pc$k_id%in%p],(n_males+2),replace=F)))
					}else{
						p<-unlist(strsplit(q,''))
						A<-list()
						for(a in 1:length(p)){A[[a]]<-grep(p[a],unique(pc$k_id))}
						g<-unique(pc$k_id)[Reduce(intersect,A)]
						list_ind<-rbind(list_ind,c(x,sample(pc$ind[pc$k_id%in%g],(n_males+2),replace=F)))
					}
				}else{
					if(q==0){
						lq<-sample(pc$ind,2,replace=F)
					}else if(q!=0 & nchar(q)==1){
						p<-unlist(strsplit(q,''))
						lq<-sample(pc$ind[pc$k_id%in%p],2,replace=F)
					}else{
						p<-unlist(strsplit(q,''))
						A<-list()
						for(a in 1:length(p)){A[[a]]<-grep(p[a],unique(pc$k_id))}
						g<-unique(pc$k_id)[Reduce(intersect,A)]
						lq<-sample(pc$ind[pc$k_id%in%g],2,replace=F)
					}
					if(d==0){
						ld<-sample(pc$ind[!(pc$ind%in%lq)],n_males,replace=F)
					}else if(d!=0 & nchar(d)==1){
						p<-unlist(strsplit(d,''))
						ld<-sample(pc$ind[pc$k_id%in%p & !(pc$ind%in%lq)],n_males,replace=F)
					}else{
						p<-unlist(strsplit(d,''))
						A<-list()
						for(a in 1:length(p)){A[[a]]<-grep(p[a],unique(pc$k_id))}
						g<-unique(pc$k_id)[Reduce(intersect,A)]
						ld<-sample(pc$ind[pc$k_id%in%g & !(pc$ind%in%lq)],n_males,replace=F)
					}	
					list_ind<-rbind(list_ind,c(x,lq,ld))
				}
			}
		}
	}else{	
	for(j in 1:n){
		if(q==d){
			if(q==0 & d==0){list_ind<-rbind(list_ind,c(x,sample(pc$ind,(n_males+2),replace=F)))
			}else if (q=='x' & d=='x'){
				g<-sample(unique(pc$k_id),1,replace=F)
				while(length(pc$ind[pc$k_id==g])<n_males+2){g<-sample(unique(pc$k_id),1,replace=F)}
				lq<-sample(pc$ind[pc$k_id==g],2,replace=F)
				ld<-sample(pc$ind[pc$k_id==g & !(pc$ind%in%lq)],n_males,replace=F)				
				list_ind<-rbind(list_ind,c(x,lq,ld))
			}else if(nchar(q)==1 & nchar(d)==1){
				p<-unlist(strsplit(q,''))
				list_ind<-rbind(list_ind,c(x,sample(pc$ind[pc$k_id%in%p],(n_males+2),replace=F)))
			}else{
				p<-unlist(strsplit(q,''))
				A<-list()
				for(a in 1:length(p)){A[[a]]<-grep(p[a],unique(pc$k_id))}
				g<-unique(pc$k_id)[Reduce(intersect,A)]
				list_ind<-rbind(list_ind,c(x,sample(pc$ind[pc$k_id%in%g],(n_males+2),replace=F)))
			}
		}else{
			if(q==0){
				lq<-sample(pc$ind,2,replace=F)
			}else if(q!=0 & nchar(q)==1){
				p<-unlist(strsplit(q,''))
				lq<-sample(pc$ind[pc$k_id%in%p],2,replace=F)
			}else{
				p<-unlist(strsplit(q,''))
				A<-list()
				for(a in 1:length(p)){A[[a]]<-grep(p[a],unique(pc$k_id))}
				g<-unique(pc$k_id)[Reduce(intersect,A)]
				lq<-sample(pc$ind[pc$k_id%in%g],2,replace=F)
			}
			if(d==0){
				ld<-sample(pc$ind[!(pc$ind%in%lq)],n_males,replace=F)
			}else if(d!=0 & nchar(d)==1){
				p<-unlist(strsplit(d,''))
				ld<-sample(pc$ind[pc$k_id%in%p & !(pc$ind%in%lq)],n_males,replace=F)
			}else{
				p<-unlist(strsplit(d,''))
				A<-list()
				for(a in 1:length(p)){A[[a]]<-grep(p[a],unique(pc$k_id))}
				g<-unique(pc$k_id)[Reduce(intersect,A)]
				ld<-sample(pc$ind[pc$k_id%in%g & !(pc$ind%in%lq)],n_males,replace=F)
			}	
		list_ind<-rbind(list_ind,c(x,lq,ld))
		}
	}		
	}
}		
list_ind<-list_ind[-1,]
list_ind<-as.data.frame(list_ind)
list_ind$ID[list_ind$scenario==0]<-'LMC x LMC'
list_ind$ID[list_ind$scenario==1]<-'LMC x LMC'
list_ind$ID[list_ind$scenario==2]<-'L x L'
list_ind$ID[list_ind$scenario==3]<-'M x M'
list_ind$ID[list_ind$scenario==4]<-'C x C'
list_ind$ID[list_ind$scenario==5]<-'L/M/C x L/M/C'
list_ind$ID[list_ind$scenario==6]<-'L/M x L/M'
list_ind$ID[list_ind$scenario==7]<-'L/C x L/C'
list_ind$ID[list_ind$scenario==8]<-'M/C x M/C'
list_ind$ID[list_ind$scenario==9]<-'L/M/C x LMC'
list_ind$ID[list_ind$scenario==10]<-'L/M x M/L'
list_ind$ID[list_ind$scenario==11]<-'L/C x C/L'
list_ind$ID[list_ind$scenario==12]<-'M/C x C/M'
list_ind$ID[list_ind$scenario==13]<-'LM x LMC'
list_ind$ID[list_ind$scenario==14]<-'LC x LMC'
list_ind$ID[list_ind$scenario==15]<-'MC x LMC'
list_ind$ID[list_ind$scenario==16]<-'LM x LM'
list_ind$ID[list_ind$scenario==17]<-'LC x LC'
list_ind$ID[list_ind$scenario==18]<-'MC x MC'
list_ind$ID[list_ind$scenario==19]<-'LM x L/M'
list_ind$ID[list_ind$scenario==20]<-'LC x L/C'
list_ind$ID[list_ind$scenario==21]<-'MC x M/C'
list_ind$ID[list_ind$scenario==22]<-'L x M'
list_ind$ID[list_ind$scenario==23]<-'L x C'
list_ind$ID[list_ind$scenario==24]<-'M x L'
list_ind$ID[list_ind$scenario==25]<-'M x C'
list_ind$ID[list_ind$scenario==26]<-'C x L'
list_ind$ID[list_ind$scenario==27]<-'C x M'

# Representation of the individuals, associated to their respective groups, based on their Q vectors
g<-list()
for(i in 1:length(unique(list_ind$scenario))){
li<-subset(list_ind,list_ind$scenario==unique(list_ind$scenario)[i])
q_scenario<-unique(unlist(li[,2:3]))
d_scenario<-unique(unlist(li[,4:18]))
pciq<-subset(pc,pc$ind%in%q_scenario)
pciq$type<-'queen'
pcid<-subset(pc,pc$ind%in%d_scenario)
pcid$type<-'drones'
pci<-rbind(pciq,pcid)
g[[i]]<-ggplot(pci)+
	geom_point(aes(x=Q1,y=Q2,col=type),alpha=0.3)+
	theme_bw()+ 
	theme(legend.title=element_blank())+
	xlim(0,1)+
	ylim(0,1)+
	xlab('')+
	ylab('')+ggtitle(li$ID)	
	}

pdf(paste0(dirin,'/choice_simul_data.pdf'),width=8,height=5)
ggplot(pc)+geom_point(aes(x=Q1,y=Q2,col=as.factor(k_id)))+theme_bw()+ theme(legend.title=element_blank())+scale_colour_manual(values=c('goldenrod1','grey34','chartreuse3','darkolivegreen1','brown','darkgreen','blue'))+xlab('')+ylab('')
ggplot(pc)+geom_point(aes(x=Q1,y=Q2,col=as.factor(k_id)))+facet_wrap(~group)+theme_bw()
dev.off()

# Write list of individuals used for each simulation as ind1 and ind2 creating the queen and the remaining individuals being the sperm pool
write.table(list_ind,paste0(dirin,'/simul_data.txt'),col.names=F,row.names=F,quote=F)
