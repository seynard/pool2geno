#---
# title: Summary statistics (genotyping error, calibration, genotype probability distribution depending on MAF, genotyping error ...) for each homogenenous group
#---

######################
# imports
######################
import pandas as pd 
import numpy as np 
import sys 
import glob
import dask.dataframe as dd
from collections import Counter
import statistics 
import os
import zipfile
import re

######################
# Variables
######################
prfx_in=sys.argv[1]
prfx_out=sys.argv[2]
arr=sys.argv[3]
arr=arr.replace("/","")
depth=sys.argv[4]
i=sys.args[5]

infile=open(prfx_out+'/Mix_'+arr+"/sim_model"+str(depth)+"_"+str(i)+".egs",'r')
firstLine=infile.readline().split(',')
unwanted={'','CHROM','POS','REF','ALT'}
firstLine=[ele for ele in firstLine if ele not in unwanted]

for j in range(0,len(firstLine)):
	print(j)
	ID=firstLine[j-1]
	col_id=ID.split('_')[2]
	simul_id=ID.split('_')[0]
	simul_id=re.sub('simul','',simul_id)
	# load all input for simulations and simulations from real data (on the whole genome and on 1000 markers)
	if '1000' in arr:
		input_geno_queen=dd.read_csv(prfx_in+'/simul_data'+simul_id+"_1000/Geno_queen*",header=0,sep=' ')
		input_geno_queen=input_geno_queen.loc[:,['CHROM','POS','queen'+col_id]]
		input_geno_queen.columns=['CHROM','POS','geno_queen']	
		input_geno_queen['CHROM']=input_geno_queen['CHROM'].replace(['NC_037638.1','NC_037639.1','NC_037640.1','NC_037641.1','NC_037642.1','NC_037643.1','NC_037644.1','NC_037645.1','NC_037646.1','NC_037647.1','NC_037648.1','NC_037649.1','NC_037650.1','NC_037651.1','NC_037652.1','NC_037653.1'],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
		input_geno_queen[['CHROM']]=input_geno_queen[['CHROM']].astype(int)
	elif 'simul_data' in arr:
		input_geno_queen=dd.read_csv(prfx_in+'/'+arr+simul_id+"/Geno_queen*",header=0,sep=' ')
		input_geno_queen=input_geno_queen.loc[:,['CHROM','POS','queen'+col_id]]
		input_geno_queen.columns=['CHROM','POS','geno_queen']	
		input_geno_queen['CHROM']=input_geno_queen['CHROM'].replace(['NC_037638.1','NC_037639.1','NC_037640.1','NC_037641.1','NC_037642.1','NC_037643.1','NC_037644.1','NC_037645.1','NC_037646.1','NC_037647.1','NC_037648.1','NC_037649.1','NC_037650.1','NC_037651.1','NC_037652.1','NC_037653.1'],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
		input_geno_queen[['CHROM']]=input_geno_queen[['CHROM']].astype(int)
	elif 'simul' in arr:
		input_geno_queen=dd.read_csv(prfx_in+'/'+arr+simul_id+"/sim_queen_geno"+str(col_id)+".txt",header=0,sep=' ')
		input_geno_queen.columns=['CHROM','POS','geno_queen']	
		input_geno_queen['CHROM']=input_geno_queen['CHROM'].replace(['NC_037638.1','NC_037639.1','NC_037640.1','NC_037641.1','NC_037642.1','NC_037643.1','NC_037644.1','NC_037645.1','NC_037646.1','NC_037647.1','NC_037648.1','NC_037649.1','NC_037650.1','NC_037651.1','NC_037652.1','NC_037653.1'],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
		input_geno_queen[['CHROM']]=input_geno_queen[['CHROM']].astype(int)
		
	freq=dd.read_csv(prfx_out+'/Mix_'+arr+"/sim_model"+str(depth)+"_"+str(i)+".freqs",header=0,sep=',')
	freq=freq.drop(['Unnamed: 0'],axis=1)
	freq.columns=['CHROM','POS','freq_drones','freq_colony']	
	BG=dd.read_csv(prfx_out+'/Mix_'+arr+"/sim_model"+str(depth)+"_"+str(i)+".bgs",header=0,sep=',')
	BG=BG.loc[:,['CHROM','POS',ID]]
	BG.columns=['CHROM','POS','BG']	
	EG=dd.read_csv(prfx_out+'/Mix_'+arr+"/sim_model"+str(depth)+"_"+str(i)+".egs",header=0,sep=',')
	EG=EG.loc[:,['CHROM','POS',ID]]
	EG.columns=['CHROM','POS','EG']	
	proba=dd.read_csv(prfx_out+'/Mix_'+arr+"/sim_model"+str(depth)+"_"+str(i)+".prob",header=0,sep=',')
	proba=proba.loc[:,['CHROM','POS',ID,'geno']]
	proba_AA=proba.loc[proba['geno'] == 'AA']
	proba_AA=proba_AA.drop(['geno'],axis=1)
	proba_AA.columns=['CHROM','POS','P_AA']
	proba_AR=proba.loc[proba['geno'] == 'AR']
	proba_AR=proba_AR.drop(['geno'],axis=1)
	proba_AR.columns=['CHROM','POS','P_AR']
	proba_RR=proba.loc[proba['geno'] == 'RR']
	proba_RR=proba_RR.drop(['geno'],axis=1)
	proba_RR.columns=['CHROM','POS','P_RR']
	dat=dd.merge(input_geno_queen,BG,on=['CHROM','POS']).merge(EG,on=['CHROM','POS']).merge(proba_AA,on=['CHROM','POS']).merge(proba_AR,on=['CHROM','POS']).merge(proba_RR,on=['CHROM','POS']).merge(freq,on=['CHROM','POS'])
	dat['col']=ID
	dat["geno_queen"]=dat["geno_queen"].astype(str)
	dat["geno_queen"]=dat["geno_queen"].mask(dat["geno_queen"]=="0.0","AA").mask(dat["geno_queen"]=="0.5","AR").mask(dat["geno_queen"]=="1.0","RR")		
	dat['geno_queen']=dat['geno_queen'].replace(['AA','AR','RR'],[0,1,2])
	dat['P_AA'] = pd.to_numeric(dat['P_AA'])
	dat['P_AR'] = pd.to_numeric(dat['P_AR'])
	dat['P_RR'] = pd.to_numeric(dat['P_RR'])
	dat['P_max']=dat[['P_AA','P_AR','P_RR']].max(axis=1)
	dat['P']=dat[['P_AA','P_AR','P_RR']].idxmax(axis=1)
	dat['P']=dat['P'].replace(['P_AA','P_AR','P_RR'],[0,1,2])
	dat['maf']=dat['freq_colony']
	dat['maf']=dat['maf'].mask((dat['maf']>0.5),1-dat['maf'])
	dat['bin_maf']=dat['maf']*20
	dat['bin_maf']=dat['bin_maf'].round()/20
	dat['simul_id']=arr
	dat['bin_proba']=dat['P_max']*20
	dat['bin_proba']=dat['bin_proba'].round()/20
	dat['BEG']=0
	dat['BEG']=dat['BEG'].mask((dat['EG']>0.33) & (dat['EG']<=0.66),1)
	dat['BEG']=dat['BEG'].mask(dat['EG']>0.66,2)
	dat['ok_bg']=dat['geno_queen']==dat['BG']
	dat['ok_eg']=dat['geno_queen']==dat['BEG']
	dat['ok_p']=dat['geno_queen']==dat['P']
	dat=dat.astype({'bin_proba':"category",'bin_maf':"category"})
	dat=dat.compute()

	# Genotyping error
	index=dat['ok_bg'].value_counts().index
	error={'index':index,'ok_bg':dat['ok_bg'].value_counts(),'ok_eg':dat['ok_eg'].value_counts(),'ok_p':dat['ok_p'].value_counts(),'col':ID,'group':i}
	error=pd.DataFrame(data=error)
	print('error')
	if j == 0:
		error.to_csv(prfx_out+'/summary_error_'+arr+str(i)+'.txt',header=True)
	else:
		error.to_csv(prfx_out+'/summary_error_'+arr+str(i)+'.txt',header=False,mode='a')
	# genotyping error as a function of MAF
	index=dat.groupby(['ok_bg','bin_maf'],dropna=False).size().index.to_frame()
	error_maf={'index':index['ok_bg'],'bin_maf':index['bin_maf'],'ok_bg':dat.groupby(['ok_bg','bin_maf']).size(),'ok_eg':dat.groupby(['ok_eg','bin_maf']).size(),'ok_p':dat.groupby(['ok_p','bin_maf']).size(),'col':ID,'group':i}
	error_maf=pd.DataFrame(data=error_maf)
	print('error_maf')
	if j == 0:
		error_maf.to_csv(prfx_out+'/summary_error_maf_'+arr+str(i)+'.txt',header=True)
	else:
		error_maf.to_csv(prfx_out+'/summary_error_maf_'+arr+str(i)+'.txt',header=False,mode='a')
	# genotyping error as a function of genotype probability
	index=dat.groupby(['ok_bg','bin_proba'],dropna=False).size().index.to_frame()
	error_proba={'index':index['ok_bg'],'bin_proba':index['bin_proba'],'ok_bg':dat.groupby(['ok_bg','bin_proba'],dropna=False).size(),'ok_eg':dat.groupby(['ok_eg','bin_proba'],dropna=False).size(),'ok_p':dat.groupby(['ok_p','bin_proba'],dropna=False).size(),'col':ID,'group':i}
	error_proba=pd.DataFrame(data=error_proba)
	print('error_proba')
	if j == 0:
		error_proba.to_csv(prfx_out+'/summary_error_proba_'+arr+str(i)+'.txt',header=True)
	else:
		error_proba.to_csv(prfx_out+'/summary_error_proba_'+arr+str(i)+'.txt',header=False,mode='a')
	# genotyping probability as a function of MAF 
	index=dat.groupby(['bin_maf','bin_proba'],dropna=False).size().index.to_frame()
	maf_proba={'bin_m':index['bin_maf'],'bin_proba':index['bin_proba'],'maf':dat.groupby(['bin_maf','bin_proba'],dropna=False).size(),'col':ID,'group':i}
	maf_proba=pd.DataFrame(data=maf_proba)
	print('maf_proba')
	if j == 0:
		maf_proba.to_csv(prfx_out+'/summary_maf_proba_'+arr+str(i)+'.txt',header=True)
	else:
		maf_proba.to_csv(prfx_out+'/summary_maf_proba_'+arr+str(i)+'.txt',header=False,mode='a')
	# three way analysis of genotyping error as function of MAF bin and genotype probability
	index=dat.groupby(['ok_bg','bin_maf','bin_proba'],dropna=False).size().index.to_frame()
	error_maf_proba1={'index':index['ok_bg'],'bin_m':index['bin_maf'],'bin_proba':index['bin_proba'],'ok_bg':dat.groupby(['ok_bg','bin_maf','bin_proba'],dropna=False).size(),'ok_eg':dat.groupby(['ok_eg','bin_maf','bin_proba'],dropna=False).size(),'ok_p':dat.groupby(['ok_p','bin_maf','bin_proba'],dropna=False).size(),'col':ID,'group':i}
	error_maf_proba=pd.DataFrame(data=error_maf_proba)
	print('error_maf_proba')
	if j == 0:
		error_maf_proba.to_csv(prfx_out+'/summary_error_maf_proba_'+arr+str(i)+'.txt',header=True)
	else:
		error_maf_proba.to_csv(prfx_out+'/summary_error_maf_proba_'+arr+str(i)+'.txt',header=False,mode='a')
	# genotype calibration
	x1=dat[['geno_queen','P_AA']]
	x1.columns=['geno_queen','P']
	x1=x1.assign(geno=0)
	x2=dat[['geno_queen','P_AR']]
	x2.columns=['geno_queen','P']
	x2=x2.assign(geno=1)
	x3=dat[['geno_queen','P_RR']]
	x3.columns=['geno_queen','P']
	x3=x3.assign(geno=2)
	x=x1.append(x2).append(x3)
	cal=x
	cal['bin_proba']=cal['P']*20
	cal['bin_proba']=cal['bin_proba'].round()/20
	cal['ok']=cal['geno_queen']==cal['geno']
	index=cal.groupby(['model','bin_proba'],dropna=False).size().index.to_frame()
	calib={'model':index['model'],'bin_proba':index['bin_proba'],'calib':cal.groupby(['model','bin_proba'])['ok'].mean(),'col':ID,'group':i}
	calib=pd.DataFrame(data=calib)
	print('calib')
	if j == 0:
		calib.to_csv(prfx_out+'/summary_calib_'+arr+str(i)+'.txt',header=True)
	else:
		calib.to_csv(prfx_out+'/summary_calib_'+arr+str(i)+'.txt',header=False,mode='a')
	# genotype probability
	prob_aa=dat[['P_AA']]
	prob_aa.columns=['P']
	prob_ar=dat[['P_AR']]
	prob_ar.columns=['P']
	prob_rr=dat[['P_RR']]
	prob_rr.columns=['P']
	prob=prob_aa.append(prob_ar).append(prob_rr)
	prob['bin']=prob['P']*20
	prob['bin']=prob['bin'].round()/20
	index=prob.groupby(['bin','model'],dropna=False).size().index.to_frame()
	prob={'model':index['model'],'bin_proba':index['bin'],'size':prob.groupby(['bin','model']).size(),'col':ID,'group':i}
	prob=pd.DataFrame(data=prob)
	print('proba')
	if j == 0:
		prob.to_csv(prfx_out+'/summary_proba_'+arr+str(i)+'.txt',header=True)
	else:
		prob.to_csv(prfx_out+'/summary_proba_'+arr+str(i)+'.txt',header=False,mode='a')
		
		
