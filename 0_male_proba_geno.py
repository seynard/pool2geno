#---
# title: Estimate queen genotype from her gametes (male offspring)
#---

######################
# imports
######################
import csv
import sys
import os
from functools import reduce
import glob
import numpy as np
import pandas as pd
import statistics
import itertools
import zipfile

######################
# Variables
######################
epsilon=float(sys.argv[3]) #sequencing error

## Define priors
# prior on queen genotype P(G_Q)
P_prior_Q_0=(1/3)
P_prior_Q_1=(1/3)
P_prior_Q_2=(1/3)

# prior on male genotype knowing queen genotype P(G_D|G_Q)
P_D_0_Q_0=1-epsilon
P_D_1_Q_0=epsilon
P_D_0_Q_1=(1/2)
P_D_1_Q_1=(1/2)
P_D_0_Q_2=epsilon
P_D_1_Q_2=1-epsilon


with open(sys.argv[2]) as male:
	header_line_male=next(male)
	header_line_male=header_line_male.strip('\n').split(' ')

col_name=list(header_line_male[4:len(header_line_male)])
col_name_new=[x[:-2] for x in col_name]
col_name=list(set(col_name_new))
col_name[:] = [item for item in col_name if item != '']

for i in col_name:
	mx=[x for x, s in enumerate(header_line_male) if i in s]
	print(i)
	w=open(sys.argv[1]+'/'+i+'QueenGenoFromMale.txt','w')
	m_g=[]
	with open(sys.argv[2]) as male:
		header_line_male=next(male)
		header_line_male=header_line_male.strip('\n').split(' ')
		for m in male:	
			m=m.strip('\n').split(' ')
			m_g=[m[y] for y in mx]
			chr=m[0]
			pos=m[1]
			m_g=[x for x in m_g if x!='NA']		
			n_male=float(len(m_g))
			n_male_0=m_g.count('0')
			n_male_1=m_g.count('2')				
			# P(G_D observed|G_Q)
			P_D_Q_0=(P_D_0_Q_0)**n_male_0*(P_D_1_Q_0)**n_male_1
			P_D_Q_1=(P_D_0_Q_1)**n_male_0*(P_D_1_Q_1)**n_male_1
			P_D_Q_2=(P_D_0_Q_2)**n_male_0*(P_D_1_Q_2)**n_male_1
			# P(G_Q|G_D)
			P_Q_0_D=P_D_Q_0*P_prior_Q_0
			P_Q_1_D=P_D_Q_1*P_prior_Q_1
			P_Q_2_D=P_D_Q_2*P_prior_Q_2
			# normalisation P(G_Q|G_D)
			P_tot=P_Q_0_D+P_Q_1_D+P_Q_2_D
			P_Q_0=P_Q_0_D/P_tot
			P_Q_1=P_Q_1_D/P_tot
			P_Q_2=P_Q_2_D/P_tot
			s=str(chr)+' '+str(pos)+' '+str(m_g)+' '+str(round(P_Q_0,2))+' '+str(round(P_Q_1,2))+' '+str(round(P_Q_2,2))
			w.write(s+'\n')
