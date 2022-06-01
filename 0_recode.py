#---
# title: Recode depth/count files so missing depth and count are 0/0
#---

######################
# imports
######################
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import sys
from scipy.stats import binom
from scipy.optimize import minimize
from multiprocessing import Pool
    
######################
# Variables
######################
prfx=sys.argv[1]

drecode=open(prfx+'/'+sys.argv[4],'w')
crecode=open(prfx+'/'+sys.argv[5],'w')

with open(prfx+'/'+sys.argv[2]) as f1, open(prfx+'/'+sys.argv[3]) as f2: 
	for d,c in zip(f1,f2):
		d=d.strip('\n').split(' ')
		c=c.strip('\n').split(' ')
		if d[0] == 'CHROM':
			drecode.write(" ".join(d)+'\n')
			crecode.write(" ".join(d)+'\n')
		else:
			indexes = [i for i, x in enumerate(c) if x == "NA"]
			d_modify = d
			c_modify = c
			replacements = '0'
			for index in indexes:		
				d_modify[index] = replacements
				c_modify[index] = replacements
			drecode.write(" ".join(d_modify)+'\n')
			crecode.write(" ".join(c_modify)+'\n')
				
drecode.close()
crecode.close()
