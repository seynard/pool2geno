#---
# title: Combine all pileup output into one file
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
######################
# Variables
######################
all_col=sys.argv[1]
to_add=sys.argv[2]
to_write=sys.argv[3]

i=pd.read_table(all_col,sep=' ',header=None)
i=i.dropna(axis='columns')
i.columns=['CHROM','POS','REF','ALT']
f=pd.read_table(to_add+'.sync',sep='\t',header=None)
f.columns=['CHROM','POS','REF',to_add.split('/')[-1]]
m=pd.merge(i,f,on=['CHROM','POS','REF'],how='outer')
m.fillna('0:0:0:0:0:0',inplace=True)
m.to_csv(to_write,sep=' ',header=True,index=False)
