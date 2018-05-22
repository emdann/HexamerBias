import pandas as pd
import numpy as np
import os
from hexVSprimed import *
from avg_deltaG_bs import *
from predictCovBs import *

# dir='/hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/'

def make_coverage_tab(dir, filename=''):
    '''
    Table of coverage per template for list of ptcounts files in dir
    '''
    files = [el for el in os.listdir() if 'ptCounts' in el]
    coverageTab = pd.DataFrame()
    for f in files:
        sample = f.split('.')[0]
        ptMat = pd.read_csv(f, compression=findCompr(f), index_col=0)
        cov = pd.DataFrame(ptMat.sum(axis=1), columns = [sample])
        coverageTab = pd.concat([coverageTab, cov], axis=1)
    if filename:
        coverageTab.to_csv(filename)
    return(coverageTab)

def make_deltaG_tab(dir, filename=''):
    '''
    Table of deltaG diagonal for list of ptcounts files in dir
    '''
files = [el for el in os.listdir() if 'ptDg' in el]
totDiagDf = pd.DataFrame(columns=['hex'])
for f in files:
sample = f.split('.')[0]
dgMat = pd.read_csv(f, compression=findCompr(f), index_col=0)
diag={}
for temp in dgMat.iterrows():
    diag[temp[0]]=temp[1][temp[0]]
    diagDf = pd.DataFrame(list(diag.items()), columns=['hex',sample])
    totDiagDf = pd.merge(totDiagDf , diagDf, on='hex', how='outer')
if filename:
    totDiagDf.to_csv(filename)
    return(totDiagDf)
