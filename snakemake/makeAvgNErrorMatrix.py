import fnmatch
import os
import pandas as pd
import numpy as np
import argparse
from hexVSprimed import *

def filter_cells(predictedDgFile, numReads, read_thresh=40000):
    '''
    Finds average of predicted Dg and standard deviation using a selection of cells
    '''
    nReads = pd.read_csv(numReads, sep='\t', dtype=int)
    commonPairs = pd.read_csv(predictedDgFile, sep=',', index_col=0)
    goodCells = [str(cell.cell) for index,cell in nReads.iterrows() if cell.numReads >= read_thresh]
    predictedDgAvg = {}
    for pair,values in commonPairs[goodCells].iterrows():
        mean,sd = np.mean(values),np.std(values)
        predictedDgAvg[pair]={'dg':mean,'sd':sd}
    return(pd.DataFrame(predictedDgAvg).T)

def make_predictedDg_matrix(df, cellAb):
    '''
    Turns file of pt - predicted Dg to matrix of template on row and primer on column
    Needs cell abundance file to fill in missing pt pairs
    '''
    # df = pd.read_csv(file, index_col=0)
    tabDic = {}
    sdDic = {}
    for pt,row in df.iterrows():
        dg,sd =  row.dg, row.sd
        t,p = pt.split('-')
        if t not in tabDic.keys():
            tabDic[t]={}
        tabDic[t][p] = dg
        if t not in sdDic.keys():
            sdDic[t]={}
        sdDic[t][p] = sd
    ptMat = pd.DataFrame(tabDic)
    ptMat = fillNsortPTmatrix(ptMat, cellAb)
    sdMat = pd.DataFrame(sdDic)
    sdMat = fillNsortPTmatrix(sdMat, cellAb)
    return((ptMat,sdMat))
