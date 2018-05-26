import fnmatch
import os
import pandas as pd
import numpy as np
import argparse
from hexVSprimed import *

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Find average predicted Delta G from common pt pairs\n By Emma Dann")
argparser.add_argument('commonPtPairs', type=str, help='File of predicted Dg for common PtPairs')
argparser.add_argument('numReads', type=str, help='File of predicted Dg for common PtPairs')
args = argparser.parse_args()

def filter_cells(predictedDgFile, numReads, read_thresh=40000):
    '''
    Finds average of predicted Dg and standard deviation using a selection of cells
    '''
    nReads = pd.read_csv(numReads, sep='\t', dtype=int)
    commonPairs = pd.read_csv(predictedDgFile, sep=',', index_col=0)
    commonPairs = np.log(commonPairs)
    goodCells = [str(cell.cell) for index,cell in nReads.iterrows() if cell.numReads >= read_thresh]
    predictedDgAvg = {}
    for pair,values in commonPairs[goodCells].iterrows():
        mean,sd = np.mean(values),np.std(values)
        predictedDgAvg[pair]={'dg':mean,'sd':sd}
    return(pd.DataFrame(predictedDgAvg).T)

def select_cell_matrix(folder, numReads, read_thresh=40000):
    '''
    Selects dG tables of cells witj=h right number of reads
    '''
    nReads = pd.read_csv(numReads, sep='\t', dtype=int)
    files = [f for f in os.listdir(folder) if 'ptDg' in f]
    mats = {}
    for file in files:
        sample = file.split('_ptDg_')[0]
        m = pd.read_csv(dir+file, index_col=0)
        flatM = np.array(m).flatten()
        logFlatM = np.log(flatM)
        mats[sample]=logFlatM

    bigArray = np.array(list(mats.values()))
    avg = np.nanmean(bigArray, axis=0).reshape((4096,4096))
    pd.DataFrame(avg, columns=tab.columns, index=tab.index).to_csv('predictedDg.csv')



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
    ptMat = fillNsortPTmatrix(ptMat, cellAb, fillna=-99999)
    sdMat = pd.DataFrame(sdDic)
    sdMat = fillNsortPTmatrix(sdMat, cellAb, fillna=0)
    return((ptMat,sdMat))

predictedDgFile = args.commonPtPairs
numReads = args.numReads
cellAbFile = numReads.split('.numReads.txt')[0] + '.cellAbundance.noN.csv'
cellAb = pd.read_csv(cellAbFile, index_col=0, compression = findCompr(cellAbFile))

# Make predicted Dg tab
pairsDg = filter_cells(predictedDgFile, numReads)
dgMat,errMat = make_predictedDg_matrix(pairsDg, cellAb)

outfile = numReads.split('.numReads.txt')[0]
dgMat.to_csv(outfile + '.dgMat.csv')
errMat.to_csv(outfile +'.errMat.csv')
