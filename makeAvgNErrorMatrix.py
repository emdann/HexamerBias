import fnmatch
import os
import pandas as pd
import numpy as np
import argparse
from hexVSprimed import *

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Find average predicted Delta G from common pt pairs\n By Emma Dann")
argparser.add_argument('commonPtPairs', type=str, help='File of predicted Dg for common PtPairs')
argparser.add_argument('noReads', type=str, help='File of predicted Dg for common PtPairs')
argparser.add_argument('cellab', type=str, help='name of output file')
argparser.add_argument('output', type=str, help='suffix for output files')
args = argparser.parse_args()

def filterCells(predictedDgFile, noReads):
    '''
    Finds average of predicted Dg and standard deviation using a selection of cells
    '''
    nReads = pd.read_csv(noReads, sep=' ', names=['noReads', 'cell'])
    commonPairs = pd.read_csv(predictedDgFile, sep=',', index_col=0)
    goodCells = [str(cell.cell) for index,cell in nReads.iterrows() if cell.noReads>=40000]
    predictedDgAvg = {}
    for pair,values in commonPairs[goodCells].iterrows():
        mean,sd = np.mean(values),np.std(values)
        predictedDgAvg[pair]={'dg':mean,'sd':sd}
    return(pd.DataFrame(predictedDgAvg).T)

def makePredictedDgMatrix(df, cellAb):
    '''
    Turns file of pt - predicted Dg (precessed in R) to matrix of template on row and primer on column
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


cellAb = pd.read_csv(args.cellab, index_col=0, compression = findCompr(args.cellab))
cellAb = cellAb.loc[[i for i in cellAb.index if 'N' not in i]]

predictedDgFile = args.commonPtPairs
noReads = args.no.reads
# Make predicted Dg tab
pairsDg = filterCells(predictedDgFile, noReads)
dgMat,errMat = makePredictedDgMatrix(pairsDg, cellAb)
dgMat.to_csv(args.output+'dgMat.csv')
errMat.to_csv(args.output+'errMat.csv')
