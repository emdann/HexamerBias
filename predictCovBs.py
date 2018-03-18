import fnmatch
import os
import pandas as pd
import numpy as np
import argparse
from hexVSprimed import *

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Coverage prediction. \n By Emma Dann")
argparser.add_argument('predDg', type=str, help='Matrix of predicted Dg')
argparser.add_argument('predSd', type=str, help='Matrix of predicted Dg standard deviation')
argparser.add_argument('abundance', type=str, help='template abundance in refgen')
argparser.add_argument('cellPtCount', type=str, help='Matrix of pt occurrencies in a cell')
args = argparser.parse_args()

def predictCov(t,DgRow):
    '''
    Computes the predicted coverage for random hexamer experiment.
    '''
    sumChi = sum(np.exp(DgRow.fillna(-9999999)))
    cov = t * (sumChi/(1 + sumChi)) # <--- check
    return(cov)

def propagateError(t,DgRow,errRow):
    '''
    Propagation of error from standard deviation of avg deltaG
    '''
    sumChi = sum(np.exp(DgRow).fillna(0))
    errChi = errRow.multiply(np.exp(DgRow)).fillna(0)
    beta = (1 - 2 * sumChi)/np.square(1 - sumChi)
    error = np.sqrt(sum(beta*errChi)*t)
    return(error)

ptMatrix = args.cellPtCount
predictedDg = args.predDg
predictedSd = args.predSd
sample = ptMatrix.split('.')[0]
cellAbundanceTab = args.abundance

tabAb = pd.read_csv(cellAbundanceTab, index_col=0, compression=findCompr(cellAbundanceTab), header=None)
genomeAb = tabAb[1]
ptMat = pd.read_csv(ptMatrix, compression=findCompr(ptMatrix), index_col=0)
dgMat = pd.read_csv(predictedDg, index_col=0, compression = findCompr(predictedDg))
errDgMat = pd.read_csv(predictedSd, compression = findCompr(predictedSd), index_col=0)

with open(sample+'.CovPred.qual.txt', 'w') as output:
    print('template','obs', 'exp','err',  sep='\t', file=output)
    for templ in dgMat.iterrows():
        t,DgRow = templ
        print(t, ptMat.loc[ptMat.index==t].fillna(0).values.sum(), predictCov(genomeAb[t],DgRow), propagateError(genomeAb[t],DgRow, errDgMat[t]), sep='\t', file=output)
