import fnmatch
import os
import pandas as pd
import numpy as np
import argparse
from hexVSprimed import *

# argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Coverage prediction. \n By Emma Dann")
# argparser.add_argument('predDg', type=str, help='Matrix of predicted Dg')
# argparser.add_argument('predSd', type=str, help='Matrix of predicted Dg standard deviation')
# argparser.add_argument('abundance', type=str, help='template abundance in refgen')
# # argparser.add_argument('cellPtCount', type=str, help='Matrix of pt occurrencies in a cell')
# args = argparser.parse_args()

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

def predictCoverage(dgMat, genomeAb, errDgMat=None):
    predictedCov = pd.DataFrame(columns=['template', 'exp', 'err'])
    for templ in dgMat.iterrows():
        t,DgRow = templ
        cov = predictCov(genomeAb[t],DgRow)
        if errDgMat:
            err = propagateError(genomeAb[t],DgRow, errDgMat[t])
        else:
            err = np.nan
        rowdf = pd.DataFrame.from_dict({'template':[t], 'exp':[cov], 'err':[err]})
        predictedCov = predictedCov.append(rowdf)
    return(predictedCov)

# ptMatrix = args.cellPtCount
# predictedDg = args.predDg
# predictedSd = args.predSd
# sample = ptMatrix.split('.')[0]
# cellAbundanceTab = args.abundance

def predict_fromFiles(predictedDg, abundanceTab, predictedSd=None):
    tabAb = pd.read_csv(abundanceTab, index_col=0, compression=findCompr(cellAbundanceTab), header=None)
    genomeAb = tabAb[1]
    dgMat = pd.read_csv(predictedDg, index_col=0, compression = findCompr(predictedDg))
    if predictedSd:
        errDgMat = pd.read_csv(predictedSd, compression = findCompr(predictedSd), index_col=0)
    return(predictCoverage(dgMat, genomeAb, errMat=errDgMat))

# with open('predictedCoverage_avgzf.txt', 'w') as output:
#     # print('template', 'exp','err',  sep='\t', file=output)
#     print('template', 'exp',  sep='\t', file=output)
#     for templ in dgMat.iterrows():
#         t,DgRow = templ
#         # print(t, predictCov(genomeAb[t],DgRow), propagateError(genomeAb[t],DgRow, errDgMat[t]), sep='\t', file=output)
#         print(t, predictCov(genomeAb[t],DgRow), sep='\t', file=output)

# ptMat.loc[ptMat.index==t].fillna(0).values.sum()
