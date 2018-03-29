import fnmatch
import os
import pandas as pd
import numpy as np
import argparse
import multiprocessing
from hexVSprimed import *

# argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Coverage prediction. \n By Emma Dann")
# argparser.add_argument('predDg', type=str, help='Matrix of predicted Dg')
# argparser.add_argument('predSd', type=str, help='Matrix of predicted Dg standard deviation')
# argparser.add_argument('abundance', type=str, help='template abundance in refgen')
# # argparser.add_argument('cellPtCount', type=str, help='Matrix of pt occurrencies in a cell')
# args = argparser.parse_args()

def predictCov(t,DgRow, primer=None):
    '''
    Computes the predicted coverage for random hexamer experiment.
    When using primer concentrations, use directly the exponential (dg tab without log)
    '''
    if primer is None:
        sumChi = sum(np.exp(DgRow.fillna(-9999999)))
    else:
        sumChi = sum(np.exp(DgRow).mul(primer.primerProb))
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

def predictCoverage(dgMat, genomeAb, errDgMat=None, primer_ppm=None):
    '''
    Takes as input the log(Dg)
    '''
    predictedCov = pd.DataFrame(columns=['template', 'exp', 'err'])
    if primer_ppm is not None:
        probs = prob_from_ppm(primer_ppm, list(genomeAb.index))
    else:
        probs=None
    for templ in dgMat.iterrows():
        t,DgRow = templ
        cov = predictCov(genomeAb[t],DgRow, primer=probs)
        if errDgMat:
            err = propagateError(genomeAb[t],DgRow, errDgMat[t])
        else:
            err = np.nan
        rowdf = pd.DataFrame.from_dict({'template':[t], 'exp':[cov], 'err':[err]})
        predictedCov = predictedCov.append(rowdf)
    return(predictedCov)

def predictCoverage_setProbs(dgMat, genomeAb, probs): # <--- This should be quicker
    '''
    Takes aspred input the log(Dg)
    '''
    predictedCov = pd.DataFrame(columns=['template', 'exp', 'err'])
    for templ in dgMat.iterrows():
        t,DgRow = templ
        cov = predictCov(genomeAb[t],DgRow, primer=probs)
        rowdf = pd.DataFrame.from_dict({'template':[t], 'exp':[cov]})
        predictedCov = predictedCov.append(rowdf)
    return(predictedCov)

### Making coverage prediction parallel
def predictCoverage_setProbs_parallel(dgMat, genomeAb, probs, cores=5): # <--- This should be quicker
    '''
    Takes aspred input the log(Dg)
    '''
    predictedCov = pd.DataFrame(columns=['template', 'exp', 'err'])
    for rowdf in workers.map(predict_row_coverage, [(t,dgMat) for t,dgMat in dgMat.iterrows()]):
        predictedCov = predictedCov.append(rowdf)
    return(predictedCov)

def predict_row_coverage(params):
    t,dgMat=params
    cov = predictCov(genomeAb[t],DgRow, primer=probs)
    rowdf = pd.DataFrame.from_dict({'template':[t], 'exp':[cov]})
    return(rowDf)

# ptMatrix = args.cellPtCount
# predictedDg = args.predDg
# predictedSd = args.predSd
# sample = ptMatrix.split('.')[0]
# cellAbundanceTab = args.abundance

def predict_fromFiles(predictedDg, abundanceTab, predictedSd=None):
    tabAb = pd.read_csv(abundanceTab, index_col=0, compression=findCompr(abundanceTab), header=None)
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
