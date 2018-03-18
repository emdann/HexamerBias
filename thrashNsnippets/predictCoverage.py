import fnmatch
import os
import pandas as pd
import numpy as np
import multiprocessing
import argparse
from hexVSprimed import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Coverage prediction. \n By Emma Dann")
argparser.add_argument('predDg', type=str, help='Matrix of predicted Dg')
argparser.add_argument('predSd', type=str, help='Matrix of predicted Dg standard deviation')
argparser.add_argument('cellPtCount', type=str, help='Matrix of pt occurrencies in a cell')
argparser.add_argument('-t', type=str, default=0, help='type (rna or bs)')
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
    sumChi = sum(np.exp(DgRow))
    errChi = errRow.multiply(np.exp(DgRow))
    beta = (1 - 2 * sumChi)/np.square(1 - sumChi)
    error = np.sqrt(sum(beta*errChi)*t)
    return(error)

def setThresh4Dg(dgMat,ptMat,thresh=1):
    '''
    Assuming that if a pt interaction is seen only once then it is an error, the predicted Dg value
    for a pt couple with conc = 1 is set to -99999
    '''
    dgMat[ptMat <= thresh ] = -0
    return(dgMat)


ptMatrix = args.cellPtCount
predictedDg = args.predDg
errDg = predictedDg.split('dgMat.csv')[0]+'errMat.csv'
sample = predictedDg.split('.')[0]
type = args.t

# for g in dgMat.iterrows():
#     temp,DgRow=g
#     t=cellAb['99'][temp]
#     print(temp,predictCov(t,DgRow), propagateError(t,DgRow,errMat.loc[temp]), sep='\t')



if type=='rna':
    cellAbundanceTab = sample + '.cellAbundance.noN.csv'
    cell = ptMatrix.split('.ptCounts')[0].split('cell')[-1]
    path = '/'.join(cellAbundanceTab.split('/')[:-1])
    tabAb = pd.read_csv(cellAbundanceTab, index_col=0, compression = findCompr(cellAbundanceTab))
    cellAb = tabAb[cell]

    ptMat = pd.read_csv(ptMatrix, compression = findCompr(ptMatrix), index_col=0)

    errDgMat = pd.read_csv(errDg, compression = findCompr(errDg), index_col=0)
    dgMat = pd.read_csv(path+predictedDg, index_col=0, compression = findCompr(predictedDg))
    # if thresh:
    #     dgMat = setThresh4Dg(dgMat,ptMat,thresh=thresh)


    list=[]
    with open(path+'predictedCov/'+sample+'.CovPred.'+str(cell)+'.qual.txt', 'w') as output:
        print('template','obs', 'exp', 'err', sep='\t') #, file=output)
        for templ in dgMat.iterrows():
            t,DgRow = templ
            print(t, ptMat.loc[ptMat.index==t].fillna(0).values.sum(), predictCov(cellAb[t],DgRow), propagateError(cellAb[t], DgRow, errDgMat.loc[t]), sep='\t', file=output)

if type=='bs':
cellAbundanceTab = '/hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv'
tabAb = pd.read_csv(cellAbundanceTab, index_col=0, compression=findCompr(cellAbundanceTab), header=None)
genomeAb = tabAb[1]
ptMat = pd.read_csv(ptMatrix, compression=findCompr(ptMatrix), index_col=0)
errDgMat = pd.read_csv(errDg, compression = findCompr(errDg), index_col=0)
path = '/'.join(ptMatrix.split('/')[:-1])
dgMat = pd.read_csv(path+predictedDg, index_col=0, compression = findCompr(predictedDg))

with open('OUD2086_distal.CovPred.qual.txt', 'w') as output:
    print('template','obs', 'exp',  sep='\t', file=output)
    for templ in dgMat.iterrows():
        t,DgRow = templ
        print(t, ptMat.loc[ptMat.index==t].fillna(0).values.sum(), predictCov(genomeAb[t],DgRow), sep='\t', file=output) #, propagateError(cellAb[t], errDgMat[t]), sep='\t')#, file=output)
