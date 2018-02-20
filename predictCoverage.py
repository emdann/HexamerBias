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

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get hexamers used in fasta file.\n By Emma Dann")
argparser.add_argument('predDg', type=str, help='Matrix of predicted Dg')
argparser.add_argument('cellPtCount', type=str, help='Matrix of pt occurrencies in a cell')
argparser.add_argument('-t', type=int, default=0, required=False, help='Threshold of pt events to predict Dg on')
args = argparser.parse_args()

def predictCov(t,DgRow):
    '''
    Computes the predicted coverage for random hexamer experiment.
    '''
    sumRow = sum(np.exp(DgRow).fillna(0))
    cov = t*sumRow
    return(cov)

def propagateError(t,ErrRow):
    '''
    Propagation of error from standard deviation of avg deltaG
    '''
    error = np.sqrt(sum(np.square(np.exp(ErrRow)).fillna(0)))
    return(error)

def setThresh4Dg(dgMat,ptMat,thresh=1):
    '''
    Assuming that if a pt interaction is seen only once then it is an error, the predicted Dg value
    for a pt couple with conc = 1 is set to -99999
    '''
    dgMat[ptMat <= thresh ] = -99999.0
    return(dgMat)

def findCompr(filename):
    if filename.endswith('gz'):
        compr='gzip'
    else:
        compr='infer'
    return(compr)

path='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/'

predictedDg = args.predDg
ptMatrix = args.cellPtCount
errMatrix = predictedDg.split('dgMat.csv')[0]+'errMat.csv'
sample = predictedDg.split('_')[0]
cellAbundanceTab = sample + '.cellAbundance.csv'
thresh=args.t


cell = ptMatrix.split('pt')[0].split('cell')[-1]

tabAb = pd.read_csv(path+cellAbundanceTab, index_col=0, compression = findCompr(cellAbundanceTab))
cellAb = tabAb[cell]
cellAb = cellAb[[i for i in cellAb.index if 'N' not in i]]

ptMat = pd.read_csv(path+ptMatrix, compression = findCompr(ptMatrix), index_col=0)
ptMat = ptMat[[i for i in ptMat.columns if 'N' not in i]]
ptMat = fillNsortPTmatrix(ptMat, cellAb)

errMat = pd.read_csv(path+errMatrix, compression = findCompr(errMatrix), index_col=0)
errMat = errMat[[i for i in errMat.columns if 'N' not in i]]
errMat = fillNsortPTmatrix(errMat, cellAb)


tab = pd.read_csv(path+predictedDg, index_col=0, compression = findCompr(predictedDg))
if thresh:
    tab=setThresh4Dg(tab,ptMat,thresh=thresh)

list=[]
with open(path+'predictedCov/'+sample+'.CovPred.'+str(cell)+'.thresh'+str(thresh)+'.qual.txt', 'w') as output:
    print('template','obs', 'exp', 'err' sep='\t') #, file=output)
for templ in tab.iterrows():
    t,DgRow = templ
    # list.append((sum(ptMat[t]), predictCov(cellAb[t],DgRow)))
    print(t, sum(ptMat[t].fillna(0)), predictCov(cellAb[t],DgRow), propagateError(cellAb[t], errMat[t]), sep='\t')#, file=output)
# obs = [x for x,y in list]
# exp = [y for x,y in list]
#
# plt.figure(1)
# plt.plot(obs,exp, 'b.')
# plt.xlabel('observed coverage')
# plt.ylabel('predicted coverage')
# plt.savefig('/hpc/hub_oudenaarden/edann/output/test_obsVSexp_cov_cell194.pdf')
