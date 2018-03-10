import fnmatch
import os
import pandas as pd
import numpy as np
import multiprocessing
from hexVSprimed import *

def predictCov(t,DgRow):
    '''
    Computes the predicted coverage for random hexamer experiment.
    '''
    sumChi = sum(np.exp(DgRow))
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

def findCompr(filename):
    if filename.endswith('gz'):
        compr='gzip'
    else:
        compr='infer'
    return(compr)
