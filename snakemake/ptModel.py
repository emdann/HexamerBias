from hexVSprimed import *
import pysam as ps
import collections
import argparse
import numpy as np
import pandas as pd
import multiprocessing
import gzip
import os

def extract_deltaG(templateRow,tempAb):
    '''
    Extract predicted p*exp(DeltaG) for row of ptCount table (one template)
    ...
    '''
    dg = templateRow/(tempAb - templateRow.values.sum())
    # dg[dg == - np.inf] = -9999
    return(dg)

def make_DgMat_per_cell(params):
    '''
    Make matrix of predicted dg for p-t couples
    Input: tab of template abundance for cell OI, matrix of pt occurrencies
    '''
    cellAb,ptMat = params
    dgMat=pd.DataFrame()
    for temp in cellAb.index:
        temprow = ptMat[ptMat.index==temp]
        temprow = temprow.fillna(0)
        tempAb=cellAb[temp]
        dg = extract_deltaG(temprow,tempAb)
        dgMat = dgMat.append(dg)
    return(dgMat)
