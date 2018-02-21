from hexVSprimed import *
import pysam as ps
import collections
import argparse
import numpy as np
import pandas as pd
import multiprocessing
import gzip
import os

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get matrix of predicted dg for primer-template complex in single cells \n By Emma Dann")
argparser.add_argument('ptmatrix', type=str, help='Bam input')
argparser.add_argument('cellabcsv', type=str, help='Fasta input')
args = argparser.parse_args()

def extractDeltaG(templateRow,tempAb):
    '''
    Extract predicted [log(p) + DeltaG] for row of ptCount table (one template)
    ...
    '''
    dg = np.log(templateRow/(tempAb - templateRow.values.sum()))
    dg[dg == - np.inf] = -99999
    return(dg)

def cellDgMat(params):
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
        dg = extractDeltaG(temprow,tempAb)
        dgMat = dgMat.append(dg)
    return(dgMat)

ptMatrix = args.ptmatrix
cellAbundanceTab = args.cellabcsv
cell = ptMatrix.split('/')[-1].split('ptCounts')[0].split('cell')[-1]
tabAb = pd.read_csv(cellAbundanceTab, index_col=0)
cellAb = tabAb[cell]
cellAb = cellAb[[i for i in cellAb.index if 'N' not in i]]

# Order pt matrix and add missing values
ptMat = pd.read_csv(ptMatrix, compression=findCompr(ptMatrix), index_col=0)
ptMat = fillNsortPTmatrix(ptMat, cellAb)

dgMat = cellDgMat((cellAb, ptMat))

outpath = '/'.join(fasta.split('/')[:-1]) + '/predictedDg/'
if not os.path.exists(outpath):
    os.makedirs(outpath)

sample = cellAbundanceTab.split('/')[-1].split('.cellAbundance')[0]

if outpath:
    dgMat.to_csv(outpath+'/'+ sample + '_cell'+ cell +'_ptDg_qual.csv')
else:
    dgMat.to_csv(outpath+ sample +'_cell'+cell+'_ptDg_qual.csv')
