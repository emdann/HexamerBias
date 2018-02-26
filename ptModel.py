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

def extract_deltaG(templateRow,tempAb):
    '''
    Extract predicted p* exp(DeltaG) for row of ptCount table (one template)
    ...
    '''
    dg = templateRow/(tempAb - templateRow.values.sum())
    # dg[dg == - np.inf] = -99999
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

ptMatrix = args.ptmatrix
# ptMatrix='ptCounts/SvdB11d1-MitoTrackerThird-Satellites-Adult.cell130.ptCounts.qualFilt.parallel.csv'
cellAbundanceTab = args.cellabcsv

sample = cellAbundanceTab.split('/')[-1].split('.cellAbundance')[0]
cell = ptMatrix.split('/')[-1].split('.ptCounts')[0].split('cell')[-1]
tabAb = pd.read_csv(cellAbundanceTab, index_col=0, compression=findCompr(cellAbundanceTab))
cellAb = tabAb[cell]
ptMat = pd.read_csv(ptMatrix, compression=findCompr(ptMatrix), index_col=0)

dgMat = make_DgMat_per_cell((cellAb, ptMat))

path = '/'.join(cellAbundanceTab.split('/')[:-1])
if path:
    outpath = path + '/predictedDg/'
else:
    outpath = './predictedDg/'

dgMat.to_csv(outpath + sample + '_cell'+ cell +'_ptDg_qual.csv')
