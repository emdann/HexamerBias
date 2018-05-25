from hexVSprimed import *
import pysam as ps
import collections
import argparse
import numpy as np
import pandas as pd
import multiprocessing
import gzip
import os

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get matrix of predicted dg for primer-template complex \n By Emma Dann")
argparser.add_argument('ptmatrix', type=str, help='ptCounts file')
argparser.add_argument('cellabcsv', type=str, help='Kmer abundance file')
argparser.add_argument('--filt', type=int, default=0, help='filter out counts under filt')
argparser.add_argument('--suff', type=str,default="", help='additional suffix to output name')
args = argparser.parse_args()

def extract_deltaG(templateRow,tempAb, totKmers, S):
    '''
    Extract predicted p*exp(DeltaG) for row of ptCount table (one template)
        templateRow: row of primer template counts for a certain template
        tempAb: genomic abundance of template
    '''
    dg = (templateRow/S)/(((tempAb - templateRow.values.sum())/totKmers)*0.00024) # Takes off primer concentration from the parameter
    # dg[dg == - np.inf] = -9999
    return(dg)

def make_DgMat_per_cell(cellAb,ptMat,S):
    '''
    Make matrix of predicted dg for p-t couples, scaled by tot number of reads
    Input:
        tab of template abundance in genome,
        matrix of pt occurrencies,
        scaling factor (no. of reads)
    '''
    totTempl = cellAb.sum()
    dgMat=pd.DataFrame()
    for temp in cellAb.index:
        temprow = ptMat[ptMat.index==temp]
        temprow = temprow.fillna(0)
        tempAb=cellAb[temp]
        dg = extract_deltaG(temprow,tempAb, totTempl, S)
        dgMat = dgMat.append(dg)
    return(dgMat)

def filter_lowCounts(ptMat,t):
    '''
    Change to zero if pt count is less than t
    '''
    filtPtMat=ptMat.copy()
    filtPtMat[filtPtMat<t]=0
    return(filtPtMat)

ptMatrix = args.ptmatrix
# ptMatrix='ptCounts/SvdB11d1-MitoTrackerThird-Satellites-Adult.cell130.ptCounts.qualFilt.parallel.csv'
cellAbundanceTab = args.cellabcsv
filt=args.filt

# if type=='bs':
sample = ptMatrix.split('/')[-1].split('.ptCounts')[0]
tabAb = pd.read_csv(cellAbundanceTab, index_col=0, compression=findCompr(cellAbundanceTab), header=None)
genomeAb = tabAb[1]
ptMat = pd.read_csv(ptMatrix, compression=findCompr(ptMatrix), index_col=0)
filtPtMat = filter_lowCounts(ptMat, filt)
dgMat = make_DgMat_per_cell(genomeAb, filtPtMat, filtPtMat.sum().sum())
path = '/'.join(ptMatrix.split('/')[:-1])
dgMat.to_csv(path +sample +'_ptDg_qual' + args.suff + '.csv')
