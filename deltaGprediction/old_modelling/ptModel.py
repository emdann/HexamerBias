from hexVSprimed import *
import sys
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/optimization')
from primerProbability import *
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
argparser.add_argument('--ppm', type=str, default=None, help='ppm of primer composition')
argparser.add_argument('--filt', type=int, default=0, help='filter out counts under filt')
argparser.add_argument('--suff', type=str,default="", help='additional suffix to output name')
args = argparser.parse_args()

def extract_deltaG(templateRow,tempAb, totKmers, S):
    '''
    Extract predicted p*exp(DeltaG) for row of ptCount table (one template)
        templateRow: row of primer template counts for a certain template
        tempAb: genomic abundance of template
    '''
    dg = (templateRow/S)/(((tempAb - templateRow.values.sum())/totKmers))
    # dg[dg == - np.inf] = -9999
    return(dg)

# def extract_deltaG_primer_conc(templateRow,tempAb, totKmers, primer_prob,S):
#     '''
#     Extract predicted p*exp(DeltaG) for row of ptCount table (one template)
#         templateRow: row of primer template counts for a certain template
#         tempAb: genomic abundance of template
#     '''
#     dg = (templateRow/S)/(((tempAb - templateRow.values.sum())/totKmers)) # Takes off primer concentration from the parameter
#     probDg  = pd.DataFrame(pd.concat([dg.T, 1/primerProb], axis=1).apply(np.prod, axis=1), columns=[templateRow.index]).T
#     # dg[dg == - np.inf] = -9999
#     return(probDg)

def divide_primer_prob(dg,primerProb):
    pp = primerProb.iloc[:,0][dg.columns]
    dgWpp = dg.divide(pp)
    return(dgWpp)

def make_DgMat_per_cell(cellAb,ptMat,S, wPrimer=False, primerProb=None):
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
    if wPrimer:
        dgWpp = divide_primer_prob(dgMat,primerProb)
    else:
        dgWpp = dgMat/0.000244
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
ppmFile = args.ppm
if ppmFile:
    wPrimer=True
else:
    wPrimer=False

# if type=='bs':
sample = ptMatrix.split('/')[-1].split('.ptCounts')[0]
tabAb = pd.read_csv(cellAbundanceTab, index_col=0, compression=findCompr(cellAbundanceTab), header=None)
genomeAb = tabAb[1]
ptMat = pd.read_csv(ptMatrix, compression=findCompr(ptMatrix), index_col=0)
filtPtMat = filter_lowCounts(ptMat, filt)

ppm = pd.read_csv(ppmFile, index_col=0, compression=findCompr(ppmFile))
# prob={'A':[0.25,0.25,0.25,0.25,0.25,0.25],'T':[0.05,0.05,0.05,0.05,0.05,0.05], 'C':[0.25,0.25,0.25,0.25,0.25,0.25], 'G':[0.45,0.45,0.45,0.45,0.45,0.45]}
# ppm = pd.DataFrame(prob).T
ppm.columns = [int(i) for i in ppm.columns]   # Make the colnames integers for the calling in prob_from_ppm
ppm = ppm.loc[['A','T', 'C', 'G'],:]

dgMat = make_DgMat_per_cell(genomeAb, filtPtMat, filtPtMat.sum().sum(), wPrimer=wPrimer, primerProb=prob_from_ppm(ppm, all_hexamers()).loc[ptMat.columns])
path = '/'.join(ptMatrix.split('/')[:-1])
dgMat.to_csv(path +sample +'_ptDg_qual' + args.suff + '.csv')
