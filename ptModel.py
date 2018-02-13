from hexVSprimed import *
import pysam as ps
import collections
import argparse
import numpy as np
import pandas as pd
import multiprocessing
import gzip

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get matrix of predicted dg for primer-template complex in single cells \n By Emma Dann")
argparser.add_argument('ptmatrix', type=str, help='Bam input')
argparser.add_argument('cellabcsv', type=str, help='Fasta input')
args = argparser.parse_args()

def extractDeltaG(templateRow,tempAb):
    '''
    Extract predicted [log(p) + DeltaG] for row of ptCount table (one template)
    ...
    '''
    dg = np.log(templateRow/tempAb)
    dg[dg == - np.inf] = -99999
    return(dg)

def cellDgMat(params):
    '''
    Make matrix of predicted dg for p-t couples
    Input: tab of template abundace for cell OI, matrix of pt occurrencies
    '''
    cellAb,ptMat = params
    dgMat=pd.DataFrame()
    for temp in cellAb.index:
        temprow = ptMat[ptMat.index==temp]
        tempAb=cellAb[temp]
        dg = extractDeltaG(temprow,tempAb)
        dgMat = dgMat.append(dg)
    return(dgMat)

ptMatrix = args.ptmatrix
cellAbundanceTab = args.cellabcsv
# ptMatrix='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/cell121ptCounts.csv.gz'
# cellAbundanceTab='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/gk2a-2.cellAbundance.noN.csv'
cell = ptMatrix.split('/')[-1].split('ptCounts')[0].split('cell')[-1]
tabAb=pd.read_csv(cellAbundanceTab, index_col=0)
cellAb = tabAb[cell]

if ptMatrix.endswith('gz'):
    compr='gzip'
elif ptMatrix.endswith('csv'):
    compr='infer'

# Order pt matrix and add missing values
ptMat = pd.read_csv(ptMatrix, compression=compr, index_col=0)
ptMat = ptMat[[i for i in ptMat.columns if 'N' not in i]]
for temp in [i for i in cellAb.index if i not in ptMat.index]:
    newRow=pd.DataFrame(0, index=[temp], columns=ptMat.columns)
    ptMat = ptMat.append(newRow)

for primer in [i for i in cellAb.index if i not in ptMat.columns]:
    newCol=pd.DataFrame(0.0, index=[primer], columns=ptMat.index).T
    ptMat = pd.concat([ptMat, newCol], axis=1)

ptMat=ptMat.sort_index(axis=1).sort_index()


dgMat = cellDgMat((cellAb, ptMat))

outpath = '/'.join(ptMatrix.split('/')[:-1])
if outpath:
    dgMat.to_csv(outpath+'/cell'+cell+'_ptDg_qual.csv')
else:
    dgMat.to_csv(outpath+'cell'+cell+'_ptDg_qual.csv')
