import pysam as ps
import os
import pandas as pd
import numpy as np
from hexVSprimed import *

def make_Dg_array(dir):
    '''
    Make array of all predicted DeltaF values for all the samples in dir.
    Transforms to log and masks the nan values.
    Returns sample names, pairs names, array of arrays
    '''
    files = [f for f in os.listdir(dir) if 'ptDg' in f]
    dgArrays = []
    for f in files:
        df = pd.read_csv(dir + '/' + f, index_col=0, compression=findCompr(f))
        array = np.ma.log(np.array(df)).reshape((16777216))
        dgArrays.append(array)
    pairNames = []
    for row in df.index:
        for p in df.columns:
            pairNames.append(row+'-'+p)
    colnames = [s.split('_bismark')[0] for s in files]
    return((colnames,pairNames, dgArrays))

def save_common_pt_pairs(dgArrays,pairNames, colnames, samplename='zebrafishBs'):
    '''
    Find instances that have at least one value != NaN and save them in csv
    '''
    with open(samplename+'.commonPtPairs.csv', 'w') as outfile:
        print(','.join(['']+colnames), file=outfile)
        for i in range(len(dgArrays[0])):
            pair = np.array([dg[i] for dg in dgArrays])
            if sum(np.isnan(pair))<5:
                print(','.join([pairNames[i]]+[str(p) for p in pair]), file=outfile)
                # commonPairs[pairNames[i]]=(n for n in pair) # <-- using generator instead of list comprehension
    return(None)

# def save_common_pt_pairs(dgArrays, pairNames, colnames, samplename='zebrafishBs', to_csv=True):
#     pairs = pd.DataFrame(commonPairs).T
#         pairs.columns = colnames
#         if to_csv:
#             pairs.to_csv(samplename + '.commonPtPairs.csv')
#         return(pairs)

dir='.'
colnames,pairNames,dgArrays = make_Dg_array(dir)

# Correlation matrix
def make_deltaG_correlation_matrix(commonPtFiles):
    '''
    Input: list of commonPtFiles from different experiments/cells
    Output: correlation matrix of all the values
    '''
    allPairs = pd.read_csv(commonPtFiles[0], index_col=0)
    for file in commonPtFiles[1:]:
        commonPt = pd.read_csv(file, index_col=0)
        allPairs = pd.concat([allPairs, commonPt], axis=1)
    corrMat = allPairs.corr()
    return(corrMat)


# Average and standard deviation
def compute_std_deltaG(dgArrays, colnames, templateDf):
    '''
    Compute std on selected colnames.
    Takes in input original ptDg table as template
    '''
    maskedArray = np.ma.array([f for f in dgArrays])
    sdMatrix = maskedArray.std(axis=0)
    stddf = pd.DataFrame(sd.reshape(4096,4096))
    stddf.columns = templateDf.columns
    stddf.index = templateDf.index
    return(stddf)

# s=0
# for a in dgArrays:
#     s+=a[-1]
#
#
# print(s/len(dgArrays))
#
# avg = np.array([])
# for i in range(len(dgArrays[0])):
#     values = [a[i] for a in dgArrays]
#     mean = np.nanmean(values)
#     avg = np.append(avg,mean)
#
#
