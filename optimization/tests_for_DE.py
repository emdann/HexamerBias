import numpy as np
import multiprocessing
import json
import sys
import argparse
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/optimization')
from primerProbability import *
from predictCovBs import *
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage')
from cov_from_density import *

scoreFile='test4_ctcf.DE.rho.txt'
matrixFile='test4_ctcf.DE.matrix.csv'

def read_score_file(scoreFile):
    scores=[]
    with open(scoreFile, 'r') as f:
        for s in f.readlines():
            scores.append(float(s.rstrip()))
    return(scores)

def read_matrix_file(matrixFile):
    matTab = pd.read_csv(matrixFile, index_col=0)
    return(matTab)

probVec = np.array(mats.iloc[1])

def test_DE_output(mats, scores, deltaGfile, abundanceFile, kmerFile):
    seqs = all_hexamers()
    dgMat = pd.read_csv(deltaGfile, index_col=0)
    abundance = pd.read_csv(abundanceFile, index_col=0, compression=findCompr(abundanceFile), header=None)
    targetKmers = pd.read_csv(kmerFile, index_col=0, header=None)
    for row in mats.iterrows():
        it,probVec = row
        newScore = new_coverage_function((np.array(probVec),seqs,dgMat,abundance,targetKmers))
        if newScore==scores[it]:
            print('Right score!')
        else:
            print("Wrong score: {new} instead of {old}".format(new=newScore, old=scores[it]))
