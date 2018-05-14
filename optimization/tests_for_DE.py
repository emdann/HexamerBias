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

scoreFile='test6_ctcf.DE.rho.txt'
matrixFile='test6_ctcf.DE.matrix.csv'
deltaGfile = '/hpc/hub_oudenaarden/edann/crypts_bs/VAN2408/CM1_tr2_R1_bismark_bt2_ptDg_qual.csv'
abundanceFile = "/hpc/hub_oudenaarden/edann/hexamers/genomes_kmers/mm10.kmerAbundance.csv"


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

def coverage_best_matrix(mats, dgMat, abundance):
    seqs = all_hexamers()
    finalMat=mats.iloc[-1]
    prob_vec=np.array(finalMat)
    ppm = from_vec_to_ppm(prob_vec)
    primer_prob = prob_from_ppm(ppm, seqs)
    coverage = predictCoverage_setProbs(dgMat, abundance[1], primer_prob)
    return(coverage)

kmerFile = '/hpc/hub_oudenaarden/edann/hexamers/strand_specific/CTCF.flank60.kmersTot.csv'
# outdir='/'.join(kmerFile.split('/')[:-1])

randomKmerFile = '/hpc/hub_oudenaarden/edann/hexamers/strand_specific/CTCF.flank60.randomize.kmersTot.csv'
