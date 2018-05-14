#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys
import argparse
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/optimization')
from primerProbability import *
from predictCovBs import *
from tests_for_DE import read_matrix_file
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage')
from cov_from_density import *

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Compute coverage fraction from optimization output")
argparser.add_argument('sample', type=str, help='sample name for optimization files')
args = argparser.parse_args()

def coverage_best_matrix(mats, dgMat, abundance):
    seqs = all_hexamers()
    finalMat=mats.iloc[-1]
    prob_vec=np.array(finalMat)
    ppm = from_vec_to_ppm(prob_vec)
    primer_prob = prob_from_ppm(ppm, seqs)
    coverage = predictCoverage_setProbs(dgMat, abundance[1], primer_prob)
    return(coverage)

matrixFile = args.sample+'.DE.matrix.csv'
deltaGfile = '/hpc/hub_oudenaarden/edann/crypts_bs/VAN2408/CM1_tr2_R1_bismark_bt2_ptDg_qual.csv'
abundanceFile = "/hpc/hub_oudenaarden/edann/hexamers/genomes_kmers/mm10.kmerAbundance.csv"

dgMat = pd.read_csv(deltaGfile, index_col=0)
abundance = pd.read_csv(abundanceFile, index_col=0, compression=findCompr(abundanceFile), header=None)

coverage_best_matrix(read_matrix_file(matrixFile), dgMat, abundance).to_csv(sample+'.bestMat.coverage.csv')
