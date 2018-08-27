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
argparser.add_argument('dg', type=str, help='DeltaG matrix file')
argparser.add_argument('genome', type=str, help='genome for abundance file (mm10, hg38, WBcel235, danRer10)')
argparser.add_argument('--ppm', type=str, default=None, help='ppm of primer composition')
args = argparser.parse_args()

deltaGfile = args.dg
ppmFile = args.ppm
abundanceFile = "/hpc/hub_oudenaarden/edann/hexamers/genomes_kmers/"+args.genome+".kmerAbundance.csv"

dgMat = pd.read_csv(deltaGfile, index_col=0)
abundance = pd.read_csv(abundanceFile, index_col=0, compression=findCompr(abundanceFile), header=None)


if ppmFile:
    ppm2 = pd.read_csv(ppmFile, index_col=0)
    ppm2.columns = [int(i) for i in ppm2.columns]   # Make the colnames integers for the calling in prob_from_ppm
    ppm = ppm2.loc[['A','T', 'C', 'G'],:]   # Right order of base pairs
    ppmName= ppmFile.rstrip(".csv")
else:
    prob_vec = np.array([[0.25,0.25,0.25, 0.25, 0.25, 0.25],
                [0.25,0.25,0.25, 0.25, 0.25, 0.25],
                [0.25,0.25,0.25, 0.25, 0.25, 0.25],
                [0.25,0.25,0.25, 0.25, 0.25, 0.25]  ])
    ppm = from_vec_to_ppm(prob_vec)
    ppmName = 'R_ppm'

primer_prob = prob_from_ppm(ppm, all_hexamers()).loc[dgMat.columns]
coverage = predictCoverage_setProbs(dgMat, abundance[1], primer_prob)
coverage.to_csv(deltaGfile.rstrip('_ptDg_qual.csv') + '_' +ppmName + '.predcoverage.csv')
