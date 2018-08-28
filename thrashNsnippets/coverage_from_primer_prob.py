import pandas as pd
import numpy as np
from primerProbability import *
from hexVSprimed import *


def compute_real_Dg(ptDgFile, log=True):
    '''
    Takes out the primer concentration factor from the Delta G parameter
    Input: file with predicted DGs at even primer concentrations
    '''
    dgTab = pd.read_csv(ptDgFile, index_col=0)
    probs = get_even_prob()
    realDg = pd.DataFrame()
    for col in dgTab.iteritems():
        primer,deltaGs = col
        print('Processing primer '+primer)
        primerProb = probs.loc[primer][0]
        dG = deltaGs/primerProb
        realDg = pd.concat([realDg, dG], axis=1)
    if log:
        realDg = np.log(realDg)
        realDg[realDg == - np.inf] = -99999
    return(realDg)

def coverage_from_primer_prob(dgTab, abundanceTab):
    '''

    '''


compute_real_Dg('/hpc/hub_oudenaarden/edann/hexamers/bootstrap/VAN1667_100_ptDg_qual.csv.gz', log=True)

tabAb = pd.read_csv(abundanceFile, index_col=0, compression=findCompr(abundanceFile), header=None)
genomeAb = tabAb[1] # Because it needs to be a series
