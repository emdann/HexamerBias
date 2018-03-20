import pandas as pd
import numpy as np
import os
from hexVSprimed import *
from avg_deltaG_bs import *
from predictCovBs import *
from make_coverage_tab import make_coverage_tab


def make_cov_prediction_tab_avg(dir):
    '''
    Make table of predicted hexamer coverage using avg of delta G tables in dir
    '''
    coverageTab = make_coverage_tab(dir)
    colnames,pairNames,dgArrays = make_Dg_array(dir)
    cov = coverageTab.sum().sort_values(ascending=False)
    f = [f for f in os.listdir(dir) if 'ptDg' in f][0]
    df = pd.read_csv(dir + '/' + f, index_col=0, compression=findCompr(f))
    tabAb = pd.read_csv(abundanceTab, index_col=0, compression=findCompr(abundanceTab), header=None)
    genomeAb = tabAb[1]
    allPredictions = pd.DataFrame({'template':df.index})
    end=2
    while end < len(dgArrays):
        samples = [s.split('_bismark')[0] for s in cov.index[0:end]]
        selectedArrays = subset_dgArrays(dgArrays, samples)
        avg = compute_avg_deltaG(selectedArrays, df)
        expCov = predictCoverage(avg, genomeAb)
        expCov = expCov.rename(columns={'exp':'top' + str(end)})
        allPredictions = pd.merge(allPredictions, expCov.drop('err', axis=1), on='template')
        end +=1
    return(allPredictions)


def make_cov_prediction_tab(dir):
    '''
    Make table of predicted hexamer coverage using delta G tables in dir
    '''
    expCov = predictCoverage(avg, genomeAb)
