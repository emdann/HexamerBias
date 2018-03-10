import fnmatch
import os
import pandas as pd
import numpy as np
import multiprocessing
from hexVSprimed import *

def makeNonInfDic(file):
    '''
    Makes a dictionary of pt pairs that have values != 0 in cell pt dg matrix
        Input: name of .csv file of predicted deltaG matrix
        Output: dictionary of {pt pair : {cellname : predicted deltaG}}
    '''
    cell = file.split('_cell')[1].split('_')[0]
    tab = pd.read_csv(file, index_col=0, compression=findCompr(file))
    cellDic={}
    for templ in tab.iterrows():
        t,p = templ
        for i in range(len(p)):
            if p[i] != 0:
                cellDic[t+'-'+p.index[i]]={}
                cellDic[t+'-'+p.index[i]][cell] = p[i]
    return(cellDic)
