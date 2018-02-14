import fnmatch
import os
import pandas as pd
import numpy as np
import multiprocessing
import argparse

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get matrix of predicted dg for primer-template complex in single cells \n By Emma Dann")
argparser.add_argument('sample', type=str, help='Pattern for ptDg files to merge')
args = argparser.parse_args()

pattern = args.sample
path='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/'

def makeNonInfDic(file):
    '''
    Makes a dictionary of pt pairs that have values != -99999 in cell pt dg matrix
        Input: name of .csv file of predicted deltaG matrix
        Output: dictionary of {pt pair : {cellname : predicted deltaG}}
    '''
    cell = file.split('_')[0][4:]
    if file.endswith('gz'):
        compr='gzip'
    elif file.endswith('csv'):
        compr='infer'
    tab = pd.read_csv(path+file, index_col=0, compression=compr)
    cellDic={}
    for templ in tab.iterrows():
        t,p = templ
        for i in range(len(p)):
            if p[i] != -99999:
                cellDic[t+'-'+p.index[i]]={}
                cellDic[t+'-'+p.index[i]][cell] = p[i]
    return(cellDic)

workers = multiprocessing.Pool(8)
finalCellDic = {}

for dic in workers.imap_unordered(makeNonInfDic, [ file for file in os.listdir(path) if fnmatch.fnmatch(file, pattern+'*ptDg_qual*')]):
    for pt,celdic in dic.items():
        if pt not in finalCellDic.keys():
            finalCellDic[pt]={}
        finalCellDic[pt] = dict(finalCellDic[pt],**celdic)

filtDic={}
for k,v in finalCellDic.items():
    if len(v)>=10:
        filtDic[k] = v

output = path + pattern + 'commonPtPairs_qual_allCells_parallel.csv'
pd.DataFrame.from_dict(filtDic).T.to_csv(output)
