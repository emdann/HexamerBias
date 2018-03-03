import fnmatch
import os
import pandas as pd
import numpy as np
import multiprocessing
import argparse

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get df of pt-pairs and predicted DeltaG values (run on 10 cores) \n By Emma Dann")
argparser.add_argument('path', type=str, help='Path for folder with ptDg files to merge')
argparser.add_argument('-n', type=int, default=10, required=False, help='Minimum number of cells to retain pair.')
args = argparser.parse_args()

def makeNonInfDic(file):
    '''
    Makes a dictionary of pt pairs that have values != 0 in cell pt dg matrix
        Input: name of .csv file of predicted deltaG matrix
        Output: dictionary of {pt pair : {cellname : predicted deltaG}}
    '''
    cell = file.split('_cell')[1].split('_')[0]
    if file.endswith('gz'):
        compr='gzip'
    elif file.endswith('csv'):
        compr='infer'
    tab = pd.read_csv(file, index_col=0, compression=compr)
    cellDic={}
    for templ in tab.iterrows():
        t,p = templ
        for i in range(len(p)):
            if p[i] != 0:
                cellDic[t+'-'+p.index[i]]={}
                cellDic[t+'-'+p.index[i]][cell] = p[i]
    return(cellDic)

path = args.path
ncells = args.n

workers = multiprocessing.Pool(10)

finalCellDic = {}
for dic in workers.imap_unordered(makeNonInfDic, [ path + '/' + file for file in os.listdir(path)]):
    for pt,celdic in dic.items():
        if pt not in finalCellDic.keys():
            finalCellDic[pt]={}
        finalCellDic[pt] = dict(finalCellDic[pt],**celdic)

filtDic={}
for k,v in finalCellDic.items():
    if len(v) >= ncells:
        filtDic[k] = v

sample = os.listdir(path)[0].split('_')[0]
output = path.strip('predictedDg') + sample + '.commonPtPairs.csv'
pd.DataFrame.from_dict(filtDic).T.to_csv(output)
