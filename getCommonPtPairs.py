import fnmatch
import os
import pandas as pd
import numpy as np


path='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/'

bigDic={}
for file in os.listdir(path):
    if fnmatch.fnmatch(file, 'cell*ptDg*'):
        cell= file.split('_')[0][4:]
        if file.endswith('gz'):
            compr='gzip'
        elif file.endswith('csv'):
            compr='infer'
        tab=pd.read_csv(path+file, index_col=0, compression=compr)
        for templ in tab.iterrows():
            t,p = templ
            for i in range(len(p)):
                if p[i]!=-99999:
                    if t+'-'+p.index[i] not in bigDic.keys():
                        bigDic[t+'-'+p.index[i]]={}
                    bigDic[t+'-'+p.index[i]][cell]=p[i]

filtDic={}
for k,v in bigDic.items():
    if len(v)>=20:
        filtDic[k]=v

output = path + 'commonPtPairs_allCells.csv'
pd.DataFrame.from_dict(filtDic).to_csv(output)
