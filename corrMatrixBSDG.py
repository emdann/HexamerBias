import numpy as np
import pandas as pd
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dir='/hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/'

files = [f for f in os.listdir(dir) if 'ptDg' in f]
mats = {}
for file in files:
    sample = file.split('_ptDg_')[0]
    m = pd.read_csv(dir+file, index_col=0)
    flatM = np.array(m).flatten()
    logFlatM = np.log(flatM)
    # logFlatM[logFlatM==-np.inf]=-99999
    mats[sample]=logFlatM

np.corrcoef(list(mats.values()))
np.std(list(mats.values()))

bigArray = np.array(list(mats.values()))
sd = np.std(bigArray, axis=0)

bigArray[bigArray==-np.inf] = nan

avg = np.nanmean(bigArray, axis=0).reshape((4096,4096))
tab = pd.read_csv(files[0], index_col=0)
pd.DataFrame(avg, columns=tab.columns, index=tab.index).to_csv('predictedDg.csv')
