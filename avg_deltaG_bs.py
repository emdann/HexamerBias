import pysam as ps
import os
import pandas as pd
import numpy as np
from hexVSprimed import *

files = [f for f in os.listdir() if 'ptDg' in f]

dgMat = pd.read_csv(files[0], index_col=0, compression=findCompr(files[0]))
for f in files[1:]:
    df = pd.read_csv(f, index_col=0, compression=findCompr(f))
    dgMat = dgMat + df

avgMat = dgMat/len(files)

longDg=np.array(mat).reshape((16777216))
longDg2=np.array(mat2).reshape((16777216))
longDg3=np.array(mat3).reshape((16777216))

np.std((np.array(mat), np.array(mat2), np.array(mat3)), axis=0, ddof=1)

dgArrays = []
for f in files:
    df = pd.read_csv(f, index_col=0)
    array = np.ma.log(np.array(df)).reshape((16777216))
    dgArrays.append(array)

avgMatrix = np.mean(dgArrays, axis=0)
sdMatrix = np.std(dgArrays, axis=0, ddof=1)

s=0
for a in dgArrays:
    s+=a[-1]


print(s/len(dgArrays))

avg = np.array([])
for i in range(len(dgArrays[0])):
    values = [a[i] for a in dgArrays]
    mean = np.nanmean(values)
    avg = np.append(avg,mean)


for f in files:
    df = pd.read_csv(f, index_col=0)
    array = np.ma.log(np.array(df))
    print(array[0][0])
