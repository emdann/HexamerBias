import sys, os
from pandas.io.parsers import read_csv
import pandas as pd

try:
    df = read_csv(sys.argv[1], sep = '\t', header = None)
    output = sys.argv[2]
except:
    sys.exit('Please, give (1) input cov2c file; (2) root for output file')


chrlist = set(df[0])

for chromo in chrlist:
    tmpdf = df[df[0] == chromo]
    #tmpdf = tmpdf[tmpdf[[3,4]].sum(axis=1)>0]
    tmpdf.to_csv('.'.join([output, chromo]), sep = '\t', index = None, header = None)



