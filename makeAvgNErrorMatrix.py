import fnmatch
import os
import pandas as pd
import numpy as np
import argparse
from hexVSprimed import *

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get hexamers used in fasta file.\n By Emma Dann")
argparser.add_argument('input', type=str, help='File of predicted Dg')
argparser.add_argument('cellab', type=str, help='name of output file')
argparser.add_argument('output', type=str, help='suffix for output files')
args = argparser.parse_args()

def makePredictedDgMatrix(file, cellAb):
    '''
    Turns file of pt - predicted Dg (precessed in R) to matrix of template on row and primer on column
    Needs cell abundance file to fill in missing pt pairs
    '''
    df = pd.read_csv(file, index_col=0)
    tabDic = {}
    sdDic = {}
    for row in df.iterrows():
        pt,dg,sd = row[1].pt, row[1].predictedDg, row[1].sd
        t,p = pt.split('-')
        if t not in tabDic.keys():
            tabDic[t]={}
        tabDic[t][p] = dg
        if t not in sdDic.keys():
            sdDic[t]={}
        sdDic[t][p] = sd
    ptMat = pd.DataFrame(tabDic)
    ptMat = fillNsortPTmatrix(ptMat, cellAb)
    sdMat = pd.DataFrame(sdDic)
    sdMat = fillNsortPTmatrix(sdMat, cellAb)
    return((ptMat,sdMat))

tabAb = pd.read_csv(args.cellab, index_col=0, compression = findCompr(args.cellab))
cellAb = tabAb[1]
cellAb = cellAb[[i for i in cellAb.index if 'N' not in i]]

predictedDgFile = args.input
# Make predicted Dg tab
dgMat,errMat = makePredictedDgMatrix(predictedDgFile, cellAb)
dgMat.to_csv(args.output+'dgMat.csv')
errMat.to_csv(args.output+'errMat.csv')
