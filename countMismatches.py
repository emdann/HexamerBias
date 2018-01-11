import itertools
import sys
import os
import pandas as pd
import collections

pileup='/hpc/hub_oudenaarden/edann/hexamers/mismatch/L1_R1_bismark_bt2_pe.sorted.filt.pileup'

def read_pileup(pileupfile):
    a = []
    with open(pileup, 'rb') as f:
        for line in f.readlines():
            line = line.strip().decode('utf8').split("\t")
            type = line[2].upper()+'->'+line[4].upper()
            pos = line[6]
            if line[4].isupper():
                strand='+'
            else:
                strand='-'
            a.append({'type':type, 'pos':pos, 'strand':strand})
    return(a)

def strandSpecificDic(a, strand):
    df = pd.DataFrame(a)
    dfStrand = df[df.strand==strand]
    d = dfStrand.T.to_dict()
    return(list(d.values()))


def make_mm_table(listDic):
    pos = [int(i['pos']) for i in listDic]
    occ = dict.fromkeys(list(range(1,max(pos)+1)))
    for pos in occ.keys():
        occ[pos]=collections.Counter()
    for line in listDic:
        occ[int(line['pos'])][line['type']]+=1
    mmOcc = pd.DataFrame(occ)
    return(mmOcc)

a = read_pileup(pileup)
tab = make_mm_table(a)
tab.to_csv('/hpc/hub_oudenaarden/edann/hexamers/mismatch/L1_mismatchfreq.csv')
tabPlus = make_mm_table(strandSpecificDic(a, '+'))
tabMinus = make_mm_table(strandSpecificDic(a, '-'))
