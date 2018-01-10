import itertools
import sys
import os
import pandas as pd
import collections

pileup='/hpc/hub_oudenaarden/edann/hexamers/mismatch/L1_R1_bismark_bt2_pe.sorted.filt.pileup'

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

def make_mm_table(listDic):
    pos = [int(i['pos']) for i in listDic]
    occ = dict.fromkeys(list(range(1,max(pos)+1)))
    for pos in occ.keys():
        occ[pos]=collections.Counter()
    for line in a:
        occ[int(line['pos'])][line['type']]+=1
    mmOcc = pd.DataFrame(occ)
    return(mmOcc)
