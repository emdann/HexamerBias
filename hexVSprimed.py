import pandas as pd
import pysam as ps
import collections
import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import argparse
import math

def mismatchKmers(params):
    seq,fqHex = params
    a = pairwise2.align.globalxs(fqHex, seq, -1000, -1000) # alignment funciton to not open gaps
    score = a[0][2]
    mmatch=6-score
    return(mmatch)

def make_occurrencies_tbl(seqDict):
    dic_oc={}
    for seqs in seqDict.values():
        primed,hex = seqs[0],seqs[1]
        if hex not in dic_oc.keys():
            dic_oc[hex]=[]
        dic_oc[hex].append(primed)
    hexs=list(dic_oc.keys())
    countOcc={}
    counts=[collections.Counter(i) for i in dic_oc.values()]
    for i in list(range(len(hexs))):
        countOcc[hexs[i]]=counts[i]
    oc_tbl=pd.DataFrame(countOcc)
    oc_tbl=oc_tbl.fillna(0)
    return(oc_tbl)

def makePWM(occTbl)

def countBases(hex):
    count={0:collections.Counter(), 1:collections.Counter(), 2:collections.Counter(), 3:collections.Counter(), 4:collections.Counter(), 5:collections.Counter()}
    for i in range(len(hex)):
        for pos in count.keys():
            count[pos][hex.index[i][pos]]+=hex[i]
    df = pd.DataFrame(count)
    return(df/df.sum())

def computeEntropy(pwm):
    H=0
    for pos in pwm.iteritems():
        h=0
        for base in ["A","T","C", "G"]:
            if pos[1][base]!=0:
                h+= - pos[1][base]*math.log(pos[1][base], 4)
        H+=h
    return(H)

# filename='L1_R2_hex_entropy.txt'
# with open(filename, 'w') as f:
#     for row in tbl.iterrows():
#         print(row[0], computeEntropy(countBases(row[1])), sep='\t', file=f)
